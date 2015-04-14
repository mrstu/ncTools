#!/usr/bin/env python

'''
netCDF.py methods.

@author: Matt Stumbaugh
'''

import os
import sys
import getopt

import glob
import string
# from datetime import datetime, timedelta
from collections import OrderedDict as OD

import numpy as np
import numpy.ma as ma
import pandas as pd

import argparse

try:
    from netCDF4 import Dataset
    from netCDF4 import num2date, date2num
except:
    print 'No netCDF4 module in path.  Proceeding without netCDF4 module support'

import PyNC

def getargs(*args):
    '''Parse command line arguments

    Current features:
    
    * Input/output directory (optional)
    * Input file names (required, some pattern matching capability)
    * Output filename (optional, but recommended)
    * Modes for selecting variables: (ULM), vic, vicQ, vicNoQ
    '''
    
    usage = "usage: %prog [options] arg"
    parser = argparse.ArgumentParser(description='Summarize a VIC flux file')

    parser.add_argument("-dd", "--dir_data", dest="dirdata", default="./")
    parser.add_argument("-dn", "--data_name", dest="dataname", nargs="+", default=None, help='i.e. "PREC_??.nc"') 
    
    parser.add_argument("-o", "--output_directory", dest="dirout", default="./")
    parser.add_argument("-p", "--output_prefix", dest="nameout", default=None)
    
#     parser.add_argument("-vn", "--variable_name", dest="varname", default="Baseflow")    

    parser.add_argument("--vic", action="store_true", default=None)    
    parser.add_argument("--vicQ", action="store_true", default=None)        
    parser.add_argument("--vicNoQ", action="store_true", default=None)        
    
    parser.add_argument("-v", "--verbose", action="store_true", default=None)    

#    parser.add_option("-f", "--suffix",
#                      dest="suffix", nargs='+', default='*.c, *.scr, *.AWK')
#     parser.add_argument("-m", "--metadata", dest="header_file", default=None, help='Place holder for attribute info')
#     parser.add_argument("-t", "--time_data", dest="timedata", default=None, help='File containing timestamps')        
#     parser.add_argument("-s", "--tstart",
#                       dest="starttime", default=None, help='Start date (yyyymmdd)')
#     parser.add_argument("-e", "--tend",
#                       dest="endtime", default=None, help='End date (yyyymmdd)')
#     parser.add_argument("-tf", "--time_frequency",
#                   default="MS", help='http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases')

    args = parser.parse_args()

    return args

def glob_data(indir,name):
    '''Limited shell completion support for name search
    
    For arg -dn: these should work: 
        
        splitlist??
        splitlist[29,30]
        splitlist29 splitlist30
        splitlist[25-30]
        splitlist0[0-9] 'splitlist[1-9][0-9]'
    '''
    masterlist=[]
#     if len(name)>1:
    print 'Input directory:', indir
#     print name
    for nm in name:
#         print nm
        masterlist.extend(glob.glob(indir+'/'+nm))
#     else:
#         
#         sdfiles = glob.glob(indir+'/'+name)
    masterlist.sort()
    print 'Infiles:', masterlist[0], 'to', masterlist[-1]
    return masterlist

def selectVars(args):
    '''List of variable names to extract from input tile files and aggregate into merged output
        
    Defaults to select ULM variables
    '''

    if args.vic:
        # VIC all
        #dt_table='Evaporation Qs Qsb SoilMoist SoilTemp SWE Sensible Latent LatentSub Ground Rnet SurfTemp Shortwave Shortwave_net Longwave Longwave_net Precip Qair PET_SatSoil PET_h20Surf PET_Short PET_Tall PET_NatVeg Transp Qsm'.split()
        dt_table='Evaporation Qs Qsb SoilMoist SoilTemp SWE Sensible Latent LatentSub Ground Rnet SurfTemp Shortwave Shortwave_net Longwave Longwave_net Precip Qair PET_SatSoil PET_h20Surf PET_Short PET_Tall PET_NatVeg Transp Qsm'.split()
    elif args.vicQ:
        # VIC runoff/baseflow only
        dt_table='Qs Qsb'.split()        
    elif args.vicNoQ:
        # VIC ex-runoff/baseflow
        dt_table='AvgSurfT Evap LWnet PotEvap Qg Qh Qle Qsm Rainf RootMoist SMLiqFrac Snowf SoilMoist SoilTemp SWE SWEVeg SWnet TVeg'.split()
    else:
        #ULM default
        dt_table='AvgSurfT Evap LWnet PotEvap Qg Qh Qle Qs Qsb Qsm Rainf SMLiqFrac Snowf SoilMoist SoilTemp SWE SWEVeg SWnet TVeg'.split()
    
    return dt_table
 
def main(*args):   
    '''Merge spatially tiled netcdf files.  
    
    Supports the following:
    
    * 1/16th degree cell size
    * ULM compressed monthly output
    * vic2nc.c compressed monthly output
    ''' 
     
    try:
        
        args = getargs(*args)
        
        '''Identify variables to merge'''
        
        dt_table = selectVars(args)
        
        '''Build list of input file paths'''
        
        datas = glob_data(args.dirdata, args.dataname)    
        
        '''Make output directory'''
        
        try:
            os.makedirs(args.dirout)
        except OSError:
            pass
        
        '''Read tiled NetCDF inputs'''
        
        ncz = OD()
        ncf = OD()
        nindex=0
        
        for ndat, dat in enumerate(datas):
    
            ncz[ndat] = PyNC.PyNC()
            ncf[ndat] = dat
            
            ncz[ndat].read_ncblock(dat)
            ncz[ndat].read_dims()
            ncz[ndat].read_numvars()
            ncz[ndat].read_scalars()
            
            nindex=nindex+len(ncz[ndat].var1d['CellID'])
            ntimes=ncz[ndat].dims['tstep']
            try:
                nlevels=ncz[ndat].dims['level']
            except KeyError:
                nlevels=ncz[ndat].dims['z']

        '''Determine unique latitude and longitude dimensions'''
                
        lats=np.zeros((nindex,),np.float)
        lons=np.zeros((nindex,),np.float)    
        offset=0
        
        for n in xrange(ndat+1):
            inds = ncz[n].var1d['CellID'].astype(int) - 1 + offset
            lats[inds] = ncz[n].var1d['nav_lat']
            lons[inds] = ncz[n].var1d['nav_lon']
            offset+=len(inds)

        ulats = np.unique(lats)
        ulons = np.unique(lons)

        ''' get xy space indices'''
        
        xres=0.0625 #ulons[1]-ulons[0]
        xll=np.min(ulons)
        yll=np.min(ulats)
        dx = ulons-xll
        dy = ulats-yll

        nx = int(np.max(dx)/xres) + 1
        ny = int(np.max(dy)/xres) + 1
        if nx != len(ulons) or ny != len(ulats):
            ''' need to add to ulats and/or ulons'''
            ulons=np.arange(xll,np.max(ulons)+xres,xres)
            ulats=np.arange(yll,np.max(ulats)+xres,xres)
            nx=len(ulons)
            ny=len(ulats)
        
        '''Create output nc'''
        
        txy = np.ones((ntimes,ny,nx),np.float) * -9999.
        
        innc=ncz[ndat].nc
     
        outnc = PyNC.PyNC()
        if not args.nameout:
            nameout=ncz[ndat].name2d[0]
        else:
            nameout=args.nameout
        outfile=os.path.join(args.dirout,nameout)
        outnc.create_NC(outfile)    
        PyNC.trans_glob_atts(innc, outnc.nc) #, exclude_atts=excludes)
    
        outnc.nc.setncattr('merge_path',outfile)
   
        if innc.dimensions.has_key('tstep'):
            # daily vic2nc.c output
            tsteps = innc.dimensions['tstep']            
            times = innc.variables['time']
            timecal=times.calendar

        units = times.units
    ##         units = times.units.replace('sec','seconds')
        try:
            units = 'seconds since %s'%times.time_origin
        except:
            units = 'seconds since %s'%times.origin
        print units

        dates = num2date(times[:], units=units, calendar=timecal)
        
        times_origin=dates[0].strftime("%Y-%m-%d")
        times_units = 'days since %s'%(times_origin)
        times_days = date2num(dates,units=times_units,calendar=timecal)
            
        atts=OD([('units',times_units), 
                 ('origin',times_origin), 
                 ('long_name','Time axis')])
        outnc.create_coord_var(outnc.nc, times_days, name='Time', dtype='i4', **atts)
        outnc.nc.createDimension('Level',nlevels) 
    #     outnc.nc.createVariable('Time','f8',('Time',)) 
        atts=OD([('units','degrees_east'), 
                 ('long_name','longitude'), 
                 ('valid_min',-180.),
                 ('valid_max',180.)])
        outnc.create_coord_var(outnc.nc, ulons, name='Longitude', dtype='f4', **atts)
        atts=OD([('units','degrees_north'), 
                 ('long_name','latitude'), 
                 ('valid_min',-90.),
                 ('valid_max',90.)])
        outnc.create_coord_var(outnc.nc, ulats, name='Latitude', dtype='f4', **atts)    

        '''Write to output file'''
    
        for name in dt_table:
            var3d=0
    
            if name in ncz[ndat].name3d:
                tzxy = np.ones((ntimes,nlevels,ny,nx),np.float) * -9999.
                var3d=1
                
            for inc, nc in ncz.iteritems():
                slats=nc.var1d['nav_lat']
                slons=nc.var1d['nav_lon']
                ax = (slons-xll)/xres 
                ay = (slats-yll)/xres 
                ix = ax.astype(int)
                iy = ay.astype(int)        
                selvar = nc.get_vector(name)
                if var3d:
                    tzxy[:,:,iy,ix]=selvar[:,:,:]# works for 3d?                
                else:
                    txy[:,iy,ix]=selvar[:,:]# works for 3d?
   
            if var3d:
                hov = outnc.nc.createVariable(name,'f4',('Time','Level','Latitude','Longitude',), zlib=True, fill_value=-9999.)
                hov.setncattr('units',innc.variables[name].__dict__['units'])
                mxy = ma.masked_values(tzxy, -9999.)
                hov[:,:,:,:] = mxy.data
                
            else:
                mxy = ma.masked_values(txy, -9999.)            
                hov = outnc.nc.createVariable(name,'f4',('Time','Latitude','Longitude',), zlib=True, fill_value=-9999.)
                hov.setncattr('units',innc.variables[name].__dict__['units'])            
                hov[:,:,:] = mxy.data
        return 0
    
    except Exception as exception:
        
        '''Generic handling of tracebacks'''
        
        print "Exception name: ",type(exception).__name__
        print exception.error
        return 1
        
if __name__=='__main__':   
    sys.exit(main())
    

#python -m mattslib.nclib.ncblk2ncgrd -dd baseflow_aug_01 -dn mean_BASEFLOW_??.nc -vn Baseflow -p test00
#python -m mattslib.nclib.ncblk2ncgrd -dd baseflow_aug_01 -dn mean_BASEFLOW_??.nc -vn Baseflow -p test00
#qsub1L.bsh 'python -m mattslib.nclib.ncblk2ncgrd -dd nc/ -dn PREC_??.nc -vn Precipitation -p test00' plotPrec.qsub
#qsub1L.bsh 'python -m mattslib.nclib.ncblk2ncgrd -dd nc/ -dn SWE_??.nc -vn SWE -p test00' plotSWE.qsub

    
    