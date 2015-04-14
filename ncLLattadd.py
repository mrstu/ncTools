#!/usr/bin/env python

'''
netCDF.py methods.

@author: Matt Stumbaugh
'''

import os
import sys
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

#from mattslib.nclib import PyNC
import PyNC
from configobj import ConfigObj

# OD([('units','degrees_east'), ('long_name','longitude'), ('valid_min',-180.), ('valid_max',180.)]),
# OD([('units','degrees_north'), ('long_name','latitude'), ('valid_min',-90.), ('valid_max',90.)]),
# OD([('units','mm/s'), ('long_name','Evaporation'), ('valid_min',-90.), ('valid_max',90.)]),
# OD([('units','degrees_north'), ('long_name','latitude'), ('valid_min',-90.), ('valid_max',90.)]),
# OD([('units','degrees_north'), ('long_name','latitude'), ('valid_min',-90.), ('valid_max',90.)]),
#         atts=OD([('units',nc.times_units), 
#                  ('origin',nc.times_origin), 
#                  ('long_name','Time axis')])
#         nc.create_coord_var(nc.nc, nc.times_days, name='Time', dtype='f8', **atts)
#      
#         atts=OD([('units','spatial_index'), 
#                  ('long_name','grid index'),
#                  ('valid_min',0)])


def getargs(*args):
    '''Parse command line arguments'''
    
    usage = "usage: %prog [options] arg"
    parser = argparse.ArgumentParser(description='Summarize a VIC flux file')
#    parser.add_option("-f", "--suffix",
#                      dest="suffix", nargs='+', default='*.c, *.scr, *.AWK')
    parser.add_argument("-dd", "--dir_data", dest="dirdata", default="./")
    parser.add_argument("-dn", "--data_name", dest="dataname", nargs="+", default=None, help='i.e. "PREC_??.nc"') 

    parser.add_argument("-o", "--output_directory", dest="dirout", default="./")
    parser.add_argument("-p", "--output_prefix", dest="nameout", default=None)

#     parser.add_argument("-vn", "--variable_name", dest="varname", default="Baseflow")    

#     parser.add_argument("-cfg", "--config_file", dest="configfile", default=None)
    parser.add_argument("--config_globals", dest="globalcfg", default=None)
    parser.add_argument("--config_vars", dest="varscfg", default=None)        
    
    parser.add_argument("-tc", "--num_time_chunks", dest="timechunks", default=-1, type=int, help='Write out data in chunks (totals steps/timechunks).') 

    parser.add_argument("--vic",
                      action="store_true", default=None)    

    parser.add_argument("--vicQ",
                      action="store_true", default=None)        

    parser.add_argument("--vicNoQ",
                      action="store_true", default=None)        

    parser.add_argument("-v", "--verbose",
                      action="store_true", default=None)    

    args = parser.parse_args()
    # dirdata, dataname, infmt
    # dirspatial, spatialname
    # dirout, nameout
    # configfile, varname
    # timedata

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
 
def main(*args):   
    '''Read a netcdf WRF file and create map for nc variables (T2,XLAT,XLON).
    
    ./readfluxnc2.py $outpath $inputs
    
    outpath is directoy where output NetCDF will be stored
    inputs is directory where input NetCDF will be read from.
    '''  
    args = getargs(*args)

#      dt_table='\n'.join('AvgSurfT Evap LWnet PotEvap Qfz Qg Qh Qle Qs Qsb Qsm Qst Rainf RootMoist SatSoil SMFrozFrac SMLiqFrac Snowf SOILDEPTH SoilMoist SoilTemp SWE SWEVeg SWnet TVeg WltSoil'.split())        
    #dt_table='AvgSurfT DelIntercept DelSoilMoist DelSurfStor DelSWE ECanop ESoil Evap LWnet PotEvap Qfz Qg Qh Qle Qs Qsb Qsm Qst Rainf RootMoist SatSoil SMFrozFrac SMLiqFrac Snowf SOILDEPTH SoilMoist SoilTemp SWE SWEVeg SWnet TVeg WltSoil'.split()
# VIC merge    
#    dt_table='Evaporation Qs Qsb SoilMoist SoilTemp SWE Sensible Latent LatentSub Ground Rnet SurfTemp Shortwave Shortwave_net Longwave Longwave_net Precip Qair PET_SatSoil PET_h20Surf PET_Short PET_Tall PET_NatVeg Transp Qsm'.split()
# ULM merge
     #dt_table='AvgSurfT Evap LWnet PotEvap Qfz Qg Qh Qle Qs Qsb Qsm Qst Rainf RootMoist SatSoil SMFrozFrac SMLiqFrac Snowf SOILDEPTH SoilMoist SoilTemp SWE SWEVeg SWnet TVeg WltSoil'.split()
    if args.vic:
        #dt_table='Evaporation Qs Qsb SoilMoist SoilTemp SWE Sensible Latent LatentSub Ground Rnet SurfTemp Shortwave Shortwave_net Longwave Longwave_net Precip Qair PET_SatSoil PET_h20Surf PET_Short PET_Tall PET_NatVeg Transp Qsm'.split()
        #dt_table='Evaporation Qs Qsb SoilMoist SoilTemp SWE Sensible Latent LatentSub Ground Rnet SurfTemp Shortwave Shortwave_net Longwave Longwave_net Precip Qair PET_SatSoil PET_h2OSurf PET_Short PET_Tall PET_NatVeg Transp Qsm TotalSoilMoist TotalRunoff'.split()
        dt_table='Evaporation Qs Qsb SoilMoist SoilTemp SWE Sensible Latent LatentSub Ground Rnet SurfTemp Shortwave ShortwaveNet Longwave LongwaveNet Precip Qair petSatSoil petH2OSurf petShort petTall petNatVeg Transp Qsm TotalSoilMoist TotalRunoff'.split()
    elif args.vicQ:
        dt_table='Qs Qsb'.split()        
    elif args.vicNoQ:
        dt_table='AvgSurfT Evap LWnet PotEvap Qg Qh Qle Qsm Rainf RootMoist SMLiqFrac Snowf SoilMoist SoilTemp SWE SWEVeg SWnet TVeg'.split()
    else:
        dt_table='AvgSurfT Evap LWnet PotEvap Qg Qh Qle Qs Qsb Qsm Rainf RootMoist SMLiqFrac Snowf SoilMoist SoilTemp SWE SWEVeg SWnet TVeg'.split()
#     dt_table='AvgSurfT Evap LWnet PotEvap Qfz Qg Qh Qle Qs Qsb Qsm Qst Rainf RootMoist SatSoil SMFrozFrac SMLiqFrac Snowf SOILDEPTH SoilMoist SoilTemp SWE SWEVeg SWnet TVeg WltSoil'.split()
#    dt_table=['SOILDEPTH']
#    dt_table='AvgSurfT Evap Qs Qsb'.split()


    datas = glob_data(args.dirdata, args.dataname)    
    try:
        os.makedirs(args.dirout)
    except OSError:
        pass
    
    ncz = OD()
    ncf = OD()
    nindex=0
    
#     for ndat, dat in enumerate(datas):
# 
#         #print 'ndat, dat', ndat, dat
#         ncz[ndat] = PyNC.PyNC()
#         ncf[ndat] = dat
#         
#         ncz[ndat].read_ncblock(dat)
#         ncz[ndat].read_dims()
#         ncz[ndat].read_numvars()
#         ncz[ndat].read_scalars()
# #         ncz[ndat].read_vector(args.varname)
# #         ncz[ndat].read_vectors()    
#     innc=ncz[ndat].nc

    print datas[0]
    ncdat=PyNC.PyNC() 
    ncdat.read_ncblock(datas[0])
    ncdat.read_dims()
    ncdat.read_numvars()
    ncdat.read_scalars()
    innc=ncdat.nc
    
    outnc = PyNC.PyNC()
    if not args.nameout:
        nameout=innc.name2d[0]
    else:
        nameout=args.nameout
    outfile=os.path.join(args.dirout,nameout)
    outnc.create_NC(outfile)    
    
    if args.globalcfg:
        print 'Reading configuration:', args.globalcfg
        cfg=ConfigObj(args.globalcfg)
        keys = cfg['GLOBAL_ATTRIBUTES'].keys()
        outnc.add_atts_global(*keys, **cfg['GLOBAL_ATTRIBUTES'])
    if args.varscfg:
        print 'Reading configuration:', args.varscfg
        cfgvars=ConfigObj(args.varscfg)

#     if args.configfile:
#         print 'Reading configuration:', args.configfile
#         options, global_atts, domain_dict, fields = v2nc.parse_config(args.configfile,**{'colskip':True})
#         print global_atts
#         keys=global_atts.keys()
#         outnc.add_atts_global(*keys, **global_atts)
    else:
        PyNC.trans_glob_atts(innc, outnc.nc) #, exclude_atts=excludes)
    
    #PyNC.trans_glob_atts(innc, outnc.nc, exclude_atts="lev")
    #outnc.nc.setncattr('merge_path',outfile)
 
#     outnc.nc.createDimension('Time',tsteps) 
#     outnc.nc.createVariable('Time','f8',('Time',))
#     outnc.nc.variables['Time'][:] = innc.variables['time'][:]    
#     for att,attv in innc.variables['time'].__dict__.iteritems():
#         outnc.nc.variables['Time'].setncattr(att,attv)

    print innc.dimensions
    if innc.dimensions.has_key('Time'):
        # daily vic2nc.c output
        tsteps = innc.dimensions['Time']            
        times = innc.variables['Time']
        timecal=times.calendar

    dates = num2date(times, units=times.units, calendar=timecal)
#         print [date.strftime("%Y-%m-%d") for date in dates]
#        times_origin=dates[0].strftime("%Y-%m-%d %H:%M:%S")
#        times_units = 'seconds since %s'%(times_origin)
    times_origin=dates[0].strftime("%Y-%m-%d")
    times_units = 'days since %s'%(times_origin)
    times_days = date2num(dates,units=times_units,calendar=timecal)
        
    atts=OD([('units',times_units), 
             ('origin',times_origin), 
             ('long_name','Time axis')])
    outnc.create_coord_var(outnc.nc, times_days, name='Time', dtype='i4', **atts)
    
    ## For some reason 
#     if innc.dims.has_key('lev'):
#     print ncdat.dims
#     outnc.nc.createDimension('Level',ncdat.dims['lev']) 
    try:
        if ncdat.dims.has_key('lev'):
            outnc.nc.createDimension('Level',ncdat.dims['lev']) 
    except:
        pass

#     try:
#         if innc.dims.has_key('lev'):
#             outnc.nc.createDimension('Level',innc.dims['lev']) 
#     except:
#         pass
        
#     outnc.nc.createVariable('Time','f8',('Time',)) 
    atts=OD([('units','degrees_east'), 
             ('long_name','longitude'), 
             ('valid_min',-180.),
             ('valid_max',180.)])
    outnc.create_coord_var(outnc.nc, innc.variables['Longitude'], name='Longitude', dtype='f8', **atts)
 
    atts=OD([('units','degrees_north'), 
             ('long_name','latitude'), 
             ('valid_min',-90.),
             ('valid_max',90.)])
    outnc.create_coord_var(outnc.nc, innc.variables['Latitude'], name='Latitude', dtype='f8', **atts)    

    timeinds=np.arange(0,len(times_days))

    if args.timechunks !=-1:    
        if len(times_days) > args.timechunks:
            splitinds=np.array_split(timeinds,args.timechunks) #chk was 20 previously for 60 splits
        else:
            splitinds=[timeinds]


    outnc.create_coord_var(outnc.nc, [1], name='crs', dtype='f4', **cfgvars['crs'])

    varavail=[str(key) for key in ncdat.nc.variables.keys()]
#     print varavail
    
    for name in dt_table:
        var4d=0

        if name in ncdat.name4d:
# #         if name in ncz[ndat].name3d:
# #             tzxy = np.ones((tsteps,nlevels,ny,nx),np.float) * -9999.
#             txyz = np.ones((ncz[ndat].dims['Time'],ncz[ndat].dims['lev'], ncz[ndat].dims['Latitude'], ncz[ndat].dims['Longitude']),np.float) * -9999.    
            var4d=1
#         else:
#             txy = np.ones((ncz[ndat].dims['Time'], ncz[ndat].dims['Latitude'], ncz[ndat].dims['Longitude']),np.float) * -9999.    
       
        if name in varavail:
            print 'name to add', name
            
            if var4d:
                print '4d'
                hov = outnc.nc.createVariable(name,'f4',('Time','Level','Latitude','Longitude',), zlib=True, fill_value=-9999.)
            else:
                print '3d'                
                hov = outnc.nc.createVariable(name,'f4',('Time','Latitude','Longitude',), zlib=True, fill_value=-9999.)

#                 hov = outnc.nc.createVariable(name,'f4',('Time','Latitude','Longitude',), zlib=True, fill_value=-9999.)
            keys = cfgvars[name].keys()
            outnc.add_atts_var(outnc.nc.variables[name], *keys, **cfgvars[name])
#             print cfgvars[name]

# #                 txy[:,:,:]=selvar[:,:,:]# works for 3d?
# #                 mxy = ma.masked_values(txy, -9999.)
# #                 hov[:,:,:] = mxy.data
# 
# #                 txy[:,:,:]=selvar[:,:,:]# works for 3d?
#             mxy = ma.masked_values(selvar[:,:,:], -9999.)
#             hov[:,:,:] = mxy.data

            if args.timechunks ==-1:  
#                 print 'No chunk', name          
                selvar = ncdat.get_vector(name)
                mxy = ma.masked_values(selvar[:,:,:], -9999.)
#                 print mxy.shape()
                hov[:,:,:] = mxy.data
                
            else:
                # Chunked output
                for chk in splitinds:
                    selvar = ncdat.get_vector(name, chunk=chk)
                    mxy = ma.masked_values(selvar[:,:,:], -9999.)                    
                    hov[chk,:,:]=mxy.data   
                    
    #                 print 'tsteps', tsteps
    #                 print selvar.shape
    #                 txy[:,iy,ix]=selvar[:,:]# works for 2d
#                     hov[chk,:,:]=txy   



                

#                 mxy = ma.masked_values(tzxy, -9999.)
#                 hov[:,:,:,:] = mxy.data


# #                 txy = np.ones((chk.size,ny,nx),np.float) * -9999.                                
# #                 tzxy[:,:,iy,ix]=selvar[:,:,:]# works for 3d?# 
# #                 mxy = ma.masked_values(txy, -9999.)                
# #                 hov[:,:,:] = mxy.data
# #                 hov = nc.nc.variables[name][:]
#                 print 'selvar', selvar
#                 hov = selvar[:,:,:]
#                 
#                 
#             selvar=nc.nc.variables[name]
#             if var3d:
#                 tzxy[:,:,iy,ix]=selvar[:,:,:]# works for 3d?
#             else:
#                 txy[:,iy,ix]=selvar[:,:]# works for 3d?
# 
#         if var3d:
#             hov = outnc.nc.createVariable(name,'f4',('Time','Level','Latitude','Longitude',), zlib=True, fill_value=-9999.)
#             keys = cfg[name].keys()
#             outnc.add_atts_var(outnc.nc.__dict__[name], *keys, **cfg[name])            
#             mxy = ma.masked_values(tzxy, -9999.)
#             hov[:,:,:,:] = mxy.data
#         else:
#             hov = outnc.nc.createVariable(name,'f4',('Time','Latitude','Longitude',), zlib=True, fill_value=-9999.)
#             keys = cfg[name].keys()
#             outnc.add_atts_var(outnc.nc.__dict__[name], *keys, **cfg[name])
#             mxy = ma.masked_values(txy, -9999.)
#             hov[:,:,:] = mxy.data

            ## Chunked output
#             for chk in splitinds:
#                 txy = np.ones((chk.size,ny,nx),np.float) * -9999.
#             
#                 for inc, nc in ncz.iteritems():
#                     slats=nc.var1d['Latitude']
#                     slons=nc.var1d['Longitude']
#                     ax = (slons-xll)/xres 
#                     ay = (slats-yll)/xres 
#                     ix = ax.astype(int)
#                     iy = ay.astype(int)        
#                     selvar = nc.get_vector(nc.name2d[0], chunk=chk)
#                     print 'tsteps', tsteps
#                     print selvar.shape
#                     txy[:,iy,ix]=selvar[:,:]# works for 2d
#                 hov[chk,:,:]=txy   

            
    #         outnc.nc.variables[name].setncattr('units','something')
    #         print innc.varg
        
if __name__=='__main__':   
    main()

#python -m mattslib.nclib.ncblk2ncgrd -dd baseflow_aug_01 -dn mean_BASEFLOW_??.nc -vn Baseflow -p test00
#python -m mattslib.nclib.ncblk2ncgrd -dd baseflow_aug_01 -dn mean_BASEFLOW_??.nc -vn Baseflow -p test00
#qsub1L.bsh 'python -m mattslib.nclib.ncblk2ncgrd -dd nc/ -dn PREC_??.nc -vn Precipitation -p test00' plotPrec.qsub
#qsub1L.bsh 'python -m mattslib.nclib.ncblk2ncgrd -dd nc/ -dn SWE_??.nc -vn SWE -p test00' plotSWE.qsub

    
    