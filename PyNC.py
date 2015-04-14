#!/usr/bin/env python

import os 
import os.path
import sys

import string

import numpy as np
import numpy.ma as ma

from collections import OrderedDict as OD

try:
    from netCDF4 import Dataset
    from netCDF4 import num2date, date2num
except:
    print 'No netCDF4 module in path.  Proceeding without netCDF4 module support'

def trans_glob_atts(parentgrp, childgrp, exclude_atts=[]):
    for att, val in parentgrp.__dict__.iteritems():
        if att in exclude_atts:
            pass
        else:
            childgrp.setncattr(att,val)
            
def trans_coord_var(parentgrp, childgrp, mapnames):
    pass

class PyNC:

    def __init__(self):
        self.nctype='NETCDF4'

    def create_NC(self,name):
        
        self.nc = Dataset(name, 'w', clobber=True, format=self.nctype) 

    def add_atts_global(self, *args, **kwargs):
        if args:
            for attribute in args:
                value = kwargs[attribute]

                if attribute in 'geospatial_lat_max geospatial_lat_min geospatial_lon_max geospatial_lon_min geospatial_lat_resolution geospatial_lon_resolution geospatial_vertical_min geospatial_vertical_max geospatial_vertical_resolution'.split():
                    value = float(kwargs[attribute])
                else:
                    value = kwargs[attribute]
                
#             for attribute, value in kwargs.iteritems():
                if hasattr(self.nc, attribute):
                    print('WARNING: Attribute {0} already \
                          exists'.format(attribute))
                    print('Renaming to g_{0} to ovoid \
                          overwriting.'.format(attribute))
                    attribute = 'g_{0}'.format(attribute)
                setattr(self.nc, attribute, value)
        else:
            for attribute, value in kwargs.iteritems():
                if hasattr(self.nc, attribute):
                    print('WARNING: Attribute {0} already \
                          exists'.format(attribute))
                    print('Renaming to g_{0} to ovoid \
                          overwriting.'.format(attribute))
                    attribute = 'g_{0}'.format(attribute)
                setattr(self.nc, attribute, value)

    def add_atts_var(self, target, *args, **kwargs):

        if args:
            for attribute in args:
                value = kwargs[attribute]

                if attribute in 'geospatial_lat_max geospatial_lat_min geospatial_lon_max geospatial_lon_min geospatial_lat_resolution geospatial_lon_resolution geospatial_vertical_min geospatial_vertical_max geospatial_vertical_resolution'.split():
                    value = float(kwargs[attribute])
                else:
                    value = kwargs[attribute]
                
#             for attribute, value in kwargs.iteritems():
                if hasattr(target, attribute):
                    print('WARNING: Attribute {0} already \
                          exists'.format(attribute))
                    print('Renaming to g_{0} to ovoid \
                          overwriting.'.format(attribute))
                    attribute = 'g_{0}'.format(attribute)
                setattr(target, attribute, value)
        else:
            for attribute, value in kwargs.iteritems():
                if hasattr(target, attribute):
                    print('WARNING: Attribute {0} already \
                          exists'.format(attribute))
                    print('Renaming to g_{0} to ovoid \
                          overwriting.'.format(attribute))
                    attribute = 'g_{0}'.format(attribute)
                setattr(target, attribute, value)
                

    def set_datenums(self, dates, timestep="days", time_calendar="standard"):
        self.times_origin=dates[0].strftime("%Y-%m-%d 0:0:0")
        self.times_units = '%s since %s'%(timestep, self.times_origin)        
        self.times_days = date2num(dates,units=self.times_units,calendar=time_calendar)

    def create_coord_var(self, rgrp, vardata, name='TIME', dtype='i4', **atts):
        
        dim = rgrp.createDimension(name,len(vardata))
        cvdim = rgrp.createVariable(name,dtype,(name,))
        if atts:
            for att, val in atts.iteritems():
                if att in 'inverse_flattening longitude_of_prime_meridian semi_major_axis'.split():
                    cvdim.setncattr(att,float(val))
                else:
                    cvdim.setncattr(att,val)
        cvdim[:] = vardata

    def add_coord_dims(self,LON,LAT,**kwargs):
    
        atts=OD([('units','degrees_east'), 
                 ('long_name','longitude'), 
                 ('valid_min',-180.),
                 ('valid_max',180.)])
        self.create_coord_var(self.nc, LON, name='LON', dtype='f4', **atts)

        atts=OD([('units','degrees_north'), 
                 ('long_name','latitude'), 
                 ('valid_min',-90.),
                 ('valid_max',90.)])
        self.create_coord_var(self.nc, LAT, name='LAT', dtype='f4', **atts)
        
        if kwargs:
            print 'Only 2D (LON,LAT) coordinate dimensions implemented!'
            raise 

    def add_var2D(self, valm, name='VARNAME', coordDims=('LAT','LON',), dtype='f', compress=True, **kwargs):

#         vmask = self.nc.createVariable('MASK','i',('LAT','LON',), zlib=True)
        ncvar = self.nc.createVariable(name, dtype, coordDims, zlib=compress)
        if kwargs and kwargs.has_key('attrs'):
            attrs = kwargs['attrs']
            for attr in attrs:
                ncvar.setncattr(attr,attrs[attr])
#         ncvar.setncattr('_FillValue',valm.fill_value)                
#         vmask.setncattr('units','0=invalid, 1=valid')
#         vmask.setncattr('long_name','binary_land_mask')
        ncvar[:] = valm

    def add_var1D(self, valm, name='VARNAME', coordDims=('Index',), dtype='f', compress=True, **kwargs):

#         vmask = self.nc.createVariable('MASK','i',('LAT','LON',), zlib=True)
        ncvar = self.nc.createVariable(name, dtype, coordDims, zlib=compress)
        if kwargs and kwargs.has_key('attrs'):
            attrs = kwargs['attrs']
            for attr in attrs:
                ncvar.setncattr(attr,attrs[attr])
#         ncvar.setncattr('_FillValue',valm.fill_value)                
#         vmask.setncattr('units','0=invalid, 1=valid')
#         vmask.setncattr('long_name','binary_land_mask')
        ncvar[:] = valm


    def readnc(self, ncfilein):
        '''Read nc file'''
        
        print 'Reading', ncfilein
        
        self.rootgrp = Dataset(ncfilein,'r', format=self.nctype)

        print 'Listing dimenstions:'
        self.dimnames=[]
        self.dims={}
#         self.dimxyt={}
# #         xyt=('lons','lats','nlayer') # this order for params.nc
#         xyt=('lons','lats','time') # this order for CSC flux normals
# #         xyt=('time','lons','lats') # this order for monthly climo (livneh refined forcings)
# #         for ndim, dimobj in enumerate(self.rootgrp.dimensions.values()):
        print self.rootgrp.variables

        for ndim, dimobj in enumerate(self.rootgrp.dimensions):
            size = len(self.rootgrp.dimensions[str(dimobj)])
            name = str(dimobj)
            print name, size, len(dimobj)
#             self.__dict__[xyt[ndim]] = dimobj
#             print dimobj
            print ndim
            if self.rootgrp.variables.has_key(name):
                self.dimnames.append(name)
                self.dims[ndim] = self.rootgrp.variables[name]

#             if self.rootgrp.variables.has_key(name):
#                 self.__dict__[xyt[ndim]] = self.rootgrp.variables[name]
#             else:
#                 self.__dict__[xyt[ndim]] = range(size)
                
#         ## Reading variables        
        self.varlist2D = []
        self.varlist3D = []
        for varobj in self.rootgrp.variables:
            varname = varobj
            varndim = self.rootgrp.variables[varobj].ndim
            if varndim == 2: 
                self.varlist2D.append(varname)
            elif varndim == 3: 
                self.varlist3D.append(varname)
        self.varlist2D.sort()
        self.varlist3D.sort()                
#         print self.varlist

    def get_datetime(self,timevar='Time'):        
        ''' >>> dates = num2date(times[:],units=times.units,calendar=times.calendar)
        >>> print 'dates corresponding to time values:\n',dates'''
        times=self.rootgrp.variables[timevar]
        dates = num2date(times[:],
                         units=times.units,
                         calendar=times.calendar) 
        return dates


    def get_lonslatstimes(self):
        dimorder=[]
        for ndim, dim in enumerate(self.dimnames):
            print dim.lower()
            if dim.lower().startswith('lat'):
                self.lats=self.rootgrp.variables[dim][:]
            elif dim.lower().startswith('lon'):
                self.lons=self.rootgrp.variables[dim][:]                        
            elif dim.lower().startswith('time'):
                self.times=self.rootgrp.variables[dim][:]                        

    #                 dimpkorder.append(ndim)


    def get_dimorder(self,vname):
        dimorder=[]
        for dimlen in self.rootgrp.variables[vname].shape:
            # match 
            for ndim, dim in enumerate(self.dimnames):
                if len(self.dims[ndim]) == dimlen:
                    print ndim, dim, dimlen
                    print self.rootgrp.variables[vname]
                    dimorder.append(dim)
        return dimorder

    def get_limits(self,vname):
        limmax = np.amax(self.rootgrp.variables[vname][:])
        limmin = np.amin(self.rootgrp.variables[vname][:])
        return [limmin, limmax]
           
    def read_ncblock(self, ncfilein):
        '''Read nc file'''
        
        print 'Reading', ncfilein
        
        self.nc = Dataset(ncfilein,'r', format=self.nctype)

    def read_dims(self): 
        '''Reading dimensions'''
        self.dims = OD()
#         print 'Listing dimenstions:'
        for ndim, dimobj in enumerate(self.nc.dimensions):
            size = len(self.nc.dimensions[str(dimobj)])
            name = str(dimobj)

            self.dims[name]=size
#             print 'Dim:', name, size
            
    def read_scalars(self): 
        '''Reading one-dimensional variables'''
        self.var1d = OD()
        for varobj in self.nc.variables:
            varname = str(varobj)
            varndim = self.nc.variables[varobj].ndim
#             print# varndim
            if varndim == 1: 
                #print 'Var-1D:', varname, varndim                
                self.var1d[varname] = self.nc.variables[varobj][:]

                
    def read_numvars(self): 
        '''Reading multi-dimensional variables'''
        self.n2d = 0
        self.n3d = 0
        self.n4d = 0
        self.name2d=[]
        self.name3d=[]        
        self.name4d=[]                
        for varobj in self.nc.variables:
            varname = str(varobj)
            varndim = self.nc.variables[varobj].ndim
            if varndim == 2:
                self.n2d+=1
                self.name2d.append(varname)
            elif varndim == 3: 
                self.n3d+=1
                self.name3d.append(varname)
            elif varndim == 4: 
                self.n4d+=1
                self.name4d.append(varname)                
#             print 'Var-%dD:'%varndim, varname

    def read_vector(self, selvarname):
        '''Reading multi-dimensional variables'''
#         self.var2d = OD()
        print selvarname
        for varobj in self.nc.variables:
            varname = str(varobj)
            varndim = self.nc.variables[varobj].ndim
            if varname == selvarname:
                self.selvar = self.nc.variables[varobj][:]
#                 print 'Select 2D-variable:', varname, varndim                

    def read_vectors(self): 
        '''Reading multi-dimensional variables'''
        self.var2d = OD()
        self.var3d = OD()
        self.var4d = OD()
        for varobj in self.nc.variables:
            varname = str(varobj)
            varndim = self.nc.variables[varobj].ndim
            if varndim == 2: 
                self.var2d[varname] = self.nc.variables[varobj][:]
            elif varndim == 3: 
                self.var3d[varname] = self.nc.variables[varobj][:]
            elif varndim == 4: 
                self.var4d[varname] = self.nc.variables[varobj][:]

    def get_vector(self, selvarname, chunk=None):
        '''Reading multi-dimensional variables'''
#         self.var2d = OD()
#         print selvarname
        for varobj in self.nc.variables:
            varname = str(varobj)
            varndim = self.nc.variables[varobj].ndim
            if varname == selvarname:
#                 print 'Select 2D-variable:', varname, varndim
#                 print self.nc.variables[varobj].shape
                if chunk != None:
                    return self.nc.variables[varobj][chunk,:]
                else:
                    return self.nc.variables[varobj]

    def read_txtLLinds(self,LLi_file='goodlatlons.txt'):
        self.cellinfo=np.loadtxt(LLi_file, dtype={'names': ('lat','lon','ilat','jlon'), 'formats': ('f4','f4','i4','i4')})

    def dump_active_cells(self,outdir="junkdir"):
        for vname in self.varlist3D:
            dmns = self.get_dimorder(vname)
            outpath = os.path.join(outdir,vname)
            try:
                os.mkdir(outpath)
            except OSError:
                outdir, 'already exists!'

            print vname
            print dmns
            print self.dims[1]
            print self.dims[2]            

            for i1 in self.dims[1]:
                for i2 in self.dims[2]:
 
                    print i1, i2
#                     fname = "%s/data_%4.5f_%4.5f"%(outpath, self.cellinfo['lat'][ncell], self.cellinfo['lon'][ncell])
#                     np.savetxt(fname, self.rootgrp.variables[vname][:, self.cellinfo['ilat'][ncell], self.cellinfo['jlon'][ncell]], fmt="%.4f")


    def dump_by_cell(self,outdir="junkdir"):
        for vname in self.varlist3D:
            outpath = os.path.join(outdir,vname)
            try:
                os.mkdir(outpath)
            except OSError:
                outdir, 'already exists!'
                
            for ncell in range(np.size(self.cellinfo)):
                fname = "%s/data_%4.5f_%4.5f"%(outpath, self.cellinfo['lat'][ncell], self.cellinfo['lon'][ncell])
                np.savetxt(fname, self.rootgrp.variables[vname][:, self.cellinfo['ilat'][ncell], self.cellinfo['jlon'][ncell]], fmt="%.4f")

    def dump_cells(self,outdir="junkdir"):
        for vname in self.varlist3D:
            outpath = os.path.join(outdir,vname)
            try:
                os.mkdir(outpath)
            except OSError:
                outdir, 'already exists!'
            fname = "%s/%s.txt"%(outpath, self.basename)
            ntimes = len(self.rootgrp.variables[vname][:,0,0])
            #ntimes=20
            #np.savetxt(fname, self.rootgrp.variables[vname][:ntimes, self.cellinfo['ilat'], self.cellinfo['jlon']].reshape(ntimes,-1), fmt="%.4f")
            print fname
#             np.savetxt(fname, self.rootgrp.variables[vname][:, self.cellinfo['ilat'], self.cellinfo['jlon']].reshape(ntimes,-1), fmt="%.4f")
#         self.rootgrp.close()
            
    def nc2txtcells(self, ncfilein):
        '''Read nc file'''
        
        print 'Reading', ncfilein
        
        self.basename = os.path.basename(ncfilein)
        print self.basename
        self.rootgrp = Dataset(ncfilein,'r', format=self.nctype)
 
        for ndim, dimobj in enumerate(self.rootgrp.dimensions):
            size = len(self.rootgrp.dimensions[str(dimobj)])
            name = str(dimobj)
            print name, size, len(dimobj)
#             self.__dict__[xyt[ndim]] = dimobj
#             print dimobj
            print ndim
                 
#         ## Reading variables        
        self.varlist2D = []
        self.varlist3D = []
        for varobj in self.rootgrp.variables:
            varname = varobj
            varndim = self.rootgrp.variables[varobj].ndim
            if varndim == 3: 
                self.varlist3D.append(varname)
        print self.varlist3D

def make_nc_dem(nc, infilename, lons, lats, valm):
    outname=infilename+'.nc'
    nc.create_NC(outname)
    nc.add_coord_dims(lons, lats)
    nc.add_var2D(valm, name='ELEV', coordDims=('LAT','LON',), dtype='f', compress=True)

# def main(*args):
#     '''
#     Adapted from mattslib.mapping.arcgrd2map
#     '''
# 
#     nc = PyNC()
#     grid1 = dem2np.Gridder()
#     
#     #outpath = 'extents_maxmin'
#     full_dem_file = args[0] #'dem_0625_pnw_kinggrd.asc'
#     griddic = grid1.gridHead(full_dem_file)
#     
#     if len(args)>1:
#         for arg in args:
#             lats_full, lons_full, topom_full, valm_full = grid1.gridGrabData(arg, griddic)
#             make_nc_dem(nc, arg, lats_full, lons_full, topom_full)
#     else:
#         lats_full, lons_full, topom_full, valm_full = grid1.gridGrabData(full_dem_file, griddic)
#         arg=full_dem_file
#         outname=arg+'.nc'
#         nc.create_NC(outname)
#         nc.add_coord_dims(lons_full, lats_full)
#         nc.add_var2D(topom_full, name='ELEV', coordDims=('LAT','LON',), dtype='f', compress=True)
#     
# if __name__=='__main__':
#     if len(sys.argv) >= 2:
#         main(*sys.argv[1:])
#     else:
#         print 'Require 2+ path arguments <arcdem_common_domain_1> <..> <arcdem_common_domain_N>'