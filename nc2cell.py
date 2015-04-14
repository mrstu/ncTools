#!/usr/bin/env python

import os 
import os.path
import sys

import numpy as np
#import numpy.ma as ma

import pandas as pd

#from collections import OrderedDict as OD

try:
    from netCDF4 import Dataset
    from netCDF4 import num2date, date2num
except:
    print 'No netCDF4 module in path.  Proceeding without netCDF4 module support'

#import PyNC 

def main(args):
    
    ncfilein = args[0]
    outdir   = args[1]
    print args
    readnc(ncfilein)

def readnc(ncfilein):
    '''Read nc file'''
    
    vars='precipitation air_temperature air_temperature_2 wind_speed'.split()
    v1=vars[0]
    print 'Reading', ncfilein
    rgrp = Dataset(ncfilein,'r', format='NETCDF4')
#    print rgrp.variables[v1]

#     coords = rgrp.variables[v1].coordinates.split()
    coords = ['time', 'lat', 'lon']
    print coords
    fillval = rgrp.variables[v1]._FillValue
    c0 = coords[0]
    c1 = coords[1]
    c2 = coords[2]
    
    ntimes=len(rgrp.variables[c0])

    dummy = np.zeros((ntimes,len(vars)))

    df = pd.DataFrame(dummy)
 
    print c0, c1, c2, fillval
     
    i,j = np.where(rgrp.variables[vars[0]][0,:,:].mask == False)
#     print i
#     print j

    for i1, i2 in zip(i[:1], j[:1]):
#     for i1, i2 in zip(i, j):


        ## populate numpy dummy then convert to DataFrame

        ## loop over all input vars
##         for nv, var in enumerate(vars):
##             dummy[:,nv] = rgrp.variables[vars[nv]][:,i1,i2]

        ## don't loop, because of custom handling
#         dummy[:,0] = rgrp.variables[vars[0]][:,i1,i2]
#         dummy[:,1] = rgrp.variables[vars[1]][:,i1,i2] - 273.15
#         dummy[:,2] = rgrp.variables[vars[2]][:,i1,i2] - 273.15
#         dummy[:,3] = rgrp.variables[vars[3]][:,i1,i2]
##         print dummy[:5,:]
#         df = pd.DataFrame(dummy)
            
        ## create dummy DataFrame and populate directly

        df.ix[:,0] = rgrp.variables[vars[0]][:,i1,i2]
        df.ix[:,1] = rgrp.variables[vars[1]][:,i1,i2] - 273.15
        df.ix[:,2] = rgrp.variables[vars[2]][:,i1,i2] - 273.15
        df.ix[:,3] = rgrp.variables[vars[3]][:,i1,i2]

        ## write out csv file
                
        oname = 'data_%4.5f_%4.5f'%(rgrp.variables[c1][i1], rgrp.variables[c2][i2]-360.)
        df.to_csv(oname, sep=" ", float_format="%4.5f", header=False, index=False)

        ## for string representation, have to re-direct stdout to file

        ## myformatter = lambda x: '[%4.1f]' % x
        ## df.to_string(formatters={'A': myformatter, 'C': myformatter})

        ## write without pandas

#     for i1, i2 in zip(i, j):
#         print rgrp.variables[c1][i1], rgrp.variables[c2][i2]
#         
#         for n in range(ntimes):
#             print n
#             print "%f %f %f %f"%(rgrp.variables[vars[0]][n,i1,i2], rgrp.variables[vars[1]][n,i1,i2], rgrp.variables[vars[2]][n,i1,i2], rgrp.variables[vars[3]][n,i1,i2])
         
     


    
if __name__=='__main__':
    print sys.argv
    if len(sys.argv)==3:
        print sys.argv[1:]
        main(sys.argv[1:])
    else:
        print 'Need args <lat_lon_ilat_ilon table> <input netcdf file> <output location>'

# if __name__=='__main__':
#     if len(sys.argv)==4:
#         print sys.argv[1:]
#         main(sys.argv[1:])
#     else:
#         print 'Need args <lat_lon_ilat_ilon table> <input netcdf file> <output location>'