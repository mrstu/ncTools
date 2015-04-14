#!/usr/bin/env python

import os 
import os.path
import sys

import numpy as np
import numpy.ma as ma

from collections import OrderedDict as OD

try:
    from netCDF4 import Dataset
    from netCDF4 import num2date, date2num
except:
    print 'No netCDF4 module in path.  Proceeding without netCDF4 module support'

import PyNC 

def main(args):
#     LLi_file = "goodlatlons.txt"
#     #LLi_file = "goodlatlons_1000.txt"
#     outdir   = "/raid3/stumbaugh/IS/CONUS/v2/nc2forcings_dump"
#     ncfilein = "CONUS_v2/macav2livneh_was_daily_CanESM2_historical_1950_1969_CONUS.nc"
#     ncfilein = "/home/raid3/stumbaugh/IS/CONUS/v2/splitmon/was/historical/CanESM2/1950-01-01_1950-01-31"
    
    ncfilein = args[0]
    outdir   = args[1]
    print args
    nc = PyNC.PyNC()
    nc.readnc(ncfilein)
    nc.nc2txtcells(ncfilein)
    #nc.dump_by_cell(outdir) #One cell per file
    nc.dump_active_cells(outdir) # All cells per file
    
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