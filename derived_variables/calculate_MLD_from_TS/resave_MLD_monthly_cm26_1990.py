#!/usr/bin/env python

# import modules
import netCDF4 as nc
import os,sys
import numpy as np
import jdcal as jc
import math
from jdcal import jd2gcal
from jdcal import gcal2jd
import time as timmee
np.set_printoptions(threshold=np.nan)

model = 'CM2.6'
exp = 'CM2.6_A_Control-1990_V03/'

path2mld = '/archive/wga/'+model+'/'+exp+'/history/'
path2save = '/archive/hfd/'+model+'/'+exp+'/monthly/MLD/'
htfilename = '/archive/Carolina.Dufour/Mosaic/'+model+'/ocean.static.nc'

files = os.listdir(path2mld)
files = np.sort(files)
count = 0
for file in files:
    if not file.endswith('ocean.nc'):
        continue
    count += 1
    year = file[0:4]
    if (int(year) < 85) or (int(year) > 96):
        continue

    print 'Loading '+file
    mldFile = nc.Dataset(path2mld+file)
    st_ocean = mldFile.variables['st_ocean'][:]
    xt_ocean = mldFile.variables['xt_ocean'][:]
    yt_ocean = mldFile.variables['yt_ocean'][:]
    time = mldFile.variables['time'][:]
    for i in range(0,np.size(time)):
        MLD = mldFile.variables['mld'][i,:,:]
        
        # Save MLD to Netcdf4 file
        outfile = path2save+'nest_1_'+year+str(i+1).zfill(2)+'15000000m.nc'
        print 'Saving MLD in: nest_1_'+year+str(i+1).zfill(2)+'15000000m.nc'
        if os.path.lexists(outfile):
            os.remove(outfile)
        netcdf_file = nc.Dataset(outfile, 'w', format='NETCDF4')
        netcdf_file.description = exp

        # dimensions
        netcdf_file.createDimension('time', None)
        netcdf_file.createDimension('st_ocean', len(st_ocean))
        netcdf_file.createDimension('yt_ocean', len(yt_ocean))
        netcdf_file.createDimension('xt_ocean', len(xt_ocean))
        
        # variables
        ti = netcdf_file.createVariable('time', 'f4', ('time',))
        ti.units = 'days since 0001-01-01 00:00:00'
        ti.calender_type = 'JULIAN'
        ti.calender = 'JULIAN'
        
        st = netcdf_file.createVariable('st_ocean', 'f4', ('st_ocean',))
        st.units = 'metres'
        st.long_name = 'tcell zstar depth'
        st.positive = 'down'
        
        yt = netcdf_file.createVariable('yt_ocean', 'f4', ('yt_ocean',))
        yt.units = 'degrees_N'
        yt.long_name = 'tcell latitude'
        
        xt = netcdf_file.createVariable('xt_ocean', 'f4', ('xt_ocean',))
        xt.units = 'degrees_E'
        xt.long_name = 'tcell longitude'
        
        MLD_var = netcdf_file.createVariable('MLD', 'f4', ('time','yt_ocean', 'xt_ocean',),fill_value=-1.e20)
        MLD_var.units = 'm'
        MLD_var.long_name = 'Mixed Layer Depth - depth of surface equal to (sigma2 at surface plus 0.03)'
        MLD_var.missing_value =  -1.e+20
        
        # data
        ti[:] = time
        st[:] = st_ocean
        yt[:] = yt_ocean
        xt[:] = xt_ocean
        MLD_var[0,:] = MLD
    
        netcdf_file.close()
sys.exit()
####################################################################


