#!/usr/bin/env python

# import modules
import netCDF4 as nc
import os,sys
import numpy as np
import jdcal as jc
import math
from wright_eos import *
from jdcal import jd2gcal
from jdcal import gcal2jd
import time as timmee
np.set_printoptions(threshold=np.nan)

model = 'CM2.6'
exp = 'CM2.6_A_Control-1990_V03/'

path2sigma2 = '/archive/wga/'+model+'/'+exp+'/history/'
path2save = '/archive/hfd/'+model+'/'+exp+'/monthly/sigma2/'
htfilename = '/archive/Carolina.Dufour/Mosaic/'+model+'/ocean.static.nc'

files = os.listdir(path2sigma2)
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
    sigma2File = nc.Dataset(path2sigma2+file)
    st_ocean = sigma2File.variables['st_ocean'][:]
    xt_ocean = sigma2File.variables['xt_ocean'][:]
    yt_ocean = sigma2File.variables['yt_ocean'][:]
    lat_n_index = (np.abs(yt_ocean - (-29.5))).argmin()
    lat_s_index = (np.abs(yt_ocean - (-81))).argmin()
    yt_ocean = yt_ocean[lat_s_index:lat_n_index]
    time = sigma2File.variables['time'][:]
    for i in range(0,np.size(time)):
        temp = sigma2File.variables['temp'][i,:,lat_s_index:lat_n_index,:]
        salt = sigma2File.variables['salt'][i,:,lat_s_index:lat_n_index,:]
        
        # Save sigma2 to Netcdf4 file
        outfile = path2save+'nest_1_'+year+str(i+1).zfill(2)+'15000000d.nc'
        print 'Saving sigma2 in: nest_1_'+year+str(i+1).zfill(2)+'15000000d.nc'
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
        
        sigma2_var = netcdf_file.createVariable('sigma2', 'f4', ('time','st_ocean','yt_ocean', 'xt_ocean',),fill_value=-1.e20)
        sigma2_var.units = 'm'
        sigma2_var.long_name = 'sigma2 (potential density referenced to 2000m)'
        sigma2_var.missing_value =  -1.e+20
        
        # data
        ti[:] = time[i]
        st[:] = st_ocean
        yt[:] = yt_ocean
        xt[:] = xt_ocean
        sigma2_var[0,:,:,:] = wright_eos(temp,salt,2)
    
        netcdf_file.close()
sys.exit()
####################################################################


