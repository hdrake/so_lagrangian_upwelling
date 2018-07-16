#!/usr/efe/evn python
# -*- coding: utf-8 -*-

import netCDF4 as nc
import os,sys
import numpy as np
from wright_eos import *
from jdcal import jd2gcal

path2files = '/archive/hfd/CM2.6/CM2.6_A_Control-1860_V03/monthly/sigma2/'
numyears = 20
firstyear = 67

files = os.listdir(path2files)
files = np.sort(files)

for file in files:
    if not file.endswith('s.nc'):
        continue
    year = file[7:11]
    month = file[11:13]
    print file
    # Size of these variables is 73 x 50 x 2700 x 3600
    # 'salt' is Practical Salinity (time,st_ocean,yt_ocean,xt_ocean)
    # 'temp' is Potential Temperature (time,st_ocean,yt_ocean,xt_ocean)
    ncSalt = nc.Dataset(path2files+file)
    ncTemp = nc.Dataset(path2files+file[:-4]+'t.nc')

    # Extract dimensions
    times = ncSalt.variables['time'][...]
    st_ocean = ncSalt.variables['st_ocean'][:]

    yt_ocean = ncSalt.variables['yt_ocean'][:]
    xt_ocean = ncSalt.variables['xt_ocean'][:]
    ntime = len(times)
    
    # Find lat indices for southern ocean:
    lat_n_index = (np.abs(yt_ocean - (-29.5))).argmin()
    lat_s_index = (np.abs(yt_ocean - (-81))).argmin()
    yt_ocean = yt_ocean[lat_s_index:lat_n_index]

    # Extract salt and temp
    salt = ncSalt.variables['salt'][0,:,lat_s_index:lat_n_index,:]
    temp = ncTemp.variables['temp'][0,:,lat_s_index:lat_n_index,:]

    path2save = '/archive/hfd/CM2.6/CM2.6_A_Control-1860_V03/monthly/sigma2/'

    # Create sigma2 netcdf file
    outfile = path2save+'nest_1_'+year+month+'15'+'000000'+'d.nc'
    print 'Saving calculation in: nest_1_'+year+month+'15'+'000000'+'d.nc'
    if os.path.lexists(outfile):
        os.remove(outfile)
    netcdf_file = nc.Dataset(outfile, 'w', format='NETCDF4')
    netcdf_file.description = 'Monthly Sigma2 for year'+year+', month '+month

    # Dimensions
    netcdf_file.createDimension('time',1)
    netcdf_file.createDimension('st_ocean', len(st_ocean))
    netcdf_file.createDimension('yt_ocean', len(yt_ocean))
    netcdf_file.createDimension('xt_ocean', len(xt_ocean))

    # Variables
    ti = netcdf_file.createVariable('time','f4',('time',))
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

    s2 = netcdf_file.createVariable('sigma2', 'f4', ('time','st_ocean','yt_ocean','xt_ocean',),fill_value=-1.e20)
    s2.units = 'kg/m^3'
    s2.long_name = 't-cell rho'

    # Data
    ti = times[0]
    st[:] = st_ocean
    yt[:] = yt_ocean
    xt[:] = xt_ocean

    sigma2 = wright_eos(temp,salt,2)
    s2[0,:,:,:] = sigma2

    netcdf_file.close()

    ncSalt.close()
    ncTemp.close()

sys.exit()


