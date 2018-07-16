#!/usr/efe/evn python
# -*- coding: utf-8 -*-

import netCDF4 as nc
import os,sys
import numpy as np
from wright_eos import *
from jdcal import jd2gcal

path2files = '/archive/rds/CM2.6/CM2.6_A_Control-1860_V03/history/'

sigmalist = ['01900101']

#sigmalist = ['01810101',
#             '01820101',
#             '01830101',
#             '01840101',
#             '01850101',
#             '01860101',
#             '01870101',
#             '01880101',
#             '01890101',
#             '01900101',
#             '01910101',
#             '01920101',
#             '01930101',
#             '01940101',
#             '01950101',
#             '01960101',
#             '01970101',
#             '01980101',
#             '01990101',
#             '02000101',
#             '02010101',
#             '02020101']

for interv in sigmalist:
    print interv
    # Size of these variables is 73 x 50 x 2700 x 3600
    # 'salt' is Practical Salinity (time,st_ocean,yt_ocean,xt_ocean)
    # 'temp' is Potential Temperature (time,st_ocean,yt_ocean,xt_ocean)
    ncSalt = nc.Dataset(path2files+interv+'.ocean_minibling_field_salt.nc')
    ncTemp = nc.Dataset(path2files+interv+'.ocean_minibling_field_temp.nc')
    
    # Extract dimensions
    times = ncSalt.variables['time'][...]
    st_ocean = ncSalt.variables['st_ocean'][...]
    yt_ocean = ncSalt.variables['yt_ocean'][...]
    xt_ocean = ncSalt.variables['xt_ocean'][...]
    ntime = len(times)
    
    # Find lat indices for southern ocean:
    lat_n_index = (np.abs(yt_ocean - (-29.5))).argmin()
    lat_s_index = (np.abs(yt_ocean - (-81))).argmin()
    yt_ocean = yt_ocean[lat_s_index:lat_n_index]
    
    # Turn times into month / days
    months = np.zeros(ntime)
    days = np.zeros(ntime)
    for t in range(ntime):
        months[t] = jd2gcal(1721423,times[t])[1]
        days[t] = jd2gcal(1721423,times[t])[2]
    
    for fiveday in range(ntime):
        # Extract salt and temp
        salt = ncSalt.variables['salt'][fiveday,:,lat_s_index:lat_n_index,:]
        temp = ncTemp.variables['temp'][fiveday,:,lat_s_index:lat_n_index,:]

        path2save = '/archive/hfd/CM2.6/CM2.6_A_Control-1860_V03/5day/sigma2/'
        # Create sigma2 netcdf file
        outfile = path2save+interv[0:4]+str(int(months[fiveday])).zfill(2)+str(int(days[fiveday])).zfill(2)+'.ocean_5day_sigma2.nc'
        print 'Saving calculation in: '+interv[0:4]+str(int(months[fiveday])).zfill(2)+str(int(days[fiveday])).zfill(2)+'.ocean_5day_sigma2.nc'
        netcdf_file = nc.Dataset(outfile, 'w', format='NETCDF4')
        netcdf_file.description = 'Sigma2 for '+interv

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
        ti = times[fiveday]
        st[:] = st_ocean
        yt[:] = yt_ocean
        xt[:] = xt_ocean

        sigma2 = wright_eos(temp,salt,2)
        s2[0,:,:,:] = sigma2

        netcdf_file.close()

    ncSalt.close()
    ncTemp.close()

sys.exit()


