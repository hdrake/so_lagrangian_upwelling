#!/usr/efe/evn python
# -*- coding: utf-8 -*-

import netCDF4 as nc
import os,sys
import numpy as np
from jdcal import jd2gcal
import matplotlib.pyplot as plt
np.set_printoptions(threshold=np.nan)

# Name of cms experiment
cms_exp = 'cm26_200year_5day'

# Model name from cms_exp
if cms_exp[0:4] == 'cm21':
    model = 'CM2-1deg'
    exp = 'CM2-1deg_A_Control-1860_V03'
    temporal_resolution = '5day'

elif cms_exp[0:4] == 'cm25':
    model = 'CM2.5'
    exp = 'CM2.5_A_Control-1860_Y03-MBLING-BLING'
    temporal_resolution = '5day'

elif cms_exp[0:4] == 'cm26':
    model = 'CM2.6'
    exp = 'CM2.6_A_Control-1860_V03'
    temporal_resolution = '5day'

elif cms_exp[0:5] == 'cesm1':
    model = 'CESM1'
    exp = 'CESM1_A_Control-1860_V03'
    temporal_resolution = 'monthly'

path_to_sigma2 = '/archive/hfd/'+model+'/'+exp+'/'+temporal_resolution+'/sigma2/'
path_to_isopycnal = '/archive/hfd/'+model+'/'+exp+'/'+temporal_resolution+'/isopycnal/'

files = os.listdir(path_to_sigma2)
files = np.sort(files)
isopyc_368 = np.zeros((np.size(files),933,3600))
isopyc_370 = np.zeros((np.size(files),933,3600))
PV_368 = np.zeros((np.size(files),933,3600))
PV_370 = np.zeros((np.size(files),933,3600))

count = 0
for file in files:
    ncFile = nc.Dataset(path_to_sigma2+file)
    sigma2 = ncFile.variables['sigma2'][0,:,:,:]-1000.0
    xt_ocean = ncFile.variables['xt_ocean'][:]
    yt_ocean = ncFile.variables['yt_ocean'][:]
    st_ocean = ncFile.variables['st_ocean'][:]
    time = ncFile.variables['time'][0]
    
    coriolis = np.zeros((933,3600))
    for ylat in range(np.size(yt_ocean)):
        coriolis[ylat,:] = 2.*(7.2921e-5)*np.sin(yt_ocean[ylat]*2*np.pi/360.)
    
    inds_368 = np.argmax(sigma2 > 36.8, axis=0)
    isopycnal_368 = st_ocean[inds_368]
    isopycnal_368 = np.ma.array(isopycnal_368, mask = isopycnal_368 < 6)

    ## 36.8 isopycnal
    above_inds_368 = np.maximum.reduce([inds_368-1,np.zeros(np.shape(inds_368))]).astype(int)
    sigma2_above_368 = np.zeros((933,3600))
    sigma2_below_368 = np.zeros((933,3600))
    z_above_368 = np.zeros((933,3600))
    z_below_368 = np.zeros((933,3600))
    for j in range(933):
        for i in range(3600):
            sigma2_above_368[j,i] = sigma2[above_inds_368[j,i],j,i]
            sigma2_below_368[j,i] = sigma2[inds_368[j,i]+1,j,i]
            z_above_368[j,i] = st_ocean[above_inds_368[j,i]]
            z_below_368[j,i] = st_ocean[inds_368[j,i]+1]

    drho_dz_368 = (sigma2_above_368-sigma2_below_368)/(z_above_368/z_below_368)
    PV_368[count,:,:] = coriolis * (drho_dz_368 / 36.8)
    PV_368[count,:,:] = np.ma.array(PV_368[count,:,:], mask = isopycnal_368 < 6)

    ## 37.0 isopycnal
    inds_370 = np.argmax(sigma2 > 37.0, axis=0)
    isopycnal_370 = st_ocean[inds_370]
    isopycnal_370 = np.ma.array(isopycnal_370, mask = isopycnal_370 < 6)

    above_inds_370 = np.maximum.reduce([inds_370-1,np.zeros(np.shape(inds_370))]).astype(int)
    sigma2_above_370 = np.zeros((933,3600))
    sigma2_below_370 = np.zeros((933,3600))
    z_above_370 = np.zeros((933,3600))
    z_below_370 = np.zeros((933,3600))
    for j in range(933):
        for i in range(3600):
            sigma2_above_370[j,i] = sigma2[above_inds_370[j,i],j,i]
            sigma2_below_370[j,i] = sigma2[inds_370[j,i]+1,j,i]
            z_above_370[j,i] = st_ocean[above_inds_370[j,i]]
            z_below_370[j,i] = st_ocean[inds_370[j,i]+1]
    
    drho_dz_370 = (sigma2_above_370-sigma2_below_370)/(z_above_370/z_below_370)
    PV_370[count,:,:] = coriolis * (drho_dz_370 / 37.0)
    PV_370[count,:,:] = np.ma.array(PV_370[count,:,:], mask = isopycnal_370 < 6)

    # Save Isopycnal depths to Netcdf4 file
    outfile = path_to_isopycnal+'nest_1_'+file[0:8]+'000000i.nc'
    print 'Saving Isopycnal depths in: '+outfile
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

    iso368_var = netcdf_file.createVariable('isopycnal368', 'f4', ('time','yt_ocean', 'xt_ocean',),fill_value=-1.e20)
    iso368_var.units = 'm'
    iso368_var.long_name = 'Isopycnal Depth - depth of 36.8 isopycnal)'
    iso368_var.missing_value =  -1.e+20

    iso370_var = netcdf_file.createVariable('isopycnal370', 'f4', ('time','yt_ocean', 'xt_ocean',),fill_value=-1.e20)
    iso370_var.units = 'm'
    iso370_var.long_name = 'Isopycnal Depth - depth of 37.0 isopycnal)'
    iso370_var.missing_value =  -1.e+20

    PV368_var = netcdf_file.createVariable('PV368', 'f4', ('time','yt_ocean', 'xt_ocean',),fill_value=-1.e20)
    PV368_var.units = 'm'
    PV368_var.long_name = 'PV on 36.8 isopycnal)'
    PV368_var.missing_value =  -1.e+20

    PV370_var = netcdf_file.createVariable('PV370', 'f4', ('time','yt_ocean', 'xt_ocean',),fill_value=-1.e20)
    PV370_var.units = 'm'
    PV370_var.long_name = 'PV on 37.0 isopycnal)'
    PV370_var.missing_value =  -1.e+20

    # data
    ti[:] = time
    st[:] = st_ocean
    yt[:] = yt_ocean
    xt[:] = xt_ocean
    iso368_var[0,:,:] = isopycnal_368
    iso370_var[0,:,:] = isopycnal_370
    PV368_var[0,:,:] = PV_368[count,:,:]
    PV370_var[0,:,:] = PV_370[count,:,:]

    isopyc_368[count,:,:] = isopycnal_368
    isopyc_370[count,:,:] = isopycnal_370
    count+=1

    netcdf_file.close()

np.save(path_to_isopycnal+'isopycnal_368.npy',isopyc_368)
np.save(path_to_isopycnal+'isopycnal_370.npy',isopyc_370)
np.save(path_to_isopycnal+'PV_368.npy',PV_368)
np.save(path_to_isopycnal+'PV_370.npy',PV_370)


