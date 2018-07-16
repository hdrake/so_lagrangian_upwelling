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


#cms_exp = 'cesm1_100year_monthly'
cms_exp = 'cm26_original_100year_monthly'

# Model name from cms_exp
if cms_exp[0:4] == 'cm21':
    model = 'CM2-1deg'
    exp = 'CM2-1deg_A_Control-1860_V03'
elif cms_exp[0:4] == 'cm25':
    model = 'CM2.5'
    exp = 'CM2.5_A_Control-1860_Y03-MBLING-BLING'
elif cms_exp[0:4] == 'cm26':
    model = 'CM2.6'
    exp = 'CM2.6_A_Control-1860_V03'
    numprocs = 6
    zfillnum = 1
elif cms_exp[0:5] == 'cesm1':
    model = 'CESM1'
    exp = 'CESM1_A_Control-1860_V03'
    numprocs = 32
    zfillnum = 2

path2MLD = '/archive/hfd/'+model+'/'+exp+'/monthly/MLD/'
htfilename = '/archive/Carolina.Dufour/Mosaic/CM2.6/ocean.static.nc'
htfile = nc.Dataset(htfilename)
ht = htfile.variables['ht'][:933,:]
htfile.close()

firstfile=1

file = 'testing_MLD1_'

nc_sigma2 = nc.Dataset('/archive/hfd/CM2.6/CM2.6_A_Control-1860_V03/5day/sigma2/01811009.ocean_5day_sigma2.nc')
#nc_sigma2 = nc.Dataset('/archive/imf/CM2.6/CM2.6_A_Control-1860_V03.ocean.sigma2_018101-018512.nc')
sigma2 = nc_sigma2.variables['sigma2'][0,:,:933,:]
sigma2_mask = sigma2.mask[:,:,:]
if firstfile:
    st_ocean = nc_sigma2.variables['st_ocean'][:]
    xt_ocean = nc_sigma2.variables['xt_ocean'][:]
    yt_ocean = nc_sigma2.variables['yt_ocean'][:933]
    time = nc_sigma2.variables['time'][:]
    firstfile = 0
MLD_sigma2 = sigma2[0,:,:]+0.03
sigma2[sigma2_mask] = 1.e+20
MLD = np.minimum.reduce([st_ocean[np.argmax(sigma2 > MLD_sigma2, axis=0)],ht])
MLD = np.ma.array(MLD,mask=sigma2_mask[0,:,:])
MLD[MLD.mask]=-1.e+20
# Save MLD to Netcdf4 file
outfile = path2MLD+file[:-4]+'m.nc'
print 'Saving MLD in: '+outfile
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


