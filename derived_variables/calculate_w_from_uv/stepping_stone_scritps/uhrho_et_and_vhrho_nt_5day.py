#!/usr/bin/env python

# import modules
import os,sys,glob
import numpy as np
import netCDF4 as nc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.cm as cm
import matplotlib.colorbar as cb
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator, LinearLocator, NullLocator
np.set_printoptions(threshold=np.nan)

##############

model = 'CM2.6'
exp = 'CM2.6_A_Control-1860_V03/'
exp_short_name = '1860'
path2save = '/archive/Henri.Drake/'+model+'/'+exp+'5day/for_CMS_SO/'
archdir_hist = '/archive/Richard.Slater/'+model+'/'+exp+'/history/'
other_dir = '/archive/wga/'+model+'/'+exp+'pp/ocean_budgets/ts/annual/5yr/'

# Get Grid Data
griddir = '/archive/Carolina.Dufour/Mosaic/'+model+'/'
grid_file = griddir+model+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
lat_t = gridFile.variables['geolat_t'][...]

first_year = 181
last_year = 200
boundary_north = -29.5
index_north = np.where(lat_t[:,0]<boundary_north)[0][-1]

dxt = np.zeros((50,935,3600))
dyt = np.zeros((50,935,3600))
dxt_uno = gridFile.variables['dxt'][:index_north,:]
dyt_uno = gridFile.variables['dyt'][:index_north,:]
for i in range(0,50):
    dxt[i,:,:]=dxt_uno
    dyt[i,:,:]=dyt_uno

files = os.listdir(path2save)
files = np.sort(files)

interval_list = ['ocean_budgets.018101-018501.',
                 'ocean_budgets.018601-019001.',
                 'ocean_budgets.019101-019501.',
                 'ocean_budgets.019601-020001.']

countu5day = 0
countfiles = 0
#### COMPUTE uhrho_et and vhrho_nt 5-day averages from u,v,rho_dzt 5-day averages

#935 lat indices for north = -29.5
up_flux = np.zeros((50,935,3600))
# u and v
for file in files:
    if not file.endswith('u.nc'):
        continue
    print 'Doing '+file
    year = file.split('.',3)[1][0:4]
    file_name_prefix = str(year)+'0101.ocean_minibling_field_'
    if (year == str(first_year).zfill(4)) and (countu5day == 0):
        # Get landmask from salt field
        data = nc.Dataset(archdir_hist + file_name_prefix + 'salt.nc')
        salt = data.variables['salt'][0,:,:index_north,:]
        lon = data.variables['xt_ocean'][:]
        lat = data.variables['yt_ocean'][:index_north]
        depth = data.variables['st_ocean'][:]
    # rho_dzt is the grid cell thickness (which varies in time) x 1035
    data = nc.Dataset(archdir_hist + file_name_prefix + 'rho_dzt.nc')
    # use filling trick to avoid spreading mask:
    rho_dzt = (data.variables['rho_dzt'][countu5day,:,:index_north+1,:])\
        .filled(fill_value=1.e20)/1035.
    # compute rho_dzu (depth, lat, lon) by taking minimum of four grid cell corners
    # Take out last latitude (most northern?) and roll longitudes over
    rho_dzt_right = np.roll(rho_dzt[:,:-1,:],-1,axis=-1)
    # Take out first latitude
    rho_dzt_up = rho_dzt[:,1:,:]
    # Take out first latitude and roll longitudes over
    rho_dzt_up_right = np.roll(rho_dzt[:,1:,:],-1,axis=-1)
    # Take minimum of thickness of tracer grid box over 4 corners to be the grid thickness
    rho_dzu = np.minimum(np.minimum(np.minimum(rho_dzt[:,:-1,:],rho_dzt_right),rho_dzt_up),rho_dzt_up_right)
    # Mask if ground or if 0 thickness
    rho_dzu = np.ma.masked_array(rho_dzu,mask=salt.mask)
    rho_dzu = rho_dzu.filled(fill_value=0)
    del rho_dzt_right,rho_dzt_up,rho_dzt_up_right

    print 'calculating vhrho_nt and uhrho_nt'
    data = nc.Dataset(archdir_hist + file_name_prefix + 'v.nc')
    v = (data.variables['v'][countu5day,:,:index_north,:]).filled(fill_value=0)
    data = nc.Dataset(archdir_hist + file_name_prefix + 'u.nc')
    u = (data.variables['u'][countu5day,:,:index_north,:]).filled(fill_value=0)
    # compute vhrho_nt and uhrho_ut:
    vhrho_net = v*rho_dzu
    vhrho_nwt = np.roll(vhrho_net,1,axis=-1)
    uhrho_net = u*rho_dzu
    uhrho_set = u[:,1:,:]*rho_dzu[:,1:,:]
    uhrho_set = np.insert(uhrho_set,0,uhrho_set[:,0,:],axis=1)
    del u,v,rho_dzu
    # this is still on upper right corner of grid cell at this stage,
    # so interpolate to top center of grid cell:
    vhrho_nt = 0.5*(vhrho_net+vhrho_nwt)
    vhrho_nt = np.ma.array(vhrho_nt,mask=salt.mask)
    vhrho_nt = np.ma.filled(vhrho_nt,fill_value=0)

    uhrho_et = 0.5*(uhrho_net+uhrho_set)
    uhrho_et = np.ma.array(uhrho_et,mask=salt.mask)
    uhrho_et = np.ma.filled(uhrho_et,fill_value=0)

    print 'Calculating w'
    # Create volume flux out (positive) and in (negative) to each tracer grid cell in each of 4 horizontal directions
    trans_nt = vhrho_nt * dxt
    trans_nt = np.ma.array(trans_nt,mask=salt.mask)
    trans_nt = np.ma.filled(trans_nt,fill_value=0)
    
    trans_st = -(vhrho_nt[:,:-1,:]*dxt[:,:-1,:])
    trans_st = np.insert(trans_st,0,trans_st[:,0,:],axis=1)
    trans_st = np.ma.array(trans_st,mask=salt.mask)
    trans_st = np.ma.filled(trans_st,fill_value=0)
    
    trans_et = uhrho_et * dyt
    trans_et = np.ma.array(trans_et,mask=salt.mask)
    trans_et = np.ma.filled(trans_et,fill_value=0)

    trans_wt = -(np.roll(trans_et,1,axis=-1))
    trans_wt = np.ma.array(trans_wt,mask=salt.mask)
    trans_wt = np.ma.filled(trans_wt,fill_value=0)
    
    # vert_flux defined as positive upwards
    vert_flux = trans_nt + trans_st + trans_et + trans_wt
    
    for i in range(0,50):
        p = 49-i
        if i == 0:
            up_flux[i,:,:]=-vert_flux[i,:,:]
        else:
            up_flux[i,:,:]=-(vert_flux[i,:,:]+up_flux[i,:,:])

    area_t = dxt*dyt
    w = up_flux/area_t

    index_south = (np.abs(lat + 81)).argmin()
    w = w[:,index_south:,:]
    st_ocean = depth
    xt_ocean = lon
    yt_ocean = lat[index_south:]

    u_File = nc.Dataset(path2save+file)
    tim = u_File.variables['time'][:]
    average_DT = u_File.variables['average_DT'][:]

    ##############################################################
    # Save w
    out_filename_w = path2save+file[:-4]+'w.nc'
    print 'Saving w in: '+file[:-4]+'w.nc'
    netcdf_file = nc.Dataset(out_filename_w,'w',format='NETCDF4')
    netcdf_file.description = exp

    # dimensions
    netcdf_file.createDimension('time', None)
    netcdf_file.createDimension('st_ocean', len(st_ocean))
    netcdf_file.createDimension('yu_ocean', len(yu_ocean))
    netcdf_file.createDimension('xu_ocean', len(xu_ocean))

    # variables
    ti = netcdf_file.createVariable('time', 'f4', ('time',))
    ti.units = 'days since 0001-01-01 00:00:00'
    ti.calender_type = 'JULIAN'
    ti.calender = 'JULIAN'

    st = netcdf_file.createVariable('st_ocean', 'f4', ('st_ocean',))
    st.units = 'metres'
    st.long_name = 'tcell zstar depth'
    st.positive = 'down'

    yu = netcdf_file.createVariable('yu_ocean', 'f4', ('yu_ocean',))
    yu.units = 'degrees_N'
    yu.long_name = 'ucell latitude'

    xu = netcdf_file.createVariable('xu_ocean', 'f4', ('xu_ocean',))
    xu.units = 'degrees_E'
    xu.long_name = 'ucell longitude'

    dt = netcdf_file.createVariable('average_DT', 'f4', ('time',))
    dt.units = 'days'
    dt.long_name = 'length of average period'

    w_var = netcdf_file.createVariable('w', 'f4', ('time','st_ocean','yu_ocean', 'xu_ocean',),fill_value=-1.e20)
    w_var.units = 'm/s'
    w_var.long_name = 'vertical velocity'
    w_var.missing_value =  -1.e+20
   
    # data
    ti[:] = tim
    st[:] = st_ocean
    yu[:] = yu_ocean
    xu[:] = xu_ocean
    #u_var[:,:,:,0] = np.swapaxes(u,2,0)
    w_var[0,:] = w
    dt[:] = average_DT

    netcdf_file.close()

    ####################################################################


    ##### Counter
    countfiles += 1
    countu5day += 1
    if countu5day == 73:
        countu5day = 0