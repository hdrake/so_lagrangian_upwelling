#!/usr/bin/env python

# import modules
import os,sys,glob
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

##############

model = 'CM2.6'
exp = 'CM2.6_A_Control-1860_V03/'
exp_short_name = '1860'
path2save = '/archive/Henri.Drake/'+model+'/'+exp+'5day/for_CMS_SO/'
archdir_hist = '/archive/Richard.Slater/'+model+'/'+exp+'/history/'
griddir = '/archive/Carolina.Dufour/Mosaic/'+model+'/'

grid_file = griddir+model+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
lat_t = gridFile.variables['geolat_t'][...]

first_year = 181
last_year = 200
boundary_north = -29.5
index_north = np.where(lat_t[:,0]<boundary_north)[0][-1]

files = os.listdir(path2save)
files = np.sort(files)
countu5day = 0
# u and v
for file in files:
    if not file.endswith('u.nc'):
        continue
    print 'Doing '+file
    year = file.split('.',3)[1][0:4])
    file_name_prefix = str(year)+'0101.ocean_minibling_field_'
    if (year == str(first_year).zfill(4)) and (countu5day == 0):
        # Get landmask from salt field
        data = nc.Dataset(archdir_hist + file_name_prefix + 'salt.nc')
        salt = data.variables['salt'][0,:,:index_north,:]
    # rho_dzt is the grid cell thickness (which varies in time) x 1035
    data = nc.Dataset(archdir_hist + file_name_prefix + 'rho_dzt.nc')
    # use filling trick to avoid spreading mask:
    rho_dzt = (data.variables['rho_dzt'][countu5day,:,:index_north+1,:])\
        .filled(fill_value=1.e20)/1035.0
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
    vhrho_nt = v*rho_dzu
    uhrho_et = u*rho_dzu
    del u,v,rho_dzu
    # this is still on upper right corner of grid cell at this stage,
    # so interpolate to top center of grid cell:
    vhrho_nt = 0.5*(vhrho_nt+np.roll(vhrho_nt,1,axis=-1))
    vhrho_nt = np.ma.array(vhrho_nt,mask=salt.mask)
    vhrho_nt = np.ma.filled(vhrho_nt,fill_value=0)

    uhrho_et = 0.5*(uhrho_nt+np.roll(uhrho_nt,1,axis=-1))
    uhrho_et = np.ma.array(uhrho_nt,mask=salt.mask)
    uhrho_et = np.ma.filled(uhrho_nt,fill_value=0)
