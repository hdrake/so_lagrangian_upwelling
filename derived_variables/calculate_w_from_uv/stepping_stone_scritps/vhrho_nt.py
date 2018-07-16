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

archdir_hist = '/archive/Richard.Slater/'+model+'/'+exp+'/history/'
griddir = '/archive/Carolina.Dufour/Mosaic/'+model+'/'

grid_file = griddir+model+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
lat_t = gridFile.variables['geolat_t'][...]

first_year = 191
last_year = 200
n_timesteps = 73
boundary_north = -28
index_north = np.where(lat_t[:,0]<boundary_north)[0][-1]

for year in range(first_year,last_year+1):
	for time_slice in range(0,n_timesteps):
		start = clock()
		file_name_prefix = str(year).zfill(4)+'0101.ocean_minibling_field_'
		print exp_short_name+', year '+str(year)+', time slice '+str(time_slice)
		
		if (year == first_year) and (time_slice == 0):
			# just get for land mask:
			data = nc.Dataset(archdir_hist + file_name_prefix + 'salt.nc')
			salt = data.variables['salt'][time_slice,:,:index_north,:]
		
		# note: rho_dzt is the grid cell thickness (which varies in time) x 1035
		data = nc.Dataset(archdir_hist + file_name_prefix + 'rho_dzt.nc')
		# use filling trick to avoid spreading mask:
		rho_dzt = (data.variables['rho_dzt'][time_slice,:,:index_north+1,:])\
			.filled(fill_value=1.e20)/1035.
		# compute rho_dzu:
		rho_dzt_right = np.roll(rho_dzt[:,:-1,:],-1,axis=-1)
		rho_dzt_up = rho_dzt[:,1:,:]
		rho_dzt_up_right = np.roll(rho_dzt[:,1:,:],-1,axis=-1)
		rho_dzu = np.minimum(np.minimum(np.minimum(rho_dzt[:,:-1,:],rho_dzt_right),
			rho_dzt_up),rho_dzt_up_right)
		rho_dzu = np.ma.masked_array(rho_dzu,mask=salt.mask)
		rho_dzu = rho_dzu.filled(fill_value=0)
		del rho_dzt_right,rho_dzt_up,rho_dzt_up_right

		#print 'calculating vhrho_nt'
		data = nc.Dataset(archdir_hist + file_name_prefix + 'v.nc')
		v = (data.variables['v'][time_slice,:,:index_north,:]).filled(fill_value=0)
		# compute vhrho_nt:
		vhrho_nt = v*rho_dzu
		del v,rho_dzu
		# this is still on upper right corner of grid cell at this stage, 
		# so interpolate to top center of grid cell:
		vhrho_nt = 0.5*(vhrho_nt+np.roll(vhrho_nt,1,axis=-1))
		vhrho_nt = np.ma.array(vhrho_nt,mask=salt.mask)
		
		vhrho_nt = np.ma.filled(vhrho_nt,fill_value=0)
