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
rickdir = '/archive/Richard.Slater/'+model+'/'+exp+'history/'
archdir_hist = '/archive/wga/'+model+'/'+exp+'pp/ocean_budgets/ts/annual/5yr/'
griddir = '/archive/Carolina.Dufour/Mosaic/'+model+'/'

grid_file = griddir+model+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
lat_t = gridFile.variables['geolat_t'][...]

first_year = 181
last_year = 200
boundary_north = -29.5
index_north = np.where(lat_t[:,0]<boundary_north)[0][-1]

interval_list = ['ocean_budgets.018101-018501.',
                 'ocean_budgets.018601-019001.',
                 'ocean_budgets.019101-019501.',
                 'ocean_budgets.019601-020001.']

# u and v
for interval_ind in range(0,4):
    for j in range(0,5):
        year = 181+j+interval_ind*5
        print 'Doing year '+str(year).zfill(4)
        file_name_prefix = str(year).zfill(4)+'0101.ocean_minibling_field_'
        if ((j==0) & (interval_ind == 0)):
            # Get landmask from salt field
            data = nc.Dataset(rickdir + file_name_prefix + 'salt.nc')
            salt = data.variables['salt'][0,:,:index_north+1,:]
        # rho_dzt is the grid cell thickness (which varies in time)
        data = nc.Dataset(archdir_hist + interval_list[interval_ind] + 'rho_dzt.nc')
        # use filling trick to avoid spreading mask:
        rho_dzt = (data.variables['rho_dzt'][j,:,:index_north+2,:])\
            .filled(fill_value=1.e20)
        data = nc.Dataset(archdir_hist + interval_list[interval_ind] + 'uhrho_et.nc')
        other_uhrho_et = data.variables['uhrho_et'][j,:,:index_north,:]
        data = nc.Dataset(archdir_hist + interval_list[interval_ind] + 'vhrho_nt.nc')
        other_vhrho_nt = data.variables['vhrho_nt'][j,:,:index_north,:]
        # compute rho_dzu (depth, lat, lon) by taking minimum of four grid cell corners
        # Take out last latitude (most northern) and roll longitudes over
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
        
        velpath = '/archive/wga/CM2.6/CM2.6_A_Control-1860_V03/pp/ocean/av/annual_1yr/ocean.'
        print 'calculating vhrho_nt and uhrho_et'
        data = nc.Dataset(velpath+str(year).zfill(4)+'.ann.nc')
        v = (data.variables['v'][0,:,:index_north+1,:]).filled(fill_value=0)
        u = (data.variables['u'][0,:,:index_north+1,:]).filled(fill_value=0)
        
        # compute vhrho_nt and uhrho_ut:
        vhrho_net = v[:,:-1,:]*rho_dzu[:,:-1,:]
        uhrho_net = u[:,:-1,:]*rho_dzu[:,:-1,:]
        uhrho_set = u[:,1:,:]*rho_dzu[:,1:,:]
        del u,v,rho_dzu
        # this is still on upper right corner of grid cell at this stage,
        # so interpolate to top center of grid cell:
        vhrho_nt = 0.5*(vhrho_net+np.roll(vhrho_net,1,axis=-1))
        vhrho_nt = np.ma.array(vhrho_nt,mask=salt[:,:-1,:].mask)
        vhrho_nt = np.ma.filled(vhrho_nt,fill_value=0)

        uhrho_et = 0.5*(uhrho_net+uhrho_set)
        uhrho_et = np.ma.array(uhrho_et,mask=salt[:,:-1,:].mask)
        uhrho_et = np.ma.filled(uhrho_et,fill_value=0)
        #print other_uhrho_et[5,800:805,1730:1735]
        #print uhrho_et[5,800:805,1730:1735]
        
        print 'Percent error of my uhrho_et calculation: '+str(np.average((np.abs(uhrho_et - other_uhrho_et))/other_uhrho_et))
        print 'Percent error of my vhrho_nt calculation: '+str(np.average((np.abs(vhrho_nt - other_vhrho_nt))/other_vhrho_nt))
