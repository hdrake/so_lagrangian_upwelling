#!/usr/bin/env python

# import modules
import netCDF4 as nc
import os,sys
import numpy as np
import math
np.set_printoptions(threshold=np.nan)

#########################################################################################################
"""
    This script takes as input individual 5-day files containing C-grid velocities:
    U-vel (with dimensions time (length 1), st_ocean, yu_ocean, xu_ocean),
    V-vel (with dimensions time (length 1), st_ocean, yu_ocean, xu_ocean),
    W-vel (with dimensions time (length 1), sw_ocean, yt_ocean, xt_ocean)
    
    and resaves the files in a different directory into 5-day files containing u,v,w on
    an A-grid with dimensions:
    time (length 1), st_ocean (+1 length of above files for bottom boundary condition), yu_ocean, xu_ocean.

"""
############################################################################################################

# Variables for file locations
user = 'Henri.Drake'
model = 'CM2-1deg'
exp = 'CM2-1deg_A_Control-1860_V03/'

# Path to velocity data:
path2files = '/archive/'+user+'/'+model+'/'+exp+'/5day/for_CMS_SO_residual/'

# Path to save velocity files
path2save   = '/archive/'+user+'/'+model+'/'+exp+'/5day/for_CMS_SO_residual_Edited/'
if not os.path.lexists(path2save):
    os.mkdir(path2save)

# Define northern boundary to save velocities:
lat_north = -29.9

files = os.listdir(path2files)
files = np.sort(files)

# Go through each different file ending with 'u.nc' to get the date prefix and then extract u,v,w from that 5-day period
# Make sure to extract the u-cell mask.
for file in files:
    if not file.endswith('u.nc'):
        continue
    year = int(file[6:10])
    month = int(file[10:12])
    day = int(file[12:14])
    u_file = path2files+file
    v_file = path2files+file[:-4]+'v.nc'
    w_file = path2files+file[:-4]+'wt.nc'
    uncFile = nc.Dataset(u_file)
    vncFile = nc.Dataset(v_file)
    wncFile = nc.Dataset(w_file)
    if ((year == 181) and (month == 1) and (day == 2)):
        st_ocean = uncFile.variables['st_ocean'][...]
        yu_ocean = uncFile.variables['yu_ocean'][...]
        xu_ocean = uncFile.variables['xu_ocean'][...]
    
        # find lat/long indices smaller domain:
        lat_n_index = len(yu_ocean)-2
        umask = uncFile.variables['u'][0,:,:lat_n_index,:].mask
        yu_ocean = yu_ocean[:lat_n_index]

    u = (uncFile.variables['u'][0,:,:lat_n_index,:]).filled(fill_value = 0)
    # Change this probably
    time = uncFile.variables['time'][0]
    uncFile.close()

    v = (vncFile.variables['v'][0,:,:lat_n_index,:]).filled(fill_value = 0)
    vncFile.close()

    wmask = wncFile.variables['wt'][0,:,:lat_n_index+1,:].mask
    w = (wncFile.variables['wt'][0,:,:lat_n_index+1,:]).filled(fill_value = 0)
    sw_ocean = wncFile.variables['sw_ocean'][...]
    xt_ocean = wncFile.variables['xt_ocean'][...]
    yt_ocean = wncFile.variables['yt_ocean'][:lat_n_index]
    wncFile.close()

    diff_bottom = sw_ocean[-1]-st_ocean[-1]

    ################################################

    # save u:

    # Name the file that contains u,v, and w
    out_filename_uvw = path2save+'ocean.'+str(year).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+'.uvw.nc'
    print 'Saving calculation in: '+out_filename_uvw

    if os.path.lexists(out_filename_uvw):
        os.remove(out_filename_uvw)

    netcdf_file = nc.Dataset(out_filename_uvw, 'w', format='NETCDF4')
    netcdf_file.description = exp

    # Assign dimensions of the netCDF file. Add one extra depth level because we are interpolating the bottom two sw_ocean levels onto the bottom st_ocean which may not be 0 anymore. So, we add an extra level with negative this non-zero value so that it asserts W-vel = 0 at the sw_ocean location (see figures from LateX technical document)
    netcdf_file.createDimension('time', None)
    netcdf_file.createDimension('st_ocean', len(st_ocean)+1)
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

    # u velocity (crucial that the fill_value for the mask is 0 to satisfy the no-slip boundary condition)
    u_var = netcdf_file.createVariable('u', 'f4', ('time','st_ocean','yu_ocean','xu_ocean'),fill_value = 0)
    u_var.units = 'm/s'
    u_var.long_name = 'zonal velocity'
    u_var.missing_value = 0

    # Add a bottom layer that will contain a boundary so that CMS has 0 velocities to interpolate with
    # Take away mask, add extra level (all masked), then reapply the mask.
    u_null_mask = u
    u_null_mask[umask == True] = -1.e20
    u_null_mask = np.insert(u_null_mask,len(st_ocean)-1,-1.e20,axis=-3)
    u_regridded = np.ma.masked_array(u_null_mask, mask=(u_null_mask == -1.e20))

    # v velocity
    v_var = netcdf_file.createVariable('v', 'f4', ('time','st_ocean','yu_ocean','xu_ocean'),fill_value = 0)
    v_var.units = 'm/s'
    v_var.long_name = 'meridional velocity'
    v_var.missing_value = 0

    # vmask (crucial that the fill_value for the mask is 0 to satisfy the no-slip boundary condition)
    vmask = umask

    # Add a bottom layer that will contain a boundary so that CMS has 0 velocities to interpolate with
    # Take away mask, add extra level (all masked), then reapply the mask.
    v_null_mask = v
    v_null_mask[vmask == True] = -1.e20
    v_null_mask = np.insert(v_null_mask,len(st_ocean)-1,-1.e20,axis=-3)
    v_regridded = np.ma.masked_array(v_null_mask, mask=(v_null_mask == -1.e20))

    # w velocity
    w_var = netcdf_file.createVariable('wt', 'f4', ('time','st_ocean','yu_ocean','xu_ocean'),fill_value = 0)
    w_var.units = 'm/s'
    w_var.long_name = 'vertical velocity'
    w_var.missing_value = 0
    
    n_xt = len(xt_ocean)
    n_yt = len(yt_ocean)
    n_sw = len(sw_ocean)
    
    # Replace mask with large values so interpolation gives large values if one of the 4 points is land [(large + 0 + 0 + 0)/4 = still large is the idea here]
    w_null_mask = w[:,:,:]
    w_null_mask[wmask == True] = -1.e20
    
    # Horizontal Interpolation
    w_sw = w_null_mask[:,:-1,:]
    w_nw = w_null_mask[:,1:,:]
    w_se = np.roll(w_sw,-1,axis=-1)
    w_ne = np.roll(w_nw,-1,axis=-1)
    w2d = 0.25*(w_sw+w_nw+w_se+w_ne)
    # Reapply the mask to any point where one of the points interpolated over was masked
    w2d[w2d < -200] = -1.e20
    
    # Vertical Interpolation
    w_regridded = np.zeros((n_sw,n_yt,n_xt))
    # Top layer (Ekman?) intepolation - I took this from CMS. Should not matter since I am not looking at particles in the top layer anyways.
    w_regridded[0,:,:] = w2d[0,:,:]*(st_ocean[0]/sw_ocean[0])
    # Tile st_ocean and sw_ocean onto all lat/lon values so that it has same dimensions as the w-field
    st_tile = np.swapaxes(np.tile(st_ocean,(n_xt,n_yt,1)),0,-1)
    sw_tile = np.swapaxes(np.tile(sw_ocean,(n_xt,n_yt,1)),0,-1)
    # Apply vertical interpolation
    w_regridded[1:,:,:] = w2d[:-1,:,:]+(w2d[1:,:,:]-w2d[:-1,:,:])*((st_tile[1:]-sw_tile[:-1])/(sw_tile[1:]-sw_tile[:-1]))
    # Reapply the mask if the layer below was masked.
    w_regridded[w_regridded < -200] = -1.e20
    
    # Insert layer at the bottom (because of trilinear interpolation method in CMS)
    w_regridded = np.insert(w_regridded,len(st_ocean)-1,-1.e20,axis=-3)
    
    # Make the cell below bottom cells have negative the velocities of the cell above. Then, halfway between the two (approx the sw_ocean bottom of the domain) will have W-vel = 0.
    for i in range(0,n_xt):
        for j in range(0,n_yt):
            b_ind = 0
            for k in range(0,n_sw):
                b_ind = k
                if (w_regridded[k,j,i] < -200):
                    if (k == 0):
                        break
                    else:
                        w_regridded[k,j,i] = -w_regridded[k-1,j,i]
                        break

    # Reapply mask
    w_regridded = np.ma.masked_array(w_regridded, mask=(w_regridded < -1.e18))

    # data
    ti[:] = time
    st[:] = np.append(st_ocean[:],[st_ocean[-1]+(diff_bottom*2.0)],axis=0)
    yu[:] = yu_ocean
    xu[:] = xu_ocean

    # Only w is regridded but the landmasks were changed to 0 for all three fields.
    u_var[0,:,:,:] = u_regridded
    v_var[0,:,:,:] = v_regridded
    w_var[0,:,:,:] = w_regridded

    netcdf_file.close()

    ####################################################################