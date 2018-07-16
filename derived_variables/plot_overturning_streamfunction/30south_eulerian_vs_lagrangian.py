

# import modules
import os,sys,glob
import numpy as np
import netCDF4 as nc

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator, LinearLocator, NullLocator

from jdcal import jd2gcal
np.set_printoptions(threshold=np.nan)

#exp = 'CM2.5_A_Control-1860_Y03-MBLING-BLING/'
#exp = 'CM2.6_A_Control-1860_V03/'
exp = 'CM2-1deg_A_Control-1860_V03/'

tempres = '5day'

if exp == 'CM2.6_A_Control-1860_V03/':
    model = 'CM2.6'
if exp == 'CM2.6_A_Control-1990_V03/':
    model = 'CM2.6'
if exp == 'CM2.5_A_Control-1860_Y03-MBLING-BLING/':
    model = 'CM2.5'
if exp == 'CM2-1deg_A_Control-1860_V03/':
    model = 'CM2-1deg'

path_to_velocities = '/archive/hfd/'+model+'/'+exp+'/'+tempres+'/for_CMS_SO'+'_Edited/'
path_to_rho_dzt = '/archive/hfd/'+model+'/'+exp+'/'+tempres+'/rho_dzt/'

# Latitude we want to release at
release_lat = -30.001

# Get dxu and dimensions from CM2.6 because this is our release grid for all experiments
rho_0 = 1035.
griddir = '/archive/Carolina.Dufour/Mosaic/CM2.6/CM2.6'
grid_file = griddir+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
dxu = gridFile.variables['dxu'][...]
gridFile.close()
CM26_file = '/archive/hfd/CM2.6/CM2.6_A_Control-1860_V03/5day/for_CMS_SO_Edited/ocean.01810107.uvw.nc'
ncFile = nc.Dataset(CM26_file)
st_ocean26 = ncFile.variables['st_ocean'][:-1]
xu_ocean26 = ncFile.variables['xu_ocean'][...]
yu_ocean26 = ncFile.variables['yu_ocean'][...]
latind26 = (np.abs(release_lat - yu_ocean26)).argmin()
dxu = dxu[latind26,0]
umask = ncFile.variables['u'][0,:-1,latind26,:].mask
ncFile.close()

# Get dxu and dimensions from CM2-1deg because this is our release grid for all experiments
rho_0 = 1035.
griddir = '/archive/Carolina.Dufour/Mosaic/CM2-1deg/CM2-1deg'
grid_file = griddir+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
dxu_cm21 = gridFile.variables['dxu'][...]
gridFile.close()
CM21_file = '/archive/hfd/CM2-1deg/CM2-1deg_A_Control-1860_V03/5day/for_CMS_SO_Edited/ocean.01810107.uvw.nc'
ncFile = nc.Dataset(CM21_file)
st_ocean21 = ncFile.variables['st_ocean'][:-1]
xu_ocean21 = ncFile.variables['xu_ocean'][...]
yu_ocean21 = ncFile.variables['yu_ocean'][...]
latind21 = (np.abs(release_lat - yu_ocean21)).argmin()
dxu_cm21 = dxu_cm21[latind21,0]
umask = ncFile.variables['u'][0,:-1,latind21,:].mask
ncFile.close()

# Release at a fixed latitude and between these two longitudes
files = os.listdir(path_to_velocities)
files = np.sort(files)

count = 1
countbad = 1

# Loop through each velocity field and release at monthly interval
# This is different for 5-day and monthly fields
file = files[86]

nx_cm26 = np.size(xu_ocean26)
ny_cm26 = np.size(yu_ocean26)
nx_cm21 = np.size(xu_ocean21)
ny_cm21 = np.size(yu_ocean21)

year = int(file[6:10])
month = int(file[10:12])
day = int(file[12:14])

if (day <= 15) and (day > 10):
    ncFile = nc.Dataset(path_to_velocities+file)
    yu_ocean = ncFile.variables['yu_ocean'][...]
    xu_ocean = ncFile.variables['xu_ocean'][...]
    st_ocean = ncFile.variables['st_ocean'][:-1]
    
    latind = (np.abs(release_lat - yu_ocean)).argmin()
    vlower = ncFile.variables['v'][0,:-1,latind,:]
    ncFile.close()
    ncFile = nc.Dataset(path_to_velocities+file[0:12]+str(day+5).zfill(2)+'.uvw.nc')
    vupper = ncFile.variables['v'][0,:-1,latind,:]
    ncFile.close()
    
    ncFile = nc.Dataset(path_to_rho_dzt+'nest_1_'+file[6:14]+'000000.rho_dzt.nc')
    yt_ocean = ncFile.variables['yt_ocean'][...]
    xt_ocean = ncFile.variables['xt_ocean'][...]
    latind_t = (np.abs(release_lat - yt_ocean)).argmin()
    rho_dzt_lower = ncFile.variables['rho_dzt'][0,:,latind_t,:]
    ncFile.close()
    ncFile = nc.Dataset(path_to_rho_dzt+'nest_1_'+file[6:12]+str(day+5).zfill(2)+'000000.rho_dzt.nc')
    rho_dzt_upper = ncFile.variables['rho_dzt'][0,:,latind_t,:]
    ncFile.close()
    
    
    v_cm21 = vlower + (vupper - vlower)*((15.-day)/(5.))
    rho_dzt_cm21 = rho_dzt_lower + (rho_dzt_upper - rho_dzt_lower)*((15.-day)/(5.))
    transport_cm21 = v_cm21*dxu_cm21*(rho_dzt_cm21/1035.)/(1.e6)
    
    vinterp_lower = np.zeros((len(st_ocean26),len(xu_ocean26)))
    vinterp_upper = np.zeros((len(st_ocean26),len(xu_ocean26)))
    rho_dzt_interp_lower = np.zeros((len(st_ocean26),len(xu_ocean26)))
    rho_dzt_interp_upper = np.zeros((len(st_ocean26),len(xu_ocean26)))
    for q in range(len(st_ocean26)):
        for j in range(len(xu_ocean26)):
            # Get longitude indices
            xu0 = (np.abs(xu_ocean-xu_ocean26[j])).argmin()
            if xu_ocean[xu0] >= xu_ocean26[j]:
                if xu0 == 0:
                    xu1 = np.size(xu_ocean)-1
                else:
                    xu1 = xu0-1
            else:
                if xu0 == np.size(xu_ocean)-1:
                    xu1 = 0
                else:
                    xu1 = xu0+1
            
            # Get longitude indices
            xt0 = (np.abs(xt_ocean-xu_ocean26[j])).argmin()
            if xt_ocean[xt0] >= xu_ocean26[j]:
                if xt0 == 0:
                    xt1 = np.size(xt_ocean)-1
                else:
                    xt1 = xu0-1
            else:
                if xt0 == np.size(xt_ocean)-1:
                    xt1 = 0
                else:
                    xt1 = xt0+1
            # Interpolate velocities and rho_dzt in space
            vinterp_upper[q,j] = vupper[q,xu0]+(vupper[q,xu1]-vupper[q,xu0])*((xu_ocean26[j]-xu_ocean[xu0])/(xu_ocean[xu1]-xu_ocean[xu0]))
            vinterp_lower[q,j] = vlower[q,xu0]+(vlower[q,xu1]-vlower[q,xu0])*((xu_ocean26[j]-xu_ocean[xu0])/(xu_ocean[xu1]-xu_ocean[xu0]))
            rho_dzt_interp_lower[q,j] = rho_dzt_lower[q,xt0]+(rho_dzt_lower[q,xt1]-rho_dzt_lower[q,xt0])*((xu_ocean26[j]-xt_ocean[xt0])/(xt_ocean[xt1]-xt_ocean[xt0]))
            rho_dzt_interp_upper[q,j] = rho_dzt_upper[q,xt0]+(rho_dzt_upper[q,xt1]-rho_dzt_upper[q,xt0])*((xu_ocean26[j]-xt_ocean[xt0])/(xt_ocean[xt1]-xt_ocean[xt0]))
    
    ######### NEED INTERPOLATION SCHEME #############
    print 'Now doing file '+file
    vinterp = vinterp_lower + (vinterp_upper - vinterp_lower)*((15.-day)/(5.))
    rho_dzt_interp = rho_dzt_interp_lower + (rho_dzt_interp_upper - rho_dzt_interp_lower)*((15.-day)/(5.))
    transport_cm26 = vinterp*dxu*(rho_dzt_interp/1035.)/(1.e6)

plt.figure(1)
plt.clf()
cm21_mask = np.isnan(transport_cm21)
cm26_mask = np.isnan(transport_cm26)
transport_cm21 = np.ma.masked_array(transport_cm21,mask = cm21_mask)
transport_cm26 = np.ma.masked_array(transport_cm26,mask = cm26_mask)
cm21_mask = transport_cm21.mask
cm26_mask = transport_cm26.mask
transport_cm21 = np.ma.masked_array(transport_cm21,mask = cm21_mask)
transport_cm26 = np.ma.masked_array(transport_cm26,mask = cm26_mask)
plt.plot(np.sum(transport_cm21,axis = 1),st_ocean21,'r')
plt.plot(np.sum(transport_cm26,axis = 1),st_ocean26,'b')
plt.plot([-5,4],[1000,1000],'c--')
plt.plot([-5,4],[3500,3500],'c--')
plt.plot([0,0],[st_ocean26[0],st_ocean26[-1]])
plt.gca().invert_yaxis()
plt.suptitle('Zonal sum of residual transport at 30 south, derived from CM2-1deg velocities in red \nand same velocities interpolated onto the CM2.6 grid in blue. The sum between \nthe cyan lines is -10.0 for interpolated fields and -13.2 for CM2-1deg fields.')
plt.ylabel('Depth (m)')
plt.xlabel('Zonal sum of transport on each model depth level (Sv)')
plt.draw()
plt.show()
plt.savefig('/work/hfd/figures/velocity_interpolation.png')

plt.figure(2)
plt.clf()
transport_cm21_neg = np.copy(transport_cm21)
transport_cm26_neg = np.copy(transport_cm26)
transport_cm21_neg = np.ma.masked_array(transport_cm21_neg,mask = cm21_mask)
transport_cm26_neg = np.ma.masked_array(transport_cm26_neg,mask = cm26_mask)
transport_cm21_neg[transport_cm21_neg > 0] = 0
transport_cm26_neg[transport_cm26_neg > 0] = 0
plt.plot(np.sum(transport_cm21_neg,axis = 1),st_ocean21,'r')
plt.plot(np.sum(transport_cm26_neg,axis = 1),st_ocean26,'b')
plt.plot([-5,4],[1000,1000],'c--')
plt.plot([-5,4],[3500,3500],'c--')
plt.plot([0,0],[st_ocean26[0],st_ocean26[-1]],'k')
plt.gca().invert_yaxis()
plt.suptitle('Zonal sum of southward transport at 30 south, derived from CM2-1deg velocities in red \nand same velocities interpolated onto the CM2.6 grid in blue. The sum between \nthe cyan lines is -70.9 for interpolated fields and -81.6 for CM2-1deg fields.')
plt.ylabel('Depth (m)')
plt.xlabel('Zonal sum of southward transport on each model depth level (Sv)')
plt.draw()
plt.show()
plt.savefig('/work/hfd/figures/southward_velocity_interpolation.png')

plt.figure(1)
plt.clf()
cm21_mask = np.isnan(transport_cm21)
cm26_mask = np.isnan(transport_cm26)
transport_cm21 = np.ma.masked_array(transport_cm21,mask = cm21_mask)
transport_cm26 = np.ma.masked_array(transport_cm26,mask = cm26_mask)
cm21_mask = transport_cm21.mask
cm26_mask = transport_cm26.mask
transport_cm21 = np.ma.masked_array(transport_cm21,mask = cm21_mask)
transport_cm26 = np.ma.masked_array(transport_cm26,mask = cm26_mask)
plt.plot(np.sum(transport_cm21,axis = 1),st_ocean21,'r')
plt.plot(np.sum(transport_cm26,axis = 1),st_ocean26,'b')
plt.plot([-5,4],[1000,1000],'c--')
plt.plot([-5,4],[3500,3500],'c--')
plt.plot([0,0],[st_ocean26[0],st_ocean26[-1]])
plt.gca().invert_yaxis()
plt.suptitle('Zonal sum of residual transport at 30 south, derived from CM2-1deg velocities in red \nand same velocities interpolated onto the CM2.6 grid in blue. The sum between \nthe cyan lines is -10.0 for interpolated fields and -13.2 for CM2-1deg fields.')
plt.ylabel('Depth (m)')
plt.xlabel('Zonal sum of transport on each model depth level (Sv)')
plt.draw()
plt.show()
plt.savefig('/work/hfd/figures/velocity_interpolation.png')