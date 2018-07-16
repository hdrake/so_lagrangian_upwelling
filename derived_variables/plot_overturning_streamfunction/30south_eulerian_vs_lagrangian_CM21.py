

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
dxu_high = gridFile.variables['dxu'][...]
gridFile.close()
CM_high_file = '/archive/hfd/CM2.6/CM2.6_A_Control-1860_V03/5day/for_CMS_SO_Edited/ocean.01810107.uvw.nc'
ncFile = nc.Dataset(CM_high_file)
st_ocean_high = ncFile.variables['st_ocean'][:-1]
xu_ocean_high = ncFile.variables['xu_ocean'][...]
yu_ocean_high = ncFile.variables['yu_ocean'][...]
latind_high = (np.abs(release_lat - yu_ocean_high)).argmin()
dxu_high = dxu_high[latind_high,0]
umask = ncFile.variables['u'][0,:-1,latind_high,:].mask
ncFile.close()

# Get dxu and dimensions from CM2-1deg because this is our release grid for all experiments
rho_0 = 1035.
griddir = '/archive/Carolina.Dufour/Mosaic/CM2-1deg/CM2-1deg'
grid_file = griddir+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
dxu_low = gridFile.variables['dxu'][...]
gridFile.close()
CM_low_file = '/archive/hfd/CM2-1deg/CM2-1deg_A_Control-1860_V03/5day/for_CMS_SO_Edited/ocean.01810107.uvw.nc'
ncFile = nc.Dataset(CM_low_file)
st_ocean_low = ncFile.variables['st_ocean'][:-1]
xu_ocean_low = ncFile.variables['xu_ocean'][...]
yu_ocean_low = ncFile.variables['yu_ocean'][...]
latind_low = (np.abs(release_lat - yu_ocean_low)).argmin()
dxu_low = dxu_low[latind_low,0]
umask_low = ncFile.variables['u'][0,:-1,latind_low,:].mask
ncFile.close()

# Release at a fixed latitude and between these two longitudes
files = os.listdir(path_to_velocities)
files = np.sort(files)

count = 1
countbad = 1

# Loop through each velocity field and release at monthly interval
# This is different for 5-day and monthly fields
file = files[86]

nx_cm_high = np.size(xu_ocean_high)
ny_cm_high = np.size(yu_ocean_high)
nx_cm_low = np.size(xu_ocean_low)
ny_cm_low = np.size(yu_ocean_low)

year = int(file[6:10])
month = int(file[10:12])
day = int(file[12:14])

zonal_sum_transport_cm_low = np.zeros(12*10)
zonal_sum_transport_cm_high = np.zeros(12*10)

count_times = 0
for file in files:
    if not file.endswith('uvw.nc'):
        continue
    year = int(file[6:10])
    month = int(file[10:12])
    day = int(file[12:14])
    if year <= 182:
        continue
    if year >= 193:
        continue
    if (day <= 15) and (day > 10):
        ncFile = nc.Dataset(path_to_velocities+file)
        yu_ocean = ncFile.variables['yu_ocean'][...]
        xu_ocean = ncFile.variables['xu_ocean'][...]
        st_ocean = ncFile.variables['st_ocean'][:-1]
        
        latind = (np.abs(release_lat - yu_ocean)).argmin()
        vlower = ncFile.variables['v'][0,:-1,latind,:].data
        ncFile.close()
        ncFile = nc.Dataset(path_to_velocities+file[0:12]+str(day+5).zfill(2)+'.uvw.nc')
        vupper = ncFile.variables['v'][0,:-1,latind,:].data
        ncFile.close()
        
        ncFile = nc.Dataset(path_to_rho_dzt+'nest_1_'+file[6:14]+'000000.rho_dzt.nc')
        yt_ocean = ncFile.variables['yt_ocean'][...]
        xt_ocean = ncFile.variables['xt_ocean'][...]
        latind_t = (np.abs(release_lat - yt_ocean)).argmin()
        rho_dzt_lower = ncFile.variables['rho_dzt'][0,:,latind_t,:].data
        rho_dzt_lower[rho_dzt_lower<-1.e-10] = 0
        ncFile.close()
        ncFile = nc.Dataset(path_to_rho_dzt+'nest_1_'+file[6:12]+str(day+5).zfill(2)+'000000.rho_dzt.nc')
        rho_dzt_upper = ncFile.variables['rho_dzt'][0,:,latind_t,:].data
        rho_dzt_upper[rho_dzt_upper<-1.e-10] = 0
        ncFile.close()
        
        
        v_cm_low = vlower + (vupper - vlower)*((15.-day)/(5.))
        rho_dzt_low = rho_dzt_lower + (rho_dzt_upper - rho_dzt_lower)*((15.-day)/(5.))
        transport_cm_low = v_cm_low*dxu_low*(rho_dzt_low/1035.)/(1.e6)
        
        vinterp_lower = np.zeros((len(st_ocean_high),len(xu_ocean_high)))
        vinterp_upper = np.zeros((len(st_ocean_high),len(xu_ocean_high)))
        rho_dzt_interp_lower = np.zeros((len(st_ocean_high),len(xu_ocean_high)))
        rho_dzt_interp_upper = np.zeros((len(st_ocean_high),len(xu_ocean_high)))
        for q in range(len(st_ocean_high)):
            for j in range(len(xu_ocean_high)):
                # Get longitude indices
                xu0 = (np.abs(xu_ocean-xu_ocean_high[j])).argmin()
                if xu_ocean[xu0] >= xu_ocean_high[j]:
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
                xt0 = (np.abs(xt_ocean-xu_ocean_high[j])).argmin()
                if xt_ocean[xt0] >= xu_ocean_high[j]:
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
                vinterp_upper[q,j] = vupper[q,xu0]+(vupper[q,xu1]-vupper[q,xu0])*((xu_ocean_high[j]-xu_ocean[xu0])/(xu_ocean[xu1]-xu_ocean[xu0]))
                vinterp_lower[q,j] = vlower[q,xu0]+(vlower[q,xu1]-vlower[q,xu0])*((xu_ocean_high[j]-xu_ocean[xu0])/(xu_ocean[xu1]-xu_ocean[xu0]))
                rho_dzt_interp_lower[q,j] = rho_dzt_lower[q,xt0]+(rho_dzt_lower[q,xt1]-rho_dzt_lower[q,xt0])*((xu_ocean_high[j]-xt_ocean[xt0])/(xt_ocean[xt1]-xt_ocean[xt0]))
                rho_dzt_interp_upper[q,j] = rho_dzt_upper[q,xt0]+(rho_dzt_upper[q,xt1]-rho_dzt_upper[q,xt0])*((xu_ocean_high[j]-xt_ocean[xt0])/(xt_ocean[xt1]-xt_ocean[xt0]))
        
        ######### NEED INTERPOLATION SCHEME #############
        print 'Now doing file '+file
        vinterp = vinterp_lower + (vinterp_upper - vinterp_lower)*((15.-day)/(5.))
        rho_dzt_interp = rho_dzt_interp_lower + (rho_dzt_interp_upper - rho_dzt_interp_lower)*((15.-day)/(5.))
        transport_cm_high = vinterp*dxu_high*(rho_dzt_interp/1035.)/(1.e6)
        zonal_sum_transport_cm_low[count_times] = np.sum(transport_cm_low[24:40,:])
        zonal_sum_transport_cm_high[count_times] = np.sum(transport_cm_high[24:40,:])
        count_times += 1

plt.figure(1)
plt.clf()
cm_low_mask = np.isnan(transport_cm_low)
cm_high_mask = np.isnan(transport_cm_high)

plt.plot(np.arange(count_times),zonal_sum_transport_cm_low,'r')
plt.plot(np.arange(count_times),zonal_sum_transport_cm_high,'b')
plt.plot([0,count_times],np.ones(2)*np.average(zonal_sum_transport_cm_low),'y--')
plt.plot([0,count_times],np.ones(2)*np.average(zonal_sum_transport_cm_high),'c--')
plt.gca().invert_yaxis()
plt.suptitle('Zonal sum of residual transport at 30 south, derived from CM2-1deg velocities in red \nand same velocities interpolated onto the CM2.6 grid in blue. The sum between \nthe cyan lines is '+str(np.average(zonal_sum_transport_cm_high))+' for interpolated fields and '+str(np.average(zonal_sum_transport_cm_low))+' for CM2-1deg fields.')
plt.ylabel('Zonally-and-depth-integrated transport at 30 south')
plt.xlabel('5-day mean timestep')
plt.draw()
plt.show()
plt.savefig('/work/hfd/figures/velocity_interpolation_time_series_CM2-1deg.png')
