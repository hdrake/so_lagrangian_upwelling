#!/usr/bin/env python
# import modules
import os,sys,glob
import numpy as np
import netCDF4 as nc

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.cm as cm
import matplotlib.colorbar as cb
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator, LinearLocator, NullLocator
from jdcal import gcal2jd
from jdcal import jd2gcal

##############

model = 'CM2.5'
exp = 'CM2.5_A_Control-1860_Y03-MBLING-BLING_rerun'
path2files = '/archive/wga/Siena/siena_201303_fix_rds-whit/CM2.5_A_Control-1860_Y03-MBLING-BLING_rerun/history/'
path2save = '/archive/hfd/CM2.5/CM2.5_A_Control-1860_Y03-MBLING-BLING_rerun/5day/for_CMS_SO/'
path2figures = '/work/hfd/figures/CM2.5/misc/'
figdir = path2figures

# Get grid Data
griddir = '/archive/Carolina.Dufour/Mosaic/'+model+'/'
grid_file = griddir+model+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
lat_t = gridFile.variables['geolat_t'][...]
lon_t = gridFile.variables['geolon_t'][...]
wet = gridFile.variables['wet'][...]

# Set parameters
first_year = 190
last_year = 200
boundary_north = -29.5
index_north = np.where(lat_t[:,0]<boundary_north)[0][-1]

# Make copies of lon thickness and lat thickness of t gridcells for each depth to make calculations easier
dxt = gridFile.variables['dxt'][:index_north,:]
dyt = gridFile.variables['dyt'][:index_north,:]
area_t = dxt*dyt
dxu = gridFile.variables['dxu'][:index_north,:]
dyu = gridFile.variables['dyu'][:index_north,:]
ht = gridFile.variables['ht'][:index_north,:]

# The 5-day field files
files = os.listdir(path2files)
files = np.sort(files)

year = 190
count_day = 1

## Online vhrho_nt and uhrho_et
budget_file_name = path2files+str(year).zfill(4)+'0101.ocean_budgets_5d.nc'
data = nc.Dataset(budget_file_name)
vhrho_nt = data.variables['vhrho_nt'][count_day,:,:index_north,:]
sw_ocean = data.variables['sw_ocean'][...]
vhrho_nt = np.ma.filled(vhrho_nt,fill_value=0)
uhrho_et = data.variables['uhrho_et'][count_day,:,:index_north,:]
uhrho_et = np.ma.filled(uhrho_et,fill_value=0)

# rate of change of volume from 5 daily rho_dzt:
data = nc.Dataset(path2files+str(year).zfill(4)+'0101.ocean_bgc_physics_field_rho_dzt.nc')
rho_dzt = data.variables['rho_dzt'][count_day-1:count_day+2,:,:index_north,:]
rho_dzt_dt = (rho_dzt[2,...]-rho_dzt[0,...])/(10.*24.*60.*60.)

### try using daily SSH instead (only works where no ice):
data = nc.Dataset(path2files+str(year).zfill(4)+'0101.ice_daily.nc')
# get from 1 day before to 1 day after this 5 day period:
SSH = data.variables['SSH'][count_day*5-1:count_day*5+6,:index_north,:]
# check where there is ice:
EXT_max = np.max(data.variables['EXT'][count_day*5-1:count_day*5+6,:index_north,:],axis=0)
EXT_max = np.tile(EXT_max,(len(sw_ocean),1,1))

# make thickness array with partial bottom cells (this will be same for all timesteps):
thickness = np.diff(np.append(0,sw_ocean))
thickness_tile = np.swapaxes(np.tile(thickness,(ht.shape[1],ht.shape[0],1)),0,2)
depth = np.cumsum(thickness)
for ii in range(ht.shape[1]):
	for jj in range(ht.shape[0]):
		last_depth_index = np.where(depth>ht[jj,ii])[0]
		if len(last_depth_index) == 0:
			continue
		else:
			thickness_tile[last_depth_index[0],jj,ii] += (ht[jj,ii] - \
				depth[last_depth_index[0]])
# mask:
thickness_tile = np.ma.array(thickness_tile,mask=rho_dzt[0,...].mask)

# modify thickness with SSH:
# average SSH at start of 5 day period:
SSH1 = np.mean(SSH[:2,...],axis=0)
# average SSH at end of 5 day period:
SSH5 = np.mean(SSH[-2:,...],axis=0)
rho_dzt_dt_offline = ((1 + SSH5/ht) - (1 + SSH1/ht))*\
	thickness_tile/(5*24.*60.*60.)*1035.

# under ice use rho_dzt_dt from 5 day averages, not daily averages from ice fields:
rho_dzt_dt_combined = np.ma.where(EXT_max==0.,rho_dzt_dt_offline,rho_dzt_dt)

print 'Calculating w from online fluxes'

lat_fluxes = (uhrho_et*dyt-np.roll(uhrho_et*dyt,1,axis=-1))/area_t + \
	(vhrho_nt*dxu-np.insert((vhrho_nt*dxu)[:,:-1,:],0,(vhrho_nt*dxu)[:,0,:],axis=1))/area_t

wrho_bt = np.zeros(vhrho_nt.shape)
for k in range(49)[::-1]:
	wrho_bt[k,...] = -lat_fluxes[k+1,...] - rho_dzt_dt_combined[k+1,...] + wrho_bt[k+1,...]

w = np.ma.array(wrho_bt/1035.,mask=(rho_dzt[0,...]).mask)

## Online w
w_file_name = path2files+str(year).zfill(4)+'0101.ocean_bgc_physics_field_w.nc'
data = nc.Dataset(w_file_name)
online_w = data.variables['w'][count_day,:,:index_north,:]

# plot:
layer = -45

plt.figure(6)
plt.clf()
plt.pcolormesh(w[layer,...]-online_w[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-1e-7,1e-7))
plt.title('diff, k='+str(layer))
plt.draw()
plt.show()

plt.figure(10)
plt.clf()
plt.pcolormesh(w[layer,...]-OG_w[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-1e-5,1e-5))
plt.title('diff, k='+str(layer))
plt.draw()
plt.show()

plt.figure(11)
plt.clf()
plt.pcolormesh(w[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-1e-5,1e-5))
plt.title('online_w, k='+str(layer))
plt.draw()
plt.show()

plt.figure(12)
plt.clf()
plt.pcolormesh(w[layer,...]-w_5d[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-1e-7,1e-7))
plt.title('diff w and online-tendency w, k='+str(layer))
plt.draw()
plt.show()

plt.figure(13)
plt.clf()
plt.pcolormesh(w[layer,...]-w_10d[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-1e-7,1e-7))
plt.title('diff w and 10-day-tendency w, k='+str(layer))
plt.draw()
plt.show()
    outfile = figdir+'diff_uhrho_et'+str(kloop)+'.png'
    plt.savefig(outfile,bbox_inches = 'tight',dpi=270)



# Compare rho_dzt_tendencies
plt.figure(7)
plt.clf()
plt.pcolormesh(rho_dzt_tendency[layer,...]-rho_dzt_tendency_10d[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-1e-6,1e-6))
plt.title('diff, k='+str(layer))
plt.draw()
plt.show()

# vhrho diff

plt.figure(8)
plt.clf()
plt.pcolormesh(vhrho_nt[layer,...]-online_vhrho_nt[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-0.0001,0.0001))
plt.title('diff vhrho_nt, k='+str(layer))
plt.draw()
plt.show()

# vhrho

plt.figure(9)
plt.clf()
plt.pcolormesh(vhrho_nt[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-1,1))
plt.title('vhrho_nt, k='+str(layer))
plt.draw()
plt.show()

plt.figure(9)
plt.clf()
plt.pcolormesh(w[layer,...]-w2[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-1e-7,1e-7))
plt.title('diff, k='+str(layer))
plt.draw()
plt.show()



sys.exit()

plt.figure(3)
plt.clf()
plt.pcolormesh(w[layer,...]-online_w[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-1e-7,1e-7))
plt.title('diff, k='+str(layer))
plt.draw()
plt.show()


plt.figure(3)
plt.clf()
plt.pcolormesh(w1[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-1e-5,1e-5))
plt.title('w2, k='+str(layer))
plt.draw()
plt.show()

plt.figure(2)
plt.clf()
plt.pcolormesh(online_w[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-1e-5,1e-5))
plt.title('online, k='+str(layer))
plt.draw()
plt.show()


sys.exit()

plt.figure(4)
plt.clf()
plt.pcolormesh((w[layer,...]-online_w[layer,...])/online_w[layer,...],cmap='RdBu')
plt.colorbar()
plt.clim((-.5,.5))
plt.title('error, k='+str(layer))
plt.draw()
plt.show()




