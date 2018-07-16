#!/usr/bin/env python

# import modules
import os,sys,glob
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.cm as cm
import matplotlib.colorbar as cb

model = 'CM2.6'
exp = 'CM2.6_A_Control-1860_V03/'
exp_short_name = '1860'
archdir = '/archive/wga/CM2.6/CM2.6_A_Control-1860_V03/pp/ocean_trans/av/annual_1yr/'

#model = 'CM2.5'
#exp = 'CM2.5_A_Control-1860_Y03/'
#exp_short_name = '1860'
#archdir = '/archive/wga/Siena/siena_201303_fix_rds-whit/'+'CM2.5_A_Control-1860_Y03-MBLING-BLING/pp/ocean_trans/ts/annual/20yr/'

#model = 'CM2-1deg'
#exp = 'CM_O1p0_C180_A02/'
#exp_short_name = '1860'
#archdir = '/archive/rds/Siena/siena_201303_fix_rds/'+'CM_O1p0_C180_A02-topaz-bling-minibling-ctrl-restart_bgc/pp/ocean_trans/av/annual_1yr/'

first_year = 192
last_year = 192

griddir = '/archive/Carolina.Dufour/Mosaic/'+model+'/'
figdir = '/work/Henri.Drake/figures/overturning/'+model+'/'+exp

grid_file = griddir+model+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
dxt = gridFile.variables['dxt'][...]
dyt = gridFile.variables['dyt'][...]
wet = gridFile.variables['wet'][...]
lat_t = gridFile.variables['geolat_t'][...]
lat_av = np.mean(lat_t,axis=1)
long_t = gridFile.variables['geolon_t'][...]

if model == 'CM2.6':
	ty_trans_rho = 0
	for year in range(first_year,last_year+1):
		#data_file = archdir + str(year).zfill(4)+'0101.ocean.nc'
		data_file = archdir + 'ocean_trans.'+str(year).zfill(4)+'0101.nc'
		print 'Opening '+data_file
		ncFile = nc.Dataset(data_file)
		ty_trans_rho = ty_trans_rho + np.mean(ncFile.variables['ty_trans_rho'][...],
			axis=0)
	ty_trans_rho = ty_trans_rho / (last_year-first_year + 1)
elif model == 'CM2.5':
	data_file = archdir + 'ocean_trans.0181-0200.ty_trans_rho.nc'
	print 'Opening '+data_file
	ncFile = nc.Dataset(data_file)
	ty_trans_rho = np.mean(ncFile.variables['ty_trans_rho'][first_year-181:last_year-181+1,...],axis=0)
elif model == 'CM2-1deg':
	ty_trans_rho = 0
	ty_trans_rho_gm = 0
	for year in range(first_year,last_year+1):
		#data_file = archdir + str(year).zfill(4)+'0101.ocean.nc'
		data_file = archdir + 'ocean_trans.'+str(year).zfill(4)+'-'+str(year).zfill(4)+\
			'.ty_trans_rho.nc'
		print 'Opening '+data_file
		ncFile = nc.Dataset(data_file)
		ty_trans_rho = ty_trans_rho + ncFile.variables['ty_trans_rho'][0,...]
		data_file = archdir + 'ocean_trans.'+str(year).zfill(4)+'-'+str(year).zfill(4)+\
			'.ty_trans_rho_gm.nc'
		print 'Opening '+data_file
		ncFile = nc.Dataset(data_file)
		ty_trans_rho_gm = ty_trans_rho_gm + ncFile.variables['ty_trans_rho_gm'][0,...]
	ty_trans_rho = ty_trans_rho / (last_year-first_year + 1)
	ty_trans_rho_gm = ty_trans_rho_gm / (last_year-first_year + 1)

potrho = ncFile.variables['potrho'][...]
# zonal sum:
ty_trans_rho = np.sum(ty_trans_rho,axis=2)
if model == 'CM2-1deg':
	ty_trans_rho_gm = np.sum(ty_trans_rho_gm,axis=2)
	ty_trans_rho_gm2 = np.zeros_like(ty_trans_rho_gm)
	ty_trans_rho_gm2[1:,:] = -np.diff(ty_trans_rho_gm,axis=0)
	net_trans = ty_trans_rho - ty_trans_rho_gm2
	trans_30 = net_trans[potrho>1036,np.where(lat_t[:,0]>-30)[0][0]]
else:
	trans_30 = ty_trans_rho[potrho>1036,np.where(lat_t[:,0]>-30)[0][0]]
net_southward = np.sum(trans_30[trans_30<0])
if model == 'CM2-1deg':
	overturning = np.flipud(np.cumsum(np.flipud(ty_trans_rho),axis=0)) - ty_trans_rho_gm
else:
	overturning = np.flipud(np.cumsum(np.flipud(ty_trans_rho),axis=0))
	
darkblue = '#000080'
blue = '#0000ff'
blue0 = '#006cd8'
lightblue2 = '#0080ff'
lightblue1 = '#76ffff'
yellow = '#ffff00'
yellow2 = '#fdd017'
orange = '#ff8000'
red = '#ff0000' 
darkred = '#670000'	
my_cmap = col.LinearSegmentedColormap.from_list('own2',
	[darkblue,blue,lightblue2,lightblue1,'w',yellow,orange,red,darkred])
cm.register_cmap(cmap=my_cmap)

contours=np.arange(-30,31,.2)
norm1 = col.Normalize(vmin=-20,vmax=20)
font = {'size':14}
tick_font=14

fig1=plt.figure(1)
plt.clf()
ax1 = fig1.add_axes([0.12,0.1,0.72,0.8])
plt.contourf(lat_av,potrho,-overturning,contours,cmap=my_cmap,norm=norm1)
plt.contour(lat_av,potrho,-overturning,
	np.append(np.arange(-50,0,2),np.arange(2,50,2)),colors='black',linewidths=1)
plt.ylim((1037.5,1034))
plt.yticks(np.arange(0,4)+1034,('1034','1035','1036','1037'))
plt.ylabel('Potential density (kgm$^{-3}$)',font)
plt.xlabel('Latitude',font)
plt.xlim((-70,-30))
plt.xticks([-70,-60,-50,-40,-30])
plt.legend(loc=2,prop=font)
plt.title('Residual overturning for '+model+', '+exp_short_name+', yrs '+\
	str(first_year)+'-'+str(last_year),font)
plt.tick_params(labelsize=tick_font)
axcb = fig1.add_axes([0.88,0.1,0.03,0.8])
cbar=cb.ColorbarBase(axcb,cmap=my_cmap,norm=norm1,orientation='vertical')
cbar.set_label('overturning (Sv)')
cbar_labels=plt.getp(cbar.ax.axes,'yticklabels')
plt.setp(cbar_labels,fontsize=tick_font)
import matplotlib
matplotlib.rcParams.update({'font.size':14})
plt.clim(vmin=-40,vmax=40) 
plt.text(-24,0.9,'Net southward = '+'%.1f' % net_southward + ' Sv')
plt.show()
outfile = os.path.join(figdir,model+'_'+exp_short_name+'_SO_psires_'+str(first_year)+\
	'-'+str(last_year)+'.png')
plt.savefig(outfile,bbox_inches='tight',dpi=200)
