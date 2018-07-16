#!/usr/bin/env python

# import modules
import os,sys,glob
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.cm as cm
import matplotlib.colorbar as cb

def interping(m1,m2,t1,t2,t):
    return m1 + (m2-m1) / (t2-t1) * (t-t1)

#model = 'CM2.6'
#exp = 'CM2.6_A_Control-1860_V03/'
#exp_short_name = '1860'
#archdir = '/archive/wga/CM2.6/CM2.6_A_Control-1860_V03/pp/ocean_trans/av/annual_1yr/'

#model = 'CM2.5'
#exp = 'CM2.5_A_Control-1860_Y03/'
#exp_short_name = '1860'
#archdir = '/archive/wga/Siena/siena_201303_fix_rds-whit/'+'CM2.5_A_Control-1860_Y03-MBLING-BLING/pp/ocean_trans/ts/annual/20yr/'

model = 'CM2-1deg'
exp = 'CM_O1p0_C180_A02/'
exp_short_name = '1860'
archdir = '/archive/rds/Siena/siena_201303_fix_rds/'+'CM_O1p0_C180_A02-topaz-bling-minibling-ctrl-restart_bgc/pp/ocean_trans/av/annual_1yr/'

first_year = 183
last_year = 192

griddir = '/archive/Carolina.Dufour/Mosaic/'+model+'/'
figdir = '/work/hfd/figures/overturning/'+model+'/'+exp

grid_file = griddir+model+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
dxt = gridFile.variables['dxt'][...]
dyt = gridFile.variables['dyt'][...]
wet = gridFile.variables['wet'][...]
lat_t = gridFile.variables['geolat_t'][...]
lat_av = np.mean(lat_t,axis=1)
long_t = gridFile.variables['geolon_t'][...]
long_av = np.mean(long_t,axis=0)

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
	ty_trans_rho = np.mean(ncFile.variables['ty_trans_rho'][1:13,...],axis=0)
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
        original_ty_trans_rho_gm = np.copy(ty_trans_rho_gm)

potrho = ncFile.variables['potrho'][...]
nrho = np.shape(ty_trans_rho)[0]
ny = np.shape(ty_trans_rho)[1]
nx = np.shape(ty_trans_rho)[2]
# How many more grid points you want in each the X and Y direction.
scale_x = 10
scale_y = 1

# Create high-res version of ty_trans_rho
high_ty_trans_rho = np.zeros((nrho,ny*scale_y,nx*scale_x))
# Fill the high-res version of ty_trans_rho with 2D-spatially-interpolated fields, dividing out by float(scale*scale) so that z
for j in range(ny):
    for i in range(nx):
        for jj in range(scale_y):
            for ii in range(scale_x):
                high_ty_trans_rho[:,j*scale_y+jj,i*scale_x+ii] = interping(interping(ty_trans_rho[:,j,i],ty_trans_rho[:,j,np.mod(i+1,nx)],
                                                        long_av[i],long_av[np.mod(i+1,nx)],
                                                        long_av[i]+(long_av[i]-long_av[np.mod(i+1,nx)])*ii/float(scale_x)),
                                              interping(ty_trans_rho[:,j,i],ty_trans_rho[:,j,np.mod(i+1,nx)],
                                                        long_av[i],long_av[np.mod(i+1,nx)],
                                                        long_av[i]+(long_av[i]-long_av[np.mod(i+1,nx)])*ii/float(scale_x)),
                                              lat_av[j],lat_av[np.mod(j+1,ny)],
                                              lat_av[j]+(lat_av[j]-lat_av[np.mod(j+1,ny)])*jj/float(scale_y))/float(scale_x*scale_y)

high_lat_av = np.zeros(ny*scale_y)
for j in range(ny):
    for jj in range(scale_y):
        high_lat_av[j*scale_y+jj] = lat_av[j]+(lat_av[j]-lat_av[np.mod(j+1,ny)])*jj/float(scale_y)

ty_trans_rho = np.sum(ty_trans_rho,axis=2)
high_ty_trans_rho = np.sum(high_ty_trans_rho,axis=2)

if model == 'CM2-1deg':
	ty_trans_rho_gm = np.sum(ty_trans_rho_gm,axis=2)
	ty_trans_rho_gm2 = np.zeros_like(ty_trans_rho_gm)
	ty_trans_rho_gm2[1:,:] = -np.diff(ty_trans_rho_gm,axis=0)
	net_trans = ty_trans_rho - ty_trans_rho_gm2
        trans_30 = net_trans[potrho>1036,np.where(lat_t[:,0]>-30)[0][0]]

        # Create high-res version of ty_trans_rho_gm
        high_ty_trans_rho_gm = np.zeros((nrho,ny*scale_y,nx*scale_x))
        # Fill the high-res version of ty_trans_rho_gm with 2D-spatially-interpolated fields, dividing out by float(scale*scale) so that z
        for j in range(ny):
            for i in range(nx):
                for jj in range(scale_y):
                    for ii in range(scale_x):
                        high_ty_trans_rho_gm[:,j*scale_y+jj,i*scale_x+ii] = interping(interping(original_ty_trans_rho_gm[:,j,i],original_ty_trans_rho_gm[:,j,np.mod(i+1,nx)],
                                                                                            long_av[i],long_av[np.mod(i+1,nx)],
                                                                                            long_av[i]+(long_av[i]-long_av[np.mod(i+1,nx)])*ii/float(scale_x)),
                                                                                  interping(original_ty_trans_rho_gm[:,j,i],original_ty_trans_rho_gm[:,j,np.mod(i+1,nx)],
                                                                                            long_av[i],long_av[np.mod(i+1,nx)],
                                                                                            long_av[i]+(long_av[i]-long_av[np.mod(i+1,nx)])*ii/float(scale_x)),
                                                                                  lat_av[j],lat_av[np.mod(j+1,ny)],
                                                                                      lat_av[j]+(lat_av[j]-lat_av[np.mod(j+1,ny)])*jj/float(scale_y))/float(scale_x*scale_y)

        high_lat_av = np.zeros(ny*scale_y)
        for j in range(ny):
            for jj in range(scale_y):
                high_lat_av[j*scale_y+jj] = lat_av[j]+(lat_av[j]-lat_av[np.mod(j+1,ny)])*jj/float(scale_y)

        high_ty_trans_rho_gm = np.sum(high_ty_trans_rho_gm,axis=2)
	high_ty_trans_rho_gm2 = np.zeros_like(high_ty_trans_rho_gm)
	high_ty_trans_rho_gm2[1:,:] = -np.diff(high_ty_trans_rho_gm,axis=0)
	high_net_trans = high_ty_trans_rho - high_ty_trans_rho_gm2
	high_trans_30 = high_net_trans[potrho>1036,np.where(high_lat_av[:]>-30)[0][0]]
else:
	trans_30 = ty_trans_rho[potrho>1036,np.where(lat_t[:,0]>-30)[0][0]]
        high_trans_30 = high_net_trans[potrho>1036,np.where(high_lat_av[:]>-30)[0][0]]
net_southward = np.sum(trans_30[trans_30<0])
high_net_southward = np.sum(high_trans_30[high_trans_30<0])
if model == 'CM2-1deg':
	overturning = np.flipud(np.cumsum(np.flipud(ty_trans_rho),axis=0)) - ty_trans_rho_gm
        high_overturning = np.flipud(np.cumsum(np.flipud(high_ty_trans_rho),axis=0)) - high_ty_trans_rho_gm
else:
	overturning = np.flipud(np.cumsum(np.flipud(ty_trans_rho),axis=0))
        high_overturning = np.flipud(np.cumsum(np.flipud(high_ty_trans_rho),axis=0))

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

fig1=plt.figure(2)
plt.clf()
ax1 = fig1.add_axes([0.12,0.1,0.72,0.8])
plt.contourf(high_lat_av,potrho,-high_overturning,contours,cmap=my_cmap,norm=norm1)
plt.contour(high_lat_av,potrho,-high_overturning,
            np.append(np.arange(-50,0,2),np.arange(2,50,2)),colors='black',linewidths=1)
plt.ylim((1037.5,1034))
plt.yticks(np.arange(0,4)+1034,('1034','1035','1036','1037'))
plt.ylabel('Potential density (kgm$^{-3}$)',font)
plt.xlabel('Latitude',font)
plt.xlim((-70,-30))
plt.xticks([-70,-60,-50,-40,-30])
plt.legend(loc=2,prop=font)
plt.title('Interpolated residual overturning for '+model+', '+exp_short_name+', yrs '+\
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
plt.text(-24,0.9,'Net southward = '+'%.1f' % high_net_southward + ' Sv')
plt.show()
outfile = os.path.join(figdir,model+'_'+exp_short_name+'_SO_interp_psires_'+str(first_year)+\
                       '-'+str(last_year)+'_10.png')
plt.savefig(outfile,bbox_inches='tight',dpi=200)
