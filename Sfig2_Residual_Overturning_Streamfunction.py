#!/usr/bin/env python

# import modules
import os,sys,glob
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.cm as cm
import matplotlib.colorbar as cb


models = ['CM2-1deg','CM2.5','CM2.6','CM2-1deg']
exps = ['CM_O1p0_C180_A02/','CM2.5_A_Control-1860_Y03/','CM2.6_A_Control-1860_V03/','CM_O1p0_C180_A02/']
exp_short_names = ['1860','1860','1860','1860']
archdirs = ['/archive/rds/Siena/siena_201303_fix_rds/CM_O1p0_C180_A02-topaz-bling-minibling-ctrl-restart_bgc/pp/ocean_trans/av/annual_1yr/',
'/archive/wga/Siena/siena_201303_fix_rds-whit/'+'CM2.5_A_Control-1860_Y03-MBLING-BLING/pp/ocean_trans/ts/annual/20yr/',
'/archive/wga/CM2.6/CM2.6_A_Control-1860_V03/pp/ocean_trans/av/annual_1yr/',
'/archive/rds/Siena/siena_201303_fix_rds/CM_O1p0_C180_A02-topaz-bling-minibling-ctrl-restart_bgc/pp/ocean_trans/av/annual_1yr/']
cms_exp_names = ['CM2-1deg including parameterized \neddy-induced velocities','CM2.5-5day','CM2.6-5day','CM2-1deg excluding parameterized \neddy induced velocities']
plot_legends = ['a)','b)','c)','d)']

first_year = 183
last_year = 192

fig=plt.figure(1,figsize = (14,14))
plt.clf()

for e in range(0,4):
    model = models[e]
    exp = exps[e]
    exp_short_name = exp_short_names[e]
    archdir = archdirs[e]
    cms_exp_name = cms_exp_names[e]
    griddir = '/archive/Carolina.Dufour/Mosaic/'+model+'/'
    figdir = '/work/Henri.Drake/figures/papers/upwelling_timescales_200m/'
    
    subp = fig.add_subplot(2,2,1+e)
    plt.title(plot_legends[e]+' '+cms_exp_name)
    
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
   	if e == 3:
   	    net_trans = ty_trans_rho
        else:
       	    net_trans = ty_trans_rho - ty_trans_rho_gm2
   	trans_30 = net_trans[potrho>1036,np.where(lat_t[:,0]>-30)[0][0]]
    else:
   	trans_30 = ty_trans_rho[potrho>1036,np.where(lat_t[:,0]>-30)[0][0]]
    net_southward = np.sum(trans_30[trans_30<0])
    if model == 'CM2-1deg':
        if e == 3:
   	    overturning = np.flipud(np.cumsum(np.flipud(ty_trans_rho),axis=0))
   	else:
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

    ctf = plt.contourf(lat_av,potrho,-overturning,contours,cmap=my_cmap,norm=norm1)
    plt.contour(lat_av,potrho,-overturning,
   	np.append(np.arange(-50,0,2),np.arange(2,50,2)),colors='black',linewidths=1)
    plt.ylim((1037.5,1034))
    plt.yticks(np.arange(0,4)+1034,('1034','1035','1036','1037'))
    plt.ylabel('Potential density (kgm$^{-3}$)',font)
    plt.xlabel('Latitude',font)
    plt.xlim((-70,-30))
    plt.xticks([-70,-60,-50,-40,-30])
    plt.legend(loc=2,prop=font)
    plt.tick_params(labelsize=tick_font)
    cbar=plt.colorbar(ctf,cmap=my_cmap,norm=norm1,orientation='vertical',ticks = np.arange(-25,30,5))
    cbar.set_label('overturning (Sv)')
    cbar_labels=plt.getp(cbar.ax.axes,'yticklabels')
    plt.setp(cbar_labels,fontsize=tick_font)
    import matplotlib
    matplotlib.rcParams.update({'font.size':15})
    plt.clim(vmin=-40,vmax=40) 

plt.suptitle('Residual overturning streamfunctions for years '+str(first_year)+'-'+str(last_year),fontsize=18)
plt.tight_layout(rect=[0,0,1,0.95],pad=0.4,h_pad = 1, w_pad= 1)
outfile = os.path.join(figdir,'Sfig2_overturning.png')
plt.savefig(outfile,bbox_inches='tight',dpi=200)
#sys.exit()
