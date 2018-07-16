#!/usr/bin/env python

# import modules
import os,sys,glob
import numpy as np
import netCDF4 as nc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.cm as cm
import matplotlib.colorbar as cb
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator, LinearLocator, NullLocator
np.set_printoptions(threshold=np.nan)

##############

model = 'CM2.6'
exp = 'CM2.6_A_Control-1860_V03/'
exp_short_name = '1860'
user = 'Henri.Drake'

path25day = '/archive/Henri.Drake/'+model+'/'+exp+'5day/for_CMS_SO/'
path2monthly = '/archive/Henri.Drake/'+model+'/'+exp+'monthly/for_CMS_SO/'
figdir = '/work/'+user+'/figures/CMS/'

all_w5000 = np.zeros((73,932,3600))
all_w4000 = np.zeros((73,932,3600))
all_w3000 = np.zeros((73,932,3600))
all_w2000 = np.zeros((73,932,3600))
all_w1500 = np.zeros((73,932,3600))
all_w1000 = np.zeros((73,932,3600))
all_w700 = np.zeros((73,932,3600))
all_w500 = np.zeros((73,932,3600))
all_w300 = np.zeros((73,932,3600))
all_w200 = np.zeros((73,932,3600))
all_w100 = np.zeros((73,932,3600))
all_w5 = np.zeros((73,932,3600))

north = -29.5

####################

griddir = '/archive/Carolina.Dufour/Mosaic/'+model+'/'
grid_file = griddir+model+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
wet = gridFile.variables['wet'][...]
lat_t = gridFile.variables['geolat_t'][...]
long_t = gridFile.variables['geolon_t'][...]

#print lat_t[0:4,0]
boundary_north = -29.5
index_north = np.where(lat_t[:,0]<boundary_north)[0][-1]

maskfile = '/archive/Richard.Slater/CM2.6/CM2.6_A_Control-1860_V03/history/01810101.ocean_minibling_field_salt.nc'
salt_file = nc.Dataset(maskfile)
salt = salt_file.variables['salt'][0,:,:index_north,:]
lonsies = salt_file.variables['yt_ocean'][:index_north]
##############          5-DAY W             ##########

files = os.listdir(path25day)
files = np.sort(files)
count = 0
for file in files:
    if not file.endswith('w.nc'):
        continue
    five_file = nc.Dataset(path25day+file)
    st_ocean = five_file.variables['st_ocean'][:]
    ind_5000 = (np.abs(st_ocean - 5000)).argmin()
    ind_4000 = (np.abs(st_ocean - 4000)).argmin()
    ind_3000 = (np.abs(st_ocean - 3000)).argmin()
    ind_2000 = (np.abs(st_ocean - 2000)).argmin()
    ind_1500 = (np.abs(st_ocean - 1500)).argmin()
    ind_1000 = (np.abs(st_ocean - 1000)).argmin()
    ind_700 = (np.abs(st_ocean - 700)).argmin()
    ind_500 = (np.abs(st_ocean - 500)).argmin()
    ind_300 = (np.abs(st_ocean - 300)).argmin()
    ind_200 = (np.abs(st_ocean - 200)).argmin()
    ind_100 = (np.abs(st_ocean - 100)).argmin()
    ind_5 = (np.abs(st_ocean - 5)).argmin()

    all_w5000[count,:,:] = five_file.variables['w'][0,ind_5000,:,:]
    all_w4000[count,:,:] = five_file.variables['w'][0,ind_4000,:,:]
    all_w3000[count,:,:] = five_file.variables['w'][0,ind_3000,:,:]
    all_w2000[count,:,:] = five_file.variables['w'][0,ind_2000,:,:]
    all_w1500[count,:,:] = five_file.variables['w'][0,ind_1500,:,:]
    all_w1000[count,:,:] = five_file.variables['w'][0,ind_1000,:,:]
    all_w700[count,:,:] = five_file.variables['w'][0,ind_700,:,:]
    all_w500[count,:,:] = five_file.variables['w'][0,ind_500,:,:]
    all_w300[count,:,:] = five_file.variables['w'][0,ind_300,:,:]
    all_w200[count,:,:] = five_file.variables['w'][0,ind_200,:,:]
    all_w100[count,:,:] = five_file.variables['w'][0,ind_100,:,:]
    all_w5[count,:,:] = five_file.variables['w'][0,ind_5,:,:]

    count = count+1
    if count == 73:
        break

w5000 = np.average(all_w5000,axis=0)
w4000 = np.average(all_w4000,axis=0)
w3000 = np.average(all_w3000,axis=0)
w2000 = np.average(all_w2000,axis=0)
w1500 = np.average(all_w1500,axis=0)
w1000 = np.average(all_w1000,axis=0)
w700 = np.average(all_w700,axis=0)
w500 = np.average(all_w500,axis=0)
w300 = np.average(all_w300,axis=0)
w200 = np.average(all_w200,axis=0)
w100 = np.average(all_w100,axis=0)
w5 = np.average(all_w5,axis=0)


all_w5000 = np.zeros((73,932,3600))
all_w4000 = np.zeros((73,932,3600))
all_w3000 = np.zeros((73,932,3600))
all_w2000 = np.zeros((73,932,3600))
all_w1500 = np.zeros((73,932,3600))
all_w1000 = np.zeros((73,932,3600))
all_w700 = np.zeros((73,932,3600))
all_w500 = np.zeros((73,932,3600))
all_w300 = np.zeros((73,932,3600))
all_w200 = np.zeros((73,932,3600))
all_w100 = np.zeros((73,932,3600))
all_w5 = np.zeros((73,932,3600))

#################           MONTHLY W           #################
files = os.listdir(path2monthly)
files = np.sort(files)
count = 0
for file in files:
    if not file.endswith('wt.nc'):
        continue
    month_file = nc.Dataset(path2monthly+file)
    st_ocean = month_file.variables['sw_ocean'][:]
    xt_ocean = month_file.variables['xt_ocean'][:]
    yt_ocean = month_file.variables['yt_ocean'][:]
    #print yt_ocean[-4:-1]
    #print yt_ocean[0:4]
    ind_5000 = (np.abs(st_ocean - 5000)).argmin()
    ind_4000 = (np.abs(st_ocean - 4000)).argmin()
    ind_3000 = (np.abs(st_ocean - 3000)).argmin()
    ind_2000 = (np.abs(st_ocean - 2000)).argmin()
    ind_1500 = (np.abs(st_ocean - 1500)).argmin()
    ind_1000 = (np.abs(st_ocean - 1000)).argmin()
    ind_700 = (np.abs(st_ocean - 700)).argmin()
    ind_500 = (np.abs(st_ocean - 500)).argmin()
    ind_300 = (np.abs(st_ocean - 300)).argmin()
    ind_200 = (np.abs(st_ocean - 200)).argmin()
    ind_100 = (np.abs(st_ocean - 100)).argmin()
    ind_5 = (np.abs(st_ocean - 5)).argmin()

    #print np.shape(month_file.variables['wt'][0,ind_5000,:-1,:])
    #print np.shape(salt[ind_5000,:-1,:])
    all_w5000[count,:,:] = month_file.variables['wt'][0,ind_5000,:-1,:]
    all_w4000[count,:,:] = month_file.variables['wt'][0,ind_4000,:-1,:]
    all_w3000[count,:,:] = month_file.variables['wt'][0,ind_3000,:-1,:]
    all_w2000[count,:,:] = month_file.variables['wt'][0,ind_2000,:-1,:]
    all_w1500[count,:,:] = month_file.variables['wt'][0,ind_1500,:-1,:]
    all_w1000[count,:,:] = month_file.variables['wt'][0,ind_1000,:-1,:]
    all_w700[count,:,:] = month_file.variables['wt'][0,ind_700,:-1,:]
    all_w500[count,:,:] = month_file.variables['wt'][0,ind_500,:-1,:]
    all_w300[count,:,:] = month_file.variables['wt'][0,ind_300,:-1,:]
    all_w200[count,:,:] = month_file.variables['wt'][0,ind_200,:-1,:]
    all_w100[count,:,:] = month_file.variables['wt'][0,ind_100,:-1,:]
    all_w5[count,:,:] = month_file.variables['wt'][0,ind_5,:-1,:]
    count = count+1
    if count == 73:
        break


wt5000 = np.average(all_w5000,axis=0)
wt4000 = np.average(all_w4000,axis=0)
wt3000 = np.average(all_w3000,axis=0)
wt2000 = np.average(all_w2000,axis=0)
wt1500 = np.average(all_w1500,axis=0)
wt1000 = np.average(all_w1000,axis=0)
wt700 = np.average(all_w700,axis=0)
wt500 = np.average(all_w500,axis=0)
wt300 = np.average(all_w300,axis=0)
wt200 = np.average(all_w200,axis=0)
wt100 = np.average(all_w100,axis=0)
wt5 = np.average(all_w5,axis=0)

print lonsies[2:5]
print yt_ocean[0:3]
print lonsies[-3:]
print yt_ocean[-3:]

w5000 = np.ma.array(w5000,mask=salt[ind_5000,2:-1,:].mask)
w4000 = np.ma.array(w4000,mask=salt[ind_4000,2:-1,:].mask)
w3000 = np.ma.array(w3000,mask=salt[ind_3000,2:-1,:].mask)
w2000 = np.ma.array(w2000,mask=salt[ind_2000,2:-1,:].mask)
w1500 = np.ma.array(w1500,mask=salt[ind_1500,2:-1,:].mask)
w1000 = np.ma.array(w1000,mask=salt[ind_1000,2:-1,:].mask)
w700 = np.ma.array(w700,mask=salt[ind_700,2:-1,:].mask)
w500 = np.ma.array(w500,mask=salt[ind_500,2:-1,:].mask)
w300 = np.ma.array(w300,mask=salt[ind_300,2:-1,:].mask)
w200 = np.ma.array(w200,mask=salt[ind_200,2:-1,:].mask)
w100 = np.ma.array(w100,mask=salt[ind_100,2:-1,:].mask)
w5 = np.ma.array(w5,mask=salt[ind_5,2:-1,:].mask)

wt5000 = np.ma.array(wt5000,mask=salt[ind_5000,2:-1,:].mask)
wt4000 = np.ma.array(wt4000,mask=salt[ind_4000,2:-1,:].mask)
wt3000 = np.ma.array(wt3000,mask=salt[ind_3000,2:-1,:].mask)
wt2000 = np.ma.array(wt2000,mask=salt[ind_2000,2:-1,:].mask)
wt1500 = np.ma.array(wt1500,mask=salt[ind_1500,2:-1,:].mask)
wt1000 = np.ma.array(wt1000,mask=salt[ind_1000,2:-1,:].mask)
wt700 = np.ma.array(wt700,mask=salt[ind_700,2:-1,:].mask)
wt500 = np.ma.array(wt500,mask=salt[ind_500,2:-1,:].mask)
wt300 = np.ma.array(wt300,mask=salt[ind_300,2:-1,:].mask)
wt200 = np.ma.array(wt200,mask=salt[ind_200,2:-1,:].mask)
wt100 = np.ma.array(wt100,mask=salt[ind_100,2:-1,:].mask)
wt5 = np.ma.array(wt5,mask=salt[ind_5,2:-1,:].mask)

print 'Saving plots.'
### Plots
north_index = (np.abs(yt_ocean - north)).argmin()
x = xt_ocean
y = yt_ocean[:north_index]
levels=LinearLocator(numticks=201).tick_values(-0.000001,0.000001)
cmap = plt.get_cmap('coolwarm')
norm = BoundaryNorm(levels,ncolors=cmap.N, clip=True)
landlat,landlon = np.where(wet == 0)

###################################################################
#                annual   w from 5-day                            #
###################################################################

fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,w5000,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 5000m from 5-day')
outfile = os.path.join(figdir,'misc/w_ann_5000_5day')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,w4000,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 4000m from 5-day')
outfile = os.path.join(figdir,'misc/w_ann_4000_5day')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,w3000,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 3000m from 5-day')
outfile = os.path.join(figdir,'misc/w_ann_3000_5day')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,w2000,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 2000m from 5-day')
outfile = os.path.join(figdir,'misc/w_ann_2000_5day')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,w1500,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 1500m from 5-day')
outfile = os.path.join(figdir,'misc/w_ann_1500_5day')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,w1000,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 1000m from 5-day')
outfile = os.path.join(figdir,'misc/w_ann_1000_5day')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,w700,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 700m from 5-day')
outfile = os.path.join(figdir,'misc/w_ann_700_5day')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,w500,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 500m from 5-day')
outfile = os.path.join(figdir,'misc/w_ann_500_5day')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,w300,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 300m from 5-day')
outfile = os.path.join(figdir,'misc/w_ann_300_5day')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

##################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,w200,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 200m from 5-day')
outfile = os.path.join(figdir,'misc/w_ann_200_5day')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,w100,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 100m from 5-day')
outfile = os.path.join(figdir,'misc/w_ann_100_5day')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,w5,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 5m from 5-day')
outfile = os.path.join(figdir,'misc/w_ann_5_5day')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
#                annual   w from monthly                          #
###################################################################

fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,wt5000,levels=levels,cmap=cmap)
plt.colorbar()
plt.axis([-280, 80, -82, -29.5])
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.title('Annual-averaged Vertical Velocity at 5000m from monthly')
outfile = os.path.join(figdir,'misc/w_ann_5000_monthly')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,wt4000,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 4000m from monthly')
outfile = os.path.join(figdir,'misc/w_ann_4000_monthly')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,wt3000,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 3000m from monthly')
outfile = os.path.join(figdir,'misc/w_ann_3000_monthly')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,wt2000,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 2000m from monthly')
outfile = os.path.join(figdir,'misc/w_ann_2000_monthly')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,wt1500,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 1500m from monthly')
outfile = os.path.join(figdir,'misc/w_ann_1500_monthly')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,wt1000,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 1000m from monthly')
outfile = os.path.join(figdir,'misc/w_ann_1000_monthly')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,wt700,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 700m from monthly')
outfile = os.path.join(figdir,'misc/w_ann_700_monthly')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,wt500,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 500m from monthly')
outfile = os.path.join(figdir,'misc/w_ann_500_monthly')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,wt300,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 300m from monthly')
outfile = os.path.join(figdir,'misc/w_ann_300_monthly')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,wt200,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 200m from monthly')
outfile = os.path.join(figdir,'misc/w_ann_200_monthly')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,wt100,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 100m from monthly')
outfile = os.path.join(figdir,'misc/w_ann_100_monthly')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

###################################################################
fig1=plt.figure(1)
plt.clf()
plt.contourf(x,y,wt5,levels=levels,cmap=cmap)
plt.colorbar()
plt.contour(long_t,lat_t,wet,colors='k',linewidths=1)
plt.axis([-280, 80, -82, -29.5])
plt.title('Annual-averaged Vertical Velocity at 5m from monthly')
outfile = os.path.join(figdir,'misc/w_ann_5_monthly')
plt.savefig(outfile,bbox_inches = 'tight',dpi = 270)

print 'Done'





