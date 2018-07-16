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

path2my5day = '/archive/Henri.Drake/'+model+'/'+exp+'5day/for_CMS_SO/'
path2og5day = '/archive/rds/CM2.6/'+exp+'history/'
figdir = '/work/'+user+'/figures/CMS/'

maskpath = path2og5day+'01810101.ocean_minibling_field_salt.nc'
maskfile = nc.Dataset(maskpath)
latsies = maskfile.variables['yt_ocean'][:]
north_ind = np.argmin(np.abs(latsies+29.6331))
latsies = latsies[2:north_ind+1]
landmask = maskfile.variables['salt'][0,:,:,:]
landmask_SO = maskfile.variables['salt'][0,:,2:north_ind+1,:]

path = path2og5day+'01810101.ocean_minibling_field_u.nc'
ncFile = nc.Dataset(path)
data = ncFile.variables['u'][0,:,:,:]
data = np.ma.array(data,mask=landmask.mask)
print '\n5day'
print 'Original u nulls'
print np.sum(np.abs(data) < 1.0e-10)
print np.size(data)
print 'Original u fill'
print np.sum(np.abs(data) > 1.0e10)
print np.size(data)

path = path2my5day+'ocean.01810102.u.nc'
ncFile = nc.Dataset(path)
data = ncFile.variables['u'][0,:,:,:]
data = np.ma.array(data,mask=landmask_SO.mask)
print 'My u nulls'
print np.sum(np.abs(data) < 1.0e-10)
print np.size(data)
print 'My u fill'
print np.sum(np.abs(data) > 1.0e10)
print np.size(data)

path = path2og5day+'01810101.ocean_minibling_field_v.nc'
ncFile = nc.Dataset(path)
data = ncFile.variables['v'][0,:,:,:]
data = np.ma.array(data,mask=landmask.mask)
print 'Original v nulls'
print np.sum(np.abs(data) < 1.0e-10)
print np.size(data)
print 'Original v fill'
print np.sum(np.abs(data) > 1.0e10)
print np.size(data)

path = path2my5day+'ocean.01810102.v.nc'
ncFile = nc.Dataset(path)
data = ncFile.variables['v'][0,:,:,:]
data = np.ma.array(data,mask=landmask_SO.mask)
print 'My v nulls'
print np.sum(np.abs(data) < 1.0e-10)
print np.size(data)
print 'My v fill'
print np.sum(np.abs(data) > 1.0e10)
print np.size(data)

path = path2my5day+'ocean.01810102.w.nc'
ncFile = nc.Dataset(path)
data = ncFile.variables['w'][0,:,:,:]
data = np.ma.array(data,mask=landmask_SO[:,:-1,:].mask)
print 'My w nulls'
print np.sum(np.abs(data) < 1.0e-10)
print np.size(data)
print 'My w fill'
print np.sum(np.abs(data) > 1.0e10)
print np.size(data)

path2myMonthly = '/archive/Henri.Drake/'+model+'/'+exp+'monthly/for_CMS_SO_20yr/'

path = path2myMonthly+'ocean.018201.u.nc'
print '\nMONTHLY'
ncFile = nc.Dataset(path)
datau = ncFile.variables['u'][0,:,:,:]
datau = np.ma.array(datau,mask=landmask_SO.mask)
print 'My u nulls'
print np.sum(np.abs(datau) < 1.0e-20)
print np.size(datau)
print 'My u fill'
print np.sum(np.abs(datau) > 1.0e10)
print np.size(datau)

path = path2myMonthly+'ocean.018201.v.nc'
ncFile = nc.Dataset(path)
datav = ncFile.variables['v'][0,:,:,:]
datav = np.ma.array(datav,mask=landmask_SO.mask)
print 'My v nulls'
print np.sum(np.abs(datav) < 1.0e-20)
print np.size(datav)
print 'My v fill'
print np.sum(np.abs(datav) > 1.0e10)
print np.size(datav)

path = path2myMonthly+'ocean.018201.wt.nc'
ncFile = nc.Dataset(path)
data = ncFile.variables['wt'][0,:,:,:]
lons = ncFile.variables['xt_ocean'][...]
lats = ncFile.variables['yt_ocean'][...]
depths = ncFile.variables['sw_ocean'][...]
print 'My wt nulls'
print np.sum(np.abs(data) < 1.0e-20)
print np.size(data)
print 'My wt fill'
print np.sum(np.abs(data) > 1.0e10)
print np.size(data)

[depthsind,latsind,lonsind] = np.where(np.abs(data) < 1.0e-20)

trimmed_datau = datau[depthsind,latsind,lonsind]
trimmed_datav = datav[depthsind,latsind,lonsind]
#print trimmed_datau
#print trimmed_datav

newdepthsind = np.zeros(len(depthsind))
newdepthsind[np.where(depthsind<49)]=depthsind[np.where(depthsind<49)]+1
newdepthsind[np.where(depthsind==49)]=49
newdepthsind = newdepthsind.astype(int)
newdepths = depths[newdepthsind]

indy = np.where(np.abs(lats[latsind]+30)<0.02)
specdepthsind = depthsind[indy]
speclonsind = lonsind[indy]

fig1=plt.figure(1)
plt.clf()
plt.plot(lons[lonsind],lats[latsind],'r.',markersize=0.05)
plt.ylabel('Latitude')
plt.xlabel('Longitude')
plt.xlim((-280,80))
plt.ylim((-80,-30))
plt.show()
plt.draw()
outfile = os.path.join(figdir,'misc/where_null_vel1.png')
plt.savefig(outfile,bbox_inches='tight',dpi=270)

fig1=plt.figure(1)
plt.clf()
plt.plot(lons[lonsind],-depths[depthsind],'r.',markersize=0.05)
plt.ylabel('Depth')
plt.xlabel('Longitude')
plt.xlim((-280,80))
plt.ylim((-5500,0))
plt.show()
plt.draw()
outfile = os.path.join(figdir,'misc/where_null_vel2.png')
plt.savefig(outfile,bbox_inches='tight',dpi=270)

fig1=plt.figure(1)
plt.clf()
plt.plot(lons[speclonsind],-depths[specdepthsind],'r.',markersize=2)
plt.ylabel('Depth')
plt.xlabel('Longitude')
plt.xlim((-280,80))
plt.ylim((-5500,0))
plt.show()
plt.draw()
outfile = os.path.join(figdir,'misc/where_null_vel3.png')
plt.savefig(outfile,bbox_inches='tight',dpi=270)
