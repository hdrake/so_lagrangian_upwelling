#!/usr/bin/env python

## DOESN'T ADJUST TIMES FOR STAGGERED RELEASES??

import numpy as np
import netCDF4 as nc
import sys, os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as col
import matplotlib.cm as cm
import matplotlib.colorbar as cb
from jdcal import jd2gcal
from jdcal import gcal2jd

# Name of CMS experiment:

exp_name = 'CM2.6_1860/'
cms_exp = 'cm26_1860_90year_5day'
#cms_exp = 'expt_cm26_1860_200year_5day'

#exp_name = 'CM2.6_1860_monthly/'
#cms_exp = 'cm26_1860_90year_monthly_rerun'

#exp_name = 'CM2.5_1860/'
#cms_exp = 'cm25_1860_90year_5day'

#exp_name = 'CESM1_2000/'
#cms_exp = 'cesm1_100year_monthly'

#exp_name = 'CM2.6_1990/'
#cms_exp = 'cm26_1990_200year_5day'

path2trajs = '/archive/hfd/CMS/expt_'+cms_exp+'/output/'
path2savedata = '/archive/akm/CMS/expt_'+cms_exp+'/output/'
path2sigma2 = '/archive/hfd/CM2.6/CM2.6_A_Control-1860_V03/5day/sigma2/'

save_directory = '/work/Adele.Morrison/figures/lagrangian_upwelling/'+exp_name
mydata_dir = '/archive/akm/CMS/lagrangian_upwelling/'+exp_name
if exp_name == 'CM2.6_1990/':
	timestep = 5
	numfiles = 4
	npt = 0
	for i in range(1,numfiles+1):
		npt_step = np.shape(np.load(path2trajs+'MLD_upwell_trajs_'+str(i)+'.npy'))[0]
		ntime = np.shape(np.load(path2trajs+'MLD_upwell_trajs_'+str(i)+'.npy'))[1]
		npt += npt_step

	# Quick summary of dataset
	traj_lons = np.zeros((npt,ntime))
	traj_lats = np.zeros((npt,ntime))
	traj_depths = np.zeros((npt,ntime))
	traj_MLD = np.zeros((npt,ntime))
	julian_dates = np.zeros((npt,ntime))
	locations = np.zeros((npt))
	release_dates = np.zeros((npt))
	transports = np.zeros((npt))
	MLD_time = np.zeros((npt))
	stuck_time = np.zeros((npt))
	north_time = np.zeros((npt))

	npt_count = 0
	for i in range(1,numfiles+1):
		print 'reading in file '+str(i)+'/'+str(numfiles)
		trajs = np.load(path2trajs+'MLD_upwell_trajs_'+str(i)+'.npy')
		trajs_meta = np.load(path2trajs+'MLD_upwell_trajs_meta_'+str(i)+'.npy')
		npt_step = np.size(trajs_meta[:,4])

		## FOR DEPTH SURFACE CRITERION
		# Extract variables
		traj_lons[npt_count:npt_count+npt_step,:] = trajs[:,:,0]
		traj_lats[npt_count:npt_count+npt_step,:] = trajs[:,:,1]
		traj_depths[npt_count:npt_count+npt_step,:] = trajs[:,:,2]
		traj_MLD[npt_count:npt_count+npt_step,:] = trajs[:,:,3]
		julian_dates[npt_count:npt_count+npt_step,:] = trajs[:,:,5]
		locations[npt_count:npt_count+npt_step] = trajs_meta[:,0]
		release_dates[npt_count:npt_count+npt_step] = trajs_meta[:,1]
		transports[npt_count:npt_count+npt_step] = trajs_meta[:,2]
		MLD_time[npt_count:npt_count+npt_step] = trajs_meta[:,3]
		stuck_time[npt_count:npt_count+npt_step] = trajs_meta[:,4]
		north_time[npt_count:npt_count+npt_step] = trajs_meta[:,5]

		npt_count += np.size(trajs_meta[:,4])
		
else:
	timestep = 30

	print 'Loading trajectory data.'
	trajs = np.load(path2trajs+'MLD_upwell_trajs.npy')
	trajs_meta = np.load(path2trajs+'MLD_upwell_trajs_meta.npy')

	## FOR DEPTH SURFACE CRITERION
	# Extract variables
	traj_lons = trajs[:,:,0]
	traj_lats = trajs[:,:,1]
	traj_depths = trajs[:,:,2]
	traj_MLD = trajs[:,:,3]
	julian_dates = trajs[:,:,5]

	# Extract variables
	locations = trajs_meta[:,0]
	release_dates = trajs_meta[:,1]
	transports = trajs_meta[:,2]
	MLD_time = trajs_meta[:,3]
	stuck_time = trajs_meta[:,4]
	north_time = trajs_meta[:,5]
	npt = np.size(stuck_time)

mindep = 1000
maxdep = 3500

# Choose depths
dep_range_inds = np.transpose(np.where(np.logical_and(traj_depths[:,0] >= mindep, 
	traj_depths[:,0] < maxdep)))[:,0]
longitude = traj_lons[dep_range_inds,:]
latitude = traj_lats[dep_range_inds,:]
depth = traj_depths[dep_range_inds,:]
locations = locations[dep_range_inds]
cutoff_time_index = MLD_time[dep_range_inds].astype(int)
north_time = north_time[dep_range_inds]
stuck_time = stuck_time[dep_range_inds]
transport = -transports[dep_range_inds]
tot_transport = -transports[dep_range_inds]
julian_dates = julian_dates[dep_range_inds]
release_dates = release_dates[dep_range_inds]
#npt = np.size(dep_range_inds)
#ntime = np.shape(traj_depths)[1]

htfilename = '/archive/Carolina.Dufour/Mosaic/CM2.6/ocean.static.nc'
htfile = nc.Dataset(htfilename)
#yt_ocean = htfile.variables['yt_ocean'][:]
#xt_ocean = htfile.variables['xt_ocean'][:]
ocean_depth = htfile.variables['ht'][:,:]

# Number of particles in trajectories
nparticles = np.shape(longitude)[0]
# Number of time steps in trajectories
ntime = np.shape(longitude)[1]
# Number of trajectory dimensions (usually 3)
ndimensions = 3

# Interpolate sigma2 to each particle's location at each timestep!
# Variable is called "sigma2" and is a function of time, st_ocean (depth), yt_ocean (lat), xt_ocean (lon)
sigma2file_example = path2sigma2+'/nest_1_01820102000000d.nc'
nc_sigma_example = nc.Dataset(sigma2file_example)
st_ocean = nc_sigma_example.variables['st_ocean'][:]
xt_ocean = nc_sigma_example.variables['xt_ocean'][:]
yt_ocean = nc_sigma_example.variables['yt_ocean'][:]
nc_sigma_example.close()

files = os.listdir(path2sigma2)
files = np.sort(files)
# find min and max files used:
julian_date_min = julian_dates.min()
first_year = int((jd2gcal(0,julian_date_min))[0])
first_month = int((jd2gcal(0,julian_date_min))[1])
first_day = int((jd2gcal(0,julian_date_min))[2])
first_sigma2_file_name = 'nest_1_'+str(first_year).zfill(4) + \
	str(first_month).zfill(2) + str(first_day).zfill(2) + '000000d.nc'
first_file_index = np.where(files == first_sigma2_file_name)[0][0]
last_year = int((jd2gcal(0,julian_dates.max()))[0])
last_month = int((jd2gcal(0,julian_dates.max()))[1])
last_day = int((jd2gcal(0,julian_dates.max()))[2])
last_sigma2_file_name = 'nest_1_'+str(last_year).zfill(4) + \
	str(last_month).zfill(2) + str(last_day).zfill(2) + '000000d.nc'
last_file_index = np.where(files == last_sigma2_file_name)[0][0]
files = files[first_file_index:last_file_index+1]

density = np.zeros_like(longitude)
model_lat_index = np.zeros_like(longitude)
model_lon_index = np.zeros_like(longitude)
model_depth_index = np.zeros_like(longitude)
from scipy.interpolate import interp1d
for p in range(nparticles):
	print 'Finding density for particle '+ str(p)+'/'+str(nparticles)
	for t in range(ntime):
		if t > cutoff_time_index[p]:
			break
		# find first year, month
		if t == 0:
			#release_year = int((jd2gcal(0,release_dates[p]))[0])
			#release_month = int((jd2gcal(0,release_dates[p]))[1])
			#release_day = int((jd2gcal(0,release_dates[p]))[2])
			#first_sigma2_file_name = 'nest_1_'+str(release_year).zfill(4) + \
			#	str(release_month).zfill(2) + str(release_day).zfill(2) + '000000d.nc'
			#file_index = np.where(files == first_sigma2_file_name)[0][0]
			file_index = int(release_dates[p]-julian_date_min-0.5)/5
		sigma2_file = path2sigma2+files[file_index]
		ncfile = nc.Dataset(sigma2_file)
			
        # Find closest mixed layer grid point to each particle at this date
		indlon =  np.argmin(np.abs(xt_ocean-longitude[p,t]))
		indlat = np.argmin(np.abs(yt_ocean-latitude[p,t]))
		inddep = np.argmin(np.abs(st_ocean-depth[p,t]))
		# interpolate in depth:
		density_profile = ncfile.variables['sigma2'][0,inddep-1:inddep+2,indlat,indlon]
		depth_interp = np.arange(int(st_ocean[inddep-1]),int(st_ocean[inddep+1]),10)
		interp_function=interp1d(st_ocean[inddep-1:inddep+2],density_profile,
			fill_value=0,kind='linear',bounds_error=False)
		density_interp = interp_function(depth[p,t]).item(0)

		density[p,t] = density_interp
		model_lat_index[p,t] = indlat
		model_lon_index[p,t] = indlon
		model_depth_index[p,t] = inddep
		
		file_index = file_index + timestep/5
		if file_index > len(files):
			file_index = int(julian_dates[p,t+1]-julian_date_min)/5
		
outfile = os.path.join(path2savedata,'sigma2_MLD.npz')
np.savez(outfile,sigma2 = density,
	model_lat_index = model_lat_index,
	model_lon_index = model_lon_index,
	model_depth_index = model_depth_index)
