#!/usr/efe/evn python
# -*- coding: utf-8 -*-

import netCDF4 as nc
import os,sys
import numpy as np
from jdcal import jd2gcal
from jdcal import gcal2jd
import time as timee

# Default parameter
firstyear = 181

# Command line parameter input
while len(sys.argv) > 1:
    option = (sys.argv[1]);             del sys.argv[1]
    if option == '-s':
        firstyear = int(str(sys.argv[1]));     del sys.argv[1]
    else:
        print sys.argv[0],': invalid option',option
        sys.exit


# Name of cms experiment
cms_exp = 'cm26_100year_monthly'


# Model name from cms_exp
if cms_exp[0:4] == 'cm21':
    model = 'CM2-1deg'
    exp = 'CM2-1deg_A_Control-1860_V03'
elif cms_exp[0:4] == 'cm25':
    model = 'CM2.5'
    exp = 'CM2.5_A_Control-1860_Y03-MBLING-BLING'
elif cms_exp[0:4] == 'cm26':
    model = 'CM2.6'
    exp = 'CM2.6_A_Control-1860_V03'
    first_release_jd = 1787176.0
    startyear = 181
    numprocs = 6
    zfillnum = 1
    total_sim_time = sum(gcal2jd(200,12,31))-sum(gcal2jd(181,1,1))+1
elif cms_exp[0:5] == 'cesm1':
    model = 'CESM1'
    exp = 'CESM1_A_Control-1860_V03'
    first_release_jd = 1745531.5
    startyear = 67
    numprocs = 32
    zfillnum = 2
    total_sim_time = sum(gcal2jd(86,12,31))-sum(gcal2jd(67,1,1))+1

path2trajs = '/archive/hfd/CMS/expt_'+cms_exp+'/output/'
path2release = '/archive/hfd/CMS/input_'+cms_exp+'/'
path2sigma2 = '/archive/hfd/'+model+'/'+exp+'/monthly/sigma2/'
monthly_add = '_monthly'
weight_by_trans = 1

monthslist = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

fast = 0
saving = 1
upwell_flag = 1
save_sigma = 0
numyears = 1
if save_sigma:
    print 'Saving sigma2 for years '+str(firstyear)+' to '+str(firstyear+1)

breakflag = 0
npt = 0
# Get dimensions of trajectory files
# Fast just looks at first file
for i in range(12):
    if breakflag:
        break
    for j in range(1,numprocs+1):
        file = 'traj_file_'+str(j).zfill(zfillnum)+'_'+monthslist[i]+'.nc'
        traj_file = path2trajs+file
        ncFile = nc.Dataset(traj_file)
        time = ncFile.variables['time'][:]
        releasedate = ncFile.variables['releasedate'][...]
        nptstep = len(releasedate)
        ntime = len(time)
        ncFile.close()
        npt += nptstep
        if fast:
            breakflag = 1
            break

# Initialize arrays
traj_lons = np.zeros((npt,ntime))
traj_lats = np.zeros((npt,ntime))
traj_depths = np.zeros((npt,ntime))
exitcodes = np.zeros(npt)
locations = np.zeros(npt)
release_dates = np.zeros(npt)
release_years = np.zeros(npt)
release_months = np.zeros(npt)
release_days = np.zeros(npt)
release_seconds_dec = np.zeros(npt)

breakflag = 0
countnpt = 0
count_location_add = 0
# import trajectory data:
for i in range(12):
    if breakflag:
        break
    count_location_add_step = 0
    for j in range(1,numprocs+1):
        file = 'traj_file_'+str(j).zfill(zfillnum)+'_'+monthslist[i]+'.nc'
        traj_file = path2trajs+file
        ncFile = nc.Dataset(traj_file)
        print 'Now extracting '+file
        traj_lon = ncFile.variables['lon'][:,:]
        traj_lat = ncFile.variables['lat'][:,:]
        traj_depth = ncFile.variables['depth'][:,:]
        time = ncFile.variables['time'][:]
        releasedate = ncFile.variables['releasedate'][...]
        exitcode = ncFile.variables['exitcode'][...]
        location = ncFile.variables['location'][...]
        nptstep = len(releasedate)
        ntime = len(time)
        ncFile.close()
        traj_lon[np.where(traj_lon>80)] = traj_lon[np.where(traj_lon>80)] - 360
        exitcodes[countnpt:countnpt+nptstep] = exitcode
        locations[countnpt:countnpt+nptstep] = location+count_location_add
        release_dates[countnpt:countnpt+nptstep] = releasedate
        traj_lons[countnpt:countnpt+nptstep,:] = traj_lon
        traj_lats[countnpt:countnpt+nptstep,:] = traj_lat
        traj_depths[countnpt:countnpt+nptstep,:] = traj_depth
        countnpt += nptstep
        count_location_add_step += nptstep
        if fast:
            breakflag = 1
            break
    count_location_add += count_location_add_step



# Convert negative times to positive times (2**32 is highest number allowed in CMS time)
maxt = 2**31
times = np.zeros(ntime)
switch = 0
multiplier = 0
for t in range(ntime):
    if time[t]<0:
        times[t] = time[t]+2*maxt
    else:
        times[t] = time[t]


# Read in the transport for each particle
if weight_by_trans:
    transports = np.zeros(npt)
    ifile = open(path2release+'releaseFile_master_'+model+monthly_add,'r')
    countlines = 0
    counti = 0
    for line in ifile:
        linesplit = line.split('  ')
        loca = linesplit[0]
        locations[counti] = loca
        if (linesplit[-1] == '') or (linesplit[-1] == '\n'):
            transports[counti] = -np.abs(float(linesplit[-2]))
            counti += 1
        elif (linesplit[-1][-1:] == '\n'):
            transports[counti] = -np.abs(float(linesplit[-1][:-1]))
            counti += 1
        else:
            transports[counti] = -np.abs(float(linesplit[-1]))
            counti += 1
        if (counti == npt):
            break
        countlines += 1

# Sort by location number!
locations_sort = locations.argsort()
locations = locations[locations_sort]
release_dates = release_dates[locations_sort]
traj_lons = traj_lons[locations_sort,:]
traj_lats = traj_lats[locations_sort,:]
traj_depths = traj_depths[locations_sort,:]
exitcodes = exitcodes[locations_sort]
transports = transports[locations_sort]

# Copy the trajectory files to save an upwelling subset of trajectories
traj_depths_OG = np.copy(traj_depths)
traj_lons_OG = np.copy(traj_lons)
traj_lats_OG = np.copy(traj_lats)

if upwell_flag:
    # Masking particles after they upwell above 1000m!
    print 'Masking upwelled particles'
    upwell_time = np.zeros(npt)
    for pt in range(0,npt):
        if np.size(np.where(np.logical_and((traj_depths_OG[pt,:] < 1000),(traj_depths_OG[pt,:] > 0)))):
            upwell_time[pt] = np.min(np.where(np.logical_and((traj_depths_OG[pt,:] < 1000),(traj_depths_OG[pt,:] > 0))))
            #traj_lats[pt,upwell_time[pt]:] = np.nan
            #traj_lons[pt,upwell_time[pt]:] = np.nan
            #traj_depths[pt,upwell_time[pt]:] = np.nan
    
    # Only consider particles that upwell
    # Meta data
    upwell1000_inds = np.transpose(np.where(upwell_time))[:,0]
    upwell_time = upwell_time[upwell1000_inds]
    exitcodes = exitcodes[upwell1000_inds]
    locations = locations[upwell1000_inds]
    release_dates = release_dates[upwell1000_inds]
    transports = transports[upwell1000_inds]
    npt = np.size(upwell1000_inds)
    
    # Trajectory data
    traj_lats = traj_lats[upwell1000_inds,:]
    traj_lats_OG = traj_lats_OG[upwell1000_inds,:]
    
    traj_lons = traj_lons[upwell1000_inds,:]
    traj_lons_OG = traj_lons_OG[upwell1000_inds,:]
    
    traj_depths = traj_depths[upwell1000_inds,:]
    traj_depths_OG = traj_depths_OG[upwell1000_inds,:]

# Calculate final julian date for each particle releasetime.
print 'Calculating time indices of final julian date for each particle release time.'
days_between_particle_release_and_start_date = np.zeros(npt)
for p in range(npt):
    release_years[p] = int((jd2gcal(0,release_dates[p]))[0])
    release_months[p] = int((jd2gcal(0,release_dates[p]))[1])
    release_days[p] = int((jd2gcal(0,release_dates[p]))[2])
    release_seconds_dec[p] = int((jd2gcal(0,release_dates[p]))[3])
    # 181,1,07 is start date for CM2.6
    # 67, 1, 1 is start date for CESM1
    days_between_particle_release_and_start_date[p] = sum(gcal2jd(release_years[p],release_months[p],release_days[p]))+release_seconds_dec[p]-first_release_jd

days_after_release = np.zeros(ntime)
for t in range(0,ntime):
    days_after_release[t] = float(times[t])/86400

# Calculate julian date of each particle at each time step,
# where gcal2jd(181,1,7) = 1787175.5 and gcal2jd(67,1,1) = 1745531.5
julian_dates = np.zeros((npt,ntime))
for p in range(npt):
    for t in range(ntime):
        julian_dates[p,t] = np.mod(days_after_release[t]+days_between_particle_release_and_start_date[p],total_sim_time) + first_release_jd

# Interpolate sigma2 to each particle's location at each timestep!
# Variable is called "sigma2" and is a function of time, st_ocean (depth), yt_ocean (lat), xt_ocean (lon)
sigma2file_example = path2sigma2+'/nest_1_'+str(startyear).zfill(4)+'0115000000d.nc'
nc_sigma_example = nc.Dataset(sigma2file_example)
st_ocean = nc_sigma_example.variables['st_ocean']
xt_ocean = nc_sigma_example.variables['xt_ocean']
yt_ocean = nc_sigma_example.variables['yt_ocean']

# Interpolate sigma2 to particle locations for every particle and at every time.
# Algorithm speed is optimized as much as I knew how.

if save_sigma:
    for yeari in range(startyear,startyear+20):
        sigma2 = np.nan*np.ones((npt,ntime))
        for monthi in range(1,13):
            if ((yeari < firstyear) or (yeari > firstyear+numyears-1)):
                continue
            # Find which two files to read in!
            sigma2file = path2sigma2+'/nest_1_'+str(yeari).zfill(4)+str(monthi).zfill(2)+'15000000d.nc'
            nc_sigma = nc.Dataset(sigma2file)
            if monthi == 12:
                if yeari == startyear+19:
                    nextyeari = startyear
                else:
                    nextyeari = yeari+1
                nextmonthi = 1
            else:
                nextmonthi = monthi+1
                nextyeari = yeari
            sigma2file_next = path2sigma2+'/nest_1_'+str(nextyeari).zfill(4)+str(nextmonthi).zfill(2)+'15000000d.nc'
            nc_sigma_next = nc.Dataset(sigma2file_next)
            sigma2_tmp = nc_sigma.variables['sigma2'][0,:,:,:].data
            sigma2_nxt = nc_sigma_next.variables['sigma2'][0,:,:,:].data
            # Calculate julian date of the sigma2files!
            juldate_tmp = sum(gcal2jd(yeari,monthi,15))
            juldate_nxt = sum(gcal2jd(nextyeari,nextmonthi,15))
            # Find particle and timesteps that fit between these two files.
            if (nextmonthi == 1) and (nextyeari == startyear):
                [indp,indt] = np.where(np.logical_or((julian_dates > juldate_tmp),(julian_dates < juldate_nxt)))
            else:
                [indp,indt] = np.where(np.logical_and((julian_dates > juldate_tmp),(julian_dates < juldate_nxt)))
            # Next, find depths / lat / lon to interpolate to
            print 'For year '+str(yeari)+' and month '+str(monthi)+', there are '+str(np.size(indp))+' particles.'
            for i in range(np.size(indp)):
                if not np.mod(i,5000):
                    print 'Now at '+str(int(i/(float(np.size(indp)))*100))+'%'
                dep = traj_depths[indp[i],indt[i]]
                lon = traj_lons[indp[i],indt[i]]
                lat = traj_lats[indp[i],indt[i]]
                if (dep>10000):
                    continue
                indlon = [int((lon-xt_ocean[0])*(1/(xt_ocean[1]-xt_ocean[0]))),
                          int((lon-xt_ocean[0])*(1/(xt_ocean[1]-xt_ocean[0])))+1]
                indlat = [int((lat-yt_ocean[0])*(1/(yt_ocean[1]-yt_ocean[0]))),
                          int((lat-yt_ocean[0])*(1/(yt_ocean[1]-yt_ocean[0])))+1]
                inddep = np.sort(np.argsort(np.abs(dep - st_ocean))[0:2])
                if indlon[1] == 3600:
                    indlon[1] = 0
                sigma2[indp[i],indt[i]] = np.interp(lat,yt_ocean[inddep],
                                                    [np.interp(lon,xt_ocean[inddep],
                                                               [np.interp(dep,st_ocean[inddep],sigma2_tmp[inddep,indlat[0],indlon[0]]),
                                                                np.interp(dep,st_ocean[inddep],sigma2_tmp[inddep,indlat[1],indlon[0]])]),
                                                     np.interp(lon,xt_ocean[inddep],
                                                               [np.interp(dep,st_ocean[inddep],sigma2_tmp[inddep,indlat[0],indlon[1]]),
                                                                np.interp(dep,st_ocean[inddep],sigma2_tmp[inddep,indlat[1],indlon[1]])])])
            # Interpolate sigma2 in time, depth, lat and lon
            nc_sigma_next.close()
            nc_sigma.close()
        if save_sigma and (yeari == firstyear):
            outfile = path2trajs+'sigma_upwell1000_'+str(yeari).zfill(4)+'.npy'
            print 'Saving '+'sigma_upwell1000_'+str(yeari).zfill(4)+'.npy'
            np.save(outfile,sigma2)
    if save_sigma:    
        sys.exit()
else:
    sigma2 = np.nan*np.ones((npt,ntime))
    for i in range(20):
        print 'Loading '+'sigma_upwell1000_'+str(startyear+i).zfill(4)+'.npy'
        outfile_tmp = path2trajs+'sigma_upwell1000_'+str(startyear+i).zfill(4)+'.npy'
        sigma2_tmp = np.load(outfile_tmp)
        sigma2[sigma2_tmp > 0] = sigma2_tmp[sigma2_tmp > 0]

# Masking stuck particles
print 'Masking stuck particles'
traj_lats_after = np.copy(traj_lats_OG[:,:-2])
traj_lats_before = np.copy(traj_lats_OG[:,1:-1])
traj_lons_after = np.copy(traj_lons_OG[:,:-2])
traj_lons_before = np.copy(traj_lons_OG[:,1:-1])
inds = np.logical_and(np.logical_and(np.abs(traj_lats_after - traj_lats_before) < 0.00001,np.abs(traj_lons_after - traj_lons_before) < 0.00001),np.abs(traj_lats_before) < 1000)
traj_lats[inds] = np.nan
traj_lons[inds] = np.nan
traj_depths[inds] = np.nan
stuck_time = np.zeros(npt)
for pt in range(0,npt):
    tmp = np.where(np.logical_and(np.logical_and(np.abs(traj_lats_after[pt,:] - traj_lats_before[pt,:]) < 0.00001,np.abs(traj_lons_after[pt,:] - traj_lons_before[pt,:]) < 0.00001),np.abs(traj_lats_before[pt,:]) < 1000))
    if np.size(tmp):
        stuck_time[pt] = tmp[0][0]+1

# Fast mask particles that go north
north_time = np.zeros(npt)
for pt in range(0,npt):
    tmp = np.where(traj_lats_OG[pt,:-1]>-29.9)
    if np.size(tmp):
        north_time[pt] = tmp[0][0]
traj_lats[traj_lats_OG>-29.9] = np.nan
traj_depths[traj_lats_OG>-29.9] = np.nan
traj_lons[traj_lats_OG>-29.9] = np.nan

print 'Number of particles released: '+str(npt)
print 'Number of particles that upwell: '+str(np.size(np.where(upwell_time)))
print 'Number of particles that get stuck: '+str(np.size(np.where(stuck_time)))
print 'Number of particles that go north of 30 South: '+str(np.size(np.where(north_time)))

if not upwell_flag:
    # Save array of trajectories
    trajs = np.zeros((npt,ntime,8))
    trajs[:,:,0] = traj_lons
    trajs[:,:,1] = traj_lats
    trajs[:,:,2] = traj_depths
    trajs[:,:,3] = traj_lons_OG
    trajs[:,:,4] = traj_lats_OG
    trajs[:,:,5] = traj_depths_OG
    trajs[:,:,6] = sigma2
    trajs[:,:,7] = julian_dates
    outfile = path2trajs+'trajs.npy'
    if saving:
        np.save(outfile,trajs)

    # Save array of traj meta data
    trajs_meta = np.zeros((npt,6))
    trajs_meta[:,0] = locations
    trajs_meta[:,1] = release_dates
    trajs_meta[:,2] = transports
    trajs_meta[:,3] = upwell_time
    trajs_meta[:,4] = stuck_time
    trajs_meta[:,5] = north_time
    outfile = path2trajs+'trajs_meta.npy'
    if saving:
        np.save(outfile,trajs_meta)
else:
    # Save array of trajectories upwelling
    trajs = np.zeros((npt,ntime,8))
    trajs[:,:,0] = traj_lons
    trajs[:,:,1] = traj_lats
    trajs[:,:,2] = traj_depths
    trajs[:,:,3] = traj_lons_OG
    trajs[:,:,4] = traj_lats_OG
    trajs[:,:,5] = traj_depths_OG
    trajs[:,:,6] = sigma2
    trajs[:,:,7] = julian_dates
    outfile = path2trajs+'trajs_upwell1000.npy'
    if saving:
        np.save(outfile,trajs)

    # Save array of traj meta data upwelling
    trajs_meta = np.zeros((npt,6))
    trajs_meta[:,0] = locations
    trajs_meta[:,1] = release_dates
    trajs_meta[:,2] = transports
    trajs_meta[:,3] = upwell_time
    trajs_meta[:,4] = stuck_time
    trajs_meta[:,5] = north_time
    outfile = path2trajs+'trajs_upwell1000_meta.npy'
    if saving:
        np.save(outfile,trajs_meta)
