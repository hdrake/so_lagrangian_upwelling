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
from jdcal import gcal2jd
from jdcal import jd2gcal
np.set_printoptions(threshold=np.nan)

##############

model = 'CM2.1'
exp = 'CM2.1_A_Control-1860_V03'
# Ice is saved in 1yr files at daily resolution
path2ice = '/archive/cod/CM2-1deg/CM2-1deg_A_Control-1860_V03/history/'
path2ice199 = '/archive/hfd/CM2-1deg/CM2-1deg_A_Control-1860_V03/5day/Missing_v_files/'
# U,V saved in 1 yr files at 5day resolution
path2files = '/archive/cod/CM2-1deg/CM2-1deg_A_Control-1860_V03/history/'
path2files199 = '/archive/cod/CM2-1deg/CM2-1deg_A_Control-1860_V03/history/tmp/'
path2save = '/archive/hfd/CM2-1deg/CM2-1deg_A_Control-1860_V03/5day/for_CMS_SO/'
path2figures = '/work/hfd/figures/CM2-1deg/misc/'
figdir = path2figures

# Get grid Data
griddir = '/archive/Carolina.Dufour/Mosaic/'+model+'/'
grid_file = griddir+model+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
lat_t = gridFile.variables['geolat_t'][...]
lon_t = gridFile.variables['geolon_t'][...]
wet = gridFile.variables['wet'][...]

# Set parameters
first_year = 181
last_year = 200
boundary_north = -27.6
index_north = np.where(lat_t[:,0]<boundary_north)[0][-1]+2

# Sw_ocean
data = nc.Dataset(path2files+'01810101.ocean_bgc_physics_field_w.nc')
sw_ocean = data.variables['sw_ocean'][:]
xt_ocean = data.variables['xt_ocean'][:]
yt_ocean = data.variables['yt_ocean'][:index_north]
data = nc.Dataset(path2files+'01810101.ocean_bgc_physics_field_u.nc')
st_ocean = data.variables['st_ocean'][:]
xu_ocean = data.variables['xu_ocean'][:]
yu_ocean = data.variables['yu_ocean'][:index_north]
ns = len(sw_ocean)
nx = len(xt_ocean)
ny = len(yt_ocean)

# Make copies of lon thickness and lat thickness of t gridcells for each depth to make calculations easier
dxt = np.zeros((ns,ny,nx))
dyt = np.zeros((ns,ny,nx))
dxu = np.zeros((ns,ny,nx))
dyu = np.zeros((ns,ny,nx))
dxt_uno = gridFile.variables['dxt'][:index_north,:]
dyt_uno = gridFile.variables['dyt'][:index_north,:]
dxu_uno = gridFile.variables['dxu'][:index_north,:]
dyu_uno = gridFile.variables['dyu'][:index_north,:]
ht = gridFile.variables['ht'][:index_north,:]
for i in range(0,ns):
    dxt[i,:,:]=dxt_uno
    dyt[i,:,:]=dyt_uno
    dxu[i,:,:]=dxu_uno
    dyu[i,:,:]=dyu_uno
area_t = dxt*dyt
area_u = dxu*dyu

# The 5-day field files
files = os.listdir(path2files)
files = np.sort(files)

# Start counters
countfiles = 0
icecount = 0

# Define arrays
up_flux = np.zeros((ns,ny,nx))
online_up_flux = np.zeros((ns,ny,nx))
#### COMPUTE uhrho_et and vhrho_nt 5-day averages from u,v,rho_dzt 5-day averages

start_year = 181
countfiles = start_year - 181

# There is one file per year with 73 5-day time steps
for file in files:
    if not file.endswith('v.nc'):
        continue
    year = int(file[0:4])
    if (year < start_year):
        continue
    for countu5day in range(0,73):
        print 'Doing '+str(year)+' and the '+str(countu5day)+'th 5-day period.'
        # Ignore first and last files because we can't get time dependency of dzt or dzu
        if ((year == start_year) and (countu5day == 0)):
            continue
        if ((year == 200) and (countu5day == 72)):
            sys.exit()
        file_name_prefix = str(year).zfill(4)+'0101.ocean_bgc_physics_field_'
        file_name_ice = str(year).zfill(4)+'0101.ice_daily.nc'
        # Get landmask of u and t field grids
        if ((year == start_year) and (countu5day == 1)):
            data = nc.Dataset(path2files+file_name_prefix+'u.nc')
            umasked = data.variables['u'][0,:,0:index_north,:]
            ulon = data.variables['xu_ocean'][:]
            ulat = data.variables['yu_ocean'][0:index_north]
            depth = data.variables['st_ocean'][:]
            first_jdate = data.variables['time'][countu5day]-5
            data = nc.Dataset(path2files+file_name_prefix+'salt.nc')
            tmasked = data.variables['salt'][0,:,0:index_north,:]
            tlon = data.variables['xt_ocean'][:]
            tlat = data.variables['yt_ocean'][0:index_north]
        # rho_dzt is the grid cell thickness (which varies in time) x 1035
        data = nc.Dataset(path2files+file_name_prefix+'rho_dzt.nc')
        if (year != 199):
            SSH_data = nc.Dataset(path2ice+file_name_ice)
            EXT_data = nc.Dataset(path2ice+file_name_ice)
        elif (year == 199):
            SSH_data = nc.Dataset(path2ice199+file_name_ice)
            EXT_data = nc.Dataset(path2ice199+file_name_ice)
        # We take an extra northern latitude because of the algorithm for calculating dzu
        rho_dzt = (data.variables['rho_dzt'][countu5day,:,:index_north+1,:]).filled(fill_value=1.e20)/1035.
        rho_dzt_time = data.variables['time'][countu5day]
        month = int(jd2gcal(1721423,rho_dzt_time)[1])
        day = int(jd2gcal(1721423,rho_dzt_time)[2])
        rho_dzt_t1 = data.variables['average_T1'][countu5day]
        rho_dzt_t2 = data.variables['average_T2'][countu5day]
        ice_time = SSH_data.variables['time'][:]
        # Compute time tendency of dzt because this will artificially change the volume contained in a grid box
        # Special cases where the before timestep is a different year
        if (countu5day == 0):
            # rho_dzt tendency terms from ocean
            file_name_prefix_before = str(year-1).zfill(4)+'0101.ocean_bgc_physics_field_'
            data_before = nc.Dataset(path2files+file_name_prefix_before+'rho_dzt.nc')
            data_after = data
            rho_dzt_before = (data_before.variables['rho_dzt'][72,:,0:index_north,:]).filled(fill_value=1.e20)/1035.
            rho_dzt_after = (data_after.variables['rho_dzt'][1,:,0:index_north,:]).filled(fill_value=1.e20)/1035.

        elif (countu5day == 72):
            # rho_dzt tendency comes from ocean
            file_name_prefix_after = str(year+1).zfill(4)+'0101.ocean_bgc_physics_field_'
            data_before = data
            data_after = nc.Dataset(path2files+file_name_prefix_after+'rho_dzt.nc')
            rho_dzt_before = (data_before.variables['rho_dzt'][71,:,0:index_north,:]).filled(fill_value=1.e20)/1035.
            rho_dzt_after = (data_after.variables['rho_dzt'][0,:,0:index_north,:]).filled(fill_value=1.e20)/1035.
        else:
            # rho_dzt tendency comes from ocean
            data_before = data
            data_after = data
            rho_dzt_before = (data_before.variables['rho_dzt'][countu5day-1,:,0:index_north,:]).filled(fill_value=1.e20)/1035.
            rho_dzt_after = (data_after.variables['rho_dzt'][countu5day+1,:,0:index_north,:]).filled(fill_value=1.e20)/1035.

        if (rho_dzt_t1 < ice_time[0]):
            t2_ind = (np.abs(ice_time - (rho_dzt_t2+0.5))).argmin()
            # SSH tendency terms from ice
            file_name_ice_before = str(year-1).zfill(4)+'0101.ice_daily.nc'
            if (year-1 != 199):
                SSH_data_before = nc.Dataset(path2ice+file_name_ice_before)
                EXT_data_before = nc.Dataset(path2ice+file_name_ice_before)
            elif (year-1 == 199):
                SSH_data_before = nc.Dataset(path2ice199+file_name_ice_before)
                EXT_data_before = nc.Dataset(path2ice199+file_name_ice_before)
            SSH = np.append([SSH_data_before.variables['SSH'][-1,:index_north,:]],SSH_data.variables['SSH'][:t2_ind+1,:index_north,:],axis=0)
            EXT_max = np.max(np.append([EXT_data_before.variables['EXT'][-1,:index_north,:]],EXT_data.variables['EXT'][:t2_ind+1,:index_north,:],axis=0),axis=0)
        elif (rho_dzt_t2 > ice_time[-1]):
            t1_ind = (np.abs(ice_time - (rho_dzt_t1-0.5))).argmin()
            # SSH tendency terms from ice
            file_name_ice_after = str(year+1).zfill(4)+'0101.ice_daily.nc'
            if (year+1 != 199):
                SSH_data_after = nc.Dataset(path2ice+file_name_ice_after)
                EXT_data_after = nc.Dataset(path2ice+file_name_ice_after)
            elif (year+1 == 199):
                SSH_data_after = nc.Dataset(path2ice199+file_name_ice_after)
                EXT_data_after = nc.Dataset(path2ice199+file_name_ice_after)
            SSH = np.append(SSH_data.variables['SSH'][t1_ind:,:index_north,:],[SSH_data_after.variables['SSH'][0,:index_north,:]],axis=0)
            EXT_max = np.max(np.append(EXT_data.variables['EXT'][t1_ind:,:index_north,:],[EXT_data.variables['EXT'][0,:index_north,:]],axis=0),axis=0)
        else:
            t1_ind = (np.abs(ice_time - (rho_dzt_t1-0.5))).argmin()
            t2_ind = (np.abs(ice_time - (rho_dzt_t2+0.5))).argmin()
            SSH = SSH_data.variables['SSH'][t1_ind:t2_ind+1,:index_north,:]
            EXT_max = np.max(EXT_data.variables['EXT'][t1_ind:t2_ind+1,:index_north,:],axis=0)

        # Compute rho_dzu (depth,lat,lon) by taking minimum of four grid cell corners
        rho_dzt_right = np.roll(rho_dzt[:,:-1,:],-1,axis=-1)
        rho_dzt_up = rho_dzt[:,1:,:]
        rho_dzt_up_right = np.roll(rho_dzt[:,1:,:],-1,axis=-1)
        # Define dzu as the minimum of the four corners
        rho_dzu = np.minimum(np.minimum(np.minimum(rho_dzt[:,:-1,:],rho_dzt_right),rho_dzt_up),rho_dzt_up_right)
        # Mask if ground or if 0 thickness
        rho_dzu = np.ma.masked_array(rho_dzu,mask=umasked.mask)
        rho_dzu = rho_dzu.filled(fill_value = 0)
        del rho_dzt_right,rho_dzt_up,rho_dzt_up_right

        print 'Calculating vhrho_nt and uhrho_nt'
        if (year != 200):
            datau = nc.Dataset(path2files+file_name_prefix+'u.nc')
            datav = nc.Dataset(path2files+file_name_prefix+'v.nc')
        elif (year == 200):
            datau = nc.Dataset(path2files199+file_name_prefix+'u.nc')
            datav = nc.Dataset(path2files+file_name_prefix+'v.nc')
        u = (datau.variables['u'][countu5day,:,0:index_north,:]).filled(fill_value=0)
        tim = datav.variables['time'][countu5day]
        average_DT = datav.variables['average_DT'][countu5day]
        v = (datav.variables['v'][countu5day,:,0:index_north,:]).filled(fill_value=0)

        # Compute vhrho_nt and uhrho_ut
        vhrho_net = v*rho_dzu
        vhrho_nwt = np.roll(vhrho_net,1,axis=-1)
        uhrho_net = u*rho_dzu*dyu
        uhrho_set = u[:,:-1,:]*rho_dzu[:,:-1,:]*dyu[:,:-1,:]
        uhrho_set = np.insert(uhrho_set,0,uhrho_set[:,0,:],axis=1)
        del u,v

        # this is still on upper right corner of grid cell at this stage,
        # so interpolate to top center of grid cell:
        vhrho_nt = (0.5*(vhrho_net+vhrho_nwt))
        vhrho_nt = np.ma.array(vhrho_nt,mask=tmasked.mask)
        vhrho_nt = np.ma.filled(vhrho_nt,fill_value=0)

        uhrho_et = (0.5*(uhrho_net+uhrho_set))/dyt
        uhrho_et = np.ma.array(uhrho_et,mask=tmasked.mask)
        uhrho_et = np.ma.filled(uhrho_et,fill_value=0)

        # Try using daily SSH instead (only works where there is no ice)
        # Get from 1 day before to 1 day after this 5 day period
        # Check where this is ice
        EXT_max = np.tile(EXT_max,(len(sw_ocean),1,1))
        rho_dzt_dt = (rho_dzt_after - rho_dzt_before)/(10.*24.*3600.)

        # Make thickness array with partial bottom cells (this will be the same for all timesteps):
        thickness = np.diff(np.append(0,sw_ocean))
        thickness_tile = np.swapaxes(np.tile(thickness,(ht.shape[1],ht.shape[0],1)),0,2)
        depth = np.cumsum(thickness)
        for ii in range(ht.shape[1]):
            for jj in range(ht.shape[0]):
                last_depth_index = np.where(depth > ht[jj,ii])[0]
                if len(last_depth_index) == 0:
                    continue
                else:
                    thickness_tile[last_depth_index[0],jj,ii] += (ht[jj,ii] - \
                                                                      depth[last_depth_index[0]])
        # Mask:
        thickness_tile = np.ma.array(thickness_tile,mask=tmasked.mask)

        # Modify thickness with SSH:
        # Average SSH at start of 5-day period
        SSH1 = np.mean(SSH[:2,...],axis=0)
        SSH5 = np.mean(SSH[-2:,...],axis=0)
        rho_dzt_dt_offline = ((1 + SSH5/ht) - (1 + SSH1/ht))*thickness_tile/(5.*24.*3600.)

        # Under ice use rho_dzt_dt from 5 day averages, not daily averages from ice fields:
        rho_dzt_dt_combined = np.ma.where(EXT_max==0.,rho_dzt_dt_offline,rho_dzt_dt)

        print 'Calculating w'
        # Create volume flux out (positive) and in (negative) to each tracer grid cell in each of 4 horizontal directions
        trans_nt = vhrho_nt * dxu

        trans_st = -(vhrho_nt[:,:-1,:]*dxu[:,:-1,:])
        #trans_st = np.insert(trans_st,0,trans_st[:,0,:],axis=1)
        trans_st = np.insert(trans_st,0,0,axis=1)

        trans_et = uhrho_et * dyt

        trans_wt = -(np.roll(trans_et,1,axis=-1))

        # Calculate change in volume of grid cell due to change in thickness (dzt)
        area_t = dxt*dyt
        delta_thickness_t = rho_dzt_dt_combined
        delta_volume_t = delta_thickness_t*area_t
        trans_from_volume_change = delta_volume_t
        trans_from_volume_change = np.ma.array(trans_from_volume_change,mask=tmasked.mask)
        trans_from_volume_change = np.ma.filled(trans_from_volume_change,fill_value=0)

        # up_flux defined as positive upwards
        out_flux = trans_nt + trans_st + trans_et + trans_wt

        # Define w to be 0 at the bottom cell
        up_flux[49,:,:] = 0
        for i in range(1,50):
            p = 49-i
            up_flux[p,:,:]= - out_flux[p+1,:,:] + up_flux[p+1,:,:] - trans_from_volume_change[p+1,:,:]

        w = up_flux/area_t
        w = np.ma.array(w,mask=tmasked.mask)
        w = np.ma.filled(w,fill_value = -1.e20)

        st_ocean = depth
        xt_ocean = tlon
        yt_ocean = tlat


        ###############################################################
        # save w
        out_filename_w = path2save+'ocean.'+str(year).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+'.wt.nc'
        print 'Saving calculation in: '+out_filename_w
        if os.path.lexists(out_filename_w):
            os.remove(out_filename_w)
        netcdf_file = nc.Dataset(out_filename_w,'w',format='NETCDF4')
        netcdf_file.description = exp
            
        # dimensions
        netcdf_file.createDimension('time', None)
        netcdf_file.createDimension('sw_ocean', len(sw_ocean))
        netcdf_file.createDimension('yt_ocean', len(yt_ocean))
        netcdf_file.createDimension('xt_ocean', len(xt_ocean))
        
        # variables
        ti = netcdf_file.createVariable('time', 'f4', ('time',))
        ti.units = 'days since 0001-01-01 00:00:00'
        ti.calender_type = 'JULIAN'
        ti.calender = 'JULIAN'
        
        sw = netcdf_file.createVariable('sw_ocean', 'f4', ('sw_ocean',))
        sw.units = 'metres'
        sw.long_name = 'tcell zstar depth'
        sw.positive = 'down'
        
        yt = netcdf_file.createVariable('yt_ocean', 'f4', ('yt_ocean',))
        yt.units = 'degrees_N'
        yt.long_name = 'tcell latitude'
        
        xt = netcdf_file.createVariable('xt_ocean', 'f4', ('xt_ocean',))
        xt.units = 'degrees_E'
        xt.long_name = 'tcell longitude'
        
        dt = netcdf_file.createVariable('average_DT', 'f4', ('time',))
        dt.units = 'days'
        dt.long_name = 'length of average period'
        
        w_var = netcdf_file.createVariable('w', 'f4', ('time','sw_ocean','yt_ocean', 'xt_ocean',),fill_value=-1.e20)
        w_var.units = 'm/s'
        w_var.long_name = 'vertical velocity'
        w_var.missing_value =  -1.e+20

        # data
        ti[:] = tim
        sw[:] = sw_ocean
        yt[:] = yt_ocean
        xt[:] = xt_ocean
        #u_var[:,:,:,0] = np.swapaxes(u,2,0)
        w_var[0,:] = w
        dt[:] = average_DT
        
        netcdf_file.close()
        ####################################################################

    countfiles +=1

