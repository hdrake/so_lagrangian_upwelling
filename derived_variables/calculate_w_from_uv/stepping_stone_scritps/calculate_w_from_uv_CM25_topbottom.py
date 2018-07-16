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
dxt = np.zeros((50,363,1440))
dyt = np.zeros((50,363,1440))
dxu = np.zeros((50,363,1440))
dyu = np.zeros((50,363,1440))
dxt_uno = gridFile.variables['dxt'][:index_north,:]
dyt_uno = gridFile.variables['dyt'][:index_north,:]
dxu_uno = gridFile.variables['dxu'][:index_north,:]
dyu_uno = gridFile.variables['dyu'][:index_north,:]
for i in range(0,50):
    dxt[i,:,:]=dxt_uno
    dyt[i,:,:]=dyt_uno
    dxu[i,:,:]=dxu_uno
    dyu[i,:,:]=dyu_uno

# The 5-day field files
files = os.listdir(path2files)
files = np.sort(files)

# Start counters
countu5day = 0
countfiles = 0

# Define arrays
up_flux = np.zeros((50,363,1440))
online_up_flux = np.zeros((50,363,1440))
#### COMPUTE uhrho_et and vhrho_nt 5-day averages from u,v,rho_dzt 5-day averages

# There is one file per year with 73 5-day time steps
for file in files:
    if not file.endswith('u.nc'):
        continue
    for countu5day in range(0,73):
        year = int(file[0:4])
        print 'Doing '+str(year)+' and the '+str(countu5day)+'th 5-day period.'
        # Ignore first and last files because we can't get time dependency of dzt or dzu
        if ((year == 190) and (countu5day == 0)):
            countu5day=1
            continue
        if ((year == 200) and (countu5day == 72)):
            sys.exit()
        file_name_prefix = str(year).zfill(4)+'0101.ocean_bgc_physics_field_'
        # Get landmask of u and t field grids
        if ((year == 190) and (countu5day == 1)):
            data = nc.Dataset(path2files+file_name_prefix+'u.nc')
            umasked = data.variables['u'][0,:,0:index_north,:]
            ulon = data.variables['xu_ocean'][:]
            ulat = data.variables['yu_ocean'][:]
            depth = data.variables['st_ocean'][:]
            first_jdate = data.variables['time'][countu5day]-5
            data = nc.Dataset(path2files+file_name_prefix+'salt.nc')
            tmasked = data.variables['salt'][0,:,0:index_north,:]
            tlon = data.variables['xt_ocean'][:]
            tlat = data.variables['yt_ocean'][:]
        # rho_dzt is the grid cell thickness (which varies in time) x 1035
        data = nc.Dataset(path2files+file_name_prefix+'rho_dzt.nc')
        # We take an extra northern latitude because of the algorithm for calculating dzu
        rho_dzt = (data.variables['rho_dzt'][countu5day,:,:index_north+1,:]).filled(fill_value=1.e20)/1035.
        # Compute time tendency of dzt because this will artificially change the volume contained in a grid box
        # Special cases where the before timestep is a different year
        if (countu5day == 0):
            file_name_prefix_before = str(year-1).zfill(4)+'0101.ocean_bgc_physics_field_'
            print file_name_prefix_before
            data_before = nc.Dataset(path2files+file_name_prefix_before+'rho_dzt.nc')
            data_after = data
            rho_dzt_before = (data_before.variables['rho_dzt'][72,:,0:index_north,:]).filled(fill_value=1.e20)/1035.
            rho_dzt_after = (data_after.variables['rho_dzt'][1,:,0:index_north,:]).filled(fill_value=1.e20)/1035.
        elif (countu5day == 72):
            file_name_prefix_after = str(year+1).zfill(4)+'0101.ocean_bgc_physics_field_'
            data_before = data
            data_after = nc.Dataset(path2files+file_name_prefix_after+'rho_dzt.nc')
            rho_dzt_before = (data_before.variables['rho_dzt'][71,:,0:index_north,:]).filled(fill_value=1.e20)/1035.
            rho_dzt_after = (data_after.variables['rho_dzt'][0,:,0:index_north,:]).filled(fill_value=1.e20)/1035.
        else:
            data_before = data
            data_after = data
            rho_dzt_before = (data_before.variables['rho_dzt'][countu5day-1,:,0:index_north,:]).filled(fill_value=1.e20)/1035.
            rho_dzt_after = (data_after.variables['rho_dzt'][countu5day+1,:,0:index_north,:]).filled(fill_value=1.e20)/1035.

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
        data = nc.Dataset(path2files+file_name_prefix+'u.nc')
        u = (data.variables['u'][countu5day,:,0:index_north,:]).filled(fill_value=0)
        data = nc.Dataset(path2files+file_name_prefix+'v.nc')
        tim = data.variables['time'][countu5day]
        average_DT = data.variables['average_DT'][countu5day]
        v = (data.variables['v'][countu5day,:,0:index_north,:]).filled(fill_value=0)

        # Compute vhrho_nt and uhrho_ut
        vhrho_net = v*rho_dzu
        vhrho_nwt = np.roll(vhrho_net,1,axis=-1)
        uhrho_net = u*rho_dzu*dyu
        uhrho_set = u[:,:-1,:]*rho_dzu[:,:-1,:]*dyu[:,:-1,:]
        uhrho_set = np.insert(uhrho_set,0,0,axis=1)
        del u,v

        # this is still on upper right corner of grid cell at this stage,
        # so interpolate to top center of grid cell:
        vhrho_nt = (0.5*(vhrho_net+vhrho_nwt))
        vhrho_nt = np.ma.array(vhrho_nt,mask=tmasked.mask)
        vhrho_nt = np.ma.filled(vhrho_nt,fill_value=0)

        uhrho_et = (0.5*(uhrho_net+uhrho_set))/dyt
        uhrho_et = np.ma.array(uhrho_et,mask=tmasked.mask)
        uhrho_et = np.ma.filled(uhrho_et,fill_value=0)

        ## Online vhrho_nt and uhrho_et
        budget_file_name = path2files+str(year).zfill(4)+'0101.ocean_budgets_5d.nc'
        data = nc.Dataset(budget_file_name)
        online_vhrho_nt = data.variables['vhrho_nt'][countu5day,:,0:index_north,:]
        online_vhrho_nt = np.ma.array(online_vhrho_nt,mask=tmasked.mask)
        online_vhrho_nt = np.ma.filled(online_vhrho_nt,fill_value=0)/1035.
        online_uhrho_et = data.variables['uhrho_et'][countu5day,:,0:index_north,:]
        online_uhrho_et = np.ma.array(online_uhrho_et,mask=tmasked.mask)
        online_uhrho_et = np.ma.filled(online_uhrho_et,fill_value=0)/1035.

        # Plots
        colors = ['DarkBlue','Blue','Cyan','Lightblue','White','Yellow','Orange','Maroon','Red']
        levels = [-40,-20,-10,-5,-1,1,5,10,20,40]
        diff_levels = [-1,-0.1,-0.01,-0.001,-0.0001,0.0001,0.001,0.01,0.1,1]
        x = tlon
        y = tlat[:index_north]
        
        if ((year == 190) and (countu5day == 0)):
            for kloop in range(0,50):
                # vhrho_nt
                fig=plt.figure(1)
                plt.clf()
                plt.contourf(x,y,vhrho_nt[kloop,:,:],levels=levels,colors=colors,extend='both')
                plt.colorbar()
                plt.contour(lon_t,lat_t,wet,colors='k',linewidths=2)
                plt.axis([-280,80,-82,-29.5])
                plt.title('Snapshot of vhrho_nt')
                outfile = figdir+'vhrho_nt'+str(kloop)+'.png'
                plt.savefig(outfile,bbox_inches = 'tight',dpi=270)
                
                # uhrho_et
                fig=plt.figure(1)
                plt.clf()
                plt.contourf(x,y,uhrho_et[kloop,:,:],levels=levels,colors=colors,extend='both')
                plt.colorbar()
                plt.contour(lon_t,lat_t,wet,colors='k',linewidths=2)
                plt.axis([-280,80,-82,-29.5])
                plt.title('Snapshot of uhrho_et')
                outfile = figdir+'uhrho_et'+str(kloop)+'.png'
                plt.savefig(outfile,bbox_inches = 'tight',dpi=270)

                # Online vhrho_nt
                fig=plt.figure(1)
                plt.clf()
                plt.contourf(x,y,online_vhrho_nt[kloop,:,:],levels=levels,colors=colors,extend='both')
                plt.colorbar()
                plt.contour(lon_t,lat_t,wet,colors='k',linewidths=2)
                plt.axis([-280,80,-82,-29.5])
                plt.title('Snapshot of online_vhrho_nt')
                outfile = figdir+'online_vhrho_nt'+str(kloop)+'.png'
                plt.savefig(outfile,bbox_inches = 'tight',dpi=270)
                
                # Online uhrho_et
                fig=plt.figure(1)
                plt.clf()
                plt.contourf(x,y,online_uhrho_et[kloop,:,:],levels=levels,colors=colors,extend='both')
                plt.colorbar()
                plt.contour(lon_t,lat_t,wet,colors='k',linewidths=2)
                plt.axis([-280,80,-82,-29.5])
                plt.title('Snapshot of online_uhrho_et')
                outfile = figdir+'online_uhrho_et'+str(kloop)+'.png'
                plt.savefig(outfile,bbox_inches = 'tight',dpi=270)

                # vhrho_nt Difference
                fig=plt.figure(1)
                plt.clf()
                plt.contourf(x,y,vhrho_nt[kloop,:,:]-online_vhrho_nt[kloop,:,:],levels=diff_levels,colors=colors,extend='both')
                plt.colorbar()
                plt.contour(lon_t,lat_t,wet,colors='k',linewidths=2)
                plt.axis([-280,80,-82,-29.5])
                plt.title('Snapshot of difference between vhrho_nt and online_vhrho_nt')
                outfile = figdir+'diff_vhrho_nt'+str(kloop)+'.png'
                plt.savefig(outfile,bbox_inches = 'tight',dpi=270)
                
                # uhrho_et Difference
                fig=plt.figure(1)
                plt.clf()
                plt.contourf(x,y,uhrho_et[kloop,:,:]-online_uhrho_et[kloop,:,:],levels=diff_levels,colors=colors,extend='both')
                plt.colorbar()
                plt.contour(lon_t,lat_t,wet,colors='k',linewidths=2)
                plt.axis([-280,80,-82,-29.5])
                plt.title('Snapshot of difference between uhrho_et and online_uhrho_et')
                outfile = figdir+'diff_uhrho_et'+str(kloop)+'.png'
                plt.savefig(outfile,bbox_inches = 'tight',dpi=270)

        #sys.exit()

        print 'Calculating w'
        # Create volume flux out (positive) and in (negative) to each tracer grid cell in each of 4 horizontal directions
        trans_nt = vhrho_nt * dxu

        trans_st = -(vhrho_nt[:,:-1,:]*dxu[:,:-1,:])
        trans_st = np.insert(trans_st,0,0,axis=1)

        trans_et = uhrho_et * dyt

        trans_wt = -(np.roll(trans_et,1,axis=-1))

        # Calculate change in volume of grid cell due to change in thickness (dzt)
        area_t = dxt*dyt
        delta_thickness_t = rho_dzt_after - rho_dzt_before
        delta_volume_t = delta_thickness_t*area_t
        trans_from_volume_change = delta_volume_t/(3600.0*24.0*10.0)
        trans_from_volume_change = np.ma.array(trans_from_volume_change,mask=tmasked.mask)
        trans_from_volume_change = np.ma.filled(trans_from_volume_change,fill_value=0)

        ## Online w
        w_file_name = path2files+str(year).zfill(4)+'0101.ocean_bgc_physics_field_w.nc'
        data = nc.Dataset(w_file_name)
        online_w = data.variables['w'][countu5day,:,:index_north,:]
        online_w = np.ma.filled(online_w,fill_value = -1.e20)

        # up_flux defined as positive upwards
        out_flux = trans_nt + trans_st + trans_et + trans_wt
        up_flux[0,:,:] = online_w[0,:,:] * area_t[0,:,:]

        for i in range(1,50):
            up_flux[i,:,:] = -(- out_flux[i,:,:] - up_flux[i-1,:,:] - trans_from_volume_change[i,:,:])

        w = up_flux/area_t
        w = np.ma.array(w,mask=tmasked.mask)
        w = np.ma.filled(w,fill_value = -1.e20)


        print 'Calculating w from online fluxes'
        # Create volume flux out (positive) and in (negative) to each tracer grid cell in each of 4 horizontal directions
        online_trans_nt = online_vhrho_nt * dxu

        online_trans_st = -(online_vhrho_nt[:,:-1,:]*dxu[:,:-1,:])
        #online_trans_st = np.insert(online_trans_st,0,online_trans_st[:,0,:],axis=1)
        online_trans_st = np.insert(online_trans_st,0,0,axis=1)

        online_trans_et = online_uhrho_et * dyt

        online_trans_wt = -(np.roll(online_trans_et,1,axis=-1))

        # up_flux defined as positive upwards
        online_out_flux = online_trans_nt + online_trans_st + online_trans_et + online_trans_wt

        # Define w to be 0 at the bottom cell
        online_up_flux[49,:,:] = 0
        for i in range(1,50):
            p = 49-i
            online_up_flux[p,:,:]= - online_out_flux[p+1,:,:] + online_up_flux[p+1,:,:] - trans_from_volume_change[p+1,:,:]

        my_online_w = online_up_flux/area_t
        my_online_w = np.ma.array(w,mask=tmasked.mask)
        my_online_w = np.ma.filled(w,fill_value = -1.e20)

        # Plots
        colors = ['Blue','Cyan','Lightblue','White','Yellow','Orange','Maroon']
        levels = [-5.e-4,-5.e-5,-5.e-6,-5.e-7,5.e-7,5.e-6,5.e-5,5.e-4]
        levels_per = [0,0.01,0.1,0.5,1]
        diff_levels = [0,5e-9,5e-8,5e-7,5e-6,5e-5]
        diff_colors = ['White','Yellow','Orange','Brown','Maroon']
        colors_per = ['White','Yellow','Orange','Maroon']
        x = tlon
        y = tlat[:index_north]

        if ((year == 190) and (countu5day == 1)):
            for d in range(0,50):
                # w
                fig=plt.figure(1)
                plt.clf()
                CS=plt.contourf(x,y,w[d,:,:],levels=levels,colors=colors,extend='both')
                CS.cmap.set_under('DarkBlue')
                CS.cmap.set_over('Red')
                plt.colorbar()
                plt.contour(lon_t,lat_t,wet,colors='k',linewidths=2)
                plt.axis([-280,80,-82,-29.5])
                plt.title('Snapshot of w at level '+str(d))
                outfile = figdir+'w_topbottom'+str(d)+'.png'
                plt.savefig(outfile,dpi=270)
                
                # online_w
                fig=plt.figure(1)
                plt.clf()
                CS=plt.contourf(x,y,online_w[d,:,:],levels=levels,colors=colors,extend='both')
                CS.cmap.set_under('DarkBlue')
                CS.cmap.set_over('Red')
                plt.colorbar()
                plt.contour(lon_t,lat_t,wet,colors='k',linewidths=2)
                plt.axis([-280,80,-82,-29.5])
                plt.title('Snapshot of online_w at level '+str(d))
                outfile = figdir+'online_w'+str(d)+'.png'
                plt.savefig(outfile,dpi=270)
                
                # w / online_w Difference
                fig=plt.figure(1)
                plt.clf()
                diffy = np.abs(w[d,:,:]-online_w[d,:,:])
                CS=plt.contourf(x,y,diffy,levels=diff_levels,colors=diff_colors,extend='max')
                CS.cmap.set_over('Red')
                plt.colorbar()
                plt.contour(lon_t,lat_t,wet,colors='k',linewidths=2)
                plt.axis([-280,80,-82,-29.5])
                plt.title('Snapshot of difference between w and online_w at level '+str(d))
                outfile = figdir+'diff_w_topbottom'+str(d)+'.png'
                plt.savefig(outfile,dpi=270)


                # w / online_w Per-Difference
                fig=plt.figure(1)
                plt.clf()
                diffy = w[d,:,:]-online_w[d,:,:]
                divizor = online_w[d,:,:]
                divizor[divizor == 0] = 1e20
                diffy_per = np.abs(diffy/divizor)
                CS=plt.contourf(x,y,diffy_per,levels=levels_per,colors=colors_per,extend='max')
                CS.cmap.set_over('Red')
                plt.colorbar()
                plt.contour(lon_t,lat_t,wet,colors='k',linewidths=2)
                plt.axis([-280,80,-82,-29.5])
                plt.title('Snapshot of difference between w and online_w at level '+str(d))
                outfile = figdir+'per_diff_w_topbottom'+str(d)+'.png'
                plt.savefig(outfile,dpi=270)

        st_ocean = depth
        xt_ocean = x
        yt_ocean = y

        ###############################################################
        # save w
        file_date = jd2gcal(1721423,first_jdate+5*(countu5day+countfiles*73))
        out_filename_w = path2save+'ocean.'+str(file_date[0]).zfill(4)+str(file_date[1]).zfill(2)+str(file_date[2]).zfill(2)+'.w.nc'
        print 'Saving calculation in: '+out_filename_w
        netcdf_file = nc.Dataset(out_filename_w,'w',format='NETCDF4')
        netcdf_file.description = exp
            
        # dimensions
        netcdf_file.createDimension('time', None)
        netcdf_file.createDimension('st_ocean', len(st_ocean))
        netcdf_file.createDimension('yt_ocean', len(yt_ocean))
        netcdf_file.createDimension('xt_ocean', len(xt_ocean))
        
        # variables
        ti = netcdf_file.createVariable('time', 'f4', ('time',))
        ti.units = 'days since 0001-01-01 00:00:00'
        ti.calender_type = 'JULIAN'
        ti.calender = 'JULIAN'
        
        st = netcdf_file.createVariable('st_ocean', 'f4', ('st_ocean',))
        st.units = 'metres'
        st.long_name = 'tcell zstar depth'
        st.positive = 'down'
        
        yt = netcdf_file.createVariable('yt_ocean', 'f4', ('yt_ocean',))
        yt.units = 'degrees_N'
        yt.long_name = 'tcell latitude'
        
        xt = netcdf_file.createVariable('xt_ocean', 'f4', ('xt_ocean',))
        xt.units = 'degrees_E'
        xt.long_name = 'tcell longitude'
        
        dt = netcdf_file.createVariable('average_DT', 'f4', ('time',))
        dt.units = 'days'
        dt.long_name = 'length of average period'
        
        w_var = netcdf_file.createVariable('w', 'f4', ('time','st_ocean','yt_ocean', 'xt_ocean',),fill_value=-1.e20)
        w_var.units = 'm/s'
        w_var.long_name = 'vertical velocity'
        w_var.missing_value =  -1.e+20
        
        # data
        ti[:] = tim
        st[:] = st_ocean
        yt[:] = yt_ocean
        xt[:] = xt_ocean
        #u_var[:,:,:,0] = np.swapaxes(u,2,0)
        w_var[0,:] = w
        dt[:] = average_DT
        
        netcdf_file.close()
        ####################################################################
    countfiles =+1

