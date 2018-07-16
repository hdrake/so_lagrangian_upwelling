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
path2save = '/archive/Henri.Drake/'+model+'/'+exp+'5day/for_CMS_SO/'
archdir_hist = '/archive/Richard.Slater/'+model+'/'+exp+'/history/'
other_dir = '/archive/wga/'+model+'/'+exp+'pp/ocean_budgets/ts/annual/5yr/'

# Get Grid Data
griddir = '/archive/Carolina.Dufour/Mosaic/'+model+'/'
grid_file = griddir+model+'_grid_spec.nc'
gridFile = nc.Dataset(grid_file)
lat_t = gridFile.variables['geolat_t'][...]

first_year = 181
last_year = 200
boundary_north = -29.5
index_north = np.where(lat_t[:,0]<boundary_north)[0][-1]

dxt = gridFile.variables['dxt'][:index_north,:]
dyt = gridFile.variables['dyt'][:index_north,:]

files = os.listdir(path2save)
files = np.sort(files)

interval_list = ['ocean_budgets.018101-018501.',
                 'ocean_budgets.018601-019001.',
                 'ocean_budgets.019101-019501.',
                 'ocean_budgets.019601-020001.']

countu5day = 0
countfiles = 0
#### COMPUTE uhrho_et and vhrho_nt 5-day averages from u,v,rho_dzt 5-day averages

#935 lat indices for north = -29.5
uhrho_et_year = np.zeros((73,50,935,3600))
vhrho_nt_year = np.zeros((73,50,935,3600))
# u and v
for file in files:
    if not file.endswith('u.nc'):
        continue
    print 'Doing '+file
    year = file.split('.',3)[1][0:4]
    file_name_prefix = str(year)+'0101.ocean_minibling_field_'
    if (year == str(first_year).zfill(4)) and (countu5day == 0):
        # Get landmask from salt field
        data = nc.Dataset(archdir_hist + file_name_prefix + 'salt.nc')
        salt = data.variables['salt'][0,:,:index_north,:]
    # rho_dzt is the grid cell thickness (which varies in time) x 1035
    data = nc.Dataset(archdir_hist + file_name_prefix + 'rho_dzt.nc')
    # use filling trick to avoid spreading mask:
    rho_dzt = (data.variables['rho_dzt'][countu5day,:,:index_north+1,:])\
        .filled(fill_value=1.e20)/1035.
    # compute rho_dzu (depth, lat, lon) by taking minimum of four grid cell corners
    # Take out last latitude (most northern?) and roll longitudes over
    rho_dzt_right = np.roll(rho_dzt[:,:-1,:],-1,axis=-1)
    # Take out first latitude
    rho_dzt_up = rho_dzt[:,1:,:]
    # Take out first latitude and roll longitudes over
    rho_dzt_up_right = np.roll(rho_dzt[:,1:,:],-1,axis=-1)
    # Take minimum of thickness of tracer grid box over 4 corners to be the grid thickness
    rho_dzu = np.minimum(np.minimum(np.minimum(rho_dzt[:,:-1,:],rho_dzt_right),rho_dzt_up),rho_dzt_up_right)
    # Mask if ground or if 0 thickness
    rho_dzu = np.ma.masked_array(rho_dzu,mask=salt.mask)
    rho_dzu = rho_dzu.filled(fill_value=0)
    del rho_dzt_right,rho_dzt_up,rho_dzt_up_right

    print 'calculating vhrho_nt and uhrho_nt'
    data = nc.Dataset(archdir_hist + file_name_prefix + 'v.nc')
    v = (data.variables['v'][countu5day,:,:index_north,:]).filled(fill_value=0)
    data = nc.Dataset(archdir_hist + file_name_prefix + 'u.nc')
    u = (data.variables['u'][countu5day,:,:index_north,:]).filled(fill_value=0)
    # compute vhrho_nt and uhrho_ut:
    vhrho_net = v*rho_dzu
    vhrho_nwt = np.roll(vhrho_net,1,axis=-1)
    uhrho_net = u*rho_dzu
    uhrho_set = u[:,1:,:]*rho_dzu[:,1:,:]
    uhrho_set = np.insert(uhrho_set,0,uhrho_set[:,0,:],axis=1)
    del u,v,rho_dzu
    # this is still on upper right corner of grid cell at this stage,
    # so interpolate to top center of grid cell:
    vhrho_nt = 0.5*(vhrho_net+vhrho_nwt)
    vhrho_nt = np.ma.array(vhrho_nt,mask=salt.mask)
    vhrho_nt = np.ma.filled(vhrho_nt,fill_value=0)

    uhrho_et = 0.5*(uhrho_net+uhrho_set)
    uhrho_et = np.ma.array(uhrho_et,mask=salt.mask)
    uhrho_et = np.ma.filled(uhrho_et,fill_value=0)

    uhrho_et_year[countu5day,:,:,:] = uhrho_et
    vhrho_nt_year[countu5day,:,:,:] = vhrho_nt

    ##### Annual Average and comparing to other calculations of vhrho_nt and uhrho_et
    if countu5day == 72:
        uhrho_et_annual = np.mean(uhrho_et_year, axis = 0)
        vhrho_nt_annual = np.mean(vhrho_nt_year, axis = 0)

        #### Archived Annual Averages
        data = nc.Dataset(other_dir + interval_list[countfiles/(73*5)] + 'uhrho_et.nc')
        other_uhrho_et = data.variables['uhrho_et'][countfiles/(73),:,:index_north,:]
        other_uhrho_et = np.ma.array(other_uhrho_et,mask=salt.mask)
        other_uhrho_et = np.ma.filled(other_uhrho_et,fill_value=0)/1035.
        lon = data.variables['xt_ocean'][:]
        lat = data.variables['yt_ocean'][:index_north]
        data = nc.Dataset(other_dir + interval_list[countfiles/(73*5)] + 'vhrho_nt.nc')
        other_vhrho_nt = data.variables['vhrho_nt'][countfiles/(73),:,:index_north,:]
        other_vhrho_nt = np.ma.array(other_vhrho_nt,mask=salt.mask)
        other_vhrho_nt = np.ma.filled(other_vhrho_nt,fill_value=0)/1035.
        
        uhnonzero_ind = np.where((uhrho_et_annual) != 0 & (other_uhrho_et != 0))
        vhnonzero_ind = np.where((vhrho_nt_annual) != 0 & (other_vhrho_nt != 0))

        print np.average(np.abs(other_uhrho_et))
        print np.average(np.abs(uhrho_et_annual))
        print np.average(np.abs(other_vhrho_nt))
        print np.average(np.abs(vhrho_nt_annual))
        print 'Non-zero including averages:'
        print np.average(np.abs(other_uhrho_et[uhnonzero_ind]))
        print np.average(np.abs(uhrho_et_annual[uhnonzero_ind]))
        print np.average(np.abs(other_vhrho_nt[vhnonzero_ind]))
        print np.average(np.abs(vhrho_nt_annual[vhnonzero_ind]))
        
        diffuh = np.abs(other_uhrho_et - uhrho_et_annual)
        diffvh = np.abs(other_vhrho_nt - vhrho_nt_annual)
        print 'Av diff: '+str(np.average(diffuh))+'  and  '+str(np.average(diffvh))
        print 'Std from diff: '+str(np.std(diffuh))+'  and  '+str(np.std(diffvh))
        print 'Av diff non-zero: '+str(np.average(diffuh[uhnonzero_ind]))+'  and  '+str(np.average(diffvh[vhnonzero_ind]))
        print 'Std from non-zero diff: '+str(np.std(diffuh[uhnonzero_ind]))+'  and  '+str(np.std(diffvh[vhnonzero_ind]))
        
        
        #raw_input("Enter")
        
        indices_other_uh = np.where(other_uhrho_et == 0)
        indices_other_vh = np.where(other_vhrho_nt == 0)
        #print 'Shapes: '+str(np.shape(other_uhrho_et))+' and '+str(np.shape(indices_other_vh))
        other_vhrho_nt[indices_other_vh] = 0
        other_uhrho_et[indices_other_uh] = 0
        percent_error_uh = np.abs((uhrho_et_annual - other_uhrho_et)/other_uhrho_et)
        percent_error_vh = np.abs((vhrho_nt_annual - other_vhrho_nt)/other_vhrho_nt)
        percent_error_uh[indices_other_uh] = 0
        percent_error_vh[indices_other_vh] = 0
        goodind_uh = np.where(percent_error_uh != 0)
        goodind_vh = np.where(percent_error_vh != 0)

        print 'Percent error of my uhrho_et calculation: '+str(np.average(percent_error_uh[goodind_uh]))
        print 'Percent error of my vhrho_nt calculation: '+str(np.average(percent_error_vh[goodind_vh]))
        print 'Max percent errors: '+str(percent_error_uh.max())+' and '+str(percent_error_vh.max())
        
        bad_indices_uh = np.where(percent_error_uh > 100.0)
        bad_indices_vh = np.where(percent_error_vh > 100.0)
        print 'This many bad points: '+str(len(bad_indices_uh[0]))
        print 'This many bad points: '+str(len(bad_indices_vh[0]))
        #raw_input("Press Enter")
        #print percent_error_uh[bad_indices_uh]
        #raw_input("Press Enter")
        #print np.abs(other_uhrho_et[bad_indices_uh] - uhrho_et_annual[bad_indices_uh])
        #raw_input("Press Enter")
        #print percent_error_vh[bad_indices_vh]
        #raw_input("Press Enter")
        
        x = lon
        y = lat
        z1 = np.mean(uhrho_et_annual, axis = 0)
        z2 = np.mean(vhrho_nt_annual, axis = 0)
        z3 = np.mean(other_uhrho_et, axis = 0)
        z4 = np.mean(other_vhrho_nt, axis = 0)
        
        zminu = min(z1.min(),z3.min())
        zmaxu = max(z1.max(),z3.max())
        zminv = min(z2.min(),z4.min())
        zmaxv = max(z2.max(),z4.max())
        
        fig1 = plt.figure(1)
        plt.clf()
        levels=LinearLocator(numticks=11).tick_values(zminu,zmaxu)
        cmap = plt.get_cmap('cool')
        norm = BoundaryNorm(levels,ncolors=cmap.N, clip=True)
        plt.contourf(x,y,z1,levels=levels,cmap=cmap)
        plt.colorbar()
        plt.axis([-280,80,-81,-30])
        plt.title('Averaged uhrho_nt from 5-day calculations')
        
        figfilepath = '/work/hfd/figures/CMS/misc/uhrho_et'
        plt.savefig(figfilepath,bbox_inches = 'tight',dpi = 270)
        
        fig2 = plt.figure(2)
        plt.clf()
        levels=LinearLocator(numticks=11).tick_values(zminv,zmaxv)
        cmap = plt.get_cmap('cool')
        plt.contourf(x,y,z2,cmap=cmap)
        plt.colorbar()
        plt.axis([-280,80,-81,-30])
        plt.title('Averaged vhrho_nt from 5-day calculations')
        
        figfilepath = '/work/hfd/figures/CMS/misc/vhrho_nt'
        plt.savefig(figfilepath,bbox_inches = 'tight',dpi = 270)
        
        fig3 = plt.figure(3)
        plt.clf()
        levels=LinearLocator(numticks=11).tick_values(zminu,zmaxu)
        cmap = plt.get_cmap('cool')
        norm = BoundaryNorm(levels,ncolors=cmap.N, clip=True)
        plt.contourf(x,y,z3,levels=levels,cmap=cmap)
        plt.colorbar()
        plt.axis([-280,80,-81,-30])
        plt.title('Archived uhrho_et')
        
        figfilepath = '/work/hfd/figures/CMS/misc/arch_uhrho_et'
        plt.savefig(figfilepath,bbox_inches = 'tight',dpi = 270)
        
        fig4 = plt.figure(4)
        plt.clf()
        levels=LinearLocator(numticks=11).tick_values(zminv,zmaxv)
        cmap = plt.get_cmap('cool')
        norm = BoundaryNorm(levels,ncolors=cmap.N, clip=True)
        plt.contourf(x,y,z4,levels=levels,cmap=cmap)
        plt.colorbar()
        plt.axis([-280,80,-81,-30])
        plt.title('Archived vhrho_nt')
        
        figfilepath = '/work/hfd/figures/CMS/misc/arch_vhrho_et'
        plt.savefig(figfilepath,bbox_inches = 'tight',dpi = 270)
        
        #######
#       fig3 = plt.figure(3)
#       plt.clf()
        
        ## Oops, can't give tuple indices
#       plt.plot(bad_indices_uh[:,1],bad_indices_vh[:,1],'r.')
#       plt.axis([-280,80,-81,-75])
#
#       figfilepath = '/work/hfd/figures/CMS/misc/vhuh_error_locations'
#plt.savefig(figfilepath,bbox_inches = 'tight',dpi = 270)

    ##### Counter
    countfiles += 1
    countu5day += 1
    if countu5day == 73:
        countu5day = 0