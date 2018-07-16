# Test upwelling_plots.py
#from upwelling_plots import *
import numpy as np
import netCDF4 as nc
import sys, os
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.colorbar as cb
from scipy.optimize import leastsq
from perceptually_uniform import *
from inversegaussian import *

mindep = 1000
maxdep = 4000

numfolders = 48
numfiles = 6
numdiv = 5
mode = 'MLD'

cms_exp_legend = ['CM2-1deg-5day (1$^{\circ}$)','CM2.5-5day (0.25$^{\circ}$)','CM2.6-5day (0.1$^{\circ}$)', 'CM2.6-Monthly (0.1$^{\circ}$)','CM2-1deg-nomeso-5day (1$^{\circ}$)']
cms_exps = ['cm21_GM_1860_500year_5day','cm25_1860_300year_5day','cm26_1860_200year_5day','cm26_1860_200year_monthly','cm21_noGM_1860_500year_5day']
models = ['CM2-1deg','CM2.5','CM2.6','CM2.6','CM2-1deg']
subplot_legend = ['a)','b)','c)','d)','e)','f)']
model_linestyles = ['dotted','dashed','solid','solid','dotted']
linestyles = ['solid','dashed','solid','dashed','solid','dashed']
loops_linestyles = ['solid','solid','solid','dashed','dashed']
model_colors = [(0,0,0),(0,0,0),(0,0,0),(1,0,0),(1,0,0)]
cmap = plasma()
facecolors = [cmap(0),cmap(0),cmap(90),cmap(90),cmap(180),cmap(180)]
labels = ['West Atlantic','East Atlantic','West Indian','East Indian','West Pacific','East Pacific','Total']
long_labels = ['West Atlantic Ocean', 'East Atlantic Ocean', 'West Indian Ocean', 'East Indian Ocean', 'West Pacific Ocean', 'East Pacific Ocean' ,'total']
mindeps = [1000,1000,1000,1000,1000]
maxdeps = [3500,3500,3500,3500,3500]
num_exp = len(cms_exps)

figuresize1 = (20,8)
figuresize2 = (10,8)
figuresize3 = (16,8)

fig1 = plt.figure(1,figsize=figuresize1)
fig1.clf()
fig2 = plt.figure(2,figsize=figuresize2)
fig2.clf()
fig3 = plt.figure(3,figsize=figuresize3)
fig3.clf()

max_ttd = np.zeros(num_exp)
west_ttd_transport = np.zeros(num_exp)
east_ttd_transport = np.zeros(num_exp)
total_ttd_transport = np.zeros(num_exp)
west_inds = np.array([0,2,4])
east_inds = np.array([1,3,5])

upwell_ttd_mode = np.zeros(num_exp)
upwell_ttd_t80 = np.zeros(num_exp)
upwell_ttd_t50 = np.zeros(num_exp)
upwell_ttd_t10 = np.zeros(num_exp)
lagrangian_transport = np.zeros((num_exp,7))
circumpolar_ttd_mode = np.zeros(num_exp)
circumpolar_num_mode = np.zeros(num_exp)
circumpolar_num_mean = np.zeros(num_exp)
ratio_of_modes = np.zeros(num_exp)

for e in range(0,num_exp):
    cms_exp = cms_exps[e]
    model = models[e]
    save_directory = '/work/hfd/figures/papers/upwelling_timescales_'+mode+'/'
    timestep = 5
    
    print 'Measuring trajectory data for '+cms_exp
    npt = 0
    for i in range(1,numfolders+1):
        print 'Measuring file '+str(i).zfill(2)
        path2trajs = '/archive/hfd/CMS/expt_'+cms_exp+'_'+str(i).zfill(2)+'/output/'
        first = 1
        for j in range(1,numfiles+1):
            for div in range(1,numdiv+1):
                filename = path2trajs+mode+'_upwell_trajs_'+str(i).zfill(2)+'_'+str(j).zfill(2)+'_'+str(div).zfill(2)+'.npz'
                if not(os.path.isfile(filename)):
                    if first == 1:
                        print 'Loading file '+str(i).zfill(2)+' failed'
                        first = 0
                else:
                    npt_step = np.shape(np.load(filename)['upwell_traj_lat'])[0]
                    ntime = np.shape(np.load(filename)['upwell_traj_lat'])[1]
                    npt += npt_step
    
    # Quick summary of dataset
    traj_lons = np.zeros((npt,ntime))
    traj_lats = np.zeros((npt,ntime))
    traj_depths = np.zeros((npt,ntime))
    traj_MLD = np.zeros((npt,ntime))
    
    locations = np.zeros((npt))
    release_dates = np.zeros((npt))
    transports = np.zeros((npt))
    d200_time = np.zeros((npt))
    
    npt_count = 0
    print 'Loading trajectory data for '+cms_exp
    # Go through 48 different releaseFiles
    for i in range(1,numfolders+1):
        print 'Loading file '+str(i).zfill(2)
        path2trajs = '/archive/hfd/CMS/expt_'+cms_exp+'_'+str(i).zfill(2)+'/output/'
        first = 1
        # Go through 6 trajectory files
        for j in range(1,numfiles+1):
            for div in range(1,numdiv+1):
                filename = path2trajs+mode+'_upwell_trajs_'+str(i).zfill(2)+'_'+str(j).zfill(2)+'_'+str(div).zfill(2)+'.npz'
                if not(os.path.isfile(filename)):
                    if first == 1:
                        print 'Loading file '+str(i).zfill(2)+' failed'
                        first = 0
                else:
                    trajs = np.load(filename)
    
                    # Extract variables
                    traj_lats_tmp = trajs['upwell_traj_lat'][:,:]
                    traj_lons_tmp = trajs['upwell_traj_lon'][:,:]
                    traj_depths_tmp = trajs['upwell_traj_depth'][:,:]
                    
                    # Extract variables
                    d200_time_tmp = trajs['upwell_time'][:]
                    release_dates_tmp = trajs['upwell_release_date'][:]
                    transports_tmp = trajs['upwell_transport'][:]
                    npt_step = np.size(d200_time_tmp)
                    
                    trajs = None
                    
                    # Add to total trajectory files
                    traj_lons[npt_count:npt_count+npt_step,:] = traj_lons_tmp
                    traj_lats[npt_count:npt_count+npt_step,:] = traj_lats_tmp
                    traj_depths[npt_count:npt_count+npt_step,:] = traj_depths_tmp

                    release_dates[npt_count:npt_count+npt_step] = release_dates_tmp
                    transports[npt_count:npt_count+npt_step] = transports_tmp
                    d200_time[npt_count:npt_count+npt_step] = d200_time_tmp
                    
                    npt_count += npt_step

    # Choose depths
    dep_range_inds = np.transpose(np.where(np.logical_and(traj_depths[:,0] >= mindep, traj_depths[:,0] < maxdep)))[:,0]
    longitude = traj_lons[dep_range_inds,:]
    traj_lons = None
    latitude = traj_lats[dep_range_inds,:]
    trajs_lats = None
    depth = traj_depths[dep_range_inds,:]
    traj_depths = None
    MLD = traj_MLD[dep_range_inds,:]
    traj_MLD = None

    locations = locations[dep_range_inds]
    cutoff_time_index = d200_time[dep_range_inds].astype(int)
    transport = -transports[dep_range_inds]/np.size(np.unique(release_dates))

    htfilename = '/archive/Carolina.Dufour/Mosaic/'+model+'/ocean.static.nc'
    htfile = nc.Dataset(htfilename)
    yt_ocean = htfile.variables['yt_ocean'][:]
    xt_ocean = htfile.variables['xt_ocean'][:]
    ocean_depth = htfile.variables['ht'][:,:]


    # Number of particles in trajectories
    nparticles = np.shape(longitude)[0]
    # Number of time steps in trajectories
    ntime = np.shape(longitude)[1]
    # Number of trajectory dimensions (usually 3)
    ndimensions = 3
    # Number of release basins
    nbasins = 6

    # Get locations (Longitude, Latidude, Depth) that particles cross mixed layer (and various depth surfaces)
    cutoff_location = np.zeros((nparticles,ndimensions))
    cutoff_MLD = np.zeros(nparticles)
    for p in range(nparticles):
        cutoff_location[p,0] = longitude[p,cutoff_time_index[p]]
        cutoff_location[p,1] = latitude[p,cutoff_time_index[p]]
        cutoff_location[p,2] = depth[p,cutoff_time_index[p]]
        cutoff_MLD[p] = MLD[p,cutoff_time_index[p]]

    # Separate particle indices based on release-basin (Atlantic, Indian, Pacific)
    basin_inds = np.zeros((nparticles,nbasins)).astype(int)
    basin_inds[np.logical_and((longitude[:,0]>=-65),(longitude[:,0]<-15)),0] = 1 # West Atlantic
    basin_inds[np.logical_and((longitude[:,0]>=-15),(longitude[:,0]<20)),1] = 2 # East Atlantic

    basin_inds[np.logical_and((longitude[:,0]>=25),(longitude[:,0]<=70)),2] = 3 # West Indian
    basin_inds[np.logical_or(np.logical_and((longitude[:,0]>=70),(longitude[:,0]<=90)),np.logical_and((longitude[:,0]>=-270),(longitude[:,0]<-240))),3] = 4 # East Indian

    basin_inds[np.logical_and((longitude[:,0]>=-210),(longitude[:,0]<-130)),4] = 5 # West Pacific
    basin_inds[np.logical_and((longitude[:,0]>=-130),(longitude[:,0]<-70)),5] = 6 # East Pacific

    # Separate transport, cutoff_time_index, and cutoff_location by basin
    total_transport = np.sum(transport)
    basin_transport = np.zeros((nparticles,nbasins))
    basin_cutoff_time_index = np.zeros((nparticles,nbasins))
    basin_cutoff_location = np.zeros((nparticles,nbasins,ndimensions))+10000
    for b in range(nbasins):
        basin_transport[basin_inds[:,b] == b+1,b] = transport[basin_inds[:,b] == b+1]
        basin_cutoff_time_index[basin_inds[:,b] == b+1,b] = cutoff_time_index[basin_inds[:,b] == b+1]
        for d in range(ndimensions):
            basin_cutoff_location[basin_inds[:,b] == b+1,b,d] = cutoff_location[basin_inds[:,b] == b+1,d]

    print 'There are '+str(int(nparticles))+' particles that are released between '+str(mindep)+' and '+str(maxdep)+' and eventually reach the mixed layer.'
    
    import matplotlib
    matplotlib.rcParams.update({'font.size': 15})
    
    #### FIGURE 1 #####
    # Bin sizes
    if e < 3:
        delta_longitude_bin = 1.0
        delta_latitude_bin = 1.0
    
        # Dependent Parameters
        longitude_bin = np.arange(-280,80+delta_longitude_bin,delta_longitude_bin)
        latitude_bin = np.arange(-80,-29,delta_latitude_bin)
    
        # Number of bins
        numlon = np.size(longitude_bin)
        numlat = np.size(latitude_bin)
    
        # Attempt Veronica's Method
        percent_visited_cumulative = np.zeros((numlat,numlon,nbasins))
        percent_visited_cumulative_basin_normalized = np.zeros((numlat,numlon,nbasins))
        nbin = [numlat,numlon]
        binlim = [[-80,-29],[-280,80]]
    
        out_of_bounds_mask = longitude > 1000
        particle_mask = np.zeros((np.shape(longitude)),dtype=bool)
        for p in range(nparticles):
            particle_mask[p,cutoff_time_index[p]:] = True
            
        masking = np.logical_and(out_of_bounds_mask,particle_mask)
        longitude = np.ma.masked_array(longitude,mask = masking)
        latitude = np.ma.masked_array(latitude,mask = masking)
    
        for b in range(nbasins):
            in_box = np.zeros((numlat,numlon))
            for p in range(0,nparticles):
                if basin_inds[p,b]:
                    H = np.histogram2d(latitude[p,:cutoff_time_index[p]],longitude[p,:cutoff_time_index[p]],nbin,binlim)
                    boxind = np.nonzero(H[0])
                    in_box[boxind[0],boxind[1]] = in_box[boxind[0],boxind[1]]+transport[p]
            percent_visited_cumulative[:,:,b] = 100.*np.float64(in_box)/np.float64(total_transport)
            percent_visited_cumulative_basin_normalized[:,:,b] = 100.*np.float64(in_box)/np.float64(np.sum(basin_transport[:,b]))
    
        # Loop the data back around so that the polar plots close nicely
        percent_visited_cumulative[:,-1,:] = percent_visited_cumulative[:,0,:]
        percent_visited_cumulative_basin_normalized[:,-1,:] = percent_visited_cumulative_basin_normalized[:,0,:]
    
        # Plot percent of particles visiting each grid column
        maxlevel = 15
        minlevel = 1.5
        deltalevel = 0.1
        levels = np.arange(minlevel,maxlevel+deltalevel,deltalevel)
        cmap = plasma_r()
        longitude_bin[-1]=-279.5
    
        # Figure
        plt.figure(1,figsize=figuresize1)
    
        # Subplots
        ax = plt.subplot(1,3,1+e)    
        
        # Plot filled contour in polar projection
        m = Basemap(projection = 'spstere',boundinglat = -30,lon_0=180,resolution='l',round='True')
        x,y = m(*np.meshgrid(longitude_bin+delta_longitude_bin/2.0,latitude_bin+delta_latitude_bin/2.0))
        z = np.sum(percent_visited_cumulative[:,:,:],axis=2)
        CS = m.contourf(x,y,z,levels=levels,extend="max",cmap=cmap)
        m.drawcoastlines()
        m.fillcontinents(color='lightgrey',lake_color='aqua')
        m.drawparallels(np.arange(-80.,81.,10.))
        ax.set_title(subplot_legend[e]+' '+cms_exp_legend[e],y=1.02)
        fig1.subplots_adjust(top = 0.85)
        cbar_ax = fig1.add_axes([0.2, 0.92, 0.6, 0.03])
        cb1 = fig1.colorbar(CS, cax = cbar_ax, orientation = 'horizontal',ticks=[1.5, 3.0,4.5,6.0,7.5,9.0,10.5,12.0,13.5,15.0])
    
        if e == 2:
            plt.suptitle('Percent of particle-transport visiting each grid column \nbetween release at 30S and the surface mixed layer.', fontsize=18)
            outfile = save_directory+'/fig1_pathways.png'
            plt.savefig(outfile,bbox_inches='tight',dpi=270)

    ####################
    #### Figure 2 ######

    # Size of bins for histograms
    delta_time_bin = int(365/timestep)

    # Bins for histograms
    time_bin = np.arange(0,ntime,delta_time_bin)

    # Plotting time of upwelling with exponent tail'
    print 'Plot probability of upwelling as a function of time, with exponential tail'
    # Create one large array that contains all of the histogram information.
    # Dimensions are: (Number of bins, number of release basins)
    time_histograms = np.zeros((np.size(time_bin)-1,nbasins+1))
    for b in range(nbasins):
        for t in range(np.size(time_bin)-1):
            time_histograms[t,b] = np.sum((basin_transport[np.logical_and(basin_cutoff_time_index[:,b] > time_bin[t],basin_cutoff_time_index[:,b] <= time_bin[t+1]),b])/(total_transport))

    time_histograms[:,nbasins] = np.sum(time_histograms[:,0:nbasins],axis=1)

    # X in timestep increments
    x = time_bin[:-1]
    x_long = np.arange(delta_time_bin, ntime*3, delta_time_bin)

    # Convert X to years and center points on middle of timesteps
    x_years = (x + delta_time_bin/2.0) * timestep/(365.)
    x_years_long = (x_long + delta_time_bin/2.0) * timestep/(365.)
        
    # Plot exponential tail
    fig2 = plt.figure(2,figsize=figuresize2)

    y_est_long_total = np.zeros((np.size(x_years_long),7))
    for b in range(0,nbasins+1):
        # Data to fit the inverse gaussian distribution to
        y_data = time_histograms[:,b]
        
        # Initial guess for Gam, Del, Sca
        Gam, Del, Sca = [100, 100, 0.5]
        p = [Gam, Del, Sca] # Initial guesses for leastsq
        y_init = inversegaussian(x_years, Gam, Del, Sca) # Initial Guess
        
        # Least squares
        plsq = leastsq(res, p, args = (y_data, x_years))
        
        if b<nbasins:
            lagrangian_transport[e,b] = plsq[0][2]*total_transport
        else:
            lagrangian_transport[e,b] = np.sum(lagrangian_transport[e,:b])

        # Fitted inverse gaussian distribution, optimizing three parameters
        if b<nbasins:
            y_est_long_total[:,b] = inversegaussian(x_years_long, plsq[0][0], plsq[0][1], plsq[0][2])
        else:
            y_est_long_total[:,b] = np.sum(y_est_long_total[:,0:nbasins],axis=1)
    
        
        # Fitted inverse gaussian distribution, optimizing three parameters
        y_est_long = y_est_long_total[:,b]
        
        # Index of max of the distribution
        argmax_dist = np.argmax(y_est_long)
        
        if b == nbasins:
            #plt.title(subplot_legend[e]+' '+cms_exp_legend[e],y=1.02, loc='left')
            total_ttd_transport[e] = west_ttd_transport[e]+east_ttd_transport[e]
        else:
            if not(b%2):
                west_ttd_transport[e] += plsq[0][2]*total_transport
            else:
                east_ttd_transport[e] += plsq[0][2]*total_transport
            #plt.plot(x_years_long, y_est_long*total_transport, color = facecolors[b], linestyle=linestyles[b], label=labels[b], linewidth = 1.5)
        
        #plt.axis([0,400,0,np.max(y_est_long_total[:,:nbasins]*total_transport)*1.05])
        #if e==0:
            #plt.legend(loc='upper right',fontsize=13)
        #plt.xlabel('Time (years)')        
        
    # Plot the model's total line
    max_ttd[e] = np.max(y_est_long*total_transport)*1.05
    upwell_ttd_mode[e] = x_years_long[argmax_dist]
    upwell_ttd_t80[e] = x_years_long[np.argmax(np.cumsum(y_est_long*total_transport)>0.8*total_transport)]
    upwell_ttd_t50[e] = x_years_long[np.argmax(np.cumsum(y_est_long*total_transport)>0.5*total_transport)]
    upwell_ttd_t10[e] = x_years_long[np.argmax(np.cumsum(y_est_long*total_transport)>0.1*total_transport)]
    
    plt.plot(x_years_long, y_est_long*total_transport, color = model_colors[e], linestyle = model_linestyles[e], label=cms_exp_legend[e], linewidth = 1.5)
    plt.legend(loc='upper right',fontsize=14)
    #plt.title(subplot_legend[num_exp]+' All experiments', y=1.02, loc='left')
    plt.xlabel('Time (years)')
    plt.axis([0,400,0,np.max(max_ttd)])
    
    if e == num_exp-1:
        plt.suptitle('Transit-Time Distributions of upwelling particle-transport', y=0.98, fontsize=18)
        outfile = save_directory+'/fig2_TTDs.png'
        plt.tight_layout(rect=[0,0,1,0.96])
        plt.savefig(outfile,bbox_inches='tight',dpi=270)

    ######################
    #### TABLE 1 #########

    ######################
    ##### FIGURE 3 #######
    if e<3:
        plt.figure(3,figsize=figuresize3)
        
        # Find how much a particle changes longitude at every timestep
        delta_longitude = (longitude[:,1:]-longitude[:,:-1])
        
        # Need to address boundary conditions (what happens if a particle goes from 79.9 degrees to -279.5 degrees?)
        delta_longitude[np.logical_and(delta_longitude < -180, delta_longitude > - 360)] = delta_longitude[np.logical_and(delta_longitude < -180, delta_longitude > - 360)] + 360
        delta_longitude[np.logical_and(delta_longitude > 180, delta_longitude < 360)] = delta_longitude[np.logical_and(delta_longitude > 180, delta_longitude < 360)] - 360
        
        # Find the number of loops by doing a cumulative sum
        num_loops_cumulative = np.cumsum(delta_longitude, axis = 1)
        
        delta_longitude = None
        
        # Count the number of loops a particle does
        num_antarctica_loops = np.zeros(nparticles)
        for p in range(nparticles):
            num_antarctica_loops[p] = num_loops_cumulative[p,min(cutoff_time_index[p],ntime-2)]/360.0
        
        # Find the average period of one loop of a particle
        acc_loop_period = (((cutoff_time_index)*timestep)/num_antarctica_loops)/365.0
        
        
        # Subplots
        ax = plt.subplot(1,2,2)
        bins = np.arange(0,8,0.25)
        y = np.histogram(num_antarctica_loops, bins=bins, weights=transport, normed = True)[0]
        bins = np.append([0],bins,axis=0)
        y = np.append([0],y,axis=0)
        plt.plot((bins[:-1]+bins[1:])/2., y, color = model_colors[e], label=cms_exp_legend[e], linestyle = model_linestyles[e], linewidth = 1.5)
        plt.legend(loc='upper right',fontsize=13)
        print cms_exp_legend[e]+' '+str(np.ma.average(num_antarctica_loops,weights=transport))
        x_axis_bins = (bins[:-1]+bins[1:])/2.
        circumpolar_num_mode[e] = x_axis_bins[np.argmax(y)]
        circumpolar_num_mean[e] = np.sum(num_antarctica_loops * transport)/total_transport
        
        ax = plt.subplot(1,2,1)
        bins = np.arange(0,255,5)
        y = np.histogram(acc_loop_period, bins=bins, weights=transport, normed = True)[0]
        bins = np.append([0],bins,axis=0)
        y = np.append([0],y,axis=0)
        plt.plot((bins[:-1]+bins[1:])/2., y, color = model_colors[e], label=cms_exp_legend[e], linestyle = model_linestyles[e], linewidth = 1.5)
        acc_loop = np.ma.masked_array(acc_loop_period, mask = np.abs(acc_loop_period)>500)
        print cms_exp_legend[e]+' '+str(np.ma.average(acc_loop,weights=transport))
        x_axis_bins = (bins[:-1]+bins[1:])/2.
        circumpolar_ttd_mode[e] = x_axis_bins[np.argmax(y)]
        
        ratio_of_modes[e] = upwell_ttd_mode[e] / circumpolar_ttd_mode[e]
        
        if e == 2:
            # Make the plots.
            plt.subplot(1,2,2)
            plt.title('b) Number of circumpolar loops before \nreaching the surface mixed layer',y=1.015)
            plt.xlabel('Number of loops')
            plt.subplot(1,2,1)
            plt.title('a) Time for one circumpolar loop',y=1.015)
            plt.xlabel('Years')
            plt.savefig(save_directory+'/fig3_loops.png')
            #plt.draw()
            #plt.show()

        
print 'Lagrangian Transport: '
print lagrangian_transport
print ' '
print 'Upwell ttd mode '+str(upwell_ttd_mode)
print 'Upwell ttd t80% '+str(upwell_ttd_t80)
print 'Upwell ttd t50% '+str(upwell_ttd_t50)
print 'Upwell ttd t10% '+str(upwell_ttd_t10)
print 'Circumpolar num mean '+str(circumpolar_num_mean)
print 'Ratio of ttd modes '+str(ratio_of_modes)

