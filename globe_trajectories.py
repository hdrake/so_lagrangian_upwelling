import sys
from perceptually_uniform import *

## Parameters
savetype='movie'
steps_per_frame = 2
first_frame = 1000
nframes = 200

## Code
import numpy as np

## Load in subsampled bathymetry from ETOPO1.
data = np.load('/Users/hdrake/Dropbox/code/papers/timescale_final/old/trajectory_movie/bathymetry/ETOPO1_30min.npz')
c = np.swapaxes(data['z'][:,:],0,1)

phi = ((data['y'][:]*np.pi*2)/360.)+np.pi/2.

theta = (data['x'][::-1]*np.pi*2)/360.

phi = np.append(phi,[np.pi],axis=0)
theta = np.append(theta,[theta[0]],axis=0)
appendc = c[0,:]
c = np.append(c, appendc[None,:], axis=0)
appendc = c[:,-1]
c = np.append(c, appendc[:,None], axis=1)

phi, theta = np.meshgrid(phi,theta)

# Scale back height of mountains for clarity
c[c>-30]=c[c>-30]/3

# Create variable dimensions
x = np.sin(phi) * np.cos(theta) * (1 + c/30000.)
y = np.sin(phi) * np.sin(theta) * (1 + c/30000.)
z = np.cos(phi) * (1 + c/30000.)

from mayavi import mlab
from perceptually_uniform import plasma

# Loop through frames (only one if 'figure')
mast_randvec = np.random.random(10000)

mlab.figure(1,size = (960,575),bgcolor = (1,1,1), fgcolor = (0.5, 0.5, 0.5))
# Create figure, specifying size, background color, and foreground color
cocean = np.copy(c)
cland = np.copy(c)
cocean[cocean>=0] = np.nan
cland[cland<0] = np.nan

mlab.clf()

# Plot the surface of our globe, with colors specified by schocastic zcol.
cbdata = mlab.mesh(x, y, z, scalars = -cocean, colormap = 'Blues', vmin = 0, vmax = 2500)
lut = cbdata.module_manager.scalar_lut_manager.lut.table.to_array()
plas = plasma_r()
for k in range(0,256):
    color_tmp = 255*np.array(plas(k))
    lut[255-k,:] = np.array(color_tmp.astype(int))
cbdata.module_manager.scalar_lut_manager.lut.table = lut[::-1]

cb = mlab.scalarbar(cbdata,nb_labels = 11, orientation = 'vertical', label_fmt = '%.0f')

cbdata.module_manager.scalar_lut_manager.scalar_bar_representation.minimum_size = np.array([1, 1])
cbdata.module_manager.scalar_lut_manager.scalar_bar_representation.moving = 0
cbdata.module_manager.scalar_lut_manager.scalar_bar_representation.position2 = np.array([ 0.1,  0.8])
cbdata.module_manager.scalar_lut_manager.scalar_bar_representation.position = np.array([ 0.05,  0.05])

cbdata.remove()

if savetype=='figure':
    mlab.clf()

# Plot the surface of our globe, with colors specified by schocastic zcol.
m = mlab.mesh(x, y, z, scalars = -cocean, colormap = 'Blues', vmin = 0, vmax = 7000)
m.module_manager.scalar_lut_manager.lut.nan_color = [0,0,0,0]

# Plot the surface of our globe, with colors specified by schocastic zcol.
m = mlab.mesh(x, y, z, scalars = cland, colormap = 'gist_earth', vmin = -1800, vmax = 1800)
m.module_manager.scalar_lut_manager.lut.nan_color = [0,0,0,0]

atext = mlab.text(0.05,0.87,'Particle depth (m)')
atext.actor.text_scale_mode = 'none'
atext.property.font_size = 22
atext.property.bold = True

for frame in range(first_frame,first_frame+nframes):
    month = int(np.floor((frame*steps_per_frame*5)/12)+1)
    year = (frame*steps_per_frame*5)/365
    trajlist = [None]*6;
    for j in range(1,7):
        # Load in particle trajectories
        path2trajs = '/Users/hdrake/Dropbox/data/upwelling_trajectories/'
        # Load trajectories
        trajs = np.load(path2trajs+'/CM2.6/traj_file_'+str(j)+'.npz')
    
        # Extract variables
        traj_lons = trajs['lon']
        traj_lats = trajs['lat']
        traj_depths = -trajs['dep']
        npt = np.shape(traj_lons)[0]
        ntime = np.shape(traj_lons)[1]
        
        # Which particle?
        partn = 0
        nplot = 1
        
        upwell_time = np.zeros((npt))
        for pt in range(npt):
            if np.size(np.where(np.logical_and((traj_depths[pt,:] > -200),(traj_depths[pt,:] < 0)))):
                upwell_time[pt] = np.min(np.where(np.logical_and((traj_depths[pt,:] > -200),(traj_depths[pt,:] < 0))))
        
        inds = np.where(np.logical_and(upwell_time>0,np.logical_and(traj_depths[:,0]<-2000,traj_depths[:,0]>-2500)))[0]
        traj_lons = traj_lons[inds,:]
        traj_lats = traj_lats[inds,:]
        traj_depths = traj_depths[inds,:]
        upwell_time = upwell_time[inds]
        npt = np.shape(traj_lons)[0]
        ntime = np.shape(traj_lons)[1]
        
        rand_vec = mast_randvec[:npt]
        inds = np.argsort(rand_vec)
        traj_lons = traj_lons[inds,:]
        traj_lats = traj_lats[inds,:]
        traj_depths = traj_depths[inds,:]
        upwell_time = upwell_time[inds]
        
        trajlist_tmp = [None]*npt
        for i in range(partn,partn+npt):
            
            # Plot particle trajectory, cutting it off when it reaches the mixed layer
            # Convert 3D spherical coordinates in degrees to 2D spherical coordinate,
            # with theta in rad and r in latitude degrees from south pole.
            
            phi_traj = ((traj_lats[i,:upwell_time[i]]*np.pi*2)/360.)+np.pi/2.
    
            theta_traj = -(traj_lons[i,:upwell_time[i]]*np.pi*2)/360.
            
            delta_r_traj = np.copy(traj_depths[i,:upwell_time[i]])
            
            # Create variable dimensions
            x_traj = np.sin(phi_traj) * np.cos(theta_traj) * (1 + delta_r_traj/40000.)
            y_traj = np.sin(phi_traj) * np.sin(theta_traj) * (1 + delta_r_traj/40000.)
            z_traj = np.cos(phi_traj) * (1 + delta_r_traj/40000.)
            
            # Plot trajectory of particle
            if frame*steps_per_frame <= upwell_time[i]:
                p = mlab.plot3d(x_traj[0:frame*steps_per_frame],
                y_traj[0:frame*steps_per_frame],
                z_traj[0:frame*steps_per_frame],
                delta_r_traj[0:frame*steps_per_frame],
                tube_radius = 0.003,colormap = 'hot',vmin = -2500, vmax = 0, opacity = 1)
                lut = p.module_manager.scalar_lut_manager.lut.table.to_array()
                plas = plasma_r()
                for k in range(0,256):
                    color_tmp = 255*np.array(plas(k))
                    lut[255-k,:] = np.array(color_tmp.astype(int))
                p.module_manager.scalar_lut_manager.lut.table = lut
            else:
                p = mlab.plot3d(x_traj[0:upwell_time[i]],
                y_traj[0:upwell_time[i]],
                z_traj[0:upwell_time[i]],
                delta_r_traj[0:upwell_time[i]],
                tube_radius = 0.003,colormap = 'hot',vmin = -2500, vmax = 0, opacity = 1)
                lut = p.module_manager.scalar_lut_manager.lut.table.to_array()
                plas = plasma_r()
                for k in range(0,256):
                    color_tmp = 255*np.array(plas(k))
                    lut[255-k,:] = np.array(color_tmp.astype(int))
                p.module_manager.scalar_lut_manager.lut.table = lut
            trajlist_tmp[i] = p
        trajlist[j-1] = trajlist_tmp
    # Show plot
    mlab.draw()
    mlab.show()
    mlab.view(-10,0,4.6,np.array([0,0,0]))
    
    month = int(np.floor(((frame*steps_per_frame*5)%365)/30.5)+1)
    year = (frame*steps_per_frame*5)/365
    tit = mlab.title(str(month).zfill(2)+' / '+str(year).zfill(2)+' since release',height=0.935,size=0.25)
    
    # Save figure
    if savetype == 'movie':
        mlab.savefig('/Users/hdrake/Dropbox/temporary/frames/'+str(frame).zfill(4)+'.png')
        tit.remove()
    else:
        atext.remove()
        cb.remove()
        tit.remove()
        mlab.savefig('global_trajectories.png')
    
    # Reset trajectories
    for j in range(len(trajlist)):
        for i in range(len(trajlist[j])):
            trajlist[j][i].remove()
    