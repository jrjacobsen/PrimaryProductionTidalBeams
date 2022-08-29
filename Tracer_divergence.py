"""
Created on Tue May 17 11:47:50 2022

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')
import sys
sys.path.append('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')
sys.path.append('/home/jjacob2/.local/bin')
from netCDF4 import Dataset as nc4
import numpy as np
import obs_depth_vs145 as dep
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean

#file directories
root = '/home/jjacob2/runs/npzd_coupled/dimensional/olig_ridge04/output02/'

#netcdf directory and file for saving
n_file = 'roms_his.nc'


#load files
ncfile = nc4(root+n_file, 'r')


def nw_load(ncfile, idx) :

    nw ={}
    nw['nute'] = ncfile.variables['dye_01'][idx['time'], :,
                                           idx['lat'],
                                           idx['lon']]
    nw['w'] = ncfile.variables['w'][idx['time'], :,
                                           idx['lat'],
                                           idx['lon']]
    nw['u'] = ncfile.variables['u'][idx['time'], :,
                                           idx['lat'],
                                           idx['ulon']]

    return nw

def grid_load(FileName,ncfile, idx, cr_dist):
    #surface bounce distance
    N = 2e-3
    om = (2*np.pi)/(3600*12.4)
    theta_g = np.pi/2 - np.arccos(om/N)
    bot1 = 2/np.tan(theta_g)
    surf = 2*bot1

    #
    grid = {}
    #depth
    romsvars = {'Vstretching' : ncfile.variables['Vstretching'][0], \
                'Vtransform' : ncfile.variables['Vtransform'][0], \
                'theta_s' :  ncfile.variables['theta_s'][0], \
                'theta_b' : ncfile.variables['theta_b'][0], \
                'N' : ncfile.variables['Cs_r'].size, \
                'h' : ncfile.variables['h'][:], \
                'hc': ncfile.variables['hc'][0]}
    grid['depth'] =dep._set_depth(FileName, romsvars, 'rho',
                                 romsvars['h'])[:,
                                                idx['lat'],
                                                idx['lon']]
    grid['depthw'] =dep._set_depth(FileName, romsvars, 'w',
                                 romsvars['h'])[:,
                                                idx['lat'],
                                                idx['lon']]
    grid['nM2'] = (ncfile.variables['ocean_time'][idx['time']])/(12.4*3600)

    #distance from ridge [km]
    dist =(ncfile.variables['x_rho'][idx['lat'],
                                              idx['lon']]/1000-600)
    #center on critical point on eastern side of ridge
    cr_cent = dist - cr_dist

    #nondimensionalize by surface bounce distance
    grid['nbnc'] = cr_cent/surf

    return grid

nM2 = ncfile.variables['ocean_time'][:]/(12.4*3600)

#slicing
idx = {'time' : slice(149,None),
       'lat' : 3,
       'lon' : slice(228,438),
       'ulon': slice(228,439)}

hr04 = {}
hr04['nw'] = nw_load(ncfile, idx)
hr04['grid'] = grid_load(n_file, ncfile, idx, 28.2796)

#gradient div dot (U N)
#shift nutrient to u and w points
pad = np.concatenate((hr04['nw']['nute'][:,:,0:1],
                      hr04['nw']['nute'][:,:,:],
                      hr04['nw']['nute'][:,:,-2:-1]), axis = 2)
nute_u = 0.5*(pad[:,:,:-1]+pad[:,:,1:])
pad = np.concatenate((hr04['nw']['nute'][:,0:1,:],
                      hr04['nw']['nute'][:,:,:],
                      hr04['nw']['nute'][:,-2:-1,:]), axis = 1)
nute_w = 0.5*(pad[:,:-1,:]+ pad[:,1:,:])

#mean product
uN = np.mean(nute_u*hr04['nw']['u']*3600, axis = 0)
wN = np.mean(nute_w*hr04['nw']['w']*3600, axis = 0)

#gradient in x
xgrad = np.diff(uN, axis = 1)/1500

#gradient in z
zgrad = np.diff(wN, axis =0) \
        /np.diff(hr04['grid']['depthw'], axis = 0)

#divergence of uN
divUN = -1*(xgrad+zgrad)

#grid set up
hr04['mbc'] = np.repeat(hr04['grid']['nbnc'][np.newaxis,:],
                                hr04['grid']['depth'].shape[0], axis = 0)

#average N
nbar = np.median(hr04['nw']['nute'], axis = 0)
nb = np.empty(nbar.shape[1])
for j in range(nbar.shape[1]) :
    nb[j] = np.interp(-75, hr04['grid']['depth'][:,j], nbar[:,j])

#plotting
fig, (ax1, ax2) = plt.subplots(2,1, sharex = True,
                               figsize = (8,6),
                               gridspec_kw = {'height_ratios':[1,4]})
ax1.plot(hr04['grid']['nbnc'], nb, 'k')
ax1.set_ylim([0,0.52])
ax1.set_ylabel('$Mdn(Tr_{75m})$ \n [$mmol \ N \ m^{-3}$]')
ax1.set_xticks(np.arange(-5,1, 0.5))
ax1.grid()

cf2 = ax2.contourf(hr04['mbc'], hr04['grid']['depth'], divUN,
                   cmap = cmocean.cm.balance,
                   vmin =-0.05, vmax =0.05,
                   levels = np.arange(-0.05,0.06,0.01),
                   extend = 'both')
ax2.plot([-5,0.5], [-75, -75], 'k--', alpha = 0.5)
ax2.patch.set_facecolor('silver')
ax2.grid()
ax2.set_ylim([-2000,0])
ax2.set_ylabel('Depth [m]')
ax2.set_xlabel('Surface Bounce [n]')
ax2.set_xticks(np.arange(-5,1, 0.5))
ax2.set_xlim([-5,0.4])
divider = make_axes_locatable(ax2)
cax = divider.append_axes('bottom', size = '7%', pad = 0.5)
cb2 =fig.colorbar(cf2,cax = cax, orientation='horizontal',
                  label = r'$-\nabla \cdot (\overline{U \ Tr})$')
cb2.ax.yaxis.set_tick_params(color = 'black')
plt.setp(plt.getp(cb2.ax.axes,'yticklabels'), color = 'black')

fig.tight_layout()
