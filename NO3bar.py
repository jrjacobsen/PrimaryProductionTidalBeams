"""
Created on Wed Jul  6 11:32:55 2022

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
root = '/home/jjacob2/runs/npzd_coupled/dimensional/'

#netcdf directory and file for saving
n_file = 'nsf_ridge04/output02/roms_his.nc'
flatname = root+'nsf_flat/output02/roms_his.nc'

#load files
ncfile = nc4(root+n_file, 'r')
ncflat = nc4(flatname,'r')


def nw_load(ncfile, idx) :

    nw ={}
    nw['nute'] = np.median(ncfile.variables['NO3'][idx['time'], :,
                                           idx['lat'],
                                           idx['lon']], axis = 0)

    return nw

def grid_load(FileName,ncfile, idx):
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

    grid['nM2'] = (ncfile.variables['ocean_time'][idx['time']])/(12.4*3600)

    #distance from ridge [km]
    dist =(ncfile.variables['x_rho'][idx['lat'],
                                              idx['lon']]/1000-600)
    #center on critical point on eastern side of ridge
    cr_cent = dist - 28.2796

    #nondimensionalize by surface bounce distance
    #surface bounce distance
    N = 2e-3
    om = (2*np.pi)/(3600*12.4)
    theta_g = np.pi/2 - np.arccos(om/N)
    bot1 = 2/np.tan(theta_g)
    surf = 2*bot1

    grid['nbnc'] = cr_cent/surf

    return grid

#slicing
idx = {'time' : slice(149, None),
       'lat' : 35,
       'lon' : slice(None,450),
       'ulon': slice(None,451)}

#ridge case
ridge = {}
ridge['nw'] = nw_load(ncfile, idx)
ridge['grid'] = grid_load(n_file, ncfile, idx)
ridge['grid']['mbc'] = np.repeat(ridge['grid']['nbnc'][np.newaxis,:],
                                 ridge['grid']['depth'].shape[0], axis = 0)

#flat case
idx = {'time' : slice(149, None),
       'lat' : 35,
       'lon' : 400,
       'ulon': slice(None,451)}
flat = {}
flat['nw'] = nw_load(ncflat, idx)
flat['grid'] = grid_load(flatname, ncflat, idx)



#plot median Nitrate
fig, (ax, ax1) = plt.subplots(1,2, sharey = True,
                       figsize = (8,5),
                       gridspec_kw = {'width_ratios' : [4,1]})
cf =ax.contourf(ridge['grid']['mbc'], ridge['grid']['depth'],
                ridge['nw']['nute'],
                cmap = cmocean.cm.dense,
                vmin = 0, vmax = 3.25,
                levels = np.arange(0,3.5,0.25))
ax.set_xlim([-5, 0.5])
ax.set_xticks(np.arange(-5,1,0.5))
ax.set_ylim([-195,-5])
ax.grid(alpha = 0.5)
ax.set_xlabel('Surface Bounce [n]')
ax.set_ylabel('Depth [m]')

ax1.plot(flat['nw']['nute'], flat['grid']['depth'], color = 'grey', linewidth = 2)
ax1.set_ylim([-195,-5])
ax1.set_xlim([-0.1,3.25])
ax1.set_xticks(np.arange(0,4,1))
ax1.set_xlabel('NO3')
ax1.grid()

fig.subplots_adjust(bottom = 0.8)
cbar_ax = fig.add_axes([0.1, -0.025, 0.67, 0.03])
cb =fig.colorbar(cf,cax = cbar_ax, orientation='horizontal',
                  label = 'NO3 [$mmol \ N \ m^{-3}$]')
cb.ax.yaxis.set_tick_params(color = 'black')
plt.setp(plt.getp(cb.ax.axes,'yticklabels'), color = 'black')

fig.tight_layout()

  
