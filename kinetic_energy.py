"""
Created on Fri Apr  8 12:33:15 2022

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')
import sys
sys.path.append('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')
sys.path.append('/home/jjacob2/python/NPZD/')
from netCDF4 import Dataset as nc4
import numpy as np
import  matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean
import obs_depth_vs145 as dep

#file directories
root = '/home/jjacob2/runs/npzd_coupled/dimensional/'
n_root = '/home/jjacob2/python/Ideal_Ridge/wave_forcing/data/'

#file name
hr02name = root+'nsf_ridge02/output02/roms_his.nc'
hr04name = root+'nsf_ridge04/output02/floats25m/roms_his.nc'
hr06name = root+'nsf_ridge06/output02/roms_his.nc'

#load files
hr02nc = nc4(hr02name, 'r')
hr04nc = nc4(hr04name, 'r')
hr06nc = nc4(hr06name, 'r')

def nw_load(ncfile, idx) :

    nw ={}
    alpha = ncfile.variables['Tcoef'][0]
    beta = ncfile.variables['Scoef'][0]
    salt = ncfile.variables['salt'][idx['time'], :,
                                           idx['lat'],
                                           idx['lon']]
    temp = ncfile.variables['temp'][idx['time'], :,
                                           idx['lat'],
                                           idx['lon']]
    r0 = ncfile.variables['R0'][0]

    #compute density T0 = 10 and S0 = 30 from .in file
    nw['rho'] = alpha*(temp - 10) + beta*(salt - 30)+r0

    #velocity
    w = ncfile.variables['w'][idx['time'], :,
                                           idx['lat'],
                                           idx['lon']]
    u = ncfile.variables['u'][idx['time'], :,
                                           idx['lat'],
                                           idx['ulon']]
    nw['u'] = (0.5*(u[:,:,1:] + u[:,:,:-1]))

    nw['w'] =0.5*(w[:,1:,:] + w[:,:-1,:])

    return nw

def grid_load(FileName,ncfile, idx, cr_dist):
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
    cr_cent = dist - cr_dist

    #nondimensionalize by surface bounce distance
    grid['nbnc'] = cr_cent/57.2

    return grid

#time
nM2 = hr04nc.variables['ocean_time'][:]/(12.4*3600)

#slicing
idx = {'time' : slice(149,None),
       'lat' : 35,
       'lon' : slice(228,438),
       'ulon': slice(228,439)}

hr02 = {}
hr02['nw'] = nw_load(hr02nc, idx)
hr02['grid'] = grid_load(hr02name, hr02nc, idx, 26.1346)

hr04 = {}
hr04['nw'] = nw_load(hr04nc, idx)
hr04['grid'] = grid_load(hr04name, hr04nc, idx, 28.2796)

hr06 = {}
hr06['nw'] = nw_load(hr06nc, idx)
hr06['grid'] = grid_load(hr06name, hr06nc, idx, 28.8889)

#grid set up
hr02['mbc'] = np.repeat(hr02['grid']['nbnc'][np.newaxis,:],
                                hr02['grid']['depth'].shape[0], axis = 0)
hr04['mbc'] = np.repeat(hr04['grid']['nbnc'][np.newaxis,:],
                                hr04['grid']['depth'].shape[0], axis = 0)
hr06['mbc'] = np.repeat(hr06['grid']['nbnc'][np.newaxis,:],
                                hr06['grid']['depth'].shape[0], axis = 0)
#KE
hr02['KE'] = np.mean(0.5*hr02['nw']['rho']*(hr02['nw']['u']**2+hr02['nw']['w']**2),
                     axis = 0)
hr04['KE'] = np.mean(0.5*hr04['nw']['rho']*(hr04['nw']['u']**2+hr04['nw']['w']**2),
                     axis = 0)
hr06['KE'] = np.mean(0.5*hr06['nw']['rho']*(hr06['nw']['u']**2+hr06['nw']['w']**2),
                     axis = 0)

#plotting
vmin = 0
vmax = 100
step = 5
fig,(ax1, ax2, ax3) = plt.subplots(3,1, sharex=True,
                                   figsize = (6,6),
                                   constrained_layout=True)
cf =ax1.contourf(hr02['mbc'], hr02['grid']['depth'], hr02['KE'],
                cmap = cmocean.cm.amp,
                vmin = vmin, vmax =vmax,
                levels = np.arange(vmin, vmax+step, step),
                extend = 'max')
ax1.contour(hr02['mbc'], hr02['grid']['depth'], hr02['nw']['rho'][149,:,:],
            colors = 'black', alpha = 0.5)
ax1.patch.set_facecolor('silver')
ax1.set_ylim([-2000,0])
ax1.set_ylabel('Depth [m]')
ax1.set_xlim([-5, 0.5])
ax1.set_xticks(np.arange(-5,1,0.5))
ax1.grid()


cf =ax2.contourf(hr04['mbc'], hr04['grid']['depth'], hr04['KE'],
                cmap = cmocean.cm.amp,
                vmin = vmin, vmax =vmax,
                levels = np.arange(vmin, vmax+step, step),
                extend = 'max')
ax2.contour(hr04['mbc'], hr04['grid']['depth'], hr04['nw']['rho'][149,:,:],
            colors = 'black', alpha = 0.5)
ax2.patch.set_facecolor('silver')
ax2.set_ylim([-2000,0])
ax2.set_ylabel('Depth [m]')
ax2.set_xlim([-5, 0.5])
ax2.set_xticks(np.arange(-5,1,0.5))
ax2.grid()

cf =ax3.contourf(hr06['mbc'], hr06['grid']['depth'], hr06['KE'],
                cmap = cmocean.cm.amp,
                vmin = vmin, vmax =vmax,
                levels = np.arange(vmin, vmax+step, step),
                extend = 'max')
ax3.contour(hr06['mbc'], hr06['grid']['depth'], hr06['nw']['rho'][149,:,:],
            colors = 'black', alpha = 0.5)
ax3.patch.set_facecolor('silver')
ax3.set_ylim([-2000,0])
ax3.set_ylabel('Depth [m]')
ax3.set_xlim([-5, 0.5])
ax3.set_xticks(np.arange(-5,1,0.5))
ax3.grid()
ax3.set_xlabel('Surface Bounce [n]')
fig.colorbar(cf, location = 'bottom', label = r'$\frac{1}{2} \ \overline{\rho (u^2 + w^2)}$ [ J ]')


plt.figure()
plt.contourf(hr04['mbc'], hr04['grid']['depth'], np.mean(hr04['nw']['w'], axis = 0))
plt.scatter(np.linspace(-1,-1.5, 100), np.linspace(-2000,0,100))


fig, ax1 = plt.subplots(1,1)
