#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 13:41:23 2022

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')
import sys
sys.path.append('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')
sys.path.append('/home/jjacob2/python/NPZD/')
from netCDF4 import Dataset as nc4
import numpy as np
from scipy.interpolate import griddata
import  matplotlib.pyplot as plt
import cmocean
import NPZD_Implicit_Lagrangian as npzd

#file directories
root = '/home/jjacob2/runs/npzd_coupled/dimensional/'
n_root = '/home/jjacob2/python/Ideal_Ridge/wave_forcing/data/'

#load files
flfile = nc4(root+'nsf_ridge04/output02/roms_flt.nc', 'r')
hsfile = nc4(root+'nsf_ridge04/output02/roms_his.nc', 'r')
ncfile = nc4(n_root+'NPZD_LagrangianFloats.nc','r')

#surface bounce distance
N = 2e-3
om = (2*np.pi)/(3600*12.4)
theta_g = np.pi/2 - np.arccos(om/N)
bot1 = 2/np.tan(theta_g)
surf = 2*bot1

#distance same units as floats
dist = hsfile.variables['x_rho'][35,:]

#NPZD parameters
param = npzd.params(hsfile)

#depth
ifltstart = 5961
#float depth slices
dslice = {}
dslice['5'] = slice(1,800)
dslice['10'] = slice(803,1602)
dslice['15'] = slice(1605,2404)
dslice['20'] = slice(2407,3206)
dslice['25'] = slice(3209,4008)
dslice['50'] = slice(4011, 4810)
dslice['75'] = slice(4813,5612)
dslice['100'] = slice(5615,6414)
dslice['125'] = slice(6417,7216)
dslice['150'] = slice(7219, 8018)
dslice['175'] = slice(8021,8820)
dslice['200'] = slice(8823,9622)

fltz = np.dstack([flfile.variables['depth'][ifltstart:,
                                            dslice['5']],
                  flfile.variables['depth'][ifltstart:,
                                            dslice['10']],
                  flfile.variables['depth'][ifltstart:,
                                            dslice['15']],
                  flfile.variables['depth'][ifltstart:,
                                            dslice['20']],
                  flfile.variables['depth'][ifltstart:,
                                            dslice['25']],
                  flfile.variables['depth'][ifltstart:,
                                            dslice['50']],
                  flfile.variables['depth'][ifltstart:,
                                            dslice['75']],
                  flfile.variables['depth'][ifltstart:,
                                            dslice['100']],
                  flfile.variables['depth'][ifltstart:,
                                            dslice['125']],
                  flfile.variables['depth'][ifltstart:,
                                            dslice['175']],
                  flfile.variables['depth'][ifltstart:,
                                            dslice['200']]])
fltz = np.moveaxis(fltz, 1, -1)
#xposition
fltx = np.dstack([flfile.variables['x'][ifltstart:,
                                            dslice['5']],
                  flfile.variables['x'][ifltstart:,
                                            dslice['10']],
                  flfile.variables['x'][ifltstart:,
                                            dslice['15']],
                  flfile.variables['x'][ifltstart:,
                                            dslice['20']],
                  flfile.variables['x'][ifltstart:,
                                            dslice['25']],
                  flfile.variables['x'][ifltstart:,
                                            dslice['50']],
                  flfile.variables['x'][ifltstart:,
                                            dslice['75']],
                  flfile.variables['x'][ifltstart:,
                                            dslice['100']],
                  flfile.variables['x'][ifltstart:,
                                            dslice['125']],
                  flfile.variables['x'][ifltstart:,
                                            dslice['175']],
                  flfile.variables['x'][ifltstart:,
                                            dslice['200']]])
fltx = np.moveaxis(fltx, 1, -1)

#light
light = np.exp(param['Kext']*fltz)

#mean position
xbar = np.mean(fltx, axis = 0)
zbar = np.mean(fltz, axis = 0)

#regular grid
dx = np.mean(np.diff(dist))
grid_x , grid_y = np.meshgrid(np.linspace(np.min(dist), np.max(dist), 800),
                              np.arange(-195,0, 5)[::-1])
glight = np.empty((fltx.shape[0], grid_x.shape[0], grid_x.shape[1]))
for t in range(fltx.shape[0]):
    glight[t,:,:]  =griddata((fltx[t,:,:].flatten(), fltz[t,:,:].flatten()),
                             light[t,:,:].flatten(),
                             (grid_x, grid_y),
                             method = 'linear')

m2name = np.array(np.arange(4,11,1),dtype = str)
#M2 slicing
iM2 = {}
iM2['4'] = slice(0,1478)
iM2['5'] = slice(1479,2966)
iM2['6'] = slice(2967,4454)
iM2['7'] = slice(4455,5942)
iM2['8'] = slice(5943,7430)
iM2['9'] = slice(7431,8918)
iM2['10']= slice(8919,10406)
delta_IR = np.empty((len(iM2),glight.shape[1], glight.shape[2]))
for i in range(glight.shape[1]) :
    for j in range(m2name.shape[0]) :
        delta_IR[j,i,:] = np.max(glight[iM2[m2name[j]],i,:], axis = 0) - \
                          np.min(glight[iM2[m2name[j]],i,:], axis = 0)

fig, (ax1, ax2)   = plt.subplots(2,1)
cf = ax1.contourf((grid_x/1000-628.2796)/surf, grid_y, np.mean(glight, axis= 0),
                   cmap = cmocean.cm.solar,
                   vmin = 0, vmax = 0.6,
                   levels = np.arange(0,0.7,0.05),
                   extend = 'max')
ax1.set_xlim([-5, 0.5])
ax1.set_xticks(np.arange(-5,1,0.5))
ax1.set_xticklabels([])
ax1.set_ylim([-200,-25])
ax1.set_ylabel('Depth [m]')
ax1.grid()
ax1.text(-4.25, -185, '$\overline{f_{IR}}$', color = 'white')

cf =ax2.contourf((grid_x/1000-628.2796)/surf,grid_y,np.mean(delta_IR, axis = 0),
                 cmap = cmocean.cm.solar,
                   vmin = 0, vmax = 0.6,
                   levels = np.arange(0,0.7,0.05),
                   extend = 'max')
ax2.set_xlim([-5, 0.5])
ax2.set_xticks(np.arange(-5,1,0.5))
ax2.set_ylim([-200,-25])
ax2.grid()
ax2.set_xlabel('Surface Bounce [n]')
ax2.set_ylabel('Depth [m]')
ax2.text(-4.45, -185, '$\overline{\Delta \ f_{IR}}$', color = 'white')
fig.subplots_adjust(right = 0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(cf, cax=cbar_ax,
             label = '$\overline{f_{\Delta IR}}$  [%]')
