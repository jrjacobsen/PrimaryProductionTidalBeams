#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 10:42:10 2022

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')
import sys
sys.path.append('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')

from netCDF4 import Dataset as nc4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import obs_depth_vs145 as dep
import cmocean

#file directories
root = '/home/jjacob2/runs/npzd_coupled/dimensional/'

#load files
flfile = nc4(root+'nsf_ridge04/output02/floats25m/roms_flt.nc', 'r')
ncfile = nc4(root+'nsf_ridge04/output02/floats25m/roms_his.nc', 'r')

#slicing
tstrt = 5961 #tstep of beginning of lagrangian tracer

#time
time = flfile.variables['ocean_time'][tstrt:]

#depth
romsvars = {'Vstretching' : ncfile.variables['Vstretching'][0], \
            'Vtransform' :ncfile.variables['Vtransform'][0], \
            'theta_s' : ncfile.variables['theta_s'][0], \
            'theta_b' : ncfile.variables['theta_b'][0], \
            'N' : ncfile.variables['Cs_r'].size, \
            'h' : ncfile.variables['h'][:], \
            'hc': ncfile.variables['hc'][0]}

z = dep._set_depth(root+'nsf_ridge04/output02/roms_his.nc',
                       romsvars, 'w',
                       romsvars['h'])[:,35,:450]

#surface bounce distance
N = 2e-3
om = (2*np.pi)/(3600*12.4)
theta_g = np.pi/2 - np.arccos(om/N)
bot1 = 2/np.tan(theta_g)
surf = 2*bot1

nb = np.repeat(((ncfile.variables['x_rho'][35,:450]/1000 - \
                 628.2796)/surf)[np.newaxis,:],
               z.shape[0], axis = 0)

#M2 slicing               
iM2 = {}
iM2['4'] = slice(0,1478)
iM2['5'] = slice(1479,2966)
iM2['6'] = slice(2967,4454)
iM2['7'] = slice(4455,5942)
iM2['8'] = slice(5943,7430)
iM2['9'] = slice(7431,8918)
iM2['10']= slice(8919,10406)

#float depth slices
dslice = {}
dslice['25'] = slice(1,801)
dslice['50'] = slice(803,1603)
dslice['75'] = slice(1605,2405)
dslice['100'] = slice(2407,3207)
dslice['125'] = slice(3209,4009)
dslice['150'] = slice(4011, 4811)
dslice['175'] = slice(4813,5613)
dslice['200'] = slice(5615,6415)

#float positons
zpos25 = flfile.variables['depth'][tstrt:, dslice['25']]
xpos25 = (flfile.variables['x'][tstrt:, dslice['25']]/1000 - 628.2796)/surf

zpos75 = flfile.variables['depth'][tstrt:, dslice['75']]
xpos75 = (flfile.variables['x'][tstrt:, dslice['75']]/1000 - 628.2796)/surf

zpos125 = flfile.variables['depth'][tstrt:, dslice['125']]
xpos125 = (flfile.variables['x'][tstrt:, dslice['125']]/1000 - 628.2796)/surf

zpos200 = flfile.variables['depth'][tstrt:, dslice['200']]
xpos200 = (flfile.variables['x'][tstrt:, dslice['200']]/1000 - 628.2796)/surf

def stack(pos, iM2):
    a = np.dstack((pos[iM2['5'], :],
                   pos[iM2['6'], :],
                   pos[iM2['7'], :],
                   pos[iM2['8'], :],
                   pos[iM2['9'], :],
                   pos[iM2['10'], :]))
    return a

#median over tidal cycle
x25 = np.median(stack(xpos25, iM2), axis = 2)
z25 = np.median(stack(zpos25, iM2), axis = 2)
x75 = np.median(stack(xpos75, iM2), axis = 2)
z75 = np.median(stack(zpos75, iM2), axis = 2)
x125 = np.median(stack(xpos125, iM2), axis = 2)
z125 = np.median(stack(zpos125, iM2), axis = 2)
x200 = np.median(stack(xpos200, iM2), axis = 2)
z200 = np.median(stack(zpos200, iM2), axis = 2)

strt = 201
stop = 450
step = 4

ocol = 'black'
icol = 'black'

fig, ax1 = plt.subplots(1,1, figsize = (6,4))
ax1.plot([0.5, 0], [0, -2000], color = 'grey', linewidth = 3)
ax1.plot([0, -0.5], [-2000, 0], color = 'grey', linewidth = 3)
ax1.plot([-0.5, -1.0], [0, -2000], color = 'grey', linewidth = 3)
ax1.plot([-1.0, -1.5], [-2000, 0], color = 'grey', linewidth = 3)
ax1.plot([-1.5, -2.0], [0, -2000], color = 'grey', linewidth = 3)
ax1.plot([-2.0, -2.5], [-2000, 0], color = 'grey', linewidth = 3)
ax1.plot([-2.5, -3.0], [0, -2000], color = 'grey', linewidth = 3)
ax1.plot([-3.0, -3.5], [-2000, 0], color = 'grey', linewidth = 3)
ax1.plot([-3.5, -4.0], [0, -2000], color = 'grey', linewidth = 3)
ax1.plot([-4.0, -4.5], [-2000, 0], color = 'grey', linewidth = 3)
ax1.plot([-4.5, -5.0], [0, -2000], color = 'grey', linewidth = 3)

ax1.plot(x25[:, strt:stop:step],
         z25[:, strt:stop:step],
         color = icol, linewidth=1.35)


ax1.plot(x75[:, strt:stop:step],
         z75[:, strt:stop:step],
         color = icol, linewidth=1.35)


ax1.plot(x125[:, strt:stop:step],
         z125[:, strt:stop:step],
         color = icol, linewidth=1.35)


ax1.plot(x200[:, strt:stop:step],
         z200[:, strt:stop:step],
         color = icol, linewidth=1.35)
ax1.set_xlim([-5,0.5])
ax1.set_xticks(np.arange(-5, 1, 0.5))
ax1.set_xlabel('Surface Bounce [n]')
ax1.set_ylim([-300,0])
ax1.set_ylabel('Depth [m]')
ax1.grid()
