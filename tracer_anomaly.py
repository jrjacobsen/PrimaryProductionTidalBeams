"""
Created on Tue May 17 14:21:32 2022

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/Ideal_Ridge/npzd_coupling/')
import sys
sys.path.append('/home/jjacob2/python/Ideal_Ridge/pub_figs/')

import numpy as np
import matplotlib.pyplot as plt
import pp_bm_tr_plot_funcs as pfun
from netCDF4 import Dataset as nc4
import cmocean
import matplotlib.ticker as ticker

#files
dat_root = '/home/jjacob2/python/Ideal_Ridge/wave_forcing/data/'
root = '/home/jjacob2/runs/npzd_coupled/dimensional/'
FileName = root+'olig_ridge04/output02/roms_his.nc'

#slicing
idx = {'time' : slice(149,410),
       'lat' : 3,
       'lon' : slice(67,734)
         }

#load nc file for grid
ncfile = nc4(FileName, 'r')

#load interpolatd data
flat = nc4(dat_root + 'InterpolatedPP_oligflat.nc', 'r')
hr02 = nc4(dat_root + 'InterpolatedPP_olighr02.nc', 'r')
hr04 = nc4(dat_root + 'InterpolatedPP_olighr04.nc', 'r')
hr06 = nc4(dat_root + 'InterpolatedPP_olighr06.nc', 'r')

#load displacement
disp = nc4(dat_root + 'MedianTidalDisplacement.nc', 'r')

#surface bounce distance
N = 2e-3
om = (2*np.pi)/(3600*12.4)
theta_g = np.pi/2 - np.arccos(om/N)
bot1 = 2/np.tan(theta_g)
sbounce = 2*bot1

#displacement variables
#distance from critical points
x02 = (disp.variables['dist'][idx['lon']]/1000-26.1346)/sbounce
x04 = (disp.variables['dist'][idx['lon']]/1000-28.2796)/sbounce
x06 = (disp.variables['dist'][idx['lon']]/1000-28.8889)/sbounce

#displacement relative to 75 m
disp02 = disp.variables['range02'][idx['lon']] - 75
disp04 = disp.variables['range04'][idx['lon']] - 75
disp06 = disp.variables['range06'][idx['lon']] - 75

#bio variables
#grid
xb02 = (hr02.variables['dist'][:] - 626.1346)/sbounce
zb02 = hr02.variables['depth'][:]
xb04 = (hr04.variables['dist'][:] - 628.2796)/sbounce
zb04 = hr04.variables['depth'][:]
xb06 = (hr06.variables['dist'][:] - 628.8889)/sbounce
zb06 = hr06.variables['depth'][:]

#average difference from flat
bioflat= flat.variables['Ntracer'][:]

bio02 = np.mean(hr02.variables['Ntracer'][:]-bioflat,
                      axis = 0)
bio04 = np.mean(hr04.variables['Ntracer'][:]-bioflat,
                      axis = 0)
bio06= np.mean(hr06.variables['Ntracer'][:]-bioflat,
                      axis = 0)

#flat reference
zflat = flat.variables['depth'][:,294]
bflatbar =  np.median(bioflat[:,:,294], axis = 0)

#example profiles
#2nd surface bounce
sb02 = np.median(hr02.variables['Ntracer'][:,:,294], axis = 0)
sb04 = np.median(hr04.variables['Ntracer'][:,:,296], axis = 0)
sb06 = np.median(hr06.variables['Ntracer'][:,:,296], axis = 0)

#2nd bottom bounce
bb02 = np.median(hr02.variables['Ntracer'][:,:,275], axis = 0)
bb04 = np.median(hr04.variables['Ntracer'][:,:,277], axis = 0)
bb06 = np.median(hr06.variables['Ntracer'][:,:,277], axis = 0)

#plot parameters
ppm = {}
ppm['xmin'] = -5
ppm['xmax'] = 0.5
ppm['xtick'] = np.arange(-5, 1, 0.5)
ppm['ymin'] = -200
ppm['ymax'] = 0
ppm['cmin'] = -1.0
ppm['cmax'] = 1.0
ppm['clev'] = np.arange(-1.0,1.1,0.1)
ppm['cmap'] = cmocean.cm.balance
ppm['lablev'] = [-2.0, -1.5, 0, 1.5, 2.0]

#plotting
fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharex = True,
                                    figsize = (7,7),
                                    gridspec_kw = {'height_ratios':[1,1,1]})
cf = ax1.contourf(xb02, zb02, bio02,
                  cmap = ppm['cmap'],
                  vmin = ppm['cmin'], vmax = ppm['cmax'],
                  levels = ppm['clev'],
                  extend = 'both')
cl = ax1.contour(xb02, zb02, bio02, levels = ppm['lablev'],
                 colors = 'white', alpha = 0.5)
ax1.set_xlim([ppm['xmin'], ppm['xmax']])
ax1.set_xticks(ppm['xtick'])
ax1.set_ylim([ppm['ymin'], ppm['ymax']])
ax1.set_ylabel('Depth [m]')
ax1.grid()

cf = ax2.contourf(xb04, zb04, bio04,
                  cmap = ppm['cmap'],
                  vmin = ppm['cmin'], vmax = ppm['cmax'],
                  levels = ppm['clev'],
                  extend = 'both')
cl = ax2.contour(xb04, zb04, bio04, levels = ppm['lablev'],
                 colors = 'white', alpha = 0.5)
ax2.set_xticks(ppm['xtick'])
ax2.set_ylim([ppm['ymin'], ppm['ymax']])
ax2.set_ylabel('Depth [m]')
ax2.grid()

cf = ax3.contourf(xb06, zb06, bio06,
                  cmap = ppm['cmap'],
                  vmin = ppm['cmin'], vmax = ppm['cmax'],
                  levels = ppm['clev'],
                  extend = 'both')
cl = ax3.contour(xb06, zb06, bio06, levels = ppm['lablev'],
                 colors = 'white', alpha = 0.5)
ax3.set_xticks(ppm['xtick'])
ax3.set_xlabel('Surface Bounce [n]')
ax3.set_ylim([ppm['ymin'], ppm['ymax']])
ax3.set_ylabel('Depth [m]')
ax3.grid()

fig.subplots_adjust(bottom = 0.8)
cbar_ax = fig.add_axes([0.1, -0.025, 0.9, 0.03])
cb =fig.colorbar(cf,cax = cbar_ax, orientation='horizontal',
                  label = 'Tracer Anomaly [$mmol \ N \ m^{-3}$]')
cb.ax.yaxis.set_tick_params(color = 'black')
plt.setp(plt.getp(cb.ax.axes,'yticklabels'), color = 'black')

fig.tight_layout()

#plot reference and example profiles
fig, ax1 = plt.subplots(1,1, figsize = (3,4))
plt.plot(sb02, zb02[:,294], 'r', linewidth = 4)
plt.plot(bb02, zb02[:,275], 'b', linewidth = 4)
plt.plot(bflatbar, zflat, color = 'grey', linewidth = 4.5)
ax1.set_ylim([-200,0])
ax1.set_yticks(np.arange(-200,50,50))
ax1.set_yticklabels([])
ax1.set_xlim([-0.1, 6])
ax1.set_xlabel('400 m Ridge ')
ax1.grid()

fig, ax1 = plt.subplots(1,1, figsize = (3,4))
plt.plot(sb04, zb04[:,296], 'r', linewidth = 4)
plt.plot(bb04, zb04[:,277], 'b', linewidth = 4)
plt.plot(bflatbar, zflat, color = 'grey', linewidth = 4.5)
ax1.set_ylim([-200,0])
ax1.set_yticks(np.arange(-200,50,50))
ax1.set_yticklabels([])
ax1.set_xlim([-0.1, 6])
ax1.set_xlabel('800 m Ridge')
ax1.grid()

fig, ax1 = plt.subplots(1,1, figsize = (3,4))
plt.plot(sb06, zb06[:,296], 'r-', linewidth = 4)
plt.plot(bb06, zb06[:,277], 'b', linewidth = 4)
plt.plot(bflatbar, zflat, color = 'grey', linewidth = 4.5)
ax1.set_ylim([-200,0])
ax1.set_yticks(np.arange(-200,50,50))
ax1.set_yticklabels([])
ax1.set_xlim([-0.1, 6])
ax1.set_xlabel('1200 m Ridge')
ax1.grid()
                          
