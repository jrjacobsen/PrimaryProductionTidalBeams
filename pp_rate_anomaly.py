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
from netCDF4 import Dataset as nc4
import cmocean


#files
dat_root = '/home/jjacob2/python/Ideal_Ridge/wave_forcing/data/'

#slicing
idx = {'time' : slice(149,410),
       'lat' : 35,
       'lon' : slice(67,734)
         }


#load interpolatd data
flat = nc4(dat_root + 'InterpolatedPP_oligflat.nc', 'r')
hr02 = nc4(dat_root + 'InterpolatedPP_olighr02.nc', 'r')
hr04 = nc4(dat_root + 'InterpolatedPP_olighr04.nc', 'r')
hr06 = nc4(dat_root + 'InterpolatedPP_olighr06.nc', 'r')


#surface bounce distance
N = 2e-3
om = (2*np.pi)/(3600*12.4)
theta_g = np.pi/2 - np.arccos(om/N)
bot1 = 2/np.tan(theta_g)
sbounce = 2*bot1


#bio variables
#grid
xb02 = (hr02.variables['dist'][:] - 626.1346)/sbounce
zb02 = hr02.variables['depth'][:]
xb04 = (hr04.variables['dist'][:] - 628.2796)/sbounce
zb04 = hr04.variables['depth'][:]
xb06 = (hr06.variables['dist'][:] - 628.8889)/sbounce
zb06 = hr06.variables['depth'][:]

#average difference from flat
bioflat= flat.variables['PP'][:]/flat.variables['BM'][:]

bio02 = np.mean(hr02.variables['PP'][:]/hr02.variables['BM'][:]-bioflat,
                      axis = 0)
bio04 = np.mean(hr04.variables['PP'][:]/hr04.variables['BM'][:]-bioflat,
                      axis = 0)
bio06= np.mean(hr06.variables['PP'][:]/hr06.variables['BM'][:]-bioflat,
                      axis = 0)

#flat reference
zflat = flat.variables['depth'][:,294]
bflatbar =  np.median(bioflat[:,:,294], axis = 0)

#example profiles
#2nd surface bounce
sb02 = np.median(hr02.variables['PP'][:,:,294]/hr02.variables['BM'][:,:,294],
                 axis = 0)
sb04 = np.median(hr04.variables['PP'][:,:,296]/hr04.variables['BM'][:,:,296],
                 axis = 0)
sb06 = np.median(hr06.variables['PP'][:,:,296]/hr06.variables['BM'][:,:,296],
                 axis = 0)

#2nd bottom bounce
bb02 = np.median(hr02.variables['PP'][:,:,275]/hr02.variables['BM'][:,:,275],
                 axis = 0)
bb04 = np.median(hr04.variables['PP'][:,:,277]/hr04.variables['BM'][:,:,277],
                 axis = 0)
bb06 = np.median(hr06.variables['PP'][:,:,277]/hr06.variables['BM'][:,:,277],
                 axis = 0)

#new color map
newcmap = cmocean.tools.crop_by_percent(cmocean.cm.delta,
                                        25, which = 'min', N = None)

#plot parameters
ppm = {}
ppm['xmin'] = -5
ppm['xmax'] = 0.5
ppm['xtick'] = np.arange(-5, 1, 0.5)
ppm['ymin'] = -200
ppm['ymax'] = 0
ppm['cmin'] = -0.02
ppm['cmax'] = 0.04
ppm['clev'] = np.arange(-0.02,0.045,0.005)
ppm['cmap'] = newcmap
ppm['lablev']= [0, 0.05, 0.1, 0.2]

#plotting
fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharex = True,
                                    figsize = (7,7),
                                    gridspec_kw = {'height_ratios':[1,1,1]})
cf = ax1.contourf(xb02, zb02, bio02,
                  cmap = ppm['cmap'],
                  vmin = ppm['cmin'], vmax = ppm['cmax'],
                  levels = ppm['clev'],
                  extend = 'max')
cl = ax1.contour(xb02, zb02, bio02,
                 levels = ppm['lablev'], colors = 'white', alpha = 0.5)
#fmt.create_dummy_axis()
#ax1.clabel(cl, cl.levels, fmt = fmt)
#ax1.clabel(cl, cl.levels, inline = True, fmt = fmt, fontsize = 8)
ax1.set_xlim([ppm['xmin'], ppm['xmax']])
ax1.set_xticks(ppm['xtick'])
ax1.set_ylim([ppm['ymin'], ppm['ymax']])
ax1.set_ylabel('Depth [m]')
ax1.grid()

cf = ax2.contourf(xb04, zb04, bio04,
                  cmap = ppm['cmap'],
                  vmin = ppm['cmin'], vmax = ppm['cmax'],
                  levels = ppm['clev'],
                  extend = 'max')
cl = ax2.contour(xb04, zb04, bio04,
                 levels = ppm['lablev'], colors = 'white', alpha = 0.5)
#fmt.create_dummy_axis()
#ax2.clabel(cl, cl.levels, fmt = fmt)
#ax2.clabel(cl, cl.levels, inline = True, fmt = fmt, fontsize = 8)
#ax2.set_xlim([ppm['xmin'], ppm['xmax']])
ax2.set_xticks(ppm['xtick'])
ax2.set_ylim([ppm['ymin'], ppm['ymax']])
ax2.set_ylabel('Depth [m]')
ax2.grid()

cf = ax3.contourf(xb06, zb06, bio06,
                  cmap = ppm['cmap'],
                  vmin = ppm['cmin'], vmax = ppm['cmax'],
                  levels = ppm['clev'],
                  extend = 'max')
cl = ax3.contour(xb06, zb06, bio06,
                 levels = ppm['lablev'], colors = 'white', alpha = 0.5)
#fmt.create_dummy_axis()
#ax3.clabel(cl, cl.levels, fmt = fmt)
#ax3.clabel(cl, cl.levels, inline = True, fmt = fmt, fontsize = 8)
#ax3.set_xlim([ppm['xmin'], ppm['xmax']])
ax3.set_xticks(ppm['xtick'])
ax3.set_xlabel('Surface Bounce [n]')
ax3.set_ylim([ppm['ymin'], ppm['ymax']])
ax3.set_ylabel('Depth [m]')
ax3.grid()

fig.subplots_adjust(bottom = 0.8)
cbar_ax = fig.add_axes([0.1, -0.025, 0.9, 0.03])
cb =fig.colorbar(cf,cax = cbar_ax, orientation='horizontal',
                  label = 'Primary Production per Biomass [$d^{-1}$]')
cb.ax.yaxis.set_tick_params(color = 'black')
plt.setp(plt.getp(cb.ax.axes,'yticklabels'), color = 'black')

fig.tight_layout()
