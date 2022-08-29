"""
Created on Thu Apr 14 16:28:11 2022

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/Ideal_Ridge/npzd_coupling/')
import sys
sys.path.append('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')

import numpy as np
import matplotlib.pyplot as plt
import pp_bm_tr_plot_funcs as pfun
from netCDF4 import Dataset as nc4

#file roots
dat_root = '/home/jjacob2/python/Ideal_Ridge/wave_forcing/data/'
nc_root = '/home/jjacob2/runs/npzd_coupled/dimensional/'

disp = nc4(dat_root + 'MedianTidalDisplacement.nc', 'r')

N = 2e-3
om = (2*np.pi)/(3600*12.4)
theta_g = np.pi/2 - np.arccos(om/N)
bot1 = 2/np.tan(theta_g)
surf = 2*bot1

#displacement variables
xslice = slice(66,732)
dist = disp.variables['dist'][:]/1000
nb02 = (disp.variables['dist'][:]/1000 - 26.1346)/surf
nb04 = (disp.variables['dist'][:]/1000 - 28.2796)/surf
nb06 = (disp.variables['dist'][:]/1000 - 28.8889)/surf
range02 = disp.variables['range02'][:]
range04 = disp.variables['range04'][:]
range06 = disp.variables['range06'][:]

#displacement
fig,(ax1, ax2) = plt.subplots(2,1, sharex = True,
                              figsize = (7,5),
                              gridspec_kw = {'height_ratios':[2,3]})
#top
ax1.plot(nb02, range02, color = 'blue', linewidth = 2,
         label = '400 m Ridge')
ax1.plot(nb04, range04, color = 'black', linewidth = 2,
         label = '800 m Ridge')
ax1.plot(nb06, range06, color = 'red', linewidth = 2,
         label = '1200 m Ridge')
ax1.set_ylim([200,300])
ax1.set_xlim([-5,0.5])
ax1.set_xticks(np.arange(-5,1,0.5))
ax1.grid()
ax1.legend(loc ='upper left')
ax1.spines['bottom'].set_visible(False)
ax1.xaxis.tick_top()
ax1.tick_params(labeltop = False)

ax1b = ax1.twiny()
ax1b.plot(dist, range06, color = 'red')
ax1b.set_ylim([100, 350])
ax1b.set_xticks(np.arange(-5*surf, 2*surf, surf))
ax1b.set_xlim([-1*surf*5+28.2796,(0.5*surf+28.2796)])
ax1b.set_xlabel('Distance from Ridge Center [km]')
ax1b.spines['bottom'].set_visible(False)
ax1b.xaxis.tick_top()


#bottom
ax2.plot(nb02, range02, color = 'blue', linewidth = 2)
ax2.plot(nb04, range04, color = 'black', linewidth = 2)
ax2.plot(nb06, range06, color = 'red', linewidth = 2)
ax2.set_ylim([0,85])
ax2.set_xlim([-5, 0.5])
ax2.set_xticks(np.arange(-5,1,0.5))
ax2.grid()
ax2.set_xlabel('Surface Bounce [n]')
ax2.set_ylabel('Displacement [m]')
ax2.spines['top'].set_visible(False)
ax2.xaxis.tick_bottom()

fig.subplots_adjust(hspace=0.1)
