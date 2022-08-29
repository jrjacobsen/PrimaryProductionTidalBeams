"""
Created on Wed Jan 12 16:04:07 2022

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/Ideal_Ridge/npzd_coupling/')
import sys
sys.path.append('/home/jjacob2/python/Ideal_Ridge/npzd_coupling/')
import matplotlib.pyplot as plt
from netCDF4 import Dataset as nc4
import numpy as np
import obs_depth_vs145 as dep
import cmocean


#file directories 

#flat case
root = '/home/jjacob2/runs/npzd_coupled/dimensional/'
FileName = root+'olig_flat/output02/roms_his.nc'

#load file
ncfile = nc4(FileName,'r')



#slicing 
idx = {'tslice' : [0,-1],
       'latslice': 3,
       'lonslice': 360,
       'target' : -75
       }

#import variables
n = ncfile.variables['NO3'][idx['tslice'],:,idx['latslice'],idx['lonslice']]
p = ncfile.variables['phytoplankton'][idx['tslice'],:,idx['latslice'],idx['lonslice']]
z = ncfile.variables['zooplankton'][idx['tslice'],:,idx['latslice'],idx['lonslice']]
d = ncfile.variables['detritus'][idx['tslice'],:,idx['latslice'],idx['lonslice']]
time = ncfile.variables['ocean_time'][idx['tslice']]/(12.4*3600)
tr = ncfile.variables['dye_01'][0,:,idx['latslice'],idx['lonslice']]

#depth
romsvars = {'Vstretching' : ncfile.variables['Vstretching'][0], \
            'Vtransform' : ncfile.variables['Vtransform'][0], \
            'theta_s' : ncfile.variables['theta_s'][0], \
            'theta_b' : ncfile.variables['theta_b'][0], \
            'N' : ncfile.variables['Cs_r'].size, \
            'h' : ncfile.variables['h'][:], \
            'hc': ncfile.variables['hc'][0]}

depth = dep._set_depth(FileName, romsvars, 'rho', romsvars['h'])[:,idx['latslice'],idx['lonslice']]

#plot variables
ymax = -210

#plotting
fig,(ax1, ax2, ax3, ax4) = plt.subplots(1, 4,figsize=(7,4))
ax1.grid()
ax1.plot(n[0,:], depth, color = 'darkviolet', linewidth = 2,
         label = 'Initial')
ax1.plot(n[-1,:], depth, color = 'darkviolet', linewidth = 3, linestyle = ':',
         label = 'Final')
ax1.plot(tr, depth, color = 'grey', linestyle ='--', linewidth = 3,
         label = 'Tracer')
ax1.set_ylim([ymax, 0])
ax1.set_yticks(np.arange(-200,50,50))
ax1.set_ylabel('Depth [m]')
ax1.set_xlim([-0.2, 6])
ax1.set_xticks(np.arange(0,7,2))
ax1.set_xlabel('Nutrient')
ax1.xaxis.set_label_position('top')

ax2.grid()
ax2.plot(p[0,:], depth, color = 'green', linewidth = 2)
ax2.plot(p[-1,:], depth, color = 'green', linewidth = 3, linestyle = ':')
ax2.set_ylim([ymax, 0])
ax2.set_yticks(np.arange(-200,50,50))
ax2.set_yticklabels([])
ax2.set_xlim([-0.01, 0.25])
ax2.set_xticks(np.arange(0,0.3,0.1))
ax2.set_xlabel('Phytoplankton')
#ax2.xaxis.set_ticks_position('top')
ax2.xaxis.set_label_position('top')

ax3.grid()
ax3.plot(z[0,:], depth, color = 'black', linewidth = 2)
ax3.plot(z[-1,:], depth, color = 'black', linewidth = 3, linestyle = ':')
ax3.set_ylim([ymax, 0])
ax3.set_yticks(np.arange(-200,50,50))
ax3.set_yticklabels([])
ax3.set_xlim([-0.001,0.05])
ax3.set_xticks(np.arange(0,0.075,0.025))
ax3.set_xlabel('Zooplankton')
ax3.xaxis.set_label_position('top')

ax4.grid()
ax4.plot(d[0,:], depth, color = 'darkorange', linewidth = 2)
ax4.plot(d[-1,:], depth, color = 'darkorange', linewidth = 3, linestyle = ':')
ax4.set_ylim([ymax, 0])
ax4.set_yticklabels([])
ax4.set_yticks(np.arange(-200,50,50))
ax4.set_xlim([-0.005, 0.05])
ax4.set_xticks(np.arange(0,0.075,0.025))
ax4.set_xlabel('Detritus')
#ax4.xaxis.set_ticks_position('top')
ax4.xaxis.set_label_position('top')
