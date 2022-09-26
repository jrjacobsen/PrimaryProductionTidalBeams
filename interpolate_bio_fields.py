"""
Created on Mon Feb 21 17:56:23 2022

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')
import sys
sys.path.append('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')

import numpy as np
from netCDF4 import Dataset as nc4
from scipy.interpolate import griddata
import obs_depth_vs145 as dep
import npzd_interpolation as npzd

#select files to load
root = '/home/jjacob2/runs/npzd_coupled/dimensional/'
FileName = root+'olig_LTflat/output02/roms_his.nc'

#save as file 
n_root = '/home/jjacob2/python/Ideal_Ridge/wave_forcing/data/'

#new file name
name_f = 'Interpolated_PP_LTflat.nc'

#slicing for grid
t_start = 149
idx = {'time' : slice(t_start,2235),
       'lat' : 3,
       'lon' : slice(1667,2343)
         }

#slicing for interp
i_idx = {'z' : slice(175,None),
         'lat' : 3,
         'lon' : slice(1667,2343)
         }

#load grid
grid = npzd.load_grid(FileName, idx)

#regular grid for interpolation
idepth = np.arange(-250, 0, 5)
ilon = np.arange(np.min(grid['dist']),
                  np.max(grid['dist']), 1.5)

#create grid 
ilon, idepth = np.meshgrid(ilon, idepth)
ilevels = {'ilon' : ilon,
          'idepth': idepth}

#loop of over space, time average,
n = np.empty((grid['time'].shape[0], idepth.shape[0], ilon.shape[1]))
p = np.empty(n.shape)
tr = np.empty(n.shape)
for t in range(n.shape[0]) :
    #select time step
    i_idx['time'] = t+t_start

    #load file
    ncdat = npzd.load_NPZD(FileName, i_idx)

    #interpolate data
    n[t,:,:] = griddata((grid['dist'][i_idx['z'],:].flatten(),
                         grid['depth'][i_idx['z'],:].flatten()),
                        ncdat['n'].flatten(),
                        (ilon, idepth))
    p[t,:,:] = griddata((grid['dist'][i_idx['z'],:].flatten(),
                         grid['depth'][i_idx['z'],:].flatten()),
                        ncdat['p'].flatten(),
                        (ilon, idepth))
    tr[t,:,:] = griddata((grid['dist'][i_idx['z'],:].flatten(),
                         grid['depth'][i_idx['z'],:].flatten()),
                        ncdat['tr'].flatten(),
                        (ilon, idepth))


#Primary production
vm = ncdat['Vm']
kext = ncdat['kext']
ks = ncdat['ks']
ipp = (vm*n*p*np.exp(kext*idepth))/(ks + n)

#save interpolated values
#create files
f = nc4(n_root+name_f, 'w', format = 'NETCDF4')

#create dimensions
f.createDimension('dist', ilon.shape[1])
f.createDimension('z', idepth.shape[0])
f.createDimension('time',grid['time'].shape[0])

#create variables
dist_f = f.createVariable('dist', 'f4', ('z','dist'))
depth_f = f.createVariable('depth', 'f4', ('z','dist'))
tr_f = f.createVariable('Ntracer', 'f4', ('time','z','dist'))
pp_f = f.createVariable('PP', 'f4', ('time','z','dist'))
bm_f = f.createVariable('BM', 'f4', ('time','z','dist'))

#store data
dist_f[:] = ilon
depth_f[:] = idepth
tr_f[:] = tr
pp_f[:] = ipp
bm_f[:] = p


f.close()
