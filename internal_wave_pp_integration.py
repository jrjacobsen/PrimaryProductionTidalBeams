"""
Created on Tue Apr 26 16:16:58 2022

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')
import sys
sys.path.append('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as nc4
import cmocean
import float_tidecycle as m2

#file directories
root = '/home/jjacob2/runs/npzd_coupled/dimensional/'
dat_root = '/home/jjacob2/python/Ideal_Ridge/wave_forcing/data/'

#load files
hsfile = nc4(root+'nsf_ridge04/output02/roms_his.nc')
flfile = nc4(root+'nsf_ridge04/output02/floats25m/roms_flt.nc', 'r')
ncfile = nc4(dat_root+'NPZD_LagrangianFloats.nc','r')
dsfile = nc4(dat_root + 'MedianTidalDisplacement.nc', 'r')

#calculated variables
dist = (ncfile.variables['dist'][:]/1000-628.2796)/57.2
depth= ncfile.variables['depth'][:]
nM2 = ncfile.variables['time'][:]/(12.4*3600)
nute = ncfile.variables['nute'][:]
phyt = ncfile.variables['phyt'][:]


#displacement
#slicing
tstrt = 5961 #tstep of beginning of lagrangian tracer

#time
time = flfile.variables['ocean_time'][tstrt:]

#4th barotropic  M2 cycle
m2_idx = slice(0,1478)

#float depth slices
dslice =  slice(1605,2405)

dispM2 = np.empty((7,800))
nuteM2 = np.empty((7,800))
for nth in range(4,11) :
    i4th, i_th = m2.get_cycle(flfile, time, tstrt,
                                        dslice, m2_idx, nth)

    #depth of floats
    x_th, z_th = m2.get_pos(flfile, time, tstrt,
                                          dslice, i_th)

    #range
    dispM2[nth-4, :] = np.max(z_th, axis = 0) - np.min(z_th, axis = 0)

    #median nutrient over M2 cycle
    nuteM2[nth-4, :] = np.median(nute[i_th[:,1], 2,:], axis = 0)


#params
par = {}
par['vm'] = hsfile.variables['Vm_NO3'][0]
par['kn'] = hsfile.variables['K_NO3'][0]
par['kex'] = hsfile.variables['K_ext'][0]

#
pp = phyt[0,2,:]*(par['vm']*nute[0,2,:])/(par['kn']+ nute[0,2,:])*np.exp(-75*par['kex'])

#select values within bounce
iwidth =5
ibounce = np.arange(247,400,38)
nute_sb = np.empty((nuteM2.shape[0],iwidth*2+1,ibounce.shape[0]))
amp_sb = np.empty((nuteM2.shape[0],iwidth*2+1,ibounce.shape[0]))
for i in range(ibounce.shape[0]) :
    iselect = np.arange(ibounce[i]-iwidth,ibounce[i]+iwidth+1,1)
    nute_sb[:,:,i] = nuteM2[:,iselect]
    amp_sb[:,:,i] = dispM2[:,iselect]/2


#median primary production
def pp_bar(p0, n0, z0, a0, par, dn, da, t) :
    #mechalis menton 
    n_up = (par['vm']*(n0+dn))/(par['kn']+(n0 + dn))

    #depth
    z = (a0+da)*np.sin((t/12.4))+z0

    #light 
    light = p0*np.exp(par['kex']*z)

    #primary production
    pprod = np.median(n_up*light)*(106/16)

    return pprod

#set initial conditions
p0 = 1
n0 = 0.01
z0 = -75
a0 = 0


#set purturbations
dn = np.linspace(0,0.5,100)
da = np.linspace(0,50,100)
t = np.linspace(0,12.4,100)

pmap = np.empty((da.shape[0], dn.shape[0]))
for i in range(da.shape[0]):
    for j in range(dn.shape[0]) :
        pmap[i,j] = pp_bar(p0, n0, z0, a0, par, dn[j], da[i], t)

#derivative of PP
dPP_dN = np.diff(pmap, axis = 1)/np.mean(np.diff(dn))
dPP_dA = np.diff(pmap, axis = 0)/np.mean(np.diff(da))

#regrid
pad = np.concatenate((dPP_dA[0:1,:],dPP_dA, dPP_dA[-2:-1,:]), axis = 0)
dPdA = 0.5*(pad[:-1,:]+pad[1:,:])

pad = np.concatenate((dPP_dN[:,0:1],dPP_dN, dPP_dN[:,-2:-1]), axis = 1)
dPdN = 0.5*(pad[:,:-1]+pad[:,1:])

gradP = np.gradient(pmap)

amp = np.repeat(da[np.newaxis ,:], 100, axis = 0)
n0 = np.repeat(dn[:, np.newaxis], 100, axis = 1)

#refernce PP
pp_ref = 0#pp_bar(p0, n0, z0, a0, par, dn[0], da[0], t)

ilat =slice(246,400)

fig, ax = plt.subplots(1,1,
                       figsize = (8,6))
cf =ax.contourf(amp, n0, pmap, cmap = cmocean.cm.algae,
             vmin = 0, vmax = 2,
             levels = np.arange(0,2.1, 0.1))
ax.scatter(dispM2[:,ilat].flatten()/2, nuteM2[:,ilat].flatten(),
           facecolors='none', edgecolors = 'lightgrey')
ax.scatter(amp_sb[:,:,-1].flatten(), nute_sb[:,:,-1].flatten(),
           facecolors='red', edgecolors = 'lightgrey')
ax.scatter(amp_sb[:,:,-2].flatten(), nute_sb[:,:,-2].flatten(),
           facecolors='blue', edgecolors = 'lightgrey')
ax.scatter(amp_sb[:,:,:-3].flatten(), nute_sb[:,:,:-3].flatten(),
           facecolors='white', edgecolors = 'lightgrey')

ax.quiver(amp[1::15, ::10], n0[1::15,::10],
          gradP[1][1::15, ::10], gradP[0][1::15, ::10],
          facecolor = 'black', edgecolor = 'white',
          angles = 'uv', pivot = 'mid',
          minshaft = 1.5,
          width = 0.01, headwidth = 4,
          headlength = 5, headaxislength = 5)
ax.set_xlabel('Amplitude  [m]')
ax.set_ylabel('$NO_3$ [$mmol \ N \ m^{-3} $]')
ax.set_ylim([0,0.5])
fig.colorbar(cf, label = r'$ \overline{PP}$ [$d^{-1} $]')
