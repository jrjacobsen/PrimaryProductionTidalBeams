"""
Created on Fri Mar 11 07:38:06 2022

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
root = '/home/jjacob2/runs/npzd_coupled/dimensional/nsf_ridge04/output02/'
n_file = 'floats25m/roms_his.nc'

#load file
ncfile = nc4(root+n_file, 'r')

#grid
#depth
romsvars = {'Vstretching' : ncfile.variables['Vstretching'][0], \
            'Vtransform' :ncfile.variables['Vtransform'][0], \
            'theta_s' : ncfile.variables['theta_s'][0], \
            'theta_b' : ncfile.variables['theta_b'][0], \
            'N' : ncfile.variables['Cs_r'].size, \
            'h' : ncfile.variables['h'][:], \
            'hc': ncfile.variables['hc'][0]}

z = dep._set_depth(root+n_file,
                       romsvars, 'rho',
                       romsvars['h'])[:,35,:450]

#surface bounce distance
N = 2e-3
om = (2*np.pi)/(3600*12.4)
theta_g = np.pi/2 - np.arccos(om/N)
bot1 = 2/np.tan(theta_g)
surf = 2*bot1

#x-origin critical point on east side of ridge
x = np.repeat(( ncfile.variables['x_rho'][35,:450]/1000 - 628.2796)[np.newaxis, :],
              z.shape[0], axis = 0)/surf

def nw_load(ncfile, idx) :

    nw ={}
    nw['nute'] = ncfile.variables['dye_01'][idx['time'], :,
                                           idx['lat'],
                                           idx['lon']]
    nw['w'] = ncfile.variables['w'][idx['time'], :,
                                           idx['lat'],
                                           idx['lon']]*3600
    nw['u'] = ncfile.variables['u'][idx['time'], :,
                                           idx['lat'],
                                           idx['ulon']]*3600

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
    grid['depthw'] =dep._set_depth(FileName, romsvars, 'w',
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
    grid['nbnc'] = cr_cent/57.2

    return grid

#slicing
idx = {'time' : slice(149, None),
       'lat' : 35,
       'lon' : slice(None,450),
       'ulon': slice(None,451)}

hr04 = {}
hr04['nw'] = nw_load(ncfile, idx)
hr04['grid'] = grid_load(n_file, ncfile, idx)
hr04['time'] = ncfile.variables['ocean_time'][idx['time']]/3600
hr04['grid']['mbc'] = np.repeat(hr04['grid']['nbnc'][np.newaxis,:],
                                hr04['grid']['depth'].shape[0], axis = 0)
hr04['kv'] = ncfile.variables['AKt'][1,2,35,1] #constant 

 #shift nutrient to u and w points
pad = np.concatenate((hr04['nw']['nute'][:,:,0:1],
                      hr04['nw']['nute'][:,:,:],
                      hr04['nw']['nute'][:,:,-2:-1]), axis = 2)
nute_u = 0.5*(pad[:,:,:-1]+pad[:,:,1:])
pad = np.concatenate((hr04['nw']['nute'][:,0:1,:],
                      hr04['nw']['nute'][:,:,:],
                      hr04['nw']['nute'][:,-2:-1,:]), axis = 1)
nute_w = 0.5*(pad[:,:-1,:]+ pad[:,1:,:])


#mean product
uN = np.mean(nute_u*hr04['nw']['u'], axis = 0)
wN = np.mean(nute_w*hr04['nw']['w'], axis = 0)

#gradient in x
xgrad = np.diff(uN, axis = 1)/1500
xgrad_t = np.diff(nute_u*hr04['nw']['u'], axis = 2)/1500

#gradient in z
zgrad = np.diff(wN, axis =0) \
        /np.diff(hr04['grid']['depthw'], axis = 0)
zgrad_t = np.empty(xgrad_t.shape)
zdiff = np.empty(xgrad_t.shape)
for t in range(xgrad_t.shape[0]) :
    zgrad_t[t,:,:] = np.diff(nute_w[t,:,:]*hr04['nw']['w'][t,:,:],
                             axis = 0) \
                                /np.diff(hr04['grid']['depthw'], axis = 0)
    #second derivative
    fh = nute_w[t,2:,:]
    f_h= nute_w[t,:-2,:]
    fc = nute_w[t,1:-1,:]
    h0 = np.diff(hr04['grid']['depthw'], axis = 0)
    h = 0.5*(h0[1:,:]+h0[:-1,:])

    d2r_dz2 = (fh - 2*fc + f_h)/h**2

    #shift to rho points
    pad = np.concatenate((d2r_dz2[0:1,:], d2r_dz2, d2r_dz2[-2:-1,:]),
                         axis = 0)
    d2rdz2 = 0.5*(pad[1:,:] + pad[:-1])

    #diffusion
    zdiff[t,:,:] = hr04['kv']*d2rdz2


#zero diffusivity through top and bottom
zdiff[:,0,:] = 0
zdiff[:,-1,:] = 0

#divergence of uN
div = xgrad+zgrad
div_tt = xgrad_t + zgrad_t

#shift to "dt" points
div_t = 0.5*(div_tt[:-1,:,:] + div_tt[1:,:,:])
xdiv_t = 0.5*(xgrad_t[:-1,:,:] + xgrad_t[1:,:,:])
zdiv_t = 0.5*(zgrad_t[:-1,:,:] + zgrad_t[1:,:,:])
zdiff_t = 0.5*(zdiff[:-1,:,:] + zdiff[1:,:,:])
time = 0.5*(hr04['time'][:-1] + hr04['time'][1:])

#average density
rhobar = np.mean(hr04['nw']['nute'], axis = 0)

#time rate of changeof density
drho_dt = np.diff(hr04['nw']['nute'], axis = 0)\
                  /(np.mean(np.diff(hr04['time'])))

#set points to extract time series
#hr04['grid']['nbnc'][381]

idist = 397, 381, 367, 346
rho = np.empty((drho_dt.shape[0], len(idist)))
drdt = np.empty((drho_dt.shape[0], len(idist)))
zdif= np.empty((drho_dt.shape[0], len(idist)))
afd = np.empty((drho_dt.shape[0], len(idist)))
xdiv =np.empty((drho_dt.shape[0], len(idist)))
zdiv = np.empty((drho_dt.shape[0], len(idist)))
for i in range(len(idist)) :
    for t in range(drho_dt.shape[0]) :
        rho[t,i] = np.interp(-75,hr04['grid']['depth'][:,idist[i]],
                              hr04['nw']['nute'][t,:,idist[i]])
        drdt[t,i] = np.interp(-75,hr04['grid']['depth'][:,idist[i]],
                              drho_dt[t,:,idist[i]])
        zdif[t,i] = np.interp(-75,hr04['grid']['depth'][:,idist[i]],
                              zdiff_t[t,:,idist[i]])
        afd[t,i] = np.interp(-75,hr04['grid']['depth'][:,idist[i]],
                              div_t[t,:,idist[i]])
        xdiv[t,i] = np.interp(-75,hr04['grid']['depth'][:,idist[i]],
                              xdiv_t[t,:,idist[i]])
        zdiv[t,i] = np.interp(-75,hr04['grid']['depth'][:,idist[i]],
                              zdiv_t[t,:,idist[i]])

##############################################
fig, (ax1, ax2) = plt.subplots(2,1, constrained_layout=True,
                              sharex = True,
                              figsize = (7,7),
                              gridspec_kw = {'height_ratios': [1,1]})

cf1 = ax1.contourf(x,z,rhobar,
             cmap = cmocean.cm.dense,
             vmin = 0, vmax = 3.5,
             levels = np.arange(0,3.75,0.25))
ax1.set_xticks(np.arange(-3.5,1.0,0.5))
ax1.set_xlim([-3.5,0.5])
ax1.set_ylim([-2000, 0])
#ax1.set_yticks(np.arange(-200,50,50))
ax1.set_ylabel('Depth \n [m]')
ax1.grid()
ax1.plot(x[0,idist[0]], -75, 'b+')
ax1.plot(x[0,idist[1]], -75, 'b+')
ax1.plot(x[0,idist[2]], -75, 'b+')
ax1.plot(x[0,idist[3]], -75, 'b+')
#ax1.text(-3.5, 10, '$\overline{ \rho }$')
ax1.patch.set_facecolor('silver')

divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size = '2%', pad = 0.05)
cb1 =fig.colorbar(cf1,cax = cax, orientation='vertical')#,
                  #ticks = [1.0,2.0,3.0])
cb1.ax.yaxis.set_tick_params(color = 'black')
plt.setp(plt.getp(cb1.ax.axes,'yticklabels'), color = 'black')


cf2 =ax2.contourf(x,z,-div*12.4,
                vmin = -1, vmax = 1,
                levels = np.arange(-1,1, 0.1),
                extend = 'both',
                cmap = cmocean.cm.balance)
ax2.set_xticks(np.arange(-3.5,1.0,0.5))
ax2.set_xlim([-3.5, 0.5])
ax2.set_xlabel('Surface Bounce [n]')
ax2.set_ylabel('Depth \n [m]')
ax2.set_ylim([-2000,0])
ax2.patch.set_facecolor('silver')
ax2.grid()
ax2.plot(x[0,idist[0]], -75, 'b+')
ax2.plot(x[0,idist[1]], -75, 'b+')
ax2.plot(x[0,idist[2]], -75, 'b+')
ax2.plot(x[0,idist[3]], -75, 'b+')
ax2.text(-3.5, 15, r'$-\nabla \cdot (\overline{\hat{u} \ NO_3})$')

divider = make_axes_locatable(ax2)
cax2 = divider.append_axes('right', size = '2%', pad = 0.05)
cb2 =fig.colorbar(cf2,cax = cax2, orientation='vertical')#,
                  #ticks = [1.0,2.0,3.0])
cb2.ax.yaxis.set_tick_params(color = 'black')
plt.setp(plt.getp(cb2.ax.axes,'yticklabels'), color = 'black')


#time series
fig, (ax3,ax4, ax5, ax6) = plt.subplots(4, 1, figsize = (6,6))
ax3.plot(time/12.4, drdt[:,2], color = 'black',
         label = r'$\frac{\partial }{\partial t}Tr$')
ax3.plot(time/12.4, -zdif[:,2], color = 'red',
         label = r'$-kv \frac{\partial^2 }{\partial z^2}Tr$')
ax3.plot(time/12.4, -afd[:,2], color = 'blue',
         label = r'$-\nabla \cdot (u \ Tr)$')
ax3.set_ylim([-0.1, 0.1])
ax3.set_yticks(np.arange(-0.1,0.15,0.05))
ax3.legend(bbox_to_anchor= (1,1))
ax3.set_ylabel('$Tr \ m^{-3} hr^{-1}$')
ax3.set_xticklabels([])
ax3.grid()


ax4.plot(time/12.4, -afd[:,2], color = 'blue',
         label = r'$-\nabla \cdot (u \ Tr)$')
ax4.plot(time/12.4, -zdiv[:,2], color = 'darkviolet',
         label = r'$-\frac{\partial}{\partial z}(w \ Tr)$')
ax4.plot(time/12.4, -xdiv[:,2], color = 'darkorange',
         label = r'$-\frac{\partial}{\partial x}(u \ Tr)$')
ax4.set_ylim([-0.1,0.1])
ax4.set_yticks(np.arange(-0.1,0.15,0.05))
ax4.legend(bbox_to_anchor= (1,1))
ax4.set_xticklabels([])
ax4.set_ylabel('$Tr \ m^{-3} hr^{-1}$')
ax4.grid()

ax5.plot(time/12.4, drdt[:,3], color = 'black')
ax5.plot(time/12.4, -zdif[:,3], color = 'red')
ax5.plot(time/12.4, -afd[:,3], color = 'blue')
ax5.set_ylim([-0.1, 0.1])
ax5.set_yticks(np.arange(-0.1,0.15,0.05))
ax5.set_ylabel('$Tr \ m^{-3} hr^{-1}$')
ax5.set_xticklabels([])
ax5.grid()


ax6.plot(time/12.4, -afd[:,3], color = 'blue')
ax6.plot(time/12.4, -zdiv[:,3], color = 'darkviolet')
ax6.plot(time/12.4, -xdiv[:,3], color = 'darkorange')
ax6.set_ylim([-0.1, 0.1])
ax6.set_yticks(np.arange(-0.1,0.15,0.05))
ax6.set_ylabel('$Tr \ m^{-3} hr^{-1}$')
ax6.set_xlabel('Time [nM2]')
ax6.grid()



