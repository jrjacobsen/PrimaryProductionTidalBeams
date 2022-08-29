"""
Created on Sun Jul 20 12:18:48 2021


Linear temperature change; constant salinity
N set to 2e-03

S-levels are set-able here

Adds grid and  both physical and NPZD variables 

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/Ideal_Ridge/npzd_coupling/')
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as nc4
import obs_depth_vs145 as dep

#root
root = '/home/jjacob2/runs/npzd_coupled/dimensional/olig_LTflat/'

#load old ini file
oldfile = nc4('/home/jjacob2/runs/npzd_coupled/dimensional/nsf_ridge04/initial/roms_ridge_hr04_NTcomb_npzd_ini.nc','r')

#load new ini file
newfile = nc4(root+'initial/roms_flat_LT_2D_olig_npzd_ini.nc','r+')

#load grid file
grid = nc4(root+'flat_grd_LT_2D.nc','r')

#load  stable solution from NPZD
F_npzd = nc4('/home/jjacob2/runs/bio_toy/stability_bio/output/olig02_7200hr/roms_his.nc')

#compute depth npzd
romsvars = {'Vstretching' : F_npzd.variables['Vstretching'][0], \
            'Vtransform' : F_npzd.variables['Vtransform'][0], \
            'theta_s' : F_npzd.variables['theta_s'][0], \
            'theta_b' : F_npzd.variables['theta_b'][0], \
            'N' : F_npzd.variables['temp'].shape[1], \
            'h' : F_npzd.variables['h'][:], \
            'hc': F_npzd.variables['hc'][0]}

npzd  = {'depth' : dep._set_depth(F_npzd, romsvars, 'rho', romsvars['h'])}


##NPZD profile with full depth
npzd['NO3'] = F_npzd.variables['NO3'][0,:,3,3]
npzd['phyt']= F_npzd.variables['phytoplankton'][0,:,3,3]
npzd['zoop']= F_npzd.variables['zooplankton'][0,:,3,3]
npzd['detr']= F_npzd.variables['detritus'][0,:,3,3]

plt.figure()
plt.plot(npzd['NO3'], npzd['depth'][:,3,3], 'b')
plt.plot(npzd['phyt'], npzd['depth'][:,3,3], 'g')
plt.plot(npzd['zoop'], npzd['depth'][:,3,3], 'r')
plt.plot(npzd['detr'], npzd['depth'][:,3,3], 'k')
plt.ylim([-200,0])
plt.xlim([0,5])

#set numer of vertical levels
n_lev = 200
dims = (oldfile.variables['gls'][:].shape[0], n_lev+1)

#repeat over other axis
x_rho0 = grid.variables['x_rho'][:]
x_v0 = grid.variables['x_v'][:]
x_psi0 = grid.variables['x_psi'][:]
x_u0 = grid.variables['x_u'][:]

y_rho0 = grid.variables['y_rho'][:]
y_u0 =  grid.variables['y_u'][:]
y_psi0 = grid.variables['y_psi'][:]
y_v0 = grid.variables['y_v'][:]

#bathymetry
h0 = grid.variables['h'][:]

#new grid info
pm0 = grid.variables['pm'][:]
pn0 = grid.variables['pn'][:]
el0 = grid.variables['el'][:]
xl0 = grid.variables['xl'][:]
f0 = grid.variables['f'][:]*0


#define dimensions of new variables
ones_srhorho = np.ones((dims[0],dims[1]-1,x_rho0.shape[0],x_rho0.shape[1]))
ones_swrho = np.ones((dims[0], dims[1], x_rho0.shape[0],x_rho0.shape[1]))
ones_srhou = np.ones((dims[0], dims[1]-1, x_u0.shape[0],x_u0.shape[1]))
ones_srhov = np.ones((dims[0], dims[1]-1, x_v0.shape[0],x_v0.shape[1]))

#set spatially varing fields to zero at appropriate dimension
#4D fields
omega0 = ones_swrho*0
gls0 = ones_swrho*0
tke0 = ones_swrho*0
u0 = ones_srhou*0
v0 = ones_srhov*0
w0  = ones_swrho*0

#3D fields
bustr0 = ones_srhou[:,0,:,:]*0
bvstr0 = ones_srhov[:,0,:,:]*0
ubar0 = ones_srhou[:,0,:,:]*0
vbar0 = ones_srhov[:,0,:,:]*0
zeta0 = ones_srhorho[:,0,:,:]*0
#set AKs and AKv to 5*10^-6 and 5*10^-5
AKs0 = ones_swrho*(1*10**(-5))
AKv0 = ones_swrho*(1*10**(-5))

#set top and bottom diffusivity to 0
AKs0[:,0,:,:] = 0
AKs0[:,-1,:,:] = 0
AKv0[:,0,:,:] = 0
AKv0[:,-1,:,:] = 0

#set grid to equal spacing
Vstretching0 = 1
Vtransform0 = 1
theta_s0 = 0
theta_b0 = 0
#generate linear change in temp 
#compute depth
romsvars = {'Vstretching' : Vstretching0, \
            'Vtransform' : Vtransform0, \
            'theta_s' : theta_s0, \
            'theta_b' : theta_b0, \
            'N' : n_lev, \
            'h' : grid.variables['h'][:], \
            'hc': newfile.variables['hc'][0]}

depth = dep._set_depth(newfile, romsvars, 'rho', romsvars['h'])

#linear temp 
lintemp = np.linspace(8.48,12.95,100)
lindepth = np.linspace(-2000, 0, 100)

#constant salt
salt0 = ones_srhorho*34

#propagate over domain
temp0 = np.empty(ones_srhorho.shape)
temp0.fill(np.nan)
for d0 in range(temp0.shape[0]) :
    for d2 in range(temp0.shape[2]) :
        for d3 in range(temp0.shape[3]):
            temp0[d0, :, d2, d3] = np.interp(depth[:,d2,d3], lindepth, lintemp)

phys = {'depth' : depth,
        'temp' : temp0}

#interpolate N, P, Z, D
#if out of range use largst value
Nutein = np.empty(phys['temp'].shape)
Nutein.fill(np.nan)
Phytin = np.empty(phys['temp'].shape)
Phytin.fill(np.nan)
Zoopin = np.empty(phys['temp'].shape)
Zoopin.fill(np.nan)
Detrin = np.empty(phys['temp'].shape)
Detrin.fill(np.nan)
for tstep in range(phys['temp'].shape[0]) :
    for ltstep in range(phys['temp'].shape[2]) :
        for lostep in range(phys['temp'].shape[3]) :
            for dstep in range(phys['temp'].shape[1]) :
                #if out of range use largst value
                if (phys['depth'][dstep,ltstep,lostep] < np.min(npzd['depth'])):
                    Nutein[tstep,dstep,ltstep,lostep] = npzd['NO3'][0]
                    Phytin[tstep,dstep,ltstep,lostep] = npzd['phyt'][0]
                    Zoopin[tstep,dstep,ltstep,lostep] = npzd['zoop'][0]
                    Detrin[tstep,dstep,ltstep,lostep] = npzd['detr'][0]
                #interpolate to physical model grid    
                else :
                    Nutein[tstep,dstep,ltstep,lostep] = np.interp(
                                 phys['depth'][dstep,ltstep,lostep],
                                 npzd['depth'][:,3,3],
                                 npzd['NO3'][:])
                    Phytin[tstep,dstep,ltstep,lostep] = np.interp(
                                 phys['depth'][dstep,ltstep,lostep],
                                 npzd['depth'][:,3,3],
                                 npzd['phyt'][:])
                    Zoopin[tstep,dstep,ltstep,lostep] = np.interp(
                                 phys['depth'][dstep,ltstep,lostep],
                                 npzd['depth'][:,3,3],
                                 npzd['zoop'][:])
                    Detrin[tstep,dstep,ltstep,lostep] = np.interp(
                                 phys['depth'][dstep,ltstep,lostep],
                                 npzd['depth'][:,3,3],
                                 npzd['detr'][:])

#create dimensions in nc file
eta_rho = newfile.createDimension('eta_rho',x_rho0.shape[0])
xi_rho = newfile.createDimension('xi_rho', x_rho0.shape[1])

eta_psi = newfile.createDimension('eta_psi',x_psi0.shape[0])
xi_psi = newfile.createDimension('xi_psi', x_psi0.shape[1])

eta_v = newfile.createDimension('eta_v', x_v0.shape[0])
xi_v = newfile.createDimension('xi_v', x_v0.shape[1])

eta_u = newfile.createDimension('eta_u',x_u0.shape[0])
xi_u = newfile.createDimension('xi_u', x_u0.shape[1])

s_rho = newfile.createDimension('s_rho',dims[1]-1)
s_w = newfile.createDimension('s_w',dims[1])

#create variables in nc files 
Vstretching = newfile.createVariable('Vstretching','i1')
Vtransform = newfile.createVariable('Vtransform','i1')
theta_s = newfile.createVariable('theta_s', 'i1')
theta_b = newfile.createVariable('theta_b', 'i1')

#dimensional variables
x_rho = newfile.createVariable('x_rho','f8',("eta_rho","xi_rho"))
x_rho.units = "m"
y_rho = newfile.createVariable('y_rho','f8',("eta_rho","xi_rho"))
y_rho.units = "m"
x_v = newfile.createVariable('x_v','f8',("eta_v","xi_v"))
x_v.units = "m"
y_v = newfile.createVariable('y_v','f8',("eta_v","xi_v"))
y_v.units = "m"
x_psi = newfile.createVariable('x_psi','f8',("eta_psi","xi_psi"))
x_psi.units = "m"
y_psi = newfile.createVariable('y_psi','f8',("eta_psi","xi_psi"))
y_psi.units = "m"
x_u = newfile.createVariable('x_u','f8',("eta_u","xi_u"))
x_u.units = "m"
y_u = newfile.createVariable('y_u','f8',("eta_u","xi_u"))
y_u.units = "m"
#grid parameters
pm = newfile.createVariable('pm','f8',("eta_rho","xi_rho"))
pm.units = "1/m"
pn = newfile.createVariable('pn','f8',("eta_rho","xi_rho"))
pn.units = "1/m"
h = newfile.createVariable('h','f8',("eta_rho","xi_rho"))
h.units = "m"
f = newfile.createVariable('f','f8',("eta_rho","xi_rho"))
f.units = "1/s"
el = newfile.createVariable('el','f8')
el.units = "m"
xl = newfile.createVariable('xl','f8')
xl.units = "m"
#4D fields on w points
AKs = newfile.createVariable('AKs','f8',("ocean_time", "s_w", "eta_rho","xi_rho"))
AKv = newfile.createVariable('AKv','f8',("ocean_time", "s_w", "eta_rho","xi_rho"))
omega = newfile.createVariable('omega','f8',("ocean_time", "s_w", "eta_rho","xi_rho"))
gls = newfile.createVariable('gls','f8',("ocean_time", "s_w", "eta_rho","xi_rho"))
tke = newfile.createVariable('tke','f8',("ocean_time", "s_w", "eta_rho","xi_rho"))
w = newfile.createVariable('w','f8',("ocean_time", "s_w", "eta_rho","xi_rho"))
#4D fields on u and v points
u = newfile.createVariable('u','f8',("ocean_time", "s_rho", "eta_u","xi_u"))
v = newfile.createVariable('v','f8',("ocean_time", "s_rho", "eta_v","xi_v"))
#3D field
bustr = newfile.createVariable('bustr','f8',("ocean_time","eta_u","xi_u"))
bvstr = newfile.createVariable('bvstr','f8',("ocean_time","eta_v","xi_v"))
ubar = newfile.createVariable('ubar','f8',("ocean_time","eta_u","xi_u"))
vbar = newfile.createVariable('vbar','f8',("ocean_time","eta_v","xi_v"))
zeta = newfile.createVariable('zeta','f8',("ocean_time","eta_rho","xi_rho"))
#4D field on rho points
temp = newfile.createVariable('temp','f8',("ocean_time", "s_rho", "eta_rho","xi_rho"))
salt = newfile.createVariable('salt','f8',("ocean_time", "s_rho", "eta_rho","xi_rho"))
#NPZD variables
NO3 = newfile.createVariable('NO3','f8',
                            ("ocean_time", "s_rho", "eta_rho","xi_rho"))
phytoplankton = newfile.createVariable('phytoplankton','f8',
                                      ("ocean_time", "s_rho", "eta_rho","xi_rho"))
zooplankton = newfile.createVariable('zooplankton','f8',
                                      ("ocean_time", "s_rho", "eta_rho","xi_rho"))
detritus = newfile.createVariable('detritus','f8',
                                      ("ocean_time", "s_rho", "eta_rho","xi_rho"))


#put variable in nc file
Vstretching[:] = Vstretching0
Vtransform[:] = Vtransform0
theta_s[:] = theta_s0
theta_b[:] = theta_b0
x_rho[:] = x_rho0
y_rho[:] = y_rho0
x_v[:] = x_v0
y_v[:] = y_v0
x_u[:] = x_u0
y_u[:] = y_u0
x_psi[:] = x_psi0
y_psi[:] = y_psi0
pm[:] = pm0
pn[:] = pn0
h[:] = h0
f[:] = f0
el[:] = el0
xl[:] = xl0
AKs[:] = AKs0
AKv[:] = AKv0
omega[:] = omega0
gls[:] = gls0
tke[:] = tke0
w[:] = w0
u[:] = u0
v[:] = v0
temp[:] = temp0
salt[:] = salt0
bustr[:] = bustr0
bvstr[:] = bvstr0
ubar[:] = ubar0
vbar[:] = vbar0
zeta[:] = zeta0

NO3[:] = Nutein
phytoplankton[:] = Phytin
zooplankton[:] = Zoopin
detritus[:] = Detrin

#check result
n = newfile.variables['NO3'][0,:,3,200]
p = newfile.variables['phytoplankton'][0,:,3,200]
z = newfile.variables['zooplankton'][0,:,3,200]
d = newfile.variables['detritus'][0,:,3,200]

newfile.close()


