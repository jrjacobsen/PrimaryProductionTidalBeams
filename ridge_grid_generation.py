"""
Created on Fri Jun 25 16:34:05 2021

@author: jjacob2
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as nc4

root = '/home/jjacob2/runs/npzd_coupled/dimensional/olig_LThr04/'


#new grid file
newfile = nc4(root+'ridge_grd_LThr04_2D.nc','r+')



#lat for coriolis
lat = 21.7

#grid
x_max = 1200*1000*5 #m
dx = 1.5*1000
y_max = 100*60
dy = 1.5*1000

#bathymetry
#ridge width, height, and bottom depth
a = 30*1000 #m
h_max = 800 #m
z_max = 2000

#new x coord's
new_x_rho = np.arange(-dx/2, x_max+dx, dx)
new_x_psi = np.arange(0, x_max+dx, dx)

#new y coord's
new_y_rho = np.arange(-dy/2, y_max+dy, dy)
new_y_psi = np.arange(0, y_max+dy, dy)

#repeat over other axis
x_rho0 = np.repeat(new_x_rho[np.newaxis, :], new_y_rho.shape[0], axis = 0)
x_v0 = np.repeat(new_x_rho[np.newaxis, :], new_y_psi.shape[0], axis = 0)
x_psi0 = np.repeat(new_x_psi[np.newaxis, :], new_y_psi.shape[0], axis = 0)
x_u0 = np.repeat(new_x_psi[np.newaxis, :], new_y_rho.shape[0], axis = 0)

y_rho0 = np.repeat(new_y_rho[:, np.newaxis], new_x_rho.shape[0], axis = 1)
y_u0 =  np.repeat(new_y_rho[:, np.newaxis], new_x_psi.shape[0], axis = 1)
y_psi0 = np.repeat(new_y_psi[:, np.newaxis], new_x_psi.shape[0], axis = 1)
y_v0 = np.repeat(new_y_psi[:, np.newaxis], new_x_rho.shape[0], axis = 1)

#calculate bathymetry
#center x dimension on mid point
x_cent = new_x_rho - new_x_rho[int(new_x_rho.shape[0]/2)]

#generate ridge
h_new = np.zeros(new_x_rho.shape)
for istep in range(new_x_rho.shape[0]) :
    if (np.abs(x_cent[istep]) < a):
        h_new[istep] = h_max*(1-(x_cent[istep]**2)/(a**2))**2

#repeat over y dim
h0 = np.repeat(np.abs(h_new - z_max)[np.newaxis, :], new_y_rho.shape[0], axis = 0)

#new grid info
pm0 = (1/dx)*np.ones(x_rho0.shape)
pn0 = (1/dy)*np.ones(x_rho0.shape)
el0 = y_max
xl0 = x_max
f0 = 0*np.ones(x_rho0.shape)*2*7.2921159*10**(-5)*np.sin(lat*np.pi/180)


#create dimensions
newfile.createDimension('eta_rho',x_rho0.shape[0])
newfile.createDimension('xi_rho', x_rho0.shape[1])

newfile.createDimension('eta_psi',x_psi0.shape[0])
newfile.createDimension('xi_psi', x_psi0.shape[1])

newfile.createDimension('eta_v', x_v0.shape[0])
newfile.createDimension('xi_v', x_v0.shape[1])

newfile.createDimension('eta_u', x_u0.shape[0])
newfile.createDimension('xi_u', x_u0.shape[1])

#create variables
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

#put data into variable
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

newfile.close()
