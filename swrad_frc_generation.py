import os
os.chdir('/home/jjacob2/python/Ideal_Ridge/npzd_coupling/')
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as nc4
import obs_depth_vs145 as dep

#root
root = '/home/jjacob2/runs/npzd_coupled/dimensional/spring_neap/diurnal_light/'

ncfile = nc4(root+'ncfiles/jrj_diurnal_slim_frc.nc','r+')
ncgrid = nc4(root+'ridge_grd_LThr04_2Dslim.nc')

#load dimension sizes
time = ncfile.variables['frc_time'][:]
xrho = ncgrid.variables['x_rho'][:]
yrho = ncgrid.variables['y_rho'][:]

#create dimensions
srf_time_dim = ncfile.createDimension('srf_time', None)
eta_rho_dim = ncfile.createDimension('eta_rho', xrho.shape[0])
xi_rho_dim = ncfile.createDimension('xi_rho', xrho.shape[1])

#create variables
srf_time = ncfile.createVariable('srf_time', 'f8', 'srf_time')
x_rho = ncfile.createVariable('x_rho', 'f8', ('eta_rho', 'xi_rho'))
y_rho = ncfile.createVariable('y_rho', 'f8', ('eta_rho', 'xi_rho'))
swrad = ncfile.createVariable('swrad', 'f8', ('srf_time', 'eta_rho', 'xi_rho'))

#put data in variables
srf_time[:] = time
x_rho[:] = xrho
y_rho[:] = yrho

#shortwave rad set to constant value used in chapter 1
swrad[:] = np.ones((time.shape[0], xrho.shape[0], yrho.shape[1]))*206.3

ncfile.close()
