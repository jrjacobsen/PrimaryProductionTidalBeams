#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 10:35:41 2022

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')
import sys
sys.path.append('/home/jjacob2/python/Ideal_Ridge/wave_forcing/')

from netCDF4 import Dataset as nc4
import numpy as np


def cycle_min(time, fltz, t_idx) :
    #cycle time
    itime = time[t_idx]

    #min over cycle
    imin = np.argmin(fltz[t_idx,:], axis = 0)

    #find index of cycle min each location
    ctime = np.empty(fltz.shape[1])
    ctime.fill(np.nan)
    for j in range(fltz.shape[1]) :
        ctime[j] = np.argwhere(time == itime[imin[j]])

    return ctime


def get_cycle(flfile, time, tstrt, dep_idx, m2_idx, nextM2) :
    """
    Parameters
    ----------
    flfile : netCDF4 file
        ROMS float file.
    time : np array float
        Time in secods.
    tstrt : integer
        Index of time of float release.
    dep_idx : slice
        Location in float file (dim 1) for desiered depth.
    m2_idx : slice
        Loaction in float file (dim 0) of 4th M2 cycle 
        relative to float release time (0.
    nextM2 : integer
        Number of next M2 cycle 
        (must be strictly greater than 4).

    Returns
    -------
    i4th : 2D array of int
        Start and end time indices of 4th M2 cycle realtive to min depth
    inth : 2d array of int
        Start and end time indices of nth M2 cycle specified by nextM2
        """

    #Number of steps in M2  cycle
    nstep = (12.4*3600)/np.mean(np.diff(time))

    #number of steps in M2 cycle
    step = (nextM2-4)*nstep

    #float depth
    fltz = flfile.variables['depth'][tstrt:,dep_idx]

    #compute min of 4th M2 cycle
    m2_4th = cycle_min(time,fltz, m2_idx)

    #array of start and end indices 4th M2
    i4th = np.array((m2_4th, m2_4th+nstep-1), dtype = int).T

    #next M2 cycle relative to minimum of 4th M2
    m2_nth = m2_4th+step

    #array of start and end indices 4th M2
    inth = np.array((m2_nth, m2_nth+nstep-1), dtype = int).T


    return i4th, inth


def get_pos(flfile, time, tstrt, dep_idx, iM2) :
    """
    Parameters
    ----------
    flfile : netCDF4 file
        ROMS float file.
    time : np array float
        Time in secods.
    tstrt : integer
        Index of time of float release.
    dep_idx : slice
        Indices in float file (dim 1) for desiered depth.
    iM2 : 2D array of int
        Start and end time indices of M2 cycle as fu
        nction of distance from ridge.

    Returns
    -------
    x : 2D array of float
        x position of floats (time, distance from ridge).
    z : 2D array of float
        z position of floats (time, distance from ridge).

    """
#Number of steps in M2  cycle
    nstep = (12.4*3600)/np.mean(np.diff(time))

    #load position for initial float depth 
    zpos = flfile.variables['depth'][tstrt:, dep_idx]
    xpos = flfile.variables['x'][tstrt:, dep_idx]

    #extract float position over M2 cycle
    z = np.empty((int(nstep-1), iM2.shape[0]))
    z.fill(np.nan)
    x = np.empty((int(nstep-1), iM2.shape[0]))
    x.fill(np.nan)
    for j in range(iM2.shape[0]) :
        #time indices
        st = iM2[j, 0]
        en = iM2[j, 1]

        #position of float over M2  cycle
        z[:,j] = zpos[st:en, j]
        x[:,j] = xpos[st:en, j]

    return x, z
