"""
Created on Tue Jun  2 13:03:41 2015

@author: pmattern

Borrowed by jjacobsen on June 26, 2020
    Added psi, u, & v depths to _set_depth. July 9, 2020
"""

import netCDF4 as nc4
import numpy as np


def _obs_depth_bilinear(depth, Xgrid, Ygrid, z_w, modify_inplace=False):
    '''
    Modelled after matlab function obs_depth.m (see below)
    but with smarter way of handling horizontal interpolation.
    '''
    ind = np.flatnonzero(depth < 0.)
    fp = np.arange(z_w.shape[0])
    for i in ind:
        x = int(np.floor(Xgrid[i]))
        y = int(np.floor(Ygrid[i]))
        a_x = np.floor(Xgrid[i])-Xgrid[i]
        a_y = np.floor(Ygrid[i])-Ygrid[i]
        z  = (z_w[:, y, x]*(1.-a_x)+z_w[:, y, x+1]*a_x)*(1.-a_y)
        z += (z_w[:, y+1, x]*(1.-a_x)+z_w[:, y+1, x+1]*a_x)*a_y
        depth[i] = np.interp(xp=z, fp=fp, x=depth[i])

    return depth


def _obs_depth(depth, Xgrid, Ygrid, z_w, modify_inplace=False):
    '''
    Modelled after matlab function obs_depth.m
    Using only the interpolation part, note that it uses floor to interpolate.
    '''

    '''
    matlab code:
%  Interpolate negative depth values to model fractional z-grid location.

ind = find(S.depth < 0);

if (~isempty(ind));
  for n=1:length(ind),
    iobs = ind(n);
    I = 1.0 + floor(S.Xgrid(iobs));
    J = 1.0 + floor(S.Ygrid(iobs));
    z = reshape(Zw(I,J,:),1,Nw);
    S.Zgrid(iobs) = interp1(z, [0:1:Nr], S.depth(iobs));
  end,
end,
    '''

    ind = np.flatnonzero(depth < 0.)
    fp = np.arange(z_w.shape[0])
    for i in ind:
        x = int(np.floor(Xgrid[i]))
        y = int(np.floor(Ygrid[i]))
        z = z_w[:, y, x]
        depth[i] = np.interp(xp=z, fp=fp, x=depth[i])

    return depth

def _set_depth(romsfile=None, romsvars=None, point_type=None, h=None, zeta=0):
    '''
    Modelled after matlab function "roms/repository/matlab/utility/set_depth.m"
    '''
    if romsvars is None or h is None:
        romsvars = dict()
        with nc4.Dataset(romsfile) as nc:
            for v in ['Vstretching', 'theta_s', 'theta_b', 'hc', 'h']:
                romsvars[v] = nc.variables[v][:]
            romsvars['N'] = nc.variables['Cs_r'].size
            h = nc.variables['h'][:]

    if point_type is None:
        point_type = 'rho'

    use_w = point_type in [5, 'w']

    s, C = _stretching(romsvars=romsvars, use_w=use_w)


    if romsvars['Vtransform'] == 1:
        #Indices to "shift" averaging grid
        N = romsvars['N']
        Lp, Mp = romsvars['h'].shape
        L = Lp - 1
        M = Mp - 1



        if point_type in [1, 'rho', 'density']:
            hr = h
            zetar = zeta

            z = np.empty((N+use_w, hr.shape[0], hr.shape[1]))
            z.fill(np.nan)

            for k in range(z.shape[0]):
                z0 = (s[k]-C[k])*romsvars['hc'] + C[k]*hr
                z[k, :, :] = z0 + zetar*(1.0 + z0/hr)

        elif point_type in [2, 'psi', 'streamfunction']:
            #average bathymetry and free surface
            hp = 0.25*(h[0:L, 0:M] + h[1:Lp, 0:M] + \
                       h[0:L, 1:Mp] + h[1:Lp, 1:Mp])
            zetap = 0.25*(zeta[:, 0:L, 0:M] + zeta[:, 1:Lp, 0:M] + \
                          zeta[:, 0:L, 1:Mp] + zeta[:, 1:Lp, 1:Mp])

            z = np.empty((N+use_w, hp.shape[0], hp.shape[1]))
            z.fill(np.nan)

            #compute depth
            for k in range(z.shape[0]):
                z0 = (s[k]-C[k])*romsvars['hc'] + C[k]*hp
                z[k, :, :] = z0 + zetap*(1.0 + z0/hp)

        elif point_type in [3, 'u']:
            #averaging
            hu = 0.5*(h[0:L, 0:Mp] + h[1:Lp, 0:Mp])
            zetau = 0.5*(zeta[:, 0:L, 0:Mp] + zeta[:, 1:Lp, 0:Mp])

            z = np.empty((N+use_w, hu.shape[0], hu.shape[1]))
            z.fill(np.nan)
            
                        #compute depth
            for k in range(z.shape[0]):
                z0 = (s[k]-C[k])*romsvars['hc'] + C[k]*hu
                z[k, :, :] = z0 + zetau*(1.0 + z0/hu)

        elif point_type in [4, 'v']:
            #averaging
            hv = 0.5*(h[0:Lp, 0:M] + h[0:Lp, 1:Mp])
            zetav = 0.5*(zeta[:, 0:Lp, 0:M] + zeta[:, 0:Lp, 1:Mp])

            z = np.empty((N+use_w, hv.shape[0], hv.shape[1]))
            z.fill(np.nan)

            #compute depth
            for k in range(z.shape[0]):
                z0 = (s[k]-C[k])*romsvars['hc'] + C[k]*hv
                z[k, :, :] = z0 + zetav*(1.0 + z0/hv)

        elif use_w:
            hr = h
            zetar = zeta

            z = np.empty((N+use_w, h.shape[0], h.shape[1]))
            z.fill(np.nan)

            z[0, :, :] = -hr
            for k in range(1, z.shape[0]):
                z0 = (s[k]-C[k])*romsvars['hc'] + C[k]*hr
                z[k, :, :] = z0 + zetar*(1.0 + z0/hr)
        else:
            raise RuntimeError('Invalid point_type "{}".'.format(point_type))

    elif (romsvars['Vtransform'] == 2) :
        if point_type in [1, 'rho', 'density']:
            hr = h
            zetar = zeta

            z = np.empty((romsvars['N']+use_w, hr.shape[0], hr.shape[1]))
            z.fill(np.nan)

            for k in range(z.shape[0]):
                z0 = (romsvars['hc']*s[k] + C[k]*hr)/(romsvars['hc']+hr)
                z[k,:,:] = zetar+(zeta+hr)*z0

    elif point_type in [2, 'psi', 'streamfunction']:
            #average bathymetry and free surface
            hp = 0.25*(h[0:L, 0:M] + h[1:Lp, 0:M] + \
                       h[0:L, 1:Mp] + h[1:Lp, 1:Mp])
            zetap = 0.25*(zeta[:, 0:L, 0:M] + zeta[:, 1:Lp, 0:M] + \
                          zeta[:, 0:L, 1:Mp] + zeta[:, 1:Lp, 1:Mp])

            z = np.empty((romsvars['N']+use_w, hp.shape[0], hp.shape[1]))
            z.fill(np.nan)

            #compute depth
            for k in range(z.shape[0]):
                z0 = (romsvars['hc']*s[k]+C[k]*hp)/(romsvars['hc']+hp)
                z[k,:,:]= zetap+(zetap+hp)*z0
         

       elif point_type in [3, 'u']:
            #averaging
            hu = 0.5*(h[0:L, 0:Mp] + h[1:Lp, 0:Mp])
            zetau = 0.5*(zeta[:, 0:L, 0:Mp] + zeta[:, 1:Lp, 0:Mp])

            z = np.empty((romsvars['N']+use_w, hu.shape[0], hu.shape[1]))
            z.fill(np.nan)

            #compute depth
            for k in range(z.shape[0]):
                z0 = (romsvars['hc']*s[k] + C[k]*hu)/(romsvars['hc']+hu)
                z[k,:,:] = zetau+(zetau+hu)*z0

        elif point_type in [4, 'v']:
            #averaging
            hv = 0.5*(h[0:Lp, 0:M] + h[0:Lp, 1:Mp])
            zetav = 0.5*(zeta[:, 0:Lp, 0:M] + zeta[:, 0:Lp, 1:Mp])

            z = np.empty((romsvars['N']+use_w, hv.shape[0], hv.shape[1]))
            z.fill(np.nan)

            #compute depth
            for k in range(z.shape[0]):
                z0 = (romsvars['hc']*s[k] + C[k]*hv)/(romsvars['hc']+hv)
                z[k,:,:] = zetav+(zetav+hv)*z0

        elif use_w:
            hr = h
            zetar = zeta

            z = np.empty((romsvars['N']+use_w, h.shape[0], h.shape[1]))
            z.fill(np.nan)

            z[0, :, :] = -hr
            for k in range(1, z.shape[0]):
                z0 = (romsvars['hc']*s[k]+C[k]*hr)/(romsvars['hc']+hr)
                z[k,:,:]=zetar+(zetar+hr)*z0

        else:
            raise RuntimeError('Invalid point_type "{}".'.format(point_type))



    return z

  def _stretching(romsvars=None, romsfile=None, use_w=False):
    '''
    Modelled after matlab function "roms/repository/matlab/utility/stretching.m"
    '''
    if romsvars is None:
        romsvars = dict()
        with nc4.Dataset(romsfile) as nc:
            for v in ['Vstretching', 'theta_s', 'theta_b', 'hc']:
                romsvars[v] = nc.variables[v][:]
            romsvars['N'] = nc.variables['Cs_r'].size

    if romsvars['Vstretching'] == 1:
        '''
        matlab code:

        Np=N+1;

        if (Vstretching == 1),

          ds=1.0/N;
          if (kgrid == 1),
            Nlev=Np;

            lev=0:N;
            s=(lev-N).*ds;
          else
            Nlev=N;
            lev=(1:N)-0.5;
            s=(lev-N).*ds;
          end
          if (theta_s > 0),
            Ptheta=sinh(theta_s.*s)./sinh(theta_s);
            Rtheta=tanh(theta_s.*(s+0.5))./(2.0*tanh(0.5*theta_s))-0.5;
            C=(1.0-theta_b).*Ptheta+theta_b.*Rtheta;
          else
            C=s;
          end
        '''
        ds = 1.0/romsvars['N']
        if use_w:
            Nlev = romsvars['N']+1
            lev = np.arange(Nlev)
        else:
            Nlev = romsvars['N']
            lev = np.arange(Nlev)+0.5
        s = (lev-romsvars['N']) * ds
        if romsvars['theta_s'] > 0:
            Ptheta = np.sinh(romsvars['theta_s']*s)/np.sinh(romsvars['theta_s'])
            Rtheta = np.tanh(romsvars['theta_s']*(s+0.5))/(2.0*np.tanh(0.5*romsvars['theta_s']))-0.5
            C = (1.0-romsvars['theta_b'])*Ptheta+romsvars['theta_b']*Rtheta
        else:
            C = s

    elif (romsvars['Vstretching'] == 4) :
        ds = 1.0/romsvars['N']
        if use_w:
            Nlev = romsvars['N']+1
            lev = np.arange(Nlev)
        else:
            Nlev = romsvars['N']
            lev = np.arange(Nlev)+0.5
        s = (lev-romsvars['N']) * ds
        if (romsvars['theta_s'] > 0) :
            Csur = (1.0 - np.cosh(romsvars['theta_s']*s))/(np.cosh(romsvars['theta_s'])-1.0)
        else:
            Csur = s**2
        if (romsvars['theta_b'] > 0) :
            Cbot = (np.exp(romsvars['theta_b']*Csur) - 1.0)/(1.0 - np.exp(-romsvars['theta_b']))
            C = Cbot
        else:
           C = Csur

    elif (romsvars['Vstretching'] == 5) :
        if use_w:
           Nlev = romsvars['N']+1.0
           lev = np.arange(Nlev)

        else :
           Nlev = romsvars['N']
           lev = np.arange(Nlev)+0.5

        s = -1*(lev*lev - 2.0*lev*romsvars['N'] + lev + romsvars['N']**2.0 \
            - romsvars['N'])/(romsvars['N']**2 -romsvars['N']) \
            - 0.01*(lev*lev - lev*romsvars['N'])/(1 - romsvars['N'])

        if (romsvars['theta_s'] > 0) :
           Csur = (1 - np.cosh(romsvars['theta_s']*s))/(np.cosh(romsvars['theta_s']) - 1.0)
        else :
           Csur = s**2
        if (romsvars['theta_b'] > 0) :
           Cbot = (np.exp(romsvars['theta_b']*Csur) - 1.0) \
                   /(1.0 - np.exp(-1*romsvars['theta_b']))
           C = Cbot
        else :
           C = Csur

    else:
        raise NotImplementedError('Vstretching != 1 | 4 | 5  not supported yet.')
    return s, C
