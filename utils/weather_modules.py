#!/home/scavallo/anaconda2/bin/python

import numpy as np
from scipy import ndimage
from pstats import *
from numpy import trapz,linspace,logspace,log,log10
from mstats import *
import os, datetime

Cp = 1004.5;
Cv = 717.5;
Rd = 287.04;
Rv = 461.6;
RvRd = Rv / Rd;
g = 9.81;
L = 2.50e6;
Lf = 3.34*10**5;
Talt = 288.1500;
Tfrez = 273.1500;
To = 300;
Po = 101325;
Pr = 1000.;
lapsesta = 6.5 / 1000;
kappa = Rd / Cp;
epsil = Rd/Rv;
pi = 3.14159265;
pid = pi/180;
R_earth = 6371200;
omeg_e = (2*pi) / (24*3600);
eo = 6.11;
missval = -9999;
eps = 2.2204e-16

def temp_to_theta(temp, pres):
    ''' Compute potential temperature '''
    ''' '''
    ''' theta: Input potential temperature (K) '''
    ''' pres:  Input pressure (Pa)'''
    ''' temp:  Output temperature (K)'''
    return temp * (100000. / pres) ** 0.286

def theta_to_temp(theta, pres):
    ''' Compute temperature '''
    ''' '''
    ''' temp:  Input temperature (K)'''
    ''' pres:  Input pressure (Pa)'''
    ''' theta: Output potential temperature (K)'''
    return theta * (pres / 100000.) ** 0.286

def td_to_mixrat(tdew, pres):
    ''' Convert from dewpoint temperature to water vapor mixing ratio '''
    ''' '''
    ''' tdew:   Input dewpoint temperature (K)'''
    ''' pres:   Input pressure (Pa)'''
    ''' mixrat: Output water vapor mixing ratio (kg/kg)'''
    pres = pres/100
    mixrat = eo / (pres * RvRd)  * np.exp( (L/Rv)*((1/Tfrez) - (1 / tdew) ) )
    return mixrat

def mixrat_to_td(qvap, pres):
    ''' Convert from water vapor mixing ratio to dewpoint temperature '''
    ''' '''
    ''' qvap: Input water vapor mixing ratio (kg/kg)'''
    ''' pres: Input pressure (Pa)'''
    ''' tdew: Output dewpoint temperature (K)'''
    pres = pres/100
    evap = qvap * pres * RvRd;
    tdew = 1/((1/Tfrez) - (Rv/L)*np.log(evap/eo))
    return tdew

def spechum_to_td(spechum, pres):
    ''' Convert from water vapor mixing ratio to dewpoint temperature '''
    ''' '''
    ''' spechum: Input specific humidity in (kg/kg)'''
    ''' pres: Input pressure (Pa)'''
    ''' tdew: Output dewpoint temperature (K)'''
    qvap = (spechum/(1-spechum))
    
    pres = pres/100
    evap = qvap * pres * RvRd;
    tdew = 1/((1/Tfrez) - (Rv/L)*np.log(evap/eo))
    return tdew
    
def claus_clap(temp):
    ''' Compute saturation vapor pressure '''
    ''' '''
    ''' temp: Input temperature (K)  '''
    ''' esat: Output satuation vapor pressure (Pa)'''
    esat = (eo * np.exp( (L / Rv) * ( 1.0/Tfrez - 1/temp) ) ) * 100.
    return esat

def claus_clap_ice(temp):
    ''' Compute saturation vapor pressure over ice '''
    ''' '''
    ''' es:  Input satuation vapor pressure of liquid (Pa)'''    
    ''' temp: Input temperature (K) '''
    ''' esi: Output satuation vapor pressure of ice (Pa)'''

    a = 273.16 / temp
    exponent = -9.09718 * (a - 1.) - 3.56654 * np.log10(a) + 0.876793 * (1. - 1./a) + np.log10(6.1071)
    esi = 10**exponent
    esi = esi*100
    
    return esi

def sat_vap(temp):
    ''' Compute saturation vapor pressure '''
    ''' '''
    ''' temp: Input temperature (K)  '''
    ''' esat: Output satuation vapor pressure (Pa)'''
    ''' '''
    ''' Updated from claus_clap(temp) to account for phase automatically '''
    
    [iinds] = np.where(temp<273.15)
    [linds] = np.where(temp>=273.15)
    esat =  np.zeros_like(temp).astype('f')
    
    nice = len(iinds)
    nliq = len(linds)

    
    tempc = temp - 273.15    
    if nliq > 1:
        esat[linds] = 6.112*np.exp(17.67*tempc[linds]/(tempc[linds]+243.12))*100.
    else:
        if nliq > 0:
            esat = 6.112*np.exp(17.67*tempc/(tempc+243.12))*100.
    if nice > 1:
        esat[iinds] = 6.112*np.exp(22.46*tempc[iinds]/(tempc[iinds]+272.62))*100.
    else:
        if nice > 0:
            esat = 6.112*np.exp(22.46*tempc/(tempc+272.62))*100. 
    #esat = (eo * np.exp( (L / Rv) * ( 1.0/Tfrez - 1/temp) ) ) * 100.
    return esat

def moist_lapse(ws, temp):
    ''' Compute moist adiabatic lapse rate '''
    ''' '''
    ''' ws:   Input saturation mixing ratio (kg kg-1)  '''
    ''' temp: Input air temperature (K)'''
    ''' '''
    ''' moist_lapse: Output moist adiabatic lapse rate '''
    
    moist_lapse = (g/Cp)*((1.0 + L*ws)/(Rd*temp)) / ( 1.0 + (ws*(L**2.0)/(Cp*Rv*temp**2.0)) )
    return moist_lapse


def satur_mix_ratio(es, pres):
    ''' Compute saturation mixing ratio '''
    ''' '''
    ''' es:   Input saturation vapor pressure (Pa)  '''
    ''' pres: Input air pressure (Pa)'''
    ''' '''
    ''' ws: Output saturation mixing ratio '''
    
    ws = 0.622 * ( es / (pres - es) )
    return ws

def VirtualTempFromMixR(tempk,mixr):
    """Virtual Temperature

    INPUTS:
    tempk: Temperature (K)
    mixr: Mixing Ratio (kg/kg)

    OUTPUTS:
    tempv: Virtual temperature (K)
    """

    return tempk*(1.0+0.6*mixr)

def Latentc(tempk):
    """Latent heat of condensation (vapourisation)

    INPUTS:
    tempk (K)

    OUTPUTS:
    L_w (J/kg)

    SOURCE:
    http://en.wikipedia.org/wiki/Latent_heat#Latent_heat_for_condensation_of_water
    """
    tempc = tempk - 273.15
    return 1000*(2500.8 - 2.36*tempc + 0.0016*tempc**2 - 0.00006*tempc**3)

def total_col(infld, pres, temp, hght):
    ''' Compute column integrated value of infld '''
    ''' '''
    ''' infld:  Input 3D field to column integrate  '''
    ''' hght:  Input 3D geopotential height field (m)  '''
    ''' temp:  Input 3D temperature field (K)  '''
    ''' pres: Input 3D air pressure (Pa)'''
    ''' '''
    ''' tcol: Output total column integrated value '''
    
    [iz,iy,ix] = np.shape(infld)
    density = pres/(287*temp)
    
    tmp = pres[0,:,:].squeeze()
    
    coltot = np.zeros_like(tmp).astype('f')
    for jj in range(0,iy):      
       for ii in range(0,ix):              
           colnow = infld[:,jj,ii]*density[:,jj,ii]
           hghtnow = hght[:,jj,ii].squeeze()
           coltot[jj,ii] = trapz(colnow[::-1],hghtnow[::-1])          

    return coltot    

def GammaW(tempk,pres,e=None):
    """Function to calculate the moist adiabatic lapse rate (deg K/Pa) based
    on the temperature, pressure, and rh of the environment.

    INPUTS:
    tempk (K)
    pres (Pa)
    RH (%)

    RETURNS:
    GammaW: The moist adiabatic lapse rate (Dec K/Pa)
    """
    
    es=sat_vap(tempk)    
    ws=satur_mix_ratio(es,pres)
    

    if e is None:
        # assume saturated
        e=es
    
    w = satur_mix_ratio(e, pres)
    
    tempv=VirtualTempFromMixR(tempk,w)
    latent=Latentc(tempk)

    A=1.0+latent*ws/(Rd*tempk)
    B=1.0+epsil*latent*latent*ws/(Cp*Rd*tempk*tempk)
    Rho=pres/(Rd*tempv)
    Gamma=(A/B)/(Cp*Rho)
    return Gamma

def dry_parcel_ascent(startpp,starttk,starttdewk,nsteps=101):
    from numpy import interp
    #--------------------------------------------------------------------
    # Lift a parcel dry adiabatically from startp to LCL.
    #
    # Inputs:
    #     startp: Pressure of parcel to lift in Pa
    #     startt: Temperature of parcel at startp in K
    #     starttdew: Dewpoint temperature of parcel at startp in K
    #
    # Returns:
    #     presdry, tempdry: pressure (Pa) and temperature (K) along dry adiabatic ascent of parcel 
    #     tempiso is in K
    #     T_lcl, P_lcl: Temperature and pressure at LCL
    #--------------------------------------------------------------------

    assert starttdewk<=starttk
    
    
    startt = starttk - 273.15
    starttdew = starttdewk - 273.15
    startp = startpp/100.

    if starttdew==startt:
        return np.array([startp]),np.array([startt]),np.array([starttdew]),

    #Pres=np.linspace(startp,600)
    Pres=logspace(log10(startp),log10(600),nsteps)
    
    # Lift the dry parcel    
    T_dry=((starttk)*(Pres/startp)**(Rd/Cp))-273.15    

    # Mixing ratio isopleth  
    starte= sat_vap(starttdewk)  
    startw = satur_mix_ratio(starte, startpp)    
    #if starttdew < 0:
    #    starte = 6.112*np.exp(22.46*startt/(startt+272.62))*100.
    #else:
    #    starte = 6.112*np.exp(17.67*startt/(startt+243.12))*100.
    
    #startw=MixRatio(starte,startp*100)
    #startw = epsil*starte/((startp*100.)-starte)
    ee=Pres*startw/(.622+startw)    
    T_iso=243.5/(17.67/np.log(ee/6.112)-1.0)
   
    
    # Solve for the intersection of these lines (LCL).
    # interp requires the x argument (argument 2)
    # to be ascending in order!
    P_lcl=interp(0,T_iso-T_dry,Pres)    
    T_lcl=interp(P_lcl,Pres[::-1],T_dry[::-1])

    presdry=np.linspace(startp,P_lcl)
    tempdry=interp(presdry,Pres[::-1],T_dry[::-1])
    tempiso=interp(presdry,Pres[::-1],T_iso[::-1])
    
    return presdry*100.,tempdry+273.15,tempiso+273.15, T_lcl+273.15, P_lcl*100.

def moist_ascent(startpp,starttk,ptop=100,nsteps=501):
    #--------------------------------------------------------------------
    # Lift a parcel moist adiabatically from startp to endp.
    # Init temp is startt in Kelvin, pressure levels are in Pa    
    #--------------------------------------------------------------------
    
    startp = startpp/100. # convert to hPa
    startt = starttk - 273.15 # convert to deg C
    
    #preswet=np.linspace(startp,ptop,101) 
    preswet=logspace(log10(startp),log10(ptop),nsteps)
    
    temp=startt
    tempwet=np.zeros(preswet.shape);tempwet[0]=startt
    for ii in range(preswet.shape[0]-1):
        delp=preswet[ii]-preswet[ii+1]
        temp=temp-100.*delp*GammaW(temp+273.15,(preswet[ii]-delp/2)*100.)
        tempwet[ii+1]=temp

    return preswet*100.,tempwet+273.15

def thetae(thta, temp, qv):
    ''' Compute equivalent potential temperature '''
    ''' '''
    ''' thta:   Input potential temperature of column (K) '''
    ''' temp:   Input temperature (K) at LCL '''
    ''' qv:     Input mixing ratio of column (kg kg-1)'''
    ''' thetae: Output equivalent potential temperature (K)'''
    thout = thta * np.exp( (L * qv) / (Cp * temp) )
    return thout

def w_to_omega(w,pres,tempk):
    ''' Compute vertical velocity on isobaric surfaces '''
    ''' '''
    ''' w:   Input vertical velocity (m s-1) '''
    ''' pres:   Input half level pressures (full field) in Pa '''
    ''' tempk:   Input temperature (K) '''

    omeg = -((pres*g) / (Rd*tempk)) * w
    return omeg

def calc_gradient(fldin, dx, dy, dz):
    '''
    Computes the horizontal gradient of any scalar given a constant
    grid spacing in the x, y, and z directions.

    fldin: Input scalar
    dx: Input x grid spacing (must be single value)
    dy: Input y grid spacing (must be single value)
    dz: Input z grid spacing (must be single value)
    '''
    dfdx, dfdy, dfdz = np.gradient(fldin, dx, dy, dz)
    return dfdx, dfdy, dfdz

def latlon_to_dlatdlon(lats,lons):
    """
    Return arrays with the spacing between latitude and longitudes

    The gradients are computed using central differences in the interior
    and first differences at the boundaries. The returned gradient hence has
    the same shape as the input array.

    Parameters
    ----------
    lats : vector of latitudes in degrees
    lons : vector of longitudes in degrees

    Returns
    -------
    dlat : array with differences in latitudes between grid points with size (lats,lons)
    dlon : array with differences in longitudes between grid points with size (lats,lons)

    Examples
    --------
    dlat,dlon = latlon_to_dlatdlon(lats,lons)

    """

    nlat = len(lats)
    nlon = len(lons)
    latarr = np.zeros((nlat,nlon))
    lonarr = np.zeros((nlat,nlon))
    dlatarr = np.zeros((nlat,nlon))
    dlonarr = np.zeros((nlat,nlon))


    for jj in range(0,nlat):
       for ii in range(0,nlon):
          latarr[jj,ii] = lats[jj]
          lonarr[jj,ii] = lons[ii]


    latrad = latarr*(pi/180)

    # use central differences on interior and first differences on endpoints

    otype = latarr.dtype.char
    if otype not in ['f', 'd', 'F', 'D']:
        otype = 'd'

    dlats = np.zeros_like(lats).astype(otype)
    dlats[1:-1] = (lats[2:] - lats[:-2])
    dlats[0] = (lats[1] - lats[0])
    dlats[-1] = (dlats[-2] - dlats[-1])

    dlons = np.zeros_like(lons).astype(otype)
    dlons[1:-1] = (lons[2:] - lons[:-2])
    dlons[0] = (lons[1] - lons[0])
    dlons[-1] = (dlons[-2] - dlons[-1])

    # Since we differenced in the reverse direction, change the sign
    dlats = -1*dlats

    for jj in range(0,nlat):
       for ii in range(0,nlon):
          dlonarr[jj,ii] = dlons[ii]
          dlatarr[jj,ii] = dlats[jj]


    return dlatarr, dlonarr

def gradient_cartesian(f, *varargs):
    """
    Return the gradient of an N-dimensional array on an evenly spaced grid.

    The gradient is computed using central differences in the interior
    and first differences at the boundaries. The returned gradient hence has
    the same shape as the input array.

    Parameters
    ----------
    f : N-dimensional array containing samples of a scalar function.
        If 2-D, must be ordered as f(y,x)
    If 3-D, must be ordered as f(z,y,x) or f(p,y,x)
    `*varargs` : scalars
          0, 1, or N scalars specifying the sample distances in each direction,
          that is: `dz`, `dy`, `dx`, ... The default distance is 1.

      If a vector is specified as the first argument of three, then the difference
         of this vector will be taken here.


    Returns
    -------
    g : ndarray
      N arrays of the same shape as `f` giving the derivative of `f` with
      respect to each dimension.

    Examples
    --------
    fldin = temperature(pressure,y,x)
    levs = pressure vector
    dy = scalar or array of grid spacing in y direction
    dx = scalar or vector of grid spacing in x direction
    >>> dfdz, dfdy = gradient_cartesian(fldin, dy, dx)
    >>> dfdz, dfdy, dfdx = gradient_cartesian(fldin, dz, dy, dx)
    >>> dfdp, dfdy, dfdx = gradient_cartesian(fldin, levs, dy, dx)
    
    
    Steven Cavallo
    University of Oklahoma
    July 2015
    """
    N = len(f.shape)  # number of dimensions
    n = len(varargs)
    argsin = list(varargs)

    if N != n:
       raise SyntaxError("dimensions of input array must match the number of remaining arguments")

    df = np.gradient(f)

    if n == 1:
        dy = argsin[0]
        dfdy = df[0]
    elif n == 2:
        dy = argsin[0]
        dx = argsin[1]
        dfdy = df[0]
        dfdx = df[1]
    elif n == 3:
        levs = argsin[0]
        dy = argsin[1]
        dx = argsin[2]
        dfdz = df[0]
        dfdy = df[1]
        dfdx = df[2]
    else:
        raise SyntaxError(
                "invalid number of arguments")

    otype = f.dtype.char
    if otype not in ['f', 'd', 'F', 'D']:
        otype = 'd'


    try:
       M = len(dx.shape)
    except:
       M = 1
    dyarr = np.zeros_like(f).astype(otype)
    dxarr = np.zeros_like(f).astype(otype)
    if M == 1:
       dyarr[:] = dy
       dxarr[:] = dx
       if N == 1:
          ny = np.shape(f)
       elif N == 2:
          ny, nx = np.shape(f)
       else:
          nz, ny, nx = np.shape(f)

    else:
       if N == 1:
          ny = np.shape(f)
          for jj in range(0,ny):
              dyarr[jj,ii] = dy[jj]
       elif N == 2:
           ny, nx = np.shape(f)
           for jj in range(0,ny):
               for ii in range(0,nx):
                   dyarr[jj,ii] = dy[jj]
                   dxarr[jj,ii] = dx[ii]
       else:
           nz, ny, nx = np.shape(f)
           for kk in range(0,nz):
               for jj in range(0,ny):
                   for ii in range(0,nx):
                       dyarr[kk,jj,ii] = dy[jj]
                       dxarr[kk,jj,ii] = dx[ii]

    if n==1:
       dfdy = dfdy/dx

       return dfdy
    elif n==2:
       dfdy = dfdy/dy
       dfdx = dfdx/dx

       return dfdy,dfdx
    elif n==3:
       dfdy = dfdy/dy
       dfdx = dfdx/dx

       nzz = np.shape(levs)
       
       if not nzz:
          nzz=0

       if nzz>1:
           zin = levs
           dz = np.zeros_like(zin).astype(otype)
           dz[1:-1] = (zin[2:] - zin[:-2])/2
           dz[0] = (zin[1] - zin[0])
           dz[-1] = (zin[-1] - zin[-2])                     
       if zin[1] < zin[0]:
           dz = dz*-1 # assume the model top is the first index and the lowest model is the last index
           dx3 = np.ones_like(f).astype(otype)
           for kk in range(0,nz):
               dx3[kk,:,:] = dz[kk]
       else:
           dx3 = np.ones_like(f).astype(otype)
           dx3[:] = dx[0]

       dfdz = dfdz/dx3
       return dfdz,dfdy,dfdx

def gradient_sphere(f, *varargs):
    """
    Return the gradient of a 2-dimensional array on a sphere given a latitude
    and longitude vector.

    The gradient is computed using central differences in the interior
    and first differences at the boundaries. The returned gradient hence has
    the same shape as the input array.

    Parameters
    ----------
    f : A 2-dimensional array containing samples of a scalar function.
    latvec: latitude vector
    lonvec: longitude vector

    Returns
    -------
    g : dfdx and dfdy arrays of the same shape as `f` giving the derivative of `f` with
        respect to each dimension.

    Examples
    --------
    temperature = temperature(pressure,latitude,longitude)
    levs = pressure vector
    lats = latitude vector
    lons = longitude vector
    >>> tempin = temperature[5,:,:]
    >>> dfdlat, dfdlon = gradient_sphere(tempin, lats, lons)

    >>> dfdp, dfdlat, dfdlon = gradient_sphere(temperature, levs, lats, lons)

    based on gradient function from /usr/lib64/python2.6/site-packages/numpy/lib/function_base.py
    """

    R_earth = 6371200.
    N = len(f.shape)  # number of dimensions
    n = len(varargs)
    argsin = list(varargs)

    if N != n:
       raise SyntaxError("dimensions of input array must match the number of remaining arguments")

    df = np.gradient(f)

    if n == 1:
        lats = argsin[0]
        dfdy = df[0]
    elif n == 2:
        lats = argsin[0]
        lons = argsin[1]
        dfdy = df[0]
        dfdx = df[1]
    elif n == 3:
        levs = argsin[0]
        lats = argsin[1]
        lons = argsin[2]
        dfdz = df[0]
        dfdy = df[1]
        dfdx = df[2]
    else:
        raise SyntaxError(
                "invalid number of arguments")

    otype = f.dtype.char
    if otype not in ['f', 'd', 'F', 'D']:
        otype = 'd'

    latarr = np.zeros_like(f).astype(otype)
    lonarr = np.zeros_like(f).astype(otype)
    if N == 1:
       nlat = np.shape(f)
       for jj in range(0,nlat):
          latarr[jj,ii] = lats[jj]
       lonarr = latarr
       lons = lats
    elif N == 2:
       nlat, nlon = np.shape(f)
       for jj in range(0,nlat):
          for ii in range(0,nlon):
             latarr[jj,ii] = lats[jj]
             lonarr[jj,ii] = lons[ii]
    else:
       nz, nlat, nlon = np.shape(f)
       for kk in range(0,nz):
          for jj in range(0,nlat):
             for ii in range(0,nlon):
                latarr[kk,jj,ii] = lats[jj]
                lonarr[kk,jj,ii] = lons[ii]

    latrad = latarr*(pi/180)

    # use central differences on interior and first differences on endpoints

    outvals = []

    dlats = np.zeros_like(lats).astype(otype)
    dlats[1:-1] = (lats[2:] - lats[:-2])
    dlats[0] = (lats[1] - lats[0])
    dlats[-1] = (dlats[-2] - dlats[-1])

    dlons = np.zeros_like(lons).astype(otype)
    dlons[1:-1] = (lons[2:] - lons[:-2])
    dlons[0] = (lons[1] - lons[0])
    dlons[-1] = (dlons[-2] - dlons[-1])

    # Since we differenced in the reverse direction, change the sign
    #dlats = -1*dlats

    dlatarr = np.tile(dlats,[nlon,1])
    dlatarr = np.reshape(dlatarr,[nlat,nlon])

    dlonarr = np.zeros_like(f).astype(otype)
    if N==2:
       for jj in range(0,nlat):
          for ii in range(0,nlon):
              dlonarr[jj,ii] = dlons[ii]
    elif N==3:
       for kk in range(0,nz):
          for jj in range(0,nlat):
             for ii in range(0,nlon):
                 dlonarr[kk,jj,ii] = dlons[ii]

    dlatsrad = dlatarr*(pi/180)
    dlonsrad = dlonarr*(pi/180)
    latrad = latarr*(pi/180)

    if n==1:
       dx1 = R_earth * dlatsrad
       dfdy = dfdy/dx1

       return dfdy
    elif n==2:
       dx1 = R_earth * dlatsrad
       dx2 = R_earth * np.cos(latrad) * dlonsrad

       dfdy = dfdy/dx1
       dfdx = dfdx/dx2

       return dfdy,dfdx
    elif n==3:
       dx1 = R_earth * dlatsrad
       dx2 = R_earth * np.cos(latrad) * dlonsrad

       dfdy = dfdy/dx1
       dfdx = dfdx/dx2

       nzz = np.shape(levs)
       if not nzz:
          nzz=0

       if nzz>1:
           zin = levs
           dz = np.zeros_like(zin).astype(otype)
           dz[1:-1] = (zin[2:] - zin[:-2])/2
           dz[0] = (zin[1] - zin[0])
           dz[-1] = (zin[-1] - zin[-2])
           if zin[0,1,1] > zin[1,1,1]:
               dz = dz*-1 # assume the model top is the first index and the lowest model is the last index

           dx3 = np.ones_like(f).astype(otype)
           for kk in range(0,nz):
               dx3[kk,:,:] = dz[kk]
       else:
           dx3 = np.ones_like(f).astype(otype)
           dx3[:] = dx[0]

       dfdz = dfdz/dx3
       return dfdz,dfdy,dfdx


def gradient_cendiff_latlon(f, *varargs):
    """
    Return the gradient of a 2-dimensional array given a latitude
    and longitude vector.
    
    Does NOT take into account sphereical geometry

    The gradient is computed using central differences in the interior
    and first differences at the boundaries. The returned gradient hence has
    the same shape as the input array.

    Parameters
    ----------
    f : A 2-dimensional array containing samples of a scalar function.
    latvec: latitude vector
    lonvec: longitude vector

    Returns
    -------
    g : dfdx and dfdy arrays of the same shape as `f` giving the derivative of `f` with
        respect to each dimension.

    Examples
    --------
    temperature = temperature(pressure,latitude,longitude)
    levs = pressure vector
    lats = latitude vector
    lons = longitude vector
    >>> tempin = temperature[5,:,:]
    >>> dfdlat, dfdlon = gradient_sphere(tempin, lats, lons)

    >>> dfdp, dfdlat, dfdlon = gradient_sphere(temperature, levs, lats, lons)

    based on gradient function from /usr/lib64/python2.6/site-packages/numpy/lib/function_base.py
    """

    R_earth = 6371200.
    N = len(f.shape)  # number of dimensions
    n = len(varargs)
    argsin = list(varargs)

    if N != n:
       raise SyntaxError("dimensions of input array must match the number of remaining argumens")

    df = np.gradient(f)

    if n == 1:
        lats = argsin[0]
        dfdy = df[0]
    elif n == 2:
        lats = argsin[0]
        lons = argsin[1]
        dfdy = df[0]
        dfdx = df[1]
    elif n == 3:
        levs = argsin[0]
        lats = argsin[1]
        lons = argsin[2]
        dfdz = df[0]
        dfdy = df[1]
        dfdx = df[2]
    else:
        raise SyntaxError(
                "invalid number of arguments")

    otype = f.dtype.char
    if otype not in ['f', 'd', 'F', 'D']:
        otype = 'd'

    latarr = np.zeros_like(f).astype(otype)
    lonarr = np.zeros_like(f).astype(otype)
    if N == 1:
       nlat = np.shape(f)
       for jj in range(0,nlat):
          latarr[jj,ii] = lats[jj]
       lonarr = latarr
       lons = lats
    elif N == 2:
       nlat, nlon = np.shape(f)
       for jj in range(0,nlat):
          for ii in range(0,nlon):
             latarr[jj,ii] = lats[jj]
             lonarr[jj,ii] = lons[ii]
    else:
       nz, nlat, nlon = np.shape(f)
       for kk in range(0,nz):
          for jj in range(0,nlat):
             for ii in range(0,nlon):
                latarr[kk,jj,ii] = lats[jj]
                lonarr[kk,jj,ii] = lons[ii]

    latrad = latarr*(pi/180)

    # use central differences on interior and first differences on endpoints

    outvals = []

    dlats = np.zeros_like(lats).astype(otype)
    dlats[1:-1] = (lats[2:] - lats[:-2])
    dlats[0] = (lats[1] - lats[0])
    dlats[-1] = (dlats[-2] - dlats[-1])

    dlons = np.zeros_like(lons).astype(otype)
    dlons[1:-1] = (lons[2:] - lons[:-2])
    dlons[0] = (lons[1] - lons[0])
    dlons[-1] = (dlons[-2] - dlons[-1])

    # Since we differenced in the reverse direction, change the sign
    #dlats = -1*dlats

    dlatarr = np.tile(dlats,[nlon,1])
    dlatarr = np.reshape(dlatarr,[nlat,nlon])

    dlonarr = np.zeros_like(f).astype(otype)
    if N==2:
       for jj in range(0,nlat):
          for ii in range(0,nlon):
             dlonarr[jj,ii] = dlons[ii]
    elif N==3:
       for kk in range(0,nz):
          for jj in range(0,nlat):
             for ii in range(0,nlon):
                 dlonarr[kk,jj,ii] = dlons[ii]

    dlatsrad = dlatarr*(pi/180)
    dlonsrad = dlonarr*(pi/180)
    latrad = latarr*(pi/180)

    if n==1:
       dx1 = dlatsrad
       dfdy = dfdy/dx1

       return dfdy
    elif n==2:
       dx1 = dlatsrad
       dx2 = dlonsrad

       dfdy = dfdy/dx1
       dfdx = dfdx/dx2

       return dfdy,dfdx
    elif n==3:
       dx1 = dlatsrad
       dx2 = dlonsrad

       dfdy = dfdy/dx1
       dfdx = dfdx/dx2

       nzz = np.shape(levs)
       if not nzz:
          nzz=0

       if nzz>1:
           zin = levs
           dz = np.zeros_like(zin).astype(otype)
           dz[1:-1] = (zin[2:] - zin[:-2])/2
           dz[0] = (zin[1] - zin[0])
           dz[-1] = (zin[-1] - zin[-2])
           if zin[0,1,1] > zin[1,1,1]:
               dz = dz*-1 # assume the model top is the first index and the lowest model is the last index

           dx3 = np.ones_like(f).astype(otype)
           for kk in range(0,nz):
               dx3[kk,:,:] = dz[kk]
       else:
           dx3 = np.ones_like(f).astype(otype)
           dx3[:] = dx[0]

       dfdz = dfdz/dx3
       return dfdz,dfdy,dfdx




def _get_gradients(u, v, dx, dy):
    #Helper function for getting convergence and vorticity from 2D arrays
    dudx, dudy = np.gradient(u, dx, dy)
    dvdx, dvdy = np.gradient(v, dx, dy)
    return dudx, dudy, dvdx, dvdy

def vertical_vorticity(u, v, dx, dy, grid_opt):
    '''
Calculate the vertical vorticity of the horizontal wind. The grid
must have a constant spacing in each direction.

u, v : 2 dimensional arrays
Arrays with the x and y components of the wind, respectively.
X must be the first dimension and y the second.

dx : scalar or array
The grid spacing in the x-direction

dy : scalar or array
The grid spacing in the y-direction

grid_opt: 1 for cartesian grid
          2 for lat/lon grid

Returns : 2 dimensional array
The vertical vorticity
'''

    if grid_opt == 1:
       dudy,dudx = gradient_cartesian(u, dy, dx)
       dvdy,dvdx = gradient_cartesian(v, dy, dx)
    else:
       dudy,dudx = gradient_sphere(u, dy, dx)
       dvdy,dvdx = gradient_sphere(v, dy, dx)

    return dvdx - dudy

def h_convergence(u, v, dx, dy):
    '''
Calculate the horizontal convergence of the horizontal wind. The grid
must have a constant spacing in each direction.

u, v : 2 dimensional arrays
Arrays with the x and y components of the wind, respectively.
X must be the first dimension and y the second.

dx : scalar
The grid spacing in the x-direction

dy : scalar
The grid spacing in the y-direction

Returns : 2 dimensional array
The horizontal convergence
'''
    dudx, dudy, dvdx, dvdy = _get_gradients(u, v, dx, dy)
    return dudx + dvdy

def geostrophic_latlon(u, v, ghgt, lats, lons):
    '''
Calculate geostrophic wind on a latitude/longitude grid

Input:
   
   u(lats,lons), v(lats,lons) : 2 dimensional u and v wind arrays 
                             dimensioned by (lats,lons).
                             Arrays correspond to the x and y 
                 components of the wind, respectively.
                 
   ghgt(lats,lons): 2 dimensional array of geopotential height
   
   lats(lats) : latitude vector
   lons(lons) : longitude array

Output:

   ug(lats,lons), vg(lats,lons): 2 dimensional geostrophic u and v 
                             wind arrays dimensioned by (lats,lons).
                             Arrays correspond to the x and y 
                 components of the geostrophic wind, 
                 respectively.
'''
    # 2D latitude array
    glats = np.zeros_like(u).astype('f')      
    for jj in range(0,len(lats)):
        for ii in range(0,len(lons)):    
            glats[jj,ii] = lats[jj]

    # Coriolis parameter
    f = 2*(7.292e-05)*np.sin(np.deg2rad(glats))    

    ghgt = ndimage.gaussian_filter(ghgt,0.75)
    
    geop = ghgt * 9.81
    dphidy,dphidx = gradient_sphere(geop, lats, lons)
    
    ug = -dphidy/f
    vg = dphidx/f
    
    return ug, vg
    
def geostrophic_cartesian(u, v, ghgt, lats, dy, dx):
    '''
Calculate geostrophic wind on a latitude/longitude grid

Input:
   
   u,v : 2 dimensional u and v wind arrays 
         corresponding to the x and y 
     components of the wind, respectively.
                 
   ghgt: 2 dimensional array of geopotential height
   
   lats : latitude array

Output:

   ug,vg : 2 dimensional geostrophic u and v 
           wind arrays corresponding to the x and y 
       components of the geostrophic wind, respectively.
'''

    # Coriolis parameter
    f = 2*(7.292e-05)*np.sin(np.deg2rad(lats))    
    ghgt = ndimage.gaussian_filter(ghgt,0.75)
    
    geop = ghgt * 9.81
    dummy,dphidy,dphidx = gradient_cartesian(geop, geop[:,0,0], dy, dx)
    
    ug = -dphidy/f
    vg = dphidx/f    
    
    return ug, vg
def hadvection_cartesian(datain, u, v, deltax, deltay):
    '''
Calculate the vertical vorticity on a Cartesian grid

Input:
   datain(y,x)    : 2 dimensional array to compute advection of
   u(y,x), v(y,x) : 2 dimensional u and v wind arrays                                 
                    Arrays correspond to the x and y 
            components of the wind, respectively.
   deltax     : horizontal x grid spacing in meters
   deltay     : horiztional y grid spacing in meters

Output: 

   dataout(y,x): Two dimensional array with the horizontal
         advection of datain
'''
    sh = np.shape(datain)
    n = len(sh)
    
    if n == 2:
        datady,datadx = gradient_cartesian(datain, deltay, deltax)
    else:
        datadz,datady,datadx = gradient_cartesian(datain, datain[:,0,0], deltay, deltax)        
    
    dataout = -u*datadx -v*datady
    
    return dataout
    
def hadvection_latlon(u, v, datain, lats, lons):
    '''
Calculate horizontal advection on a latitude/longitude grid

Input:

   datain(lats,lons): 2 dimensional array to compute advection of

   u(lats,lons), v(lats,lons) : 2 dimensional u and v wind arrays 
                             dimensioned by (lats,lons).
                             Arrays correspond to the x and y 
                 components of the wind, respectively.

   lats(lats) : latitude vector
   lons(lons) : longitude array

Output:

   dataout(lats,lons): Two dimensional array with the
                        horizontal advection of data
'''
    datady,datadx = gradient_sphere(datain, lats, lons)
    dataout = -u*datadx -v*datady
    
    return dataout
def vertical_vorticity_latlon(u, v, lats, lons, abs_opt):
    '''
Calculate the vertical vorticity on a latitude/longitude grid

Input:
   u(lats,lons), v(lats,lons) : 2 dimensional u and v wind arrays 
                                dimensioned by (lats,lons).
                                Arrays correspond to the x and y 
                    components of the wind, respectively.

   lats(lats) : latitude vector
   lons(lons) : longitude array

   abs_opt: 1 to compute absolute vorticity
            0 for relative vorticity only

Output: 

   vert_vort(lats,lons): Two dimensional array of vertical voriticity
'''
    dudy,dudx = gradient_sphere(u, lats, lons)
    dvdy,dvdx = gradient_sphere(v, lats, lons)
    
    if abs_opt == 1 :
       # 2D latitude array
       glats = np.zeros_like(u).astype('f')      
       for jj in range(0,len(lats)):
           for ii in range(0,len(lons)):    
               glats[jj,ii] = lats[jj]

       # Coriolis parameter
       f = 2*(7.292e-05)*np.sin(np.deg2rad(glats))    
       
    else:   
       f = 0.
    
    zeta = dvdx - dudy
    vert_vort = zeta + f  

    return vert_vort

def vertical_vorticity_cartesian(u, v, lats, deltax, deltay, abs_opt):
    '''
Calculate the vertical vorticity on a Cartesian grid

Input:
   u(y,x), v(y,x) : 2 dimensional u and v wind arrays                                 
                    Arrays correspond to the x and y 
            components of the wind, respectively.

   lats(lats) : 2D latitude array
   deltax     : horizontal x grid spacing in meters
   deltay     : horiztional y grid spacing in meters

   abs_opt: 1 to compute absolute vorticity
            0 for relative vorticity only

Output: 

   vert_vort(y,x): Two dimensional array of vertical voriticity
'''
    dudy,dudx = gradient_cartesian(u, deltay, deltax)
    dvdy,dvdx = gradient_cartesian(v, deltay, deltax)
        
    iy, ix = u.shape
    
    if abs_opt == 1 :
       # Coriolis parameter
       f = 2*(7.292e-05)*np.sin(np.deg2rad(lats))    
       
    else:   
       f = 0.
    
    zeta = dvdx - dudy
    vert_vort = zeta + f  

    return vert_vort

def thermal_wind_cartesian(thickness_in, lats, deltax, deltay):
    
    '''Calculate the thermal wind on a cartesian grid

    Input: thickness_in(deltax, deltay) : 2 dimensional geopotential height thickness array
    deltax     : horizontal x grid spacing in meters
    deltay     : horiztional y grid spacing in meters
    lats       : 2D latitude array used for calculating coriolis parameter
    Output: thermal_wind_u(lats,lons), thermal_wind_v(lats,lons): Two dimensional arryas of u and v wind
    '''
    sh = np.shape(thickness_in)
    n = len(sh)

    # Smooth the thicknesses
    thickness_in = ndimage.gaussian_filter(thickness_in,0.75)

    if n == 2:    
        dthickdy,dthickdx = gradient_cartesian(thickness_in, deltay, deltax)
    else:         
        dthickdz,dthickdy,dthickdx = gradient_cartesian(thickness_in, thickness_in[:,0,0], deltay, deltax)    
    
    #Calculate the Coriolis parameter

    f = 2*(7.292e-05)*np.sin(np.deg2rad(lats)) 

    thermal_wind_u = -(9.81/f)*dthickdy
    thermal_wind_v = (9.81/f)*dthickdx

    return thermal_wind_u, thermal_wind_v

def thermal_wind_sphere(thickness_in, lats, lons):
    '''
Calculate the thermal wind on a latitude/longitude grid

Input:
   thickness_in(lats,lons) : 2 dimensional geopotential height thickness array 
                             dimensioned by (lats,lons).


   lats(lats) : latitude vector
   lons(lons) : longitude array


Output: 

   thermal_wind_u(lats,lons), thermal_wind_v(lats,lons): Two dimensional arrays of 
                                                         u- and v- components of 
                             thermal wind vector
'''
    
    # Smooth the thicknesses
    thickness_in = ndimage.gaussian_filter(thickness_in,0.75)
    
    dthickdy,dthickdx = gradient_sphere(thickness_in, lats, lons)

       
    # 2D latitude array
    glats = np.zeros_like(thickness_in).astype('f')      
    for jj in range(0,len(lats)):
        for ii in range(0,len(lons)):    
            glats[jj,ii] = lats[jj]

    # Coriolis parameter
    f = 2*(7.292e-05)*np.sin(np.deg2rad(glats))    
    
           
    thermal_wind_u = -1*(9.81/f)*dthickdy
    thermal_wind_v = (9.81/f)*dthickdx

    return thermal_wind_u, thermal_wind_v

def eliassen_palm_flux_sphere(geop,theta,lats,lons,levs,normalize_option):
    """

   Computes the 3-D Eliassen-Palm flux vectors and divergence on a 
   latitude/longitude grid with any vertical coordinate 
   
   Computation is Equation 5.7 from:
   R. A. Plumb, On the Three-dimensional Propagation of Stationary Waves, J. Atmos. Sci., No. 3, 42 (1985).

   Input:    
       geop:      3D geopotential (m^2 s-2)
       theta:     3D potential temperature (K)
       lats,lons: 1D latitude and longitude vectors
       levs:      1D pressure vector (Pa)
       normalize_option: 0 = no normalizing, 1 = normalize by maximum value at each vertical level, 2 = standardize 

   Output:      
      Fx, Fy, Fz: Eliassen-Palm flux x, y, z vector components
      divF: Eliassen-Palm flux divergence
   

   Steven Cavallo
   March 2014
   University of Oklahoma
    
    """
    import utilities_modules as um
    
    iz, iy, ix = geop.shape
    
    # First need to filter out numeric nans    
    geop = um.filter_numeric_nans(geop,0,0,'low')    
    theta = um.filter_numeric_nans(theta,200,0,'low')
    theta = um.filter_numeric_nans(theta,10000,0,'high')

    theta_anom, theta_anom_std = spatial_anomaly(theta,1) 
    geop_anom , geop_anom_std = spatial_anomaly(geop,1)          
        
    latarr = np.zeros_like(geop).astype('f')    
    farr = np.zeros_like(geop).astype('f')    
    for kk in range(0,iz):
        for jj in range(0,iy):            
            latarr[kk,jj,:] = lats[jj]
            farr[kk,jj,:] = 2.0*omeg_e*np.sin(lats[jj]*(np.pi/180.0))
          
    psi = geop / farr     
    psi_anom, psi_anom_std = spatial_anomaly(psi,1) 
        
    pres = np.zeros_like(theta_anom).astype('f')   
    for kk in range(0,iz):      
        pres[kk,:,:] = levs[kk]
    
    coef = pres*np.cos(latarr*(np.pi/180.0))    
    arg1 = coef/( 2.0*np.pi*(R_earth**2)*np.cos(latarr*(np.pi/180.0))*np.cos(latarr*(np.pi/180.0)) )
    arg2 = coef/( 2.0*np.pi*(R_earth**2)*np.cos(latarr*(np.pi/180.0)) )
    arg3 = coef*(2.0*(omeg_e**2)*np.sin(latarr*(np.pi/180.0))*np.sin(latarr*(np.pi/180.0)))
   
    dthdz, yy, xx = gradient_cendiff_latlon(theta, geop/g, lats, lons) 
    xx, dthdy, dthdx = gradient_cendiff_latlon(theta_anom, geop/g, lats, lons)        
    dpsidz, dpsidy, dpsidx = gradient_cendiff_latlon(psi_anom, geop/g, lats, lons)          
    d2psidxdz, d2psidxdy, d2psidx2 = gradient_cendiff_latlon(dpsidx, geop/g, lats, lons)     
    aaa, d2psidxdz, ccc = gradient_cendiff_latlon(dpsidz, geop/g, lats, lons) 

    N2 = (g/theta)*(dthdz)
    arg4 = arg3/(N2*R_earth*np.cos(latarr*(np.pi/180.0)))        
    
    Fx = arg1*( dpsidx**2.0      - (psi_anom*d2psidx2))
    Fy = arg2*((dpsidy*dpsidx) - (psi_anom*d2psidxdy))
    Fz = arg4*((dpsidx*dpsidz) - (psi_anom*d2psidxdz))          
    
    Fx_z, Fx_y, Fx_x = gradient_sphere(Fx, geop/g, lats, lons)    
    Fy_z, Fy_y, Fy_x = gradient_sphere(Fy, geop/g, lats, lons)    
    Fz_z, Fz_y, Fz_x = gradient_sphere(Fz, geop/g, lats, lons)    
       
    divF = Fx_x + Fy_y + Fz_z
    
    if normalize_option == 1:        
        for kk in range(0,iz):  
	    Fx[kk,:,:] =  Fx[kk,:,:] / np.nanmax(np.abs(Fx[kk,:,:]))
	    Fy[kk,:,:] =  Fy[kk,:,:] / np.nanmax(np.abs(Fy[kk,:,:]))
	    Fz[kk,:,:] =  Fz[kk,:,:] / np.nanmax(np.abs(Fz[kk,:,:]))
    if normalize_option == 2:        
        for kk in range(0,iz):  
	    Fx[kk,:,:] =  Fx[kk,:,:] / np.nanstd(Fx[kk,:,:])
	    Fy[kk,:,:] =  Fy[kk,:,:] / np.nanstd(Fy[kk,:,:])
	    Fz[kk,:,:] =  Fz[kk,:,:] / np.nanstd(Fz[kk,:,:])    

    return Fx, Fy, Fz, divF


def epv_sphere(theta,pres,u,v,lats,lons):
    """

   Computes the Ertel Potential Vorticity (PV) on a latitude/longitude grid

   Input:    

       theta:       3D potential temperature array on isobaric levels
       pres:        3D pressure array
       u,v:         3D u and v components of the horizontal wind on isobaric levels
       lats,lons:   1D latitude and longitude vectors

   Output:
      
      epv: Ertel PV in potential vorticity units (PVU)
   

   Steven Cavallo
   October 2012
   University of Oklahoma
    
    """
    iz, iy, ix = theta.shape
    

    dthdp, dthdy, dthdx = gradient_sphere(theta, pres, lats, lons)
    dudp, dudy, dudx = gradient_sphere(u, pres, lats, lons)
    dvdp, dvdy, dvdx = gradient_sphere(v, pres, lats, lons)    

    avort = np.zeros_like(theta).astype('f')   
    for kk in range(0,iz):       
        avort[kk,:,:] = vertical_vorticity_latlon(u[kk,:,:].squeeze(), v[kk,:,:].squeeze(), lats, lons, 1)

    epv = (-9.81*(-dvdp*dthdx - dudp*dthdy + avort*dthdp))*10**6


    return epv
    
def epv_cartesian(theta,pres,u,v,lats,deltax,deltay):
    """

   Computes the Ertel Potential Vorticity (PV) on a Cartesian grid

   Input:    

       theta:       3D potential temperature array on isobaric levels
       pres:        1D pressure vector
       u,v:         3D u and v components of the horizontal wind on isobaric levels
       lats:        2D latitude array
       deltax, deltay: x and y horizontal grid spacing in meters

   Output:
      
      epv: Ertel PV in potential vorticity units (PVU)
   

   Steven Cavallo
   October 2012
   University of Oklahoma
    
    """
    iz, iy, ix = theta.shape
    
    dthdp, dthdy, dthdx = gradient_cartesian(theta, pres, deltay, deltax)
    dudp, dudy, dudx = gradient_cartesian(u, pres, deltay, deltax)
    dvdp, dvdy, dvdx = gradient_cartesian(v, pres, deltay, deltax)

    avort = np.zeros_like(theta).astype('f')   
    for kk in range(0,iz):       
        avort[kk,:,:] = vertical_vorticity_cartesian(u[kk,:,:].squeeze(), v[kk,:,:].squeeze(), lats, deltax, deltay, 1)

    epv = (-9.81*(-dvdp*dthdx - dudp*dthdy + avort*dthdp))*10**6


    return epv


    
def interp2pv(pv, fval, pv_surf):
    """

   Linearly interpolates a field to a PV surface

   Input:    

       pv - 3D array that contains the PV at common levels (P or theta)
       fval - 3D field to interpolate
       pv_surf - potential vorticity surface to interpolate onto

    Steven Cavallo
    October 2012
    University of Oklahoma
    
    """
    iz, iy, ix = pv.shape
    
    # Scan  from the top of the model downward.
    # The zeroth index is assumed to correspond to the top of the model.
    trop = fval[0,:,:].squeeze() 

    for jj in range(iy):
       for ii in range(ix):  

           aa = np.ravel(pv[:,jj,ii]>pv_surf)
           pvcol = pv[:,jj,ii].squeeze() 
           minpv = np.min(pvcol)

           if ( minpv >= pv_surf ):
               # If there are no PV values in the column less than what is desired to interpolate onto, then use value closest to the surface
               trop[jj,ii] = fval[-1,jj,ii]
           elif ( pv[0,jj,ii] <= pv_surf ):
               # If PV at the model top is less than what is desired to interpolate onto, then use the value at the top of the model
               trop[jj,ii] = fval[0,jj,ii]
           else:               
               for kk in range(1,iz):      
                   # linearly interpolate between the closest levels
                   if pv[kk,jj,ii] < pv_surf:
                       m = (fval[kk-1,jj,ii] - fval[kk,jj,ii]) / (pv[kk-1,jj,ii] - pv[kk,jj,ii])
                       trop[jj,ii] = m * (pv_surf - pv[kk,jj,ii]) + fval[kk,jj,ii]
                       break


    return trop

def spatial_anomaly(varin,slice_option):
    """

   Computes the spatial anomaly of varin

   Input:    

       varin:        3D array of variable to compute anomaly of
       slice_option: 1 to compute anomaly of second dimension
                     2 to compute anomaly of third dimension

   Output:
      
      varanom: Anomaly of varin
      varanom_std = Standardized anomaly of varin
   

   Steven Cavallo
   March 2014
   University of Oklahoma
    
    """
    iz, iy, ix = varin.shape
        
    mvar = np.ma.masked_array(varin,np.isnan(varin)) 
    
    tmp = np.zeros_like(varin).astype('f')  
    tmp_std = np.zeros_like(varin).astype('f')  
    
    if slice_option == 1:
        var_mean = np.mean(mvar,2)
        var_std = np.std(mvar,2)
        for kk in range(0,iz): 
            for jj in range(0,iy):     
                tmp[kk,jj,:] = varin[kk,jj,:] - var_mean[kk,jj]    
                tmp_std[kk,jj,:] = var_std[kk,jj]
    else:
        var_mean = np.mean(mvar,1)
        var_std = np.std(mvar,1)
        for kk in range(0,iz): 
            for ii in range(0,ix):     
                tmp[kk,:,ii] = varin[kk,:,ii] - var_mean[kk,ii]    
                tmp_std[kk,:,ii] = var_std[kk,ii]    
            
    varanom = tmp
    varanom_std = tmp/tmp_std
    
    return varanom, varanom_std

def rmse(predictions, targets):    
    comp = np.sqrt(((predictions - targets) ** 2))
    return np.mean(np.ma.MaskedArray(comp, np.isnan(comp))) 

def bias(predictions, targets):    
    comp = predictions - targets
    return np.mean(np.ma.MaskedArray(comp, np.isnan(comp))) 

def relative_error(predictions, targets):    
    comp = ((predictions - targets)/targets) * 100
    return np.mean(np.ma.MaskedArray(comp, np.isnan(comp))) 

def find_amplitude(datain,dataininds,thresh,min_or_max):
    """

   Finds the amplitude of datain

   Input:    

       datain - 2D data array 
       dataininds - masked array of indices denoting section of datain array to search for minimum or maximum
       thresh - threshold min/max of datain
       min_or_max - 'min' to find the amplitude from a minimum and 'max' to find the amplitude from a maximum

    Steven Cavallo
    February 2013
    University of Oklahoma
    
    """

   

    import utilities_modules as um
    import numpy as np
    datain = um.filter_numeric_nans(datain,thresh,thresh,'high')
            
    datasave = datain
    
    if min_or_max == 'min':
        mininds=np.where(datain==np.min(datain[dataininds]))   
    else:
        mininds=np.where(datain==np.max(datain[dataininds]))        
    mininds2 = np.ma.getdata(mininds)    
    
    minval = datasave[mininds2[0],mininds2[1]]       
       
    #datain = ndimage.gaussian_filter(datain,0.75)
    datain = ndimage.gaussian_filter(datain,5.0)
    
    [iy,ix] = np.shape(datain)
    minsearch = np.min(mininds2)       
    amps = np.zeros(8)
    for ii in range(0,8):
        for jj in range(minsearch-1,0,-1):
        # Search up
            if ii==0:
                try:
                    if datain[jj-1,mininds2[1]] < datain[jj,mininds2[1]]:           
                        amps[ii] = datasave[jj,mininds2[1]] - minval
                        break
                except:
                    continue
                    
        # Search NE
            if ii==1:
                try:         
                    if datain[jj-1,jj+1] < datain[jj,jj]:           
                        amps[ii] = datasave[jj,jj] - minval
                        break   
                except:
                    continue

        # Search right
            if ii==2:
                try:         
                    if datain[mininds2[0],jj+1] < datain[mininds2[0],jj]:           
                        amps[ii] = datasave[mininds2[0],jj] - minval
                        break                 
                except:
                    continue

       # Search SE
            if ii==3:
                try:         
                    if datain[jj+1,jj+1] < datain[jj,jj]:           
                        amps[ii] = datasave[jj,jj] - minval
                        break                
                except:
                    continue
         
       # Search down
            if ii==4:
                try:         
                    if datain[jj+1,mininds2[1]] < datain[jj,mininds2[1]]:           
                        amps[ii] = datasave[jj,mininds2[1]] - minval
                        break         
                except:
                    continue

       # Search SW
            if ii==5:
                try:         
                    if datain[jj+1,jj-1] < datain[jj,jj]:           
                        amps[ii] = datasave[jj,jj] - minval
                        break             
                except:
                    continue
         
       # Search left
            if ii==6:
                try:         
                    if datain[mininds2[0],jj-1] < datain[mininds2[0],jj]:           
                        amps[ii] = datasave[mininds2[0],jj] - minval
                        break         
                except:
                    continue

       # Search NW
            if ii==7:
                try:                 
                    if datain[jj-1,jj-1] < datain[jj,jj]:           
                        amps[ii] = datasave[jj,jj] - minval
                        break             
                except:
                    continue


    try:
       ampinds = np.where((amps>=1) & (amps<=50))
       #ampout = np.min(amps[ampinds])
       
       ampsort = np.sort(amps[ampinds])
       ampout = np.median(ampsort[0:2])
    except:
       ampout = float('NaN')

    
    #ampout = np.median(amps)
    
    return ampout

def readsounding(fname):
    """

    Reads a sounding text file from the University of Wyoming web page

    Useage example:    

    data,fieldnames = readsounding(wyoming_sounding_textfile)
    print fieldnames

    pres=data["pres"]
    temp=data["temp"]
    dwpt=data["dwpt"]
    hght=data["hght"]
    theta=data["thta"]
    thetae=data["thte"]
    mixr=data["mixr"]
    thetav=data["thtv"]

    Steven Cavallo
    September 2014
    University of Oklahoma
    
    """
    
    data={}
    fieldnames={}

    fid=open(fname)
    lines=fid.readlines()
    nlines=len(lines)
    ndata=nlines-34
    output={}
    fields=lines[3].split()
    units=lines[4].split()
    # First line for WRF profiles differs from the UWYO soundings
    header=lines[0]
    if header[:5]=='00000':
        # WRF profile
        data['StationNumber']='-99999'
        data['Longitude']=float(header.split()[5].strip(","))
        data['Latitude']=float(header.split()[6])
        data['SoundingDate']=header.split()[-1]
    else:
        data['StationNumber']=header[:5]
        dstr=(' ').join(header.split()[-4:])
        data['SoundingDate']=datetime.datetime.strptime(dstr,"%HZ %d %b %Y").strftime("%Y-%m-%d_%H:%M:%S")
    for ff in fields:
        output[ff.lower()]=np.zeros((nlines-34))-999.
    lhi=[1, 9,16,23,30,37,46,53,58,65,72]
    rhi=[7,14,21,28,35,42,49,56,63,70,77]
    lcounter=5
    for line,idx in zip(lines[6:],range(ndata)):
        lcounter+=1
        try: output[fields[0].lower()][idx]=float(line[lhi[0]:rhi[0]])
        except ValueError: break
        for ii in range(1,len(rhi)):
            try:
                # Debug only:
                # print fields[ii].lower(), float(line[lhi[ii]:rhi[ii]].strip())
                output[fields[ii].lower()][idx]=float(line[lhi[ii]:rhi[ii]].strip())
            except ValueError:
                pass
    for field in fields:        
        ff=field.lower()
        data[ff]=np.ma.masked_values(output[ff],-999.)
        fieldnames[ff] = field
    
    return data, fieldnames

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t.ppf((1+confidence)/2., n-1)
    return h, m-h, m+h
    

def bootstrap(data,boot_num):
    data = 1.0*np.array(data)
    n = len(data)
    boot_means = np.empty([boot_num])
    for i in np.arange(0,boot_num,1):
        current_boot = np.random.choice(data,n)
        boot_means[i] = np.mean(current_boot)
    return boot_means

def bootstrap_twosamps(data1,data2,boot_num,alpha):
    import numpy.matlib
    data1vec = np.reshape(data1, np.prod(np.size(data1)), 1)
    data2vec = np.reshape(data2, np.prod(np.size(data2)), 1)
       
    # N(datamean,var):
    data1mean = np.nanmean(data1vec)
    data1std = np.nanstd(data1vec)
    data2mean = np.nanmean(data2vec)
    data2std = np.nanstd(data2vec)
        
    n = len(data1vec)
    boot_means = np.empty([boot_num])
    for i in np.arange(0,boot_num,1):
        boot1 = data1std*np.matlib.randn((n)) + data1mean
	boot2 = data2std*np.matlib.randn((n)) + data2mean
        boot_means[i] = np.nanmean(boot1) - np.nanmean(boot2)
    
    ci_limits = [100.0*alpha/2.0,100.0*(1.0-alpha/2.0)]
    
    ci_upper = np.percentile(boot_means,ci_limits[0])
    ci_lower = np.percentile(boot_means,ci_limits[1])
    if ci_upper > ci_lower:
        ci_limits_out = [ci_lower,ci_upper]    
    else:
        ci_limits_out = [ci_upper,ci_lower]  
    return boot_means,ci_limits_out

    
