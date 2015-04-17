#!/usr/bin/python

import numpy as np
#import matplotlib.pyplot as mpl
import matplotlib as plt
from mpl_toolkits.basemap import Basemap, addcyclic 
import math

def nan2zero(data):
    ''' Convert NaNs to zero '''
    ''' '''
    ''' data: Input data array '''
    dimens = np.shape(data)
               
    # Temporarily collapse data array
    temp = np.reshape(data,np.prod(np.size(data)), 1)       
    
    # Find indices with NaNs
    inds = np.argwhere(np.isnan(temp))    
    
    # Replace NaNs with zero
    temp[inds] = 0.                 
    
    # Turn vector back into array
    data = np.reshape(temp,dimens,order='F').copy()
 
    return data

def zero2nan(data):
    ''' Convert zeros to Nans '''
    ''' '''
    ''' data: Input data array '''
    dimens = np.shape(data)
               
    # Temporarily collapse data array
    temp = np.reshape(data,np.prod(np.size(data)), 1)       
    
    # Find indices with NaNs
    inds = np.argwhere(temp==0)    
    
    # Replace zeros with NaNs
    temp[inds] = float('NaN')                 
    
    # Turn vector back into array
    data = np.reshape(temp,dimens,order='F').copy()
 
    return data

def FixNaNs(arr):
    ''' Convert NaNs to nearest neighbor value '''
    ''' '''
    ''' arr:  Input data array '''
    ''' data: Output data array '''
    dimens = np.shape(arr)
    
    # Temporarily collapse data array
    temp = np.reshape(arr,np.prod(np.size(arr)), 1) 
              
    idxs=np.nonzero(temp==temp)[0]

    if len(idxs)==0:
        return None

    ret=temp

    for i in range(0, idxs[0]):
        ret[i]=ret[idxs[0]]

    for i in range(idxs[-1]+1, ret.size):
        ret[i]=ret[idxs[-1]]

    # Turn vector back into array
    data = np.reshape(ret,dimens,order='F').copy()
    
    return data

def filter_numeric_nans(data,thresh,repl_val,high_or_low) :
    ''' Filter numerical nans above or below a specified value'''
    ''' '''
    ''' data:        (Input) array to filter '''
    ''' thresh:      (Input) threshold value to filter above or below '''
    ''' repl_val:    (Input) replacement value'''
    ''' high_or_low: (Input)''' 


    dimens = np.shape(data)    
    temp = np.reshape(data,np.prod(np.size(data)), 1)    
    if high_or_low=='high':        	
	inds = np.argwhere(temp>thresh) 	
	temp[inds] = repl_val	  
    elif high_or_low=='low':    
        inds = np.argwhere(temp<thresh) 
	temp[inds] = repl_val	  
    elif high_or_low =='both':
       	inds = np.argwhere(temp>thresh) 	
	temp[inds] = repl_val
	del inds
       	inds = np.argwhere(temp<-thresh) 	
	temp[inds] = -repl_val	                 
    else:
        inds = np.argwhere(temp>thresh) 
	temp[inds] = repl_val	  
    
    # Turn vector back into array
    data = np.reshape(temp,dimens,order='F').copy()
    
 
    return data    
    
def bold_labels(ax,fontsize=None):
    if fontsize is None:
        fontsize = 14
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')

#def draw_map_background(m, ax=mpl.gca()):
#    ''' Setup the map background '''
#    m.drawcoastlines(ax=ax, linewidth=2, color='#444444', zorder=6)
#    m.drawcountries(ax=ax, linewidth=1, color='#444444', zorder=5)
#    m.drawstates(ax=ax, linewidth=0.66, color='#444444', zorder=4)
#    m.drawmapboundary
    
def lonswap(d,subtract=0.):
        sh = np.shape(d)	
	midl = sh[1]/2	
	midl = np.round(midl)
	h=d[:,midl:].copy()-subtract
	d[:,midl:]=d[:,:midl].copy()
	d[:,:midl]=h
	return d
	    
def periodic(d,add=0.):
	return np.append( d, (d[:,0].copy()+add).reshape(-1,1) , 1)    
	
def advance_time(timestrin,timeinc):
    ''' Advances or reverses a time by timeinc'''
    ''' '''
    ''' timestrin: (Input) time string in yyyymmddhh format'''
    ''' timeinc:   (Input) number of hours to increment or decrement'''
    '''             Use a negative sign to decrement '''
    
    import datetime
        
    yyyy = timestrin[0:4]
    mm = timestrin[4:6]
    dd = timestrin[6:8]
    hh = timestrin[8:10]	
	
    date=datetime.datetime(int(yyyy),int(mm),int(dd),int(hh))
    date += datetime.timedelta(hours=timeinc)
    tt = date.timetuple()
    
    yyyy = str(tt[0])
    mm = str(tt[1])
    dd = str(tt[2])
    hh = str(tt[3])
    
    if tt[0]<1000: yy = '0'+mm 
    if tt[1]<10: mm = '0'+mm 
    if tt[2]<10: dd = '0'+dd
    if tt[3]<10: hh = '0'+hh
    
    timestrin = yyyy+mm+dd+hh        
    
    return timestrin
def get_cmap_cust():
    ''' Setup a custom colortable. '''
    cdict = {'red': ((0.00, 240/255., 220/255.),
                         (0.25, 40/255., 20/255.),
                         (0.50, 225/255., 255/255.),
                         (0.75, 150/255., 150/255.),
                         (1.00, 255/255., 255/255.)),

             'green': ((0.00, 240/255., 220/255.),
                         (0.25, 0/255., 50/255.),
                         (0.50, 255/255., 255/255.),
                         (0.75, 0/255., 35/255.),
                         (1.00, 225/255., 240/255.)),

             'blue': ((0.00, 255/255., 255/255.),
                         (0.25, 160/255., 150/255.),
                         (0.50, 255/255., 170/255.),
                         (0.75, 0/255., 35/255.),
                         (1.00, 225/255., 240/255.))}
    return plt.colors.LinearSegmentedColormap('cool2warm', cdict, 256)
def cmap_discretize(cmap, N):
   """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet. 
        N: Number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
       imshow(x, cmap=djet)
   """
   from scipy import interpolate
   
   cdict = cmap._segmentdata.copy()
   # N colors
   colors_i = np.linspace(0,1.,N)
   # N+1 indices
   indices = np.linspace(0,1.,N+1)
   for key in ('red','green','blue'):
       # Find the N colors
       D = np.array(cdict[key])
       I = interpolate.interp1d(D[:,0], D[:,1])
       colors = I(colors_i)
       # Place these colors at the correct indices.
       A = np.zeros((N+1,3), float)
       A[:,0] = indices
       A[1:,1] = colors
       A[:-1,2] = colors
       # Create a tuple for the dictionary.
       L = []
       for l in A:
           L.append(tuple(l))
       cdict[key] = tuple(L)
   # Return colormap object.
   return plt.colors.LinearSegmentedColormap('colormap',cdict,1024)

def cmap_whitezero(cmap,N,Nwhites,pos):
   """Whites out middle index of a colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet. 
        N: Number of colors.
	Nwhites: Number of whites
	pos: Position for white bar; if -1, will place in middle
	                             if  0, will place at bottom

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
       imshow(x, cmap=djet)
   """
   from scipy import interpolate
   
   if ( pos == -1 ):
      mid = np.round(N/2)
      mid = int(mid)
   else:
      mid = pos
   
   nadd = Nwhites - 1
   
   
   cdict = cmap._segmentdata.copy()
   # N colors
   colors_i = np.linspace(0,1.,N)

   # N+1 indices
   indices = np.linspace(0,1.,N+1)
   for key in ('red','green','blue'):
       # Find the N colors
       D = np.array(cdict[key])
       I = interpolate.interp1d(D[:,0], D[:,1])
       colors = I(colors_i)
       colors[mid] = 1.
       isodd = 0
       if ( np.mod(N,2) == 0 ) :           
	  colors[mid-1] = 1.
	  isodd = 1
       kk=mid-nadd
       kend=mid+nadd 
       if ( kk < kend ) :      
          while (kk <= kend) :
             colors[kk] = 1.
             kk += 1
       if (isodd == 1 ): colors[kk] = 1.   	  
       # Place these colors at the correct indices.
       A = np.zeros((N+1,3), float)
       A[:,0] = indices
       A[1:,1] = colors
       A[:-1,2] = colors       
       # Create a tuple for the dictionary.
       L = []
       for l in A:
           L.append(tuple(l))
       cdict[key] = tuple(L)
   # Return colormap object.
   return plt.colors.LinearSegmentedColormap('colormap',cdict,1024)

def earth_distm(lat1,lon1,lat2,lon2):
    """

   Calculates the distances between a point and all other points given a latitude and longitude

   Input:    
       lat1, lon1 - Coordinate pair to calculate distance from (must be single values)
       lat2, lon2 - Coordinates to calculate distance to (can be vector or array)
   Output:
       ed - distance between pairs in km
       
    Steven Cavallo
    February 2013
    University of Oklahoma
    
    """
    
    latshape = np.shape(lat2)
    latsize = np.size(lat2)
    ndims = len(latshape)
    
    R_earth = 6371200
    R_earth = R_earth / 1000
    pid = np.pi/180

    if ndims == 1:
         if lon2 < 0:
	    lon2 = lon2 + 360
	 if lon1 < 0:
	    lon1 = lon1 + 360
	    
	 X = lon1
	 Y = lat1	     
    else:

	if latsize > 1:
	   [iy,ix] = np.shape(lat2)

	   X = np.zeros_like(lat2).astype('f')   
	   Y = np.zeros_like(lon2).astype('f')   
	   for ii in range(0,ix):
               for jj in range(0,iy):   
        	   X[jj,ii] = lon1
		   Y[jj,ii] = lat1      
	else:
	   X = lon1
	   Y = lat1   
    	            
    # calculate distance
    ed = R_earth * np.arccos( np.sin(Y*pid) * np.sin(lat2*pid) + np.cos(Y*pid) * np.cos(lat2*pid) * np.cos((lon2 - X)*pid));

    return ed

def destagger_hor_wind(ustag,vstag):
   """
   u,v = destagger_hor_wind(ustag,vstag)

   destaggers horizontal wind

   Steven Cavallo
   March 2013

   """
   
   nshape = np.shape(ustag)
   nel = len(nshape)   
   if nel == 2:
      Ny,Nx = np.shape(ustag)         
      
      u = np.zeros((Ny,Nx-1))
      v = np.zeros((Ny,Nx-1))
      
      for jj in range(0,Ny):
         for ii in range(0,Nx-1):	 	     
             u[jj,ii] = (ustag[jj,ii+1] + ustag[jj,ii])/2


      for jj in range(0,Ny):
         for ii in range(0,Nx-1):			
             v[jj,ii] = (vstag[jj+1,ii] + ustag[jj,ii])/2   
   
   else:
      Nz,Ny,Nx = np.shape(ustag)   
       
      u = np.zeros((Nz,Ny,Nx-1))
      v = np.zeros((Nz,Ny,Nx-1))
      for kk in range(0,Nz):
	 for jj in range(0,Ny):
            for ii in range(0,Nx-1):	 		
        	u[kk,jj,ii] = (ustag[kk,jj,ii+1] + ustag[kk,jj,ii])/2

      for kk in range(0,Nz):
	 for jj in range(0,Ny):
            for ii in range(0,Nx-1):			
        	v[kk,jj,ii] = (vstag[kk,jj+1,ii] + ustag[kk,jj,ii])/2   

   return u,v
def destagger_vertical(varin):
   """
   varout = destagger_hor_wind(varin)

   destaggers in the vertical

   Steven Cavallo
   June 2014

   """
   
   nshape = np.shape(varin)
   nel = len(nshape)   
   if nel == 4:
      Nt,Nz,Ny,Nx = np.shape(varin)         
      
      varout = np.zeros((Nt,Nz-1,Ny,Nx))      
      
      for kk in range(0,Nz-1):         	 	     
          varout[:,kk,:,:] = (varin[:,kk+1,:,:] + varin[:,kk,:,:])/2
   
   else:
      Nz,Ny,Nx = np.shape(varin)   
       
      varout = np.zeros((Nz-1,Ny,Nx))
      for kk in range(0,Nz-1):	             	 		
          varout[kk,:,:] = (varin[kk+1,:,:] + varin[kk,:,:])/2



   return varout 
def grid_to_true_wind(lon,ug,vg,truelat1,truelat2,stdlon,proj_type):
    """

   converts from grid relative wind to true direction wind.  Based on FSL 
   mapping module.

   Input:    
       lon - 2D longitude array on grid wind points
       ug, vg - 2D arrays of grid u and v winds
       truelat1,truelat2 - true latitudes
       stdlon - standard longitude
       proj_type - projection type
   Output:
       ut,vt - output u and v true winds
       
    Steven Cavallo
    March 2013
    University of Oklahoma
    
    """
        
    
    Ny,Nx = np.shape(ug)         
      
    ut = np.zeros((Ny,Nx))
    vt = np.zeros((Ny,Nx))    
    
    pid = np.pi/180
    
    if proj_type == 'lambert':
       cone = (np.log(np.cos(truelat1*pid))-np.log(np.cos(truelat2*pid))) / (np.log(np.tan((90. - abs(truelat1)) * pid * 0.5 )) - \
            np.log(np.tan((90. - abs(truelat2)) * pid * 0.5 )) )           
    if proj_type == 'polar':
       cone = 1
       
    if ( (proj_type == 'polar') or (proj_type == 'lambert') ):
       diff = lon - stdlon
       Ny,Nx = np.shape(diff)       
       for jj in range(0,Ny):
          for ii in range(0,Nx):	 	     
	      diffnow = lon[jj,ii] - stdlon
	      if diffnow > 180:
                 diffnow = diffnow - 360
	      if diffnow < -180:
                 diffnow = diffnow + 360                 
       
	      alpha = diffnow * cone * pid * 1 * np.sign(truelat1); 
	      ut[jj,ii] = vg[jj,ii] * np.sin(alpha) + ug[jj,ii] * np.cos(alpha);
	      vt[jj,ii] = vg[jj,ii] * np.cos(alpha) - ug[jj,ii] * np.sin(alpha);    
    
    else:
       ut = ug
       vt = vg

    return ut,vt

def xsection_inds(slon, slat, elon, elat, lons, lats, m):
    '''
    Returns the indicies for creating a cross section.
    Note: The indices returned are generally south-to-north

    Parameters
    ----------
    slon : scalar
        The longitude value of the starting point
    slat : scalar
        The latitude value of the starting point
    elon : scalar
        The longitude value of the ending point
    elat : scalar
        The latitude value of the ending point
    lons : 2D array_like
        The longitude values of the underlying dataset
    lats : 2D array_like
        The latitude values of the underlying dataset
    m : Basemap Instance
        The basemap instance used for plotting data on a map

    Returns
    -------
    xinds : numpy array
        The first dimension indices of the cross section
    yinds : numpy array
        The second dimension indices of the cross section

    '''
    
    import scipy.spatial as ss
    from PIL import Image, ImageDraw
    
    x, y = m(lons, lats)
    gpoints = zip(x.ravel(), y.ravel())
    gtree = ss.cKDTree(gpoints)
    sx, sy = m(slon, slat)
    ex, ey = m(elon, elat)
    pts = np.array([(sx, sy), (ex, ey)])
    dists, inds = gtree.query(pts, k=1, distance_upper_bound=100*1000.)
    xinds, yinds = np.unravel_index(inds, x.shape)
    pts = ((yinds[0], xinds[0]), (yinds[1], xinds[1]))
    grid = np.zeros_like(lons)
    img = Image.new('L', grid.shape[::-1], 0)
    ImageDraw.Draw(img).line(pts, fill=1, width=1)
    img = np.array(img)
    xinds, yinds = np.where(img > 0)
    if slat > elat:
        xinds = xinds[::-1]
        yinds = yinds[::-1]
    return xinds, yinds

def nbrInds8_ll(iLat, iLon, nLat, nLon):
  #return list of 8-conn nbr indices: [(latNbr1, lonNbr1), (latNbr2, lonNbr2),...]
  #always have east, west neighbors. not north/south at respective pole
  #at the north pole, all pole-1 latitude cells are nbrs.
  
  iWest = (iLon-1)%nLon # -1%4=3 so don't worry about negatives
  iEast = (iLon+1)%nLon
  iSouth = iLat+1
  iNorth = iLat-1
  
  haveSouth = iSouth<nLat
  haveNorth = iNorth>-1
  
  nbrLats = [iLat, iLat]; nbrLons = [iWest, iEast]
  if (haveSouth):
    nbrLats.extend([iSouth, iSouth, iSouth])
    nbrLons.extend([iWest, iLon, iEast])
  else:
    #add on all cells on equatorward latitude. this will duplicate some cells.
    nbrLats.extend([iNorth]*nLon)
    nbrLons.extend(range(nLon))
    
  if (haveNorth):
    nbrLats.extend([iNorth, iNorth, iNorth])
    nbrLons.extend([iWest, iLon, iEast])
  else:
    #add on all cells on equatorward latitude. this will duplicate some cells.
    nbrLats.extend([iSouth]*nLon)
    nbrLons.extend(range(nLon))
    
  #to get flat array indices, call index_2dTo1d(lats, lons, nlon)
  #nbrLats = np.array(nbrLats, dtype=int); nbrLons = np.array(nbrLons, dtype=int)
  return (nbrLats, nbrLons)

def nbrInds4_ll(iLat, iLon, nLat, nLon):
  #return list of 4-conn nbr indices: [(latNbr1, lonNbr1), (latNbr2, lonNbr2),...]
  #always have east, west neighbors. not north/south at respective pole
  #at the north pole, all pole-1 latitude cells are nbrs.
  
  iWest = (iLon-1)%nLon # -1%4=3 so don't worry about negatives
  iEast = (iLon+1)%nLon
  iSouth = iLat+1
  iNorth = iLat-1
  
  haveSouth = iSouth<nLat
  haveNorth = iNorth>-1
  
  nbrLats = [iLat, iLat]; nbrLons = [iWest, iEast]
  if (haveSouth):
    nbrLats.extend([iSouth])
    nbrLons.extend([iLon])
  else:
    #add on all cells on equatorward latitude. this will duplicate some cells.
    nbrLats.extend([iNorth]*nLon)
    nbrLons.extend(range(nLon))
    
  if (haveNorth):
    nbrLats.extend([iNorth])
    nbrLons.extend([iLon])
  else:
    #add on all cells on equatorward latitude. this will duplicate some cells.
    nbrLats.extend([iSouth]*nLon)
    nbrLons.extend(range(nLon))
    
  #to get flat array indices, call index_2dTo1d(lats, lons, nlon)
  #nbrLats = np.array(nbrLats, dtype=int); nbrLons = np.array(nbrLons, dtype=int)
  return (nbrLats, nbrLons)

def nbrInds_ll(iLat, iLon, nLat, nLon):
  #return nbrInds8_ll(iLat, iLon, nLat, nLon)
  return nbrInds4_ll(iLat, iLon, nLat, nLon)

def calc_distSphere_multiple(r, lat1, lon1, lat2, lon2):
  '''
  #return the distance between 1 ll1 point and >=1 ll2 points.
  on a sphere.
  input lat/lon in degrees!
  '''
  
  d2r = np.pi/180.;
  lat1 *= d2r; lon1 *= d2r; lat2 *= d2r; lon2 *= d2r
  
  dlat = lat2-lat1
  dlon = lon2-lon1
  latTerm = np.sin(.5*dlat); latTerm = latTerm*latTerm;
  lonTerm = np.sin(.5*dlon); lonTerm = lonTerm*lonTerm*np.cos(lat1)*np.cos(lat2);
  dAngle = np.sqrt(latTerm+lonTerm)
  dist = 2.*r*np.arcsin(dAngle)
  
  return(dist)

def wrapLongitudeRange(lon):
  #given longitude value, return within range [-180, 180)
  while (lon>180):
    lon = lon-360
  while (lon<-180):
    lon = lon+360
  return lon

def gatherCellsOnLine_coords(lat, lon, latCoord, lonCoord):
  #return a unique array of cells obtained from the coordinate cells
  
  nLat = len(lat); nLon=len(lon); nPts = len(latCoord)
  latCells = np.empty(nPts, dtype=int); lonCells = np.empty(nPts, dtype=int)
  for iPt in xrange(nPts):
    iLat = np.argmin(np.absolute(lat-latCoord[iPt]))
    iLon = np.argmin(np.absolute(lon-lonCoord[iPt]))
    latCells[iPt] = iLat; lonCells[iPt] = iLon
      
  flatInds = latCells*nLon+lonCells
  u, ind = np.unique(flatInds, return_index=True)
  flatInds = u[np.argsort(ind)]

  iLon = ind%nLon
  iLat = (ind-iLon)/nLon  
  return (flatInds, nLon)

def gatherCellsOnLine(iLat0, iLon0, lat1, lon1, lat, lon):
  #the "line" is defined by center cell and latLon of point to reach.
  #return a list of cells crossed to point..
  
  doDiskConn = False
  
  nLat = len(lat); nLon=len(lon)
  
  latInds = [iLat0]; lonInds = [iLon0]
  iLat = iLat0; iLon = iLon0
  dMin = calc_distSphere_multiple(1., lat[iLat0], lon[iLon0], lat1, lon1)
  while (True):
    if (doDiskConn):
      nbrLat, nbrLon = gatherCells_region(iLat, iLon, nLat, nLon, lat, lon, 
                         6371., 30.)
    else:
      nbrLat, nbrLon = nbrInds_ll(iLat, iLon, nLat, nLon)
    d = calc_distSphere_multiple(1., lat1, lon1, lat[nbrLat], lon[nbrLon])
    iNbr = np.argmin(d)
    if (d[iNbr]<dMin):
      iLat = nbrLat[iNbr]; iLon = nbrLon[iNbr]
      latInds.append(iLat); lonInds.append(iLon)
      dMin = d[iNbr]
    else:
      break
  
  return (latInds, lonInds)

def driver_xsec_old(lat, lon, cenLat, cenLon, dLat, dLon):
  #given tpv lat/lon, create a path ll-dll -> ll -> ll+dll
  
  #mesh info ------------------------------
  #lat = data.variables['latitude'][:];
  #lon = data.variables['longitude'][:];
  
  #id cross-section cells ------------------------------
  #if cross pole in latitude, have to wrap to other side of globe
  ll0 = np.array([cenLat-dLat, cenLon-dLon])
  if (ll0[0]<-90.):
    ll0[0] = -90.-(ll0[0]+90.) #need > -np.pi/2
    ll0[1] += 180.
  ll1 = np.array([cenLat, cenLon])
  ll2 = np.array([cenLat+dLat, cenLon+dLon])
  if (ll2[0]>90.):
    ll2[0] = 90.-(ll2[0]-90.) #need < np.pi/2
    ll2[1] += 180.
  #for lat/lon mesh here, longitudes also have to be bounded
  ll0[1] = wrapLongitudeRange(ll0[1])
  ll1[1] = wrapLongitudeRange(ll1[1])
  ll2[1] = wrapLongitudeRange(ll2[1])
  
  #note that the above will fail in purpose if you input a silly dLat,
  #eg, degrees vs. radians
  print ll0, ll1, ll2
  
  #go minus to center
  iLat0 = np.argmin(np.absolute(lat-ll0[0]))
  iLon0 = np.argmin(np.absolute(lon-ll0[1]))
  latCells0, lonCells0 = gatherCellsOnLine(iLat0, iLon0, ll1[0], ll1[1], lat, lon)
  
  #go center to plus
  iLat0 = np.argmin(np.absolute(lat-ll1[0]))
  iLon0 = np.argmin(np.absolute(lon-ll1[1]))
  latCells1, lonCells1 = gatherCellsOnLine(iLat0, iLon0, ll2[0], ll2[1], lat, lon)
  
  #don't double count the tpv cell
  cellsOnLine_lat = np.array(latCells0[0:-1]+latCells1)
  cellsOnLine_lon = np.array(lonCells0[0:-1]+lonCells1)
  
  return (cellsOnLine_lat, cellsOnLine_lon)

def genLatLonCoords(lat0, lon0, dLat, dLon, nPts):
  
  latCoords = np.linspace(lat0-dLat, lat0+dLat, nPts)
  lonCoords = np.linspace(lon0-dLon, lon0+dLon, nPts)

  print latCoords
  print lonCoords
  
  #recover lats that cross pole 1x
  crossedP = latCoords<-np.pi/2
  latCoords[crossedP] = -np.pi/2-(latCoords[crossedP]+np.pi/2)
  lonCoords[crossedP] += np.pi
  
  crossedP = latCoords>np.pi/2
  latCoords[crossedP] = np.pi/2-(latCoords[crossedP]-np.pi/2)
  lonCoords[crossedP] += np.pi
  
  
  lonCoords = lonCoords%(2*np.pi)
  
  return (latCoords, lonCoords)

def driver_xsec(lat, lon, cenLat, cenLon, dLat, dLon):
  
  
  #id cross-section cells ------------------------------
  r2d = 180./np.pi; d2r = np.pi/180.; nPts = 601
  latCoords, lonCoords = genLatLonCoords(cenLat*d2r, cenLon*d2r, dLat*d2r, dLon*d2r, nPts)
  latCoords *= r2d; lonCoords *= r2d
  cellsOnLine_lat, cellsOnLine_lon = gatherCellsOnLine_coords(lat, lon, latCoords, lonCoords)
  
  return (cellsOnLine_lat, cellsOnLine_lon)

def boundary_buffer(field, npts_x, npts_y):
    ''' Replace field with NaNs near boundary defined by npts_x and npts_y '''
    ''' '''
    ''' '''
    ''' field: Input data array '''
    ''' npts_x, npts_y: number of points from boundaries to replace with NaNs '''
    
    dimens = np.shape(field)
    
    if len(dimens)==2:
        iy,ix = np.shape(field)
        field[0:npts_y,:] = float('NaN')
        field[-(npts_y+1),:] = float('NaN')	
        field[:,0:npts_x] = float('NaN')
        field[:,-(npts_x+1)] = float('NaN')		
	
    if len(dimens)==3:
        iz,iy,ix = np.shape(field)    
        field[:,0:npts_y,:] = float('NaN')
        field[:,-(npts_y+1),:] = float('NaN')	
        field[:,:,0:npts_x] = float('NaN')
        field[:,:,-(npts_x+1)] = float('NaN')		    
    
    return field 

def find_matching_indices(a, b):   
    ''' Find the indices of array a that match the values in array b'''
    ''' '''
    ''' a: Input array to match values to '''
    ''' b: Input array with values that are a subset of array a '''
    ''' '''
    ''' Output: minds (indices of a that match b)'''
    ''' '''
    ''' Steven Cavallo '''
    ''' December 2014 '''
    
    ab = np.in1d(a.ravel(), b).reshape(a.shape)
    linds = np.where(ab)

    return linds
def date_to_jd(year,month,day):
    """
    Convert a date to Julian Day.
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet',
    4th ed., Duffet-Smith and Zwart, 2011.
    Parameters
    ----------
    year : int
    Year as integer. Years preceding 1 A.D. should be 0 or negative.
    The year before 1 A.D. is 0, 10 B.C. is year -9.
    month : int
    Month as integer, Jan = 1, Feb. = 2, etc.
    day : float
    Day, may contain fractional part.
    Returns
    -------
    jd : float
    Julian Day
    Examples
    --------
    Convert 6 a.m., February 17, 1985 to Julian Day
    >>> date_to_jd(1985,2,17.25)
    2446113.75
    """
    if month == 1 or month == 2:
	yearp = year - 1
	monthp = month + 12
    else:
	yearp = year
	monthp = month
    # this checks where we are in relation to October 15, 1582, the beginning
    # of the Gregorian calendar.
    if ((year < 1582) or
	(year == 1582 and month < 10) or
	(year == 1582 and month == 10 and day < 15)):
	# before start of Gregorian calendar
	B = 0
    else:
        # after start of Gregorian calendar
        A = math.trunc(yearp / 100.)
        B = 2 - A + math.trunc(A / 4.)
    if yearp < 0:
        C = math.trunc((365.25 * yearp) - 0.75)
    else:
        C = math.trunc(365.25 * yearp)
        D = math.trunc(30.6001 * (monthp + 1))
        jd = B + C + D + day + 1720994.5
    return jd
def jd_to_date(jd):
    """
    Convert Julian Day to date.
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet',
    4th ed., Duffet-Smith and Zwart, 2011.
    Parameters
    ----------
    jd : float
    Julian Day
    Returns
    -------
    year : int
    Year as integer. Years preceding 1 A.D. should be 0 or negative.
    The year before 1 A.D. is 0, 10 B.C. is year -9.
    month : int
    Month as integer, Jan = 1, Feb. = 2, etc.
    day : float
    Day, may contain fractional part.
    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.
    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)
    """
    jd = jd + 0.5
    F, I = math.modf(jd)
    I = int(I)
    A = math.trunc((I - 1867216.25)/36524.25)
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        C = B + 1524
        D = math.trunc((C - 122.1) / 365.25)
        E = math.trunc(365.25 * D)
        G = math.trunc((C - E) / 30.6001)
        day = C - E + F - math.trunc(30.6001 * G)
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
    return year, month, day 
def smooth_onedim(x,npasses):
    """ data = smooth_onedim(x)

     Uses a 5-point moving average with filter on vector x 
     with coefficients equal to the reciprocal of the span 
     to obtain vector dat

     x = input vector
     npasses = number of smoothing passes
     
     Returns vector called data
     
     Steven Cavallo
     January 2015     
    """
    
    if npasses == 0:
       print 'Number of smoothing passes is set to zero.  Set to a value greater than zero to smooth your data.'
       data = x
       return data
    
    for tt in xrange(0,npasses):       
        if tt == 0:
	    data = np.zeros_like(x).astype('f')
	    npts = len(x)

	data[0] = x[1]
	data[1] = (x[0] + x[1] + x[2])/3    
	for ii in xrange(2,npts-2):	   
	   data[ii] = (x[ii-2] + x[ii-1] + x[ii] + x[ii+1] + x[ii+2] ) /5

	data[npts-2] = (x[npts-3] + x[npts-2] + x[npts-1])/3
	data[npts-1] = x[npts-1]

	x = data
    
    return data
    



