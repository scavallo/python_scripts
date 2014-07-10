# imports
import netCDF4

import os, datetime, pylab
import numpy as np
import matplotlib as mpl
import scipy.io
from scipy import interpolate
from scipy import ndimage
#from matplotlib.mlab import griddata
#from griddata import griddata, __version__

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um

from mstats import *

import warnings
warnings.filterwarnings("ignore")


plot_field = 'trth' # options are 'slp', 'trth', 'hght500'
# Set the default domain to be d02
#fpath_composites = '/arctic1/scavallo/wrf_output/sept2012/grid210/'
date_firstrecord = '2011122312'
imagedir = '/data1/scavallo/data/cases/dagmar/wrf_runs/2011122312/images/';
datadir = '/data1/scavallo/data/cases/dagmar/wrf_runs/2011122312/'
fname_netcdf = 'wrfout_d01_2011-12-23_12:00:00'

# General options
titlestring = 'WRF'
timeindices = [0,20]
hinc = 6 # number of hours between output times
label_fontsize = 22
cbar_fontsize = 12
map_projection = 'npstere' # 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection
plot_field = 'sfcwind'
level2plot_pres = 500.
proj_latlon = [75. , 270.]
zoom = 'false'
figname = "dagmar"
###################################

ident = fname_netcdf[7:10]
print ident
if ident == 'd01':
   ftype = 1
   fname_pres = fname_netcdf[0:6] + '_pres_' + fname_netcdf[11::]
   print fname_pres
elif ident == 'pre':
   ftype = 2
elif ident == 'pv':
   ftype = 3

fpath = datadir + fname_netcdf
f = netCDF4.Dataset(fpath,'r')

xlat = f.variables['XLAT'][0,:,:].squeeze()
xlon = f.variables['XLONG'][0,:,:].squeeze()
if ftype == 1:
   levs = f.variables['P'][0,:,100,100].squeeze() + f.variables['PB'][0,:,100,100].squeeze()
else:
   levs =f.variables['PRES'][0,:,100,100].squeeze()
levelindex = np.ravel(levs==level2plot_pres)


x_dim = len(f.dimensions['west_east'])
y_dim = len(f.dimensions['south_north'])
#ntimes = len(f.dimensions['time'])
dx = float(f.DX)
dy = float(f.DY)
cen_lat = float(f.CEN_LAT)
cen_lon = float(f.CEN_LON)
truelat1 = float(f.TRUELAT1)
truelat2 = float(f.TRUELAT2)
standlon = float(f.STAND_LON)
width_meters = dx * (x_dim - 1)
height_meters = dy * (y_dim - 1)
f.close

if plot_field == 'sfcwind':
   cint_field = 2
   min_field = 0
   max_field = 50
   cflevs_field =  np.arange(min_field, max_field, cint_field)
   cflevs_field_ticks = cflevs_field[::5]   
   cmap_now = plt.cm.hot_r
   clabeltext = 'm s-1' 
   
   min_slp = 996-(20*4)
   max_slp = 996+(20*4)
   cint_slp = 4
   cflevs_slp =  np.arange(min_slp, max_slp, cint_slp/2)
   cflevs_slp_ticks = np.arange(min_slp,max_slp,4*(cint_slp/2))  
   

X = xlon
Y = xlat

rcount = 0 
record_end = timeindices[1]
record_numnow = timeindices[0]  
datenow = date_firstrecord
while record_numnow <= record_end:   
   print record_numnow

   if plot_field == 'sfcwind':
      f = netCDF4.Dataset(fpath,'r')
      uin = f.variables['U10'][record_numnow,:,:].squeeze()
      vin = f.variables['V10'][record_numnow,:,:].squeeze()
      plotvar = np.sqrt(uin**2 + vin**2)
      f.close
      
      fpath_pres = datadir + fname_pres
      f2 = netCDF4.Dataset(fpath_pres,'r')
      slp = f2.variables['SLP'][record_numnow,:,:].squeeze()
      f2.close

   
   figname_now = figname + '_' + plot_field + datenow + '.png'
   
   fig = plt.figure(figsize=(12.,12.), dpi=128)   # New figure
   ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
   if map_projection == 'ortho':
      if zoom == 'false':   
         m = Basemap(projection='ortho', lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],
        	  resolution = 'l', area_thresh = 1000.,ax=ax1)
      else:
         m1 = Basemap(projection='ortho', lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],
        	  resolution = 'l', area_thresh = 1000.,ax=ax1)      		           

         width = m1.urcrnrx - m1.llcrnrx
         height = m1.urcrnry - m1.llcrnry

         #coef = 0.5
	 coef = 0.9
         width = width*coef
         height = height*coef
         m = Basemap(projection='ortho',lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],resolution='l',\
             llcrnrx=-0.5*width,llcrnry=-0.5*height,urcrnrx=0.5*width,urcrnry=0.5*height)		  
		  
   elif map_projection == 'lcc':
   #   m = Basemap(llcrnrlon=-125.5,llcrnrlat=15.,urcrnrlon=-30.,urcrnrlat=50.352,\
       m = Basemap(llcrnrlon=-120.0,llcrnrlat=20.,urcrnrlon=-60.0,urcrnrlat=50.0,\
               rsphere=(6378137.00,6356752.3142),\
               resolution='l',area_thresh=1000.,projection='lcc',\
               lat_1=50.,lon_0=-107.,ax=ax1)	       
   elif map_projection == 'npstere':
       if zoom == 'false':
          m = Basemap(projection='npstere',boundinglat=30,lon_0=proj_latlon[1],resolution='l')
       else:
          m = Basemap(projection='npstere',boundinglat=50,lon_0=proj_latlon[1],resolution='l')

   # draw countries, states, and differentiate land from water areas.
   m.drawcoastlines(linewidth=2, color='#444444', zorder=6)
   m.drawcountries(linewidth=1, color='#444444', zorder=5)
   m.drawstates(linewidth=0.66, color='#444444', zorder=4)
   m.drawmapboundary

   # draw lat/lon grid lines every 30 degrees.
   m.drawmeridians(np.arange(0, 360, 30))
   m.drawparallels(np.arange(-90, 90, 30))

   x, y = m(X, Y)

   #trth = um.filter_numeric_nans(trth,340,float('NaN') ,'high')
   
   # Contour tropopause potential temperature
   
   CS1 = m.contourf(x,y,plotvar,cmap=cmap_now,levels=cflevs_field, extend='both',zorder=1)
   cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both',pad=0.05)
   clabs = ['%i K' % f for f in cflevs_field_ticks]
   cbar.ax.set_yticklabels(clabs, size=14) 
   cbar.set_label(clabeltext)
   if plot_field == 'sfcwind':
      CS3 = m.contour(x, y, slp, cflevs_slp, colors='k', linewidths=1.0)        
      plt.clabel(CS3, cflevs_slp, fmt = '%i', inline=True, fontsize=10)   
   
   save_name = imagedir + figname_now
   plt.savefig(save_name, bbox_inches='tight')         
   #plt.show()   
   
   record_numnow +=1
   datenow = um.advance_time(datenow,hinc)      
   

f.close


