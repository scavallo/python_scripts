# # Cross Sections with Python
# 
# -----
# 
# Quick routine to generate cross sections with Python
# 
# Steven Cavallo
# Adapted from Patrick Marsh
# March 2013

###############################################################
# Imports
###############################################################

import os
import numpy as np
import netCDF4 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import ImageGrid
from scipy import ndimage

import utilities_modules as um
import weather_modules as wm
from mstats import *


###############################################################
# User options
###############################################################
fpath = '/data1/scavallo/data/cases/dagmar/wrf_runs/2011122312/wrfout_pres_2011-12-23_12:00:00'
level2plot = 500      # Level to plot (hPa)
timeindex = 8        # Time index to plot
#slat, slon = 55, 30   # starting latitude,longitude pair of cross section.  Use negative longitudes for Western Hemisphere. 
#elat, elon = 55, -150  # ending latitude,longitude pair of cross section.  Use negative longitudes for Western Hemisphere.
slat, slon = 50, -50   # starting latitude,longitude pair of cross section.  Use negative longitudes for Western Hemisphere. 
elat, elon = 70, 30

###############################################################
# END user options.  Be very careful when editing below!
###############################################################

nc = netCDF4.Dataset(fpath, 'r')
hgt = nc.variables['PRES'][0,:,0,0].squeeze()
levelindex = np.ravel(hgt==level2plot)
tmp = nc.variables['THETA'][timeindex,:,:,:].squeeze() # Cross section variable.  Chose temperature as an example here.
#pres = nc.variables['PRES'][timeindex,:,:,:].squeeze() # Cross section variable.  Chose temperature as an example here.
slp = nc.variables['SLP'][timeindex,:,:].squeeze() 
lons = nc.variables['XLONG'][0,:,:].squeeze()
lats = nc.variables['XLAT'][0,:,:].squeeze()

cen_lon = float(nc.CEN_LON)

nc.close

#tmp = wm.temp_to_theta(tmp, pres*100)

# Create a basemap instance. This allows for the cross section to be conducted in the appropriate projection
m = Basemap(projection='npstere',boundinglat=47,lon_0=cen_lon,resolution='l')
# Get the cross section points.  This is a function in utilities_modules.
xx, yy = um.xsection_inds(slon, slat, elon, elat, lons, lats, m)

x, y = m(lons, lats)
#clevs1 = np.arange(220, 260, 1)
#clevs2 = np.arange(190, 273, 2)
clevs1 = np.arange(240, 360, 2)
clevs2 = np.arange(250, 340, 2)

clevs2_ticks = clevs2[::5]
slplevs = np.arange(948,1000,4)
slp_ticks = slplevs[::2]

golden = (np.sqrt(5)+1.)/2.

# Plot plan view map with solid white line indicating the location of the cross section
fig = plt.figure(figsize=(8., 16./golden), dpi=128)   # New figure
ax0 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
F = plt.gcf()  # Gets the current figure


m.drawstates(color='#444444', linewidth=1.25)
m.drawcoastlines(color='#444444')
m.drawcountries(color='#444444', linewidth=1.25)


plotvar = tmp[levelindex,:,:].squeeze()
cbar0 = ax0.contourf(x, y, plotvar, levels=clevs1)
plt.colorbar(cbar0, orientation='horizontal')
cbar1 = ax0.contour(x, y, slp, levels=slplevs,colors='k', linewidths=1.0)
plt.clabel(cbar1,slp_ticks, fmt = '%i', inline=True, fontsize=10)
ax0.plot(x[xx,yy], y[xx,yy], color='w', lw=4)
ax0.plot([x[xx[0],yy[0]], x[xx[-1],yy[-1]]], [y[xx[0],yy[0]], y[xx[-1],yy[-1]]], color='w', lw=2, ls=":")



# Cross section on log-p space
fig = plt.figure(figsize=(12., 12./golden), dpi=128)   # New figure
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.82])

F = plt.gcf()  # Gets the current figure

tx = range(xx.shape[0])
ty = hgt
ttxx, ttyy = np.meshgrid(tx, ty)

plotvar = tmp[:,xx,yy]
plotvar[2:]  = ndimage.gaussian_filter(plotvar[2:],0.75)

CS1 = ax1.contourf(ttxx, ttyy, plotvar,cmap=plt.cm.jet,levels=clevs2)
cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both')
CS2 = ax1.contour(ttxx, ttyy, plotvar,levels=clevs2,colors='k', linewidths=1.0)
#plt.clabel(CS2, clevs2_ticks, fmt = '%i', inline=True, fontsize=10)
cbar.set_label('Kelvin')
labs = []
for xxx,yyy in zip(xx,yy):
    labs.append('%.1f, %.1f' % (lats[xxx,yyy], lons[xxx,yyy]))
ax1.set_xticks(np.arange(0,xx.shape[0]+1e-11, 5))
ax1.set_xticklabels(labs[0::5], rotation=60);
yticks = [1000, 925, 850, 700, 500, 300, 250, 200, 150, 100]
ax1.semilogy()
#ax1.set_ylim(1013, 200)
#ax1.set_ylim(925, 200)
ax1.set_yticks(yticks)
ax1.set_yticklabels(yticks)
#ax1.set_ylim(1013, 200)
ax1.set_ylim(950, 200)

ax0.title.set_y(1.05)
ax0.set_title('Temperature (K)', size=20)
ax1.title.set_y(1.05)
ax1.set_title('Temperature Cross Section (K)', size=20)
ax1.grid(axis='y', ls=':', lw=2)
save_name = 'theta_cross_after_smoothing2.png'
plt.savefig(save_name, bbox_inches='tight')
plt.show()





