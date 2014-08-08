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
date_string = '2011122406'
fpath = '/data1/scavallo/data/cases/dagmar/wrf_runs/2011122312/wrfout_d01_pres_2011-12-23_12:00:00'
fpath_pv = '/data1/scavallo/data/cases/dagmar/wrf_runs/2011122312/wrfout_d01_pv_2011-12-23_12:00:00'
imagedir = '/home/scavallo/chestnut_files/presentations/wwosc_2014/images/cross_sec_wrf/'
figname = "wrf_cross_section_ew_omega"
level2plot = 500      # Level to plot (hPa)
timeindex = 7        # Time index to plot
pvindex = 5
#slat, slon = 55, 30   # starting latitude,longitude pair of cross section.  Use negative longitudes for Western Hemisphere. 
#elat, elon = 55, -150  # ending latitude,longitude pair of cross section.  Use negative longitudes for Western Hemisphere.
slat, slon = 50, -50   # starting latitude,longitude pair of cross section.  Use negative longitudes for Western Hemisphere. 
elat, elon = 70, 30
plot_option = 10 # 1 for theta 
                # 2 for z
		# 3 for epv
		# 4 for v component of wind
		# 5 for wind magnitude
		# 6 for absolute angular momentum
		# 7 for heat flux
		# 8 for E-P flux
		# 9 for geopotential height
		# 10 for vertical motion (omega)		
plot_anomaly = 'false'		
plot_logscale = 'true'
plot_isentropes = 'true'
plot_legend = 'false'
pres_range = [100,1000]

# set all of the below if pick_location_from_guess = 'true'
pick_location_from_guess = 'true'
#latlon_guess = [51.50,-34.50];
#latlon_guess = [54.50,-26.50];
#latlon_guess = [59.00,339.5];
#latlon_guess = [63.00,345.0];

latlon_guess = [58.7,-14.5];
#latlon_guess = [62.5,-10.0];
#latlon_guess = [64.0,1.5];
#latlon_guess = [63.0,14.0];
#latlon_guess = [61.0,23.8];

dist_thresh = 500
dlat_dlon = [0, 50]
#dlat_dlon = [20, 0]

###############################################################
# END user options.  Be very careful when editing below!
###############################################################

nc_pv =  netCDF4.Dataset(fpath_pv, 'r')
trpr = nc_pv.variables['PRES'][timeindex,pvindex,:,:].squeeze()
nc_pv.close

nc = netCDF4.Dataset(fpath, 'r')
hgt = nc.variables['PRES'][0,:,0,0].squeeze()
levelindex = np.ravel(hgt==level2plot)
tmp = nc.variables['THETA'][timeindex,:,:,:].squeeze() # Cross section variable.  Chose temperature as an example here.
#pres = nc.variables['PRES'][timeindex,:,:,:].squeeze() # Cross section variable.  Chose temperature as an example here.
slp = nc.variables['SLP'][timeindex,:,:].squeeze() 
lons = nc.variables['XLONG'][0,:,:].squeeze()
lats = nc.variables['XLAT'][0,:,:].squeeze()

cen_lon = float(nc.CEN_LON)
dx = float(nc.DX)

if pick_location_from_guess == 'true':    
    ed = um.earth_distm(latlon_guess[0],latlon_guess[1],lats,lons)     
    latloninds = np.where( (ed < dist_thresh) ) 

    
    #a,b = np.where(slp==np.min(slp[latloninds]))
    #print slp[a,b]
    a,b = np.where(trpr==np.max(trpr[latloninds]))
    print trpr[a,b]
    
    x_ngpts = np.round(dlat_dlon[1]*(111*np.cos(lats[a,b]*np.pi/180))/(dx/1000))
    y_ngpts = np.round(dlat_dlon[0]*(111*np.cos(lats[a,b]*np.pi/180))/(dx/1000))
    print x_ngpts, y_ngpts
    astart = int(a-y_ngpts)
    aend = int(a+y_ngpts)
    bstart = int(b-x_ngpts)
    bend = int(b+x_ngpts)    
    
    slat = lats[astart,bstart] 
    elat = lats[aend,bend]
    slon = lons[astart,bstart] 
    elon = lons[aend,bend] 
       
    slat = int(np.round(slat))
    elat = int(np.round(elat))
    slon = int(np.round(slon))
    elon = int(np.round(elon))    
    
    
#tmp = wm.temp_to_theta(tmp, pres*100)

# Create a basemap instance. This allows for the cross section to be conducted in the appropriate projection
m = Basemap(projection='npstere',boundinglat=47,lon_0=cen_lon,resolution='l')
# Get the cross section points.  This is a function in utilities_modules.
xx, yy = um.xsection_inds(slon, slat, elon, elat, lons, lats, m)



if plot_option == 1 :
    plot_cross = nc.variables['THETA'][timeindex,:,xx,yy].squeeze()
    
    if plot_anomaly == 'true':
        plot_cross = theta_cross_anom
    
        cint = 5
        cbar_min = -100
        cbar_max = 100+(cint/2)
	
	cbar_labels = 'K' 
	titlestring = "Potential temperature anomaly " + date_string
	
	cmap_opt = plt.cm.RdBu_r		       
        if standardize_anomaly == 'true':            
            cint = 0.2
            cbar_min = -3 
            cbar_max = 3 + (cint/2) 

	    cbar_labels = 'Standard deviations' 
	    titlestring = "Standardized potential temperature anomaly " + date_string

	    cmap_opt = plt.cm.RdBu_r		    
    else:
        cint = 5
        cbar_min = 270
        cbar_max = 2100+(cint/2)            
        
	cbar_labels = 'K' 
	titlestring = "Potential temperature " + date_string
	
	cmap_opt = plt.cm.winter 
    
elif plot_option == 2:
    plot_cross = nc.variables['HEIGHT'][timeindex,:,xx,yy].squeeze()
    
    if plot_anomaly == 'true':
        plot_cross = z_cross_anom
        cint = 5
        cbar_min = -100
        cbar_max = 100+(cint/2)     
	
	cbar_labels = 'meters'
	titlestring = "Geopotential height anomaly " + date_string
	
	cmap_opt = plt.cm.RdBu_r		   
        if standardize_anomaly == 'true':            
            cint = 0.2
            cbar_min = -3 
            cbar_max = 3 + (cint/2)      

	    cbar_labels = 'Standard deviations' 
	    titlestring = "Standardized geopotential height anomaly " + date_string

	    cmap_opt = plt.cm.RdBu_r
    else:
        cint = 5
        cbar_min = 270
        cbar_max = 2100+(cint/2)            
	
	cbar_labels = 'meters'
	titlestring = "Geopotential height " + date_string
	
	cmap_opt = plt.cm.winter 

elif plot_option == 3:
    plot_cross = nc.variables['PV'][timeindex,:,xx,yy].squeeze()
    
    if plot_anomaly == 'true':
        plot_cross = epv_cross_anom
        cint = 0.25
        cbar_min = -3
        cbar_max = 3+(cint/2)     
	
	cbar_labels = 'PVU'
	titlestring = "EPV anomaly " + date_string
	
	cmap_opt = plt.cm.RdBu_r		   
        if standardize_anomaly == 'true':            
            cint = 0.1
            cbar_min = -2 
            cbar_max = 2 + (cint/2)      

	    cbar_labels = 'Standard deviations' 
	    titlestring = "Standardized EPV anomaly " + date_string

	    cmap_opt = plt.cm.RdBu_r
    else:
        cint = 3
        cbar_min = -45
        cbar_max = 45+(cint/2)            
	
	cbar_labels = 'PVU'
	titlestring = "EPV " + date_string
	
	cmap_opt = plt.cm.RdBu_r 	        
elif plot_option == 4:
    
    plot_cross = nc.variables['V'][timeindex,:,xx,yy].squeeze()
    titlestring = "V-component wind anomaly " + date_string

	
    if plot_anomaly == 'true':
        if cross_orient == 'east-west':
            plot_cross = v_cross_anom
        else:
	    plot_cross = u_cross_anom
	cint = 1
        cbar_min = -30
        cbar_max = 30+(cint/2)     
	
	cbar_labels = 'm s-1'
	
	
	cmap_opt = plt.cm.RdBu_r		   
        if standardize_anomaly == 'true':            
            cint = 0.2
            cbar_min = -3 
            cbar_max = 3 + (cint/2)      

	    cbar_labels = 'Standard deviations' 
	    if cross_orient == 'east-west':
	       titlestring = "Standardized V-component anomaly " + date_string
            else:
	       titlestring = "Standardized U-component anomaly " + date_string
	    cmap_opt = plt.cm.RdBu_r
    else:
        cint = 3
        cbar_min = -51
        cbar_max = 51+(cint/2)            
	
	cbar_labels = 'm s-1'
	#titlestring = "Meridional wind " + date_string
	#date_forecast[6:8] + ' January ' + date_forecast[0:4]
	titlestring = "Meridional wind " + date_string[8:10] + ' UTC ' + date_string[6:8] + ' January ' + date_string[0:4]
	
	cmap_opt = plt.cm.RdBu_r 
	        
elif plot_option == 5:
    vin = nc.variables['V'][timeindex,:,:,:].squeeze()
    uin = nc.variables['U'][timeindex,:,:,:].squeeze()
    windmag = np.sqrt(uin**2 + vin**2)
    
    plot_cross = windmag[:,xx,yy].squeeze()
    
    if plot_anomaly == 'true':        
        cint = 5
        cbar_min = -100
        cbar_max = 100+(cint/2)     
	
	cbar_labels = 'm s-1'
	titlestring = "Wind magnitude anomaly " + date_string
	
	cmap_opt = plt.cm.RdBu_r		   
        if standardize_anomaly == 'true':            
            cint = 0.2
            cbar_min = -3 
            cbar_max = 3 + (cint/2)      

	    cbar_labels = 'Standard deviations' 
	    titlestring = "Standardized V-component anomaly " + date_string

	    cmap_opt = plt.cm.RdBu_r
    else:
        cint = 4
        cbar_min = 0
        cbar_max = 80+(cint/2)            
	
	cbar_labels = 'meters'
	titlestring = "Wind magnitude " + date_string
	
	cmap_opt = plt.cm.Reds 
	
elif plot_option == 6:
    plot_cross = angular_mom
    
    if plot_anomaly == 'true':
        cint = 100
        cbar_min = 0
        cbar_max = 50000+(cint/2) 
	
	cbar_labels = 'm2 s-1'
	titlestring = "Absolute angular momentum anomaly " + date_string
	
	cmap_opt = plt.cm.winter		            
        if standardize_anomaly == 'true':
            cint = 0.2
            cbar_min = -3 
            cbar_max = 3 + (cint/2)

	    cbar_labels = 'Standard deviations' 
	    titlestring = "Standardized angular momentum anomaly " + date_string

	    cmap_opt = plt.cm.RdBu_r		     
    else:
        cint = 25000
        cbar_min = 0
        cbar_max = 500000+(cint/2)            
        
	cbar_labels = 'm2 s-1'
	titlestring = "Absolute angular momentum " + date_string
	
	cmap_opt = plt.cm.Reds 

elif plot_option == 8:
    plot_cross = u_cross # this is actually the E-P flux
    #plot_cross = ndimage.gaussian_filter(plot_cross,2.0)    
    
    if plot_anomaly == 'true':
        plot_cross = plot_cross*10**5
        plot_cross = u_cross_anom
        cint = 0.2
        cbar_min = -5
        cbar_max = 5+(cint/2) 
	
	cbar_labels = 'm2 s-1'
	titlestring = "Eliassen-Palm flux divergence anomaly " + date_string
	
	cmap_opt = plt.cm.RdBu		            
        if standardize_anomaly == 'true':
            cint = 0.1
            cbar_min = -2
            cbar_max = 2 + (cint/2)

	    cbar_labels = 'Standard deviations' 
	    titlestring = "Standardized Eliassen-Palm flux divergence " + date_string

	    cmap_opt = plt.cm.RdBu		    
    else:
        cint = -0.01
        cbar_min = -0.1
        cbar_max = 0.1+(cint/2)            
        
	cbar_labels = ''
	titlestring = "Eliassen-Palm flux divergence " + date_string
	
	cmap_opt = plt.cm.RdBu 

    mstats(plot_cross)

elif plot_option == 9:
    plot_cross = geop_cross / 9.81
    
    if plot_anomaly == 'true':
        plot_cross = geop_cross_anom / 9.81
        cint = 50
        cbar_min = -500
        cbar_max = 500+(cint/2)     
	
	cbar_labels = 'meters'
	titlestring = "Geopotential height anomaly " + date_string
	
	cmap_opt = plt.cm.RdBu_r		   
        if standardize_anomaly == 'true':            
            cint = 0.1
            cbar_min = -2 
            cbar_max = 2 + (cint/2)      

	    cbar_labels = 'Standard deviations' 
	    titlestring = "Standardized Geopotential height anomaly " + date_string

	    cmap_opt = plt.cm.RdBu_r
    else:
        cint = 3
        cbar_min = -45
        cbar_max = 45+(cint/2)            
	
	cbar_labels = 'meters'
	titlestring = "Geopotential height " + date_string
	
	cmap_opt = plt.cm.RdBu_r 	


elif plot_option == 10:
    win = nc.variables['W'][timeindex,:,:,:].squeeze()
    pin = nc.variables['PRES'][timeindex,:,:,:].squeeze()*100
    tempk = nc.variables['TEMP'][timeindex,:,:,:].squeeze()
    omeg = wm.w_to_omega(win,pin,tempk)
    plot_cross = omeg[:,xx,yy]   
    
    if plot_anomaly == 'true':
        plot_cross = w_cross_anom
        cint = 0.01
        cbar_min = -0.1
        cbar_max = 0.1+(cint/2)     
	
	cbar_labels = 'Pa s-1'
	titlestring = "Vertical velocity anomaly " + date_string
	
	cmap_opt = plt.cm.RdBu		   
        if standardize_anomaly == 'true':            
            cint = 0.1
            cbar_min = -2 
            cbar_max = 2 + (cint/2)      

	    cbar_labels = 'Standard deviations' 
	    titlestring = "Standardized vertical velocity anomaly " + date_string

	    cmap_opt = plt.cm.RdBu_r
    else:
        cint = 0.1
        cbar_min = -1
        cbar_max = 1            
	
	cbar_labels = 'Pa s-1'
	titlestring = "Vertical velocity (omega) " + date_string
	
	cmap_opt = plt.cm.RdBu 	

else:
    plot_cross = vt_cross
    
    if plot_anomaly == 'true':
        cint = 50
        cbar_min = -500
        cbar_max = 500+(cint/2) 
	
	cbar_labels = 'm K s-1'
	titlestring = "Heat flux " + date_string
	
	cmap_opt = plt.cm.RdBu_r		            
        if standardize_anomaly == 'true':
            cint = 0.2
            cbar_min = -3 
            cbar_max = 3 + (cint/2)

	    cbar_labels = 'Standard deviations' 
	    titlestring = "Standardized heat flux anomaly " + date_string

	    cmap_opt = plt.cm.RdBu		     
    else:
        cint = 30
        cbar_min = -300
        cbar_max = 300+(cint/2)            
        
	cbar_labels = 'm K s-1'
	titlestring = "Heat flux " + date_string
	
	cmap_opt = plt.cm.RdBu 

epvin = nc.variables['PV'][timeindex,:,:,:].squeeze()
epv_cross = epvin[:,xx,yy]       

thetain = nc.variables['THETA'][timeindex,:,:,:].squeeze()
theta_cross = thetain[:,xx,yy]
nc.close    

plot_cross = um.filter_numeric_nans(plot_cross,cbar_max-(cint/2),cbar_max,'high')
plot_cross = um.filter_numeric_nans(plot_cross,cbar_min,cbar_min,'low')

x, y = m(lons, lats)
clevs1 = np.arange(240, 360, 2)
#clevs2 = np.arange(250, 340, 2)
clevs2 = np.arange(cbar_min,cbar_max+(cint/2),cint)

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

cflevs_trop = [2.0, 2.0]

plotvar = plot_cross
plotvar[2:]  = ndimage.gaussian_filter(plotvar[2:],0.75)

CS1 = ax1.contourf(ttxx, ttyy, plotvar,cmap=cmap_opt,levels=clevs2)
cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both')
CS2 = plt.contour(ttxx, ttyy,epv_cross,levels=cflevs_trop,colors='k', linestyles='solid',linewidths=2.5)  
if plot_isentropes == 'true':
   cbar_min_theta = 220
   cbar_max_theta = 2100
   cint_theta = 5
   cflevs_theta =  np.arange(cbar_min_theta, cbar_max_theta, cint_theta)

   CS3 = plt.contour(ttxx,ttyy,theta_cross,levels=cflevs_theta,colors='b', linestyles='solid',linewidths=1.5)  
   if plot_legend == 'true':
      legend = ax1.legend(loc='center right', shadow=True)
      legend.set_zorder(20)

#CS2 = ax1.contour(ttxx, ttyy, plotvar,levels=clevs2,colors='k', linewidths=1.0)
#plt.clabel(CS2, clevs2_ticks, fmt = '%i', inline=True, fontsize=10)
cbar.set_label(cbar_labels)
labs = []
for xxx,yyy in zip(xx,yy):
    labs.append('%.1f, %.1f' % (lats[xxx,yyy], lons[xxx,yyy]))
ax1.set_xticks(np.arange(0,xx.shape[0]+1e-11, 5))
ax1.set_xticklabels(labs[0::5], rotation=60);

plt.ylim((pres_range[0],pres_range[1]))
ax1.set_ylim(ax1.get_ylim()[::-1])
ymin = np.min(pres_range)
ymax = np.max(pres_range)
print 'ymin is ', ymin
if ymin == 1:
   levs_label = [1, 5, 10, 20, 30, 50, 70, 100, 250, 500, 700, 850, 1000]
elif ymin == 5:
   levs_label = [5, 10, 20, 30, 50, 70, 100, 250, 500, 700, 850, 1000]   
elif ymin == 10:
   levs_label = [10, 20, 30, 50, 70, 100, 250, 500, 700, 850, 1000]   
elif ymin == 50:
   levs_label = [50, 70, 100, 250, 500, 700, 850, 1000]   
elif ymin == 100:
    levs_label = [100, 250, 500, 700, 850, 1000]  
else:
    levs_label = [100, 250, 500, 700, 850, 1000]  



if plot_logscale == 'true':
   ax1.set_yscale('log')   
   ycint = 200
   #yticks = np.arange(ymin,ymax+(ycint/10),ycint)

   ax1.set_yticks(levs_label)
   ax1.set_yticklabels(levs_label, size=10)


ax0.title.set_y(1.05)
#ax0.set_title(titlestring, size=20)
ax1.title.set_y(1.05)
#ax1.set_title(titlestring, size=20)
#ax1.grid(axis='y', ls=':', lw=2)

save_name = figname + "_" + date_string + ".png"
plt.savefig(imagedir + save_name, bbox_inches='tight')

plt.show()





