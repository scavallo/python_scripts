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
import sounding_modules as sm
from mstats import *


###############################################################
# User options
###############################################################
date_string = '2014110100'

#fpath = '/data1/scavallo/data/cases/dagmar/wrf_runs/2011122312/tpv_removed1/wrfout_d01_pres_2011-12-23_12:00:00'
#fpath_pv = '/data1/scavallo/data/cases/dagmar/wrf_runs/2011122312/tpv_removed1/wrfout_d01_pv_2011-12-23_12:00:00'
fpath = '/arctic2/kelton/arcticwrf/wrfout/pres/wrf_pres_2014-11-09_00:00:00.nc'
fpath_pv = '/arctic2/kelton/arcticwrf/wrfout/pv/wrf_pv_2014-11-09_00:00:00.nc'

#fpath = '/arctic1/scavallo/wrf_runs/dagmar/mpas_physics/wrfout_d02_pres_2011-12-24_12:00:00'
#fpath_pv = '/arctic1/scavallo/wrf_runs/dagmar/mpas_physics/wrfout_d02_pv_2011-12-24_12:00:00'
#fpath = '/arctic1/scavallo/wrf_runs/dagmar/mpas_physics/wrfout_d01_pres_2011-12-23_12:00:00'
#fpath_pv = '/arctic1/scavallo/wrf_runs/dagmar/mpas_physics/wrfout_d01_pv_2011-12-23_12:00:00'
imagedir = '/home/scavallo/scripts/python_scripts/images/'

nest = 'false'
coef = 0.8
figname = "wrf_cross_section_xxx"
figname_location = 'cross_sec_location_xxx' 
level2plot = 500      # Level to plot (hPa)
timeindex = 1       # Time index to plot
pvindex = 1

slat, slon = 30, -86.8  # Nashville
elat, elon = 70, -86.8

#slat, slon = 36.1, -100.0  # Nashville
#elat, elon = 36.1, -60.0

#slat, slon = 30, -84.4  # starting latitude,longitude pair of cross section.  Use negative longitudes for Western Hemisphere. 
#elat, elon = 70, -84.4

plot_option = 5 # 1 for theta 
                # 2 for z
		# 3 for epv
		# 4 for v component of wind
		# 5 for wind magnitude
		# 6 for absolute angular momentum
		# 7 for heat flux
		# 8 for E-P flux
		# 9 for geopotential height
		# 10 for vertical motion (omega)		
		# 11 for PVDIAB
		# 12 for T_TEND_PHYSICS
		# 13 for T_TEND_DYNAMICS
		# 14 for Sutcliffe advection
		# 15 for QCLOUD+QICE
		# 16 for CAPE
center_diagnostic_metric = 'trth' # 'trpr', 'sfcwind', 'slp', 'none'
plot_background = 'trth'		
plot_anomaly = 'false'		
standardize_anomaly = 'false'
plot_logscale = 'true'
plot_isentropes = 'true'
plot_legend = 'false'
pres_range = [50,1000]

# set all of the below if pick_location_from_guess = 'true'
pick_location_from_guess = 'true'

latlon_guess = [58.00,-90.0];
#latlon_guess = [33.7,-84.5]; # Nashville


dist_thresh = 500
#dlat_dlon = [50, 0]
dlat_dlon = [0, 30]

###############################################################
# END user options.  Be very careful when editing below!
###############################################################

nc_pv =  netCDF4.Dataset(fpath_pv, 'r')
trpr = nc_pv.variables['PRES'][timeindex,pvindex,:,:].squeeze()
trprin = nc_pv.variables['PRES'][timeindex,pvindex,:,:].squeeze()
trthin = nc_pv.variables['THETA'][timeindex,pvindex,:,:].squeeze()
nc_pv.close

nc = netCDF4.Dataset(fpath, 'r')
hgt = nc.variables['PRES'][0,:,0,0].squeeze()
levelindex = np.ravel(hgt==level2plot)
tmp = nc.variables['THETA'][timeindex,:,:,:].squeeze() # Cross section variable.  Chose temperature as an example here.
#pres = nc.variables['PRES'][timeindex,:,:,:].squeeze() # Cross section variable.  Chose temperature as an example here.
slp = nc.variables['SLP'][timeindex,:,:].squeeze() 
lons = nc.variables['XLONG'][:,:].squeeze()
lats = nc.variables['XLAT'][:,:].squeeze()
sfcwind = np.sqrt(nc.variables['U10'][timeindex,:,:].squeeze()**2 + nc.variables['V10'][timeindex,:,:].squeeze()**2)
cen_lat = nc.getncattr('CEN_LAT')
cen_lon = float(nc.CEN_LON)
dx = float(nc.DX)
dy = float(nc.DY)




if pick_location_from_guess == 'true':    
    ed = um.earth_distm(latlon_guess[0],latlon_guess[1],lats,lons)     
    latloninds = np.where( (ed < dist_thresh) ) 

    
    #if ( timeindex < 7) or (timeindex > 9):
    if center_diagnostic_metric != 'none':    
	if center_diagnostic_metric == 'slp':    
            a,b = np.where(slp==np.min(slp[latloninds]))	
            print 'Using SLP to find center of cross section with slp minimum ', slp[a,b]
	elif center_diagnostic_metric == 'trpr':
            a,b = np.where(trpr==np.nanmax(trpr[latloninds]))
	    print 'Using trop pressure to find center of cross section with trpr max ', trpr[a,b]
	elif center_diagnostic_metric == 'sfcwind':
            a,b = np.where(sfcwind==np.nanmax(sfcwind[latloninds]))    
	    print 'Using 10-m wind to find center of cross section with wind maximum ', sfcwind[a,b]
	elif center_diagnostic_metric == 'trth':
            a,b = np.where(trthin==np.nanmin(trthin[latloninds]))
	    print 'Using trop pressure to find center of cross section with trpr max ', trthin[a,b]	
	else:
            print 'Invalid entry'
	#a = a + 2
	#b = b + 3
	#if len(a) > 1:
	if 1 == 0:
	   aa = a[0]
	   bb = b[0]
	   del a, b
	   a = aa
	   b = bb
	   del aa, bb


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
	print slat,elat,slon,elon
	

    
    
#tmp = wm.temp_to_theta(tmp, pres*100)

# Create a basemap instance. This allows for the cross section to be conducted in the appropriate projection
proj_latlon = [cen_lat, cen_lon]
if nest == 'true':    
    m1 = Basemap(projection='ortho', lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],
             resolution = 'l', area_thresh = 1000.)      		           

    width = m1.urcrnrx - m1.llcrnrx
    height = m1.urcrnry - m1.llcrnry

    coef = 1.0
    width = width*coef
    height = height*coef
    m = Basemap(projection='ortho',lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],resolution='l',\
        llcrnrx=-0.5*width,llcrnry=-0.5*height,urcrnrx=0.5*width,urcrnry=0.5*height)
else:
    m = Basemap(projection='npstere',boundinglat=30,lon_0=cen_lon,resolution='l')

# Get the cross section points.  This is a function in utilities_modules.
xx, yy = um.xsection_inds(slon, slat, elon, elat, lons, lats, m)

if plot_background == 'trpr':
    cint_bg = 10
    cbar_min_bg = 300
    cbar_max_bg = 850
    cmap_opt_bg = plt.cm.gist_heat_r 
    cflevs_bg = np.arange(cbar_min_bg,cbar_max_bg+(cint_bg/2),cint_bg)
if plot_background == 'trth':
    cint_bg = 1
    cbar_min_bg = 280
    cbar_max_bg = 330
    cmap_opt_bg = plt.cm.jet 
    cflevs_bg = np.arange(cbar_min_bg,cbar_max_bg+(cint_bg/2),cint_bg)    
    
    
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

	figname = "wrf_cross_section_theta"	   	
	if dlat_dlon[0] == 0:
	    figname = "wrf_cross_section_ew_theta"
	if dlat_dlon[1] == 0:
	    figname = "wrf_cross_section_ns_theta"	    

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
	figname = "wrf_cross_section_height"	 	

	if dlat_dlon[0] == 0:
	    figname = "wrf_cross_section_ew_height"
	if dlat_dlon[1] == 0:
	    figname = "wrf_cross_section_ns_height"	    

elif plot_option == 3:
    pvin = nc.variables['PV'][timeindex,:,:,:].squeeze()

    plot_cross = pvin[:,xx,yy].squeeze()
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
        cint = 0.5
        cbar_min = 0
        cbar_max = 15+(cint/2)            
	
	cbar_labels = 'PVU'
	titlestring = "EPV " + date_string
	
	cmap_opt = plt.cm.gist_heat_r

	figname = "wrf_cross_section_epv"	   	
	if dlat_dlon[0] == 0:
	    figname = "wrf_cross_section_ew_epv"
	if dlat_dlon[1] == 0:
	    figname = "wrf_cross_section_ns_epv"	    

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

	figname = "wrf_cross_section_v"	
	if dlat_dlon[0] == 0:
	    figname = "wrf_cross_section_ew_v"
	if dlat_dlon[1] == 0:
	    figname = "wrf_cross_section_ns_v"	    

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

	    cmap_opt = plt.cm.gist_heat
    else:
        cint = 2
        cbar_min = 20
        cbar_max = 70+(cint/2)            
	
	cbar_labels = 'm s-1'
	titlestring = "Wind magnitude " + date_string
	
	cmap_opt = plt.cm.gist_heat_r 

	figname = "wrf_cross_section_windmag"	
	if dlat_dlon[0] == 0:
	    figname = "wrf_cross_section_ew_windmag"
	if dlat_dlon[1] == 0:
	    figname = "wrf_cross_section_ns_windmag"	    
 
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

	figname = "wrf_cross_section_absang"	
	if dlat_dlon[0] == 0:
	    figname = "wrf_cross_section_ew_absang"
	if dlat_dlon[1] == 0:
	    figname = "wrf_cross_section_ns_absang"	    
   
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

   
	figname = "wrf_cross_section_epflux"	
	if dlat_dlon[0] == 0:
	    figname = "wrf_cross_section_ew_epflux"
	if dlat_dlon[1] == 0:
	    figname = "wrf_cross_section_ns_epflux"	    

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

	figname = "wrf_cross_section_ghgt"	
	if dlat_dlon[0] == 0:
	    figname = "wrf_cross_section_ew_ghgt"
	if dlat_dlon[1] == 0:
	    figname = "wrf_cross_section_ns_ghgt"	    
  
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
        cint = 0.005
        cbar_min = -0.05
        cbar_max = 0.05            
	
	cbar_labels = 'Pa s-1'
	titlestring = "Vertical velocity (omega) " + date_string
	
	cmap_opt = plt.cm.RdBu 

        figname = "wrf_cross_section_omega"	    

	if dlat_dlon[0] == 0:
	    figname = "wrf_cross_section_ew_omega"
	if dlat_dlon[1] == 0:
	    figname = "wrf_cross_section_ns_omega"	    

elif plot_option == 11:
    pvdiab = nc.variables['PVDIAB'][timeindex,:,:,:].squeeze()
    plot_cross = pvdiab[:,xx,yy]   
    
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
	
	cbar_labels = 'PVU day-1'
	titlestring = "Diabatic PV tendency " + date_string
	
	cmap_opt = plt.cm.RdBu_r 	

        figname  = "wrf_cross_section_pvdiab"	 
	      
	if dlat_dlon[0] == 0 :
	    figname = "wrf_cross_section_ew_pvdiab"
	if dlat_dlon[1] == 0 :
	    figname = "wrf_cross_section_ns_pvdiab"	    
elif plot_option == 12:
    t_tend_physics = nc.variables['T_TEND_PHYSICS'][timeindex,:,:,:].squeeze()*86400
    plot_cross = t_tend_physics[:,xx,yy]   
    
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
        cint = 0.2
        cbar_min = -10
        cbar_max = 10            
	
	cbar_labels = 'K day-1'
	titlestring = "Total diabatic heating rate " + date_string
	
	cmap_opt = plt.cm.RdBu_r 	

        figname  = "wrf_cross_section_ttendf"	 
	      
	if dlat_dlon[0] == 0 :
	    figname = "wrf_cross_section_ew_ttendf"
	if dlat_dlon[1] == 0 :
	    figname = "wrf_cross_section_ns_ttendf"	     
elif plot_option == 13:
    t_tend_physics = nc.variables['T_TEND_DYNAMICS'][timeindex,:,:,:].squeeze()*86400
    plot_cross = t_tend_physics[:,xx,yy]   
    
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
        cint = 0.25
        cbar_min = -5
        cbar_max = 5            
	
	cbar_labels = 'K day-1'
	titlestring = "Total dynamics heating rate " + date_string
	
	cmap_opt = plt.cm.RdBu_r 	

        figname  = "wrf_cross_section_ttendf"	 
	      
	if dlat_dlon[0] == 0 :
	    figname = "wrf_cross_section_ew_tdynamics"
	if dlat_dlon[1] == 0 :
	    figname = "wrf_cross_section_ns_tdynamics"	     
elif plot_option == 14:
    u = nc.variables['U'][timeindex,:,:,:].squeeze()
    v = nc.variables['V'][timeindex,:,:,:].squeeze()    
    ghgt = nc.variables['HEIGHT'][timeindex,:,:,:].squeeze()
    ghgt = um.nan2zero(ghgt)   
    
    ug, vg = wm.geostrophic_cartesian(u, v, ghgt, lats, dy, dx)
    ug = um.filter_numeric_nans(ug,200,float('NaN'),'high')
    ug = um.filter_numeric_nans(ug,-200,float('NaN'),'low')
    vg = um.filter_numeric_nans(vg,200,float('NaN'),'high')
    vg = um.filter_numeric_nans(vg,-200,float('NaN'),'low')        
    
    thickness = np.zeros_like(ghgt).astype('f')
    thickness[1:-1,:,:] = (ghgt[2:,:,:] - ghgt[:-2,:,:])/2
    thickness[0,:,:] = (ghgt[1,:,:] - ghgt[0,:,:])
    thickness[-1,:,:] = (ghgt[-1,:,:] - ghgt[-2,:,:])         
    
    zeta_a = np.zeros_like(ghgt).astype('f')    
    [iz,iy,ix] = np.shape(zeta_a)
    for ii in range(iz):
        zeta_a[ii,:,:] = wm.vertical_vorticity_cartesian(ug[ii,:,:].squeeze(), vg[ii,:,:].squeeze(), lats, dx, dy, 1)    
    
    thermal_wind_u, thermal_wind_v = wm.thermal_wind_cartesian(thickness, lats, dx, dy)
    thermal_wind_v = um.nan2zero(thermal_wind_v)       
    
    vort_advT = wm.hadvection_cartesian(thermal_wind_u, thermal_wind_v, zeta_a, dx, dy)
    vort_advT = vort_advT * 86400
    #vort_advT = ndimage.gaussian_filter(vort_advT,0.75)        
    
    plot_cross = vort_advT[:,xx,yy]   
    
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
        cbar_min = -10
        cbar_max = 10            
	
	cbar_labels = 's-1 day-1'
	titlestring = "Geostrophic vorticity advection by the thermal wind " + date_string
	
	cmap_opt = plt.cm.RdBu_r 	

        figname  = "wrf_cross_section_sutcliffe"	 
	      
	if dlat_dlon[0] == 0 :
	    figname = "wrf_cross_section_ew_sutcliffe"
	if dlat_dlon[1] == 0 :
            figname = "wrf_cross_section_ns_sutcliffe"	
elif plot_option == 15:
    #qcloud = nc.variables['QCLOUD'][timeindex,:,:,:].squeeze()
    #qice = nc.variables['QICE'][timeindex,:,:,:].squeeze()
    qtot = nc.variables['H_DIABATIC'][timeindex,:,:,:].squeeze() + nc.variables['RTHCUTEN'][timeindex,:,:,:].squeeze() 
    qtot = qtot*3600*24
    plot_cross = qtot[:,xx,yy]   

    
    if plot_anomaly == 'true':
        plot_cross = qtot_cross_anom
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
        cbar_min = -6
        cbar_max = 6            
	
	cbar_labels = 'K day-1'
	titlestring = "Total diabatic heating rate " + date_string
	
	cmap_opt = plt.cm.RdBu_r 	

        figname  = "wrf_cross_section_ttendf"	 
	      
	if dlat_dlon[0] == 0 :
	    figname = "wrf_cross_section_ew_latent"
	if dlat_dlon[1] == 0 :
	    figname = "wrf_cross_section_ns_latent"	     
elif plot_option == 16:
    temp = nc.variables['TEMP'][timeindex,:,:,:].squeeze()
    pres = nc.variables['PRES'][timeindex,:,:,:].squeeze()*100    
    relh = nc.variables['RELH'][timeindex,:,:,:].squeeze()
    hght = nc.variables['HEIGHT'][timeindex,:,:,:].squeeze()
    es = wm.claus_clap(temp)
    qvs = wm.satur_mix_ratio(es, pres)
    qv = (relh/100)*qvs
    dewpt = wm.mixrat_to_td(qv, pres)
    
    temp_cross = temp[:,xx,yy]
    pres_cross = pres[:,xx,yy]
    relh_cross = relh[:,xx,yy]
    hght_cross = hght[:,xx,yy]
    dewpt_cross = dewpt[:,xx,yy]

    [iz,ix] = np.shape(temp_cross)
   
    startp = 100000.    
    capearr = np.zeros_like(temp_cross).astype('f')   
    cinarr = np.zeros_like(temp_cross).astype('f')   
    for kk in range(0,iz-1):            
        print kk   
        for jj in range(0,ix):       
	    startp = pres_cross[kk,jj]
	    startt = temp_cross[kk,jj]
	    startdp = dewpt_cross[kk,jj]
	    if ( np.isnan(startp) or np.isnan(startt) or  np.isnan(startdp)) :
	        cape = float('NaN')  	
	    else:   	
	        if startdp > startt:
		    startdp = startt
                lcl,lfc,el,cape,cin = sm.get_cape(temp_cross[:,jj],pres_cross[:,jj],dewpt_cross[:,jj],hght_cross[:,jj],startp,startt,startdp,totalcape=False)
            #print lcl,lfc,el,cape,cin
	    capearr[kk,jj] = cape
	    try:
	        cinarr[kk,jj] = cin
	    except:
	        cinarr[kk,jj] = float('NaN')   
    
    #mstats(capearr)
    plot_cross = capearr 
    plot_cross2 = cinarr  


    
    if plot_anomaly == 'true':
        plot_cross = qtot_cross_anom
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
        cint = 5
        cbar_min = 0
        cbar_max = 500           

        cint2 = 10
        cbar_min2 = 0
        cbar_max2 = 100         
	cflevs_cin = np.arange(cbar_min2,cbar_max2+(cint2/2),cint2)
	
	cbar_labels = 'J kg-1'
	titlestring = "CAPE " + date_string
	
	cmap_opt = plt.cm.gist_heat_r 	
	cmap_opt2 = plt.cm.cool_r

        figname  = "wrf_cross_section_cape"	 
	      
	if dlat_dlon[0] == 0 :
	    figname = "wrf_cross_section_ew_cape"
	if dlat_dlon[1] == 0 :
	    figname = "wrf_cross_section_ns_cape"	   

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
	figname = "wrf_cross_section_default"	   
	if dlat_dlon[0] == 0:
	    figname = "wrf_cross_section_ew_default"
	if dlat_dlon[1] == 0:
	    figname = "wrf_cross_section_ns_default"	    


if nest == 'true':
   figname = figname + '_d02'
   figname_location = figname_location + '_d02_' + date_string + ".png"
else:
   figname = figname + '_d01'
   figname_location = figname_location + '_d01_' + date_string + ".png"
   
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
slplevs = np.arange(948,1020,4)
slp_ticks = slplevs[::2]

if plot_background == 'trpr':
   plotvar_bg = trprin
if plot_background == 'trth':
   plotvar_bg = trthin
plotvar_bg = um.filter_numeric_nans(plotvar_bg,cbar_max_bg-(cint_bg/2),cbar_max_bg,'high')
plotvar_bg = um.filter_numeric_nans(plotvar_bg,cbar_min_bg,cbar_min_bg,'low')


golden = (np.sqrt(5)+1.)/2.

# Plot plan view map with solid white line indicating the location of the cross section
fig = plt.figure(figsize=(8., 16./golden), dpi=128)   # New figure
ax0 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#F = plt.gcf()  # Gets the current figure
if nest == 'true':
      proj_latlon = [cen_lat, cen_lon]
      m1 = Basemap(projection='ortho', lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],
               resolution = 'l', area_thresh = 1000.,ax=ax0)      		           

      width = m1.urcrnrx - m1.llcrnrx
      height = m1.urcrnry - m1.llcrnry

      #coef = 0.4
      #coef = 0.9
      width = width*coef
      height = height*coef
      m = Basemap(projection='ortho',lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],resolution='l',\
          llcrnrx=-0.5*width,llcrnry=-0.5*height,urcrnrx=0.5*width,urcrnry=0.5*height)
     
      npts_x = 6
      npts_y = 6
      slp = um.boundary_buffer(slp, npts_x, npts_y)
      plotvar_bg = um.boundary_buffer(plotvar_bg, npts_x, npts_y)
else:
      m = Basemap(projection='npstere',boundinglat=30,lon_0=proj_latlon[1],resolution='l')      
m.drawstates(color='#444444', linewidth=1.25)
m.drawcoastlines(color='#444444')
m.drawcountries(color='#444444', linewidth=1.25)


#plotvar = tmp[levelindex,:,:].squeeze()
cbar0 = ax0.contourf(x, y, plotvar_bg, cmap=cmap_opt_bg, levels=cflevs_bg)
plt.colorbar(cbar0,shrink=0.95, orientation='horizontal',extend='both',pad=0.05)
cbar1 = ax0.contour(x, y, slp, levels=slplevs,colors='k', linewidths=1.0)
plt.clabel(cbar1,slp_ticks, fmt = '%i', inline=True, fontsize=10)
ax0.plot(x[xx,yy], y[xx,yy], color='k', lw=4)
ax0.plot([x[xx[0],yy[0]], x[xx[-1],yy[-1]]], [y[xx[0],yy[0]], y[xx[-1],yy[-1]]], color='k', lw=2, ls=":")
plt.savefig(imagedir + figname_location, bbox_inches='tight')

# Cross section on log-p space
fig = plt.figure(figsize=(12., 12./golden), dpi=128)   # New figure
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.82])

F = plt.gcf()  # Gets the current figure

tx = range(xx.shape[0])
ty = hgt

#if ( slon<0 and elon<0):
if ( 1 == 0 ):
   if slon<elon:
       tx = tx[::-1]       
       yy = yy[::-1]
       
       
       
ttxx, ttyy = np.meshgrid(tx, ty)

cflevs_trop = [2.0, 2.0]

plotvar = plot_cross
plotvar[2:]  = ndimage.gaussian_filter(plotvar[2:],0.75)

CS1 = ax1.contourf(ttxx, ttyy, plotvar,cmap=cmap_opt,levels=clevs2)
cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both')
#if plot_option == 16:
#    CS2 = ax1.contour(ttxx, ttyy, plot_cross2,levels=cflevs_cin,colors='r', linestyles='solid',linewidths=1.5)
    
CS2 = plt.contour(ttxx, ttyy,epv_cross,levels=cflevs_trop,colors='k', linestyles='solid',linewidths=2.5)  
if plot_isentropes == 'true':
   cbar_min_theta = 220
   cbar_max_theta = 2100
   cint_theta = 2
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





