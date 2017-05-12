#!/usr/bin/python 


# imports
import netCDF4

import os, datetime, pylab
import numpy as np
import matplotlib as mpl
from scipy import ndimage
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy import ndimage
from matplotlib.colors import LogNorm

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um
from mstats import *


import warnings
warnings.filterwarnings("ignore")
###################################
# Set user options
###################################
date_firstrecord = '2010050100'

map_projection = 'npstere' # 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection 
#proj_latlon = [60. , 0.]
proj_latlon = [55. , 270.]
plot_field = 'trth' # 'trth', 'trpr', 'trtemp','pw', 'trwind', 'sfcwind', 'totqv', 'quivwind','ghgt500', 'blank' for blank map
figname = 'erainterim_reanalysis_trth'
plot_trth_contours = 'True'
plot_slp = 'False' # If true, plots slp contours
smooth_plotvar = 'False'
zoom = 'false'
coef = 0.8

record_num = [0, 3] # time index to plot.  0 means the first time, 1 the second time...etc.

hinc = 6 # number of hours between records;  only used if record_num > -1

data_gridnum = 12 # 3 = ncep
                  # 4 = ncep
		  # 11 = era interim
		  # 12 = era interim, alternative (2 PVU surface only)


label_fontsize = 16

# Now provide the path to the directory containing the .nc file. Please note,
# do NOT include the .nc file in the path.
fdir = '/data2/scavallo/era_interim/'
fname_trop = 'era_testdata_20100500_trth.nc'
fname_sfc = 'era_testdata_20100500_slp.nc'
imagedir = '/home/scavallo/scripts/python_scripts/images/'

###################################
# END user options
###################################
nplots = record_num[1] - record_num[0]
record_start = record_num[0]
record_end = record_num[1]

rcount = 0 
record_numnow = record_num[0] 
count_anal = record_numnow
datenow = date_firstrecord
print record_start, record_end, nplots
while record_start <= record_end:   
   #record_numnow = record_start
   print "rcount is %i while record_numnow is %i" %(rcount, record_numnow)

   date_string = datenow
   print date_string
   dt = datetime.datetime.strptime(date_string, '%Y%m%d%H')

   print fname_trop
   fpath = fdir + fname_trop

   print fname_sfc
   fpath_sfc = fdir + fname_sfc


   print fpath
   f = netCDF4.Dataset(fpath,'r')
   if data_gridnum == 4:
      trtempin = f.variables['TMP_P0_L109_GLL0'][record_numnow,1,::-1,:].squeeze()
      trpresin = f.variables['PRES_P0_L109_GLL0'][record_numnow,1,::-1,:].squeeze()
      #uin = f.variables['UGRD_P0_L109_GLL0'][record_numnow,1,::-1,:].squeeze()
      #vin = f.variables['VGRD_P0_L109_GLL0'][record_numnow,1,::-1,:].squeeze()
      uin = f.variables['UGRD_P0_L100_GLL0'][record_numnow,19,::-1,:].squeeze()
      vin = f.variables['VGRD_P0_L100_GLL0'][record_numnow,19,::-1,:].squeeze()
      slp = f.variables['PRMSL_P0_L101_GLL0'][record_numnow,::-1,:].squeeze()/10**2
      
      lons = f.variables['lon_0'][:]
      lats = f.variables['lat_0'][::-1] # Read in reverse direction
      levs = f.variables['lv_ISBL0'][:]  
      print levs    
    
   elif data_gridnum == 11:
      f2 = netCDF4.Dataset(fpath_sfc,'r') 
      
      
      trth = f.variables['th'][record_numnow,3,::-1,:].squeeze()      
      uin = f.variables['u'][record_numnow,3,::-1,:].squeeze()
      vin = f.variables['v'][record_numnow,3,::-1,:].squeeze()
      if plot_field == 'sfcwind':         
         uin_sfc = f2.variables['u10'][record_numnow,:,:].squeeze()
         vin_sfc = f2.variables['v10'][record_numnow,:,:].squeeze()          	 
      if plot_field == 'ghgt500':         
         levels =  f2.variables['level'][:]
	 [ind500] = np.where(levels==500)
         ghgt500 = (f2.variables['z'][record_numnow,ind500,:,:].squeeze())/9.81
	 
             
      slp = (f2.variables['msl'][record_numnow,:,:].squeeze())/100.     
      trpresin = f.variables['p'][record_numnow,3,::-1,:].squeeze()
      
      xlons = f.variables['xlong'][::-1,:]
      lons = f.variables['xlong'][0,:]
      lats = f.variables['xlat'][::-1,0] 
      xlats = f.variables['xlat'][::-1,:] 
      levs = f.variables['levels'][:]             
      
      inds = np.where(trth>400)
      uin[inds] = float('NaN');
      vin[inds] = float('NaN');  
      
      if plot_field == 'pw':
          spechum = (f2.variables['q'][record_numnow,:,:,:].squeeze()) 
	  tempin = (f2.variables['t'][record_numnow,:,:,:].squeeze()) 
	  hght = (f2.variables['z'][record_numnow,:,:,:].squeeze())/9.81 
	  pin = f2.variables['level'][:]*100         
	  qv = (spechum/(1-spechum))	

	  pres = np.zeros_like(tempin).astype('f')   
          for kk in range(0,len(levs)):      
              pres[kk,:,:] = pin[kk]     
          
          qvtot = wm.total_col(qv[::-1,:,:], pres[::-1,:,:], tempin[::-1,:,:], hght[::-1,:,:])         
      if plot_field == 'totqv':
          qvtot = (f2.variables['q'][record_numnow,::-1,:].squeeze()) 
      
      f2.close
   elif data_gridnum == 12:
      if plot_slp == 'True':
          f2 = netCDF4.Dataset(fpath_sfc,'r') 
      
      
      trth = f.variables['pt'][record_numnow,::-1,:].squeeze() 
      trpr = (f.variables['pres'][record_numnow,::-1,:].squeeze()) /100.

      uin = f.variables['u'][record_numnow,::-1,:].squeeze()
      vin = f.variables['v'][record_numnow,::-1,:].squeeze()
      if plot_field == 'sfcwind':         
         uin_sfc = f2.variables['u10'][record_numnow,:,:].squeeze()
         vin_sfc = f2.variables['v10'][record_numnow,:,:].squeeze()          	 
      if plot_slp == 'True':
          slp = (f2.variables['msl'][record_numnow,::-1,:].squeeze())/100.     
      else:
          slp = trth
      trpresin = f.variables['pres'][record_numnow,::-1,:].squeeze()
      trte = wm.theta_to_temp(trth,trpresin)

      
      #xlons = f.variables['longitude'][::-1,:]
      lons = f.variables['longitude'][:]
      lats = f.variables['latitude'][::-1] 
      #xlats = f.variables['latitude'][::-1,:] 
      xlons, xlats = np.meshgrid(lons, lats)   
      levs = [2.0]             
      
      inds = np.where(trth>400)
      uin[inds] = float('NaN');
      vin[inds] = float('NaN');  
      
      if plot_field == 'pw':
          spechum = (f2.variables['q'][record_numnow,:,:,:].squeeze()) 
	  tempin = (f2.variables['t'][record_numnow,:,:,:].squeeze()) 
	  hght = (f2.variables['z'][record_numnow,:,:,:].squeeze())/9.81 
	  pin = f2.variables['level'][:]*100         
	  qv = (spechum/(1-spechum))	

	  pres = np.zeros_like(tempin).astype('f')   
          for kk in range(0,len(levs)):      
              pres[kk,:,:] = pin[kk]     
          
          qvtot = wm.total_col(qv[::-1,:,:], pres[::-1,:,:], tempin[::-1,:,:], hght[::-1,:,:])         
      if plot_field == 'totqv':
          qvtot = (f2.variables['tcwv'][record_numnow,::-1,:].squeeze()) 
      
      if plot_slp == 'True':
          f2.close
   elif data_gridnum == -1:
      trth = f.variables['TROP_THETA'][record_numnow,::-1,:].squeeze()
      slp = (f.variables['SLP'][record_numnow,::-1,:].squeeze())/100.     
      trpresin = f.variables['TROP_PRESSURE'][record_numnow,::-1,:].squeeze()
      
      lons = f.variables['lon'][record_numnow,:].squeeze()
      lats = f.variables['lat'][record_numnow,::-1].squeeze() # Read in reverse direction
      levs = []         
   else:
      trtempin = f.variables['TMP_3_PVL_10'][record_numnow,0,::-1,:].squeeze()
      trpresin = f.variables['PRES_3_PVL_10'][record_numnow,0,::-1,:].squeeze()
      try:
         slp = f.variables['MSLET_3_MSL_10'][record_numnow,::-1,:].squeeze()/10**2   
      except:
         slp = f.variables['PRMSL_3_MSL_10'][record_numnow,::-1,:].squeeze()/10**2   
      lons = f.variables['lon_3'][:]
      lats = f.variables['lat_3'][::-1] # Read in reverse direction
      levs = f.variables['lv_ISBL3'][:]

   f.close
   

   if data_gridnum > -1:
      if ( (data_gridnum != 11) & (data_gridnum != 12) ):
         trth = wm.temp_to_theta(trtempin, trpresin)

   lonin = lons
   trth, lons = um.addcyclic(trth, lonin)
   if plot_field == 'trtemp':
       del trth
       trth, dummy = um.addcyclic(trte, lonin)
   slp, dummy = um.addcyclic(slp, lonin)

   
   if data_gridnum != 11 :
      X, Y = np.meshgrid(lons, lats)      
   else:
      xlons, dummy = um.addcyclic(xlons, lonin)
      xlats, dummy = um.addcyclic(xlats, lonin)
      X = xlons
      Y = xlats

   if map_projection == 'ortho' :
      trth_thresh = 380
      cbar_min_trth = 270   
   
   else:
      trth_thresh = 500
      #cbar_min_trth = 305
      cbar_min_trth = 270
   
   if plot_field == 'trth' :

       cbarlabel = 'Kelvin'
   
       #cbar_max_trth = trth_thresh
       #cbar_max_trth = 335
       cbar_max_trth = 380
       #cbar_max_trth = 400
       
       
       cint_trth = 1
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth, cint_trth)
       #cflevs_trth_cntrs = cflevs_trth[0:12]
       cflevs_trth_cntrs = cflevs_trth
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.jet
       #cmap_opt = plt.cm.Blues_r

       trth = ndimage.gaussian_filter(trth,0.75)
       
       trth = um.filter_numeric_nans(trth,trth_thresh+cint_trth,trth_thresh+cint_trth,'high')
         
       
       plotvar = trth
       mstats(plotvar)
          	    
       #plotvar, lons = um.addcyclic(trth[tt,:,:].squeeze(), lonin)	    	    	    	    
       #plotvar = ndimage.gaussian_filter(plotvar,0.75)   
       
       titletext1 = 'Tropopause potential temperature gradient %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))        
   if plot_field == 'trpr' :
 
       cbarlabel = 'Pressure (hPa)'    
       
       cint_trth = 20
       cbar_max_trth = 800.
       cbar_min_trth = 300.
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth, cint_trth)
       #cflevs_trth_cntrs = cflevs_trth[0:12]
       cflevs_trth_cntrs = cflevs_trth
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.gist_heat_r
       #cmap_opt = plt.cm.Blues_r

       trpr = ndimage.gaussian_filter(trpr,0.75)
       
 


       trpr = um.filter_numeric_nans(trpr,cbar_max_trth+cint_trth,cbar_max_trth+cint_trth,'high')
       #trth = um.filter_numeric_nans(trth,cbar_max_trth+cint_trth,float('NaN'),'high')       
       trpr, dummy = um.addcyclic(trpr, lonin)
       plotvar = trpr
       mstats(plotvar)
          	    
       #plotvar, lons = um.addcyclic(trth[tt,:,:].squeeze(), lonin)	    	    	    	    
       #plotvar = ndimage.gaussian_filter(plotvar,0.75)   
       
       titletext1 = 'Tropopause pressure %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))       
   if plot_field == 'trtemp' :
       trth_thresh = 220
       cbar_min_trth = 220       
       if convert_to_celsius == 'True':
           cbarlabel = r'$^{\circ}$C'
       else:
           #cbarlabel = 'Potential temperature (K)'
           cbarlabel = 'Kelvin'
   
       #cbar_max_trth = trth_thresh
       #cbar_max_trth = 335
       cbar_max_trth = 250
       #cbar_max_trth = 400
       
       if convert_to_celsius == 'True':
           cbar_max_trth = np.ceil(cbar_max_trth - 273.15)
	   cbar_min_trth = np.floor(cbar_min_trth - 273.15)
	   trth_thresh = np.ceil(trth_thresh - 273.15)
       
       cint_trth = 2
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth, cint_trth)
       cflevs_trth_cntrs = cflevs_trth[-10:]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       #cmap_opt = plt.cm.jet
       cmap_opt = plt.cm.Blues


       
       if convert_to_celsius == 'True':
           trth[:] = trth[:] - 273.15    


       trth = um.filter_numeric_nans(trth,trth_thresh-cint_trth,float('NaN'),'low')
       #trth = um.filter_numeric_nans(trth,cbar_max_trth+cint_trth,float('NaN'),'high')       

       plotvar = trth
       mstats(plotvar)

          	    
       #plotvar, lons = um.addcyclic(trth[tt,:,:].squeeze(), lonin)	    	    	    	    
       #plotvar = ndimage.gaussian_filter(plotvar,0.75)   
       
       titletext1 = 'Tropopause temperature gradient %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))        
   if ( (plot_field == 'pw') or (plot_field == 'totqv') ):
       qvtot, dummy = um.addcyclic(qvtot, lonin)   
       plotvar = qvtot
       cbarlabel = 'Precipitable water (mm)'
   
       #cbar_max_trth = 30
       cbar_max_trth = 18
       cbar_min_trth = 0
       cint_trth = 0.5
       
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth + (cint_trth/2), cint_trth)
       #cflevs_trth = np.logspace(-0.5, 1.5, num=50, endpoint=True, base=10.0)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
      
       #cflevs_trth = np.logspace(0.4, 1.2, num=100, endpoint=True, base=10.0) # Range is from 10^-1 to 10^9 meters with 100 points
       #cflevs_trth_cntrs = cflevs_trth[::2]
       #cflevs_trth_ticks = cflevs_trth[::2]
       
       #cmap_opt = plt.cm.BuGn 
       cmap_opt = plt.cm.Greens
       
       titletext1 = 'Precipitable water vapor %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00')) 
   if plot_field == 'trwind' :
       trwind = np.sqrt(uin**2 + vin**2)
       trwind, dummy = um.addcyclic(trwind, lonin)          
       plotvar = trwind
       cbarlabel = 'Tropopause wind (m s-1)'
   
       cbar_max_trth = 100
       cbar_min_trth = 20
       cint_trth = 2
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth + (cint_trth/2), cint_trth)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.hot_r
       
       titletext1 = 'Tropopause wind %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))      
   if plot_field == 'quivwind' :
       trwind = np.sqrt(uin**2 + vin**2)       
       trwind, dummy = um.addcyclic(trwind, lonin)   
       uwind, dummy = um.addcyclic(uin, lonin) 
       vwind, dummy = um.addcyclic(vin, lonin)   
       plotvar = trwind
       cbarlabel = 'Tropopause wind (m s-1)'
   
       cbar_max_trth = 50
       cbar_min_trth = 10
       cint_trth = 1
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth + (cint_trth/2), cint_trth)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       #cmap_opt = plt.cm.RdBu_r
       cmap_opt = plt.cm.gist_heat_r
       
       titletext1 = 'Tropopause wind %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00'))      

   if plot_field == 'sfcwind' :
       trwind = np.sqrt(uin**2 + vin**2)
       sfcwind = np.sqrt(uin_sfc**2 + vin_sfc**2)
       sfcwind, dummy = um.addcyclic(sfcwind, lonin)   
       trwind, dummy = um.addcyclic(trwind, lonin)
       plotvar = sfcwind
       cbarlabel = '2-meter AGL wind (m s-1)'
   
       cbar_max_trth = 35
       cbar_min_trth = 10
       cint_trth = 1
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth + (cint_trth/2), cint_trth)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.gist_heat_r
       
       titletext1 = 'Surface wind %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00')) 
   if plot_field == 'ghgt500' :
       ghgt500, dummy = um.addcyclic(ghgt500, lonin)   
       plotvar = ghgt500
       cbarlabel = '500 hPa heights (m)'
   
       cbar_max_trth = 5910
       cbar_min_trth = 4910
       cint_trth = 30
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth + (cint_trth/2), cint_trth)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.RdBu_r       
       #cmap_opt = plt.cm.jet
       contour2plot = 5420
       
       titletext1 = '500 hPa height %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00')) 
   if plot_field == 'trthgrad':
       dthdy, dthdx = wm.gradient_sphere(trth, lats, lons)
       plotvar = np.sqrt(dthdx**2 + dthdy**2)*1000
       plotvar = ndimage.gaussian_filter(plotvar,0.75)
       
       cbar_max_trth = 0.05
       cbar_min_trth = 0.01
       cint_trth = 0.001
       cflevs_trth =  np.arange(cbar_min_trth, cbar_max_trth, cint_trth)
       cflevs_trth_cntrs = cflevs_trth[::2]
       cflevs_trth_ticks = np.arange(cbar_min_trth,cbar_max_trth,4*cint_trth)
       cmap_opt = plt.cm.gist_heat_r          
       cbarlabel = 'Gradient (K km-1)'
   
       titletext1 = 'Tropopause potential temperature gradient %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00')) 
   
   base_cntr_slp = 996
   nslpconts = 20
   cint_slp = 4
   cbar_min_slp = base_cntr_slp-nslpconts*cint_slp
   cbar_max_slp = base_cntr_slp+nslpconts*cint_slp
   #cbar_max_slp = cbar_max_slp + (cint_slp/2)   
   cbar_max_slp = 1004
   cflevs_slp =  np.arange(cbar_min_slp, cbar_max_slp, cint_slp)
   cflevs_slp_ticks = np.arange(cbar_min_slp,cbar_max_slp,4*cint_slp)

   if smooth_plotvar == 'True':
       plotvar = ndimage.gaussian_filter(plotvar,0.75) 

 
   #titletext1 = 'Tropopause potential temperature valid %s at %s UTC' % (
   #        dt.strftime('%d %b %Y'), dt.strftime('%H00'))  
   #titletext1 = 'Tropopause potential temperature and SLP valid %s' % (dt.strftime('%Y%m%d%H'))     
   #titletext1 = 'Tropopause potential temperature %s at %s UTC' % (dt.strftime('%d %b %Y'), dt.strftime('%H00')) 
   
   # Set global figure properties
   golden = (np.sqrt(5)+1.)/2.   

   # Figure 1
   
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

         #coef = 0.65
	 #coef = 0.9
         width = width*coef
         height = height*coef
         m = Basemap(projection='ortho',lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],resolution='l',\
             llcrnrx=-0.5*width,llcrnry=-0.5*height,urcrnrx=0.5*width,urcrnry=0.5*height)		  
		  
   elif map_projection == 'lcc':
   #    m = Basemap(llcrnrlon=-125.5,llcrnrlat=15.,urcrnrlon=-30.,urcrnrlat=50.352,\
   #    m = Basemap(llcrnrlon=-125.5,llcrnrlat=15.,urcrnrlon=-30.,urcrnrlat=50.352,\
       m = Basemap(llcrnrlon=-120.0,llcrnrlat=20.,urcrnrlon=-60.0,urcrnrlat=50.0,\
               rsphere=(6378137.00,6356752.3142),\
               resolution='l',area_thresh=1000.,projection='lcc',\
               lat_1=50.,lon_0=-107.,ax=ax1)	       
   elif map_projection == 'npstere':
       if zoom == 'false':
          m = Basemap(projection='npstere',boundinglat=30,lon_0=proj_latlon[1],resolution='l')
       else:
          m = Basemap(projection='npstere',boundinglat=40,lon_0=proj_latlon[1],resolution='l')

   # draw countries, states, and differentiate land from water areas.
   m.drawcoastlines(linewidth=2, color='#444444', zorder=6)
   m.drawcountries(linewidth=1, color='#444444', zorder=5)
   m.drawstates(linewidth=0.66, color='#444444', zorder=4)
   m.drawmapboundary

   # draw lat/lon grid lines every 30 degrees.
   m.drawmeridians(np.arange(0, 360, 30))
   m.drawparallels(np.arange(-90, 90, 30))
   
   if plot_field == 'blank' :
       m.fillcontinents(color='tan',lake_color='lightskyblue')
       m.drawmapboundary(fill_color='lightskyblue')
       save_name = imagedir + 'blank_map' + '_' + date_string + ".png"
       plt.savefig(save_name, bbox_inches='tight')
       exit()
   x, y = m(X, Y)

   #trth = um.filter_numeric_nans(trth,340,float('NaN') ,'high')
   
   # Contour tropopause potential temperature
   
   CS1 = m.contourf(x,y,plotvar,cmap=cmap_opt,levels=cflevs_trth, extend='both',zorder=1)
   if plot_field == 'trtemp' :
       CS1.cmap.set_under('w')
   if plot_field == 'trth' :
       CS1.cmap.set_over('w')
   cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both',pad=0.05)
   
   #clabs = ['%i' % f for f in cflevs_trth]
   
   labels = [item.get_text() for item in cbar.ax.get_xticklabels()]
   cbar.ax.set_xticklabels(labels, size=20)
   cbar.set_label(cbarlabel,size=20)

   if plot_field == 'quivwind' :
       ugrid,newlons = shiftgrid(180.,uwind,lons,start=False)
       vgrid,newlons = shiftgrid(180.,vwind,lons,start=False)
       # transform vectors to projection grid.
       
       uproj,vproj,xx,yy = m.transform_vector(ugrid[::-1,:],vgrid[::-1,:],newlons,lats[::-1],61,61,returnxy=True,masked=True)
       # now plot.
       Q = m.quiver(xx,yy,uproj,vproj,scale=700,color='g',linewidths=(1,), edgecolors=('g'), headaxislength=2)
     
   if plot_trth_contours == 'True':
       CS2 = m.contour(x, y, trth, cflevs_trth_cntrs, colors='k', linewidths=0.25)
   
   if plot_slp == 'True':
      CS3 = m.contour(x, y, slp, cflevs_slp, colors='k', linewidths=1.5)
      #plt.clabel(CS3, inline=1, fontsize=10, fmt='%i')   
      plt.clabel(CS3, cflevs_slp, fmt = '%i', inline=True, fontsize=10)


   if plot_field == 'ghgt500':
      CS3 = m.contour(x, y, plotvar, [contour2plot, contour2plot], colors='green', linewidths=6.0)
    
     

   #ax1.set_title(titletext1)
   #plt.title('%s' % (titletext1), fontsize=label_fontsize,bbox=dict(facecolor='white', alpha=0.65),x=0.5,y=.95,weight = 'demibold',style='oblique', \
   #		stretch='normal', family='sans-serif')
   save_name = imagedir + figname + '_' + date_string + ".png"
   plt.savefig(save_name, bbox_inches='tight')

   record_start += 1
   rcount +=1
   record_numnow +=1

   datenow = um.advance_time(datenow,hinc)      

if nplots==1:
    plt.show()
