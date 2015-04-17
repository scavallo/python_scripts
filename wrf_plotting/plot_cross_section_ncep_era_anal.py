# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

#!/usr/bin/python

# imports
import netCDF4

import os, datetime
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um
import sounding_modules as sm
from scipy import ndimage

from mstats import *
###################################
# Set user options
###################################
#date_string = '2011122412'
#date_string = '2013123112' # yyyymmddhh
date_string_start = '2014070400'
data_gridnum = 3 # 1 for ncep analyses
                  # 2 for ncep real-time forecasts
		  # 3 for gfs final analysis
		  # 11 for era interim
plot_option = 11 # 1 for theta 
                # 2 for z
		# 3 for epv
		# 4 for v component of wind
		# 5 for wind magnitude
		# 6 for absolute angular momentum
		# 7 for heat flux
		# 8 for E-P flux
		# 9 for geopotential height
		# 10 for vertical motion (omega)
		# 11 for CAPE
		# 12 for temperature advection
		# 13 for temperature lapse rate
cross_method = 1 # 2 for Nick's, 1 otherwise		
plot_anomaly = 'false'
standardize_anomaly = 'true'
average_crosssec = 'false'
plot_logscale = 'true'
plot_isentropes = 'true'
plot_legend = 'false'
plot_background = 'trth'	
#record_num = 0
record_numbers = [0,2]
hinc = 6 # number of hours between record numbers

cenLat = 41.0
cenLon = 265.

#cenLat = 43.0
#cenLon = 295.0


#cenLat = 59.5
#cenLon = 0.

#cenLat = 45.0
#cenLon = 360.- 91.

#cenLat = 48.0
#cenLon = 360.-75.0

#cenLat = 48.5 cenLon = 299.5 
#cenLat = 51.5 cenLon = 309.5 
#cenLat = 56.5  cenLon = 309.5
#cenLat = 59.0 cenLon = 309.5
#cenLat = 61.5 cenLon = 309.5

dLat = 0.
dLon = 15.


#cenLat = 75.
#cenLon = 272.
#dLat = 30.
#dLon = 0.

#cenLat = 65.
#cenLon = 272.
#dLat = 0.
#dLon = 40.

#lat_range = [45.0,45.0]
#lon_range = [265,310]
#lon_range = [255,290]


#lat_range = [68.0,68.0]
#lon_range = [200.0,350.0]
#first_guess_latlon = [68.0,262.0]

#lat_range = [61.5,61.5]
#lat_range = [56.5,56.5]
#lon_range = [250.0,359.5]
#lat_range = [40.0,85.0]
#lon_range = [272.0,272.0]

#lat_range = [62.5,62.5]
#lon_range = [250.0,359.5]
#lat_range = [30.0,89.0]
#lon_range = [269.0,269.0]
#lon_range = [262.0,262.0]

find_vortex = 'false'
#first_guess_latlon = [44.0,280.0]
first_guess_latlon = [50.0,360.0-70.0]

#lat_range = [75.0,75.0]
#lon_range = [0.5,179.5]
#first_guess_latlon = [75.0,65.0]

boxsize = [5.0,5.0]
pres_range = [150,1000]
label_fontsize = 18



# Now provide the path to the directory containing the .nc file. Please note,
# do NOT include the .nc file in the path.
#imagedir = '/home/scavallo/chestnut_files/papers/polar_vortex/images/'
imagedir = '/home/scavallo/scripts/python_scripts/images/'
#imagedir = '/home/scavallo/chestnut_files/presentations/nwc_colloquium_jan2015/images/'
#imagedir = '/home/scavallo/Documents/web/wxdisc/20150114/'
#imagedir = '/home/scavallo/Documents/web/wxdisc/20150114/gfs_cross_0E_titles/'
#fname = 'gfs_analyses_all_trop_2014010700.nc'

#fpath = '/data1/scavallo/data/'
fpath = '/data1/scavallo/data/cases/july2014/'
#fpath = '/data1/scavallo/data/cases/decjan_2013/gfs_fnl/'
#fpath = '/data2/scavallo/era_interim/'
#fpath = '/data1/scavallo/data/cases/decjan_2013/gfs_fnl/'
#fpath = '/data1/scavallo/data/cases/jan2015/'

#fname = 'erainterim_pressure_2011122400.nc'
#fname = 'erainterim_pressure_2014010700.nc'
#fname = 'erainterim_pressure_2013123100.nc'
#fpath = '/data1/scavallo/data/save_realtime/tmp/'
#fname = 'gfs_analyses_2014022600_startf012.nc'
#fname = 'gfs_forecasts_2014022500_036_048.nc'

fname = 'gfs_analyses_2014040400_2014040418.nc'
#fname = 'gfs_analyses_2015020712_2015020800.nc'
#fname = 'gfs_analyses_2014010200_2014010700.nc'
#fname = 'gfs_analyses_2014010700_2014011200.nc'
#fname = 'gfs_analysis_2015010700_2015011200.nc'
#fname = 'gfs_analysis_2015010300_2015011200.nc'
###################################
# END user options
###################################


lat_range = [cenLat-dLat, cenLat+dLat]
lon_range = [cenLon-dLon,cenLon+dLon]

if cross_method == 1:
   if lat_range[1] >=90:
      lat_range[1] = 89.5

print lat_range
print lon_range

date_string = date_string_start
# Open the netcdf file and read select variables
dt = datetime.datetime.strptime(date_string, '%Y%m%d%H')
fpath = os.path.join(fpath, fname)
print fpath


if (lat_range[0] == lat_range[1]):
   cross_orient = 'east-west'
else:
   cross_orient = 'north-south'   


iternum = record_numbers[0]
iterend = record_numbers[1]
while iternum <= iterend:

    record_num = iternum
    

    f = netCDF4.Dataset(fpath,'r')
    if data_gridnum == 1:
       lons = f.variables['lon'][record_num,:].squeeze()
       lats = f.variables['lat'][record_num,::-1].squeeze() # Read in reverse direction
       levs = f.variables['PRES'][record_num,:,0,0].squeeze()
    elif data_gridnum == 2:
       lons = f.variables['lon_0'][:]
       lats = f.variables['lat_0'][::-1]
       levs = f.variables['lv_ISBL0'][:]
       levs = levs/100
    elif data_gridnum == 3:
       lons = f.variables['lon_3'][:]
       lats = f.variables['lat_3'][::-1]
       
       levs = f.variables['lv_ISBL3'][:]	   
       levs2 = f.variables['lv_ISBL9'][:]
	   
      
    else:
       lons = f.variables['longitude'][:].squeeze()
       lats = f.variables['latitude'][::-1].squeeze() # Read in reverse direction
       nloninds = np.ravel(lons<0)
       if cross_method == 1:
	  lons[nloninds] = lons[nloninds] + 360.
       levs = f.variables['level'][:]        

    if record_num > -1:
       if data_gridnum == 1:
	  u_plot =  f.variables['U'][record_num,:,::-1,:].squeeze()
	  v_plot =  f.variables['V'][record_num,:,::-1,:].squeeze() 
	  w_plot =  f.variables['W'][record_num,:,::-1,:].squeeze() 
	  theta = f.variables['THETA'][record_num,:,::-1,:].squeeze()
	  epv = f.variables['EPV'][record_num,:,::-1,:].squeeze()
	  trth = f.variables['TROP_THETA'][record_numnow,::-1,:].squeeze()

       elif data_gridnum == 2:
	  u_plot =  f.variables['UGRD_P0_L100_GLL0'][record_num,:,::-1,:].squeeze()
	  v_plot =  f.variables['VGRD_P0_L100_GLL0'][record_num,:,::-1,:].squeeze()
	  w_plot =  f.variables['VVEL_P0_L100_GLL0'][record_num,:,::-1,:].squeeze()
	  t_plot =  f.variables['TMP_P0_L100_GLL0'][record_num,:,::-1,:].squeeze()      
	  geop =  f.variables['HGT_P0_L100_GLL0'][record_num,:,::-1,:].squeeze()    
	  temptrop = f.variables['TMP_P0_L109_GLL0'][record_num,5,::-1,:].squeeze()
	  prestrop = f.variables['PRES_P0_L109_GLL0'][record_num,5,::-1,:].squeeze()  
	  p_in = np.zeros_like(t_plot).astype('f')        
	  z_dim = len(f.dimensions['lv_ISBL0'])
	  levs2 = f.variables['lv_ISBL5'][:]	  
	  levs2 = levs2/100 
	  p_in2 = np.zeros_like(w_plot).astype('f')    
	  z_dim2 = len(f.dimensions['lv_ISBL5'])

	  for kk in range(0,z_dim):      
             p_in[kk,:,:] = levs[kk] 

	  for kk in range(0,z_dim2):      
             p_in2[kk,:,:] = levs2[kk] 	  

	  theta = t_plot * (1000. / p_in) ** 0.286 
	  trth = temptrop * (100000. / prestrop) ** 0.286

	  epv = wm.epv_sphere(theta,p_in*100,u_plot,v_plot,lats,lons)
	  #epv = epv/10**6
       elif data_gridnum == 3:
	  u_plot =  f.variables['U_GRD_3_ISBL_10'][record_num,:,::-1,:].squeeze()
	  v_plot =  f.variables['V_GRD_3_ISBL_10'][record_num,:,::-1,:].squeeze()
	  w_plot =  f.variables['V_VEL_3_ISBL_10'][record_num,:,::-1,:].squeeze()
	  t_plot =  f.variables['TMP_3_ISBL_10'][record_num,:,::-1,:].squeeze()            
	  geop =  f.variables['HGT_3_ISBL_10'][record_num,:,::-1,:].squeeze()      
	  p_in = np.zeros_like(t_plot).astype('f')   
	  z_dim = len(f.dimensions['lv_ISBL3'])
	  for kk in range(0,z_dim):      
             p_in[kk,:,:] = levs[kk]  

	  theta = t_plot * (1000. / p_in) ** 0.286 
	  epv = wm.epv_sphere(theta,p_in*100,u_plot,v_plot,lats,lons)
	  #epv = epv/10**6      

	  trtempin = f.variables['TMP_3_PVL_10'][record_num,0,::-1,:].squeeze()
	  trpresin = f.variables['PRES_3_PVL_10'][record_num,0,::-1,:].squeeze()      
	  trth = wm.temp_to_theta(trtempin,trpresin)  
       else:
	  u_plot =  f.variables['u'][record_num,:,::-1,:].squeeze()
	  v_plot =  f.variables['v'][record_num,:,::-1,:].squeeze() 
	  w_plot =  f.variables['w'][record_num,:,::-1,:].squeeze() 
	  t_in = f.variables['t'][record_num,:,::-1,:].squeeze()

	  epv = f.variables['pv'][record_num,:,::-1,:].squeeze() 
	  geop = f.variables['z'][record_num,:,::-1,:].squeeze() 
	  epv = epv*10**6
	  p_in = np.zeros_like(t_in).astype('f')   
	  z_dim = len(f.dimensions['level'])
	  for kk in range(0,z_dim):      
             p_in[kk,:,:] = levs[kk]  

	  theta = t_in * (1000. / p_in) ** 0.286

	  erapvfile = fname[0:10] + '_pv_' + fname[20:]
	  fpath2 = os.path.join(fpath, erapvfile)
	  f2 = netCDF4.Dataset(fpath2,'r')
	  trth = f2.variables['th'][record_num,:,::-1,:].squeeze()
	  f2.close
    else:
       if data_gridnum == 1:
	  u_plot =  f.variables['U'][:,::-1,:].squeeze()
	  v_plot =  f.variables['V'][:,::-1,:].squeeze() 
	  w_plot =  f.variables['W'][:,::-1,:].squeeze() 
	  theta = f.variables['THETA'][:,::-1,:].squeeze()
	  epv = f.variables['EPV'][:,::-1,:].squeeze()
       else:      
	  u_plot =  f.variables['u'][:,::-1,:].squeeze()
	  v_plot =  f.variables['v'][:,::-1,:].squeeze() 
	  w_plot = f.variables['w'][:,::-1,:].squeeze()
	  t_in = f.variables['t'][:,::-1,:].squeeze()      
	  epv = f.variables['pv'][:,::-1,:].squeeze() 
	  geop = f.variables['z'][:,::-1,:].squeeze() 
	  p_in = np.zeros_like(t_in).astype('f')   
	  z_dim = len(f.dimensions['level'])
	  for kk in range(0,z_dim):      
             p_in[kk,:,:] = levs[kk]  

	  theta = t_in * (1000. / p_in) ** 0.286   

	  erapvfile = fname[0:10] + '_pv_' + fname[20:]
	  fpath2 = os.path.join(fpath, erapvfile)
	  f2 = netCDF4.Dataset(fpath2,'r')
	  trth = f2.variables['th'][:,::-1,:].squeeze()
	  f2.close      

    f.close

    nz, ny, nx = theta.shape



    pres = np.zeros_like(theta).astype('f')   
    for kk in range(0,len(levs)):      
       pres[kk,:,:] = levs[kk]

    tmp = theta[0,:,:].squeeze()
    latarr = np.zeros_like(tmp).astype('f')
    lonarr = np.zeros_like(tmp).astype('f')
    for ii in range(0,len(lons)):
       for jj in range(0,len(lats)):
	  latarr[jj,ii] = lats[jj]
	  lonarr[jj,ii] = lons[ii]

    wind_plot = np.sqrt(u_plot**2 + v_plot**2)

    if plot_option == 6:
       plot_anomaly = 'true'


    if plot_option == 8:
	Fx,Fy,Fz,divF = wm.eliassen_palm_flux_sphere(geop,theta,lats,lons,levs*100)    
	divF = ndimage.gaussian_filter(divF,1.75)  
	u_plot = divF
	finds = np.where( (epv<2.0) | (pres>200) )

	Fx[finds] = float('NaN') 
	Fy[finds] = float('NaN')      
	Fz[finds] = float('NaN')    


	Fx[:] = 0.005*Fx
	Fy[:] = 0.005*Fy 

	Fx_anom, Fx_anom_std = wm.spatial_anomaly(Fx,1)
	Fy_anom, Fy_anom_std = wm.spatial_anomaly(Fy,1)
	Fz_anom, Fz_anom_std = wm.spatial_anomaly(Fz,1)    

	#Fx = Fx_anom / Fx_anom_std
	#Fy = Fy_anom / Fy_anom_std
	#Fz = Fz_anom / Fz_anom_std

    heat_flux = np.zeros_like(theta).astype('f')
    heat_flux_anom = np.zeros_like(theta).astype('f')
    if plot_anomaly == 'true':

       epv_anom, epv_anom_std = wm.spatial_anomaly(epv,1)
       u_anom, u_anom_std = wm.spatial_anomaly(u_plot,1)
       v_anom, v_anom_std = wm.spatial_anomaly(v_plot,1)
       w_anom, w_anom_std = wm.spatial_anomaly(w_plot,1)
       theta_anom, theta_anom_std = wm.spatial_anomaly(theta,1)   
       wind_anom, wind_anom_std = wm.spatial_anomaly(wind_plot,1) 
       geop_anom, geop_anom_std = wm.spatial_anomaly(geop,1) 

       heat_flux = theta_anom * v_anom
       heat_flux_anom, heat_flux_std = wm.spatial_anomaly(heat_flux,1)               

       if standardize_anomaly == 'true':
	  theta_anom = theta_anom_std      
	  epv_anom = epv_anom_std
	  u_anom = u_anom_std
	  v_anom = v_anom_std
	  geop_anom = geop_anom_std
	  wind_anom = wind_anom_std
	  heat_flux = heat_flux_std



    if cross_orient == 'east-west':
       orient = 'ew'
       latinds =  np.ravel(lats==lat_range[0])
       loninds = np.ravel((lons<=lon_range[1])&(lons>=lon_range[0]))
       x_cross = lons[loninds]

       epv_cross = epv[:,latinds,loninds].squeeze()
       theta_cross = theta[:,latinds,loninds].squeeze()
       geop_cross = geop[:,latinds,loninds].squeeze()   
       u_cross = u_plot[:,latinds,loninds].squeeze()
       v_cross = v_plot[:,latinds,loninds].squeeze()
       w_cross = w_plot[:,latinds,loninds].squeeze()
       wind_cross = wind_plot[:,latinds,loninds].squeeze()
       vt_cross = heat_flux[:,latinds,loninds].squeeze()
       if plot_option == 8:
	   Fx_cross = Fx[:,latinds,loninds].squeeze()
	   Fz_cross = Fz[:,latinds,loninds].squeeze()

       minlon = first_guess_latlon[1] - boxsize[1]
       maxlon = first_guess_latlon[1] + boxsize[1]   
       latloninds = np.where( (x_cross>=minlon) & (x_cross<=maxlon) )

       
       if data_gridnum == 2:
           levelindex = np.ravel(levs==425)
       else:
           levelindex = np.ravel(levs==400)

       
       vnow = v_cross[levelindex,latloninds]
       vnow = np.abs(vnow)
       thetanow = theta_cross[levelindex,latloninds]


       lonsnow = x_cross[latloninds]   
       
       angular_mom  = np.zeros_like(v_cross).astype('f') 
       if 1 == 0:
	   mininds=np.where(vnow==np.min(vnow[:]))     
	   #mininds=np.where(thetanow==np.min(thetanow[:]))     

	   first_guess_latlon[1] = lonsnow[mininds[1]]


	   vmin = np.zeros_like(levs).astype('f') 
	   thetamin = np.zeros_like(levs).astype('f') 
	   lonsmin = np.zeros_like(levs).astype('f')   
	   rdist = np.zeros_like(theta_cross).astype('f') 
	   if find_vortex == 'true':       
	      for kk in range(0,nz):
		  vnow = v_cross[kk,latloninds] 	   
		  vnow = np.abs(vnow)    
		  thetanow = theta_cross[kk,latloninds] 	
		  lonsnow = x_cross[latloninds]                  

		  mininds=np.where(vnow==np.min(vnow[:]))            
		  #mininds=np.where(thetanow==np.min(thetanow[:]))   

		  try:
        	     vmin[kk] = vnow[mininds]
		     lonsmin[kk] = lonsnow[mininds[1]]
		     thetamin[kk] = thetanow[mininds[kk]]	   
		     first_guess_latlon[1] = lonsmin[kk]
		  except:
		     vmin[kk] = float('NaN')
		     thetamin[kk] = float('NaN')
		     lonsmin[kk] = float('NaN')


		  minlon = first_guess_latlon[1] - boxsize[1]
		  maxlon = first_guess_latlon[1] + boxsize[1]   
		  latloninds = np.where( (x_cross>=minlon) & (x_cross<=maxlon) )	   


		  rdist[kk,:] = (x_cross - lonsmin[kk])*111000*np.cos(lat_range[0]*(np.pi/180))


	   omeg_e = (2*np.pi) / (24*3600);
	   f = 2*omeg_e*np.cos(lat_range[0]*(np.pi/180))
	   angular_mom = rdist*v_cross + ( (f*rdist**2)/2 )   

       if plot_anomaly == 'true':
	  epv_cross_anom = epv_anom[:,latinds,loninds].squeeze()
	  theta_cross_anom = theta_anom[:,latinds,loninds].squeeze()          
	  u_cross_anom = u_anom[:,latinds,loninds].squeeze()  
	  v_cross_anom = v_anom[:,latinds,loninds].squeeze()  
	  w_cross_anom = w_anom[:,latinds,loninds].squeeze()  
	  geop_cross_anom = geop_anom[:,latinds,loninds].squeeze()  
	  vt_cross_anom = vt_cross
	  wind_cross_anom = wind_anom[:,latinds,loninds].squeeze()   


    else:
       orient = 'ns'
       if average_crosssec == 'false':
	  loninds =  np.ravel(lons==lon_range[0])
       else:
	  loninds = np.ravel((lons<=lon_range[1])&(lons>=lon_range[0]))

       if cross_method == 2:  
	  print cenLat, cenLon, dLat, dLon

	  latinds, loninds = um.driver_xsec_old(lats, lons, cenLat, cenLon, dLat, dLon)
	  print lats[latinds]
	  print lons[loninds]

	  xArray = um.calc_distSphere_multiple(6371., lats[latinds[0]], lons[loninds[0]], lats[latinds], lons[loninds])      
	  x_cross, levs = np.meshgrid(xArray, levs)# indexing='ij');
	  f = np.zeros_like(x_cross).astype('f') 
       else:   
	  latinds = np.ravel((lats<=lat_range[1])&(lats>=lat_range[0]))
	  x_cross = lats[latinds]

       geop_cross = geop[:,latinds,loninds].squeeze()
       epv_cross = epv[:,latinds,loninds].squeeze()
       theta_cross = theta[:,latinds,loninds].squeeze()
       u_cross = u_plot[:,latinds,loninds].squeeze()
       v_cross = v_plot[:,latinds,loninds].squeeze()
       w_cross = w_plot[:,latinds,loninds].squeeze()
       uv_cross = wind_plot[:,latinds,loninds].squeeze()
       vt_cross = heat_flux[:,latinds,loninds].squeeze()
       if plot_option == 8:
	   Fx_cross = Fy[:,latinds,loninds].squeeze()
	   Fz_cross = Fz[:,latinds,loninds].squeeze()


       minlat = first_guess_latlon[0] - boxsize[0]
       maxlat = first_guess_latlon[0] + boxsize[0]   
       latloninds = np.where( (x_cross>=minlat) & (x_cross<=maxlat) )

       thetamin = np.zeros_like(levs).astype('f') 
       lonsmin = np.zeros_like(levs).astype('f') 
       rdist = np.zeros_like(theta_cross).astype('f')     
       if find_vortex == 'true':
	  for kk in range(0,nz):
	      thetanow = theta_cross[kk,latloninds]
	      lonsnow = x_cross[latloninds]                            
	      if levs[kk] >= 0.5:                 
        	  mininds=np.where(thetanow==np.min(thetanow[:]))            	   
        	  thetamin[kk] = thetanow[mininds]
		  lonsmin[kk] = lonsnow[mininds[1]]
		  #thetamin[kk] = 0.
		  #lonsmin[kk] = 0.	   
		  first_guess_latlon[0] = lonsmin[kk]


		  minlat = first_guess_latlon[0] - boxsize[0]
        	  maxlat = first_guess_latlon[0] + boxsize[0]   
        	  latloninds = np.where( (x_cross>=minlat) & (x_cross<=maxlat) )	   
	      else:
        	  thetamin[kk] = float('NaN')
		  lonsmin[kk] = float('NaN')	   

	      rdist[kk,:] = (x_cross - lonsmin[kk])*111000



       omeg_e = (2*np.pi) / (24*3600);
       nzz, nyy = theta_cross.shape
       for jj in range(0,nyy): 

	  if cross_method == 2:
              f[:,jj] = 2*omeg_e*np.cos(x_cross[:,jj]*(np.pi/180))
	  else:
              f = 2*omeg_e*np.cos(x_cross[jj]*(np.pi/180))      

	  angular_mom = rdist*u_cross + ( (f*rdist**2)/2 )          

       if average_crosssec == 'true':      
	  epv_cross = np.mean(z_cross,axis=2)
	  theta_cross = np.mean(theta_cross,axis=2)           
	  u_cross = np.mean(u_cross,axis=2)
	  v_cross = np.mean(v_cross,axis=2)
	  w_cross = np.mean(w_cross,axis=2)
	  geop_cross = np.mean(geop_cross,axis=2)
	  wind_cross = np.mean(wind_cross,axis=2)       

       if plot_anomaly == 'true':
	  epv_cross_anom = epv_anom[:,latinds,loninds].squeeze()
	  theta_cross_anom = theta_anom[:,latinds,loninds].squeeze()
	  u_cross_anom = u_anom[:,latinds,loninds].squeeze()
	  v_cross_anom = v_anom[:,latinds,loninds].squeeze()
	  w_cross_anom = w_anom[:,latinds,loninds].squeeze()
	  geop_cross_anom = geop_anom[:,latinds,loninds].squeeze()
	  wind_cross_anom = wind_anom[:,latinds,loninds].squeeze()
	  vt_cross_anom = vt_cross

    #u_cross = u_cross * 1.94384
    #v_cross = v_cross * 1.94384
    #uv_cross = uv_cross * 1.94384

    mom = u_cross*v_cross
    absmom = np.abs(mom)


    angular_mom = np.abs(angular_mom)
    angular_mom = angular_mom/10**3


    if plot_option == 1 :
	plot_cross = theta_cross

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
	figname = "erainterim_cross_section_analysis_theta_" + orient
    elif plot_option == 2:
	plot_cross = z_cross

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
	figname = "erainterim_cross_section_analysis_z_" + orient
    elif plot_option == 3:
	plot_cross = epv_cross

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
	figname = "erainterim_cross_section_analysis_epv_" +orient	        
    elif plot_option == 4:
	if cross_orient == 'east-west':
            plot_cross = v_cross
	    titlestring = "V-component wind anomaly " + date_string
	else:
            plot_cross = u_cross    
	    titlestring = "U-component wind anomaly " + date_string

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
            cint = 5
            cbar_min = -70
            cbar_max = 70+(cint/2)            

	    cbar_labels = 'm s-1'
	    #titlestring = "Meridional wind " + date_string
	    #date_forecast[6:8] + ' January ' + date_forecast[0:4]
	    titlestring = "Meridional wind " + date_string[8:10] + ' UTC ' + date_string[6:8] + ' January ' + date_string[0:4]

	    cmap_opt = plt.cm.RdBu_r 
	figname = "erainterim_cross_section_analysis_vwind_" + orient	        	        
    elif plot_option == 5:
	plot_cross = wind_cross

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
            cint = 1
            cbar_min = 5
            cbar_max = 50+(cint/2)            

	    cbar_labels = 'm s-1'
	    titlestring = "Wind magnitude " + date_string

	    cmap_opt = plt.cm.Reds 
	figname = "erainterim_cross_section_analysis_windmag_" + orient	        	
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
	figname = "erainterim_cross_section_analysis_absangmom_" + orient	        
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


	figname = "erainterim_cross_section_analysis_epflux_" + orient	        
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

	figname = "erainterim_cross_section_analysis_ghgt_" + orient	        
    elif plot_option == 10:
	plot_cross = w_cross
	mstats(plot_cross)

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
            cint = 0.05
            cbar_min = -1.0
            cbar_max = 1.0+(cint/2)            

	    cbar_labels = 'Pa s-1'
	    titlestring = "Vertical velocity " + date_string

	    cmap_opt = plt.cm.RdBu 	
	figname = "erainterim_cross_section_analysis_omega_" + orient	        

    elif plot_option == 11:                 
        f = netCDF4.Dataset(fpath,'r')
	temp_cross = t_plot[:,latinds,loninds].squeeze()
	pres_cross = p_in[:,latinds,loninds].squeeze()*100    
	relh_cross = np.zeros_like(temp_cross).astype('f')   
	relh_cross[2:,:] = f.variables['R_H_3_ISBL_10'][record_num,1:,latinds,loninds].squeeze()
	hght_cross = f.variables['HGT_3_ISBL_10'][record_num,:,latinds,loninds].squeeze()
	f.close

	es = wm.claus_clap(temp_cross)
	qvs = wm.satur_mix_ratio(es, pres_cross)
	qv = (relh_cross/100)*qvs
	dewpt_cross = wm.mixrat_to_td(qv, pres_cross)


        #dewpt_cross = np.zeros_like(temp_cross).astype('f')   
	plot_cross = np.zeros_like(temp_cross).astype('f')   

	[iz,ix] = np.shape(temp_cross)

	startp = 100000.    
	capearr = np.zeros_like(temp_cross).astype('f')   
	cinarr = np.zeros_like(temp_cross).astype('f')  
	
	pres_cross2 = pres_cross[::-1,:] 
	temp_cross2 = temp_cross[::-1,:] 
        dewpt_cross2 = dewpt_cross[::-1,:] 
	hght_cross2 = hght_cross[::-1,:] 
	for kk in range(0,iz-1):            
	#for kk in range(iz-1,0,-1):
            print kk   
            for jj in range(0,ix):       
		startp = pres_cross2[kk,jj]
		startt = temp_cross2[kk,jj]
		startdp = dewpt_cross2[kk,jj]
		if ( np.isnan(startp) or np.isnan(startt) or  np.isnan(startdp)) :
	            cape = float('NaN')  	
		else:   	
	            if startdp > startt:
			startdp = startt
                    lcl,lfc,el,cape,cin = sm.get_cape(temp_cross2[:,jj],pres_cross2[:,jj],dewpt_cross2[:,jj],hght_cross2[:,jj],startp,startt,startdp,totalcape=False)
        	#print lcl,lfc,el,cape,cin
		capearr[kk,jj] = cape
		try:
	            cinarr[kk,jj] = cin
		except:
	            cinarr[kk,jj] = float('NaN')   

	mstats(capearr)
	
	plot_cross = capearr[::-1,:] 
	plot_cross2 = cinarr[::-1,:]  
	 
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
            cint = 100
            cbar_min = 500
            cbar_max = 3550           

            cint2 = 5
            cbar_min2 = 0
            cbar_max2 = 100         
	    cflevs_cin = np.arange(cbar_min2,cbar_max2+(cint2/2),cint2)

	    cbar_labels = 'J kg-1'
	    titlestring = "CAPE " + date_string

	    cmap_opt = plt.cm.gist_heat_r 	
	    cmap_opt2 = plt.cm.cool_r

            figname  = "era_cross_section_cape_" + orient	 

    elif plot_option == 12:
        #f = netCDF4.Dataset(fpath,'r')
	#temp_cross = t_plot[:,latinds,loninds].squeeze()
	#u_cross = u_plot[:,latinds,loninds].squeeze()
	#v_cross = u_plot[:,latinds,loninds].squeeze()
	#f.close    
        
	tadv_plot =  np.zeros_like(t_plot).astype('f') 
	for kk in range(0,nz-1):              	    
	    print kk
	    tadv_plot[kk,:,:] = wm.hadvection_latlon(u_plot[kk,:,:].squeeze(), v_plot[kk,:,:].squeeze(), t_plot[kk,:,:].squeeze(), lats, lons)	 

	plot_cross = tadv_plot[:,latinds,loninds].squeeze() 
	plot_cross = plot_cross * 3600 * 24
	mstats(t_plot)
	mstats(plot_cross)


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
            cint = 0.5
            cbar_min = -10.0
            cbar_max = 10.0+(cint/2)            

	    cbar_labels = 'K day-1'
	    titlestring = "Temperature advection " + date_string

	    cmap_opt = plt.cm.RdBu_r 	
	figname = "erainterim_cross_section_analysis_tempadv_" + orient	        	 
    elif plot_option == 13:
        f = netCDF4.Dataset(fpath,'r')
	hght_plot = f.variables['HGT_3_ISBL_10'][record_num,:,:,:].squeeze()
	f.close 	  

        tlapse_plot, dfdlat, dfdlon = wm.gradient_sphere(t_plot, hght_plot, lats, lons)        	 
	plot_cross = tlapse_plot[:,latinds,loninds].squeeze() 
	plot_cross = plot_cross * 1000	
	mstats(plot_cross)



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
            cbar_min = 5
            cbar_max = 11+(cint/2)            

	    cbar_labels = 'K km-1'
	    titlestring = "Temperature lapse rate " + date_string

	    #cmap_opt = plt.cm.RdBu_r 
	    cmap_opt = plt.cm.gist_heat_r 		
	figname = "erainterim_cross_section_analysis_templapse_" + orient	        	 

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
	figname = "erainterim_cross_section_analysis_vtcross_" + orient	        


    # Set global figure properties
    golden = (np.sqrt(5)+1.)/2.
    figprops = dict(figsize=(8., 16./golden), dpi=128)


    if ( pres_range[0] > 10):
       cbar_min_theta = 220
       cbar_max_theta = 2100
       cint_theta = 5

       cbar_min_theta = 220
       #cbar_max_theta = 331
       cbar_max_theta = 401
       cint_theta = 2   
       cflevs_theta =  np.arange(cbar_min_theta, cbar_max_theta, cint_theta)
    else:
       arr1 = np.arange(250,360,5)   
              
       arr2 = np.arange(360,500,10)
       arr3 = np.arange(550,2000,50)
       arr4 = np.arange(250,331,2)
       #cflevs_theta = np.concatenate((arr1,arr2,arr3), axis=0)
       cflevs_theta = arr1
       cbar_min_theta = cflevs_theta[0]
       cbar_max_theta = cflevs_theta[::-1]   


    cflevs = np.arange(cbar_min, cbar_max, cint)


    #cflevs_std =  np.arange(cbar_min_std, cbar_max_std, cint_std)
    #cflevs_epv =  np.arange(cbar_min_epv, cbar_max_epv, cint_epv)
    #cflevs_epv_ticks = np.arange(cbar_min_epv,cbar_max_epv,4*cint_epv)

    cflevs_trop = [2.0,2.0]
    #cflevs_strat = [95.0,95.0] 
    cflevs_strat = [5000.0,5000.0] 

    #cflevs_winds = np.arange(50,200,5)
    #cflevs_temp = np.arange(200,501,3)
    #cflevs_theta = np.arange(cbar_min_theta,cbar_max_theta,cint_theta)
    #cflevs_theta_anom = np.arange(cbar_min_theta_anom,cbar_max_theta_anom,cint_theta_anom)   



    if data_gridnum == 1:
       levs_hPa = levs/100
    else:
       levs_hPa = levs


    if plot_background == 'trpr':   
	lonin = lons
	trth, lons = um.addcyclic(trth, lonin)
	plotvar_bg = trprin

	cint_bg = 10
	cbar_min_bg = 300
	cbar_max_bg = 850
	cmap_opt_bg = plt.cm.gist_heat_r 
	cflevs_bg = np.arange(cbar_min_bg,cbar_max_bg+(cint_bg/2),cint_bg)
    if plot_background == 'trth':
	lonin = lons

	trth, lons = um.addcyclic(trth, lonin)        

	plotvar_bg = trth


	latarr = np.zeros_like(trth).astype('f')
	lonarr = np.zeros_like(trth).astype('f')
	for ii in range(0,len(lons)):
	   for jj in range(0,len(lats)):
               latarr[jj,ii] = lats[jj]
               lonarr[jj,ii] = lons[ii]

	cint_bg = 1
	cbar_min_bg = 270
	cbar_max_bg = 380
	cmap_opt_bg = plt.cm.jet 
	cflevs_bg = np.arange(cbar_min_bg,cbar_max_bg+(cint_bg/2),cint_bg)       

    plotvar_bg = um.filter_numeric_nans(plotvar_bg,cbar_max_bg-(cint_bg/2),cbar_max_bg,'high')
    plotvar_bg = um.filter_numeric_nans(plotvar_bg,cbar_min_bg,cbar_min_bg,'low')


    golden = (np.sqrt(5)+1.)/2.
    #m = Basemap(projection='npstere',boundinglat=30,lon_0=270.,resolution='l')

    # Get the cross section points.  This is a function in utilities_modules.

    if cross_method == 1:
	m = Basemap(projection='npstere',boundinglat=30,lon_0=270.,resolution='l')

	
	xx, yy = um.xsection_inds(lon_range[0], lat_range[0], lon_range[1], lat_range[1], lonarr, latarr, m)
	#x, y = m(lonarr, latarr)
	# Plot plan view map with solid white line indicating the location of the cross section
	fig = plt.figure(figsize=(8., 16./golden), dpi=128)   # New figure
	ax0 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
	F = plt.gcf()  # Gets the current figure


	m.drawstates(color='#444444', linewidth=1.25)
	m.drawcoastlines(linewidth=1.25)
	#m.fillcontinents(color='peru',lake_color='aqua', zorder=1)
	m.drawcountries(color='#444444', linewidth=1.25)
	m.drawmapboundary
	m.drawparallels(np.arange(-80.,81.,10.))
	m.drawmeridians(np.arange(-180.,181.,10.))

        x, y = m(lonarr, latarr)
	if plot_background == 'trth':        
	    
            cmap0 = m.contourf(x, y, plotvar_bg, cmap=cmap_opt_bg, levels=cflevs_bg, extend='both',zorder=1)
            cbar0 = plt.colorbar(cmap0,shrink=0.95, orientation='horizontal',extend='both',pad=0.05)    
	    cbar0.set_label('Kelvin',size=20)

	    ax0.plot(x[xx,yy], y[xx,yy], color='k', lw=6)

	else:
            m.fillcontinents(color='peru',lake_color='aqua', zorder=1)
	    m.drawmapboundary(fill_color='aqua')
	    ax0.plot(x[xx,yy], y[xx,yy], color='w', lw=4)

	#ax0.plot(x[xx,yy], y[xx,yy], color='w', lw=4)
	ax0.title.set_y(1.00)
	ax0.set_title('Area of Cross Section', size=20)
    else:
	d = um.calc_distSphere_multiple(6371., lats[latinds[0]], lons[loninds[0]], lats[latinds], lons[loninds])    
	latMin = min(lats[latinds])

	m = Basemap(projection='npstere',boundinglat=latMin-10,lon_0=0,resolution='l')
	x,y = m(lons[loninds],lats[latinds])

	plt.figure()
	m.drawcoastlines(linewidth=0.5)
	m.drawmapboundary()
	m.plot(x,y,'k', linewidth=2.0)
	m.plot(x[0], y[0], 'go')
	m.plot(x[-1], y[-1], 'rx')

	m.drawparallels([70, 80]) # draw parallels
	m.drawmeridians(np.arange(0.,360.,60.)) # draw meridians

	#space markers along specified distance
	props = dict(boxstyle='round', facecolor='white', alpha=.6)
	dMarker = 500.;
	dLast = 0.
	for iMarker in xrange(1,len(d)):
	  dLast += d[iMarker]-d[iMarker-1]
	  if (dLast>=dMarker):
	    labl = '{0}'.format(int(d[iMarker]))
	    plt.text(x[iMarker], y[iMarker], labl, color='b', size=12, bbox=props)
	    dLast = 0.;



    plt.savefig(imagedir + 'area_trth_' + date_string + '.png' , bbox_inches='tight')


    fig = plt.figure(figsize=(8., 16./golden), dpi=128)   # New figure
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])      


    #plot_cross = um.filter_numeric_nans(plot_cross,cbar_max-(cint/2),cbar_max-(cint/2),'high')
    
    plot_cross = um.filter_numeric_nans(plot_cross,cbar_min,cbar_min,'low')    
    plot_cross = um.filter_numeric_nans(plot_cross,np.max(cflevs),np.max(cflevs),'high')   

   
    mstats(levs_hPa)
    mstats(plot_cross)
    
    
    xpts = [x_cross[0],x_cross[1],x_cross[2]]
    ypts = [levs_hPa[0],levs_hPa[1],levs_hPa[2]]
    if plot_legend == 'true':
	CSa = plt.plot(xpts,ypts,color='b',linewidth=1.5, label='Potential temperature <= 360 K')
    if plot_option != 8:

	#plot_cross = um.filter_numeric_nans(plot_cross,cbar_min,cbar_min,'low')    
	#plot_cross = um.filter_numeric_nans(plot_cross,np.max(cflevs),np.max(cflevs),'high')

	try:
            CS1 = plt.contourf(x_cross,levs_hPa,plot_cross,cmap=cmap_opt,levels=cflevs)
	except:
            #pres_range[0] = 100.	    
            CS1 = plt.contourf(x_cross,levs2,plot_cross,cmap=cmap_opt,levels=cflevs)    	
    else:
	CS1 = plt.contourf(x_cross,levs_hPa,plot_cross,cmap=cmap_opt,levels=cflevs)
	Fx_cross[:,::2] = float('NaN')
	Fz_cross[:,::2] = float('NaN')
	Fx_cross[:,0:4:] = float('NaN')
	Fz_cross[:,0:4:] = float('NaN')    
	CSQ = plt.quiver(x_cross,levs_hPa,Fx_cross,Fz_cross,color='k', units='x', linewidths=(1,), edgecolors=('k'), headaxislength=4 )



    CS1.cmap.set_under('w')
    #if find_vortex == 'true':
    #    CS2 = plt.scatter(lonsmin,levs_hPa,c='k',s=100,zorder=10)

    cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both',pad=.10)
    
    cbar.set_label(cbar_labels)

    if plot_isentropes == 'true':
       CS3 = plt.contour(x_cross,levs_hPa,theta_cross,levels=cflevs_theta,colors='b', linestyles='solid',linewidths=1.5)  
       if plot_legend == 'true':
	  legend = ax1.legend(loc='center right', shadow=True)
	  legend.set_zorder(20)

    CS4 = plt.contour(x_cross,levs_hPa,epv_cross,levels=cflevs_trop,colors='k', linestyles='solid',linewidths=2.5)  
    #CS5 = plt.contour(x_cross,levs_hPa,epv_cross,levels=cflevs_strat,colors='k', linestyles='solid',linewidths=2.5)  

    ymin = np.min(pres_range)
    ymax = np.max(pres_range)

    if cross_method == 1:
       plt.ylim((pres_range[0],pres_range[1]))
       ax1.set_ylim(ax1.get_ylim()[::-1])
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
	   levs_label = [150, 250, 500, 700, 850, 1000]  

       #if cross_orient == 'east-west':
       #   plt.xlim((lon_range[0],lon_range[1]))
       #else:
       #   plt.xlim((lat_range[0],lat_range[1]))

       if plot_logscale == 'true':
	  ax1.set_yscale('log')   
	  ycint = 200
	  #yticks = np.arange(ymin,ymax+(ycint/10),ycint)

	  ax1.set_yticks(levs_label)
	  ax1.set_yticklabels(levs_label, size=10)
    else:
       ax1.set_yscale('log')  
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

       ax1.set_yticks(levs_label)
       ax1.set_yticklabels(levs_label, size=14)   
       plt.ylim((pres_range[0],pres_range[1]))
       plt.gca().invert_yaxis()   

    if cross_orient == 'east-west':
       plt.xlabel('Longitude (Degrees)',fontsize=label_fontsize)
    else:
       plt.xlabel('Latitude (Degrees North)',fontsize=label_fontsize)
    plt.ylabel('Pressure (hPa)',fontsize=label_fontsize)


    ax1.set_title(titlestring, fontsize=22)

    save_name = figname + "_" + date_string + ".png"
    plt.savefig(imagedir + save_name, bbox_inches='tight')
    if record_numbers[0] == record_numbers[1]:
        plt.show()
    
    date_string = um.advance_time(date_string,hinc) 
    iternum += 1


