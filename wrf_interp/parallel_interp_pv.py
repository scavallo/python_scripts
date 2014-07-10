# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

#from IPython.parallel import Client
import numpy 
from scipy.interpolate import splprep, splev
from netCDF4 import Dataset
import datetime
import os
import gc
import subprocess

#subprocess.Popen("ipcluster start -n 4", shell=True) 
#rc = Client()
#lview = rc.load_balanced_view()
#dview = rc[:]
#lview.block = True
#dview.block = True

# <codecell>

def wrf_unstagger( grid, dim ):
    """ Unstagger a staggered WRF grid in the X, Y, or Z (U, V, or W) direction.
    ---------------------
    grid (numpy.ndarray): The 2D, 3D, 4D, or 5D array to be unstaggered.
    dim (str): A string specifying what dimension to unstagger. Must be
    X, Y, Z, U, V or W.
    ---------------------
    returns:
        numpy.ndarray unstaggered grid (dim-1)
    ---------------------
    EXAMPLE:
        import numpy as np
        
        arr = np.random.randint( low=1, high=10, size=( 9,10,9 ) ) ## create a random array staggered in the Y direction
        arr_unstaggered = wrf_unstagger( arr, 'Y' )
        """
    nd = len( grid.shape )
    if dim == 'X' or dim == 'U':
        if nd == 5:
            gridout = ( grid[ :, :, :, :, :-1 ] + grid[ :, :, :, :, 1: ] ) / 2
        elif nd == 4:
            gridout = ( grid[ :, :, :, :-1 ] + grid[ :, :, :, 1: ] ) / 2
        elif nd == 3:
            gridout = ( grid[ :, :, :-1 ] + grid[ :, :, 1: ] ) / 2
        elif nd == 2:
            gridout = ( grid[ :, :-1 ] + grid[ :, 1: ] ) / 2
        else: pass
    if dim == 'Y' or dim == 'V':
        if nd == 5:
            gridout = ( grid[ :, :, :, :-1, : ] + grid[ :, :, :, 1:, : ] ) / 2
        elif nd == 4:
            gridout = ( grid[ :, :, :-1, : ] + grid[ :, :, 1:, : ] ) / 2
        elif nd == 3:
            gridout = ( grid[ :, :-1, : ] + grid[ :, 1:, : ] ) / 2
        elif nd == 2:
            gridout = ( grid[ :-1, : ] + grid[ 1:, : ] ) / 2
        else: pass
    if dim == 'Z' or dim == 'W':
        if nd == 5:
            gridout = ( grid[ :, :, :-1, :, : ] + grid[ :, :, 1:, :, : ] ) / 2
        elif nd == 4:
            gridout = ( grid[ :, :-1, :, : ] + grid[ :, 1:, :, : ] ) / 2
        elif nd == 3:
            gridout = ( grid[ :-1, :, : ] + grid[ 1:, :, : ] ) / 2
        else: pass
    return gridout

def wrf_copy_attributes( infilename, outfilename ):
    ## open the files
    infile = Dataset( infilename )
    outfile = Dataset( outfilename, 'w', format='NETCDF4' )
    
    ## create dimensions
    level = outfile.createDimension( 'bottom_top', None )
    time = outfile.createDimension( 'time', None )
    lon = outfile.createDimension( 'south_north', infile.getncattr('SOUTH-NORTH_PATCH_END_UNSTAG') )
    lat = outfile.createDimension( 'west_east', infile.getncattr('WEST-EAST_PATCH_END_UNSTAG') )
    
    ## copy the global attributes to the new file
    inattrs = infile.ncattrs()
    for attr in inattrs:
        outfile.setncattr( attr, infile.getncattr( attr ) )
    infile.close()
    outfile.close()

def wrf_copy_static_fields( infilename, outfilename ):
    
    infile = Dataset( infilename )
    outfile = Dataset( outfilename, 'a' )

    XLONG = outfile.createVariable( 'XLONG', 'f4', ('time', 'south_north', 'west_east') )
    XLAT = outfile.createVariable( 'XLAT', 'f4', ('time', 'south_north', 'west_east') )
    SINALPHA = outfile.createVariable( 'SINALPHA', 'f4', ('time', 'south_north', 'west_east') )
    COSALPHA = outfile.createVariable( 'COSALPHA', 'f4', ('time', 'south_north', 'west_east') )
    MAPFAC_M = outfile.createVariable( 'MAPFAC_M', 'f4', ('time', 'south_north', 'west_east') )
    

    XLONG[:] = infile.variables['XLONG'][:]
    XLAT[:] = infile.variables['XLAT'][:]
    SINALPHA[:] = infile.variables['SINALPHA'][:]
    COSALPHA[:] = infile.variables['COSALPHA'][:]
    MAPFAC_M[:] = infile.variables['MAPFAC_M'][:]
    
    #F = outfile.createVariable( 'F', 'f4', ('south_north', 'west_east') )
    XTIME = outfile.createVariable( 'XTIME', 'f4', ('time') )

    #F[:] = infile.variables['F'][0][:]
    #MU = infile.variables['MU'][0][:]
    #MUB = infile.variables['MUB'][0][:]

    XTIME[:] = infile.variables['XTIME'][:]
    infile.close()
    outfile.close()

def wrf_copy_sfc_fields( infilename, outfilename ):
    
    infile = Dataset( infilename )
    outfile = Dataset( outfilename, 'a' )
    
    T2 = outfile.createVariable( 'T2', 'f4', ('time', 'south_north', 'west_east') )
    TH2 = outfile.createVariable( 'TH2', 'f4', ('time', 'south_north', 'west_east') )
    Q2 = outfile.createVariable( 'Q2', 'f4', ('time', 'south_north', 'west_east') )
    PSFC = outfile.createVariable( 'PSFC', 'f4', ('time', 'south_north', 'west_east') )
    U10 = outfile.createVariable( 'U10', 'f4', ('time', 'south_north', 'west_east') )
    V10 = outfile.createVariable( 'V10', 'f4', ('time', 'south_north', 'west_east') )
    SNOWH = outfile.createVariable( 'SNOWH', 'f4', ('time', 'south_north', 'west_east') )
    SEAICE = outfile.createVariable( 'SEAICE', 'f4', ('time', 'south_north', 'west_east') )
    RAINNC = outfile.createVariable( 'RAINNC', 'f4', ('time', 'south_north', 'west_east') )
    SNOWNC = outfile.createVariable( 'SNOWNC', 'f4', ('time', 'south_north', 'west_east') )
    #REFL_10CM = outfile.createVariable( 'REFL_10CM', 'f4', ('time', 'south_north', 'west_east') )
    
    T2[:] = infile.variables['T2'][:]
    TH2[:] = infile.variables['TH2'][:]
    Q2[:] = infile.variables['Q2'][:]
    PSFC[:] = infile.variables['PSFC'][:] * .01
    U10[:] = infile.variables['U10'][:]
    V10[:] = infile.variables['V10'][:]
    SNOWH[:] = infile.variables['SNOWH'][:]
    SEAICE[:] = infile.variables['SEAICE'][:]
    RAINNC[:] = infile.variables['RAINNC'][:]
    SNOWNC[:] = infile.variables['SNOWNC'][:]
    #REFL_10CM[:] = infile.variables['REFL_10CM'][ :, 0, :, : ]
    
    infile.close()
    outfile.close()
def wrf_pressure( P, PB ):
    """ Calculate the pressure in hPa given Perturbation Pressure (P) and Base State Pressure (PB) in Pa.
    ---------------------
    P (numpy.ndarray): ndarray of Perturbation Pressure from WRF 
    PB (numpy.ndarray): ndarray of Base State Pressure from WRF
    ---------------------
    returns:
        numpy.ndarray of pressure in hPa the same shape as P and PB 
        """
    assert P.shape == PB.shape, 'Arrays are different shapes. They must be the same shape.'
    return ( P + PB ) * .01

def wrf_height( PH, PHB ):
    """ Calculate the geopotential height given the Perturbation Geopotential (PH) and the Base State
        Geopotential (PHB). PH and PHB must be the same shape. 
    ---------------------
    PH (numpy.ndarray): ndarray of Perturbation Geopotential from WRF
    PHB (numpy.ndarray): ndarray of Base State Geopotential from WRF
    ---------------------
    returns:
        numpy.ndarray of geopotential height in meters in sane shape as PH and PHB
        """
    assert PH.shape == PHB.shape, 'Arrays are different shapes. They must be the same shape.'
    return ( PH + PHB ) / 9.81

def wrf_theta( T ):
    """ Calculate the potential temperature given the Perturbation Potential Temperature (T) in degrees Kelvin.
    ---------------------
    T (numpy.ndarray): ndarray of Perturbation Potential Temperature from WRF
    ---------------------
    returns:
        numpy.ndarray of potential temperature in degrees Kelvin, same shape as T
        """
    return T + 300

def wrf_temp( THETA, PRES ):
    """ Calculate the 'normal' temperature in degrees Kelvin given the 
    Potential Temperature (THETA in Kelvin) and Pressure (PRES in hPa or mb).
    PRES and THETA must be the same shape.
    ---------------------
    THETA (numpy.ndarray): ndarray of potential temperature in degrees Kelvin
    PRES (numpy.ndarray): ndarray of pressure in hPa or mb same shape as THETA
    ---------------------
    returns:
        numpy.ndarray of 'normal' temperature in degrees Kelvin same shape as THETA and PRES
        """
    assert THETA.shape == PRES.shape, 'Arrays are different shapes. They must be the same shape.'
    K = 0.2854
    return THETA * ( PRES / 1000 )**K

def wrf_rh( TEMP, PRES, QVAPOR ):
    """ Calculate relative humidity given the Temperature in Kelvin (TEMP), 
    Pressure in hPa or mb (PRES), and Water Vapor Mixing Ratio (QVAPOR).
    TEMP, PRES, and QVAPOR must be the same shape.
    ---------------------
    TEMP (numpy.ndarray): ndarray of 'normal' temperature in degrees Kelvin
    PRES (numpy.ndarray): ndarray of pressure in hPa or mb same shape as TEMP
    QVAPOR (numpy.ndarray): ndarray of Water Vapor Mixing Ratio from WRF, same shape as PRES and TEMP
    ---------------------
    returns:
        numpy.ndarray of relative humidity values (%) same shape as TEMP, PRES, and QVAPOR
        """
    assert TEMP.shape == PRES.shape == QVAPOR.shape, 'Arrays are different shapes. They must be the same shape.'
    e_s = ( 2.53 * pow( 10, 11 ) ) * numpy.exp( -5420/TEMP )
    w_s = 6.22 * ( e_s / PRES )
    return ( QVAPOR * 1000 / w_s ) * 100

def wrf_vort( U, V, dx ):
    """Calculate the relative vorticity given the U and V vector components in m/s
    and the grid spacing dx in meters.
    U and V must be the same shape.
    ---------------------
    U (numpy.ndarray): ndarray of U vector values in m/s
    V (numpy.ndarray): ndarray of V vector values in m/s
    dx (float or int): float or integer of U and V grispacing in meters
    ---------------------
    returns:
        numpy.ndarray of vorticity values s^-1 same shape as U and V
    """
    assert U.shape == V.shape, 'Arrays are different shapes. They must be the same shape.'
    dy = dx
    du = numpy.gradient( U )
    dv = numpy.gradient( V )
    return ( dv[-1]/dx - du[-2]/dy )

def wrf_absvort( U, V, F, dx ):
    """Calculate the absolute vorticity given the U and V vector components in m/s,
    the Coriolis sine latitude term (F) in s^-1, and gridspacing dx in meters. U, V, and F
    must be the same shape.
    ---------------------
    U (numpy.ndarray): ndarray of U vector values in m/s
    V (numpy.ndarray): ndarray of V vector values in m/s
    F (numpy.ndarray): ndarray of Coriolis sine latitude values in s^-1
    dx (float or int): float or integer of U and V grispacing in meters
    ---------------------
    returns:
        numpy.ndarray of absolute vorticity values s^-1 same shape as U and V
    
    """
    assert U.shape == V.shape == F.shape, 'Arrays are different shapes. They must be the same shape.'
    return wrf_vort( U, V, dx ) + F

def wrf_pv( U, V, F, THETA, PRES, MAPFAC_M, dx ):
    """Calculate the potential vorticity given the U and V vector components in m/s,
    the Coriolis sine latitude term (F) in s^-1, THETA potential temperature in degrees
    Kelvin, PRES pressure in hPa or mb, the map scale factor on mass grid and the gridspacing 
    dx in meters. U, V, F, THETA, and PRES must be 4D arrays.
    ---------------------
    U (numpy.ndarray): ndarray of U vector values in m/s
    V (numpy.ndarray): ndarray of V vector values in m/s
    F (numpy.ndarray): ndarray of Coriolis sine latitude values in s^-1
    THETA (numpy.ndarray): ndarray of potential temperature in degrees Kelvin
    PRES (numpy.ndarray): ndarray of pressure in hPa or mb same shape as THETA
    MAPFAC_M (numpy.ndarray): 2D of map scale factor on mass grid. 
    dx (float or int): float or integer of U and V grispacing in meters
    ---------------------
    returns:
        numpy.ndarray of potential vorticity values in ( K * m^2 * kg^-1 * s^-1 ) * 10^6 
        ( or 1 PVU * 10^6).
    """
    assert U.shape == V.shape == F.shape == THETA.shape == PRES.shape, 'Arrays are different shapes. They must be the same shape.'
    ## pres in hPa needs to convert to Pa
    PRES = PRES * 100
    dx = dx * MAPFAC_M
    dy = dx
    grav = 9.8

    dVt,dVp,dVy,dVx = numpy.gradient( V )
    dUt,dUp,dUy,dUx = numpy.gradient( U )
    dTt,dTp,dTy,dTx = numpy.gradient( THETA )
    dPt,dp,dPy,dPx = numpy.gradient( PRES )
    return ( -grav * ( -dVp/dp * dTx/dx + dUp/dp * dTy/dy + ( dVx/dx - dUy/dy + F ) * dTp/dp ) ) * pow(10, 6)
    
    

def wrf_to_pv( grid, surface, interplevels ):
    """Linearly interpolates a grid to potential voericity levels defined by 
        interplevels. This uses the numpy.interp function that requires
        the function coordinate values to be monotonically increasing,
        and this function does not check for this. See numpy.interp docs
        for more information. Grid and surface must be the same shape.
    -----------------------
    grid (numpy.ndarray): 4D array of values to be interpolated onto a vertical surface.  Must be monotonically increasing.
    surface (numpy.ndarray): 4D array of the verical coordinate values of the surface to be interpolated to. Must be monotonically increasing.
    interplevels (numpy.ndarray): 1D array of vertical coordinate values that are desired to be interpolated to.
    -----------------------
    returns:
        numpy.ndarray of values of shape (grid.shape[0], len( interplevels ), grid.shape[2], grid.shape[3] )
    """
    shape = ( grid.shape[0], interplevels.shape[0], grid.shape[2], grid.shape[3] )
    outgrid = numpy.ones( shape )
    for time in numpy.arange( grid.shape[0] ):
        for idx, val in numpy.ndenumerate( grid[0][0] ):
            column = surface[ time, :, idx[0], idx[1] ]
            column_GRID = grid[ time, :, idx[0], idx[1] ]
            pvidx = numpy.where( ( column < interplevels[0] ) & ( 0 < column ) )[0]
            if len( pvidx ) == 0:
                pvidx = numpy.where( ( column < interplevels[1] ) & ( 0 < column ) ) [0]
            else:
                pvidx = pvidx
            try:
                value = numpy.interp( interplevels, column[ pvidx[-1] : ], column_GRID[ pvidx[-1] : ], left=numpy.nan, right=numpy.nan )
                outgrid[ time, :, idx[0], idx[1] ] = value[:]
            except:
                value = numpy.interp( interplevels, column, column_GRID, left=numpy.nan, right=numpy.nan )
                outgrid[ time, :, idx[0], idx[1] ] = value[:]
    return outgrid

# <codecell>

#dview.execute('import numpy')
#dview['Dataset'] = Dataset
#dview['wrf_to_pv'] = wrf_to_pv

## open the files
#infile = Dataset('wrf_pres_out.nc')
#outfile = Dataset( 'wrf_pv_out.nc', 'w', format='NETCDF4' )
infilename = 'wrf_pres_out.nc'
outfilename = 'wrf_pv_out.nc'

fname_done = 'pv_done'
cmd1 = 'rm -f ' + fname_done
os.system(cmd1)

## copy global attributes an initialize an output file
print 'Copying global attributes'
wrf_copy_attributes( infilename, outfilename )

## copy the static fields to the output file
print 'Copying static fields'
wrf_copy_static_fields( infilename, outfilename )

## copy over the surface fields
print 'Copying surface fields'
wrf_copy_sfc_fields( infilename, outfilename )

infile = Dataset( infilename )
outfile = Dataset( outfilename, 'a')

print 'Reading in 3D fields'
HGHT_in = infile.variables['HEIGHT'][:]
PRES_in = infile.variables['PRES'][:]

THETA_in = infile.variables['THETA'][:]

QVAPOR_in = infile.variables['QVAPOR'][:]
QCLOUD_in = infile.variables['QCLOUD'][:]
QICE_in = infile.variables['QICE'][:]

U_in = infile.variables['U'][:]
V_in = infile.variables['V'][:]
W_in = infile.variables['W'][:]
#F_in = infile.variables['F'][:]


SLP = outfile.createVariable( 'SLP', 'f4', ('time', 'south_north', 'west_east') )
SLP[:] = infile.variables['SLP'][:]

gc.collect()
T_TEND_DYNAMICS_in = infile.variables['T_TEND_DYNAMICS'][:]
T_TEND_PHYSICS_in = infile.variables['T_TEND_PHYSICS'][:]

#RU_TENDF_in = infile.variables['RU_TENDF'][:]
#RV_TENDF_in = infile.variables['RV_TENDF'][:]
#RW_TENDF_in = infile.variables['RW_TENDF'][:]
T_TEND_SMALL_in = infile.variables['T_TEND_SMALL'][:]
T_TEND_LARGE_in = infile.variables['T_TEND_LARGE'][:]

RTHRATEN_in = infile.variables['RTHRATEN'][:]
RTHRATLW_in = infile.variables['RTHRATLW'][:]
RTHRATSW_in = infile.variables['RTHRATSW'][:]
RTHBLTEN_in = infile.variables['RTHBLTEN'][:]
RTHCUTEN_in = infile.variables['RTHCUTEN'][:]
RTHMIX_in = infile.variables['RTHMIX'][:]

H_DIABATIC_in = infile.variables['H_DIABATIC'][:]
PV_in = infile.variables['PV'][:]
PVDIAB_in = infile.variables['PVDIAB'][:]
PVRATEN_in = infile.variables['PVRATEN'][:]
PVRATLW_in = infile.variables['PVRATLW'][:]
PVRATSW_in = infile.variables['PVRATSW'][:]
PVBLTEN_in = infile.variables['PVBLTEN'][:]
PVCUTEN_in = infile.variables['PVCUTEN'][:]
PVCOND_in = infile.variables['PVCOND'][:]
PVMIX_in = infile.variables['PVMIX'][:]

gc.collect()

names = [ 'PV', 'HEIGHT', 'PRES', 'THETA', 'U', 'V', 'W','QVAPOR', 'QCLOUD', 'QICE',
         'T_TEND_DYNAMICS', 'T_TEND_PHYSICS','T_TEND_SMALL', 'T_TEND_LARGE', 'RTHRATEN', 'RTHRATLW', 'RTHRATSW', 
         'RTHBLTEN', 'RTHCUTEN', 'H_DIABATIC', 'RTHMIX','PVDIAB', 
	 'PVRATEN', 'PVRATLW', 'PVRATSW', 'PVCOND', 'PVBLTEN','PVCUTEN','PVMIX']
grids = [ PV_in, HGHT_in, PRES_in, THETA_in, U_in, V_in, W_in, QVAPOR_in, QCLOUD_in, QICE_in,
         T_TEND_DYNAMICS_in, T_TEND_PHYSICS_in, T_TEND_SMALL_in, T_TEND_LARGE_in,RTHRATEN_in, RTHRATLW_in, RTHRATSW_in, 
	 RTHBLTEN_in, RTHCUTEN_in, H_DIABATIC_in, RTHMIX_in, PVDIAB_in, 
	 PVRATEN_in, PVRATLW_in, PVRATSW_in, PVCOND_in, PVBLTEN_in, PVCUTEN_in, PVMIX_in ]

#names2 = ['RTHBLTEN', 'RTHCUTEN', 'H_DIABATIC', 'RTHMIX','PVDIAB', 
#	 'PVRATEN', 'PVRATLW', 'PVRATSW', 'PVCOND', 'PVBLTEN','PVCUTEN','PVMIX']
#grids2 = [RTHBLTEN_in, RTHCUTEN_in, H_DIABATIC_in, RTHMIX_in, PVDIAB_in, 
#	 PVRATEN_in, PVRATLW_in, PVRATSW_in, PVCOND_in, PVBLTEN_in, PVCUTEN_in, PVMIX_in ]


surface=PV_in
#mydict=dict( surface=PV_in )
#dview.push( mydict )

#interplevels = numpy.asarray( [ 1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0 ] )
interplevels = numpy.asarray( [ 0.5,0.7,0.9,1.2,1.5,2.0,2.7,3.6,4.8,6.4,8.5,11.4,15.1,20.0,26.6,35.3,47.0,62.3,82.8,110.0 ])

#mydict=dict( interplevels=interplevels )
#dview.push( mydict )

print 'interpolating'
func = lambda grid : wrf_to_pv( grid, surface, interplevels )
start_time = datetime.datetime.now()
#data = lview.map( func, grids, block=True )
data = map( func, grids )
end_time = datetime.datetime.now()
secs = end_time - start_time
print 'interpolation took', secs
gc.collect()
count = 0
for name in names:
    grid = outfile.createVariable( name, 'f4', ('time', 'bottom_top', 'south_north', 'west_east') )
    grid[:] = data[ count ][:]
    count += 1

infile.close() 
outfile.close()

with open(fname_done, 'a'):
    os.utime(fname_done, None)
