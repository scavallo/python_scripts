# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

#from IPython.parallel import Client
import numpy
from netCDF4 import Dataset
import datetime
import os
import gc
import subprocess

#rc = Client()
#lview = rc.load_balanced_view()
#dview = rc[:]
#lview.block = True
#dview.block = True
#gc.collect()

# <codecell>

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
    REFL_10CM = outfile.createVariable( 'REFL_10CM', 'f4', ('time', 'south_north', 'west_east') )
    
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
    REFL_10CM[:] = infile.variables['REFL_10CM'][ :, 0, :, : ]
    
    infile.close()
    outfile.close()

def wrf_slp( infilename, outfilename ):
    
    infile = Dataset( infilename )
    outfile = Dataset( outfilename, 'a' )
    
    stemps = infile.variables['T2'][:]+6.5*infile.variables['HGT'][:]/1000.
    SLP_in = infile.variables['PSFC'][:]*numpy.exp(9.81/(287.0*stemps)*infile.variables['HGT'][:])*0.01 + (6.7 * infile.variables['HGT'][:] / 1000)
    SLP = outfile.createVariable( 'SLP', 'f4', ('time', 'south_north', 'west_east') )
    SLP[:] = SLP_in    
    
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
    
    #XLONG = outfile.createVariable( 'XLONG', 'f4', ('south_north', 'west_east') )
    #XLAT = outfile.createVariable( 'XLAT', 'f4', ('south_north', 'west_east') )

    #SINALPHA = outfile.createVariable( 'SINALPHA', 'f4', ('south_north', 'west_east') )
    #COSALPHA = outfile.createVariable( 'COSALPHA', 'f4', ('south_north', 'west_east') )
    #MAPFAC_M = outfile.createVariable( 'MAPFAC_M', 'f4', ('south_north', 'west_east') )
    F = outfile.createVariable( 'F', 'f4', ('south_north', 'west_east') )
    XTIME = outfile.createVariable( 'XTIME', 'f4', ('time') )

    #XLONG[:] = infile.variables['XLONG'][0][:]
    #XLAT[:] = infile.variables['XLAT'][0][:]
    F[:] = infile.variables['F'][0][:]
    MU = infile.variables['MU'][0][:]
    MUB = infile.variables['MUB'][0][:]

    #MAPFAC_M[:] = infile.variables['MAPFAC_M'][0][:]
    #SINALPHA[:] = infile.variables['SINALPHA'][0][:]
    #COSALPHA[:] = infile.variables['COSALPHA'][0][:]
    XTIME[:] = infile.variables['XTIME'][:]
    infile.close()
    outfile.close()


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
    e_0 = 6.1173 ## mb
    t_0 = 273.16 ## K
    Rv = 461.50 ## J K-1 Kg-1
    Lv_0 = 2.501 * 10**6 ## J Kg-1
    K1 = Lv_0 / Rv ## K 
    K2 = 1 / t_0 ## K-1
    K3 = 1 / TEMP ## K-1
    ## Clausius Clapeyron Equation
    e_s = e_0 * numpy.exp( K1 * ( K2 - K3 ) )
    w_s = ( 0.622 * e_s ) / ( PRES - e_s )
    return ( QVAPOR / w_s ) * 100

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
    assert U.shape == V.shape, 'Arrays are different shapes. They must be the same shape.'
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
    assert U.shape == V.shape == THETA.shape == PRES.shape, 'Arrays are different shapes. They must be the same shape.'
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

def wrf_to_pres( grid, surface, interplevels ):
    shape = ( grid.shape[0], interplevels.shape[0], grid.shape[2], grid.shape[3] )
    outgrid = numpy.ones( shape )
    for time in numpy.arange( grid.shape[0] ):
        for idx, val in numpy.ndenumerate( grid[0][0] ):
            column = surface[ time, :, idx[0], idx[1] ]
            column_GRID = grid[ time, :, idx[0], idx[1] ]
            
            value = numpy.interp( interplevels, column, column_GRID, left=numpy.nan, right=numpy.nan )
            outgrid[ time, :, idx[0], idx[1] ] = value[:]
    return outgrid

# <codecell>

#dview.execute('import numpy')
#dview['Dataset'] = Dataset
#dview['wrf_to_pres'] = wrf_to_pres

print 'Starting!'

fname_done = 'pres_done'
cmd1 = 'rm -f ' + fname_done
os.system(cmd1)

## define file names
infilename = 'wrf_eta_in.nc'
outfilename = 'wrf_pres_out.nc'

## copy global attributes an initialize an output file
wrf_copy_attributes( infilename, outfilename )

## copy the static fields to the output file
wrf_copy_static_fields( infilename, outfilename )

## copy over the surface fields
wrf_copy_sfc_fields( infilename, outfilename )

## compute sea level pressure and write
print 'Computing SLP'
wrf_slp( infilename, outfilename )

gc.collect()

infile = Dataset( infilename )
outfile = Dataset( outfilename, 'a' )

dx = float(infile.DX)
print 'dx is ', dx

print 'I am now moving on to metr fields'
PRES_in = wrf_pressure( infile.variables['P'][:], infile.variables['PB'][:] )
print 'One'
HGHT_in = wrf_unstagger( wrf_height( infile.variables['PH'][:], infile.variables['PHB'][:] ), 'Z' )
print 'Two'
THETA_in = wrf_theta( infile.variables['T'][:] )
TEMP_in = wrf_temp( THETA_in, PRES_in )
RH_in = wrf_rh( TEMP_in, PRES_in, infile.variables['QVAPOR'][:] )
print 'Three'
U_in = wrf_unstagger( infile.variables['U'][:], 'U' )
V_in = wrf_unstagger( infile.variables['V'][:], 'V' )
W_in = wrf_unstagger( infile.variables['W'][:], 'W' )


print 'now comes the computational goodies'
VORT_in = wrf_vort( U_in, V_in, dx )
ABSVORT_in = wrf_absvort( U_in, V_in, infile.variables['F'][0][:], dx )
PV_in = wrf_pv( U_in, V_in, infile.variables['F'][0][:], THETA_in, PRES_in, infile.variables['MAPFAC_M'][0][:], dx )
gc.collect()


QVAPOR_in = infile.variables['QVAPOR'][:]
QCLOUD_in = infile.variables['QCLOUD'][:]
QICE_in = infile.variables['QICE'][:]
QRAIN_in = infile.variables['QRAIN'][:]
QSNOW_in = infile.variables['QSNOW'][:]

print 'Diabatic Tendencies'
T_TEND_in = infile.variables['T_TEND'][:] / ( infile.variables['MU'][0][:] + infile.variables['MUB'][0][:] )
T_TENDF_in = infile.variables['T_TENDF'][:] / ( infile.variables['MU'][0][:] + infile.variables['MUB'][0][:] )
T_TEND_SMALL_in = infile.variables['T_TEND_SMALL'][:]
T_TEND_LARGE_in = infile.variables['T_TEND_LARGE'][:]


#RU_TENDF_in = wrf_unstagger( infile.variables['RU_TENDF'][:], 'U' )
#RV_TENDF_in = wrf_unstagger( infile.variables['RV_TENDF'][:], 'V' )
#RW_TENDF_in = wrf_unstagger( infile.variables['RW_TENDF'][:], 'Z' )

gc.collect()

RTHRATEN_in = infile.variables['RTHRATEN'][:] / ( infile.variables['MU'][0][:] + infile.variables['MUB'][0][:] )
RTHRATLW_in = infile.variables['RTHRATLW'][:]
RTHRATSW_in = infile.variables['RTHRATSW'][:]
RTHBLTEN_in = infile.variables['RTHBLTEN'][:] / ( infile.variables['MU'][0][:] + infile.variables['MUB'][0][:] )
RTHCUTEN_in = infile.variables['RTHCUTEN'][:] / ( infile.variables['MU'][0][:] + infile.variables['MUB'][0][:] )
H_DIABATIC_in = infile.variables['H_DIABATIC'][:]

# Total physics tendency
T_TEND_PHYSICS_in = T_TENDF_in + H_DIABATIC_in
T_TEND_DYNAMICS_in = T_TEND_SMALL_in + T_TEND_LARGE_in

# Unaccounted for tendency
RTHMIX_in = T_TEND_PHYSICS_in - (RTHRATEN_in + RTHCUTEN_in + RTHBLTEN_in + H_DIABATIC_in)
gc.collect()
## convert from PVU/s to PVU/day by * 86400
PVRATEN_in = wrf_pv( U_in, V_in, infile.variables['F'][0][:], RTHRATEN_in, PRES_in, infile.variables['MAPFAC_M'][0][:], dx ) * 86400
gc.collect()
PVRATLW_in = wrf_pv( U_in, V_in, infile.variables['F'][0][:], RTHRATLW_in, PRES_in, infile.variables['MAPFAC_M'][0][:], dx ) * 86400
gc.collect()
PVRATSW_in = wrf_pv( U_in, V_in, infile.variables['F'][0][:], RTHRATSW_in, PRES_in, infile.variables['MAPFAC_M'][0][:], dx ) * 86400
gc.collect()
PVCOND_in = wrf_pv( U_in, V_in, infile.variables['F'][0][:], H_DIABATIC_in, PRES_in, infile.variables['MAPFAC_M'][0][:], dx ) * 86400
gc.collect()
PVBLTEN_in = wrf_pv( U_in, V_in, infile.variables['F'][0][:], RTHBLTEN_in, PRES_in, infile.variables['MAPFAC_M'][0][:], dx ) * 86400
gc.collect()
PVCUTEN_in = wrf_pv( U_in, V_in, infile.variables['F'][0][:], RTHCUTEN_in, PRES_in, infile.variables['MAPFAC_M'][0][:], dx ) * 86400
gc.collect()
PVDIAB_in = wrf_pv( U_in, V_in, infile.variables['F'][0][:], T_TEND_PHYSICS_in, PRES_in, infile.variables['MAPFAC_M'][0][:], dx ) * 86400
gc.collect()
PVMIX_in = wrf_pv( U_in, V_in, infile.variables['F'][0][:], RTHMIX_in, PRES_in, infile.variables['MAPFAC_M'][0][:], dx ) * 86400
gc.collect()


#PV_DIAB_in = PVRATEN_in + PVCOND_in + PVBLTEN_in + PVCUTEN_in
print 'All done!'
infile.close()
names = [ 'HEIGHT', 'PRES', 'THETA', 'TEMP', 'QVAPOR', 'QCLOUD', 'QICE', 'QRAIN', 'QSNOW', 'ABSVORT', 'PV', 'RELH', 'U', 'V', 'W', 
         'T_TEND_DYNAMICS', 'T_TEND_PHYSICS','T_TEND_SMALL', 'T_TEND_LARGE']
grids = [ HGHT_in, PRES_in, THETA_in, TEMP_in, QVAPOR_in, QCLOUD_in, QICE_in, QRAIN_in, QSNOW_in, ABSVORT_in, PV_in, RH_in, U_in, V_in, W_in,
         T_TEND_DYNAMICS_in, T_TEND_PHYSICS_in, T_TEND_SMALL_in, T_TEND_LARGE_in]
surface=PRES_in


names2 = ['RTHRATEN', 'RTHRATLW', 'RTHRATSW', 
          'RTHBLTEN', 'RTHCUTEN', 'H_DIABATIC', 'RTHMIX','PVDIAB', 
	  'PVRATEN', 'PVRATLW', 'PVRATSW', 'PVCOND', 'PVBLTEN','PVCUTEN','PVMIX']
grids2 = [ RTHRATEN_in, RTHRATLW_in, RTHRATSW_in, 
           RTHBLTEN_in, RTHCUTEN_in, H_DIABATIC_in, RTHMIX_in, PVDIAB_in,
	   PVRATEN_in, PVRATLW_in, PVRATSW_in, PVCOND_in, PVBLTEN_in, PVCUTEN_in, PVMIX_in ]


#mydict=dict( surface=PRES_in )
#dview.push( mydict )

interplevels = numpy.asarray( [ 1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100 ,70, 50, 30 ] )
#mydict=dict( interplevels=interplevels )
#dview.push( mydict )


## most re-order the arrays becuase the interpolation expects values to be
## monotonically increasing
#def func( grid, name ):
#    interp = wrf_to_pres( grid[:, ::-1, :, :], surface[:, ::-1, :, :], interplevels )
#    newgrid = outfile.createVariable( name, 'f4', ('time', 'bottom_top', 'south_north', 'west_east') )
#    newgrid[:] = interp[:]
func = lambda grid : wrf_to_pres( grid[:, ::-1, :, :], surface[:, ::-1, :, :], interplevels )
start_time = datetime.datetime.now()
print "Starting first round of interpolation..."
#data = lview.map( func, grids, block=True )
data = map( func, grids )
end_time = datetime.datetime.now()
secs = end_time - start_time
print 'interpolation took', secs

gc.collect()

count = 0
for name in names:
    print 'Writing ', name
    grid = outfile.createVariable( name, 'f4', ('time', 'bottom_top', 'south_north', 'west_east') )
    grid[:] = data[ count ][:]
    count += 1 

del data, HGHT_in, PRES_in, THETA_in, TEMP_in, QVAPOR_in, QCLOUD_in, QICE_in, QRAIN_in, QSNOW_in, ABSVORT_in, PV_in, RH_in, U_in, V_in, W_in, T_TEND_DYNAMICS_in, T_TEND_PHYSICS_in, T_TEND_SMALL_in, T_TEND_LARGE_in 
gc.collect()

start_time = datetime.datetime.now()
print "Starting second round of interpolation..."
#data = lview.map( func, grids, block=True )
data2 = map( func, grids2 )
end_time = datetime.datetime.now()
secs = end_time - start_time
print 'interpolation took', secs



gc.collect()
count = 0
for name in names2:
    print 'Writing ', name
    grid = outfile.createVariable( name, 'f4', ('time', 'bottom_top', 'south_north', 'west_east') )
    grid[:] = data2[ count ][:]
    count += 1     
outfile.close()


with open(fname_done, 'a'):
    os.utime(fname_done, None)
