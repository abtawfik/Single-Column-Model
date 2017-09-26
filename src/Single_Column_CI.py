#!/bin/python

### Handling Arrays and Reading Grib package
import time
import numpy as np
from netCDF4 import Dataset

### Plotting packages
import sys
import math
from HandyMet import evaluate_ci_plevels



#################################################
#####
#####  Directories and time variables
#####
#################################################
directory   =  '/glade/scratch/abtawfik/NARR_Data/'
iyear       =  sys.argv[1]
imon        =  sys.argv[2]
#iyear       =  '2002'
#imon        =  '06'
print iyear, imon
print type(iyear),type(imon)
#sys.exit()


def create_grid(inputV, gridSize):
    return  np.repeat(inputV[:,np.newaxis], gridSize, 1).ravel()


def read_with_UTC( inputFile, varname, offset, myaxis, dimensions ):
    return np.roll(inputFile.variables[varname][:], offset, axis=myaxis).reshape(dimensions)


def replace_missing( input_variable, badvalue, new_missing ):
    output_variable  =  input_variable
    output_variable[ (output_variable > badvalue) | (output_variable < -badvalue) ]  =  new_missing
    return output_variable



#################################################
#####
#####  Read in Temperature and Pressure Level
#####
#################################################
input_file  =  Dataset(directory + "air." + iyear + imon + ".nc" , "r")
timestamp   =  input_file.variables['time'][:]
lats        =  input_file.variables['lat'][:]
lons        =  input_file.variables['lon'][:]
pres        =  input_file.variables['level'][:]
pres        =  pres * 1e2

badvalue    =  1e10
missing     =  -99999.
ntim        =  timestamp.shape[0]
nlev        =  pres.shape[0]
nlat        =  lats.shape[0]
nlon        =  lats.shape[1]
nhr         =  8
ndays       =  ntim/nhr
days        =  np.arange(0,ndays)

temp        =  replace_missing( read_with_UTC(input_file,'air',-2,0,(ndays,nhr,nlev,nlat,nlon)), badvalue, missing )
print "  Time   ===>>    ",ntim
print "  Lats   ===>>    ",nlat
print "  Lons   ===>>    ",nlon
print "  Levs   ===>>    ",nlev



#################################################
#####
#####  Read in 2-m Temperature
#####  and latitude and longitude grid from T
#####
#################################################
input_file  =  Dataset(directory + "air.2m." + iyear + imon + ".nc" , "r")
t2m         =  replace_missing( read_with_UTC( input_file, 'air', -2, 0, (ndays,nhr,nlat,nlon)), badvalue, missing )


#################################################
#####
#####  Read in 2-m Specific Humidity
#####
#################################################
input_file  =  Dataset(directory + "shum.2m." + iyear + imon + ".nc" , "r")
q2m         =  replace_missing( read_with_UTC( input_file, 'shum', -2, 0, (ndays,nhr,nlat,nlon)), badvalue, missing )


#################################################
#####
#####  Read in surface pressure [Pa]
#####
#################################################
input_file  =  Dataset(directory + "lhtfl." + iyear + imon + ".nc" , "r")
rnet        =  replace_missing( read_with_UTC( input_file, 'lhtfl', -2, 0, (ndays,nhr,nlat,nlon)), badvalue, missing )

input_file  =  Dataset(directory + "shtfl." + iyear + imon + ".nc" , "r")
rnet        =  replace_missing( rnet + read_with_UTC( input_file, 'shtfl', -2, 0, (ndays,nhr,nlat,nlon)), badvalue, missing )
rnet        =  rnet * -1.0


#################################################
#####
#####  Read in surface pressure [Pa]
#####
#################################################
input_file  =  Dataset(directory + "pres.sfc." + iyear + imon + ".nc" , "r")
psfc        =  replace_missing( read_with_UTC( input_file, 'pres', -2, 0, (ndays,nhr,nlat,nlon)), badvalue, missing )



#################################################
#####
#####  Read in Specific Humidity
#####
#################################################
input_file  =  Dataset(directory + "shum." + iyear + imon + ".nc" , "r")
shum        =  replace_missing( read_with_UTC( input_file, 'shum', -2, 0, (ndays,nhr,nlev,nlat,nlon)), badvalue, missing )



#################################################
#####
#####  Read in Vertical Velocity [Pa/s]
#####
#################################################
Omega_On  =  True
if Omega_On:
    input_file  =  Dataset(directory + "omega." + iyear + imon + ".nc" , "r")
    omega       =  replace_missing( read_with_UTC( input_file, 'omega', -2, 0, (ndays,nhr,nlev,nlat,nlon)), badvalue, missing )
else:
    omega = np.zeros( shum.shape )





#################################################
#####
#####  Assign some variables
#####
#################################################
ef    =  np.linspace(0.05,0.95,10)
dt    =  600.                        # 10-minute intervals
min5  =  np.arange(0,21,10.0/60.0)




#################################################
#####
#####  Read in Polynomial fit of net radiation
#####
#################################################
timeOfday  =  np.arange(0,23,3) * 1.0
NETRAD     =  np.zeros( (ndays,min5.shape[0],nlat,nlon) )
for dd in np.arange(0,ndays):
    for yy in np.arange(0,nlat):
        for xx in np.arange(0,nlon):
            polyno              =  np.poly1d( np.polyfit(timeOfday,rnet[dd,:,yy,xx],6) )
            NETRAD[dd,:,yy,xx]  =  polyno(min5)


#inputFile  =  Dataset('NETRAD.'+iyear+imon+'.nc', 'r')
#NETRAD     =  inputFile.variables['NETRAD'][:]



#################################################
#####
#####  Calculate the LoCo CI-Flux equation terms
#####
#################################################
tbm      =  np.zeros( (ndays,ef.shape[0],nlat,nlon) )
bclp     =  np.zeros( (ndays,ef.shape[0],nlat,nlon) )
timeofci =  np.zeros( (ndays,ef.shape[0],nlat,nlon) )
qbcl     =  np.zeros( (ndays,ef.shape[0],nlat,nlon) )

#d0       =  10
#d1       =  11
#i0       =  4
#i1       =  5
#i2       =  285
#i3       =  286



print("Start the CI evaluation")
start = time.time()
tbm,bclp,timeofci,qbcl  =  evaluate_ci_plevels(temp[:,0,:,:,:], shum[:,0,:,:,:], pres, omega, 1, t2m[:,0,:,:], q2m[:,0,:,:], psfc[:,0,:,:], NETRAD, ef, dt, missing)
#tbm[d0:d1,:,i0:i1,i2:i3],bclp[d0:d1,:,i0:i1,i2:i3],timeofci[d0:d1,:,i0:i1,i2:i3], qbcl[d0:d1,:,i0:i1,i2:i3]  =  \
#evaluate_ci_plevels(temp[d0:d1,0,:,i0:i1,i2:i3], shum[d0:d1,0,:,i0:i1,i2:i3], pres, 1, t2m[d0:d1,0,i0:i1,i2:i3], q2m[d0:d1,0,i0:i1,i2:i3], psfc[d0:d1,0,i0:i1,i2:i3], NETRAD[d0:d1,:,i0:i1,i2:i3], ef, dt, missing)
#tbm[d0:d1,:,i0:i1,i2:i3],bclp[d0:d1,:,i0:i1,i2:i3],timeofci[d0:d1,:,i0:i1,i2:i3], qbcl[d0:d1,:,i0:i1,i2:i3]  =  \
#ECI77(temp[d0:d1,0,:,i0:i1,i2:i3], shum[d0:d1,0,:,i0:i1,i2:i3], pres, 1, t2m[d0:d1,0,i0:i1,i2:i3], q2m[d0:d1,0,i0:i1,i2:i3], psfc[d0:d1,0,i0:i1,i2:i3], NETRAD[d0:d1,:,i0:i1,i2:i3], ef, dt, missing)
end = time.time()
print(end - start)



################################################################################
# Open a new NetCDF file to write the data to. For format, you can choose from #
################################################################################
if Omega_On:
    ncid = Dataset('Single_Column_CI.OMEGA.'+iyear+imon+'.nc', 'w')
else:
    ncid = Dataset('Single_Column_CI.no_Omega.'+iyear+imon+'.nc', 'w')
    
#####################################################
# Assign the dimension data to the new NetCDF file. #
#####################################################
ncid.createDimension('time', None )
ncid.createDimension('day' , ndays)
ncid.createDimension('ef'  , ef.shape[0])
ncid.createDimension('y'   , nlat )
ncid.createDimension('x'   , nlon )

#####################################################
# Assign the dimension data to the new NetCDF file.
#####################################################
ncid.createVariable('time'  , np.float64, ('time' ))[:]  =  timestamp
ncid.createVariable('day'   , np.float64, ('day'  ))[:]  =  days
ncid.createVariable('ef'    , np.float64, ('ef'   ))[:]  =  ef
ncid.createVariable('lat2d' , np.float32, ('y','x'))[:]  =  lats
ncid.createVariable('lon2d' , np.float32, ('y','x'))[:]  =  lons

#####################################################
# Ok, time to create our variable
#####################################################
ncid.createVariable('TBM'      , np.float32, ('day','ef','y','x'))[:]  =  tbm
ncid.createVariable('BCLP'     , np.float32, ('day','ef','y','x'))[:]  =  bclp
ncid.createVariable('QBCL'     , np.float32, ('day','ef','y','x'))[:]  =  qbcl
#ncid.createVariable('CAPE'     , np.float32, ('day','ef','y','x'))[:]  =  cape
ncid.createVariable('TimeOfCI' , np.float32, ('day','ef','y','x'))[:]  =  timeofci

ncid.close()  # close the new file

























exit()

def read_with_UTC( inputFile, varname, offset, myaxis, dimensions ):
    return np.roll(inputFile.variables[varname][:], offset, axis=myaxis).reshape(dimensions)


#################################################
#####
#####  Read in Temperature and Pressure Level
#####
#################################################
input_file  =  Dataset(directory + "air." + iyear + imon + ".nc" , "r")
timestamp   =  input_file.variables['time'][:]
lats        =  input_file.variables['lat'][:]
lons        =  input_file.variables['lon'][:]
pres        =  input_file.variables['level'][:]
pres        =  pres * 1e2

missing     =  9.96921e36
ntim        =  timestamp.shape[0]
nlev        =  pres.shape[0]
nlat        =  lats.shape[0]
nlon        =  lats.shape[1]
nhr         =  8
ndays       =  ntim/nhr
days        =  np.arange(0,ndays)

temp        =  read_with_UTC( input_file, 'air', -2, 0, (ndays,nhr,nlev,nlat,nlon) )
print "  Time   ===>>    ",ntim
print "  Lats   ===>>    ",nlat
print "  Lons   ===>>    ",nlon
print "  Levs   ===>>    ",nlev



#################################################
#####
#####  Read in 2-m Temperature
#####  and latitude and longitude grid from T
#####
#################################################
input_file  =  Dataset(directory + "air.2m." + iyear + imon + ".nc" , "r")
t2m         =  read_with_UTC( input_file, 'air', -2, 0, (ndays,nhr,nlat,nlon) )


#################################################
#####
#####  Read in 2-m Specific Humidity
#####
#################################################
input_file  =  Dataset(directory + "shum.2m." + iyear + imon + ".nc" , "r")
q2m         =  read_with_UTC( input_file, 'shum', -2, 0, (ndays,nhr,nlat,nlon) )
h2m         =  2.0


#################################################
#####
#####  Read in accumulated surface evaporation
#####
#################################################
#input_file  =  Dataset(directory + "evap." + iyear + imon + ".nc" , "r")
#evap        =  input_file.variables['evap'][:]
#evap        =  evap0
#evap[0:ntim-1,:,:]  =  evap0[1:,:,:]
#evap[  ntim-1,:,:]  =  missing
#del evap0


#################################################
#####
#####  Read in accumulated potential evaporation
#####
#################################################
#input_file  =  Dataset(directory + "pevap." + iyear + imon + ".nc" , "r")
#pevap       =  input_file.variables['pevap'][:]
#pevap       =  pevap0
#pevap[0:ntim-1,:,:]  =  pevap0[1:,:,:]
#pevap[  ntim-1,:,:]  =  missing
#del pevap0

#################################################
#####
#####  Read in surface pressure [Pa]
#####
#################################################
input_file  =  Dataset(directory + "lhtfl." + iyear + imon + ".nc" , "r")
rnet        =  read_with_UTC( input_file, 'lhtfl', -2, 0, (ndays,nhr,nlat,nlon) )

input_file  =  Dataset(directory + "shtfl." + iyear + imon + ".nc" , "r")
rnet        =  rnet + read_with_UTC( input_file, 'shtfl', -2, 0, (ndays,nhr,nlat,nlon) )


#################################################
#####
#####  Read in surface pressure [Pa]
#####
#################################################
input_file  =  Dataset(directory + "pres.sfc." + iyear + imon + ".nc" , "r")
psfc        =  read_with_UTC( input_file, 'pres', -2, 0, (ndays,nhr,nlat,nlon) )


#################################################
#####
#####  Read in Specific Humidity
#####
#################################################
input_file  =  Dataset(directory + "shum." + iyear + imon + ".nc" , "r")
shum        =  read_with_UTC( input_file, 'shum', -2, 0, (ndays,nhr,nlev,nlat,nlon) )


#################################################
#####
#####  Read in Height
#####
#################################################
input_file  =  Dataset(directory + "hgt." + iyear + imon + ".nc" , "r")
hgt         =  read_with_UTC( input_file, 'hgt', -2, 0, (ndays,nhr,nlev,nlat,nlon) )



#yy = 106
#xx = 206
#dd = 0
#polyno     =  np.poly1d( np.polyfit(timeOfday,rnet[dd,:,yy,xx],6) )
#NETRAD     =  np.zeros( (ndays,min5.shape[0],nlat,nlon) )
#for dd in np.arange(0,ndays):
#    NETRAD[dd,:,:,:]  =  np.repeat(polyno(min5), nlat*nlon).reshape(min5.shape[0],nlat,nlon)
##exit()
#for dd in np.arange(0,ndays):
#    for yy in np.arange(0,nlat):
#        for xx in np.arange(0,nlon):
#            polyno              =  np.poly1d( np.polyfit(timeOfday,rnet[dd,:,yy,xx],6) )
#            NETRAD[dd,:,yy,xx]  =  polyno(min5)
#ncid = 
#ncid.createDimension('hour126', None )
#ncid.createDimension('day'    , ndays)
#ncid.createDimension('y'      , nlat )
#ncid.createDimension('x'      , nlon )
#ncid.createVariable ('day'     , np.float64, ('day'     ))[:]  =  days
#ncid.createVariable ('hour126' , np.float64, ('hour126' ))[:]  =  min5
#ncid.createVariable ('lat2d'   , np.float32, ('y','x'   ))[:]  =  lats
#ncid.createVariable ('lon2d'   , np.float32, ('y','x'   ))[:]  =  lons
#ncid.createVariable ('NETRAD'  , np.float32, ('day','hour126','y','x'))[:]  =  NETRAD
#ncid.close()  # close the new file

