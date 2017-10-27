#!/bin/python
'''
All the tests in this file are for the modules in Useful_Met_Functions.f90
'''


################################################
### Handling Arrays and Reading Grib package
################################################
import numpy as np
import pandas as pd
import pytest


#############################################################
### These are helper Fortran functions to speed things up
#############################################################
from HandyMet import *



###################################################################################
###
### Tests that need to be added 
###   -- Check PBLH and PBLP
###   -- make sure missing value handling is correct; insert random missing value
###   -- check for negative values maybe??
###
###################################################################################


'''
Test the height above ground calculations
---function name:  calculate_height_above_ground---
'''
def test_height_above_ground():
    missing      =  -99999.0
    temperature  =  [ 287., 289. ]
    pressure     =  [ 100000., 98000. ]
    height       =  calculate_height_above_ground(temperature, pressure, missing)
    assert height[0] == 0
    assert height[1] == pytest.approx(170.2454, 0.001)


'''
Test potential temperature calculations
Do this for 1-D and scalar subroutines from HandyMet
---function name:  potentialTemperature ---
---function name:  potentialTemperature1D---
'''
def test_potential_temperature():
    missing      =  -99999.0
    temperature  =  [   288. ,  288. ,   239.  ]
    pressure     =  [ 100000., 98000.,  20000. ]
    theta        =  potentialtemperature  (temperature, pressure, missing)
    theta1D      =  potentialtemperature1d(temperature, pressure, missing)
    assert theta[0] == pytest.approx(288     , 0.001)
    assert theta[1] == pytest.approx(289.6672, 0.001)
    assert theta[2] == pytest.approx(378.5329, 0.001)

    for zz,tt in enumerate(temperature):
        theta1D      =  potentialtemperature1d(temperature[zz], pressure[zz], missing)
        if zz == 0:  assert theta1D == pytest.approx(288     , 0.001)
        if zz == 1:  assert theta1D == pytest.approx(289.6672, 0.001)
        if zz == 2:  assert theta1D == pytest.approx(378.5329, 0.001)



'''
Test getting temperature from potential temperature
---function name:  calculate_temperature---
'''
def test_temperature():
    missing      =  -99999.0
    theta        =  [   288. , 289.6672 , 378.5329 ]
    pressure     =  [ 100000.,  98000.  ,  20000.  ]
    temperature  =  calculate_temperature(theta, pressure, missing)
    assert temperature[0] == pytest.approx(288 , 0.001)
    assert temperature[1] == pytest.approx(288 , 0.001)
    assert temperature[2] == pytest.approx(239 , 0.001)



'''
Test linear average between multiple levels
---linear_layer_avg---
'''
def test_linear_avg():
    missing          =  -99999.0
    even_values      =  np.array( np.arange(0,22,2) )
    expected_values  =  np.array( np.arange(1,21,2) )
    mid_point        =  linear_layer_avg(even_values, missing)
    nlev             =  even_values.shape[0]
    for ii in range(0,expected_values.shape[0]):
        assert mid_point[ii] == expected_values[ii]
    assert mid_point[nlev-1] == missing



'''
Test the weighted average function in three ways 
Don't perform weighted average for sufficiently coarse vertical profiles
---weighted_avg---
'''
def test_weighted_avg():
    missing       =  -99999.0
    pressure      =  np.array( [100000, 99900, 99700, 99500, 99000] )
    theta         =  np.array( [ 2  ,  4 ,  6 ,  8 ,  10] )
    mid_point     =  linear_layer_avg(theta, missing)
    dp            =  depthpressure(pressure, missing)
    numpy_wgt     =  np.average  (mid_point[:4], weights=dp[:4])
    handy_wgt, p  =  weighted_avg(theta, pressure, missing) 
    assert handy_wgt[0] == pytest.approx(numpy_wgt , 0.001)



'''
Test difference between pressure levels
---depthPressure---
'''
def test_depth_pressure():
    nlev      =  27
    missing   =  -99999.0
    pressure  =  np.linspace(1000,350,nlev)
    dp        =  depthpressure(pressure, missing)
    assert all(dp[0:nlev-1] == pytest.approx(25 , 0.001))



'''
Test column density calculation that uses hydrostatic equation
rho * dz  =  -dp / grav
---columnDensity---
---columnDensityNoMix---
---columndensity_cummulative---
---columndensitynomix_cummulative---
'''
def test_column_density():
    nlev           =  10
    missing        =  -99999.0
    dp             =  np.ones(nlev) * 9.81
    qhum           =  np.ones(nlev)
    density        =  columndensitynomix(       dp, missing )
    ratio_density  =  columndensity     ( qhum, dp, missing )
    assert all(density       == pytest.approx(1.0 , 0.001))
    assert all(ratio_density == pytest.approx(1.0 , 0.001))

    sumdensity               =  columndensitynomix_cummulative(       dp, missing )
    sumdensity_mixing_ratio  =  columndensity_cummulative     ( qhum, dp, missing )
    for zz,rho in enumerate(sumdensity):
        assert rho == zz+1
    for zz,rho in enumerate(sumdensity_mixing_ratio):
        assert rho == zz+1



'''
Test the cummulative sum from Fortran Module
---cummulative sum---
'''
def test_cumsum():
    missing         =  -99999.0
    dummy_data      =  np.linspace(0,10,11)
    cummulative     =  cummulative_sum( dummy_data, missing )
    expected_value  =  np.cumsum(dummy_data)
    for ii,cumsum in enumerate(cummulative):
        assert cumsum == expected_value[ii]



'''
Test Hi-resoluation real-world PBL from
Craig Enhanced Sounding Data
'''
def test_hi_res_PBL():
    # Read in the hi-res data
    test_files = [ './test_data/Hi_res_test_profile.20150711.173000.csv', 
                   './test_data/Hi_res_test_profile.20150711.213000.csv', 
                   './test_data/Hi_res_test_profile.20150826.223000.csv',
                   './test_data/Hi_res_test_profile.csv']
    nprofiles  = len(test_files)
    bounds     = [ (91000,94000), (84000,87000), (83000,86000), (77000,80500) ]
    
    for t in range(0,nprofiles):
        df              =  pd.read_csv(test_files[t])
        pressure        =  df["Pressure"]
        temperature     =  df["Temperature"]
        height          =  df["Height"]
        nlev            =  df.shape[0]
        missing         =  -99999.
        theta           =  potentialtemperature(temperature, pressure, missing)
        pblp,pblt,pblh  =  pbl_gradient        (missing, theta, pressure, height)
        assert pblp >= bounds[t][0] and pblp <= bounds[t][1]





'''
Test Stable boundary layer
'''
def test_stable_PBL():
    # common variables across tests
    nlev            =  27
    missing         =  -99999.
    theta           =  np.linspace(300,430,27)
    pressure        =  np.linspace(1000,350,27)
    pressure        =  pressure * 1e2
    theta[0]        =  285.0
    temperature     =  calculate_temperature        (theta      ,pressure,missing)
    height          =  calculate_height_above_ground(temperature,pressure,missing)
    pblp,pblt,pblh  =  pbl_gradient        (missing, theta, pressure, height)
    assert pblt >= 294 and pblt <= 296



'''
Test unstable boundary layer
'''
def test_unstable_PBL():
    # common variables across tests
    nlev            =  27
    missing         =  -99999.
    theta           =  np.linspace(300,430,27)
    pressure        =  np.linspace(1000,350,27)
    pressure        =  pressure * 1e2
    theta[0:5]      =  322.5
   
    temperature     =  calculate_temperature        (theta      ,pressure,missing)
    height          =  calculate_height_above_ground(temperature,pressure,missing)
    pblp,pblt,pblh  =  pbl_gradient        (missing, theta, pressure, height)
    print(pblp)
    print(pblt)
    assert pblt >= 320 and pblt <= 325



'''
Isothermal profile
'''






'''
Test to see if the change in integrated column heat is equal 
to the addition of surface sensible heat
'''
'''
def test_change_in_column_heat():
    # desired heat change
    sensible_heat  =  200.0
    dt             =  60.0 * 10.0
    change_in_heat =  sensible_heat * dt
    
    
    # Assume a linear decrease in temperature with pressure (makes math easier)
    pressure        =  np.linspace(1000,350,nlev)
    pressure        =  pressure * 1e2
    temperature     =  9.8e-4 * pressure  +  200.0
    
    
    assert pblt == 295

'''


'''
Test missing values
def test_PBL_missing_value():
    # common variables across tests
    nlev            =  30
    theta           =  np.zeros( nlev )
    missing         =  -99999.
    theta[0:3]      =  missing
    theta[3:]       =  np.linspace(300,430,27)
    pressure        =  np.linspace(1000,350,nlev)
    pressure        =  pressure * 1e2

    # test with buoyancy
    theta[3]        =  322.5
    theta           =  packit( theta[1:], theta[0], nlev-1, pressure[0], pressure[1:], missing ) 
    temperature     =  calculate_temperature        (theta      ,pressure,missing)
    height          =  calculate_height_above_ground(temperature,pressure,missing)
    pblp,pblt,pblh  =  pblheat                      (missing    ,theta   ,pressure,height)
    assert pblt == 322.5
'''
