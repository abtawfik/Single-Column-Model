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
### This is a collection of micro and macro tests that attempt to test each
### each section of the code.  The tests are sometimes arithmetic tests that use
### very simple examples to make sure the function is behaving properly and physics
### checks in other instants that verify that we are within bounds of reasonable 
### values
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
        assert rho == pytest.approx(zz+1, 0.0001)
    for zz,rho in enumerate(sumdensity_mixing_ratio):
        assert rho == pytest.approx(zz+1, 0.0001)


'''
Test use some dummy data to perform this integration
Solves this equation --> Energy = cp/g integrate T dp from top to bottom
---column_energy---
'''
def test_column_energy():
    missing        =  -99999.0
    temperature    =  [ 7  ,  5  ,   3  ,   1  ]
    pressure       =  [3/12, 1/12, -2/12, -8/12]
    total_energy   =  column_energy( temperature, pressure, missing )
    assert total_energy == pytest.approx(307.5535168, 0.001)



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
                   './test_data/Hi_res_test_profile.csv'                ,
                   './test_data/Hi_res_test_profile.20150815.113000.csv',
                   './test_data/Hi_res_test_profile.20150815.133000.csv',
                   './test_data/Hi_res_test_profile.20150815.153000.csv',
                   './test_data/Hi_res_test_profile.20150815.193000.csv',
                   './test_data/Hi_res_test_profile.20150815.203000.csv',
                   './test_data/Hi_res_test_profile.20150816.142900.csv',
                   './test_data/Hi_res_test_profile.20150816.152900.csv',
                   './test_data/Hi_res_test_profile.20150816.203800.csv',
                   './test_data/Hi_res_test_profile.20150827.113300.csv',
                   './test_data/Hi_res_test_profile.20150827.132800.csv']
    nprofiles  = len(test_files)
    #########----   7/11 173000    7/11 213000    8/26 223000   No time stamp   8/15 113000    8/15 133000    8/15 153000    8/15 193000    8/15 203000
    bounds     = [ (91000,94000), (84000,87000), (83000,86000), (77000,80500), (96500,99000), (97000,99000), (90500,93000), (81000,83000), (82500,85000),
    ########        8/16 142900    8/16 152900    8/16 203800    8/27 113300    8/27 132800
                   (95000,97000), (95000,97000), (77500,81000), (96000,98000), (96000,98000) ]
    
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
Test whether the correct boundary layer integrated energy is returned
This takes dummy data and prescribes a pressure level of the PBL
---return_pbl_energy
'''
def test_return_pbl_energy():
    missing          =  -99999.0
    temperature      =  [ 7  ,  5  ,   3  ,   1  ]
    pressure         =  [3/12, 1/12, -2/12, -8/12]
    pblp             =  2/12
    pbl_energy,ipbl  =  return_pbl_energy(temperature,pressure,pblp,missing)
    assert pbl_energy == pytest.approx(102.5178389, 0.001)

    pblp             =  1/12
    pbl_energy,ipbl  =  return_pbl_energy(temperature,pressure,pblp,missing)
    assert pbl_energy == pytest.approx(102.5178389, 0.001)

    pblp             =  0
    pbl_energy,ipbl  =  return_pbl_energy(temperature,pressure,pblp,missing)
    assert pbl_energy == pytest.approx(102.5178389, 0.001)

    pblp             =  -3/12
    pbl_energy,ipbl  =  return_pbl_energy(temperature,pressure,pblp,missing)
    assert pbl_energy == pytest.approx(205.0356779, 0.001)

    pblp             =  -2/12
    pbl_energy,ipbl  =  return_pbl_energy(temperature,pressure,pblp,missing)
    assert pbl_energy == pytest.approx(205.0356779, 0.001)

    pblp             =  -11/12
    pbl_energy,ipbl  =  return_pbl_energy(temperature,pressure,pblp,missing)
    assert pbl_energy == pytest.approx(307.5535168, 0.001)

    ## edge case where pblp is greater than all pressure levels
    pblp             =  6/12
    pbl_energy,ipbl  =  return_pbl_energy(temperature,pressure,pblp,missing)
    assert pbl_energy == pytest.approx(102.5178389, 0.001)









'''
Test the lowest level inject moisture in PBL subroutine
---inject_moisture---
and the higher level wrapper that returns mixing ratio
---inject_and_mix
'''
def test_inject_moisture():
    ############################################################################################
    #### Test the low-level subroutine to ensure that moisture increases by exactly the flux 
    #### amount in the first layer
    ############################################################################################
    nlev         =  8
    missing      =  -99999.
    mixLevel     =  900.
    tracer       =  np.ones( nlev )
    pressure     =  np.linspace(1000,300,nlev)
    layer_depth  =  max(pressure) - mixLevel 
    flux_in      =  1.0
    dp           =  depthpressure(pressure, missing)
    new_density  =  inject_moisture(tracer,flux_in,mixLevel,pressure,dp,layer_depth,missing)
    assert     new_density[0 ] == 2.0
    assert all(new_density[1:] == 1.0)

    ############################################################################################
    #### Test the full moisture flux injection call
    ############################################################################################
    nlev         =  4
    missing      =  -99999.
    mixLevel     =  99901.9
    mixing_ratio =  np.ones( nlev )/10.
    pressure     =  np.array( [100000, 99901.9, 95000, 80000] )
    flux_in      =  1.0
    dp           =  depthpressure(pressure, missing)
    mixedprofile =  inject_lh_into_pbl(mixLevel,flux_in,pressure,mixing_ratio,missing)
    assert mixedprofile[0] == pytest.approx(0.2, 0.0001)
    assert all(mixedprofile[1:] == 0.1)


    ############################################################################################
    #### Test the individual components that build up to create the moisture injection
    #### it is tested by comparing the sums of the moisture profile before and after injection
    #### the column should also increase exactly the amount of the flux injection
    ############################################################################################
    nlev         =  4
    missing      =  -99999.
    mixLevel     =  98000.
    mixing_ratio =  np.ones( nlev )/10.
    pressure     =  np.array( [100000, 99901.9, 95000, 80000] )
    layer_depth  =  max(pressure) - mixLevel 
    flux_in      =  1.0
    dp           =  depthpressure(pressure, missing)
    rho          =  columndensity     ( mixing_ratio, dp, missing )
    density      =  columndensitynomix(               dp, missing )
    trues        =  np.where( np.ones(nlev) == 1, True, False )
    num_in_layer =  countplus( np.where(pressure>mixLevel,True,False), 1)
    ilayer_top   =  maxindex(1,np.where(pressure>mixLevel,True,False))
    layer_depth  =  sumit   (1,trues[0:ilayer_top],dp[0:ilayer_top])
    new_tracer   =  inject_moisture(density, flux_in, mixLevel, pressure, dp, layer_depth, missing)
    sum_old      =  np.sum( density[:nlev-1] )
    sum_new      =  np.sum( new_tracer[:nlev-1] )
    difference   =  sum_new - sum_old
    assert difference == flux_in





'''
Test whether the pbl theta calculation returns the known
prescribed boundary layer theta -- it should be 300K in this case
Will likely be some numerical errors because I'm doing several transformations
from theta -> temperature -> pbl_energy -> back to boundary layer theta
---get_pbl_theta_using_energy
'''
def test_get_pbl_theta():
    nlev            =  300
    R               =  287.
    grav            =  9.81
    cp              =  1005.7
    dry_lapse       =  -grav/cp
    alpha           =  R/cp
    a_plus_one      =  alpha + 1
    missing         =  -99999.
    sensible        =  0.0
    dt              =  60. * 10.
    
    theta           =  np.zeros(nlev)
    theta           =  300. + theta
    height          =  np.linspace(0,5,nlev) * 1e3
    temperature     =  np.zeros(nlev)
    pressure        =  np.zeros(nlev)
    temperature[0]  =  300.
    pressure[0]     =  100000.0
    for zz in range(1,nlev):
        temperature[zz] = temperature[zz-1]  +  dry_lapse*(height[zz]-height[zz-1])
        dz              = height[zz]-height[zz-1]
        Tavg            = 0.5*(temperature[zz]+temperature[zz-1])
        pressure   [zz] = pressure[zz-1] * np.exp( -(grav*dz)/(R*Tavg) )

    # Try a medium height boundary layer
    pblp            =  93000.
    theta           =  np.where( pressure >= pblp, theta, theta + 1.0 )
    temp            =  calculate_temperature(theta,pressure,missing)
    ipbl            =  maxindex     ( 1, np.where(pressure>=pblp, True,False) )
    pbl_energy      =  column_energy(temp[:ipbl], pressure[:ipbl], missing)
    pbl_theta       =  get_pbl_theta_using_energy(pressure, pbl_energy, ipbl, missing)
    absolute_error  =  np.abs( pbl_theta - theta[0])
    percent_error   =  np.abs( absolute_error / theta[0] ) * 1e2
    assert absolute_error <= 0.005
    assert percent_error  <= 0.005






'''
Test whether this method of finding the theta intersection to
quantify the boundary layer top is accurate.  It only works
on resolved levels
NOTE: This is only used internally in the get_new_pbl_by_adding_sh function
      the more accurate and between model boundary layer calculation is
      handled by the pbl_gradient function. This is a helper function 
---get_new_pbl_using_theta
'''
def test_get_new_pbl_using_theta():
    pass


'''
Test the updating of the theta profile using theta_pbl
Essentially the function fills in theta_pbl for layers below the mixed layer
and then converts from theta to temperature
---reconstruct_temperature
'''
def test_reconstruct_temperature():
    pass


'''
Test whether this method of finding the theta intersection to
quantify the boundary layer top is accurate.  It only works
on resolved levels
NOTE: This is only used internally in the get_new_pbl_by_adding_sh function
      the more accurate and between model boundary layer calculation is
      handled by the pbl_gradient function. This is a helper function 
---get_new_pbl_by_adding_sh
'''
def test_get_new_pbl_by_adding_sh():
    nlev            =  300
    R               =  287.
    grav            =  9.81
    cp              =  1005.7
    dry_lapse       =  -grav/cp
    missing         =  -99999.
    sensible        =  200.0
    dt              =  60. * 10.

    height          =  np.linspace(0,5,nlev) * 1e3
    temperature     =  np.zeros(nlev)
    pressure        =  np.zeros(nlev)
    qhum            =  np.zeros(nlev)
    temperature[0]  =  300.0
    pressure[0]     =  100000.0
    for zz in range(1,nlev):
        temperature[zz] = temperature[zz-1]  +  dry_lapse*(height[zz]-height[zz-1])
        dz              = height[zz]-height[zz-1]
        Tavg            = 0.5*(temperature[zz]+temperature[zz-1])
        pressure   [zz] = pressure[zz-1] * np.exp( -(grav*dz)/(R*Tavg) )
 
    # Try a few initial boundary layer depths Medium, High, and two shallow ones
    #    pblps  =  np.array( [93000, 80000, 99808, 99810] )
    pblps  =  np.array( [ pp+10 for pp in pressure[1:]] )
    for pblp in pblps:
        theta                  =  np.zeros(nlev)
        theta                  =  300. + theta
        theta                  =  np.where( pressure >= pblp, theta, theta + 1.0 )
        temp                   =  calculate_temperature(theta,pressure,missing)
        constructed_temp       =  get_new_pbl_by_adding_sh(temp,pressure,qhum,pblp,sensible,dt,missing)
        initial_column_energy  =  column_energy( temp            , pressure, missing )
        new_column_energy      =  column_energy( constructed_temp, pressure, missing )
        absolute_error         =  (new_column_energy - initial_column_energy) - (sensible * dt)
        percent_error          =  np.abs( absolute_error / (sensible * dt) ) * 1e2
        print('  pblp={}      abs =  {}            percent = {}        '.format(pblp,absolute_error,percent_error))
        print('  Old Energy={}      New Energy={}   Difference should be close to {} J/m2'.format(initial_column_energy,new_column_energy,sensible*dt))
        assert percent_error  <= 2.0






'''
Test the totality of adding surface fluxes to the column.
This test focuses on conservation. Specifically whether the total column
thermal energy and moisture content change exactly the amount delivered
through the sensible and latent heat fluxes respectively.
---test_add_surface_fluxes
'''
def test_add_surface_fluxes():
    nlev            =  300
    R               =  287.
    grav            =  9.81
    cp              =  1005.7
    dry_lapse       =  -grav/cp
    missing         =  -99999.
    sensible        =  200.0
    dt              =  60. * 10.
    latent          =  0.5

    height          =  np.linspace(0,5,nlev) * 1e3
    qhum            =  np.ones (nlev)
    qhum            =  qhum * 2e-3
    temperature     =  np.zeros(nlev)
    pressure        =  np.zeros(nlev)
    temperature[0]  =  300.0
    pressure[0]     =  100000.0
    for zz in range(1,nlev):
        temperature[zz] = temperature[zz-1]  +  dry_lapse*(height[zz]-height[zz-1])
        dz              = height[zz]-height[zz-1]
        Tavg            = 0.5*(temperature[zz]+temperature[zz-1])
        pressure   [zz] = pressure[zz-1] * np.exp( -(grav*dz)/(R*Tavg) )
 
    # Try a few initial boundary layer depths Medium, High, and two shallow ones
#    pblps  =  np.array( [93000, 80000, 99808, 99810] )
    pblps  =  np.array( [ pp+10 for pp in pressure[1:]] )
    for pblp in pblps:
        theta                  =  np.zeros(nlev)
        theta                  =  300. + theta
        theta                  =  np.where( pressure >= pblp, theta, theta + 1.0 )
        temp                   =  calculate_temperature(theta,pressure,missing)
        newT,newQ,newPBLP      =  add_surface_fluxes(temp,height,pressure,qhum,pblp,sensible,latent,dt,missing)
        initial_column_energy  =  column_energy( temp , pressure, missing )
        new_column_energy      =  column_energy( newT , pressure, missing )
        absolute_error         =  (new_column_energy - initial_column_energy) - (sensible * dt)
        percent_error          =  np.abs( absolute_error / (sensible * dt) ) * 1e2
        print('  pblp={}      abs =  {}            percent = {}        '.format(pblp,absolute_error,percent_error))
        print('  Old Energy={}      New Energy={}   Difference should be close to {} J/m2'.format(initial_column_energy,new_column_energy,sensible*dt))
        assert percent_error  <= 2.0












'''
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import pylab
from matplotlib import style
from matplotlib import animation
import matplotlib.cm as cm
from pylab import savefig, figure
style.use('fivethirtyeight')
'''

