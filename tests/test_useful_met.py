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
    mixedprofile =  inject_and_mix(mixLevel,flux_in,pressure,mixing_ratio,dp,missing)
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
Test to see if when given an atmospheric profile of temperature,
the code can return a PBL energy [J/m2] and then reconstruct the original
temperature profile.  This shows the integration works and no info is lost
during the transform.
'''
def test_pbl_energy_temp_reconstruction():
    nlev            =  300
    R               =  287.
    grav            =  9.81
    cp              =  1005.7
    dry_lapse       =  -grav/cp
    missing         =  -99999.
    height          =  np.linspace(0,5,nlev) * 1e3
    temperature     =  np.zeros(nlev)
    pressure        =  np.zeros(nlev)
    temperature[0]  =  300.0
    pressure[0]     =  100000.0
    for zz in range(1,nlev):
        temperature[zz] = temperature[zz-1]  +  dry_lapse*(height[zz]-height[zz-1])
        dz              = height[zz]-height[zz-1]
        Tavg            = 0.5*(temperature[zz]+temperature[zz-1])
        pressure   [zz] = pressure[zz-1] * np.exp( -(grav*dz)/(R*Tavg) )

    pblp              =  93000.
    ipbl              =  maxindex(1,np.where(pressure >= pblp,True,False))
    pbl_energy        =  column_energy(temperature[0:ipbl],pressure[0:ipbl],missing)
    constructed_temp  =  reconstruct_pbl_temperature(temperature,height,pressure,pbl_energy,ipbl,missing)
    absolute_error    =  constructed_temp - temperature
    percent_error     =  np.abs( absolute_error / temperature ) * 1e2
    assert all(absolute_error <= 0.1 )
    assert all(percent_error  <= 0.05)


    pblp              =  99904.7995543
    ipbl              =  maxindex(1,np.where(pressure >= pblp,True,False))
    pbl_energy        =  column_energy(temperature[0:ipbl],pressure[0:ipbl],missing)
    constructed_temp  =  reconstruct_pbl_temperature(temperature,height,pressure,pbl_energy,ipbl,missing)
    absolute_error    =  constructed_temp - temperature
    percent_error     =  np.abs( absolute_error / temperature ) * 1e2
    assert all(absolute_error <= 0.1 )
    assert all(percent_error  <= 0.05)





'''
Test Hi-resoluation real-world PBL from
Craig Enhanced Sounding Data
'''
def test_hi_res_add_sensible_heat():
    missing         =  -99999.
    theta           =  np.linspace(0,10,11) + 1
    dt              =  1.0
    sensible_heat   =  1.0
    pbl_depth       =  9.81/1005.7
    newT            =  add_sensible_heat (theta[0], sensible_heat, pbl_depth, dt, missing )
    assert newT == 2 





'''
nlev         =  4
missing      =  -99999.
mixLevel     =  99901.9
mixing_ratio =  np.ones( nlev )/10.
pressure     =  np.array( [100000, 99901.9, 95000, 80000] )
flux_in      =  1.0
dp           =  depthpressure(pressure, missing)
mixedprofile =  inject_and_mix(mixLevel,flux_in,pressure,mixing_ratio,dp,missing)


nlev         =  8
missing      =  -99999.
mixLevel     =  900.
layer_depth  =  100.
tracer       =  np.ones( nlev )
pressure     =  np.linspace(1000,300,nlev)
flux_in      =  1.0
dp           =  depthpressure(pressure, missing)
new_density  =  inject_moisture(tracer,flux_in,mixLevel,pressure,dp,layer_depth,missing)
'''


'''
Code not test yet:
---- total_density
---- packIt
---- assign_layer
---- avg_over_layer
---- add_sensible_heat
---- idealGas_P
---- virtualTemp   
---- Latent_heat_of_condensation
---- midLevel
---- trueMidLevel
---- surfacePressure
---- virtualPotentialTemp
---- virtualTemperature
---- layerDepth
---- logMiss
---- inject_and_mix
----   

'''





'''
missing         =  -99999.
theta           =  np.linspace(0,10,11) + 1
dt              =  1.0
sensible_heat   =  1.0
pbl_depth       =  9.81/1005.7
newT            =  add_sensible_heat (theta[0], sensible_heat, pbl_depth, dt, missing )


missing         =  -99999.
temperature     =  np.linspace(0,10,11)
height          =  np.linspace(0,10,11)
pressure        =  temperature*1000 + (90000.)
pressure        =  pressure[::-1]
pblh            =  3.0
dt              =  1.0
sensible_heat   =  20.0
pbl_depth       =  0.1950879984
newT            =  add_sensible_heat (temperature, sensible_heat, pbl_depth, height, pblh, dt, missing )
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



'''
import numpy as np
import numpy.ma as ma
import pandas as pd
import sys
from HandyMet import *

from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import pylab
from matplotlib import style
from matplotlib import animation
import matplotlib.cm as cm
from pylab import savefig, figure
style.use('fivethirtyeight')

dt              =  3 * 3600.0
sensible_heat   =  700.0

test_file       =  './test_data/Hi_res_test_profile.20150815.203000.csv'
df              =  pd.read_csv(test_file)
pressure        =  df["Pressure"]
temperature     =  df["Temperature"]
height          =  df["Height"]
nlev            =  df.shape[0]
missing         =  -99999.
theta           =  potentialtemperature(temperature, pressure, missing)
pblp,pblt,pblh  =  pbl_gradient(missing,theta,pressure,height)
ilower          =  maxindex(1, np.where( pressure>=pblp, True, False) )
pbl_energy      =  column_energy(temperature[:ilower],pressure[:ilower],missing)
new_pbl_energy  =  add_sensible_heat_energy(pbl_energy,sensible_heat,dt,missing)
dp              =  depthpressure(pressure, missing)
new_temperature = reconstruct_pbl_temperature(temperature,height,dp,pbl_energy,ilower,missing)
add_temperature = reconstruct_pbl_temperature(temperature,height,dp,new_pbl_energy,ilower,missing)

theta           =  potentialtemperature(    temperature, pressure, missing)
comp_theta      =  potentialtemperature(new_temperature, pressure, missing)
add_theta       =  potentialtemperature(add_temperature, pressure, missing)

avg_T  =  np.avg(temperature[:ilower])
newT   =  np.array( [ avg_T - (g_cp * (height[zz]-(0.5 * (height[ilower] + height[0])))) if zz < ilower else temperature[zz] for zz in range(0,len(height)) ] )

f, ax1 = plt.subplots(1,1)
ax1.plot(comp_theta, height/1e3, '-k')
ax1.plot(     theta, height/1e3, 'ob')
ax1.plot( add_theta, height/1e3, 'or')
ax1.set_ylim([0,3])
ax1.set_xlim([280,320])
plt.show()


f, ax1 = plt.subplots(1,1)
ax1.plot(new_temperature, height/1e3, '-k')
ax1.plot(    temperature, height/1e3, 'ob')
ax1.plot(add_temperature, height/1e3, 'or')
ax1.plot( newT          , height/1e3, 'ow')
ax1.set_ylim([0,3])
ax1.set_xlim([280,310])
plt.show()


dtdp         =  smoothed_lapse_rate(theta,pressure,missing)
f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
ax1.plot(theta         , pressure/1e2         , '-k')
ax1.plot(smoothed_theta, smoothed_pressure/1e2, '--r')
ax1.set_ylim([1000,750])
ax1.set_xlim([300 ,312])

ax2.plot(dtdp, pressure/1e2, 'ob')
ax2.set_ylim([1000   , 750  ])
ax2.set_xlim([-0.004,0.004])
plt.gca().invert_yaxis()

ax3.plot(dtdp, smoothed_pressure/1e2, 'ob')
ax3.set_ylim([1000   , 750  ])
ax3.set_xlim([-0.004,0.004])

plt.show()
del ax1
del ax2
del ax3
'''

