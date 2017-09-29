#!/bin/python
'''

All the tests in this file are for the modules in Useful_Met_Functions.f90

'''


################################################
### Handling Arrays and Reading Grib package
################################################
import numpy as np



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
Test theta first about the pblHeat subroutine
'''
def test_pbl_heat_theta():

    # common variables across tests
    nlev            =  27
    missing         =  -99999.
    theta           =  np.linspace(300,430,27)
    pressure        =  np.linspace(1000,350,27)
    pressure        =  pressure * 1e2

    # test of stable boundary layer
    theta[0]        =  285.0
    temperature     =  calculate_temperature        (theta      ,pressure,missing)
    height          =  calculate_height_above_ground(temperature,pressure,missing)
    pblp,pblt,pblh  =  pblheat                      (missing    ,theta   ,pressure,height)
    assert pblt == 295

    # test with buoyancy
    theta[0]        =  322.5
    temperature     =  calculate_temperature        (theta      ,pressure,missing)
    height          =  calculate_height_above_ground(temperature,pressure,missing)
    pblp,pblt,pblh  =  pblheat                      (missing    ,theta   ,pressure,height)
    assert pblt == 322.5

