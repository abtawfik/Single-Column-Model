#!/bin/python

### Handling Arrays and Reading Grib package
import numpy as np
import numpy.ma as ma
import pandas as pd

### These are helper Fortran functions to speed things up
from HandyMet import *


### Plotting packages
from mpl_toolkits.basemap import Basemap , addcyclic
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import pylab
from matplotlib import style
from matplotlib import animation
import matplotlib.cm as cm
from pylab import savefig, figure
style.use('fivethirtyeight')




def create_grid(inputV, gridSize):
    return  np.repeat(inputV[:,np.newaxis], gridSize, 1).ravel()


def intialize( temp, pres, qhum, missing ):
    height               =  calculate_height_above_ground(temp, pres, missing)
    Theta                =  potentialtemperature         (temp, pres, missing)
    pblp,pbl_theta,pblh  =  pblheat                      (missing, Theta, pres, height)
    newTheta             =  assign_layer                 (Theta, pbl_theta, pblh, height, missing)
    Theta                =  newTheta
    avgRHO               =  avg_over_layer( total_density(pres,temp,qhum,missing), pblh, height, missing )
    return pblp,pbl_theta,pblh,avgRHO,Theta,height

def relative_humidity(t, p, q, missing):
    qsat  =  saturationhumidity( t, p, missing )
    relh  =  np.where( q == 0, 0.001, q/qsat * 1e2 )
    return relh


def advance_in_time( temp0, pres0, qhum0, pblh0, pblp0, latent, sensible, dt, missing ):
    pblh                 =  pblh0
    pblp                 =  pblp0
    temp                 =  temp0
    qhum                 =  qhum0
    pres                 =  pres0
#    print('           1')
    Theta                =  potentialtemperature         (temp, pres, missing)
#    print('           2')
    height               =  calculate_height_above_ground(temp, pres, missing)
#    print('           3')
    avgRHO               =  avg_over_layer( total_density(pres,temp,qhum,missing), pblh, height, missing )
#    print('           4')
#    print(pblp.shape)
#    print(sensible.shape)
#    print(pres[0].shape)
#    print('')
#    print('')
#    print('')
    newTheta             =  Theta
    Theta                =  add_sensible_heat(newTheta, sensible, pres[0]-pblp, height, pblh, dt, missing)
#    Theta                =  add_sensible_heat(newTheta, sensible, avgRHO, height, pblh, dt, missing)
    print(sensible,pres[0]-pblp,latent/(sensible+latent),dt * (sensible / ((pres[0]-pblp)*(1005.7/9.81))))
    temp                 =  calculate_temperature(Theta, pres, missing)
#    print('           6    *****')
    pblp,pbl_theta,pblh  =  pblheat                      (missing, Theta, pres, height)
#    print('           6    =====')
    newTheta             =  assign_layer (Theta, pbl_theta, pblh, height, missing)
    Theta                =  newTheta
    pbl_depth            =  pres[0] - pblp
#    print('           7')
    temp                 =  calculate_temperature(Theta, pres, missing)
    dpress               =  layerdepth( pres, temp[0], qhum[0], pres[0], height[0], missing )
    evap                 =  latent * dt / latent_heat_of_condensation(temp[0], missing)
#    print('           8')
    newQ                 =  inject_and_mix(pblp, evap, pres, qhum, dpress, missing)
    qhum                 =  newQ
#    print('           9')
    tbm,tdef,bclp,qbcl   =  hcfcalc(missing, temp, pres, qhum, height, pblp)
    pressure_deficit     =  pblp - bclp
#    print('           10')
    return pressure_deficit,temp,qhum,Theta,height,pblp,pblh







####################################################
#####  Prepare profile data
####################################################
#filename   =  'LAMONT.JUNE_06_2002_1200Z.csv'
filename   =  'LAMONT.June_03_2006_1200Z.csv'
df         =  pd.read_csv(filename)
nlev       =  df.shape[0]
nlat       =  1 
nlon       =  1
nday       =  1
missing    =  -99999.0


### Break down of components
T = np.zeros( (nlev-1) )
Q = np.zeros( (nlev-1) )
Z = np.zeros( (nlev-1) )
P = np.zeros( (nlev-1) )

T2M  = np.zeros( 1 )
Q2M  = np.zeros( 1 )
H2M  = np.zeros( 1 )
Psfc = np.zeros( 1 )

T[:]     =  df[  'Temp'  ][1:]  + 273.16
P[:]     =  df['Pressure'][1:]  * 1e2
Z[:]     =  df[ 'Height' ][1:]  
Q[:]     =  df[  'Shum'  ][1:]  * 1e-3

T2M [0]  =  df[  'Temp'  ][0 ]  + 273.16
Psfc[0]  =  df['Pressure'][0 ]  * 1e2
H2M [0]  =  df[ 'Height' ][0 ]  
Q2M [0]  =  df[  'Shum'  ][0 ]  * 1e-3

df_flux =  pd.read_csv('LH_and_SH_2002_June.csv')
df_flux[df_flux <= -999] = np.nan
df_flux =  df_flux.shift(-6)

lh      =  ma.array(df_flux['LH'], fill_value=missing)
sh      =  ma.array(df_flux['SH'], fill_value=missing)
nhr     =  lh[:24].shape[0]
rnet    =  lh[:nhr] + sh[:nhr]


#exit()


ef         =  np.linspace(0.05,0.95,10)
numEF      =  ef.shape[0]
#dt         =  600.   
dt         =  300.   

min5       =  np.arange(0,nhr-1,5.0/60.0)
#min5       =  np.arange(0,nhr-1,10.0/60.0)
timeOfday  =  np.arange(0,nhr  ,1) * 1.0
polyno     =  np.poly1d( np.polyfit(timeOfday,rnet,10) )
NETRAD     =  polyno(min5)



nlev1  =  nlev + 1
tpack  =  create_grid(packit(T,T2M ,nlev,Psfc,P,missing), numEF).reshape(nlev,numEF)
qpack  =  create_grid(packit(Q,Q2M ,nlev,Psfc,P,missing), numEF).reshape(nlev,numEF)
#hpack  =  create_grid(packit(Z,H2M ,nlev,Psfc,P,missing), numEF).reshape(nlev,numEF)
ppack  =  create_grid(packit(P,Psfc,nlev,Psfc,P,missing), numEF).reshape(nlev,numEF)
#hpack  =  np.zeros( (nlev,numEF) )
#hpack  =  np.array([ calculate_height_above_ground(tpack[:,ee],ppack[:,ee],missing) for ee,evapf in enumerate(ef) ])
#hpack  =  np.transpose(np.array([ calculate_height_above_ground(tpack[:,ee],ppack[:,ee],missing) for ee,evapf in enumerate(ef) ]))

pblp       =  np.zeros(       numEF  )
pblh       =  np.zeros(       numEF  )
pbl_theta  =  np.zeros(       numEF  )
avgRHO     =  np.zeros(       numEF  )
Theta      =  np.zeros( (nlev,numEF) )
hpack      =  np.zeros( (nlev,numEF) )
relh       =  np.zeros( (nlev,numEF) )
for ee,evapf in enumerate(ef):
    pblp[ee],pbl_theta[ee],pblh[ee],avgRHO[ee],Theta[:,ee],hpack[:,ee]  =  intialize( tpack[:,ee], ppack[:,ee], qpack[:,ee], missing )
    relh[:,ee]   =  relative_humidity(tpack[:,ee], ppack[:,ee], qpack[:,ee], missing)

#hpack  =  np.transpose([ calculate_height_above_ground(tpack[:,ee],ppack[:,ee],missing) for ee,evapf in enumerate(ef) ])
#Theta  =  np.transpose([ potentialtemperature         (tpack[:,ee],ppack[:,ee],missing) for ee,evapf in enumerate(ef) ])
#pblp,pbl_theta,pblh  =  pblheat(missing,tpack[:,0],ppack[:,0],hpack[:,0])
#pblp       =  create_grid( [pblp    ] , numEF )
#pbl_theta  =  create_grid( pbl_theta, numEF )
#pblh       =  create_grid( pblh     , numEF )


grav   =  9.81
cp     =  1005.7
cp_g   =  cp/grav

omega  =  0.0
ilev   =  35



firstP = plt.subplot(1,2,1)
line1  = plt.plot(Theta [:ilev,0], hpack[:ilev,0]*1e-3, 'k'  , lw=2.5)
line2  = plt.plot(Theta [:ilev,0], hpack[:ilev,0]*1e-3, '--r', lw=2.5)
firstP.set_ylim([-0.3,9])
firstP.set_xlim([290,340])
secondP = plt.subplot(1,2,2)
#line3 = plt.plot(qpack [:ilev,0]*1e3, hpack[:ilev,0], 'k'  , lw=2.5)
#line4 = plt.plot(qpack [:ilev,0]*1e3, hpack[:ilev,0], '--r', lw=2.5)
line3 = plt.plot(relh [:ilev,0], hpack[:ilev,0]*1e-3, 'k'  , lw=2.5)
line4 = plt.plot(relh [:ilev,0], hpack[:ilev,0]*1e-3, '--r', lw=2.5)
secondP.set_ylim([-0.3,9])
secondP.set_xlim([0,100])
plt.show()


line_colors  =  [ '#F08080', '#7a4d1f', '#a2662a', '#d08c49', '#e5bf9a', '#99ccff', '#1e90ff', '#0073e6', '#004d99', 'k'  ]

triggered  =  np.zeros( numEF )
for tt,rnet in enumerate(NETRAD):
    print('   RNET:  {}   {} '.format(tt,rnet))
    if rnet >= 0:
        for ee,evapf in enumerate(ef):
            latent_heat    =  evapf * rnet
            sensible_heat  =  (1.0 - evapf) * rnet

            if triggered[ee] == 0:
                pressure_deficit, tpack[:,ee], qpack[:,ee], Theta[:,ee], hpack[:,ee], pblp[ee], pblh[ee]  =  \
                advance_in_time( tpack[:,ee], ppack[:,ee], qpack[:,ee], pblh[ee], pblp[ee], latent_heat, sensible_heat, dt, missing )
                relh[:,ee]   =  relative_humidity(tpack[:,ee], ppack[:,ee], qpack[:,ee], missing)            

            if pressure_deficit <= 0:
                triggered[ee]  =  1

            if np.mod(tt,1) == 0:
                first = plt.subplot(2,2,1)
                line1 = plt.plot(Theta [:ilev,ee], hpack[:ilev,ee]*1e-3, line_colors[ee]  , lw=2.5)
                first.set_ylim([-0.3,9])
                first.set_xlim([290,340])
                second = plt.subplot(2,2,2)
#                line3 = plt.plot(qpack [:ilev,ee]*1e3, hpack[:ilev,ee], line_colors[ee]  , lw=2.5)
                line3 = plt.plot(relh [:ilev,ee], hpack[:ilev,ee]*1e-3, line_colors[ee]  , lw=2.5)
                second.set_ylim([-0.3,9])
                second.set_xlim([0,100])

                third = plt.subplot(2,2,(3,4))
                line4 = plt.plot(min5[:tt+1], (1.0 - evapf) * NETRAD[:tt+1], line_colors[ee]  , lw=2.5)
                third.set_xlim([0,24])
                third.set_ylim([-10,620])

        if np.mod(tt,1) == 0:
            outfile  =  './animation/profile_with_rad{}.png'.format(tt)
            savefig( outfile, bbox_inches='tight')
#            plt.show()






print("   ALL DONE  ")


