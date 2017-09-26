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





####################################################
#####  Prepare profile data
####################################################
filename  =  'LAMONT.JUNE_06_2002_1200Z.csv'
df        =  pd.read_csv(filename)
nlev      =  df.shape[0]
nlat      =  1 
nlon      =  1
nday      =  1

#dt       =  300.0     ## 5min timestep
#nhr      =  24 * 12   ## int(dt/60)
dt       =  3600.0     ## 1 hour timestep
nhr      =  24         ## int(dt/60)


missing  =  -99999.0

T = np.zeros( (nday,nlev-1,nlat,nlon) )
Q = np.zeros( (nday,nlev-1,nlat,nlon) )
Z = np.zeros( (nday,nlev-1,nlat,nlon) )
P = np.zeros( (nday,nlev-1,nlat,nlon) )

T2M  = np.zeros( (nday,nlat,nlon) )
Q2M  = np.zeros( (nday,nlat,nlon) )
H2M  = np.zeros( (nday,nlat,nlon) )
Psfc = np.zeros( (nday,nlat,nlon) )

T[0,:,0,0]  =  df[  'Temp'  ][1:] + 273.16
P[0,:,0,0]  =  df['Pressure'][1:]  * 1e2
Z[0,:,0,0]  =  df[ 'Height' ][1:]  
Q[0,:,0,0]  =  df[  'Shum'  ][1:]  * 1e-3

T2M [0,0,0]  =  df[  'Temp'  ][0] + 273.16
Psfc[0,0,0]  =  df['Pressure'][0]  * 1e2
H2M [0,0,0]  =  df[ 'Height' ][0]  
Q2M [0,0,0]  =  df[  'Shum'  ][0]  * 1e-3


####################################################
#####  Prepare surface flux data
####################################################
df_flux =  pd.read_csv('LH_and_SH_2002_June.csv')
df_flux[df_flux <= -999] = np.nan
df_flux =  df_flux.shift(-6)
ef      =  np.linspace(0.05,0.95,10)

lh      =  ma.array(df_flux['LH'], fill_value=missing)
sh      =  ma.array(df_flux['SH'], fill_value=missing)
nrows   =  lh.shape[0]

rnet           =  np.zeros( (nday,nhr,nlat,nlon) )
rnet[0,:,0,0]  =  lh[:nhr] + sh[:nhr]

exit()

#line = plt.plot(rnet[0,:,0,0])
#plt.show()
#df_flux.plot.line()


### Complete run
tbm,bclp,timeofci = evaluate_ci(T,Q,P,Z,3,T2M,Q2M,Psfc,H2M,rnet,ef,dt,missing)



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
ef      =  np.linspace(0.05,0.95,10)

lh      =  ma.array(df_flux['LH'], fill_value=missing)
sh      =  ma.array(df_flux['SH'], fill_value=missing)
nrows   =  lh.shape[0]
#rnet    =  create_grid(lh[:24] + sh[:24],  12)
rnet    =  lh[:nhr] + sh[:nhr]


nlev1 = nlev + 1
tpack = packit(T,T2M ,nlev,Psfc,P,missing)
hpack = packit(Z,H2M ,nlev,Psfc,P,missing)
qpack = packit(Q,Q2M ,nlev,Psfc,P,missing)
ppack = packit(P,Psfc,nlev,Psfc,P,missing)

Theta = potentialtemperature(tpack, ppack, missing )
itime = 3
grav  = 9.81
cp    = 1005.7
cp_g  = cp/grav

omega = 0.1
ilev  = 35


#for zz,th in enumerate(Theta[:ilev]):
#    print('    {}      {}       {}       {}  '.format(zz,ppack[zz]/1e2,th,qpack[zz]*1e3))


plt.subplot(1,2,1)
line1 = plt.plot(Theta [:ilev], Z[:ilev], 'k'  , lw=2.5)
line2 = plt.plot(Theta [:ilev], Z[:ilev], '--r', lw=2.5)

plt.subplot(1,2,2)
line3 = plt.plot(qpack [:ilev]*1e3, Z[:ilev], 'k'  , lw=2.5)
line4 = plt.plot(qpack [:ilev]*1e3, Z[:ilev], '--r', lw=2.5)
plt.show()



print itime, nhr
for tt in np.arange(itime,nhr):

    print('   RNET:  {}   {} '.format(tt,rnet[tt]))
    if rnet[tt] >= 0:
        latent_heat    =  ef[0] * rnet[tt]
        sensible_heat  =  (1.0 - ef[0]) * rnet[tt]

        ##
        pblp,pbl_theta,pblh  =  pblheat     (missing, tpack    , ppack, hpack         )
        newTheta             =  assign_layer(Theta  , pbl_theta, pblh , hpack, missing)
        pbl_depth            =  ppack[0] - pblp
        avgRHO               =  avg_over_layer   ( total_density(ppack,tpack,qpack,missing), pblh, hpack, missing)
        Theta                =  add_sensible_heat( newTheta, sensible_heat, avgRHO, hpack, pblh, dt, missing)
        tpack                =  calculate_temperature(Theta, ppack, missing)

#        newTheta             =  assign_layer(Theta  , pbl_theta, pblh , hpack, missing)
#        Theta                =  newTheta
#        pbl_depth            =  ppack[0] - pblp
#        tpack                =  calculate_temperature(Theta, ppack, missing)

        dpress               =  layerdepth(ppack, tpack[0], qpack[0], ppack[0], hpack[0], missing)
        evap                 =  latent_heat * dt / latent_heat_of_condensation(tpack[0], missing)
        newQhum              =  inject_and_mix( pblp, evap, ppack, qpack, dpress, missing)
        qpack                =  newQhum
        tbm,tdef,bclp        =  hcfcalc(missing, tpack, ppack, qpack, hpack, pblp)
        pressure_deficit     =  (pblp + omega*dt) - bclp

        print tt, latent_heat, sensible_heat, pblp/1e2, pblh/1e3, evap, tdef, pressure_deficit/1e2
        print("  ")
        print("  ")
        print("  ")
        print("  ")


        if np.mod(tt,1) == 0:
            plt.subplot(1,2,1)
            line1 = plt.plot(Theta [:ilev], Z[:ilev], 'k'  , lw=2.5)
            line2 = plt.plot(Theta [:ilev], Z[:ilev], '--r', lw=2.5)

            plt.subplot(1,2,2)
            line3 = plt.plot(qpack [:ilev]*1e3, Z[:ilev], 'k'  , lw=2.5)
            line4 = plt.plot(qpack [:ilev]*1e3, Z[:ilev], '--r', lw=2.5)
            plt.show()







print("   ALL DONE  ")


