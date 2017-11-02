
!-----------------------------------------------------------------------------
!
! Primary Subroutine: At a high-level this subroutine takes real atmospheric
! profiles and sees whether convection initiates under various evaporative 
! fraction regimes.  This subroutine begins with an early morning sounding
! and then evolves forward in time with a given net radiation that is split by 
! the desired evaporative fraction.  Overall the subroutine can explore the entire
! possible space of evaporative fraction for a given day to see whether CI occurred
! - Is convection more likely over dry or 'wet' surface flux properties?
! - Applying it spatially can give an answer to the question of positive or
!   negative flux-CI feedbacks over certain regions and times of year
!
!-----------------------------------------------------------------------------
subroutine Evaluate_CI_plevels( T   , Q   , P        , omega , itime                     , &
                                t2m , q2m , psfc     , rnet  , ef      , dt    , nhr8    , &
                                nlev, nlat, nlon     , nhr   , nday    , num_EF, missing , &
                                TBM , BCLP, TimeOfCI , QBCL   )   !!! , CAPE                      )

   implicit none
!
! Input/Output Variables
!
   integer, intent(in  )                                       ::  nday        ! *** # days
   integer, intent(in  )                                       ::  nhr         ! *** # number of hours per day
   integer, intent(in  )                                       ::  nhr8        ! *** # number of hours per day Original
   integer, intent(in  )                                       ::  nlev        ! *** # of atmospheric levels
   integer, intent(in  )                                       ::  nlat        ! *** # latitude
   integer, intent(in  )                                       ::  nlon        ! *** # longitude
   integer, intent(in  )                                       ::  num_EF      ! *** # number of evaporative fraction breakdowns 
   integer, intent(in  )                                       ::  itime       ! *** Start time of morning sounding

   real(4), intent(in  )                                       ::  missing     ! *** Missing values
   real(4), intent(in  )                                       ::  dt          ! *** timestep in seconds [s/timestep]
   real(4), intent(in  ), dimension(nday,nhr8,nlev,nlat,nlon)  ::  omega       ! *** vertical velocity [Pa/s]
   real(4), intent(in  ), dimension(nday     ,nlev,nlat,nlon)  ::  T  , Q      ! *** Temp and Humidity over level in SI
   real(4), intent(in  ), dimension(nday          ,nlat,nlon)  ::  t2m, q2m    ! *** 2-m temp and humidity, Height in SI

   real(4), intent(in  ), dimension(          nlev          )  ::  P           ! *** Pressure in (level) SI
   real(4), intent(in  ), dimension(nday          ,nlat,nlon)  ::  psfc        ! *** 2-m pressure and height in SI
   real(4), intent(in  ), dimension(nday,nhr      ,nlat,nlon)  ::  rnet        ! *** Net Radiation time series [W/m2]
   real(4), intent(in  ), dimension(     num_EF             )  ::  ef          ! *** Evaporative Fraction levels

   real(4), intent(out ), dimension(nday,num_EF   ,nlat,nlon)  ::  TBM         ! *** Buoyant mixing theta at time of CI [K]
   real(4), intent(out ), dimension(nday,num_EF   ,nlat,nlon)  ::  BCLP        ! *** Cloud base pressure at time of CI [K]
   real(4), intent(out ), dimension(nday,num_EF   ,nlat,nlon)  ::  QBCL        ! *** Humidity at BCL [kg/kg]
!   real(4), intent(out ), dimension(nday,num_EF   ,nlat,nlon)  ::  CAPE        ! *** CAPE from NARR [J/kg]
   real(4), intent(out ), dimension(nday,num_EF   ,nlat,nlon)  ::  TimeOfCI    ! *** Time of day of CI [0-24 hour]


!
! Local variables
!
!   integer, parameter         ::  itime = 3  
!   real(4), parameter         ::  omega = 0.0
   real(4), dimension(nlev+1) ::  ppack, tpack, hpack, qpack, Theta, newTheta, newQhum, dpress, density, OmegaPack
   real(4)                    ::  pbl_depth, pblp, pbl_theta, pblh, latent_heat, sensible_heat
   real(4)                    ::  evap, Lc, tbm_out, tdef_out, bclp_out, pressure_deficit
   real(4)                    ::  avgRHO, qbcl_out, avgOMEGA
   integer                    ::  xx, yy, dd, tt, ee, zz
   integer                    ::  nlev1



!-----------------------------------------------------------------------------


      !-----------------------------------------------------------------------------
      !-- Initialize output variables
      !-----------------------------------------------------------------------------
      TBM       =  missing
      BCLP      =  missing
      QBCL      =  missing
!      CAPE      =  missing
      TimeOfCI  =  missing
      nlev1     =  nlev + 1


      !-----------------------------------------------------------------------------
      !-- Loop over time, lat, and lon
      !-----------------------------------------------------------------------------
      lati_loop: do yy = 1,nlat 
      long_loop: do xx = 1,nlon
      ef_loop:   do ee = 1,num_EF
      day_loop:  do dd = 1,nday


         !--------------------------------------
         !-- Append 2-m quantiy to the profile 
         !--------------------------------------
         call packIt(T    (dd,:,yy,xx), t2m(dd  ,yy,xx) , nlev   , nlev1, &
                     psfc (dd  ,yy,xx), P               , missing, tpack  )

         call packIt(Q    (dd,:,yy,xx), q2m(dd  ,yy,xx) , nlev   , nlev1, &
                     psfc (dd  ,yy,xx), P               , missing, qpack  )

         call packIt(P                , psfc(dd  ,yy,xx), nlev   , nlev1, &
                     psfc (dd  ,yy,xx), P               , missing, ppack  )

         !-----------------------------------------------------------------------------
         !-- Calculate the potential temperature [This is mutable]
         !-----------------------------------------------------------------------------
         call potentialTemperature(tpack, ppack, nlev1, missing, Theta)

         !-----------------------------------------------------------------------------
         !-- Calculate the height above ground [m]
         !-----------------------------------------------------------------------------
         call calculate_height_above_ground (tpack, ppack, nlev1, missing, hpack)

         !-----------------------------------------------------------------------------
         !-- Calcuate the boundary layer top variables but with no extra water vapor
         !-- added to the column
         !-----------------------------------------------------------------------------
         call pblHeat ( nlev1, missing, Theta, ppack, hpack, pblp, pbl_theta, pblh )

         !-----------------------------------------------------------------------------
         !-- Assign PBL quantities across the depth of the PBL
         !-----------------------------------------------------------------------------
         call assign_layer (Theta, pbl_theta, pblh, hpack, nlev1, missing, newTheta )
         Theta      =  newTheta
         pbl_depth  =  ppack(1) - pblp

         !-----------------------------------------------------------------------------------
         !-- Calculate density for column -> used for sensible heat calculation
         !-----------------------------------------------------------------------------------
         call total_density    ( ppack  , tpack, qpack, nlev1, missing, density )
         call avg_over_layer   ( density, pblh , hpack, nlev1, missing, avgRHO  )

         !-----------------------------------------------------------------------------------
         !-- Return corresponding temperature given the previous heat added
         !-----------------------------------------------------------------------------------
         call calculate_temperature(Theta, ppack, nlev1, missing, tpack)




         !-----------------------------------------------------------------------------
         !-- Loop over time of day until Initiation occurs
         !-- Running a little mini-simulation!!!!!  
         !-----------------------------------------------------------------------------
         time_of_day: do tt = itime,nhr

            !-----------------------------------------------------------------------------
            !-- Don't do anything if there is no heat or moisture being added
            !-- Otherwise partition the evaporative fraction between LH and SH
            !-----------------------------------------------------------------------------
            if( rnet(dd,tt,yy,xx).le.0 ) cycle time_of_day
            latent_heat    =  ef(ee)         * rnet(dd,tt,yy,xx)
            sensible_heat  =  (1.0 - ef(ee)) * rnet(dd,tt,yy,xx)


!    write(*,*) '--->>   ',tt, pblh/1e3,pblp/1e2,latent_heat,sensible_heat
            !-----------------------------------------------------------------------------
            !-- Calculate the height above ground [m]
            !-----------------------------------------------------------------------------
            call calculate_height_above_ground (tpack, ppack, nlev1, missing, hpack)


            !-----------------------------------------------------------------------------
            !-- Calcuate the boundary layer top variables but with no extra water vapor
            !-- added to the column
            !-----------------------------------------------------------------------------
            !call pblHeat ( nlev1, missing, tpack, ppack, hpack, pblp, pbl_theta, pblh )
            !write(*,*) '--->>   ',tt, pblh/1e3,pblp/1e2,latent_heat,sensible_heat


            !-----------------------------------------------------------------------------
            !-- Assign PBL quantities across the depth of the PBL
            !-----------------------------------------------------------------------------
            !call assign_layer (Theta, pbl_theta, pblh, hpack, nlev1, missing, newTheta )
            !Theta      =  newTheta
            !pbl_depth  =  ppack(1) - pblp


            !-----------------------------------------------------------------------------------
            !-- Calculate density for column -> used for sensible heat calculation
            !-----------------------------------------------------------------------------------
            call total_density    ( ppack  , tpack, qpack, nlev1, missing, density )
            call avg_over_layer   ( density, pblh , hpack, nlev1, missing, avgRHO  )


            !-----------------------------------------------------------------------------------
            !-- Adds some sensible heat flux to the boundary layer to increase the temperature
            !-----------------------------------------------------------------------------------
            newTheta   =  Theta
            pbl_depth  =  ppack(1) - pblp
!            call add_sensible_heat( newTheta, sensible_heat, avgRHO, hpack, pblh, dt, nlev1, missing, Theta ) 
            call add_sensible_heat( newTheta, sensible_heat, pbl_depth, hpack, pblh, dt, nlev1, missing, Theta ) 
            if( newTheta(1).gt.400  .or. Theta(1).gt.400 ) then
               write(*,*)  "!!!! MAJOR ERROR !!!!"
               write(*,*) '       indices             ',dd,yy,xx,ee,tt
               write(*,*)  pblh, avgRHO, sensible_heat, dt
               write(*,*)  ( (dt * sensible_heat) / (avgRHO*1005.7*pblh) )
               write(*,*)  newTheta(1), Theta(1)
               write(*,*)  pblp/1e2, pbl_theta
               do zz=1,nlev1
                  write(*,*)  hpack(zz), ppack(zz), Theta(zz), newTheta(1) - newTheta(zz), Theta(1) - Theta(zz)!density(zz)
               end do
               return
            end if


            !-----------------------------------------------------------------------------------
            !-- Return corresponding temperature given the previous heat added
            !-----------------------------------------------------------------------------------
            call calculate_temperature(Theta, ppack, nlev1, missing, tpack)





!  do zz=1,nlev1
!      write(*,*)  hpack(zz), ppack(zz), Theta(zz), Theta(1) - Theta(zz)   !density(zz)
!  end do


            !-----------------------------------------------------------------------------
            !-- Calcuate the boundary layer top variables but with no extra water vapor
            !-- added to the column
            !-----------------------------------------------------------------------------
            call pblHeat ( nlev1, missing, Theta, ppack, hpack, pblp, pbl_theta, pblh )


!  write(*,*) '---##   ',tt, pblh/1e3,pblp/1e2
!  write(*,*) "  "
!  write(*,*) "  "
!  write(*,*) "  "
!  write(*,*) "  "
!  write(*,*) "  "
!  write(*,*) "  "
!  write(*,*) "  "

            !-----------------------------------------------------------------------------
            !-- Assign PBL quantities across the depth of the PBL
            !-----------------------------------------------------------------------------
            call assign_layer (Theta, pbl_theta, pblh, hpack, nlev1, missing, newTheta )
            Theta      =  newTheta
            pbl_depth  =  ppack(1) - pblp
            
            !-----------------------------------------------------------------------------------
            !-- Return corresponding temperature given the previous heat added
            !-----------------------------------------------------------------------------------
            call calculate_temperature(Theta, ppack, nlev1, missing, tpack)






            !-----------------------------------------------------------------------------------
            !-- Adds some latent heat flux to the boundary layer to increase the speciific humidity
            !-----------------------------------------------------------------------------------
            !-----------------------------------------------------------------------------
            !----- Get the depth of each layer
            !-----------------------------------------------------------------------------
            call layerDepth( ppack, tpack(1), qpack(1), ppack(1), hpack(1), nlev1, missing, dpress )

            !-----------------------------------------------------------------------------
            !----- Convert from latent heat flux [W/m2] to evapotranspiration [kg/m2/timestep]
            !-----------------------------------------------------------------------------
            call Latent_heat_of_condensation(tpack(1), missing, Lc)
            evap  =  latent_heat * dt / Lc 

            !-----------------------------------------------------------------------------
            !----- inject into the boundary layer -- NOTE returns humidity as kg/m2
            !-----------------------------------------------------------------------------
            call inject_and_mix( pblp, evap, ppack, qpack, dpress, nlev1, missing, newQhum )
            qpack  =  newQhum

            !-----------------------------------------------------------------------------
            !----- Calculate HCF variables in new system and check if CI occurred
            !-----------------------------------------------------------------------------
            !call hcfcalc ( nlev1, missing, tpack, ppack, qpack, hpack, pblp, tbm_out, tdef_out, bclp_out )
            call hcfcalc ( nlev1, missing, tpack, ppack, qpack, hpack, pblp, tbm_out, tdef_out, bclp_out, qbcl_out )


            !-----------------------------------------------------------------------------
            !----- Average vertical velocity at and below the Boundary Layer
            !-----------------------------------------------------------------------------
            call packIt(omega(dd,(tt/18)+1,:,yy,xx) , 0.0 , nlev   , nlev1      , &
                        psfc (dd            ,yy,xx) , P   , missing, OmegaPack    )

            call avg_over_layer( OmegaPack, pblh , hpack, nlev1, missing, avgOMEGA  )


            !-----------------------------------------------------------------------------
            !----- Check for CI
            !-----------------------------------------------------------------------------
            pressure_deficit  =  (pblp + avgOMEGA*dt) - bclp_out 
            if( pressure_deficit.le.0 ) then

               !!!write(*,*) "  We have achieved CI  ", dd,tt, tbm_out, bclp_out/1e2, pressure_deficit/1e2
               TBM     (dd,ee,yy,xx)  =  tbm_out
               BCLP    (dd,ee,yy,xx)  =  bclp_out
               TimeOfCI(dd,ee,yy,xx)  =  tt * 1.0
               QBCL    (dd,ee,yy,xx)  =  qbcl_out
!               CAPE    (dd,ee,yy,xx)  =  cape_in(dd,(tt/18)+1,yy,xx)

               !!!!add CAPE calculation later to figure out some sort of intensity
!               call calculate_temperature(tbm_out, ppack, nlev1, missing, tbcl_out)
!               call calc_cape( nlev1, tpack, qpack, ppack, tbcl_out, qbcl_out, bclp_out, CAPE(dd,ee,yy,xx), missing )

               exit time_of_day
            end if

!!! Remember tpack and qpack need to be re-assigned  qpack = newQhum
!!! Add an omega feature to calculating CI
!!! make sure to define all these variables
!!! define what ef is -- maybe pass it in?
!!! exit this time_of_day loop if CI occurs
!!! clean up code?

         end do time_of_day


      end do day_loop
      end do ef_loop
      end do long_loop
      end do lati_loop



end subroutine Evaluate_CI_plevels














!----------------------------------------------------------------------------------------------
!
! function: Calculate the height above the ground using the hypsometric equation
!
!----------------------------------------------------------------------------------------------
subroutine calculate_height_above_ground (temperature, pressure, nlev, missing, height)

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  temperature(nlev), pressure(nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  height(nlev)
      real(4), parameter    ::  grav = 9.81, Rd = 287.04

      real(4)               ::  Tavg
      integer               ::  zz

      Tavg       =  missing
      height     =  missing
      height(1)  =  0.0
      do zz = 2,nlev

         if( temperature(zz-1).ne.missing  .and.  temperature(zz).ne.missing ) then 
            Tavg  =  0.5 * (temperature(zz-1) + temperature(zz))
         end if

         if( height  (zz-1).ne.missing  .and.  Tavg        .ne.missing  .and.  &
             pressure(zz-1).ne.missing  .and.  pressure(zz).ne.missing         ) then 
            height(zz)  =  height(zz-1)  +  ((Rd * Tavg)/grav  *  log(pressure(zz-1)/pressure(zz)) )
         end if
      end do

end subroutine calculate_height_above_ground






!---------------------------------------------------------------------------------
! function: Center finite difference
!---------------------------------------------------------------------------------
subroutine centerDiff(x, r, nlev, missing, dXdZ)
     integer, intent(in )  ::  nlev
     real(4), intent(in )  ::  x   (nlev)
     real(4), intent(in )  ::  r   (nlev)
     real(4), intent(in )  ::  missing
     real(4), intent(out)  ::  dXdZ(nlev)

     !*********************************
     !****  Initialize
     !*********************************
     dXdz  =  missing

     !*********************************
     !****  Get the non-endpoints 
     !*********************************
     where( x(3:nlev).ne.missing  .and.  x(1:nlev-2).ne.missing  .and.  & 
            r(3:nlev).ne.missing  .and.  r(1:nlev-2).ne.missing  .and.  &
            (r(3:nlev) - r(1:nlev-2)).ne.0                              )
            dXdZ(2:nlev-1)  =  (x(3:nlev) - x(1:nlev-2)) / (r(3:nlev) - r(1:nlev-2))
     endwhere

     !*********************************
     !****  Get the endpoints 
     !*********************************
     if( (r(2)    - r(1)     ).ne.0 )   dXdZ(1)     =  (x(2)    - x(1)     ) / (r(2)    - r(1)     )
     if( (r(nlev) - r(nlev-1)).ne.0 )   dXdZ(nlev)  =  (x(nlev) - x(nlev-1)) / (r(nlev) - r(nlev-1))
end subroutine centerDiff






!----------------------------------------------------------------------------------------------
!
! function: Approximate a smoothed profile lapse rate
!
!----------------------------------------------------------------------------------------------
subroutine smoothed_lapse_rate (theta, pressure, nlev, missing, dTdP)

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  theta(nlev), pressure(nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  dTdP(nlev)

      real(4), dimension(nlev) :: smoothed_theta, smoothed_pressure
      real(4), dimension(nlev) :: temporary_profile, allmissing
      logical, dimension(nlev) :: notmissing

      dTdP               =  missing
      allmissing         =  missing
      notmissing         =  .false.
      temporary_profile  =  missing
      call weighted_avg(theta         , pressure         , nlev, missing, smoothed_theta, smoothed_pressure)
      call centerDiff  (smoothed_theta, smoothed_pressure, nlev, missing, temporary_profile)
      dTdP  =  temporary_profile

      !-----------------------------------------------------------------------------
      !-- Collapse the missing values along the level dimension so that mid-points
      !-- can be calculated using levels on either side of the missing level.
      !-- This would avoid just ignoring the level
      !-- Use the PACK function
      !-----------------------------------------------------------------------------
      ! notmissing  =  .not.( temporary_profile.eq.missing )
      ! dTdP        =  pack ( temporary_profile, notmissing, allmissing)

end subroutine smoothed_lapse_rate






!----------------------------------------------------------------------------------------------
!
! function: Calculate the weighted average using the pressure layer depth as the weight
!           This function is used to smooth high-resolution profiles.  Here Hi-resolution
!           refers to vertical profiles that have layers less than 10hPa thick
!           Reason:  The reason for smoothing the profile is to remove high frequency perturbations
!                    in d(theta)/d(p) which will produce errors when calculating boundary layer height
!                    using the gradient technique
!
!----------------------------------------------------------------------------------------------
subroutine weighted_avg (theta, pressure, nlev, missing, smoothed_theta, smoothed_pressure)

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  theta(nlev), pressure(nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  smoothed_theta(nlev), smoothed_pressure(nlev)

      integer               ::  zz, i0, i1
      real(4)               ::  pressure_layer(nlev), depth, mid_theta(nlev)
      real(4), parameter    ::  threshold_depth = 1e3


      !==========================================================================
      !===  Missing value
      !==========================================================================
      smoothed_theta     =  missing
      smoothed_pressure  =  missing

      !==========================================================================
      !===  Calculate the depth of each layer and mid-point of each theta layer
      !==========================================================================
      call depthPressure    (pressure, nlev, missing, pressure_layer)


      !==========================================================================
      !===  If all layers are greater than the threshold depth then do nothing
      !===  Basically, it is smooth enough
      !==========================================================================
      if( all(pressure_layer.ge.threshold_depth) ) then
         smoothed_theta     =  theta
         smoothed_pressure  =  pressure
         return
      end if

      !==========================================================================
      !===  Calculate the mid-point of each theta layer
      !==========================================================================
      call linear_layer_avg (theta   , nlev, missing, mid_theta )

      !==========================================================================
      !===  Loop over all layers and average over layers when 
      !==========================================================================
      i0  =  1
      do zz = 2,nlev
         i1  =  zz
         depth  =  sum( pressure_layer(i0:i1), mask = pressure_layer(i0:i1).ne.missing )
         if( depth.ge.threshold_depth ) then
            smoothed_theta   (i0:i1)  =  sum( mid_theta   (i0:i1) * pressure_layer(i0:i1) / depth , &
                                              mask = pressure_layer(i0:i1).ne.missing           )
            smoothed_pressure(i0:i1)  =  sum( pressure(i0:i1) * pressure_layer(i0:i1) / depth , &
                                              mask = pressure_layer(i0:i1).ne.missing           )
            i0 = i1
         end if
      end do

end subroutine weighted_avg








!----------------------------------------------------------------------------------------------
!
! function: Add sensible heat to the across the mixed layer
!
!----------------------------------------------------------------------------------------------
subroutine add_sensible_heat (initialProfile, sensible_heat, pbl_depth, height, pblh, dt, nlev, missing, newProfile )

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  initialProfile(nlev), sensible_heat , dt
      real(4), intent(in )  ::  pbl_depth           , height(nlev), pblh
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  newProfile(nlev)
      real(4), parameter    ::  grav = 9.81, cp = 1005.7, cp_g = cp/grav

      newProfile = initialProfile
      where( initialProfile.ne.missing  .and.  initialProfile(1).ne.missing  .and.  height.le.pblh )
          newProfile  =  initialProfile(1) + (dt * (sensible_heat / (pbl_depth*cp_g)))
      endwhere

end subroutine add_sensible_heat






!---------------------------------------------------------------------------------
!
! function: Average over a certain layer Non-linearly
!           Follows equation from Rogers and Yao textbook on Cloud Physics
!
!---------------------------------------------------------------------------------
subroutine middle_layer_avg (scalar, pressure, nlev, missing, middle )

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  scalar(nlev), pressure(nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  middle(nlev)

      integer :: nlev1

      !-------------------------------------------------------------------
      !--- Calculate layer averages (mid-point)
      !-------------------------------------------------------------------
      nlev1  = nlev - 1
      middle = missing
      where( scalar  (2:nlev).ne.missing  .and.  scalar  (:nlev1).ne.missing  .and.  &
             pressure(2:nlev).ne.missing  .and.  pressure(:nlev1).ne.missing         )

         middle(:nlev1)  =  ( (scalar  (2:nlev)*log(pressure(2:nlev))  +  scalar(:nlev1)*log(pressure(:nlev1)))  /  &
                           log(pressure(2:nlev)*    pressure(:nlev1)) )

      end where

end subroutine middle_layer_avg



!---------------------------------------------------------------------------------
!
! function: Average over a certain layer -- Note this is purely linear
!
!---------------------------------------------------------------------------------
subroutine linear_layer_avg (scalar, nlev, missing, middle )

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  scalar(nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  middle(nlev)

      integer :: nlev1

      !-------------------------------------------------------------------
      !--- Calculate layer averages (mid-point)
      !-------------------------------------------------------------------
      nlev1  = nlev - 1
      middle = missing
      where( scalar(2:nlev).ne.missing  .and.  scalar(:nlev1).ne.missing )
         middle(:nlev1)  =  0.5 * (scalar  (2:nlev) + scalar(:nlev1))
      end where

end subroutine linear_layer_avg











!---------------------------------------------------------------------------------
!
! function: Average over a certain layer
!
!---------------------------------------------------------------------------------
subroutine avg_over_layer (incoming, targetDepth, depth, nlev, missing, outgoing )

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  incoming(nlev), depth(nlev)
      real(4), intent(in )  ::  targetDepth
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  outgoing

      real(4)               ::  ngood !# of non-missing values within the desired layer

      outgoing  =  missing
      ngood     =  real(count(incoming.ne.missing .and. depth.le.targetDepth))
      if( ngood.gt.0 ) then
         outgoing  =  sum(incoming, mask=(incoming.ne.missing .and. depth.le.targetDepth)) / ngood 
      else
         outgoing  =  incoming(1)
      end if

end subroutine avg_over_layer




!---------------------------------------------------------------------------------
!
! function: Assign a value over a certain layer
!
!---------------------------------------------------------------------------------
subroutine assign_layer (incoming, value2assign, targetDepth, depth, nlev, missing, outgoing )

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  incoming(nlev), depth  (nlev)
      real(4), intent(in )  ::  targetDepth   , value2assign
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  outgoing(nlev)

      integer               ::  iTarget


      !-----------------------------------------------------------------------------
      !-- Initialize output variables
      !-----------------------------------------------------------------------------
      outgoing  =  missing

      !-----------------------------------------------------------------------------
      !-- Initialize output variables
      !-----------------------------------------------------------------------------
      call maxIndex ( nlev, 1, (incoming.ne.missing .and. depth.le.targetDepth), iTarget )
      outgoing(:iTarget   )  =  value2assign
      outgoing( iTarget+1:)  =  incoming(iTarget+1:) 


end subroutine assign_layer









!---------------------------------------------------------------------------------
! function: put lowest level variable in profile and move all non-missing in the 1st indices
!---------------------------------------------------------------------------------
subroutine packIt(profile, height2m, nlev, nlev1, psfc, press, missing, packed)
      integer, intent(in )  ::  nlev, nlev1
      real(4), intent(in )  ::  missing
      real(4), intent(in )  ::  profile(nlev)
      real(4), intent(in )  ::  height2m
      real(4), intent(in )  ::  press(nlev)
      real(4), intent(in )  ::  psfc
      real(4), intent(out)  ::  packed(nlev1)

      logical              ::  notmissing(nlev1)
      real(4)              ::  allmissing(nlev1)
      real(4)              ::  working   (nlev1)
      real(4)              ::  presswork (nlev1)

      working(2:)    =  profile
      working(1 )    =  height2m
      presswork(2:)  =  press
      presswork(1 )  =  psfc
      allmissing     =  missing
      notmissing     =  .false.

      !-----------------------------------------------------------------------------
      !-- Make sure there are no higher pressure levels than that given by surface pressure
      !-- Important for models that have atmo levels below the surface
      !-----------------------------------------------------------------------------
      where( presswork(2:).ge.psfc ) presswork(2:) = missing

      !-----------------------------------------------------------------------------
      !-- Collapse the missing values along the level dimension so that mid-points
      !-- can be calculated using levels on either side of the missing level.
      !-- This would avoid just ignoring the level
      !-- Use the PACK function
      !-----------------------------------------------------------------------------
      notmissing  =  .not.( presswork.eq.missing .or. working.eq.missing )
      packed      =  pack (working, notmissing, allmissing)
end subroutine packIt




!---------------------------------------------------------------------------------
!
! function: Returns the density ignoring solids and liquids present
!
!---------------------------------------------------------------------------------
subroutine total_density (press, temp, mixing_ratio, nsize, missing, density )

      integer, intent(in )  ::  nsize
      real(4), intent(in )  ::  mixing_ratio(nsize)
      real(4), intent(in )  ::  press  (nsize), temp(nsize)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  density(nsize)

      real(4), parameter        ::  Rd=287.04, ep=0.622
      real(4), dimension(nsize) ::  denom

      density  =  missing
      denom    =  missing
      where( temp.ne.missing  .and.  press.ne.missing  .and.  mixing_ratio.ne.missing )
         denom   =   Rd * (temp) * ((1. + (mixing_ratio/ep))/(1. + mixing_ratio))
      end where
      where( denom.ne.missing  .and.  press.ne.missing  .and.  denom.ne.0 )  density = (press) / denom

end subroutine total_density




!---------------------------------------------------------------------------------
!
! function: Returns the density ignoring solids and liquids present
!
!---------------------------------------------------------------------------------
subroutine totaldensity (press, temp, mixing_ratio, nsize, nlev, missing, density )

      integer, intent(in )  ::  nsize, nlev
      real(4), intent(in )  ::  mixing_ratio(nlev,nsize)
      real(4), intent(in )  ::  press  (nlev,nsize), temp(nlev,nsize)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  density(nlev,nsize)

      real(4), parameter             ::  Rd=287.04, ep=0.622
      real(4), dimension(nlev,nsize) ::  denom

      density  =  missing
      denom    =  missing
      where( temp.ne.missing  .and.  press.ne.missing  .and.  mixing_ratio.ne.missing )
         denom   =   Rd * (temp) * ((1. + (mixing_ratio/ep))/(1. + mixing_ratio))
      end where
      where( denom.ne.missing  .and.  press.ne.missing  .and.  denom.ne.0 )  density = (press) / denom

end subroutine totaldensity





!---------------------------------------------------------------------------------
!
! function: Returns the air pressure using the ideal gas law
!
!---------------------------------------------------------------------------------
subroutine idealGas_P (temp, density, mixQ, nsize, nlev, missing, pressure )

      integer, intent(in )  ::  nsize, nlev
      real(4), intent(in )  ::  temp   (nlev,nsize)
      real(4), intent(in )  ::  density(nlev,nsize), mixQ(nlev,nsize)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  pressure(nlev,nsize)
      real(4), parameter    ::  Rd=287.04
      real(4)               ::  Tvirt(nlev,nsize)

      pressure  =  missing
      Tvirt     =  missing
      call virtualTemp(temp, mixQ, nsize, nlev, missing, Tvirt)
      where( temp .ne.missing  .and.  density.ne.missing  .and.  &
             Tvirt.ne.missing  .and.  density.ne.0               )
          pressure = Rd * Tvirt * density
      endwhere

end subroutine idealGas_P




!---------------------------------------------------------------------------------
! function: Calculates virtual temperature
!---------------------------------------------------------------------------------
subroutine virtualTemp(temp, qhum, nsize, nlev, missing, Tvirt)
     integer, intent(in )  ::  nlev, nsize
     real(4), intent(in )  ::  temp (nlev,nsize)
     real(4), intent(in )  ::  qhum (nlev,nsize)
     real(4), intent(in )  ::  missing
     real(4), intent(out)  ::  Tvirt(nlev,nsize)
 
     Tvirt  =  missing
     where( temp.ne.missing  .and.  qhum.ne.missing )
            Tvirt  =  temp * (1 + 0.61 * qhum)
     endwhere
end subroutine virtualTemp










!---------------------------------------------------------------------------------
!         cpv    = 
!   Calculate latent heat of condensation as a temperature dependent variable
!
!---------------------------------------------------------------------------------
subroutine Latent_heat_of_condensation(T, missing, Lc)

      real(4), intent(in )  ::  missing
      real(4), intent(in )  ::  T 
      real(4), intent(out)  ::  Lc

      real(4), parameter :: cvl   = 4218.     ,  cvv = 1463. !J K-1 kg-1
      real(4), parameter :: Ttrip = 273.16    ,  Rv = 461    !J K-1 kg-1
      real(4), parameter :: Eov   = 2375073.24

      Lc = missing
      if( T.ne.missing )  Lc  =  Eov + Rv*T + (cvv-cvl) * (T-Ttrip)

end subroutine Latent_heat_of_condensation










!---------------------------------------------------------------------------------
! function: Calculate the mid-level of a certain variable using log-pressure averaging
!---------------------------------------------------------------------------------
subroutine midLevel(inputVariable, press, nlev, missing, middle)

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  inputVariable(nlev)
      real(4), intent(in )  ::  press(nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  middle(nlev)

      middle  =  inputVariable
      where( inputVariable(1:nlev-1).ne.missing  .and.  press(1:nlev-1).ne.missing  .and. &
             inputVariable(2:nlev  ).ne.missing  .and.  press(2:nlev  ).ne.missing        )
             middle(2:nlev)    =  ((inputVariable(2:nlev  )*log(press(2:nlev  ))  + &
                                    inputVariable(1:nlev-1)*log(press(1:nlev-1))) / &
                                    log(press(2:nlev)* press(1:nlev-1)))
      endwhere
end subroutine midLevel



!---------------------------------------------------------------------------------
! function: Calculate the mid-level of a certain variable using log-pressure averaging
!---------------------------------------------------------------------------------
subroutine trueMidLevel(inputVariable, press, nlev, missing, middle)

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  inputVariable(nlev)
      real(4), intent(in )  ::  press(nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  middle(nlev)

      middle  =  missing
      where( inputVariable(1:nlev-1).ne.missing  .and.  press(1:nlev-1).ne.missing  .and. &
             inputVariable(2:nlev  ).ne.missing  .and.  press(2:nlev  ).ne.missing        )
             middle(1:nlev-1)    =  ((inputVariable(2:nlev  )*log(press(2:nlev  ))  + &
                                      inputVariable(1:nlev-1)*log(press(1:nlev-1))) / &
                                      log(press(2:nlev)* press(1:nlev-1)))
      endwhere
end subroutine trueMidLevel








!---------------------------------------------------------------------------------
!
! function: Calculate column density to convert a mixing ratio to a column density
!
!---------------------------------------------------------------------------------
subroutine columnDensity( mixing_ratio, dpress, nlev, missing, column_rho )
      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  mixing_ratio(nlev)
      real(4), intent(in )  ::  dpress      (nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  column_rho  (nlev)
      real(4), parameter    ::  grav = 9.81

      column_rho = missing
      where( dpress.ne.missing .and. mixing_ratio.ne.missing )
              column_rho  =  mixing_ratio  * dpress/grav
      endwhere
end subroutine columnDensity



!---------------------------------------------------------------------------------
!
! function: Sum over a particular variable;  This is a cummulative sum
!
!---------------------------------------------------------------------------------
subroutine cummulative_sum( incoming, nlev, missing, outgoing )
      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  incoming(nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  outgoing(nlev)
      integer               ::  zz

      outgoing  =  incoming
      do zz = 2,nlev
         if( outgoing(zz).ne.missing .and. outgoing(zz-1).ne.missing ) then
             outgoing(zz)  =  outgoing(zz-1) + outgoing(zz)
         end if
      end do

end subroutine cummulative_sum





!---------------------------------------------------------------------------------
!
! function: Calculate column density to convert a mixing ratio to a column density
!           
!---------------------------------------------------------------------------------
subroutine columnDensityNoMix( dpress, nlev, missing, column_rho )
      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  dpress    (nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  column_rho(nlev)
      real(4), parameter    ::  grav = 9.81

      column_rho = missing
      where( dpress.ne.missing ) column_rho  =  dpress/grav

end subroutine columnDensityNoMix









!---------------------------------------------------------------------------------
!
! function: Calculate column density to convert a mixing ratio to a column density
!
!---------------------------------------------------------------------------------
subroutine columnDensityNoMix_Cummulative( dpress, nlev, missing, column_rho )
      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  dpress     (nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  column_rho (nlev)

      real(4)               ::  column_rho0(nlev)

      call columnDensityNoMix( dpress     , nlev, missing, column_rho0 )
      call cummulative_sum   ( column_rho0, nlev, missing, column_rho  )

end subroutine columnDensityNoMix_Cummulative









!---------------------------------------------------------------------------------
!
! function: Calculate column density to convert a mixing ratio to a column density
!
!---------------------------------------------------------------------------------
subroutine columnDensity_Cummulative( mixing_ratio, dpress, nlev, missing, column_rho )
      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  mixing_ratio(nlev)
      real(4), intent(in )  ::  dpress      (nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  column_rho  (nlev)

      real(4)               ::  column_rho0(nlev)

      call columnDensity  ( mixing_ratio, dpress, nlev, missing, column_rho0 )
      call cummulative_sum( column_rho0, nlev, missing, column_rho )

end subroutine columnDensity_Cummulative







!---------------------------------------------------------------------------------
! function: Calculate potential temperature
!---------------------------------------------------------------------------------
subroutine potentialTemperature(temperature, pressure, nlev, missing, theta)
     integer, intent(in )  ::  nlev
     real(4), intent(in )  ::  temperature(nlev)
     real(4), intent(in )  ::  pressure   (nlev)
     real(4), intent(in )  ::  missing
     real(4), intent(out)  ::  theta(nlev)
     real(4), parameter    ::  p_ref = 1e5, R_cp=287.04/1005.7

     theta  =  missing
     where(temperature.ne.missing  .and.  pressure.ne.missing)  
          theta  =  temperature * ((p_ref/pressure)**(R_cp))
     endwhere
end subroutine potentialTemperature




!---------------------------------------------------------------------------------
! function: Calculate potential temperature
!---------------------------------------------------------------------------------
subroutine calculate_temperature(theta, pressure, nlev, missing, temperature)
     integer, intent(in )  ::  nlev
     real(4), intent(in )  ::  pressure(nlev)
     real(4), intent(in )  ::  missing
     real(4), intent(in )  ::  theta(nlev)
     real(4), intent(out)  ::  temperature(nlev)
     real(4), parameter    ::  p_ref = 1e5, R_cp=287.04/1005.7

     temperature  =  missing
     where(theta.ne.missing  .and.  pressure.ne.missing)  
          temperature  =  theta / ((p_ref/pressure)**(R_cp))
     endwhere

end subroutine calculate_temperature








!---------------------------------------------------------------------------------
! function: Calculate potential temperature
!---------------------------------------------------------------------------------
subroutine potentialTemperature1D(temperature, pressure, missing, theta)
     real(4), intent(in )  ::  temperature
     real(4), intent(in )  ::  pressure   
     real(4), intent(in )  ::  missing
     real(4), intent(out)  ::  theta
     real(4), parameter    ::  p_ref = 1e5, R_cp=287.04/1005.7

     theta  =  missing
     if(temperature.ne.missing  .and.  pressure.ne.missing)  then
          theta  =  temperature * ((p_ref/pressure)**(R_cp))
     end if
end subroutine potentialTemperature1D




!---------------------------------------------------------------------------------
! function: Calculate near surface pressure using ideal gas law
!---------------------------------------------------------------------------------
subroutine surfacePressure(t, p, q, h, psfc)
     real(4), intent(in )  ::  t
     real(4), intent(in )  ::  p
     real(4), intent(in )  ::  q
     real(4), intent(in )  ::  h
     real(4), intent(out)  ::  psfc
     real(4), parameter    ::  grav=9.81, Rd=287.04, ep=0.622
     psfc = (p / (Rd * t * ((1. + (q/ep)) / (1. + q)) )) * grav * h
end subroutine surfacePressure



!---------------------------------------------------------------------------------
! function: Calculates the difference between two layers (except the 1st)
!---------------------------------------------------------------------------------
subroutine depthPressure(press, nlev, missing, pressureDiff)
     integer, intent(in )  ::  nlev
     real(4), intent(in )  ::  press(nlev)
     real(4), intent(in )  ::  missing
     real(4), intent(out)  ::  pressureDiff(nlev)
     pressureDiff  =  missing
     where( press(1:nlev-1).ne.missing  .and.  press(2:nlev).ne.missing )
            pressureDiff(1:nlev-1)  =  press(1:nlev-1) - press(2:nlev)
     endwhere
end subroutine depthPressure




!---------------------------------------------------------------------------------
! function: Calculates virtual potential temperature
!---------------------------------------------------------------------------------
subroutine virtualPotentialTemp(temp, qhum, press, nlev, missing, thetaV)
     integer, intent(in )  ::  nlev
     real(4), intent(in )  ::  temp  (nlev)
     real(4), intent(in )  ::  qhum  (nlev)
     real(4), intent(in )  ::  press (nlev)
     real(4), intent(in )  ::  missing
     real(4), intent(out)  ::  thetaV(nlev)
     
     real(4)               ::  Tvirt(nlev)
     call virtualTemperature  (temp, qhum, nlev, missing, Tvirt)
     call potentialTemperature(Tvirt, press, nlev, missing, thetaV)
end subroutine virtualPotentialTemp


!---------------------------------------------------------------------------------
! function: Calculates virtual temperature
!---------------------------------------------------------------------------------
subroutine virtualTemperature(temp, qhum, nlev, missing, Tvirt)
     integer, intent(in )  ::  nlev
     real(4), intent(in )  ::  temp  (nlev)
     real(4), intent(in )  ::  qhum  (nlev)
     real(4), intent(in )  ::  missing
     real(4), intent(out)  ::  Tvirt(nlev)
 
     Tvirt  =  missing
     where( temp.ne.missing  .and.  qhum.ne.missing )
            Tvirt  =  temp * (1 + 0.61 * qhum)
     endwhere
end subroutine virtualTemperature



!---------------------------------------------------------------------------------
! function: Calculates the difference between two layers includes surface layer
!---------------------------------------------------------------------------------
subroutine layerDepth(press, t, q, p, h, nlev, missing, Depth)
     integer, intent(in )  ::  nlev
     real(4), intent(in )  ::  press(nlev)
     real(4), intent(in )  ::  t
     real(4), intent(in )  ::  p
     real(4), intent(in )  ::  q
     real(4), intent(in )  ::  h
     real(4), intent(in )  ::  missing
     real(4), intent(out)  ::  Depth(nlev)

     call depthPressure(press, nlev, missing, Depth)
     Depth(1)  =  press(2) - press(1)
     if( Depth(1).le.0 ) then
         Depth(1)  =  1.    !*** set to 1 Pa because the h2m is likley zero
     else
         call surfacePressure(t, p, q, h, Depth(1))
     end if
end subroutine layerDepth



!---------------------------------------------------------------------------------
! function: Calculates the log of a variable while ignoring the missing values
!---------------------------------------------------------------------------------
subroutine logMiss(variableToBeLogged, dimsize, missing, logVal)
     integer, intent(in )  ::  dimsize
     real(4), intent(in )  ::  missing
     real(4), intent(in )  ::  variableToBeLogged(dimsize)
     real(4), intent(out)  ::  logVal(dimsize)
     logVal  =  missing
     where( variableToBeLogged.ne.missing ) logVal = log(variableToBeLogged)
end subroutine logMiss





!-----------------------------------------------------------------------------
! function: Inject latent heat flux across the mixed layer [kg/m2]
!-----------------------------------------------------------------------------
subroutine inject_moisture( tracer_density, injection_amount, mixLevel, press, dp, layer_depth, nlev, missing, new_density )

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  missing
      real(4), intent(in )  ::  injection_amount      !*** amount of tracer to be injected [kg/m2]
      real(4), intent(in )  ::  mixLevel              !*** Pressure of level to mix down from [Pa]
      real(4), intent(in )  ::  tracer_density(nlev)  !*** tracer to be mixed assuming units of [kg/kg]
      real(4), intent(in )  ::  dp(nlev)              !*** depth of each pressure level [Pa]
      real(4), intent(in )  ::  press (nlev)          !*** pressure levels [Pa]
      real(4), intent(in )  ::  layer_depth           !*** depth of boundary layer psfc - pblp [Pa]
      real(4), intent(out)  ::  new_density(nlev)     !*** return mixed tracer profile [kg/kg]

     new_density  =  tracer_density
     where( press.gt.mixLevel .and. press.ne.missing .and. new_density.ne.missing )
         new_density  =  new_density + (injection_amount * (dp/layer_depth))
     end where

end subroutine inject_moisture




!-----------------------------------------------------------------------------
! function: Inject latent heat flux across the mixed layer [kg/m2]
!-----------------------------------------------------------------------------
subroutine inject_and_mix( mixLevel, injection, press, tracer, dpress, nlev, missing, MixedProfile )

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  missing
      real(4), intent(in )  ::  injection           !*** amount of tracer to be injected [kg/m2]
      real(4), intent(in )  ::  mixLevel            !*** Pressure of level to mix down from [Pa]
      real(4), intent(in )  ::  tracer(nlev)        !*** tracer to be mixed assuming units of [kg/kg]
      real(4), intent(in )  ::  dpress(nlev)        !*** depth of each pressure level [Pa]
      real(4), intent(in )  ::  press (nlev)        !*** pressure levels [Pa]
      real(4), intent(out)  ::  MixedProfile(nlev)  !*** return mixed tracer profile [kg/kg]

      integer  ::  zz, ilayer_top
      real(4)  ::  ilevel(nlev), layer_depth
      real(4)  ::  tracer_density(nlev), air_density(nlev), tracer_density_w_flux(nlev)
      real(4)  ::  tracer_mix          , air_mix



      !-----------------------------------------------------
      !--------  Initialize output variable
      !--------  don't do anything if there is no flux
      !-----------------------------------------------------
      MixedProfile  =  tracer
      if( injection.le.0 ) return


      !-----------------------------------------------------
      !--------  Return column density
      !-----------------------------------------------------
      call columnDensity     (tracer,dpress,nlev,missing,tracer_density)
      call columnDensityNoMix(       dpress,nlev,missing,   air_density)


      !-----------------------------------------------------
      !--------  index of levels
      !-----------------------------------------------------
      do zz = 1,nlev
         ilevel(zz)  =  real(zz)
      end do


      !-----------------------------------------------------
      !--------  Are there levels within the mixed layer?
      !-----------------------------------------------------
      call maxIndex ( nlev, 1, (press.gt.mixLevel  .and.  press.ne.missing)        , ilayer_top   )
      call sumIt    ( nlev, 1, (press.gt.mixLevel  .and.  press.ne.missing), dpress, layer_depth  )


      !*********************************************************************
      !**** Inject only when latent heat is greater than zero
      !**** There are two cases involved:
      !****
      !**** 1) When the boundary layer is above the 1st resolved level
      !****    in this case inject across all levels below the PBL that is
      !****    scaled by the depth of each layer
      !****
      !**** 2) When the boundary layer is WITHIN the 1st resolved level
      !****    in this case inject across the first two levels to avoid the 
      !****    saturating the level between 2-meter and level one right after
      !****    sunrise where you have some ET but weak to no mixing
      !****
      !*********************************************************************
      !--------------------------
      !---  Inject moisture
      !--------------------------
      call inject_moisture( tracer_density, injection, mixLevel, press, dpress, layer_depth, nlev, missing, &
                            tracer_density_w_flux )


      !--------------------------------------------
      !---  Sum over layer and return mixing ratio
      !--------------------------------------------
      tracer_mix  =  sum( tracer_density_w_flux(1:ilayer_top) )
      air_mix     =  sum(           air_density(1:ilayer_top) )
      if( tracer_mix.ne.missing  .and.  air_mix.ne.missing  .and.  air_mix.ne.0 ) then
         MixedProfile(:ilayer_top)  =  tracer_mix / air_mix
      end if


end subroutine inject_and_mix



!---------------------------------------------------------------------------------
!
! function: Calculate integrated column energy [J/m2]
!
!---------------------------------------------------------------------------------
subroutine column_energy (temperature, pressure, nlev, missing, internal_energy )

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  temperature(nlev), pressure(nlev)
      real(4), intent(in )  ::  missing
      real(4), intent(out)  ::  internal_energy

      integer               ::  nlev1
      real(4), parameter    ::  grav = 9.81, cp = 1005.7, cp_g = cp/grav

      !-------------------------------------------------------------------
      !--- Calculate layer averages (mid-point)
      !-------------------------------------------------------------------
      internal_energy  =  missing
      nlev1            =  nlev - 1

      internal_energy  =  sum( 0.5  * (temperature(:nlev1) + temperature(2:nlev)) *                                 & 
                               cp_g * (pressure   (:nlev1) - pressure   (2:nlev))  ,                                &
                               mask= (temperature(:nlev1).ne.missing .and. temperature(2:nlev).ne.missing  .and.    &
                                      pressure   (:nlev1).ne.missing .and. pressure   (2:nlev).ne.missing         ) )


end subroutine column_energy







!-----------------------------------------------------------------------------
! function: Calculate saturation specific humidity [kg/kg]
!-----------------------------------------------------------------------------
subroutine saturationHumidity( t, p, nlev, missing, qsat )
      integer, intent(in ) :: nlev
      real(4), intent(in ) :: t(nlev)
      real(4), intent(in ) :: p(nlev)
      real(4), intent(in ) :: missing
      real(4), intent(out) :: qsat(nlev)

      real(4), parameter   ::  by100 = 1e2
      real(4), parameter   ::  t0=273.15, ep=0.622, es0=6.11, a=17.269, b=35.86
      real(4), parameter   ::  onemep=1.0 - ep

      qsat  =  missing
      where( t.ne.missing .and. p.ne.missing )
         qsat  =  by100*0.01 *(ep* (es0*exp((a*( t-t0))/( t-b))) ) /  &
                   ((p/1e2)-onemep*(es0*exp((a*( t-t0))/( t-b))))
      endwhere
      where( qsat.ne.missing ) 
         qsat  =  qsat/(1.+qsat)
      endwhere
end subroutine saturationHumidity



!-----------------------------------------------------------------------------
! function: Calculate saturation specific humidity [kg/kg]
!-----------------------------------------------------------------------------
subroutine specificHumidity( t, p, relh, nlev, missing, qhum )
      integer, intent(in ) :: nlev
      real(4), intent(in ) :: t   (nlev)
      real(4), intent(in ) :: p   (nlev)
      real(4), intent(in ) :: relh(nlev)
      real(4), intent(in ) :: missing
      real(4), intent(out) :: qhum(nlev)

      real(4), parameter   :: by100 = 1e2
      real(4)              :: qsat(nlev), mixRatio(nlev)

      qhum  =  missing
      call saturationHumidity( t, p, nlev, missing, qsat )
      where( qsat    .ne.missing .and. relh.ne.missing )  mixRatio  =  qsat * (relh/by100)
      where( mixRatio.ne.missing )                        qhum      =  mixRatio / (1.0 + mixRatio)

end subroutine specificHumidity






!-----------------------------------------------------------------------------
! function: Fill above a certain value
!-----------------------------------------------------------------------------
subroutine maskAbove( applyMaskTo, maskValue, limitingVariable, fillValue, nlev, maskedAbove )
      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  applyMaskTo     (nlev)
      real(4), intent(in )  ::  maskValue       
      real(4), intent(in )  ::  limitingVariable(nlev)
      real(4), intent(in )  ::  fillValue       
      real(4), intent(out)  ::  maskedAbove     (nlev)
      maskedAbove  =  applyMaskTo
      where( limitingVariable.ge.fillValue )  maskedAbove  =  maskValue
end subroutine maskAbove






!-----------------------------------------------------------------------------
! function: use the fortran count function
!-----------------------------------------------------------------------------
subroutine sumIt ( nlev, dims, conditional, var_to_sum, This_Sum )
      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  var_to_sum (nlev)
      logical, intent(in )  ::  conditional(nlev)
      integer, intent(in )  ::  dims
      real(4), intent(out)  ::  This_Sum
     
      This_Sum  =  sum( var_to_sum, mask = conditional, dim = dims )
end subroutine sumIt







!-----------------------------------------------------------------------------
! function: use the fortran count function
!-----------------------------------------------------------------------------
subroutine countPlus ( nlev, conditional, dims, ThisCount )
      integer, intent(in )  ::  nlev
      logical, intent(in )  ::  conditional(nlev)
      integer, intent(in )  ::  dims
      integer, intent(out)  ::  ThisCount

      ThisCount  =  count( MASK = conditional, DIM = dims )
end subroutine countPlus






!-----------------------------------------------------------------------------
! function: use the fortran maxloc function
!-----------------------------------------------------------------------------
subroutine minIndex ( nlev, dims, mymask, ThisIndex )
      integer, intent(in )  ::  nlev
      logical, intent(in )  ::  mymask(nlev)
      integer, intent(in )  ::  dims
      integer, intent(out)  ::  ThisIndex

      real(4)               ::  ilevels(nlev)
      integer               ::  zz

      !-----------------------------------------------------
      !--------  index of levels
      !-----------------------------------------------------
      ThisIndex  =  -1
      do zz = 1,nlev
         ilevels(zz)  =  real(zz)
      end do
      ThisIndex  =  minloc( ilevels, DIM = dims, MASK = mymask)

end subroutine minIndex



!-----------------------------------------------------------------------------
! function: use the fortran maxloc function
!-----------------------------------------------------------------------------
subroutine maxIndex ( nlev, dims, mymask, ThisIndex )
      integer, intent(in )  ::  nlev
      logical, intent(in )  ::  mymask(nlev)
      integer, intent(in )  ::  dims
      integer, intent(out)  ::  ThisIndex

      real(4)               ::  ilevels(nlev)
      integer               ::  zz

      !-----------------------------------------------------
      !--------  index of levels
      !-----------------------------------------------------
      ThisIndex  =  -1
      do zz = 1,nlev
         ilevels(zz)  =  real(zz)
      end do
      ThisIndex  =  maxloc( ilevels, DIM = dims, MASK = mymask)

end subroutine maxIndex






!-----------------------------------------------------------------------------
!--- Subroutine:  Check to see if first level is satured; (Foggy scenario)
!--- If so then check the second and third layers to see if fog will dissipate
!--- If the 2nd and/or 3rd are not saturated then recalc CONVECTIVE saturation
!--- transition level
!-----------------------------------------------------------------------------
subroutine checkForFog( ilower, iupper, psfc, pblp, press, qdef, nlev, missing, i_unsat, i_sat )

      integer, intent(out)  ::  i_sat
      integer, intent(out)  ::  i_unsat
      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  psfc
      real(4), intent(in )  ::  pblp
      real(4), intent(in )  ::  press(nlev)
      real(4), intent(in )  ::  qdef (nlev)
      real(4), intent(in )  ::  missing
      integer, intent(in )  ::  iupper
      integer, intent(in )  ::  ilower

      integer                 ::  zz, ilevels(nlev), cc



      !---------------------------------------------------------------
      !--- Create the loop index for minloc calculation
      !---------------------------------------------------------------
      do zz = 1,nlev
         ilevels(zz)  =  zz
      end do

      !---------------------------------------------------------------
      !--- Loop over levels and see if the fog layer can be eroded
      !---------------------------------------------------------------
      i_unsat = ilower
      i_sat   = iupper
      cc      = 0
      do zz=2,nlev-1
         !**** make sure initiation level is below 150 hPa above the ground to ensure it is actually fog
         if( (psfc-press(zz))/1e2.gt.150   .or.  (press(zz).le.pblp .and. (psfc-pblp)/1e2.gt.150) ) then
            i_sat   =  0
            i_unsat = -1
            exit
         end if

         !**** If within the 150 hPa layer above the ground then try to erode the fog layer first
         !**** to determine the convective initiation layer
         i_sat  =  minloc( ilevels(zz:), DIM = 1, MASK = qdef(zz:).ne.missing .and. qdef(zz:).le.0 )
         cc     =  cc + 1

         !**** If still saturated then cycle
         if( qdef(zz).le.0 ) cycle
         i_sat  =  i_sat + cc
         exit
      end do
      i_unsat   =  i_sat - 1


end subroutine checkForFog





!-----------------------------------------------------------------------------
!--- Subroutine: Find level of transition from positive to negative number
!-----------------------------------------------------------------------------
subroutine transition( profile, nlev, missing, i_unsat, i_sat )

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  missing
      real(4), intent(in )  ::  profile(nlev)
      integer, intent(out)  ::  i_unsat, i_sat
      integer               ::  zz,ilevels(nlev)
      integer               ::  notmissing, lessthan0



      !--------------------------------------------
      !--- Initialize
      !--------------------------------------------
      i_sat   = 0
      i_unsat = 0


      !---------------------------------------------------------------
      !--- -- If all levels are saturated then just exit
      !--- Create the loop index for minloc calculation
      !--- Find new maximum level that isn't missing
      !---------------------------------------------------------------
      notmissing  =  count(profile.ne.missing)
      lessthan0   =  count(profile.ne.missing .and. profile.le.0)
      if( notmissing.eq.lessthan0 ) then
         i_unsat = i_sat - 1
         return
      end if



      !--------------------------------------------
      !--- Loop to find sign change
      !--------------------------------------------
      do zz=1,nlev
         if( profile(zz).ne.missing .and. profile(zz).le.0 ) then
            i_sat = zz
            exit
         end if
      end do
      i_unsat = i_sat - 1



      !---------------------------------------------------------------
      !--- If all NOT saturated 
      !--- Then set to highest non-missing level
      !---------------------------------------------------------------
      lessthan0   =  count(profile.ne.missing .and. profile.gt.0)
      if( notmissing.eq.lessthan0 ) then
         do zz = 1,nlev
            ilevels(zz)  =  zz
         end do
         i_sat = maxloc( ilevels, MASK = profile.ne.missing, DIM = 1 )
         i_unsat = i_sat - 1
      end if


end subroutine transition




!-----------------------------------------------------------------------------
!--- Subroutine:  Get indicies of transition and 
!--- Check to see if first level is satured; (Foggy scenario)
!-----------------------------------------------------------------------------
subroutine transitionCheck( profile, psfc, pblp, press, nlev, missing, i_unsat, i_sat )

      integer, intent(in )  ::  nlev
      real(4), intent(in )  ::  missing
      real(4), intent(in )  ::  profile(nlev)
      real(4), intent(in )  ::  press(nlev)
      real(4), intent(in )  ::  pblp, psfc
      integer, intent(out)  ::  i_unsat, i_sat


      !-----------------------------------------------------------------------------
      !--- Check to see if first level is satured; (Foggy scenario)
      !--- If so then check the second and third layers to see if fog will dissipate
      !--- If the 2nd and/or 3rd are not saturated then recalc CONVECTIVE saturation
      !--- transition level
      !-----------------------------------------------------------------------------
      call transition ( profile, nlev, missing, i_unsat, i_sat )
      if( i_sat.le.1 ) then
         call checkForFog( i_sat, i_unsat, psfc, pblp, press, profile, nlev, missing, i_unsat, i_sat )
      end if

end subroutine transitionCheck



!-----------------------------------------------------------------------------
!---   Subroutine:  Calculate slope of each variable to find the      
!---   y-intercept for each variable;                    
!---   Meaning locate the two data points surrounding    
!---   the sign change in qdef and linearly interpolate  
!---   to find the "zero point"                          
!-----------------------------------------------------------------------------
subroutine findIntercept( nlev, iupper, ilower, y, x, missing, interceptValue )
     integer, intent(in )  ::  nlev
     integer, intent(in )  ::  iupper, ilower
     real(4), intent(in )  ::  y(nlev), x(nlev)
     real(4), intent(in )  ::  missing
     real(4), intent(out)  ::  interceptValue
     real(4)               ::  y0, y1, x0, x1

     !-----------------------------------------------------------------------------
     !--- Get the upper and lower bounds for each variable to be calc'd at the BCL
     !-----------------------------------------------------------------------------
     x0  =  x(ilower)
     x1  =  x(iupper)
     y0  =  y(ilower)
     y1  =  y(iupper)
     if(x1.ne.missing .and. x0.ne.missing .and.  &
        y1.ne.missing .and. y0.ne.missing        )  then

           interceptValue   =   ( x1 - ((x1-x0)/(y1-y0))*y1 )
     end if

end subroutine findIntercept






!---------------------------------------------------------------------------------
!
! subroutines: Estimates the PBL height by using the lowest level potential temperature
!              to determine where it intersects the environmental profile of Theta
!              The lowest level in the profile is assumed to be the source of buoyancy
!
!---------------------------------------------------------------------------------
subroutine pblHeat ( nlev , missing , pot_k, press, height, PBLP, PBLT, PBLH )

   implicit none
!
! Input/Output Variables
!
   integer, intent(in   )                   ::  nlev        ! *** # of atmospheric levels
   real(4), intent(in   )                   ::  missing     ! *** Missing values
   real(4), intent(in   ), dimension(nlev)  ::  pot_k       ! *** Potential Temperature (level), [K]
   real(4), intent(in   ), dimension(nlev)  ::  press       ! *** Pressure (level) [Pa]
   real(4), intent(in   ), dimension(nlev)  ::  height      ! *** Height (level) [m]
   real(4), intent(out  )                   ::  PBLP        ! *** pressure of boundary layer height [Pa]
   real(4), intent(out  )                   ::  PBLT        ! *** temperature of boundary layer height [K]
   real(4), intent(out  )                   ::  PBLH        ! *** height of boundary layer height [m]
!
! Local variables
!
   real(4), parameter         ::  p_ref = 1e5 , Lv=2.5e6 , cp=1005.7, R_cp=287.04/1005.7
   real(4)                    ::  pot2m
   real(4), dimension(nlev)   ::  pot_diff

   real(4)                    ::  p_up, p_lo, t_up, t_lo, h_up, h_lo
   integer                    ::  i_nobuoy, i_buoy, num_nobuoy, num_buoy



!-----------------------------------------------------------------------------

      !-----------------------------------------------------------------------------
      !-- Initialize output variables
      !-----------------------------------------------------------------------------
      PBLP  =  missing
      PBLT  =  missing
      PBLH  =  missing


      !-----------------------------------------------------------------------------
      !-- Check input arrays to make sure all are NOT missing
      !-----------------------------------------------------------------------------
      if( all(press.eq.missing)  .or.  all(height.eq.missing)  .or.  all(pot_k.eq.missing) )  return


      !-----------------------------------------------------------------------------
      !-- Initialize middle level variables
      !-----------------------------------------------------------------------------
      pot_diff   =  missing
      pot2m      =  pot_k(1)


      !---------------------------------------------------
      !--------
      !--------  Find pressure level and mixed specific
      !--------  humidity deficit given a potential temperature
      !--------
      !---------------------------------------------------
      !---------------------------------------------------
      !-- Calculate difference between reference pot. temp. (K)
      !-- Note the 0.0001 being added to pot2m
      !-- this is to ensure there is a level above the lowest
      !-- level.  It does not influence any other calculations
      !---------------------------------------------------
      where( pot_k.ne.missing .and. press.ne.missing )   pot_diff  =  pot2m  -  pot_k


      !-----------------------------------------------------------------------------
      !----- Find the point where the sign first turns negative from the ground up
      !-----------------------------------------------------------------------------
      !*** Highest buoyant level ***
      call countPlus ( nlev, (pot_diff.ne.missing .and. pot_diff.ge.0), 1, num_buoy )
      if( num_buoy.gt.0 ) then
         call maxIndex ( nlev, 1, (pot_diff.ne.missing .and. pot_diff.ge.0), i_buoy )
      end if

      !*** Lowest negatively buoyant level ***
      call countPlus ( nlev, (pot_diff.ne.missing .and. pot_diff.lt.0), 1, num_nobuoy )
      if( num_nobuoy.gt.0 ) then
         call minIndex ( nlev, 1, (pot_diff.ne.missing .and. pot_diff.lt.0), i_nobuoy )
      end if


      !-----------------------------------------------------------------------------
      !--- All levels are not buoyant then define PBL as average between the 1st and 2nd layer
      !--- This is done for nocturnal boundary layers which are ill-defined
      !-----------------------------------------------------------------------------
      if( num_buoy.eq.0  .or. (i_buoy.eq.1  .and.  i_nobuoy.eq.2) ) then

          i_nobuoy   =  2
          i_buoy     =  1
          
          p_up  =  press (i_nobuoy)
          p_lo  =  press (i_buoy)

          t_up  =  pot_k (i_nobuoy)
          t_lo  =  pot_k (i_buoy)

          h_up  =  height(i_nobuoy)
          h_lo  =  height(i_buoy)

          !-----------------------------------------------------------------------------
          !--- Average the adjacent to get level in between (linear assumption of course)
          !-----------------------------------------------------------------------------
          PBLP  =  (p_up + p_lo) / 2.0
          PBLT  =  (t_up + t_lo) / 2.0 
          PBLH  =  (h_up + h_lo) / 2.0 


          !-----------------------------------------------------------------------------
          !--- Make sure the PBL is not crazy shallow according to buoyancy alone
          !--- If it is less than 100m than perform calc for next layers
          !--- There is still no great way to avoid the very stable boundary layer conditions
          !--- in this case, especially with such low-res veritical data (as in NARR)
          !-----------------------------------------------------------------------------
          if(PBLH.lt.100) then
             i_nobuoy   =  3
             i_buoy     =  2
             p_up  =  press (i_nobuoy)
             p_lo  =  press (i_buoy)
             t_up  =  pot_k (i_nobuoy)
             t_lo  =  pot_k (i_buoy)
             h_up  =  height(i_nobuoy)
             h_lo  =  height(i_buoy)
             !-----------------------------------------------------------------------------
             !--- Average the adjacent to get level in between (linear assumption of course)
             !-----------------------------------------------------------------------------
             PBLP  =  (p_up + p_lo) / 2.0
             PBLT  =  (t_up + t_lo) / 2.0 
             PBLH  =  (h_up + h_lo) / 2.0 
          end if

      else

          !-----------------------------------------------------------------------------
          !--- Average the adjacent to get level in between (linear assumption of course)
          !-----------------------------------------------------------------------------
          call findIntercept( nlev, i_nobuoy, i_buoy, pot_diff, press , missing, PBLP )
          call findIntercept( nlev, i_nobuoy, i_buoy, pot_diff, pot_k , missing, PBLT )
          call findIntercept( nlev, i_nobuoy, i_buoy, pot_diff, height, missing, PBLH )

          !-----------------------------------------------------------------------------
          !--- Make sure the PBL is not crazy shallow according to buoyancy alone
          !--- If it is less than 100m than perform calc for next layers
          !--- There is still no great way to avoid the very stable boundary layer conditions
          !--- in this case, especially with such low-res veritical data (as in NARR)
          !-----------------------------------------------------------------------------
          if(PBLH.lt.100) then
             p_up  =  press (i_nobuoy)
             p_lo  =  press (i_buoy)
             t_up  =  pot_k (i_nobuoy)
             t_lo  =  pot_k (i_buoy)
             h_up  =  height(i_nobuoy)
             h_lo  =  height(i_buoy)
             !-----------------------------------------------------------------------------
             !--- Average the adjacent to get level in between (linear assumption of course)
             !-----------------------------------------------------------------------------
             PBLP  =  (p_up + p_lo) / 2.0
             PBLT  =  (t_up + t_lo) / 2.0 
             PBLH  =  (h_up + h_lo) / 2.0 
          end if


      end if



      return

end subroutine pblHeat



!---------------------------------------------------------------------------------
!
! subroutines:  calculates buoyant condensation level and basic variables (THETA_BM; TDEF)
!
!---------------------------------------------------------------------------------
subroutine hcfcalc ( nlev, missing, tmp_k, press, qhum, hgt, pblp, TBM, TDEF, BCLP, QBCL )

   implicit none
!
! Input/Output Variables
!
   integer, intent(in   )                  ::  nlev        ! *** # of atmospheric levels
   real(4), intent(in   )                  ::  missing     ! *** Missing values
   real(4), intent(in   ), dimension(nlev) ::  tmp_k       ! *** Temperature (level), [K]
   real(4), intent(in   ), dimension(nlev) ::  hgt         ! *** Geometric Height above ground (level) [m]
   real(4), intent(in   ), dimension(nlev) ::  qhum        ! *** Specific Humidity (level) [kg/kg]
   real(4), intent(in   ), dimension(nlev) ::  press       ! *** Pressure (level) [Pa]
   real(4), intent(in   )                  ::  pblp        ! *** pressure of boundary layer [Pa]
   real(4), intent(out  )                  ::  TBM         ! *** buoyant mixing pot. temp (convective threshold) [K]
   real(4), intent(out  )                  ::  TDEF        ! *** potential temperature deficit need to initiate [K]
   real(4), intent(out  )                  ::  BCLP        ! *** pressure of the buoyant condensation level [Pa]
   real(4), intent(out  )                  ::  QBCL        ! *** specific humidity at the buoyant condensation level [kg/kg]
!
! Local variables
!
   real(4), parameter        ::  p_ref = 1e5 , Lv=2.5e6 , cp=1005.7, R_cp=287.04/1005.7
   real(4), dimension(nlev)  ::  rhoh
   real(4), dimension(nlev)  ::  qdef, qmix, qsat, dpress, logp, pot_k

   real(4)                   ::  pot2m
   integer                   ::  i_unsat, i_sat

!-----------------------------------------------------------------------------

      !-----------------------------------------------------------------------------
      !-- Initialize output variables
      !-----------------------------------------------------------------------------
      TBM   =   missing
      BCLP  =   missing
      TDEF  =   missing
      QBCL  =   missing 

      
      !-----------------------------------------------------------------------------
      !-- Check input arrays to make sure all are NOT missing
      !-----------------------------------------------------------------------------
      if( all(tmp_k.eq.missing)  .or.  all(press.eq.missing) .or.  &
          all(hgt  .eq.missing)  .or.  all(qhum .eq.missing) )  return


      !-----------------------------------------------------------------------------
      !-- Initialize middle level variables
      !-----------------------------------------------------------------------------
      qdef    =  missing
      dpress  =  missing
      logp    =  missing
      qmix    =  missing
      qsat    =  missing


      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate column potential temperature (K)
      !--------
      !-----------------------------------------------------------------------------
      call potentialTemperature(tmp_k,press,nlev,missing,pot_k)



      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate pressure difference of each layer
      !--------
      !-----------------------------------------------------------------------------
!      call layerDepth(press, tmp_k(1), qhum(1), press(1), hgt(1), nlev, missing, dpress)
      call depthPressure(press, nlev, missing, dpress)


      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate log pressure to linearize it for slope calculation
      !--------
      !-----------------------------------------------------------------------------
      call logMiss(press, nlev, missing, logp)



      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate column density [kg/m2]
      !--------
      !-----------------------------------------------------------------------------
      call columnDensityNoMix_Cummulative( dpress, nlev, missing, rhoh )



      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate column density of water vapor [kg/m2] AKA precipitable water
      !--------
      !-----------------------------------------------------------------------------
      call columnDensity_Cummulative( qhum, dpress, nlev, missing, qmix )



      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate saturation specific humidity at each level [kg/kg]
      !--------
      !-----------------------------------------------------------------------------
      call saturationHumidity( tmp_k, press, nlev, missing, qsat )



      !-----------------------------------------------------------------------------
      !--------
      !--------  Calculate specific humidity deficit [kg/kg]
      !--------
      !-----------------------------------------------------------------------------
      where( qmix.ne.missing .and. rhoh.ne.missing )
         qmix  =  qmix / rhoh
      endwhere


      !-----------------------------------------------------------------------------
      !--- Check to make sure qdef is always negative when outside of the tropo
      !--- Assume a tropopause height of 10 km; so BCL cannot be higher
      !-----------------------------------------------------------------------------
      where( qmix.ne.missing .and. qsat.ne.missing )  qdef = qsat - qmix
      ! where( hgt.ge.10000 )  qdef = -1.0



      !---------------------------------------------------------------------------------
      !----- Find the point where the sign first turns negative from the ground up
      !---------------------------------------------------------------------------------
      call transitionCheck( qdef, press(1), pblp, press, nlev, missing, i_unsat, i_sat )


      !-----------------------------------------------------------------------------
      !--- If all layers below 150 hPa above the ground are still saturated then call it all saturated
      !--- And use the 1st level stats and call it "convective" because fog is unlikely to be
      !--- deeper than 100 hPa above the ground
      !-----------------------------------------------------------------------------
      if( i_unsat.le.0 .or. i_sat.le.0 ) then
         pot2m  =    missing
         BCLP   =    missing
         TBM    =    missing
         TDEF   =    missing
         return
      end if



      !-----------------------------------------------------------------------------
      !--- Get the upper and lower bounds for each variable to be calc'd at the BCL
      !--- Calculate output variables; BCL height, BCL pressure,
      !--- Buoyant Mixing Potential Temp, and Potential Temperature Deficit
      !-----------------------------------------------------------------------------
      call findIntercept( nlev, i_sat, i_unsat, qdef, logp , missing, BCLP )
      call findIntercept( nlev, i_sat, i_unsat, qdef, tmp_k, missing, TBM  )
      call findIntercept( nlev, i_sat, i_unsat, qdef, qmix , missing, QBCL )
      BCLP  =  exp(BCLP)
      call potentialTemperature1D( TBM, BCLP, missing, TBM)


      !-----------------------------------------------------------------------------
      !--- Virtual Potential Temperature (K) calculation using the mixed humidty,
      !--- *** THIS is an assumption!!  only influences TDEF but an important
      !--- effect because if pot2m is close to TBM then a slight change in qbcl
      !--- can mean the difference between initiation (TDEF=0) or not
      !--- This should only be an issue over very shallow pbls ...
      !-----------------------------------------------------------------------------
      call potentialTemperature1D(tmp_k(1),press(1),missing,pot2m)
      TDEF   =  TBM  - pot2m

      return


end subroutine hcfcalc








!------------------------------------------------------------------------   
!------------------------------------------------------------------------    

subroutine int2p(PPIN,XXIN,NPIN,PPOUT,XXOUT,NPOUT,LINLOG,XMSG)
    implicit none

! routine to interpolate from one set of pressure levels                     
! .   to another set  using linear or ln(p) interpolation                    
!                                                                            
! NCL: xout = int2p (pin,xin,pout,linlog)                                    
!                                                                            
! This code was originally written for a specific purpose.                   
!    Several features were added for incorporation into NCL's                
!    function suite including linear extrapolation.                          
!                                                                            
! nomenclature:                                                              
!                                                                            
!    ppin   - input pressure levels. The pin can be                          
!             be in ascending or descending order                            
!    xxin   - data at corresponding input pressure levels                    
!    npin   - number of input pressure levels >= 2                           
!    ppout  - output pressure levels (input by user)                         
!             same (ascending or descending) order as pin                    
!    xxout  - data at corresponding output pressure levels                   
!    npout  - number of output pressure levels                               
!    linlog - if abs(linlog)=1 use linear interp in pressure                 
!             if abs(linlog) anything but =1 linear interp in                
!                 ln(pressure)                                               
!             If the value is negative then the routine will                 
!             extrapolate. Be wary of results in this case.                  
!    xmsg   - missing data code. if none, set to some number                 
!             which will not be encountered (e.g., 1.e+36)                   
!    ier    - error code                                                     

      ! input types                                                          
      integer, intent(in ) :: NPIN,NPOUT,LINLOG
      real(4), intent(in ) :: PPIN(NPIN),XXIN(NPIN),PPOUT(NPOUT),XMSG
      ! output                                                               
      real(4), intent(out) :: XXOUT(NPOUT)

      ! local                                                                
      integer :: NP,NL,NLMAX,NLSAVE,NP1,NO1,N1,N2,LOGLIN,NLSTRT
      real(4) :: SLOPE,PA,PB,PC

      ! automatic arrays                                                                                                                  
      real(4) :: PIN(NPIN),XIN(NPIN),P(NPIN),X(NPIN)
      real(4) :: POUT(NPOUT),XOUT(NPOUT)

      ! error code
      integer :: IER

!------------------------------------------------------------------------    


      ! error check: enough points: pressures consistency?                   
      LOGLIN = ABS(LINLOG)
      IER = 0
      IF (NPOUT.GT.0) THEN
          DO NP = 1,NPOUT
              XXOUT(NP) = XMSG
          END DO
      END IF

      IF (NPIN.LT.2 .OR. NPOUT.LT.1) IER = IER + 1

      IF (IER.NE.0) THEN
          RETURN
      END IF

      ! should input arrays be reordered: want p(1) > p(2) > p(3) etc        
      ! so that it will match order for which code was originally designed   
      ! copy to local arrays                                                 
      NP1 = 0
      NO1 = 0
      IF (PPIN(1).LT.PPIN(2)) THEN
          NP1 = NPIN + 1
          NO1 = NPOUT + 1
      END IF

      DO NP = 1,NPIN
          PIN(NP) = PPIN(ABS(NP1-NP))
          XIN(NP) = XXIN(ABS(NP1-NP))
      END DO

      DO NP = 1,NPOUT
          POUT(NP) = PPOUT(ABS(NO1-NP))
      END DO
      !                                                                      
      ! eliminate levels with missing data. This can easily                  
      ! .   happen with observational data.                                  
      !                                                                      
      NL = 0
      DO NP = 1,NPIN
          IF (XIN(NP).NE.XMSG .AND. PIN(NP).NE.XMSG) THEN
              NL = NL + 1
              P(NL) = PIN(NP)
              X(NL) = XIN(NP)
          END IF
      END DO
      NLMAX = NL
      ! all missing data                                                                                                                  
      IF (NLMAX.LT.2) THEN
          IER = IER + 1000
          RETURN
      END IF


      ! ===============> pressure in decreasing order <================      
      ! perform the interpolation  [pin(1)>pin(2)>...>pin(npin)]             
      !                                                      ( p ,x)         
      ! ------------------------- p(nl+1), x(nl+1)   example (200,5)         
      ! .                                                                    
      ! ------------------------- pout(np), xout(np)         (250,?)         
      ! .                                                                    
      ! ------------------------- p(nl)  , x(nl)             (300,10)        
      ! exact p-level matches                                                
      NLSTRT = 1
      NLSAVE = 1
      DO NP = 1,NPOUT
          XOUT(NP) = XMSG
          DO NL = NLSTRT,NLMAX
              IF (POUT(NP).EQ.P(NL)) THEN
                  XOUT(NP) = X(NL)
                  NLSAVE = NL + 1
                  GO TO 10
              END IF
          END DO
   10     NLSTRT = NLSAVE
      END DO

      IF (LOGLIN.EQ.1) THEN
          DO NP = 1,NPOUT
              DO NL = 1,NLMAX - 1
                  IF (POUT(NP).LT.P(NL) .AND. POUT(NP).GT.P(NL+1)) THEN
                      SLOPE = (X(NL)-X(NL+1))/ (P(NL)-P(NL+1))
                      XOUT(NP) = X(NL+1) + SLOPE* (POUT(NP)-P(NL+1))
                  END IF
              END DO
          END DO
      ELSE
          DO NP = 1,NPOUT
              DO NL = 1,NLMAX - 1
                  IF (POUT(NP).LT.P(NL) .AND. POUT(NP).GT.P(NL+1)) THEN
                      PA = LOG(P(NL))
                      PB = LOG(POUT(NP))
                      ! special case: In case someome inadvertently enter p=0
                      if (p(nl+1).gt.0.d0) then
                          PC = LOG(P(NL+1))
                      else
                          PC = LOG(1.d-4)
                      end if

                      SLOPE = (X(NL)-X(NL+1))/ (PA-PC)
                      XOUT(NP) = X(NL+1) + SLOPE* (PB-PC)
                  END IF
              END DO
          END DO
      END IF

      ! extrapolate?                                   
      ! . use the 'last' valid slope for extrapolating 
      IF (LINLOG.LT.0) THEN
          DO NP = 1,NPOUT
              DO NL = 1,NLMAX
                  IF (POUT(NP).GT.P(1)) THEN
                      IF (LOGLIN.EQ.1) THEN
                          SLOPE = (X(2)-X(1))/ (P(2)-P(1))
                          XOUT(NP) = X(1) + SLOPE* (POUT(NP)-P(1))
                      ELSE
                          PA = LOG(P(2))
                          PB = LOG(POUT(NP))
                          PC = LOG(P(1))
                          SLOPE = (X(2)-X(1))/ (PA-PC)
                          XOUT(NP) = X(1) + SLOPE* (PB-PC)
                      END IF
                  ELSE IF (POUT(NP).LT.P(NLMAX)) THEN
                      N1 = NLMAX
                      N2 = NLMAX - 1
                      IF (LOGLIN.EQ.1) THEN
                          SLOPE = (X(N1)-X(N2))/ (P(N1)-P(N2))
                          XOUT(NP) = X(N1) + SLOPE* (POUT(NP)-P(N1))
                      ELSE
                          PA = LOG(P(N1))
                          PB = LOG(POUT(NP))
                          PC = LOG(P(N2))
                          SLOPE = (X(N1)-X(N2))/ (PA-PC)
                          XOUT(NP) = X(N1) + SLOPE* (PB-PC)
                      END IF
                  END IF
              END DO
          END DO
      END IF

      ! place results in the return array;    
      ! .   reverse to original order         

      DO NP = 1,NPOUT
          XXOUT(NP) = XOUT(ABS(NO1-NP))
      END DO

end subroutine int2p












!---------------------------------------------------------------------------------
!
! Description:
!
! Calculate the convective available potential energy (CAPE) given launch level
! quantities of temperature, specific humidity, and pressure AND the actual profiles
! of T, q, and P.  Launch level variables must be contained within the T, q, P profiles.
!
! Author and Revision History: 
! A.B. Tawfik on May 2016
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!
! subroutines:  Calculates CAPE given a launch and equilibrium pressure level
!
!---------------------------------------------------------------------------------
 subroutine calc_cape( nlev, tlev_in, qlev_in, plev_in, tlaunch, qlaunch, plaunch, CAPE, missing )

   implicit none

!
! Input/Output Variables
!
   integer, intent(in ) :: nlev           !** # of vertical levels
   real(4), intent(in ) :: missing        !** missing value
   real(4), intent(in ) :: qlev_in(nlev)  !** specific humidity [kg/kg]
   real(4), intent(in ) :: plev_in(nlev)  !** pressure [Pa]
   real(4), intent(in ) :: tlev_in(nlev)  !** temperature [K]
   real(4), intent(in ) :: tlaunch        !** launch level temperature [K]
   real(4), intent(in ) :: plaunch        !** launch level pressure [Pa]
   real(4), intent(in ) :: qlaunch        !** launch level spc. humidity [kg/kg]
   real(4), intent(out) :: CAPE           !** CAPE from Launch-to-Equilibrium [J/kg]
!
! Local variables
!
   real(4), parameter         ::  Rd=287.04, cp=1005.7
   real(4), parameter         ::  Lv=2.5e6 , Rv=461.5, grav= 9.81
   real(4), parameter         ::  epsilon=0.622

   !---- constants for vapor pressure calculation Arden Buck 1996
   real(4), parameter         ::  c  = 6.1121, a0 = 18.678, b0 = 234.5
   real(4), parameter         ::  c0 = 257.14, C2K = 273.15  

   logical, dimension(nlev )  ::  notmissing
   real(4), dimension(nlev )  ::  allmissing
   real(4), dimension(nlev )  ::  vapor_pressure, rlev_in

   real(4), dimension(nlev+1) ::  dz
   real(4), dimension(nlev+1) ::  ilev
   real(4), dimension(nlev+1) ::  numerator
   real(4), dimension(nlev+1) ::  denomenator
   real(4), dimension(nlev+1) ::  moist_lapse
   real(4), dimension(nlev+1) ::  temp_difference
   real(4), dimension(nlev+1) ::  tmid, rmid, qmid, pmid
   real(4), dimension(nlev+1) ::  virtual_temperature   
   real(4), dimension(nlev+1) ::  tlev, rlev, plev, qlev
   real(4), dimension(nlev+1) ::  parcel_tmid, parcel_temp
   real(4), dimension(nlev+1) ::  wsat, esat

   integer                    ::  zz, nlev1
   integer                    ::  above_LNB, below_LNB
   real(4)                    ::  rlaunch, vplaunch
   real(4)                    ::  top_level_temp_diff
   real(4)                    ::  t_not_buoyant, t_buoyant
   real(4)                    ::  p_not_buoyant, p_buoyant
   real(4)                    ::  LNB            !** Level of Neutral Buoyancy (LNB) [Pa]

   

!---------------------------------------------------------------------------------

      !--------------------------------
      !-- Initialize variables
      !--------------------------------
      nlev1                =  nlev + 1
      dz                   =  missing
      rlev_in              =  missing
      rlev                 =  missing
      tlev                 =  missing
      plev                 =  missing
      qlev                 =  missing
      rmid                 =  missing
      tmid                 =  missing
      qmid                 =  missing
      moist_lapse          =  missing
      vapor_pressure       =  missing
      virtual_temperature  =  missing
      temp_difference      =  missing 
      parcel_temp          =  missing
      parcel_tmid          =  missing 
      allmissing           =  missing
      notmissing           =  .false.

      CAPE                 =  missing
      LNB                  =  missing


      !-------------------------------------------------------------------
      !--- Make sure launch level pressure is within the profile bounds
      !-------------------------------------------------------------------
      if( plaunch.lt.minval(plev_in, mask = plev_in.ne.missing) .or. plaunch.gt.maxval(plev_in, mask = plev_in.ne.missing) ) then
         write(*,*) " ERROR:  Launch level is not within the bounds "
         write(*,*) "         Launch   pressure --->  ",plaunch,qlaunch,tlaunch
         do zz=1,nlev
            write(*,*) " Pressure Levels ---> ",zz,plev_in(zz),qlev_in(zz)*1e3,tlev_in(zz)
         end do
         stop
      end if


 

      !--------------------------------
      !-- define the level index
      !-- this is used for locating
      !-- indices below and above a
      !-- desired pressure level
      !--------------------------------
      do zz=1,nlev1
         ilev(zz)  =  real(zz)
      end do



      !-------------------------------------------------------------------
      !--- convert specific humidity to vapor mixing ratio (rlev)
      !--- For both column variables and launch variables
      !-------------------------------------------------------------------
      where( qlev_in       .ne.missing .and. plev_in.ne.missing )   &
      vapor_pressure  =  (qlev_in*plev_in) / (epsilon + qlev_in - (qlev_in*epsilon))

      where( vapor_pressure.ne.missing .and. plev_in.ne.missing )   &
      rlev_in         =  epsilon * ( vapor_pressure / (plev_in - vapor_pressure) )

      vplaunch =  (qlaunch*plaunch) / (epsilon + qlaunch - (qlaunch*epsilon))
      rlaunch  =  epsilon * ( vplaunch / (plaunch - vplaunch) )


      !------------------------------------------------------------------------------------------
      !-- Remove all levels below the launch level to make calculations easier and code cleaner
      !-- Then place launch-level variables in the first index signifying the lowest level
      !-- where the CAPE calculation begins
      !------------------------------------------------------------------------------------------
      call packIt(tlev_in, tlaunch , nlev   , nlev1, &
                  plaunch, plev_in , missing, tlev  )

      call packIt(rlev_in, rlaunch , nlev   , nlev1, &
                  plaunch, plev_in , missing, rlev  )

      call packIt(qlev_in, qlaunch , nlev   , nlev1, &
                  plaunch, plev_in , missing, qlev  )

      call packIt(plev_in, plaunch , nlev   , nlev1, &
                  plaunch, plev_in , missing, plev  )


!      tlev(2:)  =  tlev_in
!      plev(2:)  =  plev_in
!      qlev(2:)  =  qlev_in
!      where( plev(2:).gt.plaunch )   plev(2:)  =  missing
!      notmissing  =  .not.( plev(2:).eq.missing )
!      tlev(2:)    =  pack ( tlev_in, notmissing, allmissing )
!      rlev(2:)    =  pack ( rlev   , notmissing, allmissing )
!      plev(2:)    =  pack ( plev_in, notmissing, allmissing )
!      qlev(2:)    =  pack ( qlev_in, notmissing, allmissing )
!      tlev(1)     =  tlaunch
!      rlev(1)     =  rlaunch
!      plev(1)     =  plaunch
!      qlev(1)     =  qlaunch


      !-------------------------------------------------------------------
      !--- Check level is present or not
      !-------------------------------------------------------------------
      if( all(tlev.eq.missing)  .or.  all(plev.eq.missing)  .or. all(qlev.eq.missing) ) then
         CAPE = 0.0
         LNB  = 0.0
         return
      end if

      !-------------------------------------------------------------------
      !--- Calculate layer averages (mid-point)
      !-------------------------------------------------------------------
      where( tlev(2:nlev1).ne.missing  .and.  tlev(:nlev).ne.missing  .and.  &
             plev(2:nlev1).ne.missing  .and.  plev(:nlev).ne.missing         )
         tmid(:nlev)  =  ( (tlev(2:nlev1)*log(plev(2:nlev1))  +  tlev(:nlev)*log(plev(:nlev)))  /  &
                            log(plev(2:nlev1) * plev(:nlev)) )
      end where

      where( rlev(2:nlev1).ne.missing  .and.  rlev(:nlev).ne.missing  .and.  &
             plev(2:nlev1).ne.missing  .and.  plev(:nlev).ne.missing         )
         rmid(:nlev)  =  ( (rlev(2:nlev1)*log(plev(2:nlev1))  +  rlev(:nlev)*log(plev(:nlev)))  /  &
                            log(plev(2:nlev1) * plev(:nlev)) )
      end where

      where( qlev(2:nlev1).ne.missing  .and.  qlev(:nlev).ne.missing  .and.  &
             plev(2:nlev1).ne.missing  .and.  plev(:nlev).ne.missing         )
         qmid(:nlev)  =  ( (qlev(2:nlev1)*log(plev(2:nlev1))  +  qlev(:nlev)*log(plev(:nlev)))  /  &
                            log(plev(2:nlev1) * plev(:nlev)) )
      end where

      where( plev(2:nlev1).ne.missing  .and.  plev(:nlev).ne.missing         )
         pmid(:nlev)  =  ( (plev(2:nlev1)*log(plev(2:nlev1))  +  plev(:nlev)*log(plev(:nlev)))  /  &
                            log(plev(2:nlev1) * plev(:nlev)) )
      end where


      !-------------------------------------------------------------------
      !--- Calculate layer average virtual temperature [K]
      !-------------------------------------------------------------------
      where( tmid.ne.missing  .and.  rmid.ne.missing )
          virtual_temperature  =  tmid * (1.0 + rmid/epsilon) * (1.0 + rmid)
      end where


      !------------------------------------------------------------------------
      !--- calculate saturation mixing ratio for the moist lapse rate formula
      !--- Following Goff-Gratch (1946) which is accurate range from 
      !--- -50 to 100 C given in the Smithsonian Meteorological Tables
      !------------------------------------------------------------------------
      where( tmid.ne.missing  .and.  pmid.ne.missing )
         esat  =  c * exp( (a0 - (tmid-C2K)/b0) * ((tmid-C2K)/(c0+tmid-C2K)) ) * 1e2
         wsat  =  epsilon * ( esat / (pmid-esat) )
      end where



      !----------------------------------------------------------------
      !--- Get moist adiabatic lapse rate [K/m] and
      !--- the depth of the layer from lower to upper level
      !--- 1)  Get thickness of each layer using hypsometric equation
      !---     Wallace and Hobbs (1977) Atmospheric Science: An Introductory Survey. Academic Press, 55-57.
      !--- 2)  Then calculate moist adiabatic lapse rate using AMS standard
      !----------------------------------------------------------------
      where( virtual_temperature(:nlev).ne.missing   .and.                            &
             plev               (:nlev).ne.missing   .and.  plev(2:nlev1).ne.missing  )
          dz(:nlev)  =  (Rd * virtual_temperature(:nlev) / grav) * log(plev(:nlev)/plev(2:nlev1))
      end where

      where( tmid.ne.missing  .and.  wsat.ne.missing ) 
          numerator    =  1.0  +  (Lv    * wsat          )/(     Rd * tmid   )
          denomenator  =  1.0  +  (Lv**2 * wsat * epsilon)/(cp * Rd * tmid**2)
          moist_lapse  =  grav/cp  *  (numerator/denomenator)
      end where



      !----------------------------------------------------
      !--- Return the temperature change at each level
      !--- due to lifting a parcel moist adiabatically
      !--- note that this is just the temperature change 
      !--- of a parcel and NOT the parcel temperature yet
      !----------------------------------------------------
      !----------------------------------------------------
      !--- Get the parcel temperature [K]
      !----------------------------------------------------
      parcel_temp(1)  =  tlaunch
      do zz=2,nlev1
         if( dz(zz-1).ne.missing .and. moist_lapse(zz-1).ne.missing ) then
            parcel_temp(zz)  =  parcel_temp(zz-1) - moist_lapse(zz-1)*dz(zz-1)
         end if
      end do
 
      where( parcel_temp(2:nlev1).ne.missing  .and.  parcel_temp(:nlev).ne.missing  .and.  &
             plev       (2:nlev1).ne.missing  .and.  plev       (:nlev).ne.missing         )
         parcel_tmid(:nlev)  =  ( (parcel_temp(2:nlev1)*log(plev(2:nlev1))  +  &
                                   parcel_temp( :nlev )*log(plev( :nlev)))  /  &
                                   log(plev(2:nlev1) * plev(:nlev)) )
      end where


      !---------------------------------------------------------------------
      !--- Get the parcel temperature minus environmental temperature [K]
      !---------------------------------------------------------------------
      where( parcel_tmid.ne.missing  .and.  tmid.ne.missing )   temp_difference  =  parcel_tmid  -  tmid



      !------------------------------------------------------------------------
      !--- Find the level of neutral buoyancy (LNB) in [Pa]
      !--- above_LNB = level index corresponding to first negatively buoyant layer
      !--- below_LNB = level index corresponding to last positively buoyant layer
      !------------------------------------------------------------------------
      above_LNB  =  minloc( ilev, dim = 1, mask = temp_difference.ne.missing  .and.  temp_difference.lt.0 )
      below_LNB  =  maxloc( ilev, dim = 1, mask = temp_difference.ne.missing  .and.  temp_difference.ge.0 )



      !-----------------------------------------------------------------------------
      !--- Linearly interpolate between the levels directly below and above where the
      !--- LNB is located (the LNB is located somewhere inbetween these levels so
      !--- so find the zero point)
      !--- Note that the LNB serves as an upper limit of the potential cloud top
      !-----------------------------------------------------------------------------
      t_not_buoyant  =  temp_difference(above_LNB)
      t_buoyant      =  temp_difference(below_LNB)

      p_not_buoyant  =  pmid(above_LNB)
      p_buoyant      =  pmid(below_LNB)

      LNB            =  ( p_not_buoyant - ((p_not_buoyant-p_buoyant) / (t_not_buoyant-t_buoyant))*t_not_buoyant )


      !------------------------------------------------------------------------
      !--- Now calculate CAPE by integrating until parcel temperature is less
      !--- than the environmental temperature (meaning parcel is no longer
      !--- buoyant)
      !------------------------------------------------------------------------
      CAPE  =  sum( Rd * temp_difference(:nlev) * log(plev(:nlev)/plev(2:nlev1)),                                &
                    mask = plev           (:nlev).ne.missing  .and.  plev           (2:nlev1).ne.missing  .and.  &
                           temp_difference(:nlev).ne.missing  .and.  temp_difference(:nlev  ).ge.0               )  !*** temp diff assures buoyancy

      top_level_temp_diff  =  temp_difference(below_LNB)*log(plev(below_LNB))  /  log(plev(below_LNB)*LNB)
      CAPE                 =  CAPE  +  (Rd * top_level_temp_diff * log(plev(below_LNB)/LNB))



      return


 end subroutine calc_cape
!---------------------------------------------------------------------------------







!-----------------------------------------------------------------------------
!
! Primary Subroutine: At a high-level this subroutine takes real atmospheric
! profiles and sees whether convection initiates under various evaporative 
! fraction regimes.  This subroutine begins with an early morning sounding
! and then evolves forward in time with a given net radiation that is split by 
! the desired evaporative fraction.  Overall the subroutine can explore the entire
! possible space of evaporative fraction for a given day to see whether CI occurred
! - Is convection more likely over dry or 'wet' surface flux properties?
! - Applying it spatially can give an answer to the question of positive or
!   negative flux-CI feedbacks over certain regions and times of year
!
!-----------------------------------------------------------------------------
subroutine Evaluate_CI( T   , Q   , P        , itime                  , &
                        t2m , q2m , psfc     , rnet  , ef    , dt     , &
                        nlev, nlat, nlon     , nhr , nday  , num_EF, missing, &
                        TBM , BCLP, TimeOfCI , CAPE                           )

   implicit none
!
! Input/Output Variables
!
   integer, intent(in  )                                       ::  nday        ! *** # days
   integer, intent(in  )                                       ::  nhr         ! *** # number of hours per day
   integer, intent(in  )                                       ::  nlev        ! *** # of atmospheric levels
   integer, intent(in  )                                       ::  nlat        ! *** # latitude
   integer, intent(in  )                                       ::  nlon        ! *** # longitude
   integer, intent(in  )                                       ::  num_EF      ! *** # number of evaporative fraction breakdowns 
   integer, intent(in  )                                       ::  itime       ! *** Start time of morning sounding

   real(4), intent(in  )                                       ::  missing     ! *** Missing values
   real(4), intent(in  )                                       ::  dt          ! *** timestep in seconds [s/timestep]
   real(4), intent(in  ), dimension(nday     ,nlev,nlat,nlon)  ::  T  , Q      ! *** Temp and Humidity over level in SI
   real(4), intent(in  ), dimension(nday          ,nlat,nlon)  ::  t2m, q2m    ! *** 2-m temp and humidity, Height in SI

   real(4), intent(in  ), dimension(nday     ,nlev,nlat,nlon)  ::  P           ! *** Pressure and Height in (level) SI
   real(4), intent(in  ), dimension(nday          ,nlat,nlon)  ::  psfc        ! *** 2-m pressure and height in SI
   real(4), intent(in  ), dimension(nday,nhr      ,nlat,nlon)  ::  rnet        ! *** Net Radiation time series [W/m2]
   real(4), intent(in  ), dimension(     num_EF             )  ::  ef          ! *** Evaporative Fraction levels

   real(4), intent(out ), dimension(nday,num_EF   ,nlat,nlon)  ::  TBM         ! *** Buoyant mixing theta at time of CI [K]
   real(4), intent(out ), dimension(nday,num_EF   ,nlat,nlon)  ::  BCLP        ! *** Cloud base pressure at time of CI [K]
   real(4), intent(out ), dimension(nday,num_EF   ,nlat,nlon)  ::  TimeOfCI    ! *** Time of day of CI [0-24 hour]
   real(4), intent(out ), dimension(nday,num_EF   ,nlat,nlon)  ::  CAPE        ! *** CAPE when BCL is reached [J/kg]


!
! Local variables
!
!   integer, parameter         ::  itime = 3  
   real(4), parameter         ::  omega = 0.0
   real(4), dimension(nlev+1) ::  ppack, tpack, hpack, qpack, Theta, newTheta, newQhum, dpress, density
   real(4)                    ::  pbl_depth, pblp, pbl_theta, pblh, latent_heat, sensible_heat
   real(4)                    ::  evap, Lc, tbm_out, tdef_out, bclp_out, pressure_deficit
   real(4)                    ::  avgRHO, qbcl_out
   integer                    ::  xx, yy, dd, tt, ee
   integer                    ::  nlev1



!-----------------------------------------------------------------------------

      !-----------------------------------------------------------------------------
      !-- Initialize output variables
      !-----------------------------------------------------------------------------
      TBM       =  missing
      BCLP      =  missing
      CAPE      =  missing
      TimeOfCI  =  missing
      nlev1     =  nlev + 1


      !-----------------------------------------------------------------------------
      !-- Loop over time, lat, and lon
      !-----------------------------------------------------------------------------
      lati_loop: do yy = 1,nlat 
      long_loop: do xx = 1,nlon
      ef_loop:   do ee = 1,num_EF
      day_loop:  do dd = 1,nday


         !--------------------------------------
         !-- Append 2-m quantiy to the profile 
         !--------------------------------------
         call packIt(T   (dd,:,yy,xx), t2m(dd  ,yy,xx), nlev   , nlev1, &
                     psfc(dd  ,yy,xx), P  (dd,:,yy,xx), missing, tpack  )

         call packIt(Q   (dd,:,yy,xx), q2m(dd  ,yy,xx), nlev   , nlev1, &
                     psfc(dd  ,yy,xx), P  (dd,:,yy,xx), missing, qpack  )

         call packIt(P   (dd,:,yy,xx), psfc(dd  ,yy,xx), nlev   , nlev1, &
                     psfc(dd  ,yy,xx), P   (dd,:,yy,xx), missing, ppack  )


         !-----------------------------------------------------------------------------
         !-- Calculate the potential temperature [This is mutable]
         !-----------------------------------------------------------------------------
         call potentialTemperature(tpack, ppack, nlev1, missing, Theta)


         !-----------------------------------------------------------------------------
         !-- Check if CI occurs yet??  No don't need to because I want to erode the 
         !-- fog layer if one exists
         !-----------------------------------------------------------------------------



         !-----------------------------------------------------------------------------
         !-- Loop over time of day until Initiation occurs
         !-- Running a little mini-simulation!!!!!  
         !-----------------------------------------------------------------------------
         time_of_day: do tt = itime,nhr

            !-----------------------------------------------------------------------------
            !-- Don't do anything if there is no heat or moisture being added
            !-- Otherwise partition the evaporative fraction between LH and SH
            !-----------------------------------------------------------------------------
            if( rnet(dd,tt,yy,xx).le.0 ) cycle time_of_day
            latent_heat    =  ef(ee)         * rnet(dd,tt,yy,xx)
            sensible_heat  =  (1.0 - ef(ee)) * rnet(dd,tt,yy,xx)


            !-----------------------------------------------------------------------------
            !-- Calculate the height above ground [m]
            !-----------------------------------------------------------------------------
            call calculate_height_above_ground (tpack, ppack, nlev1, missing, hpack)


            !-----------------------------------------------------------------------------
            !-- Calcuate the boundary layer top variables but with no extra water vapor
            !-- added to the column
            !-----------------------------------------------------------------------------
            call pblHeat ( nlev1, missing, tpack, ppack, hpack, pblp, pbl_theta, pblh )


            !-----------------------------------------------------------------------------
            !-- Assign PBL quantities across the depth of the PBL
            !-----------------------------------------------------------------------------
            call assign_layer (Theta, pbl_theta, pblh, hpack, nlev1, missing, newTheta )
            Theta      =  newTheta
            pbl_depth  =  ppack(1) - pblp


            !-----------------------------------------------------------------------------------
            !-- Calculate density for column -> used for sensible heat calculation
            !-----------------------------------------------------------------------------------
            call total_density    ( ppack  , tpack, qpack, nlev1, missing, density )
            call avg_over_layer   ( density, pblh , hpack, nlev1, missing, avgRHO  )

            !-----------------------------------------------------------------------------------
            !-- Adds some sensible heat flux to the boundary layer to increase the temperature
            !-----------------------------------------------------------------------------------
            call add_sensible_heat( newTheta, sensible_heat, avgRHO, hpack, pblh, dt, nlev1, missing, Theta ) 


            !-----------------------------------------------------------------------------------
            !-- Return corresponding temperature given the previous heat added
            !-----------------------------------------------------------------------------------
            call calculate_temperature(Theta, ppack, nlev1, missing, tpack)





            !-----------------------------------------------------------------------------------
            !-- Adds some latent heat flux to the boundary layer to increase the speciific humidity
            !-----------------------------------------------------------------------------------
            !-----------------------------------------------------------------------------
            !----- Get the depth of each layer
            !-----------------------------------------------------------------------------
            call layerDepth( ppack, tpack(1), qpack(1), ppack(1), hpack(1), nlev1, missing, dpress )

            !-----------------------------------------------------------------------------
            !----- Convert from latent heat flux [W/m2] to evapotranspiration [kg/m2/timestep]
            !-----------------------------------------------------------------------------
            call Latent_heat_of_condensation(tpack(1), missing, Lc)
            evap  =  latent_heat * dt / Lc 

            !-----------------------------------------------------------------------------
            !----- inject into the boundary layer -- NOTE returns humidity as kg/m2
            !-----------------------------------------------------------------------------
            call inject_and_mix( pblp, evap, ppack, qpack, dpress, nlev1, missing, newQhum )
            qpack  =  newQhum

            !-----------------------------------------------------------------------------
            !----- Calculate HCF variables in new system and check if CI occurred
            !-----------------------------------------------------------------------------
            !call hcfcalc ( nlev1, missing, tpack, ppack, qpack, hpack, pblp, tbm_out, tdef_out, bclp_out )
            call hcfcalc ( nlev1, missing, tpack, ppack, qpack, hpack, pblp, tbm_out, tdef_out, bclp_out, qbcl_out )

            !-----------------------------------------------------------------------------
            !----- Check for CI
            !-----------------------------------------------------------------------------
            pressure_deficit  =  (pblp + omega*dt) - bclp_out 
            if( pressure_deficit.le.0 ) then

               !! write(*,*) "  We have achieved CI  ", dd,tt, tbm_out, bclp_out/1e2, pressure_deficit/1e2
               TBM     (dd,ee,yy,xx)  =  tbm_out
               BCLP    (dd,ee,yy,xx)  =  bclp_out
               TimeOfCI(dd,ee,yy,xx)  =  tt * 1.0

               !!!!add CAPE calculation later to figure out some sort of intensity
               !call calculate_temperature(tbm_out, ppack, nlev1, missing, tbcl_out)
               !call calc_cape( nlev1, tpack, qpack, ppack, tbcl_out, qbcl_out, bclp_out, CAPE(dd,ee,yy,xx), missing )

               exit time_of_day
            end if
!!! Remember tpack and qpack need to be re-assigned  qpack = newQhum
!!! Add an omega feature to calculating CI
!!! make sure to define all these variables
!!! define what ef is -- maybe pass it in?
!!! exit this time_of_day loop if CI occurs
!!! clean up code?

         end do time_of_day


      end do day_loop
      end do ef_loop
      end do long_loop
      end do lati_loop



end subroutine Evaluate_CI




!---------------------------------------------------------------------------------
!
! subroutines:  calculates buoyant condensation level and basic variables (THETA_BM; TDEF)
!
!---------------------------------------------------------------------------------
subroutine hcfloop ( nlev, nlat, nlon, nday, missing, T, P, Q, t2m, q2m, psfc, TBM, TDEF, BCLP, QBCL )

   implicit none
!
! Input/Output Variables
!
   integer, intent(in  )                                  ::  nlev        ! *** # of atmospheric levels
   integer, intent(in  )                                  ::  nlat        ! *** # of horizontal grid size
   integer, intent(in  )                                  ::  nlon        ! *** # of horizontal grid size
   integer, intent(in  )                                  ::  nday        ! *** # of horizontal grid size
   real(4), intent(in  )                                  ::  missing     ! *** Missing values
   real(4), intent(in  ), dimension(     nlev          )  ::  P           ! *** Pressure (level) [Pa]
   real(4), intent(in  ), dimension(nday,nlev,nlat,nlon)  ::  T           ! *** Temperature (level), [K]
   real(4), intent(in  ), dimension(nday,nlev,nlat,nlon)  ::  Q           ! *** Specific Humidity (level) [kg/kg]

   real(4), intent(in  ), dimension(nday     ,nlat,nlon)  ::  t2m         ! *** 2-m temperature [K]
   real(4), intent(in  ), dimension(nday     ,nlat,nlon)  ::  q2m         ! *** 2-m specific humidity [kg/kg]
   real(4), intent(in  ), dimension(nday     ,nlat,nlon)  ::  psfc        ! *** surface pressure [Pa]

   real(4), intent(out ), dimension(nday     ,nlat,nlon)  ::  TBM         ! *** buoyant mixing pot. temp (convective threshold) [K]
   real(4), intent(out ), dimension(nday     ,nlat,nlon)  ::  TDEF        ! *** pot. temp deficit need to initiate [K]
   real(4), intent(out ), dimension(nday     ,nlat,nlon)  ::  BCLP        ! *** pressure of buoyant condensation level [Pa]
   real(4), intent(out ), dimension(nday     ,nlat,nlon)  ::  QBCL        ! *** Specific Humidity at bcl [kg/kg]


!
! Local variables
!
   integer                    ::  xx, yy, dd, nlev1
   real(4), dimension(nlev+1) ::  ppack, tpack, hpack, qpack, Theta
   real(4)                    ::  pblp, pbl_theta, pblh

!-----------------------------------------------------------------------------

      !-----------------------------------------------------------------------------
      !-- Initialize output variables
      !-----------------------------------------------------------------------------
      nlev1 =  nlev + 1
      TBM   =  missing
      TDEF  =  missing
      BCLP  =  missing
      QBCL  =  missing
      

      !-----------------------------------------------------------------------------
      !-- Loop over time, lat, and lon
      !-----------------------------------------------------------------------------
      lati_loop: do yy = 1,nlat 
      long_loop: do xx = 1,nlon
      days_loop: do dd = 1,nday


         !--------------------------------------
         !-- Append 2-m quantiy to the profile 
         !--------------------------------------
         call packIt(T    (dd,:,yy,xx), t2m(dd  ,yy,xx) , nlev   , nlev1, &
                     psfc (dd  ,yy,xx), P               , missing, tpack  )

         call packIt(Q    (dd,:,yy,xx), q2m(dd  ,yy,xx) , nlev   , nlev1, &
                     psfc (dd  ,yy,xx), P               , missing, qpack  )

         call packIt(P                , psfc(dd  ,yy,xx), nlev   , nlev1, &
                     psfc (dd  ,yy,xx), P               , missing, ppack  )


         !-----------------------------------------------------------------------------
         !-- Calculate the height above ground [m]
         !-----------------------------------------------------------------------------
         call calculate_height_above_ground (tpack, ppack, nlev1, missing, hpack)


         !-----------------------------------------------------------------------------
         !-- Calcuate the boundary layer top variables but with no extra water vapor
         !-- added to the column
         !-----------------------------------------------------------------------------
         call pblHeat ( nlev1, missing, Theta, ppack, hpack, pblp, pbl_theta, pblh )


         !-----------------------------------------------------------------------------
         !-- Calcuate the HCF variables but with no extra water vapor
         !-- added to the column
         !-----------------------------------------------------------------------------
!        call hcfcalc ( nlev , missing, tmp(:,xx), press(:,xx), qhum(:,xx), hgt(:,xx), PBLP(xx), TBM(xx), TDEF(xx), BCLP(xx) )
         call hcfcalc ( nlev1, missing, tpack, ppack, qpack, hpack, pblp, &
                        TBM(dd,yy,xx), TDEF(dd,yy,xx), BCLP(dd,yy,xx), QBCL(dd,yy,xx) )


      end do days_loop
      end do long_loop
      end do lati_loop



end subroutine hcfloop







!---------------------------------------------------------------------------------
!
! subroutines: Estimates the PBL height by using the gradient in potential temperature
!              1st smooth the profiles if the resolution is high
!              2nd calculate the vertical gradient of d(theta)/d(pressure)
!              3rd find the level where the gradient is <= -0.0002 K/Pa
!
!---------------------------------------------------------------------------------
subroutine pbl_gradient ( nlev, missing, theta, pressure, height, PBLP, PBLT, PBLH )

   implicit none
!
! Input/Output Variables
!
   integer, intent(in   )                   ::  nlev        ! *** # of atmospheric levels
   real(4), intent(in   )                   ::  missing     ! *** Missing values
   real(4), intent(in   ), dimension(nlev)  ::  theta       ! *** Potential Temperature (level), [K]
   real(4), intent(in   ), dimension(nlev)  ::  pressure    ! *** Pressure (level) [Pa]
   real(4), intent(in   ), dimension(nlev)  ::  height      ! *** Height (level) [m]
   real(4), intent(out  )                   ::  PBLP        ! *** pressure of boundary layer height [Pa]
   real(4), intent(out  )                   ::  PBLT        ! *** temperature of boundary layer height [K]
   real(4), intent(out  )                   ::  PBLH        ! *** height of boundary layer height [m]
!
! Local variables
!
   real(4), parameter  ::  threshold_gradient = -0.0003

   integer             ::  iupper, ilower
   real(4)             ::  dTdP(nlev)


!-----------------------------------------------------------------------------

      !-----------------------------------------------------------------------------
      !-- Initialize output variables
      !-----------------------------------------------------------------------------
      PBLP  =  missing
      PBLT  =  missing
      PBLH  =  missing


      !-----------------------------------------------------------------------------
      !-- calculate a smoothed lapse rate only if profile is hi-res 
      !-- e.g. if depths are less than 10 hPa  otherwise just use the same profile
      !-----------------------------------------------------------------------------
      call smoothed_lapse_rate (theta, pressure, nlev, missing, dTdP)


      !-----------------------------------------------------------------------------
      !-- If all levels do not exceed the threshold 
      !-- or threshold is only exceeded above an unrealistic 
      !-----------------------------------------------------------------------------
      ! if( all(dTdP.gt.threshold_gradient) ) then
      !    
      ! end if

      !-----------------------------------------------------------------------------
      !-- Return indices where the threshold exceeds
      !-----------------------------------------------------------------------------
      call minIndex (nlev , 1, (dTdP.ne.missing .and. dTdP.le.threshold_gradient), ilower )
      iupper  =  ilower + 1


      !-----------------------------------------------------------------------------
      !-- Calculate boundary layer top properties -->  Just get linear average
      !-----------------------------------------------------------------------------
      PBLH  =  sum(height  (ilower:iupper)) / 2.0
      PBLT  =  sum(theta   (ilower:iupper)) / 2.0
      PBLP  =  sum(pressure(ilower:iupper)) / 2.0
      return

end subroutine pbl_gradient
