      SUBROUTINE checkdefs
!
!svn $Id: checkdefs.F 873 2017-10-05 20:27:10Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2018 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine checks activated C-preprocessing options for        !
!  consistency.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_strings
!
      USE strings_mod, ONLY : uppercase
!
      implicit none
!
!  Local variable declarations.
!
      integer :: iatms = 0
      integer :: ibbl = 0
      integer :: ibiology = 0
      integer :: idriver = 0
      integer :: itrcHadv = 0
      integer :: itrcVadv = 0
      integer :: itrcHadvtl = 0
      integer :: itrcVadvtl = 0
      integer :: ivelHadv = 0
      integer :: ivelVadv = 0
      integer :: ivmix = 0
      integer :: nearshore = 0
      integer :: sbed_load = 0
      integer :: sbed_type = 0
      integer :: is, lstr, ng
!
!-----------------------------------------------------------------------
!  Report activated C-preprocessing options.
!-----------------------------------------------------------------------
!
      Coptions=' '
      IF (Master) WRITE (stdout,10)
  10  FORMAT (/,' Activated C-preprocessing Options:',/)
  20  FORMAT (1x,a,t22,a)
!
      IF (Master) THEN
        WRITE (stdout,20) TRIM(ADJUSTL(MyAppCPP)), TRIM(ADJUSTL(title))
      END IF
      is=LEN_TRIM(Coptions)+1
      lstr=LEN_TRIM(MyAppCPP)
      Coptions(is:is+lstr)=TRIM(ADJUSTL(MyAppCPP))
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is)=','
!
      IF (Master) WRITE (stdout,20) 'ANA_BPFLUX',                       &
     &   'Analytical bottom passive tracers fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BPFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_BSFLUX',                       &
     &   'Analytical kinematic bottom salinity flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BSFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_BTFLUX',                       &
     &   'Analytical kinematic bottom temperature flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BTFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_GRID',                         &
     &   'Analytical grid set-up.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_GRID,'
!
      IF (ANY(Lnudging(:))) THEN
        IF (Master) WRITE (stdout,20) 'ANA_NUDGCOEF',                   &
     &     'Analytical spatially varying nudging time-scales.'
        is=LEN_TRIM(Coptions)+1
        Coptions(is:is+14)=' ANA_NUDGCOEF,'
      END IF
!
      IF (Master) WRITE (stdout,20) 'ANA_SMFLUX',                       &
     &   'Analytical kinematic surface momentum flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SMFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_SPFLUX',                       &
     &   'Analytical surface passive tracer fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SPFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_SSFLUX',                       &
     &   'Analytical kinematic surface salinity flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SSFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_STFLUX',                       &
     &   'Analytical kinematic surface temperature flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_STFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_WWAVE',                        &
     &   'Analytical wind induced waves.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_WWAVE,'
!
      IF (Master) WRITE (stdout,20) 'ASSUMED_SHAPE',                    &
     &   'Using assumed-shape arrays.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' ASSUMED_SHAPE,'
!
      IF (Master) WRITE (stdout,20) 'DJ_GRADPS',                        &
     &   'Parabolic Splines density Jacobian (Shchepetkin, 2002).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' DJ_GRADPS,'
!
      IF (Master) WRITE (stdout,20) 'DOUBLE_PRECISION',                 &
     &   'Double precision arithmetic.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' DOUBLE_PRECISION,'
! DDMITRY
! DD END
!
      IF (Master) WRITE (stdout,20) 'GLS_MIXING',                       &
     &   'Generic Length-Scale turbulence closure.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' GLS_MIXING,'
      ivmix=ivmix+1
!
      IF (Master) WRITE (stdout,20) 'KANTHA_CLAYSON',                   &
     &   'Kantha and Clayson stability function formulation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' KANTHA_CLAYSON,'
!
      IF (Master) WRITE (stdout,20) 'NONLINEAR',                        &
     &   'Nonlinear Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' NONLINEAR,'
!
      IF (Master) WRITE (stdout,20) '!NONLIN_EOS',                      &
     &   'Linear Equation of State for seawater.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' !NONLIN_EOS,'
!
      IF (Master) WRITE (stdout,20) 'N2S2_HORAVG',                      &
     &   'Horizontal smoothing of buoyancy and shear.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' N2S2_HORAVG,'
!
      IF (Master) WRITE (stdout,20) 'OUT_DOUBLE',                       &
     &   'Double precision output fields in NetCDF files.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' OUT_DOUBLE,'
!
      IF (Master) WRITE (stdout,20) 'POWER_LAW',                        &
     &   'Power-law shape time-averaging barotropic filter.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' POWER_LAW,'
!
      IF (Master) WRITE (stdout,20) 'PROFILE',                          &
     &   'Time profiling activated .'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' PROFILE,'
!
      IF (Master) WRITE (stdout,20) 'K_C4ADVECTION',                    &
     &   'Fourth-order centered differences advection of TKE fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' K_C4ADVECTION,'
!
      IF (Master) WRITE (stdout,20) 'RI_SPLINES',                       &
     &   'Parabolic Spline Reconstruction for Richardson Number.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' RI_SPLINES,'
!
      IF (Master) WRITE (stdout,20) '!RST_SINGLE',                      &
     &   'Double precision fields in restart NetCDF file.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' !RST_SINGLE,'
!
      IF (Master) WRITE (stdout,20) 'SALINITY',                         &
     &   'Using salinity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SALINITY,'
!
      IF (Master) WRITE (stdout,20) 'SEDIMENT',                         &
     &   'Cohesive and noncohesive sediments.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SEDIMENT,'
!
      IF (Master) WRITE (stdout,20) 'SED_DENS',                         &
     &   'Include density of suspended sediment in equation of state.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SED_DENS,'
!
      IF (Master) WRITE (stdout,20) 'SUSPLOAD',                         &
     &   'Activate suspended sediment transport.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SUSPLOAD,'
      IF (Master) WRITE (stdout,20) 'NONCOHESIVE_BED2',                 &
     &   'Activate non-cohesive (modified) sediment bed.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NONCOHESIVE_BED2,'
      sbed_type=sbed_type+1
      IF (Master) WRITE (stdout,20) 'SED_OPA',                          &
     &   'Activate sediment-oil aggregation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SED_OPA,'
      IF (Master) WRITE (stdout,20) 'FLOC_TURB_DISS',                   &
     &   'Activate dissipation for flocculation based on turbulence.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' FLOC_TURB_DISS,'
      IF (Master) WRITE (stdout,20) 'SED_TAU_CD_CONST',                 &
     &   'Activate constant critical stress for deposition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SED_TAU_CD_CONST,'
!
      IF (Master) WRITE (stdout,20) 'SOLVE3D',                          &
     &   'Solving 3D Primitive Equations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SOLVE3D,'
!
      IF (Master) WRITE (stdout,20) 'SPLINES_VDIFF',                    &
     &   'Parabolic Spline Reconstruction for Vertical Diffusion.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' SPLINES_VDIFF,'
!
      IF (Master) WRITE (stdout,20) 'SPLINES_VVISC',                    &
     &   'Parabolic Spline Reconstruction for Vertical Viscosity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' SPLINES_VVISC,'
!
       IF (Master) WRITE (stdout,20) 'SSW_BBL',                         &
     &   'Styles and Glenn Bottom Boundary Layer - modified.'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+8)=' SSW_BBL,'
       ibbl=ibbl+1
!
       IF (Master) WRITE (stdout,20) 'SSW_CALC_UB',                     &
     &   'Internal computation of bottom orbital velocity.'
       is=LEN_TRIM(Coptions)+1
       Coptions(is:is+14)=' SSW_CALC_UB,'
      IF (Master) WRITE (stdout,20) 'TS_HSIMT',                         &

     &   'the 3rd HSIMT-TVD scheme for tracers (Wu and Zhu, 2010).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' TS_HSIMT,'
      itrcHadv=itrcHadv+1
      itrcVadv=itrcVadv+1
!
      IF (Master) WRITE (stdout,20) 'UV_ADV',                           &
     &   'Advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_ADV,'
!
      IF (Master) WRITE (stdout,20) 'UV_COR',                           &
     &   'Coriolis term.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_COR,'
!
      IF (Master) WRITE (stdout,20) 'UV_U3HADVECTION',                  &
     &   'Third-order upstream horizontal advection of 3D momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' UV_U3HADVECTION,'
      ivelHadv=ivelHadv+1
!
      IF (Master) WRITE (stdout,20) 'UV_SADVECTION',                    &
     &   'Parabolic splines vertical advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' UV_SADVECTION,'
      ivelVadv=ivelVadv+1
!
      IF (Master) WRITE (stdout,20) 'VAR_RHO_2D',                       &
     &   'Variable density barotropic mode.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' VAR_RHO_2D,'
!
!-----------------------------------------------------------------------
!  Stop if unsupported C-preprocessing options or report issues with
!  particular options.
!-----------------------------------------------------------------------
!
      CALL checkadj
!
!-----------------------------------------------------------------------
!  Check C-preprocessing options.
!-----------------------------------------------------------------------
!
!  Stop if more than one vertical closure scheme is selected.
!
      IF (Master.and.(ivmix.gt.1)) THEN
        WRITE (stdout,30)
  30    FORMAT (/,' CHECKDEFS - only one vertical closure scheme',      &
     &            ' is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more that one bottom stress formulation is selected.
!
      IF (Master.and.(ibbl.gt.1)) THEN
        WRITE (stdout,40)
  40    FORMAT (/,' CHECKDEFS - only one bottom stress formulation is', &
     &            ' allowed.')
        exit_flag=5
      END IF
!
!  Stop if no bottom stress formulation is selected.
!
      IF (Master.and.(ibbl.eq.0)) THEN
        WRITE (stdout,50)
  50    FORMAT (/,' CHECKDEFS - no bottom stress formulation is',       &
     &            ' selected.')
        exit_flag=5
      END IF
!
!  Stop if more than one biological module is selected.
!
      IF (Master.and.(ibiology.gt.1)) THEN
        WRITE (stdout,60)
  60    FORMAT (/,' CHECKDEFS - only one biology MODULE is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more that one model driver is selected.
!
      IF (Master.and.(idriver.gt.1)) THEN
        WRITE (stdout,70)
  70    FORMAT (/,' CHECKDEFS - only one model example is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more than one advection scheme has been activated.
!
      IF (Master.and.(ivelHadv.gt.1)) THEN
        WRITE (stdout,140) 'horizontal','momentum','ivelHadv =',ivelHadv
        exit_flag=5
      END IF
      IF (Master.and.(ivelVadv.gt.1)) THEN
        WRITE (stdout,140) 'vertical','momentum','ivelVadv =',ivelVadv
        exit_flag=5
      END IF
      IF (Master.and.(itrcHadv.gt.1)) THEN
        WRITE (stdout,140) 'horizontal','tracers','itrcHadv =',itrcHadv
        exit_flag=5
      END IF
      IF (Master.and.(itrcVadv.gt.1)) THEN
        WRITE (stdout,140) 'vertical','tracers','itrcVadv =',itrcVadv
        exit_flag=5
      END IF
 140  FORMAT (/,' CHECKDEFS - only one ',a,' advection scheme',         &
     &        /,13x,'is allowed for ',a,', ',a,1x,i1)
      IF (Master.and.(sbed_type.gt.1)) THEN
        WRITE (stdout,170)
 170    FORMAT (/,' CHECKDEFS - More than one sediment bed type   ',    &
     &            ' selected.')
        exit_flag=5
      END IF
      IF (Master.and.(sbed_load.gt.1)) THEN
        WRITE (stdout,180)
 180    FORMAT (/,' CHECKDEFS - More than one sediment bed load   ',    &
     &            ' selected.')
        exit_flag=5
      END IF
      RETURN
      END SUBROUTINE checkdefs
