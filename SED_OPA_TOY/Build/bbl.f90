      MODULE bbl_mod
!
!svn $Id: bbl.F 830 2017-01-24 21:21:11Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2018 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes bottom momentum stress via a bottom boundary  !
!  layer formulation.                                                  !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: bblm
      CONTAINS
      SUBROUTINE bblm (ng, tile)
!
!svn $Id: ssw_bbl.h 854 2017-07-18 23:28:45Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2018 The ROMS/TOMS Group        Chris Sherwood   !
!    Licensed under a MIT/X style license               Rich Signell   !
!    See License_ROMS.txt                             John C. Warner   !
!=======================================================================
!                                                                      !
!  This routine compute bottom stresses for the case when the wave     !
!  solution in the wave boundary layer is based on a  2-layer eddy     !
!  viscosity that is linear increasing above Zo and constant above     !
!  Z1.                                                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Styles, R. and S.M. glenn,  2000: Modeling stratified wave and      !
!    current bottom boundary layers in the continental shelf, JGR,     !
!    105, 24119-24139.                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_bbl
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbed
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL wclock_on (ng, iNLM, 37, 44, "ROMS/Nonlinear/ssw_bbl.h")
      CALL bblm_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                nrhs(ng),                                         &
     &                GRID(ng) % h,                                     &
     &                GRID(ng) % z_r,                                   &
     &                GRID(ng) % z_w,                                   &
     &                GRID(ng) % angler,                                &
     &                GRID(ng) % ZoBot,                                 &
     &                FORCES(ng) % Hwave,                               &
     &                FORCES(ng) % Dwave,                               &
     &                FORCES(ng) % Pwave_bot,                           &
     &                SEDBED(ng) % bottom,                              &
     &                OCEAN(ng) % rho,                                  &
     &                OCEAN(ng) % u,                                    &
     &                OCEAN(ng) % v,                                    &
     &                OCEAN(ng) % zeta,                                 &
     &                BBL(ng) % Iconv,                                  &
     &                BBL(ng) % Ubot,                                   &
     &                BBL(ng) % Vbot,                                   &
     &                BBL(ng) % Ur,                                     &
     &                BBL(ng) % Vr,                                     &
     &                BBL(ng) % bustrc,                                 &
     &                BBL(ng) % bvstrc,                                 &
     &                BBL(ng) % bustrw,                                 &
     &                BBL(ng) % bvstrw,                                 &
     &                BBL(ng) % bustrcwmax,                             &
     &                BBL(ng) % bvstrcwmax,                             &
     &                FORCES(ng) % bustr,                               &
     &                FORCES(ng) % bvstr)
      CALL wclock_off (ng, iNLM, 37, 87, "ROMS/Nonlinear/ssw_bbl.h")
      RETURN
      END SUBROUTINE bblm
!
!***********************************************************************
      SUBROUTINE bblm_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs,                                       &
     &                      h, z_r, z_w, angler, ZoBot,                 &
     &                      Hwave,                                      &
     &                      Dwave, Pwave_bot,                           &
     &                      bottom, rho, u, v,                          &
     &                      zeta,                                       &
     &                      Iconv,                                      &
     &                      Ubot, Vbot, Ur, Vr,                         &
     &                      bustrc, bvstrc,                             &
     &                      bustrw, bvstrw,                             &
     &                      bustrcwmax, bvstrcwmax,                     &
     &                      bustr, bvstr)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_sediment
!
      USE bc_2d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
      integer, intent(inout) :: Iconv(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: ZoBot(LBi:,LBj:)
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Pwave_bot(LBi:,LBj:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(out) :: Ubot(LBi:,LBj:)
      real(r8), intent(out) :: Vbot(LBi:,LBj:)
      real(r8), intent(out) :: Ur(LBi:,LBj:)
      real(r8), intent(out) :: Vr(LBi:,LBj:)
      real(r8), intent(out) :: bustrc(LBi:,LBj:)
      real(r8), intent(out) :: bvstrc(LBi:,LBj:)
      real(r8), intent(out) :: bustrw(LBi:,LBj:)
      real(r8), intent(out) :: bvstrw(LBi:,LBj:)
      real(r8), intent(out) :: bustrcwmax(LBi:,LBj:)
      real(r8), intent(out) :: bvstrcwmax(LBi:,LBj:)
      real(r8), intent(out) :: bustr(LBi:,LBj:)
      real(r8), intent(out) :: bvstr(LBi:,LBj:)
!
!  Local variable declarations.
!
      logical :: ITERATE
      integer :: Iter, i, j, k
      real(r8), parameter :: eps = 1.0E-10_r8
      real(r8) :: Kbh, Kbh2, Kdh
      real(r8) :: taucr, wsedr, tstar, coef_st
      real(r8) :: coef_b1, coef_b2, coef_b3, d0
      real(r8) :: dolam, dolam1, doeta1, doeta2, fdo_etaano
      real(r8) :: lamorb, lamanorb
      real(r8) :: m_ubr, m_wr, m_ucr, m_zr, m_phicw, m_kb
      real(r8) :: m_ustrc, m_ustrwm, m_ustrr, m_fwc, m_zoa
      real(r8) :: zo
      real(r8) :: Kb, Kdelta, Ustr
      real(r8) :: anglec, anglew
      real(r8) :: cff, cff1, cff2, cff3, og, fac, fac1, fac2
      real(r8) :: sg_ab, sg_abokb, sg_a1, sg_b1, sg_chi, sg_c1, d50
      real(r8) :: sg_epsilon, ssw_eta, sg_fofa, sg_fofb, sg_fofc, sg_fwm
      real(r8) :: sg_kbs, ssw_lambda, sg_mu, sg_phicw, sg_ro, sg_row
      real(r8) :: sg_shdnrm, sg_shld, sg_shldcr, sg_scf, rhos, sg_star
      real(r8) :: sg_ub, sg_ubokur, sg_ubouc, sg_ubouwm, sg_ur
      real(r8) :: sg_ustarc, sg_ustarcw, sg_ustarwm, sg_znot, sg_znotp
      real(r8) :: sg_zr, sg_zrozn, sg_z1, sg_z1ozn, sg_z2, twopi, z1, z2
      real(r8) :: zoMIN, zoMAX
      real(r8) :: coef_fd
      real(r8), parameter :: absolute_zoMIN = 5.0d-5  ! in Harris-Wiberg
      real(r8), parameter ::  Cd_fd = 0.5_r8
      real(r8), parameter :: K1 = 0.6666666666_r8     ! Coefficients for
      real(r8), parameter :: K2 = 0.3555555555_r8     ! explicit
      real(r8), parameter :: K3 = 0.1608465608_r8     ! wavenumber
      real(r8), parameter :: K4 = 0.0632098765_r8     ! calculation
      real(r8), parameter :: K5 = 0.0217540484_r8     ! (Dean and
      real(r8), parameter :: K6 = 0.0065407983_r8     !  Dalrymple, 1991)
      real(r8), parameter :: coef_a1=0.095_r8         ! Coefficients for
      real(r8), parameter :: coef_a2=0.442_r8         ! ripple predictor
      real(r8), parameter :: coef_a3=2.280_r8         ! (Wiberg-Harris)
      real(r8), parameter :: ar = 0.267_r8         ! Nielsen (1992)
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ab
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Fwave_bot
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauc
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taucwmax
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ur_sg
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vr_sg
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ub
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ucur
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Umag
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vcur
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Zr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: phic
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: phicw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rheight
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rlength
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: u100
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: znot
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: znotc
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoN
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoST
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoBF
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoDEF
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoBIO
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Set currents above the bed.
!-----------------------------------------------------------------------
!
      twopi=2.0_r8*pi
      DO j=JstrV-1,Jend+1
        DO i=IstrU-1,Iend+1
          Zr(i,j)=z_r(i,j,1)-z_w(i,j,0)
          Ur_sg(i,j)=u(i,j,1,nrhs)
          Vr_sg(i,j)=v(i,j,1,nrhs)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute bottom stresses.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!
!  Set bed wave orbital velocity and excursion amplitude.  Use data
!  from wave models (SWAN) or use Dean and Dalrymple (1991) 6th-degree
!  polynomial to approximate wave number on shoaling water.
          Fwave_bot(i,j)=twopi/MAX(Pwave_bot(i,j),1.0_r8)
          Kdh=(h(i,j)+zeta(i,j,nrhs))*Fwave_bot(i,j)**2/g
          Kbh2=Kdh*Kdh+                                                 &
     &         Kdh/(1.0_r8+Kdh*(K1+Kdh*(K2+Kdh*(K3+Kdh*(K4+             &
     &              Kdh*(K5+K6*Kdh))))))
          Kbh=SQRT(Kbh2)
          Ab(i,j)=0.5_r8*Hwave(i,j)/SINH(Kbh)+eps
          Ub(i,j)=Fwave_bot(i,j)*Ab(i,j)+eps
!
!  Compute bottom current magnitude at RHO-points.
!
          Ucur(i,j)=0.5_r8*(Ur_sg(i,j)+Ur_sg(i+1,j))
          Vcur(i,j)=0.5_r8*(Vr_sg(i,j)+Vr_sg(i,j+1))
          Umag(i,j)=SQRT(Ucur(i,j)*Ucur(i,j)+Vcur(i,j)*Vcur(i,j)+eps)
!
!  Compute angle between currents and waves (radians)
!
          IF (Ucur(i,j).eq.0.0_r8) THEN
            phic(i,j)=0.5_r8*pi*SIGN(1.0_r8,Vcur(i,j))
          ELSE
            phic(i,j)=ATAN2(Vcur(i,j),Ucur(i,j))
          ENDIF
          phicw(i,j)=1.5_r8*pi-Dwave(i,j)-phic(i,j)-angler(i,j)
        END DO
      END DO
!
!  Loop over RHO points.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!
!  Load sediment properties, stresses, and roughness from previous time
!  step (stresses are in m2 s-2).
!
          d50=bottom(i,j,isd50)
          rhos=bottom(i,j,idens)/(rho(i,j,1)+1000.0_r8)
          wsedr=bottom(i,j,iwsed)
          Taucr=bottom(i,j,itauc)
          Tauc(i,j)=SQRT(bustrc(i,j)**2+bvstrc(i,j)**2)
          Tauw(i,j)=SQRT(bustrw(i,j)**2+bvstrw(i,j)**2)
          Taucwmax(i,j)=SQRT( bustrcwmax(i,j)**2+bvstrcwmax(i,j)**2)
!
          rheight(i,j)=bottom(i,j,irhgt)
          rlength(i,j)=bottom(i,j,irlen)
          zoMAX=0.9_r8*Zr(i,j)
          zoMIN=MAX(absolute_zoMIN,2.5_r8*d50/30.0_r8)
!
!  Initialize arrays.
!
          zoN(i,j)=MIN(MAX(2.5_r8*d50/30.0_r8, zoMIN ),zoMAX)
          zoST(i,j)=0.0_r8
          zoBF(i,j)=0.0_r8
          zoBIO(i,j)=0.0_r8
          zoDEF(i,j)=ZoBot(i,j)
          IF (zoDEF(i,j).lt.absolute_zoMIN) THEN
            zoDEF(i,j)=absolute_zoMIN
            IF (Master) THEN
              PRINT *, ' Warning: default zo < 0.05 mm, replaced with:',&
     &                 zoDEF
            END IF
          END IF
          zo=zoDEF(i,j)
!
!  Compute stresses.
!
!  Default stress calcs for pure currents
!
          zo=MIN(MAX(zo,zoMIN),zoMAX)
!
          cff1=vonKar/LOG(Zr(i,j)/zo)
          cff2=MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
          Tauc(i,j)=cff2*Umag(i,j)*Umag(i,j)
          Tauw(i,j)=0.0_r8
          Taucwmax(i,j)=Tauc(i,j)
          znot(i,j)=zo
          znotc(i,j)=zo
!
          IF ((Umag(i,j).le.eps).and.(Ub(i,j).gt.eps)) THEN
!
!  Pure waves - use wave friction factor approach from Madsen
!  (1994, eqns 32-33).
!
            sg_abokb=Ab(i,j)/(30.0_r8*zo)
            sg_fwm=0.3_r8
            IF ((sg_abokb.gt.0.2_r8).and.(sg_abokb.le.100.0_r8)) THEN
              sg_fwm=EXP(-8.82_r8+7.02_r8*sg_abokb**(-0.078_r8))
            ELSE IF (sg_abokb.gt.100.0_r8)THEN
              sg_fwm=EXP(-7.30_r8+5.61_r8*sg_abokb**(-0.109_r8))
            END IF
            Tauc(i,j)= 0.0_r8
            Tauw(i,j)= 0.5_r8*sg_fwm*Ub(i,j)*Ub(i,j)
            Taucwmax(i,j)=Tauw(i,j)
            znot(i,j)=zo
            znotc(i,j)=zo
          ELSE IF ((Umag(i,j).gt.0.0_r8).and.(Ub(i,j).gt.eps).and.      &
     &             ((Zr(i,j)/zo).le.1.0_r8)) THEN
!
!  Waves and currents, but zr <= zo.
!
            IF (Master) THEN
              PRINT *,' Warning: w-c calcs ignored because zr <= zo'
            END IF
          ELSE IF ((Umag(i,j).gt.0.0_r8).and.(Ub(i,j).gt.eps).and.      &
     &             ((Zr(i,j)/zo).gt.1.0_r8)) THEN
!
!  Waves and currents, zr > zo.
!
            m_ubr=Ub(i,j)
            m_wr=Fwave_bot(i,j)
            m_ucr=Umag(i,j)
            m_zr=Zr(i,j)
            m_phicw=phicw(i,j)
            m_kb=30.0_r8*zo
            CALL madsen94 (m_ubr, m_wr, m_ucr,                          &
     &                     m_zr, m_phicw, m_kb,                         &
     &                     m_ustrc, m_ustrwm, m_ustrr, m_fwc, m_zoa)
            Tauc(i,j)=m_ustrc*m_ustrc
            Tauw(i,j)=m_ustrwm*m_ustrwm
            Taucwmax(i,j)=m_ustrr*m_ustrr
            znotc(i,j)=min( m_zoa, zoMAX )
            u100(i,j)=(m_ustrc/vonKar)*LOG(1.0_r8/m_zoa)
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute kinematic bottom stress components due current and wind-
!  induced waves.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          anglec=Ur_sg(i,j)/(0.5*(Umag(i-1,j)+Umag(i,j)))
          bustr(i,j)=0.5_r8*(Tauc(i-1,j)+Tauc(i,j))*anglec
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          anglec=Vr_sg(i,j)/(0.5_r8*(Umag(i,j-1)+Umag(i,j)))
          bvstr(i,j)=0.5_r8*(Tauc(i,j-1)+Tauc(i,j))*anglec
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          anglec=Ucur(i,j)/Umag(i,j)
          anglew=COS(1.5_r8*pi-Dwave(i,j)-angler(i,j))
          bustrc(i,j)=Tauc(i,j)*anglec
          bustrw(i,j)=Tauw(i,j)*anglew
          bustrcwmax(i,j)=Taucwmax(i,j)*anglew
          Ubot(i,j)=Ub(i,j)*anglew
          Ur(i,j)=Ucur(i,j)
!
          anglec=Vcur(i,j)/Umag(i,j)
          anglew=SIN(1.5_r8*pi-Dwave(i,j)-angler(i,j))
          bvstrc(i,j)=Tauc(i,j)*anglec
          bvstrw(i,j)=Tauw(i,j)*anglew
          bvstrcwmax(i,j)=Taucwmax(i,j)*anglew
          Vbot(i,j)=Ub(i,j)*anglew
          Vr(i,j)=Vcur(i,j)
!
          bottom(i,j,ibwav)=Ab(i,j)
          bottom(i,j,irhgt)=rheight(i,j)
          bottom(i,j,irlen)=rlength(i,j)
          bottom(i,j,izdef)=zoDEF(i,j)
          bottom(i,j,izapp)=znotc(i,j)
          bottom(i,j,izNik)=zoN(i,j)
          bottom(i,j,izbio)=zoBIO(i,j)
          bottom(i,j,izbfm)=zoBF(i,j)
          bottom(i,j,izbld)=zoST(i,j)
          bottom(i,j,izwbl)=znot(i,j)
        END DO
      END DO
!
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstr)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustrc)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstrc)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustrw)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstrw)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustrcwmax)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstrcwmax)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Ubot)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Vbot)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Ur)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Vr)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,ibwav))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,irhgt))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,irlen))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izdef))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izapp))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izNik))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izbio))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izbfm))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izbld))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izwbl))
      RETURN
      END SUBROUTINE bblm_tile
      SUBROUTINE madsen94 (ubr, wr, ucr, zr, phiwc, kN,                 &
     &                     ustrc, ustrwm, ustrr, fwc, zoa)
!
!=======================================================================
!                                                                      !
!  Grant-Madsen model from Madsen (1994).                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ubr     Rep. wave-orbital velocity amplitude outside WBL (m/s).  !
!     wr      Rep. angular wave frequency,  2* pi/T (rad/s).           !
!     ucr     Current velocity at height zr (m/s).                     !
!     zr      Reference height for current velocity (m).               !
!     phiwc   Angle between currents and waves at zr (radians).        !
!     kN      Bottom roughness height, like Nikuradse k, (m).          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ustrc   Current friction velocity, u*c (m/s).                    !
!     ustrwm  Wave maximum friction velocity, u*wm (m/s).              !
!     ustrr   Wave-current combined friction velocity, u*r (m/s).      !
!     fwc     Wave friction factor (nondimensional).                   !
!     zoa     Apparent bottom roughness (m).                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      real(r8), intent(in) ::  ubr, wr, ucr, zr, phiwc, kN
      real(r8), intent(out) ::  ustrc, ustrwm, ustrr, fwc, zoa
!
!  Local variable declarations.
!
      integer, parameter :: MAXIT = 20
      integer :: i, nit
      integer :: iverbose = 1
      real(r8) :: bigsqr, cosphiwc, cukw, diff, lndw, lnln, lnzr
      real(r8) :: phicwc, zo
      real(r8) :: dval = 99.99_r8
      real(r8), dimension(MAXIT) :: Cmu
      real(r8), dimension(MAXIT) :: dwc
      real(r8), dimension(MAXIT) :: fwci
      real(r8), dimension(MAXIT) :: rmu
      real(r8), dimension(MAXIT) :: ustrci
      real(r8), dimension(MAXIT) :: ustrr2
      real(r8), dimension(MAXIT) :: ustrwm2
!
!-----------------------------------------------------------------------
!  Compute bottom friction velocities and roughness.
!-----------------------------------------------------------------------
!
!  Set special default values.
!
      ustrc=dval
      ustrwm=dval
      ustrr=dval
      fwc=0.4_r8
      zoa=kN/30.0_r8
      phicwc=phiwc
      zo = kN/30.0_r8
      IF (ubr.le.0.01_r8) THEN
        IF (ucr.le. 0.01_r8) THEN          ! no waves or currents
          ustrc=0.0_r8
          ustrwm=0.0_r8
          ustrr=0.0_r8
          RETURN
        END IF
        ustrc=ucr*vonKar/LOG(zr/zo)        ! no waves
        ustrwm=0.0_r8
        ustrr=ustrc
        RETURN
      END IF
!
!  Iterate to compute friction velocities, roughness, and wave friction
!  factor.  Notice that the computation of the wave friction factor
!  has been inlined for efficiency.
!
      cosphiwc=ABS(COS(phiwc))
      rmu(1)=0.0_r8
      Cmu(1)=1.0_r8
      cukw=Cmu(1)*ubr/(kN*wr)
      IF ((cukw.gt.0.2_r8).and.(cukw.le.100.0_r8)) THEN       ! Eq 32/33
        fwci(1)=Cmu(1)*EXP(7.02_r8*cukw**(-0.078_r8)-8.82_r8)
      ELSE IF ((cukw.gt.100.).and.(cukw.le.10000.0_r8)) THEN
        fwci(1)=Cmu(1)*EXP(5.61_r8*cukw**(-0.109_r8)-7.30_r8)
      ELSE IF (cukw.gt.10000.0_r8 ) THEN
        fwci(1)=Cmu(1)*EXP(5.61_r8*10000.0_r8**(-0.109_r8)-7.30_r8)
      ELSE
        fwci(1)=Cmu(1)*0.43_r8
      END IF
      ustrwm2(1)=0.5_r8*fwci(1)*ubr*ubr                       ! Eq 29
      ustrr2(1)=Cmu(1)*ustrwm2(1)                             ! Eq 26
      ustrr=SQRT(ustrr2(1))
      IF (cukw.ge.8.0_r8) THEN
        dwc(1)=2.0_r8*vonKar*ustrr/wr                         ! Eq 36
      ELSE
        dwc(1)=kN
      END IF
      lnzr=LOG(zr/dwc(1))
      lndw=LOG(dwc(1)/zo)
      lnln=lnzr/lndw
      bigsqr=-1.0_r8+SQRT(1.0_r8+((4.0_r8*vonKar*lndw)/                 &
     &                            (lnzr*lnzr))*ucr/ustrr)
      ustrci(1)=0.5_r8*ustrr*lnln*bigsqr
!
      i=1
      diff=1.0_r8
      DO WHILE ((i.lt.MAXIT).and.(diff.gt.0.000005_r8))
        i=i+1
        rmu(i)=ustrci(i-1)*ustrci(i-1)/ustrwm2(i-1)
        Cmu(i)=SQRT(1.0_r8+                                             &
     &              2.0_r8*rmu(i)*cosphiwc+rmu(i)*rmu(i))     ! Eq 27
        cukw=Cmu(i)*ubr/(kN*wr)
        IF ((cukw.gt.0.2_r8).and.(cukw.le.100.0_r8)) THEN     ! Eq 32/33
          fwci(i)=Cmu(i)*EXP(7.02_r8*cukw**(-0.078_r8)-8.82_r8)
        ELSE IF ((cukw.gt.100.).and.(cukw.le.10000.0_r8)) THEN
          fwci(i)=Cmu(i)*EXP(5.61_r8*cukw**(-0.109_r8)-7.30_r8)
        ELSE IF (cukw.gt.10000.0_r8 ) THEN
          fwci(i)=Cmu(i)*EXP(5.61_r8*10000.0_r8**(-0.109_r8)-7.30_r8)
        ELSE
          fwci(i)=Cmu(i)*0.43_r8
        END IF
        ustrwm2(i)=0.5_r8*fwci(i)*ubr*ubr                     ! Eq 29
        ustrr2(i)=Cmu(i)*ustrwm2(i)                           ! Eq 26
        ustrr=SQRT(ustrr2(i))
        IF (cukw.ge.8.0_r8) THEN
          dwc(i)=2.0_r8*vonKar*ustrr/wr                       ! Eq 36
        ELSE
          dwc(i)=kN
        END IF
        lnzr=LOG(zr/dwc(i))
        lndw=LOG(dwc(i)/zo)
        lnln=lnzr/lndw
        bigsqr=-1.0_r8+SQRT(1.0_r8+((4.0_r8*vonKar*lndw)/               &
     &                              (lnzr*lnzr))*ucr/ustrr)
        ustrci(i)=0.5_r8*ustrr*lnln*bigsqr                    ! Eq 38
        diff=ABS((fwci(i)-fwci(i-1))/fwci(i))
      END DO
      ustrwm=SQRT(ustrwm2(i))
      ustrc=ustrci(i)
      ustrr=SQRT(ustrr2(i))
      phicwc=phiwc
      zoa=EXP(LOG(dwc(i))-(ustrc/ustrr)*LOG(dwc(i)/zo))       ! Eq 11
      fwc=fwci(i)
      RETURN
      END SUBROUTINE madsen94
      END MODULE bbl_mod
