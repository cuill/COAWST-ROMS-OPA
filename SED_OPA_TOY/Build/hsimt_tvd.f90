      MODULE hsimt_tvd_mod
!
!======================= Hui Wu, Tarandeep Kalra, and John C. Warner====
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!This routine computes anti-diffusive tracer flux based on HSIMT-TVD   !
!by Wu and Zhu (2010). This routine is for personal test only currently!
!                                                                      !
!  On Output: FX, FE                                                   !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Hui Wu and Jianrong Zhu (2010), Advection scheme with 3rd         !
!    high-order spatial interpolation at the middle temporal level     !
!    and its application to saltwater intrusion in the Changjiang      !
!    Estuary, Ocean Modelling 33, 33-51.                               !
!   Please contact Hui Wu (hwusklec@gmail.com) if have any questions   !
!                                                                      !
!=======================================================================
!
      implicit none
      PUBLIC :: hsimt_tvd_tile
      CONTAINS
!
!***********************************************************************
      SUBROUTINE hsimt_tvd_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           pm, pn,                                &
     &                           Huon_k, Hvom_k, oHz_k, t_k,            &
     &                           FX, FE)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Huon_k(LBi:,LBj:)
      real(r8), intent(in) :: Hvom_k(LBi:,LBj:)
      real(r8), intent(in) :: oHz_k(IminS:,JminS:)
      real(r8), intent(in) :: t_k(LBi:,LBj:)
      real(r8), intent(out) :: FX(IminS:,JminS:)
      real(r8), intent(out) :: FE(IminS:,JminS:)
!
!  Local variable declarations.
!
      integer  :: i, is, j, k, ii, jj
      real(r8) :: cc1, cc2, cc3
      real(r8) :: sw_xi, rl, rkal, a1, b1, betal, rt, rkar, betar
      real(r8) :: sw_eta, rd, rkad, betad, ru, rkau, betau
      real(r8) :: cff, cff1, cff2, epson
      real(r8), dimension(IminS:ImaxS) :: kax, kax_inverse
      real(r8), dimension(IminS:ImaxS) :: grad_x
      real(r8), dimension(JminS:JmaxS) :: grad_y
      real(r8), dimension(JminS:JmaxS) :: kay, kay_inverse
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
!************Declare some constants locally***************************
      cc1=0.25_r8
      cc2=0.5_r8
      cc3=1.0_r8/12.0_r8
      epson=1.0E-12_r8
!
      DO j=Jstr,Jend
        DO i=Istr-1,Iend+2
          grad_x(i)=(t_k(i,j)-t_k(i-1,j))
          cff=0.125_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))*         &
     &        (oHz_k(i-1,j)+oHz_k(i,j))
          kax(i)=(1.0_r8-abs(Huon_k(i,j)*dt(ng)*cff))
        END DO
        IF (.not.EWperiodic(ng)) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            IF (Huon_k(Istr,j).ge.0.0_r8) THEN
              grad_x(Istr-1)=0.0_r8
              kax(Istr-1)=0.0_r8
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            IF (Huon_k(Iend+1,j).lt.0.0_r8) THEN
              grad_x(Iend+2)=0.0_r8
              kax(Iend+2)=0.0_r8
            END IF
          END IF
        END IF
        DO i=Istr,Iend+1
          IF (kax(i).le.epson) THEN
            kax_inverse(i)=0.0_r8
          ELSE
            kax_inverse(i)=1.0_r8/MAX(kax(i),epson)
          END IF
          IF (Huon_k(i,j).ge.0.0_r8) THEN
            IF (abs(grad_x(i)).le.epson) THEN
              rl=0.0_r8
              rkal=0.0_r8
            ELSE
              rl=grad_x(i-1)/(grad_x(i))
              rkal=kax(i-1)*kax_inverse(i)
!             cff1=Huon_k(i-1,j)*0.25_r8*(oHz_k(i-2,j)+oHz_k(i-1,j))*(pn(i-2,j)+pn(i-1,j))
!             cff2=Huon_k(i,j)*0.25_r8*(oHz_k(i-1,j)+oHz_k(i,j))*(pn(i-1,j)+pn(i,j))
!             rkal=kax(i-1)*kax_inverse(i)*cff1/cff2
            END IF
            a1= cc1*kax(i)+cc2-cc3*kax_inverse(i)
            b1=-cc1*kax(i)+cc2+cc3*kax_inverse(i)
            betal=a1+b1*rl
            cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rl*rkal,betal))*    &
     &                                  grad_x(i)*kax(i)
            sw_xi=t_k(i-1,j)+cff
          ELSE
            IF (abs(grad_x(i)).le.epson) THEN
              rt=0.0_r8
              rkar=0.0_r8
            ELSE
              rt=grad_x(i+1)/(grad_x(i))
              rkar=kax(i+1)*kax_inverse(i)
!             cff1=Huon_k(i+1,j)*0.25_r8*(oHz_k(i,j)+oHz_k(i+1,j))*(pn(i,j)+pn(i+1,j))
!             cff2=Huon_k(i,j)*0.25_r8*(oHz_k(i-1,j)+oHz_k(i,j))*(pn(i-1,j)+pn(i,j))
!             rkar=kax(i+1)*kax_inverse(i)*cff1/cff2
            END IF
            a1= cc1*kax(i)+cc2-cc3*kax_inverse(i)
            b1=-cc1*kax(i)+cc2+cc3*kax_inverse(i)
            betar=a1+b1*rt
            cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rt*rkar,betar))*    &
     &                                grad_x(i)*kax(i)
            sw_xi=t_k(i,j)-cff
          END IF
          FX(i,j)=sw_xi*huon_k(i,j)
        END DO
      END DO
!
      DO i=Istr,Iend
        DO j=Jstr-1,Jend+2
          grad_y(j)=(t_k(i,j)-t_k(i,j-1))
          cff=0.125_r8*(pn(i,j)+pn(i,j-1))*(pm(i,j)+pm(i,j-1))*         &
     &        (oHz_k(i,j)+oHz_k(i,j-1))
          kay(j)=(1.0_r8-abs(Hvom_k(i,j)*dt(ng)*cff))
        END DO
        IF (.not.NSperiodic(ng)) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            IF (Hvom_k(i,Jstr).ge.0.0_r8) THEN
              grad_y(Jstr-1)=0.0_r8
              kay(Jstr-1)=0.0_r8
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            IF (Hvom_k(i,Jend+1).lt.0.0_r8) THEN
              grad_y(Jend+2)=0.0_r8
              kay(Jend+2)=0.0_r8
            END IF
          END IF
        END IF
        DO j=Jstr,Jend+1
          IF (kay(j).le.epson) THEN
            kay_inverse(j)=0.0_r8
          ELSE
            kay_inverse(j)=1.0_r8/MAX(kay(j),epson)
          END IF
          IF (Hvom_k(i,j).ge.0.0_r8) THEN
            IF (abs(grad_y(j)).le.epson) THEN
              rd=0.0_r8
              rkad=0.0_r8
            ELSE
              rd=grad_y(j-1)/grad_y(j)
              rkad=kay(j-1)*kay_inverse(j)
!             cff1=Hvom_k(i,j-1)*0.25_r8*(oHz_k(i,j-2)+oHz_k(i,j-1))*(pm(i,j-2)+pm(i,j-1))
!             cff2=Hvom_k(i,j)*0.25_r8*(oHz_k(i,j-1)+oHz_k(i,j))*(pm(i,j-1)+pm(i,j))
!             rkad=kay(j-1)*kay_inverse(j)*cff1/cff2
            END IF
            a1= cc1*kay(j)+cc2-cc3*kay_inverse(j)
            b1=-cc1*kay(j)+cc2+cc3*kay_inverse(j)
            betad=a1+b1*rd
            cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rd*rkad,betad))*    &
     &                              grad_y(j)*kay(j)
            sw_eta=t_k(i,j-1)+cff
          ELSE
            IF (abs(grad_y(j)).le.epson) THEN
              ru=0.0_r8
              rkau=0.0_r8
            ELSE
              ru=grad_y(j+1)/(grad_y(j))
              rkau=kay(j+1)*kay_inverse(j)
!             cff1=Hvom_k(i,j+1)*0.25_r8*(oHz_k(i,j+1)+oHz_k(i,j))*(pm(i,j+1)+pm(i,j))
!             cff2=Hvom_k(i,j)*0.25_r8*(oHz_k(i,j-1)+oHz_k(i,j))*(pm(i,j-1)+pm(i,j))
!             rkau=kay(j+1)*kay_inverse(j)*cff1/cff2
            END IF
            a1= cc1*kay(j)+cc2-cc3*kay_inverse(j)
            b1=-cc1*kay(j)+cc2+cc3*kay_inverse(j)
            betau=a1+b1*ru
            cff=0.5*max(0.0_r8,min(2.0_r8,2.0_r8*ru*rkau,betau))*       &
     &                            grad_y(j)*kay(j)
            sw_eta=t_k(i,j)-cff
          END IF
          FE(i,j)=sw_eta*hvom_k(i,j)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE hsimt_tvd_tile
      END MODULE hsimt_tvd_mod
 
