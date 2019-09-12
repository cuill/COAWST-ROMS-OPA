      MODULE step3d_t_mod
!
!svn $Id: step3d_t.F 854 2017-07-18 23:28:45Z arango $
!=======================================================================
!  Copyright (c) 2002-2018 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine time-steps tracer equations.  Notice that advective    !
!  and diffusive terms are time-stepped differently. It applies the    !
!  corrector time-step for horizontal/vertical advection,  vertical    !
!  diffusion, nudging if necessary, and lateral boundary conditions.   !
!                                                                      !
!  Notice that at input the tracer arrays have:                        !
!                                                                      !
!    t(:,:,:,nnew,:)   m Tunits  n+1     horizontal/vertical diffusion !
!                                        terms plus source/sink terms  !
!                                        (biology, sediment), if any   !
!                                                                      !
!    t(:,:,:,3   ,:)   Tunits    n+1/2   advective terms and vertical  !
!                                        diffusion predictor step      !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: step3d_t
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE step3d_t (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
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
      CALL wclock_on (ng, iNLM, 35, 60, "ROMS/Nonlinear/step3d_t.F")
      CALL step3d_t_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nrhs(ng), nstp(ng), nnew(ng),                 &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % Hvom,                              &
     &                    GRID(ng) % z_r,                               &
     &                    MIXING(ng) % Akt,                             &
     &                    OCEAN(ng) % W,                                &
     &                    OCEAN(ng) % t)
      CALL wclock_off (ng, iNLM, 35, 116, "ROMS/Nonlinear/step3d_t.F")
      RETURN
      END SUBROUTINE step3d_t
!
!***********************************************************************
      SUBROUTINE step3d_t_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nrhs, nstp, nnew,                       &
     &                          pm, pn,                                 &
     &                          Hz, Huon, Hvom,                         &
     &                          z_r,                                    &
     &                          Akt,                                    &
     &                          W,                                      &
     &                          t)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_ncparam
      USE mod_scalars
      USE mod_sources
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
      USE hsimt_tvd_mod
      USE t3dbc_mod, ONLY : t3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
!
!  Local variable declarations.
!
      integer :: Isrc, Jsrc
      integer :: i, ibt, ic, is, itrc, j, k, ltrc
      real(r8), parameter :: eps = 1.0E-16_r8
      real(r8) :: cff, cff1, cff2, cff3
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: oHz
      real(r8):: rd, rkad, a1, b1, betad, sw, ru, rkau
      real(r8):: betau, epson, cc1, cc2, cc3
      real(r8), dimension(0:N(ng)) :: kaz, kaz_inverse, grad_k
      real(r8) :: my_maxbio(15)
!
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
!  Time-step horizontal advection term.
!-----------------------------------------------------------------------
!
!  Compute inverse thickness.
!
      DO k=1,N(ng)
        DO j=Jstrm2,Jendp2
          DO i=Istrm2,Iendp2
            oHz(i,j,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
      END DO
!
!  The MPDATA algorithm requires a three-point footprint, so exchange
!  boundary data on t(:,:,:,nnew,:) so other processes computed earlier
!  (horizontal diffusion, biology, or sediment) are accounted.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nnew,itrc))
        END DO
      END IF
!
!  Compute horizontal tracer advection fluxes.
!
      T_LOOP : DO itrc=1,NT(ng)
        K_LOOP : DO k=1,N(ng)
!   Third-order HSIMT-TVD horizontal advective fluxes.
          CALL hsimt_tvd_tile (ng, tile,                                &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          pm, pn,                                 &
     &                          Huon(:,:,k), Hvom(:,:,k),               &
     &                          oHz(:,:,k), t(:,:,k,3,itrc),            &
     &                          FX, FE)
!
!  Apply tracers point sources to the horizontal advection terms,
!  if any.
!
          IF (LuvSrc(ng).and.ANY(LtracerSrc(:,ng))) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
                IF (((IstrUm2.le.Isrc).and.(Isrc.le.Iendp3)).and.       &
     &              ((JstrVm2.le.Jsrc).and.(Jsrc.le.Jendp2i))) THEN
                  IF (LtracerSrc(itrc,ng)) THEN
                    FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                    &
     &                      SOURCES(ng)%Tsrc(is,k,itrc)
                  END IF
                END IF
              ELSE IF (INT(SOURCES(ng)%Dsrc(is)).eq.1) THEN
                IF (((IstrUm2.le.Isrc).and.(Isrc.le.Iendp2i)).and.      &
     &              ((JstrVm2.le.Jsrc).and.(Jsrc.le.Jendp3))) THEN
                  IF (LtracerSrc(itrc,ng)) THEN
                    FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                    &
     &                      SOURCES(ng)%Tsrc(is,k,itrc)
                  END IF
                END IF
              END IF
            END DO
          END IF
!
!  Time-step horizontal advection term.  Advective fluxes have units
!  of Tunits m3/s.  The new tracer has units of m Tunits.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=dt(ng)*pm(i,j)*pn(i,j)
              cff1=cff*(FX(i+1,j)-FX(i,j))
              cff2=cff*(FE(i,j+1)-FE(i,j))
              cff3=cff1+cff2
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff3
            END DO
          END DO
        END DO K_LOOP
      END DO T_LOOP
!
!-----------------------------------------------------------------------
!  Time-step vertical advection term.
!-----------------------------------------------------------------------
      cc1=0.25_r8
      cc2=0.5_r8
      cc3=1.0_r8/12.0_r8
      epson=1.0E-12_r8
      DO j=Jstr,Jend
        DO itrc=1,NT(ng)
          kaz(0)=0.0_r8
          kaz(N(ng))=0.0_r8
          kaz_inverse(0)=0.0_r8
          kaz_inverse(N(ng))=0.0_r8
          grad_k(0)=0.0_r8
          grad_k(N(ng))=0.0_r8
          DO i=Istr,Iend
            DO k=1,N(ng)-1
              kaz(k)=1.0_r8-abs(w(i,j,k)*pm(i,j)*pn(i,j)*dt(ng)/        &
     &                         (z_r(i,j,k+1)-z_r(i,j,k)))
              kaz_inverse(k)=1.0_r8/(kaz(k))
              grad_k(k)=t(i,j,k+1,3,itrc)-t(i,j,k,3,itrc)
            END DO
            DO k=1,N(ng)-1
              IF (k.eq.1.and.W(i,j,k).ge.0.0_r8) THEN
                FC(i,k)=W(i,j,k)*t(i,j,k,3,itrc)
              ELSE IF (k.eq.N(ng)-1.and.W(i,j,k).lt.0.0_r8) THEN
                FC(i,k)=W(i,j,k)*t(i,j,k+1,3,itrc)
              ELSE
                IF (W(i,j,k).ge.0) THEN
                  IF (abs(grad_k(k)).le.epson) THEN
                    rd=0.0_r8
                    rkad=0.0_r8
                  ELSE
                    rd=grad_k(k-1)/grad_k(k)
                    rkad=kaz(k-1)*kaz_inverse(k)
!                   rkad=kaz(k-1)*kaz_inverse(k)*W(i,j,k-1)/W(i,j,k)
                  END IF
                  a1= cc1*kaz(k)+cc2-cc3*kaz_inverse(k)
                  b1=-cc1*kaz(k)+cc2+cc3*kaz_inverse(k)
                  betad=a1+b1*rd
                  sw=t(i,j,k,3,itrc)+                                   &
     &               0.5_r8*MAX(0.0_r8,                                 &
     &                          MIN(2.0_r8,2.0_r8*rd*rkad,betad))*      &
     &                grad_k(k)*kaz(k)
                ELSE
                  IF (abs(grad_k(k)).le.epson) THEN
                    ru=0.0_r8
                    rkau=0.0_r8
                  ELSE
                    ru=grad_k(k+1)/grad_k(k)
                    rkau=kaz(k+1)*kaz_inverse(k)
!                   rkau=kaz(k+1)*kaz_inverse(k)*W(i,j,k+1)/W(i,j,k)
                  END IF
                  a1= cc1*kaz(k)+cc2-cc3*kaz_inverse(k)
                  b1=-cc1*kaz(k)+cc2+cc3*kaz_inverse(k)
                  betau=a1+b1*ru
                  sw=t(i,j,k+1,3,itrc)-                                 &
     &               0.5_r8*MAX(0.0_r8,                                 &
     &                          MIN(2.0_r8,2.0_r8*ru*rkau,betau))*      &
     &               grad_k(k)*kaz(k)
                END IF 
                FC(i,k)=W(i,j,k)*sw
              END IF
            END DO
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Time-step vertical advection term (Tunits).
!
          DO i=Istr,Iend
            CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
          END DO
!
!  Apply mass point sources (volume vertical influx), if any.
!
          IF (LwSrc(ng).and.ANY(LtracerSrc(:,ng))) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (LtracerSrc(itrc,ng).and.                              &
     &            ((IstrUm2.le.Isrc).and.(Isrc.le.Iendp2i)).and.        &
     &            (j.eq.Jsrc)) THEN
                DO k=1,N(ng)-1
                  FC(Isrc,k)=FC(Isrc,k)+0.5_r8*                         &
     &                       (SOURCES(ng)%Qsrc(is,k  )*                 &
     &                        SOURCES(ng)%Tsrc(is,k  ,itrc)+            &
     &                        SOURCES(ng)%Qsrc(is,k+1)*                 &
     &                        SOURCES(ng)%Tsrc(is,k+1,itrc))
                END DO
              END IF
            END DO
          END IF
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*oHz(i,j,k)
            END DO
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Time-step vertical diffusion term.
!-----------------------------------------------------------------------
!
        DO itrc=1,NT(ng)
          ltrc=MIN(NAT,itrc)
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hz(i,j,k  )-                                 &
     &                dt(ng)*Akt(i,j,k-1,ltrc)*oHz(i,j,k  )
              CF(i,k)=cff1*Hz(i,j,k+1)-                                 &
     &                dt(ng)*Akt(i,j,k+1,ltrc)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hz(i,j,k)+Hz(i,j,k+1))+                     &
     &                dt(ng)*Akt(i,j,k,ltrc)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(t(i,j,k+1,nnew,itrc)-t(i,j,k,nnew,itrc)-     &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)*Akt(i,j,k,ltrc)
              cff1=dt(ng)*oHz(i,j,k)*(DC(i,k)-DC(i,k-1))
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff1
            END DO
          END DO
        END DO
      END DO
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions and, if appropriate, nudge
!  to tracer data and apply Land/Sea mask.
!-----------------------------------------------------------------------
!
      ic=0
!  Initialize tracer counter index. The "tclm" array is only allocated
!  to the NTCLM fields that need to be processed. This is done to
!  reduce memory.
!
!
      DO itrc=1,NT(ng)
!
!  Set compact reduced memory tracer index for nudging coefficients and
!  climatology arrays.
!
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
        END IF
!
!  Set lateral boundary conditions.
!
        CALL t3dbc_tile (ng, tile, itrc, ic,                            &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   t)
!
!  Nudge towards tracer climatology.
!
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+                  &
     &                             dt(ng)*                              &
     &                             CLIMA(ng)%Tnudgcof(i,j,k,ic)*        &
     &                             (CLIMA(ng)%tclm(i,j,k,ic)-           &
     &                              t(i,j,k,nnew,itrc))
              END DO
            END DO
          END DO
        END IF
!
!  Apply periodic boundary conditions.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nnew,itrc))
        END IF
      END DO
      RETURN
      END SUBROUTINE step3d_t_tile
      END MODULE step3d_t_mod
