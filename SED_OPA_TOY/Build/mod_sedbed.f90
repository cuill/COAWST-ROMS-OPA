      MODULE mod_sedbed
!
!svn $Id: sedbed_mod.h 830 2017-01-24 21:21:11Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2018 The ROMS/TOMS Group        John C. Warner   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Sediment Model Kernel Variables:                                    !
!                                                                      !
!  bed            Sediment properties in each bed layer:               !
!                   bed(:,:,:,ithck) => layer thickness                !
!                   bed(:,:,:,iaged) => layer age                      !
!                   bed(:,:,:,iporo) => layer porosity                 !
!                   bed(:,:,:,idiff) => layer bio-diffusivity          !
!  bed_frac       Sediment fraction of each size class in each bed     !
!                   layer(nondimensional: 0-1.0).  Sum of              !
!                   bed_frac = 1.0.                                    !
!  bed_mass       Sediment mass of each size class in each bed layer   !
!                   (kg/m2).
!  bottom         Exposed sediment layer properties:                   !
!                   bottom(:,:,isd50) => mean grain diameter           !
!                   bottom(:,:,idens) => mean grain density            !
!                   bottom(:,:,iwsed) => mean settling velocity        !
!                   bottom(:,:,itauc) => mean critical erosion stress  !
!                   bottom(:,:,irlen) => ripple length                 !
!                   bottom(:,:,irhgt) => ripple height                 !
!                   bottom(:,:,ibwav) => bed wave excursion amplitude  !
!                   bottom(:,:,izNik) => Nikuradse bottom roughness    !
!                   bottom(:,:,izbio) => biological bottom roughness   !
!                   bottom(:,:,izbfm) => bed form bottom roughness     !
!                   bottom(:,:,izbld) => bed load bottom roughness     !
!                   bottom(:,:,izapp) => apparent bottom roughness     !
!                   bottom(:,:,izwbl) => wave bottom roughness         !
!                   bottom(:,:,izdef) => default bottom roughness      !
!                   bottom(:,:,iactv) => active layer thickness        !
!                   bottom(:,:,ishgt) => saltation height              !
!  ero_flux       Flux from erosion.                                   !
!  settling_flux  Flux from settling.                                  !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      TYPE T_SEDBED
!
!  Nonlinear model state.
!
        real(r8), pointer :: bed(:,:,:,:)
        real(r8), pointer :: bed_frac(:,:,:,:)
        real(r8), pointer :: bed_mass(:,:,:,:,:)
        real(r8), pointer :: bottom(:,:,:)
        real(r8), pointer :: ero_flux(:,:,:)
        real(r8), pointer :: settling_flux(:,:,:)
      END TYPE T_SEDBED
      TYPE (T_SEDBED), allocatable :: SEDBED(:)
      CONTAINS
      SUBROUTINE allocate_sedbed (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_sediment
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate structure variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( SEDBED(Ngrids) )
!
!  Nonlinear model state.
!
      allocate ( SEDBED(ng) % bed(LBi:UBi,LBj:UBj,Nbed,MBEDP) )
      allocate ( SEDBED(ng) % bed_frac(LBi:UBi,LBj:UBj,Nbed,NST) )
      allocate ( SEDBED(ng) % bed_mass(LBi:UBi,LBj:UBj,Nbed,2,NST) )
      allocate ( SEDBED(ng) % bottom(LBi:UBi,LBj:UBj,MBOTP) )
      allocate ( SEDBED(ng) % ero_flux(LBi:UBi,LBj:UBj,NST) )
      allocate ( SEDBED(ng) % settling_flux(LBi:UBi,LBj:UBj,NST) )
      RETURN
      END SUBROUTINE allocate_sedbed
      SUBROUTINE initialize_sedbed (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize structure variables in the module using     !
!  first touch distribution policy. In shared-memory configuration,    !
!  this operation actually performs the propagation of the  shared     !
!  arrays  across the cluster,  unless another policy is specified     !
!  to  override the default.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_sediment
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, itrc, j, k
      real(r8), parameter :: IniVal = 0.0_r8
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
!  Set array initialization range.
!
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
!
!-----------------------------------------------------------------------
!  Initialize sediment structure variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO itrc=1,MBEDP
            DO k=1,Nbed
              DO i=Imin,Imax
                SEDBED(ng) % bed(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
          DO itrc=1,NST
            DO k=1,Nbed
              DO i=Imin,Imax
                SEDBED(ng) % bed_frac(i,j,k,itrc) = IniVal
                SEDBED(ng) % bed_mass(i,j,k,1,itrc) = IniVal
                SEDBED(ng) % bed_mass(i,j,k,2,itrc) = IniVal
              END DO
            END DO
          END DO
          DO itrc=1,MBOTP
            DO i=Imin,Imax
              SEDBED(ng) % bottom(i,j,itrc) = IniVal
            END DO
          END DO
          DO itrc=1,NST
            DO i=Imin,Imax
              SEDBED(ng) % ero_flux(i,j,itrc) = IniVal
              SEDBED(ng) % settling_flux(i,j,itrc) = IniVal
            END DO
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE initialize_sedbed
      END MODULE mod_sedbed
