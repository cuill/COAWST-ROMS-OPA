      MODULE mod_oil_Eulvar
! Subroutines for allocating and initializing
! 3D arrays for keeping oil variables
! in Eularian coordinates for
! providing input fields to the Biological and sediments modules
!
!=======================================================================
!                                                                      !
!=======================================================================
!                                                                      !
! Allocate and initialize arrays                                       !
!                                                                      !
!  COil - Oil Content/concentr in a grid cell, kg/m3 3D array          !
!         Coil(:,:,:,1:Nocmp) => Oil content/conc by oil components    !
!         Coil = Mass oil/Vgrid by oil components                      !
!                                                                      !
!  Doil - Average Oil droplet size in 3D grid cells                    !
!                                                                      !
!  NFLT3D - Number of oil Lagr. floats in 3D grid cells                !
!                                                                      !
!                                                                      !
!=======================================================================
      USE mod_kinds
      USE mod_floats
!
      implicit none
!
      TYPE T_OIL3D
!
! Oil Conc in grid cells, 3D array by components
! Mean droplet size
!
      real(r8), pointer :: Coil(:,:,:,:)
      real(r8), pointer :: Doil(:,:,:)
      integer, pointer :: NFLT3D(:,:,:)
!
      END TYPE T_OIL3D
      TYPE (T_OIL3D), allocatable :: OIL3D(:)
      CONTAINS
!
      SUBROUTINE allocate_oil_Eulvar (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates 3D arrays for Oil variables                  !
!  mapped onto  Eularian coordinates from Lagrangian                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_parallel  !! 
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
      integer,parameter ::Nocmp = 3
!-----------------------------------------------------------------------
!  Allocate structure variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( OIL3D(Ngrids) )
!
      allocate ( OIL3D(ng) % Coil(LBi:UBi,LBj:UBj,N(ng),Nocmp) )
      allocate ( OIL3D(ng) % Doil(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OIL3D(ng) % NFLT3D(LBi:UBi,LBj:UBj,N(ng)) ) 
!
! Debug
      IF (MyRank.eq.MyMaster) THEN
        print*,'Allocating: My Rank =', MyRank
        print*,'  LBi=',LBi,' UBi=',UBi,' LBj=',LBj,' UBj=',UBj
        print*,'Size I NFLT3D=',size( OIL3d(ng) % NFLT3D,1 )
        print*,'Size J NFLT3D=',size( OIL3d(ng) % NFLT3D,2 )
        print*,'Size K NFLT3D=',size( OIL3d(ng) % NFLT3D,3 )
      ENDIF
! end debug
!
      RETURN
      END SUBROUTINE allocate_oil_Eulvar
      SUBROUTINE initialize_oil_Eulvar (ng, tile, model)
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
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
      integer,parameter ::Nocmp = 3
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
!  Initialize oil structure variables, only nonlinear model state
!-----------------------------------------------------------------------
!
!
!      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO itrc=1,Nocmp
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                OIL3D(ng) % Coil(i,j,k,itrc) = IniVal
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              OIL3D(ng) % Doil(i,j,k) = IniVal
              OIL3D(ng) % NFLT3D(i,j,k) = -1
            ENDDO
          ENDDO
        ENDDO
!     ENDIF
      RETURN
      END SUBROUTINE initialize_oil_Eulvar
      END MODULE mod_oil_Eulvar
