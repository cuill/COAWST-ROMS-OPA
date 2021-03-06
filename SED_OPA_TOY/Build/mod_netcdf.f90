      MODULE mod_netcdf
!
!svn $Id: mod_netcdf.F 873 2017-10-05 20:27:10Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2018 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains all NetCDF variables definitions. It also      !
!  contains several variables and routines to facilitate  generic      !
!  manipulations of NetCDF data.  Usually, several NetCDF library      !
!  calls are required to inquire and read a dimension or variable.     !
!  These routines provide a single interface for such operations.      !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE netcdf
!
      implicit none
!
      INTERFACE netcdf_get_fvar
        MODULE PROCEDURE netcdf_get_fvar_0d
        MODULE PROCEDURE netcdf_get_fvar_1d
        MODULE PROCEDURE netcdf_get_fvar_2d
        MODULE PROCEDURE netcdf_get_fvar_3d
        MODULE PROCEDURE netcdf_get_fvar_4d
      END INTERFACE netcdf_get_fvar
      INTERFACE netcdf_get_lvar
        MODULE PROCEDURE netcdf_get_lvar_0d
        MODULE PROCEDURE netcdf_get_lvar_1d
      END INTERFACE netcdf_get_lvar
      INTERFACE netcdf_get_ivar
        MODULE PROCEDURE netcdf_get_ivar_0d
        MODULE PROCEDURE netcdf_get_ivar_1d
        MODULE PROCEDURE netcdf_get_ivar_2d
      END INTERFACE netcdf_get_ivar
      INTERFACE netcdf_get_svar
        MODULE PROCEDURE netcdf_get_svar_0d
        MODULE PROCEDURE netcdf_get_svar_1d
        MODULE PROCEDURE netcdf_get_svar_2d
        MODULE PROCEDURE netcdf_get_svar_3d
      END INTERFACE netcdf_get_svar
      INTERFACE netcdf_put_fvar
        MODULE PROCEDURE netcdf_put_fvar_0d
        MODULE PROCEDURE netcdf_put_fvar_1d
        MODULE PROCEDURE netcdf_put_fvar_2d
        MODULE PROCEDURE netcdf_put_fvar_3d
        MODULE PROCEDURE netcdf_put_fvar_4d
      END INTERFACE netcdf_put_fvar
      INTERFACE netcdf_put_ivar
        MODULE PROCEDURE netcdf_put_ivar_0d
        MODULE PROCEDURE netcdf_put_ivar_1d
        MODULE PROCEDURE netcdf_put_ivar_2d
      END INTERFACE netcdf_put_ivar
      INTERFACE netcdf_put_lvar
        MODULE PROCEDURE netcdf_put_lvar_0d
        MODULE PROCEDURE netcdf_put_lvar_1d
        MODULE PROCEDURE netcdf_put_lvar_2d
      END INTERFACE netcdf_put_lvar
      INTERFACE netcdf_put_svar
        MODULE PROCEDURE netcdf_put_svar_0d
        MODULE PROCEDURE netcdf_put_svar_1d
        MODULE PROCEDURE netcdf_put_svar_2d
        MODULE PROCEDURE netcdf_put_svar_3d
      END INTERFACE netcdf_put_svar
!
!  Local dimension parameters.
!
      integer, parameter :: Mdims = 50  ! maximum number of dimensions
      integer, parameter :: Mvars = 900 ! maximum number of variables
      integer, parameter :: NvarD = 5   ! number of variable dimensions
      integer, parameter :: NvarA = 50  ! number of variable attributes
!
!  Generic information about current NetCDF for all dimensions and
!  all variables.
!
      integer :: n_dim                  ! number of dimensions
      integer :: n_var                  ! number of variables
      integer :: n_gatt                 ! number of global attributes
      integer :: ncformat               ! file format number (version)
      integer :: rec_id                 ! unlimited dimension ID
      integer :: rec_size               ! unlimited dimension value
      integer :: dim_id(Mdims)          ! dimensions ID
      integer :: dim_size(Mdims)        ! dimensions value
      integer :: var_id(Mvars)          ! variables ID
      integer :: var_natt(Mvars)        ! variables number of attributes
      integer :: var_flag(Mvars)        ! Variables water points flag
      integer :: var_type(Mvars)        ! variables external data type
      integer :: var_ndim(Mvars)        ! variables number of dimensions
      integer :: var_dim(NvarD,Mvars)   ! variables dimensions ID
!
      character (len=40) :: dim_name(Mdims)      ! dimensions name
      character (len=40) :: var_name(Mvars)      ! variables name
!
!  Generic information about requested current variable.
!
      integer :: n_vdim                 ! number of variable dimensions
      integer :: n_vatt                 ! number of variable attributes
      integer :: var_kind               ! external data type
      integer :: var_Dids(NvarD)        ! dimensions ID
      integer :: var_Dsize(NvarD)       ! dimensions values
      integer :: var_Aint(NvarA)        ! attribute integer values
      real(r8) :: var_Afloat(NvarA)     ! attribute float values
!
      character (len=40)   :: var_Aname(NvarA)   ! Attribute names
      character (len=40)   :: var_Dname(NvarD)   ! dimension names
      character (len=1024) :: var_Achar(NvarA)   ! Attribute char values
!
!  External data representation for floating-point variables.
!
      integer, parameter :: NF_FOUT = nf90_double
      integer, parameter :: NF_FRST = nf90_double
      integer, parameter :: NF_TYPE = nf90_double
!
!  Netcdf file basic creation mode flag. Its value is further modified
!  in "netcdf_create" using the IOR function to include additional
!  flags, like "nf90_clobber" to overwrite existing file.
!
      integer :: CMODE = nf90_64bit_offset ! NetCDF 64-bit offset format
!                                            (large file support)
!
!  Set NetCDF-4/HDF5 deflate (file compression) parameters.
!
      integer :: shuffle = 1
      integer :: deflate = 1
      integer :: deflate_level = 1
!
      CONTAINS
!
      SUBROUTINE netcdf_get_dim (ng, model, ncname, ncid, DimName,      &
     &                           DimSize, DimID)
!
!=======================================================================
!                                                                      !
!  This routine inquires a NetCDF file dimensions names and values.    !
!  All the dimension information is stored in the module variables     !
!  declared above.  In addition, if a particular dimension name is     !
!  provided, this routine returns the requested information in the     !
!  optional arguments.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     DimName      Requested dimension name (string, OPTIONAL)         !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     DimSize      Size of requested dimension (integer, OPTIONAL)     !
!     DimID        ID of requested dimension (integer, OPTIONAL)       !
!                                                                      !
!  Other information stored in this module:                            !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     rec_size     Size of unlimited dimension                         !
!     dim_name     Dimensions name (1:n_dim)                           !
!     dim_id       Dimensions ID (1:n_dim)                             !
!     dim_size     Dimensions value (1:n_dim)                          !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in), optional :: DimName
      integer, intent(out), optional :: DimSize
      integer, intent(out), optional :: DimID
!
!  Local variable declarations.
!
      logical :: foundit
      integer :: my_ncid, i, j, status
      integer :: myID, myValue
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF dimensions (names and values).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      n_dim=0
      n_var=0
      n_gatt=0
      rec_id=-1
      rec_size=0
      dim_id=0
      dim_size=0
      DO i=1,Mdims
        DO j=1,LEN(dim_name(1))
          dim_name(i)(j:j)=' '
        END DO
      END DO
!
!  Open file for reading.
!
      IF (.not.PRESENT(ncid)) THEN
         CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
         IF (FoundError(exit_flag, NoError, 280,                        &
     &                  "ROMS/Modules/mod_netcdf.F" //                                     &
     &                  ", netcdf_get_dim")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Inquire file.
!
      IF (InpThread) THEN
        status=nf90_inquire(my_ncid, n_dim, n_var, n_gatt, rec_id)
        ncformat=-1
        IF ((status.eq.nf90_noerr).and.(n_dim.le.Mdims)) THEN
!
!  Inquire about dimensions: names, ID, and size.
!
          rec_size=-1
          DO i=1,n_dim
            dim_id(i)=i
            status=nf90_inquire_dimension(my_ncid, dim_id(i),           &
     &                                    dim_name(i), dim_size(i))
            IF (FoundError(status, nf90_noerr, 312,                     &
     &                     "ROMS/Modules/mod_netcdf.F" //                                  &
     &                     ", netcdf_get_dim")) THEN
              WRITE (stdout,10) dim_id(i), TRIM(ncname),                &
     &                          TRIM(SourceFile)
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            IF (dim_id(i).eq.rec_id) THEN
              rec_size=dim_size(i)
            END IF
          END DO
        ELSE
          IF (n_dim.gt.Mdims) THEN
            WRITE (stdout,20) ' Mdims = ', Mdims, n_dim
            exit_flag=2
            ioerror=0
          END IF
          IF (FoundError(status, nf90_noerr, 334,                       &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_dim")) THEN
            WRITE (stdout,30) TRIM(ncname), TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF
!
!  Load requested information.
!
      IF (exit_flag.eq.NoError) THEN
        foundit=.FALSE.
        IF (PRESENT(DimName)) THEN
          DO i=1,n_dim
            IF (TRIM(dim_name(i)).eq.TRIM(DimName)) THEN
              foundit=.TRUE.
              myID=dim_id(i)
              myValue=dim_size(i)
            END IF
          END DO
          IF (foundit) THEN
            IF (PRESENT(DimSize)) THEN
              DimSize=myValue
            END IF
            IF (PRESENT(DimID)) THEN
              DimID=myID
            END IF
          ELSE
            WRITE (stdout,40) TRIM(DimName), TRIM(ncname)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_DIM - error while reading dimension ID:',  &
     &        2x,i3,/,18x,'in input file:',2x,a,/,18x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_DIM - too small dimension parameter,',a,   &
     &        2i5,/,18x,'change file  mod_netcdf.F  and recompile')
  30  FORMAT (/,' NETCDF_GET_DIM - unable to inquire about contents',   &
     &          ' of input NetCDF file:',2x,a,/,18x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_GET_DIM - requested dimension: ',a,/18x,       &
     &        'not found in input file:',2x,a,/,18x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_dim
      SUBROUTINE netcdf_check_dim (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine inquires a NetCDF file dimensions names and values.    !
!  It checks the file dimensions against model dimension parameters    !
!  for consitency.  All the dimensions information is stored in the    !
!  module variables declared above.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!                                                                      !
!  On output the following information is stored in this module:       !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     rec_size     Size of unlimited dimension                         !
!     dim_name     Dimensions name (1:n_dim)                           !
!     dim_id       Dimensions ID (1:n_dim)                             !
!     dim_size     Dimensions value (1:n_dim)                          !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: i, status
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF dimensions (names and values).
!-----------------------------------------------------------------------
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_get_dim (ng, model, ncname)
      ELSE
        CALL netcdf_get_dim (ng, model, ncname,                         &
     &                       ncid = ncid)
      END IF
      IF (FoundError(exit_flag, NoError, 466,                           &
     &               "ROMS/Modules/mod_netcdf.F" //                                        &
     &               ", netcdf_check_dim")) RETURN
!
!-----------------------------------------------------------------------
!  Check dimensions for consistency.
!-----------------------------------------------------------------------
!
      DO i=1,n_dim
        SELECT CASE (TRIM(ADJUSTL(dim_name(i))))
          CASE ('xi_psi')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_psi) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_psi,        &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_psi')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_psi) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_psi,       &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('xi_rho')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_rho) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_rho,        &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_rho')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_rho) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_rho,       &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('xi_u')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_u) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_u,          &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_u')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_u) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_u,         &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('xi_v')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_v) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_v,          &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_v')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_v) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_v,         &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('IorJ')
            IF (dim_size(i).ne.IOBOUNDS(ng)%IorJ) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%IorJ,          &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('s_rho')
            IF (dim_size(i).ne.N(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i), N(ng),         &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
        END SELECT
      END DO
 10   FORMAT (/,' NETCDF_CHECK_DIM - inconsistent size of dimension: ', &
     &        a,2x,2i5,/,20x,'in file: ',a)
      RETURN
      END SUBROUTINE netcdf_check_dim
      SUBROUTINE netcdf_check_var (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine inquires the NetCDF file variables and check if the    !
!  values of few of them are consitent with the parameters provided    !
!  in input scripts. All the variables information is stored in the    !
!  module variables declared above.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!                                                                      !
!  On output the following information is stored in this module:       !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     var_name     Variables name (1:n_var)                            !
!     var_id       Variables ID (1:n_var)                              !
!     var_natt     Variables number of attributes (1:n_var)            !
!     var_flag     Variables flag [1=full field, -1=water points only] !
!     var_type     Variables external data type (1:n_var)              !
!     var_ndim     Variables number of dimensions (1:n_var)            !
!     var_dim      Variables dimensions ID (:,1:n_var)                 !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: IDmod, Npts, i, ib, ic, j, j1, j2, status
      integer :: IvarS
      real(r8), parameter :: RoundOff = 1.0e-7_r8
      real(r8) :: FvarS, FvarV(50), VarVal
      character (len=80) :: text
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF variables.
!-----------------------------------------------------------------------
!
!  Limit model identifier. The profiling is limited to iNLM, iTLM, iRPM,
!  and iADM.
!
      IF ((model.lt.1).or.(model.gt.4)) THEN
        IDmod=iNLM
      ELSE
        IDmod=model
      END IF
!
! Inquire about all variables.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_inq_var (ng, IDmod, ncname)
      ELSE
        CALL netcdf_inq_var (ng, IDmod, ncname,                         &
     &                       ncid = ncid)
      END IF
      IF (FoundError(exit_flag, NoError, 680,                           &
     &               "ROMS/Modules/mod_netcdf.F" //                                        &
     &               ", netcdf_check_var")) RETURN
!
!-----------------------------------------------------------------------
!  Check several important variables for consistency.
!-----------------------------------------------------------------------
!
      DO i=1,n_var
        SELECT CASE (TRIM(ADJUSTL(var_name(i))))
          CASE ('Vtransform')
            IF (.not.PRESENT(ncid)) THEN
              CALL netcdf_get_ivar (ng, IDmod, ncname, var_name(i),     &
     &                              IvarS)
            ELSE
              CALL netcdf_get_ivar (ng, IDmod, ncname, var_name(i),     &
     &                              IvarS,                              &
     &                              ncid = ncid)
            END IF
            IF (FoundError(exit_flag, NoError, 701,                     &
     &                     "ROMS/Modules/mod_netcdf.F" //                                  &
     &                     ", netcdf_check_var")) RETURN
            IF (IvarS.ne.Vtransform(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(var_name(i)),          &
     &                                      IvarS, Vtransform(ng),      &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('Vstretching')
            IF (.not.PRESENT(ncid)) THEN
              CALL netcdf_get_ivar (ng, IDmod, ncname, var_name(i),     &
     &                              IvarS)
            ELSE
              CALL netcdf_get_ivar (ng, IDmod, ncname, var_name(i),     &
     &                              IvarS,                              &
     &                              ncid = ncid)
            END IF
            IF (FoundError(exit_flag, NoError, 722,                     &
     &                     "ROMS/Modules/mod_netcdf.F" //                                  &
     &                     ", netcdf_check_var")) RETURN
            IF (IvarS.ne.Vstretching(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(var_name(i)),          &
     &                                      IvarS, Vstretching(ng),     &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('hc')
            IF (.not.PRESENT(ncid)) THEN
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS)
            ELSE
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS,                              &
     &                              ncid = ncid)
            END IF
            IF (FoundError(exit_flag, NoError, 743,                     &
     &                     "ROMS/Modules/mod_netcdf.F" //                                  &
     &                     ", netcdf_check_var")) RETURN
            IF (ABS(hc(ng)-FvarS).gt.RoundOff) THEN
              IF (Master) WRITE (stdout,20) TRIM(var_name(i)),          &
     &                                      FvarS, hc(ng),              &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('theta_s')
            IF (.not.PRESENT(ncid)) THEN
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS)
            ELSE
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS,                              &
     &                              ncid = ncid)
            END IF
            IF (FoundError(exit_flag, NoError, 763,                     &
     &                     "ROMS/Modules/mod_netcdf.F" //                                  &
     &                     ", netcdf_check_var")) RETURN
            IF (ABS(theta_s(ng)-FvarS).gt.RoundOff) THEN
              IF (Master) WRITE (stdout,20) TRIM(var_name(i)),          &
     &                                      FvarS, theta_s(ng),         &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('theta_b')
            IF (.not.PRESENT(ncid)) THEN
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS)
            ELSE
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS,                              &
     &                              ncid = ncid)
            END IF
            IF (FoundError(exit_flag, NoError, 783,                     &
     &                     "ROMS/Modules/mod_netcdf.F" //                                  &
     &                     ", netcdf_check_var")) RETURN
            IF (ABS(theta_b(ng)-FvarS).gt.RoundOff) THEN
              IF (Master) WRITE (stdout,20) TRIM(var_name(i)),          &
     &                                      FvarS, theta_b(ng),         &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('Tcline')
            IF (.not.PRESENT(ncid)) THEN
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS)
            ELSE
              CALL netcdf_get_fvar (ng, IDmod, ncname, var_name(i),     &
     &                              FvarS,                              &
     &                              ncid = ncid)
            END IF
            IF (FoundError(exit_flag, NoError, 803,                     &
     &                     "ROMS/Modules/mod_netcdf.F" //                                  &
     &                     ", netcdf_check_var")) RETURN
            IF (ABS(Tcline(ng)-FvarS).gt.RoundOff) THEN
              IF (Master) WRITE (stdout,20) TRIM(var_name(i)),          &
     &                                      FvarS, Tcline(ng),          &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
        END SELECT
      END DO
 10   FORMAT (/,' NETCDF_CHECK_VAR - inconsistent value of variable: ', &
     &        a,2x,2i5,/,20x,'in file: ',a)
 20   FORMAT (/,' NETCDF_CHECK_VAR - inconsistent value of variable: ', &
     &        a,2x,2(1pe14.6),/,20x,'in file: ',a)
      RETURN
      END SUBROUTINE netcdf_check_var
      SUBROUTINE netcdf_inq_var (ng, model, ncname, ncid, myVarName,    &
     &                           SearchVar, VarID, nVarDim, nVarAtt)
!
!=======================================================================
!                                                                      !
!  This routine inquires a NetCDF file dimensions names and values.    !
!  All the dimension information is stored in the module variables     !
!  declared above.  In addition,  if a particular variable name is     !
!  provided, this routine returns the requested information in the     !
!  optional arguments.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     myVarName    Requested variable name (string, OPTIONAL)          !
!     SearchVar    Switch used when searching a variable over          !
!                    multiple NetCDF files (logical, OPTIONAL)         !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     VarID        ID of requested variable (integer, OPTIONAL)        !
!     nVarDim      Number of variable dimensions (integer, OPTIONAL)   !
!     nVarAtt      Number of variable attributes (integer, OPTIONAL)   !
!                                                                      !
!  Other information stored in this module:                            !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     var_name     Variables name (1:n_var)                            !
!     var_id       Variables ID (1:n_var)                              !
!     var_natt     Variables number of attributes (1:n_var)            !
!     var_flag     Variables flag [1=full field, -1=water points only] !
!     var_type     Variables external data type (1:n_var)              !
!     var_ndim     Variables number of dimensions (1:n_var)            !
!     var_dim      Variables dimensions ID (:,1:n_var)                 !
!                                                                      !
!  If the OPTIONAL argument myVarName is provided, the following       !
!  information for requested variable is also stored:                  !
!                                                                      !
!     n_vatt       Number of variable attributes                       !
!     n_vdim       Number of variable dimensions                       !
!     var_kind     Variable external data type                         !
!     var_Aname    Variable attribute names (1:n_vatt)                 !
!     var_Achar    Variable string attribute values (1:n_vatt)         !
!     var_Afloat   Variable float attribute values (1:n_vatt)          !
!     var_Aint     Variable integer attribute values (1:n_vatt)        !
!     var_Dids     Variable dimensions ID (1:n_vdim)                   !
!     var_Dname    Variable dimensions name (1:n_vdim)                 !
!     var_Dsize    Variable dimensions size (1:n_vdim)                 !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in), optional :: myVarName
      logical, intent(out), optional :: SearchVar
      integer, intent(in), optional :: ncid
      integer, intent(out), optional :: VarID
      integer, intent(out), optional :: nVarDim
      integer, intent(out), optional :: nVarAtt
!
!  Local variable declarations.
!
      logical :: foundit, WriteError
      integer :: i, j, status
      integer :: my_Alen, my_Atype, my_id, my_ncid
      character (len=1024) :: text
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF dimensions (names and values).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      n_dim=0
      n_var=0
      n_gatt=0
      rec_id=-1
      var_id=0
      var_natt=0
      var_flag=0
      var_type=0
      var_ndim=0
      var_dim=0
      DO i=1,Mvars
        DO j=1,LEN(var_name(1))
          var_name(i)(j:j)=' '
        END DO
      END DO
!
!  Open file for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 1398,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_inq_var")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Inquire NetCDF file.
!
      IF (InpThread) THEN
        status=nf90_inquire(my_ncid, n_dim, n_var, n_gatt, rec_id)
        IF ((status.eq.nf90_noerr).and.(n_var.le.Mvars)) THEN
!
!  Inquire about variables: name, ID, dimensions, data type, and number
!  of attributes.
!
          DO i=1,n_var
            var_id(i)=i
            var_flag(i)=1
            status=nf90_inquire_variable(my_ncid, var_id(i),            &
     &                                   var_name(i), var_type(i),      &
     &                                   var_ndim(i), var_dim(:,i),     &
     &                                   var_natt(i))
            IF (status.eq.nf90_noerr) THEN
              DO j=1,MIN(NvarA,var_natt(i))
                status=nf90_inq_attname(my_ncid, var_id(i), j,          &
     &                                  var_Aname(j))
                IF (status.eq.nf90_noerr) THEN
                  IF (TRIM(var_Aname(j)).eq.'water_points'.and.         &
     &                (var_ndim(i).gt.0)) THEN
                    var_flag(i)=-1
                  END IF
                ELSE
                  WRITE (stdout,10) j, TRIM(var_name(i)),               &
     &                              TRIM(ncname), TRIM(SourceFile)
                  exit_flag=2
                  ioerror=status
                  EXIT
                END IF
              END DO
            ELSE
              WRITE (stdout,20) var_id(i), TRIM(ncname),                &
     &                          TRIM(SourceFile)
              exit_flag=2
              ioerror=status
              EXIT
            END IF
          END DO
        ELSE
          IF (n_var.gt.Mvars) THEN
            WRITE (stdout,30) 'Mvars = ', Mvars, n_var
            exit_flag=2
          END IF
          IF (FoundError(status, nf90_noerr, 1457,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_inq_var")) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF
!
!  Load requested requested variable information.
!
      IF (exit_flag.eq.NoError) THEN
        foundit=.FALSE.
        IF (PRESENT(myVarName)) THEN
          var_Dids=-1
          var_Dsize=0
          var_Aint=0
          var_Afloat=0.0_r8
          DO i=1,NvarA
            DO j=1,LEN(var_Aname(1))
              var_Aname(i)(j:j)=' '
            END DO
            DO j=1,LEN(var_Achar(1))
              var_Achar(i)(j:j)=' '
            END DO
          END DO
          DO i=1,NvarD
            DO j=1,LEN(var_Dname(1))
              var_Dname(i)(j:j)=' '
            END DO
          END DO
!
          DO i=1,n_var
            IF (TRIM(var_name(i)).eq.TRIM(myVarName)) THEN
              foundit=.TRUE.
              my_id=var_id(i)
              n_vatt=var_natt(i)
              n_vdim=var_ndim(i)
              DO j=1,n_vdim
                var_Dids(j)=var_dim(j,i)
              END DO
              var_kind=var_type(i)
            END IF
          END DO
          IF (foundit) THEN
            IF (PRESENT(VarID)) THEN
              VarID=my_id
            END IF
            IF (PRESENT(nVarDim)) THEN
              nVarDim=n_vdim
            END IF
            IF (PRESENT(nVarAtt)) THEN
              nVarAtt=n_vatt
            END IF
          END IF
!
!  If founded requested variable, inquire about is dimensions and
!  attributes.
!
          IF (foundit.and.InpThread) THEN
            DO i=1,n_vdim
              status=nf90_inquire_dimension(my_ncid, var_Dids(i),       &
     &                                      var_Dname(i), var_Dsize(i))
              IF (FoundError(status, nf90_noerr, 1543,                  &
     &                       "ROMS/Modules/mod_netcdf.F" //                                &
     &                       ", netcdf_inq_var")) THEN
                WRITE (stdout,50) i, TRIM(myVarName), TRIM(ncname),     &
     &                            TRIM(SourceFile)
                PRINT *, trim(nf90_strerror(status))
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END DO
            IF (status.eq.nf90_noerr) THEN
              DO i=1,MIN(NvarA, n_vatt)
                status=nf90_inq_attname(my_ncid, my_id, i, var_Aname(i))
                IF (status.eq.nf90_noerr) THEN
                  status=nf90_inquire_attribute(my_ncid, my_id,         &
     &                                          TRIM(var_Aname(i)),     &
     &                                          xtype = my_Atype,       &
     &                                          len = my_Alen)
                  IF (status.eq.nf90_noerr) THEN
                    IF ((my_Alen.eq.1).and.                             &
     &                  (my_Atype.eq.NF90_INT)) THEN
                      status=nf90_get_att(my_ncid, my_id,               &
     &                                    TRIM(var_Aname(i)),           &
     &                                    var_Aint(i))
                      IF (FoundError(status, nf90_noerr, 1568,          &
     &                               "ROMS/Modules/mod_netcdf.F" //                        &
     &                               ", netcdf_inq_var")) THEN
                        WRITE (stdout,60) 'integer',                    &
     &                                    TRIM(var_Aname(i)),           &
     &                                    TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
                        PRINT *, trim(nf90_strerror(status))
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
                    ELSE IF ((my_Alen.eq.1).and.                        &
     &                       ((my_Atype.eq.NF90_FLOAT ).or.             &
     &                        (my_Atype.eq.NF90_DOUBLE))) THEN
                      status=nf90_get_att(my_ncid, my_id,               &
     &                                    TRIM(var_Aname(i)),           &
     &                                    var_Afloat(i))
                      IF (FoundError(status, nf90_noerr, 1587,          &
     &                               "ROMS/Modules/mod_netcdf.F" //                        &
     &                               ", netcdf_inq_var")) THEN
                        WRITE (stdout,60) 'float',                      &
     &                                    TRIM(var_Aname(i)),           &
     &                                    TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
                        PRINT *, trim(nf90_strerror(status))
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
                    ELSE IF (my_Atype.eq.NF90_CHAR) THEN
                      status=nf90_get_att(my_ncid, my_id,               &
     &                                    TRIM(var_Aname(i)),           &
     &                                    text(1:my_Alen))
                      IF (FoundError(status, nf90_noerr, 1604,          &
     &                               "ROMS/Modules/mod_netcdf.F" //                        &
     &                               ", netcdf_inq_var")) THEN
                        WRITE (stdout,60) 'string',                     &
     &                                    TRIM(var_Aname(i)),           &
     &                                    TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
                        PRINT *, trim(nf90_strerror(status))
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
                      var_Achar(i)=text(1:my_Alen)
                    END IF
                  ELSE
                    WRITE (stdout,10) i, TRIM(myVarName), TRIM(ncname), &
     &                                TRIM(SourceFile)
                    exit_flag=2
                    ioerror=status
                    EXIT
                  END IF
                ELSE
                  WRITE (stdout,70) i, TRIM(myVarName), TRIM(ncname),   &
     &                              TRIM(SourceFile)
                  exit_flag=4
                  ioerror=status
                  EXIT
                END IF
              END DO
            END IF
          END IF
!
!  Ignore error message if requested variable not found when searching
!  over multiple input NetCDF files.
!
          IF (PRESENT(SearchVar)) THEN
            SearchVar=foundit
            WriteError=.FALSE.
          ELSE
            WriteError=.TRUE.
          END IF
          IF (.not.foundit.and.WriteError) THEN
            IF (Master) WRITE (stdout,80) TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring attribute ',   &
     &        i1,' for variable: ',a,/,18x,'in input file:',2x,a,       &
     &        /,18x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring variable ID:', &
     &        2x,i3,/,18x,'in input file:',2x,a,/,18x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_INQ_VAR - too small dimension parameter, ',a,  &
     &        2i5,/,18x,'change file  mod_netcdf.F  and recompile')
  40  FORMAT (/,' NETCDF_INQ_VAR - unable to inquire about contents',   &
     &          ' of input NetCDF file:',2x,a,/,18x,'call from:',2x,a)
  50  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring dimension ',   &
     &        i1,' for variable:',2x,a,/,18x,'in input file:',2x,a,     &
     &        /,18x,'call from:',2x,a)
  60  FORMAT (/,' NETCDF_INQ_VAR - error while reading ',a,             &
     &        'attribute:',1x,a,' for variable ',a,/,18x,               &
     &        'in input file:',2x,a,/,18x,'call from:',2x,a)
  70  FORMAT (/,' NETCDF_INQ_VAR - unable to inquire name of ',         &
     &        'attribute ',i1,' for variable ',a,/,18x,                 &
     &        'in input file:',2x,a,/,18x,'call from:',2x,a)
  80  FORMAT (/,' NETCDF_INQ_VAR - requested variable:',2x,a,/18x,      &
     &        'not found in input file:',2x,a,/,18x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_inq_var
      SUBROUTINE netcdf_inq_varid (ng, model, ncname, myVarName,        &
     &                             ncid, VarID)
!
!=======================================================================
!                                                                      !
!  This routine inquires ID of requested NetCDF variable.              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Requested variable name (string)                    !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     VarID        ID of requested variable (integer)                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      integer, intent(out) :: VarID
!
!  Local variable declarations.
!
      integer :: status
!
!-----------------------------------------------------------------------
!  Inquire ID of requested variable.
!-----------------------------------------------------------------------
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncid, TRIM(myVarName), VarID)
        IF (FoundError(status, nf90_noerr, 1761,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_inq_varid")) THEN
          WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          PRINT *, trim(nf90_strerror(status))
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
  10  FORMAT (/,' NETCDF_INQ_VARID - error while inquiring ID for ',    &
     &        'variable:',2x,a,/,20x,'in input file:',2x,a,/,           &
     &        20x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_inq_varid
      SUBROUTINE netcdf_get_fatt (ng, model, ncname, varid, AttName,    &
     &                            AttValue, foundit, ncid)
!
!=======================================================================
!                                                                      !
!  This routine gets requested variable floating-point attribute(s).   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     varid        Variable ID for variable attribute or               !
!                    NF90_GLOBAL for a global attribute (integer)      !
!     AttName      Attribute name to read (string array)               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     AttValue     Attribute value (real array)                        !
!     foundit      Switch (T/F) activated when the requested           !
!                    attribute is found (logical array)                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, varid
      integer, intent(in), optional :: ncid
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: AttName(:)
      logical, intent(out) :: foundit(:)
      real(r8), intent(out) :: AttValue(:)
!
!  Local variable declarations.
!
      integer :: i, j, my_natts, my_ncid, nvatts, status
      character (len=40) :: my_Aname
      character (len=40) :: my_Vname
!
!-----------------------------------------------------------------------
!  Inquire ID of requested variable.
!-----------------------------------------------------------------------
!
!  If appropriate, open file for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 1856,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_fatt")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Get number of variable attributes to process.
!
      nvatts=UBOUND(AttValue, DIM=1)
!
!  Inquire about requested attribute value.
!
      DO i=1,nvatts
        foundit(i)=.FALSE.
        AttValue(i)=0.0_r8
      END DO
      IF (InpThread) THEN
        IF (varid.eq.nf90_global) THEN
          status=nf90_inquire(my_ncid,                                  &
     &                        nAttributes = my_natts)
        ELSE
          status=nf90_inquire_variable(my_ncid, varid,                  &
     &                                 name = my_Vname,                 &
     &                                 nAtts = my_natts)
        END IF
        IF (status.eq.nf90_noerr) THEN
          DO j=1,my_natts
            status=nf90_inq_attname(my_ncid, varid, j, my_Aname)
            IF (status.eq.nf90_noerr) THEN
              DO i=1,nvatts
                IF (TRIM(my_Aname).eq.TRIM(AttName(i))) THEN
                  status=nf90_get_att(my_ncid, varid, TRIM(AttName(i)), &
     &                                AttValue(i))
                  IF (FoundError(status, nf90_noerr, 1891,              &
     &                           "ROMS/Modules/mod_netcdf.F" //                            &
     &                           ", netcdf_get_fatt")) THEN
                    IF (Master) WRITE (stdout,10) TRIM(AttName(i)),     &
     &                                            TRIM(my_Vname),       &
     &                                            TRIM(ncname),         &
     &                                            TRIM(SourceFile)
                    PRINT *, trim(nf90_strerror(status))
                    exit_flag=2
                    ioerror=status
                  END IF
                  foundit(i)=.TRUE.
                  EXIT
                END IF
              END DO
            ELSE
              IF (Master) WRITE (stdout,20) j,                          &
     &                                      TRIM(my_Vname),             &
     &                                      TRIM(ncname),               &
     &                                      TRIM(SourceFile)
              exit_flag=2
              ioerror=status
              EXIT
            END IF
          END DO
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(my_Vname),                 &
     &                                  TRIM(ncname),                   &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  If applicable, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FATT - error while reading attribute:',    &
     &        1x,a,'for variable',1x,a,/,19x,'in input file:',2x,a,/,   &
     &        19x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_FATT - error while inquiring attribute:',  &
     &        1x,i2.2,'for variable',1x,a,/,19x,'in input file:',2x,a,  &
     &        /,19x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_FATT - error while inquiring number of ',  &
     &        'attributes for variable:',1x,a,/,19x,'in input file:',   &
     &        2x,a,/,19x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_fatt
      SUBROUTINE netcdf_get_satt (ng, model, ncname, varid, AttName,    &
     &                            AttValue, foundit, ncid)
!
!=======================================================================
!                                                                      !
!  This routine gets requested variable string attribute(s).           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     varid        Variable ID for variable attribute or               !
!                    NF90_GLOBAL for a global attribute (integer)      !
!     AttName      Attribute name to read (string array)               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     AttValue     Attribute value (string array)                      !
!     foundit      Switch (T/F) activated when the requested           !
!                    attribute is found (logical array)                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, varid
      integer, intent(in), optional :: ncid
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: AttName(:)
      logical, intent(out) :: foundit(:)
      character (len=*), intent(out) :: AttValue(:)
!
!  Local variable declarations.
!
      integer :: i, j, my_natts, my_ncid, nvatts, status
      character (len=40) :: my_Aname
      character (len=40) :: my_Vname
!
!-----------------------------------------------------------------------
!  Inquire ID of requested variable.
!-----------------------------------------------------------------------
!
!  If appropriate, open file for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 2040,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_satt")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Get number of variable attributes to process.
!
      nvatts=UBOUND(AttValue, DIM=1)
!
!  Inquire about requested attribute value.
!
      DO i=1,nvatts
        foundit(i)=.FALSE.
        AttValue(i)=' '
      END DO
      IF (InpThread) THEN
        IF (varid.eq.nf90_global) THEN
          status=nf90_inquire(my_ncid,                                  &
     &                        nAttributes = my_natts)
        ELSE
          status=nf90_inquire_variable(my_ncid, varid,                  &
     &                                 name = my_Vname,                 &
     &                                 nAtts = my_natts)
        END IF
        IF (status.eq.nf90_noerr) THEN
          DO j=1,my_natts
            status=nf90_inq_attname(my_ncid, varid, j, my_Aname)
            IF (status.eq.nf90_noerr) THEN
              DO i=1,nvatts
                IF (TRIM(my_Aname).eq.TRIM(AttName(i))) THEN
                  status=nf90_get_att(my_ncid, varid, TRIM(AttName(i)), &
     &                                AttValue(i))
                  IF (FoundError(status, nf90_noerr, 2075,              &
     &                           "ROMS/Modules/mod_netcdf.F" //                            &
                                 ", netcdf_get_satt")) THEN
                    IF (Master) WRITE (stdout,10) TRIM(AttName(i)),     &
     &                                            TRIM(my_Vname),       &
     &                                            TRIM(ncname),         &
     &                                            TRIM(SourceFile)
                    exit_flag=2
                    ioerror=status
                  END IF
                  foundit(i)=.TRUE.
                  EXIT
                END IF
              END DO
            ELSE
              IF (Master) WRITE (stdout,20) j,                          &
     &                                      TRIM(my_Vname),             &
     &                                      TRIM(ncname),               &
     &                                      TRIM(SourceFile)
              exit_flag=2
              ioerror=status
              EXIT
            END IF
          END DO
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(my_Vname),                 &
     &                                  TRIM(ncname),                   &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  If applicable, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_SATT - error while reading attribute:',    &
     &        1x,a,'for variable',1x,a,/,19x,'in input file:',2x,a,/,   &
     &        19x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_SATT - error while inquiring attribute:',  &
     &        1x,i2.2,'for variable',1x,a,/,19x,'in input file:',2x,a,  &
     &        /,19x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_SATT - error while inquiring number of ',  &
     &        'attributes for variable:',1x,a,/,19x,'in input file:',   &
     &        2x,a,/,19x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_satt
      SUBROUTINE netcdf_get_fvar_0d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total,                &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point scalar variable from    !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (real)                         !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(1))    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val
      real(r8), intent(out) :: A
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid
      real(r8), dimension(1) :: my_A
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 2233,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_fvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, my_A, start, total)
            A=my_A(1)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, 2251,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_fvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          PRINT *, trim(nf90_strerror(status))
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable. Notice that
!  the same read value is assigned since a scalar variable was
!  processed.
!
      IF (PRESENT(min_val)) THEN
        min_val=A
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=A
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_0D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_fvar_0d
      SUBROUTINE netcdf_get_fvar_1d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total,                &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 1D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(0:))   !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(:,1))  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val
      real(r8), intent(out) :: A(:)
!
!  Local variable declarations.
!
      logical, dimension(3) :: foundit
      integer :: i, my_ncid, status, varid
      integer, dimension(1) :: Asize
      real(r8) :: Afactor, Aoffset, Aspval
      real(r8), parameter :: Aepsilon = 1.0E-8_r8
      real(r8), dimension(3) :: AttValue
      character (len=12), dimension(3) :: AttName
!
!-----------------------------------------------------------------------
!  Read in a floating-point 1D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=1
        DO i=1,SIZE(total)              ! this logic is for the case
          Asize(1)=Asize(1)*total(i)    ! of reading multidimensional
        END DO                          ! data into a compact 1D array
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
      END IF
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 2406,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_fvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, 2423,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_fvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '
      CALL netcdf_get_fatt (ng, model, ncname, varid, AttName,          &
     &                      AttValue, foundit,                          &
     &                      ncid = ncid)
      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=AttValue(1)
        END IF
        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=AttValue(2)
        END IF
        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF
        DO i=1,Asize(1)                       ! zero out missing values
          IF ((foundit(3).and.(ABS(A(i)-Aspval).lt.Aepsilon)).or.       &
     &        (.not.foundit(3).and.(ABS(A(i)).ge.ABS(Aspval)))) THEN
            A(i)=0.0_r8
          END IF
        END DO
        IF (foundit(1)) THEN                  ! scale data
          DO i=1,Asize(1)
            A(i)=Afactor*A(i)
          END DO
        END IF
        IF (foundit(2)) THEN                  ! add data offset
          DO i=1,Asize(1)
            A(i)=A(i)+Aoffset
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_1D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_fvar_1d
      SUBROUTINE netcdf_get_fvar_2d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total,                &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 2D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 2D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(0:,:)) !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(0:,0:))!
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(:,:,1))!
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val
      real(r8), intent(out) :: A(:,:)
!
!  Local variable declarations.
!
      logical, dimension(3) :: foundit
      integer :: i, j, my_ncid, status, varid
      integer, dimension(2) :: Asize
      real(r8) :: Afactor, Aoffset, Aspval
      real(r8), parameter :: Aepsilon = 1.0E-8_r8
      real(r8), dimension(3) :: AttValue
      character (len=12), dimension(3) :: AttName
!
!-----------------------------------------------------------------------
!  Read in a floating-point 2D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=total(1)
        Asize(2)=total(2)
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
        Asize(2)=UBOUND(A, DIM=2)
      END IF
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 2636,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_fvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, 2653,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_fvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '
      CALL netcdf_get_fatt (ng, model, ncname, varid, AttName,          &
     &                      AttValue, foundit,                          &
     &                      ncid = ncid)
      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=AttValue(1)
        END IF
        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=AttValue(2)
        END IF
        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF
        DO j=1,Asize(2)                       ! zero out missing values
          DO i=1,Asize(1)
            IF ((foundit(3).and.(ABS(A(i,j)-Aspval).lt.Aepsilon)).or.   &
     &          (.not.foundit(3).and.(ABS(A(i,j)).ge.ABS(Aspval)))) THEN
              A(i,j)=0.0_r8
            END IF
          END DO
        END DO
        IF (foundit(1)) THEN                  ! scale data
          DO j=1,Asize(2)
            DO i=1,Asize(1)
              A(i,j)=Afactor*A(i,j)
            END DO
          END DO
        END IF
        IF (foundit(2)) THEN                  ! add data offset
          DO j=1,Asize(2)
            DO i=1,Asize(1)
              A(i,j)=A(i,j)+Aoffset
            END DO
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_2D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_fvar_2d
      SUBROUTINE netcdf_get_fvar_3d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total,                &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 3D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 3D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val
      real(r8), intent(out) :: A(:,:,:)
!
!  Local variable declarations.
!
      logical, dimension(3) :: foundit
      integer :: i, j, k, my_ncid, status, varid
      integer, dimension(3) :: Asize
      real(r8) :: Afactor, Aoffset, Aspval
      real(r8), parameter :: Aepsilon = 1.0E-8_r8
      real(r8), dimension(3) :: AttValue
      character (len=12), dimension(3) :: AttName
!
!-----------------------------------------------------------------------
!  Read in a floating-point 2D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=total(1)
        Asize(2)=total(2)
        Asize(3)=total(3)
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
        Asize(2)=UBOUND(A, DIM=2)
        Asize(3)=UBOUND(A, DIM=3)
      END IF
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 2870,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_fvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, 2887,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_fvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '
      CALL netcdf_get_fatt (ng, model, ncname, varid, AttName,          &
     &                      AttValue, foundit,                          &
     &                      ncid = ncid)
      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=AttValue(1)
        END IF
        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=AttValue(2)
        END IF
        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF
        DO k=1,Asize(3)                       ! zero out missing values
          DO j=1,Asize(2)
            DO i=1,Asize(1)
              IF ((foundit(3).and.                                      &
     &             (ABS(A(i,j,k)-Aspval).lt.Aepsilon)).or.              &
     &            (.not.foundit(3).and.                                 &
     &             (ABS(A(i,j,k)).ge.ABS(Aspval)))) THEN
                A(i,j,k)=0.0_r8
              END IF
            END DO
          END DO
        END DO
        IF (foundit(1)) THEN                  ! scale data
          DO k=1,Asize(3)
            DO j=1,Asize(2)
              DO i=1,Asize(1)
                A(i,j,k)=Afactor*A(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (foundit(2)) THEN                  ! add data offset
          DO k=1,Asize(3)
            DO j=1,Asize(2)
              DO i=1,Asize(1)
                A(i,j,k)=A(i,j,k)+Aoffset
              END DO
            END DO
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_3D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_3D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_fvar_3d
      SUBROUTINE netcdf_get_fvar_4d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total,                &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 4D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 4D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val
      real(r8), intent(out) :: A(:,:,:,:)
!
!  Local variable declarations.
!
      logical, dimension(3) :: foundit
      integer :: i, j, k, l, my_ncid, status, varid
      integer, dimension(4) :: Asize
      real(r8) :: Afactor, Aoffset, Aspval
      real(r8), parameter :: Aepsilon = 1.0E-8_r8
      real(r8), dimension(3) :: AttValue
      character (len=12), dimension(3) :: AttName
!
!-----------------------------------------------------------------------
!  Read in a floating-point 2D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=total(1)
        Asize(2)=total(2)
        Asize(3)=total(3)
        Asize(4)=total(4)
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
        Asize(2)=UBOUND(A, DIM=2)
        Asize(3)=UBOUND(A, DIM=3)
        Asize(4)=UBOUND(A, DIM=4)
      END IF
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 3115,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_fvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, 3132,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_fvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '
      CALL netcdf_get_fatt (ng, model, ncname, varid, AttName,          &
     &                      AttValue, foundit,                          &
     &                      ncid = ncid)
      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=AttValue(1)
        END IF
        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=AttValue(2)
        END IF
        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF
        DO l=1,Asize(4)                       ! zero out missing values
          DO k=1,Asize(3)
            DO j=1,Asize(2)
              DO i=1,Asize(1)
                IF ((foundit(3).and.                                    &
     &               (ABS(A(i,j,k,l)-Aspval).lt.Aepsilon)).or.          &
     &              (.not.foundit(3).and.                               &
     &               (ABS(A(i,j,k,l)).ge.ABS(Aspval)))) THEN
                  A(i,j,k,l)=0.0_r8
                END IF
              END DO
            END DO
          END DO
        END DO
        IF (foundit(1)) THEN                  ! scale data
          DO l=1,Asize(4)
            DO k=1,Asize(3)
              DO j=1,Asize(2)
                DO i=1,Asize(1)
                  A(i,j,k,l)=Afactor*A(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (foundit(2)) THEN                  ! add data offset
          DO l=1,Asize(4)
            DO k=1,Asize(3)
              DO j=1,Asize(2)
                DO i=1,Asize(1)
                  A(i,j,k,l)=A(i,j,k,l)+Aoffset
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_4D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_4D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_fvar_4d
      SUBROUTINE netcdf_get_lvar_0d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested logical scalar variable from specified !
!  NetCDF file.  The variable can be stored as an interger (0 or 1) or !
!  as a character ('T' or 'F').   Reading a character variable is very !
!  inefficient in parallel I/O.                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (logical)                      !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_lvar (ng, iNLM, 'file.nc', 'VarName', lvar)       !
!    CALL netcdf_get_lvar (ng, iNLM, 'file.nc', 'VarName', lvar(1))    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      logical, intent(out) :: A
!
!  Local variable declarations.
!
      integer :: my_ncid, my_type, status, varid
      integer :: AI
      integer, dimension(1) :: my_AI
      character (len=1) :: Achar
!
!-----------------------------------------------------------------------
!  Read in an integer scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 3344,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_lvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          status=nf90_inquire_variable(my_ncid, varid,                  &
     &                                 xtype = my_type)
          IF (status.eq.nf90_noerr) THEN
            IF (my_type.eq.nf90_int) THEN
              IF (PRESENT(start).and.PRESENT(total)) THEN
                status=nf90_get_var(my_ncid, varid, my_AI, start, total)
                AI=my_AI(1)
              ELSE
                status=nf90_get_var(my_ncid, varid, AI)
              END IF
              IF (status.eq.nf90_noerr) THEN
                IF (AI.eq.0) THEN
                  A=.FALSE.
                ELSE
                  A=.TRUE.
                END IF
              END IF
            ELSE IF (my_type.eq.nf90_char) THEN
              IF (PRESENT(start).and.PRESENT(total)) THEN
                status=nf90_get_var(my_ncid, varid, Achar, start, total)
              ELSE
                status=nf90_get_var(my_ncid, varid, Achar)
              END IF
              IF (status.eq.nf90_noerr) THEN
                A=.FALSE.
                IF ((Achar.eq.'t').or.(Achar.eq.'T')) THEN
                  A=.TRUE.
                END IF
              END IF
            END IF
            IF (FoundError(status, nf90_noerr, 3386,                    &
     &                     "ROMS/Modules/mod_netcdf.F" //                                  &
     &                     ", netcdf_get_lvar")) THEN
              WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),          &
     &                        TRIM(SourceFile)
              PRINT *, trim(nf90_strerror(status))
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_LVAR_0D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_LVAR_0D - error while inquiring type for ',&
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_LVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_lvar_0d
      SUBROUTINE netcdf_get_lvar_1d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested logical 1D-array variable from         !
!  specified  NetCDF file.  The  variable  can be stored as an         !
!  interger (0 or 1) or  as a character ('T' or 'F').  Reading         !
!  a character variable is very inefficient in parallel I/O.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (logical)                    !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_lvar (ng, iNLM, 'file.nc', 'VarName', lvar)       !
!    CALL netcdf_get_lvar (ng, iNLM, 'file.nc', 'VarName', lvar(0:))   !
!    CALL netcdf_get_lvar (ng, iNLM, 'file.nc', 'VarName', lvar(:,1))  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      logical, intent(out) :: A(:)
!
!  Local variable declarations.
!
      integer :: i, my_ncid, my_type, status, varid
      integer, dimension(SIZE(A,1)) :: AI
      character (len=1), dimension(SIZE(A,1)) :: Achar
!
!-----------------------------------------------------------------------
!  Read in an integer scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 3518,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_lvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          status=nf90_inquire_variable(my_ncid, varid,                  &
     &                                 xtype = my_type)
          IF (status.eq.nf90_noerr) THEN
            IF (my_type.eq.nf90_int) THEN
              IF (PRESENT(start).and.PRESENT(total)) THEN
                status=nf90_get_var(my_ncid, varid, AI, start, total)
              ELSE
                status=nf90_get_var(my_ncid, varid, AI)
              END IF
              IF (status.eq.nf90_noerr) THEN
                DO i=1,SIZE(A,1)
                  IF (AI(i).eq.0) THEN
                    A(i)=.FALSE.
                  ELSE
                    A(i)=.TRUE.
                  END IF
                END DO
              END IF
            ELSE IF (my_type.eq.nf90_char) THEN
              IF (PRESENT(start).and.PRESENT(total)) THEN
                status=nf90_get_var(my_ncid, varid, Achar, start, total)
              ELSE
                status=nf90_get_var(my_ncid, varid, Achar)
              END IF
              IF (status.eq.nf90_noerr) THEN
                DO i=1,SIZE(A,1)
                  A(i)=.FALSE.
                  IF ((Achar(i).eq.'t').or.(Achar(i).eq.'T')) THEN
                    A(i)=.TRUE.
                  END IF
                END DO
              END IF
            END IF
            IF (FoundError(status, nf90_noerr, 3563,                    &
     &                     "ROMS/Modules/mod_netcdf.F" //                                  &
     &                     ", netcdf_get_lvar")) THEN
              WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),          &
     &                        TRIM(SourceFile)
              PRINT *, trim(nf90_strerror(status))
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_LVAR_1D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_LVAR_1D - error while inquiring type for ',&
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_LVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_lvar_1d
      SUBROUTINE netcdf_get_ivar_0d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested integer scalar variable from specified !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (integer)                      !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(1))    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      integer, intent(out) :: A
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid
      integer, dimension(1) :: my_A
!
!-----------------------------------------------------------------------
!  Read in an integer scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 3691,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_ivar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, my_A, start, total)
            A=my_A(1)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, 3709,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_ivar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_IVAR_0D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_IVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_ivar_0d
      SUBROUTINE netcdf_get_ivar_1d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested integer 1D-array variable from         !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (integer)                    !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(0:))   !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(:,1))  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      integer, intent(out) :: A(:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid
!
!-----------------------------------------------------------------------
!  Read in an integer 1D-array variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 3827,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_ivar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, 3844,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_ivar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_IVAR_1D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_IVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_ivar_1d
      SUBROUTINE netcdf_get_ivar_2d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested integer 2D-array variable from         !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 2D-array variable (integer)                    !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(0:,:)) !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(0:,0:))!
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(:,:,1))!
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      integer, intent(out) :: A(:,:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid
!
!-----------------------------------------------------------------------
!  Read in an integer 2D-array variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 3963,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_ivar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, 3980,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_ivar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_IVAR_2D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_IVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_ivar_2d
      SUBROUTINE netcdf_get_svar_0d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string scalar variable from specified  !
!  NetCDF file. The CDL of the scalar variable has one-dimension in    !
!  the NetCDF file for the number of characters:                       !
!                                                                      !
!     char string(Nchars)                                         CDL  !
!                                                                      !
!     character (len=Nchars) :: string                            F90  !
!                                                                      !
!  to read a scalar string use:                                        !
!                                                                      !
!     start = (/1/)                                                    !
!     total = (/Nchars/)                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (string)                       !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_svar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_svar (ng, iNLM, 'file.nc', 'VarName', ivar(1))    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      character (len=*), intent(out) :: A
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid
!
!-----------------------------------------------------------------------
!  Read in a string scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 4107,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_svar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, 4124,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_svar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_SVAR_0D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_SVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_svar_0d
      SUBROUTINE netcdf_get_svar_1d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string 1D-array variable or array      !
!  element from specified NetCDF file. The CDL of the 1D-array         !
!  variable has two-dimensions in the NetCDF file, and the first       !
!  dimension is the number of characters:                              !
!                                                                      !
!     char string(dim1, Nchars)                                   CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1)                      F90  !
!                                                                      !
!  to read a single array element at location (i) use:                 !
!                                                                      !
!     start = (/1, i/)                                                 !
!     total = (/Nchars, 1/)                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable or array element (string)    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      character (len=*), intent(out) :: A(:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid
!
!-----------------------------------------------------------------------
!  Read in a string 1D-array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 4247,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_svar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, 4264,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_svar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_SVAR_1D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_SVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_svar_1d
      SUBROUTINE netcdf_get_svar_2d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string 2D-array variable or array      !
!  element from specified NetCDF file. The CDL of the 2D-array         !
!  variable has three-dimensions in the NetCDF file, and the first     !
!  dimension is the number of characters:                              !
!                                                                      !
!     char string(dim2, dim1, Nchars)                             CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1, dim2)                F90  !
!                                                                      !
!  to read a single array element at location (i,j) use:               !
!                                                                      !
!     start = (/1, i, j/)                                              !
!     total = (/Nchars, 1, 1/)                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (3D vector integer, OPTIONAL)        !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (3D vector integer, OPTIONAL)           !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 2D-array variable or array element (string)    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      character (len=*), intent(out) :: A(:,:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid
!
!-----------------------------------------------------------------------
!  Read in a string 2D-array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 4387,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_svar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, 4404,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_svar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_SVAR_2D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_SVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_svar_2d
      SUBROUTINE netcdf_get_svar_3d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string 3D-array variable or array      !
!  element from specified NetCDF file. The CDL of the 3D-array         !
!  variable has four-dimensions in the NetCDF file, and the first      !
!  dimension is the number of characters:                              !
!                                                                      !
!     char string(dim3, dim2, dim1, Nchars)                       CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1, dim2, dim3)          F90  !
!                                                                      !
!  to write a single array element at location (i,j,k) use:            !
!                                                                      !
!     start = (/1, i, j, k/)                                           !
!     total = (/Nchars, 1, 1, 1/)                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (4D vector      !
!                    integer, OPTIONAL)                                !
!     total        Number of data values to be read along each         !
!                    dimension (3D vector integer, OPTIONAL)           !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 3D-array variable or element (string)          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
      character (len=*), intent(out) :: A(:,:,:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid
!
!-----------------------------------------------------------------------
!  Read in a string 3D-array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 4526,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_get_svar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  Read in variable.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
        IF (status.eq.nf90_noerr) THEN
          IF (PRESENT(start).and.PRESENT(total)) THEN
            status=nf90_get_var(my_ncid, varid, A, start, total)
          ELSE
            status=nf90_get_var(my_ncid, varid, A)
          END IF
          IF (FoundError(status, nf90_noerr, 4543,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_get_svar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_SVAR_3D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_SVAR_3D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_get_svar_3d
      SUBROUTINE netcdf_put_fvar_0d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point scalar variable into a NetCDF  !
!  file.  If the NetCDF ID is not provided, it opens the file, writes  !
!  data, and then closes the file.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      real(r8), intent(in) :: A
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
      real(r8), dimension(1) :: my_A
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 4660,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_fvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 4672,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_fvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
            status=nf90_put_var(my_ncid, my_varid, A)
          ELSE
            my_A(1)=A
            status=nf90_put_var(my_ncid, my_varid, my_A, start, total)
          END IF
          IF (FoundError(status, nf90_noerr, 4694,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_fvar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_0D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_fvar_0d
      SUBROUTINE netcdf_put_fvar_1d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point 1D-array variable into a file. !
!  If the NetCDF ID is not provided,  it opens the file,  writes data, !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      real(r8), intent(in) :: A(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 4803,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_fvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 4815,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_fvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, 4832,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_fvar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_1D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_fvar_1d
      SUBROUTINE netcdf_put_fvar_2d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point 2D-array variable into a file. !
!  If the NetCDF ID is not provided,  it opens the file,  writes data, !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      real(r8), intent(in) :: A(:,:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 4941,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_fvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 4953,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_fvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, 4970,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_fvar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_2D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_fvar_2d
      SUBROUTINE netcdf_put_fvar_3d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point 3D-array variable into a file. !
!  If the NetCDF ID is not provided,  it opens the file,  writes data, !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      real(r8), intent(in) :: A(:,:,:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 5079,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_fvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 5091,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_fvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, 5108,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_fvar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_3D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_3D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_fvar_3d
      SUBROUTINE netcdf_put_fvar_4d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point 4D-array variable into a file. !
!  If the NetCDF ID is not provided,  it opens the file,  writes data, !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      real(r8), intent(in) :: A(:,:,:,:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 5217,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_fvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 5229,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_fvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, 5246,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_fvar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_4D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_4D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_fvar_4d
      SUBROUTINE netcdf_put_ivar_0d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a integer scalar variable into a NetCDF file.   !
!  If the NetCDF ID is not provided, it opens the file, writes data,   !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (integer)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      integer, intent(in) :: A
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
      integer, dimension(1) :: my_A
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 5357,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_ivar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 5369,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_ivar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
            status=nf90_put_var(my_ncid, my_varid, A)
          ELSE
            my_A(1)=A
            status=nf90_put_var(my_ncid, my_varid, my_A, start, total)
          END IF
          IF (FoundError(status, nf90_noerr, 5391,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_ivar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_IVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_IVAR_0D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_ivar_0d
      SUBROUTINE netcdf_put_ivar_1d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a integer 1D-array variable into a file. If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (integer)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      integer, intent(in) :: A(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 5500,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_ivar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 5512,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_ivar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, 5529,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_ivar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_IVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_IVAR_1D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_ivar_1d
      SUBROUTINE netcdf_put_ivar_2d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a integer 2D-array variable into a file. If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (integer)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      integer, intent(in) :: A(:,:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 5638,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_ivar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 5650,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_ivar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, 5667,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_ivar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_IVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_IVAR_2D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_ivar_2d
      SUBROUTINE netcdf_put_lvar_0d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a logical scalar variable into a NetCDF file.   !
!  If the NetCDF ID is not provided, it opens the file, writes data,   !
!  and then closes the file.                                           !
!                                                                      !
!  The input logical data is converted to integer such that .FALSE.    !
!  is interpreted as zero, and .TRUE. is interpreted as one.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (logical)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      logical, intent(in) :: A
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
      integer :: AI
      integer, dimension(1) :: my_AI
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 5782,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_lvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 5794,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_lvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Convert logical data to integer: .FALSE. is interpreted as zero, and
!  .TRUE. is interpreted as one.
!
        IF (A) THEN
          AI=1
        ELSE
          AI=0
        END IF
!
!  Write out logical data as integers.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
            status=nf90_put_var(my_ncid, my_varid, AI)
          ELSE
            my_AI(1)=AI
            status=nf90_put_var(my_ncid, my_varid, my_AI, start, total)
          END IF
          IF (FoundError(status, nf90_noerr, 5825,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_lvar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_LVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_LVAR_0D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_lvar_0d
      SUBROUTINE netcdf_put_lvar_1d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a logical 1D-array variable into a file. If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  The input logical data is converted to integer such that .FALSE.    !
!  is interpreted as zero, and .TRUE. is interpreted as one.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (logical)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      logical, intent(in) :: A(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: i, my_ncid, my_varid, status
      integer, dimension(SIZE(A,1)) :: AI
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 5939,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_lvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 5951,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_lvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Convert logical data to integer: .FALSE. is interpreted as zero, and
!  .TRUE. is interpreted as one.
!
        DO i=1,SIZE(A,1)
          IF (A(i)) THEN
            AI(i)=1
          ELSE
            AI(i)=0
          END IF
        END DO
!
!  Write out logical data as integers.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, AI, start, total)
          IF (FoundError(status, nf90_noerr, 5979,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_lvar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_LVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_LVAR_1D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_lvar_1d
      SUBROUTINE netcdf_put_lvar_2d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a logical 2D-array variable into a file. If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  The input logical data is converted to integer such that .FALSE.    !
!  is interpreted as zero, and .TRUE. is interpreted as one.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (logical)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      logical, intent(in) :: A(:,:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: i, j, my_ncid, my_varid, status
      integer, dimension(SIZE(A,1),SIZE(A,2)) :: AI
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 6093,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_lvar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 6105,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_lvar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Convert logical data to integer: .FALSE. is interpreted as zero, and
!  .TRUE. is interpreted as one.
!
        DO j=1,SIZE(A,2)
          DO i=1,SIZE(A,1)
            IF (A(i,j)) THEN
              AI(i,j)=1
            ELSE
              AI(i,j)=0
            END IF
          END DO
        END DO
!
!  Write out logical data as integers.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, AI, start, total)
          IF (FoundError(status, nf90_noerr, 6135,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_lvar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_LVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_LVAR_2D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_lvar_2d
      SUBROUTINE netcdf_put_svar_0d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a string scalar variable into a file.  If       !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file. The CDL of the scalar variable has        !
!  one-dimension in the NetCDF file for the number of characters:      !
!                                                                      !
!     char string(Nchars)                                         CDL  !
!                                                                      !
!     character (len=Nchars) :: string                            F90  !
!                                                                      !
!  to write a scalar string use:                                       !
!                                                                      !
!     start = (/1/)                                                    !
!     total = (/Nchars/)                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (string)                !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (1D vector   !
!                    integer)                                          !
!     total        Number of data values to be written along each      !
!                    dimension (1D vector integer)                     !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)    ! dummy variables
      integer, intent(in), optional :: ncid, varid
      character (len=*), intent(in) :: A
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
!
!-----------------------------------------------------------------------
!  Write out a scalar string.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 6254,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_svar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 6266,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_svar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, 6283,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_svar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_SVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_SVAR_0D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_svar_0d
      SUBROUTINE netcdf_put_svar_1d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a string 1D-array variable into a file.  If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file. The CDL of the 1D-array variable has      !
!  two-dimensions in the NetCDF file, and the first dimension is       !
!  the number of characters:                                           !
!                                                                      !
!     char string(dim1, Nchars)                                   CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1)                      F90  !
!                                                                      !
!  to write a single array element at location (i) use:                !
!                                                                      !
!     start = (/1, i/)                                                 !
!     total = (/Nchars, 1/)                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (1D string array)       !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (2D vector   !
!                    integer)                                          !
!     total        Number of data values to be written along each      !
!                    dimension (2D vector integer)                     !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      character (len=*), intent(in) :: A(:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
      integer, dimension(2) :: mystart, mytotal
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 6406,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_svar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 6418,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_svar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, 6435,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_svar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_SVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_SVAR_1D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_svar_1d
      SUBROUTINE netcdf_put_svar_2d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a string 2D-array variable into a file.  If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file. The CDL of the 2D-array variable has      !
!  three-dimensions in the NetCDF file, and the first dimension is     !
!  the number of characters:                                           !
!                                                                      !
!     char string(dim2, dim1, Nchars)                             CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1, dim2)                F90  !
!                                                                      !
!  to write a single array element at location (i,j) use:              !
!                                                                      !
!     start = (/1, i, j/)                                              !
!     total = (/Nchars, 1, 1/)                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (2D string array)       !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (3D vector   !
!                    integer)                                          !
!     total        Number of data values to be written along each      !
!                    dimension (3D vector integer)                     !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      character (len=*), intent(in) :: A(:,:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
!
!-----------------------------------------------------------------------
!  Write out a string array element or full 2D array.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 6554,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_svar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 6566,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_svar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, 6582,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_svar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_SVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_SVAR_2D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_svar_2d
      SUBROUTINE netcdf_put_svar_3d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a string 3D-array variable into a file.  If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file. The CDL of the 3D-array variable has      !
!  four-dimensions in the NetCDF file, and the first dimension is      !
!  the number of characters:                                           !
!                                                                      !
!     char string(dim3, dim2, dim1, Nchars)                       CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1, dim2, dim3)          F90  !
!                                                                      !
!  to write a single array element at location (i,j,k) use:            !
!                                                                      !
!     start = (/1, i, j, k/)                                           !
!     total = (/Nchars, 1, 1, 1/)                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (3D string array)       !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (4D vector   !
!                    integer)                                          !
!     total        Number of data values to be written along each      !
!                    dimension (4D vector integer)                     !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
      integer, intent(in), optional :: ncid, varid
      character (len=*), intent(in) :: A(:,:,:)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status
!
!-----------------------------------------------------------------------
!  Write out a string array element or full 3D array.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_open (ng, model, TRIM(ncname), 1, my_ncid)
        IF (FoundError(exit_flag, NoError, 6700,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_put_svar")) RETURN
      ELSE
        my_ncid=ncid
      END IF
!
!  If variable ID is not provided, inquire its value.
!
      IF (OutThread) THEN
        IF (.not.PRESENT(varid)) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
          IF (FoundError(status, nf90_noerr, 6712,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_svar")) THEN
            WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_varid=varid
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (FoundError(status, nf90_noerr, 6728,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_put_svar")) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(ncid)) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' NETCDF_PUT_SVAR_3D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_SVAR_3D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_put_svar_3d
      SUBROUTINE netcdf_close (ng, model, ncid, ncname, Lupdate)
!
!=======================================================================
!                                                                      !
!  This routine closes requested NetCDF file. If appropriate, it       !
!  also performs additional processing, like updating the global       !
!  attributes, before closing the file.                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncid         NetCDF file ID (integer)                            !
!     ncname       NetCDF file name (string, OPTIONAL)                 !
!     Lupdate      Update global attribute (logical, OPTIONAl)         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(inout) :: ncid
      logical, intent(in), optional :: Lupdate
      character (len=*), intent(in), optional :: ncname
!
!  Local variable declarations.
!
      integer :: i, status
      character (len=200) :: my_ncname
!
!-----------------------------------------------------------------------
!  If open, close requested NetCDF file.
!-----------------------------------------------------------------------
!
      IF (OutThread.and.(ncid.ne.-1)) THEN
        DO i=1,LEN(my_ncname)
          my_ncname(i:i)=' '
        END DO
        IF (.not.PRESENT(ncname)) THEN
!
!  Get filename, if any. It will be nice if there is a function in
!  the NetCDF library to do this. Fortunately, the filename is
!  written as a global attribute.
!
          status=nf90_get_att(ncid, nf90_global, 'file', my_ncname)
        ELSE
          my_ncname=TRIM(ncname)
        END IF
!
!  Close requested NetCDF file.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_close(ncid)
          IF (FoundError(status, nf90_noerr, 6896,                      &
     &                   "ROMS/Modules/mod_netcdf.F" //                                    &
     &                   ", netcdf_close")) THEN
            WRITE (stdout,20) ncid, TRIM(my_ncname), TRIM(SourceFile)
            PRINT *, trim(nf90_strerror(status))
            exit_flag=3
            ioerror=status
          END IF
          ncid=-1
        END IF
      END IF
!
  10  FORMAT (/,' NETCDF_CLOSE - error while writing global ',          &
     &        'attribute:',2x,a,/,16x,'file:',2x,a,/,16x,               &
     &        'call from:',2x,a)
  20  FORMAT (/,' NETCDF_CLOSE - error during closing of file, ',       &
     &        'ncid = ',i3,/,16x,'file:',2x,a,/,16x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_close
      SUBROUTINE netcdf_create (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine creates a new NetCDF file. Currently, it only creates  !
!  file for serial or parallel I/O access.                             !
!                                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(out) :: ncid
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: my_cmode, status
!
!-----------------------------------------------------------------------
!  Create requested NetCDF file.
!-----------------------------------------------------------------------
!
      my_cmode=IOR(nf90_clobber, CMODE)
      IF (OutThread) THEN
        status=nf90_create(TRIM(ncname), my_cmode, ncid)
        IF (FoundError(status, nf90_noerr, 7005,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_create")) THEN
          WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
          PRINT *, trim(nf90_strerror(status))
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
  10  FORMAT (/,' NETCDF_CREATE - unable to create output NetCDF ',     &
     &        'file:',/,17x,a,/,17x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_create
      SUBROUTINE netcdf_enddef (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine ends definition in an  open NetCDF  dataset.  The      !
!  changes made in define mode  are checked and committed to disk      !
!  if no errors occurred. The dataset is then placed in data mode,     !
!  so variable data can be read or written.                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!  Note:   There is not need to call this function when processing     !
!          NetCDF-4 format type file.                                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: status
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_enddef(ncid)
        IF (FoundError(status, nf90_noerr, 7092,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_enddef")) THEN
          IF (Master) WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
          PRINT *, trim(nf90_strerror(status))
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
  10  FORMAT (/,' NETCDF_ENDDEF - unable to end definition mode for',   &
     &        ' file:',/,17x,a,/,17x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_enddef
      SUBROUTINE netcdf_open (ng, model, ncname, omode, ncid)
!
!=======================================================================
!                                                                      !
!  This routine opens an existing NetCDF file for access. Currently,   !
!  it only opens file for serial or parallel I/O access.               !
!                                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!     omode        Open mode flag:                                     !
!                    omode = 0,  read-only access (NF90_NOWRITE)       !
!                    omode = 1,  read and write access (NF90_WRITE)    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, omode
      integer, intent(out) :: ncid
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: status
!
!-----------------------------------------------------------------------
!  Create requested NetCDF file.
!-----------------------------------------------------------------------
!
      IF (InpThread) THEN
        SELECT CASE (omode)
          CASE (0)
            status=nf90_open(TRIM(ncname), nf90_nowrite, ncid)
          CASE (1)
            status=nf90_open(TRIM(ncname), nf90_write, ncid)
          CASE DEFAULT
            status=nf90_open(TRIM(ncname), nf90_nowrite, ncid)
        END SELECT
        IF (FoundError(status, nf90_noerr, 7220,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_open")) THEN
          WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
  10  FORMAT (/,' NETCDF_OPEN - unable to open existing NetCDF ',       &
     &        'file:',/,15x,a,/,15x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_open
      SUBROUTINE netcdf_redef (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine puts an open NetCDF dataset into define mode, so       !
!  dimensions, variables, and attributes can be added or renamed       !
!  an attributes can be deleted.                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!  Note:   There is not need to call this function when processing     !
!          NetCDF-4 format type files.                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: status
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_redef(ncid)
        IF (FoundError(status, nf90_noerr, 7300,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_redef")) THEN
          IF (Master) WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
          PRINT *, trim(nf90_strerror(status))
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
  10  FORMAT (/,' NETCDF_REDEF - unable to put in definition mode',     &
     &        ' file:',/,16x,a,/,16x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_redef
      SUBROUTINE netcdf_sync (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine synchronize to disk requested NetCDF file with         !
!  in-memory buffer to make data available to other  processes         !
!  immediately after it is written.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!  NOTE:                                                               !
!                                                                      !
!  Nowadays, it is recommended to have the writer and readers open     !
!  the NetCDF file with the NF90_SHARE flag to improve performance.    !
!  See NetCDF user manual for more details.                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: status
!
!-----------------------------------------------------------------------
!  Create requested NetCDF file.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncid)
        IF (FoundError(status, nf90_noerr, 7377,                        &
     &                 "ROMS/Modules/mod_netcdf.F" //                                      &
     &                 ", netcdf_sync")) THEN
          WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
          PRINT *, trim(nf90_strerror(status))
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
  10  FORMAT (/,' NETCDF_SYNC - unable to synchronize to disk file:',   &
     &        /,15x,a,/,15x,'call from:',2x,a)
      RETURN
      END SUBROUTINE netcdf_sync
      END MODULE mod_netcdf
