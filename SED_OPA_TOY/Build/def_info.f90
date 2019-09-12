      SUBROUTINE def_info (ng, model, ncid, ncname, DimIDs)
!
!svn $Id: def_info.F 855 2017-07-29 03:10:57Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2018 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine defines information variables in requested NetCDF      !
!  file.                                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer).                           !
!     model    Calling model identifier (integer).                     !
!     ncid     NetCDF file ID (integer).                               !
!     ncname   NetCDF file name (character).                           !
!     DimIDs   NetCDF dimensions IDs (integer vector):                 !
!                DimIDs( 1) => XI-dimension at RHO-points.             !
!                DimIDs( 2) => XI-dimension at U-points.               !
!                DimIDs( 3) => XI-dimension at V-points.               !
!                DimIDs( 4) => XI-dimension at PSI-points.             !
!                DimIDs( 5) => ETA-dimension at RHO-points.            !
!                DimIDs( 6) => ETA-dimension at U-points.              !
!                DimIDs( 7) => ETA-dimension at V-points.              !
!                DimIDs( 8) => ETA-dimension at PSI-points.            !
!                DimIDs( 9) => S-dimension at RHO-points.              !
!                DimIDs(10) => S-dimension at W-points.                !
!                DimIDs(11) => Number of tracers dimension.            !
!                DimIDs(12) => Unlimited time record dimension.        !
!                DimIDs(13) => Number of stations dimension.           !
!                DimIDs(14) => Boundary dimension.                     !
!                DimIDs(15) => Number of floats dimension.             !
!                DimIDs(16) => Number sediment bed layers dimension.   !
!                DimIDs(17) => Dimension 2D water RHO-points.          !
!                DimIDs(18) => Dimension 2D water U-points.            !
!                DimIDs(19) => Dimension 2D water V-points.            !
!                DimIDs(20) => Dimension 3D water RHO-points.          !
!                DimIDs(21) => Dimension 3D water U-points.            !
!                DimIDs(23) => Dimension 3D water W-points.            !
!                DimIDs(24) => Dimension sediment bed water points.    !
!                DimIDs(25) => Number of EcoSim phytoplankton groups.  !
!                DimIDs(26) => Number of EcoSim bateria groups.        !
!                DimIDs(27) => Number of EcoSim DOM groups.            !
!                DimIDs(28) => Number of EcoSim fecal groups.          !
!                DimIDs(29) => Number of state variables.              !
!                DimIDs(30) => Number of 3D variables time levels (2). !
!                DimIDs(31) => Number of 2D variables time levels (3). !
!                DimIDs(32) => Number of sediment tracers.             !
!                DimIDs(33) => Number of wave directions.              !
!                DimIDs(35) => Number of vegetation properties.        !
!                DimIDs(36) => Number of oil components                ! DDMITRY
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_strings
!     
      USE def_var_mod, ONLY : def_var
      USE strings_mod, ONLY : FoundError, join_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
      integer, intent(in) :: DimIDs(32)
      character (*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 25
      integer :: brydim, i, ie, is, j, lstr, varid
      integer :: srdim, stadim, status, swdim, trcdim, usrdim
      integer :: seddim
      integer :: p2dgrd(2), tbrydim(2)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: def_dim
      real(r8) :: Aval(6)
      character (len=11 ) :: bryatt
      character (len=11 ) :: frcatt
      character (len=11 ) :: clmatt
      character (len=50 ) :: tiling
      character (len=80 ) :: type
      character (len=120) :: Vinfo(Natt)
      character (len=4096) :: string
!
!-----------------------------------------------------------------------
!  Set dimension variables.
!-----------------------------------------------------------------------
!
      p2dgrd(1)=DimIDs(4)
      p2dgrd(2)=DimIDs(8)
      t2dgrd(1)=DimIDs(1)
      t2dgrd(2)=DimIDs(5)
      u2dgrd(1)=DimIDs(2)
      u2dgrd(2)=DimIDs(6)
      v2dgrd(1)=DimIDs(3)
      v2dgrd(2)=DimIDs(7)
      srdim=DimIDs(9)
      swdim=DimIDs(10)
      trcdim=DimIDs(11)
      seddim=DimIDs(32)
      stadim=DimIDs(13)
      brydim=DimIDs(14)
      tbrydim(1)=DimIDs(11)
      tbrydim(2)=DimIDs(14)
!
!  Set dimension for generic user parameters.
!
      IF ((Nuser.gt.0).and.(ncid.ne.GST(ng)%ncid)) THEN
        status=def_dim(ng, model, ncid, ncname, 'Nuser', 25, usrdim)
        IF (FoundError(exit_flag, NoError, 193,                         &
     &                 "ROMS/Utility/def_info.F")) RETURN
      END IF
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
      DO i=1,6
        Aval(i)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  Define history global attribute.
!
        IF (LEN_TRIM(date_str).gt.0) THEN
          WRITE (history,'(a,1x,a,", ",a)') 'ROMS/TOMS, Version',       &
     &                                      TRIM( version),             &
     &                                      TRIM(date_str)
        ELSE
          WRITE (history,'(a,1x,a)') 'ROMS/TOMS, Version',              &
     &                               TRIM(version)
        END IF
!
!  Set tile decomposition global attribute.
!
        WRITE (tiling,10) NtileI(ng), NtileJ(ng)
!
!  Define file name global attribute.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'file',                &
     &                        TRIM(ncname))
          IF (FoundError(status, nf90_noerr, 234,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define NetCDF format type.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'format',              &
     &                        'netCDF-3 64bit offset file')
          IF (FoundError(status, nf90_noerr, 254,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'format', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define file climate and forecast metadata convention global
!  attribute.
!
        type='CF-1.4, SGRID-0.3'
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'Conventions',         &
     &                        TRIM(type))
          IF (FoundError(status, nf90_noerr, 271,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'Conventions', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define file type global attribute.
!
        IF (ncid.eq.ADM(ng)%ncid) THEN
          type='ROMS/TOMS adjoint history file'
        ELSE IF (ncid.eq.AVG(ng)%ncid) THEN
          type='ROMS/TOMS nonlinear model averages file'
        ELSE IF (ncid.eq.DIA(ng)%ncid) THEN
          type='ROMS/TOMS diagnostics file'
        ELSE IF (ncid.eq.FLT(ng)%ncid) THEN
          type='ROMS/TOMS floats file'
        ELSE IF (ncid.eq.ERR(ng)%ncid) THEN
          type='ROMS/TOMS posterior analysis error covariance matrix'
        ELSE IF (ncid.eq.GST(ng)%ncid) THEN
          type='ROMS/TOMS GST check pointing restart file'
        ELSE IF (ncid.eq.HAR(ng)%ncid) THEN
          type='ROMS/TOMS Least-squared Detiding Harmonics file'
        ELSE IF (ncid.eq.HSS(ng)%ncid) THEN
          type='ROMS/TOMS 4DVAR Hessian eigenvectors file'
        ELSE IF (ncid.eq.HIS(ng)%ncid) THEN
          type='ROMS/TOMS history file'
        ELSE IF (ncid.eq.ITL(ng)%ncid) THEN
          type='ROMS/TOMS tangent linear model initial file'
        ELSE IF (ncid.eq.LCZ(ng)%ncid) THEN
          type='ROMS/TOMS 4DVAR Lanczos vectors file'
        ELSE IF (ncid.eq.LZE(ng)%ncid) THEN
          type='ROMS/TOMS 4DVAR Evolved Lanczos vectors file'
        ELSE IF (ncid.eq.NRM(1,ng)%ncid) THEN
          type='ROMS/TOMS initial conditions error covariance norm file'
        ELSE IF (ncid.eq.NRM(2,ng)%ncid) THEN
          type='ROMS/TOMS model error covariance norm file'
        ELSE IF (ncid.eq.NRM(3,ng)%ncid) THEN
         type='ROMS/TOMS boundary conditions error covariance norm file'
        ELSE IF (ncid.eq.NRM(4,ng)%ncid) THEN
          type='ROMS/TOMS surface forcing error covariance norm file'
        ELSE IF (ncid.eq.QCK(ng)%ncid) THEN
          type='ROMS/TOMS quicksave file'
        ELSE IF (ncid.eq.RST(ng)%ncid) THEN
          type='ROMS/TOMS restart file'
        ELSE IF (ncid.eq.STA(ng)%ncid) THEN
          type='ROMS/TOMS station file'
        ELSE IF (ncid.eq.TLF(ng)%ncid) THEN
          type='ROMS/TOMS tangent linear impulse forcing file'
        ELSE IF (ncid.eq.TLM(ng)%ncid) THEN
          type='ROMS/TOMS tangent linear history file'
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'type',                &
     &                        TRIM(type))
          IF (FoundError(status, nf90_noerr, 335,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'type', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define other global attributes to NetCDF file.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'title',               &
     &                        TRIM(title))
          IF (FoundError(status, nf90_noerr, 362,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'title', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'var_info',            &
     &                        TRIM(varname))
          IF (FoundError(status, nf90_noerr, 373,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'var_info', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'rst_file',            &
     &                        TRIM(RST(ng)%name))
          IF (FoundError(status, nf90_noerr, 411,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'rst_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (LdefHIS(ng)) THEN
            IF (ndefHIS(ng).gt.0) THEN
              status=nf90_put_att(ncid, nf90_global, 'his_base',        &
     &                            TRIM(HIS(ng)%base))
            ELSE
              status=nf90_put_att(ncid, nf90_global, 'his_file',        &
     &                            TRIM(HIS(ng)%name))
            END IF
            IF (FoundError(status, nf90_noerr, 428,                     &
     &                     "ROMS/Utility/def_info.F")) THEN
              IF (Master) WRITE (stdout,20) 'his_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'ini_file',            &
     &                        TRIM(INI(ng)%name))
          IF (FoundError(status, nf90_noerr, 581,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'ini_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (ObcData(ng)) THEN
          DO i=1,nBCfiles(ng)
            IF (exit_flag.eq.NoError) THEN
              CALL join_string (BRY(i,ng)%files, BRY(i,ng)%Nfiles,      &
     &                          string, lstr)
              WRITE (bryatt,50) i
              status=nf90_put_att(ncid, nf90_global, bryatt,            &
     &                            string(1:lstr))
              IF (FoundError(status, nf90_noerr, 703,                   &
     &                     "ROMS/Utility/def_info.F")) THEN
                IF (Master) WRITE (stdout,20) TRIM(bryatt), TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END DO
        END IF
        IF (Lclimatology(ng)) THEN
          DO i=1,nCLMfiles(ng)
            IF (exit_flag.eq.NoError) THEN
              CALL join_string (CLM(i,ng)%files, CLM(i,ng)%Nfiles,      &
     &                          string, lstr)
              WRITE (clmatt,40) i
              status=nf90_put_att(ncid, nf90_global, clmatt,            &
     &                            string(1:lstr))
              IF (FoundError(status, nf90_noerr, 723,                   &
     &                     "ROMS/Utility/def_info.F")) THEN
                IF (Master) WRITE (stdout,20) TRIM(clmatt), TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END DO
        END IF
!
!  NLM Lateral boundary conditions.
!
        IF (exit_flag.eq.NoError) THEN
          CALL lbc_putatt (ng, ncid, ncname, 'NLM_LBC', LBC, status)
          IF (FoundError(status, nf90_noerr, 848,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'NLM_LBC', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  SVN repository information.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'svn_url',             &
     &                        TRIM(svn_url))
          IF (FoundError(status, nf90_noerr, 877,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'svn_url', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'svn_rev',             &
     &                        TRIM(svn_rev))
          IF (FoundError(status, nf90_noerr, 889,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'svn_rev', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Local root directory, cpp header directory and file, and analytical
!  directory
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'code_dir',            &
     &                        TRIM(Rdir))
          IF (FoundError(status, nf90_noerr, 906,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'code_dir', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'header_dir',          &
     &                        TRIM(Hdir))
          IF (FoundError(status, nf90_noerr, 919,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'header_dir', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'header_file',         &
     &                        TRIM(Hfile))
          IF (FoundError(status, nf90_noerr, 932,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'header_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Attributes describing platform and compiler
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'os',                  &
     &                        TRIM(my_os))
          IF (FoundError(status, nf90_noerr, 948,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'os', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'cpu',                 &
     &                        TRIM(my_cpu))
          IF (FoundError(status, nf90_noerr, 959,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'cpu', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_system',     &
     &                        TRIM(my_fort))
          IF (FoundError(status, nf90_noerr, 970,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'compiler_system',            &
     &                                    TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_command',    &
     &                        TRIM(my_fc))
          IF (FoundError(status, nf90_noerr, 982,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'compiler_command',           &
     &                                    TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_flags',      &
     &                        TRIM(my_fflags))
          IF (FoundError(status, nf90_noerr, 994,                       &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'compiler_flags', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Tiling and history attributes.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'tiling',              &
     &                        TRIM(tiling))
          IF (FoundError(status, nf90_noerr, 1007,                      &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'tiling', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'history',             &
     &                        TRIM(history))
          IF (FoundError(status, nf90_noerr, 1018,                      &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'history', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Analytical header files used.
!
        IF (exit_flag.eq.NoError) THEN
          CALL join_string (ANANAME, SIZE(ANANAME), string, lstr)
          IF (lstr.gt.0) THEN
            status=nf90_put_att(ncid, nf90_global, 'ana_file',          &
     &                          string(1:lstr))
            IF (FoundError(status, nf90_noerr, 1033,                    &
     &                     "ROMS/Utility/def_info.F")) THEN
              IF (Master) WRITE (stdout,20) 'ana_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  Activated CPP options.
!
        IF (exit_flag.eq.NoError) THEN
          lstr=LEN_TRIM(Coptions)-1
          status=nf90_put_att(ncid, nf90_global, 'CPP_options',         &
     &                        TRIM(Coptions(1:lstr)))
          IF (FoundError(status, nf90_noerr, 1070,                      &
     &                   "ROMS/Utility/def_info.F")) THEN
            IF (Master) WRITE (stdout,20) 'CPP_options', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
      IF (FoundError(exit_flag, NoError, 1087,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!-----------------------------------------------------------------------
!  Define running parameters.
!-----------------------------------------------------------------------
!
!  Time stepping parameters.
!
      Vinfo( 1)='ntimes'
      Vinfo( 2)='number of long time-steps'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1109,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='ndtfast'
      Vinfo( 2)='number of short time-steps'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1117,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='dt'
      Vinfo( 2)='size of long time-steps'
      Vinfo( 3)='second'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1126,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='dtfast'
      Vinfo( 2)='size of short time-steps'
      Vinfo( 3)='second'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1135,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='dstart'
      Vinfo( 2)='time stamp assigned to model initilization'
      WRITE (Vinfo( 3),'(a,a)') 'days since ', TRIM(Rclock%string)
      Vinfo( 4)=TRIM(Rclock%calendar)
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1145,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='nHIS'
      Vinfo( 2)='number of time-steps between history records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1179,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='ndefHIS'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of history files'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1188,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='nRST'
      Vinfo( 2)='number of time-steps between restart records'
      IF (LcycleRST(ng)) THEN
        Vinfo(13)='only latest two records are maintained'
      END IF
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1199,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      Vinfo( 1)='Falpha'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1461,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Fbeta'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1469,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Fgamma'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1477,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  Horizontal mixing coefficients.
!
!
!  Background vertical mixing coefficients.
!
      Vinfo( 1)='Akt_bak'
      Vinfo( 2)='background vertical mixing coefficient for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1681,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Akv_bak'
      Vinfo( 2)='background vertical mixing coefficient for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1690,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Akk_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for turbulent energy'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1701,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Akp_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for length scale'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1711,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  Drag coefficients.
!
      Vinfo( 1)='rdrg'
      Vinfo( 2)='linear drag coefficient'
      Vinfo( 3)='meter second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1773,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='rdrg2'
      Vinfo( 2)='quadratic drag coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo ,ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1781,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Zob'
      Vinfo( 2)='bottom roughness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1791,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Zos'
      Vinfo( 2)='surface roughness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1800,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  Generic length-scale parameters.
!
      Vinfo( 1)='gls_p'
      Vinfo( 2)='stability exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1812,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='gls_m'
      Vinfo( 2)='turbulent kinetic energy exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1820,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='gls_n'
      Vinfo( 2)='turbulent length scale exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1828,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='gls_cmu0'
      Vinfo( 2)='stability coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1836,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='gls_c1'
      Vinfo( 2)='shear production coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1844,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='gls_c2'
      Vinfo( 2)='dissipation coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1852,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='gls_c3m'
      Vinfo( 2)='buoyancy production coefficient (minus)'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1860,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='gls_c3p'
      Vinfo( 2)='buoyancy production coefficient (plus)'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1868,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='gls_sigk'
      Vinfo( 2)='constant Schmidt number for TKE'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1876,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='gls_sigp'
      Vinfo( 2)='constant Schmidt number for PSI'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1884,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='gls_Kmin'
      Vinfo( 2)='minimum value of specific turbulent kinetic energy'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1892,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='gls_Pmin'
      Vinfo( 2)='minimum Value of dissipation'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1900,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Charnok_alpha'
      Vinfo( 2)='Charnok factor for surface roughness'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1908,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Zos_hsig_alpha'
      Vinfo( 2)='wave amplitude factor for surface roughness'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1916,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='sz_alpha'
      Vinfo( 2)='surface flux from wave dissipation'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1924,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='CrgBan_cw'
      Vinfo( 2)='surface flux due to Craig and Banner wave breaking'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1932,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  Nudging inverse time scales used in various tasks.
!
      Vinfo( 1)='Znudg'
      Vinfo( 2)='free-surface nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1954,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='M2nudg'
      Vinfo( 2)='2D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1963,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='M3nudg'
      Vinfo( 2)='3D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1973,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Tnudg'
      Vinfo( 2)='Tracers nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1982,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Tnudg_SSS'
      Vinfo( 2)='SSS nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1991,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        Vinfo( 1)='FSobc_in'
        Vinfo( 2)='free-surface inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2017,                        &
     &                 "ROMS/Utility/def_info.F")) RETURN
        Vinfo( 1)='FSobc_out'
        Vinfo( 2)='free-surface outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2026,                        &
     &                 "ROMS/Utility/def_info.F")) RETURN
        Vinfo( 1)='M2obc_in'
        Vinfo( 2)='2D momentum inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2035,                        &
     &                 "ROMS/Utility/def_info.F")) RETURN
        Vinfo( 1)='M2obc_out'
        Vinfo( 2)='2D momentum outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2044,                        &
     &                 "ROMS/Utility/def_info.F")) RETURN
        Vinfo( 1)='Tobc_in'
        Vinfo( 2)='tracers inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 2, tbrydim, Aval, Vinfo, ncname,                 &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2054,                        &
     &                 "ROMS/Utility/def_info.F")) RETURN
        Vinfo( 1)='Tobc_out'
        Vinfo( 2)='tracers outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 2, tbrydim, Aval, Vinfo, ncname,                 &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2063,                        &
     &                 "ROMS/Utility/def_info.F")) RETURN
        Vinfo( 1)='M3obc_in'
        Vinfo( 2)='3D momentum inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2072,                        &
     &                 "ROMS/Utility/def_info.F")) RETURN
        Vinfo( 1)='M3obc_out'
        Vinfo( 2)='3D momentum outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2081,                        &
     &                 "ROMS/Utility/def_info.F")) RETURN
      END IF
!
!  Equation of State parameters.
!
      Vinfo( 1)='rho0'
      Vinfo( 2)='mean density used in Boussinesq approximation'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2095,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='R0'
      Vinfo( 2)='background density used in linear equation of state'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2117,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Tcoef'
      Vinfo( 2)='thermal expansion coefficient'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2126,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Scoef'
      Vinfo( 2)='Saline contraction coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2134,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  Various parameters.
!
!
!  Slipperiness parameters.
!
      Vinfo( 1)='gamma2'
      Vinfo( 2)='slipperiness parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2167,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
! Logical switches to activate horizontal momentum transport
! point Sources/Sinks (like river runoff transport) and mass point
! Sources/Sinks (like volume vertical influx).
!
      Vinfo( 1)='LuvSrc'
      Vinfo( 2)='momentum point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2181,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='LwSrc'
      Vinfo( 2)='mass point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2191,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  Logical switches indicating which tracer variables are processed
!  during point Sources/Sinks.
!
      Vinfo( 1)='LtracerSrc'
      Vinfo( 2)='tracer point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2206,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  Logical switches to process climatology fields.
!
      Vinfo( 1)='LsshCLM'
      Vinfo( 2)='sea surface height climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2219,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Lm2CLM'
      Vinfo( 2)='2D momentum climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2229,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Lm3CLM'
      Vinfo( 2)='3D momentum climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2240,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='LtracerCLM'
      Vinfo( 2)='tracer climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2250,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  Logical switches for nudging of climatology fields.
!
      Vinfo( 1)='LnudgeM2CLM'
      Vinfo( 2)='2D momentum climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2263,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
      Vinfo( 1)='LnudgeM3CLM'
      Vinfo( 2)='3D momentum climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2274,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
      Vinfo( 1)='LnudgeTCLM'
      Vinfo( 2)='tracer climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2284,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  Define sediment model parameters.
!
      Vinfo( 1)='minlayer_thick'
      Vinfo( 2)='depositional bed layer minimum thickness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 24,                            &
     &               "ROMS/Nonlinear/Sediment/sediment_def.h")) RETURN
      Vinfo( 1)='newlayer_thick'
      Vinfo( 2)='depositional bed thickness criteria to crate a '//     &
     &          'new layer'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 34,                            &
     &               "ROMS/Nonlinear/Sediment/sediment_def.h")) RETURN
      Vinfo( 1)='Sd50'
      Vinfo( 2)='median sediment grain diameter used in '//             &
     &          'uniform initial conditions'
      Vinfo( 3)='millimeter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 55,                            &
     &               "ROMS/Nonlinear/Sediment/sediment_def.h")) RETURN
      Vinfo( 1)='Srho'
      Vinfo( 2)='sediment grain density used in '//                     &
     &          'uniform initial conditions'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 65,                            &
     &               "ROMS/Nonlinear/Sediment/sediment_def.h")) RETURN
      Vinfo( 1)='Csed'
      Vinfo( 2)='sediment concentration used in '//                     &
     &          'uniform initial conditions'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 75,                            &
     &               "ROMS/Nonlinear/Sediment/sediment_def.h")) RETURN
      Vinfo( 1)='Wsed'
      Vinfo( 2)='sediment particle settling velocity'
      Vinfo( 3)='meter second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 85,                            &
     &               "ROMS/Nonlinear/Sediment/sediment_def.h")) RETURN
      Vinfo( 1)='Erate'
      Vinfo( 2)='sediment surface layer erosing rate'
      Vinfo( 3)='kilogram meter-2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 94,                            &
     &               "ROMS/Nonlinear/Sediment/sediment_def.h")) RETURN
      Vinfo( 1)='tau_ce'
      Vinfo( 2)='sediment critical shear for erosion'
      Vinfo( 3)='Newton meter-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 103,                           &
     &               "ROMS/Nonlinear/Sediment/sediment_def.h")) RETURN
      Vinfo( 1)='tau_cd'
      Vinfo( 2)='sediment critical shear for deposition'
      Vinfo( 3)='Newton meter-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 112,                           &
     &               "ROMS/Nonlinear/Sediment/sediment_def.h")) RETURN
      Vinfo( 1)='poros'
      Vinfo( 2)='sediment porosity'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 120,                           &
     &               "ROMS/Nonlinear/Sediment/sediment_def.h")) RETURN
! DDMITRY
! Number of oil components 
!      Vinfo( 1)='Noil_cmp'
!      Vinfo( 2)='Number of oil components'
!      Vinfo(14)='Noil_cmp, scalar'
!      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
!     &               1, (/ncmpdim/), Aval, Vinfo, ncname,               &
!     &               SetParAccess = .FALSE.)
!      IF (FoundError(exit_flag, NoError, 3040,                         &
!     &               "ROMS/Utility/def_info.F")) RETURN
! END DD
!
!-----------------------------------------------------------------------
!  Define grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch: Spherical or Cartesian. Writing characters in
!  parallel I/O is extremely inefficient.  It is better to write
!  this as an integer switch: 0=Cartesian, 1=spherical.
!
      Vinfo( 1)='spherical'
      Vinfo( 2)='grid type logical switch'
      Vinfo( 9)='Cartesian'
      Vinfo(10)='spherical'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3064,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  Domain Length.
!
      Vinfo( 1)='xl'
      Vinfo( 2)='domain length in the XI-direction'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3075,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='el'
      Vinfo( 2)='domain length in the ETA-direction'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3084,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  S-coordinate parameters.
!
      Vinfo( 1)='Vtransform'
      Vinfo( 2)='vertical terrain-following transformation equation'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3095,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Vstretching'
      Vinfo( 2)='vertical terrain-following stretching function'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3103,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='theta_s'
      Vinfo( 2)='S-coordinate surface control parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3111,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='theta_b'
      Vinfo( 2)='S-coordinate bottom control parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3119,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='Tcline'
      Vinfo( 2)='S-coordinate surface/bottom layer width'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3128,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
      Vinfo( 1)='hc'
      Vinfo( 2)='S-coordinate parameter, critical depth'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3137,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  SGRID conventions for staggered data on structured grids.
!
      Vinfo( 1)='grid'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3146,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  S-coordinate non-dimensional independent variable at RHO-points.
!
      Vinfo( 1)='s_rho'
      Vinfo( 2)='S-coordinate at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='s_rho, scalar'
      IF (Vtransform(ng).eq.1) THEN
        Vinfo(21)='ocean_s_coordinate_g1'
      ELSE IF (Vtransform(ng).eq.2) THEN
        Vinfo(21)='ocean_s_coordinate_g2'
      END IF
      Vinfo(23)='s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
      vinfo(25)='up'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/srdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3172,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  S-coordinate non-dimensional independent variable at W-points.
!
      Vinfo( 1)='s_w'
      Vinfo( 2)='S-coordinate at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='s_w, scalar'
      Vinfo(21)='ocean_s_coordinate'
      IF (Vtransform(ng).eq.1) THEN
        Vinfo(21)='ocean_s_coordinate_g1'
      ELSE IF (Vtransform(ng).eq.2) THEN
        Vinfo(21)='ocean_s_coordinate_g2'
      END IF
      Vinfo(23)='s: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
      vinfo(25)='up'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/swdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3199,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  S-coordinate non-dimensional stretching curves at RHO-points.
!
      Vinfo( 1)='Cs_r'
      Vinfo( 2)='S-coordinate stretching curves at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='Cs_r, scalar'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/srdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3214,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  S-coordinate non-dimensional stretching curves at W-points.
!
      Vinfo( 1)='Cs_w'
      Vinfo( 2)='S-coordinate stretching curves at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='Cs_w, scalar'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/swdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3229,                          &
     &               "ROMS/Utility/def_info.F")) RETURN
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        Vinfo( 1)='user'
        Vinfo( 2)='user generic parameters'
        Vinfo(14)='user, scalar'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/usrdim/), Aval, Vinfo, ncname,              &
     &               SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 3242,                        &
     &                 "ROMS/Utility/def_info.F")) RETURN
      END IF
      IF (ncid.ne.FLT(ng)%ncid) THEN
!
!  Bathymetry.
!
        Vinfo( 1)='h'
        Vinfo( 2)='bathymetry at RHO-points'
        Vinfo( 3)='meter'
        Vinfo(14)='bath, scalar'
        Vinfo(22)='coordinates'
        Aval(5)=REAL(r2dvar,r8)
        IF (ncid.eq.STA(ng)%ncid) THEN
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   1, (/stadim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3285,                      &
     &                   "ROMS/Utility/def_info.F")) RETURN
        ELSE
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3290,                      &
     &                   "ROMS/Utility/def_info.F")) RETURN
        END IF
!
!  Coriolis Parameter.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)='f'
          Vinfo( 2)='Coriolis parameter at RHO-points'
          Vinfo( 3)='second-1'
          Vinfo(14)='coriolis, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3323,                      &
     &                   "ROMS/Utility/def_info.F")) RETURN
        END IF
!
!  Curvilinear coordinate metrics.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)='pm'
          Vinfo( 2)='curvilinear coordinate metric in XI'
          Vinfo( 3)='meter-1'
          Vinfo(14)='pm, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3338,                      &
     &                   "ROMS/Utility/def_info.F")) RETURN
          Vinfo( 1)='pn'
          Vinfo( 2)='curvilinear coordinate metric in ETA'
          Vinfo( 3)='meter-1'
          Vinfo(14)='pn, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3349,                      &
     &                   "ROMS/Utility/def_info.F")) RETURN
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_rho'
          Vinfo( 2)='longitude of RHO-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_rho, scalar'
          Vinfo(21)='longitude'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3364,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3369,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
          Vinfo( 1)='lat_rho'
          Vinfo( 2)='latitude of RHO-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_rho, scalar'
          Vinfo(21)='latitude'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo,  ncname)
            IF (FoundError(exit_flag, NoError, 3381,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3386,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_rho'
          Vinfo( 2)='x-locations of RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_rho, scalar'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3397,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3402,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
          Vinfo( 1)='y_rho'
          Vinfo( 2)='y-locations of RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_rho, scalar'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3413,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3418,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_u'
          Vinfo( 2)='longitude of U-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_u, scalar'
          Vinfo(21)='longitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3434,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
          Vinfo( 1)='lat_u'
          Vinfo( 2)='latitude of U-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_u, scalar'
          Vinfo(21)='latitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3446,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_u'
          Vinfo( 2)='x-locations of U-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_u, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3457,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
          Vinfo( 1)='y_u'
          Vinfo( 2)='y-locations of U-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_u, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3468,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_v'
          Vinfo( 2)='longitude of V-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_v, scalar'
          Vinfo(21)='longitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3484,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
          Vinfo( 1)='lat_v'
          Vinfo( 2)='latitude of V-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_v, scalar'
          Vinfo(21)='latitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3496,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_v'
          Vinfo( 2)='x-locations of V-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_v, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3507,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
          Vinfo( 1)='y_v'
          Vinfo( 2)='y-locations of V-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_v, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3518,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_psi'
          Vinfo( 2)='longitude of PSI-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_psi, scalar'
          Vinfo(21)='longitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3534,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
          Vinfo( 1)='lat_psi'
          Vinfo( 2)='latitude of PSI-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_psi, scalar'
          Vinfo(21)='latitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3546,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_psi'
          Vinfo( 2)='x-locations of PSI-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_psi, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3557,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
          Vinfo( 1)='y_psi'
          Vinfo( 2)='y-locations of PSI-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_psi, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3568,                    &
     &                     "ROMS/Utility/def_info.F")) RETURN
          END IF
        END IF
      END IF
  10  FORMAT (i3.3,'x',i3.3)
  20  FORMAT (/,' DEF_INFO - error while creating global attribute: ',  &
     &        a,/,12x,a)
  30  FORMAT ('frc_file_',i2.2)
  40  FORMAT ('clm_file_',i2.2)
  50  FORMAT ('bry_file_',i2.2)
      RETURN
      END SUBROUTINE def_info
