      SUBROUTINE read_SedPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in cohesive and non-cohesive sediment model      !
!  parameters.                                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
      USE mod_sediment
      USE mod_sedflocs
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval
      integer :: iTrcStr, iTrcEnd
      integer :: i, ifield, igrid, itracer, itrc, ng, nline, status
      integer :: decode_line, load_i, load_l, load_lbc, load_r
      logical, dimension(Ngrids) :: Lbed
      logical, dimension(Ngrids) :: Lbottom
      logical, dimension(NCS,Ngrids) :: Lmud
      logical, dimension(NNS,Ngrids) :: Lsand
      real(r8), dimension(Ngrids) :: Rbed
      real(r8), dimension(NCS,Ngrids) :: Rmud
      real(r8), dimension(NNS,Ngrids) :: Rsand
      real(r8), dimension(100) :: Rval
      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(200) :: Cval
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      igrid=1                            ! nested grid counter
      itracer=0                          ! LBC tracer counter
      iTrcStr=1                          ! first LBC tracer to process
      iTrcEnd=NST                        ! last  LBC tracer to process
      nline=0                            ! LBC multi-line counter
!
!-----------------------------------------------------------------------
!  Read in cohesive and non-cohesive model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Lsediment')
              Npts=load_l(Nval, Cval, Ngrids, Lsediment)
            CASE ('NEWLAYER_THICK')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                newlayer_thick(ng)=Rbed(ng)
              END DO
            CASE ('MINLAYER_THICK')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                minlayer_thick(ng)=Rbed(ng)
              END DO
            CASE ('BEDLOAD_COEFF')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                bedload_coeff(ng)=Rbed(ng)
              END DO
            CASE ('LBC(isTvar)')
              IF (itracer.lt.NST) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idsed(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      idsed(iTrcStr), idsed(iTrcEnd),             &
     &                      Vname(1,idTvar(idsed(itracer))), LBC)
            CASE ('MUD_SD50')
              IF (.not.allocated(Sd50)) allocate (Sd50(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Sd50(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_CSED')
              IF (.not.allocated(Csed)) allocate (Csed(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud )
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Csed(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_SRHO')
              IF (.not.allocated(Srho)) allocate (Srho(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Srho(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_WSED')
              IF (.not.allocated(Wsed)) allocate (Wsed(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Wsed(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_ERATE')
              IF (.not.allocated(Erate)) allocate (Erate(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Erate(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TAU_CE')
              IF (.not.allocated(tau_ce)) allocate (tau_ce(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  tau_ce(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TAU_CD')
              IF (.not.allocated(tau_cd)) allocate (tau_cd(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  tau_cd(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_POROS')
              IF (.not.allocated(poros)) allocate (poros(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  poros(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TNU2')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  nl_tnu2(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TNU4')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  nl_tnu4(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('ad_MUD_TNU2')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  ad_tnu2(i,ng)=Rmud(itrc,ng)
                  tl_tnu2(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('ad_MUD_TNU4')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  ad_tnu4(i,ng)=Rmud(itrc,ng)
                  nl_tnu4(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_Sponge')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  LtracerSponge(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_AKT_BAK')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  Akt_bak(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_AKT_fac')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  ad_Akt_fac(i,ng)=Rmud(itrc,ng)
                  tl_Akt_fac(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TNUDG')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  Tnudg(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_MORPH_FAC')
              IF (.not.allocated(morph_fac)) THEN
                allocate (morph_fac(NST,Ngrids))
              END IF
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  morph_fac(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_Ltsrc', 'MUD_Ltracer')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  LtracerSrc(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_Ltclm')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  LtracerCLM(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_Tnudge')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  LnudgeTCLM(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idmud)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idTvar(idsed(itrc))
                  Hout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iMfrac)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idfrac(itrc)
                  Hout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iMmass)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idBmas(itrc)
                  Hout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idmud)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idTvar(idsed(itrc))
                  Qout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iSmud)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsurT(idsed(itrc))
                  Qout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iMfrac)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idfrac(itrc)
                  Qout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iMmass)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idBmas(itrc)
                  Hout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('SAND_SD50')
              IF (.not.allocated(Sd50)) allocate (Sd50(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Sd50(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_CSED')
              IF (.not.allocated(Csed)) allocate (Csed(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand )
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Csed(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_SRHO')
              IF (.not.allocated(Srho)) allocate (Srho(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Srho(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_WSED')
              IF (.not.allocated(Wsed)) allocate (Wsed(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Wsed(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_ERATE')
              IF (.not.allocated(Erate)) allocate (Erate(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Erate(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TAU_CE')
              IF (.not.allocated(tau_ce)) allocate (tau_ce(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  tau_ce(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TAU_CD')
              IF (.not.allocated(tau_cd)) allocate (tau_cd(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  tau_cd(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_POROS')
              IF (.not.allocated(poros)) allocate (poros(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  poros(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TNU2')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  nl_tnu2(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TNU4')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  nl_tnu4(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('ad_SAND_TNU2')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  ad_tnu2(i,ng)=Rsand(itrc,ng)
                  tl_tnu2(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('ad_SAND_TNU4')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  ad_tnu4(i,ng)=Rsand(itrc,ng)
                  tl_tnu4(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_Sponge')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  LtracerSponge(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_AKT_BAK')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  Akt_bak(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_AKT_fac')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  ad_Akt_fac(i,ng)=Rsand(itrc,ng)
                  tl_Akt_fac(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TNUDG')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  Tnudg(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_MORPH_FAC')
              IF (.not.allocated(morph_fac)) THEN
                allocate (morph_fac(NST,Ngrids))
              END IF
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  morph_fac(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_Ltsrc', 'SAND_Ltracer')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  LtracerSrc(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_Ltclm')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  LtracerCLM(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_Tnudge')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  LnudgeTCLM(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idsand)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idTvar(idsed(NCS+itrc))
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iSfrac)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idfrac(NCS+itrc)
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iSmass)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idBmas(NCS+itrc)
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idsand)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idTvar(idsed(NCS+itrc))
                  Qout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iSsand)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsurT(idsed(NCS+itrc))
                  Qout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iSfrac)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idfrac(NCS+itrc)
                  Qout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iSmass)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idBmas(NCS+itrc)
                  Qout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(ithck)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(ithck)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
            CASE ('Hout(iaged)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(iaged)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
            CASE ('Hout(iporo)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(iporo)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
            CASE ('Hout(idiff)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(idiff)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
            CASE ('Hout(isd50)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(isd50)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idens)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idens)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
             END DO
            CASE ('Hout(iwsed)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(iwsed)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(itauc)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(itauc)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(irlen)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(irlen)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(irhgt)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(irhgt)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(ibwav)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(ibwav)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izdef)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izdef)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izapp)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izapp)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izNik)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izNik)
             DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izbio)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbio)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izbfm)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbfm)
              DO ng=1,Ngrids
               Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izbld)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbld)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izwbl)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izwbl)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(iactv)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(iactv)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(ishgt)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(ishgt)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(imaxD)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(imaxD)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idnet)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idnet)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(ithck)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(ithck)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbed(ng)
              END DO
            CASE ('Qout(iaged)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(iaged)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbed(ng)
              END DO
            CASE ('Qout(iporo)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(iporo)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbed(ng)
              END DO
            CASE ('Qout(idiff)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(idiff)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbed(ng)
              END DO
            CASE ('Qout(isd50)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(isd50)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(idens)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idens)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
             END DO
            CASE ('Qout(iwsed)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(iwsed)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(itauc)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(itauc)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(irlen)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(irlen)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(irhgt)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(irhgt)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(ibwav)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(ibwav)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izdef)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izdef)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izapp)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izapp)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izNik)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izNik)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izbio)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbio)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izbfm)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbfm)
              DO ng=1,Ngrids
               Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izbld)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbld)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izwbl)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izwbl)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(iactv)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(iactv)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(ishgt)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(ishgt)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,40) line
      exit_flag=4
      RETURN
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lsediment(ng)) THEN
            WRITE (out,50) ng
            WRITE (out,60)
            DO itrc=1,NST
              WRITE (out,70) itrc, Sd50(itrc,ng), Csed(itrc,ng),        &
     &                       Srho(itrc,ng), Wsed(itrc,ng),              &
     &                       Erate(itrc,ng), poros(itrc,ng)
            END DO
            WRITE (out,80)
            DO itrc=1,NST
              i=idsed(itrc)
              WRITE (out,70) itrc, tau_ce(itrc,ng), tau_cd(itrc,ng),    &
     &                       nl_tnu2(i,ng), nl_tnu4(i,ng),              &
     &                       Akt_bak(i,ng), Tnudg(i,ng)
            END DO
            WRITE (out,90)
            DO itrc=1,NST
              WRITE (out,70) itrc,  morph_fac(itrc,ng)
            END DO
            WRITE (out,100) newlayer_thick(ng)
            WRITE (out,110) minlayer_thick(ng)
            WRITE (out,120) bedload_coeff(ng)
            DO itrc=1,NST
              i=idsed(itrc)
              IF (LtracerSponge(i,ng)) THEN
                WRITE (out,150) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning ON  sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning OFF sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NST
              i=idsed(itrc)
              IF (LtracerSrc(i,ng)) THEN
                WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc', i,      &
     &              'Turning ON  point sources/Sink on tracer ', i,     &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc', i,      &
     &              'Turning OFF point sources/Sink on tracer ', i,     &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NST
              i=idsed(itrc)
              IF (LtracerCLM(i,ng)) THEN
                WRITE (out,150) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning ON  processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning OFF processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NST
              i=idsed(itrc)
              IF (LnudgeTCLM(i,ng)) THEN
                WRITE (out,150) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning ON  nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning OFF nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            IF ((nHIS(ng).gt.0).and.ANY(Hout(:,ng))) THEN
              WRITE (out,'(1x)')
              DO itrc=1,NST
                i=idTvar(idsed(itrc))
                IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),             &
     &              'Hout(idTvar)',                                     &
     &              'Write out sediment', itrc, TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idfrac(itrc)
                IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),             &
     &              'Hout(idfrac)',                                     &
     &              'Write out bed fraction, sediment ', itrc,          &
     &              TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idBmas(itrc)
                IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),             &
     &              'Hout(idfrac)',                                     &
     &              'Write out mass, sediment ', itrc,                  &
     &              TRIM(Vname(1,i))
              END DO
              DO itrc=1,MBEDP
                i=idSbed(itrc)
                IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),             &
     &              'Hout(idSbed)',                                     &
     &              'Write out BED property ', itrc, TRIM(Vname(1,i))
              END DO
              DO itrc=1,MBOTP
                i=idBott(itrc)
                IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),             &
     &             'Hout(idBott)',                                      &
     &             'Write out BOTTOM property ', itrc, TRIM(Vname(1,i))
              END DO
            END IF
            IF ((nQCK(ng).gt.0).and.ANY(Qout(:,ng))) THEN
              WRITE (out,'(1x)')
              DO itrc=1,NST
                i=idTvar(idsed(itrc))
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &              'Qout(idTvar)',                                     &
     &              'Write out sediment', itrc, TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idsurT(idsed(itrc))
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &              'Qout(idTvar)',                                     &
     &              'Write out surface sediment', itrc, TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idfrac(itrc)
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &              'Qout(idfrac)',                                     &
     &              'Write out bed fraction, sediment ', itrc,          &
     &              TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idBmas(itrc)
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &              'Qout(idfrac)',                                     &
     &              'Write out mass, sediment ', itrc,                  &
     &              TRIM(Vname(1,i))
              END DO
              DO itrc=1,MBEDP
                i=idSbed(itrc)
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &              'Qout(idSbed)',                                     &
     &              'Write out BED property ', itrc, TRIM(Vname(1,i))
              END DO
              DO itrc=1,MBOTP
                i=idBott(itrc)
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &             'Qout(idBott)',                                      &
     &             'Write out BOTTOM property ', itrc, TRIM(Vname(1,i))
              END DO
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Scale relevant input parameters
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO i=1,NST
          Sd50(i,ng)=Sd50(i,ng)*0.001_r8
          Wsed(i,ng)=Wsed(i,ng)*0.001_r8
          tau_ce(i,ng)=tau_ce(i,ng)/rho0
          tau_cd(i,ng)=tau_cd(i,ng)/rho0
          nl_tnu4(idsed(i),ng)=SQRT(ABS(nl_tnu4(idsed(i),ng)))
          IF (Tnudg(idsed(i),ng).gt.0.0_r8) THEN
            Tnudg(idsed(i),ng)=1.0_r8/(Tnudg(idsed(i),ng)*86400.0_r8)
          ELSE
            Tnudg(idsed(i),ng)=0.0_r8
          END IF
        END DO
      END DO
  30  FORMAT (/,' READ_SedPar - variable info not yet loaded, ', a)
  40  FORMAT (/,' READ_SedPar - Error while processing line: ',/,a)
  50  FORMAT (/,/,' Sediment Parameters, Grid: ',i2.2,                  &
     &        /,  ' =============================',/)
  60  FORMAT (/,1x,'Size',5x,'Sd50',8x,'Csed',8x,'Srho',8x,'Wsed',      &
     &        8x,'Erate',7x,'poros',/,1x,'Class',4x,'(mm)',7x,          &
     &        '(kg/m3)',5x,'(kg/m3)',5x,'(mm/s)',5x,'(kg/m2/s)',4x,     &
     &        '(nondim)',/)
  70  FORMAT (2x,i2,2x,6(1x,1p,e11.4))
  80  FORMAT (/,9x,'tau_ce',6x,'tau_cd',6x,'nl_tnu2',5x,'nl_tnu4',5x,   &
     &        'Akt_bak',6x,'Tnudg',/,9x,'(N/m2)',6x,'(N/m2)',6x,        &
     &        '(m2/s)',6x,'(m4/s)',7x,'(m2/s)',6x,'(day)',/)
  90  FORMAT (/,9x,'morph_fac',/,9x,'(nondim)',/)
 100  FORMAT (/,' New bed layer formed when deposition exceeds ',e12.5, &
     &        ' (m).')
 110  FORMAT (' Two first layers are combined when 2nd layer smaller ', &
     &         'than ',e12.5,' (m).')
 120  FORMAT (' Rate coefficient for bed load transport = ',e12.5,/)
 130  FORMAT (' Transition for mixed sediment =',e12.5,/)
 140  FORMAT (' Transition for cohesive sediment =',e12.5,/)
 150  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)
 160  FORMAT (10x,l1,2x,a,t29,a,i2.2,':',1x,a)
      RETURN
      END SUBROUTINE read_SedPar
