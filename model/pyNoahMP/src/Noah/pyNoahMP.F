subroutine initialize(MMINLU,SNOW,SNOWH,CANWAT,ISLTYP,IVGTYP,XLAT, &
       TSLB,SMOIS,SH2O,DZS,FNDSOILW,FNDSNOWH,&
       TSK,isnowxy,tvxy,tgxy,canicexy,TMN,XICE,&
       canliqxy ,eahxy    ,tahxy    ,cmxy     ,chxy     ,                     &
       fwetxy   ,sneqvoxy ,alboldxy ,qsnowxy  ,wslakexy ,zwtxy    ,waxy     , &
       wtxy     ,tsnoxy   ,zsnsoxy  ,snicexy  ,snliqxy  ,lfmassxy ,rtmassxy , &
       stmassxy ,woodxy   ,stblcpxy ,fastcpxy ,xsaixy   ,lai      ,           &
       grainxy  ,gddxy,croptype ,cropcat  ,t2mvxy   ,t2mbxy   ,chstarxy,            &
       restart,                 &
       allowed_to_read , iopt_run,  iopt_crop,                        &
       sf_urban_physics,                         &
       smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, areaxy, dx, dy, msftx, msfty,&
       wtddt    ,stepwtd  ,dt       ,qrfsxy     ,qspringsxy  , qslatxy    ,  &
       fdepthxy ,ht     ,riverbedxy ,eqzwt     ,rivercondxy ,pexpxy       ,  &
       rechclim,                                                             &
       gecros_state,nsoil,nhru, &
       soilparm_file,genparm_file,mptable_file,&
       dzb,dzg,dzr,isurban,low_density_residential_table,high_density_residential_table,&
       high_intensity_industrial_table,num_urban_ndm,urban_map_zrd,urban_map_zwd,&
       urban_map_gd,urban_map_zd,urban_map_zdf,urban_map_bd,urban_map_wd,&
       urban_map_gbd,urban_map_fbd,num_urban_hi,tsurface0_urb,tlayer0_urb,&
       tdeep0_urb,tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d,xxxr_urb2d,xxxb_urb2d,&
       xxxg_urb2d,xxxc_urb2d,drelr_urb2d,drelb_urb2d,drelg_urb2d,flxhumr_urb2d,&
       flxhumb_urb2d,flxhumg_urb2d,cmcr_urb2d,tgr_urb2d,trl_urb3d,tbl_urb3d,&
       tgl_urb3d,tgrl_urb3d,smr_urb3d,sh_urb2d,lh_urb2d,g_urb2d,rn_urb2d,ts_urb2d,& 
       trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,tlev_urb3d,qlev_urb3d,tw1lev_urb3d,&
       tw2lev_urb3d,tglev_urb3d,tflev_urb3d,lf_ac_urb3d,sf_ac_urb3d,cm_ac_urb3d,&
       sfvent_urb3d,lfvent_urb3d,sfwin1_urb3d,sfwin2_urb3d,sfw1_urb3d,sfw2_urb3d,&
       sfr_urb3d,sfg_urb3d,hi_urb2d,lp_urb2d,lb_urb2d,hgt_urb2d,mh_urb2d,stdh_urb2d,&
       lf_urb2d,a_u_bep,a_v_bep,a_t_bep,a_q_bep,a_e_bep,b_u_bep,b_v_bep,b_t_bep,&
       b_q_bep,b_e_bep,vl_bep,dlg_bep,sf_bep,dl_u_bep,frc_urb2d,utype_urb2d,&
       urbtable_file&
       )

 USE module_sf_noahmpdrv, only: noahmp_init
 USE module_sf_urban, only: urban_param_init,urban_var_init

 implicit none

  character(len=256),intent(in) :: soilparm_file,genparm_file,mptable_file,urbtable_file
    INTEGER, INTENT(INOUT)    ::     nhru
    INTEGER, INTENT(INOUT)       ::     NSOIL, iopt_run, iopt_crop
    LOGICAL      ::     restart,allowed_to_read
    INTEGER, INTENT(INOUT)       ::     sf_urban_physics
    REAL,    DIMENSION( NSOIL), INTENT(INOUT)    ::     DZS  ! Thickness of the soil layers [m]
    REAL,    INTENT(INOUT) ::     DX, DY
    REAL,    DIMENSION(nhru) ,  INTENT(INOUT) ::MSFTX,MSFTY
    REAL,    DIMENSION(nhru, NSOIL) ,INTENT(INOUT)::SMOIS,SH2O,TSLB
    REAL,    DIMENSION(nhru),INTENT(INOUT) :: SNOW,SNOWH,CANWAT
  INTEGER, DIMENSION(nhru),                      &
         &   INTENT(INOUT)       ::     ISLTYP,  &
                                     IVGTYP

    LOGICAL, INTENT(INOUT)       ::     FNDSOILW,FNDSNOWH

    character(len=*), intent(in) :: MMINLU
    REAL, DIMENSION(nhru), INTENT(INOUT) :: XLAT         !latitude
    REAL, DIMENSION(nhru), INTENT(INOUT) :: TSK         !skin temperature (k)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: TMN         !deep soil temperature (k)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: XICE         !sea ice fraction
    INTEGER, DIMENSION(nhru), INTENT(INOUT) :: isnowxy     !actual no. of snow layers
    REAL, DIMENSION(nhru,-2:NSOIL), INTENT(INOUT) :: zsnsoxy  !snow layer depth [m]
    REAL, DIMENSION(nhru,-2:0), INTENT(INOUT) :: tsnoxy !snow temperature [K]
    REAL, DIMENSION(nhru,-2:0), INTENT(INOUT) :: snicexy  !snow layer ice [mm]
    REAL, DIMENSION(nhru,-2:0), INTENT(INOUT) :: snliqxy  !snow layer liquid water [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: tvxy        !vegetation canopy temperature
    REAL, DIMENSION(nhru), INTENT(INOUT) :: tgxy        !ground surface temperature
    REAL, DIMENSION(nhru), INTENT(INOUT) :: canicexy !canopy-intercepted ice (mm)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: canliqxy !canopy-intercepted liquid water (mm)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: eahxy       !canopy air vapor pressure (pa)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: tahxy       !canopy air temperature (k)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: cmxy        !momentum drag coefficient
    REAL, DIMENSION(nhru), INTENT(INOUT) :: chxy        !sensible heat exchange coefficient
    REAL, DIMENSION(nhru), INTENT(INOUT) :: fwetxy      !wetted or snowed fraction of the canopy (-)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: sneqvoxy    !snow mass at last time step(mm h2o)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: alboldxy    !snow albedo at last time step (-)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: qsnowxy     !snowfall on the ground [mm/s]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: wslakexy    !lake water storage [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: zwtxy       !water table depth [m]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: waxy        !water in the "aquifer" [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: wtxy        !groundwater storage [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: lfmassxy    !leaf mass [g/m2]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: rtmassxy    !mass of fine roots [g/m2]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: stmassxy    !stem mass [g/m2]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: woodxy      !mass of wood (incl. woody roots) [g/m2]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: grainxy     !mass of grain [g/m2] !XING
    REAL, DIMENSION(nhru), INTENT(INOUT) :: gddxy       !growing degree days !XING
    REAL, DIMENSION(nhru), INTENT(INOUT) :: stblcpxy    !stable carbon in deep soil [g/m2]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: fastcpxy    !short-lived carbon, shallow soil [g/m2]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: xsaixy      !stem area index
    REAL, DIMENSION(nhru), INTENT(INOUT) :: lai         !leaf area index
    INTEGER, DIMENSION(nhru), INTENT(INOUT) :: cropcat
    REAL   , DIMENSION(nhru,5), INTENT(IN ) :: croptype
! IOPT_RUN = 5 option

    REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT)  :: smoiseq !equilibrium soil moisture content [m3m-3]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: smcwtdxy !deep soil moisture content [m3m-3]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: deeprechxy !deep recharge [m]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: rechxy !accumulated recharge [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: qrfsxy !accumulated flux from groundwater to rivers [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: qspringsxy !accumulated seeping water [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: qslatxy !accumulated lateral flow [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: areaxy !grid cell area [m2]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: FDEPTHXY !efolding depth for transmissivity (m)
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: HT !terrain height (m)
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: RIVERBEDXY !riverbed depth (m)
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: EQZWT !equilibrium water table depth (m)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: RIVERCONDXY !river conductance
    REAL, DIMENSION(nhru), INTENT(INOUT) :: PEXPXY !factor for river conductance
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: rechclim
    REAL, DIMENSION(nhru,60), INTENT(INOUT) :: gecros_state ! Optional gecros crop
    INTEGER,  INTENT(INOUT)  :: STEPWTD
    REAL, INTENT(INOUT)  :: DT, WTDDT
!jref:start
    REAL, DIMENSION(nhru), INTENT(INOUT) :: t2mvxy        !2m temperature vegetation part (k)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: t2mbxy        !2m temperature bare ground part (k)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: chstarxy        !dummy
!jref:end
! Urban canopy model
    REAL, DIMENSION(nsoil), INTENT(INOUT) :: DZR
    REAL, DIMENSION(nsoil), INTENT(INOUT) :: DZB
    REAL, DIMENSION(nsoil), INTENT(INOUT) :: DZG
   INTEGER, INTENT(INOUT) :: ISURBAN
   INTEGER, INTENT(INOUT) :: LOW_DENSITY_RESIDENTIAL_TABLE,&
                             HIGH_DENSITY_RESIDENTIAL_TABLE, &
                             HIGH_INTENSITY_INDUSTRIAL_TABLE
   INTEGER, INTENT(INOUT) :: num_urban_ndm
   INTEGER, INTENT(INOUT) :: urban_map_zrd
   INTEGER, INTENT(INOUT) :: urban_map_zwd
   INTEGER, INTENT(INOUT) :: urban_map_gd
   INTEGER, INTENT(INOUT) :: urban_map_zd
   INTEGER, INTENT(INOUT) :: urban_map_zdf
   INTEGER, INTENT(INOUT) :: urban_map_bd
   INTEGER, INTENT(INOUT) :: urban_map_wd
   INTEGER, INTENT(INOUT) :: urban_map_gbd
   INTEGER, INTENT(INOUT) :: urban_map_fbd
   INTEGER, INTENT(INOUT) :: num_urban_hi                                !multi-layer urban
   REAL, DIMENSION(nhru), INTENT(INOUT)                    :: TSURFACE0_URB
   REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT) :: TLAYER0_URB
   REAL, DIMENSION(nhru), INTENT(INOUT)                    :: TDEEP0_URB
   !INTEGER, DIMENSION(nhru), INTENT(INOUT)                 :: IVGTYP
   !LOGICAL , INTENT(INOUT) :: restart
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TR_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TB_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TG_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TC_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: QC_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: XXXR_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: XXXB_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: XXXG_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: XXXC_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: DRELR_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: DRELB_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: DRELG_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: FLXHUMR_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: FLXHUMB_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: FLXHUMG_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CMCR_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TGR_URB2D
   REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT) :: TRL_URB3D
   REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT) :: TBL_URB3D
   REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT) :: TGL_URB3D
   REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT) :: TGRL_URB3D
   REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT) :: SMR_URB3D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: SH_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: LH_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: G_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: RN_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TS_URB2D
! multi-layer UCM variables
   REAL, DIMENSION(nhru,1:urban_map_zrd), INTENT(INOUT) :: TRB_URB4D
   REAL, DIMENSION(nhru,1:urban_map_zwd), INTENT(INOUT) :: TW1_URB4D
   REAL, DIMENSION(nhru,1:urban_map_zwd), INTENT(INOUT) :: TW2_URB4D
   REAL, DIMENSION(nhru,1:urban_map_gd), INTENT(INOUT) :: TGB_URB4D
   REAL, DIMENSION(nhru,1:urban_map_bd), INTENT(INOUT) :: TLEV_URB3D
   REAL, DIMENSION(nhru,1:urban_map_bd), INTENT(INOUT) :: QLEV_URB3D
   REAL, DIMENSION(nhru,1:urban_map_wd), INTENT(INOUT) :: TW1LEV_URB3D
   REAL, DIMENSION(nhru,1:urban_map_wd), INTENT(INOUT) :: TW2LEV_URB3D
   REAL, DIMENSION(nhru,1:urban_map_gbd), INTENT(INOUT) :: TGLEV_URB3D
   REAL, DIMENSION(nhru,1:urban_map_fbd), INTENT(INOUT) :: TFLEV_URB3D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: LF_AC_URB3D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: SF_AC_URB3D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CM_AC_URB3D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: SFVENT_URB3D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: LFVENT_URB3D
   REAL, DIMENSION(nhru,1:urban_map_wd), INTENT(INOUT) :: SFWIN1_URB3D
   REAL, DIMENSION(nhru,1:urban_map_wd), INTENT(INOUT) :: SFWIN2_URB3D
   REAL, DIMENSION(nhru,1:urban_map_zd), INTENT(INOUT) :: SFW1_URB3D
   REAL, DIMENSION(nhru,1:urban_map_zd), INTENT(INOUT) :: SFW2_URB3D
   REAL, DIMENSION(nhru,1:urban_map_zdf), INTENT(INOUT) :: SFR_URB3D
   REAL, DIMENSION(nhru,1:num_urban_ndm), INTENT(INOUT) :: SFG_URB3D
   REAL, DIMENSION(nhru,1:num_urban_hi), INTENT(INOUT) :: HI_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: LP_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: LB_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: HGT_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: MH_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: STDH_URB2D
   REAL, DIMENSION(nhru,4), INTENT(INOUT) :: LF_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: A_U_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: A_V_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: A_T_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: A_Q_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: A_E_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: B_U_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: B_V_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: B_T_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: B_Q_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: B_E_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: VL_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: DLG_BEP
   REAL, DIMENSION(nhru),INTENT(INOUT) :: SF_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: DL_U_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: FRC_URB2D
   INTEGER, DIMENSION(nhru), INTENT(INOUT) :: UTYPE_URB2D

integer :: ids,ide, jds,jde, kds,kde, ims,ime, jms,jme, kms,kme, its,ite
integer :: jts,jte, kts,kte
logical :: debug_flag
debug_flag = .False.

ids = 1
ide = nhru+1
jds = 1
jde = 2
kds = 1
kde = 2
ims = 1
ime = nhru
jms = 1
jme = 1
kms = 1
kme = 2
its = 1
ite = nhru
jts = 1
jte = 1
kts = 1
kte = 2

 call NOAHMP_INIT ( MMINLU, SNOW , SNOWH , CANWAT , ISLTYP ,   IVGTYP,XLAT, &
       TSLB , SMOIS , SH2O , DZS , FNDSOILW , FNDSNOWH ,             &
       TSK, isnowxy , tvxy     ,tgxy     ,canicexy ,         TMN,     XICE,   &
       canliqxy ,eahxy    ,tahxy    ,cmxy     ,chxy     ,                     &
       fwetxy   ,sneqvoxy ,alboldxy ,qsnowxy  ,wslakexy ,zwtxy    ,waxy     , &
       wtxy     ,tsnoxy   ,zsnsoxy  ,snicexy  ,snliqxy  ,lfmassxy ,rtmassxy , &
       stmassxy ,woodxy   ,stblcpxy ,fastcpxy ,xsaixy   ,lai      ,           &
       grainxy  ,gddxy    ,                                                   &
       croptype ,cropcat  ,                      &
       t2mvxy   ,t2mbxy   ,chstarxy,            &
       NSOIL, restart,                 &
       allowed_to_read , iopt_run,  iopt_crop,                        &
       sf_urban_physics,                         &
       ids,ide, jds,jde, kds,kde,                &
       ims,ime, jms,jme, kms,kme,                &
       its,ite, jts,jte, kts,kte,                &
       smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, areaxy, dx, dy, msftx, msfty,&
       wtddt    ,stepwtd  ,dt       ,qrfsxy     ,qspringsxy  , qslatxy    ,  &
       fdepthxy ,ht     ,riverbedxy ,eqzwt     ,rivercondxy ,pexpxy       ,  &
       rechclim,                                                             &
       gecros_state, &
       soilparm_file,genparm_file,mptable_file)

         IF(SF_URBAN_PHYSICS > 0 ) THEN  !urban

                 CALL urban_param_init(DZR,DZB,DZG,nsoil,sf_urban_physics,urbtable_file)
                 CALL urban_var_init(ISURBAN,TSK,TSLB,TMN,IVGTYP,                        & !urban
                      ims,ime,jms,jme,kms,kme,nsoil,                          & !urban
                      LOW_DENSITY_RESIDENTIAL_TABLE,                                    &
		      HIGH_DENSITY_RESIDENTIAL_TABLE,                                   &
		      HIGH_INTENSITY_INDUSTRIAL_TABLE,                                  &
                      .false.,sf_urban_physics,                                          & !urban
                      XXXR_URB2D,    XXXB_URB2D,    XXXG_URB2D,XXXC_URB2D,              & !urban
                        TR_URB2D,      TB_URB2D,      TG_URB2D,  TC_URB2D, QC_URB2D,    & !urban
                       TRL_URB3D,     TBL_URB3D,     TGL_URB3D,                         & !urban
                        SH_URB2D,      LH_URB2D,       G_URB2D,  RN_URB2D, TS_URB2D,    & !urban
                   num_urban_ndm,  urban_map_zrd,  urban_map_zwd, urban_map_gd,          & !I multi-layer urban
                    urban_map_zd,  urban_map_zdf,   urban_map_bd, urban_map_wd,          & !I multi-layer urban
                   urban_map_gbd,  urban_map_fbd,                                        & !I multi-layer urban
                    num_urban_hi,                                                       & !urban
                       TRB_URB4D,     TW1_URB4D,     TW2_URB4D, TGB_URB4D,              & !urban
                      TLEV_URB3D,    QLEV_URB3D,                                        & !urban
                    TW1LEV_URB3D,  TW2LEV_URB3D,                                        & !urban
                     TGLEV_URB3D,   TFLEV_URB3D,                                        & !urban
                     SF_AC_URB3D,   LF_AC_URB3D,   CM_AC_URB3D,                         & !urban
                    SFVENT_URB3D,  LFVENT_URB3D,                                        & !urban
                    SFWIN1_URB3D,  SFWIN2_URB3D,                                        & !urban
                      SFW1_URB3D,    SFW2_URB3D,     SFR_URB3D, SFG_URB3D,              & !urban
                        LP_URB2D,      HI_URB2D,      LB_URB2D,                         & !urban
                       HGT_URB2D,      MH_URB2D,    STDH_URB2D,                         & !urban
                        LF_URB2D,                                                       & !urban
                      CMCR_URB2D,     TGR_URB2D,    TGRL_URB3D, SMR_URB3D,              & !urban
                     DRELR_URB2D,   DRELB_URB2D,   DRELG_URB2D,                         & !urban
                   FLXHUMR_URB2D, FLXHUMB_URB2D, FLXHUMG_URB2D,                         & !urban
                         A_U_BEP,       A_V_BEP,       A_T_BEP,   A_Q_BEP,              & !multi-layer urban
                         A_E_BEP,       B_U_BEP,       B_V_BEP,                         & !multi-layer urban
                         B_T_BEP,       B_Q_BEP,       B_E_BEP,   DLG_BEP,              & !multi-layer urban
                        DL_U_BEP,        SF_BEP,        VL_BEP,                         & !multi-layer urban
                       FRC_URB2D,   UTYPE_URB2D)                                          !urban
          if (debug_flag)then
          print*,'DZR',DZR
          print*,'DZB',DZB
          print*,'DZG',DZG
          print*,'ISURBAN',ISURBAN
          print*,'TSK',TSK
          print*,'TSLB',TSLB
          print*,'TMN',TMN
          print*,'IVGTYP',IVGTYP
          print*,'LOW_DENSITY_RESIDENTIAL_TABLE',LOW_DENSITY_RESIDENTIAL_TABLE
          print*,'HIGH_DENSITY_RESIDENTIAL_TABLE',HIGH_DENSITY_RESIDENTIAL_TABLE
          print*,'HIGH_INTENSITY_INDUSTRIAL_TABLE',HIGH_INTENSITY_INDUSTRIAL_TABLE
          print*,'XXXR_URB2D',XXXR_URB2D
          print*,'XXXB_URB2D',XXXB_URB2D
          print*,'XXXG_URB2D',XXXG_URB2D
          print*,'XXXC_URB2D',XXXC_URB2D
          print*,'TR_URB2D',TR_URB2D
          print*,'TB_URB2D',TB_URB2D
          print*,'TG_URB2D',TG_URB2D
          print*,'TC_URB2D',TC_URB2D
          print*,'QC_URB2D',QC_URB2D
          print*,'TRL_URB3D',TRL_URB3D
          print*,'TBL_URB3D',TBL_URB3D
          print*,'TGL_URB3D',TGL_URB3D
          print*,'SH_URB2D',SH_URB2D
          print*,'LH_URB2D',LH_URB2D
          print*,'G_URB2D',G_URB2D
          print*,'RN_URB2D',RN_URB2D
          print*,'TS_URB2D',TS_URB2D
          print*,'num_urban_ndm',num_urban_ndm
          print*,'urban_map_zrd',urban_map_zrd
          print*,'urban_map_zwd',urban_map_zwd
          print*,'urban_map_gd',urban_map_gd
          print*,'urban_map_zd',urban_map_zd
          print*,'urban_map_zdf',urban_map_zdf
          print*,'urban_map_bd',urban_map_bd
          print*,'urban_map_wd',urban_map_wd
          print*,'urban_map_gbd',urban_map_gbd
          print*,'urban_map_fbd',urban_map_fbd
          print*,'num_urban_hi',num_urban_hi
          print*,'TRB_URB4D',TRB_URB4D
          print*,'TW1_URB4D',TW1_URB4D
          print*,'TW2_URB4D',TW2_URB4D
          print*,'TGB_URB4D',TGB_URB4D
          print*,'TLEV_URB3D',TLEV_URB3D
          print*,'QLEV_URB3D',QLEV_URB3D
          print*,'TW1LEV_URB3D',TW1LEV_URB3D
          print*,'TW2LEV_URB3D',TW2LEV_URB3D
          print*,'TGLEV_URB3D',TGLEV_URB3D
          print*,'TFLEV_URB3D',TFLEV_URB3D
          print*,'SF_AC_URB3D',SF_AC_URB3D
          print*,'LF_AC_URB3D',LF_AC_URB3D
          print*,'CM_AC_URB3D',CM_AC_URB3D
          print*,'SFVENT_URB3D',SFVENT_URB3D
          print*,'LFVENT_URB3D',LFVENT_URB3D
          print*,'SFWIN1_URB3D',SFWIN1_URB3D
          print*,'SFWIN2_URB3D',SFWIN2_URB3D
          print*,'SFW1_URB3D',SFW1_URB3D
          print*,'SFW2_URB3D',SFW2_URB3D
          print*,'SFR_URB3D',SFR_URB3D
          print*,'SFG_URB3D',SFG_URB3D
          print*,'LP_URB2D',LP_URB2D
          print*,'HI_URB2D',HI_URB2D      
          print*,'LB_URB2D',LB_URB2D
          print*,'HGT_URB2D',HGT_URB2D
          print*,'MH_URB2D',MH_URB2D
          print*,'STDH_URB2D',STDH_URB2D
          print*,'LF_URB2D',LF_URB2D
          print*,'CMCR_URB2D',CMCR_URB2D
          print*,'TGR_URB2D',TGR_URB2D
          print*,'TGRL_URB3D',TGRL_URB3D
          print*,'SMR_URB3D',SMR_URB3D
          print*,'DRELR_URB2D',DRELR_URB2D
          print*,'DRELB_URB2D',DRELB_URB2D
          print*,'DRELG_URB2D',DRELG_URB2D
          print*,'FLXHUMR_URB2D',FLXHUMR_URB2D
          print*,'FLXHUMB_URB2D',FLXHUMB_URB2D
          print*,'FLXHUMG_URB2D',FLXHUMG_URB2D
          print*,'A_U_BEP',A_U_BEP
          print*,'A_V_BEP',A_V_BEP
          print*,'A_T_BEP',A_T_BEP
          print*,'A_Q_BEP',A_Q_BEP
          print*,'A_E_BEP',A_E_BEP
          print*,'B_U_BEP',B_U_BEP
          print*,'B_V_BEP',B_V_BEP
          print*,'B_T_BEP',B_T_BEP
          print*,'B_Q_BEP',B_Q_BEP
          print*,'B_E_BEP',B_E_BEP
          print*,'DLG_BEP',DLG_BEP
          print*,'DL_U_BEP',DL_U_BEP  
          print*,'SF_BEP',SF_BEP
          print*,'VL_BEP',VL_BEP
          print*,'FRC_URB2D',FRC_URB2D
          print*,'UTYPE_URB2D',UTYPE_URB2D
          endif
          ENDIF

if (debug_flag)then
print*,'OUT'
print*,'LLANDUSE ',MMINLU
print*,'SNOW ',SNOW
print*,'SNOWH ',SNOWH
print*,'CANWAT ',CANWAT
print*,'ISLTYP ',ISLTYP
print*,'IVGTYP ',IVGTYP
print*,'XLAT ',XLAT
print*,'TSLB ',TSLB
print*,'SMOIS ',SMOIS
print*,'SH2O ',SH2O
print*,'DZS ',DZS
print*,'FNDSOILW ',FNDSOILW
print*,'FNDSNOWH ',FNDSNOWH
print*,'TSK ',TSK
print*,'ISNOWXY ',ISNOWXY
print*,'TVXY ',TVXY
print*,'TGXY ',TGXY
print*,'CANICEXY ',CANICEXY
print*,'TMN ',TMN
print*,'XICE ',XICE
print*,'CANLIQXY',CANLIQXY
print*,'EAHXY',EAHXY
print*,'TAHXY',TAHXY
print*,'CMXY',CMXY
print*,'CHXY',CHXY
print*,'FWETXY',FWETXY
print*,'SNEQVOXY',SNEQVOXY
print*,'ALBOLDXY',ALBOLDXY
print*,'QSNOWXY',QSNOWXY
print*,'WSLAKEXY',WSLAKEXY
print*,'ZWTXY',ZWTXY
print*,'WAXY',WAXY
print*,'WOODXY',WOODXY
print*,'STBLCPXY',STBLCPXY
print*,'FASTCPXY',FASTCPXY
print*,'XSAIXY',XSAIXY
print*,'LAI',LAI
print*,'GRAINXY',GRAINXY
print*,'GDDXY',GDDXY
print*,'CROPTYPE',CROPTYPE
print*,'CROPCAT',CROPCAT
print*,'T2MVXY',T2MVXY
print*,'T2MBXY',T2MBXY
print*,'CHSTARXY',CHSTARXY
print*,'NSOIL',NSOIL
print*,'runoff_option',iopt_run
print*,'crop_option',iopt_crop
print*,'sf_urban_physics',sf_urban_physics
print*,'ids',ids
print*,'ide',ide
print*,'jds',jds
print*,'jde',jde
print*,'kds',kds
print*,'kde',kde
print*,'ims',ims
print*,'ime',ime
print*,'jms',jms
print*,'jme',jme
print*,'kms',kms
print*,'kme',kme
print*,'its',its
print*,'ite',ite
print*,'jts',jts
print*,'jte',jte
print*,'kts',kts
print*,'kte',kte
print*,'smoiseq',smoiseq
print*,'smcwtdxy',smcwtdxy
print*,'rechxy',rechxy
print*,'deeprechxy',deeprechxy
print*,'areaxy',areaxy
print*,'dx',dx
print*,'dy',dy
print*,'msftx',msftx
print*,'msfty',msfty
print*,'wtddt',wtddt
print*,'stepwtd',stepwtd
print*,'dt',dt
print*,'qrfsxy',qrfsxy
print*,'qspringsxy',qspringsxy
print*,'qslatxy',qslatxy
print*,'fdepthxy',fdepthxy
print*,'terrain',ht
print*,'riverbedxy',riverbedxy
print*,'eqzwt',eqzwt
print*,'rivercondxy',rivercondxy
print*,'pexpxy',pexpxy
print*,'rechclim',rechclim
print*,'gecros_state',gecros_state
endif

end subroutine


subroutine update(zlvl,dt,lwdn,swdn,u_ml,v_ml,q_ml,t_ml,prcp,psfc,&
                  nowdate,xlat,xlong,coszin,julian,nhru,&
                  ITIMESTEP,        YR,   &
                  DZS,    NSOIL,       DX,&
                  IVGTYP,   ISLTYP,    VEGFRA,   VEGMAX,      TMN, &
                  CROPCAT,&
                  PLANTING,  HARVEST,SEASON_GDD,&
                  IDVEG, IOPT_CRS,  IOPT_BTR, IOPT_RUN, IOPT_SFC, IOPT_FRZ,  &
                  IOPT_INF, IOPT_RAD,  IOPT_ALB, IOPT_SNF,IOPT_TBOT, IOPT_STC, &
                  IOPT_GLA, IOPT_RSF, IOPT_CROP,            &
                  IZ0TLND, SF_URBAN_PHYSICS,&
                  TSK,      HFX,      QFX,        LH,   GRDFLX,    SMSTAV, &
                  SMSTOT,SFCRUNOFF, UDRUNOFF,    ALBEDO,    SNOWC,     SMOIS, &
                  SH2O,     TSLB,     SNOW,     SNOWH,   CANWAT,    ACSNOM, &
                  ACSNOW,    EMISS,     QSFC,&
                  Z0,      ZNT,&
                  ISNOWXY,     TVXY,     TGXY,  CANICEXY, CANLIQXY,     EAHXY, &
                  TAHXY,     CMXY,     CHXY,    FWETXY, SNEQVOXY,  ALBOLDXY, &
                  QSNOWXY, WSLAKEXY,    ZWTXY,      WAXY,     WTXY,    TSNOXY, &
                  ZSNSOXY,  SNICEXY,  SNLIQXY,  LFMASSXY, RTMASSXY,  STMASSXY, &
                  WOODXY, STBLCPXY, FASTCPXY,    XLAIXY,   XSAIXY,   TAUSSXY, &
                  SMOISEQ, SMCWTDXY,DEEPRECHXY,   RECHXY,  GRAINXY,GDDXY,PGSXY,  &
                  GECROS_STATE,                                                &
                  T2MVXY,   T2MBXY,    Q2MVXY,   Q2MBXY,                      &
                  TRADXY,    NEEXY,    GPPXY,     NPPXY,   FVEGXY,   RUNSFXY, &
                  RUNSBXY,   ECANXY,   EDIRXY,   ETRANXY,    FSAXY,    FIRAXY, &
                  APARXY,    PSNXY,    SAVXY,     SAGXY,  RSSUNXY,   RSSHAXY, &
                  BGAPXY,   WGAPXY,    TGVXY,     TGBXY,    CHVXY,     CHBXY, &
                  SHGXY,    SHCXY,    SHBXY,     EVGXY,    EVBXY,     GHVXY, &
                  GHBXY,    IRGXY,    IRCXY,     IRBXY,     TRXY,     EVCXY, &
                  CHLEAFXY,   CHUCXY,   CHV2XY,    CHB2XY, RS, &
                  PAHXY,&
                  BEXP_3D,SMCDRY_3D,SMCWLT_3D,SMCREF_3D,SMCMAX_3D,DKSAT_3D,&
                  DWSAT_3D,PSISAT_3D,QUARTZ_3D,&
                  HDIVXY,&
                  SFCHEADRTXY,&
                  dzb,dzg,dzr,isurban,&
                  num_urban_ndm,urban_map_zrd,urban_map_zwd,&
                  urban_map_gd,urban_map_zd,urban_map_zdf,urban_map_bd,urban_map_wd,&
                  urban_map_gbd,urban_map_fbd,num_urban_hi,&
                  tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d,uc_urb2d,xxxr_urb2d,xxxb_urb2d,&
                  xxxg_urb2d,xxxc_urb2d,drelr_urb2d,drelb_urb2d,drelg_urb2d,flxhumr_urb2d,&
                  flxhumb_urb2d,flxhumg_urb2d,cmcr_urb2d,tgr_urb2d,trl_urb3d,tbl_urb3d,&
                  tgl_urb3d,tgrl_urb3d,smr_urb3d,sh_urb2d,lh_urb2d,g_urb2d,rn_urb2d,ts_urb2d,& 
                  trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,tlev_urb3d,qlev_urb3d,tw1lev_urb3d,&
                  tw2lev_urb3d,tglev_urb3d,tflev_urb3d,lf_ac_urb3d,sf_ac_urb3d,cm_ac_urb3d,&
                  sfvent_urb3d,lfvent_urb3d,sfwin1_urb3d,sfwin2_urb3d,sfw1_urb3d,sfw2_urb3d,&
                  sfr_urb3d,sfg_urb3d,hi_urb2d,lp_urb2d,lb_urb2d,hgt_urb2d,mh_urb2d,stdh_urb2d,&
                  lf_urb2d,a_u_bep,a_v_bep,a_t_bep,a_q_bep,a_e_bep,b_u_bep,b_v_bep,b_t_bep,&
                  b_q_bep,b_e_bep,vl_bep,dlg_bep,sf_bep,dl_u_bep,frc_urb2d,utype_urb2d,&
                  psim_urb2d,psih_urb2d,u10_urb2d,v10_urb2d,gz1oz0_urb2d,&
                  akms_urb2d,th2_urb2d,q2_urb2d,ust_urb2d,&
                  cmr_sfcdif,chr_sfcdif,cmc_sfcdif,chc_sfcdif,cmgr_sfcdif,chgr_sfcdif,&
                  chs,chs2,cqs2&
                 )

 USE module_sf_noahmpdrv, only: noahmplsm,noahmp_urban

 implicit none

 !I/O
 integer,intent(in) :: nhru
 real,dimension(nhru),intent(in) :: lwdn,swdn,u_ml,v_ml,q_ml,t_ml,prcp,psfc
 real,dimension(nhru),intent(in) :: xlat,xlong
 real,dimension(nhru),intent(inout) :: coszin
 real,intent(in) :: zlvl,dt
 real,intent(in) :: julian
 CHARACTER(LEN=19), INTENT(IN)  :: NOWDATE
! IN only

    INTEGER,                                         INTENT(IN   ) ::  ITIMESTEP ! timestep number
    INTEGER,                                         INTENT(IN   ) ::  YR        ! 4-digit year
    !REAL,                                            INTENT(IN   ) ::  JULIAN    ! Julian day
    !REAL,    DIMENSION(nhru), INTENT(INOUT) ::  COSZIN    ! cosine zenith angle
    !REAL,    DIMENSION(nhru), INTENT(IN   ) ::  XLAT      ! latitude [rad]
    !REAL,    DIMENSION(nhru), INTENT(IN   ) ::  XLONG     ! latitude [rad]
    !REAL,    DIMENSION(nhru,2), INTENT(IN   ) ::  DZ8W      ! thickness of atmo layers [m]
    !REAL,                                            INTENT(IN   ) ::  DT        ! timestep [s]
    REAL,    DIMENSION(1:nsoil),                     INTENT(IN   ) ::  DZS       ! thickness of soil layers [m]
    INTEGER,                                         INTENT(IN   ) ::  NSOIL     ! number of soil layers
    REAL,                                            INTENT(IN   ) ::  DX        ! horizontal grid spacing [m]
    INTEGER, DIMENSION(nhru), INTENT(IN   ) ::  IVGTYP    ! vegetation type
    INTEGER, DIMENSION(nhru), INTENT(IN   ) ::  ISLTYP    ! soil type
    REAL,    DIMENSION(nhru), INTENT(IN   ) ::  VEGFRA    ! vegetation fraction []
    REAL,    DIMENSION(nhru), INTENT(IN   ) ::  VEGMAX    ! annual max vegetation fraction []
    REAL,    DIMENSION(nhru), INTENT(IN   ) ::  TMN       ! deep soil temperature [K]
    !REAL,    DIMENSION(nhru), INTENT(IN   ) ::  XLAND     ! =2 ocean; =1 land/seaice
    !REAL,    DIMENSION(nhru), INTENT(IN   ) ::  XICE      ! fraction of grid that is seaice
    !REAL,                                            INTENT(IN   ) ::  XICE_THRES! fraction of grid determining seaice
    INTEGER,                                         INTENT(IN   ) ::  IDVEG     ! dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1      
    INTEGER,                                         INTENT(IN   ) ::  IOPT_CRS  ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_BTR  ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_RUN  ! runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_SFC  ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_FRZ  ! supercooled liquid water (1-> NY06; 2->Koren99)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_INF  ! frozen soil permeability (1-> NY06; 2->Koren99)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_RAD  ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_ALB  ! snow surface albedo (1->BATS; 2->CLASS)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_SNF  ! rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_TBOT ! lower boundary of soil temperature (1->zero-flux; 2->Noah)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_STC  ! snow/soil temperature time scheme
    INTEGER,                                         INTENT(IN   ) ::  IOPT_GLA  ! glacier option (1->phase change; 2->simple)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_RSF  ! surface resistance (1->Sakaguchi/Zeng; 2->Seller; 3->mod Sellers; 4->1+snow)
    INTEGER,                                         INTENT(IN   ) ::  IOPT_CROP ! crop model option (0->none; 1->Liu et al.; 2->Gecros)
    INTEGER,                                         INTENT(IN   ) ::  IZ0TLND   ! option of Chen adjustment of Czil (not used)
    INTEGER,                                         INTENT(IN   ) ::  sf_urban_physics   ! urban physics option
    !REAL,    DIMENSION(nhru,8), INTENT(IN   ) ::  SOILCOMP  ! soil sand and clay percentage
    !REAL,    DIMENSION(nhru), INTENT(IN   ) ::  SOILCL1   ! soil texture in layer 1
    !REAL,    DIMENSION(nhru), INTENT(IN   ) ::  SOILCL2   ! soil texture in layer 2
    !REAL,    DIMENSION(nhru), INTENT(IN   ) ::  SOILCL3   ! soil texture in layer 3
    !REAL,    DIMENSION(nhru), INTENT(IN   ) ::  SOILCL4   ! soil texture in layer 4
    !REAL,    DIMENSION(nhru,2), INTENT(IN   ) ::  T3D       ! 3D atmospheric temperature valid at mid-levels [K]
    !REAL,    DIMENSION(nhru,2), INTENT(IN   ) ::  QV3D      ! 3D water vapor mixing ratio [kg/kg_dry]
    !REAL,    DIMENSION(nhru,2), INTENT(IN   ) ::  U_PHY     ! 3D U wind component [m/s]
    !REAL,    DIMENSION(nhru,2), INTENT(IN   ) ::  V_PHY     ! 3D V wind component [m/s]
    !REAL,    DIMENSION(nhru), INTENT(IN   ) ::  SWDOWN    ! solar down at surface [W m-2]
    !REAL,    DIMENSION(nhru), INTENT(IN   ) ::  GLW       ! longwave down at surface [W m-2]
    !REAL,    DIMENSION(nhru,2), INTENT(IN   ) ::  P8W3D     ! 3D pressure, valid at interface [Pa]
    !REAL,    DIMENSION(nhru), INTENT(IN   ) ::  PRECIP_IN ! total input precipitation [mm]
    !REAL,    DIMENSION(nhru), INTENT(IN   ) ::  SR        ! frozen precipitation ratio [-]

!Optional Detailed Precipitation Partitioning Inputs
    !REAL,    DIMENSION(nhru), INTENT(IN   ), OPTIONAL ::  MP_RAINC  ! convective precipitation entering land model [mm] ! MB/AN : v3.7
    !REAL,    DIMENSION(nhru), INTENT(IN   ), OPTIONAL ::  MP_RAINNC ! large-scale precipitation entering land model [mm]! MB/AN : v3.7
    !REAL,    DIMENSION(nhru), INTENT(IN   ), OPTIONAL ::  MP_SHCV   ! shallow conv precip entering land model [mm]      ! MB/AN : v3.7
    !REAL,    DIMENSION(nhru), INTENT(IN   ), OPTIONAL ::  MP_SNOW   ! snow precipitation entering land model [mm]       ! MB/AN : v3.7
    !REAL,    DIMENSION(nhru), INTENT(IN   ), OPTIONAL ::  MP_GRAUP  ! graupel precipitation entering land model [mm]    ! MB/AN : v3.7
    !REAL,    DIMENSION(nhru), INTENT(IN   ), OPTIONAL ::  MP_HAIL   ! hail precipitation entering land model [mm]       ! MB/AN : v3.7

! Crop Model
    INTEGER, DIMENSION(nhru), INTENT(IN   ) ::  CROPCAT   ! crop catagory
    REAL,    DIMENSION(nhru), INTENT(IN   ) ::  PLANTING  ! planting date
    REAL,    DIMENSION(nhru), INTENT(IN   ) ::  HARVEST   ! harvest date
    REAL,    DIMENSION(nhru), INTENT(IN   ) ::  SEASON_GDD! growing season GDD
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  GRAINXY   ! mass of grain XING [g/m2]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  GDDXY     ! growing degree days XING (based on 10C) 
 INTEGER,    DIMENSION(nhru), INTENT(INOUT) ::  PGSXY

! gecros model
    REAL,    DIMENSION(nhru,60), INTENT(INOUT) :: gecros_state !  gecros crop

!#ifdef WRF_HYDRO
!    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  sfcheadrt,INFXSRT,soldrain   ! for WRF-Hydro
!#endif
! placeholders for 3D soil
    REAL,    DIMENSION( nhru,nsoil ), INTENT(INOUT) ::  BEXP_3D   ! C-H B exponent
    REAL,    DIMENSION( nhru,nsoil ), INTENT(INOUT) ::  SMCDRY_3D ! Soil Moisture Limit: Dry
    REAL,    DIMENSION( nhru,nsoil ), INTENT(INOUT) ::  SMCWLT_3D ! Soil Moisture Limit: Wilt
    REAL,    DIMENSION( nhru,nsoil ), INTENT(INOUT) ::  SMCREF_3D ! Soil Moisture Limit: Reference
    REAL,    DIMENSION( nhru,nsoil ), INTENT(INOUT) ::  SMCMAX_3D ! Soil Moisture Limit: Max
    REAL,    DIMENSION( nhru,nsoil ), INTENT(INOUT) ::  DKSAT_3D  ! Saturated Soil Conductivity
    REAL,    DIMENSION( nhru,nsoil ), INTENT(INOUT) ::  DWSAT_3D  ! Saturated Soil Diffusivity
    REAL,    DIMENSION( nhru,nsoil ), INTENT(INOUT) ::  PSISAT_3D ! Saturated Matric Potential
    REAL,    DIMENSION( nhru,nsoil ), INTENT(INOUT) ::  QUARTZ_3D ! Soil quartz content
    !REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  REFDK_2D  ! Reference Soil Conductivity
    !REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)          ::  REFKDT_2D ! Soil Infiltration Parameter

! INOUT (with generic LSM equivalent)

    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  TSK       ! surface radiative temperature [K]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  HFX       ! sensible heat flux [W m-2]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  QFX       ! latent heat flux [kg s-1 m-2]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  LH        ! latent heat flux [W m-2]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  GRDFLX    ! ground/snow heat flux [W m-2]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SMSTAV    ! soil moisture avail. [not used]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SMSTOT    ! total soil water [mm][not used]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SFCRUNOFF ! accumulated surface runoff [m]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  UDRUNOFF  ! accumulated sub-surface runoff [m]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  ALBEDO    ! total grid albedo []
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SNOWC     ! snow cover fraction []
    REAL,    DIMENSION( nhru, 1:nsoil), INTENT(INOUT) ::  SMOIS     ! volumetric soil moisture [m3/m3]
    REAL,    DIMENSION( nhru, 1:nsoil), INTENT(INOUT) ::  SH2O      ! volumetric liquid soil moisture [m3/m3]
    REAL,    DIMENSION( nhru, 1:nsoil), INTENT(INOUT) ::  TSLB      ! soil temperature [K]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SNOW      ! snow water equivalent [mm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SNOWH     ! physical snow depth [m]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  CANWAT    ! total canopy water + ice [mm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  ACSNOM    ! accumulated snow melt leaving pack
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  ACSNOW    ! accumulated snow on grid
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  EMISS     ! surface bulk emissivity
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  QSFC      ! bulk surface specific humidity
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  Z0        ! combined z0 sent to coupled model
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  ZNT       ! combined z0 sent to coupled model
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  RS        ! Total stomatal resistance (s/m)

    INTEGER, DIMENSION(nhru), INTENT(INOUT) ::  ISNOWXY   ! actual no. of snow layers
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  TVXY      ! vegetation leaf temperature
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  TGXY      ! bulk ground surface temperature
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  CANICEXY  ! canopy-intercepted ice (mm)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  CANLIQXY  ! canopy-intercepted liquid water (mm)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  EAHXY     ! canopy air vapor pressure (pa)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  TAHXY     ! canopy air temperature (k)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  CMXY      ! bulk momentum drag coefficient
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  CHXY      ! bulk sensible heat exchange coefficient
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  FWETXY    ! wetted or snowed fraction of the canopy (-)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SNEQVOXY  ! snow mass at last time step(mm h2o)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  ALBOLDXY  ! snow albedo at last time step (-)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  QSNOWXY   ! snowfall on the ground [mm/s]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  WSLAKEXY  ! lake water storage [mm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  ZWTXY     ! water table depth [m]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  WAXY      ! water in the "aquifer" [mm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  WTXY      ! groundwater storage [mm]
    REAL,    DIMENSION(nhru,-2:0), INTENT(INOUT) ::  TSNOXY    ! snow temperature [K]
    REAL,    DIMENSION(nhru,-2:NSOIL), INTENT(INOUT) ::  ZSNSOXY   ! snow layer depth [m]
    REAL,    DIMENSION(nhru,-2:0), INTENT(INOUT) ::  SNICEXY   ! snow layer ice [mm]
    REAL,    DIMENSION(nhru,-2:0), INTENT(INOUT) ::  SNLIQXY   ! snow layer liquid water [mm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  LFMASSXY  ! leaf mass [g/m2]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  RTMASSXY  ! mass of fine roots [g/m2]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  STMASSXY  ! stem mass [g/m2]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  WOODXY    ! mass of wood (incl. woody roots) [g/m2]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  STBLCPXY  ! stable carbon in deep soil [g/m2]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  FASTCPXY  ! short-lived carbon, shallow soil [g/m2]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  XLAIXY    ! leaf area index
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  XSAIXY    ! stem area index
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  TAUSSXY   ! snow age factor
    REAL,    DIMENSION(nhru,1:nsoil), INTENT(INOUT) ::  SMOISEQ   ! eq volumetric soil moisture [m3/m3]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SMCWTDXY  ! soil moisture content in the layer to the water table when deep
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  DEEPRECHXY ! recharge to the water table when deep
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  RECHXY    ! recharge to the water table (diagnostic) 

! OUT (with no Noah LSM equivalent)

    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  T2MVXY    ! 2m temperature of vegetation part
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  T2MBXY    ! 2m temperature of bare ground part
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  Q2MVXY    ! 2m mixing ratio of vegetation part
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  Q2MBXY    ! 2m mixing ratio of bare ground part
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  TRADXY    ! surface radiative temperature (k)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  NEEXY     ! net ecosys exchange (g/m2/s CO2)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  GPPXY     ! gross primary assimilation [g/m2/s C]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  NPPXY     ! net primary productivity [g/m2/s C]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  FVEGXY    ! Noah-MP vegetation fraction [-]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  RUNSFXY   ! surface runoff [mm/s]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  RUNSBXY   ! subsurface runoff [mm/s]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  ECANXY    ! evaporation of intercepted water (mm/s)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  EDIRXY    ! soil surface evaporation rate (mm/s]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  ETRANXY   ! transpiration rate (mm/s)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  FSAXY     ! total absorbed solar radiation (w/m2)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  FIRAXY    ! total net longwave rad (w/m2) [+ to atm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  APARXY    ! photosyn active energy by canopy (w/m2)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  PSNXY     ! total photosynthesis (umol co2/m2/s) [+]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SAVXY     ! solar rad absorbed by veg. (w/m2)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SAGXY     ! solar rad absorbed by ground (w/m2)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  RSSUNXY   ! sunlit leaf stomatal resistance (s/m)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  RSSHAXY   ! shaded leaf stomatal resistance (s/m)
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  BGAPXY    ! between gap fraction
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  WGAPXY    ! within gap fraction
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  TGVXY     ! under canopy ground temperature[K]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  TGBXY     ! bare ground temperature [K]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  CHVXY     ! sensible heat exchange coefficient vegetated
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  CHBXY     ! sensible heat exchange coefficient bare-ground
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SHGXY     ! veg ground sen. heat [w/m2]   [+ to atm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SHCXY     ! canopy sen. heat [w/m2]   [+ to atm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  SHBXY     ! bare sensible heat [w/m2]     [+ to atm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  EVGXY     ! veg ground evap. heat [w/m2]  [+ to atm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  EVBXY     ! bare soil evaporation [w/m2]  [+ to atm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  GHVXY     ! veg ground heat flux [w/m2]  [+ to soil]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  GHBXY     ! bare ground heat flux [w/m2] [+ to soil]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  IRGXY     ! veg ground net LW rad. [w/m2] [+ to atm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  IRCXY     ! canopy net LW rad. [w/m2] [+ to atm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  IRBXY     ! bare net longwave rad. [w/m2] [+ to atm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  TRXY      ! transpiration [w/m2]  [+ to atm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  EVCXY     ! canopy evaporation heat [w/m2]  [+ to atm]
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  CHLEAFXY  ! leaf exchange coefficient 
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  CHUCXY    ! under canopy exchange coefficient 
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  CHV2XY    ! veg 2m exchange coefficient 
    REAL,    DIMENSION(nhru), INTENT(INOUT) ::  CHB2XY    ! bare 2m exchange coefficient 
 !Added
    REAL,    DIMENSION(nhru) :: PAHXY
    REAL,    DIMENSION(nhru,1:nsoil), INTENT(INOUT) ::  HDIVXY     ! divergence
    REAL,    DIMENSION(nhru) :: SFCHEADRTXY
! Urban canopy model
    REAL, DIMENSION(nsoil), INTENT(INOUT) :: DZR
    REAL, DIMENSION(nsoil), INTENT(INOUT) :: DZB
    REAL, DIMENSION(nsoil), INTENT(INOUT) :: DZG
   INTEGER, INTENT(INOUT) :: ISURBAN
   INTEGER, INTENT(INOUT) :: num_urban_ndm
   INTEGER, INTENT(INOUT) :: urban_map_zrd
   INTEGER, INTENT(INOUT) :: urban_map_zwd
   INTEGER, INTENT(INOUT) :: urban_map_gd
   INTEGER, INTENT(INOUT) :: urban_map_zd
   INTEGER, INTENT(INOUT) :: urban_map_zdf
   INTEGER, INTENT(INOUT) :: urban_map_bd
   INTEGER, INTENT(INOUT) :: urban_map_wd
   INTEGER, INTENT(INOUT) :: urban_map_gbd
   INTEGER, INTENT(INOUT) :: urban_map_fbd
   INTEGER, INTENT(INOUT) :: num_urban_hi                                !multi-layer urban
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TR_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TB_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TG_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TC_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: QC_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: UC_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: XXXR_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: XXXB_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: XXXG_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: XXXC_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: DRELR_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: DRELB_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: DRELG_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: FLXHUMR_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: FLXHUMB_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: FLXHUMG_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CMCR_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TGR_URB2D
   REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT) :: TRL_URB3D
   REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT) :: TBL_URB3D
   REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT) :: TGL_URB3D
   REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT) :: TGRL_URB3D
   REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT) :: SMR_URB3D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: SH_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: LH_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: G_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: RN_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TS_URB2D
! Unknown location UCM variables
   REAL, DIMENSION(nhru), INTENT(INOUT) :: PSIM_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: PSIH_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: U10_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: V10_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: GZ1OZ0_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: AKMS_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: TH2_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: Q2_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: UST_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CMR_SFCDIF
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CHR_SFCDIF
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CMC_SFCDIF
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CHC_SFCDIF
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CMGR_SFCDIF
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CHGR_SFCDIF
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CHS
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CHS2
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CQS2
! multi-layer UCM variables
   REAL, DIMENSION(nhru,1:urban_map_zrd), INTENT(INOUT) :: TRB_URB4D
   REAL, DIMENSION(nhru,1:urban_map_zwd), INTENT(INOUT) :: TW1_URB4D
   REAL, DIMENSION(nhru,1:urban_map_zwd), INTENT(INOUT) :: TW2_URB4D
   REAL, DIMENSION(nhru,1:urban_map_gd), INTENT(INOUT) :: TGB_URB4D
   REAL, DIMENSION(nhru,1:urban_map_bd), INTENT(INOUT) :: TLEV_URB3D
   REAL, DIMENSION(nhru,1:urban_map_bd), INTENT(INOUT) :: QLEV_URB3D
   REAL, DIMENSION(nhru,1:urban_map_wd), INTENT(INOUT) :: TW1LEV_URB3D
   REAL, DIMENSION(nhru,1:urban_map_wd), INTENT(INOUT) :: TW2LEV_URB3D
   REAL, DIMENSION(nhru,1:urban_map_gbd), INTENT(INOUT) :: TGLEV_URB3D
   REAL, DIMENSION(nhru,1:urban_map_fbd), INTENT(INOUT) :: TFLEV_URB3D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: LF_AC_URB3D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: SF_AC_URB3D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: CM_AC_URB3D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: SFVENT_URB3D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: LFVENT_URB3D
   REAL, DIMENSION(nhru,1:urban_map_wd), INTENT(INOUT) :: SFWIN1_URB3D
   REAL, DIMENSION(nhru,1:urban_map_wd), INTENT(INOUT) :: SFWIN2_URB3D
   REAL, DIMENSION(nhru,1:urban_map_zd), INTENT(INOUT) :: SFW1_URB3D
   REAL, DIMENSION(nhru,1:urban_map_zd), INTENT(INOUT) :: SFW2_URB3D
   REAL, DIMENSION(nhru,1:urban_map_zdf), INTENT(INOUT) :: SFR_URB3D
   REAL, DIMENSION(nhru,1:num_urban_ndm), INTENT(INOUT) :: SFG_URB3D
   REAL, DIMENSION(nhru,1:num_urban_hi), INTENT(INOUT) :: HI_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: LP_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: LB_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: HGT_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: MH_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: STDH_URB2D
   REAL, DIMENSION(nhru,4), INTENT(INOUT) :: LF_URB2D
   REAL, DIMENSION(nhru), INTENT(INOUT) :: A_U_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: A_V_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: A_T_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: A_Q_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: A_E_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: B_U_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: B_V_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: B_T_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: B_Q_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: B_E_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: VL_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: DLG_BEP
   REAL, DIMENSION(nhru),INTENT(INOUT) :: SF_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: DL_U_BEP
   REAL, DIMENSION(nhru), INTENT(INOUT) :: FRC_URB2D
   INTEGER, DIMENSION(nhru), INTENT(INOUT) :: UTYPE_URB2D
 !Local
  integer :: i
  REAL, DIMENSION(nhru,2,1)  ::  DZ8W      ! thickness of atmo layers [m]
  REAL, DIMENSION(nhru,2,1)  ::  T_PHY     ! 3D atmospheric temperature valid at mid-levels [K]
  REAL, DIMENSION(nhru,2,1)  ::  QV_CURR   ! 3D water vapor mixing ratio [kg/kg_dry]
  REAL, DIMENSION(nhru,2,1)  ::  U_PHY     ! 3D U wind component [m/s]
  REAL, DIMENSION(nhru,2,1)  ::  V_PHY     ! 3D V wind component [m/s]
  REAL, DIMENSION(nhru,1)    ::  SWDOWN    ! solar down at surface [W m-2]
  REAL, DIMENSION(nhru,1)    ::  GLW       ! longwave down at surface [W m-2]
  REAL, DIMENSION(nhru,2,1)  ::  P8W       ! 3D pressure, valid at interface [Pa]
  REAL, DIMENSION(nhru,1)    ::  RAINBL, RAINBL_tmp    ! precipitation entering land model [mm]
  REAL, DIMENSION(nhru,1)    ::  SNOWBL    ! snow entering land model [mm]
  REAL, DIMENSION(nhru,1)    ::  SR        ! frozen precip ratio entering land model [-]
  REAL, DIMENSION(nhru,1)    ::  RAINCV    ! convective precip forcing [mm]
  REAL, DIMENSION(nhru,1)    ::  RAINNCV   ! non-convective precip forcing [mm]
  REAL, DIMENSION(nhru,1)    ::  RAINSHV   ! shallow conv. precip forcing [mm]
  REAL, DIMENSION(nhru,1)    ::  SNOWNCV   ! non-covective snow forcing (subset of rainncv) [mm]
  REAL, DIMENSION(nhru,1)    ::  GRAUPELNCV! non-convective graupel forcing (subset of rainncv) [mm]
  REAL, DIMENSION(nhru,1)    ::  HAILNCV   ! non-convective hail forcing (subset of rainncv) [mm]
  REAL, DIMENSION(nhru) :: HRANG
  REAL :: DECLIN,GMT
  INTEGER :: JULDAY
  REAL,    DIMENSION(nhru) ::  XLAND     ! =2 ocean; =1 land/seaice
  REAL,    DIMENSION(nhru) ::  XICE      ! fraction of grid that is seaice
  REAL ::XICE_THRES! fraction of grid determining seaice
integer :: ids,ide, jds,jde, kds,kde, ims,ime, jms,jme, kms,kme, its,ite
integer :: jts,jte, kts,kte
logical :: debug_flag
  INTEGER ::  IOPT_SOIL ! soil configuration option (Fixed for HB, parameters set outside)
  INTEGER ::  IOPT_PEDO ! soil pedotransfer function option (Fixed for HB, parameters set outside)
  REAL,    DIMENSION(nhru,8) ::  SOILCOMP  ! soil sand and clay percentage
  REAL,    DIMENSION(nhru) ::  SOILCL1   ! soil texture in layer 1
  REAL,    DIMENSION(nhru) ::  SOILCL2   ! soil texture in layer 2
  REAL,    DIMENSION(nhru) ::  SOILCL3   ! soil texture in layer 3
  REAL,    DIMENSION(nhru) ::  SOILCL4   ! soil texture in layer 4
  !UCM variables
   INTEGER ::  num_roof_layers = 4
   INTEGER ::  num_road_layers = 4
   INTEGER ::  num_wall_layers = 4
  REAL,DIMENSION(nhru) :: theta_urban,u_urban,v_urban,dz_urban,rho_urban,p_urban,ust
!Hard coded for HB (on purpose)
IOPT_SOIL = 1
IOPT_PEDO = 1

debug_flag = .False.
ids = 1
ide = nhru+1
jds = 1
jde = 2
kds = 1
kde = 2
ims = 1
ime = nhru
jms = 1
jme = 1
kms = 1
kme = 2
its = 1
ite = nhru
jts = 1
jte = 1
kts = 1
kte = 2
  XLAND(:) = 1
  XICE(:) = 0
  XICE_THRES = 0.5

 !Prepare the forcings
 QV_CURR(:,1,1) = q_ml/(1.0 - q_ml)  ! Assuming input forcing are specific hum.;WRF wants mixing ratio at driver level
 P8W(:,1,1) = psfc
 P8W(:,2,1)     = P8W(:,1,1)      ! WRF uses lowest two layers
 T_PHY(:,1,1) = t_ml
 T_PHY(:,2,1)   = T_PHY(:,1,1)    ! Only pressure is needed in two layer but fill the rest
 U_PHY(:,1,1) = u_ml
 U_PHY(:,2,1)   = U_PHY(:,1,1)    ! 
 V_PHY(:,1,1) = v_ml
 V_PHY(:,2,1)   = V_PHY(:,1,1)    ! 
 QV_CURR(:,2,1) = QV_CURR(:,1,1)  ! 
 RAINBL_tmp(:,1) = prcp
 RAINBL = RAINBL_tmp * DT       ! RAINBL in WRF is [mm]
 SNOWBL = 0.0 ! Could add if info available
 !SNOWBL = SNOWBL * DTBL           ! 
 SR         = 0.0                 ! Will only use component if opt_snf=4
 RAINCV     = 0.0
 RAINNCV    = RAINBL
 RAINSHV    = 0.0
 SNOWNCV    = SNOWBL
 GRAUPELNCV = 0.0
 HAILNCV    = 0.0
 DZ8W = 2*ZLVL                    ! 2* to be consistent with WRF model level
 SWDOWN(:,1) = swdn
 GLW(:,1) = lwdn

   ! create a few fields that are IN in WRF - coszen, julian,yr

    do i = 1,nhru
      !IF(SF_URBAN_PHYSICS == 0) THEN
      !  CALL CALC_DECLIN(NOWDATE,XLAT(I), XLONG(I),COSZEN(I),JULIAN)
      !ELSE
        CALL CALC_DECLIN(NOWDATE,XLAT(I), XLONG(I),COSZIN(I),JULIAN, &
                         HRANG(I), DECLIN, GMT, JULDAY)
      !ENDIF
    enddo

if (debug_flag)then
    print*,'NOAHMPLSM IN'
    print*,'ITIMESTEP',ITIMESTEP
    print*,'YR',YR
    print*,'JULIAN',JULIAN
    print*,'COSZEN',COSZIN
    print*,'XLAT',XLAT
    print*,'XLONG',XLONG
    print*,'DZ8W',DZ8W
    print*,'DTBL',DT
    print*,'DZS',DZS
    print*,'NSOIL',NSOIL
    print*,'DX',DX
    print*,'IVGTYP',IVGTYP
    print*,'ISLTYP',ISLTYP
    print*,'VEGFRA',VEGFRA
    print*,'GVFMAX',VEGMAX
    print*,'TMN',TMN
    print*,'XLAND',XLAND
    print*,'XICE',XICE
    !print*,'XICE_THRESHOLD',XICE_THRESHOLD
    print*,'CROPCAT',CROPCAT
    print*,'PLANTING',PLANTING
    print*,'HARVEST',HARVEST
    print*,'SEASON_GDD',SEASON_GDD
    print*,'IDVEG',IDVEG
    print*,'IOPT_CRS',IOPT_CRS
    print*,'IOPT_BTR',IOPT_RUN
    print*,'IOPT_SFC',IOPT_SFC
    print*,'IOPT_FRZ',IOPT_FRZ
    print*,'IOPT_INF',IOPT_INF
    print*,'IOPT_RAD',IOPT_RAD
    print*,'IOPT_ALB',IOPT_ALB
    print*,'IOPT_SNF',IOPT_SNF
    print*,'IOPT_TBOT',IOPT_TBOT
    print*,'IOPT_STC',IOPT_STC
    print*,'IOPT_GLA',IOPT_GLA
    print*,'IOPT_RSF',IOPT_RSF
    !print*,'IOPT_SOIL',IOPT_SOIL
    !print*,'IOPT_PEDO',IOPT_PEDO
    print*,'IOPT_CROP',IOPT_CROP
    print*,'IZ0TLND',IZ0TLND
    print*,'sf_urban_physics',sf_urban_physics
    !print*,'SOILCOMP',SOILCOMP
    !print*,'SOILCL1',SOILCL1
    !print*,'SOILCL2',SOILCL2
    !print*,'SOILCL3',SOILCL3
    !print*,'SOILCL4',SOILCL4
    print*,'T_PHY',T_PHY
    print*,'QV_CURR',QV_CURR
    print*,'U_PHY',U_PHY
    print*,'V_PHY',V_PHY
    print*,'SWDOWN',SWDOWN
    print*,'GLW',GLW
    print*,'P8W',P8W
    print*,'RAINBL',RAINBL
    print*,'SR',SR
    print*,'TSK',TSK
    print*,'HFX',HFX
    print*,'QFX',QFX
    print*,'LH',LH
    print*,'GRDFLX',GRDFLX
    print*,'SMSTAV',SMSTAV
    print*,'SMSTOT',SMSTOT
    print*,'SFCRUNOFF',SFCRUNOFF
    print*,'UDRUNOFF',UDRUNOFF
    print*,'ALBEDO',ALBEDO
    print*,'SNOWC',SNOWC
    print*,'SMOIS',SMOIS
    print*,'SH2O',SH2O
    print*,'TSLB',TSLB
    print*,'SNOW',SNOW
    print*,'SNOWH',SNOWH
    print*,'CANWAT',CANWAT
    print*,'ACSNOM',ACSNOM
    print*,'ACSNOW',ACSNOW
    print*,'EMISS',EMISS
    print*,'QSFC',QSFC
    print*,'Z0',Z0
    print*,'ZNT',ZNT
    print*,'ISNOWXY',ISNOWXY
    print*,'TVXY',TVXY
    print*,'TGXY',TGXY
    print*,'CANICEXY',CANICEXY
    print*,'CANLIQXY',CANLIQXY
    print*,'EAHXY',EAHXY
    print*,'TAHXY',TAHXY
    print*,'CMXY',CMXY
    print*,'CHXY',CHXY
    print*,'FWETXY',FWETXY
    print*,'SNEQVOXY',SNEQVOXY
    print*,'ALBOLDXY',ALBOLDXY
    print*,'QSNOWXY',QSNOWXY
    print*,'WSLAKEXY',WSLAKEXY
    print*,'WAXY',WAXY
    print*,'WTXY',WTXY
    print*,'TSNOXY',TSNOXY
    print*,'ZSNSOXY',ZSNSOXY
    print*,'SNICEXY',SNICEXY
    print*,'SNLIQXY',SNLIQXY
    print*,'LFMASSXY',LFMASSXY
    print*,'RTMASSXY',RTMASSXY
    print*,'STMASSXY',STMASSXY
    print*,'WOODXY',WOODXY
    print*,'STBLCPXY',STBLCPXY
    print*,'FASTCPXY',FASTCPXY
    print*,'XLAIXY',XLAIXY
    print*,'XSAIXY',XSAIXY
    print*,'TAUSSXY',TAUSSXY
    print*,'SMOISEQ',SMOISEQ
    print*,'SMCWTDXY',SMCWTDXY
    print*,'DEEPRECHXY',DEEPRECHXY
    print*,'RECHXY',RECHXY
    print*,'GRAINXY',GRAINXY
    print*,'GDDXY',GDDXY
    print*,'PGSXY',PGSXY
    print*,'GECROS_STATE',GECROS_STATE
    print*,'T2MVXY',T2MVXY
    print*,'T2MBXY',T2MBXY
    print*,'Q2MVXY',Q2MVXY
    print*,'Q2MBXY',Q2MBXY
    print*,'TRADXY',TRADXY
    print*,'NEEXY',NEEXY
    print*,'GPPXY',GPPXY
    print*,'NPPXY',NPPXY
    print*,'FVEGXY',FVEGXY
    print*,'RUNSFXY',RUNSFXY
    print*,'RUNSBXY',RUNSBXY
    print*,'ECANXY',ECANXY
    print*,'EDIRXY',EDIRXY
    print*,'ETRANXY',ETRANXY
    print*,'FSAXY',FSAXY
    print*,'FIRAXY',FIRAXY
    print*,'APARXY',APARXY
    print*,'PSNXY',PSNXY
    print*,'SAVXY',SAVXY
    print*,'SAGXY',SAGXY
    print*,'RSSUNXY',RSSUNXY
    print*,'RSSHAXY',RSSHAXY
    print*,'BGAPXY',BGAPXY
    print*,'WGAPXY',WGAPXY
    print*,'TGVXY',TGVXY
    print*,'TGBXY',TGBXY
    print*,'CHVXY',CHVXY
    print*,'CHBXY',CHBXY
    print*,'SHGXY',SHGXY
    print*,'SHCXY',SHCXY
    print*,'SHBXY',SHBXY
    print*,'EVGXY',EVGXY
    print*,'EVBXY',EVBXY
    print*,'GHVXY',GHVXY
    print*,'GHBXY',GHBXY
    print*,'IRGXY',IRGXY
    print*,'IRCXY',IRCXY
    print*,'IRBXY',IRBXY
    print*,'TRXY',TRXY
    print*,'EVCXY',EVCXY
    print*,'CHLEAFXY',CHLEAFXY
    print*,'CHUCXY',CHUCXY
    print*,'CHV2XY',CHV2XY
    print*,'CHB2XY',CHB2XY
    print*,'RS',RS
    print*,'BEXP_3D',BEXP_3D
    print*,'SMCDRY_3D',SMCDRY_3D
    print*,'SMCWLT_3D',SMCWLT_3D
    print*,'SMCREF_3D',SMCREF_3D
    print*,'SMCMAX_3D',SMCMAX_3D
    print*,'DKSAT_3D',DKSAT_3D
    print*,'DWSAT_3D',DWSAT_3D
    print*,'PSISAT_3D',PSISAT_3D
    print*,'QUARTZ_3D',QUARTZ_3D
    stop
endif
    !IF (ITIME == 1 .AND. .NOT. RESTART_FLAG ) THEN
    !  EAHXY = (P8W(:,1,:)*QV_CURR(:,1,:))/(0.622+QV_CURR(:,1,:)) ! Initial guess only.
    !  TAHXY = T_PHY(:,1,:)                                       ! Initial guess only.
    !  CHXY = 0.1
    !  CMXY = 0.1
    !ENDIF

    !print*,'BEXP_3D',shape(BEXP_3D)
    !print*,'SMCDRY_3D',SMCDRY_3D
    !print*,'SMCWLT_3D',SMCWLT_3D
    !print*,'SMCREF_3D',SMCREF_3D
    !print*,'SMCMAX_3D',SMCMAX_3D
    !print*,'DKSAT_3D',DKSAT_3D
    !print*,'DWSAT_3D',DWSAT_3D
    !print*,'PSISAT_3D',PSISAT_3D
    !print*,'QUARTZ_3D',QUARTZ_3D

  call noahmplsm(ITIMESTEP,        YR,   JULIAN,   COSZIN,XLAT,XLONG, & ! IN : Time/Space-related
                  DZ8W,       DT,       DZS,    NSOIL,       DX,            & ! IN : Model configuration 
	        IVGTYP,   ISLTYP,    VEGFRA,   VEGMAX,      TMN,            & ! IN : Vegetation/Soil characteristics
		 XLAND,     XICE,XICE_THRES,  CROPCAT,                      & ! IN : Vegetation/Soil characteristics
	       PLANTING,  HARVEST,SEASON_GDD,                               &
                 IDVEG, IOPT_CRS,  IOPT_BTR, IOPT_RUN, IOPT_SFC, IOPT_FRZ,  & ! IN : User options
              IOPT_INF, IOPT_RAD,  IOPT_ALB, IOPT_SNF,IOPT_TBOT, IOPT_STC,  & ! IN : User options
              IOPT_GLA, IOPT_RSF, IOPT_SOIL,IOPT_PEDO,IOPT_CROP,            & ! IN : User options
              IZ0TLND, SF_URBAN_PHYSICS,                                    & ! IN : User options
	      SOILCOMP,  SOILCL1,  SOILCL2,   SOILCL3,  SOILCL4,            & ! IN : User options
                   T_PHY,     QV_CURR,     U_PHY,    V_PHY,   SWDOWN,      GLW,  & ! IN : Forcing
		 P8W,RAINBL,        SR,                                & ! IN : Forcing
                   TSK,      HFX,      QFX,        LH,   GRDFLX,    SMSTAV, & ! IN/OUT LSM eqv
                SMSTOT,SFCRUNOFF, UDRUNOFF,    ALBEDO,    SNOWC,     SMOIS, & ! IN/OUT LSM eqv
		  SH2O,     TSLB,     SNOW,     SNOWH,   CANWAT,    ACSNOM, & ! IN/OUT LSM eqv
		ACSNOW,    EMISS,     QSFC,                                 & ! IN/OUT LSM eqv
 		    Z0,      ZNT,                                           & ! IN/OUT LSM eqv
               ISNOWXY,     TVXY,     TGXY,  CANICEXY, CANLIQXY,     EAHXY, & ! IN/OUT Noah MP only
	         TAHXY,     CMXY,     CHXY,    FWETXY, SNEQVOXY,  ALBOLDXY, & ! IN/OUT Noah MP only
               QSNOWXY, WSLAKEXY,    ZWTXY,      WAXY,     WTXY,    TSNOXY, & ! IN/OUT Noah MP only
	       ZSNSOXY,  SNICEXY,  SNLIQXY,  LFMASSXY, RTMASSXY,  STMASSXY, & ! IN/OUT Noah MP only
	        WOODXY, STBLCPXY, FASTCPXY,    XLAIXY,   XSAIXY,   TAUSSXY, & ! IN/OUT Noah MP only
	       SMOISEQ, SMCWTDXY,DEEPRECHXY,   RECHXY,  GRAINXY,    GDDXY,PGSXY,  & ! IN/OUT Noah MP only
               GECROS_STATE,                                                & ! IN/OUT gecros model
	        T2MVXY,   T2MBXY,    Q2MVXY,   Q2MBXY,                      & ! OUT Noah MP only
	        TRADXY,    NEEXY,    GPPXY,     NPPXY,   FVEGXY,   RUNSFXY, & ! OUT Noah MP only
	       RUNSBXY,   ECANXY,   EDIRXY,   ETRANXY,    FSAXY,    FIRAXY, & ! OUT Noah MP only
	        APARXY,    PSNXY,    SAVXY,     SAGXY,  RSSUNXY,   RSSHAXY, & ! OUT Noah MP only
		BGAPXY,   WGAPXY,    TGVXY,     TGBXY,    CHVXY,     CHBXY, & ! OUT Noah MP only
		 SHGXY,    SHCXY,    SHBXY,     EVGXY,    EVBXY,     GHVXY, & ! OUT Noah MP only
		 GHBXY,    IRGXY,    IRCXY,     IRBXY,     TRXY,     EVCXY, & ! OUT Noah MP only
              CHLEAFXY,   CHUCXY,   CHV2XY,    CHB2XY, RS,                  & ! OUT Noah MP only
                 BEXP_3D,SMCDRY_3D,SMCWLT_3D,SMCREF_3D,SMCMAX_3D,          & ! placeholders to activate 3D soil
		 DKSAT_3D,DWSAT_3D,PSISAT_3D,QUARTZ_3D,                     &
!		 REFDK_2D,REFKDT_2D,                                        &
!#ifdef WRF_HYDRO
!               sfcheadrt,INFXSRT,soldrain,                                  &
!#endif
               ids,ide,  jds,jde,  kds,kde,                    &
               ims,ime,  jms,jme,  kms,kme,                    &
               its,ite,  jts,jte,  kts,kte,                    &
!               MP_RAINC, MP_RAINNC, MP_SHCV, MP_SNOW, MP_GRAUP, MP_HAIL     )
               PAHXY,HDIVXY,SFCHEADRTXY)

if (debug_flag)then
    print*,'NOAHMPLSM OUT'
    print*,'ITIMESTEP',ITIMESTEP
    print*,'YR',YR
    print*,'JULIAN',JULIAN
    print*,'COSZEN',COSZIN
    print*,'XLAT',XLAT
    print*,'XLONG',XLONG
    print*,'DZ8W',DZ8W
    print*,'DTBL',DT
    print*,'DZS',DZS
    print*,'NSOIL',NSOIL
    print*,'DX',DX
    print*,'IVGTYP',IVGTYP
    print*,'ISLTYP',ISLTYP
    print*,'VEGFRA',VEGFRA
    print*,'GVFMAX',VEGMAX
    print*,'TMN',TMN
    print*,'XLAND',XLAND
    print*,'XICE',XICE
    !print*,'XICE_THRESHOLD',XICE_THRESHOLD
    print*,'CROPCAT',CROPCAT
    print*,'PLANTING',PLANTING
    print*,'HARVEST',HARVEST
    print*,'SEASON_GDD',SEASON_GDD
    print*,'IDVEG',IDVEG
    print*,'IOPT_CRS',IOPT_CRS
    print*,'IOPT_BTR',IOPT_RUN
    print*,'IOPT_SFC',IOPT_SFC
    print*,'IOPT_FRZ',IOPT_FRZ
    print*,'IOPT_INF',IOPT_INF
    print*,'IOPT_RAD',IOPT_RAD
    print*,'IOPT_ALB',IOPT_ALB
    print*,'IOPT_SNF',IOPT_SNF
    print*,'IOPT_TBOT',IOPT_TBOT
    print*,'IOPT_STC',IOPT_STC
    print*,'IOPT_GLA',IOPT_GLA
    print*,'IOPT_RSF',IOPT_RSF
    !print*,'IOPT_SOIL',IOPT_SOIL
    !print*,'IOPT_PEDO',IOPT_PEDO
    print*,'IOPT_CROP',IOPT_CROP
    print*,'IZ0TLND',IZ0TLND
    print*,'sf_urban_physics',sf_urban_physics
    !print*,'SOILCOMP',SOILCOMP
    !print*,'SOILCL1',SOILCL1
    !print*,'SOILCL2',SOILCL2
    !print*,'SOILCL3',SOILCL3
    !print*,'SOILCL4',SOILCL4
    print*,'T_PHY',T_PHY
    print*,'QV_CURR',QV_CURR
    print*,'U_PHY',U_PHY
    print*,'V_PHY',V_PHY
    print*,'SWDOWN',SWDOWN
    print*,'GLW',GLW
    print*,'P8W',P8W
    print*,'RAINBL',RAINBL
    print*,'SR',SR
    print*,'TSK',TSK
    print*,'HFX',HFX
    print*,'QFX',QFX
    print*,'LH',LH
    print*,'GRDFLX',GRDFLX
    print*,'SMSTAV',SMSTAV
    print*,'SMSTOT',SMSTOT
    print*,'SFCRUNOFF',SFCRUNOFF
    print*,'UDRUNOFF',UDRUNOFF
    print*,'ALBEDO',ALBEDO
    print*,'SNOWC',SNOWC
    print*,'SMOIS',SMOIS
    print*,'SH2O',SH2O
    print*,'TSLB',TSLB
    print*,'SNOW',SNOW
    print*,'SNOWH',SNOWH
    print*,'CANWAT',CANWAT
    print*,'ACSNOM',ACSNOM
    print*,'ACSNOW',ACSNOW
    print*,'EMISS',EMISS
    print*,'QSFC',QSFC
    print*,'Z0',Z0
    print*,'ZNT',ZNT
    print*,'ISNOWXY',ISNOWXY
    print*,'TVXY',TVXY
    print*,'TGXY',TGXY
    print*,'CANICEXY',CANICEXY
    print*,'CANLIQXY',CANLIQXY
    print*,'EAHXY',EAHXY
    print*,'TAHXY',TAHXY
    print*,'CMXY',CMXY
    print*,'CHXY',CHXY
    print*,'FWETXY',FWETXY
    print*,'SNEQVOXY',SNEQVOXY
    print*,'ALBOLDXY',ALBOLDXY
    print*,'QSNOWXY',QSNOWXY
    print*,'WSLAKEXY',WSLAKEXY
    print*,'WAXY',WAXY
    print*,'WTXY',WTXY
    print*,'TSNOXY',TSNOXY
    print*,'ZSNSOXY',ZSNSOXY
    print*,'SNICEXY',SNICEXY
    print*,'SNLIQXY',SNLIQXY
    print*,'LFMASSXY',LFMASSXY
    print*,'RTMASSXY',RTMASSXY
    print*,'STMASSXY',STMASSXY
    print*,'WOODXY',WOODXY
    print*,'STBLCPXY',STBLCPXY
    print*,'FASTCPXY',FASTCPXY
    print*,'XLAIXY',XLAIXY
    print*,'XSAIXY',XSAIXY
    print*,'TAUSSXY',TAUSSXY
    print*,'SMOISEQ',SMOISEQ
    print*,'SMCWTDXY',SMCWTDXY
    print*,'DEEPRECHXY',DEEPRECHXY
    print*,'RECHXY',RECHXY
    print*,'GRAINXY',GRAINXY
    print*,'GDDXY',GDDXY
    print*,'PGSXY',PGSXY
    print*,'GECROS_STATE',GECROS_STATE
    print*,'T2MVXY',T2MVXY
    print*,'T2MBXY',T2MBXY
    print*,'Q2MVXY',Q2MVXY
    print*,'Q2MBXY',Q2MBXY
    print*,'TRADXY',TRADXY
    print*,'NEEXY',NEEXY
    print*,'GPPXY',GPPXY
    print*,'NPPXY',NPPXY
    print*,'FVEGXY',FVEGXY
    print*,'RUNSFXY',RUNSFXY
    print*,'RUNSBXY',RUNSBXY
    print*,'ECANXY',ECANXY
    print*,'EDIRXY',EDIRXY
    print*,'ETRANXY',ETRANXY
    print*,'FSAXY',FSAXY
    print*,'FIRAXY',FIRAXY
    print*,'APARXY',APARXY
    print*,'PSNXY',PSNXY
    print*,'SAVXY',SAVXY
    print*,'SAGXY',SAGXY
    print*,'RSSUNXY',RSSUNXY
    print*,'RSSHAXY',RSSHAXY
    print*,'BGAPXY',BGAPXY
    print*,'WGAPXY',WGAPXY
    print*,'TGVXY',TGVXY
    print*,'TGBXY',TGBXY
    print*,'CHVXY',CHVXY
    print*,'CHBXY',CHBXY
    print*,'SHGXY',SHGXY
    print*,'SHCXY',SHCXY
    print*,'SHBXY',SHBXY
    print*,'EVGXY',EVGXY
    print*,'EVBXY',EVBXY
    print*,'GHVXY',GHVXY
    print*,'GHBXY',GHBXY
    print*,'IRGXY',IRGXY
    print*,'IRCXY',IRCXY
    print*,'IRBXY',IRBXY
    print*,'TRXY',TRXY
    print*,'EVCXY',EVCXY
    print*,'CHLEAFXY',CHLEAFXY
    print*,'CHUCXY',CHUCXY
    print*,'CHV2XY',CHV2XY
    print*,'CHB2XY',CHB2XY
    print*,'RS',RS
    print*,'BEXP_3D',BEXP_3D
    print*,'SMCDRY_3D',SMCDRY_3D
    print*,'SMCWLT_3D',SMCWLT_3D
    print*,'SMCREF_3D',SMCREF_3D
    print*,'SMCMAX_3D',SMCMAX_3D
    print*,'DKSAT_3D',DKSAT_3D
    print*,'DWSAT_3D',DWSAT_3D
    print*,'PSISAT_3D',PSISAT_3D
    print*,'QUARTZ_3D',QUARTZ_3D
    stop
endif

if (debug_flag)then
           print*,'URBAN IN'
           print*,'sf_urban_physics',sf_urban_physics 
           print*,'nsoil',nsoil
           print*,'ivgtyp',ivgtyp
           print*,'itimestep',itimestep
           print*,'dtbl',dt
           print*,'coszen',coszin
           print*,'xlat',xlat
           print*,'t_phy',t_phy
           print*,'qv_curr',qv_curr
           print*,'u_phy',u_phy
           print*,'v_phy',v_phy
           print*,'swdown',swdown
           print*,'glw',glw
           print*,'p8w',p8w
           print*,'rainbl',rainbl
           print*,'dz8w',dz8w
           print*,'znt',znt
           print*,'tsk',tsk
           print*,'hfx',hfx
           print*,'qfx',qfx
           print*,'lh',lh
           print*,'grdflx',grdflx
           print*,'albedo',albedo
           print*,'emiss',emiss
           print*,'qsfc',qsfc
           print*,'ids',ids
           print*,'ide',ide
           print*,'jds',jds
           print*,'jde',jde
           print*,'kds',kds
           print*,'kde',kde
           print*,'ims',ims
           print*,'ime',ime
           print*,'jms',jms
           print*,'jme',jme
           print*,'kms',kms
           print*,'kme',kme
           print*,'its',its
           print*,'ite',ite
           print*,'jts',jts
           print*,'jte',jte
           print*,'kts',kts
           print*,'kte',kte
           print*,'cmr_sfcdif',cmr_sfcdif
           print*,'chr_sfcdif',chr_sfcdif
           print*,'cmc_sfcdif',cmc_sfcdif
           print*,'chc_sfcdif',chc_sfcdif
           print*,'cmgr_sfcdif',cmgr_sfcdif
           print*,'chgr_sfcdif',chgr_sfcdif
          print*,'TR_URB2D',TR_URB2D
          print*,'TB_URB2D',TB_URB2D
          print*,'TG_URB2D',TG_URB2D
          print*,'TC_URB2D',TC_URB2D
          print*,'QC_URB2D',QC_URB2D
          print*,'UC_URB2D',UC_URB2D
          print*,'XXXR_URB2D',XXXR_URB2D
          print*,'XXXB_URB2D',XXXB_URB2D
          print*,'XXXG_URB2D',XXXG_URB2D
          print*,'XXXC_URB2D',XXXC_URB2D
          print*,'TRL_URB3D',TRL_URB3D
          print*,'TBL_URB3D',TBL_URB3D
          print*,'TGL_URB3D',TGL_URB3D
          print*,'SH_URB2D',SH_URB2D
          print*,'LH_URB2D',LH_URB2D
          print*,'G_URB2D',G_URB2D
          print*,'RN_URB2D',RN_URB2D
          print*,'TS_URB2D',TS_URB2D
          print*,'PSIM_URB2D',PSIM_URB2D
          print*,'PSIH_URB2D',PSIH_URB2D
          print*,'U10_URB2D',U10_URB2D
          print*,'V10_URB2D',V10_URB2D
          print*,'GZ1OZ0_URB2D',GZ1OZ0_URB2D
          print*,'AKMS_URB2D',AKMS_URB2D
          print*,'TH2_URB2D',TH2_URB2D
          print*,'Q2_URB2D',Q2_URB2D
          print*,'UST_URB2D',UST_URB2D
          print*,'declin',declin
          print*,'hrang',hrang
          print*,'num_roof_layers',num_roof_layers
          print*,'num_wall_layers',num_wall_layers
          print*,'num_road_layers',num_road_layers
          print*,'DZR',DZR
          print*,'DZB',DZB
          print*,'DZG',DZG
          print*,'CMCR_URB2D',CMCR_URB2D
          print*,'TGR_URB2D',TGR_URB2D
          print*,'TGRL_URB3D',TGRL_URB3D
          print*,'SMR_URB3D',SMR_URB3D
          print*,'DRELR_URB2D',DRELR_URB2D
          print*,'DRELB_URB2D',DRELB_URB2D
          print*,'DRELG_URB2D',DRELG_URB2D
          print*,'FLXHUMR_URB2D',FLXHUMR_URB2D
          print*,'FLXHUMB_URB2D',FLXHUMB_URB2D
          print*,'FLXHUMG_URB2D',FLXHUMG_URB2D
          print*,'julday',julday
          print*,'yr',yr
          print*,'FRC_URB2D',FRC_URB2D
          print*,'UTYPE_URB2D',UTYPE_URB2D
          print*,'chs',chs
          print*,'chs2',chs2
          print*,'cqs2',cqs2
          print*,'num_urban_ndm',num_urban_ndm
          print*,'urban_map_zrd',urban_map_zrd
          print*,'urban_map_zwd',urban_map_zwd
          print*,'urban_map_gd',urban_map_gd
          print*,'urban_map_zd',urban_map_zd
          print*,'urban_map_zdf',urban_map_zdf
          print*,'urban_map_bd',urban_map_bd
          print*,'urban_map_wd',urban_map_wd
          print*,'urban_map_gbd',urban_map_gbd
          print*,'urban_map_fbd',urban_map_fbd
          print*,'num_urban_hi',num_urban_hi
          print*,'TRB_URB4D',TRB_URB4D
          print*,'TW1_URB4D',TW1_URB4D
          print*,'TW2_URB4D',TW2_URB4D
          print*,'TGB_URB4D',TGB_URB4D
          print*,'TLEV_URB3D',TLEV_URB3D
          print*,'QLEV_URB3D',QLEV_URB3D
          print*,'TW1LEV_URB3D',TW1LEV_URB3D
          print*,'TW2LEV_URB3D',TW2LEV_URB3D
          print*,'TGLEV_URB3D',TGLEV_URB3D
          print*,'TFLEV_URB3D',TFLEV_URB3D
          print*,'SF_AC_URB3D',SF_AC_URB3D
          print*,'LF_AC_URB3D',LF_AC_URB3D
          print*,'CM_AC_URB3D',CM_AC_URB3D
          print*,'SFVENT_URB3D',SFVENT_URB3D
          print*,'LFVENT_URB3D',LFVENT_URB3D
          print*,'SFWIN1_URB3D',SFWIN1_URB3D
          print*,'SFWIN2_URB3D',SFWIN2_URB3D
          print*,'SFW1_URB3D',SFW1_URB3D
          print*,'SFW2_URB3D',SFW2_URB3D
          print*,'SFR_URB3D',SFR_URB3D
          print*,'SFG_URB3D',SFG_URB3D
          print*,'LP_URB2D',LP_URB2D
          print*,'HI_URB2D',HI_URB2D      
          print*,'LB_URB2D',LB_URB2D
          print*,'HGT_URB2D',HGT_URB2D
          print*,'MH_URB2D',MH_URB2D
          print*,'STDH_URB2D',STDH_URB2D
          print*,'LF_URB2D',LF_URB2D
          print*,'theta_urban',theta_urban
          print*,'rho_urban',rho_urban
          print*,'p_urban',p_urban
          print*,'ust',ust
          print*,'gmt',gmt
          print*,'julday',julday
          print*,'xlong',xlong
          print*,'xlat',xlat
          print*,'A_U_BEP',A_U_BEP
          print*,'A_V_BEP',A_V_BEP
          print*,'A_T_BEP',A_T_BEP
          print*,'A_Q_BEP',A_Q_BEP
          print*,'A_E_BEP',A_E_BEP
          print*,'B_U_BEP',B_U_BEP
          print*,'B_V_BEP',B_V_BEP
          print*,'B_T_BEP',B_T_BEP
          print*,'B_Q_BEP',B_Q_BEP
          print*,'B_E_BEP',B_E_BEP
          print*,'DLG_BEP',DLG_BEP
          print*,'DL_U_BEP',DL_U_BEP  
          print*,'SF_BEP',SF_BEP
          print*,'VL_BEP',VL_BEP
endif
	   call noahmp_urban (sf_urban_physics,     NSOIL,     IVGTYP,ITIMESTEP,  & ! IN : Model configuration 
                               DT,         COSZIN,           XLAT,                        & ! IN : Time/Space-related
                              T_PHY,        QV_CURR,          U_PHY,      V_PHY,   SWDOWN,  & ! IN : Forcing
		                GLW,            P8W,         RAINBL,       DZ8W,      ZNT,  & ! IN : Forcing
                                TSK,            HFX,            QFX,         LH,   GRDFLX,  & ! IN/OUT : LSM 
		             ALBEDO,          EMISS,           QSFC,                        & ! IN/OUT : LSM 
                            ids,ide,        jds,jde,        kds,kde,                        &
                            ims,ime,        jms,jme,        kms,kme,                        &
                            its,ite,        jts,jte,        kts,kte,                        &
                         cmr_sfcdif,     chr_sfcdif,     cmc_sfcdif,                        &
	                 chc_sfcdif,    cmgr_sfcdif,    chgr_sfcdif,                        &
                           tr_urb2d,       tb_urb2d,       tg_urb2d,                        & !H urban
	                   tc_urb2d,       qc_urb2d,       uc_urb2d,                        & !H urban
                         xxxr_urb2d,     xxxb_urb2d,     xxxg_urb2d, xxxc_urb2d,            & !H urban
                          trl_urb3d,      tbl_urb3d,      tgl_urb3d,                        & !H urban
                           sh_urb2d,       lh_urb2d,        g_urb2d,   rn_urb2d,  ts_urb2d, & !H urban
                         psim_urb2d,     psih_urb2d,      u10_urb2d,  v10_urb2d,            & !O urban
                       GZ1OZ0_urb2d,     AKMS_URB2D,                                        & !O urban
                          th2_urb2d,       q2_urb2d,      ust_urb2d,                        & !O urban
                             declin,          hrang,                                        & !I urban
                    num_roof_layers,num_wall_layers,num_road_layers,                        & !I urban
                                dzr,            dzb,            dzg,                        & !I urban
                         cmcr_urb2d,      tgr_urb2d,     tgrl_urb3d,  smr_urb3d,            & !H urban
                        drelr_urb2d,    drelb_urb2d,    drelg_urb2d,                        & !H urban
                      flxhumr_urb2d,  flxhumb_urb2d,  flxhumg_urb2d,                        & !H urban
                             julday,             yr,                                        & !H urban
                          frc_urb2d,    utype_urb2d,                                        & !I urban
                                chs,           chs2,           cqs2,                        & !H
                      num_urban_ndm,  urban_map_zrd,  urban_map_zwd, urban_map_gd,          & !I multi-layer urban
                       urban_map_zd,  urban_map_zdf,   urban_map_bd, urban_map_wd,          & !I multi-layer urban
                      urban_map_gbd,  urban_map_fbd,                                        & !I multi-layer urban
                       num_urban_hi,                                                        & !I multi-layer urban
                          trb_urb4d,      tw1_urb4d,      tw2_urb4d,  tgb_urb4d,            & !H multi-layer urban
                         tlev_urb3d,     qlev_urb3d,                                        & !H multi-layer urban
                       tw1lev_urb3d,   tw2lev_urb3d,                                        & !H multi-layer urban
                        tglev_urb3d,    tflev_urb3d,                                        & !H multi-layer urban
                        sf_ac_urb3d,    lf_ac_urb3d,    cm_ac_urb3d,                        & !H multi-layer urban
                       sfvent_urb3d,   lfvent_urb3d,                                        & !H multi-layer urban
                       sfwin1_urb3d,   sfwin2_urb3d,                                        & !H multi-layer urban
                         sfw1_urb3d,     sfw2_urb3d,      sfr_urb3d,  sfg_urb3d,            & !H multi-layer urban
                           lp_urb2d,       hi_urb2d,       lb_urb2d,  hgt_urb2d,            & !H multi-layer urban
                           mh_urb2d,     stdh_urb2d,       lf_urb2d,                        & !SLUCM
                        theta_urban,      rho_urban,        p_urban,        ust,            & !I multi-layer urban
                                gmt,         julday,          XLONG,       XLAT,            & !I multi-layer urban
                            a_u_bep,        a_v_bep,        a_t_bep,    a_q_bep,            & !O multi-layer urban
                            a_e_bep,        b_u_bep,        b_v_bep,                        & !O multi-layer urban
                            b_t_bep,        b_q_bep,        b_e_bep,    dlg_bep,            & !O multi-layer urban
                           dl_u_bep,         sf_bep,         vl_bep)                          !O multi-layer urban
if (debug_flag)then
           print*,'URBAN OUT'
           print*,'sf_urban_physics',sf_urban_physics 
           print*,'nsoil',nsoil
           print*,'ivgtyp',ivgtyp
           print*,'itimestep',itimestep
           print*,'dtbl',dt
           print*,'coszen',coszin
           print*,'xlat',xlat
           print*,'t_phy',t_phy
           print*,'qv_curr',qv_curr
           print*,'u_phy',u_phy
           print*,'v_phy',v_phy
           print*,'swdown',swdown
           print*,'glw',glw
           print*,'p8w',p8w
           print*,'rainbl',rainbl
           print*,'dz8w',dz8w
           print*,'znt',znt
           print*,'tsk',tsk
           print*,'hfx',hfx
           print*,'qfx',qfx
           print*,'lh',lh
           print*,'grdflx',grdflx
           print*,'albedo',albedo
           print*,'emiss',emiss
           print*,'qsfc',qsfc
           print*,'ids',ids
           print*,'ide',ide
           print*,'jds',jds
           print*,'jde',jde
           print*,'kds',kds
           print*,'kde',kde
           print*,'ims',ims
           print*,'ime',ime
           print*,'jms',jms
           print*,'jme',jme
           print*,'kms',kms
           print*,'kme',kme
           print*,'its',its
           print*,'ite',ite
           print*,'jts',jts
           print*,'jte',jte
           print*,'kts',kts
           print*,'kte',kte
           print*,'cmr_sfcdif',cmr_sfcdif
           print*,'chr_sfcdif',chr_sfcdif
           print*,'cmc_sfcdif',cmc_sfcdif
           print*,'chc_sfcdif',chc_sfcdif
           print*,'cmgr_sfcdif',cmgr_sfcdif
           print*,'chgr_sfcdif',chgr_sfcdif
          print*,'TR_URB2D',TR_URB2D
          print*,'TB_URB2D',TB_URB2D
          print*,'TG_URB2D',TG_URB2D
          print*,'TC_URB2D',TC_URB2D
          print*,'QC_URB2D',QC_URB2D
          print*,'UC_URB2D',UC_URB2D
          print*,'XXXR_URB2D',XXXR_URB2D
          print*,'XXXB_URB2D',XXXB_URB2D
          print*,'XXXG_URB2D',XXXG_URB2D
          print*,'XXXC_URB2D',XXXC_URB2D
          print*,'TRL_URB3D',TRL_URB3D
          print*,'TBL_URB3D',TBL_URB3D
          print*,'TGL_URB3D',TGL_URB3D
          print*,'SH_URB2D',SH_URB2D
          print*,'LH_URB2D',LH_URB2D
          print*,'G_URB2D',G_URB2D
          print*,'RN_URB2D',RN_URB2D
          print*,'TS_URB2D',TS_URB2D
          print*,'PSIM_URB2D',PSIM_URB2D
          print*,'PSIH_URB2D',PSIH_URB2D
          print*,'U10_URB2D',U10_URB2D
          print*,'V10_URB2D',V10_URB2D
          print*,'GZ1OZ0_URB2D',GZ1OZ0_URB2D
          print*,'AKMS_URB2D',AKMS_URB2D
          print*,'TH2_URB2D',TH2_URB2D
          print*,'Q2_URB2D',Q2_URB2D
          print*,'UST_URB2D',UST_URB2D
          print*,'declin',declin
          print*,'hrang',hrang
          print*,'num_roof_layers',num_roof_layers
          print*,'num_wall_layers',num_wall_layers
          print*,'num_road_layers',num_road_layers
          print*,'DZR',DZR
          print*,'DZB',DZB
          print*,'DZG',DZG
          print*,'CMCR_URB2D',CMCR_URB2D
          print*,'TGR_URB2D',TGR_URB2D
          print*,'TGRL_URB3D',TGRL_URB3D
          print*,'SMR_URB3D',SMR_URB3D
          print*,'DRELR_URB2D',DRELR_URB2D
          print*,'DRELB_URB2D',DRELB_URB2D
          print*,'DRELG_URB2D',DRELG_URB2D
          print*,'FLXHUMR_URB2D',FLXHUMR_URB2D
          print*,'FLXHUMB_URB2D',FLXHUMB_URB2D
          print*,'FLXHUMG_URB2D',FLXHUMG_URB2D
          print*,'julday',julday
          print*,'yr',yr
          print*,'FRC_URB2D',FRC_URB2D
          print*,'UTYPE_URB2D',UTYPE_URB2D
          print*,'chs',chs
          print*,'chs2',chs2
          print*,'cqs2',cqs2
          print*,'num_urban_ndm',num_urban_ndm
          print*,'urban_map_zrd',urban_map_zrd
          print*,'urban_map_zwd',urban_map_zwd
          print*,'urban_map_gd',urban_map_gd
          print*,'urban_map_zd',urban_map_zd
          print*,'urban_map_zdf',urban_map_zdf
          print*,'urban_map_bd',urban_map_bd
          print*,'urban_map_wd',urban_map_wd
          print*,'urban_map_gbd',urban_map_gbd
          print*,'urban_map_fbd',urban_map_fbd
          print*,'num_urban_hi',num_urban_hi
          print*,'TRB_URB4D',TRB_URB4D
          print*,'TW1_URB4D',TW1_URB4D
          print*,'TW2_URB4D',TW2_URB4D
          print*,'TGB_URB4D',TGB_URB4D
          print*,'TLEV_URB3D',TLEV_URB3D
          print*,'QLEV_URB3D',QLEV_URB3D
          print*,'TW1LEV_URB3D',TW1LEV_URB3D
          print*,'TW2LEV_URB3D',TW2LEV_URB3D
          print*,'TGLEV_URB3D',TGLEV_URB3D
          print*,'TFLEV_URB3D',TFLEV_URB3D
          print*,'SF_AC_URB3D',SF_AC_URB3D
          print*,'LF_AC_URB3D',LF_AC_URB3D
          print*,'CM_AC_URB3D',CM_AC_URB3D
          print*,'SFVENT_URB3D',SFVENT_URB3D
          print*,'LFVENT_URB3D',LFVENT_URB3D
          print*,'SFWIN1_URB3D',SFWIN1_URB3D
          print*,'SFWIN2_URB3D',SFWIN2_URB3D
          print*,'SFW1_URB3D',SFW1_URB3D
          print*,'SFW2_URB3D',SFW2_URB3D
          print*,'SFR_URB3D',SFR_URB3D
          print*,'SFG_URB3D',SFG_URB3D
          print*,'LP_URB2D',LP_URB2D
          print*,'HI_URB2D',HI_URB2D      
          print*,'LB_URB2D',LB_URB2D
          print*,'HGT_URB2D',HGT_URB2D
          print*,'MH_URB2D',MH_URB2D
          print*,'STDH_URB2D',STDH_URB2D
          print*,'LF_URB2D',LF_URB2D
          print*,'theta_urban',theta_urban
          print*,'rho_urban',rho_urban
          print*,'p_urban',p_urban
          print*,'ust',ust
          print*,'gmt',gmt
          print*,'julday',julday
          print*,'xlong',xlong
          print*,'xlat',xlat
          print*,'A_U_BEP',A_U_BEP
          print*,'A_V_BEP',A_V_BEP
          print*,'A_T_BEP',A_T_BEP
          print*,'A_Q_BEP',A_Q_BEP
          print*,'A_E_BEP',A_E_BEP
          print*,'B_U_BEP',B_U_BEP
          print*,'B_V_BEP',B_V_BEP
          print*,'B_T_BEP',B_T_BEP
          print*,'B_Q_BEP',B_Q_BEP
          print*,'B_E_BEP',B_E_BEP
          print*,'DLG_BEP',DLG_BEP
          print*,'DL_U_BEP',DL_U_BEP  
          print*,'SF_BEP',SF_BEP
          print*,'VL_BEP',VL_BEP
endif

end subroutine

!------------------------------------------------------------------------------------

SUBROUTINE CALC_DECLIN ( NOWDATE, LATITUDE, LONGITUDE, COSZ, JULIAN, &
                         HRANG_OUT, DECLIN_OUT, GMT_OUT, JULDAY_OUT )

  USE MODULE_DATE_UTILITIES
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------

! !ARGUMENTS:
   CHARACTER(LEN=19), INTENT(IN)  :: NOWDATE    ! YYYY-MM-DD_HH:MM:SS
   REAL,              INTENT(IN)  :: LATITUDE
   REAL,              INTENT(IN)  :: LONGITUDE
   REAL,              INTENT(OUT) :: COSZ
   REAL,              INTENT(OUT) :: JULIAN
   !REAL,    OPTIONAL, INTENT(OUT) :: HRANG_OUT
   !REAL,    OPTIONAL, INTENT(OUT) :: DECLIN_OUT
   !REAL,    OPTIONAL, INTENT(OUT) :: GMT_OUT
   !INTEGER, OPTIONAL, INTENT(OUT) :: JULDAY_OUT
   REAL,              INTENT(OUT) :: HRANG_OUT
   REAL,              INTENT(OUT) :: DECLIN_OUT
   REAL,              INTENT(OUT) :: GMT_OUT
   INTEGER,           INTENT(OUT) :: JULDAY_OUT
   REAL                           :: OBECL
   REAL                           :: SINOB
   REAL                           :: SXLONG
   REAL                           :: ARG
   REAL                           :: TLOCTIM
   REAL                           :: HRANG
   REAL                           :: DECLIN
   REAL                           :: GMT
   INTEGER                        :: JULDAY
   INTEGER                        :: IDAY
   INTEGER                        :: IHOUR
   INTEGER                        :: IMINUTE
   INTEGER                        :: ISECOND

   REAL, PARAMETER :: DEGRAD = 3.14159265/180.
   REAL, PARAMETER :: DPD    = 360./365.

   CALL GETH_IDTS(NOWDATE(1:10), NOWDATE(1:4)//"-01-01", IDAY)
   READ(NOWDATE(12:13), *) IHOUR
   READ(NOWDATE(15:16), *) IMINUTE
   READ(NOWDATE(18:19), *) ISECOND
   GMT = REAL(IHOUR) + IMINUTE/60.0 + ISECOND/3600.0
   JULIAN = REAL(IDAY) + GMT/24.
   JULDAY = IDAY

!
! FOR SHORT WAVE RADIATION

   DECLIN=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.

   OBECL=23.5*DEGRAD
   SINOB=SIN(OBECL)

!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:

   IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)*DEGRAD
   IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)*DEGRAD
   ARG=SINOB*SIN(SXLONG)
   DECLIN=ASIN(ARG)

   TLOCTIM = REAL(IHOUR) + REAL(IMINUTE)/60.0 + REAL(ISECOND)/3600.0 + LONGITUDE/15.0 ! LOCAL TIME IN HOURS
   TLOCTIM = AMOD(TLOCTIM+24.0, 24.0)
   HRANG=15.*(TLOCTIM-12.)*DEGRAD
   COSZ=SIN(LATITUDE*DEGRAD)*SIN(DECLIN)+COS(LATITUDE*DEGRAD)*COS(DECLIN)*COS(HRANG)
   
   !IF (PRESENT(HRANG_OUT ))  HRANG_OUT = HRANG
   !IF (PRESENT(DECLIN_OUT)) DECLIN_OUT = DECLIN
   !IF (PRESENT(GMT_OUT   ))    GMT_OUT = GMT
   !IF (PRESENT(JULDAY_OUT)) JULDAY_OUT = JULDAY
   HRANG_OUT = HRANG
   DECLIN_OUT = DECLIN
   GMT_OUT = GMT
   JULDAY_OUT = JULDAY

 END SUBROUTINE CALC_DECLIN
