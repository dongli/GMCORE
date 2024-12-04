module zm_conv

  !-----------------------------------------------------------------------------
  ! Purpose:
  !
  ! Interface from Zhang-McFarlane convection scheme, includes evaporation of
  ! convective precipitation from the ZM scheme.
  !
  ! Apr 2006: RBN: Code added to perform a dilute ascent for closure of the CM
  ! mass flux based on an entraining plume a la Raymond and Blythe (1992)
  !
  ! Author: Byron Boville, from code in tphysbc
  !
  !-----------------------------------------------------------------------------

  use shr_kind_mod   , only: r8 => shr_kind_r8
  use spmd_utils     , only: masterproc
  use ppgrid         , only: pcols, pver, pverp
  use cloud_fraction , only: cldfrc_fice
  use physconst      , only: cpair, epsilo, gravit, latice, latvap, tmelt, rair, &
                             cpwv, cpliq, rh2o
  use cam_abortutils , only: endrun
  use cam_logfile    , only: iulog
  use zm_microphysics, only: zm_mphy, zm_aero_t, zm_conv_t

  implicit none

  save

  private

  public zm_convi                 ! ZM schemea
  public zm_convr                 ! ZM schemea
  public zm_conv_evap             ! Evaporation of precip from ZM schemea
  public convtran                 ! Convective transport
  public momtran                  ! Convective momentum transport

  real(r8), parameter :: capelmt    = 70.0_r8  ! Threshold value of CAPE for deep convection.
  real(r8), parameter :: tiedke_add = 0.5_r8
  real(r8), parameter :: c1         = 6.112_r8
  real(r8), parameter :: c2         = 17.67_r8
  real(r8), parameter :: c3         = 243.5_r8
  real(r8), parameter :: dcon       = 25.e-6_r8
  real(r8), parameter :: mucon      = 5.3_r8

  real(r8) rl                     ! Latent heat of vaporization
  real(r8) cpres                  ! Specific heat at constant pressure in j/kg-degk.
  real(r8) ke                     ! Tunable evaporation efficiency set from namelist input zmconv_ke
  real(r8) ke_lnd
  real(r8) c0_lnd                 ! Set from namelist input zmconv_c0_lnd
  real(r8) c0_ocn                 ! Set from namelist input zmconv_c0_ocn
  integer  num_cin                ! Set from namelist input zmconv_num_cin
                                  ! The number of negative buoyancy regions that are allowed
                                  ! Before the convection top and CAPE calculations are completed.
  logical  zm_org
  real(r8) tau                    ! Convective time scale
  real(r8) tfreez
  real(r8) eps1
  real(r8) momcu
  real(r8) momcd

  logical  zmconv_microp

  logical  no_deep_pbl            ! default = .false.
                                  ! no_deep_pbl = .true. eliminates deep convection entirely within PBL

  ! moved from moistconvection.F90
  real(r8) rgrav                  ! reciprocal of grav
  real(r8) rgas                   ! gas constant for dry air
  real(r8) grav
  real(r8) cp

  integer  limcnv                 ! top interface level limit for convection

contains

  subroutine zm_convi(limcnv_in, zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke, zmconv_ke_lnd, &
                      zmconv_momcu, zmconv_momcd, zmconv_num_cin, zmconv_org, &
                      zmconv_microp_in, no_deep_pbl_in)

    integer , intent(in)           :: limcnv_in
    integer , intent(in)           :: zmconv_num_cin
    real(r8), intent(in)           :: zmconv_c0_lnd
    real(r8), intent(in)           :: zmconv_c0_ocn
    real(r8), intent(in)           :: zmconv_ke
    real(r8), intent(in)           :: zmconv_ke_lnd
    real(r8), intent(in)           :: zmconv_momcu
    real(r8), intent(in)           :: zmconv_momcd
    logical , intent(in)           :: zmconv_org
    logical , intent(in)           :: zmconv_microp_in
    logical , intent(in), optional :: no_deep_pbl_in

    limcnv  = limcnv_in
    tfreez  = tmelt
    eps1    = epsilo
    rl      = latvap
    cpres   = cpair
    rgrav   = 1.0_r8 / gravit
    rgas    = rair
    grav    = gravit
    cp      = cpres
    c0_lnd  = zmconv_c0_lnd
    c0_ocn  = zmconv_c0_ocn
    num_cin = zmconv_num_cin
    ke      = zmconv_ke
    ke_lnd  = zmconv_ke_lnd
    zm_org  = zmconv_org
    momcu   = zmconv_momcu
    momcd   = zmconv_momcd

    zmconv_microp = zmconv_microp_in

    if (present(no_deep_pbl_in))  then
      no_deep_pbl = no_deep_pbl_in
    else
      no_deep_pbl = .false.
    end if

    tau = 3600.0_r8

    if (masterproc) then
      write(iulog, *) 'tuning parameters zm_convi: tau', tau
      write(iulog, *) 'tuning parameters zm_convi: c0_lnd', c0_lnd, ', c0_ocn', c0_ocn
      write(iulog, *) 'tuning parameters zm_convi: num_cin', num_cin
      write(iulog, *) 'tuning parameters zm_convi: ke', ke
      write(iulog, *) 'tuning parameters zm_convi: no_deep_pbl', no_deep_pbl
    end if

    if (masterproc) write(iulog, *) '**** ZM: DILUTE Buoyancy Calculation ****'

  end subroutine zm_convi

  subroutine zm_convr( &
    lchnk            , &
    ncol             , &
    t                , &
    qh               , & ! Specific humidity
    prec             , &
    jctop            , &
    jcbot            , &
    pblh             , &
    zm               , &
    geos             , &
    zi               , &
    qtnd             , &
    heat             , &
    pap              , &
    paph             , &
    dpp              , &
    delt             , &
    mcon             , &
    cme              , &
    cape             , & ! Convective available potential energy
    tpert            , &
    dlf              , &
    pflx             , &
    zdu              , &
    rprd             , &
    mu               , & ! Upward cloud mass flux (positive up) specified at interface
    md               , & ! Downward cloud mass flux (positive up)
    du               , & ! Detrainment in updraft specified in mid-layer
    eu               , & ! Detrainment in updraft
    ed               , & ! Detrainment in downdraft
    dp               , &
    dsubcld          , &
    jt               , &
    maxg             , &
    ideep            , &
    ql               , & ! Cloud liquid water
    rliq             , &
    landfrac         , &
    org              , &
    orgt             , &
    org2d            , &
    dif              , &
    dnlf             , &
    dnif             , &
    conv             , &
    aero             , &
    rice             )
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! Main driver for zhang-mcfarlane convection scheme
    !
    ! Method:
    ! performs deep convective adjustment based on mass-flux closure
    ! algorithm.
    !
    ! Author:guang jun zhang, m.lazare, n.mcfarlane. CAM Contact: P. Rasch
    !
    ! This is contributed code not fully standardized by the CAM core group.
    ! All variables have been typed, where most are identified in comments
    ! The current procedure will be reimplemented in a subsequent version
    ! of the CAM where it will include a more straightforward formulation
    ! and will make use of the standard CAM nomenclature
    !
    !---------------------------------------------------------------------------

    use phys_control, only: cam_physpkg_is

    integer , intent(in   ) :: lchnk
    integer , intent(in   ) :: ncol
    real(r8), intent(in   ), dimension(pcols,pver  ) :: t 
    real(r8), intent(in   ), dimension(pcols,pver  ) :: qh      ! Specific humidity.
    real(r8), intent(in   ), dimension(pcols,pver  ) :: pap
    real(r8), intent(in   ), dimension(pcols,pver+1) :: paph
    real(r8), intent(in   ), dimension(pcols,pver  ) :: dpp     ! Local sigma half-level thickness (i.e. dshj).
    real(r8), intent(in   ), dimension(pcols,pver  ) :: zm
    real(r8), intent(in   ), dimension(pcols       ) :: geos
    real(r8), intent(in   ), dimension(pcols,pver+1) :: zi
    real(r8), intent(in   ), dimension(pcols       ) :: pblh
    real(r8), intent(in   ), dimension(pcols       ) :: tpert
    real(r8), intent(in   ), dimension(pcols       ) :: landfrac
    real(r8), intent(  out), dimension(pcols,pver  ) :: qtnd    ! Specific humidity tendency (kg/kg/s)
    real(r8), intent(  out), dimension(pcols,pver  ) :: heat    ! Heating rate (dry static energy tendency, W/kg)
    real(r8), intent(  out), dimension(pcols,pverp ) :: mcon
    real(r8), intent(  out), dimension(pcols,pver  ) :: dlf     ! Scattrd version of the detraining cld h2o tend
    real(r8), intent(  out), dimension(pcols,pverp ) :: pflx    ! Scattered precip flux at each level
    real(r8), intent(  out), dimension(pcols,pver  ) :: cme
    real(r8), intent(  out), dimension(pcols       ) :: cape    ! Convective available potential energy.
    real(r8), intent(  out), dimension(pcols,pver  ) :: zdu
    real(r8), intent(  out), dimension(pcols,pver  ) :: rprd    ! Rain production rate
    real(r8), intent(  out), dimension(pcols,pver  ) :: dif     ! Detrained convective cloud ice mixing ratio.
    real(r8), intent(  out), dimension(pcols,pver  ) :: dnlf    ! Detrained convective cloud water num concen.
    real(r8), intent(  out), dimension(pcols,pver  ) :: dnif    ! Detrained convective cloud ice num concen.
    real(r8), intent(  out), dimension(pcols,pver  ) :: mu      ! Upward cloud mass flux (positive up) specified at interface
    real(r8), intent(  out), dimension(pcols,pver  ) :: eu
    real(r8), intent(  out), dimension(pcols,pver  ) :: du
    real(r8), intent(  out), dimension(pcols,pver  ) :: md      ! Downward cloud mass flux (positive up)
    real(r8), intent(  out), dimension(pcols,pver  ) :: ed
    real(r8), intent(  out), dimension(pcols,pver  ) :: dp      ! Layer thickness in mbs (between upper/lower interface).
    real(r8), intent(  out), dimension(pcols       ) :: dsubcld ! Layer thickness in mbs between LCL and maxi.
    real(r8), intent(  out), dimension(pcols       ) :: jctop   ! Row of top-of-deep-convection indices passed out.
    real(r8), intent(  out), dimension(pcols       ) :: jcbot   ! Row of base of cloud indices passed out.
    real(r8), intent(  out), dimension(pcols       ) :: prec
    real(r8), intent(  out), dimension(pcols       ) :: rliq    ! Reserved liquid (not yet in cldliq) for energy integrals
    real(r8), intent(  out), dimension(pcols       ) :: rice    ! Reserved ice (not yet in cldce) for energy integrals
    integer , intent(  out), dimension(pcols       ) :: ideep   ! Column indices of gathered points
    integer , intent(inout), dimension(pcols       ) :: jt      ! Top  level index of deep cumulus convection.
    integer , intent(inout), dimension(pcols       ) :: maxg    ! Gathered values of maxi.
    real(r8), intent(inout), dimension(pcols,pver  ) :: ql      ! Grid slice of cloud liquid water.
    type(zm_conv_t), intent(inout) :: conv
    type(zm_aero_t), intent(inout) :: aero         ! Aerosol object. intent(inout) because the
                                                   ! gathered arrays are set here
                                                   ! before passing object
                                                   ! to microphysics

    type(zm_conv_t) :: loc_conv

    real(r8), pointer, dimension(:,:) :: org           ! Only used if zm_org is true
    real(r8), pointer, dimension(:,:) :: orgt          ! Only used if zm_org is true
    real(r8), pointer, dimension(:,:) :: org2d         ! Only used if zm_org is true

    real(r8), dimension(pcols       ) :: zs
    real(r8), dimension(pcols,pver  ) :: dlg       ! Gathered version of the detraining cloud liquid water tendency
    real(r8), dimension(pcols,pverp ) :: pflxg     ! Gathered precipitation flux at each level
    real(r8), dimension(pcols,pver  ) :: cug       ! Gathered condensation rate
    real(r8), dimension(pcols,pver  ) :: evpg      ! Gathered evaperation rate of rain in downdraft
    real(r8), dimension(pcols       ) :: orgavg
    real(r8), dimension(pcols       ) :: dptot
    real(r8), dimension(pcols       ) :: mumax
    real(r8), dimension(pcols       ) :: pblt      ! Row of PBL top indices
    real(r8), dimension(pcols,pver  ) :: q         ! Grid slice of mixing ratio
    real(r8), dimension(pcols,pver  ) :: p         ! Grid slice of ambient mid-layer pressure in mbs
    real(r8), dimension(pcols,pver  ) :: z         ! Grid slice of ambient mid-layer height in metres
    real(r8), dimension(pcols,pver  ) :: s         ! Grid slice of scaled dry static energy (t+gz/cp)
    real(r8), dimension(pcols,pver  ) :: tp        ! Grid slice of parcel temperature
    real(r8), dimension(pcols,pver+1) :: zf        ! Grid slice of ambient interface height in metres
    real(r8), dimension(pcols,pver+1) :: pf        ! Grid slice of ambient interface pressure in mbs
    real(r8), dimension(pcols,pver  ) :: qstp      ! Grid slice of parcel temp. saturation mixing ratio
    real(r8), dimension(pcols       ) :: tl        ! Row of parcel temperature at LCL
    integer , dimension(pcols       ) :: lcl       ! Base level index of deep cumulus convection
    integer , dimension(pcols       ) :: lel       ! Index of highest theoretical convective plume
    integer , dimension(pcols       ) :: lon       ! Index of onset level for deep convection
    integer , dimension(pcols       ) :: maxi      ! Index of level with largest moist static energy
    real(r8), dimension(pcols,pver  ) :: qg        ! Grid slice of gathered values of q
    real(r8), dimension(pcols,pver  ) :: tg        ! Grid slice of temperature at interface
    real(r8), dimension(pcols,pver  ) :: pg        ! Grid slice of gathered values of p
    real(r8), dimension(pcols,pver  ) :: zg        ! Grid slice of gathered values of z
    real(r8), dimension(pcols,pver  ) :: sg        ! Grid slice of gathered values of s
    real(r8), dimension(pcols,pver  ) :: tpg       ! Grid slice of gathered values of tp
    real(r8), dimension(pcols,pver+1) :: zfg       ! Grid slice of gathered values of zf
    real(r8), dimension(pcols,pver  ) :: qstpg     ! Grid slice of gathered values of qstp
    real(r8), dimension(pcols,pver  ) :: ug        ! Grid slice of gathered values of u
    real(r8), dimension(pcols,pver  ) :: vg        ! Grid slice of gathered values of v
    real(r8), dimension(pcols,pver  ) :: cmeg
    real(r8), dimension(pcols,pver  ) :: rprdg     ! Gathered rain production rate
    real(r8), dimension(pcols       ) :: capeg     ! Gathered convective available potential energy
    real(r8), dimension(pcols       ) :: tlg       ! Grid slice of gathered values of tl
    real(r8), dimension(pcols       ) :: landfracg ! Grid slice of gathered land fraction
    integer , dimension(pcols       ) :: lclg      ! Gathered values of LCL
    integer , dimension(pcols       ) :: lelg
    real(r8), dimension(pcols,pver  ) :: dqdt      ! Mixing ratio tendency at gathered points
    real(r8), dimension(pcols,pver  ) :: dsdt      ! Dry static energy tendency at gathered points
    real(r8), dimension(pcols,pver  ) :: sd        ! Grid slice of dry static energy in downdraft
    real(r8), dimension(pcols,pver  ) :: qd        ! Grid slice of mixing ratio in downdraft
    real(r8), dimension(pcols,pver  ) :: mc        ! Net upward (scaled by mb) cloud mass flux
    real(r8), dimension(pcols,pver  ) :: qhat      ! Grid slice of upper interface mixing ratio
    real(r8), dimension(pcols,pver  ) :: qu        ! Grid slice of mixing ratio in updraft
    real(r8), dimension(pcols,pver  ) :: su        ! Grid slice of dry static energy in updraft
    real(r8), dimension(pcols,pver  ) :: qs        ! Grid slice of saturation mixing ratio
    real(r8), dimension(pcols,pver  ) :: shat      ! Grid slice of upper interface dry static energy
    real(r8), dimension(pcols,pver  ) :: hmn       ! Moist static energy
    real(r8), dimension(pcols,pver  ) :: hsat      ! Saturated moist static energy
    real(r8), dimension(pcols,pver  ) :: qlg
    real(r8), dimension(pcols,pver  ) :: dudt      ! U-wind tendency at gathered points
    real(r8), dimension(pcols,pver  ) :: dvdt      ! V-wind tendency at gathered points
    real(r8), dimension(pcols,pver  ) :: qldeg     ! Cloud liquid water mixing ratio for detrainment (kg/kg)
    real(r8), dimension(pcols       ) :: mb        ! Cloud base mass flux
    integer , dimension(pcols       ) :: jlcl
    integer , dimension(pcols       ) :: j0        ! Detrainment initiation level index
    integer , dimension(pcols       ) :: jd        ! Downdraft initiation level index

    logical doliq
    integer lengath
    integer i, k, kk, l, m
    integer msg                                        ! Number of missing moisture levels at the top of model.
    real(r8) delt                                      ! Length of model time-step in seconds.
    real(r8) precip
    real(r8) qdifr
    real(r8) sdifr
    real(r8) negadq

    ! Set internal variable "msg" (convection limit) to "limcnv-1"
    msg = limcnv - 1

    ! Initialize necessary arrays.
    if (zm_org) then
      orgt = 0
    end if

    qtnd = 0
    heat = 0
    mcon = 0
    rliq(:ncol) = 0
    rice(:ncol) = 0

    if (zmconv_microp) then
      allocate(loc_conv%frz       (pcols,pver))
      allocate(loc_conv%sprd      (pcols,pver))
      allocate(loc_conv%wu        (pcols,pver))
      allocate(loc_conv%qi        (pcols,pver))
      allocate(loc_conv%qliq      (pcols,pver))
      allocate(loc_conv%qice      (pcols,pver))
      allocate(loc_conv%qrain     (pcols,pver))
      allocate(loc_conv%qsnow     (pcols,pver))
      allocate(loc_conv%di        (pcols,pver))
      allocate(loc_conv%dnl       (pcols,pver))
      allocate(loc_conv%dni       (pcols,pver))
      allocate(loc_conv%qnl       (pcols,pver))
      allocate(loc_conv%qni       (pcols,pver))
      allocate(loc_conv%qnr       (pcols,pver))
      allocate(loc_conv%qns       (pcols,pver))
      allocate(loc_conv%qide      (pcols,pver))
      allocate(loc_conv%qncde     (pcols,pver))
      allocate(loc_conv%qnide     (pcols,pver))
      allocate(loc_conv%autolm    (pcols,pver))
      allocate(loc_conv%accrlm    (pcols,pver))
      allocate(loc_conv%bergnm    (pcols,pver))
      allocate(loc_conv%fhtimm    (pcols,pver))
      allocate(loc_conv%fhtctm    (pcols,pver))
      allocate(loc_conv%fhmlm     (pcols,pver))
      allocate(loc_conv%hmpim     (pcols,pver))
      allocate(loc_conv%accslm    (pcols,pver))
      allocate(loc_conv%dlfm      (pcols,pver))
      allocate(loc_conv%cmel      (pcols,pver))
      allocate(loc_conv%autoln    (pcols,pver))
      allocate(loc_conv%accrln    (pcols,pver))
      allocate(loc_conv%bergnn    (pcols,pver))
      allocate(loc_conv%fhtimn    (pcols,pver))
      allocate(loc_conv%fhtctn    (pcols,pver))
      allocate(loc_conv%fhmln     (pcols,pver))
      allocate(loc_conv%accsln    (pcols,pver))
      allocate(loc_conv%activn    (pcols,pver))
      allocate(loc_conv%dlfn      (pcols,pver))
      allocate(loc_conv%autoim    (pcols,pver))
      allocate(loc_conv%accsim    (pcols,pver))
      allocate(loc_conv%difm      (pcols,pver))
      allocate(loc_conv%cmei      (pcols,pver))
      allocate(loc_conv%nuclin    (pcols,pver))
      allocate(loc_conv%autoin    (pcols,pver))
      allocate(loc_conv%accsin    (pcols,pver))
      allocate(loc_conv%hmpin     (pcols,pver))
      allocate(loc_conv%difn      (pcols,pver))
      allocate(loc_conv%trspcm    (pcols,pver))
      allocate(loc_conv%trspcn    (pcols,pver))
      allocate(loc_conv%trspim    (pcols,pver))
      allocate(loc_conv%trspin    (pcols,pver))
      allocate(loc_conv%lambdadpcu(pcols,pver))
      allocate(loc_conv%mudpcu    (pcols,pver))
      allocate(loc_conv%dcape     (pcols     ))
    end if

    !
    ! Initialize convective tendencies
    !
    prec(:ncol) = 0
    do k = 1, pver
      do i = 1, ncol
        dqdt (i,k) = 0
        dsdt (i,k) = 0
        dudt (i,k) = 0
        dvdt (i,k) = 0
        pflx (i,k) = 0
        pflxg(i,k) = 0
        cme  (i,k) = 0
        rprd (i,k) = 0
        zdu  (i,k) = 0
        ql   (i,k) = 0
        qlg  (i,k) = 0
        dlf  (i,k) = 0
        dlg  (i,k) = 0
        qldeg(i,k) = 0
        dif  (i,k) = 0
        dnlf (i,k) = 0
        dnif (i,k) = 0
      end do
    end do

    if (zmconv_microp) then
      do k = 1, pver
        do i = 1, ncol
          loc_conv%qliq   (i,k) = 0
          loc_conv%qice   (i,k) = 0
          loc_conv%di     (i,k) = 0
          loc_conv%qrain  (i,k) = 0
          loc_conv%qsnow  (i,k) = 0
          loc_conv%dnl    (i,k) = 0
          loc_conv%dni    (i,k) = 0
          loc_conv%wu     (i,k) = 0
          loc_conv%qnl    (i,k) = 0
          loc_conv%qni    (i,k) = 0
          loc_conv%qnr    (i,k) = 0
          loc_conv%qns    (i,k) = 0
          loc_conv%frz    (i,k) = 0
          loc_conv%sprd   (i,k) = 0
          loc_conv%qide   (i,k) = 0
          loc_conv%qncde  (i,k) = 0
          loc_conv%qnide  (i,k) = 0
          loc_conv%autolm (i,k) = 0
          loc_conv%accrlm (i,k) = 0
          loc_conv%bergnm (i,k) = 0
          loc_conv%fhtimm (i,k) = 0
          loc_conv%fhtctm (i,k) = 0
          loc_conv%fhmlm  (i,k) = 0
          loc_conv%hmpim  (i,k) = 0
          loc_conv%accslm (i,k) = 0
          loc_conv%dlfm   (i,k) = 0
          loc_conv%autoln (i,k) = 0
          loc_conv%accrln (i,k) = 0
          loc_conv%bergnn (i,k) = 0
          loc_conv%fhtimn (i,k) = 0
          loc_conv%fhtctn (i,k) = 0
          loc_conv%fhmln  (i,k) = 0
          loc_conv%accsln (i,k) = 0
          loc_conv%activn (i,k) = 0
          loc_conv%dlfn   (i,k) = 0
          loc_conv%cmel   (i,k) = 0
          loc_conv%autoim (i,k) = 0
          loc_conv%accsim (i,k) = 0
          loc_conv%difm   (i,k) = 0
          loc_conv%cmei   (i,k) = 0
          loc_conv%nuclin (i,k) = 0
          loc_conv%autoin (i,k) = 0
          loc_conv%accsin (i,k) = 0
          loc_conv%hmpin  (i,k) = 0
          loc_conv%difn   (i,k) = 0
          loc_conv%trspcm (i,k) = 0
          loc_conv%trspcn (i,k) = 0
          loc_conv%trspim (i,k) = 0
          loc_conv%trspin (i,k) = 0
          conv%qi         (i,k) = 0
          conv%frz        (i,k) = 0
          conv%sprd       (i,k) = 0
          conv%qi         (i,k) = 0
          conv%qliq       (i,k) = 0
          conv%qice       (i,k) = 0
          conv%qnl        (i,k) = 0
          conv%qni        (i,k) = 0
          conv%qnr        (i,k) = 0
          conv%qns        (i,k) = 0
          conv%qrain      (i,k) = 0
          conv%qsnow      (i,k) = 0
          conv%wu         (i,k) = 0
          conv%autolm     (i,k) = 0
          conv%accrlm     (i,k) = 0
          conv%bergnm     (i,k) = 0
          conv%fhtimm     (i,k) = 0
          conv%fhtctm     (i,k) = 0
          conv%fhmlm      (i,k) = 0
          conv%hmpim      (i,k) = 0
          conv%accslm     (i,k) = 0
          conv%dlfm       (i,k) = 0
          conv%autoln     (i,k) = 0
          conv%accrln     (i,k) = 0
          conv%bergnn     (i,k) = 0
          conv%fhtimn     (i,k) = 0
          conv%fhtctn     (i,k) = 0
          conv%fhmln      (i,k) = 0
          conv%accsln     (i,k) = 0
          conv%activn     (i,k) = 0
          conv%dlfn       (i,k) = 0
          conv%cmel       (i,k) = 0
          conv%autoim     (i,k) = 0
          conv%accsim     (i,k) = 0
          conv%difm       (i,k) = 0
          conv%cmei       (i,k) = 0
          conv%nuclin     (i,k) = 0
          conv%autoin     (i,k) = 0
          conv%accsin     (i,k) = 0
          conv%hmpin      (i,k) = 0
          conv%difn       (i,k) = 0
          conv%trspcm     (i,k) = 0
          conv%trspcn     (i,k) = 0
          conv%trspim     (i,k) = 0
          conv%trspin     (i,k) = 0
        end do
      end do
      conv%lambdadpcu     = (mucon + 1.0_r8) / dcon
      conv%mudpcu         = mucon
      loc_conv%lambdadpcu = conv%lambdadpcu
      loc_conv%mudpcu     = conv%mudpcu
    end if

    do i = 1, ncol
      pflx (i,pverp) = 0
      pflxg(i,pverp) = 0
    end do

    do i = 1, ncol
      pblt   (i) = pver
      dsubcld(i) = 0
      jctop  (i) = pver
      jcbot  (i) = 1
    end do

    if (zmconv_microp) then
      do i = 1, ncol
        conv%dcape    (i) = 0
        loc_conv%dcape(i) = 0
      end do
    end if

    if (zm_org) then
      ! Compute vertical average here
      orgavg = 0
      dptot  = 0

      do k = 1, pver
        do i = 1, ncol
          if (org(i,k) > 0) then
            orgavg(i) = orgavg(i) + dpp(i,k) * org(i,k)
            dptot (i) = dptot (i) + dpp(i,k)
          end if
        end do
      end do

      do i = 1, ncol
        if (dptot(i) > 0) then
          orgavg(i) = orgavg(i) / dptot(i)
        end if
      end do

      do k = 1, pver
        do i = 1, ncol
          org2d(i,k) = orgavg(i)
        end do
      end do
    end if

    !
    ! Calculate local pressure (mbs) and height (m) for both interface
    ! and mid-layer locations.
    !
    do i = 1, ncol
      zs(i) = geos(i) * rgrav
      pf(i,pver+1) = paph(i,pver+1) * 0.01_r8
      zf(i,pver+1) = zi(i,pver+1) + zs(i)
    end do
    do k = 1, pver
      do i = 1, ncol
        p (i,k) = pap (i,k) * 0.01_r8
        pf(i,k) = paph(i,k) * 0.01_r8
        z (i,k) = zm(i,k) + zs(i)
        zf(i,k) = zi(i,k) + zs(i)
      end do
    end do

    do k = pver - 1, msg + 1, -1
      do i = 1, ncol
        if (abs(z(i,k) - zs(i) - pblh(i)) < (zf(i,k) - zf(i,k+1)) * 0.5_r8) pblt(i) = k
      end do
    end do
    !
    ! Store incoming specific humidity field for subsequent calculation
    ! of precipitation (through change in storage).
    ! define dry static energy (normalized by cp).
    !
    do k = 1, pver
      do i = 1, ncol
        q   (i,k) = qh(i,k)
        s   (i,k) = t(i,k) + (grav / cpres) * z(i,k)
        tp  (i,k) = 0
        shat(i,k) = s(i,k)
        qhat(i,k) = q(i,k)
      end do
    end do

    do i = 1, ncol
      capeg  (i) = 0
      lclg   (i) = 1
      lelg   (i) = pver
      maxg   (i) = 1
      tlg    (i) = 400
      dsubcld(i) = 0
    end do

    if (cam_physpkg_is('cam3')) then

      ! For cam3 physics package, call non-dilute
      call buoyan( &
        lchnk    , &
        ncol     , &
        q        , &
        t        , &
        p        , &
        z        , &
        pf       , &
        tp       , &
        qstp     , &
        tl       , &
        rl       , &
        cape     , &
        pblt     , &
        lcl      , &
        lel      , &
        lon      , &
        maxi     , &
        rgas     , &
        grav     , &
        cpres    , &
        msg      , &
        tpert    )
    else

      !  Evaluate Tparcel, qs(Tparcel), buoyancy and CAPE,
      !     lcl, lel, parcel launch level at index maxi()=hmax
      call buoyan_dilute( &
        lchnk           , &
        ncol            , &
        q               , &
        t               , &
        p               , &
        z               , &
        pf              , &
        tp              , &
        qstp            , &
        tl              , &
        rl              , &
        cape            , &
        pblt            , &
        lcl             , &
        lel             , &
        lon             , &
        maxi            , &
        rgas            , &
        grav            , &
        cpres           , &
        msg             , &
        tpert           , &
        org2d           , &
        landfrac        )
    end if

    !
    ! Determine whether grid points will undergo some deep convection
    ! (ideep=1) or not (ideep=0), based on values of cape,lcl,lel
    ! (require cape.gt. 0 and lel<lcl as minimum conditions).
    !
    lengath = 0
    ideep   = 0
    do i = 1, ncol
      if (cape(i) > capelmt) then
        lengath = lengath + 1
        ideep(lengath) = i
      end if
    end do

    if (lengath == 0) return
    !
    ! Obtain gathered arrays necessary for ensuing calculations.
    !
    do k = 1, pver
      do i = 1, lengath
        dp   (i,k) = 0.01_r8 * dpp(ideep(i),k)
        qg   (i,k) =             q(ideep(i),k)
        tg   (i,k) =             t(ideep(i),k)
        pg   (i,k) =             p(ideep(i),k)
        zg   (i,k) =             z(ideep(i),k)
        sg   (i,k) =             s(ideep(i),k)
        tpg  (i,k) =            tp(ideep(i),k)
        zfg  (i,k) =            zf(ideep(i),k)
        qstpg(i,k) =          qstp(ideep(i),k)
        ug   (i,k) = 0
        vg   (i,k) = 0
      end do
    end do

    if (zmconv_microp) then
      if (aero%scheme == 'modal') then
        do m = 1, aero%nmodes
          do k = 1, pver
            do i = 1, lengath
              aero%numg_a(i,k,m) = aero%num_a(m)%val(ideep(i),k)
              aero%dgnumg(i,k,m) = aero%dgnum(m)%val(ideep(i),k)
            end do
          end do
          do l = 1, aero%nspec(m)
            do k = 1, pver
              do i = 1, lengath
                aero%mmrg_a(i,k,l,m) = aero%mmr_a(l,m)%val(ideep(i),k)
              end do
            end do
          end do
        end do
      else if (aero%scheme == 'bulk') then
        do m = 1, aero%nbulk
          do k = 1, pver
            do i = 1, lengath
              aero%mmrg_bulk(i,k,m) = aero%mmr_bulk(m)%val(ideep(i),k)
            end do
          end do
        end do
      end if
    end if

    do i = 1, lengath
      zfg(i,pver+1) = zf(ideep(i),pver+1)
    end do
    do i = 1, lengath
      capeg    (i) = cape    (ideep(i))
      lclg     (i) = lcl     (ideep(i))
      lelg     (i) = lel     (ideep(i))
      maxg     (i) = maxi    (ideep(i))
      tlg      (i) = tl      (ideep(i))
      landfracg(i) = landfrac(ideep(i))
    end do
    !
    ! Calculate sub-cloud layer pressure "thickness" for use in
    ! closure and tendency routines.
    !
    do k = msg + 1, pver
      do i = 1, lengath
        if (k >= maxg(i)) then
          dsubcld(i) = dsubcld(i) + dp(i,k)
        end if
      end do
    end do
    !
    ! Define array of factors (alpha) which defines interfacial
    ! values, as well as interfacial values for (q,s) used in
    ! subsequent routines.
    !
    do k = msg + 2, pver
      do i = 1, lengath
        sdifr = 0
        qdifr = 0
        if (sg(i,k) > 0 .or. sg(i,k-1) > 0) &
          sdifr = abs((sg(i,k) - sg(i,k-1)) / max(sg(i,k-1), sg(i,k)))
        if (qg(i,k) > 0 .or. qg(i,k-1) > 0) &
          qdifr = abs((qg(i,k) - qg(i,k-1)) / max(qg(i,k-1), qg(i,k)))
        if (sdifr > 1.E-6_r8) then
          shat(i,k) = log(sg(i,k-1) / sg(i,k)) * sg(i,k-1) * sg(i,k) / (sg(i,k-1) - sg(i,k))
        else
          shat(i,k) = 0.5_r8 * (sg(i,k) + sg(i,k-1))
        end if
        if (qdifr > 1.E-6_r8) then
          qhat(i,k) = log(qg(i,k-1) / qg(i,k)) * qg(i,k-1) * qg(i,k) / (qg(i,k-1) - qg(i,k))
        else
          qhat(i,k) = 0.5_r8 * (qg(i,k) + qg(i,k-1))
        end if
      end do
    end do
    !
    ! Obtain cloud properties.
    !
    call cldprp( &
      lchnk    , &
      qg       , &
      tg       , &
      ug       , &
      vg       , &
      pg       , &
      zg       , &
      sg       , &
      mu       , &
      eu       , &
      du       , &
      md       , &
      ed       , &
      sd       , &
      qd       , &
      mc       , &
      qu       , &
      su       , &
      zfg      , &
      qs       , &
      hmn      , &
      hsat     , &
      shat     , &
      qlg      , &
      cmeg     , &
      maxg     , &
      lelg     , &
      jt       , &
      jlcl     , &
      maxg     , &
      j0       , &
      jd       , &
      rl       , &
      lengath  , &
      rgas     , &
      grav     , &
      cpres    , &
      msg      , &
      pflxg    , &
      evpg     , &
      cug      , &
      rprdg    , &
      limcnv   , &
      landfracg, &
      qldeg    , &
      aero     , &
      loc_conv , &
      qhat     )

    if (zmconv_microp) then
      do i = 1, lengath
        capeg(i) = capeg(i)+ loc_conv%dcape(i)
      end do
    end if

    !
    ! Convert detrainment from units of "1/m" to "1/mb".
    !
    do k = msg + 1, pver
      do i = 1, lengath
        du   (i,k) = du   (i,k) * (zfg(i,k) - zfg(i,k+1)) / dp(i,k)
        eu   (i,k) = eu   (i,k) * (zfg(i,k) - zfg(i,k+1)) / dp(i,k)
        ed   (i,k) = ed   (i,k) * (zfg(i,k) - zfg(i,k+1)) / dp(i,k)
        cug  (i,k) = cug  (i,k) * (zfg(i,k) - zfg(i,k+1)) / dp(i,k)
        cmeg (i,k) = cmeg (i,k) * (zfg(i,k) - zfg(i,k+1)) / dp(i,k)
        rprdg(i,k) = rprdg(i,k) * (zfg(i,k) - zfg(i,k+1)) / dp(i,k)
        evpg (i,k) = evpg (i,k) * (zfg(i,k) - zfg(i,k+1)) / dp(i,k)
      end do
    end do

    if (zmconv_microp) then
      do k = msg + 1, pver
        do i = 1, lengath
          loc_conv%sprd(i,k) = loc_conv%sprd(i,k) * (zfg(i,k) - zfg(i,k+1)) / dp(i,k)
          loc_conv%frz (i,k) = loc_conv%frz (i,k) * (zfg(i,k) - zfg(i,k+1)) / dp(i,k)
        end do
      end do
    end if

    call closure( &
      lchnk     , &
      qg        , &
      tg        , &
      pg        , &
      zg        , &
      sg        , &
      tpg       , &
      qs        , &
      qu        , &
      su        , &
      mc        , &
      du        , &
      mu        , &
      md        , &
      qd        , &
      sd        , &
      qhat      , &
      shat      , &
      dp        , &
      qstpg     , &
      zfg       , &
      qlg       , &
      dsubcld   , &
      mb        , &
      capeg     , &
      tlg       , &
      lclg      , &
      lelg      , &
      jt        , &
      maxg      , &
      1         , &
      lengath   , &
      rgas      , &
      grav      , &
      cpres     , &
      rl        , &
      msg       , &
      capelmt   )
    !
    ! Limit cloud base mass flux to theoretical upper bound.
    !
    do i = 1, lengath
      mumax(i) = 0
    end do
    do k = msg + 2, pver
      do i = 1, lengath
        mumax(i) = max(mumax(i), mu(i,k) / dp(i,k))
      end do
    end do

    do i = 1, lengath
      if (mumax(i) > 0) then
        mb(i) = min(mb(i), 0.5_r8 / (delt * mumax(i)))
      else
        mb(i) = 0
      end if
    end do
    ! If no_deep_pbl = .true., don't allow convection entirely
    ! within PBL (suggestion of Bjorn Stevens, 8-2000)
    if (no_deep_pbl) then
      do i = 1, lengath
        if (zm(ideep(i),jt(i)) < pblh(ideep(i))) mb(i) = 0
      end do
    end if

    if (zmconv_microp) then
      do k = msg + 1, pver
        do i = 1, lengath
          loc_conv%sprd(i,k) = loc_conv%sprd(i,k) * mb(i)
          loc_conv%frz (i,k) = loc_conv%frz (i,k) * mb(i)
        end do
      end do
    end if

    do k = msg + 1, pver
      do i = 1, lengath
        mu   (i,k  ) = mu   (i,k  ) * mb(i)
        md   (i,k  ) = md   (i,k  ) * mb(i)
        mc   (i,k  ) = mc   (i,k  ) * mb(i)
        du   (i,k  ) = du   (i,k  ) * mb(i)
        eu   (i,k  ) = eu   (i,k  ) * mb(i)
        ed   (i,k  ) = ed   (i,k  ) * mb(i)
        cmeg (i,k  ) = cmeg (i,k  ) * mb(i)
        rprdg(i,k  ) = rprdg(i,k  ) * mb(i)
        cug  (i,k  ) = cug  (i,k  ) * mb(i)
        evpg (i,k  ) = evpg (i,k  ) * mb(i)
        pflxg(i,k+1) = pflxg(i,k+1) * mb(i) * 100 / grav

        if (zmconv_microp .and. mb(i) == 0) then
          qlg            (i,k) = 0
          loc_conv%qliq  (i,k) = 0
          loc_conv%qice  (i,k) = 0
          loc_conv%qrain (i,k) = 0
          loc_conv%qsnow (i,k) = 0
          loc_conv%wu    (i,k) = 0
          loc_conv%qnl   (i,k) = 0
          loc_conv%qni   (i,k) = 0
          loc_conv%qnr   (i,k) = 0
          loc_conv%qns   (i,k) = 0
          loc_conv%autolm(i,k) = 0
          loc_conv%accrlm(i,k) = 0
          loc_conv%bergnm(i,k) = 0
          loc_conv%fhtimm(i,k) = 0
          loc_conv%fhtctm(i,k) = 0
          loc_conv%fhmlm (i,k) = 0
          loc_conv%hmpim (i,k) = 0
          loc_conv%accslm(i,k) = 0
          loc_conv%dlfm  (i,k) = 0
          loc_conv%autoln(i,k) = 0
          loc_conv%accrln(i,k) = 0
          loc_conv%bergnn(i,k) = 0
          loc_conv%fhtimn(i,k) = 0
          loc_conv%fhtctn(i,k) = 0
          loc_conv%fhmln (i,k) = 0
          loc_conv%accsln(i,k) = 0
          loc_conv%activn(i,k) = 0
          loc_conv%dlfn  (i,k) = 0
          loc_conv%cmel  (i,k) = 0
          loc_conv%autoim(i,k) = 0
          loc_conv%accsim(i,k) = 0
          loc_conv%difm  (i,k) = 0
          loc_conv%cmei  (i,k) = 0
          loc_conv%nuclin(i,k) = 0
          loc_conv%autoin(i,k) = 0
          loc_conv%accsin(i,k) = 0
          loc_conv%hmpin (i,k) = 0
          loc_conv%difn  (i,k) = 0
          loc_conv%trspcm(i,k) = 0
          loc_conv%trspcn(i,k) = 0
          loc_conv%trspim(i,k) = 0
          loc_conv%trspin(i,k) = 0
        end if
      end do
    end do
    !
    ! Compute temperature and moisture changes due to convection.
    !
    call q1q2_pjr( &
      lchnk      , &
      dqdt       , &
      dsdt       , &
      qg         , &
      qs         , &
      qu         , &
      su         , &
      du         , &
      qhat       , &
      shat       , &
      dp         , &
      mu         , &
      md         , &
      sd         , &
      qd         , &
      qldeg      , &
      dsubcld    , &
      jt         , &
      maxg       , &
      1          , &
      lengath    , &
      cpres      , &
      rl         , &
      msg        , &
      dlg        , &
      evpg       , &
      cug        , &
      loc_conv   )
    !
    ! Gather back temperature and mixing ratio.
    !
    if (zmconv_microp) then
      do k = msg + 1, pver
        do i = 1, lengath
          if (dqdt(i,k) * 2 * delt + qg(i,k) < 0.0_r8) then
            negadq = (dqdt(i,k) + 0.5_r8 * qg(i,k) / delt) / 0.9999_r8
            dqdt(i,k) = dqdt(i,k) - negadq
            do kk = k, jt(i), -1
              if (negadq < 0) then
                if (rprdg(i,kk)> -negadq * dp(i,k) / dp(i,kk)) then
                  dsdt(i,k) = dsdt(i,k) + negadq * rl / cpres
                  if (rprdg(i,kk) > loc_conv%sprd(i,kk)) then
                    if (rprdg(i,kk) - loc_conv%sprd(i,kk) < -negadq * dp(i,k) / dp(i,kk)) then
                      dsdt(i,k) = dsdt(i,k) + (negadq + (rprdg(i,kk) - loc_conv%sprd(i,kk)) * dp(i,kk) / dp(i,k)) * latice / cpres
                      loc_conv%sprd(i,kk) = negadq * dp(i,k) / dp(i,kk) + rprdg(i,kk)
                    end if
                  else
                    loc_conv%sprd(i,kk) = loc_conv%sprd(i,kk) + negadq * dp(i,k) / dp(i,kk)
                    dsdt(i,k) = dsdt(i,k) + negadq * latice / cpres
                  end if
                  rprdg(i,kk) = rprdg(i,kk) + negadq * dp(i,k) / dp(i,kk)
                  negadq = 0
                else
                  negadq = rprdg(i,kk) * dp(i,kk) / dp(i,k) + negadq
                  dsdt(i,k) = dsdt(i,k) - rprdg(i,kk) * rl / cpres * dp(i,kk) / dp(i,k)
                  if (rprdg(i,kk)>loc_conv%sprd(i,kk)) then
                    dsdt(i,k) = dsdt(i,k) - loc_conv%sprd(i,kk) * latice / cpres * dp(i,kk) / dp(i,k)
                    loc_conv%sprd(i,kk) = 0
                  else
                    dsdt(i,k) = dsdt(i,k) - rprdg(i,kk) * latice / cpres * dp(i,kk) / dp(i,k)
                    loc_conv%sprd(i,kk)= loc_conv%sprd(i,kk) - rprdg(i,kk)
                  end if
                  rprdg(i,kk) = 0
                end if

                if (dlg(i,kk) > loc_conv%di(i,kk)) then
                  doliq = .true.
                else
                  doliq = .false.
                end if

                if (negadq < 0) then
                  if (doliq) then
                    if (dlg(i,kk) > -negadq * dp(i,k) / dp(i,kk)) then
                      dsdt(i,k) = dsdt(i,k) + negadq * rl / cpres
                      loc_conv%dnl(i,kk) = loc_conv%dnl(i,kk) * (1 + negadq * dp(i,k) / dp(i,kk) / dlg(i,kk))
                      dlg(i,kk)  = dlg(i,kk) + negadq * dp(i,k) / dp(i,kk)
                      negadq = 0
                    else
                      negadq = negadq + dlg(i,kk) * dp(i,kk) / dp(i,k)
                      dsdt(i,k) = dsdt(i,k) - dlg(i,kk) * dp(i,kk) / dp(i,k) * rl / cpres
                      dlg(i,kk) = 0
                      loc_conv%dnl(i,kk) = 0
                    end if
                  else
                    if (loc_conv%di(i,kk) > -negadq * dp(i,k) / dp(i,kk)) then
                      dsdt(i,k) = dsdt(i,k) + negadq * (rl + latice) / cpres
                      loc_conv%dni(i,kk) = loc_conv%dni(i,kk) * (1 + negadq * dp(i,k) / dp(i,kk) / loc_conv%di(i,kk))
                      loc_conv%di(i,kk)  = loc_conv%di(i,kk) + negadq * dp(i,k) / dp(i,kk)
                      negadq = 0
                    else
                      negadq = negadq + loc_conv%di(i,kk) * dp(i,kk) / dp(i,k)
                      dsdt(i,k) = dsdt(i,k) - loc_conv%di(i,kk) * dp(i,kk) / dp(i,k) * (rl + latice) / cpres
                      loc_conv%di(i,kk) = 0
                      loc_conv%dni(i,kk) = 0
                    end if
                    doliq = .false.
                  end if
                end if
                if (negadq < 0 .and. doliq) then
                  if (dlg(i,kk) > -negadq * dp(i,k) / dp(i,kk)) then
                    dsdt(i,k) = dsdt(i,k) + negadq * rl / cpres
                    loc_conv%dnl(i,kk) = loc_conv%dnl(i,kk) * (1 + negadq * dp(i,k) / dp(i,kk) / dlg(i,kk))
                    dlg(i,kk)  = dlg(i,kk) + negadq * dp(i,k) / dp(i,kk)
                    negadq = 0
                  else
                    negadq = negadq + dlg(i,kk) * dp(i,kk) / dp(i,k)
                    dsdt(i,k) = dsdt(i,k) - dlg(i,kk) * dp(i,kk) / dp(i,k) * rl / cpres
                    dlg(i,kk) = 0
                    loc_conv%dnl(i,kk) = 0
                  end if
                end if
              end if
            end do
            if (negadq < 0) then
              dqdt(i,k) = dqdt(i,k) + negadq
            end if
          end if
        end do
      end do
    end if

    do k = msg + 1, pver
      do i = 1, lengath
        !
        ! q is updated to compute net precip.
        !
        q(ideep(i),k) = qh(ideep(i),k) + 2._r8*delt*dqdt(i,k)
        qtnd(ideep(i),k) = dqdt (i,k)
        cme (ideep(i),k) = cmeg (i,k)
        rprd(ideep(i),k) = rprdg(i,k)
        zdu (ideep(i),k) = du   (i,k)
        mcon(ideep(i),k) = mc   (i,k)
        heat(ideep(i),k) = dsdt (i,k) * cpres
        dlf (ideep(i),k) = dlg  (i,k)
        pflx(ideep(i),k) = pflxg(i,k)
        ql  (ideep(i),k) = qlg  (i,k)
      end do
    end do

    if (zmconv_microp) then
      do k = msg + 1, pver
        do i = 1, lengath
          dif            (ideep(i),k) = loc_conv%di        (i,k)
          dnlf           (ideep(i),k) = loc_conv%dnl       (i,k)
          dnif           (ideep(i),k) = loc_conv%dni       (i,k)
          conv%qi        (ideep(i),k) = loc_conv%qice      (i,k)
          conv%frz       (ideep(i),k) = loc_conv%frz       (i,k) * latice / cpres
          conv%sprd      (ideep(i),k) = loc_conv%sprd      (i,k)
          conv%wu        (ideep(i),k) = loc_conv%wu        (i,k)
          conv%qliq      (ideep(i),k) = loc_conv%qliq      (i,k)
          conv%qice      (ideep(i),k) = loc_conv%qice      (i,k)
          conv%qrain     (ideep(i),k) = loc_conv%qrain     (i,k)
          conv%qsnow     (ideep(i),k) = loc_conv%qsnow     (i,k)
          conv%qnl       (ideep(i),k) = loc_conv%qnl       (i,k)
          conv%qni       (ideep(i),k) = loc_conv%qni       (i,k)
          conv%qnr       (ideep(i),k) = loc_conv%qnr       (i,k)
          conv%qns       (ideep(i),k) = loc_conv%qns       (i,k)
          conv%autolm    (ideep(i),k) = loc_conv%autolm    (i,k)
          conv%accrlm    (ideep(i),k) = loc_conv%accrlm    (i,k)
          conv%bergnm    (ideep(i),k) = loc_conv%bergnm    (i,k)
          conv%fhtimm    (ideep(i),k) = loc_conv%fhtimm    (i,k)
          conv%fhtctm    (ideep(i),k) = loc_conv%fhtctm    (i,k)
          conv%fhmlm     (ideep(i),k) = loc_conv%fhmlm     (i,k)
          conv%hmpim     (ideep(i),k) = loc_conv%hmpim     (i,k)
          conv%accslm    (ideep(i),k) = loc_conv%accslm    (i,k)
          conv%dlfm      (ideep(i),k) = loc_conv%dlfm      (i,k)
          conv%autoln    (ideep(i),k) = loc_conv%autoln    (i,k)
          conv%accrln    (ideep(i),k) = loc_conv%accrln    (i,k)
          conv%bergnn    (ideep(i),k) = loc_conv%bergnn    (i,k)
          conv%fhtimn    (ideep(i),k) = loc_conv%fhtimn    (i,k)
          conv%fhtctn    (ideep(i),k) = loc_conv%fhtctn    (i,k)
          conv%fhmln     (ideep(i),k) = loc_conv%fhmln     (i,k)
          conv%accsln    (ideep(i),k) = loc_conv%accsln    (i,k)
          conv%activn    (ideep(i),k) = loc_conv%activn    (i,k)
          conv%dlfn      (ideep(i),k) = loc_conv%dlfn      (i,k)
          conv%cmel      (ideep(i),k) = loc_conv%cmel      (i,k)
          conv%autoim    (ideep(i),k) = loc_conv%autoim    (i,k)
          conv%accsim    (ideep(i),k) = loc_conv%accsim    (i,k)
          conv%difm      (ideep(i),k) = loc_conv%difm      (i,k)
          conv%cmei      (ideep(i),k) = loc_conv%cmei      (i,k)
          conv%nuclin    (ideep(i),k) = loc_conv%nuclin    (i,k)
          conv%autoin    (ideep(i),k) = loc_conv%autoin    (i,k)
          conv%accsin    (ideep(i),k) = loc_conv%accsin    (i,k)
          conv%hmpin     (ideep(i),k) = loc_conv%hmpin     (i,k)
          conv%difn      (ideep(i),k) = loc_conv%difn      (i,k)
          conv%trspcm    (ideep(i),k) = loc_conv%trspcm    (i,k)
          conv%trspcn    (ideep(i),k) = loc_conv%trspcn    (i,k)
          conv%trspim    (ideep(i),k) = loc_conv%trspim    (i,k)
          conv%trspin    (ideep(i),k) = loc_conv%trspin    (i,k)
          conv%lambdadpcu(ideep(i),k) = loc_conv%lambdadpcu(i,k)
          conv%mudpcu    (ideep(i),k) = loc_conv%mudpcu    (i,k)
        end do
      end do
      do k = msg + 1, pver
        do i = 1, ncol
          ! Convert it from units of "kg/kg" to "g/m3"
          if (k < pver) then
            conv%qice (i,k) = 0.5_r8 * (conv%qice (i,k) + conv%qice (i,k+1))
            conv%qliq (i,k) = 0.5_r8 * (conv%qliq (i,k) + conv%qliq (i,k+1))
            conv%qrain(i,k) = 0.5_r8 * (conv%qrain(i,k) + conv%qrain(i,k+1))
            conv%qsnow(i,k) = 0.5_r8 * (conv%qsnow(i,k) + conv%qsnow(i,k+1))
            conv%qni  (i,k) = 0.5_r8 * (conv%qni  (i,k) + conv%qni  (i,k+1))
            conv%qnl  (i,k) = 0.5_r8 * (conv%qnl  (i,k) + conv%qnl  (i,k+1))
            conv%qnr  (i,k) = 0.5_r8 * (conv%qnr  (i,k) + conv%qnr  (i,k+1))
            conv%qns  (i,k) = 0.5_r8 * (conv%qns  (i,k) + conv%qns  (i,k+1))
            conv%wu   (i,k) = 0.5_r8 * (conv%wu   (i,k) + conv%wu   (i,k+1))
          end if
          if (t(i,k) > 273.15_r8 .and. t(i,k-1) <= 273.15_r8) then
            conv%qice (i,k-1) = conv%qice (i,k-1) + conv%qice (i,k)
            conv%qice (i,k  ) = 0
            conv%qni  (i,k-1) = conv%qni  (i,k-1) + conv%qni  (i,k)
            conv%qni  (i,k  ) = 0
            conv%qsnow(i,k-1) = conv%qsnow(i,k-1) + conv%qsnow(i,k)
            conv%qsnow(i,k  ) = 0
            conv%qns  (i,k-1) = conv%qns  (i,k-1) + conv%qns  (i,k)
            conv%qns  (i,k  ) = 0
          end if
          conv%qice (i,k) = conv%qice (i,k) * pap(i,k) / t(i,k) / rgas * 1000
          conv%qliq (i,k) = conv%qliq (i,k) * pap(i,k) / t(i,k) / rgas * 1000
          conv%qrain(i,k) = conv%qrain(i,k) * pap(i,k) / t(i,k) / rgas * 1000
          conv%qsnow(i,k) = conv%qsnow(i,k) * pap(i,k) / t(i,k) / rgas * 1000
          conv%qni  (i,k) = conv%qni  (i,k) * pap(i,k) / t(i,k) / rgas
          conv%qnl  (i,k) = conv%qnl  (i,k) * pap(i,k) / t(i,k) / rgas
          conv%qnr  (i,k) = conv%qnr  (i,k) * pap(i,k) / t(i,k) / rgas
          conv%qns  (i,k) = conv%qns  (i,k) * pap(i,k) / t(i,k) / rgas
        end do
      end do
    end if

    do i = 1, lengath
      jctop(ideep(i)) = jt(i)
      jcbot(ideep(i)) = maxg(i)
      pflx (ideep(i),pverp) = pflxg(i,pverp)
    end do

    if (zmconv_microp) then
      do i = 1, lengath
        conv%dcape(ideep(i)) = loc_conv%dcape(i)
      end do
    end if

    ! Compute precipitation by integrating change in water vapor minus detrained cloud water
    do k = pver, msg + 1, -1
      do i = 1, ncol
        prec(i) = prec(i) - dpp(i,k)* (q(i,k)-qh(i,k)) - dpp(i,k)*(dlf(i,k)+dif(i,k))*2._r8*delt
      end do
    end do

    ! Obtain final precipitation rate in m/s.
    do i = 1, ncol
      prec(i) = rgrav * max(prec(i), 0.0_r8) / (2.0_r8 * delt) / 1000.0_r8
    end do

    ! Compute reserved liquid (not yet in cldliq) for energy integrals.
    ! Treat rliq as flux out bottom, to be added back later.
    do k = 1, pver
      do i = 1, ncol
        rliq(i) = rliq(i) + (dlf(i,k) + dif(i,k)) * dpp(i,k) / gravit
        rice(i) = rice(i) + dif(i,k) * dpp(i,k) / gravit
      end do
    end do
    rliq(:ncol) = rliq(:ncol) / 1000.0_r8
    rice(:ncol) = rice(:ncol) / 1000.0_r8

    if (zmconv_microp) then
      deallocate(loc_conv%frz       )
      deallocate(loc_conv%sprd      )
      deallocate(loc_conv%wu        )
      deallocate(loc_conv%qi        )
      deallocate(loc_conv%qliq      )
      deallocate(loc_conv%qice      )
      deallocate(loc_conv%qrain     )
      deallocate(loc_conv%qsnow     )
      deallocate(loc_conv%di        )
      deallocate(loc_conv%dnl       )
      deallocate(loc_conv%dni       )
      deallocate(loc_conv%qnl       )
      deallocate(loc_conv%qni       )
      deallocate(loc_conv%qnr       )
      deallocate(loc_conv%qns       )
      deallocate(loc_conv%qide      )
      deallocate(loc_conv%qncde     )
      deallocate(loc_conv%qnide     )
      deallocate(loc_conv%autolm    )
      deallocate(loc_conv%accrlm    )
      deallocate(loc_conv%bergnm    )
      deallocate(loc_conv%fhtimm    )
      deallocate(loc_conv%fhtctm    )
      deallocate(loc_conv%fhmlm     )
      deallocate(loc_conv%hmpim     )
      deallocate(loc_conv%accslm    )
      deallocate(loc_conv%dlfm      )
      deallocate(loc_conv%cmel      )
      deallocate(loc_conv%autoln    )
      deallocate(loc_conv%accrln    )
      deallocate(loc_conv%bergnn    )
      deallocate(loc_conv%fhtimn    )
      deallocate(loc_conv%fhtctn    )
      deallocate(loc_conv%fhmln     )
      deallocate(loc_conv%accsln    )
      deallocate(loc_conv%activn    )
      deallocate(loc_conv%dlfn      )
      deallocate(loc_conv%autoim    )
      deallocate(loc_conv%accsim    )
      deallocate(loc_conv%difm      )
      deallocate(loc_conv%cmei      )
      deallocate(loc_conv%nuclin    )
      deallocate(loc_conv%autoin    )
      deallocate(loc_conv%accsin    )
      deallocate(loc_conv%hmpin     )
      deallocate(loc_conv%difn      )
      deallocate(loc_conv%trspcm    )
      deallocate(loc_conv%trspcn    )
      deallocate(loc_conv%trspim    )
      deallocate(loc_conv%trspin    )
      deallocate(loc_conv%lambdadpcu)
      deallocate(loc_conv%mudpcu    )
      deallocate(loc_conv%dcape     )
    end if

  end subroutine zm_convr

  subroutine zm_conv_evap( &
    ncol                 , &
    lchnk                , &
    t                    , &
    pmid                 , &
    pdel                 , &
    q                    , &
    landfrac             , &
    tend_s               , &
    tend_s_snwprd        , &
    tend_s_snwevmlt      , &
    tend_q               , &
    prdprec              , &
    cldfrc               , &
    deltat               , &
    prec                 , &
    snow                 , &
    ntprprd              , &
    ntsnprd              , &
    flxprec              , &
    flxsnow              , &
    prdsnow              )


    ! --------------------------------------------------------------------------
    ! - Compute tendencies due to evaporation of rain from ZM scheme
    ! - Compute the total precipitation and snow fluxes at the surface.
    ! - Add in the latent heat of fusion for snow formation and melt, since it
    !   not dealt with in the Zhang-MacFarlane parameterization.
    ! - Evaporate some of the precip directly into the environment using a
    !   Sundqvist type algorithm
    ! --------------------------------------------------------------------------

    use wv_saturation, only: qsat
    use phys_grid    , only: get_rlat_all_p

    integer , intent(in   ) :: ncol, lchnk
    real(r8), intent(in   ), dimension(pcols,pver ) :: t               ! Temperature (K)
    real(r8), intent(in   ), dimension(pcols,pver ) :: pmid            ! Midpoint pressure (Pa)
    real(r8), intent(in   ), dimension(pcols,pver ) :: pdel            ! Layer thickness (Pa)
    real(r8), intent(in   ), dimension(pcols,pver ) :: q               ! Water vapor (kg/kg)
    real(r8), intent(in   ), dimension(pcols      ) :: landfrac
    real(r8), intent(inout), dimension(pcols,pver ) :: tend_s          ! Heating rate (J/kg/s)
    real(r8), intent(inout), dimension(pcols,pver ) :: tend_q          ! Water vapor tendency (kg/kg/s)
    real(r8), intent(  out), dimension(pcols,pver ) :: tend_s_snwprd   ! Heating rate of snow production
    real(r8), intent(  out), dimension(pcols,pver ) :: tend_s_snwevmlt ! Heating rate of evap/melting of snow
    real(r8), intent(in   ), dimension(pcols,pver ) :: prdprec         ! Precipitation production (kg/ks/s)
    real(r8), intent(in   ), dimension(pcols,pver ) :: cldfrc          ! Cloud fraction
    real(r8), intent(in   ) :: deltat                                  ! Time step
    real(r8), intent(inout), dimension(pcols      ) :: prec            ! Convective-scale preciptn rate
    real(r8), intent(  out), dimension(pcols      ) :: snow            ! Convective-scale snowfall rate
    real(r8), intent(in   ), allocatable, optional  :: prdsnow(:,:)    ! Snow production (kg/ks/s)
    real(r8), intent(  out), dimension(pcols,pverp) :: flxprec         ! Convective-scale flux of precip at interfaces (kg/m2/s)
    real(r8), intent(  out), dimension(pcols,pverp) :: flxsnow         ! Convective-scale flux of snow   at interfaces (kg/m2/s)
    real(r8), intent(  out), dimension(pcols,pver ) :: ntprprd         ! net precip production in layer
    real(r8), intent(  out), dimension(pcols,pver ) :: ntsnprd         ! net snow production in layer

    real(r8), dimension(pcols,pver) :: es         ! Saturation vapor pressure
    real(r8), dimension(pcols,pver) :: fice       ! Ice fraction in precip production
    real(r8), dimension(pcols,pver) :: fsnow_conv ! Snow fraction in precip production
    real(r8), dimension(pcols,pver) :: qs         ! Saturation specific humidity
    real(r8), dimension(pcols     ) :: evpvint    ! Vertical integral of evaporation
    real(r8), dimension(pcols     ) :: evpprec    ! Evaporation of precipitation (kg/kg/s)
    real(r8), dimension(pcols     ) :: evpsnow    ! Evaporation of snowfall (kg/kg/s)
    real(r8), dimension(pcols     ) :: snowmlt    ! Snow melt tendency in layer
    real(r8), dimension(pcols     ) :: flxsntm    ! Flux of snow into layer, after melting
    real(r8), dimension(pcols     ) :: rlat

    real(r8) work1
    real(r8) work2
    real(r8) kemask
    real(r8) evplimit
    real(r8) dum
    real(r8) omsm

    integer i, k
    logical old_snow

    ! If prdsnow is passed in and allocated, then use it in the calculation, otherwise
    ! use the old snow calculation
    old_snow = .true.
    if (present(prdsnow)) then
      if (allocated(prdsnow)) then
        old_snow = .false.
      end if
    end if

    ! Convert input precip to kg/m2/s
    prec(:ncol) = prec(:ncol) * 1000

    ! Determine saturation vapor pressure
    call qsat(t(:ncol,:pver), pmid(:ncol,:pver), es(:ncol,:pver), qs(:ncol,:pver))

    ! Determine ice fraction in rain production (use cloud water parameterization fraction at present)
    call cldfrc_fice(ncol, t, fice, fsnow_conv)

    ! Zero the flux integrals on the top boundary
    flxprec(:ncol,1) = 0
    flxsnow(:ncol,1) = 0
    evpvint(:ncol)   = 0
    omsm = 0.9999_r8

    do k = 1, pver
      do i = 1, ncol
        ! Melt snow falling into layer, if necessary.
        if (old_snow) then
          if (t(i,k) > tmelt) then
            flxsntm(i) = 0
            snowmlt(i) = flxsnow(i,k) * gravit / pdel(i,k)
          else
            flxsntm(i) = flxsnow(i,k)
            snowmlt(i) = 0
          end if
        else
          ! make sure melting snow doesn't reduce temperature below threshold
          if (t(i,k) > tmelt) then
            dum = -latice / cpres * flxsnow(i,k) * gravit / pdel(i,k) * deltat
            if (t(i,k) + dum <= tmelt) then
              dum = (t(i,k) - tmelt) * cpres / latice / deltat
              dum = dum / (flxsnow(i,k) * gravit / pdel(i,k))
              dum = max(0.0_r8, dum)
              dum = min(1.0_r8, dum)
            else
              dum = 1
            end if
            dum = dum * omsm
            flxsntm(i) = flxsnow(i,k) * (1 - dum)
            snowmlt(i) = dum * flxsnow(i,k) * gravit / pdel(i,k)
          else
            flxsntm(i) = flxsnow(i,k)
            snowmlt(i) = 0
          end if
        end if

        ! Relative humidity depression must be > 0 for evaporation
        evplimit = max(1.0_r8 - q(i,k) / qs(i,k), 0.0_r8)

        if (zm_org) then
          kemask = ke * (1 - landfrac(i)) + ke_lnd * landfrac(i)
        else
          kemask = ke
        end if

        ! Total evaporation depends on flux in the top of the layer flux prec
        ! is the net production above layer minus evaporation into environmet
        evpprec(i) = kemask * (1 - cldfrc(i,k)) * evplimit * sqrt(flxprec(i,k))
        !**********************************************************
        ! evpprec(i) = 0.    ! turn off evaporation for now
        !**********************************************************

        ! Don't let evaporation supersaturate layer (approx). Layer may already be saturated.
        ! Currently does not include heating/cooling change to qs
        evplimit   = max(0.0_r8, (qs(i,k) - q(i,k)) / deltat)

        ! Don't evaporate more than is falling into the layer - do not evaporate rain formed
        ! in this layer but if precip production is negative, remove from the available precip
        ! Negative precip production occurs because of evaporation in downdrafts.
        ! evplimit   = flxprec(i,k) * gravit / pdel(i,k) + min(prdprec(i,k), 0.)
        evplimit   = min(evplimit, flxprec(i,k) * gravit / pdel(i,k))

        ! Total evaporation cannot exceed input precipitation
        evplimit   = min(evplimit, (prec(i) - evpvint(i)) * gravit / pdel(i,k))

        evpprec(i) = min(evplimit, evpprec(i))
        if (.not. old_snow) then
          evpprec(i) = max(0.0_r8, evpprec(i))
          evpprec(i) = evpprec(i) * omsm
        end if

        ! Evaporation of snow depends on snow fraction of total precipitation in the top after melting
        if (flxprec(i,k) > 0) then
          ! evpsnow(i) = evpprec(i) * flxsntm(i) / flxprec(i,k)
          ! Prevent roundoff problems
          work1 = min(max(0.0_r8, flxsntm(i) / flxprec(i,k)), 1.0_r8)
          evpsnow(i) = evpprec(i) * work1
        else
          evpsnow(i) = 0
        end if

        ! Vertically integrated evaporation
        evpvint(i) = evpvint(i) + evpprec(i) * pdel(i,k)/gravit

        ! Net precip production is production - evaporation
        ntprprd(i,k) = prdprec(i,k) - evpprec(i)
        ! net snow production is precip production * ice fraction - evaporation - melting
        ! pjrworks ntsnprd(i,k) = prdprec(i,k)*fice(i,k) - evpsnow(i) - snowmlt(i)
        ! pjrwrks2 ntsnprd(i,k) = prdprec(i,k)*fsnow_conv(i,k) - evpsnow(i) - snowmlt(i)
        ! the small amount added to flxprec in the work1 expression has been increased from
        ! 1e-36 to 8.64e-11 (1e-5 mm/day).  This causes the temperature based partitioning
        ! scheme to be used for small flxprec amounts.  This is to address error growth problems.

        if (old_snow) then
#ifdef PERGRO
          work1 = min(max(0.0_r8, flxsnow(i,k) / (flxprec(i,k) + 8.64e-11_r8)), 1.0_r8)
#else
          if (flxprec(i,k) > 0) then
            work1 = min(max(0.0_r8, flxsnow(i,k) / flxprec(i,k)), 1.0_r8)
          else
            work1 = 0
          end if
#endif
          work2 = max(fsnow_conv(i,k), work1)
          if (snowmlt(i) > 0) work2 = 0
          ntsnprd(i,k) = prdprec(i,k) * work2 - evpsnow(i) - snowmlt(i)
          tend_s_snwprd  (i,k) = prdprec(i,k) * work2 * latice
          tend_s_snwevmlt(i,k) = -(evpsnow(i) + snowmlt(i)) * latice
        else
          ntsnprd(i,k) = prdsnow(i,k) - min(flxsnow(i,k) * gravit / pdel(i,k), evpsnow(i) + snowmlt(i))
          tend_s_snwprd  (i,k) = prdsnow(i,k) * latice
          tend_s_snwevmlt(i,k) = -min(flxsnow(i,k) * gravit / pdel(i,k), evpsnow(i) + snowmlt(i)) * latice
        end if

        ! Precipitation fluxes
        flxprec(i,k+1) = flxprec(i,k) + ntprprd(i,k) * pdel(i,k) / gravit
        flxsnow(i,k+1) = flxsnow(i,k) + ntsnprd(i,k) * pdel(i,k) / gravit

        ! Protect against rounding error
        flxprec(i,k+1) = max(flxprec(i,k+1), 0.0_r8)
        flxsnow(i,k+1) = max(flxsnow(i,k+1), 0.0_r8)
        ! More protection (pjr)
        ! flxsnow(i,k+1) = min(flxsnow(i,k+1), flxprec(i,k+1))

        ! Heating (cooling) and moistening due to evaporation
        ! - latent heat of vaporization for precip production has already been accounted for
        ! - snow is contained in prec
        if (old_snow) then
          tend_s(i,k) = -evpprec(i) * latvap + ntsnprd(i,k) * latice
        else
          tend_s(i,k) = -evpprec(i) * latvap + tend_s_snwevmlt(i,k)
        end if
        tend_q(i,k) = evpprec(i)
      end do
    end do

    ! Set output precipitation rates (m/s)
    prec(:ncol) = flxprec(:ncol,pver+1) / 1000.0_r8
    snow(:ncol) = flxsnow(:ncol,pver+1) / 1000.0_r8

  end subroutine zm_conv_evap

  subroutine convtran( &
    lchnk            , &
    doconvtran       , &
    q                , &
    ncnst            , &
    mu               , &
    md               , &
    du               , &
    eu               , &
    ed               , &
    dp               , &
    dsubcld          , &
    jt               , &
    mx               , &
    ideep            , &
    il1g             , &
    il2g             , &
    nstep            , &
    fracis           , &
    dqdt             , &
    dpdry            , &
    dt               )
    ! --------------------------------------------------------------------------
    !
    ! Purpose:
    ! Convective transport of trace species
    !
    ! Mixing ratios may be with respect to either dry or moist air
    !
    ! Author: P. Rasch
    !
    ! --------------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use constituents, only: cnst_get_type_byind
    use ppgrid

    integer , intent(in) :: lchnk                                   ! Chunk identifier
    integer , intent(in) :: ncnst                                   ! Number of tracers to transport
    logical , intent(in), dimension(           ncnst) :: doconvtran ! Flag for doing convective transport
    real(r8), intent(in), dimension(pcols,pver,ncnst) :: q          ! Tracer array including moisture
    real(r8), intent(in), dimension(pcols,pver      ) :: mu         ! Mass flux up
    real(r8), intent(in), dimension(pcols,pver      ) :: md         ! Mass flux down
    real(r8), intent(in), dimension(pcols,pver      ) :: du         ! Mass detraining from updraft
    real(r8), intent(in), dimension(pcols,pver      ) :: eu         ! Mass entraining from updraft
    real(r8), intent(in), dimension(pcols,pver      ) :: ed         ! Mass entraining from downdraft
    real(r8), intent(in), dimension(pcols,pver      ) :: dp         ! Delta pressure between interfaces
    real(r8), intent(in), dimension(pcols           ) :: dsubcld    ! Delta pressure from cloud base to sfc
    real(r8), intent(in), dimension(pcols,pver,ncnst) :: fracis     ! Fraction of tracer that is insoluble
    integer , intent(in), dimension(pcols) :: jt                    ! Index of cloud top for each column
    integer , intent(in), dimension(pcols) :: mx                    ! Index of cloud top for each column
    integer , intent(in), dimension(pcols) :: ideep                 ! Gathering array
    integer , intent(in) :: il1g                                    ! Gathered min lon indices over which to operate
    integer , intent(in) :: il2g                                    ! Gathered max lon indices over which to operate
    integer , intent(in) :: nstep                                   ! Time step index
    real(r8), intent(in), dimension(pcols,pver) :: dpdry            ! Delta pressure between interfaces
    real(r8), intent(in) :: dt                                      ! 2 delta t (model time increment)

    real(r8), intent(out), dimension(pcols,pver,ncnst) :: dqdt  ! Tracer tendency array

    integer i                 ! Work index
    integer k                 ! Work index
    integer kbm               ! Highest altitude index of cloud base
    integer kk                ! Work index
    integer kkp1              ! Work index
    integer km1               ! Work index
    integer kp1               ! Work index
    integer ktm               ! Highest altitude index of cloud top
    integer m                 ! Work index

    real(r8) cabv                 ! Mix ratio of constituent above
    real(r8) cbel                 ! Mix ratio of constituent below
    real(r8) cdifr                ! Normalized diff between cabv and cbel
    real(r8) chat  (pcols,pver)   ! Mix ratio in env at interfaces
    real(r8) cond  (pcols,pver)   ! Mix ratio in downdraft at interfaces
    real(r8) const (pcols,pver)   ! Gathered tracer array
    real(r8) fisg  (pcols,pver)   ! gathered insoluble fraction of tracer
    real(r8) conu  (pcols,pver)   ! Mix ratio in updraft at interfaces
    real(r8) dcondt(pcols,pver)   ! Gathered tend array
    real(r8) small                ! A small number
    real(r8) mbsth                ! Threshold for mass fluxes
    real(r8) mupdudp              ! A work variable
    real(r8) minc                 ! A work variable
    real(r8) maxc                 ! A work variable
    real(r8) fluxin               ! A work variable
    real(r8) fluxout              ! A work variable
    real(r8) netflux              ! A work variable

    real(r8) dutmp(pcols,pver)    ! Mass detraining from updraft
    real(r8) eutmp(pcols,pver)    ! Mass entraining from updraft
    real(r8) edtmp(pcols,pver)    ! Mass entraining from downdraft
    real(r8) dptmp(pcols,pver)    ! Delta pressure between interfaces
    real(r8) total(pcols)
    real(r8) negadt,qtmp

    small = 1.0e-36_r8
    ! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
    mbsth = 1.0e-15_r8

    ! Find the highest level top and bottom levels of convection
    ktm = pver
    kbm = pver
    do i = il1g, il2g
      ktm = min(ktm, jt(i))
      kbm = min(kbm, mx(i))
    end do

    ! Loop ever each constituent
    do m = 2, ncnst
      if (doconvtran(m)) then
        if (cnst_get_type_byind(m) == 'dry') then
          do k = 1, pver
            do i = il1g, il2g
              dptmp(i,k) = dpdry(i,k)
              dutmp(i,k) = du(i,k) * dp(i,k) / dpdry(i,k)
              eutmp(i,k) = eu(i,k) * dp(i,k) / dpdry(i,k)
              edtmp(i,k) = ed(i,k) * dp(i,k) / dpdry(i,k)
            end do
          end do
        else
          do k = 1, pver
            do i = il1g, il2g
              dptmp(i,k) = dp(i,k)
              dutmp(i,k) = du(i,k)
              eutmp(i,k) = eu(i,k)
              edtmp(i,k) = ed(i,k)
            end do
          end do
        end if
        ! Gather up the constituent and set tend to zero
        do k = 1, pver
          do i = il1g, il2g
            const(i,k) = q(ideep(i),k,m)
            fisg (i,k) = fracis(ideep(i),k,m)
          end do
        end do

        ! From now on work only with gathered data

        ! Interpolate environment tracer values to interfaces
        do k = 1, pver
          km1 = max(1, k - 1)
          do i = il1g, il2g
            minc = min(const(i,km1), const(i,k))
            maxc = max(const(i,km1), const(i,k))
            if (minc < 0) then
              cdifr = 0
            else
              cdifr = abs(const(i,k) - const(i,km1)) / max(maxc, small)
            end if

            ! If the two layers differ significantly use a geometric averaging
            ! procedure
            if (cdifr > 1.0e-6_r8) then
              cabv = max(const(i,km1), maxc * 1.0e-12_r8)
              cbel = max(const(i,k  ), maxc * 1.0e-12_r8)
              chat(i,k) = log(cabv / cbel) / (cabv - cbel) * cabv * cbel
            else             ! Small diff, so just arithmetic mean
              chat(i,k) = 0.5_r8 * (const(i,k) + const(i,km1))
            end if

            ! Provisional up and down draft values
            conu(i,k) = chat(i,k)
            cond(i,k) = chat(i,k)
            dcondt(i,k) = 0
          end do
        end do

        ! Do levels adjacent to top and bottom
        k   = 2
        km1 = 1
        kk  = pver
        do i = il1g, il2g
          mupdudp = mu(i,kk) + dutmp(i,kk) * dptmp(i,kk)
          if (mupdudp > mbsth) then
            conu(i,kk) = (+eutmp(i,kk ) * fisg(i,kk ) * const(i,kk ) * dptmp(i,kk )) / mupdudp
          end if
          if (md(i,k) < -mbsth) then
            cond(i,k ) = (-edtmp(i,km1) * fisg(i,km1) * const(i,km1) * dptmp(i,km1)) / md(i,k)
          end if
        end do
        ! Updraft from bottom to top
        do kk = pver - 1, 1, -1
          kkp1 = min(pver, kk + 1)
          do i = il1g, il2g
            mupdudp = mu(i,kk) + dutmp(i,kk) * dptmp(i,kk)
            if (mupdudp > mbsth) then
              conu(i,kk) = (mu(i,kkp1) * conu(i,kkp1) + eutmp(i,kk) * fisg(i,kk) * const(i,kk) * dptmp(i,kk)) / mupdudp
            end if
          end do
        end do
        ! Downdraft from top to bottom
        do k = 3, pver
          km1 = max(1, k - 1)
          do i = il1g, il2g
            if (md(i,k) < -mbsth) then
              cond(i,k) = (md(i,km1) * cond(i,km1) - edtmp(i,km1) * fisg(i,km1) * const(i,km1) * dptmp(i,km1)) / md(i,k)
            end if
          end do
        end do

        do k = ktm, pver
          km1 = max(1   , k - 1)
          kp1 = min(pver, k + 1)
          do i = il1g, il2g
            ! version 1 hard to check for roundoff errors
            !               dcondt(i,k) =
            !     $                  +(+mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))
            !     $                    -mu(i,k)*   (conu(i,k)-chat(i,k))
            !     $                    +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))
            !     $                    -md(i,k)*   (cond(i,k)-chat(i,k))
            !     $                   )/dp(i,k)

            ! version 2 hard to limit fluxes
            !               fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k)
            !     $                 -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
            !               fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)
            !     $                 -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))

            ! version 3 limit fluxes outside convection to mass in appropriate layer
            ! these limiters are probably only safe for positive definite quantitities
            ! it assumes that mu and md already satify a courant number limit of 1
            fluxin =  mu(i,kp1)*conu(i,kp1)+ mu(i,k)*min(chat(i,k),const(i,km1)) &
                   -(md(i,k)  *cond(i,k) + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
            fluxout = mu(i,k)*conu(i,k) + mu(i,kp1)*min(chat(i,kp1),const(i,k)) &
                    -(md(i,kp1)*cond(i,kp1) + md(i,k)*min(chat(i,k),const(i,k)))

            netflux = fluxin - fluxout
            if (abs(netflux) < max(fluxin, fluxout) * 1.0e-12_r8) then
              netflux = 0
            end if
            dcondt(i,k) = netflux / dptmp(i,k)
          end do
        end do

        do k = kbm, pver
          km1 = max(1, k - 1)
          do i = il1g, il2g
            if (k == mx(i)) then
              ! version 1
              !                  dcondt(i,k) = (1./dsubcld(i))*
              !     $              (-mu(i,k)*(conu(i,k)-chat(i,k))
              !     $               -md(i,k)*(cond(i,k)-chat(i,k))
              !     $              )

              ! version 2
              !                  fluxin =  mu(i,k)*chat(i,k) - md(i,k)*cond(i,k)
              !                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*chat(i,k)
              ! version 3
              fluxin  = mu(i,k) * min(chat(i,k), const(i,km1)) - md(i,k) * cond(i,k)
              fluxout = mu(i,k) * conu(i,k) - md(i,k) * min(chat(i,k), const(i,k))

              netflux = fluxin - fluxout
              if (abs(netflux) < max(fluxin, fluxout) * 1.0e-12_r8) then
                netflux = 0
              end if
              dcondt(i,k) = netflux/dptmp(i,k)
            else if (k > mx(i)) then
              dcondt(i,k) = 0
            end if
          end do
        end do

        if (zmconv_microp) then
          do i = il1g, il2g
            do k = jt(i), mx(i)
              if (dcondt(i,k) * dt + const(i,k) < 0) then
                negadt = dcondt(i,k) + const(i,k) / dt
                dcondt(i,k) = -const(i,k) / dt
                do kk = k + 1, mx(i)
                  if (negadt < 0 .and. dcondt(i,kk) * dt + const(i,kk) > 0) then
                    qtmp = dcondt(i,kk) + negadt * dptmp(i,k) / dptmp(i,kk)
                    if (qtmp*dt + const(i,kk) > 0) then
                      dcondt(i,kk)= qtmp
                      negadt = 0
                    else
                      negadt = negadt + (const(i,kk) / dt + dcondt(i,kk)) * dptmp(i,kk) / dptmp(i,k)
                      dcondt(i,kk) = -const(i,kk) / dt
                    end if
                  end if
                end do
                do kk = k - 1, jt(i), -1
                  if (negadt < 0 .and. dcondt(i,kk) * dt + const(i,kk) > 0) then
                    qtmp = dcondt(i,kk) + negadt * dptmp(i,k) / dptmp(i,kk)
                    if (qtmp * dt + const(i,kk) > 0) then
                      dcondt(i,kk) = qtmp
                      negadt = 0
                    else
                      negadt = negadt + (const(i,kk) / dt + dcondt(i,kk)) * dptmp(i,kk) / dptmp(i,k)
                      dcondt(i,kk) = -const(i,kk) / dt
                    end if
                  end if
                end do

                if (negadt < 0) then
                  dcondt(i,k) = dcondt(i,k) + negadt
                end if
              end if
            end do
          end do
        end if

        ! Initialize to zero everywhere, then scatter tendency back to full array
        dqdt(:,:,m) = 0
        do k = 1, pver
          kp1 = min(pver, k + 1)
          do i = il1g, il2g
            dqdt(ideep(i),k,m) = dcondt(i,k)
          end do
        end do
      end if ! for doconvtran
    end do

  end subroutine convtran

  subroutine momtran( &
    lchnk           , &
    ncol            , &
    domomtran       , &
    q               , &
    ncnst           , &
    mu              , &
    md              , &
    du              , &
    eu              , &
    ed              , &
    dp              , &
    dsubcld         , &
    jt              , &
    mx              , &
    ideep           , &
    il1g            , &
    il2g            , &
    nstep           , &
    dqdt            , &
    pguall          , &
    pgdall          , &
    icwu            , &
    icwd            , &
    dt              , &
    seten             &
  )
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Convective transport of momentum
    !
    ! Mixing ratios may be with respect to either dry or moist air
    !
    ! Method:
    ! Based on the convtran subroutine by P. Rasch
    ! <Also include any applicable external references.>
    !
    ! Author: J. Richter and P. Rasch
    !
    !-----------------------------------------------------------------------
    use constituents, only: cnst_get_type_byind

    integer , intent(in ) :: lchnk
    integer , intent(in ) :: ncol
    integer , intent(in ) :: ncnst
    logical , intent(in ), dimension(           ncnst) :: domomtran ! flag for doing convective transport
    real(r8), intent(in ), dimension(pcols,pver,ncnst) :: q         ! Wind array
    real(r8), intent(in ), dimension(pcols,pver      ) :: mu        ! Mass flux up
    real(r8), intent(in ), dimension(pcols,pver      ) :: md        ! Mass flux down
    real(r8), intent(in ), dimension(pcols,pver      ) :: du        ! Mass detraining from updraft
    real(r8), intent(in ), dimension(pcols,pver      ) :: eu        ! Mass entraining from updraft
    real(r8), intent(in ), dimension(pcols,pver      ) :: ed        ! Mass entraining from downdraft
    real(r8), intent(in ), dimension(pcols,pver      ) :: dp        ! Delta pressure between interfaces
    real(r8), intent(in ), dimension(pcols           ) :: dsubcld   ! Delta pressure from cloud base to sfc
    real(r8), intent(in ) :: dt                                     ! Time step in seconds
    integer , intent(in ), dimension(pcols           ) :: jt        ! Index of cloud top for each column
    integer , intent(in ), dimension(pcols           ) :: mx        ! Index of cloud top for each column
    integer , intent(in ), dimension(pcols           ) :: ideep     ! Gathering array
    integer , intent(in ) :: il1g                                   ! Gathered min lon indices over which to operate
    integer , intent(in ) :: il2g                                   ! Gathered max lon indices over which to operate
    integer , intent(in ) :: nstep                                  ! Time step index
    real(r8), intent(out), dimension(pcols,pver,ncnst) :: dqdt      ! Tracer tendency array

    integer i, k, kk, m
    integer kkp1, kkm1, km1, kp1
    integer ktm               ! Highest altitude index of cloud top
    integer kbm               ! Highest altitude index of cloud base

    real(r8) cabv                 ! Mix ratio of constituent above
    real(r8) cbel                 ! Mix ratio of constituent below
    real(r8) cdifr                ! Normalized diff between cabv and cbel
    real(r8) chat(pcols,pver)     ! Mix ratio in env at interfaces
    real(r8) cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
    real(r8) const(pcols,pver)    ! Gathered wind array
    real(r8) conu(pcols,pver)     ! Mix ratio in updraft at interfaces
    real(r8) dcondt(pcols,pver)   ! Gathered tend array
    real(r8) mbsth                ! Threshold for mass fluxes
    real(r8) mupdudp              ! A work variable
    real(r8) minc                 ! A work variable
    real(r8) maxc                 ! A work variable
    real(r8) fluxin               ! A work variable
    real(r8) fluxout              ! A work variable
    real(r8) netflux              ! A work variable

    real(r8) mududp(pcols,pver) ! working variable
    real(r8) mddudp(pcols,pver)     ! working variable

    real(r8) pgu(pcols,pver)      ! Pressure gradient term for updraft
    real(r8) pgd(pcols,pver)      ! Pressure gradient term for downdraft

    real(r8),intent(out) ::  pguall(pcols,pver,ncnst)      ! Apparent force from  updraft PG
    real(r8),intent(out) ::  pgdall(pcols,pver,ncnst)      ! Apparent force from  downdraft PG

    real(r8),intent(out) ::  icwu(pcols,pver,ncnst)      ! In-cloud winds in updraft
    real(r8),intent(out) ::  icwd(pcols,pver,ncnst)      ! In-cloud winds in downdraft

    real(r8),intent(out) ::  seten(pcols,pver) ! Dry static energy tendency
    real(r8)                 gseten(pcols,pver) ! Gathered dry static energy tendency

    real(r8) mflux(pcols,pverp,ncnst)   ! Gathered momentum flux

    real(r8) wind0(pcols,pver,ncnst)    ! Gathered wind before time step
    real(r8) windf(pcols,pver,ncnst)    ! Gathered wind after time step
    real(r8) fkeb, fket, ketend_cons, ketend, utop, ubot, vtop, vbot, gset2

    ! Initialize outgoing fields
    pguall     = 0
    pgdall     = 0
    ! Initialize in-cloud winds to environmental wind
    icwu(:ncol,:,:) = q(:ncol,:,:)
    icwd(:ncol,:,:) = q(:ncol,:,:)
    ! Initialize momentum flux and  final winds
    mflux      = 0
    wind0      = 0
    windf      = 0
    ! Initialize dry static energy
    seten      = 0
    gseten     = 0

    ! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
    mbsth = 1.e-15_r8

    ! Find the highest level top and bottom levels of convection
    ktm = pver
    kbm = pver
    do i = il1g, il2g
      ktm = min(ktm, jt(i))
      kbm = min(kbm, mx(i))
    end do

    ! Loop ever each wind component
    do m = 1, ncnst ! Start at m = 1 to transport momentum
      if (domomtran(m)) then
        ! Gather up the winds and set tend to zero
        do k = 1, pver
          do i = il1g, il2g
            const(i,k  ) = q(ideep(i),k,m)
            wind0(i,k,m) = const(i,k)
          end do
        end do

        ! From now on work only with gathered data
        ! Interpolate winds to interfaces
        do k = 1, pver
          km1 = max(1, k - 1)
          do i = il1g, il2g
            ! Use arithmetic mean
            chat(i,k) = 0.5_r8 * (const(i,k) + const(i,km1))
            ! Provisional up and down draft values
            conu(i,k) = chat(i,k)
            cond(i,k) = chat(i,k)
            ! Provisional tends
            dcondt(i,k) = 0
          end do
        end do
        !
        ! Pressure Perturbation Term
        !
        ! Top boundary: assume mu is zero
        k = 1
        pgu(:il2g,k) = 0
        pgd(:il2g,k) = 0
        do k = 2, pver - 1
          km1 = max(1   , k - 1)
          kp1 = min(pver, k + 1)
          do i = il1g, il2g
            ! Interior points
            mududp(i,k) = mu(i,k  ) * (const(i,k  ) - const(i,km1)) / dp(i,km1) &
                        + mu(i,kp1) * (const(i,kp1) - const(i,k  )) / dp(i,k  )
            mddudp(i,k) = md(i,k  ) * (const(i,k  ) - const(i,km1)) / dp(i,km1) &
                        + md(i,kp1) * (const(i,kp1) - const(i,k  )) / dp(i,k  )
            pgu(i,k) = - momcu * 0.5_r8 * mududp(i,k)
            pgd(i,k) = - momcd * 0.5_r8 * mddudp(i,k)
          end do
        end do
        ! Bottom boundary
        k = pver
        km1 = max(1, k - 1)
        do i=il1g, il2g
          mududp(i,k) =   mu(i,k) * (const(i,k) - const(i,km1)) / dp(i,km1)
          pgu   (i,k) = - momcu *  mududp(i,k)
          mddudp(i,k) =   md(i,k) * (const(i,k) - const(i,km1)) / dp(i,km1)
          pgd   (i,k) = - momcd * mddudp(i,k)
        end do
        !
        ! In-cloud velocity calculations
        !
        ! Do levels adjacent to top and bottom
        k = 2
        km1 = 1
        kk = pver
        kkm1 = max(1, kk - 1)
        do i = il1g, il2g
          mupdudp = mu(i,kk) + du(i,kk) * dp(i,kk)
          if (mupdudp > mbsth) then
            conu(i,kk) = ( eu(i,kk ) * const(i,kk ) * dp(i,kk ) + pgu(i,kk) * dp(i,kk)) / mupdudp
          end if
          if (md(i,k) < -mbsth) then
            cond(i,k ) = (-ed(i,km1) * const(i,km1) * dp(i,km1)) - pgd(i,km1) * dp(i,km1) / md(i,k)
          end if
        end do

        ! Updraft from bottom to top
        do kk = pver - 1, 1, -1
          kkm1 = max(1   , kk - 1)
          kkp1 = min(pver, kk + 1)
          do i = il1g, il2g
            mupdudp = mu(i,kk) + du(i,kk) * dp(i,kk)
            if (mupdudp > mbsth) then
              conu(i,kk) = (mu(i,kkp1) * conu(i,kkp1) + eu(i,kk) * &
                            const(i,kk) * dp(i,kk) + pgu(i,kk) * dp(i,kk)) / mupdudp
            end if
          end do
        end do

        ! Downdraft from top to bottom
        do k = 3, pver
          km1 = max(1, k - 1)
          do i = il1g, il2g
            if (md(i,k) < -mbsth) then
              cond(i,k) = (md(i,km1) * cond(i,km1) - ed(i,km1) * const(i,km1) &
                          * dp(i,km1)-pgd(i,km1) * dp(i,km1)) / md(i,k)
            end if
          end do
        end do

        do k = ktm,pver
          km1 = max(1   , k - 1)
          kp1 = min(pver, k + 1)
          do i = il1g, il2g
            ! Version 1 hard to check for roundoff errors
            dcondt(i,k) = &
              (mu(i,kp1) * (conu(i,kp1) - chat(i,kp1)) - &
               mu(i,k  ) * (conu(i,k  ) - chat(i,k  )) + &
               md(i,kp1) * (cond(i,kp1) - chat(i,kp1)) - &
               md(i,k  ) * (cond(i,k  ) - chat(i,k  )) &
              ) / dp(i,k)
          end do
        end do

        ! dcont for bottom layer
        do k = kbm, pver
          km1 = max(1, k - 1)
          do i = il1g, il2g
            if (k == mx(i)) then
              ! Version 1
              dcondt(i,k) = (1.0_r8 / dp(i,k)) *    &
                (-mu(i,k) * (conu(i,k) - chat(i,k)) &
                 -md(i,k) * (cond(i,k) - chat(i,k)))
            end if
          end do
        end do

        ! Initialize to zero everywhere, then scatter tendency back to full array
        dqdt(:,:,m) = 0
        do k = 1, pver
          do i = il1g, il2g
            dqdt  (ideep(i),k,m) = dcondt(i,k)
            ! Output apparent force on the mean flow from pressure gradient
            pguall(ideep(i),k,m) = -pgu(i,k)
            pgdall(ideep(i),k,m) = -pgd(i,k)
            icwu  (ideep(i),k,m) = conu(i,k)
            icwd  (ideep(i),k,m) = cond(i,k)
          end do
        end do
        ! Calculate momentum flux in units of mb*m/s2
        do k = ktm,pver
          do i = il1g, il2g
            mflux(i,k,m) = -mu(i,k) * (conu(i,k) - chat(i,k)) &
                           -md(i,k) * (cond(i,k) - chat(i,k))
          end do
        end do
        ! Calculate winds at the end of the time step
        do k = ktm,pver
          do i = il1g, il2g
            km1 = max(1, k - 1)
            kp1 = k + 1
            windf(i,k,m) = const(i,k)    -   (mflux(i,kp1,m) - mflux(i,k,m)) * dt /dp(i,k)
          end do
        end do
      end if
    end do

    ! Need to add an energy fix to account for the dissipation of kinetic energy
    ! Formulation follows from Boville and Bretherton (2003)
    ! formulation by PJR

    do k = ktm, pver
      km1 = max(1   , k - 1)
      kp1 = min(pver, k + 1)
      do i = il1g, il2g
        ! Calculate the KE fluxes at top and bottom of layer based on a discrete
        ! approximation to b&b eq(35) F_KE = u*F_u + v*F_v at interface
        utop = (wind0(i,k  ,1) + wind0(i,km1,1)) * 0.5_r8
        vtop = (wind0(i,k  ,2) + wind0(i,km1,2)) * 0.5_r8
        ubot = (wind0(i,kp1,1) + wind0(i,k  ,1)) * 0.5_r8
        vbot = (wind0(i,kp1,2) + wind0(i,k  ,2)) * 0.5_r8
        fket = utop*mflux(i,k  ,1) + vtop * mflux(i,k  ,2) ! top of layer
        fkeb = ubot*mflux(i,k+1,1) + vbot * mflux(i,k+1,2) ! bot of layer
        ! Divergence of these fluxes should give a conservative redistribution of KE
        ketend_cons = (fket-fkeb)/dp(i,k)
        ! Tendency in kinetic energy resulting from the momentum transport
        ketend = ((windf(i,k,1)**2 + windf(i,k,2)**2) - (wind0(i,k,1)**2 + wind0(i,k,2)**2)) * 0.5_r8 / dt
        ! The difference should be the dissipation
        gset2 = ketend_cons - ketend
        gseten(i,k) = gset2
      end do
    end do

    ! Scatter dry static energy to full array
    do k = 1, pver
      do i = il1g, il2g
        seten(ideep(i),k) = gseten(i,k)
      end do
    end do

  end subroutine momtran

  subroutine buoyan( &
    lchnk   , &
    ncol    , &
    q       , &
    t       , &
    p       , &
    z       , &
    pf      , &
    tp      , &
    qstp    , &
    tl      , &
    rl      , &
    cape    , &
    pblt    , &
    lcl     , &
    lel     , &
    lon     , &
    mx      , &
    rd      , &
    grav    , &
    cp      , &
    msg     , &
    tpert     &
  )

    integer , intent(in ) :: lchnk
    integer , intent(in ) :: ncol
    integer , intent(in ) :: msg
    real(r8), intent(in ) :: rd
    real(r8), intent(in ) :: grav
    real(r8), intent(in ) :: cp
    real(r8), intent(in ) :: rl
    real(r8), intent(in ), dimension(pcols,pver  ) :: q
    real(r8), intent(in ), dimension(pcols,pver  ) :: t
    real(r8), intent(in ), dimension(pcols,pver  ) :: p
    real(r8), intent(in ), dimension(pcols,pver  ) :: z
    real(r8), intent(in ), dimension(pcols,pver+1) :: pf
    real(r8), intent(in ), dimension(pcols       ) :: pblt
    real(r8), intent(in ), dimension(pcols       ) :: tpert   ! Perturbation temperature by pbl processes
    real(r8), intent(out), dimension(pcols,pver  ) :: tp      ! Parcel temperature
    real(r8), intent(out), dimension(pcols,pver  ) :: qstp    ! Saturation mixing ratio of parcel
    real(r8), intent(out), dimension(pcols       ) :: tl      ! Parcel temperature at LCL
    real(r8), intent(out), dimension(pcols       ) :: cape

    integer, dimension(pcols) :: lcl
    integer, dimension(pcols) :: lel
    integer, dimension(pcols) :: lon ! Level of onset of deep convection
    integer, dimension(pcols) :: mx  ! Level of max moist static energy

    logical , dimension(pcols        ) :: plge600
    integer , dimension(pcols        ) :: knt
    integer , dimension(pcols,num_cin) :: lelten
    real(r8), dimension(pcols,num_cin) :: capeten  ! provisional value of cape
    real(r8), dimension(pcols,pver   ) :: tv  
    real(r8), dimension(pcols,pver   ) :: tpv 
    real(r8), dimension(pcols,pver   ) :: buoy
    real(r8), dimension(pcols        ) :: a1
    real(r8), dimension(pcols        ) :: a2
    real(r8), dimension(pcols        ) :: estp
    real(r8), dimension(pcols        ) :: pl
    real(r8), dimension(pcols        ) :: plexp
    real(r8), dimension(pcols        ) :: hmax
    real(r8), dimension(pcols        ) :: hmn
    real(r8), dimension(pcols        ) :: y

    real(r8) e
    integer i, k, n

#ifdef PERGRO
    real(r8) rhd
#endif

    do n = 1, num_cin
      do i = 1, ncol
        lelten (i,n) = pver
        capeten(i,n) = 0
      end do
    end do

    do i = 1, ncol
      lon (i) = pver
      knt (i) = 0
      lel (i) = pver
      mx  (i) = lon(i)
      cape(i) = 0
      hmax(i) = 0
    end do

    tp  (:ncol,:) = t(:ncol,:)
    qstp(:ncol,:) = q(:ncol,:)

    ! Initialize tv and buoy for output.
    tv  (:ncol,:) = t (:ncol,:) * (1 + 1.608_r8 * q(:ncol,:)) / (1 + q(:ncol,:))
    tpv (:ncol,:) = tv(:ncol,:)
    buoy(:ncol,:) = 0

    !
    ! set "launching" level(mx) to be at maximum moist static energy.
    ! search for this level stops at planetary boundary layer top.
    !
#ifdef PERGRO
    do k = pver, msg + 1,-1
      do i = 1, ncol
        hmn(i) = cp * t(i,k) + grav * z(i,k) + rl * q(i,k)
        !
        ! Reset max moist static energy level when relative difference exceeds 1.e-4
        !
        rhd = (hmn(i) - hmax(i))/(hmn(i) + hmax(i))
        if (k >= nint(pblt(i)) .and. k <= lon(i) .and. rhd > -1.0e-4_r8) then
          hmax(i) = hmn(i)
          mx  (i) = k
        end if
      end do
    end do
#else
    do k = pver, msg + 1,-1
      do i = 1, ncol
        hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
        if (k >= nint(pblt(i)) .and. k <= lon(i) .and. hmn(i) > hmax(i)) then
          hmax(i) = hmn(i)
          mx  (i) = k
        end if
      end do
    end do
#endif

    do i = 1, ncol
      lcl(i) = mx(i)
      e = p(i,mx(i)) * q(i,mx(i)) / (eps1 + q(i,mx(i)))
      tl(i) = 2840.0_r8 / (3.5_r8 * log(t(i,mx(i))) - log(e) - 4.805_r8) + 55
      if (tl(i) < t(i,mx(i))) then
        plexp(i) = (1.0_r8 / (0.2854_r8 * (1 - 0.28_r8 * q(i,mx(i)))))
        pl   (i) = p(i,mx(i)) * (tl(i) / t(i,mx(i)))**plexp(i)
      else
        tl   (i) = t(i,mx(i))
        pl   (i) = p(i,mx(i))
      end if
    end do
    !
    ! Calculate lifting condensation level (LCL).
    !
    do k = pver, msg + 2, -1
      do i = 1, ncol
        if (k <= mx(i) .and. (p(i,k) > pl(i) .and. p(i,k-1) <= pl(i))) then
            lcl(i) = k - 1
        end if
      end do
    end do
    !
    ! If lcl is above the nominal level of non-divergence (600 mbs),
    ! no deep convection is permitted (ensuing calculations
    ! skipped and cape retains initialized value of zero).
    !
    do i = 1, ncol
      plge600(i) = pl(i) >= 600
    end do
    !
    ! initialize parcel properties in sub-cloud layer below lcl.
    !
    do k = pver, msg + 1,-1
      do i = 1, ncol
        if (k > lcl(i) .and. k <= mx(i) .and. plge600(i)) then
          tv  (i,k) = t(i,k) * (1 + 1.608_r8 * q(i,k)) / (1 + q(i,k))
          qstp(i,k) = q(i,mx(i))
          tp  (i,k) = t(i,mx(i)) * (p(i,k) / p(i,mx(i)))**(0.2854_r8 * (1 - 0.28_r8 * q(i,mx(i))))
          !
          ! buoyancy is increased by 0.5 k as in tiedtke
          !
          !-jjh          tpv (i,k)=tp(i,k)*(1.+1.608*q(i,mx(i)))/
          !-jjh     1                     (1.+q(i,mx(i)))
          tpv (i,k) = (tp(i,k) + tpert(i)) * (1 + 1.608_r8 * q(i,mx(i))) / (1 + q(i,mx(i)))
          buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add
        end if
      end do
    end do

    !
    ! Define parcel properties at LCL (i.e. level immediately above pl).
    !
    do k = pver, msg + 1,-1
      do i = 1, ncol
        if (k == lcl(i) .and. plge600(i)) then
          tv  (i,k) = t(i,k) * (1 + 1.608_r8 * q(i,k)) / (1 + q(i,k))
          qstp(i,k) = q(i,mx(i))
          tp  (i,k) = tl(i) * (p(i,k) / pl(i))**(0.2854_r8 * (1 - 0.28_r8 * qstp(i,k)))
          !              estp(i)  =exp(21.656_r8 - 5418._r8/tp(i,k))
          ! use of different formulas for es has about 1 g/kg difference
          ! in qs at t= 300k, and 0.02 g/kg at t=263k, with the formula
          ! above giving larger qs.
          call qsat_hPa(tp(i,k), p(i,k), estp(i), qstp(i,k))
          a1(i) = cp / rl + qstp(i,k) * (1 + qstp(i,k) / eps1) * rl * eps1 / (rd * tp(i,k)**2)
          a2(i) = 0.5_r8 * (qstp(i,k)* (1 + 2.0_r8 / eps1 * qstp(i,k)) * &
                  (1 + qstp(i,k) / eps1) * eps1**2 * rl * rl / &
                  (rd**2 * tp(i,k)**4) - qstp(i,k) * &
                  (1 + qstp(i,k) / eps1) * 2 * eps1 * rl / (rd * tp(i,k)**3))
          a1(i) = 1.0_r8 / a1(i)
          a2(i) = -a2(i) * a1(i)**3
          y (i) = q(i,mx(i)) - qstp(i,k)
          tp(i,k) = tp(i,k) + a1(i) * y(i) + a2(i) * y(i)**2
          call qsat_hPa(tp(i,k), p(i,k), estp(i), qstp(i,k))
          !
          ! buoyancy is increased by 0.5 k in cape calculation.
          ! dec. 9, 1994
          !-jjh          tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/(1.+q(i,mx(i)))
          !
          tpv (i,k) = (tp(i,k) + tpert(i)) * (1 + 1.608_r8 * qstp(i,k)) / (1 + q(i,mx(i)))
          buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add
        end if
      end do
    end do
    !
    ! Main buoyancy calculation.
    !
    do k = pver - 1, msg + 1, -1
      do i = 1, ncol
        if (k < lcl(i) .and. plge600(i)) then
          tv  (i,k) = t(i,k) * (1 + 1.608_r8 * q(i,k))/ (1 + q(i,k))
          qstp(i,k) = qstp(i,k+1)
          tp  (i,k) = tp(i,k+1) * (p(i,k) / p(i,k+1))**(0.2854_r8 * (1 - 0.28_r8 * qstp(i,k)))
          call qsat_hPa(tp(i,k), p(i,k), estp(i), qstp(i,k))
          a1(i) = cp / rl + qstp(i,k) * (1 + qstp(i,k) / eps1) * rl * eps1 / (rd * tp(i,k)**2)
          a2(i) = 0.5_r8 * (qstp(i,k) * (1 + 2.0_r8 / eps1 * qstp(i,k)) * &
                (1 + qstp(i,k) / eps1) * eps1**2 * rl * rl / &
                (rd**2 * tp(i,k)**4) - qstp(i,k) * &
                (1 + qstp(i,k) / eps1) * 2 * eps1 * rl / (rd * tp(i,k)**3))
          a1(i) = 1.0_r8 / a1(i)
          a2(i) = -a2(i) * a1(i)**3
          y (i) = qstp(i,k+1) - qstp(i,k)
          tp(i,k) = tp(i,k) + a1(i) * y(i) + a2(i) * y(i)**2
          call qsat_hPa(tp(i,k), p(i,k), estp(i), qstp(i,k))
          tpv (i,k) = (tp(i,k) + tpert(i)) * (1 + 1.608_r8 * qstp(i,k)) / (1 + q(i,mx(i)))
          buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add
        end if
      end do
    end do

    do k = msg + 2, pver
      do i = 1, ncol
        if (k < lcl(i) .and. plge600(i)) then
          if (buoy(i,k+1) > 0.0_r8 .and. buoy(i,k) <= 0) then
            knt(i) = min(5,knt(i) + 1)
            lelten(i,knt(i)) = k
          end if
        end if
      end do
    end do
    !
    ! Calculate convective available potential energy (CAPE).
    !
    do n = 1, 5
      do k = msg + 1, pver
        do i = 1, ncol
          if (plge600(i) .and. k <= mx(i) .and. k > lelten(i,n)) then
            capeten(i,n) = capeten(i,n) + rd * buoy(i,k) * log(pf(i,k+1) / pf(i,k))
          end if
        end do
      end do
    end do
    !
    ! Find maximum cape from all possible tentative capes from
    ! one sounding, and use it as the final cape, april 26, 1995
    !
    do n = 1,5
      do i = 1, ncol
        if (capeten(i,n) > cape(i)) then
          cape(i) = capeten(i,n)
          lel (i) = lelten (i,n)
        end if
      end do
    end do
    !
    ! Put lower bound on cape for diagnostic purposes.
    !
    do i = 1, ncol
      cape(i) = max(cape(i), 0.0_r8)
    end do

  end subroutine buoyan

  subroutine cldprp( &
    lchnk          , &
    q              , &
    t              , &
    u              , &
    v              , &
    p              , &
    z              , &
    s              , &
    mu             , &
    eu             , &
    du             , &
    md             , &
    ed             , &
    sd             , &
    qd             , &
    mc             , &
    qu             , &
    su             , &
    zf             , &
    qst            , &
    hmn            , &
    hsat           , &
    shat           , &
    ql             , &
    cmeg           , &
    jb             , &
    lel            , &
    jt             , &
    jlcl           , &
    mx             , &
    j0             , &
    jd             , &
    rl             , &
    il2g           , &
    rd             , &
    grav           , &
    cp             , &
    msg            , &
    pflx           , &
    evp            , &
    cu             , &
    rprd           , &
    limcnv         , &
    landfrac       , &
    qcde           , &
    aero           , &
    loc_conv,qhat  )

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! <Say what the routine does>
    !
    ! Method:
    ! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
    !             original version cldprop.
    !
    ! Author: See above, modified by P. Rasch
    ! This is contributed code not fully standardized by the CCM core group.
    !
    ! this code is very much rougher than virtually anything else in the CCM
    ! there are debug statements left strewn about and code segments disabled
    ! these are to facilitate future development. We expect to release a
    ! cleaner code in a future release
    !
    ! the documentation has been enhanced to the degree that we are able
    !
    !-----------------------------------------------------------------------

    integer , intent(in ) :: lchnk
    integer , intent(in ) :: limcnv                                 ! Convection limiting level
    integer , intent(in ) :: il2g
    integer , intent(in ) :: msg                                    ! missing moisture vals (always 0)
    real(r8), intent(in ) :: rd
    real(r8), intent(in ) :: grav
    real(r8), intent(in ) :: cp
    real(r8), intent(in ) :: rl                                     ! Latent heat of vap
    real(r8), intent(in ), dimension(pcols,pver ) :: q
    real(r8), intent(in ), dimension(pcols,pver ) :: t
    real(r8), intent(in ), dimension(pcols,pver ) :: p
    real(r8), intent(in ), dimension(pcols,pver ) :: z
    real(r8), intent(in ), dimension(pcols,pver ) :: s
    real(r8), intent(in ), dimension(pcols,pverp) :: zf
    real(r8), intent(in ), dimension(pcols,pver ) :: u
    real(r8), intent(in ), dimension(pcols,pver ) :: v
    real(r8), intent(in ), dimension(pcols      ) :: landfrac
    integer , intent(in ), dimension(pcols      ) :: jb             ! updraft base level
    integer , intent(in ), dimension(pcols      ) :: lel            ! updraft launch level
    integer , intent(out), dimension(pcols      ) :: jt             ! updraft plume top
    integer , intent(out), dimension(pcols      ) :: jlcl           ! updraft lifting cond level
    integer , intent(in ), dimension(pcols      ) :: mx             ! updraft base level (same is jb)
    integer , intent(out), dimension(pcols      ) :: j0             ! level where updraft begins detraining
    integer , intent(out), dimension(pcols      ) :: jd             ! level of downdraft
    real(r8), intent(in ), dimension(pcols,pver ) :: shat           ! Interface values of dry stat energy
    real(r8), intent(in ), dimension(pcols,pver ) :: qhat           ! Grid slice of upper interface mixing ratio
    real(r8), intent(out), dimension(pcols,pver ) :: rprd           ! Production rate of precipitation at that layer
    real(r8), intent(out), dimension(pcols,pver ) :: du             ! detrainement rate of updraft
    real(r8), intent(out), dimension(pcols,pver ) :: ed             ! Entrainment rate of downdraft
    real(r8), intent(out), dimension(pcols,pver ) :: eu             ! Entrainment rate of updraft
    real(r8), intent(out), dimension(pcols,pver ) :: hmn            ! Moist static energy of environment
    real(r8), intent(out), dimension(pcols,pver ) :: hsat           ! Saturated moist static energy of environment
    real(r8), intent(out), dimension(pcols,pver ) :: mc             ! Net mass flux
    real(r8), intent(out), dimension(pcols,pver ) :: md             ! Downdraft mass flux
    real(r8), intent(out), dimension(pcols,pver ) :: mu             ! Updraft mass flux
    real(r8), intent(out), dimension(pcols,pverp) :: pflx           ! Precipitation flux
    real(r8), intent(out), dimension(pcols,pver ) :: qd             ! Specific humidity of downdraft
    real(r8), intent(out), dimension(pcols,pver ) :: ql             ! Liquid water of updraft
    real(r8), intent(out), dimension(pcols,pver ) :: qst            ! Saturation mixing ratio of env.
    real(r8), intent(out), dimension(pcols,pver ) :: qu             ! Specific humidity of updraft
    real(r8), intent(out), dimension(pcols,pver ) :: sd             ! Normalized dry stat energy of downdraft
    real(r8), intent(out), dimension(pcols,pver ) :: su             ! Normalized dry stat energy of updraft
    real(r8), intent(out), dimension(pcols,pver ) :: qcde           ! Cloud water mixing ratio for detrainment (kg/kg)
    real(r8), intent(out), dimension(pcols,pver ) :: cmeg
    real(r8), intent(out), dimension(pcols,pver ) :: evp
    real(r8), intent(out), dimension(pcols,pver ) :: cu
    type(zm_aero_t), intent(in) :: aero

    type(zm_conv_t) loc_conv

    real(r8), dimension(pcols,pver) :: gamma
    real(r8), dimension(pcols,pver) :: dz
    real(r8), dimension(pcols,pver) :: iprm
    real(r8), dimension(pcols,pver) :: hu
    real(r8), dimension(pcols,pver) :: hd
    real(r8), dimension(pcols,pver) :: eps
    real(r8), dimension(pcols,pver) :: f
    real(r8), dimension(pcols,pver) :: k1
    real(r8), dimension(pcols,pver) :: i2
    real(r8), dimension(pcols,pver) :: ihat
    real(r8), dimension(pcols,pver) :: i3
    real(r8), dimension(pcols,pver) :: idag
    real(r8), dimension(pcols,pver) :: i4
    real(r8), dimension(pcols,pver) :: qsthat
    real(r8), dimension(pcols,pver) :: hsthat
    real(r8), dimension(pcols,pver) :: gamhat
    real(r8), dimension(pcols,pver) :: qds
    real(r8), dimension(pcols     ) :: c0mask
    real(r8), dimension(pcols     ) :: hmin
    real(r8), dimension(pcols     ) :: expdif
    real(r8), dimension(pcols     ) :: expnum
    real(r8), dimension(pcols     ) :: ftemp
    real(r8), dimension(pcols     ) :: eps0
    real(r8), dimension(pcols     ) :: rmue
    real(r8), dimension(pcols     ) :: zuef
    real(r8), dimension(pcols     ) :: zdef
    real(r8), dimension(pcols     ) :: epsm
    real(r8), dimension(pcols     ) :: ratmjb
    real(r8), dimension(pcols     ) :: est
    real(r8), dimension(pcols     ) :: totpcp
    real(r8), dimension(pcols     ) :: totevp
    real(r8), dimension(pcols     ) :: alfa
    real(r8), dimension(pcols,pver) :: fice
    real(r8), dimension(pcols,pver) :: tug
    real(r8), dimension(pcols,pver) :: tvuo
    real(r8), dimension(pcols,pver) :: tvu
    real(r8), dimension(pcols     ) :: totfrz
    real(r8), dimension(pcols,pver) :: frz
    integer , dimension(pcols     ) :: jto
    integer , dimension(pcols     ) :: tmplel

    real(r8) ql1
    real(r8) tu
    real(r8) estu
    real(r8) qstu
    real(r8) small
    real(r8) mdt
    integer iter, itnum, i, k, m
    integer khighest
    integer klowest
    integer kount
    logical doit(pcols)
    logical done(pcols)

    if (zmconv_microp) then
      loc_conv%autolm(:il2g,:) = 0
      loc_conv%accrlm(:il2g,:) = 0
      loc_conv%bergnm(:il2g,:) = 0
      loc_conv%fhtimm(:il2g,:) = 0
      loc_conv%fhtctm(:il2g,:) = 0
      loc_conv%fhmlm (:il2g,:) = 0
      loc_conv%hmpim (:il2g,:) = 0
      loc_conv%accslm(:il2g,:) = 0
      loc_conv%dlfm  (:il2g,:) = 0

      loc_conv%autoln(:il2g,:) = 0
      loc_conv%accrln(:il2g,:) = 0
      loc_conv%bergnn(:il2g,:) = 0
      loc_conv%fhtimn(:il2g,:) = 0
      loc_conv%fhtctn(:il2g,:) = 0
      loc_conv%fhmln (:il2g,:) = 0
      loc_conv%accsln(:il2g,:) = 0
      loc_conv%activn(:il2g,:) = 0
      loc_conv%dlfn  (:il2g,:) = 0

      loc_conv%autoim(:il2g,:) = 0
      loc_conv%accsim(:il2g,:) = 0
      loc_conv%difm  (:il2g,:) = 0

      loc_conv%nuclin(:il2g,:) = 0
      loc_conv%autoin(:il2g,:) = 0
      loc_conv%accsin(:il2g,:) = 0
      loc_conv%hmpin (:il2g,:) = 0
      loc_conv%difn  (:il2g,:) = 0

      loc_conv%trspcm(:il2g,:) = 0
      loc_conv%trspcn(:il2g,:) = 0
      loc_conv%trspim(:il2g,:) = 0
      loc_conv%trspin(:il2g,:) = 0

      loc_conv%dcape (:il2g)   = 0
    end if
    do i = 1, il2g
      ftemp (i) = 0
      expnum(i) = 0
      expdif(i) = 0
      c0mask(i) = c0_ocn * (1 - landfrac(i)) + c0_lnd * landfrac(i)
    end do
    !
    !jr Change from msg+1 to 1 to prevent blowup
    !
    do k = 1, pver
      do i = 1, il2g
        dz(i,k) = zf(i,k) - zf(i,k+1)
      end do
    end do

    !
    ! Initialize many output and work variables to zero
    !
    pflx(:il2g,1) = 0
    do k = 1, pver
      do i = 1, il2g
        k1  (i,k) = 0
        i2  (i,k) = 0
        i3  (i,k) = 0
        i4  (i,k) = 0
        mu  (i,k) = 0
        f   (i,k) = 0
        eps (i,k) = 0
        eu  (i,k) = 0
        du  (i,k) = 0
        ql  (i,k) = 0
        cu  (i,k) = 0
        evp (i,k) = 0
        cmeg(i,k) = 0
        qds (i,k) = q(i,k)
        md  (i,k) = 0
        ed  (i,k) = 0
        sd  (i,k) = s(i,k)
        qd  (i,k) = q(i,k)
        mc  (i,k) = 0
        qu  (i,k) = q(i,k)
        su  (i,k) = s(i,k)
        call qsat_hPa(t(i,k), p(i,k), est(i), qst(i,k))
        if (p(i,k) - est(i) <= 0) then
          qst(i,k) = 1
        end if

        gamma(i,k) = qst(i,k) * (1 + qst(i,k) / eps1) * eps1*rl/(rd*t(i,k)**2) * rl / cp
        hmn  (i,k) = cp * t(i,k) + grav * z(i,k) + rl * q(i,k)
        hsat (i,k) = cp * t(i,k) + grav * z(i,k) + rl * qst(i,k)
        hu   (i,k) = hmn(i,k)
        hd   (i,k) = hmn(i,k)
        rprd (i,k) = 0
        fice (i,k) = 0
        tug  (i,k) = 0
        qcde (i,k) = 0
        tvuo (i,k) = (shat(i,k) - grav / cp * zf(i,k)) * (1 + 0.608_r8 * qhat(i,k))
        tvu  (i,k) = tvuo(i,k)
        frz  (i,k) = 0
      end do
    end do

    if (zmconv_microp) then
      do k = 1, pver
        do i = 1, il2g
          loc_conv%sprd (i,k) = 0
          loc_conv%wu   (i,k) = 0
          loc_conv%cmel (i,k) = 0
          loc_conv%cmei (i,k) = 0
          loc_conv%qliq (i,k) = 0
          loc_conv%qice (i,k) = 0
          loc_conv%qnl  (i,k) = 0
          loc_conv%qni  (i,k) = 0
          loc_conv%qide (i,k) = 0
          loc_conv%qncde(i,k) = 0
          loc_conv%qnide(i,k) = 0
          loc_conv%qnr  (i,k) = 0
          loc_conv%qns  (i,k) = 0
          loc_conv%qrain(i,k) = 0
          loc_conv%qsnow(i,k) = 0
          loc_conv%frz  (i,k) = 0
        end do
      end do
    end if
    !
    !jr Set to zero things which make this routine blow up
    !
    do k = 1, msg
      do i = 1, il2g
        rprd(i,k) = 0
      end do
    end do
    !
    ! Interpolate the layer values of qst, hsat and gamma to
    ! layer interfaces
    !
    do k = 1, msg + 1
      do i = 1, il2g
        hsthat(i,k) = hsat(i,k)
        qsthat(i,k) = qst(i,k)
        gamhat(i,k) = gamma(i,k)
      end do
    end do
    do i = 1, il2g
      totpcp(i) = 0
      totevp(i) = 0
    end do
    do k = msg + 2, pver
      do i = 1, il2g
        if (abs(qst(i,k-1) - qst(i,k)) > 1.0e-6_r8) then
          qsthat(i,k) = log(qst(i,k-1)/qst(i,k))*qst(i,k-1)*qst(i,k)/ (qst(i,k-1)-qst(i,k))
        else
          qsthat(i,k) = qst(i,k)
        end if
        hsthat(i,k) = cp*shat(i,k) + rl*qsthat(i,k)
        if (abs(gamma(i,k-1)-gamma(i,k)) > 1.E-6_r8) then
          gamhat(i,k) = log(gamma(i,k-1) / gamma(i,k)) * gamma(i,k-1) * gamma(i,k) / (gamma(i,k-1)-gamma(i,k))
        else
          gamhat(i,k) = gamma(i,k)
        end if
      end do
    end do
    !
    ! Initialize cloud top to highest plume top.
    !jr changed hard-wired 4 to limcnv+1 (not to exceed pver)
    !
    jt(:) = pver
    do i = 1, il2g
      jt(i) = max(lel(i),limcnv+1)
      jt(i) = min(jt(i),pver)
      jd(i) = pver
      jlcl(i) = lel(i)
      hmin(i) = 1.0e6_r8
    end do
    !
    ! Find the level of minimum hsat, where detrainment starts
    !
    do k = msg + 1, pver
      do i = 1, il2g
        if (hsat(i,k) <= hmin(i) .and. k >= jt(i) .and. k <= jb(i)) then
          hmin(i) = hsat(i,k)
          j0(i) = k
        end if
      end do
    end do
    do i = 1, il2g
      j0(i) = min(j0(i),jb(i)-2)
      j0(i) = max(j0(i),jt(i)+2)
      !
      ! Fix from Guang Zhang to address out of bounds array reference
      !
      j0(i) = min(j0(i), pver)
    end do
    !
    ! Initialize certain arrays inside cloud
    !
    do k = msg + 1, pver
      do i = 1, il2g
        if (k >= jt(i) .and. k <= jb(i)) then
          hu(i,k) = hmn(i,mx(i)) + cp * tiedke_add
          su(i,k) = s(i,mx(i)) + tiedke_add
        end if
      end do
    end do
    !
    ! *********************************************************
    ! Compute taylor series for approximate eps(z) below
    ! *********************************************************
    !
    do k = pver - 1,msg + 1,-1
      do i = 1, il2g
        if (k < jb(i) .and. k >= jt(i)) then
          k1  (i,k) = k1(i,k+1) + (hmn(i,mx(i)) - hmn(i,k)) * dz(i,k)
          ihat(i,k) = 0.5_r8 * (k1(i,k+1) + k1(i,k))
          i2  (i,k) = i2(i,k+1) + ihat(i,k) * dz(i,k)
          idag(i,k) = 0.5_r8 * (i2(i,k+1) + i2(i,k))
          i3  (i,k) = i3(i,k+1) + idag(i,k) * dz(i,k)
          iprm(i,k) = 0.5_r8 * (i3(i,k+1) + i3(i,k))
          i4  (i,k) = i4(i,k+1) + iprm(i,k) * dz(i,k)
        end if
      end do
    end do
    !
    ! Re-initialize hmin array for ensuing calculation.
    !
    do i = 1, il2g
      hmin(i) = 1.0e6_r8
    end do
    do k = msg + 1, pver
      do i = 1, il2g
        if (k >= j0(i) .and. k <= jb(i) .and. hmn(i,k) <= hmin(i)) then
          hmin  (i) = hmn(i,k)
          expdif(i) = hmn(i,mx(i)) - hmin(i)
        end if
      end do
    end do
    !
    ! *********************************************************
    ! Compute approximate eps(z) using above taylor series
    ! *********************************************************
    !
    do k = msg + 2, pver
      do i = 1, il2g
        expnum(i) = 0
        ftemp (i) = 0
        if (k < jt(i) .or. k >= jb(i)) then
          k1(i,k) = 0
          expnum(i) = 0
        else
          expnum(i) = hmn(i,mx(i)) - (hsat(i,k-1) * (zf(i,k) - z(i,k)) + &
                      hsat(i,k) * (z(i,k-1) - zf(i,k))) / (z(i,k-1) - z(i,k))
        end if
        if ((expdif(i) > 100 .and. expnum(i) > 0) .and. k1(i,k) > expnum(i) * dz(i,k)) then
          ftemp(i) = expnum(i) / k1(i,k)
          f(i,k) = ftemp(i) + i2(i,k) / k1(i,k) * ftemp(i)**2 + &
                  (2 * i2(i,k)**2 - k1(i,k) * i3(i,k)) / k1(i,k)**2 * &
                  ftemp(i)**3 + (-5 * k1(i,k) * i2(i,k) * i3(i,k) + &
                  5 * i2(i,k)**3 + k1(i,k)**2 * i4(i,k)) / &
                  k1(i,k)**3 * ftemp(i)**4
          f(i,k) = max(f(i,k), 0.0_r8)
          f(i,k) = min(f(i,k), 0.0002_r8)
        end if
      end do
    end do
    do i = 1, il2g
      if (j0(i) < jb(i)) then
        if (f(i,j0(i)) < 1.0e-6_r8 .and. f(i,j0(i)+1) > f(i,j0(i))) j0(i) = j0(i) + 1
      end if
    end do
    do k = msg + 2,pver
      do i = 1, il2g
        if (k >= jt(i) .and. k <= j0(i)) then
          f(i,k) = max(f(i,k), f(i,k-1))
        end if
      end do
    end do
    do i = 1, il2g
      eps0(i) = f(i,j0(i))
      eps(i,jb(i)) = eps0(i)
    end do
    !
    ! This is set to match the Rasch and Kristjansson paper
    !
    do k = pver, msg + 1,-1
      do i = 1, il2g
        if (k >= j0(i) .and. k <= jb(i)) then
          eps(i,k) = f(i,j0(i))
        end if
      end do
    end do
    do k = pver, msg + 1,-1
      do i = 1, il2g
        if (k < j0(i) .and. k >= jt(i)) eps(i,k) = f(i,k)
      end do
    end do
    if (zmconv_microp) then
      itnum = 2
    else
      itnum = 1
    end if

    do iter = 1, itnum
      if (zmconv_microp) then
        do k = pver, msg + 1,-1
          do i = 1, il2g
            cu(i,k) = 0
            ql(i,k) = 0
            loc_conv%qliq(i,k) = 0
            loc_conv%qice(i,k) = 0
            loc_conv%frz (i,k) = 0
          end do
        end do
        do i = 1, il2g
          totpcp(i) = 0
          hu(i,jb(i)) = hmn(i,jb(i)) + cp * tiedke_add
        end do
      end if
      !
      ! Specify the updraft mass flux mu, entrainment eu, detrainment du
      ! and moist static energy hu.
      ! here and below mu, eu,du, md and ed are all normalized by mb
      !
      do i = 1, il2g
        if (eps0(i) > 0) then
          mu(i,jb(i)) = 1
          eu(i,jb(i)) = mu(i,jb(i)) / dz(i,jb(i))
        end if
        if (zmconv_microp) then
          tmplel(i) = lel(i)
        else
          tmplel(i) = jt(i)
        end if
      end do
      do k = pver, msg + 1,-1
        do i = 1, il2g
          if (eps0(i) > 0.0_r8 .and. (k >= tmplel(i) .and. k < jb(i))) then
            zuef(i) = zf(i,k) - zf(i,jb(i))
            rmue(i) = (1.0_r8 / eps0(i))* (exp(eps(i,k+1) * zuef(i)) - 1) / zuef(i)
            mu(i,k) = (1.0_r8 / eps0(i))* (exp(eps(i,k  ) * zuef(i)) - 1) / zuef(i)
            eu(i,k) = (rmue(i) - mu(i,k+1)) / dz(i,k)
            du(i,k) = (rmue(i) - mu(i,k  )) / dz(i,k)
          end if
        end do
      end do

      khighest = pverp
      klowest = 1
      do i = 1, il2g
        khighest = min(khighest, lel(i))
        klowest  = max(klowest , jb (i))
      end do
      do k = klowest - 1, khighest, -1
        do i = 1, il2g
          if (k <= jb(i)-1 .and. k >= lel(i) .and. eps0(i) > 0) then
            if (mu(i,k) < 0.02_r8) then
              hu(i,k) = hmn(i,k)
              mu(i,k) = 0
              eu(i,k) = 0
              du(i,k) = mu(i,k+1) / dz(i,k)
            else
              if (zmconv_microp) then
                hu(i,k) = (mu(i,k+1) * hu(i,k+1) + dz(i,k) * (eu(i,k) * hmn(i,k) +   &
                          latice * frz(i,k))) / (mu(i,k) + dz(i,k) * du(i,k))
              else
                hu(i,k) = mu(i,k+1) / mu(i,k) * hu(i,k+1) + &
                          dz(i,k) / mu(i,k) * (eu(i,k) * hmn(i,k) - du(i,k) * hsat(i,k))
              end if
            end if
          end if
        end do
      end do
      !
      ! Reset cloud top index beginning from two layers above the
      ! cloud base (i.e. if cloud is only one layer thick, top is not reset
      !
      do i = 1, il2g
        doit  (i) = .true.
        totfrz(i)= 0
        do k = pver, msg + 1, -1
          totfrz(i) = totfrz(i) + frz(i,k) * dz(i,k)
        end do
      end do
      do k = klowest - 2, khighest - 1, -1
        do i = 1, il2g
          if (doit(i) .and. k <= jb(i)-2 .and. k >= lel(i)-1) then
            if (hu(i,k) <= hsthat(i,k) .and. hu(i,k+1) > hsthat(i,k+1) .and. mu(i,k) >= 0.02_r8) then
              if (hu(i,k) - hsthat(i,k) < -2000) then
                jt(i) = k + 1
                doit(i) = .false.
              else
                jt(i) = k
                doit(i) = .false.
              end if
            else if ( (hu(i,k) > hu(i,jb(i)) .and. totfrz(i)<=0.0_r8) .or. mu(i,k) < 0.02_r8) then
              jt(i) = k + 1
              doit(i) = .false.
            end if
          end if
        end do
      end do

      if (iter == 1)  jto(:) = jt(:)

      do k = pver, msg + 1,-1
        do i = 1, il2g
          if (k >= lel(i) .and. k <= jt(i) .and. eps0(i) > 0) then
            mu(i,k) = 0
            eu(i,k) = 0
            du(i,k) = 0
            hu(i,k) = hmn(i,k)
          end if
          if (k == jt(i) .and. eps0(i) > 0) then
            du(i,k) = mu(i,k+1) / dz(i,k)
            eu(i,k) = 0
            mu(i,k) = 0
          end if
        end do
      end do

      do i = 1, il2g
        done(i) = .false.
      end do
      kount = 0
      do k = pver, msg + 2, -1
        do i = 1, il2g
          if (k == jb(i) .and. eps0(i) > 0) then
            qu(i,k) = q(i,mx(i))
            su(i,k) = (hu(i,k) - rl * qu(i,k)) / cp
          end if
          if ((.not. done(i) .and. k > jt(i) .and. k < jb(i)) .and. eps0(i) > 0) then
            su(i,k) = mu(i,k+1) / mu(i,k) * su(i,k+1) + &
                      dz(i,k  ) / mu(i,k) * (eu(i,k) - du(i,k)) * s(i,k)
            qu(i,k) = mu(i,k+1) / mu(i,k) * qu(i,k+1) + &
                      dz(i,k  ) / mu(i,k) * (eu(i,k) * q(i,k) - du(i,k) * qst(i,k))
            tu = su(i,k) - grav / cp * zf(i,k)
            call qsat_hPa(tu, (p(i,k) + p(i,k-1)) * 0.5_r8, estu, qstu)
            if (qu(i,k) >= qstu) then
              jlcl(i) = k
              kount = kount + 1
              done(i) = .true.
            end if
          end if
        end do
        if (kount >= il2g) exit
      end do
      do k = msg + 2,pver
        do i = 1, il2g
          if ((k > jt(i) .and. k <= jlcl(i)) .and. eps0(i) > 0) then
            su(i,k) = shat(i,k) + (hu(i,k) - hsthat(i,k)) / (cp * (1 + gamhat(i,k)))
            qu(i,k) = qsthat(i,k) + gamhat(i,k) * (hu(i,k) - hsthat(i,k)) / &
                      (rl * (1 + gamhat(i,k)))
          end if
        end do
      end do

      ! Compute condensation in updraft
      if (zmconv_microp) then
        tmplel(:il2g) = jlcl(:il2g) + 1
      else
        tmplel(:il2g) = jb(:il2g)
      end if

      do k = pver, msg + 2, -1
        do i = 1, il2g
          if (k >= jt(i) .and. k < tmplel(i) .and. eps0(i) > 0) then
            if (zmconv_microp) then
              cu(i,k) = ((mu(i,k) * su(i,k) - mu(i,k+1) * su(i,k+1)) / &
                      dz(i,k) - eu(i,k) * s(i,k) + du(i,k) * su(i,k)) / (rl / cp)  &
                      - latice * frz(i,k) / rl
            else
              cu(i,k) = ((mu(i,k) * su(i,k) - mu(i,k+1) * su(i,k+1)) / &
                      dz(i,k) - (eu(i,k) - du(i,k)) * s(i,k)) / (rl / cp)
            end if
            if (k == jt(i)) cu(i,k) = 0
            cu(i,k) = max(0.0_r8, cu(i,k))
          end if
        end do
      end do

      if (zmconv_microp) then
        tug(:il2g,:) = t(:il2g,:)
        fice(:,:)    = 0
        do k = pver, msg + 2, -1
          do i = 1, il2g
            tug(i,k) = su(i,k) - grav / cp * zf(i,k)
          end do
        end do
        do k = 1, pver - 1
          do i = 1, il2g
            if (tug(i,k+1) > 273.15_r8) then
              ! If warmer than tmax then water phase
              fice(i,k) = 0
            else if (tug(i,k+1) < 233.15_r8) then
              ! If colder than tmin then ice phase
              fice(i,k) = 1
            else
              ! Otherwise mixed phase, with ice fraction decreasing linearly
              ! from tmin to tmax
              fice(i,k) = (273.15_r8 - tug(i,k+1)) / 40.0_r8
            end if
          end do
        end do
        do k = 1, pver
          do i = 1, il2g
            loc_conv%cmei(i,k) = cu(i,k) * fice(i,k)
            loc_conv%cmel(i,k) = cu(i,k) * (1 - fice(i,k))
          end do
        end do
        call  zm_mphy(         &
          su                 , &
          qu                 , &
          mu                 , &
          du                 , &
          eu                 , &
          loc_conv%cmel      , &
          loc_conv%cmei      , &
          zf                 , &
          p                  , &
          t                  , &
          q                  , &
          eps0               , &
          jb                 , &
          jt                 , &
          jlcl               , &
          msg                , &
          il2g               , &
          grav               , &
          cp                 , &
          rd                 , &
          aero               , &
          gamhat             , &
          loc_conv%qliq      , &
          loc_conv%qice      , &
          loc_conv%qnl       , &
          loc_conv%qni       , &
          qcde               , &
          loc_conv%qide      , &
          loc_conv%qncde     , &
          loc_conv%qnide     , &
          rprd               , &
          loc_conv%sprd      , &
          frz                , &
          loc_conv%wu        , &
          loc_conv%qrain     , &
          loc_conv%qsnow     , &
          loc_conv%qnr       , &
          loc_conv%qns       , &
          loc_conv%autolm    , &
          loc_conv%accrlm    , &
          loc_conv%bergnm    , &
          loc_conv%fhtimm    , &
          loc_conv%fhtctm    , &
          loc_conv%fhmlm     , &
          loc_conv%hmpim     , &
          loc_conv%accslm    , &
          loc_conv%dlfm      , &
          loc_conv%autoln    , &
          loc_conv%accrln    , &
          loc_conv%bergnn    , &
          loc_conv%fhtimn    , &
          loc_conv%fhtctn    , &
          loc_conv%fhmln     , &
          loc_conv%accsln    , &
          loc_conv%activn    , &
          loc_conv%dlfn      , &
          loc_conv%autoim    , &
          loc_conv%accsim    , &
          loc_conv%difm      , &
          loc_conv%nuclin    , &
          loc_conv%autoin    , &
          loc_conv%accsin    , &
          loc_conv%hmpin     , &
          loc_conv%difn      , &
          loc_conv%trspcm    , &
          loc_conv%trspcn    , &
          loc_conv%trspim    , &
          loc_conv%trspin    , &
          loc_conv%lambdadpcu, &
          loc_conv%mudpcu      &
        )
        do k = pver, msg + 2, -1
          do i = 1, il2g
            ql(i,k) = loc_conv%qliq(i,k) + loc_conv%qice(i,k)
            loc_conv%frz(i,k) = frz(i,k)
          end do
        end do
        do i = 1, il2g
          if (iter == 2 .and. jt(i)> jto(i)) then
            do k = jt(i), jto(i), -1
              loc_conv%frz(i,k) = 0
              cu(i,k) = 0
            end do
          end if
        end do
        do k = pver, msg + 2, -1
          do i = 1, il2g
            if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0 .and. mu(i,k) >= 0) then
              totpcp(i) = totpcp(i) + dz(i,k) * (cu(i,k) - du(i,k) * (qcde(i,k+1) + loc_conv%qide(i,k+1)))
            end if
          end do
        end do
        do k = msg + 2, pver
          do i = 1, il2g
            if ((k > jt(i) .and. k <= jlcl(i)) .and. eps0(i) > 0) then
              if (iter == 1) tvuo(i,k) = (su(i,k) - grav / cp * zf(i,k)) * (1 + 0.608_r8 * qu(i,k))
              if (iter == 2 .and. k > max(jt(i),jto(i))) then
                tvu(i,k) = (su(i,k) - grav / cp * zf(i,k)) * (1 + 0.608_r8 * qu(i,k))
                loc_conv%dcape(i) = loc_conv%dcape(i) + rd * (tvu(i,k) - tvuo(i,k)) * log(p(i,k) / p(i,k-1))
              end if
            end if
          end do
        end do
      else ! No convective microphysics

        ! compute condensed liquid, rain production rate
        ! accumulate total precipitation (condensation - detrainment of liquid)
        ! Note ql1 = ql(k) + rprd(k)*dz(k)/mu(k)
        ! The differencing is somewhat strange (e.g. du(i,k)*ql(i,k+1)) but is
        ! consistently applied.
        !    mu, ql are interface quantities
        !    cu, du, eu, rprd are midpoint quantites

        do k = pver, msg + 2, -1
          do i = 1, il2g
            rprd(i,k) = 0
            if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0 .and. mu(i,k) >= 0) then
              if (mu(i,k) > 0) then
                ql1 = 1.0_r8 / mu(i,k) * (mu(i,k+1) * ql(i,k+1) - &
                      dz(i,k) * du(i,k) * ql(i,k+1) + dz(i,k) * cu(i,k))
                ql(i,k) = ql1 / (1 + dz(i,k) * c0mask(i))
              else
                ql(i,k) = 0
              end if
              totpcp(i) = totpcp(i) + dz(i,k) * (cu(i,k) - du(i,k) * ql(i,k+1))
              rprd(i,k) = c0mask(i) * mu(i,k) * ql(i,k)
              qcde(i,k) = ql(i,k)

              if (zmconv_microp) then
                loc_conv%qide (i,k) = 0
                loc_conv%qncde(i,k) = 0
                loc_conv%qnide(i,k) = 0
                loc_conv%sprd (i,k) = 0
              end if
            end if
          end do
        end do
      end if  ! zmconv_microp
    end do   !iter
    !
    ! specify downdraft properties (no downdrafts if jd.ge.jb).
    ! scale down downward mass flux profile so that net flux
    ! (up-down) at cloud base in not negative.
    !
    do i = 1, il2g
      ! In normal downdraft strength run alfa=0.2.  In test4 alfa=0.1
      alfa(i) = 0.1_r8
      jt(i) = min(jt(i),jb(i)-1)
      jd(i) = max(j0(i),jt(i)+1)
      jd(i) = min(jd(i),jb(i))
      hd(i,jd(i)) = hmn(i,jd(i)-1)
      if (jd(i) < jb(i) .and. eps0(i) > 0) then
        epsm(i) = eps0(i)
        md(i,jd(i)) = -alfa(i) * epsm(i) / eps0(i)
      end if
    end do
    do k = msg + 1, pver
      do i = 1, il2g
         if ((k > jd(i) .and. k <= jb(i)) .and. eps0(i) > 0.0_r8) then
            zdef(i) = zf(i,jd(i)) - zf(i,k)
            md(i,k) = -alfa(i)/ (2._r8*eps0(i))*(exp(2._r8*epsm(i)*zdef(i))-1.0_r8)/zdef(i)
         end if
      end do
    end do

    do k = msg + 1, pver
      do i = 1, il2g
        if (k >= jt(i) .and. k <= jb(i) .and. eps0(i) > 0 .and. jd(i) < jb(i)) then
          ratmjb(i) = min(abs(mu(i,jb(i)) / md(i,jb(i))), 1.0_r8)
          md(i,k) = md(i,k) * ratmjb(i)
        end if
      end do
    end do

    small = 1.0e-20_r8
    do k = msg + 1, pver
      do i = 1, il2g
        if (k >= jt(i) .and. k <= pver .and. eps0(i) > 0.0_r8) then
          ed(i,k-1) = (md(i,k-1) - md(i,k)) / dz(i,k-1)
          mdt = min(md(i,k), -small)
          hd(i,k) = (md(i,k-1) * hd(i,k-1) - dz(i,k-1) * ed(i,k-1) * hmn(i,k-1)) / mdt
        end if
      end do
    end do
    !
    ! Calculate updraft and downdraft properties.
    !
    do k = msg + 2, pver
      do i = 1, il2g
        if (k >= jd(i) .and. k <= jb(i) .and. eps0(i) > 0 .and. jd(i) < jb(i)) then
          qds(i,k) = qsthat(i,k) + gamhat(i,k)*(hd(i,k)-hsthat(i,k))/ (rl*(1.0_r8 + gamhat(i,k)))
        end if
      end do
    end do

    do i = 1, il2g
      qd(i,jd(i)) = qds(i,jd(i))
      sd(i,jd(i)) = (hd(i,jd(i)) - rl*qd(i,jd(i)))/cp
    end do

    do k = msg + 2,pver
      do i = 1, il2g
        if (k >= jd(i) .and. k < jb(i) .and. eps0(i) > 0.0_r8) then
          qd(i,k+1) = qds(i,k+1)
          evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
          evp(i,k) = max(evp(i,k),0.0_r8)
          mdt = min(md(i,k+1),-small)
          if (zmconv_microp) then
            evp(i,k) = min(evp(i,k),rprd(i,k))
          end if
          sd(i,k+1) = ((rl/cp*evp(i,k)-ed(i,k)*s(i,k))*dz(i,k) + md(i,k)*sd(i,k))/mdt
          totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
        end if
      end do
    end do
    do i = 1, il2g
      totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) - md(i,jb(i))*qd(i,jb(i))
    end do
    if (.false.) then
      do i = 1, il2g
        k = jb(i)
        if (eps0(i) > 0.0_r8) then
          evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k))/dz(i,k)
          evp(i,k) = max(evp(i,k),0.0_r8)
          totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
        end if
      end do
    end if

    do i = 1, il2g
      totpcp(i) = max(totpcp(i), 0.0_r8)
      totevp(i) = max(totevp(i), 0.0_r8)
    end do

    do k = msg + 2, pver
      do i = 1, il2g
        if (totevp(i) > 0 .and. totpcp(i) > 0) then
          md(i,k)  = md (i,k)*min(1.0_r8, totpcp(i)/(totevp(i)+totpcp(i)))
          ed(i,k)  = ed (i,k)*min(1.0_r8, totpcp(i)/(totevp(i)+totpcp(i)))
          evp(i,k) = evp(i,k)*min(1.0_r8, totpcp(i)/(totevp(i)+totpcp(i)))
        else
          md (i,k) = 0
          ed (i,k) = 0
          evp(i,k) = 0
        end if
        ! cmeg is the cloud water condensed - rain water evaporated
        ! rprd is the cloud water converted to rain - (rain evaporated)
        cmeg(i,k) = cu(i,k) - evp(i,k)
        rprd(i,k) = rprd(i,k)-evp(i,k)
      end do
    end do

    ! Compute the net precipitation flux across interfaces
    pflx(:il2g,1) = 0
    do k = 2, pverp
      do i = 1, il2g
        pflx(i,k) = pflx(i,k-1) + rprd(i,k-1) * dz(i,k-1)
      end do
    end do

    do k = msg + 1, pver
      do i = 1, il2g
        mc(i,k) = mu(i,k) + md(i,k)
      end do
    end do
  
  end subroutine cldprp

  subroutine closure( &
    lchnk           , &
    q               , &
    t               , &
    p               , &
    z               , &
    s               , &
    tp              , &
    qs              , &
    qu              , &
    su              , &
    mc              , &
    du              , &
    mu              , &
    md              , &
    qd              , &
    sd              , &
    qhat            , &
    shat            , &
    dp              , &
    qstp            , &
    zf              , &
    ql              , &
    dsubcld         , &
    mb              , &
    cape            , &
    tl              , &
    lcl             , &
    lel             , &
    jt              , &
    mx              , &
    il1g            , &
    il2g            , &
    rd              , &
    grav            , &
    cp              , &
    rl              , &
    msg             , &
    capelmt         )
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! <Say what the routine does>
    !
    ! Method:
    ! <Describe the algorithm(s) used in the routine.>
    ! <Also include any applicable external references.>
    !
    ! Author: G. Zhang and collaborators. CCM contact:P. Rasch
    ! This is contributed code not fully standardized by the CCM core group.
    !
    ! this code is very much rougher than virtually anything else in the CCM
    ! We expect to release cleaner code in a future release
    !
    ! the documentation has been enhanced to the degree that we are able
    !
    !-----------------------------------------------------------------------

    integer , intent(in   ) :: lchnk
    real(r8), intent(inout), dimension(pcols,pver  ) :: q       ! Specific humidity
    real(r8), intent(inout), dimension(pcols,pver  ) :: t       ! Temperature
    real(r8), intent(inout), dimension(pcols,pver  ) :: p       ! Pressure (mb)
    real(r8), intent(inout), dimension(pcols       ) :: mb      ! Cloud base mass flux
    real(r8), intent(in   ), dimension(pcols,pver  ) :: z       ! Height (m)
    real(r8), intent(in   ), dimension(pcols,pver  ) :: s       ! Normalized dry static energy
    real(r8), intent(in   ), dimension(pcols,pver  ) :: tp      ! Parcel temperature
    real(r8), intent(in   ), dimension(pcols,pver  ) :: qs      ! Saturated specific humidity
    real(r8), intent(in   ), dimension(pcols,pver  ) :: qu      ! Updraft specific humidity
    real(r8), intent(in   ), dimension(pcols,pver  ) :: su      ! Normalized dry stat energy of updraft
    real(r8), intent(in   ), dimension(pcols,pver  ) :: mc      ! Net convective mass flux
    real(r8), intent(in   ), dimension(pcols,pver  ) :: du      ! Detrainment from updraft
    real(r8), intent(in   ), dimension(pcols,pver  ) :: mu      ! Mass flux of updraft
    real(r8), intent(in   ), dimension(pcols,pver  ) :: md      ! Mass flux of downdraft
    real(r8), intent(in   ), dimension(pcols,pver  ) :: qd      ! Specific humidity of downdraft
    real(r8), intent(in   ), dimension(pcols,pver  ) :: sd      ! Dry static energy of downdraft
    real(r8), intent(in   ), dimension(pcols,pver  ) :: qhat    ! Environment specific humidity at interfaces
    real(r8), intent(in   ), dimension(pcols,pver  ) :: shat    ! Environment normalized dry static energy at interfaces
    real(r8), intent(in   ), dimension(pcols,pver  ) :: dp      ! Pressure thickness of layers
    real(r8), intent(in   ), dimension(pcols,pver  ) :: qstp    ! Specific humidity of parcel
    real(r8), intent(in   ), dimension(pcols,pver+1) :: zf      ! Height of interface levels
    real(r8), intent(in   ), dimension(pcols,pver  ) :: ql      ! Liquid water mixing ratio
    real(r8), intent(in   ), dimension(pcols       ) :: cape    ! available pot. energy of column
    real(r8), intent(in   ), dimension(pcols       ) :: tl
    real(r8), intent(in   ), dimension(pcols       ) :: dsubcld ! thickness of subcloud layer
    integer , intent(in   ), dimension(pcols       ) :: lcl     ! index of lcl
    integer , intent(in   ), dimension(pcols       ) :: lel     ! index of launch leve
    integer , intent(in   ), dimension(pcols       ) :: jt      ! top of updraft
    integer , intent(in   ), dimension(pcols       ) :: mx      ! base of updraft

    real(r8), dimension(pcols,pver) :: dtpdt
    real(r8), dimension(pcols,pver) :: dqsdtp
    real(r8), dimension(pcols,pver) :: dtmdt
    real(r8), dimension(pcols,pver) :: dqmdt
    real(r8), dimension(pcols,pver) :: dboydt
    real(r8), dimension(pcols,pver) :: thetavp
    real(r8), dimension(pcols,pver) :: thetavm
    real(r8), dimension(pcols     ) :: dtbdt
    real(r8), dimension(pcols     ) :: dqbdt
    real(r8), dimension(pcols     ) :: dtldt
    real(r8), dimension(pcols     ) :: dadt
    real(r8) beta
    real(r8) capelmt
    real(r8) cp
    real(r8) debdt
    real(r8) dltaa
    real(r8) eb
    real(r8) grav
    integer i
    integer il1g
    integer il2g
    integer k, kmin, kmax
    integer msg
    real(r8) rd
    real(r8) rl

    ! Change of subcloud layer properties due to convection is
    ! related to cumulus updrafts and downdrafts.
    ! mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
    ! to define betau, betad and f(z).
    ! note that this implies all time derivatives are in effect
    ! time derivatives per unit cloud-base mass flux, i.e. they
    ! have units of 1/mb instead of 1/sec.
    do i = il1g, il2g
      mb(i) = 0
      eb = p(i,mx(i)) * q(i,mx(i)) / (eps1 + q(i,mx(i)))
      dtbdt(i) = (1.0_r8 / dsubcld(i)) * (mu(i,mx(i)) * (shat(i,mx(i)) - su(i,mx(i))) + &
                 md(i,mx(i)) * (shat(i,mx(i)) - sd(i,mx(i))))
      dqbdt(i) = (1.0_r8 / dsubcld(i)) * (mu(i,mx(i)) * (qhat(i,mx(i)) - qu(i,mx(i))) + &
                 md(i,mx(i)) * (qhat(i,mx(i)) - qd(i,mx(i))))
      debdt = eps1*p(i,mx(i)) / (eps1 + q(i,mx(i)))**2 * dqbdt(i)
      dtldt(i) = -2840 * (3.5_r8 / t(i,mx(i)) * dtbdt(i) - debdt/eb) / &
                 (3.5_r8 * log(t(i,mx(i))) - log(eb) - 4.805_r8)**2
    end do
    ! dtmdt and dqmdt are cumulus heating and drying.
    do k = msg + 1, pver
      do i = il1g, il2g
        dtmdt(i,k) = 0
        dqmdt(i,k) = 0
      end do
    end do

    do k = msg + 1, pver - 1
      do i = il1g, il2g
        if (k == jt(i)) then
          dtmdt(i,k) = (1.0_r8 / dp(i,k)) * (mu(i,k+1) * (su(i,k+1) - shat(i,k+1) - &
                        rl / cp * ql(i,k+1)) + md(i,k+1) * (sd(i,k+1) - shat(i,k+1)))
          dqmdt(i,k) = (1.0_r8 / dp(i,k)) * (mu(i,k+1) * (qu(i,k+1) - qhat(i,k+1) + &
                                  ql(i,k+1)) + md(i,k+1) * (qd(i,k+1) - qhat(i,k+1)))
        end if
      end do
    end do

    beta = 0
    do k = msg + 1, pver - 1
      do i = il1g, il2g
        if (k > jt(i) .and. k < mx(i)) then
          dtmdt(i,k) = (mc(i,k) * (shat(i,k) - s(i,k)) + mc(i,k+1) * (s(i,k) - shat(i,k+1))) / &
                      dp(i,k) - rl / cp * du(i,k) * (beta * ql(i,k) + (1-beta) * ql(i,k+1))
          dqmdt(i,k) = (mu(i,k+1) * (qu(i,k+1) - qhat(i,k+1) + cp / rl * (su(i,k+1) - s(i,k))) - &
                        mu(i,k) * (qu(i,k) - qhat(i,k) + cp / rl * (su(i,k) - s(i,k))) + md(i,k+1) * &
                       (qd(i,k+1) - qhat(i,k+1) + cp / rl * (sd(i,k+1) - s(i,k))) - md(i,k) * &
                       (qd(i,k) - qhat(i,k) + cp / rl * (sd(i,k) - s(i,k)))) / dp(i,k) + &
                        du(i,k) * (beta * ql(i,k) + (1 - beta) * ql(i,k+1))
        end if
      end do
    end do

    do k = msg + 1, pver
      do i = il1g, il2g
        if (k >= lel(i) .and. k <= lcl(i)) then
          thetavp(i,k) = tp  (i,k) * (1000.0_r8 / p(i,k))**(rd/cp) * (1 + 1.608_r8 * qstp(i,k) - q(i,mx(i)))
          thetavm(i,k) = t   (i,k) * (1000.0_r8 / p(i,k))**(rd/cp) * (1 + 0.608_r8 * q(i,k))
          dqsdtp (i,k) = qstp(i,k) * (1 + qstp(i,k) / eps1) * eps1 * rl / (rd * tp(i,k)**2)
          !
          ! dtpdt is the parcel temperature change due to change of
          ! subcloud layer properties during convection.
          !
          dtpdt(i,k) = tp(i,k) / (1 + rl / cp * (dqsdtp(i,k) - qstp(i,k) / tp(i,k))) * &
                      (dtbdt(i) / t(i,mx(i)) + rl / cp * (dqbdt(i) / tl(i) - q(i,mx(i)) / &
                      tl(i)**2 * dtldt(i)))
          !
          ! dboydt is the integrand of cape change.
          !
          dboydt(i,k) = ((dtpdt(i,k) / tp(i,k) + 1.0_r8 / (1 + 1.608_r8 * qstp(i,k) - q(i,mx(i)))* &
                        (1.608_r8 * dqsdtp(i,k) * dtpdt(i,k) - dqbdt(i))) - (dtmdt(i,k) / t(i,k) + 0.608_r8 / &
                        (1 + 0.608_r8 * q(i,k)) * dqmdt(i,k))) * grav * thetavp(i,k) / thetavm(i,k)
        end if
      end do
    end do

    do k = msg + 1, pver
      do i = il1g, il2g
        if (k > lcl(i) .and. k < mx(i)) then
          thetavp(i,k) = tp(i,k) * (1000.0_r8 / p(i,k))**(rd / cp) * (1 + 0.608_r8 * q(i,mx(i)))
          thetavm(i,k) = t (i,k) * (1000.0_r8 / p(i,k))**(rd / cp) * (1 + 0.608_r8 * q(i,k))
          !
          ! dboydt is the integrand of cape change.
          !
          dboydt(i,k) = (dtbdt(i) / t(i,mx(i)) + 0.608_r8 / (1 + 0.608_r8 * q(i,mx(i))) * dqbdt(i) - &
                        dtmdt(i,k) / t(i,k) - 0.608_r8/ (1 + 0.608_r8 * q(i,k)) * dqmdt(i,k)) * &
                        grav * thetavp(i,k) / thetavm(i,k)
        end if
      end do
    end do
    !
    ! Buoyant energy change is set to 2/3*excess cape per 3 hours
    !
    dadt(il1g:il2g)  = 0
    kmin = minval(lel(il1g:il2g))
    kmax = maxval(mx(il1g:il2g)) - 1
    do k = kmin, kmax
      do i = il1g, il2g
        if (k >= lel(i) .and. k <= mx(i) - 1) then
          dadt(i) = dadt(i) + dboydt(i,k)* (zf(i,k)-zf(i,k+1))
        end if
      end do
    end do
    do i = il1g, il2g
      dltaa = -1.0_r8 * (cape(i) - capelmt)
      if (dadt(i) /= 0) mb(i) = max(dltaa / tau / dadt(i), 0.0_r8)
    end do

  end subroutine closure

  subroutine q1q2_pjr( &
    lchnk            , &
    dqdt             , &
    dsdt             , &
    q                , &
    qs               , &
    qu               , &
    su               , &
    du               , &
    qhat             , &
    shat             , &
    dp               , &
    mu               , &
    md               , &
    sd               , &
    qd               , &
    ql               , &
    dsubcld          , &
    jt               , &
    mx               , &
    il1g             , &
    il2g             , &
    cp               , &
    rl               , &
    msg              , &
    dl               , &
    evp              , &
    cu               , &
    loc_conv         )

    real(r8), intent(in ) :: cp
    integer , intent(in ) :: lchnk
    integer , intent(in ) :: il1g
    integer , intent(in ) :: il2g
    integer , intent(in ) :: msg
    real(r8), intent(in ), dimension(pcols,pver) :: q
    real(r8), intent(in ), dimension(pcols,pver) :: qs
    real(r8), intent(in ), dimension(pcols,pver) :: qu
    real(r8), intent(in ), dimension(pcols,pver) :: su
    real(r8), intent(in ), dimension(pcols,pver) :: du
    real(r8), intent(in ), dimension(pcols,pver) :: qhat
    real(r8), intent(in ), dimension(pcols,pver) :: shat
    real(r8), intent(in ), dimension(pcols,pver) :: dp
    real(r8), intent(in ), dimension(pcols,pver) :: mu
    real(r8), intent(in ), dimension(pcols,pver) :: md
    real(r8), intent(in ), dimension(pcols,pver) :: sd
    real(r8), intent(in ), dimension(pcols,pver) :: qd
    real(r8), intent(in ), dimension(pcols,pver) :: ql
    real(r8), intent(in ), dimension(pcols,pver) :: evp
    real(r8), intent(in ), dimension(pcols,pver) :: cu
    real(r8), intent(in ), dimension(pcols     ) :: dsubcld
    real(r8), intent(out), dimension(pcols,pver) :: dqdt
    real(r8), intent(out), dimension(pcols,pver) :: dsdt
    real(r8), intent(out), dimension(pcols,pver) :: dl
    type(zm_conv_t) loc_conv

    integer kbm
    integer ktm
    integer jt(pcols)
    integer mx(pcols)

    integer i, k

    real(r8) emc
    real(r8) rl

    do k = msg + 1, pver
      do i = il1g, il2g
        dsdt(i,k) = 0
        dqdt(i,k) = 0
        dl  (i,k) = 0
      end do
    end do

    if (zmconv_microp) then
      do k = msg + 1, pver
        do i = il1g, il2g
          loc_conv%di (i,k) = 0
          loc_conv%dnl(i,k) = 0
          loc_conv%dni(i,k) = 0
        end do
      end do
    end if

    ! Find the highest level top and bottom levels of convection
    ktm = pver
    kbm = pver
    do i = il1g, il2g
      ktm = min(ktm, jt(i))
      kbm = min(kbm, mx(i))
    end do

    do k = ktm,pver-1
      do i = il1g, il2g
        emc = -cu (i,k) & ! Condensation in updraft
              +evp(i,k)   ! Evaporating rain in downdraft
        dsdt(i,k) = -rl / cp * emc &
                    + (+mu(i,k+1) * (su(i,k+1) - shat(i,k+1)) &
                       -mu(i,k  ) * (su(i,k  ) - shat(i,k  )) &
                       +md(i,k+1) * (sd(i,k+1) - shat(i,k+1)) &
                       -md(i,k  ) * (sd(i,k  ) - shat(i,k  )) &
                      ) / dp(i,k)

        if (zmconv_microp) dsdt(i,k) = dsdt(i,k) + latice / cp * loc_conv%frz(i,k)

        dqdt(i,k) = emc + &
                  (+mu(i,k+1) * (qu(i,k+1) - qhat(i,k+1)) &
                   -mu(i,k  ) * (qu(i,k  ) - qhat(i,k  )) &
                   +md(i,k+1) * (qd(i,k+1) - qhat(i,k+1)) &
                   -md(i,k  ) * (qd(i,k  ) - qhat(i,k  )) &
                  ) / dp(i,k)

        dl(i,k) = du(i,k) * ql(i,k+1)

        if (zmconv_microp) then
          loc_conv%di (i,k) = du(i,k) * loc_conv%qide (i,k+1)
          loc_conv%dnl(i,k) = du(i,k) * loc_conv%qncde(i,k+1)
          loc_conv%dni(i,k) = du(i,k) * loc_conv%qnide(i,k+1)
        end if
      end do
    end do

    do k = kbm, pver
      do i = il1g, il2g
        if (k == mx(i)) then
          dsdt(i,k) = (1.0_r8 / dsubcld(i)) * &
                      (-mu(i,k)* (su(i,k)-shat(i,k)) &
                       -md(i,k)* (sd(i,k)-shat(i,k)) &
                      )
          dqdt(i,k) = (1.0_r8 / dsubcld(i)) * &
                      (-mu(i,k)*(qu(i,k)-qhat(i,k)) &
                       -md(i,k)*(qd(i,k)-qhat(i,k)) &
                      )
        else if (k > mx(i)) then
          dsdt(i,k) = dsdt(i,k-1)
          dqdt(i,k) = dqdt(i,k-1)
        end if
      end do
    end do

  end subroutine q1q2_pjr

  subroutine buoyan_dilute( &
    lchnk                 , &
    ncol                  , &
    q                     , &
    t                     , &
    p                     , &
    z                     , &
    pf                    , &
    tp                    , &
    qstp                  , &
    tl                    , &
    rl                    , &
    cape                  , &
    pblt                  , &
    lcl                   , &
    lel                   , &
    lon                   , &
    mx                    , &
    rd                    , &
    grav                  , &
    cp                    , &
    msg                   , &
    tpert                 , &
    org                   , &
    landfrac              )
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Calculates CAPE the lifting condensation level and the convective top
    ! where buoyancy is first -ve.
    !
    ! Method: Calculates the parcel temperature based on a simple constant
    ! entraining plume model. CAPE is integrated from buoyancy.
    ! 09/09/04 - Simplest approach using an assumed entrainment rate for
    !            testing (dmpdp).
    ! 08/04/05 - Swap to convert dmpdz to dmpdp
    !
    ! SCAM Logical Switches - DILUTE:RBN - Now Disabled
    ! ---------------------
    ! switch(1) = .T. - Uses the dilute parcel calculation to obtain tendencies.
    ! switch(2) = .T. - Includes entropy/q changes due to condensate loss and freezing.
    ! switch(3) = .T. - Adds the PBL Tpert for the parcel temperature at all levels.
    !
    ! References:
    ! Raymond and Blythe (1992) JAS
    !
    ! Author:
    ! Richard Neale - September 2004
    !
    !-----------------------------------------------------------------------

    integer , intent(in ) :: lchnk
    integer , intent(in ) :: ncol
    real(r8), intent(in ), dimension(pcols,pver  ) :: q         ! Specific humidity
    real(r8), intent(in ), dimension(pcols,pver  ) :: t         ! Temperature
    real(r8), intent(in ), dimension(pcols,pver  ) :: p         ! Pressure
    real(r8), intent(in ), dimension(pcols,pver  ) :: z         ! Height
    real(r8), intent(in ), dimension(pcols,pver+1) :: pf        ! Pressure at interfaces
    real(r8), intent(in ), dimension(pcols       ) :: pblt      ! Index of pbl depth
    real(r8), intent(in ), dimension(pcols       ) :: tpert     ! Perturbation temperature by pbl processes
    real(r8), intent(in ), dimension(pcols       ) :: landfrac
    real(r8), intent(out), dimension(pcols,pver  ) :: tp        ! Parcel temperature
    real(r8), intent(out), dimension(pcols,pver  ) :: qstp      ! Saturation mixing ratio of parcel (only above lcl, just q below)
    real(r8), intent(out), dimension(pcols       ) :: tl        ! Parcel temperature at lcl
    real(r8), intent(out), dimension(pcols       ) :: cape      ! Convective available potential energy
    integer , intent(out), dimension(pcols       ) :: lcl
    integer , intent(out), dimension(pcols       ) :: lel
    integer , intent(out), dimension(pcols       ) :: lon       ! Level of onset of deep convection
    integer , intent(out), dimension(pcols       ) :: mx        ! Level of max moist static energy

    real(r8), pointer :: org(:,:)                               ! Organization parameter

    real(r8), dimension(pcols,5   ) :: capeten     ! Provisional value of cape
    real(r8), dimension(pcols,pver) :: tv
    real(r8), dimension(pcols,pver) :: tpv
    real(r8), dimension(pcols,pver) :: buoy

    real(r8), dimension(pcols     ) :: a1
    real(r8), dimension(pcols     ) :: a2
    real(r8), dimension(pcols     ) :: estp
    real(r8), dimension(pcols     ) :: pl
    real(r8), dimension(pcols     ) :: plexp
    real(r8), dimension(pcols     ) :: hmax
    real(r8), dimension(pcols     ) :: hmn
    real(r8), dimension(pcols     ) :: y

    logical , dimension(pcols     ) :: plge600
    integer , dimension(pcols     ) :: knt
    integer , dimension(pcols,5   ) :: lelten

    real(r8) cp
    real(r8) e
    real(r8) grav
    real(r8) rd
    real(r8) rl
#ifdef PERGRO
    real(r8) rhd
#endif
    integer i, k, n, msg

    do n = 1, 5
      do i = 1, ncol
        lelten (i,n) = pver
        capeten(i,n) = 0
      end do
    end do

    do i = 1, ncol
      lon (i) = pver
      knt (i) = 0
      lel (i) = pver
      mx  (i) = lon(i)
      cape(i) = 0
      hmax(i) = 0
    end do

    tp  (:ncol,:) = t(:ncol,:)
    qstp(:ncol,:) = q(:ncol,:)

    tv  (:ncol,:) = t (:ncol,:) * (1 + 1.608_r8 * q(:ncol,:)) / (1 + q(:ncol,:))
    tpv (:ncol,:) = tv(:ncol,:)
    buoy(:ncol,:) = 0

    !
    ! Set "launching" level(mx) to be at maximum moist static energy.
    ! Search for this level stops at planetary boundary layer top.
    !
#ifdef PERGRO
    do k = pver, msg + 1,-1
      do i = 1, ncol
        hmn(i) = cp*t(i,k) + grav * z(i,k) + rl * q(i,k)
        !
        ! Reset max moist static energy level when relative difference exceeds 1.e-4
        !
        rhd = (hmn(i) - hmax(i)) / (hmn(i) + hmax(i))
        if (k >= nint(pblt(i)) .and. k <= lon(i) .and. rhd > -1.0e-4_r8) then
          hmax(i) = hmn(i)
          mx  (i) = k
        end if
      end do
    end do
#else
    do k = pver, msg + 1,-1
      do i = 1, ncol
        hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
        if (k >= nint(pblt(i)) .and. k <= lon(i) .and. hmn(i) > hmax(i)) then
          hmax(i) = hmn(i)
          mx  (i) = k
        end if
      end do
    end do
#endif

    ! LCL dilute calculation - initialize to mx(i)
    ! Determine lcl in parcel_dilute and get pl,tl after parcel_dilute
    ! Original code actually sets LCL as level above wher condensate forms.
    ! Therefore in parcel_dilute lcl(i) will be at first level where qsmix < qtmix.

    do i = 1, ncol ! Initialise LCL variables.
      lcl(i) = mx(i)
      tl (i) = t(i,mx(i))
      pl (i) = p(i,mx(i))
    end do

    !
    ! Main buoyancy calculation.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! DILUTE PLUME CALCULATION USING ENTRAINING PLUME !!!
    !!!   RBN 9/9/04                                    !!!

    call parcel_dilute( &
      lchnk           , &
      ncol            , & 
      msg             , &
      mx              , &
      p               , &
      t               , &
      q               , &
      tpert           , &
      tp              , &
      tpv             , &
      qstp            , &
      pl              , &
      tl              , &
      lcl             , &
      org             , &
      landfrac        )


    ! If lcl is above the nominal level of non-divergence (600 mbs),
    ! no deep convection is permitted (ensuing calculations
    ! skipped and cape retains initialized value of zero).
    do i = 1, ncol
      plge600(i) = pl(i) >= 600 ! Just change to always allow buoy calculation.
    end do
    !
    ! Main buoyancy calculation.
    !
    do k = pver, msg + 1,-1
      do i = 1, ncol
        if (k <= mx(i) .and. plge600(i)) then   ! Define buoy from launch level to cloud top.
          tv(i,k) = t(i,k) * (1 + 1.608_r8 * q(i,k)) / (1 + q(i,k))
          buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add  ! +0.5K or not?
        else
          qstp(i,k) = q (i,k)
          tp  (i,k) = t (i,k)
          tpv (i,k) = tv(i,k)
        end if
      end do
    end do

    do k = msg + 2,pver
      do i = 1, ncol
        if (k < lcl(i) .and. plge600(i)) then
          if (buoy(i,k+1) > 0 .and. buoy(i,k) <= 0) then
            knt(i) = min(num_cin,knt(i) + 1)
            lelten(i,knt(i)) = k
          end if
        end if
      end do
    end do
    !
    ! Calculate convective available potential energy (cape).
    !
    do n = 1, num_cin
      do k = msg + 1, pver
        do i = 1, ncol
          if (plge600(i) .and. k <= mx(i) .and. k > lelten(i,n)) then
            capeten(i,n) = capeten(i,n) + rd * buoy(i,k) * log(pf(i,k+1) / pf(i,k))
          end if
        end do
      end do
    end do
    !
    ! Find maximum cape from all possible tentative capes from one sounding,
    ! and use it as the final cape, april 26, 1995
    !
    do n = 1, num_cin
      do i = 1, ncol
        if (capeten(i,n) > cape(i)) then
          cape(i) = capeten(i,n)
          lel (i) = lelten (i,n)
        end if
      end do
    end do
    !
    ! Put lower bound on cape for diagnostic purposes.
    !
    do i = 1, ncol
      cape(i) = max(cape(i), 0.0_r8)
    end do

  end subroutine buoyan_dilute

  subroutine parcel_dilute( &
    lchnk                 , &
    ncol                  , &
    msg                   , &
    klaunch               , &
    p                     , &
    t                     , &
    q                     , &
    tpert                 , &
    tp                    , &
    tpv                   , &
    qstp                  , &
    pl                    , &
    tl                    , &
    lcl                   , &
    org                   , &
    landfrac              )

    ! Routine  to determine
    !   1. Tp   - Parcel temperature
    !   2. qstp - Saturated mixing ratio at the parcel temperature.

    integer , intent(in   ) :: lchnk
    integer , intent(in   ) :: ncol
    integer , intent(in   ) :: msg
    integer , intent(in   ), dimension(pcols     ) :: klaunch(pcols)
    real(r8), intent(in   ), dimension(pcols,pver) :: p
    real(r8), intent(in   ), dimension(pcols,pver) :: t
    real(r8), intent(in   ), dimension(pcols,pver) :: q
    real(r8), intent(in   ), dimension(pcols     ) :: tpert     ! PBL temperature perturbation.
    real(r8), intent(in   ), dimension(pcols     ) :: landfrac
    real(r8), intent(inout), dimension(pcols,pver) :: tp        ! Parcel temperature
    real(r8), intent(inout), dimension(pcols,pver) :: qstp      ! Parcel water vapour (sat value above lcl).
    real(r8), intent(inout), dimension(pcols     ) :: tl        ! Actual temperature of LCL.
    real(r8), intent(inout), dimension(pcols     ) :: pl        ! Actual pressure of LCL.
    integer , intent(inout), dimension(pcols     ) :: lcl       ! Lifting condesation level (first model level with saturation).
    real(r8), intent(  out), dimension(pcols,pver) :: tpv       ! Define tpv within this routine.
    real(r8), pointer, dimension(:,:) :: org

    ! Have to be careful as s is also dry static energy.

    ! If we are to retain the fact that CAM loops over grid-points in the internal
    ! loop then we need to dimension sp,atp,mp,xsh2o with ncol.

    real(r8), dimension(pcols,pver) :: tmix       ! Tempertaure of the entraining parcel.
    real(r8), dimension(pcols,pver) :: qtmix      ! Total water of the entraining parcel.
    real(r8), dimension(pcols,pver) :: qsmix      ! Saturated mixing ratio at the tmix.
    real(r8), dimension(pcols,pver) :: smix       ! Entropy of the entraining parcel.
    real(r8), dimension(pcols,pver) :: xsh2o      ! Precipitate lost from parcel.
    real(r8), dimension(pcols,pver) :: ds_xsh2o   ! Entropy change due to loss of condensate.
    real(r8), dimension(pcols,pver) :: ds_freeze  ! Entropy change sue to freezing of precip.
    real(r8), dimension(pcols,pver) :: dmpdz2d    ! variable detrainment rate
    real(r8), dimension(pcols     ) :: mp         ! Parcel mass flux.
    real(r8), dimension(pcols     ) :: qtp        ! Parcel total water.
    real(r8), dimension(pcols     ) :: sp         ! Parcel entropy.
    real(r8), dimension(pcols     ) :: sp0        ! Parcel launch entropy.
    real(r8), dimension(pcols     ) :: qtp0       ! Parcel launch total water.
    real(r8), dimension(pcols     ) :: mp0        ! Parcel launch relative mass flux.
    real(r8) lwmax                                ! Maximum condesate that can be held in cloud before rainout.
    real(r8) dmpdp                                ! Parcel fractional mass entrainment rate (/mb).
    real(r8) dmpdz                                ! Parcel fractional mass entrainment rate (/m)
    real(r8) dpdz,dzdp                            ! Hydrstatic relation and inverse of.
    real(r8) senv                                 ! Environmental entropy at each grid point.
    real(r8) qtenv                                ! Environmental total water "   "   ".
    real(r8) penv                                 ! Environmental total pressure "   "   ".
    real(r8) tenv                                 ! Environmental total temperature "   "   ".
    real(r8) new_s                                ! Hold value for entropy after condensation/freezing adjustments.
    real(r8) new_q                                ! Hold value for total water after condensation/freezing adjustments.
    real(r8) dp                                   ! Layer thickness (center to center)
    real(r8) tfguess                              ! First guess for entropy inversion - crucial for efficiency!
    real(r8) tscool                               ! Super cooled temperature offset (in degC) (eg -35).
    real(r8) qxsk, qxskp1                         ! LCL excess water (k, k+1)
    real(r8) dsdp, dqtdp, dqxsdp                  ! LCL s, qt, p gradients (k, k+1)
    real(r8) slcl,qtlcl,qslcl                     ! LCL s, qt, qs values.
    real(r8) org2rkm, org2Tpert
    real(r8) dmpdz_lnd, dmpdz_mask
    integer rcall                                 ! Number of ientropy call for errors recording
    integer nit_lheat                             ! Number of iterations for condensation/freezing loop.
    integer i, k, ii

    !======================================================================
    !    SUMMARY
    !
    !  9/9/04 - Assumes parcel is initiated from level of maxh (klaunch)
    !           and entrains at each level with a specified entrainment rate.
    !
    ! 15/9/04 - Calculates lcl(i) based on k where qsmix is first < qtmix.
    !
    !======================================================================
    !
    ! Set some values that may be changed frequently.
    !
    if (zm_org) then
      org2rkm   = 10
      org2Tpert = 0
    end if
    nit_lheat =  2
    dmpdz     = -1.e-3_r8    ! Entrainment rate. (-ve for /m)
    dmpdz_lnd = -1.e-3_r8
    lwmax     =  1.e-3_r8    ! Need to put formula in for this.
    tscool    = 0            ! Temp at which water loading freezes in the cloud.
    qtmix     = 0
    smix      = 0
    qtenv     = 0
    senv      = 0
    tenv      = 0
    penv      = 0
    qtp0      = 0
    sp0       = 0
    mp0       = 0
    qtp       = 0
    sp        = 0
    mp        = 0
    new_q     = 0
    new_s     = 0

    do k = pver, msg+1, -1
      do i = 1, ncol
        ! Initialize parcel values at launch level.
        if (k == klaunch(i)) then
          qtp0  (i  ) = q(i,k)   ! Parcel launch total water (assuming subsaturated) - OK????.
          sp0   (i  ) = entropy(t(i,k),p(i,k),qtp0(i))  ! Parcel launch entropy.
          mp0   (i  ) = 1       ! Parcel launch relative mass (i.e. 1 parcel stays 1 parcel for dmpdp=0, undilute).
          smix  (i,k) = sp0(i)
          qtmix (i,k) = qtp0(i)
          tfguess = t(i,k)
          rcall = 1
          call ientropy(rcall, i, lchnk, smix(i,k), p(i,k), qtmix(i,k), tmix(i,k), qsmix(i,k), tfguess)
        end if
        ! Entraining levels
        if (k < klaunch(i)) then
          ! Set environmental values for this level.
          dp    = p(i,k) - p(i,k+1)          ! In -ve mb as p decreasing with height - difference between center of layers.
          qtenv = 0.5_r8*(q(i,k) + q(i,k+1)) ! Total water of environment.
          tenv  = 0.5_r8*(t(i,k) + t(i,k+1))
          penv  = 0.5_r8*(p(i,k) + p(i,k+1))
          senv  = entropy(tenv, penv, qtenv) ! Entropy of environment.
          ! Determine fractional entrainment rate /pa given value /m.
          dpdz = -(penv * grav) / (rgas * tenv) ! in mb/m since  p in mb.
          dzdp = 1.0_r8 / dpdz                  ! in m/mb
          if (zm_org) then
            dmpdz_mask = landfrac(i) * dmpdz_lnd + (1 - landfrac(i)) * dmpdz
            dmpdp = (dmpdz_mask / (1 + org(i,k) * org2rkm)) * dzdp              ! /mb Fractional entrainment
          else
            dmpdp = dmpdz * dzdp
          end if
          ! Sum entrainment to current level
          ! entrains q,s out of intervening dp layers, in which linear variation is assumed
          ! so really it entrains the mean of the 2 stored values.
          sp (i) = sp (i) - dmpdp * dp * senv
          qtp(i) = qtp(i) - dmpdp * dp * qtenv
          mp (i) = mp (i) - dmpdp * dp
          ! Entrain s and qt to next level.
          smix (i,k) = (sp0 (i) + sp (i)) / (mp0(i) + mp(i))
          qtmix(i,k) = (qtp0(i) + qtp(i)) / (mp0(i) + mp(i))
          ! Invert entropy from s and q to determine T and saturation-capped q of mixture.
          ! t(i,k) used as a first guess so that it converges faster.
          tfguess = tmix(i,k+1)
          rcall = 2
          call ientropy(rcall, i, lchnk, smix(i,k), p(i,k), qtmix(i,k), tmix(i,k), qsmix(i,k), tfguess)
          ! Determine if this is lcl of this column if qsmix <= qtmix.
          ! FIRST LEVEL where this happens on ascending.
          if (qsmix(i,k) <= qtmix(i,k) .and. qsmix(i,k+1) > qtmix(i,k+1)) then
            lcl(i)  = k
            qxsk    = qtmix(i,k) - qsmix(i,k)
            qxskp1  = qtmix(i,k+1) - qsmix(i,k+1)
            dqxsdp  = (qxsk - qxskp1)/dp
            pl(i)   = p(i,k+1) - qxskp1/dqxsdp    ! pressure level of actual lcl.
            dsdp    = (smix(i,k)  - smix(i,k+1))/dp
            dqtdp   = (qtmix(i,k) - qtmix(i,k+1))/dp
            slcl    = smix(i,k+1)  +  dsdp* (pl(i)-p(i,k+1))
            qtlcl   = qtmix(i,k+1) +  dqtdp*(pl(i)-p(i,k+1))
            tfguess = tmix(i,k)
            rcall   = 3
            call ientropy(rcall,i,lchnk,slcl,pl(i),qtlcl,tl(i),qslcl,tfguess)
          end if
        end if
      end do
    end do

    ! Could stop now and test with this as it will provide some estimate of buoyancy
    ! without the effects of freezing/condensation taken into account for tmix.

    ! So we now have a profile of entropy and total water of the entraining parcel
    ! Varying with height from the launch level klaunch parcel=environment. To the
    ! top allowed level for the existence of convection.

    ! Now we have to adjust these values such that the water held in vaopor is < or
    ! = to qsmix. Therefore, we assume that the cloud holds a certain amount of
    ! condensate (lwmax) and the rest is rained out (xsh2o). This, obviously
    ! provides latent heating to the mixed parcel and so this has to be added back
    ! to it. But does this also increase qsmix as well? Also freezing processes

    xsh2o     = 0
    ds_xsh2o  = 0
    ds_freeze = 0

    ! Iterate solution twice for accuracy
    do k = pver, msg + 1, -1
      do i = 1, ncol
        ! Initialize variables at k=klaunch
        if (k == klaunch(i)) then
          ! Set parcel values at launch level assume no liquid water.
          tp  (i,k) = tmix(i,k)
          qstp(i,k) = q   (i,k)
          if (zm_org) then
            tpv(i,k) = (tp(i,k) + (org2Tpert * org(i,k) + tpert(i))) * (1 + 1.608_r8 * qstp(i,k)) / (1 + qstp(i,k))
          else
            tpv(i,k) = (tp(i,k) + tpert(i)) * (1 + 1.608_r8 * qstp(i,k)) / (1 + qstp(i,k))
          end if
        end if
        if (k < klaunch(i)) then
          ! Initiaite loop if switch(2) = .T. - RBN:DILUTE - TAKEN OUT BUT COULD BE RETURNED LATER.
          ! Iterate nit_lheat times for s,qt changes.
          do ii = 0, nit_lheat - 1
            ! Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).
            xsh2o(i,k) = max(0.0_r8, qtmix(i,k) - qsmix(i,k) - lwmax)
            ! Contribution to ds from precip loss of condensate (Accumulated change from smix).(-ve)
            ds_xsh2o(i,k) = ds_xsh2o(i,k+1) - cpliq * log(tmix(i,k) / tfreez) * max(0.0_r8, (xsh2o(i,k) - xsh2o(i,k+1)))
            !
            ! Entropy of freezing: latice times amount of water involved divided by T.
            !
            if (tmix(i,k) <= tfreez+tscool .and. ds_freeze(i,k+1) == 0) then ! One off freezing of condensate.
              ds_freeze(i,k) = (latice / tmix(i,k)) * max(0.0_r8, qtmix(i,k) - qsmix(i,k) - xsh2o(i,k)) ! Gain of LH
            end if
            if (tmix(i,k) <= tfreez+tscool .and. ds_freeze(i,k+1) /= 0) then ! Continual freezing of additional condensate.
               ds_freeze(i,k) = ds_freeze(i,k+1) + (latice / tmix(i,k)) * max(0.0_r8, qsmix(i,k+1) - qsmix(i,k))
            end if
            ! Adjust entropy and accordingly to sum of ds (be careful of signs).
            new_s = smix(i,k) + ds_xsh2o(i,k) + ds_freeze(i,k)
            ! Adjust liquid water and accordingly to xsh2o.
            new_q = qtmix(i,k) - xsh2o(i,k)
            ! Invert entropy to get updated Tmix and qsmix of parcel.
            tfguess = tmix(i,k)
            rcall = 4
            call ientropy(rcall, i, lchnk, new_s, p(i,k), new_q, tmix(i,k), qsmix(i,k), tfguess)
          end do  ! Iteration loop for freezing processes.

          ! tp  - Parcel temp is temp of mixture.
          ! tpv - Parcel v. temp should be density temp with new_q total water.
          tp(i,k) = tmix(i,k)

          ! tpv = tprho in the presence of condensate (i.e. when new_q > qsmix)
          if (new_q > qsmix(i,k)) then  ! Super-saturated so condensate present - reduces buoyancy.
            qstp(i,k) = qsmix(i,k)
          else                          ! Just saturated/sub-saturated - no condensate virtual effects.
            qstp(i,k) = new_q
          end if

          if (zm_org) then
            tpv(i,k) = (tp(i,k)+(org2Tpert*org(i,k)+tpert(i)))* (1.0_r8+1.608_r8*qstp(i,k)) / (1.0_r8+ new_q)
          else
            tpv(i,k) = (tp(i,k)+tpert(i))* (1.0_r8+1.608_r8*qstp(i,k)) / (1.0_r8+ new_q)
          end if
        end if ! k < klaunch
      end do
    end do

  end subroutine parcel_dilute

  real(r8) function entropy(t, p, qtot)

    ! from Raymond and Blyth 1992

    real(r8), intent(in) :: t     ! Temperature (K)
    real(r8), intent(in) :: p     ! Pressure (hPa)
    real(r8), intent(in) :: qtot  ! Total water mixing ratio (kg/kg)

    real(r8) qv, qst, e, est, L
    real(r8), parameter :: pref = 1000.0_r8 ! Reference pressure (hPa)

    L = rl - (cpliq - cpwv) * (t - tfreez)

    call qsat_hPa(t, p, est, qst)

    qv = min(qtot, qst)              ! Partition qtot into vapor part only.
    e  = qv * p / (eps1 + qv)

    entropy = (cpres + qtot * cpliq) * log(t / tfreez) - rgas * log((p - e) / pref) + &
              L * qv / t - qv * rh2o * log(qv / qst)

  end function entropy

  subroutine ientropy(rcall, icol, lchnk, s, p, qt, t, qst, tfg)

    ! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg).
    ! Inverts entropy, pressure and total water qt
    ! for T and saturated vapor mixing ratio

    use phys_grid, only: get_rlon_p, get_rlat_p

    integer , intent(in ) :: rcall
    integer , intent(in ) :: icol
    integer , intent(in ) :: lchnk
    real(r8), intent(in ) :: s
    real(r8), intent(in ) :: p
    real(r8), intent(in ) :: tfg
    real(r8), intent(in ) :: qt
    real(r8), intent(out) :: t
    real(r8), intent(out) :: qst

    real(r8) est, this_lat, this_lon
    real(r8) a, b, c, d, ebr, fa, fb, fc, pbr, qbr, rbr, sbr, tol1, xm, tol
    integer i
    logical converged

    ! Max number of iteration loops.
    integer , parameter :: loopmax = 100
    real(r8), parameter :: eps = 3.0e-8_r8

    converged = .false.

    ! Invert the entropy equation -- use Brent's method
    ! Brent, R. P. Ch. 3-4 in Algorithms for Minimization Without Derivatives. Englewood Cliffs, NJ: Prentice-Hall, 1973.

    t = tfg         ! Better first guess based on Tprofile from conv.

    a = tfg - 10    ! Low bracket
    b = tfg + 10    ! High bracket

    fa = entropy(a, p, qt) - s
    fb = entropy(b, p, qt) - s

    c   = b
    fc  = fb
    tol = 0.001_r8

    converge: do i = 0, loopmax
      if ((fb > 0.0_r8 .and. fc > 0.0_r8) .or. &
          (fb < 0.0_r8 .and. fc < 0.0_r8)) then
        c   = a
        fc  = fa
        d   = b - a
        ebr = d
      end if
      if (abs(fc) < abs(fb)) then
        a  = b
        b  = c
        c  = a
        fa = fb
        fb = fc
        fc = fa
      end if

      tol1      = 2.0_r8 * eps * abs(b) + 0.5_r8 * tol
      xm        = 0.5_r8 * (c - b)
      converged = (abs(xm) <= tol1 .or. fb == 0)
      if (converged) exit converge

      if (abs(ebr) >= tol1 .and. abs(fa) > abs(fb)) then
        sbr = fb / fa
        if (a == c) then
          pbr = 2 * xm * sbr
          qbr = -sbr
        else
          qbr = fa / fc
          rbr = fb / fc
          pbr = sbr * (2 * xm * qbr * (qbr - rbr) - (b - a) * (rbr - 1))
          qbr = (qbr - 1) * (rbr - 1) * (sbr - 1)
        end if
        if (pbr > 0) qbr = -qbr
        pbr = abs(pbr)
        if (2 * pbr < min(3 * xm * qbr - abs(tol1 * qbr), abs(ebr*qbr))) then
          ebr = d
          d   = pbr / qbr
        else
          d   = xm
          ebr = d
        end if
      else
        d   = xm
        ebr = d
      end if
      a  = b
      fa = fb
      b  = b + merge(d, sign(tol1, xm), abs(d) > tol1)
      fb = entropy(b, p, qt) - s
    end do converge

    T = b
    call qsat_hPa(T, p, est, qst)

    if (.not. converged) then
      this_lat = get_rlat_p(lchnk, icol)*57.296_r8
      this_lon = get_rlon_p(lchnk, icol)*57.296_r8
      write(iulog, *) '*** ZM_CONV: IENTROPY: Failed and about to exit, info follows ****'
      write(iulog, '(A,I1,I4,I4,7(A,F6.2))') &
        'ZM_CONV: IENTROPY. Details: call#,lchnk,icol= ', rcall, lchnk, icol, &
        ' lat: ', this_lat, ' lon: ', this_lon, &
        ' P(mb)= ', p, ' Tfg(K)= ', Tfg, ' qt(g/kg) = ', 1000 * qt, &
        ' qst(g/kg) = ', 1000 * qst, ', s(J/kg) = ', s
      call endrun('**** ZM_CONV IENTROPY: Tmix did not converge ****')
    end if

  end subroutine ientropy

  ! Wrapper for qsat_water that does translation between Pa and hPa
  ! qsat_water uses Pa internally, so get it right, need to pass in Pa.
  ! Afterward, set es back to hPa.
  elemental subroutine qsat_hPa(t, p, es, qm)

    use wv_saturation, only: qsat_water

    real(r8), intent(in ) :: t   ! Temperature (K)
    real(r8), intent(in ) :: p   ! Pressure (hPa)
    real(r8), intent(out) :: es  ! Saturation vapor pressure (hPa)
    real(r8), intent(out) :: qm  ! Saturation mass mixing ratio
                                 ! (vapor mass over dry mass, kg/kg)

    call qsat_water(t, p * 100, es, qm)

    es = es * 0.01_r8

  end subroutine qsat_hPa

end module zm_conv
