! ==============================================================================
! This file is part of GoMars since 2023.
!
! GoMars is a Martian general circulation model developed in Institute of
! Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
!
! GMCORE is a dynamical core for atmospheric model used in GoMars.
!
! GoMars and GMCORE are distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module gomars_v1_types_mod

  use physics_types_mod
  use tracer_mod
  use gomars_v1_const_mod
  use gomars_v1_rad_mod

  implicit none

  private

  public gomars_v1_state_type
  public gomars_v1_tend_type
  public physics_use_wet_tracers

  type, extends(physics_state_type) :: gomars_v1_state_type
    ! Surface pressure at time level n (Pa)
    real(r8), allocatable, dimension(:      ) :: ps_old
    ! Pressure of planetary boundary layer top (Pa)
    real(r8), allocatable, dimension(:      ) :: ptop_pbl
    ! Stratospheric temperature (K)
    real(r8), allocatable, dimension(:      ) :: tstrat
    ! Surface CO2 ice (?)
    real(r8), allocatable, dimension(:      ) :: co2ice_sfc
    !
    real(r8), allocatable, dimension(:      ) :: latheat
    !
    real(r8), allocatable, dimension(:,:    ) :: atmcond
    ! Tracer mass on the surface (kg m-2)
    real(r8), allocatable, dimension(:,    :) :: tm_sfc
    ! Visible extinction efficiency for dust
    real(r8), allocatable, dimension(  :,:  ) :: qxvdst
    ! Visible scattering efficiency for dust
    real(r8), allocatable, dimension(  :,:  ) :: qsvdst
    ! Visible asymmetry parameter for dust
    real(r8), allocatable, dimension(  :,:  ) :: gvdst
    ! Infrared extinction efficiency for dust
    real(r8), allocatable, dimension(  :,:  ) :: qxidst
    ! Infrared scattering efficiency for dust
    real(r8), allocatable, dimension(  :,:  ) :: qsidst
    ! Infrared asymmetry parameter for dust
    real(r8), allocatable, dimension(  :,:  ) :: gidst
    !
    real(r8), allocatable, dimension(  :    ) :: qextrefdst
    ! Visible extinction efficiency for ice clouds
    real(r8), allocatable, dimension(  :,:  ) :: qxvcld
    ! Visible scattering efficiency for ice clouds
    real(r8), allocatable, dimension(  :,:  ) :: qsvcld
    ! Visible asymmetry parameter for ice clouds
    real(r8), allocatable, dimension(  :,:  ) :: gvcld
    ! Infrared extinction efficiency for ice clouds
    real(r8), allocatable, dimension(  :,:  ) :: qxicld
    ! Infrared scattering efficiency for ice clouds
    real(r8), allocatable, dimension(  :,:  ) :: qsicld
    ! Infrared asymmetry parameter for ice clouds
    real(r8), allocatable, dimension(  :,:  ) :: gicld
    !
    real(r8), allocatable, dimension(  :    ) :: qextrefcld
    !
    real(r8), allocatable, dimension(  :,:,:) :: wbarv
    !
    real(r8), allocatable, dimension(  :,:,:) :: cosbv
    !
    real(r8), allocatable, dimension(  :,:,:) :: wbari
    !
    real(r8), allocatable, dimension(  :,:,:) :: cosbi
    !
    real(r8), allocatable, dimension(  :    ) :: fluxupv
    !
    real(r8), allocatable, dimension(  :    ) :: fluxdnv
    !
    real(r8), allocatable, dimension(  :    ) :: fmnetv
    !
    real(r8), allocatable, dimension(  :    ) :: fluxupi
    !
    real(r8), allocatable, dimension(  :    ) :: fluxdni
    !
    real(r8), allocatable, dimension(  :    ) :: fmneti
    !
    real(r8), allocatable, dimension(:      ) :: fuptopv
    !
    real(r8), allocatable, dimension(:      ) :: fdntopv
    !
    real(r8), allocatable, dimension(:      ) :: fupsfcv
    !
    real(r8), allocatable, dimension(:      ) :: fdnsfcv
    !
    real(r8), allocatable, dimension(:      ) :: fuptopi
    !
    real(r8), allocatable, dimension(:      ) :: fupsfci
    !
    real(r8), allocatable, dimension(:      ) :: fdnsfci
    ! [restart] Downward diffuse visible (solar) flux at the surface (W m-2)
    real(r8), allocatable, dimension(:      ) :: vsdif_sfc_dn
    ! Downward visible flux at the surface (W m-2)
    real(r8), allocatable, dimension(:      ) :: vsflx_sfc_dn
    ! Downward IR flux at the surface (W m-2)
    real(r8), allocatable, dimension(:      ) :: irflx_sfc_dn
    ! Reference optical depth for dust in 0.67 um
    real(r8), allocatable, dimension(  :    ) :: taurefdst
    !
    real(r8), allocatable, dimension(  :    ) :: taurefcld
    !
    real(r8), allocatable, dimension(  :    ) :: taucum
    !
    real(r8), allocatable, dimension(    :,:) :: taugsurf
    !
    real(r8), allocatable, dimension(:      ) :: tausurf
    !
    real(r8), allocatable, dimension(:,    :) :: taudst
    !
    real(r8), allocatable, dimension(:,    :) :: taucld
    !
    real(r8), allocatable, dimension(  :,:,:) :: dtauv
    !
    real(r8), allocatable, dimension(  :,:,:) :: tauv
    !
    real(r8), allocatable, dimension(  :,:,:) :: taucumv
    !
    real(r8), allocatable, dimension(  :,:,:) :: dtaui
    !
    real(r8), allocatable, dimension(  :,:,:) :: taui
    !
    real(r8), allocatable, dimension(  :,:,:) :: taucumi
    ! Delta-Eddington optical depth (???)
    real(r8), allocatable, dimension(:,  :,:) :: detau
    !
    real(r8), allocatable, dimension(  :    ) :: suntot
    !
    real(r8), allocatable, dimension(  :    ) :: irtot
    ! Downward solar flux at the surface (W m-2)
    real(r8), allocatable, dimension(    :  ) :: solar_sfc_dn   ! solar
    ! Total absorption of solar energy by the atmosphere (?)
    real(r8), allocatable, dimension(:      ) :: ssun
    !
    real(r8), allocatable, dimension(:      ) :: fluxsfc
    ! Heat rate on full levels (including TOA) due to radiation (K s-1)
    real(r8), allocatable, dimension(:,:    ) :: ht_rad         ! qrad
    ! Heat rate at the surface due to PBL (K s-1)
    real(r8), allocatable, dimension(:      ) :: ht_pbl
    ! Exchange of heat between surface and air (positive downward) (W m-2)
    real(r8), allocatable, dimension(:      ) :: ht_sfc         ! fa
    ! Surface albedo from data
    real(r8), allocatable, dimension(:      ) :: alsp
    ! Surface albedo considering surface CO2 ice
    real(r8), allocatable, dimension(:      ) :: als
    ! Polar cap flag
    logical , allocatable, dimension(:      ) :: polarcap
    ! cpd * rho * ustar * cdh
    real(r8), allocatable, dimension(:      ) :: rhouch
    ! Squared wind shear (s-1)
    real(r8), allocatable, dimension(:,:    ) :: shr2
    ! Gradient Richardson number
    real(r8), allocatable, dimension(:,:    ) :: ri
    ! Eddy mixing coefficient for momentum (m2 s-1)
    real(r8), allocatable, dimension(:,:    ) :: km
    ! Eddy mixing coefficient for heat (m2 s-1)
    real(r8), allocatable, dimension(:,:    ) :: kh
    ! Thermal inertia (J m-2 K-1 s-1/2)
    real(r8), allocatable, dimension(:,:    ) :: zin
    real(r8), allocatable, dimension(:,:    ) :: rhosoil
    real(r8), allocatable, dimension(:,:    ) :: cpsoil
    !
    real(r8), allocatable, dimension(:,:    ) :: scond
    !
    real(r8), allocatable, dimension(:,:    ) :: stemp
    ! Water ice upward sublimation amount per unit area at the surface (kg m-2)
    real(r8), allocatable, dimension(:      ) :: h2osub_sfc     ! subflux
    ! [restart] Water ice at the surface (???)
    real(r8), allocatable, dimension(:      ) :: h2oice_sfc     ! gndice
    !
    real(r8), allocatable, dimension(:      ) :: dmadt
    ! Wind stress dust lifting flux (kg m-2 s-1)
    real(r8), allocatable, dimension(:      ) :: dstflx_wsl
    ! Dust devil lifting flux (kg m-2 s-1)
    real(r8), allocatable, dimension(:      ) :: dstflx_ddl
    ! Square of Mars distance from sun
    real(r8) :: rsdist
    ! Temperature at all levels with two extra levels at the top (K)
    real(r8), allocatable, dimension(    :  ) :: tl
    ! Temperature at radiation levels (K)
    real(r8), allocatable, dimension(    :  ) :: tlev_rad
    ! Temperature at radiation levels (K)
    real(r8), allocatable, dimension(    :  ) :: tmid_rad
    ! Potential temperature at all levels with two extra levels at the top (K)
    real(r8), allocatable, dimension(    :  ) :: teta
    ! Pressure at all levels with two extra levels at the top (Pa)
    real(r8), allocatable, dimension(    :  ) :: pl
    real(r8), allocatable, dimension(    :  ) :: plogadj
    ! Pressure at radiation levels (Pa)
    real(r8), allocatable, dimension(    :  ) :: plev_rad
    ! Pressure at radiation levels (Pa)
    real(r8), allocatable, dimension(    :  ) :: pmid_rad
    !
    real(r8), allocatable, dimension(    :  ) :: aadj
    !
    real(r8), allocatable, dimension(    :  ) :: badj
    ! Exner pressure (1)
    real(r8), allocatable, dimension(    :  ) :: om
    ! Water mixing ratio on the full and half levels (kg kg-1)
    real(r8), allocatable, dimension(    :  ) :: qh2o
  contains
    procedure :: init  => gomars_v1_state_init
    procedure :: clear => gomars_v1_state_clear
    final gomars_v1_state_final
  end type gomars_v1_state_type

  type, extends(physics_tend_type) :: gomars_v1_tend_type
    ! Wind speed tendency due to planetary boundary layer (m s-1 s-1)
    real(r8), allocatable, dimension(:,:) :: dudt_pbl
  contains
    procedure :: init  => gomars_v1_tend_init
    procedure :: clear => gomars_v1_tend_clear
    final gomars_v1_tend_final
  end type gomars_v1_tend_type

contains

  subroutine gomars_v1_state_init(this, mesh)

    class(gomars_v1_state_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%clear()

    allocate(this%ps_old        (mesh%ncol                   )); this%ps_old        = 0
    allocate(this%ptop_pbl      (mesh%ncol                   )); this%ptop_pbl      = 0
    allocate(this%tstrat        (mesh%ncol                   )); this%tstrat        = 0
    allocate(this%co2ice_sfc    (mesh%ncol                   )); this%co2ice_sfc    = 0
    allocate(this%latheat       (mesh%ncol                   )); this%latheat       = 0
    allocate(this%atmcond       (mesh%ncol,mesh%nlev         )); this%atmcond       = 0
    allocate(this%tm_sfc        (mesh%ncol,ntracers          )); this%tm_sfc        = 0
    allocate(this%qxvdst        (2*mesh%nlev+4,nspectv       )); this%qxvdst        = 0
    allocate(this%qsvdst        (2*mesh%nlev+4,nspectv       )); this%qsvdst        = 0
    allocate(this%gvdst         (2*mesh%nlev+4,nspectv       )); this%gvdst         = 0
    allocate(this%qxidst        (2*mesh%nlev+4,nspecti       )); this%qxidst        = 0
    allocate(this%qsidst        (2*mesh%nlev+4,nspecti       )); this%qsidst        = 0
    allocate(this%gidst         (2*mesh%nlev+4,nspecti       )); this%gidst         = 0
    allocate(this%qextrefdst    (2*mesh%nlev+4               )); this%qextrefdst    = 0
    allocate(this%qxvcld        (2*mesh%nlev+4,nspectv       )); this%qxvcld        = 0
    allocate(this%qsvcld        (2*mesh%nlev+4,nspectv       )); this%qsvcld        = 0
    allocate(this%gvcld         (2*mesh%nlev+4,nspectv       )); this%gvcld         = 0
    allocate(this%qxicld        (2*mesh%nlev+4,nspecti       )); this%qxicld        = 0
    allocate(this%qsicld        (2*mesh%nlev+4,nspecti       )); this%qsicld        = 0
    allocate(this%gicld         (2*mesh%nlev+4,nspecti       )); this%gicld         = 0
    allocate(this%qextrefcld    (2*mesh%nlev+4               )); this%qextrefcld    = 0
    allocate(this%wbarv         (nlayrad,nspectv,ngauss      )); this%wbarv         = 0
    allocate(this%cosbv         (nlayrad,nspectv,ngauss      )); this%cosbv         = 0
    allocate(this%wbari         (nlayrad,nspecti,ngauss      )); this%wbari         = 0
    allocate(this%cosbi         (nlayrad,nspecti,ngauss      )); this%cosbi         = 0
    allocate(this%fluxupv       (nlayrad                     )); this%fluxupv       = 0
    allocate(this%fluxdnv       (nlayrad                     )); this%fluxdnv       = 0
    allocate(this%fmnetv        (nlayrad                     )); this%fmnetv        = 0
    allocate(this%fluxupi       (nlayrad                     )); this%fluxupi       = 0
    allocate(this%fluxdni       (nlayrad                     )); this%fluxdni       = 0
    allocate(this%fmneti        (nlayrad                     )); this%fmneti        = 0
    allocate(this%fuptopv       (mesh%ncol                   )); this%fuptopv       = 0
    allocate(this%fdntopv       (mesh%ncol                   )); this%fdntopv       = 0
    allocate(this%fupsfcv       (mesh%ncol                   )); this%fupsfcv       = 0
    allocate(this%fdnsfcv       (mesh%ncol                   )); this%fdnsfcv       = 0
    allocate(this%fuptopi       (mesh%ncol                   )); this%fuptopi       = 0
    allocate(this%fupsfci       (mesh%ncol                   )); this%fupsfci       = 0
    allocate(this%fdnsfci       (mesh%ncol                   )); this%fdnsfci       = 0
    allocate(this%vsdif_sfc_dn  (mesh%ncol                   )); this%vsdif_sfc_dn  = 0
    allocate(this%vsflx_sfc_dn  (mesh%ncol                   )); this%vsflx_sfc_dn  = 0
    allocate(this%irflx_sfc_dn  (mesh%ncol                   )); this%irflx_sfc_dn  = 0
    allocate(this%taurefdst     (2*mesh%nlev+4               )); this%taurefdst     = 0
    allocate(this%taurefcld     (2*mesh%nlev+4               )); this%taurefcld     = 0
    allocate(this%taucum        (2*mesh%nlev+3               )); this%taucum        = 0
    allocate(this%taugsurf      (        nspectv,ngauss-1    )); this%taugsurf      = 0
    allocate(this%tausurf       (mesh%ncol                   )); this%tausurf       = 0
    allocate(this%taudst        (mesh%ncol,                 2)); this%taudst        = 0
    allocate(this%taucld        (mesh%ncol,                 2)); this%taucld        = 0
    allocate(this%dtauv         (nlayrad,nspectv,ngauss      )); this%dtauv         = 0
    allocate(this%tauv          (nlevrad,nspectv,ngauss      )); this%tauv          = 0
    allocate(this%taucumv       (2*mesh%nlev+3,nspectv,ngauss)); this%taucumv       = 0
    allocate(this%dtaui         (nlayrad,nspecti,ngauss      )); this%dtaui         = 0
    allocate(this%taui          (nlevrad,nspecti,ngauss      )); this%taui          = 0
    allocate(this%taucumi       (2*mesh%nlev+3,nspecti,ngauss)); this%taucumi       = 0
    allocate(this%detau         (mesh%ncol,nspectv,ngauss    )); this%detau         = 0
    allocate(this%suntot        (2*mesh%nlev+3               )); this%suntot        = 0
    allocate(this%irtot         (2*mesh%nlev+3               )); this%irtot         = 0
    allocate(this%solar_sfc_dn  (nspectv                     )); this%solar_sfc_dn  = 0
    allocate(this%ssun          (mesh%ncol                   )); this%ssun          = 0
    allocate(this%fluxsfc       (mesh%ncol                   )); this%fluxsfc       = 0
    allocate(this%ht_rad        (mesh%ncol,0:mesh%nlev       )); this%ht_rad        = 0
    allocate(this%ht_pbl        (mesh%ncol                   )); this%ht_pbl        = 0
    allocate(this%ht_sfc        (mesh%ncol                   )); this%ht_sfc        = 0
    allocate(this%alsp          (mesh%ncol                   )); this%alsp          = 0
    allocate(this%als           (mesh%ncol                   )); this%als           = 0
    allocate(this%polarcap      (mesh%ncol                   )); this%polarcap      = .false.
    allocate(this%rhouch        (mesh%ncol                   )); this%rhouch        = 0
    allocate(this%shr2          (mesh%ncol,mesh%nlev+1       )); this%shr2          = 0
    allocate(this%ri            (mesh%ncol,mesh%nlev+1       )); this%ri            = 0
    allocate(this%km            (mesh%ncol,mesh%nlev+1       )); this%km            = 0
    allocate(this%kh            (mesh%ncol,mesh%nlev+1       )); this%kh            = 0
    allocate(this%zin           (mesh%ncol,nsoil             )); this%zin           = 0
    allocate(this%rhosoil       (mesh%ncol,nsoil             )); this%rhosoil       = 0
    allocate(this%cpsoil        (mesh%ncol,nsoil             )); this%cpsoil        = 0
    allocate(this%scond         (mesh%ncol,2*nsoil+1         )); this%scond         = 0
    allocate(this%stemp         (mesh%ncol,2*nsoil+1         )); this%stemp         = 0
    allocate(this%h2osub_sfc    (mesh%ncol                   )); this%h2osub_sfc    = 0
    allocate(this%h2oice_sfc    (mesh%ncol                   )); this%h2oice_sfc    = 0
    allocate(this%dmadt         (mesh%ncol                   )); this%dmadt         = 0
    allocate(this%dstflx_wsl    (mesh%ncol                   )); this%dstflx_wsl    = 0
    allocate(this%dstflx_ddl    (mesh%ncol                   )); this%dstflx_ddl    = 0
    allocate(this%tl            (2*mesh%nlev+3               )); this%tl            = 0
    allocate(this%tlev_rad      (2*mesh%nlev+3               )); this%tlev_rad      = 0
    allocate(this%tmid_rad      (2*mesh%nlev+3               )); this%tmid_rad      = 0
    allocate(this%teta          (2*mesh%nlev+3               )); this%teta          = 0
    allocate(this%pl            (2*mesh%nlev+3               )); this%pl            = 0
    allocate(this%plogadj       (2*mesh%nlev+3               )); this%plogadj       = 0
    allocate(this%plev_rad      (2*mesh%nlev+3               )); this%plev_rad      = 0
    allocate(this%pmid_rad      (2*mesh%nlev+3               )); this%pmid_rad      = 0
    allocate(this%aadj          (2*mesh%nlev+3               )); this%aadj          = 0
    allocate(this%badj          (2*mesh%nlev+3               )); this%badj          = 0
    allocate(this%om            (2*mesh%nlev+3               )); this%om            = 0
    allocate(this%qh2o          (2*mesh%nlev+3               )); this%qh2o          = 0

    call this%physics_state_init(mesh)

  end subroutine gomars_v1_state_init

  subroutine gomars_v1_state_clear(this)

    class(gomars_v1_state_type), intent(inout) :: this

    if (allocated(this%ps_old       )) deallocate(this%ps_old       )
    if (allocated(this%ptop_pbl     )) deallocate(this%ptop_pbl     )
    if (allocated(this%tstrat       )) deallocate(this%tstrat       )
    if (allocated(this%co2ice_sfc   )) deallocate(this%co2ice_sfc   )
    if (allocated(this%latheat      )) deallocate(this%latheat      )
    if (allocated(this%atmcond      )) deallocate(this%atmcond      )
    if (allocated(this%tm_sfc       )) deallocate(this%tm_sfc       )
    if (allocated(this%qxvdst       )) deallocate(this%qxvdst       )
    if (allocated(this%qsvdst       )) deallocate(this%qsvdst       )
    if (allocated(this%gvdst        )) deallocate(this%gvdst        )
    if (allocated(this%qxidst       )) deallocate(this%qxidst       )
    if (allocated(this%qsidst       )) deallocate(this%qsidst       )
    if (allocated(this%gidst        )) deallocate(this%gidst        )
    if (allocated(this%qextrefdst   )) deallocate(this%qextrefdst   )
    if (allocated(this%qxvcld       )) deallocate(this%qxvcld       )
    if (allocated(this%qsvcld       )) deallocate(this%qsvcld       )
    if (allocated(this%gvcld        )) deallocate(this%gvcld        )
    if (allocated(this%qxicld       )) deallocate(this%qxicld       )
    if (allocated(this%qsicld       )) deallocate(this%qsicld       )
    if (allocated(this%gicld        )) deallocate(this%gicld        )
    if (allocated(this%qextrefcld   )) deallocate(this%qextrefcld   )
    if (allocated(this%wbarv        )) deallocate(this%wbarv        )
    if (allocated(this%cosbv        )) deallocate(this%cosbv        )
    if (allocated(this%wbari        )) deallocate(this%wbari        )
    if (allocated(this%cosbi        )) deallocate(this%cosbi        )
    if (allocated(this%fluxupv      )) deallocate(this%fluxupv      )
    if (allocated(this%fluxdnv      )) deallocate(this%fluxdnv      )
    if (allocated(this%fmnetv       )) deallocate(this%fmnetv       )
    if (allocated(this%fluxupi      )) deallocate(this%fluxupi      )
    if (allocated(this%fluxdni      )) deallocate(this%fluxdni      )
    if (allocated(this%fmneti       )) deallocate(this%fmneti       )
    if (allocated(this%fuptopv      )) deallocate(this%fuptopv      )
    if (allocated(this%fdntopv      )) deallocate(this%fdntopv      )
    if (allocated(this%fupsfcv      )) deallocate(this%fupsfcv      )
    if (allocated(this%fdnsfcv      )) deallocate(this%fdnsfcv      )
    if (allocated(this%fuptopi      )) deallocate(this%fuptopi      )
    if (allocated(this%fupsfci      )) deallocate(this%fupsfci      )
    if (allocated(this%fdnsfci      )) deallocate(this%fdnsfci      )
    if (allocated(this%taurefcld    )) deallocate(this%taurefcld    )
    if (allocated(this%vsdif_sfc_dn )) deallocate(this%vsdif_sfc_dn )
    if (allocated(this%vsflx_sfc_dn )) deallocate(this%vsflx_sfc_dn )
    if (allocated(this%irflx_sfc_dn )) deallocate(this%irflx_sfc_dn )
    if (allocated(this%taurefdst    )) deallocate(this%taurefdst    )
    if (allocated(this%taurefcld    )) deallocate(this%taurefcld    )
    if (allocated(this%taucum       )) deallocate(this%taucum       )
    if (allocated(this%taugsurf     )) deallocate(this%taugsurf     )
    if (allocated(this%tausurf      )) deallocate(this%tausurf      )
    if (allocated(this%taudst       )) deallocate(this%taudst       )
    if (allocated(this%taucld       )) deallocate(this%taucld       )
    if (allocated(this%dtauv        )) deallocate(this%dtauv        )
    if (allocated(this%tauv         )) deallocate(this%tauv         )
    if (allocated(this%taucumv      )) deallocate(this%taucumv      )
    if (allocated(this%dtaui        )) deallocate(this%dtaui        )
    if (allocated(this%taui         )) deallocate(this%taui         )
    if (allocated(this%taucumi      )) deallocate(this%taucumi      )
    if (allocated(this%detau        )) deallocate(this%detau        )
    if (allocated(this%suntot       )) deallocate(this%suntot       )
    if (allocated(this%irtot        )) deallocate(this%irtot        )
    if (allocated(this%solar_sfc_dn )) deallocate(this%solar_sfc_dn )
    if (allocated(this%ssun         )) deallocate(this%ssun         )
    if (allocated(this%fluxsfc      )) deallocate(this%fluxsfc      )
    if (allocated(this%ht_rad       )) deallocate(this%ht_rad       )
    if (allocated(this%ht_pbl       )) deallocate(this%ht_pbl       )
    if (allocated(this%ht_sfc       )) deallocate(this%ht_sfc       )
    if (allocated(this%alsp         )) deallocate(this%alsp         )
    if (allocated(this%als          )) deallocate(this%als          )
    if (allocated(this%polarcap     )) deallocate(this%polarcap     )
    if (allocated(this%rhouch       )) deallocate(this%rhouch       )
    if (allocated(this%shr2         )) deallocate(this%shr2         )
    if (allocated(this%ri           )) deallocate(this%ri           )
    if (allocated(this%km           )) deallocate(this%km           )
    if (allocated(this%kh           )) deallocate(this%kh           )
    if (allocated(this%zin          )) deallocate(this%zin          )
    if (allocated(this%rhosoil      )) deallocate(this%rhosoil      )
    if (allocated(this%cpsoil       )) deallocate(this%cpsoil       )
    if (allocated(this%scond        )) deallocate(this%scond        )
    if (allocated(this%stemp        )) deallocate(this%stemp        )
    if (allocated(this%h2osub_sfc   )) deallocate(this%h2osub_sfc   )
    if (allocated(this%h2oice_sfc   )) deallocate(this%h2oice_sfc   )
    if (allocated(this%dmadt        )) deallocate(this%dmadt        )
    if (allocated(this%dstflx_wsl   )) deallocate(this%dstflx_wsl   )
    if (allocated(this%dstflx_ddl   )) deallocate(this%dstflx_ddl   )
    if (allocated(this%tl           )) deallocate(this%tl           )
    if (allocated(this%tlev_rad     )) deallocate(this%tlev_rad     )
    if (allocated(this%tmid_rad     )) deallocate(this%tmid_rad     )
    if (allocated(this%teta         )) deallocate(this%teta         )
    if (allocated(this%pl           )) deallocate(this%pl           )
    if (allocated(this%plogadj      )) deallocate(this%plogadj      )
    if (allocated(this%plev_rad     )) deallocate(this%plev_rad     )
    if (allocated(this%pmid_rad     )) deallocate(this%pmid_rad     )
    if (allocated(this%aadj         )) deallocate(this%aadj         )
    if (allocated(this%badj         )) deallocate(this%badj         )
    if (allocated(this%om           )) deallocate(this%om           )
    if (allocated(this%qh2o         )) deallocate(this%qh2o         )

    call this%physics_state_clear()

  end subroutine gomars_v1_state_clear

  subroutine gomars_v1_state_final(this)

    type(gomars_v1_state_type), intent(inout) :: this

    call this%clear()

  end subroutine gomars_v1_state_final

  subroutine gomars_v1_tend_init(this, mesh)

    class(gomars_v1_tend_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%clear()

    call this%physics_tend_init(mesh)

    allocate(this%dudt_pbl(mesh%ncol,2*mesh%nlev+1))

  end subroutine gomars_v1_tend_init

  subroutine gomars_v1_tend_clear(this)

    class(gomars_v1_tend_type), intent(inout) :: this

    if (allocated(this%dudt_pbl)) deallocate(this%dudt_pbl)

  end subroutine gomars_v1_tend_clear

  subroutine gomars_v1_tend_final(this)

    type(gomars_v1_tend_type), intent(inout) :: this

    call this%clear()

  end subroutine gomars_v1_tend_final

end module gomars_v1_types_mod