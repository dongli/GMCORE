! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module physics_types_mod

  use flogger
  use const_mod
  use namelist_mod
  use physics_mesh_mod
  use tracer_types_mod

  implicit none

  private

  public physics_mesh_type
  public physics_state_type
  public physics_tend_type
  public physics_use_wet_tracers

  type, abstract :: physics_state_type
    type(physics_mesh_type), pointer :: mesh => null()
    ! U-wind speed at time level n (m s-1)
    real(r8), allocatable, dimension(:,:  ) :: u_old
    ! U-wind speed (m s-1)
    real(r8), allocatable, dimension(:,:  ) :: u
    ! V-wind speed at time level n (m s-1)
    real(r8), allocatable, dimension(:,:  ) :: v_old
    ! V-wind speed (m s-1)
    real(r8), allocatable, dimension(:,:  ) :: v
    ! Air temperature at time level n (K)
    real(r8), allocatable, dimension(:,:  ) :: t_old
    ! Air temperature (K)
    real(r8), allocatable, dimension(:,:  ) :: t
    ! Air temperature on half levels (K)
    real(r8), allocatable, dimension(:,:  ) :: t_lev
    ! Ground temperature (K)
    real(r8), allocatable, dimension(:    ) :: tg
    ! Lowest model level temperature (K)
    real(r8), pointer    , dimension(:    ) :: t_bot
    ! Potential temperature on full levels at time level n (K)
    real(r8), allocatable, dimension(:,:  ) :: pt_old
    ! Potential temperature on full levels (K)
    real(r8), allocatable, dimension(:,:  ) :: pt
    ! Potential temperature on half levels (K)
    real(r8), allocatable, dimension(:,:  ) :: pt_lev
    ! Tracer wet mixing ratio (kg kg-1)
    real(r8), allocatable, dimension(:,:,:) :: q
    ! Tracer wet mixing ratio at time level n (kg kg-1)
    real(r8), allocatable, dimension(:,:,:) :: q_old
    ! Tracer mass down flux on the surface (kg m-2 s-1)
    real(r8), allocatable, dimension(:,:  ) :: tmflx_sfc_dn
    ! Full pressure (hydrostatic) on full level (Pa)
    real(r8), allocatable, dimension(:,:  ) :: p
    ! Full pressure (hydrostatic) on half level (Pa)
    real(r8), allocatable, dimension(:,:  ) :: p_lev
    ! Exner function of full pressure (hydrostatic) on full level
    real(r8), allocatable, dimension(:,:  ) :: pk
    ! Exner function of full pressure (hydrostatic) on half levels
    real(r8), allocatable, dimension(:,:  ) :: pk_lev
    ! Logarithm of pressure on full levels (Pa)
    real(r8), allocatable, dimension(:,:  ) :: lnp
    ! Logarithm of pressure on half levels (Pa)
    real(r8), allocatable, dimension(:,:  ) :: lnp_lev
    ! Full pressure thickness (Pa)
    real(r8), allocatable, dimension(:,:  ) :: dp
    ! Dry air pressure thickness (Pa)
    real(r8), allocatable, dimension(:,:  ) :: dp_dry
    ! Vertical pressure velocity (Pa s-1)
    real(r8), allocatable, dimension(:,:  ) :: omg
    ! Height on full levels (m)
    real(r8), allocatable, dimension(:,:  ) :: z
    ! Height on half levels (m)
    real(r8), allocatable, dimension(:,:  ) :: z_lev
    ! Height thickness on full levels (m)
    real(r8), allocatable, dimension(:,:  ) :: dz
    ! Height thickness on half levels (m)
    real(r8), allocatable, dimension(:,:  ) :: dz_lev
    ! Air density on full levels (kg m-3)
    real(r8), allocatable, dimension(:,:  ) :: rho
    ! Air density on half levels (kg m-3)
    real(r8), allocatable, dimension(:,:  ) :: rho_lev
    ! Surface pressure (hydrostatic) (Pa)
    real(r8), allocatable, dimension(:    ) :: ps
    ! Surface temperature (K)
    real(r8), allocatable, dimension(:    ) :: ts
    ! Wind speed on lowest model level (m s-1)
    real(r8), allocatable, dimension(:    ) :: wsp_bot
    ! Background roughness length (m)
    real(r8), allocatable, dimension(:    ) :: z0
    ! u* in similarity theory (m s-1)
    real(r8), allocatable, dimension(:    ) :: ustar
    ! Temperature scale (K)
    real(r8), allocatable, dimension(:    ) :: tstar
    ! w* in similarity theory (m s-1)
    real(r8), allocatable, dimension(:    ) :: wstar
    ! Surface drag coefficient for momentum
    real(r8), allocatable, dimension(:    ) :: cdm
    ! Surface drag coefficient for heat
    real(r8), allocatable, dimension(:    ) :: cdh
    ! Wind stress (N m-2)
    real(r8), allocatable, dimension(:    ) :: taux
    real(r8), allocatable, dimension(:    ) :: tauy
    ! Similarity function for momentum
    real(r8), allocatable, dimension(:    ) :: psim
    ! Similarity function for heat
    real(r8), allocatable, dimension(:    ) :: psih
    ! Integrated function or stability function for momentum
    real(r8), allocatable, dimension(:    ) :: fm
    ! Integrated function or stability function for heat
    real(r8), allocatable, dimension(:    ) :: fh
    ! PBL height (m)
    real(r8), allocatable, dimension(:    ) :: pblh
    ! PBL top vertical index
    integer , allocatable, dimension(:    ) :: pblk
    ! Upward heat flux at surface (W m-2)
    real(r8), allocatable, dimension(:    ) :: hflx
    ! Surface albedo
    real(r8), allocatable, dimension(:    ) :: alb
    ! Cosine of solar zenith angle
    real(r8), allocatable, dimension(:    ) :: cosz
    ! Downward solar shortwave flux on the top of atmosphere (W m-2)
    real(r8), allocatable, dimension(:    ) :: fdns_toa
    ! Direct downward solar shortwave flux on the surface (W m-2)
    real(r8), allocatable, dimension(:    ) :: fdns_dir
    ! Diffusive downward solar shortwave flux on the surface (W m-2)
    real(r8), allocatable, dimension(:    ) :: fdns_dif
    ! Downward solar shortwave flux on the surface (W m-2)
    real(r8), allocatable, dimension(:    ) :: fdns
    ! Downward longwave flux on the surface (W m-2)
    real(r8), allocatable, dimension(:    ) :: fdnl
    ! Topography (m)
    real(r8), allocatable, dimension(:    ) :: zs
    ! Total precipitation (kg m-2)
    real(r8), allocatable, dimension(:    ) :: prect
    ! Convective precipitation (kg m-2)
    real(r8), allocatable, dimension(:    ) :: precc
    ! Large-scale precipitation (kg m-2)
    real(r8), allocatable, dimension(:    ) :: precl
  contains
    procedure physics_state_init
    procedure physics_state_clear
  end type physics_state_type

  type, abstract :: physics_tend_type
    type(physics_mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:  ) :: dudt
    real(r8), allocatable, dimension(:,:  ) :: dvdt
    real(r8), allocatable, dimension(:,:  ) :: dtdt
    real(r8), allocatable, dimension(:,:  ) :: dptdt
    real(r8), allocatable, dimension(:,:,:) :: dqdt
    logical :: updated_u  = .false.
    logical :: updated_v  = .false.
    logical :: updated_t  = .false.
    logical :: updated_pt = .false.
    logical, allocatable :: updated_q(:)
  contains
    procedure physics_tend_init
    procedure physics_tend_clear
    procedure physics_tend_reset
  end type physics_tend_type

  logical, allocatable :: physics_use_wet_tracers(:)

contains

  subroutine physics_state_init(this, mesh)

    class(physics_state_type), intent(inout), target :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%physics_state_clear()

    this%mesh => mesh

    if (ntracers < 1) call log_error('ntracers is less than 1!', __FILE__, __LINE__)

    allocate(this%u_old       (mesh%ncol,mesh%nlev         )); this%u_old        = 0
    allocate(this%u           (mesh%ncol,mesh%nlev         )); this%u            = 0
    allocate(this%v_old       (mesh%ncol,mesh%nlev         )); this%v_old        = 0
    allocate(this%v           (mesh%ncol,mesh%nlev         )); this%v            = 0
    allocate(this%t_old       (mesh%ncol,mesh%nlev         )); this%t_old        = 0
    allocate(this%t           (mesh%ncol,mesh%nlev         )); this%t            = 0
    allocate(this%t_lev       (mesh%ncol,mesh%nlev+1       )); this%t_lev        = 0
    allocate(this%tg          (mesh%ncol                   )); this%tg           = 0
    allocate(this%pt_old      (mesh%ncol,mesh%nlev         )); this%pt_old       = 0
    allocate(this%pt          (mesh%ncol,mesh%nlev         )); this%pt           = 0
    allocate(this%pt_lev      (mesh%ncol,mesh%nlev+1       )); this%pt_lev       = 0
    allocate(this%q_old       (mesh%ncol,mesh%nlev,ntracers)); this%q_old        = 0
    allocate(this%q           (mesh%ncol,mesh%nlev,ntracers)); this%q            = 0
    allocate(this%tmflx_sfc_dn(mesh%ncol,          ntracers)); this%tmflx_sfc_dn = 0
    allocate(this%p           (mesh%ncol,mesh%nlev         )); this%p            = 0
    allocate(this%p_lev       (mesh%ncol,mesh%nlev+1       )); this%p_lev        = 0
    allocate(this%pk          (mesh%ncol,mesh%nlev         )); this%pk           = 0
    allocate(this%pk_lev      (mesh%ncol,mesh%nlev+1       )); this%pk_lev       = 0
    allocate(this%lnp         (mesh%ncol,mesh%nlev         )); this%lnp          = 0
    allocate(this%lnp_lev     (mesh%ncol,mesh%nlev+1       )); this%lnp_lev      = 0
    allocate(this%dp          (mesh%ncol,mesh%nlev         )); this%dp           = 0
    allocate(this%dp_dry      (mesh%ncol,mesh%nlev         )); this%dp_dry       = 0
    allocate(this%omg         (mesh%ncol,mesh%nlev         )); this%omg          = 0
    allocate(this%z           (mesh%ncol,mesh%nlev         )); this%z            = 0
    allocate(this%z_lev       (mesh%ncol,mesh%nlev+1       )); this%z_lev        = 0
    allocate(this%dz          (mesh%ncol,mesh%nlev         )); this%dz           = 0
    allocate(this%dz_lev      (mesh%ncol,mesh%nlev+1       )); this%dz_lev       = 0
    allocate(this%rho         (mesh%ncol,mesh%nlev         )); this%rho          = 0
    allocate(this%rho_lev     (mesh%ncol,mesh%nlev+1       )); this%rho_lev      = 0
    allocate(this%ps          (mesh%ncol                   )); this%ps           = 0
    allocate(this%ts          (mesh%ncol                   )); this%ts           = 0
    allocate(this%wsp_bot     (mesh%ncol                   )); this%wsp_bot      = 0
    allocate(this%z0          (mesh%ncol                   )); this%z0           = 0
    allocate(this%ustar       (mesh%ncol                   )); this%ustar        = 0
    allocate(this%tstar       (mesh%ncol                   )); this%tstar        = 0
    allocate(this%wstar       (mesh%ncol                   )); this%wstar        = 0
    allocate(this%cdm         (mesh%ncol                   )); this%cdm          = 0
    allocate(this%cdh         (mesh%ncol                   )); this%cdh          = 0
    allocate(this%taux        (mesh%ncol                   )); this%taux         = 0
    allocate(this%tauy        (mesh%ncol                   )); this%tauy         = 0
    allocate(this%psim        (mesh%ncol                   )); this%psim         = 0
    allocate(this%psih        (mesh%ncol                   )); this%psih         = 0
    allocate(this%fm          (mesh%ncol                   )); this%fm           = 0
    allocate(this%fh          (mesh%ncol                   )); this%fh           = 0
    allocate(this%pblh        (mesh%ncol                   )); this%pblh         = 0
    allocate(this%pblk        (mesh%ncol                   )); this%pblk         = 0
    allocate(this%hflx        (mesh%ncol                   )); this%hflx         = 0
    allocate(this%alb         (mesh%ncol                   )); this%alb          = 0
    allocate(this%cosz        (mesh%ncol                   )); this%cosz         = 0
    allocate(this%fdns_toa    (mesh%ncol                   )); this%fdns_toa     = 0
    allocate(this%fdns_dir    (mesh%ncol                   )); this%fdns_dir     = 0
    allocate(this%fdns_dif    (mesh%ncol                   )); this%fdns_dif     = 0
    allocate(this%fdns        (mesh%ncol                   )); this%fdns         = 0
    allocate(this%fdnl        (mesh%ncol                   )); this%fdnl         = 0
    allocate(this%zs          (mesh%ncol                   )); this%zs           = 0
    allocate(this%prect       (mesh%ncol                   )); this%prect        = 0
    allocate(this%precc       (mesh%ncol                   )); this%precc        = 0
    allocate(this%precl       (mesh%ncol                   )); this%precl        = 0

    this%t_bot => this%t(:,mesh%nlev)

  end subroutine physics_state_init

  subroutine physics_state_clear(this)

    class(physics_state_type), intent(inout) :: this

    this%mesh => null()

    if (allocated(this%u_old        )) deallocate(this%u_old        )
    if (allocated(this%u            )) deallocate(this%u            )
    if (allocated(this%v_old        )) deallocate(this%v_old        )
    if (allocated(this%v            )) deallocate(this%v            )
    if (allocated(this%t_old        )) deallocate(this%t_old        )
    if (allocated(this%t            )) deallocate(this%t            )
    if (allocated(this%t_lev        )) deallocate(this%t_lev        )
    if (allocated(this%tg           )) deallocate(this%tg           )
    if (allocated(this%pt_old       )) deallocate(this%pt_old       )
    if (allocated(this%pt           )) deallocate(this%pt           )
    if (allocated(this%pt_lev       )) deallocate(this%pt_lev       )
    if (allocated(this%q_old        )) deallocate(this%q_old        )
    if (allocated(this%q            )) deallocate(this%q            )
    if (allocated(this%tmflx_sfc_dn )) deallocate(this%tmflx_sfc_dn )
    if (allocated(this%p            )) deallocate(this%p            )
    if (allocated(this%p_lev        )) deallocate(this%p_lev        )
    if (allocated(this%pk           )) deallocate(this%pk           )
    if (allocated(this%pk_lev       )) deallocate(this%pk_lev       )
    if (allocated(this%lnp          )) deallocate(this%lnp          )
    if (allocated(this%lnp_lev      )) deallocate(this%lnp_lev      )
    if (allocated(this%dp           )) deallocate(this%dp           )
    if (allocated(this%dp_dry       )) deallocate(this%dp_dry       )
    if (allocated(this%omg          )) deallocate(this%omg          )
    if (allocated(this%z            )) deallocate(this%z            )
    if (allocated(this%z_lev        )) deallocate(this%z_lev        )
    if (allocated(this%dz           )) deallocate(this%dz           )
    if (allocated(this%rho          )) deallocate(this%rho          )
    if (allocated(this%rho_lev      )) deallocate(this%rho_lev      )
    if (allocated(this%ps           )) deallocate(this%ps           )
    if (allocated(this%ts           )) deallocate(this%ts           )
    if (allocated(this%wsp_bot      )) deallocate(this%wsp_bot      )
    if (allocated(this%z0           )) deallocate(this%z0           )
    if (allocated(this%ustar        )) deallocate(this%ustar        )
    if (allocated(this%tstar        )) deallocate(this%tstar        )
    if (allocated(this%wstar        )) deallocate(this%wstar        )
    if (allocated(this%cdm          )) deallocate(this%cdm          )
    if (allocated(this%cdh          )) deallocate(this%cdh          )
    if (allocated(this%taux         )) deallocate(this%taux         )
    if (allocated(this%tauy         )) deallocate(this%tauy         )
    if (allocated(this%psim         )) deallocate(this%psim         )
    if (allocated(this%psih         )) deallocate(this%psih         )
    if (allocated(this%fm           )) deallocate(this%fm           )
    if (allocated(this%fh           )) deallocate(this%fh           )
    if (allocated(this%pblh         )) deallocate(this%pblh         )
    if (allocated(this%pblk         )) deallocate(this%pblk         )
    if (allocated(this%hflx         )) deallocate(this%hflx         )
    if (allocated(this%alb          )) deallocate(this%alb          )
    if (allocated(this%cosz         )) deallocate(this%cosz         )
    if (allocated(this%fdns_toa     )) deallocate(this%fdns_toa     )
    if (allocated(this%fdns_dir     )) deallocate(this%fdns_dir     )
    if (allocated(this%fdns_dif     )) deallocate(this%fdns_dif     )
    if (allocated(this%fdns         )) deallocate(this%fdns         )
    if (allocated(this%fdnl         )) deallocate(this%fdnl         )
    if (allocated(this%zs           )) deallocate(this%zs           )
    if (allocated(this%prect        )) deallocate(this%prect        )
    if (allocated(this%precc        )) deallocate(this%precc        )
    if (allocated(this%precl        )) deallocate(this%precl        )

  end subroutine physics_state_clear

  subroutine physics_tend_init(this, mesh)

    class(physics_tend_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%physics_tend_clear()

    this%mesh => mesh

    allocate(this%dudt (mesh%ncol,mesh%nlev         ))
    allocate(this%dvdt (mesh%ncol,mesh%nlev         ))
    allocate(this%dtdt (mesh%ncol,mesh%nlev         ))
    allocate(this%dptdt(mesh%ncol,mesh%nlev         ))
    allocate(this%dqdt (mesh%ncol,mesh%nlev,ntracers))
    allocate(this%updated_q(ntracers))

  end subroutine physics_tend_init

  subroutine physics_tend_clear(this)

    class(physics_tend_type), intent(inout) :: this

    this%mesh => null()

    if (allocated(this%dudt     )) deallocate(this%dudt     )
    if (allocated(this%dvdt     )) deallocate(this%dvdt     )
    if (allocated(this%dtdt     )) deallocate(this%dtdt     )
    if (allocated(this%dptdt    )) deallocate(this%dptdt    )
    if (allocated(this%dqdt     )) deallocate(this%dqdt     )
    if (allocated(this%updated_q)) deallocate(this%updated_q)

  end subroutine physics_tend_clear

  subroutine physics_tend_reset(this)

    class(physics_tend_type), intent(inout) :: this

    this%dudt  = 0; this%updated_u  = .false.
    this%dvdt  = 0; this%updated_v  = .false.
    this%dtdt  = 0; this%updated_t  = .false.
    this%dptdt = 0; this%updated_pt = .false.
    this%dqdt  = 0; this%updated_q  = .false.

  end subroutine physics_tend_reset

end module physics_types_mod
