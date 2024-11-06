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

module gomars_v2_types_mod

  use fiona
  use process_mod
  use physics_types_mod
  use gomars_v2_const_mod
  use gomars_v2_namelist_mod
  use gomars_v2_spectra_mod
  use gomars_v2_rad_kcoef_mod
  use gomars_v2_tracers_mod

  implicit none

  private

  public gomars_v2_state_type
  public gomars_v2_tend_type

  type, extends(physics_state_type) :: gomars_v2_state_type
    ! Surface latent heat flux (W m-2)
    real(r8), allocatable, dimension(    :    ) :: lhflx
    ! Atmosphere CO2 condensation (???)
    real(r8), allocatable, dimension(    :,:  ) :: atmcond
    ! Dust particle median radius
    real(r8), allocatable, dimension(    :,:  ) :: ro_dst
    ! Cloud particle median radius
    real(r8), allocatable, dimension(    :,:  ) :: ro_cld
    ! delta-Eddington optical thickness on the surface (1)
    real(r8), allocatable, dimension(:,:,:    ) :: detau
    real(r8), allocatable, dimension(:,:,:,:  ) :: tau_gas_vs
    real(r8), allocatable, dimension(:,:,:,:  ) :: tau_dst_vs
    real(r8), allocatable, dimension(:,:,:,:  ) :: tau_cld_vs
    ! Top or stratosphere temperature (K)
    real(r8), allocatable, dimension(      :  ) :: t_top
    real(r8), allocatable, dimension(      :  ) :: co2ice
    ! Boundary values of tracers (kg m-2???)
    real(r8), allocatable, dimension(    :,  :) :: q_sfc
    ! Intermediate variable ρ u* cp ch (???)
    real(r8), allocatable, dimension(    :    ) :: rhouch
    ! Eddy mixing coefficients in vertical (???)
    real(r8), allocatable, dimension(    :,:  ) :: eddy_km
    real(r8), allocatable, dimension(    :,:  ) :: eddy_kh
    ! Saved square of wind shear (s-2)
    real(r8), allocatable, dimension(    :,:  ) :: shr2
    ! Heat exchange between atmosphere and surface (W m-2)
    real(r8), allocatable, dimension(    :    ) :: heat_sfc
    ! Soil conductivity (W m-1 K-1)
    real(r8), allocatable, dimension(    :,:  ) :: soil_cond
    real(r8), allocatable, dimension(    :,:  ) :: soil_cond_lev
    ! Soil density (kg m-3)
    real(r8), allocatable, dimension(    :,:  ) :: soil_rho
    ! Soil specific heat capacity (J kg-1 K-1)
    real(r8), allocatable, dimension(    :,:  ) :: soil_cp
    ! Soil temperature (K)
    real(r8), allocatable, dimension(    :,:  ) :: soil_t
    ! The fact that Mars is now a desert planet without oceans and lakes means
    ! that the thermal inertia of the surface is small.
    ! Surface thermal interia
    real(r8), allocatable, dimension(:  ) :: tin
    real(r8), allocatable, dimension(:,:) :: soil_tin
    ! Ground ice indicator for initializing soil model
    real(r8), allocatable, dimension(:  ) :: gnd_ice
  contains
    procedure :: init  => gomars_v2_state_init
    procedure :: clear => gomars_v2_state_clear
    final gomars_v2_state_final
  end type gomars_v2_state_type

  type, extends(physics_tend_type) :: gomars_v2_tend_type
    real(r8), allocatable :: dpsdt(:)
  contains
    procedure :: init  => gomars_v2_tend_init
    procedure :: clear => gomars_v2_tend_clear
    final gomars_v2_tend_final
  end type gomars_v2_tend_type

contains

  subroutine gomars_v2_state_init(this, mesh)

    class(gomars_v2_state_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%clear()

    allocate(this%lhflx        (                 mesh%ncol            ))
    allocate(this%atmcond      (                 mesh%ncol,mesh%nlev  ))
    allocate(this%ro_dst       (                 mesh%ncol,mesh%nlev  ))
    allocate(this%ro_cld       (                 mesh%ncol,mesh%nlev  ))
    allocate(this%detau        (spec_vs%n,ngauss,mesh%ncol            ))
    allocate(this%tau_gas_vs   (spec_vs%n,ngauss,mesh%ncol,mesh%nlev  ))
    allocate(this%tau_dst_vs   (spec_vs%n,ngauss,mesh%ncol,mesh%nlev  ))
    allocate(this%tau_cld_vs   (spec_vs%n,ngauss,mesh%ncol,mesh%nlev  ))
    allocate(this%t_top        (                 mesh%ncol            ))
    allocate(this%co2ice       (                 mesh%ncol            ))
    allocate(this%q_sfc        (                 mesh%ncol,           ntracers))
    allocate(this%rhouch       (                 mesh%ncol            ))
    allocate(this%eddy_km      (                 mesh%ncol,mesh%nlev+1))
    allocate(this%eddy_kh      (                 mesh%ncol,mesh%nlev+1))
    allocate(this%shr2         (                 mesh%ncol,mesh%nlev+1))
    allocate(this%heat_sfc     (                 mesh%ncol            ))
    allocate(this%soil_cond    (                 mesh%ncol,nlev_soil  ))
    allocate(this%soil_cond_lev(                 mesh%ncol,nlev_soil+1))
    allocate(this%soil_rho     (                 mesh%ncol,nlev_soil  ))
    allocate(this%soil_cp      (                 mesh%ncol,nlev_soil  ))
    allocate(this%soil_t       (                 mesh%ncol,nlev_soil  ))
    allocate(this%tin          (                 mesh%ncol            ))
    allocate(this%soil_tin     (                 mesh%ncol,nlev_soil  ))
    allocate(this%gnd_ice      (                 mesh%ncol            ))

    call this%physics_state_init(mesh)

  end subroutine gomars_v2_state_init

  subroutine gomars_v2_state_clear(this)

    class(gomars_v2_state_type), intent(inout) :: this

    if (allocated(this%lhflx        )) deallocate(this%lhflx        )
    if (allocated(this%atmcond      )) deallocate(this%atmcond      )
    if (allocated(this%ro_dst       )) deallocate(this%ro_dst       )
    if (allocated(this%ro_cld       )) deallocate(this%ro_cld       )
    if (allocated(this%detau        )) deallocate(this%detau        )
    if (allocated(this%tau_gas_vs   )) deallocate(this%tau_gas_vs   )
    if (allocated(this%tau_dst_vs   )) deallocate(this%tau_dst_vs   )
    if (allocated(this%tau_cld_vs   )) deallocate(this%tau_cld_vs   )
    if (allocated(this%t_top        )) deallocate(this%t_top        )
    if (allocated(this%co2ice       )) deallocate(this%co2ice       )
    if (allocated(this%q_sfc        )) deallocate(this%q_sfc        )
    if (allocated(this%rhouch       )) deallocate(this%rhouch       )
    if (allocated(this%eddy_km      )) deallocate(this%eddy_km      )
    if (allocated(this%eddy_kh      )) deallocate(this%eddy_kh      )
    if (allocated(this%shr2         )) deallocate(this%shr2         )
    if (allocated(this%heat_sfc     )) deallocate(this%heat_sfc     )
    if (allocated(this%soil_cond    )) deallocate(this%soil_cond    )
    if (allocated(this%soil_cond_lev)) deallocate(this%soil_cond_lev)
    if (allocated(this%soil_rho     )) deallocate(this%soil_rho     )
    if (allocated(this%soil_cp      )) deallocate(this%soil_cp      )
    if (allocated(this%soil_t       )) deallocate(this%soil_t       )
    if (allocated(this%tin          )) deallocate(this%tin          )
    if (allocated(this%soil_tin     )) deallocate(this%soil_tin     )
    if (allocated(this%gnd_ice      )) deallocate(this%gnd_ice      )

    call this%physics_state_clear()

  end subroutine gomars_v2_state_clear

  subroutine gomars_v2_state_final(this)

    type(gomars_v2_state_type), intent(inout) :: this

    call this%clear()

  end subroutine gomars_v2_state_final

  subroutine gomars_v2_tend_init(this, mesh)

    class(gomars_v2_tend_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%clear()

    allocate(this%dpsdt(mesh%ncol))

    call this%physics_tend_init(mesh)

  end subroutine gomars_v2_tend_init

  subroutine gomars_v2_tend_clear(this)

    class(gomars_v2_tend_type), intent(inout) :: this

    if (allocated(this%dpsdt)) deallocate(this%dpsdt)

    call this%physics_tend_clear()

  end subroutine gomars_v2_tend_clear

  subroutine gomars_v2_tend_final(this)

    type(gomars_v2_tend_type), intent(inout) :: this

    call this%clear()

  end subroutine gomars_v2_tend_final

end module gomars_v2_types_mod
