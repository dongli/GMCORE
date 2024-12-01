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
    ! Stratospheric temperature (K)
    real(r8), allocatable, dimension(:      ) :: tstrat
    ! Surface temperature (K)
    real(r8), allocatable, dimension(:      ) :: gt
    ! Surface CO2 ice (?)
    real(r8), allocatable, dimension(:      ) :: co2ice
    !
    real(r8), allocatable, dimension(:      ) :: latheat
    !
    real(r8), allocatable, dimension(:,:    ) :: atmcond
    ! Surface boundary condition for tracers
    real(r8), allocatable, dimension(:,    :) :: qcond
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
    real(r8), allocatable, dimension(  :    ) :: fluxupv
    !
    real(r8), allocatable, dimension(  :    ) :: fluxdnv
    !
    real(r8), allocatable, dimension(  :    ) :: fmnetv
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
    ! Downward diffuse visible (solar) flux at the surface
    real(r8), allocatable, dimension(:      ) :: dndiffv
    ! Downward visible flux at the surface
    real(r8), allocatable, dimension(:      ) :: dnvflux
    ! Downward IR flux at the surface
    real(r8), allocatable, dimension(:      ) :: dnirflux
    !
    real(r8), allocatable, dimension(:,    :) :: srfdnflx
    !
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
    ! Delta-Eddington optical depth (???)
    real(r8), allocatable, dimension(:,  :,:) :: detau
    !
    real(r8), allocatable, dimension(  :    ) :: suntot
    ! Solar flux at the current Mars distance
    real(r8), allocatable, dimension(    :  ) :: solar
    ! Total absorption of solar energy by the atmosphere (?)
    real(r8), allocatable, dimension(:      ) :: ssun
    ! Exchange of heat between surface and air (W m-2)
    real(r8), allocatable, dimension(:      ) :: fa
    ! Surface albedo
    real(r8), allocatable, dimension(:      ) :: alsp
    !
    real(r8), allocatable, dimension(:      ) :: surfalb
    ! North pole cap flag
    logical , allocatable, dimension(:      ) :: npcflag
    !
    real(r8), allocatable, dimension(:      ) :: rhouch
    ! Thermal inertia (J m-2 K-1 s-1/2)
    real(r8), allocatable, dimension(:,:    ) :: zin
    real(r8), allocatable, dimension(:,:    ) :: rhosoil
    real(r8), allocatable, dimension(:,:    ) :: cpsoil
    !
    real(r8), allocatable, dimension(:,:    ) :: scond
    !
    real(r8), allocatable, dimension(:,:    ) :: stemp
    !
    real(r8), allocatable, dimension(:      ) :: subflux
    !
    real(r8), allocatable, dimension(:      ) :: gndice
    !
    real(r8), allocatable, dimension(:      ) :: dmadt
    ! Square of Mars distance from sun
    real(r8) :: rsdist
    ! Temperature at all levels with two extra levels at the top (K)
    real(r8), allocatable, dimension(    :  ) :: tl
    ! Potential temperature at all levels with two extra levels at the top (K)
    real(r8), allocatable, dimension(    :  ) :: teta
    ! Pressure at all levels with two extra levels at the top (Pa)
    real(r8), allocatable, dimension(    :  ) :: pl
    real(r8), allocatable, dimension(    :  ) :: plogadj
    ! Zonal wind at all levels with two extra levels at the top (m s-1)
    real(r8), allocatable, dimension(    :  ) :: upi
    ! Meridional wind at all levels with two extra levels at the top (m s-1)
    real(r8), allocatable, dimension(    :  ) :: vpi
    ! Tracer mixing ratio at all levels with two extra levels at the top (kg kg-1)
    real(r8), allocatable, dimension(    :,:) :: qpi
    real(r8), allocatable, dimension(    :  ) :: qpig
    real(r8), allocatable, dimension(    :  ) :: aadj
    real(r8), allocatable, dimension(    :  ) :: badj
    real(r8), allocatable, dimension(    :  ) :: om
    !
    real(r8), allocatable, dimension(    :  ) :: qh2o
  contains
    procedure :: init  => gomars_v1_state_init
    procedure :: clear => gomars_v1_state_clear
    final gomars_v1_state_final
  end type gomars_v1_state_type

  type, extends(physics_tend_type) :: gomars_v1_tend_type
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

    allocate(this%tstrat    (mesh%ncol                   )); this%tstrat     = 0
    allocate(this%gt        (mesh%ncol                   )); this%gt         = 0
    allocate(this%co2ice    (mesh%ncol                   )); this%co2ice     = 0
    allocate(this%latheat   (mesh%ncol                   )); this%latheat    = 0
    allocate(this%atmcond   (mesh%ncol,mesh%nlev         )); this%atmcond    = 0
    allocate(this%qcond     (mesh%ncol,ntracers          )); this%qcond      = 0

    allocate(this%qxvdst    (2*mesh%nlev+4,nspectv       )); this%qxvdst     = 0
    allocate(this%qsvdst    (2*mesh%nlev+4,nspectv       )); this%qsvdst     = 0
    allocate(this%gvdst     (2*mesh%nlev+4,nspectv       )); this%gvdst      = 0
    allocate(this%qxidst    (2*mesh%nlev+4,nspecti       )); this%qxidst     = 0
    allocate(this%qsidst    (2*mesh%nlev+4,nspecti       )); this%qsidst     = 0
    allocate(this%gidst     (2*mesh%nlev+4,nspecti       )); this%gidst      = 0
    allocate(this%qextrefdst(2*mesh%nlev+4               )); this%qextrefdst = 0
    allocate(this%qxvcld    (2*mesh%nlev+4,nspectv       )); this%qxvcld     = 0
    allocate(this%qsvcld    (2*mesh%nlev+4,nspectv       )); this%qsvcld     = 0
    allocate(this%gvcld     (2*mesh%nlev+4,nspectv       )); this%gvcld      = 0
    allocate(this%qxicld    (2*mesh%nlev+4,nspecti       )); this%qxicld     = 0
    allocate(this%qsicld    (2*mesh%nlev+4,nspecti       )); this%qsicld     = 0
    allocate(this%gicld     (2*mesh%nlev+4,nspecti       )); this%gicld      = 0
    allocate(this%qextrefcld(2*mesh%nlev+4               )); this%qextrefcld = 0
    allocate(this%wbarv     (nlayrad,nspectv,ngauss      )); this%wbarv      = 0
    allocate(this%cosbv     (nlayrad,nspectv,ngauss      )); this%cosbv      = 0
    allocate(this%fluxupv   (nlayrad                     )); this%fluxupv    = 0
    allocate(this%fluxdnv   (nlayrad                     )); this%fluxdnv    = 0
    allocate(this%fmnetv    (nlayrad                     )); this%fmnetv     = 0
    allocate(this%fuptopv   (mesh%ncol                   )); this%fuptopv    = 0
    allocate(this%fdntopv   (mesh%ncol                   )); this%fdntopv    = 0
    allocate(this%fupsfcv   (mesh%ncol                   )); this%fupsfcv    = 0
    allocate(this%fdnsfcv   (mesh%ncol                   )); this%fdnsfcv    = 0
    allocate(this%fuptopi   (mesh%ncol                   )); this%fuptopi    = 0
    allocate(this%fupsfci   (mesh%ncol                   )); this%fupsfci    = 0
    allocate(this%fdnsfci   (mesh%ncol                   )); this%fdnsfci    = 0
    allocate(this%dndiffv   (mesh%ncol                   )); this%dndiffv    = 0
    allocate(this%dnvflux   (mesh%ncol                   )); this%dnvflux    = 0
    allocate(this%dnirflux  (mesh%ncol                   )); this%dnirflux   = 0
    allocate(this%srfdnflx  (mesh%ncol,          ntracers)); this%srfdnflx   = 0
    allocate(this%taurefdst (2*mesh%nlev+4               )); this%taurefdst  = 0
    allocate(this%taurefcld (2*mesh%nlev+4               )); this%taurefcld  = 0
    allocate(this%taucum    (2*mesh%nlev+3               )); this%taucum     = 0
    allocate(this%taugsurf  (        nspectv,ngauss-1    )); this%taugsurf   = 0
    allocate(this%tausurf   (mesh%ncol                   )); this%tausurf    = 0
    allocate(this%taudst    (mesh%ncol,                 2)); this%taudst     = 0
    allocate(this%taucld    (mesh%ncol,                 2)); this%taucld     = 0
    allocate(this%dtauv     (nlayrad,nspectv,ngauss      )); this%dtauv      = 0
    allocate(this%tauv      (nlevrad,nspectv,ngauss      )); this%tauv       = 0
    allocate(this%taucumv   (2*mesh%nlev+3,nspectv,ngauss)); this%taucumv    = 0
    allocate(this%detau     (mesh%ncol,nspectv,ngauss    )); this%detau      = 0
    allocate(this%suntot    (2*mesh%nlev+3               )); this%suntot     = 0
    allocate(this%solar     (nspectv                     )); this%solar      = 0
    allocate(this%ssun      (mesh%ncol                   )); this%ssun       = 0
    allocate(this%fa        (mesh%ncol                   )); this%fa         = 0
    allocate(this%alsp      (mesh%ncol                   )); this%alsp       = 0
    allocate(this%surfalb   (mesh%ncol                   )); this%surfalb    = 0
    allocate(this%npcflag   (mesh%ncol                   )); this%npcflag    = .false.
    allocate(this%rhouch    (mesh%ncol                   )); this%rhouch     = 0
    allocate(this%zin       (mesh%ncol,nl                )); this%zin        = 0
    allocate(this%rhosoil   (mesh%ncol,nl                )); this%rhosoil    = 0
    allocate(this%cpsoil    (mesh%ncol,nl                )); this%cpsoil     = 0
    allocate(this%scond     (mesh%ncol,2*nl+1            )); this%scond      = 0
    allocate(this%stemp     (mesh%ncol,2*nl+1            )); this%stemp      = 0
    allocate(this%subflux   (mesh%ncol                   )); this%subflux    = 0
    allocate(this%gndice    (mesh%ncol                   )); this%gndice     = 0
    allocate(this%dmadt     (mesh%ncol                   )); this%dmadt      = 0
    allocate(this%tl        (2*mesh%nlev+3               )); this%tl         = 0
    allocate(this%teta      (2*mesh%nlev+3               )); this%teta       = 0
    allocate(this%pl        (2*mesh%nlev+3               )); this%pl         = 0
    allocate(this%plogadj   (2*mesh%nlev+3               )); this%plogadj    = 0
    allocate(this%upi       (2*mesh%nlev+3               )); this%upi        = 0
    allocate(this%vpi       (2*mesh%nlev+3               )); this%vpi        = 0
    allocate(this%qpi       (2*mesh%nlev+3,ntracers      )); this%qpi        = 0
    allocate(this%qpig      (              ntracers      )); this%qpig       = 0
    allocate(this%aadj      (2*mesh%nlev+3               )); this%aadj       = 0
    allocate(this%badj      (2*mesh%nlev+3               )); this%badj       = 0
    allocate(this%om        (2*mesh%nlev+3               )); this%om         = 0
    allocate(this%qh2o      (2*mesh%nlev+3               )); this%qh2o       = 0

    call this%physics_state_init(mesh)

  end subroutine gomars_v1_state_init

  subroutine gomars_v1_state_clear(this)

    class(gomars_v1_state_type), intent(inout) :: this

    if (allocated(this%tstrat    )) deallocate(this%tstrat    )
    if (allocated(this%gt        )) deallocate(this%gt        )
    if (allocated(this%co2ice    )) deallocate(this%co2ice    )
    if (allocated(this%latheat   )) deallocate(this%latheat   )
    if (allocated(this%atmcond   )) deallocate(this%atmcond   )
    if (allocated(this%qcond     )) deallocate(this%qcond     )
    if (allocated(this%qxvdst    )) deallocate(this%qxvdst    )
    if (allocated(this%qsvdst    )) deallocate(this%qsvdst    )
    if (allocated(this%gvdst     )) deallocate(this%gvdst     )
    if (allocated(this%qxidst    )) deallocate(this%qxidst    )
    if (allocated(this%qsidst    )) deallocate(this%qsidst    )
    if (allocated(this%gidst     )) deallocate(this%gidst     )
    if (allocated(this%qextrefdst)) deallocate(this%qextrefdst)
    if (allocated(this%qxvcld    )) deallocate(this%qxvcld    )
    if (allocated(this%qsvcld    )) deallocate(this%qsvcld    )
    if (allocated(this%gvcld     )) deallocate(this%gvcld     )
    if (allocated(this%qxicld    )) deallocate(this%qxicld    )
    if (allocated(this%qsicld    )) deallocate(this%qsicld    )
    if (allocated(this%gicld     )) deallocate(this%gicld     )
    if (allocated(this%qextrefcld)) deallocate(this%qextrefcld)
    if (allocated(this%wbarv     )) deallocate(this%wbarv     )
    if (allocated(this%cosbv     )) deallocate(this%cosbv     )
    if (allocated(this%fluxupv   )) deallocate(this%fluxupv   )
    if (allocated(this%fluxdnv   )) deallocate(this%fluxdnv   )
    if (allocated(this%fmnetv    )) deallocate(this%fmnetv    )
    if (allocated(this%fuptopv   )) deallocate(this%fuptopv   )
    if (allocated(this%fdntopv   )) deallocate(this%fdntopv   )
    if (allocated(this%fupsfcv   )) deallocate(this%fupsfcv   )
    if (allocated(this%fdnsfcv   )) deallocate(this%fdnsfcv   )
    if (allocated(this%fuptopi   )) deallocate(this%fuptopi   )
    if (allocated(this%fupsfci   )) deallocate(this%fupsfci   )
    if (allocated(this%fdnsfci   )) deallocate(this%fdnsfci   )
    if (allocated(this%taurefcld )) deallocate(this%taurefcld )
    if (allocated(this%dndiffv   )) deallocate(this%dndiffv   )
    if (allocated(this%dnvflux   )) deallocate(this%dnvflux   )
    if (allocated(this%dnirflux  )) deallocate(this%dnirflux  )
    if (allocated(this%srfdnflx  )) deallocate(this%srfdnflx  )
    if (allocated(this%taurefdst )) deallocate(this%taurefdst )
    if (allocated(this%taurefcld )) deallocate(this%taurefcld )
    if (allocated(this%taucum    )) deallocate(this%taucum    )
    if (allocated(this%taugsurf  )) deallocate(this%taugsurf  )
    if (allocated(this%tausurf   )) deallocate(this%tausurf   )
    if (allocated(this%taudst    )) deallocate(this%taudst    )
    if (allocated(this%taucld    )) deallocate(this%taucld    )
    if (allocated(this%dtauv     )) deallocate(this%dtauv     )
    if (allocated(this%tauv      )) deallocate(this%tauv      )
    if (allocated(this%taucumv   )) deallocate(this%taucumv   )
    if (allocated(this%detau     )) deallocate(this%detau     )
    if (allocated(this%suntot    )) deallocate(this%suntot    )
    if (allocated(this%solar     )) deallocate(this%solar     )
    if (allocated(this%ssun      )) deallocate(this%ssun      )
    if (allocated(this%fa        )) deallocate(this%fa        )
    if (allocated(this%alsp      )) deallocate(this%alsp      )
    if (allocated(this%surfalb   )) deallocate(this%surfalb   )
    if (allocated(this%npcflag   )) deallocate(this%npcflag   )
    if (allocated(this%rhouch    )) deallocate(this%rhouch    )
    if (allocated(this%zin       )) deallocate(this%zin       )
    if (allocated(this%rhosoil   )) deallocate(this%rhosoil   )
    if (allocated(this%cpsoil    )) deallocate(this%cpsoil    )
    if (allocated(this%scond     )) deallocate(this%scond     )
    if (allocated(this%stemp     )) deallocate(this%stemp     )
    if (allocated(this%subflux   )) deallocate(this%subflux   )
    if (allocated(this%gndice    )) deallocate(this%gndice    )
    if (allocated(this%dmadt     )) deallocate(this%dmadt     )
    if (allocated(this%tl        )) deallocate(this%tl        )
    if (allocated(this%teta      )) deallocate(this%teta      )
    if (allocated(this%pl        )) deallocate(this%pl        )
    if (allocated(this%plogadj   )) deallocate(this%plogadj   )
    if (allocated(this%upi       )) deallocate(this%upi       )
    if (allocated(this%vpi       )) deallocate(this%vpi       )
    if (allocated(this%qpi       )) deallocate(this%qpi       )
    if (allocated(this%aadj      )) deallocate(this%aadj      )
    if (allocated(this%badj      )) deallocate(this%badj      )
    if (allocated(this%om        )) deallocate(this%om        )
    if (allocated(this%qh2o      )) deallocate(this%qh2o      )

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

  end subroutine gomars_v1_tend_init

  subroutine gomars_v1_tend_clear(this)

    class(gomars_v1_tend_type), intent(inout) :: this

  end subroutine gomars_v1_tend_clear

  subroutine gomars_v1_tend_final(this)

    type(gomars_v1_tend_type), intent(inout) :: this

    call this%clear()

  end subroutine gomars_v1_tend_final

end module gomars_v1_types_mod