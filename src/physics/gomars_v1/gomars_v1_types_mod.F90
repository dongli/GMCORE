module gomars_v1_types_mod

  use physics_types_mod
  use gomars_v1_const_mod
  use nasa_rad_mod

  implicit none

  private

  public gomars_v1_state_type
  public gomars_v1_tend_type

  type, extends(physics_state_type) :: gomars_v1_state_type
    real(r8), allocatable, dimension(:      ) :: latheat
    real(r8), allocatable, dimension(:,:    ) :: atmcond
    real(r8), allocatable, dimension(  :,:  ) :: qxvdst
    real(r8), allocatable, dimension(  :,:  ) :: qsvdst
    real(r8), allocatable, dimension(  :,:  ) :: gvdst
    real(r8), allocatable, dimension(  :,:  ) :: qxidst
    real(r8), allocatable, dimension(  :,:  ) :: qsidst
    real(r8), allocatable, dimension(  :,:  ) :: gidst
    real(r8), allocatable, dimension(  :    ) :: qextrefdst
    real(r8), allocatable, dimension(  :,:  ) :: qxvcld
    real(r8), allocatable, dimension(  :,:  ) :: qsvcld
    real(r8), allocatable, dimension(  :,:  ) :: gvcld
    real(r8), allocatable, dimension(  :,:  ) :: qxicld
    real(r8), allocatable, dimension(  :,:  ) :: qsicld
    real(r8), allocatable, dimension(  :,:  ) :: gicld
    real(r8), allocatable, dimension(  :    ) :: qextrefcld
    real(r8), allocatable, dimension(  :    ) :: taurefcld
    ! Downward diffuse visible (solar) flux at the surface
    real(r8), allocatable, dimension(:      ) :: dndiffv
    ! Downward IR flux at the surface
    real(r8), allocatable, dimension(:      ) :: dnirflux
    ! Delta-Eddington optical depth (???)
    real(r8), allocatable, dimension(:,  :,:) :: detau
    ! Solar flux at the current Mars distance
    real(r8), allocatable, dimension(    :  ) :: solar
    ! Mars distance from sun
    real(r8) :: rsdist
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

    allocate(this%latheat(mesh%ncol          )); this%latheat = 0
    allocate(this%atmcond(mesh%ncol,mesh%nlev)); this%atmcond = 0

    allocate(this%qxvdst    (mesh%nlev+1,l_nspectv)); this%qxvdst     = 0
    allocate(this%qsvdst    (mesh%nlev+1,l_nspectv)); this%qsvdst     = 0
    allocate(this%gvdst     (mesh%nlev+1,l_nspectv)); this%gvdst      = 0
    allocate(this%qxidst    (mesh%nlev+1,l_nspecti)); this%qxidst     = 0
    allocate(this%qsidst    (mesh%nlev+1,l_nspecti)); this%qsidst     = 0
    allocate(this%gidst     (mesh%nlev+1,l_nspecti)); this%gidst      = 0
    allocate(this%qextrefdst(mesh%nlev+1          )); this%qextrefdst = 0

    allocate(this%qxvcld    (mesh%nlev+1,l_nspectv)); this%qxvcld     = 0
    allocate(this%qsvcld    (mesh%nlev+1,l_nspectv)); this%qsvcld     = 0
    allocate(this%gvcld     (mesh%nlev+1,l_nspectv)); this%gvcld      = 0
    allocate(this%qxicld    (mesh%nlev+1,l_nspecti)); this%qxicld     = 0
    allocate(this%qsicld    (mesh%nlev+1,l_nspecti)); this%qsicld     = 0
    allocate(this%gicld     (mesh%nlev+1,l_nspecti)); this%gicld      = 0
    allocate(this%qextrefcld(mesh%nlev+1          )); this%qextrefcld = 0
    allocate(this%taurefcld (mesh%nlev+1          )); this%taurefcld  = 0

    allocate(this%dndiffv   (mesh%ncol)); this%dndiffv  = 0
    allocate(this%dnirflux  (mesh%ncol)); this%dnirflux = 0

    allocate(this%detau(mesh%ncol,l_nspectv,l_ngauss)); this%detau = 0

    allocate(this%solar(l_nspectv)); this%solar = 0

    call this%physics_state_init(mesh)

  end subroutine gomars_v1_state_init

  subroutine gomars_v1_state_clear(this)

    class(gomars_v1_state_type), intent(inout) :: this

    if (allocated(this%latheat   )) deallocate(this%latheat   )
    if (allocated(this%atmcond   )) deallocate(this%atmcond   )
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
    if (allocated(this%taurefcld )) deallocate(this%taurefcld )
    if (allocated(this%dndiffv   )) deallocate(this%dndiffv   )
    if (allocated(this%dnirflux  )) deallocate(this%dnirflux  )
    if (allocated(this%detau     )) deallocate(this%detau     )
    if (allocated(this%solar     )) deallocate(this%solar     )

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