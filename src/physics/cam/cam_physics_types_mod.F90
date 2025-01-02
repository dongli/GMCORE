! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module cam_physics_types_mod

  use const_mod
  use tracer_mod
  use physics_types_mod

  implicit none

  private

  public cam_state_type
  public cam_tend_type

  type, extends(physics_state_type) :: cam_state_type
    ! Wind speed at 10 m (m s-1)
    real(r8), allocatable, dimension(:  ) :: wsp10
    ! Upward constituent flux at surface (kg m-2 s-1)
    real(r8), allocatable, dimension(:,:) :: qflx
    ! Albedo for direct short wave radiation
    real(r8), allocatable, dimension(:  ) :: asdir
    ! Albedo for diffuse short wave radiation
    real(r8), allocatable, dimension(:  ) :: asdif
    ! Albedo for direct long wave radiation
    real(r8), allocatable, dimension(:  ) :: aldir
    ! Albedo for diffuse long wave radiation
    real(r8), allocatable, dimension(:  ) :: aldif
  contains
    procedure :: init  => cam_state_init
    procedure :: clear => cam_state_clear
    final cam_state_final
  end type cam_state_type

  type, extends(physics_tend_type) :: cam_tend_type
  contains
    procedure :: init  => cam_tend_init
    procedure :: clear => cam_tend_clear
    procedure :: reset => cam_tend_reset
    final cam_tend_final
  end type cam_tend_type

contains

  subroutine cam_state_init(this, mesh)

    class(cam_state_type), intent(inout), target :: this
    type(physics_mesh_type), intent(in) :: mesh

    call this%clear()

    allocate(this%wsp10(mesh%ncol         ))
    allocate(this%qflx (mesh%ncol,ntracers))
    allocate(this%asdir(mesh%ncol         ))
    allocate(this%asdif(mesh%ncol         ))
    allocate(this%aldir(mesh%ncol         ))
    allocate(this%aldif(mesh%ncol         ))

    call this%physics_state_init(mesh)

  end subroutine cam_state_init

  subroutine cam_state_clear(this)

    class(cam_state_type), intent(inout) :: this

    if (allocated(this%wsp10)) deallocate(this%wsp10)
    if (allocated(this%qflx )) deallocate(this%qflx )
    if (allocated(this%asdir)) deallocate(this%asdir)
    if (allocated(this%asdif)) deallocate(this%asdif)
    if (allocated(this%aldir)) deallocate(this%aldir)
    if (allocated(this%aldif)) deallocate(this%aldif)

    call this%physics_state_clear()

  end subroutine cam_state_clear

  subroutine cam_state_final(this)

    type(cam_state_type), intent(inout) :: this

    call this%clear()

  end subroutine cam_state_final

  subroutine cam_tend_init(this, mesh)

    class(cam_tend_type), intent(inout), target :: this
    type(physics_mesh_type), intent(in) :: mesh

    call this%clear()

    call this%physics_tend_init(mesh)

  end subroutine cam_tend_init

  subroutine cam_tend_clear(this)

    class(cam_tend_type), intent(inout) :: this

    call this%physics_tend_clear()

  end subroutine cam_tend_clear

  subroutine cam_tend_reset(this)

    class(cam_tend_type), intent(inout) :: this

    call this%physics_tend_reset()

  end subroutine cam_tend_reset

  subroutine cam_tend_final(this)

    type(cam_tend_type), intent(inout) :: this

    call this%clear()

  end subroutine cam_tend_final

end module cam_physics_types_mod
