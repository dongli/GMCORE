module testbed_types_mod

  use const_mod
  use physics_types_mod

  implicit none

  private

  public physics_mesh_type
  public testbed_state_type
  public testbed_tend_type

  type, extends(physics_state_type) :: testbed_state_type
  contains
    procedure :: init  => testbed_state_init
    procedure :: clear => testbed_state_clear
    final testbed_state_final
  end type testbed_state_type

  type, extends(physics_tend_type) :: testbed_tend_type
  contains
    procedure :: init  => testbed_tend_init
    procedure :: clear => testbed_tend_clear
    final testbed_tend_final
  end type testbed_tend_type

contains

  subroutine testbed_state_init(this, mesh)

    class(testbed_state_type), intent(inout), target :: this
    type(physics_mesh_type), intent(in) :: mesh

    call this%clear()

    call this%physics_state_init(mesh)

  end subroutine testbed_state_init

  subroutine testbed_state_clear(this)

    class(testbed_state_type), intent(inout) :: this

    call this%physics_state_clear()

  end subroutine testbed_state_clear

  subroutine testbed_state_final(this)

    type(testbed_state_type), intent(inout) :: this

    call this%clear()

  end subroutine testbed_state_final

  subroutine testbed_tend_init(this, mesh)

    class(testbed_tend_type), intent(inout), target :: this
    type(physics_mesh_type), intent(in) :: mesh

    call this%clear()

    call this%physics_tend_init(mesh)

  end subroutine testbed_tend_init

  subroutine testbed_tend_clear(this)

    class(testbed_tend_type), intent(inout) :: this

    call this%physics_tend_clear()

  end subroutine testbed_tend_clear

  subroutine testbed_tend_final(this)

    type(testbed_tend_type), intent(inout) :: this

    call this%clear()

  end subroutine testbed_tend_final

end module testbed_types_mod