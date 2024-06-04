module dyn_comp

  use cam_physics_objects_mod

  implicit none

  private

  public dyn_init
  public dyn_final
  public dyn_register
  public dyn_import_t
  public dyn_export_t

  type pstate_ptr_type
    type(cam_state_type), pointer :: ptr => null()
  end type pstate_ptr_type

  type dyn_import_t
    type(pstate_ptr_type), allocatable :: state(:)
  contains
    final :: dyn_import_final
  end type dyn_import_t

  type ptend_ptr_type
    type(cam_tend_type), pointer :: ptr => null()
  end type ptend_ptr_type

  type dyn_export_t
    type(ptend_ptr_type), allocatable :: tend(:)
  contains
    final :: dyn_export_final
  end type dyn_export_t

contains

  subroutine dyn_init(dyn_in, dyn_out)

    type(dyn_import_t), intent(out) :: dyn_in
    type(dyn_export_t), intent(out) :: dyn_out

    integer iblk

    allocate(dyn_in%state(size(objects)))
    allocate(dyn_out%tend(size(objects)))
    do iblk = 1, size(objects)
      dyn_in%state(iblk)%ptr => objects(iblk)%state
      dyn_out%tend(iblk)%ptr => objects(iblk)%tend
    end do

  end subroutine dyn_init

  subroutine dyn_final()

  end subroutine dyn_final

  subroutine dyn_register()

  end subroutine dyn_register

  subroutine dyn_import_final(this)

    type(dyn_import_t), intent(inout) :: this

    if (allocated(this%state)) deallocate(this%state)

  end subroutine dyn_import_final

  subroutine dyn_export_final(this)

    type(dyn_export_t), intent(inout) :: this

    if (allocated(this%tend)) deallocate(this%tend)

  end subroutine dyn_export_final

end module dyn_comp