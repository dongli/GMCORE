module testbed_driver_mod

  use const_mod
  use testbed_objects_mod
  use testbed_output_mod

  implicit none

  public testbed_init_stage1
  public testbed_final
  public testbed_run
  public testbed_add_output
  public testbed_output
  public objects

  real(r8) dt

contains

  subroutine testbed_init_stage1(namelist_path, dt_phys)

    character(*), intent(in) :: namelist_path
    real(r8), intent(in) :: dt_phys

    call testbed_final()

    dt = dt_phys

  end subroutine testbed_init_stage1

  subroutine testbed_final()

  end subroutine testbed_final

  subroutine testbed_run()

  end subroutine testbed_run

end module testbed_driver_mod
