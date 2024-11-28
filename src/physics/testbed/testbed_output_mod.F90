module testbed_output_mod

  use fiona
  use testbed_objects_mod

  implicit none

  private

  public testbed_add_output
  public testbed_output

  character(8), parameter ::    cell_dims_2d(3) = ['lon ', 'lat ',         'time']
  character(8), parameter ::    cell_dims_3d(4) = ['lon ', 'lat ', 'lev ', 'time']

contains

  subroutine testbed_add_output(tag, dtype)

    character(*), intent(in) :: tag
    character(*), intent(in) :: dtype

  end subroutine testbed_add_output

  subroutine testbed_output(tag, iblk)

    character(*), intent(in) :: tag
    integer, intent(in) :: iblk

  end subroutine testbed_output

end module testbed_output_mod