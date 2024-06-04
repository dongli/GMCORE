! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module simple_physics_output_mod

  use fiona
  use simple_physics_objects_mod

  implicit none

  private

  public simple_physics_add_output
  public simple_physics_output

  character(8), parameter ::    cell_dims_2d(3) = ['lon ', 'lat ',         'time']
  character(8), parameter ::    cell_dims_3d(4) = ['lon ', 'lat ', 'lev ', 'time']

contains

  subroutine simple_physics_add_output(tag, dtype)

    character(*), intent(in) :: tag
    character(*), intent(in) :: dtype

    call fiona_add_var(tag, 'precl', long_name='Large scale precipitation', units='m s-1', dim_names=cell_dims_2d)

  end subroutine simple_physics_add_output

  subroutine simple_physics_output(tag, iblk)

    character(*), intent(in) :: tag
    integer, intent(in) :: iblk

    associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
    call fiona_output(tag, 'precl', reshape(state%precl, mesh%cell_count_2d(1:2)), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    end associate

  end subroutine simple_physics_output

end module simple_physics_output_mod
