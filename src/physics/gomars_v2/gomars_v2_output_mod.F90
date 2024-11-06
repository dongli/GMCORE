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

module gomars_v2_output_mod

  use fiona
  use gomars_v2_objects_mod

  implicit none

  private

  public gomars_v2_add_output
  public gomars_v2_output

contains

  subroutine gomars_v2_add_output(tag, dtype)

    character(*), intent(in) :: tag
    character(*), intent(in) :: dtype

    character(3) :: dims(2) = ['lon', 'lat']

    call fiona_add_var(tag, 'tin' , long_name='Surface thermal inertia'     , units='', dim_names=dims, dtype=dtype)
    call fiona_add_var(tag, 'gnd_ice' , long_name='Ground ice indicator'     , units='', dim_names=dims, dtype=dtype)

  end subroutine gomars_v2_add_output

  subroutine gomars_v2_output(tag, iblk)

    character(*), intent(in) :: tag
    integer, intent(in) :: iblk

    associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
    call fiona_output(tag, 'tin' , reshape(state%tin, mesh%cell_count_2d(1:2)), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'gnd_ice', reshape(state%gnd_ice, mesh%cell_count_2d(1:2)), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    end associate

  end subroutine gomars_v2_output

end module gomars_v2_output_mod
