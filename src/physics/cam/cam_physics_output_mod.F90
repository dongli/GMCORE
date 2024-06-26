module cam_physics_output_mod

  use fiona
  use cam_physics_objects_mod

  implicit none

  private

  public cam_physics_add_output
  public cam_physics_output

contains

  subroutine cam_physics_add_output(tag, dtype)

    character(*), intent(in) :: tag
    character(*), intent(in) :: dtype

    character(4) :: dims(3) = ['lon ', 'lat ', 'time']

    call fiona_add_var(tag, 'wsp10', long_name='Wind speed at 10 m', units='m s-1', dim_names=dims, dtype=dtype)
    call fiona_add_var(tag, 'qflx_qv', long_name='Water vapor flux at surface', units='', dim_names=dims, dtype=dtype)
    call fiona_add_var(tag, 'asdir', long_name='Albedo for direct short wave radiation', units='', dim_names=dims, dtype=dtype)
    call fiona_add_var(tag, 'asdif', long_name='Albedo for diffuse short wave radiation', units='', dim_names=dims, dtype=dtype)
    call fiona_add_var(tag, 'aldir', long_name='Albedo for direct long wave radiation', units='', dim_names=dims, dtype=dtype)
    call fiona_add_var(tag, 'aldif', long_name='Albedo for diffuse long wave radiation', units='', dim_names=dims, dtype=dtype)

  end subroutine cam_physics_add_output

  subroutine cam_physics_output(tag, iblk)

    character(*), intent(in) :: tag
    integer, intent(in) :: iblk

    associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
    call fiona_output(tag, 'wsp10', reshape(state%wsp10, mesh%cell_count_2d(1:2)), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'qflx_qv', reshape(state%qflx(:,1), mesh%cell_count_2d(1:2)), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'asdir', reshape(state%asdir, mesh%cell_count_2d(1:2)), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'asdif', reshape(state%asdif, mesh%cell_count_2d(1:2)), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'aldir', reshape(state%aldir, mesh%cell_count_2d(1:2)), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'aldif', reshape(state%aldif, mesh%cell_count_2d(1:2)), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    end associate

  end subroutine cam_physics_output

end module cam_physics_output_mod
