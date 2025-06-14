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

module gomars_v2_driver_mod

  use datetime
  use tracer_mod
  use gomars_v2_const_mod
  use gomars_v2_namelist_mod
  use gomars_v2_tracers_mod
  use gomars_v2_types_mod
  use gomars_v2_output_mod
  use gomars_v2_objects_mod
  use gomars_v2_orbit_mod
  use gomars_v2_solar_mod
  use gomars_v2_rad_mod
  use gomars_v2_lsm_mod
  use gomars_v2_pbl_mod
  use gomars_v2_mp_mod

  implicit none

  private

  public gomars_v2_init_stage2
  public gomars_v2_init_stage3
  public gomars_v2_final
  public gomars_v2_run
  public gomars_v2_d2p
  public gomars_v2_p2d
  public gomars_v2_add_output
  public gomars_v2_output
  public objects

  real(r8) dt

contains

  subroutine gomars_v2_init_stage2(namelist_path, mesh, dt_adv, dt_phys, input_ngroup, model_root)

    character(*), intent(in) :: namelist_path
    type(physics_mesh_type), intent(in), target :: mesh(:)
    real(r8), intent(in) :: dt_adv
    real(r8), intent(in) :: dt_phys
    integer , intent(in) :: input_ngroup
    character(*), intent(in), optional :: model_root

    call gomars_v2_final()
    call gomars_v2_parse_namelist(namelist_path, model_root)
    call gomars_v2_tracers_init(dt_adv)
    call gomars_v2_objects_init(mesh)
    call gomars_v2_read_static_data(input_ngroup)
    call gomars_v2_lsm_init()
    call gomars_v2_rad_init(mesh(1)%nlev)
    call mars_orbit_init()

  end subroutine gomars_v2_init_stage2

  subroutine gomars_v2_init_stage3()

    integer iblk, icol

    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
      do icol = 1, mesh%ncol
        state%t_top(icol) = state%t(icol,1)
        state%tg   (icol) = state%t(icol,mesh%nlev)
      end do
      end associate
    end do

  end subroutine gomars_v2_init_stage3

  subroutine gomars_v2_final()

    call gomars_v2_rad_final()
    call gomars_v2_objects_final()
    call gomars_v2_lsm_final()
    call gomars_v2_pbl_final()

  end subroutine gomars_v2_final

  subroutine gomars_v2_run(time)

    type(datetime_type), intent(in) :: time

    integer iblk, icol
    real(r8) ls

    ls = time%solar_longitude()

    call update_solar_decl_angle(ls)
    call update_solar_flux(ls)

    do iblk = 1, size(objects)
      associate (mesh   => objects(iblk)%mesh  , &
                 state  => objects(iblk)%state , &
                 tend   => objects(iblk)%tend  )
      do icol = 1, objects(iblk)%mesh%ncol
        state%cosz(icol) = solar_cos_zenith_angle(mesh%lon(icol), mesh%lat(icol), time%time_of_day())
      end do
      call gomars_v2_pbl_run(state, tend, dt)
      call gomars_v2_lsm_run(state, tend)
      call gomars_v2_rad_run(state, tend)
      end associate
    end do

  end subroutine gomars_v2_run

  subroutine gomars_v2_read_static_data(input_ngroups)

    use latlon_interp_mod

    integer , intent(in) :: input_ngroups

    real(r8), allocatable :: lon(:)
    real(r8), allocatable :: lat(:)
    real(r8), allocatable :: array(:,:)
    integer iblk

    if (allocated(objects)) then
      associate (mesh => objects(1)%mesh, state => objects(1)%state)
      ! Surface albedo
      call fiona_open_dataset('alb', file_path=albedo_file, mpi_comm=proc%comm_model, ngroups=input_ngroups)
      call fiona_set_dim('alb', 'lon', span=[0, 360], cyclic=.true.)
      call fiona_set_dim('alb', 'lat', span=[-90, 90])
      call fiona_start_input('alb')
      call fiona_input_range('alb', 'lon', lon, coord_range=[mesh%min_lon, mesh%max_lon]); lon = lon * rad
      call fiona_input_range('alb', 'lat', lat, coord_range=[mesh%min_lat, mesh%max_lat]); lat = lat * rad
      call fiona_input_range('alb', 'albedo', array, coord_range_1=[mesh%min_lon, mesh%max_lon], coord_range_2=[min_lat, max_lat])
      call fiona_end_input('alb')
      call latlon_interp_bilinear_column(lon, lat, array, mesh%lon, mesh%lat, state%alb)
      deallocate(lon, lat, array)

      ! Surface thermal inertia
      call fiona_open_dataset('tin', file_path=thermal_inertia_file, mpi_comm=proc%comm_model, ngroups=input_ngroups)
      call fiona_set_dim('tin', 'lon', span=[0, 360], cyclic=.true.)
      call fiona_set_dim('tin', 'lat', span=[-90, 90])
      call fiona_start_input('tin')
      call fiona_input_range('tin', 'lon', lon, coord_range=[min_lon, max_lon]); lon = lon * rad
      call fiona_input_range('tin', 'lat', lat, coord_range=[min_lat, max_lat]); lat = lat * rad
      call fiona_input_range('tin', 'thin', array, coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
      call fiona_end_input('tin')
      call latlon_interp_bilinear_column(lon, lat, array, mesh%lon, mesh%lat, state%tin)
      deallocate(lon, lat, array)

      ! Ground ice indicator
      call fiona_open_dataset('gnd_ice', file_path=gnd_ice_file, mpi_comm=proc%comm_model, ngroups=input_ngroups)
      call fiona_set_dim('gnd_ice', 'lon', span=[0, 360], cyclic=.true.)
      call fiona_set_dim('gnd_ice', 'lat', span=[-90, 90])
      call fiona_start_input('gnd_ice')
      call fiona_input_range('gnd_ice', 'lon', lon, coord_range=[min_lon, max_lon]); lon = lon * rad
      call fiona_input_range('gnd_ice', 'lat', lat, coord_range=[min_lat, max_lat]); lat = lat * rad
      call fiona_input_range('gnd_ice', 'gnd_ice', array, coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
      call fiona_end_input('gnd_ice')
      call latlon_interp_bilinear_column(lon, lat, array, mesh%lon, mesh%lat, state%gnd_ice)
      deallocate(lon, lat, array)
      end associate
    end if

  end subroutine gomars_v2_read_static_data

  subroutine gomars_v2_d2p()

    ! Calculate ts.

  end subroutine gomars_v2_d2p

  subroutine gomars_v2_p2d()

  end subroutine gomars_v2_p2d

end module gomars_v2_driver_mod
