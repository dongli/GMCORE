! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   This module writes and reads restart files in NetCDF format.
!
!   When using Intel OneAPI 2021 or so, the sliced array arguments cause
!   segmentation fault, so I used allocated tmp array to avoid this problem.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module restart_mod

  use mpi
  use container
  use fiona
  use string
  use flogger
  use datetime
  use const_mod
  use namelist_mod, start_time_array => start_time, end_time_array => end_time
  use time_mod, old => old_time_idx
  use block_mod
  use tracer_mod
  use latlon_field_types_mod
  use latlon_parallel_mod

  implicit none

  private

  public restart_init_stage2
  public restart_write
  public restart_read

  character(8), parameter ::    cell_dims_2d(3) = ['lon ', 'lat ',         'time']
  character(8), parameter ::    cell_dims_3d(4) = ['lon ', 'lat ', 'lev ', 'time']
  character(8), parameter ::     vtx_dims_2d(3) = ['ilon', 'ilat',         'time']
  character(8), parameter ::     vtx_dims_3d(4) = ['ilon', 'ilat', 'lev ', 'time']
  character(8), parameter ::     lon_dims_2d(3) = ['ilon', 'lat ',         'time']
  character(8), parameter ::     lon_dims_3d(4) = ['ilon', 'lat ', 'lev ', 'time']
  character(8), parameter ::     lat_dims_2d(3) = ['lon ', 'ilat',         'time']
  character(8), parameter ::     lat_dims_3d(4) = ['lon ', 'ilat', 'lev ', 'time']
  character(8), parameter ::     lev_dims_3d(4) = ['lon ', 'lat ', 'ilev', 'time']

contains

  subroutine restart_init_stage2()

    character(10) time_value, time_units
    real(r8) seconds

    if (restart_interval == 'N/A') then
      if (proc%is_root()) call log_warning('Parameter restart_interval is not set, so no restart file outputted.')
      return
    end if
    if (case_name == 'N/A') call log_error('Parameter case_name is not set!')

    time_value = split_string(restart_interval, ' ', 1)
    time_units = split_string(restart_interval, ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days', 'sol')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case default
      call log_error('Invalid restart interval ' // trim(restart_interval) // '!')
    end select

    call time_add_alert('restart_write', seconds=seconds)

  end subroutine restart_init_stage2

  subroutine restart_write(itime)

    integer, intent(in) :: itime

    integer iblk, i
    real(8) time1, time2
    class(*), pointer :: field

    if (proc%is_root()) then
      call log_notice('Write restart.')
      time1 = MPI_WTIME()
    end if

    call fiona_create_dataset('r0', desc=case_desc, file_prefix=trim(case_name) // '.' // trim(curr_time_str), &
      mpi_comm=proc%comm, ngroup=output_ngroups)

    call fiona_add_att('r0', 'start_time', start_time%isoformat())
    call fiona_add_att('r0', 'time_step_size', dt_dyn)
    call fiona_add_att('r0', 'restart_interval', restart_interval)
    call fiona_add_dim('r0', 'time'     , add_var=.true.)
    call fiona_add_dim('r0', 'lon'      , size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('r0', 'lat'      , size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('r0', 'lev'      , size=global_mesh%full_nlev)
    call fiona_add_dim('r0', 'ilon'     , size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('r0', 'ilat'     , size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('r0', 'ilev'     , size=global_mesh%half_nlev)
    call fiona_add_var('r0', 'time_step', long_name='', units='', dim_names=['time'], dtype='i4')

    call add_fields('r0', blocks(1)%dstate(1)%fields)
    call add_fields('r0', blocks(1)%static   %fields, static=.true.)
    call add_fields('r0', blocks(1)%aux      %fields)
    field => tracers(1)%q
    call add_field ('r0', field)
    do i = 1, nbatches
      call add_fields('r0', blocks(1)%adv_batches(i)%fields)
    end do

    call fiona_start_output('r0', dble(elapsed_seconds), new_file=.true.)
    call fiona_output('r0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('r0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('r0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('r0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat))
    call fiona_output('r0', 'time_step', time_step)

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh         , &
                 dstate => blocks(iblk)%dstate(itime), &
                 static => blocks(iblk)%static       , &
                 aux    => blocks(iblk)%aux          , &
                 q       => tracers(iblk)%q          )
      call write_fields('r0', mesh, dstate%fields)
      call write_fields('r0', mesh, static%fields)
      call write_fields('r0', mesh, aux   %fields)
      call write_field ('r0', mesh, q)
      do i = 1, nbatches
        associate (batch => blocks(iblk)%adv_batches(i))
        if (batch%step /= -1) then
          call log_error('Restart advection batch ' // trim(batch%name) // ', but its step is not -1!', pid=proc%id)
        end if
        call write_fields('r0', mesh, batch%fields)
        end associate
      end do
      end associate
    end do

    call fiona_end_output('r0')

    if (proc%is_root()) then
      time2 = MPI_WTIME()
      call log_notice('Done write restart cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine restart_write

  subroutine restart_read()

    type(datetime_type) time
    integer iblk, i
    real(r8) time_value, time1, time2
    character(50) time_units
    class(*), pointer :: field

    character(30) start_time_str

    if (restart_file == 'N/A') then
      call log_error('Parameter restart_file is needed to restart!')
    end if

    if (proc%is_root()) then
      call log_notice('Read restart file ' // trim(restart_file) // '.')
      time1 = MPI_WTIME()
    end if

    call fiona_open_dataset('r0', file_path=restart_file, mpi_comm=proc%comm, ngroup=input_ngroups)
    call fiona_start_input('r0')

    call fiona_get_att('r0', 'start_time', start_time_str)
    call time_fast_forward(start_time_str, change_end_time=.true.)
    call fiona_input('r0', 'time', time_value)
    call fiona_get_att('r0', 'time', 'units', time_units)
    call fiona_input('r0', 'time_step', time_step)
    call time_fast_forward(time_value, time_units, change_end_time=.false.)

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh       , &
                 dstate => blocks(iblk)%dstate(old), &
                 static => blocks(iblk)%static     , &
                 aux    => blocks(iblk)%aux        , &
                 q      => tracers(iblk)%q         )
      call read_fields('r0', mesh, dstate%fields)
      call read_fields('r0', mesh, static%fields)
      call read_fields('r0', mesh, aux   %fields)
      field => tracers(old)%q
      call read_field ('r0', mesh, q)
      do i = 1, nbatches
        associate (batch => blocks(iblk)%adv_batches(i))
        call read_fields('r0', mesh, batch%fields)
        batch%step = -1
        end associate
      end do
      end associate
    end do

    call fiona_end_input('r0')

    if (proc%is_root()) then
      time2 = MPI_WTIME()
      call log_notice('Restart to ' // trim(curr_time_str) // ' cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine restart_read

  subroutine add_fields(dtag, fields, static)

    character(*), intent(in) :: dtag
    type(array_type), intent(in) :: fields
    logical, intent(in), optional :: static

    class(*), pointer :: field
    integer i

    do i = 1, fields%size
      field => fields%value_at(i)
      call add_field(dtag, field, static)
    end do

  end subroutine add_fields

  subroutine add_field(dtag, field, static)

    character(*), intent(in) :: dtag
    class(*), pointer :: field
    logical, intent(in), optional :: static

    integer i, n

    select type (field)
    type is (latlon_field2d_type)
      if (.not. field%restart) return
      n = merge(2, 3, present(static))
      select case (field%loc)
      case ('cell')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names=cell_dims_2d(:n), dtype='r8')
      case ('lon')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names= lon_dims_2d(:n), dtype='r8')
      case ('lat')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names= lat_dims_2d(:n), dtype='r8')
      case ('vtx')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names= vtx_dims_2d(:n), dtype='r8')
      end select
    type is (latlon_field3d_type)
      if (.not. field%restart) return
      n = merge(3, 4, present(static))
      select case (field%loc)
      case ('cell')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names=cell_dims_3d(:n), dtype='r8')
      case ('lon')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names= lon_dims_3d(:n), dtype='r8')
      case ('lat')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names= lat_dims_3d(:n), dtype='r8')
      case ('vtx')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names= vtx_dims_3d(:n), dtype='r8')
      case ('lev')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names= lev_dims_3d(:n), dtype='r8')
      end select
    type is (latlon_field4d_type)
      if (.not. field%restart) return
      do i = 1, field%dim4_size
        call fiona_add_var(dtag, trim(field%name) // '_' // trim(field%var4_names(i)), long_name=trim(field%long_name) // ' of ' // trim(field%var4_names(i)), &
          units=field%units, dim_names=cell_dims_3d, dtype='r8')
      end do
    end select

  end subroutine add_field

  subroutine write_fields(dtag, mesh, fields)

    character(*), intent(in) :: dtag
    type(latlon_mesh_type), intent(in) :: mesh
    type(array_type), intent(in) :: fields

    class(*), pointer :: field
    integer i

    do i = 1, fields%size
      field => fields%value_at(i)
      call write_field(dtag, mesh, field)
    end do

  end subroutine write_fields

  subroutine write_field(dtag, mesh, field)

    character(*), intent(in) :: dtag
    type(latlon_mesh_type), intent(in) :: mesh
    class(*), intent(in) :: field

    integer is, ie, js, je, ks, ke, i
    integer start2d(2), count2d(2), start3d(3), count3d(3)

    select type (field)
    type is (latlon_field2d_type)
      if (.not. field%restart) return
      select case (field%loc)
      case ('cell')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
      case ('lon')
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%full_jds; je = mesh%full_jde
      case ('lat')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%half_jds; je = mesh%half_jde
      case ('vtx')
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%half_jds; je = mesh%half_jde
      end select
      start2d = [is,js]
      count2d = [ie-is+1,je-js+1]
      call fiona_output(dtag, field%name, field%d(is:ie,js:je), start=start2d, count=count2d)
    type is (latlon_field3d_type)
      if (.not. field%restart) return
      select case (field%loc)
      case ('cell')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
      case ('lon')
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
      case ('lat')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%half_jds; je = mesh%half_jde
        ks = mesh%full_kds; ke = mesh%full_kde
      case ('lev')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%half_kds; ke = mesh%half_kde
      case ('vtx')
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%half_jds; je = mesh%half_jde
        ks = mesh%full_kds; ke = mesh%full_kde
      end select
      start3d = [is,js,ks]
      count3d = [ie-is+1,je-js+1,ke-ks+1]
      call fiona_output(dtag, field%name, field%d(is:ie,js:je,ks:ke), start=start3d, count=count3d)
    type is (latlon_field4d_type)
      if (.not. field%restart) return
      select case (field%loc)
      case ('cell')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
      end select
      start3d = [is,js,ks]
      count3d = [ie-is+1,je-js+1,ke-ks+1]
      do i = 1, field%dim4_size
        call fiona_output(dtag, trim(field%name) // '_' // trim(field%var4_names(i)), &
          field%d(is:ie,js:je,ks:ke,i), start=start3d, count=count3d)
      end do
    end select

  end subroutine write_field

  subroutine read_fields(dtag, mesh, fields)

    character(*), intent(in) :: dtag
    type(latlon_mesh_type), intent(in) :: mesh
    type(array_type), intent(in) :: fields

    class(*), pointer :: field
    integer i

    do i = 1, fields%size
      field => fields%value_at(i)
      call read_field(dtag, mesh, field)
    end do

  end subroutine read_fields

  subroutine read_field(dtag, mesh, field)

    character(*), intent(in) :: dtag
    type(latlon_mesh_type), intent(in) :: mesh
    class(*), intent(in) :: field

    integer is, ie, js, je, ks, ke, i
    integer start2d(2), count2d(2), start3d(3), count3d(3)

    select type (field)
    type is (latlon_field2d_type)
      if (.not. field%restart) return
      select case (field%loc)
      case ('cell')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
      case ('lon')
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%full_jds; je = mesh%full_jde
      case ('lat')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%half_jds; je = mesh%half_jde
      case ('vtx')
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%half_jds; je = mesh%half_jde
      end select
      start2d = [is,js]
      count2d = [ie-is+1,je-js+1]
      call fiona_input(dtag, field%name, field%d(is:ie,js:je), start=start2d, count=count2d)
      call fill_halo(field)
    type is (latlon_field3d_type)
      if (.not. field%restart) return
      select case (field%loc)
      case ('cell')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
      case ('lon')
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
      case ('lat')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%half_jds; je = mesh%half_jde
        ks = mesh%full_kds; ke = mesh%full_kde
      case ('lev')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%half_kds; ke = mesh%half_kde
      case ('vtx')
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%half_jds; je = mesh%half_jde
        ks = mesh%full_kds; ke = mesh%full_kde
      end select
      start3d = [is,js,ks]
      count3d = [ie-is+1,je-js+1,ke-ks+1]
      call fiona_input(dtag, field%name, field%d(is:ie,js:je,ks:ke), start=start3d, count=count3d)
      call fill_halo(field)
    type is (latlon_field4d_type)
      if (.not. field%restart) return
      select case (field%loc)
      case ('cell')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
      end select
      start3d = [is,js,ks]
      count3d = [ie-is+1,je-js+1,ke-ks+1]
      do i = 1, field%dim4_size
        call fiona_input(dtag, trim(field%name) // '_' // trim(field%var4_names(i)), &
          field%d(is:ie,js:je,ks:ke,i), start=start3d, count=count3d)
        call fill_halo(field, i)
      end do
    end select

  end subroutine read_field

end module restart_mod
