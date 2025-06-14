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
  character(8), parameter :: lev_lon_dims_3d(4) = ['ilon', 'lat ', 'ilev', 'time']
  character(8), parameter :: lev_lat_dims_3d(4) = ['lon ', 'ilat', 'ilev', 'time']

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
      mpi_comm=proc%comm_model, ngroups=output_ngroups)

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

    call fiona_add_dim('r0', 'adv_batch', size=nbatches)
    call fiona_add_dim('r0', 'ntracers' , size=ntracers)
    if (nbatches > 0) then
      call fiona_add_dim('r0', 'strlen_scheme'     , size=strlen_scheme   )
      call fiona_add_dim('r0', 'strlen_loc'        , size=strlen_loc      )
      call fiona_add_dim('r0', 'strlen_name'       , size=strlen_name     )
      call fiona_add_dim('r0', 'strlen_long_name'  , size=strlen_long_name)
      call fiona_add_dim('r0', 'strlen_units'      , size=strlen_units    )
      call fiona_add_var('r0', 'adv_batch_scheme_h', long_name='Tracer advection horizontal scheme', units='' , dim_names=['strlen_scheme   ', 'adv_batch       '], dtype='s' )
      call fiona_add_var('r0', 'adv_batch_scheme_v', long_name='Tracer advection vertical scheme'  , units='' , dim_names=['strlen_scheme   ', 'adv_batch       '], dtype='s' )
      call fiona_add_var('r0', 'adv_batch_loc'     , long_name='Tracer advection location'         , units='' , dim_names=['strlen_loc      ', 'adv_batch       '], dtype='s' )
      call fiona_add_var('r0', 'adv_batch_names'   , long_name='Tracer advection batch name'       , units='' , dim_names=['strlen_name     ', 'adv_batch       '], dtype='s' )
      call fiona_add_var('r0', 'adv_batch_dt'      , long_name='Tracer advection time step'        , units='s', dim_names=[                    'adv_batch       '], dtype='r8')
      call fiona_add_var('r0', 'adv_batch_ntracers', long_name='Number of tracers in each batch'   , units='' , dim_names=[                    'adv_batch       '], dtype='i4')
      call fiona_add_var('r0', 'adv_batch_idx'     , long_name='Tracer indices for each batch'     , units='' , dim_names=['ntracers        ', 'adv_batch       '], dtype='i4')
      call fiona_add_var('r0', 'tracer_names'      , long_name='Tracer names'                      , units='' , dim_names=['strlen_name     ', 'ntracers        '], dtype='s' )
      call fiona_add_var('r0', 'tracer_long_names' , long_name='Tracer long names'                 , units='' , dim_names=['strlen_long_name', 'ntracers        '], dtype='s' )
      call fiona_add_var('r0', 'tracer_units'      , long_name='Tracer units'                      , units='' , dim_names=['strlen_units    ', 'ntracers        '], dtype='s' )
      call fiona_add_var('r0', 'tracer_types'      , long_name='Tracer types'                      , units='' , dim_names=[                    'ntracers        '], dtype='i4')
    end if

    call add_fields('r0', blocks(1)%dstate(1)%fields)
    call add_fields('r0', blocks(1)%static   %fields, static=.true.)
    call add_fields('r0', blocks(1)%aux      %fields)
    field => tracers(1)%q
    call add_field ('r0', field)
    field => tracers(1)%qm
    call add_field ('r0', field)
    field => tracers(1)%qm_lev
    call add_field ('r0', field)
    do i = 1, nbatches
      call add_fields('r0', blocks(1)%adv_batches(i)%fields)
    end do

    call fiona_start_output('r0', dble(elapsed_seconds), new_file=.true.)
    call fiona_output('r0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%full_nlon), only_root=.true.)
    call fiona_output('r0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%full_nlat), only_root=.true.)
    call fiona_output('r0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon), only_root=.true.)
    call fiona_output('r0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat), only_root=.true.)
    call fiona_output('r0', 'time_step', time_step)

    if (nbatches > 0) then
      do i = 1, nbatches
        associate (batch => blocks(1)%adv_batches(i))
        call fiona_output('r0', 'adv_batch_scheme_h', batch%scheme_h, start=[1,i], count=[strlen_scheme ,1])
        call fiona_output('r0', 'adv_batch_scheme_v', batch%scheme_v, start=[1,i], count=[strlen_scheme ,1])
        call fiona_output('r0', 'adv_batch_loc'     , batch%loc     , start=[1,i], count=[strlen_loc    ,1])
        call fiona_output('r0', 'adv_batch_names'   , batch%name    , start=[1,i], count=[strlen_name   ,1])
        call fiona_output('r0', 'adv_batch_dt'      , batch%dt      )
        call fiona_output('r0', 'adv_batch_ntracers', batch%ntracers)
        call fiona_output('r0', 'adv_batch_idx'     , batch%idx     , start=[1,i], count=[batch%ntracers,1])
        end associate
      end do
      do i = 1, ntracers
        call fiona_output('r0', 'tracer_names'      , tracer_names     (i), start=[1,i], count=[strlen_name     ,1])
        call fiona_output('r0', 'tracer_long_names' , tracer_long_names(i), start=[1,i], count=[strlen_long_name,1])
        call fiona_output('r0', 'tracer_units'      , tracer_units     (i), start=[1,i], count=[strlen_units    ,1])
        call fiona_output('r0', 'tracer_types'      , tracer_types     (i), start=[  i])
      end do
    end if

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh         , &
                 dstate => blocks(iblk)%dstate(itime), &
                 static => blocks(iblk)%static       , &
                 aux    => blocks(iblk)%aux          , &
                 q      => tracers(iblk)%q           , &
                 qm     => tracers(iblk)%qm          , &
                 qm_lev => tracers(iblk)%qm_lev      )
      call write_fields('r0', mesh, dstate%fields)
      call write_fields('r0', mesh, static%fields)
      call write_fields('r0', mesh, aux   %fields)
      call write_field ('r0', mesh, q            )
      call write_field ('r0', mesh, qm           )
      call write_field ('r0', mesh, qm_lev       )
      do i = 1, nbatches
        associate (batch => blocks(iblk)%adv_batches(i))
        if (batch%step /= 0) then
          if (proc%is_root()) call log_error('Restart advection batch ' // trim(batch%name) // ', but its step is not 0!')
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
    integer iblk, i, j, n
    real(r8) time_value, time1, time2, dt
    character(50) time_units
    class(*), pointer :: field
    character(30) start_time_str
    character(strlen_scheme) scheme_h, scheme_v
    character(strlen_loc   ) loc
    character(strlen_name  ) name
    integer, allocatable :: idx(:)

    if (restart_file == 'N/A') then
      call log_error('Parameter restart_file is needed to restart!', pid=proc%id_model)
    end if

    if (proc%is_root()) then
      call log_notice('Read restart file ' // trim(restart_file) // '.')
      time1 = MPI_WTIME()
    end if

    call fiona_open_dataset('r0', file_path=restart_file, mpi_comm=proc%comm_model, ngroups=input_ngroups)
    call fiona_start_input('r0')

    call fiona_get_att('r0', 'start_time', start_time_str)
    call time_fast_forward(start_time_str, change_end_time=.true.)
    call fiona_input('r0', 'time', time_value)
    call fiona_get_att('r0', 'time', 'units', time_units)
    call fiona_input('r0', 'time_step', time_step)
    call time_fast_forward(time_value, time_units, change_end_time=.false.)

    select case (split_string(time_units, ' ', 1))
    case ('days', 'sol')
      elapsed_seconds = time_value * 86400.0
    case ('hours')
      elapsed_seconds = time_value * 3600.0
    case ('minutes')
      elapsed_seconds = time_value * 60.0
    case ('seconds')
      elapsed_seconds = time_value
    case default
      if (proc%is_root()) call log_error('Invalid time_units ' // trim(time_units) // '!', __FILE__, __LINE__)
    end select

    call fiona_get_dim('r0', 'adv_batch', nbatches)
    call fiona_get_dim('r0', 'ntracers' , ntracers)
    if (nbatches > 0) then
      if (allocated(tracer_names     )) deallocate(tracer_names     ); allocate(tracer_names     (ntracers))
      if (allocated(tracer_long_names)) deallocate(tracer_long_names); allocate(tracer_long_names(ntracers))
      if (allocated(tracer_units     )) deallocate(tracer_units     ); allocate(tracer_units     (ntracers))
      if (allocated(tracer_types     )) deallocate(tracer_types     ); allocate(tracer_types     (ntracers))
      do i = 1, ntracers
        call fiona_input('r0', 'tracer_names'     , tracer_names     (i), start=[1,i], count=[strlen_name     ,1])
        call fiona_input('r0', 'tracer_long_names', tracer_long_names(i), start=[1,i], count=[strlen_long_name,1])
        call fiona_input('r0', 'tracer_units'     , tracer_units     (i), start=[1,i], count=[strlen_units    ,1])
        call fiona_input('r0', 'tracer_types'     , tracer_types     (i), start=[  i])
        call tracer_catalog(i)
      end do
      allocate(idx(ntracers))
      do iblk = 1, size(blocks)
        allocate(blocks(iblk)%adv_batches(nbatches))
        do i = 1, nbatches
          associate (batch => blocks(iblk)%adv_batches(i))
          call fiona_input('r0', 'adv_batch_scheme_h', scheme_h, start=[1,i], count=[strlen_scheme,1])
          call fiona_input('r0', 'adv_batch_scheme_v', scheme_v, start=[1,i], count=[strlen_scheme,1])
          call fiona_input('r0', 'adv_batch_loc'     , loc     , start=[1,i], count=[strlen_loc   ,1])
          call fiona_input('r0', 'adv_batch_names'   , name    , start=[1,i], count=[strlen_name  ,1])
          call fiona_input('r0', 'adv_batch_dt'      , dt      )
          call fiona_input('r0', 'adv_batch_ntracers', n       )
          call fiona_input('r0', 'adv_batch_idx'     , idx     , start=[1,i], count=[ntracers     ,1])
          call batch%init(                                      &
            blocks(iblk)%big_filter                           , &
            blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, &
            blocks(iblk)%mesh, blocks(iblk)%halo              , &
            trim(scheme_h)//':'//trim(scheme_v), loc, name, dt, &
            dynamic=.false., passive=.true., idx=idx(1:n)     , &
            bg=blocks(iblk)%adv_batch_bg)
          end associate
        end do
      end do
      deallocate(idx)

      if (proc%is_root()) then
        call log_notice('There are ' // to_str(size(blocks(1)%adv_batches)) // ' advection batches.')
        do i = 1, size(blocks(1)%adv_batches)
          write(*, *) '- ', trim(blocks(1)%adv_batches(i)%name), ' dt_adv = ', to_str(int(blocks(1)%adv_batches(i)%dt))
          do j = 1, blocks(1)%adv_batches(i)%ntracers
            write(*, *) '  * ', trim(tracer_names(blocks(1)%adv_batches(i)%idx(j)))
          end do
        end do
      end if

      call tracer_init_stage2()
    end if

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh       , &
                 dstate => blocks(iblk)%dstate(old), &
                 static => blocks(iblk)%static     , &
                 aux    => blocks(iblk)%aux        , &
                 q      => tracers(iblk)%q         , &
                 qm     => tracers(iblk)%qm        , &
                 qm_lev => tracers(iblk)%qm_lev    )
      call read_fields('r0', mesh, dstate%fields)
      call read_fields('r0', mesh, static%fields)
      call read_fields('r0', mesh, aux   %fields)
      field => tracers(iblk)%q
      call read_field ('r0', mesh, q     )
      field => tracers(iblk)%qm
      call read_field ('r0', mesh, qm    )
      field => tracers(iblk)%qm_lev
      call read_field ('r0', mesh, qm_lev)
      ! FIXME: We need to aqcuire tracer advection batches information from restart file, and register them.
      do i = 1, nbatches
        associate (batch => blocks(iblk)%adv_batches(i))
        call read_fields('r0', mesh, batch%fields)
        batch%step = 0
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
      case ('lev_lon')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names=lev_lon_dims_3d(:n), dtype='r8')
      case ('lev_lat')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names=lev_lat_dims_3d(:n), dtype='r8')
      case default
        call log_error('Invalid field location ' // trim(field%loc) // '!', pid=proc%id_model)
      end select
    type is (latlon_field4d_type)
      if (.not. field%restart) return
      do i = 1, field%dim4_size
        call fiona_add_var(dtag, field%var4_names(i), long_name=trim(field%long_name) // ' of ' // trim(field%var4_names(i)), &
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
      case ('lev_lon')
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%half_kds; ke = mesh%half_kde
      case ('lev_lat')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%half_jds; je = mesh%half_jde
        ks = mesh%half_kds; ke = mesh%half_kde
      case default
        call log_error('Invalid field location ' // trim(field%loc) // '!', pid=proc%id_model)
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
        call fiona_output(dtag, field%var4_names(i), field%d(is:ie,js:je,ks:ke,i), start=start3d, count=count3d)
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
    class(*), intent(inout) :: field

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
      case ('lev_lon')
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%half_kds; ke = mesh%half_kde
      case ('lev_lat')
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%half_jds; je = mesh%half_jde
        ks = mesh%half_kds; ke = mesh%half_kde
      case default
        call log_error('Invalid field location ' // trim(field%loc) // '!', pid=proc%id_model)
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
        call fiona_input(dtag, field%var4_names(i), field%d(is:ie,js:je,ks:ke,i), start=start3d, count=count3d)
        call fill_halo(field, i)
      end do
    end select

  end subroutine read_field

end module restart_mod
