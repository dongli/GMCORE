! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module history_mod

  use mpi
  use container
  use fiona
  use flogger
  use string
  use const_mod
  use namelist_mod
  use time_mod
  use latlon_parallel_mod
  use block_mod
  use tracer_mod
  use operators_mod
  use physics_mod
  use regrid_mod

  implicit none

  private

  public history_init_stage1
  public history_init_stage2
  public history_init_stage3
  public history_final
  public history_write_h0
  public history_write_h1
  public history_write_h2

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

  subroutine history_init_stage1()

    call fiona_init()

  end subroutine history_init_stage1

  subroutine history_init_stage2()

  end subroutine history_init_stage2

  subroutine history_init_stage3()

    character(10) time_value, time_units
    real(r8) seconds, months

    if (history_interval(1) == 'N/A') call log_error('Parameter history_interval is not set!')
    if (case_name == 'N/A') call log_error('Parameter case_name is not set!')

    time_value = split_string(history_interval(1), ' ', 1)
    time_units = split_string(history_interval(1), ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days', 'sol')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      if (proc%is_root()) call log_error('Invalid history interval ' // trim(history_interval(1)) // '!', __FILE__, __LINE__)
    end select

    call time_add_alert('history_write', seconds=seconds)

    if (trim(output_h0_new_file) == '') then
      call time_add_alert('h0_new_file', seconds=seconds)
      if (proc%is_root()) call log_notice('Output data every ' // trim(history_interval(1)) // '.')
    else if (output_h0_new_file == 'one_file') then
      if (proc%is_root()) call log_notice('Output data in one file.')
    else
      time_value = split_string(output_h0_new_file, ' ', 1)
      time_units = split_string(output_h0_new_file, ' ', 2)
      if (time_units == 'months') then
        read(time_value, *) months
        if (proc%is_root()) call log_notice('Output data every ' // trim(time_value) // ' months.')
        call time_add_alert('h0_new_file', months=months)
      else
        read(time_value, *) seconds
        select case (time_units)
        case ('days')
          if (proc%is_root()) call log_notice('Output data every ' // trim(time_value) // ' days.')
          seconds = seconds * 86400
        case ('sol')
          if (proc%is_root()) call log_notice('Output data every ' // trim(time_value) // ' sol.')
          seconds = seconds * 86400
        case ('hours')
          if (proc%is_root()) call log_notice('Output data every ' // trim(time_value) // ' hours.')
          seconds = seconds * 3600
        case ('minutes')
          if (proc%is_root()) call log_notice('Output data every ' // trim(time_value) // ' minutes.')
          seconds = seconds * 60
        case ('seconds')
          if (proc%is_root()) call log_notice('Output data every ' // trim(time_value) // ' seconds.')
          seconds = seconds
        case default
          if (proc%is_root()) call log_error('Invalid output_h0_new_file ' // trim(output_h0_new_file) // '!', __FILE__, __LINE__)
        end select
        call time_add_alert('h0_new_file', seconds=seconds)
      end if
    end if

    call fiona_set_time(time_units, start_time_str)

    call history_setup_h0()
    call history_setup_h1()
    call history_setup_h2()

  end subroutine history_init_stage3

  subroutine history_final()

    call fiona_final()

  end subroutine history_final

  subroutine history_setup_h0()

    class(*), pointer :: field
    integer i

    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name), &
      mpi_comm=proc%comm_io, ngroups=output_ngroups, async=use_async_io)
    ! Global attributes
    call fiona_add_att('h0', 'planet', planet)
    call namelist_add_atts('h0')
    ! Dimensions
    if (output_h0_new_file == history_interval(1)) then
      call fiona_add_dim('h0', 'time', size=1, add_var=.true.)
    else
      call fiona_add_dim('h0', 'time', add_var=.true.)
    end if
    call fiona_add_dim('h0', 'lon' , size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lat' , size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lev' , size=global_mesh%full_nlev)
    call fiona_add_dim('h0', 'ilon', size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilat', size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilev', size=global_mesh%half_nlev)
    ! Variables
    call fiona_add_var('h0', 'area', long_name='Cell area', units='m2', dim_names=['lat'], dtype=output_h0_dtype)
    select case (planet)
    case ('mars')
      call fiona_add_var('h0', 'Ls', long_name='Solar longitude', units='deg', dim_names=['time'], dtype=output_h0_dtype)
    end select
    if (.not. advection) then
      call fiona_add_var('h0', 'tm' , long_name='Total mass'               , units='', dim_names=['time'])
      call fiona_add_var('h0', 'te' , long_name='Total energy'             , units='', dim_names=['time'])
      call fiona_add_var('h0', 'tpe', long_name='Total potential enstrophy', units='', dim_names=['time'])
    end if

    call add_fields('h0', blocks(1)%dstate(1)   %fields)
    call add_fields('h0', blocks(1)%dtend       %fields)
    call add_fields('h0', blocks(1)%static      %fields, static=.true.)
    call add_fields('h0', blocks(1)%aux         %fields)
    call add_fields('h0', blocks(1)%adv_batch_bg%fields)
    call add_fields('h0', blocks(1)%adv_batch_pt%fields)
    call add_fields('h0', blocks(1)%adv_batch_nh%fields)
    if (allocated(blocks(1)%adv_batches)) then
      do i = 1, size(blocks(1)%adv_batches)
        call add_fields('h0', blocks(1)%adv_batches(i)%fields)
      end do
    end if
    field => tracers(1)%q
    call add_field ('h0', field)

    call physics_add_output('h0')

  end subroutine history_setup_h0

  subroutine history_write_h0(itime)

    integer, intent(in) :: itime

    logical, save :: first_call = .true.
    real(8) time1, time2
    integer iblk, i

    if (proc%is_root()) then
      call log_notice('Write h0 file.')
      time1 = MPI_WTIME()
    end if

    if (.not. time_has_alert('h0_new_file')) then
      call fiona_start_output('h0', dble(elapsed_seconds), new_file=first_call)
    else
      call fiona_start_output('h0', dble(elapsed_seconds), new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    end if

    call fiona_output('h0', 'area', global_mesh%area_cell(1:global_mesh%full_nlat))
    select case (planet)
    case ('mars')
      call fiona_output('h0', 'Ls', curr_time%solar_longitude() * deg)
    end select

    if (.not. advection) then
      call fiona_output('h0', 'tm' , blocks(1)%dstate(itime)%tm )
      call fiona_output('h0', 'te' , blocks(1)%dstate(itime)%te )
      call fiona_output('h0', 'tpe', blocks(1)%dstate(itime)%tpe)
    end if

    if (first_call .or. time_has_alert('h0_new_file')) then
      call fiona_output('h0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
      call fiona_output('h0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
      call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon))
      call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat))
      do iblk = 1, size(blocks)
        call write_fields('h0', blocks(iblk)%mesh, blocks(iblk)%static%fields)
      end do
    end if

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh         , &
                 dstate => blocks(iblk)%dstate(itime), &
                 dtend  => blocks(iblk)%dtend        , &
                 aux    => blocks(iblk)%aux          , &
                 q      => tracers(iblk)%q           )
      if (.not. use_div_damp .and. .not. advection) then
        call calc_div(blocks(iblk), dstate)
      end if
      call write_fields('h0', mesh, dstate                   %fields)
      call write_fields('h0', mesh, dtend                    %fields)
      call write_fields('h0', mesh, aux                      %fields)
      call write_fields('h0', mesh, blocks(iblk)%adv_batch_bg%fields)
      call write_fields('h0', mesh, blocks(iblk)%adv_batch_pt%fields)
      call write_fields('h0', mesh, blocks(iblk)%adv_batch_nh%fields)
      if (allocated(blocks(iblk)%adv_batches)) then
        do i = 1, size(blocks(iblk)%adv_batches)
          call write_fields('h0', mesh, blocks(iblk)%adv_batches(i)%fields)
        end do
      end if
      call write_field ('h0', mesh, q)
      end associate
      call physics_output('h0', iblk)
    end do

    call fiona_end_output('h0', keep_dataset=.true.)

    if (proc%is_root()) then
      time2 = MPI_WTIME()
      call log_notice('Done write h0 file cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

    first_call = .false.

  end subroutine history_write_h0

  subroutine history_setup_h1()

    integer i

    call fiona_create_dataset('h1', desc=case_desc, file_prefix=trim(case_name), &
      mpi_comm=proc%comm_io, ngroups=output_ngroups)
    call namelist_add_atts('h1')
    call fiona_add_dim('h1', 'time' , add_var=.true.)
    call fiona_add_dim('h1', 'lon'  , size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lat'  , size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lev'  , size=global_mesh%full_nlev)
    call fiona_add_dim('h1', 'ilon' , size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilat' , size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilev' , size=global_mesh%half_nlev)

    call add_fields('h1', blocks(1)%dstate(1)%fields)
    call add_fields('h1', blocks(1)%dtend    %fields)
    call add_fields('h1', blocks(1)%static   %fields)
    call add_fields('h1', blocks(1)%aux      %fields)
    call add_fields('h1', blocks(1)%adv_batch_bg%fields)
    call add_fields('h1', blocks(1)%adv_batch_pt%fields)
    call add_fields('h1', blocks(1)%adv_batch_nh%fields)
    if (physics_suite /= 'N/A' .or. physics_suite /= '') then
      do i = 1, ntracers
        call fiona_add_var('h1', 'd' // trim(tracer_names(i)) // '_phys', long_name='Physics tendency of ' // tracer_long_names(i), &
          units=tracer_units(i) // ' s-1', dim_names=cell_dims_3d, dtype=output_h0_dtype)
      end do
    end if

    ! Filter parameters
    call fiona_add_var('h1', 'fw', long_name='Filter width', units='', dim_names=['lat'])

  end subroutine history_setup_h1

  subroutine history_write_h1(itime)

    integer, intent(in) :: itime

    logical, save :: first_call = .true.
    real(8) time1, time2
    integer iblk

    if (proc%is_root()) then
      call log_notice('Write h1 file.')
      time1 = MPI_WTIME()
    end if

    if (.not. time_has_alert('h0_new_file')) then
      call fiona_start_output('h1', dble(elapsed_seconds), new_file=first_call)
    else
      call fiona_start_output('h1', dble(elapsed_seconds), new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    end if

    first_call = .false.

    call fiona_output('h1', 'lon' , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('h1', 'lat' , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('h1', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('h1', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat))

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh         , &
                 dstate => blocks(iblk)%dstate(itime), &
                 dtend  => blocks(iblk)%dtend        , &
                 static => blocks(iblk)%static       , &
                 aux    => blocks(iblk)%aux          )
      call write_fields('h1', mesh, dstate%fields)
      call write_fields('h1', mesh, dtend %fields)
      call write_fields('h1', mesh, static%fields)
      call write_fields('h1', mesh, aux   %fields)
      call write_fields('h1', mesh, blocks(iblk)%adv_batch_bg%fields)
      call write_fields('h1', mesh, blocks(iblk)%adv_batch_pt%fields)
      call write_fields('h1', mesh, blocks(iblk)%adv_batch_nh%fields)

      ! Filter parameters
      call fiona_output('h1', 'fw', blocks(iblk)%big_filter%width_lon(mesh%full_jds:mesh%full_jde), start=[mesh%full_jds], count=[mesh%full_nlat])
      end associate
    end do

    call fiona_end_output('h1', keep_dataset=.true.)

    if (proc%is_root()) then
      time2 = MPI_WTIME()
      call log_notice('Done write h1 file cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine history_write_h1

  subroutine history_setup_h2

    integer i

    if (.not. regrid_initialized) return

    call fiona_create_dataset('h2', desc=case_desc, file_prefix=trim(case_name), &
      mpi_comm=proc%comm_io, ngroups=output_ngroups)
    ! Global attributes
    call fiona_add_att('h2', 'planet', planet)
    ! Dimensions
    call fiona_add_dim('h2', 'time' , add_var=.true.)
    call fiona_add_dim('h2', 'lon'  , size=regrid_global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h2', 'lat'  , size=regrid_global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h2', 'lev'  , long_name='Pressure level', units='Pa', size=regrid_global_mesh%full_nlev, add_var=.true., decomp=.false.)

    do i = 1, regrids(1)%nfields
      call fiona_add_var('h2', regrids(1)%fields(i)%name, regrids(1)%fields(i)%long_name, regrids(1)%fields(i)%units, &
        dim_names=cell_dims_3d, dtype=output_h0_dtype, missing_value=real(inf))
    end do

  end subroutine history_setup_h2

  subroutine history_write_h2()

    integer is, ie, js, je, ks, ke, i
    integer start(3), count(3)
    real(8) time1, time2

    if (proc%is_root()) then
      call log_notice('Write h2 file.')
      time1 = MPI_WTIME()
    end if

    if (.not. time_has_alert('h0_new_file')) then
      call fiona_start_output('h2', dble(elapsed_seconds), new_file=time_step==0)
    else if (time_is_alerted('h0_new_file')) then
      call fiona_start_output('h2', dble(elapsed_seconds), new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    else
      call fiona_start_output('h2', dble(elapsed_seconds), new_file=.false.)
    end if
    call fiona_output('h2', 'lon', regrid_global_mesh%full_lon_deg(1:regrid_global_mesh%full_nlon))
    call fiona_output('h2', 'lat', regrid_global_mesh%full_lat_deg(1:regrid_global_mesh%full_nlat))
    call fiona_output('h2', 'lev', regrid_global_mesh%full_lev(1:regrid_global_mesh%full_nlev))

    is = regrids(1)%mesh%full_ids; ie = regrids(1)%mesh%full_ide
    js = regrids(1)%mesh%full_jds; je = regrids(1)%mesh%full_jde
    ks = regrids(1)%mesh%full_kds; ke = regrids(1)%mesh%full_kde
    start = [is,js,ks]
    count = [regrids(1)%mesh%full_nlon,regrids(1)%mesh%full_nlat,regrids(1)%mesh%full_nlev]
    do i = 1, regrids(1)%nfields
      call fiona_output('h2', regrids(1)%fields(i)%name, regrids(1)%fields(i)%d(is:ie,js:je,ks:ke), start=start, count=count)
    end do

    call fiona_end_output('h2', keep_dataset=.true.)

    if (proc%is_root()) then
      time2 = MPI_WTIME()
      call log_notice('Done write h2 file cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine history_write_h2

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
      if (field%output /= dtag) return
      n = merge(2, 3, present(static))
      select case (field%loc)
      case ('cell')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names=cell_dims_2d(:n), dtype=output_h0_dtype)
      case ('lon')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names= lon_dims_2d(:n), dtype=output_h0_dtype)
      case ('lat')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names= lat_dims_2d(:n), dtype=output_h0_dtype)
      case ('vtx')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names= vtx_dims_2d(:n), dtype=output_h0_dtype)
      end select
    type is (latlon_field3d_type)
      if (field%output /= dtag) return
      n = merge(3, 4, present(static))
      select case (field%loc)
      case ('cell')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names=   cell_dims_3d(:n), dtype=output_h0_dtype)
      case ('lon')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names=    lon_dims_3d(:n), dtype=output_h0_dtype)
      case ('lat')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names=    lat_dims_3d(:n), dtype=output_h0_dtype)
      case ('vtx')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names=    vtx_dims_3d(:n), dtype=output_h0_dtype)
      case ('lev')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names=    lev_dims_3d(:n), dtype=output_h0_dtype)
      case ('lev_lon')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names=lev_lon_dims_3d(:n), dtype=output_h0_dtype)
      case ('lev_lat')
        call fiona_add_var(dtag, field%name, long_name=field%long_name, units=field%units, dim_names=lev_lat_dims_3d(:n), dtype=output_h0_dtype)
      case default
        call log_error('Invalid field location ' // trim(field%loc) // '!', __FILE__, __LINE__)
      end select
    type is (latlon_field4d_type)
      if (field%output /= dtag) return
      do i = 1, field%dim4_size
        call fiona_add_var(dtag, trim(field%name) // '_' // trim(field%var4_names(i)), long_name=trim(field%long_name) // ' of ' // trim(field%var4_names(i)), &
          units=field%units, dim_names=cell_dims_3d, dtype=output_h0_dtype)
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
      if (field%output /= dtag) return
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
      if (field%output /= dtag) return
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
      if (field%output /= dtag) return
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

end module history_mod
