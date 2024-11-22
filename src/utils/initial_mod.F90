! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module initial_mod

  use fiona
  use string
  use flogger
  use datetime
  use const_mod
  use namelist_mod
  use formula_mod
  use time_mod
  use block_mod
  use tracer_mod
  use latlon_parallel_mod
  use operators_mod
  use prepare_mod
  use process_mod, only: proc

  implicit none

  private

  public initial_write
  public initial_read_init
  public initial_read

contains

  subroutine initial_write(initial_file, initial_time)

    character(*), intent(in) :: initial_file
    character(*), intent(in) :: initial_time

    character(4) cell_dims(4), cell_dims_2d(3)
    character(4) lon_dims(4)
    character(4) lat_dims(4)
    character(4) lev_dims(4)
    integer iblk, is, ie, js, je, ks, ke
    integer start(3), count(3)
    real(8) time1, time2

    cell_dims   (1) =  'lon';  cell_dims   (2) =  'lat';  cell_dims   (3) =  'lev';     cell_dims(4) = 'time'
     lon_dims   (1) = 'ilon';   lon_dims   (2) =  'lat';   lon_dims   (3) =  'lev';      lon_dims(4) = 'time'
     lat_dims   (1) =  'lon';   lat_dims   (2) = 'ilat';   lat_dims   (3) =  'lev';      lat_dims(4) = 'time'
     lev_dims   (1) =  'lon';   lev_dims   (2) =  'lat';   lev_dims   (3) = 'ilev';      lev_dims(4) = 'time'
    cell_dims_2d(1) =  'lon';  cell_dims_2d(2) =  'lat';  cell_dims_2d(3) = 'time'

    if (proc%is_root()) then
      call log_notice('Write ' // trim(initial_file) // '.')
      call cpu_time(time1)
    end if

    call fiona_create_dataset('i0', file_path=initial_file, start_time=initial_time, time_units='hours', mpi_comm=proc%comm_model, ngroups=output_ngroups)
    call fiona_add_dim('i0', 'time', add_var=.true.)
    call fiona_add_dim('i0',  'lon', size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('i0',  'lat', size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('i0', 'ilon', size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('i0', 'ilat', size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    if (baroclinic) then
      call fiona_add_dim('i0',  'lev', size=global_mesh%full_nlev, add_var=.true., decomp=.false.)
      call fiona_add_dim('i0', 'ilev', size=global_mesh%half_nlev, add_var=.true., decomp=.false.)
      call fiona_add_var('i0', 'pt'         , long_name='Potential temperature'                     , units='K'       , dtype=output_i0_dtype, dim_names=cell_dims)
      call fiona_add_var('i0', 'qv'         , long_name='Water vapor dry mixing ratio'              , units='kg kg-1' , dtype=output_i0_dtype, dim_names=cell_dims)
      call fiona_add_var('i0', 'mgs'        , long_name='Surface dry-air weight'                    , units='Pa'      , dtype=output_i0_dtype, dim_names=cell_dims_2d)
      call fiona_add_var('i0', 'u'          , long_name='U wind component'                          , units='m s-1'   , dtype=output_i0_dtype, dim_names=lon_dims)
      call fiona_add_var('i0', 'v'          , long_name='V wind component'                          , units='m s-1'   , dtype=output_i0_dtype, dim_names=lat_dims)
      call fiona_add_var('i0', 'z'          , long_name='Height'                                    , units='m'       , dtype=output_i0_dtype, dim_names=lev_dims)
      call fiona_add_var('i0', 'ref_ps_smth', long_name='Smoothed reference surface pressure'       , units='Pa'      , dtype=output_i0_dtype, dim_names=cell_dims_2d)
      call fiona_add_var('i0', 'ref_ps_perb', long_name='Perturbed reference surface pressure'      , units='Pa'      , dtype=output_i0_dtype, dim_names=cell_dims_2d)
    end if
    call fiona_add_var('i0', 'zs'           , long_name='Surface height'                            , units='m'       , dtype=output_i0_dtype, dim_names=cell_dims_2d)
    call fiona_add_var('i0', 'zs_std'       , long_name='Surface height subgrid standard deviation' , units='m2'      , dtype=output_i0_dtype, dim_names=cell_dims_2d)
    call fiona_add_var('i0', 'dzsdx'        , long_name='Zonal zs gradient'                         , units=''        , dtype=output_i0_dtype, dim_names=lon_dims(1:2))
    call fiona_add_var('i0', 'dzsdy'        , long_name='Meridional zs gradient'                    , units=''        , dtype=output_i0_dtype, dim_names=lat_dims(1:2))
    call fiona_add_var('i0', 'landmask'     , long_name='Land mask'                                 , units='1'       , dtype=output_i0_dtype, dim_names=cell_dims_2d)

    call fiona_start_output('i0', 0.0d0)
    call fiona_output('i0',  'lon', global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('i0',  'lat', global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('i0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('i0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat))
    if (baroclinic) then
      call fiona_output('i0',  'lev', global_mesh%full_lev(1:global_mesh%full_nlev))
      call fiona_output('i0', 'ilev', global_mesh%half_lev(1:global_mesh%half_nlev))
    end if

    do iblk = 1, size(blocks)
      associate (mesh     => blocks(iblk)%mesh            , &
                 gzs      => blocks(iblk)%static%gzs      , &
                 zs_std   => blocks(iblk)%static%zs_std   , &
                 landmask => blocks(iblk)%static%landmask , &
                 dzsdx    => blocks(iblk)%static%dzsdx    , &
                 dzsdy    => blocks(iblk)%static%dzsdy    , &
                 mgs      => blocks(iblk)%dstate(1)%mgs   , &
                 u_lon    => blocks(iblk)%dstate(1)%u_lon , &
                 v_lat    => blocks(iblk)%dstate(1)%v_lat , &
                 pt       => blocks(iblk)%dstate(1)%pt    , &
                 t        => blocks(iblk)%aux%t           , &
                 gz_lev   => blocks(iblk)%dstate(1)%gz_lev, &
                 q        => tracers(iblk)%q              , &
                 qm       => tracers(iblk)%qm             )
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]

      if (baroclinic) then
        call fiona_output('i0', 'pt', pt%d(is:ie,js:je,ks:ke), start=start, count=count)
        if (idx_qv > 0) then
          call fiona_output('i0', 'qv', q%d(is:ie,js:je,ks:ke,idx_qv), start=start, count=count)
        end if
        if (idx_qc > 0 .and. fiona_has_var('i0', 'qc')) then
          call fiona_output('i0', 'qc', q%d(is:ie,js:je,ks:ke,idx_qc), start=start, count=count)
        end if
        if (idx_qi > 0) then
          call fiona_output('i0', 'qi', q%d(is:ie,js:je,ks:ke,idx_qi), start=start, count=count)
        end if
        if (idx_nc > 0) then
          call fiona_output('i0', 'nc', q%d(is:ie,js:je,ks:ke,idx_nc), start=start, count=count)
        end if
        if (idx_ni > 0) then
          call fiona_output('i0', 'ni', q%d(is:ie,js:je,ks:ke,idx_ni), start=start, count=count)
        end if
        call fiona_output('i0', 'mgs', mgs%d(is:ie,js:je), start=start, count=count)
      end if
      call fiona_output('i0', 'zs'      , gzs     %d(is:ie,js:je) / g, start=start, count=count)
      call fiona_output('i0', 'zs_std'  , zs_std  %d(is:ie,js:je)    , start=start, count=count)
      call fiona_output('i0', 'landmask', landmask%d(is:ie,js:je)    , start=start, count=count)

      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]

      call fiona_output('i0', 'dzsdx', dzsdx%d(is:ie,js:je)      , start=start, count=count)
      call fiona_output('i0', 'u_lon', u_lon%d(is:ie,js:je,ks:ke), start=start, count=count)

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%half_jds; je = mesh%half_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]

      call fiona_output('i0', 'dzsdy', dzsdy%d(is:ie,js:je)      , start=start, count=count)
      call fiona_output('i0', 'v_lat', v_lat%d(is:ie,js:je,ks:ke), start=start, count=count)

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%half_kds; ke = mesh%half_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%half_nlev]

      call fiona_output('i0', 'z', gz_lev%d(is:ie,js:je,ks:ke) / g, start=start, count=count)
      end associate
    end do
    call fiona_end_output('i0')

    if (proc%is_root()) then
      call cpu_time(time2)
      call log_notice('Done write initial data cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine initial_write

  subroutine initial_read_init(initial_file_)

    character(*), intent(in), optional :: initial_file_

    if (present(initial_file_)) then
      call fiona_open_dataset('i0', file_path=initial_file_, mpi_comm=proc%comm_model, ngroups=input_ngroups)
      if (proc%is_root()) call log_notice('Read initial data from ' // trim(initial_file_) // '.')
    else
      call fiona_open_dataset('i0', file_path=initial_file, mpi_comm=proc%comm_model, ngroups=input_ngroups)
      if (proc%is_root()) call log_notice('Read initial data from ' // trim(initial_file) // '.')
    end if

    call fiona_start_input('i0')

    if (fiona_has_var('i0', 'qv')) then
      if (proc%is_root()) call log_notice('Add tracer qv.')
      call tracer_add('moist', dt_adv, 'qv', 'Water vapor', 'kg kg-1')
    end if
    if (fiona_has_var('i0', 'qc')) then
      if (proc%is_root()) call log_notice('Add tracer qc.')
      call tracer_add('moist', dt_adv, 'qc', 'Cloud water', 'kg kg-1')
    end if
    if (fiona_has_var('i0', 'qi')) then
      if (proc%is_root()) call log_notice('Add tracer qi.')
      call tracer_add('moist', dt_adv, 'qi', 'Cloud ice'  , 'kg kg-1')
    end if
    if (fiona_has_var('i0', 'nc')) then
      if (proc%is_root()) call log_notice('Add tracer nc.')
      call tracer_add('moist', dt_adv, 'nc', 'Cloud number concentration', 'kg-1')
    end if
    if (fiona_has_var('i0', 'ni')) then
      if (proc%is_root()) call log_notice('Add tracer ni.')
      call tracer_add('moist', dt_adv, 'ni', 'Cloud ice number concentration', 'kg-1')
    end if
    if (fiona_has_var('i0', 'qr')) then
      if (proc%is_root()) call log_notice('Add tracer qr.')
      call tracer_add('moist', dt_adv, 'qr', 'Rain water' , 'kg kg-1')
    end if
    if (fiona_has_var('i0', 'qs')) then
      if (proc%is_root()) call log_notice('Add tracer qs.')
      call tracer_add('moist', dt_adv, 'qs', 'Snow water' , 'kg kg-1')
    end if

  end subroutine initial_read_init

  subroutine initial_read(initial_file_)

    character(*), intent(in), optional :: initial_file_

    integer iblk, is, ie, js, je, ks, ke, i, j, k
    integer start(3), count(3)
    real(8) time1, time2
    logical :: input_zs       = .false.
    logical :: input_mgs      = .false.
    logical :: input_pt       = .false.
    logical :: input_t        = .false.
    logical :: input_qv       = .false.
    logical :: input_qc       = .false.
    logical :: input_qi       = .false.
    logical :: input_nc       = .false.
    logical :: input_ni       = .false.
    logical :: input_u_lon    = .false.
    logical :: input_u        = .false.
    logical :: input_v_lat    = .false.
    logical :: input_v        = .false.
    logical :: input_gz_lev   = .false.

    if (proc%is_root()) call cpu_time(time1)
    
    do iblk = 1, size(blocks)
      associate (block    => blocks(iblk)                 , &
                 mesh     => blocks(iblk)%mesh            , &
                 dstate   => blocks(iblk)%dstate(1)       , &
                 gzs      => blocks(iblk)%static%gzs      , &
                 zs_std   => blocks(iblk)%static%zs_std   , &
                 landmask => blocks(iblk)%static%landmask , &
                 dzsdx    => blocks(iblk)%static%dzsdx    , &
                 dzsdy    => blocks(iblk)%static%dzsdy    , &
                 mgs      => blocks(iblk)%dstate(1)%mgs   , &
                 phs      => blocks(iblk)%dstate(1)%phs   , &
                 ph       => blocks(iblk)%dstate(1)%ph    , &
                 ph_lev   => blocks(iblk)%dstate(1)%ph_lev, &
                 u_lon    => blocks(iblk)%dstate(1)%u_lon , &
                 v_lat    => blocks(iblk)%dstate(1)%v_lat , &
                 u        => blocks(iblk)%aux%u           , &
                 v        => blocks(iblk)%aux%v           , &
                 pt       => blocks(iblk)%dstate(1)%pt    , &
                 t        => blocks(iblk)%aux%t           , &
                 gz       => blocks(iblk)%dstate(1)%gz    , &
                 gz_lev   => blocks(iblk)%dstate(1)%gz_lev, &
                 q        => tracers(iblk)%q              , &
                 qm       => tracers(iblk)%qm             )
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]

      if (fiona_has_var('i0', 'zs')) then
        call fiona_input('i0', 'zs', gzs%d(is:ie,js:je), start=start, count=count); gzs%d = gzs%d * g
        call fill_halo(gzs)
        input_zs = .true.
      end if
      if (fiona_has_var('i0', 'zs_std')) then
        call fiona_input('i0', 'zs_std', zs_std%d(is:ie,js:je), start=start, count=count)
        call fill_halo(zs_std)
      end if
      if (fiona_has_var('i0', 'landmask')) then
        call fiona_input('i0', 'landmask', landmask%d(is:ie,js:je), start=start, count=count)
        call fill_halo(landmask)
      end if
      if (baroclinic) then
        if (fiona_has_var('i0', 'mgs')) then
          call fiona_input('i0', 'mgs', mgs%d(is:ie,js:je), start=start, count=count)
          call fill_halo(mgs)
          input_mgs = .true.
        end if
        if (fiona_has_var('i0', 'pt')) then
          call fiona_input('i0', 'pt', pt%d(is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(pt)
          input_pt = .true.
        else if (fiona_has_var('i0', 't')) then
          call fiona_input('i0', 't', t%d(is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(t)
          input_t = .true.
        end if
        if (idx_qv > 0 .and. fiona_has_var('i0', 'qv')) then
          call fiona_input('i0', 'qv', q%d(is:ie,js:je,ks:ke,idx_qv), start=start, count=count)
          call fill_halo(q, idx_qv)
          input_qv = .true.
        end if
        if (fiona_has_var('i0', 'u')) then
          call fiona_input('i0', 'u', u%d(is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(u)
          input_u = .true.
        end if
        if (fiona_has_var('i0', 'v')) then
          call fiona_input('i0', 'v', v%d(is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(v)
          input_v = .true.
        end if
      else
        call fiona_input('i0', 'z' , gz%d(is:ie,js:je,ks:ke), start=start, count=count); gz%d = gz%d * g
        call fill_halo(gz)
      end if

      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]

      if (fiona_has_var('i0', 'u_lon')) then
        call fiona_input('i0', 'u_lon', u_lon%d(is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(u_lon)
        input_u_lon = .true.
      end if

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%half_jds; je = mesh%half_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]

      if (fiona_has_var('i0', 'v_lat')) then
        call fiona_input('i0', 'v_lat', v_lat%d(is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(v_lat)
        input_v_lat = .true.
      end if

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%half_kds; ke = mesh%half_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%half_nlev]

      if (nonhydrostatic) then
        if (fiona_has_var('i0', 'gz_lev')) then
          call fiona_input('i0', 'z', gz_lev%d(is:ie,js:je,ks:ke), start=start, count=count); gz_lev%d = gz_lev%d * g
          call fill_halo(gz_lev)
          input_gz_lev = .true.
        end if
      end if

      ! Calculate horizontal wind components on the edges.
      if (input_u .and. input_v .and. .not. (input_u_lon .and. input_v_lat)) then
        call dstate%a2c(u, v)
        call fill_halo(u_lon)
        call fill_halo(v_lat)
      end if
      if (input_mgs) then
        call calc_mg (block, dstate)
        call calc_dmg(block, dstate)
      else
        stop 111
      end if
      ! Convert wet mixing ratio to dry mixing ratio.
      call tracer_calc_qm(block)
      if (input_qv) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              q%d(i,j,k,idx_qv) = q%d(i,j,k,idx_qv) / (1 - qm%d(i,j,k))
            end do
          end do
        end do
        call fill_halo(q, idx_qv)
      end if
      if (input_qc) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              q%d(i,j,k,idx_qc) = q%d(i,j,k,idx_qc) / (1 - qm%d(i,j,k))
            end do
          end do
        end do
        call fill_halo(q, idx_qc)
      end if
      if (input_qi) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              q%d(i,j,k,idx_qi) = q%d(i,j,k,idx_qi) / (1 - qm%d(i,j,k))
            end do
          end do
        end do
        call fill_halo(q, idx_qi)
      end if
      if (input_nc) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              q%d(i,j,k,idx_nc) = q%d(i,j,k,idx_nc) / (1 - qm%d(i,j,k))
            end do
          end do
        end do
        call fill_halo(q, idx_nc)
      end if
      if (input_ni) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              q%d(i,j,k,idx_ni) = q%d(i,j,k,idx_ni) / (1 - qm%d(i,j,k))
            end do
          end do
        end do
        call fill_halo(q, idx_ni)
      end if
      call tracer_calc_qm(block)
      call calc_ph(block, dstate)
      ! Calculate modified potential temperature.
      if (input_t .and. .not. input_pt) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              pt%d(i,j,k) = modified_potential_temperature(t%d(i,j,k), ph%d(i,j,k), q%d(i,j,k,idx_qv))
            end do
          end do
        end do
        call fill_halo(pt)
      end if
      end associate
    end do
    call fiona_end_input('i0')

    ! Prepare topography if not read.
    if (.not. input_zs) then
      call prepare_topo()
    end if

    if (proc%is_root()) then
      call cpu_time(time2)
      call log_notice('Done read initial file cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine initial_read

end module initial_mod
