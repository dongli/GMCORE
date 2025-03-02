! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module latlon_bkg_mod

  use flogger
  use namelist_mod
  use const_mod
  use block_mod
  use vert_coord_mod
  use formula_mod
  use process_mod
  use latlon_parallel_mod
  use era5_reader_mod
  use cam_reader_mod
  use openmars_reader_mod
  use latlon_interp_mod
  use vert_interp_mod
  use tracer_mod
  use operators_mod

  implicit none

  private

  public latlon_bkg_read
  public latlon_bkg_final
  public latlon_bkg_regrid_phs
  public latlon_bkg_calc_ph
  public latlon_bkg_regrid_wet_qv
  public latlon_bkg_regrid_wet_qc
  public latlon_bkg_regrid_wet_qi
  public latlon_bkg_regrid_wet_nc
  public latlon_bkg_regrid_wet_ni
  public latlon_bkg_regrid_wet_qr
  public latlon_bkg_regrid_wet_qs
  public latlon_bkg_calc_dry_qv
  public latlon_bkg_calc_dry_qc
  public latlon_bkg_calc_dry_qi
  public latlon_bkg_calc_dry_nc
  public latlon_bkg_calc_dry_ni
  public latlon_bkg_calc_dry_qr
  public latlon_bkg_calc_dry_qs
  public latlon_bkg_regrid_gz
  public latlon_bkg_calc_t
  public latlon_bkg_regrid_t
  public latlon_bkg_regrid_u
  public latlon_bkg_regrid_v
  public latlon_bkg_calc_mgs
  public latlon_bkg_calc_mg
  public latlon_bkg_calc_tv
  public latlon_bkg_calc_pt

contains

  subroutine latlon_bkg_read(min_lon, max_lon, min_lat, max_lat)

    real(r8), intent(in) :: min_lon, max_lon, min_lat, max_lat

    select case (bkg_type)
    case ('era5')
      call era5_reader_run(bkg_file, min_lon, max_lon, min_lat, max_lat)
    case ('cam')
      call cam_reader_run(bkg_file, min_lon, max_lon, min_lat, max_lat)
    case ('openmars')
      call openmars_reader_run(bkg_file)
    case default
      if (proc%is_root()) call log_error('Unknown bkg_type ' // trim(bkg_type) // '!')
    end select

  end subroutine latlon_bkg_read

  subroutine latlon_bkg_final()

    call era5_reader_final()
    call cam_reader_final()
    call openmars_reader_final()

  end subroutine latlon_bkg_final

  subroutine latlon_bkg_regrid_phs()

    real(r8), allocatable, dimension(:,:) :: ps0, zs0
    real(r8), allocatable, dimension(:,:,:) :: t0, p0
    integer iblk, i, j, k, k0

    logical do_hydrostatic_correct

    if (proc%is_root()) call log_notice('Regrid surface hydrostatic pressure.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh         , &
                 gzs  => blocks(iblk)%static%gzs   , &
                 phs  => blocks(iblk)%dstate(1)%phs)
      allocate(ps0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))
      allocate(zs0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))

      select case (bkg_type)
      case ('era5')
        allocate(t0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        allocate(p0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_ps, mesh, ps0)
        call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_zs, mesh, zs0)
        do k = 1, era5_nlev
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_tv(:,:,k), mesh, t0(:,:,k))
          p0(:,:,k) = era5_lev(k)
        end do
        do_hydrostatic_correct = .true.
      case ('cam')
        call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_ps, mesh, phs%d)
        do_hydrostatic_correct = .false.
      case ('openmars')
        call latlon_interp_bilinear_cell(openmars_lon, openmars_lat, openmars_ps, mesh, phs%d)
        do_hydrostatic_correct = .false.
      end select
      ! According to pressure-height formula based on hydrostatic assumption.
      if (do_hydrostatic_correct) then
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            do k0 = era5_nlev, 1, -1
              if (p0(i,j,k0) <= ps0(i,j)) exit
            end do
            phs%d(i,j) = ps0(i,j) * (1.0_r8 + lapse_rate * (gzs%d(i,j) / g - zs0(i,j)) / t0(i,j,k0))**(-1.0_r8 / lapse_rate / rd_o_g)
          end do
        end do
      end if
      call fill_halo(phs)
      if (allocated(t0)) deallocate(t0)
      if (allocated(p0)) deallocate(p0)
      deallocate(ps0, zs0)
      end associate
    end do

  end subroutine latlon_bkg_regrid_phs

  subroutine latlon_bkg_calc_ph()

    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Calculate hydrostatic pressure on each level.')

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh            , &
                 phs    => blocks(iblk)%dstate(1)%phs   , & ! in
                 ph     => blocks(iblk)%dstate(1)%ph    , & ! out
                 ph_lev => blocks(iblk)%dstate(1)%ph_lev)   ! out
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ph%d(i,j,k) = vert_coord_calc_mg(k, phs%d(i,j))
          end do
        end do
      end do
      do k = mesh%half_kds, mesh%half_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ph_lev%d(i,j,k) = vert_coord_calc_mg_lev(k, phs%d(i,j))
          end do
        end do
      end do
      end associate
    end do

  end subroutine latlon_bkg_calc_ph

  subroutine latlon_bkg_regrid_gz()

    real(r8), allocatable, dimension(:,:,:) :: z0, p0
    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Regrid geopotential height.')

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh            , &
                 gzs    => blocks(iblk)%static%gzs      , & ! in
                 ph_lev => blocks(iblk)%dstate(1)%ph_lev, & ! in
                 gz_lev => blocks(iblk)%dstate(1)%gz_lev, & ! out
                 gz     => blocks(iblk)%dstate(1)%gz    )   ! out
      select case (bkg_type)
      case ('era5')
        allocate(z0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        allocate(p0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_z(:,:,k), mesh, z0(:,:,k))
          p0(:,:,k) = era5_lev(k)
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_log_linear(p0(i,j,:), z0(i,j,:), ph_lev%d(i,j,1:mesh%full_nlev), gz_lev%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
            gz_lev%d(i,j,1:mesh%full_nlev) = gz_lev%d(i,j,1:mesh%full_nlev) * g
            gz_lev%d(i,j,mesh%half_kde) = gzs%d(i,j)
          end do
        end do
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              gz%d(i,j,k) = 0.5_r8 * (gz_lev%d(i,j,k) + gz_lev%d(i,j,k+1))
            end do
          end do
        end do
        deallocate(z0, p0)
      case default
        if (proc%is_root()) call log_error('Unsupported bkg_type ' // trim(bkg_type) // '!', __FILE__, __LINE__)
      end select
      call fill_halo(gz_lev)
      call fill_halo(gz)
      end associate
    end do

  end subroutine latlon_bkg_regrid_gz

  subroutine latlon_bkg_calc_t()

    real(r8) tv
    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Calculate temperature.')

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh            , &
                 ph_lev => blocks(iblk)%dstate(1)%ph_lev, & ! in
                 gz_lev => blocks(iblk)%dstate(1)%gz_lev, & ! in
                 q      => tracers(iblk)%q              , & ! in
                 qm     => tracers(iblk)%qm             , & ! in
                 t      => blocks(iblk)%aux%t           )   ! out
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            tv = (gz_lev%d(i,j,k) - gz_lev%d(i,j,k+1)) / rd / log(ph_lev%d(i,j,k+1) / ph_lev%d(i,j,k))
            if (idx_qv == 0) then
              t%d(i,j,k) = tv
            else
              t%d(i,j,k) = temperature_from_virtual_temperature_and_dry_mixing_ratio(tv, q%d(i,j,k,idx_qv), qm%d(i,j,k))
            end if
          end do
        end do
      end do
      end associate
    end do

  end subroutine latlon_bkg_calc_t

  subroutine latlon_bkg_regrid_t()

    real(r8), allocatable, dimension(:,:,:) :: t0, p0
    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Regrid temperature.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh        , &
                 ph   => blocks(iblk)%dstate(1)%ph, & ! in
                 t    => blocks(iblk)%aux%t       )   ! out
        select case (bkg_type)
        case ('era5')
          allocate(t0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
          do k = 1, era5_nlev
            call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_t(:,:,k), mesh, t0(:,:,k))
          end do
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              call vert_interp_log_linear(era5_lev(:), t0(i,j,:), ph%d(i,j,1:mesh%full_nlev), t%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
            end do
          end do
          deallocate(t0)
        case ('cam')
          allocate(t0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
          allocate(p0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
          do k = 1, cam_nlev
            call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_t(:,:,k), mesh, t0(:,:,k))
            call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_p(:,:,k), mesh, p0(:,:,k))
          end do
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              call vert_interp_log_linear(p0(i,j,:), t0(i,j,:), ph%d(i,j,1:mesh%full_nlev), t%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
            end do
          end do
          deallocate(t0, p0)
        case ('openmars')
          allocate(t0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,openmars_nlev))
          allocate(p0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,openmars_nlev))
          do k = 1, openmars_nlev
            call latlon_interp_bilinear_cell(openmars_lon, openmars_lat, openmars_t(:,:,k), mesh, t0(:,:,k))
            call latlon_interp_bilinear_cell(openmars_lon, openmars_lat, openmars_p(:,:,k), mesh, p0(:,:,k))
          end do
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              call vert_interp_log_linear(p0(i,j,:), t0(i,j,:), ph%d(i,j,1:mesh%full_nlev), t%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
            end do
          end do
          deallocate(t0, p0)
        end select
        call fill_halo(t)
      end associate
    end do

  end subroutine latlon_bkg_regrid_t

  subroutine latlon_bkg_regrid_u()

    real(r8), allocatable, dimension(:,:,:) :: u0, p0
    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Regrid u wind component.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh        , &
                 ph   => blocks(iblk)%dstate(1)%ph, &
                 u    => blocks(iblk)%dstate(1)%u_lon)
      select case (bkg_type)
      case ('era5')
        allocate(u0(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_lon_edge(era5_lon, era5_lat, era5_u(:,:,k), mesh, u0(:,:,k), extrap=.true., zero_pole=.true.)
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
            call vert_interp_linear(era5_lev, u0(i,j,:), ph%d(i,j,1:mesh%full_nlev), u%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(u0)
      case ('cam')
        allocate(u0(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        allocate(p0(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_lon_edge(cam_lon, cam_slat, cam_us(:,:,k), mesh, u0(:,:,k), zero_pole=.true.)
          call latlon_interp_bilinear_lon_edge(cam_lon, cam_lat , cam_p (:,:,k), mesh, p0(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
            call vert_interp_linear(p0(i,j,:), u0(i,j,:), ph%d(i,j,1:mesh%full_nlev), u%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(u0, p0)
      case ('openmars')
        allocate(u0(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,openmars_nlev))
        allocate(p0(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,openmars_nlev))
        do k = 1, openmars_nlev
          call latlon_interp_bilinear_lon_edge(openmars_lon, openmars_lat, openmars_u(:,:,k), mesh, u0(:,:,k), zero_pole=.true.)
          call latlon_interp_bilinear_lon_edge(openmars_lon, openmars_lat, openmars_p(:,:,k), mesh, p0(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
            call vert_interp_linear(p0(i,j,:), u0(i,j,:), ph%d(i,j,1:mesh%full_nlev), u%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(u0, p0)
      end select
      call fill_halo(u)
      end associate
    end do

  end subroutine latlon_bkg_regrid_u

  subroutine latlon_bkg_regrid_v()

    real(r8), allocatable, dimension(:,:,:) :: v0, p0
    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Regrid v wind component.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh        , &
                 ph   => blocks(iblk)%dstate(1)%ph, &
                 v    => blocks(iblk)%dstate(1)%v_lat)
      select case (bkg_type)
      case ('era5')
        allocate(v0(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_lat_edge(era5_lon, era5_lat, era5_v(:,:,k), mesh, v0(:,:,k), extrap=.true.)
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(era5_lev, v0(i,j,:), ph%d(i,j,1:mesh%full_nlev), v%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(v0)
      case ('cam')
        allocate(v0(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,cam_nlev))
        allocate(p0(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_lat_edge(cam_slon, cam_lat, cam_vs(:,:,k), mesh, v0(:,:,k))
          call latlon_interp_bilinear_lat_edge(cam_lon , cam_lat, cam_p (:,:,k), mesh, p0(:,:,k))
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p0(i,j,:), v0(i,j,:), ph%d(i,j,1:mesh%full_nlev), v%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(v0, p0)
      case ('openmars')
        allocate(v0(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,openmars_nlev))
        allocate(p0(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,openmars_nlev))
        do k = 1, openmars_nlev
          call latlon_interp_bilinear_lat_edge(openmars_lon, openmars_lat, openmars_v(:,:,k), mesh, v0(:,:,k))
          call latlon_interp_bilinear_lat_edge(openmars_lon, openmars_lat, openmars_p(:,:,k), mesh, p0(:,:,k))
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p0(i,j,:), v0(i,j,:), ph%d(i,j,1:mesh%full_nlev), v%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(v0, p0)
      end select
      call fill_halo(v)
      end associate
    end do

  end subroutine latlon_bkg_regrid_v

  subroutine latlon_bkg_regrid_wet_qv()

    real(r8), allocatable, dimension(:,:,:) :: q0, p0
    integer iblk, i, j, k

    if (idx_qv == 0) return

    if (proc%is_root()) call log_notice('Regrid water vapor wet mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%filter_mesh , &
                 ph   => blocks(iblk)%dstate(1)%ph, & ! in
                 q    => tracers(iblk)%q          )   ! out
      select case (bkg_type)
      case ('era5')
        allocate(q0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_qv(:,:,k), mesh, q0(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(era5_lev, q0(i,j,:), ph%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qv), allow_extrap=.true.)
          end do
        end do
        deallocate(q0)
      case ('cam')
        allocate(q0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        allocate(p0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_q(:,:,k), mesh, q0(:,:,k))
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_p(:,:,k), mesh, p0(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p0(i,j,:), q0(i,j,:), ph%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qv), allow_extrap=.true.)
          end do
        end do
        deallocate(q0, p0)
      end select
      where (q%d(:,:,:,idx_qv) < 0) q%d(:,:,:,idx_qv) = 0
      call fill_halo(q, idx_qv)
      end associate
    end do

  end subroutine latlon_bkg_regrid_wet_qv

  subroutine latlon_bkg_regrid_wet_qc()

    real(r8), allocatable, dimension(:,:,:) :: q0, p0
    integer iblk, i, j, k

    if (idx_qc == 0) return

    if (proc%is_root()) call log_notice('Regrid cloud liquid wet mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%filter_mesh , &
                 ph   => blocks(iblk)%dstate(1)%ph, & ! in
                 q    => tracers(iblk)%q          )   ! out
      select case (bkg_type)
      case ('era5')
        allocate(q0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_qc(:,:,k), mesh, q0(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(era5_lev, q0(i,j,:), ph%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qc), allow_extrap=.true.)
          end do
        end do
        deallocate(q0)
      case ('cam')
        allocate(q0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        allocate(p0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_cldliq(:,:,k), mesh, q0(:,:,k))
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_p(:,:,k), mesh, p0(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p0(i,j,:), q0(i,j,:), ph%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qc), allow_extrap=.true.)
          end do
        end do
        deallocate(q0, p0)
      end select
      where (q%d(:,:,:,idx_qc) < 0) q%d(:,:,:,idx_qc) = 0
      call fill_halo(q, idx_qc)
      end associate
    end do

  end subroutine latlon_bkg_regrid_wet_qc

  subroutine latlon_bkg_regrid_wet_qi()

    real(r8), allocatable, dimension(:,:,:) :: q0, p0
    integer iblk, i, j, k

    if (idx_qi == 0) return

    if (proc%is_root()) call log_notice('Regrid cloud ice wet mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%filter_mesh , &
                 ph   => blocks(iblk)%dstate(1)%ph, & ! in
                 q    => tracers(iblk)%q          )   ! out
      select case (bkg_type)
      case ('era5')
        allocate(q0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_qi(:,:,k), mesh, q0(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(era5_lev, q0(i,j,:), ph%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qi), allow_extrap=.true.)
          end do
        end do
        deallocate(q0)
      case ('cam')
        allocate(q0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        allocate(p0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_cldice(:,:,k), mesh, q0(:,:,k))
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_p(:,:,k), mesh, p0(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p0(i,j,:), q0(i,j,:), ph%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qi), allow_extrap=.true.)
          end do
        end do
        deallocate(q0, p0)
      end select
      where (q%d(:,:,:,idx_qi) < 0) q%d(:,:,:,idx_qi) = 0
      call fill_halo(q, idx_qi)
      end associate
    end do

  end subroutine latlon_bkg_regrid_wet_qi

  subroutine latlon_bkg_regrid_wet_nc()

    real(r8), allocatable, dimension(:,:,:) :: n0, p0
    integer iblk, i, j, k

    if (idx_nc == 0) return

    if (proc%is_root()) call log_notice('Regrid cloud liquid number wet mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%filter_mesh , &
                 ph   => blocks(iblk)%dstate(1)%ph, & ! in
                 q    => tracers(iblk)%q          )   ! out
      select case (bkg_type)
      case ('cam')
        allocate(n0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        allocate(p0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_numliq(:,:,k), mesh, n0(:,:,k))
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_p(:,:,k), mesh, p0(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p0(i,j,:), n0(i,j,:), ph%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_nc), allow_extrap=.true.)
          end do
        end do
        deallocate(n0, p0)
      end select
      where (q%d(:,:,:,idx_nc) < 0) q%d(:,:,:,idx_nc) = 0
      call fill_halo(q, idx_nc)
      end associate
    end do

  end subroutine latlon_bkg_regrid_wet_nc

  subroutine latlon_bkg_regrid_wet_ni()

    real(r8), allocatable, dimension(:,:,:) :: n0, p0
    integer iblk, i, j, k

    if (idx_ni == 0) return

    if (proc%is_root()) call log_notice('Regrid cloud ice number wet mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%filter_mesh , &
                 ph   => blocks(iblk)%dstate(1)%ph, & ! in
                 q    => tracers(iblk)%q          )   ! out
      select case (bkg_type)
      case ('cam')
        allocate(n0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        allocate(p0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_numice(:,:,k), mesh, n0(:,:,k))
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_p(:,:,k), mesh, p0(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p0(i,j,:), n0(i,j,:), ph%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_ni), allow_extrap=.true.)
          end do
        end do
        deallocate(n0, p0)
      end select
      where (q%d(:,:,:,idx_ni) < 0) q%d(:,:,:,idx_ni) = 0
      call fill_halo(q, idx_ni)
      end associate
    end do

  end subroutine latlon_bkg_regrid_wet_ni

  subroutine latlon_bkg_regrid_wet_qr()

    real(r8), allocatable, dimension(:,:,:) :: q0, p0
    integer iblk, i, j, k

    if (idx_qr == 0) return

    if (proc%is_root()) call log_notice('Regrid rain water wet mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%filter_mesh , &
                 ph   => blocks(iblk)%dstate(1)%ph, & ! in
                 q    => tracers(iblk)%q          )   ! out
      select case (bkg_type)
      case ('era5')
        allocate(q0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_qr(:,:,k), mesh, q0(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(era5_lev, q0(i,j,:), ph%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qr), allow_extrap=.true.)
          end do
        end do
        deallocate(q0)
      end select
      where (q%d(:,:,:,idx_qr) < 0) q%d(:,:,:,idx_qr) = 0
      call fill_halo(q, idx_qr)
      end associate
    end do

  end subroutine latlon_bkg_regrid_wet_qr

  subroutine latlon_bkg_regrid_wet_qs()

    real(r8), allocatable, dimension(:,:,:) :: q0, p0
    integer iblk, i, j, k

    if (idx_qs == 0) return

    if (proc%is_root()) call log_notice('Regrid snow water wet mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%filter_mesh , &
                 ph   => blocks(iblk)%dstate(1)%ph, & ! in
                 q    => tracers(iblk)%q          )   ! out
      select case (bkg_type)
      case ('era5')
        allocate(q0(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_qs(:,:,k), mesh, q0(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(era5_lev, q0(i,j,:), ph%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qs), allow_extrap=.true.)
          end do
        end do
        deallocate(q0)
      end select
      where (q%d(:,:,:,idx_qs) < 0) q%d(:,:,:,idx_qs) = 0
      call fill_halo(q, idx_qs)
      end associate
    end do

  end subroutine latlon_bkg_regrid_wet_qs

  subroutine latlon_bkg_calc_mgs()

    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Calculate surface dry-air pressure.')

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh            , &
                 phs    => blocks(iblk)%dstate(1)%phs   , & ! in
                 ph_lev => blocks(iblk)%dstate(1)%ph_lev, & ! in
                 qm     => tracers(iblk)%qm             , & ! in
                 mgs    => blocks(iblk)%dstate(1)%mgs   )   ! out
      call mgs%copy(phs)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            mgs%d(i,j) = mgs%d(i,j) - (ph_lev%d(i,j,k+1) - ph_lev%d(i,j,k)) * qm%d(i,j,k)
          end do
        end do
      end do
      call fill_halo(mgs)
      end associate
    end do

  end subroutine latlon_bkg_calc_mgs

  subroutine latlon_bkg_calc_mg()

    integer iblk

    if (proc%is_root()) call log_notice('Calculate dry-air pressure.')

    do iblk = 1, size(blocks)
      call calc_mg (blocks(iblk), blocks(iblk)%dstate(1))
      call calc_dmg(blocks(iblk), blocks(iblk)%dstate(1))
    end do

  end subroutine latlon_bkg_calc_mg

  subroutine latlon_bkg_calc_tv()

    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Calculate virtual temperature.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh  , &
                 t    => blocks(iblk)%aux%t , & ! in
                 q    => tracers(iblk)%q    , & ! in
                 qm   => tracers(iblk)%qm   , & ! in
                 tv   => blocks(iblk)%aux%tv)   ! out
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            if (idx_qv == 0) then
              tv%d(i,j,k) = t%d(i,j,k)
            else
              tv%d(i,j,k) = virtual_temperature_from_dry_mixing_ratio(t%d(i,j,k), q%d(i,j,k,idx_qv), qm%d(i,j,k))
            end if
          end do
        end do
      end do
      end associate
    end do

  end subroutine latlon_bkg_calc_tv

  subroutine latlon_bkg_calc_pt()

    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Calculate modified potential temperature.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh        , &
                 ph   => blocks(iblk)%dstate(1)%ph, & ! in
                 t    => blocks(iblk)%aux%t       , & ! in
                 q    => tracers(iblk)%q          , & ! in
                 pt   => blocks(iblk)%dstate(1)%pt)   ! out
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            if (idx_qv == 0) then
              pt%d(i,j,k) = modified_potential_temperature(t%d(i,j,k), ph%d(i,j,k), 0.0_r8)
            else
              pt%d(i,j,k) = modified_potential_temperature(t%d(i,j,k), ph%d(i,j,k), q%d(i,j,k,idx_qv))
            end if
          end do
        end do
      end do
      call fill_halo(pt)
      end associate
    end do

  end subroutine latlon_bkg_calc_pt

  subroutine latlon_bkg_calc_dry_qv()

    integer iblk, i, j, k

    if (idx_qv == 0) return

    if (proc%is_root()) call log_notice('Calculate water vapor dry mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh, &
                 qm   => tracers(iblk)%qm , & ! in
                 q    => tracers(iblk)%q  )   ! out
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              ! NOTE: qm is currently the total wet mixing ratio of water substances.
              q%d(i,j,k,idx_qv) = q%d(i,j,k,idx_qv) / (1 - qm%d(i,j,k))
            end do
          end do
        end do
      end associate
    end do

  end subroutine latlon_bkg_calc_dry_qv

  subroutine latlon_bkg_calc_dry_qc()

    integer iblk, i, j, k

    if (idx_qc == 0) return

    if (proc%is_root()) call log_notice('Calculate cloud liquid dry mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh, &
                 qm   => tracers(iblk)%qm , & ! in
                 q    => tracers(iblk)%q  )   ! out
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              ! NOTE: qm is currently the total wet mixing ratio of water substances.
              q%d(i,j,k,idx_qc) = q%d(i,j,k,idx_qc) / (1 - qm%d(i,j,k))
            end do
          end do
        end do
      end associate
    end do

  end subroutine latlon_bkg_calc_dry_qc

  subroutine latlon_bkg_calc_dry_qi()

    integer iblk, i, j, k

    if (idx_qi == 0) return

    if (proc%is_root()) call log_notice('Calculate cloud ice dry mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh, &
                 qm   => tracers(iblk)%qm , & ! in
                 q    => tracers(iblk)%q  )   ! out
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              ! NOTE: qm is currently the total wet mixing ratio of water substances.
              q%d(i,j,k,idx_qi) = q%d(i,j,k,idx_qi) / (1 - qm%d(i,j,k))
            end do
          end do
        end do
      end associate
    end do

  end subroutine latlon_bkg_calc_dry_qi

  subroutine latlon_bkg_calc_dry_nc()

    integer iblk, i, j, k

    if (idx_nc == 0) return

    if (proc%is_root()) call log_notice('Calculate cloud liquid number dry mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh, &
                 qm   => tracers(iblk)%qm , & ! in
                 q    => tracers(iblk)%q  )   ! out
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              ! NOTE: qm is currently the total wet mixing ratio of water substances.
              q%d(i,j,k,idx_nc) = q%d(i,j,k,idx_nc) / (1 - qm%d(i,j,k))
            end do
          end do
        end do
      end associate
    end do

  end subroutine latlon_bkg_calc_dry_nc

  subroutine latlon_bkg_calc_dry_ni()

    integer iblk, i, j, k

    if (idx_ni == 0) return

    if (proc%is_root()) call log_notice('Calculate cloud ice number dry mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh, &
                 qm   => tracers(iblk)%qm , & ! in
                 q    => tracers(iblk)%q  )   ! out
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              ! NOTE: qm is currently the total wet mixing ratio of water substances.
              q%d(i,j,k,idx_ni) = q%d(i,j,k,idx_ni) / (1 - qm%d(i,j,k))
            end do
          end do
        end do
      end associate
    end do

  end subroutine latlon_bkg_calc_dry_ni

  subroutine latlon_bkg_calc_dry_qr()

    integer iblk, i, j, k

    if (idx_qr == 0) return

    if (proc%is_root()) call log_notice('Calculate rain water dry mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh, &
                 qm   => tracers(iblk)%qm , & ! in
                 q    => tracers(iblk)%q  )   ! out
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              ! NOTE: qm is currently the total wet mixing ratio of water substances.
              q%d(i,j,k,idx_qr) = q%d(i,j,k,idx_qr) / (1 - qm%d(i,j,k))
            end do
          end do
        end do
      end associate
    end do

  end subroutine latlon_bkg_calc_dry_qr

  subroutine latlon_bkg_calc_dry_qs()

    integer iblk, i, j, k

    if (idx_qs == 0) return

    if (proc%is_root()) call log_notice('Calculate snow water dry mixing ratio.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh, &
                 qm   => tracers(iblk)%qm , & ! in
                 q    => tracers(iblk)%q  )   ! out
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              ! NOTE: qm is currently the total wet mixing ratio of water substances.
              q%d(i,j,k,idx_qs) = q%d(i,j,k,idx_qs) / (1 - qm%d(i,j,k))
            end do
          end do
        end do
      end associate
    end do

  end subroutine latlon_bkg_calc_dry_qs

end module latlon_bkg_mod
