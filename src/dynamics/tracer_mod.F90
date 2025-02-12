! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module tracer_mod

  use flogger
  use string
  use const_mod, only: r8
  use namelist_mod
  use block_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use tracer_types_mod
  use interp_mod

  implicit none

  private

  public tracer_init_stage1
  public tracer_init_stage2
  public tracer_final
  public tracer_add
  public tracer_get_idx
  public tracer_calc_qm
  public ntracers
  public ntracers_water
  public nbatches
  public idx_qv
  public idx_qc, idx_nc
  public idx_qi, idx_ni
  public idx_qr, idx_nr
  public idx_qs, idx_ns
  public idx_qg
  public idx_qh
  public idx_qo3
  public idx_qso2
  public batch_names
  public batch_dts
  public tracer_batches
  public tracer_names
  public tracer_long_names
  public tracer_units
  public tracer_types
  public tracers_type
  public is_water_tracer
  public tracers

contains

  subroutine tracer_init_stage1()

    integer iblk

    call tracer_final()

    allocate(batch_names      (10 )); batch_names       = 'N/A'
    allocate(batch_dts        (10 )); batch_dts         = 0
    allocate(tracer_batches   (100)); tracer_batches    = 'N/A'
    allocate(tracer_names     (100)); tracer_names      = 'N/A'
    allocate(tracer_long_names(100)); tracer_long_names = 'N/A'
    allocate(tracer_units     (100)); tracer_units      = 'kg kg-1'
    allocate(tracer_types     (100)); tracer_types      = 0
    allocate(is_water_tracer  (100)); is_water_tracer   = .false.

    allocate(tracers(size(blocks)))
    do iblk = 1, size(blocks)
      call tracers(iblk)%init_stage1(blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, blocks(iblk)%mesh, blocks(iblk)%halo)
    end do

  end subroutine tracer_init_stage1

  subroutine tracer_init_stage2()

    integer iblk, i

    if (ntracers == 0) return

    ! Allocate tracer arrays for each block.
    do iblk = 1, size(blocks)
      call tracers(iblk)%init_stage2(blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, blocks(iblk)%mesh, blocks(iblk)%halo)
    end do

    do iblk = 1, size(blocks)
      ! Allocate physics tendency arrays for dynamics.
      call blocks(iblk)%aux%init_phys(blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, blocks(iblk)%mesh, blocks(iblk)%halo)
    end do

  end subroutine tracer_init_stage2

  subroutine tracer_final()

    nbatches       = 0
    ntracers       = 0
    ntracers_water = 0
    idx_qv         = 0
    idx_qc         = 0
    idx_nc         = 0
    idx_qi         = 0
    idx_ni         = 0
    idx_qr         = 0
    idx_nr         = 0
    idx_qs         = 0
    idx_ns         = 0
    idx_qg         = 0
    idx_qh         = 0
    idx_qo3        = 0

    if (allocated(batch_names      )) deallocate(batch_names      )
    if (allocated(batch_dts        )) deallocate(batch_dts        )
    if (allocated(tracer_batches   )) deallocate(tracer_batches   )
    if (allocated(tracer_names     )) deallocate(tracer_names     )
    if (allocated(tracer_long_names)) deallocate(tracer_long_names)
    if (allocated(tracer_units     )) deallocate(tracer_units     )
    if (allocated(tracer_types     )) deallocate(tracer_types     )
    if (allocated(is_water_tracer  )) deallocate(is_water_tracer  )
    if (allocated(tracers          )) deallocate(tracers          )

  end subroutine tracer_final

  subroutine tracer_add(batch_name, dt, name, long_name, units, type)

    character(*), intent(in) :: batch_name
    real(r8), intent(in) :: dt
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in), optional :: units
    integer, intent(in), optional :: type

    integer i
    logical found

    found = .false.
    do i = 1, nbatches
      if (batch_name == batch_names(i)) then
        found = .true.
        exit
      end if
    end do
    if (.not. found) then
      nbatches = nbatches + 1
      batch_names(nbatches) = batch_name
      batch_dts(nbatches) = dt
    end if

    ntracers = ntracers + 1
    tracer_batches(ntracers) = batch_name
    tracer_names(ntracers) = name
    tracer_long_names(ntracers) = long_name
    if (present(units)) tracer_units(ntracers) = units
    if (present(type)) tracer_types(ntracers) = type

    ! Set tracer indices.
    select case (name)
    case ('qv', 'Q')
      idx_qv    = ntracers; is_water_tracer(ntracers) = .true.; ntracers_water = ntracers_water + 1
    case ('qc', 'CLDLIQ')
      idx_qc    = ntracers; is_water_tracer(ntracers) = .true.; ntracers_water = ntracers_water + 1
    case ('nc', 'NUMLIQ')
      idx_nc    = ntracers
    case ('qi', 'CLDICE')
      idx_qi    = ntracers; is_water_tracer(ntracers) = .true.; ntracers_water = ntracers_water + 1
    case ('ni', 'NUMICE')
      idx_ni    = ntracers
    case ('qr', 'RAINQM')
      idx_qr    = ntracers; is_water_tracer(ntracers) = .true.; ntracers_water = ntracers_water + 1
    case ('qs', 'SNOWQM')
      idx_qs    = ntracers; is_water_tracer(ntracers) = .true.; ntracers_water = ntracers_water + 1
    case ('qg')
      idx_qg    = ntracers; ntracers_water = ntracers_water + 1
    case ('qh')
      idx_qh    = ntracers; ntracers_water = ntracers_water + 1
    case ('qo3')
      idx_qo3   = ntracers
    case ('qso2', 'SO2')
      idx_qso2  = ntracers
    end select

  end subroutine tracer_add

  pure integer function tracer_get_idx(name) result(res)

    character(*), intent(in) :: name

    integer i

    do i = 1, ntracers
      if (name == tracer_names(i)) then
        res = i
        return
      end if
    end do
    res = 0

  end function tracer_get_idx

  subroutine tracer_calc_qm(block)

    type(block_type), intent(in) :: block

    integer i, j, k, m, is, ie, js, je

    if (.not. allocated(tracers)) return
    if (.not. tracers(block%id)%qm%initialized) return

    associate (mesh   => block%mesh              , &
               q      => tracers(block%id)%q     , & ! in
               qm     => tracers(block%id)%qm    , & ! out
               qm_lev => tracers(block%id)%qm_lev)   ! out
    is = mesh%full_ids - 1
    ie = mesh%full_ide + 1
    js = mesh%full_jds - merge(0, 1, mesh%has_south_pole())
    je = mesh%full_jde + merge(0, 1, mesh%has_north_pole())
    qm%d = 0
    do m = 1, ntracers
      if (is_water_tracer(m)) then
        call wait_halo(q, m)
        do k = mesh%full_kds, mesh%full_kde
          do j = js, je
            do i = is, ie
              qm%d(i,j,k) = qm%d(i,j,k) + q%d(i,j,k,m)
            end do
          end do
        end do
      end if
    end do
    if (nonhydrostatic) call interp_run(qm, qm_lev, extrap=.false.)
    end associate

  end subroutine tracer_calc_qm

end module tracer_mod
