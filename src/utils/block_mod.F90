module block_mod

  use mpi
  use container
  use flogger
  use namelist_mod
  use latlon_mesh_mod, mesh_type => latlon_mesh_type
  use latlon_halo_mod, halo_type => latlon_halo_type
  use dynamics_types_mod
  use physics_types_mod
  use adv_batch_mod
  use filter_types_mod
  use accum_mod

  implicit none

  private

  public block_type
  public blocks
  public global_mesh
  public mesh_type
  public halo_type
  public static_type
  public dstate_type
  public dtend_type
  public accum_type

  type block_type
    integer id
    type(mesh_type) filter_mesh
    type(mesh_type) mesh
    type(static_type) static
    type(dstate_type), allocatable :: dstate(:)
    type(dtend_type) dtend
    type(aux_array_type) aux
    type(adv_batch_type) adv_batch_bg
    type(adv_batch_type) adv_batch_pt
    type(adv_batch_type) adv_batch_nh
    type(adv_batch_type), allocatable :: adv_batches(:)
    type(filter_type) big_filter
    type(filter_type) small_filter
    type(halo_type), allocatable :: filter_halo(:)
    type(halo_type), allocatable :: halo(:)
    type(array_type) accum_list
  contains
    procedure :: init_stage1 => block_init_stage1
    procedure :: init_stage2 => block_init_stage2
    procedure :: clear => block_clear
    procedure :: accum => block_accum
    final :: block_final
  end type block_type

  type(block_type), allocatable, target :: blocks(:)

contains

  subroutine block_init_stage1(this, id, ids, ide, jds, jde)

    class(block_type), intent(inout) :: this
    integer, intent(in) :: id
    integer, intent(in) :: ids
    integer, intent(in) :: ide
    integer, intent(in) :: jds
    integer, intent(in) :: jde

    type(accum_type) accum
    integer cell_dims(3)

    this%id = id

    call this%filter_mesh%init_from_parent(global_mesh, this%id, ids, ide, jds, jde)
    call this%mesh%init_from_parent(global_mesh, this%id, ids, ide, jds, jde)
    call this%big_filter%init(this%filter_mesh, 'big_filter')
    call this%small_filter%init(this%filter_mesh, 'small_filter')

    cell_dims = [this%mesh%full_nlon,this%mesh%full_nlat,this%mesh%full_nlev]
    call accum%init(                                &
      name='daily_avg_t'                          , &
      units='K'                                   , &
      long_name='Daily averaged temperature'      , &
      from='dstate'                               , &
      var_name='t'                                , &
      freq=freq_daily                             , &
      stat=stat_avg                               , &
      array_shape=cell_dims                       , &
      active=.false.)
    call this%accum_list%append(accum)
    call accum%init(                                &
      name='daily_max_t'                          , &
      units='K'                                   , &
      long_name='Daily maximum temperature'       , &
      from='dstate'                               , &
      var_name='t'                                , &
      freq=freq_daily                             , &
      stat=stat_max                               , &
      array_shape=cell_dims                       , &
      active=.false.)
    call this%accum_list%append(accum)
    call accum%init(                                &
      name='daily_min_t'                          , &
      units='K'                                   , &
      long_name='Daily minimum temperature'       , &
      from='dstate'                               , &
      var_name='t'                                , &
      freq=freq_daily                             , &
      stat=stat_min                               , &
      array_shape=cell_dims                       , &
      active=.false.)
    call this%accum_list%append(accum)
    call accum%init(                                &
      name='daily_avg_u'                          , &
      units='m s-1'                               , &
      long_name='Daily averaged u-component speed', &
      from='dstate'                               , &
      var_name='u'                                , &
      freq=freq_daily                             , &
      stat=stat_avg                               , &
      array_shape=cell_dims                       , &
      active=.false.)
    call this%accum_list%append(accum)
    call accum%init(                                &
      name='daily_avg_v'                          , &
      units='m s-1'                               , &
      long_name='Daily averaged v-component speed', &
      from='dstate'                               , &
      var_name='v'                                , &
      freq=freq_daily                             , &
      stat=stat_avg                               , &
      array_shape=cell_dims                       , &
      active=.false.)
    call this%accum_list%append(accum)

  end subroutine block_init_stage1

  subroutine block_init_stage2(this)

    class(block_type), intent(inout) :: this

    integer i

    call this%filter_mesh%reinit()

    if (.not. allocated(this%dstate)) then
      select case (trim(time_scheme))
      case ('euler', 'rk2')
        allocate(this%dstate(2))
      case ('pc2', 'wrfrk3')
        allocate(this%dstate(3))
      case ('N/A')
        allocate(this%dstate(1))
      case default
        if (this%id == 0) call log_error('Unknown time scheme ' // trim(time_scheme))
      end select
      if (allocated(this%dstate)) then
        do i = 1, size(this%dstate)
          call this%dstate(i)%init(this%filter_mesh, this%filter_halo, this%mesh, this%halo)
        end do
      end if
      call this%dtend%init(this%filter_mesh, this%filter_halo, this%mesh, this%halo)
      call this%aux%init(this%filter_mesh, this%filter_halo, this%mesh, this%halo)
      call this%static%init_stage1(this%filter_mesh, this%filter_halo, this%mesh, this%halo)
    end if

  end subroutine block_init_stage2

  subroutine block_clear(this)

    class(block_type), intent(inout) :: this

    integer i

    call this%filter_mesh %clear()
    call this%mesh        %clear()
    call this%big_filter  %clear()
    call this%small_filter%clear()
    call this%dtend       %clear()
    call this%static      %clear()
    call this%aux         %clear()
    call this%adv_batch_bg%clear()
    call this%adv_batch_pt%clear()
    call this%adv_batch_nh%clear()
    call this%accum_list  %clear()
    if (allocated(this%dstate     )) then
      do i = 1, size(this%dstate)
        call this%dstate(i)     %clear()
      end do
      deallocate(this%dstate)
    end if
    if (allocated(this%adv_batches)) then
      do i = 1, size(this%adv_batches)
        call this%adv_batches(i)%clear()
      end do
      deallocate(this%adv_batches)
    end if
    if (allocated(this%filter_halo)) then
      do i = 1, size(this%filter_halo)
        call this%filter_halo(i)%clear()
      end do
      deallocate(this%filter_halo)
    end if
    if (allocated(this%halo       )) then
      do i = 1, size(this%halo)
        call this%halo(i)       %clear()
      end do
      deallocate(this%halo)
    end if

  end subroutine block_clear

  subroutine block_accum(this, itime)

    class(block_type), intent(inout) :: this
    integer, intent(in) :: itime

    class(*), pointer :: accum
    integer i, is, ie, js, je, ks, ke

    is = this%mesh%full_ids; ie = this%mesh%full_ide
    js = this%mesh%full_jds; je = this%mesh%full_jde
    ks = this%mesh%full_kds; ke = this%mesh%full_kde

    do i = 1, this%accum_list%size
      accum => this%accum_list%value_at(i)
      select type (accum)
      type is (accum_type)
        if (.not. accum%active) cycle
        select case (accum%from)
        case ('dstate')
          select case (accum%var_name)
          case ('t')
            call accum%accum_run_3d(this%aux%t%d(is:ie,js:je,ks:ke))
          case ('u')
            call accum%accum_run_3d(this%aux%u%d(is:ie,js:je,ks:ke))
          case ('v')
            call accum%accum_run_3d(this%aux%v%d(is:ie,js:je,ks:ke))
          end select
        end select
      end select
    end do

  end subroutine block_accum

  subroutine block_final(this)

    type(block_type), intent(inout) :: this

    call this%clear()

  end subroutine block_final

end module block_mod
