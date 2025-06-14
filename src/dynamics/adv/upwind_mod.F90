! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module upwind_mod

  use const_mod
  use namelist_mod
  use adv_batch_mod
  use latlon_field_types_mod
  use latlon_parallel_mod
  use perf_mod

  implicit none

  private

  public upwind_calc_tracer_hflx
  public upwind_calc_tracer_vflx
  public upwind1
  public upwind3
  public upwind5

contains

  subroutine upwind_calc_mass_hflx(batch, m, mfx, mfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: m
    type(latlon_field3d_type), intent(inout) :: mfx
    type(latlon_field3d_type), intent(inout) :: mfy
    real(r8), intent(in), optional :: dt

    integer ks, ke, i, j, k

    call perf_start('upwind_calc_mass_hflx')

    associate (mesh => m%mesh , &
               u    => batch%u, & ! in
               v    => batch%v)   ! in
    select case (batch%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, batch%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, batch%loc == 'cell')
      select case (upwind_order_h)
      case (3)
        do k = ks, ke
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              mfx%d(i,j,k) = upwind3(sign(1.0_r8, u%d(i,j,k)), upwind_wgt, m%d(i-1:i+2,j,k)) * u%d(i,j,k)
            end do
          end do
          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              mfy%d(i,j,k) = upwind3(sign(1.0_r8, v%d(i,j,k)), upwind_wgt, m%d(i,j-1:j+2,k)) * v%d(i,j,k)
            end do
          end do
        end do
      case (5)
        do k = ks, ke
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              mfx%d(i,j,k) = upwind5(sign(1.0_r8, u%d(i,j,k)), upwind_wgt, m%d(i-2:i+3,j,k)) * u%d(i,j,k)
            end do
          end do
          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              mfy%d(i,j,k) = upwind5(sign(1.0_r8, v%d(i,j,k)), upwind_wgt, m%d(i,j-2:j+3,k)) * v%d(i,j,k)
            end do
          end do
        end do
      end select
    end select
    call fill_halo(mfx, south_halo=.false., north_halo=.false., east_halo =.false.)
    call fill_halo(mfy, west_halo =.false., east_halo =.false., north_halo=.false.)
    end associate

    call perf_stop('upwind_calc_mass_hflx')

  end subroutine upwind_calc_mass_hflx

  subroutine upwind_calc_tracer_hflx(batch, q, qmfx, qmfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfx
    type(latlon_field3d_type), intent(inout) :: qmfy
    real(r8), intent(in), optional :: dt

    integer ks, ke, i, j, k

    call perf_start('upwind_calc_tracer_hflx')

    associate (mesh => q%mesh   , &
               mfx  => batch%mfx, & ! in
               mfy  => batch%mfy)   ! in
    select case (batch%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, batch%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, batch%loc == 'cell')
      select case (upwind_order_h)
      case (0)
        do k = ks, ke
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              qmfx%d(i,j,k) = 0.5_r8 * (q%d(i,j,k) + q%d(i+1,j,k)) * mfx%d(i,j,k)
            end do
          end do
          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfy%d(i,j,k) = 0.5_r8 * (q%d(i,j,k) + q%d(i,j+1,k)) * mfy%d(i,j,k)
            end do
          end do
        end do
      case (1)
        do k = ks, ke
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              qmfx%d(i,j,k) = upwind1(sign(1.0_r8, mfx%d(i,j,k)), upwind_wgt, q%d(i:i+1,j,k)) * mfx%d(i,j,k)
            end do
          end do
          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfy%d(i,j,k) = upwind1(sign(1.0_r8, mfy%d(i,j,k)), upwind_wgt, q%d(i,j:j+1,k)) * mfy%d(i,j,k)
            end do
          end do
        end do
      case (3)
        do k = ks, ke
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              qmfx%d(i,j,k) = upwind3(sign(1.0_r8, mfx%d(i,j,k)), upwind_wgt, q%d(i-1:i+2,j,k)) * mfx%d(i,j,k)
            end do
          end do
          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfy%d(i,j,k) = upwind3(sign(1.0_r8, mfy%d(i,j,k)), upwind_wgt, q%d(i,j-1:j+2,k)) * mfy%d(i,j,k)
            end do
          end do
        end do
      case (5)
        do k = ks, ke
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              qmfx%d(i,j,k) = upwind5(sign(1.0_r8, mfx%d(i,j,k)), upwind_wgt, q%d(i-2:i+3,j,k)) * mfx%d(i,j,k)
            end do
          end do
          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfy%d(i,j,k) = upwind5(sign(1.0_r8, mfy%d(i,j,k)), upwind_wgt, q%d(i,j-2:j+3,k)) * mfy%d(i,j,k)
            end do
          end do
        end do
      end select
    end select
    call fill_halo(qmfx, south_halo=.false., north_halo=.false., east_halo =.false.)
    call fill_halo(qmfy, west_halo =.false., east_halo =.false., north_halo=.false.)
    end associate

    call perf_stop('upwind_calc_tracer_hflx')

  end subroutine upwind_calc_tracer_hflx

  subroutine upwind_calc_tracer_vflx(batch, q, qmfz, dt)

    type(adv_batch_type), intent(inout) :: batch
    type(latlon_field3d_type), intent(in) :: q
    type(latlon_field3d_type), intent(inout) :: qmfz
    real(r8), intent(in), optional :: dt

    integer i, j, k

    call perf_start('upwind_calc_tracer_vflx')

    associate (mesh => q%mesh   , &
               mfz  => batch%mfz)   ! in
    select case (batch%loc)
    case ('cell')
      select case (upwind_order_v)
      case (0)
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfz%d(i,j,k) = 0.5_r8 * (q%d(i,j,k-1) + q%d(i,j,k)) * mfz%d(i,j,k)
            end do
          end do
        end do
      case (1)
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfz%d(i,j,k) = upwind1(sign(1.0_r8, mfz%d(i,j,k)), upwind_wgt, q%d(i,j,k-1:k)) * mfz%d(i,j,k)
            end do
          end do
        end do
      case (3)
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfz%d(i,j,k) = upwind3(sign(1.0_r8, mfz%d(i,j,k)), upwind_wgt, q%d(i,j,k-2:k+1)) * mfz%d(i,j,k)
            end do
          end do
        end do
      case (5)
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfz%d(i,j,k) = upwind5(sign(1.0_r8, mfz%d(i,j,k)), upwind_wgt, q%d(i,j,k-3:k+2)) * mfz%d(i,j,k)
            end do
          end do
        end do
      end select
    case ('lev')
      select case (upwind_order_v)
      case (0)
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfz%d(i,j,k) = 0.5_r8 * (q%d(i,j,k) + q%d(i,j,k+1)) * mfz%d(i,j,k)
            end do
          end do
        end do
      case (1)
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfz%d(i,j,k) = upwind1(sign(1.0_r8, mfz%d(i,j,k)), upwind_wgt, q%d(i,j,k:k+1)) * mfz%d(i,j,k)
            end do
          end do
        end do
      case (3)
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfz%d(i,j,k) = upwind3(sign(1.0_r8, mfz%d(i,j,k)), upwind_wgt, q%d(i,j,k-1:k+2)) * mfz%d(i,j,k)
            end do
          end do
        end do
      case (5)
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfz%d(i,j,k) = upwind5(sign(1.0_r8, mfz%d(i,j,k)), upwind_wgt, q%d(i,j,k-2:k+3)) * mfz%d(i,j,k)
            end do
          end do
        end do
      end select
    end select
    end associate

    call perf_stop('upwind_calc_tracer_vflx')

  end subroutine upwind_calc_tracer_vflx

  pure real(r8) function upwind1(dir, wgt, f) result(res)

    real(r8), intent(in) :: dir
    real(r8), intent(in) :: f(0:1)
    real(r8), intent(in) :: wgt

    real(r8), parameter :: c11 =  0.5_r8
    real(r8), parameter :: c12 = -0.5_r8

    res = c11 * (f(1) + f(0)) + c12 * (f(1) - f(0)) * wgt * dir

  end function upwind1

  pure real(r8) function upwind3(dir, wgt, f) result(res)

    real(r8), intent(in) :: dir
    real(r8), intent(in) :: wgt
    real(r8), intent(in) :: f(-1:2)

    real(r8), parameter :: c31 =  7.0_r8 / 12.0_r8
    real(r8), parameter :: c32 = -1.0_r8 / 12.0_r8
    real(r8), parameter :: c33 =  1.0_r8 / 12.0_r8

    res = c31 * (f(1) + f( 0)) &
        + c32 * (f(2) + f(-1)) &
        + c33 * (f(2) - f(-1) - 3 * (f(1) - f(0))) * wgt * dir

  end function upwind3

  pure real(r8) function upwind5(dir, wgt, f) result(res)

    real(r8), intent(in) :: dir
    real(r8), intent(in) :: wgt
    real(r8), intent(in) :: f(-2:3)

    real(r8), parameter :: c41 = 37.0_r8 / 60.0_r8
    real(r8), parameter :: c42 = -2.0_r8 / 15.0_r8
    real(r8), parameter :: c43 =  1.0_r8 / 60.0_r8
    real(r8), parameter :: c44 = -1.0_r8 / 60.0_r8

    res = c41 * (f(1) + f( 0)) &
        + c42 * (f(2) + f(-1)) &
        + c43 * (f(3) + f(-2)) &
        + c44 * (f(3) - f(-2) - 5 * (f(2) - f(-1)) + 10 * (f(1) - f(0))) * wgt * dir

  end function upwind5

end module upwind_mod
