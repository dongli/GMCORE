module weno_mod

  use const_mod
  use namelist_mod
  use adv_batch_mod
  use latlon_field_types_mod
  use latlon_parallel_mod
  use perf_mod

  implicit none

  private

  public weno_calc_tracer_hflx
  public weno_calc_tracer_vflx
  public weno3
  public weno5

contains

  subroutine weno_calc_tracer_hflx(batch, q, qmfx, qmfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfx
    type(latlon_field3d_type), intent(inout) :: qmfy
    real(r8), intent(in), optional :: dt

    integer ks, ke, i, j, k

    call perf_start('weno_calc_tracer_hflx')

    associate (mesh => q%mesh   , &
               mfx  => batch%mfx, & ! in
               mfy  => batch%mfy)   ! in
    select case (batch%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, batch%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, batch%loc == 'cell')
      select case (weno_order_h)
      case (3)
        do k = ks, ke
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              qmfx%d(i,j,k) = weno3(sign(1.0_r8, mfx%d(i,j,k)), q%d(i-1:i+2,j,k)) * mfx%d(i,j,k)
            end do
          end do
          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfy%d(i,j,k) = weno3(sign(1.0_r8, mfy%d(i,j,k)), q%d(i,j-1:j+2,k)) * mfy%d(i,j,k)
            end do
          end do
        end do
      case (5)
        do k = ks, ke
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              qmfx%d(i,j,k) = weno5(sign(1.0_r8, mfx%d(i,j,k)), q%d(i-2:i+3,j,k)) * mfx%d(i,j,k)
            end do
          end do
          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfy%d(i,j,k) = weno5(sign(1.0_r8, mfy%d(i,j,k)), q%d(i,j-2:j+3,k)) * mfy%d(i,j,k)
            end do
          end do
        end do
      end select
    end select
    call fill_halo(qmfx, south_halo=.false., north_halo=.false., east_halo =.false.)
    call fill_halo(qmfy, west_halo =.false., east_halo =.false., north_halo=.false.)
    end associate

    call perf_stop('weno_calc_tracer_hflx')

  end subroutine weno_calc_tracer_hflx

  subroutine weno_calc_tracer_vflx(batch, q, qmfz, dt)

    type(adv_batch_type), intent(inout) :: batch
    type(latlon_field3d_type), intent(in) :: q
    type(latlon_field3d_type), intent(inout) :: qmfz
    real(r8), intent(in), optional :: dt

    integer i, j, k

    call perf_start('weno_calc_tracer_vflx')

    associate (mesh => q%mesh   , &
               mfz  => batch%mfz)   ! in
    select case (batch%loc)
    case ('cell')
      select case (weno_order_v)
      case (3)
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfz%d(i,j,k) = weno3(sign(1.0_r8, mfz%d(i,j,k)), q%d(i,j,k-2:k+1)) * mfz%d(i,j,k)
            end do
          end do
        end do
      case (5)
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfz%d(i,j,k) = weno5(sign(1.0_r8, mfz%d(i,j,k)), q%d(i,j,k-3:k+2)) * mfz%d(i,j,k)
            end do
          end do
        end do
      end select
    case ('lev')
      select case (weno_order_v)
      case (3)
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfz%d(i,j,k) = weno3(sign(1.0_r8, mfz%d(i,j,k)), q%d(i,j,k-1:k+2)) * mfz%d(i,j,k)
            end do
          end do
        end do
      case (5)
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qmfz%d(i,j,k) = weno5(sign(1.0_r8, mfz%d(i,j,k)), q%d(i,j,k-2:k+3)) * mfz%d(i,j,k)
            end do
          end do
        end do
      end select
    end select
    end associate

    call perf_stop('weno_calc_tracer_vflx')

  end subroutine weno_calc_tracer_vflx

  pure real(r8) function weno3(dir, f) result(res)

    real(r8), intent(in) :: dir
    real(r8), intent(in) :: f(4)

    real(r8), parameter :: c11 = -0.5_r8, c12 =  1.5_r8
    real(r8), parameter :: c21 =  0.5_r8, c22 =  0.5_r8
    real(r8), parameter :: weno_coef(2) = [1.0_r8 / 4.0_r8, 3.0_r8 / 4.0_r8]
    real(r8), parameter :: eps = 1.0e-16_r8

    real(r8) fs(2), beta(2), alpha(2)

    ! o-----o--x--o-----o
    ! 1     2     3     4
    select case (int(dir))
    case (1)
      fs(1) = c11 * f(1) + c12 * f(2); beta(1) = (f(2) - f(1))**2
      fs(2) = c21 * f(2) + c22 * f(3); beta(2) = (f(3) - f(2))**2
    case (-1)
      fs(1) = c11 * f(4) + c12 * f(3); beta(1) = (f(3) - f(4))**2
      fs(2) = c21 * f(3) + c22 * f(2); beta(2) = (f(2) - f(3))**2
    end select

    alpha = weno_coef / (eps + beta)**2
    res = dot_product(fs, alpha / sum(alpha))

  end function weno3

  pure real(r8) function weno5(dir, f) result(res)

    real(r8), intent(in) :: dir
    real(r8), intent(in) :: f(6)

    real(r8), parameter :: c11 =  1.0_r8 / 3.0_r8, c12 = -7.0_r8 / 6.0_r8, c13 = 11.0_r8 / 6.0_r8
    real(r8), parameter :: c21 = -1.0_r8 / 6.0_r8, c22 =  5.0_r8 / 6.0_r8, c23 =  1.0_r8 / 3.0_r8
    real(r8), parameter :: c31 =  1.0_r8 / 3.0_r8, c32 =  5.0_r8 / 6.0_r8, c33 = -1.0_r8 / 6.0_r8
    real(r8), parameter :: b1 = 13.0_r8 / 12.0_r8, b2 = 0.25_r8
    real(r8), parameter :: weno_coef(3) = [0.1_r8, 0.6_r8, 0.3_r8]
    real(r8), parameter :: eps = 1.0e-16_r8

    real(r8) fs(3), beta(3), alpha(3)

    ! o-----o-----o--x--o-----o-----o
    ! 1     2     3     4     5     6
    select case (int(dir))
    case (1)
      fs(1) = c11 * f(1) + c12 * f(2) + c13 * f(3)
      fs(2) = c21 * f(2) + c22 * f(3) + c23 * f(4)
      fs(3) = c31 * f(3) + c32 * f(4) + c33 * f(5)
      beta(1) = b1 * (f(1) - 2 * f(2) + f(3))**2 + b2 * (    f(1) - 4 * f(2) + 3 * f(3))**2
      beta(2) = b1 * (f(2) - 2 * f(3) + f(4))**2 + b2 * (    f(2)            -     f(4))**2
      beta(3) = b1 * (f(3) - 2 * f(4) + f(5))**2 + b2 * (3 * f(3) - 4 * f(4) +     f(5))**2
    case (-1)
      fs(1) = c11 * f(6) + c12 * f(5) + c13 * f(4)
      fs(2) = c21 * f(5) + c22 * f(4) + c23 * f(3)
      fs(3) = c31 * f(4) + c32 * f(3) + c33 * f(2)
      beta(1) = b1 * (f(6) - 2 * f(5) + f(4))**2 + b2 * (    f(6) - 4 * f(5) + 3 * f(4))**2
      beta(2) = b1 * (f(5) - 2 * f(4) + f(3))**2 + b2 * (    f(5)            -     f(3))**2
      beta(3) = b1 * (f(4) - 2 * f(3) + f(2))**2 + b2 * (3 * f(4) - 4 * f(3) +     f(2))**2
    end select

    alpha = weno_coef / (eps + beta)**2
    res = dot_product(fs, alpha / sum(alpha))

  end function weno5

end module weno_mod
