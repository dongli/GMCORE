! This vertical coordinate calculation is copied from WRF v4.4.1.
!
! History:
! 
!   - 2021-04-02: First version.

module hybrid_coord_wrf_mod

  use flogger
  use const_mod
  use namelist_mod, only: nlev, tiso, dzmax, dzbot, dzstretch_s, dzstretch_u, eta_b
  use latlon_mesh_mod, only: global_mesh
  use latlon_parallel_types_mod, only: proc

  implicit none

  private

  public hybrid_coord_wrf
  public wrf_compute_eta

contains

  subroutine hybrid_coord_wrf(p0, ptop, hyai, hybi)

    real(r8), intent(in) :: p0
    real(r8), intent(in) :: ptop
    real(r8), intent(out) :: hyai(:)
    real(r8), intent(out) :: hybi(:)

    real(r8) b0, b1, b2, b3, eta
    integer k

    call wrf_compute_eta(nlev, ptop, global_mesh%half_lev(1:nlev+1))

    b0 = 2 * eta_b**2 / (1 - eta_b)**3
    b1 = -eta_b * (4 + eta_b + eta_b**2) / (1 - eta_b)**3
    b2 = 2 * (1 + eta_b + eta_b**2) / (1 - eta_b)**3
    b3 = -(1 + eta_b) / (1 - eta_b)**3

    do k = 1, nlev + 1
      eta = global_mesh%half_lev(k)
      if (eta >= eta_b) hybi(k) = b0 + b1 * eta + b2 * eta**2 + b3 * eta**3
      hyai(k) = eta - hybi(k)
    end do

  end subroutine hybrid_coord_wrf

  subroutine wrf_compute_eta(nlev, ptop, eta_lev)

    integer, intent(in) :: nlev
    real(r8), intent(in) :: ptop
    real(r8), intent(out) :: eta_lev(nlev+1)

    real(r8) ztop, pbot, dz
    real(r8) z_lev(nlev+1), p_lev(nlev+1)
    integer k, k0

    ztop = rd * tiso / g * log(p0 / ptop)
    pbot = p0 * exp(-g * dzbot / (rd * tiso))

    eta_lev(1) = 1.0_r8
    eta_lev(2) = (pbot - ptop) / (p0 - ptop)
    z_lev  (1) = 0.0_r8
    z_lev  (2) = dzbot
    p_lev  (1) = p0
    p_lev  (2) = pbot

    ! The stretching transitions from dzstretch_s to dzstretch_u by the time the thickness reaches dzmax/2.
    dz = dzbot
    k0 = 0
    do k = 3, nlev + 1
      dz = dz * (dzstretch_u + (dzstretch_s - dzstretch_u) * max((dzmax / 2 - dz) / (dzmax / 2), 0.0_r8))
      if ((ztop - z_lev(k-1)) / (nlev - k + 1) < dz) then
        k0 = k - 1
        exit
      end if
      z_lev(k) = z_lev(k-1) + dz
      p_lev(k) = p0 * exp(-g * z_lev(k) / (rd * tiso))
      eta_lev(k) = (p_lev(k) - ptop) / (p0 - ptop)
      if (k == nlev + 1) then
        call log_error('Too few vertical levels for given parameters. Increase nlev!', __FILE__, __LINE__, pid=proc%id_model)
      end if
    end do

    dz = (ztop - z_lev(k0)) / (nlev - k0 + 1)
    if (dz > 1.5 * dzmax) then
      call log_error('Upper levels may be too coarse!', __FILE__, __LINE__, pid=proc%id_model)
    end if

    do k = k0 + 1, nlev + 1
      z_lev(k) = z_lev(k-1) + dz
      p_lev(k) = p0 * exp(-g * z_lev(k) / (rd * tiso))
      eta_lev(k) = (p_lev(k) - ptop) / (p0 - ptop)
    end do
    eta_lev(nlev+1) = 0.0_r8

    ! Reverse order to up-bottom.
    eta_lev = eta_lev(nlev+1:1:-1)

  end subroutine wrf_compute_eta

end module hybrid_coord_wrf_mod
