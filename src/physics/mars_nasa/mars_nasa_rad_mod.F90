! ==============================================================================
! This file is part of GoMars since 2023.
!
! GoMars is a Martian general circulation model developed in Institute of
! Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
!
! GMCORE is a dynamical core for atmospheric model used in GoMars.
!
! GoMars and GMCORE are distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module mars_nasa_rad_mod

  use fiona
  use mars_nasa_const_mod
  use mars_nasa_namelist_mod
  use mars_nasa_physics_types_mod
  use mars_nasa_tracers_mod
  use mars_nasa_spectra_mod
  use mars_nasa_optics_mod
  use mars_nasa_solar_mod
  use mars_nasa_rad_kcoef_mod

  implicit none

  private

  public mars_nasa_rad_init
  public mars_nasa_rad_final
  public mars_nasa_rad_run

  ! Number of radiation levels
  integer nlev_rad
  ! Maximum optical depth
  integer, parameter :: taumax = 35

  ! Reference pressure for Rayleigh scattering (Pa)
  real(r8), parameter :: ray_p0 = 9.423e6_r8
  ! Rayleigh scattering optical depth at each vsible spectral interval
  real(r8), allocatable, dimension(:  ) :: tauray_vs
  ! Integrated Planck function lookup table at each IR spectral interval (W m-2 cm-1)
  real(r8), allocatable, dimension(:,:) :: planck_ir

contains

  subroutine mars_nasa_rad_init(nlev)

    ! Number of model full levels
    integer, intent(in) :: nlev

    real(r8) a, b, t
    integer i, j

    call mars_nasa_rad_final()

    nlev_rad = 2 * nlev + 3

    call mars_nasa_spectra_init()
    call mars_nasa_optics_init(nlev_rad)
    call mars_nasa_solar_init()
    call mars_nasa_rad_kcoef_init()

    allocate(tauray_vs(spec_vs%n     ))
    allocate(planck_ir(spec_ir%n,8501))

    do i = 1, spec_vs%n
      tauray_vs(i) = 8.7_r8 / g * (1.527_r8 * (1 + 0.013_r8 / spec_vs%wl(i)**2) / spec_vs%wl(i)**4) / (ray_p0 / 100)
    end do

    ! For each IR wavelength interval, compute the integral of B(T), the
    ! Planck function, divided by the wavelength interval, in cm-1.  The
    ! integration is in MKS units (W m^-2 cm^-1).
    do i = 1, spec_ir%n
      a = 0.5d-2 * (spec_ir%bwn(i+1) - spec_ir%bwn(i)) / (spec_ir%bwn(i) * spec_ir%bwn(i+1))
      b = 0.5d-2 * (spec_ir%bwn(i+1) + spec_ir%bwn(i)) / (spec_ir%bwn(i) * spec_ir%bwn(i+1))
      do j = 500, 9000
        t = j * 0.1_r8
        planck_ir(i,j-499) = integrate_planck_function(a, b, t) * a / (pi * spec_ir%dwn(i))
      end do
    end do

  end subroutine mars_nasa_rad_init

  subroutine mars_nasa_rad_final()

    call mars_nasa_spectra_final()
    call mars_nasa_optics_final()
    call mars_nasa_solar_final()
    call mars_nasa_rad_kcoef_final()

    if (allocated(planck_ir)) deallocate(planck_ir)
    if (allocated(tauray_vs)) deallocate(tauray_vs)

  end subroutine mars_nasa_rad_final

  subroutine mars_nasa_rad_run(state, tend)

    type(mars_nasa_state_type), intent(inout) :: state
    type(mars_nasa_tend_type), intent(inout) :: tend

    real(r8) t_rad(nlev_rad)
    real(r8) p_rad(nlev_rad)
    real(r8) qh2o_rad(nlev_rad)
    real(r8) tau_vs(spec_vs%n,ngauss)
    real(r8) tau_vs_sfc(spec_vs%n,ngauss)
    integer icol, k

    associate (mesh  => state%mesh , &
               t_top => state%t_top, &
               tg    => state%tg   , &
               t     => state%t    , &
               p     => state%p    , &
               p_lev => state%p_lev, &
               q     => state%q    )
    call calc_fdns_dir(state)
    state%fdns = state%fdns_dir + state%fdns_dif
    do icol = 1, mesh%ncol
      ! Set radiation levels (t, p and qh2o).
      p_rad(1) = mesh%ptop / 2.0_r8
      p_rad(2) = mesh%ptop / 2.0_r8
      do k = 1, mesh%nlev
        p_rad(2*(k-1)+4) = p(icol,k)
      end do
      do k = 1, mesh%nlev + 1
        p_rad(2*(k-1)+3) = p_lev(icol,k)
      end do
      t_rad(1) = t_top(icol)
      t_rad(2) = t_top(icol)
      t_rad(3) = t_top(icol)
      do k = 1, mesh%nlev
        t_rad(2*(k-1)+4) = t(icol,k)
      end do
      do k = 5, nlev_rad - 2, 2
        t_rad(k) = t_rad(k+1) + (t_rad(k-1) - t_rad(k+1)) * &
                   log(p_rad(k  ) / p_rad(k+1)) / &
                   log(p_rad(k-1) / p_rad(k+1))
      end do
      t_rad(nlev_rad) = tg(icol)
      if (active_water) then
        do k = 1, mesh%nlev
          qh2o_rad(2*(k-1)+4) = m_co2 / m_h2o * q(icol,k,idx_m_vap)
          qh2o_rad(2*(k-1)+5) = qh2o_rad(2*(k-1)+4)
        end do
      else
        do k = 1, nlev_rad
          qh2o_rad(k) = 1.0e-7_r8
        end do
      end if
      ! Interpolate K coefficient and calculate gas optical depth.
      tau_vs_sfc = 0
      do k = 2, nlev_rad
        call get_kcoef_vs(t_rad(k), p_rad(k), qh2o_rad(k), tau_vs(:,:))
        tau_vs = cmk * (p_rad(k) - p_rad(k-1)) * tau_vs
        tau_vs_sfc = tau_vs_sfc + tau_vs
        ! print *, k, t_rad(k), p_rad(k), qh2o_rad(k), tau_vs(1,1)
      end do
      ! stop 999
    end do
    end associate

  end subroutine mars_nasa_rad_run

  pure real(r8) function integrate_planck_function(s1, s2, t) result(res)

    real(r8), intent(in) :: s1
    real(r8), intent(in) :: s2
    real(r8), intent(in) :: t   ! Temperature (K)

    ! C1 and C2 values from Goody and Yung (2nd edition) MKS units
    ! These values lead to a "sigma" (sigma*T^4) of 5.67032E-8 W m^-2 K^-4
    real(8), parameter :: c1 = 3.741832d-16      ! W m-2
    real(8), parameter :: c2 = 1.438786d-2       ! m K
    real(8), parameter :: x(12) = [ & ! Quadrature points
      -0.981560634246719D0,  -0.904117256370475D0, &
      -0.769902674194305D0,  -0.587317954286617D0, &
      -0.367831498998180D0,  -0.125233408511469D0, &
       0.125233408511469D0,   0.367831498998180D0, &
       0.587317954286617D0,   0.769902674194305D0, &
       0.904117256370475D0,   0.981560634246719D0  &
    ]
    real(8), parameter :: w(12) = [ & ! Quadrature weights
       0.047175336386512D0,   0.106939325995318D0, &
       0.160078328543346D0,   0.203167426723066D0, &
       0.233492536538355D0,   0.249147045813403D0, &
       0.249147045813403D0,   0.233492536538355D0, &
       0.203167426723066D0,   0.160078328543346D0, &
       0.106939325995318D0,   0.047175336386512D0  &
    ]
    real(8) wl
    integer i

    res = 0
    do i = 1, 12
      wl = s1 * x(i) + s2
      res = res + w(i) * c1 / (wl**5 * (exp(c2 / wl / t) - 1)) ! Planck function
    end do

  end function integrate_planck_function

  subroutine calc_fdns_dir(state)

    ! Calculate direct solar flux at surface.

    type(mars_nasa_state_type), intent(inout) :: state

    integer icol, is, ig
    real(r8) c

    associate (mesh     => state%mesh    , &
               cosz     => state%cosz    , & ! in
               detau    => state%detau   , & ! in
               fdns_dir => state%fdns_dir)
    do icol = 1, mesh%ncol
      fdns_dir(icol) = 0
      if (cosz(icol) > 0) then
        do is = 1, spec_vs%n
          c = cosz(icol) * fsol_spec_mars(is)
          do ig = 1, ngauss - 1
            if (detau(is,ig,icol) <= 5) then
              fdns_dir(icol) = fdns_dir(icol) + c * exp(-detau(is,ig,icol) / cosz(icol)) * gwgt(ig) * (1 - f0_vs(is))
            end if
          end do
          ig = ngauss
          if (detau(is,ig,icol) <= 5) then
            fdns_dir(icol) = fdns_dir(icol) + c * exp(-detau(is,ig,icol) / cosz(icol)) * f0_vs(is)
          end if
        end do
      end if
    end do
    end associate

  end subroutine calc_fdns_dir

end module mars_nasa_rad_mod
