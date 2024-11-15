subroutine opt_dst(q, pl, qxv, qsv, gv, qxi, qsi, gi, qextref, tauref, taudst)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_tracers_mod

  implicit none

  real(r8), intent(in ) :: q      (nlev,ntracers)
  real(r8), intent(in ) :: pl     (2*nlev+3)
  real(r8), intent(out) :: qxv    (2*nlev+4,nspectv)
  real(r8), intent(out) :: qsv    (2*nlev+4,nspectv)
  real(r8), intent(out) :: gv     (2*nlev+4,nspectv)
  real(r8), intent(out) :: qxi    (2*nlev+4,nspecti)
  real(r8), intent(out) :: qsi    (2*nlev+4,nspecti)
  real(r8), intent(out) :: gi     (2*nlev+4,nspecti)
  real(r8), intent(out) :: qextref(2*nlev+4)
  real(r8), intent(out) :: tauref (2*nlev+4)
  real(r8), intent(out) :: taudst(2)

  integer n, i, k, l
  real(r8) dev2, cst
  real(r8) Mo, No, Rs, Ao
  real(r8) surf(nbin_rt)

  n = 2 * nlev + 3

  do k = 1, n + 1
    qextref(k) = 1
    tauref (k) = 0
  end do

  do i = 1, nspectv
    do k = 1, n + 1
      qxv(k,i) = qextref(k)
      qsv(k,i) = qextref(k) * 0.99_r8
      gv (k,i) = 0
    end do
  end do
  do i = 1, nspecti
    do k = 1, n + 1
      qxi(k,i) = qextref(k)
      qsi(k,i) = qextref(k) * 0.99_r8
      gi (k,i) = 0
    end do
  end do

  dev2 = 1.0_r8 / (sqrt(2.0_r8) * dev_dt)
  cst = 0.75_r8 / (pi * dpden_dt)

  taudst = 0

  do l = 1, nlev
    if (q(l,iMa_dt) > 1.0e-8_r8 .and. q(l,iNb_dt) > 1) then
      ! Calculate the cross-section mean radius (Rs) of the log-normal distribution.
      Mo = q(l,iMa_dt)
      No = q(l,iNb_dt) + 1
      Rs = (Mo / No * cst)**athird * exp(-0.5_r8 * dev_dt**2)
      ! Calculate the total cross sectional area (Ao) of water ice particles.
      Ao = No * pi * Rs**2
      ! Define the cross-section weighted distribution, i.e. surface/size bin. Change Rs to Reff.
      Rs = Rs * exp(1.5_r8 * dev_dt**2)
      Rs = 1.0_r8 / min(max(Rs, 1.0e-7_r8), 50.0e-6_r8)
      do i = 1, nbin_rt
        ! surf(i) = 0.5_r8 * (erf(log()))
      end do
    end if
  end do

end subroutine opt_dst