subroutine opt_dst(q, plev, qxv, qsv, gv, qxi, qsi, gi, qextrefdst, taurefdst, taudst)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_tracers_mod
  use gomars_v1_mp_mod

  implicit none

  real(r8), intent(in ) :: q         (nlev,ntracers)
  real(r8), intent(in ) :: plev      (2*nlev+3)
  real(r8), intent(out) :: qxv       (2*nlev+4,nspectv)
  real(r8), intent(out) :: qsv       (2*nlev+4,nspectv)
  real(r8), intent(out) :: gv        (2*nlev+4,nspectv)
  real(r8), intent(out) :: qxi       (2*nlev+4,nspecti)
  real(r8), intent(out) :: qsi       (2*nlev+4,nspecti)
  real(r8), intent(out) :: gi        (2*nlev+4,nspecti)
  real(r8), intent(out) :: qextrefdst(2*nlev+4)
  real(r8), intent(out) :: taurefdst (2*nlev+4)
  real(r8), intent(out) :: taudst    (2)

  integer n, k, l, is, ib
  real(r8) dev2, cst
  real(r8) Mo, No, Rs, Ao
  real(r8) surf(nbin_rt)

  n = 2 * nlev + 3

  do l = 1, n + 1
    qextrefdst(l) = 1
    taurefdst (l) = 0
  end do

  do is = 1, nspectv
    do l = 1, n + 1
      qxv(l,is) = qextrefdst(l)
      qsv(l,is) = qextrefdst(l) * 0.99_r8
      gv (l,is) = 0
    end do
  end do
  do is = 1, nspecti
    do l = 1, n + 1
      qxi(l,is) = qextrefdst(l)
      qsi(l,is) = qextrefdst(l) * 0.99_r8
      gi (l,is) = 0
    end do
  end do

  dev2 = 1.0_r8 / (sqrt2 * dev_dst)
  cst = 0.75_r8 / (pi * rho_dst)

  taudst = 0

  do k = 1, nlev
    if (q(k,iMa_dst) > 1.0e-8_r8 .and. q(k,iNb_dst) > 1) then
      ! Calculate the cross-section mean radius (Rs) of the log-normal distribution.
      Mo = q(k,iMa_dst)
      No = q(k,iNb_dst) + 1
      Rs = (Mo / No * cst)**athird * exp(-0.5_r8 * dev_dst**2)
      ! Calculate the total cross sectional area (Ao) of water ice particles.
      Ao = No * pi * Rs**2
      ! Define the cross-section weighted distribution, i.e. surface/size bin. Change Rs to Reff.
      Rs = 1.0_r8 / min(max(Rs * exp(1.5_r8 * dev_dst**2), 1.0e-7_r8), 50.0e-6_r8)
      do ib = 1, nbin_rt
        surf(ib) = 0.5_r8 * (erf(log(radb_rt(ib+1) * Rs) * dev2) - erf(log(radb_rt(ib) * Rs) * dev2))
      end do
      ! Calculate the averaged values of the optical properties for the whole distribution.
      do l = 2 * k + 2, 2 * k + 3
        do is = 1, nspectv
          qxv(l,is) = 0
          qsv(l,is) = 0
          do ib = 1, nbin_rt
            qxv(l,is) = qxv(l,is) + qextv_dst (ib,is) * surf(ib)
            qsv(l,is) = qsv(l,is) + qscatv_dst(ib,is) * surf(ib)
            gv (l,is) = gv (l,is) + gv_dst    (ib,is) * surf(ib)
          end do
          qsv(l,is) = min(qsv(l,is), 0.99999_r8 * qxv(l,is))
        end do
        do is = 1, nspecti
          qxi(l,is) = 0
          qsi(l,is) = 0
          do ib = 1, nbin_rt
            qxi(l,is) = qxi(l,is) + qexti_dst (ib,is) * surf(ib)
            qsi(l,is) = qsi(l,is) + qscati_dst(ib,is) * surf(ib)
            gi (l,is) = gi (l,is) + gi_dst    (ib,is) * surf(ib)
          end do
          qsi(l,is) = min(qsi(l,is), 0.99999_r8 * qxi(l,is))
        end do
        qextrefdst(l) = qxv(l,nrefv)
        taurefdst (l) = Ao * qextrefdst(l) * (plev(l) - plev(l-1)) / g
        ! For diagnostics: dust opacity at reference wavelengths (vis and ir).
        taudst(1) = taudst(1) + taurefdst(l)
        taudst(2) = taudst(2) + taurefdst(l) / qextrefdst(l) * (qxi(l,4) - qsi(l,4))
      end do
    end if
  end do

end subroutine opt_dst