subroutine opt_dst(q, pl, qxv, qsv, gv, qxi, qsi, gi, qextrefdst, taurefdst, taudst)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_tracers_mod
  use gomars_v1_mp_mod

  implicit none

  real(r8), intent(in ) :: q         (nlev,ntracers)
  real(r8), intent(in ) :: pl        (2*nlev+3)
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

  do k = 1, n + 1
    qextrefdst(k) = 1
    taurefdst (k) = 0
  end do

  do is = 1, nspectv
    do k = 1, n + 1
      qxv(k,is) = qextrefdst(k)
      qsv(k,is) = qextrefdst(k) * 0.99_r8
      gv (k,is) = 0
    end do
  end do
  do is = 1, nspecti
    do k = 1, n + 1
      qxi(k,is) = qextrefdst(k)
      qsi(k,is) = qextrefdst(k) * 0.99_r8
      gi (k,is) = 0
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
      do ib = 1, nbin_rt
        surf(ib) = 0.5_r8 * (erf(log(radb_rt(ib+1) * Rs) * dev2) - erf(log(radb_rt(ib) * Rs) * dev2))
      end do
      ! Calculate the averaged values of the optical properties for the whole distribution.
      do k = 2 * l + 2, 2 * l + 3
        do is = 1, nspectv
          qxv(k,is) = 0
          qsv(k,is) = 0
          do ib = 1, nbin_rt
            qxv(k,is) = qxv(k,is) + qextv_dst (ib,is) * surf(ib)
            qsv(k,is) = qsv(k,is) + qscatv_dst(ib,is) * surf(ib)
            gv (k,is) = gv (k,is) + gv_dst    (ib,is) * surf(ib)
          end do
          qsv(k,is) = min(qsv(k,is), 0.99999_r8 * qxv(k,is))
        end do
        do is = 1, nspecti
          qxi(k,is) = 0
          qsi(k,is) = 0
          do ib = 1, nbin_rt
            qxi(k,is) = qxi(k,is) + qexti_dst (ib,is) * surf(ib)
            qsi(k,is) = qsi(k,is) + qscati_dst(ib,is) * surf(ib)
            gi (k,is) = gi (k,is) + gi_dst    (ib,is) * surf(ib)
          end do
          qsi(k,is) = min(qsi(k,is), 0.99999_r8 * qxi(k,is))
        end do
        qextrefdst(k) = qxv(k,nrefv)
        taurefdst (k) = Ao * qextrefdst(k) * (pl(k) - pl(k-1)) / g
        ! For diagnostics: dust opacity at reference wavelengths (vis and ir).
        taudst (1) = taudst(1) + taurefdst(k)
        taudst (2) = taudst(2) + taurefdst(k) / qextrefdst(k) * (qxi(k,4) - qsi(k,4))
      end do
    end if
  end do

end subroutine opt_dst