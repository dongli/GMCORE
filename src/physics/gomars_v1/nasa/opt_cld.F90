subroutine opt_cld(q, pl, qxv, qsv, gv, qxi, qsi, gi, qextrefcld, taurefcld, taucld)

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
  real(r8), intent(out) :: qextrefcld(2*nlev+4)
  real(r8), intent(out) :: taurefcld (2*nlev+4)
  real(r8), intent(out) :: taucld    (2)

  integer n, k, l, i, s, b
  integer irap
  real(r8) dev2, cst
  real(r8) mantletocore
  real(r8) Mo, No, Rs, Ao
  real(r8) surf(nbin_rt)

  n = 2 * nlev + 3

  do k = 1, n + 1
    qextrefcld(k) = 1
    taurefcld (k) = 0
  end do

  do s = 1, nspectv
    do k = 1, n + 1
      qxv(k,s) = qextrefcld(k)
      qsv(k,s) = qextrefcld(k) * 0.99_r8
      gv (k,s) = 0
    end do
  end do

  do s = 1, nspecti
    do k = 1, n + 1
      qxi(k,s) = qextrefcld(k)
      qsi(k,s) = qextrefcld(k) * 0.99_r8
      gi (k,s) = 0
    end do
  end do

  dev2 = 1.0_r8 / (sqrt(2.0_r8) * dev_ice)

  taucld = 0

  do l = 1, nlev
    if (q(l,iMa_cld) + q(l,iMa_cor) > 1.0e-7 .and. q(l,iNb_cld) > 1) then
      ! Determine the ratio of dust core over that of ice mantle.
      mantletocore = (q(l,iMa_cor) / dpden_dt / (q(l,iMa_cld) / dpden_ice + q(l,iMa_cor) / dpden_dt))**athird
      ! Find the index to which corresponds the optical properties of the core to mantle ratio.
      ! Those properties were determined off-line using the Toon and Ackerman coated spheres code.
      irap = nratio
      do i = 1, nratio
        if (mantletocore < cor_ratio(i) .and. i /= 1) then
          irap = i - 1
          exit
        else if (mantletocore == cor_ratio(i)) then
          irap = i
          exit
        else if (mantletocore < cor_ratio(i) .and. i == 1) then
          irap = i
          exit
        end if
      end do
      ! Calculate the cross-section mean radius (Rs) of the log-normal distribution.
      Mo = q(l,iMa_cld) + q(l,iMa_cor)
      No = q(l,iNb_cld)
      cst = 0.75_r8 / (pi * (q(l,iMa_cld) / dpden_ice + q(l,iMa_cor) / dpden_dt) / Mo)
      Rs = (Mo / No * cst)**athird * exp(-0.5_r8 * dev_ice**2)
      ! Calculate the total cross sectional area (Ao) of water ice particles.
      Ao = No * pi * Rs**2
      ! Define the cross-section weighted distribution, i.e. surface/size bin. Change Rs to Reff.
      Rs = Rs * exp(1.5_r8 * dev_ice**2)
      Rs = 1.0_r8 / min(max(Rs, 1.0e-7_r8), 100.0e-6_r8)
      do b = 1, nbin_rt
        surf(b) = 0.5_r8 * (erf(log(radb_rt(b+1) * Rs) * dev2) - erf(log(radb_rt(b) * Rs) * dev2))
      end do
      ! Calculate the averaged values of the optical properties for the whole distribution.
      do k = 2 * l + 2, 2 * l + 3
        do s = 1, nspectv
          qxv(k,s) = 0
          qsv(k,s) = 0
          do b = 1, nbin_rt
            qxv(k,s) = qxv(k,s) + surf(b) * qextv_cld (irap,b,s)
            qsv(k,s) = qsv(k,s) + surf(b) * qscatv_cld(irap,b,s)
            gv (k,s) = gv (k,s) + surf(b) * gv_cld    (irap,b,s)
          end do
          qsv(k,s) = min(qsv(k,s), 0.99999_r8 * qxv(k,s))
        end do
        do s = 1, nspecti
          qxi(k,s) = 0
          qsi(k,s) = 0
          do b = 1, nbin_rt
            qxi(k,s) = qxi(k,s) + surf(b) * qexti_cld (irap,b,s)
            qsi(k,s) = qsi(k,s) + surf(b) * qscati_cld(irap,b,s)
            gi (k,s) = gi (k,s) + surf(b) * gi_cld    (irap,b,s)
          end do
          qsi(k,s) = min(qsi(k,s), 0.99999_r8 * qxi(k,s))
        end do
        qextrefcld(k) = qxv(k,nrefv)
        taurefcld (k) = Ao * qextrefcld(k) * (pl(k) - pl(k-1)) / g
        ! For diagnostics: cloud opacity at reference wavelengths (vis and ir).
        taucld(1) = taucld(1) + taurefcld(k)
        taucld(2) = taucld(2) + taurefcld(k) / qextrefcld(k) * (qxi(k,4) - qsi(k,4))
      end do
    end if
  end do

end subroutine opt_cld