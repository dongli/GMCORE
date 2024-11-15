subroutine ini_optdst(qextv, qscatv, grefv, qexti, qscati, grefi, qxv, qsv, gv, qxi, qsi, gi, qextref)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in) :: qextv (nspectv)
  real(r8), intent(in) :: qscatv(nspectv)
  real(r8), intent(in) :: grefv (nspectv)
  real(r8), intent(in) :: qexti (nspecti)
  real(r8), intent(in) :: qscati(nspecti)
  real(r8), intent(in) :: grefi (nspecti)
  real(r8), intent(out) :: qxv(2*nlev+4,nspectv)
  real(r8), intent(out) :: qsv(2*nlev+4,nspectv)
  real(r8), intent(out) :: gv (2*nlev+4,nspectv)
  real(r8), intent(out) :: qxi(2*nlev+4,nspecti)
  real(r8), intent(out) :: qsi(2*nlev+4,nspecti)
  real(r8), intent(out) :: gi (2*nlev+4,nspecti)
  real(r8), intent(out) :: qextref(2*nlev+4)

  integer i, k, n

  n = 2 * nlev + 4

  do k = 1, n
    qextref(k) = qextv(nrefv)
  end do

  do i = 1, nspectv
    do k = 1, n
      qxv(k,i) = qextv (i)
      qsv(k,i) = qscatv(i)
      gv (k,i) = grefv (i)
    end do
  end do

  do i = 1, nspecti
    do k = 1, n
      qxi(k,i) = qexti (i)
      qsi(k,i) = qscati(i)
      gi (k,i) = grefi (i)
    end do
  end do

end subroutine ini_optdst