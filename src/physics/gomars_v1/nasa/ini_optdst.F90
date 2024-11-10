subroutine ini_optdst(l_nspectv, l_nspecti, l_nrefv, nlev, qextv, qscatv, grefv, qexti, qscati, grefi, qxv, qsv, gv, qxi, qsi, gi, qextref)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod

  implicit none

  integer, intent(in) :: l_nspectv
  integer, intent(in) :: l_nspecti
  integer, intent(in) :: l_nrefv
  integer, intent(in) :: nlev
  real(r8), intent(in) :: qextv (l_nspectv)
  real(r8), intent(in) :: qscatv(l_nspectv)
  real(r8), intent(in) :: grefv (l_nspectv)
  real(r8), intent(in) :: qexti (l_nspecti)
  real(r8), intent(in) :: qscati(l_nspecti)
  real(r8), intent(in) :: grefi (l_nspecti)
  real(r8), intent(out) :: qxv(2*nlev+4,l_nspectv)
  real(r8), intent(out) :: qsv(2*nlev+4,l_nspectv)
  real(r8), intent(out) :: gv (2*nlev+4,l_nspectv)
  real(r8), intent(out) :: qxi(2*nlev+4,l_nspecti)
  real(r8), intent(out) :: qsi(2*nlev+4,l_nspecti)
  real(r8), intent(out) :: gi (2*nlev+4,l_nspecti)
  real(r8), intent(out) :: qextref(2*nlev+4)

  integer i, k, n

  n = 2 * nlev + 4

  do k = 1, n
    qextref(k) = qextv(l_nrefv)
  end do

  do i = 1, l_nspectv
    do k = 1, n
      qxv(k,i) = qextv (i)
      qsv(k,i) = qscatv(i)
      gv (k,i) = grefv (i)
    end do
  end do

  do i = 1, l_nspecti
    do k = 1, n
      qxi(k,i) = qexti (i)
      qsi(k,i) = qscati(i)
      gi (k,i) = grefi (i)
    end do
  end do

end subroutine ini_optdst