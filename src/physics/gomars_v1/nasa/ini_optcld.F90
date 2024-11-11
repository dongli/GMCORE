subroutine ini_optcld(l_nspectv, l_nspecti, qxv, qsv, gv, qxi, qsi, gi, qextrefcld, taurefcld)
  
  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod

  implicit none

  integer, intent(in) :: l_nspectv
  integer, intent(in) :: l_nspecti
  real(r8), intent(out) :: qxv(2*nlev+4,l_nspectv)
  real(r8), intent(out) :: qsv(2*nlev+4,l_nspectv)
  real(r8), intent(out) :: gv (2*nlev+4,l_nspectv)
  real(r8), intent(out) :: qxi(2*nlev+4,l_nspecti)
  real(r8), intent(out) :: qsi(2*nlev+4,l_nspecti)
  real(r8), intent(out) :: gi (2*nlev+4,l_nspecti)
  real(r8), intent(out) :: qextrefcld(2*nlev+4)
  real(r8), intent(out) :: taurefcld (2*nlev+4)

  integer i, k, n
  
  n = 2 * nlev + 4
  
  do k = 1, n
    qextrefcld(k) = 1
    taurefcld (k) = 0
  end do
  
  do i = 1, l_nspectv
    do k = 1, n
      qxv(k,i) = 0
      qsv(k,i) = 0
      gv (k,i) = 0
    end do
  end do
  
  do i = 1, l_nspecti
    do k = 1, n
      qxi(k,i) = 0
      qsi(k,i) = 0
      gi (k,i) = 0
    end do
  end do

end subroutine ini_optcld