subroutine microphys(dt_mp, ps, p, dp, t, taux, tauy, ht_pbl, p_pbl_top)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in ) :: dt_mp
  real(r8), intent(in ) :: ps
  real(r8), intent(in ) :: p       (nlev)
  real(r8), intent(in ) :: dp      (nlev)
  real(r8), intent(in ) :: t       (nlev)
  real(r8), intent(in ) :: taux
  real(r8), intent(in ) :: tauy
  real(r8), intent(in ) :: ht_pbl
  real(r8), intent(in ) :: p_pbl_top



end subroutine microphys

subroutine dust_update(dt_mp, ps, p, dp, t, taux, tauy, ht_pbl, p_pbl_top)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_namelist_mod

  implicit none

  real(r8), intent(in ) :: dt_mp
  real(r8), intent(in ) :: ps
  real(r8), intent(in ) :: p       (nlev)
  real(r8), intent(in ) :: dp      (nlev)
  real(r8), intent(in ) :: t       (nlev)
  real(r8), intent(in ) :: taux
  real(r8), intent(in ) :: tauy
  real(r8), intent(in ) :: ht_pbl
  real(r8), intent(in ) :: p_pbl_top

  if (wsl_newman) then

  else if (wsl_kmh) then

  end if

end subroutine dust_update