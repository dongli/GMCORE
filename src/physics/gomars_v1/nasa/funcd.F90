subroutine funcd( &
  astar, downir, rhouch, rhoucht, &
  scond, stemp, sthick, &
  tg, ps, q_vap_sfc, gndice, npcflag, &
  f, df)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_namelist_mod

  implicit none

  real(r8), intent(in ) :: astar
  real(r8), intent(in ) :: downir
  real(r8), intent(in ) :: rhouch
  real(r8), intent(in ) :: rhoucht
  real(r8), intent(in ) :: scond
  real(r8), intent(in ) :: stemp
  real(r8), intent(in ) :: sthick
  real(r8), intent(in ) :: tg
  real(r8), intent(in ) :: ps
  real(r8), intent(in ) :: q_vap_sfc
  real(r8), intent(in ) :: gndice
  logical , intent(in ) :: npcflag
  real(r8), intent(out) :: f
  real(r8), intent(out) :: df

  real(r8) qg

  f = astar + downir + rhoucht - rhouch * tg + 2 * scond * (stemp - tg) / sthick - stbo * tg**4
  df = -rhouch - 2 * scond / sthick - 4 * stbo * tg**3

  if (latent_heat) then
    if (gndice > 0 .or. npcflag) then
      qg = (18.0_r8 / 44.0_r8) * 611 * exp(22.5_r8 * (1 - (t0_trip / tg))) / ps
      f = f + rhouch * (2.8e6_r8 / cpd) * (q_vap_sfc - qg)
      df = df - 6146.1_r8 * rhouch * (2.8e6_r8 / cpd) * qg / tg**2
    end if
  end if

end subroutine funcd