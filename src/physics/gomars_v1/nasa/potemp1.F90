subroutine potemp1( &
  tstrat          , & 
  lnp             , &
  lnp_lev         , &
  pk_lev          , &
  pt              , &
  pt_lev          , &
  t_lev           )

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in ) :: tstrat
  real(r8), intent(in ) :: lnp    (nlev  )
  real(r8), intent(in ) :: lnp_lev(nlev+1)
  real(r8), intent(in ) :: pk_lev (nlev+1)
  real(r8), intent(in ) :: pt     (nlev  )
  real(r8), intent(out) :: pt_lev (nlev+1)
  real(r8), intent(out) :: t_lev  (nlev+1)

  integer k
  real(r8) dlnp1, dlnp2, slope
 
  do k = 2, nlev
    dlnp1 =          lnp_lev(k) - lnp(k-1)
    dlnp2 = lnp(k) - lnp_lev(k)
    slope = (pt(k) - pt(k-1)) / (dlnp1 + dlnp2)
    pt_lev(k) = pt(k-1) + slope * dlnp1
    t_lev (k) = pt_lev(k) * pk_lev(k)
  end do

  ! Model top level
  k = 1
  dlnp1 =          lnp_lev(k) - lnpstrat
  dlnp2 = lnp(k) - lnp_lev(k)
  slope = (pt(k) - tstrat / pstratk) / (dlnp1 + dlnp2)
  pt_lev(k) = tstrat / pstratk + slope * dlnp1
  t_lev (k) = pt_lev(k) * pk_lev(k)

  ! Model bottom level
  k = nlev + 1
  dlnp1 =            lnp_lev(k-1) - lnp(k-2)
  dlnp2 = lnp(k-1) - lnp_lev(k-1)
  slope = (pt(k-1) - pt(k-2)) / (dlnp1 + dlnp2)
  pt_lev(k) = pt(k-1) + slope * (lnp_lev(k) - lnp_lev(k-1))

end subroutine potemp1