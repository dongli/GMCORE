! ==============================================================================
! This file is part of GoMars since 2023.
!
! GoMars is a Martian general circulation model developed in Institute of
! Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
!
! GMCORE is a dynamical core for atmospheric model used in GoMars.
!
! GoMars and GMCORE are distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   CO2 in the atmosphere and polar caps acts as two reservoirs.
!
! ==============================================================================

module gomars_v2_lsm_mod

  use gomars_v2_const_mod
  use gomars_v2_namelist_mod
  use gomars_v2_types_mod
  use gomars_v2_tracers_mod
  use gomars_v2_objects_mod

  implicit none

  private

  public gomars_v2_lsm_init
  public gomars_v2_lsm_final
  public gomars_v2_lsm_run

  real(r8), allocatable, dimension(:) :: soil_z
  real(r8), allocatable, dimension(:) :: soil_z_lev
  real(r8), allocatable, dimension(:) :: soil_dz
  real(r8), allocatable, dimension(:) :: soil_dz_lev

  real(r8), parameter :: factl = 0.25_r8
  real(r8), parameter :: factm = 1.2_r8
  ! Diurnal skin depth (m) corresponding to a surface thermal inertia of 336 SI units
  real(r8), parameter :: skind = 0.06_r8
  ! Ground ice depth in the northern hemisphere (m).
  real(r8), parameter :: gidn  = 0.0545_r8
  ! Ground ice depth in the southern hemisphere (m).
  real(r8), parameter :: gids  = 0.0805_r8

contains

  subroutine gomars_v2_lsm_init()

    integer iblk, icol, k

    call gomars_v2_lsm_final()

    allocate(soil_z     (nlev_soil  ))
    allocate(soil_z_lev (nlev_soil+1))
    allocate(soil_dz    (nlev_soil  ))
    allocate(soil_dz_lev(nlev_soil+1))

    ! Setup soil levels.
    soil_dz(1) = factl * skind
    do k = 2, nlev_soil
      soil_dz(k) = factm * soil_dz(k-1)
    end do
    soil_z_lev(1) = 0
    do k = 2, nlev_soil+1
      soil_z_lev(k) = soil_z_lev(k-1) + soil_dz(k-1)
    end do
    do k = 1, nlev_soil
      soil_z(k) = 0.5_r8 * (soil_z_lev(k) + soil_z_lev(k+1))
    end do
    do k = 2, nlev_soil
      soil_dz_lev(k) = 0.5_r8 * (soil_dz(k-1) + soil_dz(k))
    end do

    ! Initialize soil conductivity.
    do iblk = 1, size(objects)
      associate (mesh          => objects(iblk)%mesh               , &
                 gnd_ice       => objects(iblk)%state%gnd_ice      , &
                 tin           => objects(iblk)%state%tin          , &
                 soil_tin      => objects(iblk)%state%soil_tin     , &
                 soil_rho      => objects(iblk)%state%soil_rho     , &
                 soil_cp       => objects(iblk)%state%soil_cp      , &
                 soil_cond     => objects(iblk)%state%soil_cond    , &
                 soil_cond_lev => objects(iblk)%state%soil_cond_lev)
      do icol = 1, objects(iblk)%mesh%ncol
        ! Set soil parameters.
        if (gnd_ice(icol) > 0.5_r8) then
          if (mesh%lat(icol) > 0) then
            do k = 1, nlev_soil
              soil_tin(icol,k) = 2236.995_r8
            end do
          else
            do k = 1, nlev_soil
              soil_tin(icol,k) = 1100.0_r8
            end do
          end if
          do k = 1, nlev_soil
            soil_rho (icol,k) = 1781.99_r8
            soil_cp  (icol,k) = 1404.09_r8
          end do
        else
          do k = 1, nlev_soil
            soil_rho (icol,k) = 1481.39_r8 ! 1500.0_r8
            soil_cp  (icol,k) = 840.0_r8   ! 627.9_r8
          end do
        end if
        do k = 1, nlev_soil
          soil_cond(icol,k) = tin(icol)**2 / (soil_rho(icol,k) * soil_cp(icol,k))
        end do
        soil_cond_lev(icol,1) = soil_cond(icol,1)
        do k = 2, nlev_soil - 1
          soil_cond_lev(icol,k) = 0.5_r8 * (soil_cond(icol,k-1) + soil_cond(icol,k))
        end do
      end do
      end associate
    end do

  end subroutine gomars_v2_lsm_init

  subroutine gomars_v2_lsm_final()

    if (allocated(soil_z     )) deallocate(soil_z     )
    if (allocated(soil_z_lev )) deallocate(soil_z_lev )
    if (allocated(soil_dz    )) deallocate(soil_dz    )
    if (allocated(soil_dz_lev)) deallocate(soil_dz_lev)

  end subroutine gomars_v2_lsm_final

  subroutine gomars_v2_lsm_run(state, tend)

    type(gomars_v2_state_type ), intent(inout), target :: state
    type(gomars_v2_tend_type  ), intent(inout) :: tend

    integer icol
    real(r8), pointer :: qice(:)
    real(r8) tsat

    associate (mesh      => state%mesh     , &
               t_bot     => state%t_bot    , & ! in
               ps        => state%ps       , & ! in
               co2ice    => state%co2ice   , & ! inout
               q_sfc     => state%q_sfc    , & ! in
               fdnl      => state%fdnl     , & ! in
               fdns      => state%fdns     , & ! in
               rhouch    => state%rhouch   , & ! in
               soil_cond => state%soil_cond, & ! in
               soil_t    => state%soil_t   , & ! in
               heat_sfc  => state%heat_sfc , & ! in
               gnd_ice   => state%gnd_ice  , & ! in
               alb       => state%alb      , & ! out
               tg        => state%tg       , & ! out
               dpsdt     => tend%dpsdt     )   ! out
    qice => state%q(:,mesh%nlev,idx_m_vap)
    do icol = 1, state%mesh%ncol
      ! Calculate surface albedo.
      if (co2ice(icol) > 0) then
        alb(icol) = merge(alb_ice_np, alb_ice_sp, mesh%lat(icol) > 0)
      else if (albedo_feedback .and. q_sfc(icol,idx_m_vap) > ice_thresh_kgm2) then
        alb(icol) = ice_albedo
      else
        alb(icol) = state%alb(icol)
      end if
      ! Calculate CO2 frost point at this surface pressure (hPa).
      tsat = 3182.48_r8 / (23.3494_r8 - log(ps(icol) * 0.01))
      if (co2ice(icol) <= 0) then ! No CO2 on the ground
        co2ice(icol) = 0
        ! Calculate surface temperature (K) from surface energy balance.
        call calc_tg(alb(icol), fdnl(icol), fdns(icol), rhouch(icol), t_bot(icol), &
          soil_cond(icol,1), soil_t(icol,1), soil_dz(icol), ps(icol), qice(icol), &
          gnd_ice(icol), tg(icol))
      end if
    end do
    end associate
    stop 999

  end subroutine gomars_v2_lsm_run

  subroutine calc_tg(alb, fdnl, fdns, rhouch, t_bot, soil_cond, soil_t, soil_dz, ps, qice, gnd_ice, tg)

    real(r8), intent(in) :: alb
    real(r8), intent(in) :: fdnl
    real(r8), intent(in) :: fdns
    real(r8), intent(in) :: rhouch
    real(r8), intent(in) :: t_bot
    real(r8), intent(in) :: soil_cond
    real(r8), intent(in) :: soil_t
    real(r8), intent(in) :: soil_dz
    real(r8), intent(in) :: ps
    real(r8), intent(in) :: qice
    real(r8), intent(in) :: gnd_ice
    real(r8), intent(out) :: tg

    integer, parameter :: max_iter = 30
    real(r8) fl, fh, f, df, tl, th, dt, dt_save
    integer iter

    !
    ! f(Tsfc) =  Fl↓ + (1 - α)Fs↓ - Fconv - Fcond - εσTsfc⁴
    !
    ! Upward heat exchange with the atmosphere (W m⁻²)
    ! Fconv = ρ u* cp ch (Tsfc - Tbot)
    !
    ! Downward heat exchange with the ground (W m⁻²)
    ! Fcond = -K (Tsoil - Tsfc) / (Δz / 2)
    !
    ! f'(T) =  -4εσT³ - 2 K / Δz - ρ u* cp ch
    !
    tl = 50.0_r8
    th = 350.0_r8
    call surface_energy_balance(tl, fl, df)
    call surface_energy_balance(th, fh, df)
    if (fl == 0) then
      tg = tl
      return
    else if (fh == 0) then
      tg = th
      return
    else if (fl > 0) then
      f = tl; tl = th; th = f
    end if
    tg = 0.5 * (tl + th)
    dt = th - tl
    dt_save = dt
    call surface_energy_balance(tg, f, df)
    do iter = 1, max_iter
      if (((tg - th) * df - f) * ((tg - tl) * df - f) >= 0 .or. &
          abs(2 * f) > abs(dt_save * df)) then
        dt_save = dt
        dt = 0.5_r8 * (th - tl)
        tg = tl + dt
      else
        dt_save = dt
        dt = f / df
        tg = tg - dt
      end if
      if (dt == 0 .or. abs(dt) < 1.0E-3) return
      call surface_energy_balance(tg, f, df)
      if (f < 0) then
        tl = tg
      else
        th = tg
      end if
    end do

    ! Iteration failed to converge.
    stop 'calc_t_sfc: Iteration failed to converge.'

  contains

    subroutine surface_energy_balance(t_sfc, f, df)

      real(r8), intent(in) :: t_sfc ! Surface temperature (K)
      real(r8), intent(out) :: f    ! Energy residual
      real(r8), intent(out) :: df   ! Energy residual derivative respect to t_sfc

      f = (1 - alb) * fdns                           & ! input energy
        + emis_gnd_15um * fdnl                       & ! input energy
        + rhouch * (t_bot - t_sfc)                   & ! input energy
        + 2 * soil_cond * (soil_t - t_sfc) / soil_dz & ! input energy
        - stbo * t_sfc**4                              ! output energy
      df = -rhouch - 2 * soil_cond / soil_dz - 4 * stbo * t_sfc**3

      if (water_ice_latent_heat) then
        ! Latent heat of water ice
      end if

    end subroutine surface_energy_balance

  end subroutine calc_tg

end module gomars_v2_lsm_mod