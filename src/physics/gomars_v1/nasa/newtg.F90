subroutine newtg( &
  als, dnvflux, downir, rhouch, rhoucht, &
  scond, stemp, sthick, ps, q_vap_sfc, gndice, &
  polarcap, subflux, tg)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! Use modified Newton-Raphson method to find the ground temperature.

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in   ) :: als
  real(r8), intent(in   ) :: dnvflux
  real(r8), intent(in   ) :: downir
  real(r8), intent(in   ) :: rhouch
  real(r8), intent(in   ) :: rhoucht
  real(r8), intent(in   ) :: scond
  real(r8), intent(in   ) :: stemp
  real(r8), intent(in   ) :: sthick
  real(r8), intent(in   ) :: ps
  real(r8), intent(in   ) :: q_vap_sfc
  real(r8), intent(inout) :: gndice
  logical , intent(in   ) :: polarcap
  real(r8), intent(  out) :: subflux
  real(r8), intent(  out) :: tg

  integer , parameter :: max_iter = 30
  real(r8), parameter :: t1 = 50
  real(r8), parameter :: t2 = 350
  integer iter
  logical done
  real(r8) astar
  real(r8) wflux
  real(r8) qg
  real(r8) f, fl, fh, df
  real(r8) tl, th
  real(r8) dx_old, dx

  done = .false.
  astar = (1 - als) * dnvflux

  call funcd(astar, downir, rhouch, rhoucht, scond, stemp, sthick, t1, ps, q_vap_sfc, gndice, polarcap, fl, df)
  call funcd(astar, downir, rhouch, rhoucht, scond, stemp, sthick, t2, ps, q_vap_sfc, gndice, polarcap, fh, df)

  if (fl == 0) then
    tg = t1
    done = .true.
  else if (fh == 0) then
    tg = t2
    done = .true.
  else if (fl < 0) then
    tl = t1
    th = t2
  else
    tl = t2
    th = t1
  end if

  if (.not. done) then
    tg = 0.5_r8 * (t1 + t2)
    dx_old = abs(t2 - t1)
    dx = dx_old

    call funcd(astar, downir, rhouch, rhoucht, scond, stemp, sthick, tg, ps, q_vap_sfc, gndice, polarcap, f, df)

    do iter = 1, max_iter
      if (((tg - th) * df - f) * ((tg - tl) * df - f) >= 0 .or. abs(2 * f) > abs(dx_old * df)) then
        dx_old = dx
        dx = 0.5 * (th - tl)
        tg = tl + dx
        if (tl == tg) then
          done = .true.
          exit
        end if
      else
        dx_old = dx
        dx = f / df
        if (tg == tg - dx) then
          done = .true.
          exit
        else
          tg = tg - dx
        end if
      end if
      if (abs(dx) < 1.0e-3) then
        done = .true.
        exit
      end if
      call funcd(astar, downir, rhouch, rhoucht, scond, stemp, sthick, tg, ps, q_vap_sfc, gndice, polarcap, f, df)
      if (f < 0) then
        tl = tg
      else
        th = tg
      end if
    end do
  end if

  if (done) then
    wflux = 0
    if (gndice > 0 .and. .not. polarcap) then
      qg = 611 * exp(22.5_r8 * (1 - (t0_trip / tg))) / ps / mwratio
      wflux = -rhouch * (q_vap_sfc - qg) / cpd
      if (wflux * dt >= gndice) then
        wflux = gndice / dt
        gndice = 0
      else
        gndice = gndice - wflux * dt
      end if
    else if (polarcap) then
      qg = 611 * exp(22.5_r8 * (1 - (t0_trip / tg))) / ps / mwratio
      wflux = -rhouch * (q_vap_sfc - qg) / cpd
      gndice = gndice - wflux * dt
    end if
  else
    write(*, *) 'Subroutine newtg did not converge!'
    stop 1
  end if

  subflux = subflux + wflux * dt

end subroutine newtg