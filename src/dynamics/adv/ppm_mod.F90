! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module ppm_mod

  use const_mod
  use namelist_mod
  use limiter_mod

  implicit none

  private

  public ppm1
  public ppm2
  public ppm3

contains

  subroutine ppm1(fm2, fm1, f, fp1, fp2, fl, df, f6)

    real(r8), intent(in ) :: fm2
    real(r8), intent(in ) :: fm1
    real(r8), intent(in ) :: f
    real(r8), intent(in ) :: fp1
    real(r8), intent(in ) :: fp2
    real(r8), intent(out) :: fl
    real(r8), intent(out) :: df
    real(r8), intent(out) :: f6

    real(r8) dfl, dfr, fr

    ! Calculate values at left and right cell interfaces.
    dfl = slope(fm2, fm1, f  )
    df  = slope(fm1, f  , fp1)
    dfr = slope(f  , fp1, fp2)
    ! Why (B2) in Lin (2004) divide (dfl - df) and (df - dfr) by 3?
    fl = 0.5_r8 * (fm1 + f) + (dfl - df) / 6.0_r8
    fr = 0.5_r8 * (fp1 + f) + (df - dfr) / 6.0_r8
    ! Why (B3) and (B4) in Lin (2004) multiply df by 2?
    fl = f - sign(min(abs(df), abs(fl - f)), df)
    fr = f + sign(min(abs(df), abs(fr - f)), df)
    f6 = 6 * f - 3 * (fl + fr)
    df = fr - fl

  end subroutine ppm1

  pure real(r8) function ppm2(cfl, fm2, fm1, f, fp1, fp2) result(res)

    real(r8), intent(in) :: cfl
    real(r8), intent(in) :: fm2
    real(r8), intent(in) :: fm1
    real(r8), intent(in) :: f
    real(r8), intent(in) :: fp1
    real(r8), intent(in) :: fp2

    real(r8), parameter :: c1 = -1.0_r8 / 12.0_r8
    real(r8), parameter :: c2 =  7.0_r8 / 12.0_r8

    real(r8) fl, fr, tau

    fl = c1 * (fm2 + fp1) + c2 * (fm1 + f)
    fr = c1 * (fm1 + fp2) + c2 * (f + fp1)

    tau = (4 * fl + 2 * fr - 6 * f) / (6 * fl + 6 * fr - 12 * f)
    if (tau * (1 - tau) > 0) then
      fl = f
      fr = f
    end if

    if ((fl - fm1) * (f - fl) < 0) then
      fl = min(max(f, fm1), max(fl, min(f, fm1)))
    end if
    if ((fr - f) * (fp1 - fr) < 0) then
      fr = min(max(fp1, f), max(fr, min(fp1, f)))
    end if

    if (cfl >= 0) then
      res = (cfl**2 - 2 * cfl + 1) * fr - (2 * cfl**2 - 3 * cfl) * f + (cfl**2 - cfl) * fl
    else
      res = (cfl**2 + 2 * cfl + 1) * fl - (2 * cfl**2 + 3 * cfl) * f + (cfl**2 + cfl) * fr
    end if

    ! tau = (2 * fl + fr - 3 * f) / (3 * fl + 3 * fr - 6 * f)
    ! if (tau * (1 - tau) > 0) then
    !   res = f
    ! end if

  end function ppm2

  pure real(r8) function ppm3(cfl, fm2, fm1, f, fp1, fp2) result(res)

    real(r8), intent(in) :: cfl
    real(r8), intent(in) :: fm2
    real(r8), intent(in) :: fm1
    real(r8), intent(in) :: f
    real(r8), intent(in) :: fp1
    real(r8), intent(in) :: fp2

    real(r8) dfl, df, dfr, fl, fr, tau

    ! Calculate values at left and right cell interfaces.
    dfl = slope(fm2, fm1, f  )
    df  = slope(fm1, f  , fp1)
    dfr = slope(f  , fp1, fp2)
    ! Why (B2) in Lin (2004) divide (dfl - df) and (df - dfr) by 3?
    fl = 0.5_r8 * (fm1 + f) + (dfl - df) / 6.0_r8
    fr = 0.5_r8 * (fp1 + f) + (df - dfr) / 6.0_r8
    ! Why (B3) and (B4) in Lin (2004) multiply df by 2?
    fl = f - sign(min(abs(df), abs(fl - f)), df)
    fr = f + sign(min(abs(df), abs(fr - f)), df)

    if (cfl >= 0) then
      res = (cfl**2 - 2 * cfl + 1) * fr - (2 * cfl**2 - 3 * cfl) * f + (cfl**2 - cfl) * fl
    else
      res = (cfl**2 + 2 * cfl + 1) * fl - (2 * cfl**2 + 3 * cfl) * f + (cfl**2 + cfl) * fr
    end if

    ! tau = (2 * fl + fr - 3 * f) / (3 * fl + 3 * fr - 6 * f)
    ! if (tau * (1 - tau) > 0) then
    !   res = f
    ! end if

  end function ppm3

end module ppm_mod
