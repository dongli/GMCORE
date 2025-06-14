! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module hybrid_coord_mod

  use flogger
  use const_mod, only: r8
  use namelist_mod, p0 => hybrid_coord_p0
  use hybrid_coord_ncep_mod
  use hybrid_coord_wrf_mod
  use hybrid_coord_test_mod
  use hybrid_coord_ecmwf_mod
  use mars_vert_coord_mod
  use latlon_mesh_mod
  use latlon_parallel_types_mod

  implicit none

  private

  public hybrid_coord_init
  public hybrid_coord_final
  public hybrid_coord_calc_mg
  public hybrid_coord_calc_mg_lev
  public hybrid_coord_calc_dmgdt
  public hybrid_coord_calc_dmgdt_lev
  public hybrid_coord_calc_ddmgdt
  public hyai, hybi, hyam, hybm

  real(r8), allocatable, target, dimension(:) :: hyai
  real(r8), allocatable, target, dimension(:) :: hybi
  real(r8), allocatable, target, dimension(:) :: hyam
  real(r8), allocatable, target, dimension(:) :: hybm

  namelist /hybrid_coord/ &
    hyai, hybi, hyam, hybm

  real(r8) :: local_ptop = 0

contains

  subroutine hybrid_coord_init(namelist_file, template)

    character(*), intent(in), optional :: namelist_file
    character(*), intent(in), optional :: template

    integer k, ierr

    if (allocated(hyai)) deallocate(hyai)
    if (allocated(hybi)) deallocate(hybi)
    if (allocated(hyam)) deallocate(hyam)
    if (allocated(hybm)) deallocate(hybm)

    allocate(hyai(nlev+1)); hyai = 0
    allocate(hybi(nlev+1)); hybi = 0
    allocate(hyam(nlev  )); hyam = 0
    allocate(hybm(nlev  )); hybm = 0

    if (present(namelist_file)) then
      open(10, file=namelist_file, status='old')
      read(10, nml=hybrid_coord, iostat=ierr)
      close(10)
    end if

    if (present(template)) then
      select case (template)
      case ('test_l15')
        call hybrid_coord_test_l15(p0, ptop, hyai, hybi)
      case ('test_l15_2')
        call hybrid_coord_test_l15_2(p0, ptop, hyai, hybi)
      case ('test_l26')
        call hybrid_coord_test_l26(p0, ptop, hyai, hybi)
      case ('test_l30')
        call hybrid_coord_test_l30(p0, ptop, hyai, hybi)
      case ('ecmwf_l50')
        call hybrid_coord_ecmwf_l50(p0, ptop, hyai, hybi)
      case ('ecmwf_l90')
        call hybrid_coord_ecmwf_l90(p0, ptop, hyai, hybi)
      case ('wrf_l32')
        call hybrid_coord_wrf_l32(p0, ptop, hyai, hybi)
        local_ptop = ptop
      case ('wrf_l64')
        call hybrid_coord_wrf_l64(p0, ptop, hyai, hybi)
        local_ptop = ptop
      case ('supercell_l40')
        call hybrid_coord_supercell_l40(p0, ptop, hyai, hybi)
      case ('schar_l40')
        call hybrid_coord_schar_l40(p0, ptop, hyai, hybi)
      case ('schar_l60')
        call hybrid_coord_schar_l60(p0, ptop, hyai, hybi)
      case ('dcmip21_l60')
        call hybrid_coord_dcmip21_l60(p0, ptop, hyai, hybi)
      case ('dcmip31_l10')
        call hybrid_coord_dcmip31_l10(p0, ptop, hyai, hybi)
      case ('modon_l5')
        call hybrid_coord_modon_l5(p0, ptop, hyai, hybi)
      case ('waccm_l70')
        call hybrid_coord_waccm_l70(p0, ptop, hyai, hybi)
      case ('emars28')
        call mars_vert_coord_emars28(p0, ptop, hyai, hybi)
      case ('dcmip12')
        call hybrid_coord_dcmip12(p0, ptop, hyai, hybi)
      case ('cam_l32')
        call hybrid_coord_cam_l32(p0, ptop, hyai, hybi)
      case ('cam_l30')
        call hybrid_coord_cam_l30(p0, ptop, hyai, hybi)
      case ('ncep')
        call hybrid_coord_ncep(p0, ptop, hyai, hybi)
      case ('wrf')
        call hybrid_coord_wrf(p0, ptop, hyai, hybi)
        local_ptop = ptop
      case default
        if (baroclinic .and. proc%is_root()) then
          call log_error('Hybrid vertical coordinate template "' // trim(template) // '" is invalid!')
        end if
      end select
    end if

    if ((baroclinic .or. advection) .and. global_mesh%full_nlev > 1 .and. (all(hyai == 0) .and. all(hybi == 0))) then
      if (proc%is_root()) call log_error('Hybrid coordinate parameters are not set!', __FILE__, __LINE__)
    end if

    if (all(hyam == 0)) then
      do k = 1, nlev
        hyam(k) = 0.5d0 * (hyai(k) + hyai(k+1))
        hybm(k) = 0.5d0 * (hybi(k) + hybi(k+1))
      end do
    end if

    do k = 1, nlev
      global_mesh%full_lev(k) = hyam(k) + hybm(k)
    end do
    do k = 1, nlev + 1
      global_mesh%half_lev(k) = hyai(k) + hybi(k)
    end do

  end subroutine hybrid_coord_init

  subroutine hybrid_coord_final()

    deallocate(hyai)
    deallocate(hybi)
    deallocate(hyam)
    deallocate(hybm)

  end subroutine hybrid_coord_final

  pure real(r8) function hybrid_coord_calc_mg(k, mgs, ref_ps_perb) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: mgs
    real(r8), intent(in), optional :: ref_ps_perb

    res = hyam(k) * (p0 - local_ptop) + hybm(k) * (mgs - local_ptop) + local_ptop

  end function hybrid_coord_calc_mg

  pure real(r8) function hybrid_coord_calc_mg_lev(k, mgs, ref_ps_perb) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: mgs
    real(r8), intent(in), optional :: ref_ps_perb

    res = hyai(k) * (p0 - local_ptop) + hybi(k) * (mgs - local_ptop) + local_ptop

  end function hybrid_coord_calc_mg_lev

  pure real(r8) function hybrid_coord_calc_dmgdt(k, dmgsdt) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: dmgsdt

    res = hybm(k) * dmgsdt

  end function hybrid_coord_calc_dmgdt

  pure real(r8) function hybrid_coord_calc_dmgdt_lev(k, dmgsdt) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: dmgsdt

    res = hybi(k) * dmgsdt

  end function hybrid_coord_calc_dmgdt_lev

  pure real(r8) function hybrid_coord_calc_ddmgdt(k, dmgsdt) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: dmgsdt

    res = (hybi(k+1) - hybi(k)) * dmgsdt

  end function hybrid_coord_calc_ddmgdt

end module hybrid_coord_mod
