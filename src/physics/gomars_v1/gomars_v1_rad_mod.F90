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

module gomars_v1_rad_mod

  use gomars_v1_const_mod
  use gomars_v1_orbit_mod

  implicit none

  private

  public gomars_v1_rad_init
  public gomars_v1_rad_final
  public update_solar
  public qexti, qscati, gi, wi, fzeroi
  public qextv, qscatv, gv, wv, fzerov
  public solar_1au, solar, gweight
  public co2v, co2i, tgasref, pgasref, pfgasref, tauray
  public wrefco2, wrefh2o
  public wnov, dwnv, wavev
  public wnoi, dwni, wavei
  public planckir

  real(r8), public, parameter :: ubari = 0.5_r8
  ! If optical depth is less than this value, place the Gauss-point into the
  ! "zero" channel.
  real(r8), public, parameter :: tlimits = 1.0e-3_r8
  real(r8), public, parameter :: tlimiti = 5.0e-3_r8
  real(r8), public, parameter :: maxexp  = 35.0_r8
  real(r8), public, parameter :: taumax  = 35.0_r8

  real(r8), allocatable, dimension(:) :: wnoi
  real(r8), allocatable, dimension(:) :: dwni
  real(r8), allocatable, dimension(:) :: wavei
  real(r8), allocatable, dimension(:) :: wnov
  real(r8), allocatable, dimension(:) :: dwnv
  real(r8), allocatable, dimension(:) :: wavev
  real(r8), allocatable, dimension(:) :: solar_1au
  real(r8), allocatable, dimension(:) :: solar
  real(r8), allocatable, dimension(:) :: tauray

  real(r8), allocatable, dimension(:,:,:,:,:) :: co2i
  real(r8), allocatable, dimension(:,:,:,:,:) :: co2v

  real(r8), allocatable, dimension(:) :: fzeroi
  real(r8), allocatable, dimension(:) :: fzerov
  real(r8), allocatable, dimension(:) :: pgasref
  real(r8), allocatable, dimension(:) :: tgasref

  ! Put detau into physics state object.

  real(r8) qextref, ptop
  real(r8), allocatable, dimension(:) :: qextv
  real(r8), allocatable, dimension(:) :: qscatv
  real(r8), allocatable, dimension(:) :: wv
  real(r8), allocatable, dimension(:) :: gv
  real(r8), allocatable, dimension(:) :: qexti
  real(r8), allocatable, dimension(:) :: qscati
  real(r8), allocatable, dimension(:) :: wi
  real(r8), allocatable, dimension(:) :: gi

  real(r8), allocatable, dimension(:,:) :: planckir

  real(r8), allocatable, dimension(:) :: pfgasref

  real(r8), parameter :: gweight(ngauss) = [   &
    4.8083554740D-02, 1.0563099137D-01,        &
    1.4901065679D-01, 1.7227479710D-01,        &
    1.7227479710D-01, 1.4901065679D-01,        &
    1.0563099137D-01, 4.8083554740D-02,        &
    2.5307134073D-03, 5.5595258613D-03,        &
    7.8426661469D-03, 9.0670945845D-03,        &
    9.0670945845D-03, 7.8426661469D-03,        &
    5.5595258613D-03, 2.5307134073D-03,  0.0D0 &
  ]

  real(r8), parameter :: wrefco2(nrefh2o) = [      &
    9.999999D-1, 9.99999D-1, 9.9999D-1, 9.999D-1,   &
    9.99D-1, 9.9D-1, 9.0D-1, 8.0D-1, 7.0D-1, 6.0D-1 &
  ]

  real(r8), parameter :: wrefh2o(nrefh2o) = [      &
    1.0D-7, 1.0D-6, 1.0D-5, 1.0D-4, 1.0D-3, 1.0D-2, &
    1.0D-1, 2.0D-1, 3.0D-1, 4.0D-1                  &
  ]

contains

  subroutine gomars_v1_rad_init()

    call gomars_v1_rad_final()

    allocate(wnoi     (nspecti))
    allocate(dwni     (nspecti))
    allocate(wavei    (nspecti))
    allocate(wnov     (nspectv))
    allocate(dwnv     (nspectv))
    allocate(wavev    (nspectv))
    allocate(solar_1au(nspectv))
    allocate(solar    (nspectv))
    allocate(tauray   (nspectv))

    allocate(co2i(ntref,npint,nrefh2o,nspecti,ngauss))
    allocate(co2v(ntref,npint,nrefh2o,nspectv,ngauss))

    allocate(fzeroi(nspecti))
    allocate(fzerov(nspectv))

    allocate(pgasref(npref))
    allocate(tgasref(ntref))

    allocate(qextv (nspectv))
    allocate(qscatv(nspectv))
    allocate(wv    (nspectv))
    allocate(gv    (nspectv))
    allocate(qexti (nspecti))
    allocate(qscati(nspecti))
    allocate(wi    (nspecti))
    allocate(gi    (nspecti))

    allocate(planckir(nspecti,8501))

    allocate(pfgasref(npint))

    call setspv(wnov, dwnv, wavev, solar_1au, tauray)
    call setspi(wnoi, dwni, wavei, planckir)
    call setrad(tgasref, pfgasref, co2v, co2i, qextv, qscatv, wv, gv, &
                qexti, qscati, wi, gi, fzeroi, fzerov)

    ptop = 10**pfgasref(1)

  end subroutine gomars_v1_rad_init

  subroutine gomars_v1_rad_final()

    if (allocated(wnoi     )) deallocate(wnoi     )
    if (allocated(dwni     )) deallocate(dwni     )
    if (allocated(wavei    )) deallocate(wavei    )
    if (allocated(wnov     )) deallocate(wnov     )
    if (allocated(dwnv     )) deallocate(dwnv     )
    if (allocated(wavev    )) deallocate(wavev    )
    if (allocated(solar_1au)) deallocate(solar_1au)
    if (allocated(solar    )) deallocate(solar    )
    if (allocated(tauray   )) deallocate(tauray   )
    if (allocated(co2i     )) deallocate(co2i     )
    if (allocated(co2v     )) deallocate(co2v     )
    if (allocated(fzeroi   )) deallocate(fzeroi   )
    if (allocated(fzerov   )) deallocate(fzerov   )
    if (allocated(pgasref  )) deallocate(pgasref  )
    if (allocated(tgasref  )) deallocate(tgasref  )
    if (allocated(qextv    )) deallocate(qextv    )
    if (allocated(qscatv   )) deallocate(qscatv   )
    if (allocated(wv       )) deallocate(wv       )
    if (allocated(gv       )) deallocate(gv       )
    if (allocated(qexti    )) deallocate(qexti    )
    if (allocated(qscati   )) deallocate(qscati   )
    if (allocated(wi       )) deallocate(wi       )
    if (allocated(gi       )) deallocate(gi       )
    if (allocated(planckir )) deallocate(planckir )
    if (allocated(pfgasref )) deallocate(pfgasref )

  end subroutine gomars_v1_rad_final

  subroutine update_solar(ls)

    real(r8), intent(in) :: ls

    integer is
    real(r8) rsdist

    rsdist = solar_dist(ls)**2
    
    ! Calculate solar flux at the current Mars distance.
    do is = 1, nspectv
      solar(is) = solar_1au(is) * rsdist
    end do

  end subroutine update_solar

end module gomars_v1_rad_mod