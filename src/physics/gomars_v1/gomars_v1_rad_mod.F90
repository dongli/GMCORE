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

  implicit none

  private

  public gomars_v1_rad_init
  public gomars_v1_rad_final
  public qexti, qscati, gi, wi, fzeroi
  public qextv, qscatv, gv, wv, fzerov
  public solar_1au, solar, gweight

  real(r8), parameter :: ubari = 0.5_r8
  real(r8), parameter :: tlimits = 1.0e-3_r8
  real(r8), parameter :: tlimiti = 5.0e-3_r8

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

  real(r8), allocatable, dimension(:) :: tauref
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

  real(r8), parameter :: wrefco2(l_refh2o) = [      &
    9.999999D-1, 9.99999D-1, 9.9999D-1, 9.999D-1,   &
    9.99D-1, 9.9D-1, 9.0D-1, 8.0D-1, 7.0D-1, 6.0D-1 &
  ]

  real(r8), parameter :: wrefh2o(l_refh2o) = [      &
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

    allocate(co2i(ntref,l_pint,l_refh2o,nspecti,ngauss))
    allocate(co2v(ntref,l_pint,l_refh2o,nspectv,ngauss))

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

    allocate(pfgasref(l_pint))

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

end module gomars_v1_rad_mod