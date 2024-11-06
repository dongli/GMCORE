module nasa_rad_mod

  use gomars_v1_const_mod

  implicit none

  integer, parameter :: l_nspecti = 5
  integer, parameter :: l_nspectv = 7
  integer, parameter :: l_ngauss  = 17
  integer, parameter :: l_npref   = 11
  integer, parameter :: l_ntref   = 7
  integer, parameter :: l_taumax  = 35
  integer, parameter :: l_pint    = 51
  integer, parameter :: l_refh2o  = 10
  integer, parameter :: l_nrefi   = 4
  integer, parameter :: l_nrefv   = 6

  real(r8), parameter :: ubari = 0.5_r8
  real(r8), parameter :: tlimits = 1.0e-3_r8
  real(r8), parameter :: tlimiti = 5.0e-3_r8

  real(r8), allocatable, dimension(:) :: wnoi
  real(r8), allocatable, dimension(:) :: dwni
  real(r8), allocatable, dimension(:) :: wavei
  real(r8), allocatable, dimension(:) :: wnov
  real(r8), allocatable, dimension(:) :: dwnv
  real(r8), allocatable, dimension(:) :: wavev
  real(r8), allocatable, dimension(:) :: solarf
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

  real(r8), allocatable, dimension(:) :: tauref(:)
  real(r8), allocatable, dimension(:) :: pfgasref(:)

  real(r8), parameter :: gweight(l_ngauss) = [ &
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

  subroutine nasa_rad_init()

    call nasa_rad_final()

    allocate(wnoi  (l_nspecti))
    allocate(dwni  (l_nspecti))
    allocate(wavei (l_nspecti))
    allocate(wnov  (l_nspectv))
    allocate(dwnv  (l_nspectv))
    allocate(wavev (l_nspectv))
    allocate(solarf(l_nspectv))
    allocate(tauray(l_nspectv))

    allocate(co2i(l_ntref,l_pint,l_refh2o,l_nspecti,l_ngauss))
    allocate(co2v(l_ntref,l_pint,l_refh2o,l_nspectv,l_ngauss))

    allocate(fzeroi(l_nspecti))
    allocate(fzerov(l_nspectv))

    allocate(pgasref(l_npref))
    allocate(tgasref(l_ntref))

    allocate(qextv (l_nspectv))
    allocate(qscatv(l_nspectv))
    allocate(wv    (l_nspectv))
    allocate(gv    (l_nspectv))
    allocate(qexti (l_nspecti))
    allocate(qscati(l_nspecti))
    allocate(wi    (l_nspecti))
    allocate(gi    (l_nspecti))

    allocate(planckir(l_nspecti,8501))

    allocate(pfgasref(l_pint))

    call setspv(wnov, dwnv, wavev, solarf, tauray)
    call setspi(wnoi, dwni, wavei)
    call setrad(tgasref, pfgasref, co2v, co2i, qextv, qscatv, wv, gv, &
                qexti, qscati, wi, gi, fzeroi, fzerov)

    ptop = 10**pfgasref(1)

  end subroutine nasa_rad_init

  subroutine nasa_rad_final()

    if (allocated(wnoi    )) deallocate(wnoi    )
    if (allocated(dwni    )) deallocate(dwni    )
    if (allocated(wavei   )) deallocate(wavei   )
    if (allocated(wnov    )) deallocate(wnov    )
    if (allocated(dwnv    )) deallocate(dwnv    )
    if (allocated(wavev   )) deallocate(wavev   )
    if (allocated(solarf  )) deallocate(solarf  )
    if (allocated(tauray  )) deallocate(tauray  )
    if (allocated(co2i    )) deallocate(co2i    )
    if (allocated(co2v    )) deallocate(co2v    )
    if (allocated(fzeroi  )) deallocate(fzeroi  )
    if (allocated(fzerov  )) deallocate(fzerov  )
    if (allocated(pgasref )) deallocate(pgasref )
    if (allocated(tgasref )) deallocate(tgasref )
    if (allocated(qextv   )) deallocate(qextv   )
    if (allocated(qscatv  )) deallocate(qscatv  )
    if (allocated(wv      )) deallocate(wv      )
    if (allocated(gv      )) deallocate(gv      )
    if (allocated(qexti   )) deallocate(qexti   )
    if (allocated(qscati  )) deallocate(qscati  )
    if (allocated(wi      )) deallocate(wi      )
    if (allocated(gi      )) deallocate(gi      )
    if (allocated(planckir)) deallocate(planckir)
    if (allocated(pfgasref)) deallocate(pfgasref)

  end subroutine nasa_rad_final

  subroutine setspv(wnov, dwnv, wavev, solarf, tauray)

    ! Legacy Mars GCM v24
    ! Mars Climate Modeling Center
    ! NASA Ames Research Center
  
    ! PURPOSE:
    !    Set up the spectral intervals in the visible (solar).  Based
    ! on Chris McKay's SETSPV code.
    !
    ! INPUT PARAMETERS
    ! L_NSPECTV  - Number of spectral intervals in the visible
    !
    ! OUTPUT PARAMETERS
    ! WNOV       - Array of wavenumbers at the spectral interval
    !              center for the visible.  Array is NSPECTV
    !              elements long.
    ! DWNV       - Array of "delta wavenumber", i.e., the width,
    !              in wavenumbers (cm^-1) of each visible spectral
    !              interval.  NSPECTV elements long.
    ! WAVEV      - Array (NSPECTV elements long) of the wavelenght
    !              (in microns) at the center of each visible spectral
    !              interval.
    ! SOLARF     - Array (NSPECTV elements) of solar flux (W/M^2) in
    !              each spectral interval.  Values are for 1 AU, and
    !              are scaled to the Mars distance elsewhere.
    ! TAURAY     - Array (NSPECTV elements) of the wavelength independent
    !              part of Rayleigh Scattering.  The pressure dependent 
    !              part is computed elsewhere (OPTCV).
   
    real(8) wnov  (l_nspectv)
    real(8) dwnv  (l_nspectv)
    real(8) wavev (l_nspectv)
    real(8) solarf(l_nspectv)
    real(8) tauray(l_nspectv)
  
    real(8) wl
    integer n, m
  
    real(8), parameter :: p0 = 9.423d+4 ! Rayleigh scattering reference pressure in hPa.
   
    ! BWNV - Bin wavenumber of the edges of the visible spectral bins
    ! units are inverse centimeters.  Dimension needs to be changed
    ! if the number of visible bins changes.

    ! Bin wavenumber - wavenumber [cm^(-1)] at the edges of the visible
    ! spectral bins.  Go from smaller to larger wavenumbers, the same as
    ! in the IR.
  
    ! 2222.22D0    ->   4.50 microns
    ! 3087.37D0    ->   3.24 microns
    ! 4030.63D0    ->   2.48 microns
    ! 5370.57D0    ->   1.86 microns
    ! 7651.11D0    ->   1.31 microns
    ! 12500.00D0   ->   0.80 microns
    ! 25000.00D0   ->   0.40 microns
    ! 41666.67D0   ->   0.24 microns
  
    real(8), parameter :: bwnv(l_nspectv+1) = [ &
      2222.22d0, 3087.37d0, 4030.63d0,  &
      5370.57d0, 7651.11d0, 12500.00d0, &
      25000.00d0, 41666.67d0            &
    ]
  
    ! Solar flux within each spectral interval, at 1AU (W/M^2)
    ! Sum equals 1356 W/m^2 (values from Wehrli, 1985)
  
    real(8), parameter :: solar(l_nspectv) = [ &
      12.7, 24.2, 54.6, 145.9, 354.9, 657.5, 106.3 ]
  
    ! Set up mean wavenumbers and wavenumber deltas.  Units of 
    ! wavenumbers is cm^(-1); units of wavelengths is microns.
  
    do m = 1, l_nspectv
      wnov(m)  = 0.5 * (bwnv(m+1) + bwnv(m))
      dwnv(m)  = bwnv(m+1) - bwnv(m)
      wavev(m) = 1.0e+4 / wnov(m)
    end do
  
    do n = 1, l_nspectv
      solarf(n) = solar(n)
    end do
  
    ! Set up the wavelength independent part of Rayleigh Scattering.
    ! The pressure dependent part will be computed elsewhere (OPTCV).
    ! WAVEV is in microns.  There is no Rayleigh scattering in the IR.
  
    do n = 1, l_nspectv
      wl = wavev(n)
      tauray(n) = (8.7_r8 / g) * (1.527_r8 * (1.0 + 0.013 / wl**2)/ wl**4) / p0
    end do
  
  end subroutine setspv

  subroutine setspi(wnoi, dwni, wavei)

    ! Legacy Mars GCM v24
    ! Mars Climate Modeling Center
    ! NASA Ames Research Center
  
    ! PURPOSE:
    !    Set up the spectral intervals in the infrared.  Based on
    ! Chris McKay's SETSPI code.
    !
    ! INPUT PARAMETERS
    ! L_NSPECTI  - Number of spectral intervals in the INFRARED
    !
    ! OUTPUT PARAMETERS
    ! WNOI       - Array of wavenumbers at the spectral interval
    !              centers for the infrared.  Array is NSPECTI
    !              elements long.
    ! DWNI       - Array of "delta wavenumber", i.e., the width,
    !              in wavenumbers (cm^-1) of each IR spectral
    !              interval.  NSPECTI elements long.
    ! WAVEI      - Array (NSPECTI elements long) of the wavelenght
    !              (in microns) at the center of each IR spectral
    !              interval.
  
    real(8) wnoi (l_nspecti)
    real(8) dwni (l_nspecti)
    real(8) wavei(l_nspecti)
  
    real(8) a, b, ans, y, bpa, bma, t
    real(8) wn1, wn2
    integer n, nw, nt, m
  
    ! C1 and C2 values from Goody and Yung (2nd edition)  MKS units
    ! These values lead to a "sigma" (sigma*T^4) of 5.67032E-8 W m^-2 K^-4
    real(8), parameter :: c1 = 3.741832d-16 ! W m^-2
    real(8), parameter :: c2 = 1.438786d-2  ! m K
    
    real(8) :: x(12) = [-0.981560634246719D0,  -0.904117256370475D0, &
                        -0.769902674194305D0,  -0.587317954286617D0, &
                        -0.367831498998180D0,  -0.125233408511469D0, &
                         0.125233408511469D0,   0.367831498998180D0, &
                         0.587317954286617D0,   0.769902674194305D0, &
                         0.904117256370475D0,   0.981560634246719D0]
  
    real(8) :: w(12) = [ 0.047175336386512D0,   0.106939325995318D0, &
                         0.160078328543346D0,   0.203167426723066D0, &
                         0.233492536538355D0,   0.249147045813403D0, &
                         0.249147045813403D0,   0.233492536538355D0, &
                         0.203167426723066D0,   0.160078328543346D0, &
                         0.106939325995318D0,   0.047175336386512D0]
  
    ! BWNI - Bin wavenumber of the edges of the IR spectral bins
    ! units are inverse centimeters.  Dimension needs to be changed
    ! if the number of IR bins changes.
  
    ! Bin wavenumber - wavenumber [cm^(-1)] at the edges of the IR
    ! spectral bins.
  
    real(8), parameter :: bwni(l_nspecti+1) = [ &
        10.000d0, & ! 1000.0 microns
       166.667d0, & !   60.0 microns
       416.667d0, & !   24.0 microns
       833.333d0, & !   12.0 microns
      1250.000d0, & !    8.0 microns
      2222.222d0  & !    4.5 microns
    ]
  
    ! Set up mean wavenumbers and wavenumber deltas.  Units of 
    ! wavenumbers is cm^(-1); units of wavelengths is microns.
    do m = 1, l_nspecti
      wnoi (m) = 0.5 * (bwni(m+1) + bwni(m))
      dwni (m) = bwni(m+1) - bwni(m)
      wavei(m) = 1.0e+4 / wnoi(m)
    end do
  
    !  For each IR wavelength interval, compute the integral of B(T), the
    !  Planck function, divided by the wavelength interval, in cm-1.  The
    !  integration is in MKS units, the final answer is the same as the
    !  original planck.f; W m^-2 wavenumber^-1, where wavenumber is in CM^-1.
    do nw = 1, l_nspecti
      a = 1.0d-2 / bwni(nw+1)
      b = 1.0d-2 / bwni(nw)
      bpa = (b + a) / 2.0
      bma = (b - a) / 2.0
      do nt = 500, 9000
        t   = dble(nt) / 1.0d+1
        ans = 0.0d0
        do m = 1, 12
          y   = bma * x(m) + bpa
          ans = ans + w(m) * c1 / (y**5 * (exp(c2 / (y * t)) - 1.0d0))
        end do
        planckir(nw,nt-499) = ans * bma / (pi * dwni(nw))
      end do
    end do
  
  end subroutine setspi  

  subroutine setrad(tgasref, pfgasref, co2v, co2i, qextv, qscatv, wv, gv, &
    qexti, qscati, wi, gi, fzeroi, fzerov)

    ! Legacy Mars GCM v24
    ! Mars Climate Modeling Center
    ! NASA Ames Research Center

    ! PURPOSE:
    !    Set up values used by the radiation code, such as the CO2 gas
    ! absorption coefficients.  True constants are defined, and the 
    ! time-independent quantities used by the radiation code are 
    ! calculated. 
    !
    ! INPUT PARAMETERS
    ! DTAU(L,M)      - Dust optical depth of layer L, and for aerosol 
    !                  species M.
    ! ptrop          - Pressure of the tropopause (mb)
    !
    ! OUTPUT PARAMETERS
    !
    ! AEROSOL RADIATIVE OPTICAL CONSTANTS
    ! Values are at the wavelenght interval center
    !
    ! MIE SCATTERING - Size distribution weighted
    ! Qextv    - Extinction efficiency - in the visible.
    ! Qscatv   - Scattering efficiency - in the visible.
    ! WV       - Single scattering albedo - in the visible.
    ! GV       - Asymmetry parameter - in the visible.
    !
    ! Qexti    - Extinction efficiency - in the infrared.
    ! Qscati   - Scattering efficiency - in the infrared.
    ! WI       - Single scattering albedo - in the infrared.
    ! GI       - Asymmetry parameter - in the infrared.

    real(8) co2i(l_ntref,l_pint,l_refh2o,l_nspecti,l_ngauss)
    real(8) co2v(l_ntref,l_pint,l_refh2o,l_nspectv,l_ngauss)
    real(8) pfgasref(l_pint)
    real(8) pgasref(l_npref)
    real(8) tgasref(l_ntref)
    real(8) qextv(l_nspectv)
    real(8) qscatv(l_nspectv)
    real(8) wv(l_nspectv)
    real(8) gv(l_nspectv)
    real(8) qexti(l_nspecti)
    real(8) qscati(l_nspecti)
    real(8) wi(l_nspecti)
    real(8) gi(l_nspecti)
    real(8) fzeroi(l_nspecti)
    real(8) fzerov(l_nspectv)

    integer n, ns, nt, np, nw, ng

    ! Visible dust properties:  M. Wolff Planck-weighted values (T=6000K)
    ! Log-normal size distribution:  Reff = 1.5 microns, Veff = 0.5

    ! Qext - M. Wolff values (order is increasing waveNUMBER)
    ! VISULAL WAVELENGTHS.
    real(8), parameter :: qev1(l_nspectv) = [ &
      1.834d0, 2.296d0, 2.672d0,        &
      2.829d0, 2.698d0, 2.452d0, 2.261d0]

    ! Qscat - M. Wolff values
    ! VISIBLE wavelengths
    real(8), parameter :: qsv1(l_nspectv) = [ &
      1.695d0, 2.031d0, 2.583d0,        &
      2.744d0, 2.626d0, 2.225d0, 1.525d0]

    ! G - M. Wolff values
    ! VISIBLE wavelengths

    real(8), parameter :: gv1(l_nspectv) = [ &
      0.551d0, 0.640d0, 0.661d0,        &
      0.678d0, 0.690d0, 0.743d0, 0.868d0]

    ! M. Wolff Planck-weighted values (T=215K)
    ! INFRARED wavelengths.  (The order is increasing waveNUMBER.)
    ! Qext for the IR
    real(8), parameter :: qei1(l_nspecti) = [ &
      0.008d0, 0.262d0, 0.491d0,        &
      1.017d0, 0.444d0                  ]

    ! Qsca for M. Wolff      INFRARED wavelengths
    real(8), parameter :: qsi1(l_nspecti) = [ &
      0.001d0, 0.037d0, 0.122d0,        &
      0.351d0, 0.336d0                  ]

    ! g for M. Wolff values  INFRARED wavelengths
    real(8), parameter :: gi1(L_NSPECTI)  = [ &
      0.004D0, 0.030D0, 0.095D0,        &
      0.214D0, 0.316D0                  ]

    ! Set the reference pressure and temperature arrays.  These are
    ! the pressures and temperatures at which we have k-coefficients.

    pgasref( 1) = 1.0e-6
    pgasref( 2) = 1.0e-5
    pgasref( 3) = 1.0e-4
    pgasref( 4) = 1.0e-3
    pgasref( 5) = 1.0e-2
    pgasref( 6) = 1.0e-1
    pgasref( 7) = 1.0
    pgasref( 8) = 1.0e+1
    pgasref( 9) = 1.0e+2
    pgasref(10) = 1.0e+3
    pgasref(11) = 1.0e+4

    tgasref(1)  =  50.0
    tgasref(2)  = 100.0
    tgasref(3)  = 150.0
    tgasref(4)  = 200.0
    tgasref(5)  = 250.0
    tgasref(6)  = 300.0
    tgasref(7)  = 350.0

    ! Fill the (VISIBLE) arrays Qextv, Qscatv, WV, GV
    do n = 1, l_nspectv
      qextv (n) = qev1(n)
      qscatv(n) = qsv1(n)
      if (qscatv(n) >= qextv(n)) then
        qscatv(n) = 0.99999 * qextv(n)
      end if
      wv(n) = qscatv(n) / qextv(n)
      gv(n) = gv1(n)
    end do

    ! Fill the (INFRARED) arrays Qexti, Qscati, WI, GI
    do n = 1, l_nspecti
      qexti (n) = qei1(n)
      qscati(n) = qsi1(n)
      if (qscati(n) >= qexti(n)) then
        qscati(n) = 0.99999 * qexti(n)
      end if
      wi(n) = qscati(n) / qexti(n)
      gi(n) = gi1(n)
    end do

    ! Interpolate CO2 k coefficients to the finer pressure grid.
    call laginterp(pgasref, pfgasref, co2i, co2v, fzeroi, fzerov)

  end subroutine setrad

  subroutine laginterp(pgref, pint, co2i, co2v, fzeroi, fzerov)

    !  Lagrange interpolation (linear in log pressure) of the CO2 
    !  k-coefficients in the pressure domain.  Subsequent use of these
    !  values will use a simple linear interpolation in pressure.
    
    real(8) co2i8(l_ntref,l_npref,l_refh2o,l_nspecti,l_ngauss)
    real(8) co2v8(l_ntref,l_npref,l_refh2o,l_nspectv,l_ngauss)
    real(8) pgref(l_npref)
    
    real(8) co2i(l_ntref,l_pint,l_refh2o,l_nspecti,l_ngauss)
    real(8) co2v(l_ntref,l_pint,l_refh2o,l_nspectv,l_ngauss)
    
    real(8) fzeroi(l_nspecti)
    real(8) fzerov(l_nspectv)
    
    real(8) x, xi(4), yi(4), ans
    real(8) pint(l_pint), pref(l_npref), p
    integer n, nt, np, nh, ng, nw, m, i
    
    real(8), parameter :: pin(l_pint) =  [    &
      -6.0d0, -5.8d0, -5.6d0, -5.4d0, -5.2d0, &
      -5.0d0, -4.8d0, -4.6d0, -4.4d0, -4.2d0, &
      -4.0d0, -3.8d0, -3.6d0, -3.4d0, -3.2d0, &
      -3.0d0, -2.8d0, -2.6d0, -2.4d0, -2.2d0, &
      -2.0d0, -1.8d0, -1.6d0, -1.4d0, -1.2d0, &
      -1.0d0, -0.8d0, -0.6d0, -0.4d0, -0.2d0, &
       0.0d0,  0.2d0,  0.4d0,  0.6d0,  0.8d0, &
       1.0d0,  1.2d0,  1.4d0,  1.6d0,  1.8d0, &
       2.0d0,  2.2d0,  2.4d0,  2.6d0,  2.8d0, &
       3.0d0,  3.2d0,  3.4d0,  3.6d0,  3.8d0, &
       4.0d0                                ]
    
    !  Fill pint for output from this subroutine
    do n = 1, l_pint
      pint(n) = pin(n)
    end do

    !  Take log of the reference pressures
    do n = 1, l_npref
      pref(n) = log10(pgref(n))
    end do
    
    ! Get CO2 k coefficients
    open(20, file='data/CO2H2O_V_12_95_INTEL', form='unformatted')
    read(20) co2v8
    read(20) fzerov
    close(20)
    
    open(20, file='data/CO2H2O_IR_12_95_INTEL', form='unformatted')
    read(20) co2i8
    read(20) fzeroi
    close(20)
    
    ! Take Log10 of the values - we interpolate the log10 of the values,
    ! not the values themselves.   Smallest value is 1.0E-200.
    do nt = 1, l_ntref
      do np = 1, l_npref
        do nh = 1, l_refh2o
          do ng = 1, l_ngauss
            do nw = 1, l_nspectv
              if (co2v8(nt,np,nh,nw,ng) > 1.0d-200) then
                co2v8(nt,np,nh,nw,ng) = log10(co2v8(nt,np,nh,nw,ng))
              else
                co2v8(nt,np,nh,nw,ng) = -200.0
              end if
            end do
            do nw = 1, l_nspecti
              if (co2i8(nt,np,nh,nw,ng) > 1.0d-200) then
                co2i8(nt,np,nh,nw,ng) = log10(co2i8(nt,np,nh,nw,ng))
              else
                co2i8(nt,np,nh,nw,ng) = -200.0
              end if
            end do
          end do
        end do
      end do
    end do
    !  Interpolate the values:  first the IR
    do nt = 1, l_ntref
      do nh = 1, l_refh2o
        do nw = 1, l_nspecti
          do ng = 1, l_ngauss
            ! First, the initial interval (P=1e-6 to 1e-5)
            n = 1 
            do m = 1, 5
              x     = pint(m)
              xi(1) = pref(n)
              xi(2) = pref(n+1)
              xi(3) = pref(n+2)
              xi(4) = pref(n+3)
              yi(1) = co2i8(nt,n,nh,nw,ng)
              yi(2) = co2i8(nt,n+1,nh,nw,ng)
              yi(3) = co2i8(nt,n+2,nh,nw,ng)
              yi(4) = co2i8(nt,n+3,nh,nw,ng)
              call lagrange(x, xi, yi, ans)
              co2i(nt,m,nh,nw,ng) = 10.0**ans
            end do 
            do n = 2, l_npref - 2
              do m = 1, 5
                i     = (n - 1) * 5 + m
                x     = pint(i)
                xi(1) = pref(n-1)
                xi(2) = pref(n)
                xi(3) = pref(n+1)
                xi(4) = pref(n+2)
                yi(1) = co2i8(nt,n-1,nh,nw,ng)
                yi(2) = co2i8(nt,n,nh,nw,ng)
                yi(3) = co2i8(nt,n+1,nh,nw,ng)
                yi(4) = co2i8(nt,n+2,nh,nw,ng)
                call lagrange(x, xi, yi, ans)
                co2i(nt,i,nh,nw,ng) = 10.0**ans
              end do 
            end do
            !  Now, get the last interval (P=1e+3 to 1e+4)
            n = L_NPREF-1
            do m = 1, 5
              i     = (n - 1) * 5 + m
              x     = pint(i)
              xi(1) = pref(n-2)
              xi(2) = pref(n-1)
              xi(3) = pref(n)
              xi(4) = pref(n+1)
              yi(1) = co2i8(nt,n-2,nh,nw,ng)
              yi(2) = co2i8(nt,n-1,nh,nw,ng)
              yi(3) = co2i8(nt,n,nh,nw,ng)
              yi(4) = co2i8(nt,n+1,nh,nw,ng)
              call lagrange(x, xi, yi, ans)
              co2i(nt,i,nh,nw,ng) = 10.0**ans
            end do  
            !  Fill the last pressure point
            co2i(nt,l_pint,nh,nw,ng) = 10.0**co2i8(nt,l_npref,nh,nw,ng)
          end do
        end do
      end do
    end do
    ! Interpolate the values:  now the visible
    do nt=1,L_NTREF
      do nh=1,L_REFH2O
        do nw=1,L_NSPECTV
          do ng=1,L_NGAUSS
            ! First, the initial interval (P=1e-6 to 1e-5)
            n = 1 
            do m = 1, 5
              x     = pint(m)
              xi(1) = pref(n)
              xi(2) = pref(n+1)
              xi(3) = pref(n+2)
              xi(4) = pref(n+3)
              yi(1) = co2v8(nt,n,nh,nw,ng)
              yi(2) = co2v8(nt,n+1,nh,nw,ng)
              yi(3) = co2v8(nt,n+2,nh,nw,ng)
              yi(4) = co2v8(nt,n+3,nh,nw,ng)
              call lagrange(x, xi, yi, ans)
              co2v(nt,m,nh,nw,ng) = 10.0**ans
            end do 
            do n = 2, l_npref - 2
              do m = 1, 5
                i     = (n - 1) * 5 + m
                x     = pint(i)
                xi(1) = pref(n-1)
                xi(2) = pref(n)
                xi(3) = pref(n+1)
                xi(4) = pref(n+2)
                yi(1) = co2v8(nt,n-1,nh,nw,ng)
                yi(2) = co2v8(nt,n,nh,nw,ng)
                yi(3) = co2v8(nt,n+1,nh,nw,ng)
                yi(4) = co2v8(nt,n+2,nh,nw,ng)
                call lagrange(x, xi, yi, ans)
                co2v(nt,i,nh,nw,ng) = 10.0**ans
              end do 
            end do
            !  Now, get the last interval (P=1e+3 to 1e+4)
            n = l_npref - 1
            do m = 1, 5
              i     = (n - 1) * 5 + m
              x     = pint(i)
              xi(1) = pref(n-2)
              xi(2) = pref(n-1)
              xi(3) = pref(n)
              xi(4) = pref(n+1)
              yi(1) = co2v8(nt,n-2,nh,nw,ng)
              yi(2) = co2v8(nt,n-1,nh,nw,ng)
              yi(3) = co2v8(nt,n,nh,nw,ng)
              yi(4) = co2v8(nt,n+1,nh,nw,ng)
              call lagrange(x, xi, yi, ans)
              co2v(nt,i,nh,nw,ng) = 10.0**ans
            end do  
            !  Fill the last pressure point
            co2v(nt,l_pint,nh,nw,ng) = 10.0**co2v8(nt,l_npref,nh,nw,ng)
          end do
        end do
      end do
    end do
    
  end subroutine laginterp

  subroutine lagrange(x, xi, yi, ans)

    !  Lagrange interpolation - Polynomial interpolation at point x 
    !  xi(1) <= x <= xi(4).  Yi(n) is the functional value at XI(n).
    
    real(8) x, xi(4), yi(4), ans
    real(8) fm1, fm2, fm3, fm4
    
    fm1 = x - xi(1)
    fm2 = x - xi(2)
    fm3 = x - xi(3)
    fm4 = x - xi(4)
    
    !  Get the "answer" at the requested X
    ans = fm2 * fm3 * fm4 * yi(1) / ((xi(1) - xi(2)) * (xi(1) - xi(3)) * (xi(1) - xi(4))) + &
          fm1 * fm3 * fm4 * yi(2) / ((xi(2) - xi(1)) * (xi(2) - xi(3)) * (xi(2) - xi(4))) + &
          fm1 * fm2 * fm4 * yi(3) / ((xi(3) - xi(1)) * (xi(3) - xi(2)) * (xi(3) - xi(4))) + &
          fm1 * fm2 * fm3 * yi(4) / ((xi(4) - xi(1)) * (xi(4) - xi(2)) * (xi(4) - xi(3))) 
    
  end subroutine lagrange

end module nasa_rad_mod