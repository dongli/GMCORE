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

module gomars_v1_mp_mod

  use gomars_v1_const_mod
  use gomars_v1_tracers_mod

  implicit none

  private

  public mp_init
  public mp_final
  public rho_aer, std_aer
  public rad_rt, radb_rt
  public qextv_dst, qscatv_dst, gv_dst
  public qexti_dst, qscati_dst, gi_dst
  public qextv_cld, qscatv_cld, gv_cld
  public qexti_cld, qscati_cld, gi_cld
  public cor_ratio

  real(r8), parameter :: rmin  = 0.1e-6_r8
  real(r8), parameter :: rmax  = 10.0e-6_r8
  real(r8), parameter :: rbmin = 0.0001e-6_r8
  real(r8), parameter :: rbmax = 1.0e-2_r8
  real(r8), parameter :: cor_ratio(nratio) = [   &
    0.10_r8, 0.20_r8, 0.25_r8, 0.30_r8, 0.35_r8, &
    0.40_r8, 0.45_r8, 0.50_r8, 0.55_r8, 0.60_r8, &
    0.65_r8, 0.70_r8, 0.80_r8, 0.90_r8, 0.99_r8  &
  ]

  real(r8), allocatable, dimension(:    ) :: rho_aer
  real(r8), allocatable, dimension(:    ) :: std_aer
  real(r8), allocatable, dimension(:    ) :: rad_rt
  real(r8), allocatable, dimension(:    ) :: radb_rt
  real(r8), allocatable, dimension(:,:,:) :: qextv_cld
  real(r8), allocatable, dimension(:,:,:) :: qscatv_cld
  real(r8), allocatable, dimension(:,:,:) :: gv_cld
  real(r8), allocatable, dimension(:,:,:) :: qexti_cld
  real(r8), allocatable, dimension(:,:,:) :: qscati_cld
  real(r8), allocatable, dimension(:,:,:) :: gi_cld
  real(r8), allocatable, dimension(:,:  ) :: qextv_dst
  real(r8), allocatable, dimension(:,:  ) :: qscatv_dst
  real(r8), allocatable, dimension(:,:  ) :: gv_dst
  real(r8), allocatable, dimension(:,:  ) :: qexti_dst
  real(r8), allocatable, dimension(:,:  ) :: qscati_dst
  real(r8), allocatable, dimension(:,:  ) :: gi_dst

contains

  subroutine mp_init()

    integer i, j, l
    real(r8) vrat_rt
    real(r8) factor

    call mp_final()

    allocate(rho_aer(ntracers))
    allocate(std_aer(ntracers))

    rho_aer(iMa_dst) = rho_dst
    rho_aer(iNb_dst) = rho_dst
    rho_aer(iMa_cld) = rho_ice
    rho_aer(iNb_cld) = rho_ice
    rho_aer(iMa_cor) = rho_ice

    std_aer(iMa_dst) = dev_dst
    std_aer(iNb_dst) = dev_dst
    std_aer(iMa_cld) = dev_ice
    std_aer(iNb_cld) = dev_ice
    std_aer(iMa_cor) = dev_ice

    allocate(rad_rt (nbin_rt  ))
    allocate(radb_rt(nbin_rt+1))

    rad_rt(1) = 1.0e-7_r8
    rad_rt(nbin_rt) = 50.0e-6_r8
    radb_rt(1) = rbmin
    radb_rt(nbin_rt+1) = rbmax

    vrat_rt = log(rad_rt(nbin_rt) / rad_rt(1)) / (nbin_rt - 1.0_r8) * 3
    vrat_rt = exp(vrat_rt)

    do i = 1, nbin_rt - 1
      rad_rt(i+1) = rad_rt(i) * vrat_rt**athird
      radb_rt(i+1) = ((2 * vrat_rt) / (vrat_rt + 1))**athird * radb_rt(i)
    end do

    allocate(qextv_cld (nratio,nbin_rt,nspectv))
    allocate(qscatv_cld(nratio,nbin_rt,nspectv))
    allocate(gv_cld    (nratio,nbin_rt,nspectv))
    allocate(qexti_cld (nratio,nbin_rt,nspecti))
    allocate(qscati_cld(nratio,nbin_rt,nspecti))
    allocate(gi_cld    (nratio,nbin_rt,nspecti))
    allocate(qextv_dst (nbin_rt,nspectv))
    allocate(qscatv_dst(nbin_rt,nspectv))
    allocate(gv_dst    (nbin_rt,nspectv))
    allocate(qexti_dst (nbin_rt,nspecti))
    allocate(qscati_dst(nbin_rt,nspecti))
    allocate(gi_dst    (nbin_rt,nspecti))

    open(60, file='data/waterCoated_vis_JD_12bands.dat')
    open(61, file='data/waterCoated_ir_JD_12bands.dat')
    open(62, file='data/Dust_vis_wolff2010_JD_12bands.dat')
    open(63, file='data/Dust_ir_wolff2010_JD_12bands.dat')

    do j = 1, nratio
      do i = 1, nbin_rt
        read(60, '(7(e12.7,x))') (qextv_cld (j,i,l), l=1,nspectv)
        read(60, '(7(e12.7,x))') (qscatv_cld(j,i,l), l=1,nspectv)
        read(60, '(7(e12.7,x))') (gv_cld    (j,i,l), l=1,nspectv)
        read(61, '(5(e12.7,x))') (qexti_cld (j,i,l), l=1,nspecti)
        read(61, '(5(e12.7,x))') (qscati_cld(j,i,l), l=1,nspecti)
        read(61, '(5(e12.7,x))') (gi_cld    (j,i,l), l=1,nspecti)
      end do
    end do

    do i = 1, nbin_rt
      read(62, '(7(e11.5,x))') (qextv_dst (i,l), l=1,nspectv)
      read(62, '(7(e11.5,x))') (qscatv_dst(i,l), l=1,nspectv)
      read(62, '(7(e11.5,x))') (gv_dst    (i,l), l=1,nspectv)
      read(63, '(5(e11.5,x))') (qexti_dst (i,l), l=1,nspecti)
      read(63, '(5(e11.5,x))') (qscati_dst(i,l), l=1,nspecti)
      read(63, '(5(e11.5,x))') (gi_dst    (i,l), l=1,nspecti)
      ! factor = qextv_dst(i,6) / (vistoir*qexti_dst(i,4))
      factor = 1
      do l = 1, nspecti
        qexti_dst (i,l) = qexti_dst (i,l) * factor 
        qscati_dst(i,l) = qscati_dst(i,l) * factor 
      end do
    end do

    close(60)
    close(61)
    close(62)
    close(63)

  end subroutine mp_init

  subroutine mp_final()

    if (allocated(rho_aer   )) deallocate(rho_aer   )
    if (allocated(std_aer   )) deallocate(std_aer   )
    if (allocated(rad_rt    )) deallocate(rad_rt    )
    if (allocated(radb_rt   )) deallocate(radb_rt   )
    if (allocated(qextv_cld )) deallocate(qextv_cld )
    if (allocated(qscatv_cld)) deallocate(qscatv_cld)
    if (allocated(gv_cld    )) deallocate(gv_cld    )
    if (allocated(qexti_cld )) deallocate(qexti_cld )
    if (allocated(qscati_cld)) deallocate(qscati_cld)
    if (allocated(gi_cld    )) deallocate(gi_cld    )
    if (allocated(qextv_dst )) deallocate(qextv_dst )
    if (allocated(qscatv_dst)) deallocate(qscatv_dst)
    if (allocated(gv_dst    )) deallocate(gv_dst    )
    if (allocated(qexti_dst )) deallocate(qexti_dst )
    if (allocated(qscati_dst)) deallocate(qscati_dst)
    if (allocated(gi_dst    )) deallocate(gi_dst    )

  end subroutine mp_final

end module gomars_v1_mp_mod