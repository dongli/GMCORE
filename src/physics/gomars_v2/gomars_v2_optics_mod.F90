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

module gomars_v2_optics_mod

  use fiona
  use const_mod
  use gomars_v2_namelist_mod
  use gomars_v2_spectra_mod

  implicit none

  private

  public gomars_v2_optics_init
  public gomars_v2_optics_final

  ! Symbol notation from Petty (2006):
  ! Absorption coefficient (m-1): βa
  ! Scattering coefficient (m-1): βs
  ! Extinction coefficient (m-1): βe = βa + βs
  ! Mass extinction coefficient (m2 kg-1): ke = βe / ρ
  ! Volume extinction coefficient or extinction cross-section (m2): 𝜎e = ke m = βe / N
  ! Scattering cross-section (m2): 𝜎s
  ! Extinction efficiency (1): Qe = 𝜎e / A = 𝜎e / (𝜋 r2)
  ! Single scattering alqedo (1): ω = βs / βe
  ! Optical depth (1): τ(s1,s2) = ∫ βe(s) ds from s1 to s2
  ! Transmittance (1): t(s1,s2) = exp(-τ(s1,s2))
  ! Scattering asymmetry parameter: g = ∫ p(θ) cos(θ) dΩ / ∫ p(θ) dΩ

  integer, public :: nbin_dst = 0
  integer, public :: nratio_cld = 0
  real(r8), public, allocatable, dimension(    :) :: dst_r_bnds
  real(r8), public, allocatable, dimension(    :) :: dst_qe0_vs, dst_qe0_ir
  real(r8), public, allocatable, dimension(    :) :: dst_qs0_vs, dst_qs0_ir
  real(r8), public, allocatable, dimension(    :) :: dst_gs0_vs, dst_gs0_ir
  real(r8), public, allocatable, dimension(    :) :: dst_ws0_vs, dst_ws0_ir
  real(r8), public, allocatable, dimension(  :,:) :: dst_qe_vs , dst_qe_ir
  real(r8), public, allocatable, dimension(  :,:) :: dst_qs_vs , dst_qs_ir
  real(r8), public, allocatable, dimension(  :,:) :: dst_gs_vs , dst_gs_ir
  real(r8), public, allocatable, dimension(  :,:) :: dst_ws_vs , dst_ws_ir
  real(r8), public, allocatable, dimension(:    ) :: cld_ratio
  real(r8), public, allocatable, dimension(:,:,:) :: cld_qe_vs , cld_qe_ir
  real(r8), public, allocatable, dimension(:,:,:) :: cld_qs_vs , cld_qs_ir
  real(r8), public, allocatable, dimension(:,:,:) :: cld_gs_vs , cld_gs_ir
  real(r8), public, allocatable, dimension(:,:,:) :: cld_ws_vs , cld_ws_ir

contains

  subroutine gomars_v2_optics_init(nlev_rad)

    integer, intent(in) :: nlev_rad

    integer i, n

    call gomars_v2_optics_final()

    call fiona_open_dataset('dust_optics', dust_optics_file)
    call fiona_get_dim('dust_optics', 'bin', size=nbin_dst)
    call fiona_get_dim('dust_optics', 'spec_vis', size=n)
    if (n /= spec_vs%n) then
      stop 999
    end if
    call fiona_get_dim('dust_optics', 'spec_ir', size=n)
    if (n /= spec_ir%n) then
      stop 999
    end if
    allocate(dst_r_bnds (nbin_dst+1))
    allocate(dst_qe0_vs (         spec_vs%n))
    allocate(dst_qs0_vs (         spec_vs%n))
    allocate(dst_gs0_vs (         spec_vs%n))
    allocate(dst_ws0_vs (         spec_vs%n))
    allocate(dst_qe_vs  (nbin_dst,spec_vs%n))
    allocate(dst_qs_vs  (nbin_dst,spec_vs%n))
    allocate(dst_gs_vs  (nbin_dst,spec_vs%n))
    allocate(dst_ws_vs  (nbin_dst,spec_vs%n))
    allocate(dst_qe0_ir (         spec_ir%n))
    allocate(dst_qs0_ir (         spec_ir%n))
    allocate(dst_gs0_ir (         spec_ir%n))
    allocate(dst_ws0_ir (         spec_ir%n))
    allocate(dst_qe_ir  (nbin_dst,spec_ir%n))
    allocate(dst_qs_ir  (nbin_dst,spec_ir%n))
    allocate(dst_gs_ir  (nbin_dst,spec_ir%n))
    allocate(dst_ws_ir  (nbin_dst,spec_ir%n))
    call fiona_start_input('dust_optics')
    call fiona_input('dust_optics', 'dust_r_bnds' , dst_r_bnds)
    call fiona_input('dust_optics', 'dust_qe0_vis', dst_qe0_vs)
    call fiona_input('dust_optics', 'dust_qs0_vis', dst_qs0_vs)
    call fiona_input('dust_optics', 'dust_gs0_vis', dst_gs0_vs)
    call fiona_input('dust_optics', 'dust_qe_vis' , dst_qe_vs )
    call fiona_input('dust_optics', 'dust_qs_vis' , dst_qs_vs )
    call fiona_input('dust_optics', 'dust_gs_vis' , dst_gs_vs )
    call fiona_input('dust_optics', 'dust_qe0_ir' , dst_qe0_ir)
    call fiona_input('dust_optics', 'dust_qs0_ir' , dst_qs0_ir)
    call fiona_input('dust_optics', 'dust_gs0_ir' , dst_gs0_ir)
    call fiona_input('dust_optics', 'dust_qe_ir'  , dst_qe_ir )
    call fiona_input('dust_optics', 'dust_qs_ir'  , dst_qs_ir )
    call fiona_input('dust_optics', 'dust_gs_ir'  , dst_gs_ir )
    call fiona_end_input('dust_optics')

    do i = 1, spec_vs%n
      if (dst_qs0_vs(i) >= dst_qe0_vs(i)) then
        dst_qs0_vs(i) = 0.99999_r8 * dst_qe0_vs(i)
      end if
      dst_ws0_vs(i) = dst_qs0_vs(i) / dst_qe0_vs(i)
    end do
    do i = 1, spec_ir%n
      if (dst_qs0_ir(i) >= dst_qe0_ir(i)) then
        dst_qs0_ir(i) = 0.99999_r8 * dst_qe0_ir(i)
      end if
      dst_ws0_ir(i) = dst_qs0_ir(i) / dst_qe0_ir(i)
    end do

    call fiona_open_dataset('cld_optics', cld_optics_file)
    call fiona_get_dim('cld_optics', 'ratio', size=nratio_cld)
    call fiona_get_dim('cld_optics', 'spec_vis', size=n)
    if (n /= spec_vs%n) then
      stop 999
    end if
    call fiona_get_dim('cld_optics', 'spec_ir', size=n)
    if (n /= spec_ir%n) then
      stop 999
    end if
    allocate(cld_ratio (nratio_cld))
    allocate(cld_qe_vs(nratio_cld,nbin_dst,spec_vs%n))
    allocate(cld_qs_vs(nratio_cld,nbin_dst,spec_vs%n))
    allocate(cld_gs_vs(nratio_cld,nbin_dst,spec_vs%n))
    allocate(cld_ws_vs(nratio_cld,nbin_dst,spec_vs%n))
    allocate(cld_qe_ir(nratio_cld,nbin_dst,spec_ir%n))
    allocate(cld_qs_ir(nratio_cld,nbin_dst,spec_ir%n))
    allocate(cld_gs_ir(nratio_cld,nbin_dst,spec_ir%n))
    allocate(cld_ws_ir(nratio_cld,nbin_dst,spec_ir%n))
    call fiona_start_input('cld_optics')
    call fiona_input('cld_optics', 'cld_ratio' , cld_ratio )
    call fiona_input('cld_optics', 'cld_qe_vis', cld_qe_vs)
    call fiona_input('cld_optics', 'cld_qs_vis', cld_qs_vs)
    call fiona_input('cld_optics', 'cld_gs_vis', cld_gs_vs)
    call fiona_input('cld_optics', 'cld_qe_ir' , cld_qe_ir)
    call fiona_input('cld_optics', 'cld_qs_ir' , cld_qs_ir)
    call fiona_input('cld_optics', 'cld_gs_ir' , cld_gs_ir)
    call fiona_end_input('cld_optics')

  end subroutine gomars_v2_optics_init

  subroutine gomars_v2_optics_final()

    if (allocated(dst_r_bnds)) deallocate(dst_r_bnds)
    if (allocated(dst_qe0_vs)) deallocate(dst_qe0_vs)
    if (allocated(dst_qs0_vs)) deallocate(dst_qs0_vs)
    if (allocated(dst_gs0_vs)) deallocate(dst_gs0_vs)
    if (allocated(dst_ws0_vs)) deallocate(dst_ws0_vs)
    if (allocated(dst_qe0_ir)) deallocate(dst_qe0_ir)
    if (allocated(dst_qs0_ir)) deallocate(dst_qs0_ir)
    if (allocated(dst_gs0_ir)) deallocate(dst_gs0_ir)
    if (allocated(dst_ws0_ir)) deallocate(dst_ws0_ir)
    if (allocated(dst_qe_vs )) deallocate(dst_qe_vs )
    if (allocated(dst_qs_vs )) deallocate(dst_qs_vs )
    if (allocated(dst_gs_vs )) deallocate(dst_gs_vs )
    if (allocated(dst_ws_vs )) deallocate(dst_ws_vs )
    if (allocated(dst_qe_ir )) deallocate(dst_qe_ir )
    if (allocated(dst_qs_ir )) deallocate(dst_qs_ir )
    if (allocated(dst_gs_ir )) deallocate(dst_gs_ir )
    if (allocated(dst_ws_ir )) deallocate(dst_ws_ir )
    if (allocated(cld_ratio )) deallocate(cld_ratio )
    if (allocated(cld_qe_vs )) deallocate(cld_qe_vs )
    if (allocated(cld_qs_vs )) deallocate(cld_qs_vs )
    if (allocated(cld_gs_vs )) deallocate(cld_gs_vs )
    if (allocated(cld_ws_vs )) deallocate(cld_ws_vs )
    if (allocated(cld_qe_ir )) deallocate(cld_qe_ir )
    if (allocated(cld_qs_ir )) deallocate(cld_qs_ir )
    if (allocated(cld_gs_ir )) deallocate(cld_gs_ir )
    if (allocated(cld_ws_ir )) deallocate(cld_ws_ir )

  end subroutine gomars_v2_optics_final

end module gomars_v2_optics_mod
