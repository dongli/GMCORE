
module radheat

  !-----------------------------------------------------------------------
  !
  ! Purpose:  Provide an interface to convert shortwave and longwave
  !           radiative heating terms into net heating.
  !
  !           This module provides a hook to allow incorporating additional
  !           radiative terms (eUV heating and nonLTE longwave cooling).
  !
  ! Original version: B.A. Boville
  !-----------------------------------------------------------------------

  use shr_kind_mod  , only: r8 => shr_kind_r8
  use ppgrid        , only: pcols, pver
  use physics_types , only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc

  implicit none
  private
  save

  public radheat_readnl
  public radheat_init
  public radheat_timestep_init
  public radheat_tend
  public radheat_disable_waccm

contains

  subroutine radheat_readnl(nlfile)

    character(*), intent(in) :: nlfile

    ! No options for this version of radheat; this is just a stub.

  end subroutine radheat_readnl

  subroutine radheat_init(pref_mid)

    use pmgrid        , only: plev
    use physics_buffer, only: physics_buffer_desc

    real(r8), intent(in) :: pref_mid(plev)

  end subroutine radheat_init

  subroutine radheat_timestep_init (state, pbuf2d)

    use physics_types , only: physics_state
    use ppgrid        , only: begchunk, endchunk
    use physics_buffer, only: physics_buffer_desc

    type(physics_state), intent(in) :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

  end subroutine radheat_timestep_init

  subroutine radheat_tend(state, pbuf,  ptend, qrl, qrs, fsns, &
                          fsnt, flns, flnt, asdir, net_flx)

    !-----------------------------------------------------------------------
    ! Compute net radiative heating from qrs and qrl, and the associated net
    ! boundary flux.
    !-----------------------------------------------------------------------

    type(physics_state), intent(in) :: state              ! Physics state variables

    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_ptend), intent(out) :: ptend             ! Indivdual parameterization tendencie
    real(r8),            intent(in)  :: qrl(pcols,pver)   ! Longwave heating
    real(r8),            intent(in)  :: qrs(pcols,pver)   ! Shortwave heating
    real(r8),            intent(in)  :: fsns(pcols)       ! Surface solar absorbed flux
    real(r8),            intent(in)  :: fsnt(pcols)       ! Net column abs solar flux at model top
    real(r8),            intent(in)  :: flns(pcols)       ! Srf longwave cooling (up-down) flux
    real(r8),            intent(in)  :: flnt(pcols)       ! Net outgoing lw flux at model top
    real(r8),            intent(in)  :: asdir(pcols)      ! Shortwave, direct albedo
    real(r8),            intent(out) :: net_flx(pcols)

    integer i, k, ncol

    ncol = state%ncol

    call physics_ptend_init(ptend,state%psetcols, 'radheat', ls=.true.)

    ptend%s(:ncol,:) = (qrs(:ncol,:) + qrl(:ncol,:))

    do i = 1, ncol
      net_flx(i) = fsnt(i) - fsns(i) - flnt(i) + flns(i)
    end do

  end subroutine radheat_tend

  subroutine radheat_disable_waccm()
  end subroutine radheat_disable_waccm

end module radheat
