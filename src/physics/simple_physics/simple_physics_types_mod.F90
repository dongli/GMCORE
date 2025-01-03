! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module simple_physics_types_mod

  use const_mod
  use tracer_mod
  use physics_types_mod

  implicit none

  type, extends(physics_state_type) :: simple_state_type
    real(r8), pointer, dimension(:,:) :: qv, qv_old
    real(r8), pointer, dimension(:,:) :: qc
    real(r8), pointer, dimension(:,:) :: qr
  contains
    procedure :: init  => simple_state_init
    procedure :: clear => simple_state_clear
    final simple_state_final
  end type simple_state_type

  type, extends(physics_tend_type) :: simple_tend_type
    real(r8), pointer, dimension(:,:) :: dqvdt
    real(r8), pointer, dimension(:,:) :: dqcdt
    real(r8), pointer, dimension(:,:) :: dqrdt
  contains
    procedure :: init  => simple_tend_init
    procedure :: clear => simple_tend_clear
    procedure :: reset => simple_tend_reset
    final simple_tend_final
  end type simple_tend_type

contains

  subroutine simple_state_init(this, mesh)

    class(simple_state_type), intent(inout), target :: this
    type(physics_mesh_type), intent(in) :: mesh

    call this%clear()

    call this%physics_state_init(mesh)

    if (idx_qv > 0) then
      this%qv_old => this%q_old(:,:,idx_qv)
      this%qv     => this%q    (:,:,idx_qv)
    end if
    if (idx_qc > 0) this%qc => this%q(:,:,idx_qc)
    if (idx_qr > 0) this%qr => this%q(:,:,idx_qr)

  end subroutine simple_state_init

  subroutine simple_state_clear(this)

    class(simple_state_type), intent(inout) :: this

    call this%physics_state_clear()

    nullify(this%qv_old)
    nullify(this%qv    )
    nullify(this%qc    )
    nullify(this%qr    )

  end subroutine simple_state_clear

  subroutine simple_state_final(this)

    type(simple_state_type), intent(inout) :: this

    call this%clear()

  end subroutine simple_state_final

  subroutine simple_tend_init(this, mesh)

    class(simple_tend_type), intent(inout), target :: this
    type(physics_mesh_type), intent(in) :: mesh

    call this%clear()

    call this%physics_tend_init(mesh)

    if (idx_qv > 0) this%dqvdt => this%dqdt(:,:,idx_qv)
    if (idx_qc > 0) this%dqcdt => this%dqdt(:,:,idx_qc)
    if (idx_qr > 0) this%dqrdt => this%dqdt(:,:,idx_qr)

  end subroutine simple_tend_init

  subroutine simple_tend_clear(this)

    class(simple_tend_type), intent(inout) :: this

    call this%physics_tend_clear()

    nullify(this%dqvdt)
    nullify(this%dqcdt)
    nullify(this%dqrdt)

  end subroutine simple_tend_clear

  subroutine simple_tend_reset(this)

    class(simple_tend_type), intent(inout) :: this

    call this%physics_tend_reset()

  end subroutine simple_tend_reset

  subroutine simple_tend_final(this)

    type(simple_tend_type), intent(inout) :: this

    call this%clear()

  end subroutine simple_tend_final

end module simple_physics_types_mod
