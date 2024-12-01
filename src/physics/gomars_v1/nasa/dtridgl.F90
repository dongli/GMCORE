subroutine dtridgl(n, a, b, c, d, x)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod, only: r8

  implicit none

  integer, parameter :: nmax = 201

  integer , intent(in ) :: n
  real(r8), intent(in ) :: a(n)
  real(r8), intent(in ) :: b(n)
  real(r8), intent(in ) :: c(n)
  real(r8), intent(in ) :: d(n)
  real(r8), intent(out) :: x(n)

  integer i
  real(r8) as(nmax), ds(nmax), cs

  as(n) = a(n) / b(n)
  ds(n) = d(n) / b(n)

  do i = 2, n
    cs = 1.0_r8 / (b(n+1-i) - c(n+1-i) * as(n+2-i))
    as(n+1-i) = a(n+1-i) * cs
    ds(n+1-i) = (d(n+1-i) + c(n+1-i) * ds(n+2-i)) * cs
  end do

  x(1) = ds(1)
  do i = 2, n
    x(i) = ds(i) - as(i) * x(i-1)
  end do

end subroutine dtridgl