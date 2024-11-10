subroutine lagrange(x, xi, yi, ans)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  !  Lagrange interpolation - Polynomial interpolation at point x 
  !  xi(1) <= x <= xi(4).  Yi(n) is the functional value at XI(n).
  
  use gomars_v1_const_mod

  implicit none

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