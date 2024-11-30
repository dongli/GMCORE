subroutine laginterp(pgref, pint, co2i, co2v, fzeroi, fzerov)

  !  Lagrange interpolation (linear in log pressure) of the CO2 
  !  k-coefficients in the pressure domain.  Subsequent use of these
  !  values will use a simple linear interpolation in pressure.
  
  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in ) :: pgref (npref)
  real(r8), intent(out) :: pint  (npint)
  real(r8), intent(out) :: co2i  (ntref,npint,nrefh2o,nspecti,ngauss)
  real(r8), intent(out) :: co2v  (ntref,npint,nrefh2o,nspectv,ngauss)
  real(r8), intent(out) :: fzeroi(nspecti)
  real(r8), intent(out) :: fzerov(nspectv)
  
  real(r8) co2i8(ntref,npref,nrefh2o,nspecti,ngauss)
  real(r8) co2v8(ntref,npref,nrefh2o,nspectv,ngauss)
  real(r8) x, xi(4), yi(4), ans
  real(r8) pref(npref), p
  integer n, nt, np, nh, ng, nw, m, i
  
  real(8) pin(npint)

  pint = [                                  &
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
     4.0d0                                  &
  ]


  !  Fill pint for output from this subroutine
  do n = 1, npint
    pint(n) = pin(n)
  end do

  !  Take log of the reference pressures
  do n = 1, npref
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
  do nt = 1, ntref
    do np = 1, npref
      do nh = 1, nrefh2o
        do ng = 1, ngauss
          do nw = 1, nspectv
            if (co2v8(nt,np,nh,nw,ng) > 1.0d-200) then
              co2v8(nt,np,nh,nw,ng) = log10(co2v8(nt,np,nh,nw,ng))
            else
              co2v8(nt,np,nh,nw,ng) = -200.0
            end if
          end do
          do nw = 1, nspecti
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
  do nt = 1, ntref
    do nh = 1, nrefh2o
      do nw = 1, nspecti
        do ng = 1, ngauss
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
          do n = 2, npref - 2
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
          n = npref - 1
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
          co2i(nt,npint,nh,nw,ng) = 10.0**co2i8(nt,npref,nh,nw,ng)
        end do
      end do
    end do
  end do
  ! Interpolate the values:  now the visible
  do nt = 1, ntref
    do nh = 1, nrefh2o
      do nw = 1, nspectv
        do ng = 1, ngauss
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
          do n = 2, npref - 2
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
          n = npref - 1
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
          co2v(nt,npint,nh,nw,ng) = 10.0**co2v8(nt,npref,nh,nw,ng)
        end do
      end do
    end do
  end do
  
end subroutine laginterp