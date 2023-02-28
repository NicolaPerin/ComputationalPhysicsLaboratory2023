      program deriv
!c
!c  numerical derivative: left, right, symmetric in SINGLE PRECISION
!c
      real :: h(8) 
      real :: x, exact
      integer :: i, N=8 
      data h/0.5, 0.2, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001 /
!c
      print*, " h, f'_ds, error, f'_sin, error, f'_simm, error "
      x = 1.0                   ! inizialize variables
      exact = cos(x)
       do i=1,N
      deriv_ds = (sin(x+h(i))-sin(x)) / h(i)
      deriv_sin = (sin(x)-sin(x-h(i))) / h(i)
      deriv_simm = (sin(x+h(i))-sin(x-h(i))) / (2*h(i))
         print*,  h(i), deriv_ds, deriv_ds - exact, deriv_sin, deriv_sin - exact, &
  & deriv_simm, deriv_simm - exact
      end do                 
      stop
      end program deriv
      
