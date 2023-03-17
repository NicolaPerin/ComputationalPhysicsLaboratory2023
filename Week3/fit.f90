!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  fit.f: Least square fit to decay spectrum                           c
!								       c
!  from: "Projects in Computational Physics" by Landau and Paez  c 
!	       copyrighted by John Wiley and Sons, New York            c      
!    			                                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Program fit
      Implicit none
!
! declarations
      Integer :: i
      Integer, parameter :: dp=selected_real_kind(13)
      Real(kind=dp) :: s, sx, sy, sxx, sxy, delta, inter, slope
      Real(kind=dp),dimension(12) :: x, y, d
!
! input value y - exponential fit y > 0      
      y(1:12) =  (/32.0_dp, 17.0_dp, 21.0_dp, 7.0_dp, 8.0_dp, 6.0_dp, 5.0_dp,&
                   2.0_dp, 2.0_dp, 0.1_dp, 4.0_dp, 1.0_dp/)
!
! input values x
      Do i=1, 12
         x(i)=i*10-5
      end do
!
! input value delta y - estimate
      Do i=1, 12
         d(i)=1.0
      end do
!
! take logs of y values for exponential fit
      Do i=1, 12
         y(i)=log(y(i))
      end do
!
! calculate all the sums
      Do i=1, 12
         s   = s   +         1 / (d(i)*d(i))
         sx  = sx  +      x(i) / (d(i)*d(i))
         sy  = sy  +      y(i) / (d(i)*d(i))
         sxx = sxx + x(i)*x(i) / (d(i)*d(i)) 
         sxy = sxy + x(i)*y(i) / (d(i)*d(i))
      end do
!
! calculate the coefficients
      delta= s*sxx-sx*sx
      slope=  (s*sxy-sx*sy) / delta		
      inter=(sxx*sy-sx*sxy) / delta
      print*, "intercept=", inter
      print*, "slope=", slope
      print*, "correlation=", -sx/sqrt(sxx*s)		
!
      End program fit
