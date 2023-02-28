!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     exp-bad.f : a BAD ALGORITHM to calculate e^-x
!                 as a FINITE sum of the series
!                 (to compare with exp-good.f 
!                 and with the machine intrinsic function)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      program expbad
!
! variable declaration:
!     accuracy limit,
!     min, max in x, step in x
!     numerator (up), denominator (down)
!
      implicit none
      real :: down, min, max, step, sum, up, x
      integer :: i, n
      min = 1.e-10
      max = 10.0
      step =0.1
!
! execute
!     
      x=0.0
      do  
         if(x > max)exit
         x=x+step
         sum = 1
         do  n=1, 10000
            up   = 1
            down = 1
            do  i=1,n
               up   = -up*x
               down = down*i
            enddo
            sum = sum + up/down
            if((abs(up/down) < min) .and. (sum /= 0)) then
               print*,x,sum,exp(-x),up/down
               go to 10
            endif
        enddo
 10    continue
      enddo
      end program expbad


















