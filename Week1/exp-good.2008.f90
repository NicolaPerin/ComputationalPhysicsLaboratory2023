!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     exp-good.f : a GOOD ALGORITHM to calculate e^-x
!                  as a FINTE sum of a series
!                  (to compare with exp-bad.f 
!                  and with the machine intrinsic function)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program expgood
!
! variable declaration:
!     accuracy limit,
!     min, max in x, step in x
!
      implicit none
      real ::  element, min, max, step, sum, x
      integer :: n
      min = 1.e-10
      max = 10.0
      step =0.1
      open(unit=7,file="exp-good.dat",status="replace",action="write")
               write(unit=7,fmt=*) "x, n, sum, exp(-x), abs(sum-exp(-x))/sum" 
!
! execute
!
        x = 0.0
      do 
        x = x + step
        if (x > max) exit
         sum     = 1
         element = 1
         do  n=1, 10000
            element = element*(-x)/n
            sum = sum + element
            if((abs(element/sum) < min) .and. (sum /= 0)) then
               write(unit=7,fmt=*) x, n, sum, exp(-x), abs(sum-exp(-x))/sum 
               go to 10
            endif
         enddo
 10   continue
       end do
       close(7)
!      stop "data saved in exp-good.dat"        
        end program expgood


















