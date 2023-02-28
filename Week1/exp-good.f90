!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     exp-good.f : a GOOD ALGORITHM to calculate e^-x
!                  as a FINITE sum of a series
!                  (to compare with exp-bad.f 
!                  and with the machine intrinsic function)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program expgood
  !
  ! variable declaration:
  !     x
  !     accuracy limit: min
  !
  implicit none
  real ::  element, sum, x, min = 1.e-10
  integer :: n
  open(unit=7,file="exp-good.dat",position="append",action="write")
  write(unit=7,fmt=*) "x, n, sum, exp(-x), abs(sum-exp(-x))/sum" 
  !
  ! execute
  !
  write(*,*)' enter x:'
  read(*,*) x
  sum     = 1
  element = 1
  do  n=1, 10000
     element = element*(-x)/n
     sum = sum + element
     if((abs(element/sum) < min) .and. (sum /= 0)) then
        write(*,*) x, n, sum, exp(-x), abs(sum-exp(-x))/sum 
        write(unit=7,fmt=*) x, n, sum, exp(-x), abs(sum-exp(-x))/sum 
        go to 10
     endif
  enddo
10 continue
  close(7)
  !      stop "data saved in exp-good.dat"        
end program expgood


















