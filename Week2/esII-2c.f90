  program rantest_coppie
!
! test program, call to intrinsic f90 random number generator
!  generate random numbers in [0,1[ ; then,
!  generate random integers between n_min and n_max.
!
  implicit none
  real, dimension(2) :: rnd 
  integer :: n, i,j

  print*,' n?'
  read(*,*)n
       do j = 1,n
          call random_number(rnd)
          write(2,*)rnd
       end do  


!   HINT

!    x(n)   di n. random

    do i=1,n,2
   write(1,*)x(2*i-1),x(2*i)
   end do

end program rantest_coppie

