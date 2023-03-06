  program rantest_isto
!
! test program, call to intrinsic f90 random number generator
!  generate random numbers in [0,1[ ; then,
!  generate random integers between n_min and n_max.
!
  implicit none
  real :: rnd , delta_r=0.1
  integer, dimension(10) :: histog=0
  integer :: n, i,j

  print*,' n?'
  read(*,*)n
       do j = 1,n
          call random_number(rnd)
          i = int(rnd/delta_r) + 1
          histog(i) = histog(i) + 1
       end do  

  do i=1,10
     write(1,*)(i-0.5)*delta_r, histog(i)
  end do
end program rantest_isto

