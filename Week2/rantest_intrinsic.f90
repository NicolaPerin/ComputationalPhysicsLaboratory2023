program rantest_intrinsic
!
! test program, call to intrinsic f90 random number generator
!  generate random numbers in [0,1[ ; then,
!  generate random integers between n_min and n_max.
!
  implicit none
  real :: rnd
  real, dimension (:), allocatable :: x
  integer :: L,i,n_min,n_max,ran_int

  ! generates ONE random number in [0,1[
     call random_number(rnd)
     print *,' A real random number in [0,1[ is:',rnd

  ! generates L random numbers in [0,1[
     print*,' How many random numbers do you want to generate in [0,1[ ?'
     print*,' Insert the length of the sequence >'
     read(*,*)L        ! length of sequence
    do i = 1,L
     call random_number(rnd)
     print *,rnd
  end do

  !  generates integer random numbers between n_min and n_max
     print*,' Generate ',L,' integer random numbers in [n_min,n_max[ ;'
     print*,' insert n_min, n_max >'
     read(*,*),n_min,n_max
    do i = 1,L
     call random_number(rnd)
     ran_int = (n_max - n_min + 1)*rnd + n_min
     print *,ran_int
  end do

  !  use array x to generate and store L random numbers with a unique call
  print*,' Generate other ',L,' real random numbers in [0,1[:'
  allocate(x(L))
  call random_number(x)
  print*, x
  deallocate(x)

end program rantest_intrinsic

