program rantest_intrinsi_with_seed
!
  ! test program, call to intrinsic f90 (ifort) random number generator
  !  generate random numbers in [0,1[ ; then,
  !  generate random integers between n_min and n_max.
  !  dimension of seed depends on the machine architecture and compiler.
  !  
!
  implicit none
  real :: rnd
  real, dimension (:), allocatable :: x
  integer, dimension(:), allocatable :: seed, seed_old
  integer :: L,i,sizer,n_min,n_max,ran_int
  call random_seed(sizer)
  allocate(seed(sizer))
  allocate(seed_old(sizer))
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'
  read(*,*)seed
!  illustrate use of put and get  
  call random_seed(put=seed)
  call random_seed(get=seed_old)
!  print the value of seed
  print *, "The seed you inserted is: ", seed_old

! generates ONE random number in [0,1[                                             
  call random_number(rnd)
     print *,' A real random number in [0,1[ is:',rnd

! generates L random numbers in [0,1[                                              
     print*,' How many random numbers do you want to generate in [0,1[ ?'
     print*,' Insert the length of the sequence >'
     read(*,*)L        ! length of sequence                                          
    do i = 1,L
     call random_number(rnd)
     print*, rnd
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
  print "(5f10.6)", x
  call random_seed(get=seed_old)
  print *, "now the value of seed is: ", seed_old
  
  deallocate(x)
  deallocate(seed)
  deallocate(seed_old)

end program rantest_intrinsic_with_seed

