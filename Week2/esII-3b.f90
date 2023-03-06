program rantest_intrinsic_with_seed
!
  ! test program, call to intrinsic f90 (ifort) random number generator
  !  generate random numbers in [0,1[ 
  !  dimension of seed depends on the machine architecture and compiler.
  !  
!
  implicit none
  real :: rnd
  real, dimension (3) :: sum=0.
  integer, dimension(:), allocatable :: seed, seed_old
  integer :: L,i,sizer,n_min,n_max,ran_int,k
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

! generates L random numbers in [0,1[                                              
     print*,' How many random numbers do you want to generate in [0,1[ ?'
     print*,' Insert the length of the sequence >'
     read(*,*)L        ! length of sequence                                          
    do i = 1,L
     call random_number(rnd)
     sum(1) = sum(1) + rnd      ! <x>, k=1
     sum(2) = sum(2) + rnd**3   ! <x>, k=3
     sum(3) = sum(3) + rnd**7   ! <x>, k=7
     write(1,*) i,  sum(1)/i - 1/2.,  sum(2)/i - 1/4.,  sum(3)/i - 1/8.
  end do
  
  deallocate(seed)
  deallocate(seed_old)

end program rantest_intrinsic_with_seed

