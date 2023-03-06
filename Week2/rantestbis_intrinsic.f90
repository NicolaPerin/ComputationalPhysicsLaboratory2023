program rantestbis_intrinsic
  ! test program, call to intrinsic random number f90 generator
  ! iluustrate the use of "put" and "get" in  call random_seed
  implicit none
  integer, dimension(:), allocatable :: seed, seed_old
  integer :: L,i,sizer,n_min,n_max,ran_int
  real, dimension(3) :: harvest

  call random_seed(sizer)
  allocate(seed(sizer))
  allocate(seed_old(sizer))
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'
  read(*,*)seed

  call random_seed(put=seed)
  call random_seed(get=seed_old)
  print*, "Old starting value: ",seed_old

  call random_number(harvest)
  print*,"3 random numbers: ",harvest

  do i=1,3
     call random_seed(get=seed_old)
     print*,"Present values of seed: ",seed_old
     call random_number(harvest)
     print*,"Other 3 random numbers: ",harvest
     call random_number(harvest)
     print*,"and other 3 random numbers: ",harvest
  end do

  deallocate(seed)
  deallocate(seed_old)

end program rantestbis_intrinsic
