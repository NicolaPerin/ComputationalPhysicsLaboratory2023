program rantest_es3_bruteforce
!
! use intrinsic f90 random number generator and test quality
! through the calculation of one momentum;
! sequences of different lenghts are generated with different seeds
! time-consuming but more "randomized" features
!
implicit none
real, dimension (:), allocatable :: rnd
real :: somma
integer :: N, i, k

print*,' Insert the maximum length of the sequence >'
read(*,*)N

print*,' Insert the order of momentum >'
read(*,*)k

do i=1,N
allocate (rnd(i))

call random_number(rnd)  ! generate a new sequence of "i" random numbers
                         ! (seed changes automatically)
somma = sum(rnd**k)
write(1,*)i,somma/i, abs(somma/i - 1./(k+1))
! somma/i is the PARTIAL sum of the sequence for the momentum k
deallocate(rnd)

end do  ! i

close(1)

end program rantest_es3_bruteforce

