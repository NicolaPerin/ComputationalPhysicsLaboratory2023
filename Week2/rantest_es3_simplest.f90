program rantest_es3_simplest
!
! use intrinsic f90 random number generator and test quality
! through the calculation of one momentum and one sequence (= one seed)
!
implicit none
integer :: N, i, k
real :: sum
real, dimension (:), allocatable :: rnd

print*,' Insert how many random numbers >'
read(*,*)N
allocate (rnd(N))
call random_number(rnd)

print*,' Insert the order of momentum >'
read(*,*)k

sum = 0.

open (unit=1,file='momentumk.dat')

do i=1,N
sum = sum + rnd(i)**k
write(1,*)i,sum/i, abs(sum/i - 1./(k+1))
! sum/i is the PARTIAL sum of the sequence for the momentum k
end do

close(1)

deallocate(rnd)

end program rantest_es3_simplest

