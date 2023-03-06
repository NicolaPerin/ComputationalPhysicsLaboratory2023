program rantest_es3_simple
!
! use intrinsic f90 random number generator and test quality
! through the calculation of several momenta, but always with one seed
!
implicit none
real, dimension (:), allocatable :: rnd, sum
integer :: N, i, k, kmax, klabel
character(len = 15) :: file_name

print*,' Insert how many random numbers >'
read(*,*)N
allocate (rnd(N))
call random_number(rnd)

print*,' Insert the maximum order of momentum >'
read(*,*)kmax
allocate(sum(kmax))

sum = 0.

do k = 1, kmax  ! loop for the different momenta
klabel = k + 10  ! shift by 10 the unit label
                 ! to avoid conflict with units 5 and 6 (stdin / stdout)
!     assign number.dat to file_name using write statement
write(unit=file_name,fmt="(i2.2,a)") k,".dat"
!     // is concatenation operator
file_name = "momentum"//file_name
open (klabel,file=file_name,action="write")

do i=1,N
sum(k) = sum(k) + rnd(i)**k
write(klabel,*)i,sum(k)/i, abs(sum(k)/i - 1./(k+1))
! sum(k)/i is the PARTIAL sum of the sequence for the momentum k
end do

close(klabel)
end do ! k

deallocate(rnd)

end program rantest_es3_simple

