program prec_mac_double
implicit none

integer, parameter :: dp=selected_real_kind(13)
real(kind=dp) :: eps = 1._dp, one = 1._dp, zero = 0._dp
integer :: n = 0

do n = 1, 100
  eps = eps / 2._dp
  one = eps + 1._dp
  print*, n, eps + 1._dp , one, eps
  if ( eps + 1._dp == 1._dp)   print*,"at iteration",n," 1+eps=1"
  if ( one == 1._dp) print*,"at iteration",n," one=1 (defined: one=1+eps)"
  if ( eps + 1._dp == 1._dp .and. one == 1._dp) then
  print*," stop at iteration",n; exit
  end if
end do

print*," machine precision:", epsilon(1._dp)

end program prec_mac_double

