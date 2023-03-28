program differenze
  implicit none
  integer, parameter :: dp=kind(1.d0)
integer :: i
  real(dp) :: diff, z, z1

  z = 1.e-16_dp
  z1 = 0.9e-16_dp
  do  i = 1, 10
    print*,'z,z1,abs(z-z1)',z,z1,abs(z-z1)
    z = z/10
    z1 = z1/10
  end do

end program differenze
