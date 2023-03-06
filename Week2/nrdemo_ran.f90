module ran_module
  implicit none
  public :: ran_func

contains
  FUNCTION ran_func(idum) result(ran)
    !	IMPLICIT NONE
    INTEGER, PARAMETER :: K4B=selected_int_kind(9)
    INTEGER(kind=K4B) , intent(inout) :: idum
    REAL :: ran
    ! "minimal" random number generator;
    ! returns a uniform random deviate between 0.0 and 1.0 (not endpoints).
    ! Fully portable, scalar generator;
    ! has the "traditional" (NOT Fortran 90) calling sequence with
    ! a randm deviate as the returned function value:
    ! call with IDUM a NEGATIVE integer to initialize;
    ! thereafter, do not alter IDUM except to reinitialize.
    ! The period of this generator is about 3.1 * 10^18 
    INTEGER(kind=K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
    REAL, SAVE :: am
    INTEGER(kind=K4B), SAVE :: ix=-1,iy=-1,k
    if (idum <= 0 .or. iy < 0) then                ! initialize
       am=nearest(1.0,-1.0)/IM
       iy=ior(ieor(888889999,abs(idum)),1)
       ix=ieor(777755555,abs(idum))
       idum=abs(idum)+1                       ! set idum positive
    end if
    ix=ieor(ix,ishft(ix,13))  ! Marsaglia shift sequence, period 2^32-1
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    k=iy/IQ                   ! Park-Miller sequence, period 2^31-2
    iy=IA*(iy-k*IQ)-IR*k
    if (iy < 0) iy=iy+IM
    ran=am*ior(iand(IM,ieor(ix,iy)),1)  ! combine the two generators with
  END FUNCTION ran_func               ! masking to ensure nonzero value

end module ran_module

program demo

  use ran_module
  implicit none
  integer :: i,idum
  real :: x

  print*, "idum (<0) = "
  read*,idum
  x =ran_func(idum)
  print*,"Random number: ",x

  do i=1,10
     x = ran_func(idum)
     print*,"Random number: ",x
  end do

end program demo


