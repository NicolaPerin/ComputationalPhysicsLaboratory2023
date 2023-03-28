!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       pi.f90: Calculate pi using Monte Carlo
!C       ("acceptance/rejection" in a circle)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Program pi

  Implicit none
  integer, dimension(:), allocatable :: seed
  real, dimension(2) :: rnd
  Real :: area, x, y
  Integer :: i, max, pigr, sizer
  call random_seed(sizer)
  allocate(seed(sizer))
    print*,' max number of points >'
  read*, max
  print*,' enter seed (or type /) >'
  read*, seed
  call random_seed(put=seed)
  ! open file to save data, initialize the counter of random number points in the circle
  Open(7, File='pigr.dat', Status='Replace')
  pigr=0.
  ! random numbers generated between -1 and 1
  ! i.e., within a square of side 2 and centered in (0,0):
  ! points falling inside the circle (with x*x+y*y < or = 1) are counted for pigr, 
  ! which is 4*area(circle)/area(square)
  ! i.e. 4*number of points within the circle / total number of points withing the square
  Do i=1, max
     call random_number(rnd)
     x = rnd(1)*2-1
     y = rnd(2)*2-1
     If ((x*x + y*y) <= 1) then
        pigr = pigr+1
     Endif
     area = 4.0 * pigr/Real(i)
     if (mod(i,10)==0) write(7,*) i, abs(acos(-1.)-area)  !write every 10 points
  end do
  Close(7)
  Stop 'data saved in pigr.dat '
End program pi



















