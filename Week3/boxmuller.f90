! boxmuller.90
! uses the Box-Muller algorithm to generate 
! a random variate with a gaussian distribution (sigma = 1)
!
program boxmuller
  implicit none
  real :: rnd,delta
  real, dimension(:), allocatable :: histog
  integer :: npts,i,ibin,maxbin,m

  print*,' input npts, maxbin >'
  read*, npts,maxbin
  allocate(histog(-maxbin/2:maxbin/2))
  histog = 0
  delta = 10./maxbin
  do i = 1, npts
     call gasdev(rnd)
     ibin = nint(rnd/delta)
     if (abs(ibin) < maxbin/2) histog(ibin) = histog(ibin) + 1
  end do

  open(1,file='gasdev.dat',status='replace')
  do ibin = -maxbin/2 , maxbin/2
     write(1,*)ibin*delta, histog(ibin)/real(npts)/delta
  end do
  close(1)
  deallocate(histog)
  stop

contains

  SUBROUTINE gasdev(rnd)
    IMPLICIT NONE
    REAL, INTENT(OUT) :: rnd
    REAL :: r2,x,y
    REAL, SAVE :: g
    LOGICAL, SAVE :: gaus_stored=.false.
    if (gaus_stored) then
       rnd=g
       gaus_stored=.false.
    else
       do
          call random_number(x)
          call random_number(y)
          x=2.*x-1.
          y=2.*y-1.
          r2=x**2+y**2
          if (r2 > 0. .and. r2 < 1.) exit
       end do
       r2=sqrt(-2.*log(r2)/r2)
       rnd=x*r2
       g=y*r2
       gaus_stored=.true.
    end if
  END SUBROUTINE gasdev

end program boxmuller

