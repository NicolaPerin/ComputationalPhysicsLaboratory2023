program test_expdev
  implicit none
  real :: lambda,delta,x
  integer :: i,n,nbin,ibin, sizer
  integer, dimension(:), allocatable :: histo, seed
  print*, " Generates random numbers distributed as exp(-x*lambda)"
  call random_seed(sizer)
  allocate(seed(sizer))
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  read(*,*)seed
  call random_seed(put=seed)
  print *," length of the sequence >"
  read *, n
  print *," exponential decay factor (lambda)>"	
  read *, lambda	
  print *," Collecting numbers generated up to 2/lambda (disregard the others)"
  print *," and normalizing the distribution in [0,+infinity[ "
  print *," Insert number of bins in the histogram>"
  read *, nbin
  delta = 2./lambda/nbin
  allocate (histo(nbin))
  histo = 0
  do i = 1,n
     call expdev(x)
     ibin = int (x/lambda/delta) + 1
     if (ibin <= nbin)histo(ibin) = histo(ibin) + 1
  end do
  open (unit=7,file="expdev.dat",status="replace",action="write")
  do ibin= 1 ,nbin
     write(unit=7,fmt=*)(ibin-0.5)*delta,histo(ibin)/float(n)/delta
  end do

contains

  subroutine expdev(x)
    REAL, intent (out) :: x
    REAL :: r
    do
       call random_number(r)
       if(r > 0) exit
    end do
    x = -log(r)
  END subroutine expdev

end program test_expdev
