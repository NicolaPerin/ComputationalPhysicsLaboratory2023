!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! gauss_metropolis.f90
!
! METROPOLIS generation of random numbers with a Gaussian distribution
! P(x) = exp(-x**2/(2*sigma**2))/sqrt(2*pi*sigma**2)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program gauss_metropolis
  implicit none
  integer, parameter :: dp=selected_real_kind(13)
  integer :: i,n,ibin,maxbin,m
  real(kind=dp):: sigma,rnd,delta,x0,deltaisto
  real(kind=dp):: x,x1,x2,x3,x4,xp,expx,expxp,w,acc
  real, dimension(:), allocatable :: istog
  character(len=13), save :: format1 = "(a7,2x,2f9.5)"

  print*,' insert n, sigma, x0, delta, max number of bin in the histogram  >'
  read*, n, sigma,x0,delta,maxbin
  allocate(istog(-maxbin/2:maxbin/2))
  istog = 0.
  deltaisto = 10.*sigma/maxbin  ! histogram over a range of 10*sigma
  acc = 0.0_dp
  x = x0
  x1 = 0.0_dp
  x2 = 0.0_dp
  x3 = 0.0_dp
  x4 = 0.0_dp

  do i=1,n
     x1 = x1 + x
     x2 = x2 + x**2
     x3 = x3 + x**3
     x4 = x4 + x**4
     !ccccccccccccccccccccccccccccccc
     expx = - x**2 /(2*sigma**2)    !
     call random_number(rnd)        !
     xp = x + delta * (rnd-0.5_dp)  !
     expxp = - xp**2 /(2*sigma**2)  !   metropolis
     w = exp (expxp-expx)           !   algorithm
     call random_number(rnd)        !
     if (w > rnd) then              !
        x = xp                      !
     !ccccccccccccccccccccccccccccccc
        acc=acc+1.0_dp                
     endif
     ibin = nint(x/deltaisto)
     if (abs(ibin) < maxbin/2) istog(ibin) = istog(ibin) + 1
  enddo

  write(unit=*,fmt=*)"# n, x0, delta = ",n,x0,delta
  write(unit=*,fmt=*)"# acceptance ratio = ",acc/n
  write(unit=*,fmt=*)"# Results (simulation vs. exact results):"
  write(unit=*,fmt=format1)"# <x>  = ",x1/n,0.0_dp
  write(unit=*,fmt=format1)"# <x^2>= ",x2/n,sigma**2
  write(unit=*,fmt=format1)"# var2 = ",x2/n-(x1/n)**2,sigma**2
  write(unit=*,fmt=format1)"# <x^3>= ",x3/n,0.0_dp
  write(unit=*,fmt=format1)"# <x^4>= ",x4/n,3.0_dp*sigma**4

  open(1,file='gauss_metropolis.dat',status='replace')
  write(unit=1,fmt=*)"# n, x0, delta = ",n,x0,delta
  do ibin = -maxbin/2 , maxbin/2
     write(1,*)ibin*deltaisto, istog(ibin)/real(n)/deltaisto
  end do

  close(1)
  deallocate(istog)

end program gauss_metropolis
