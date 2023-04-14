!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! metropolis_sampling.f90
!
! METROPOLIS sampling of several physical observables for the
! hamiltonian:         h = -1/2 \nabla^2 + 1/2 x^2),
! comparison exact expected results with numerical results
! on psi^2(x), with  psi(x) = exp(-x^2/(4\sigma^2))
! \sigma=1 => psi^2(x) = costant * standard gaussian
!  P(x) =  exp(-x**2/(2*sigma**2))/sqrt(2*pi*sigma**2)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program metropolis_sampling
  implicit none
  integer, parameter :: dp=selected_real_kind(13)
  integer :: i,n
  real(kind=dp):: sigma,etot,ekin,epot,rnd
  real(kind=dp):: x,x1,x2,x3,x4,xp,delta,expx,expxp,w,acc
  character(len=13), save :: format1 = "(a7,2x,2f9.5)"

  acc = 0.0_dp
  x1 = 0.0_dp
  x2 = 0.0_dp
  x3 = 0.0_dp
  x4 = 0.0_dp
  ekin = 0.0_dp
  epot = 0.0_dp
  print*, "n, sigma, x0, delta"
  read*,   n,sigma,x,delta

  do i=1,n
     ekin = ekin - 0.5_dp * ((x/(2*sigma**2))**2 - 1/(2*sigma**2))
     epot = epot + 0.5_dp * x**2
     etot = ekin + epot
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
  enddo

  write(unit=*,fmt=*)"acceptance ratio = ",acc/n
  write(unit=*,fmt=*)"# Results (simulation vs. exact results):"
  write(unit=*,fmt=format1)"etot = ",etot/n,1.0_dp/(8.0_dp*sigma**2)&
       +0.5_dp*sigma**2
  write(unit=*,fmt=format1)"ekin = ",ekin/n,1.0_dp/(8.0_dp*sigma**2)
  write(unit=*,fmt=format1)"epot = ",epot/n,0.5_dp*sigma**2
  write(unit=*,fmt=format1)"<x>  = ",x1/n,0.0_dp
  write(unit=*,fmt=format1)"<x^2>= ",x2/n,sigma**2
  write(unit=*,fmt=format1)"<x^3>= ",x3/n,0.0_dp
  write(unit=*,fmt=format1)"<x^4>= ",x4/n,3.0_dp*sigma**4

end program metropolis_sampling
