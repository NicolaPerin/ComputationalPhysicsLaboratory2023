!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! boltzmann_metropolis.f90
!
! Metropolis algorithm used as importance-sampling:
! generation of microstates with Boltzmann distribution,
! here for a classical particle in 1D.
! The interesting quantity is the probability P(E)dE for a particle
! to have energy between E and E+dE (here E can label a microstate,
! a part from the sign +/- of the velocity)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module common
  implicit none
  public :: initial, Metropolis, data, probability, averages
  real, public :: E,T,del_E,beta,dvmax,vel,accept
  integer, public, dimension(:), allocatable :: seed
  integer, public :: nbin,nmcs,sizer
  real, public, dimension(:), allocatable :: P
contains

  subroutine initial(nequil,vcum,ecum,e2cum)
    real, intent(out) :: vcum,ecum,e2cum
    integer, intent(out) :: nequil
    print*," number of MC steps >"
    read *, nmcs
    print*," absolute temperature >"
    read *, T
    print*," initial velocity >"
    read *, vel
    print*," maximum variation of the velocity (hint: 4*sqrt(T)=",4*sqrt(T),") >"
    read *, dvmax
    call random_seed(sizer)
    allocate(seed(sizer))
    print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
    read *, seed
    call random_seed(put=seed)
    beta   = 1/T
    nequil = 0.1 * nmcs   ! WARNING : VERIFY this choice !
    E      = 0.5 * vel * vel
    del_E  = T/20         ! a reasonable width of the bin for the histogram of P(E)
    nbin   = int(4*T / del_E)  ! max. number of bins
    print *,"# T       :",T
    print *,"# <E0>    :",E
    print *,"# <v0>    :",vel
    print *,"# dvmax   :",dvmax
    print *,"# nMCsteps:",nmcs
    print *,"# deltaE  :",del_E
    print *,"# nbin    :",nbin
    open(unit=9,file="boltzmann.dat",status="replace",action="write")
    write(unit=9,fmt=*)"# T       :",T
    write(unit=9,fmt=*)"# <E0>    :",E
    write(unit=9,fmt=*)"# <v0>    :",vel
    write(unit=9,fmt=*)"# dvmax   :",dvmax
    write(unit=9,fmt=*)"# nMCsteps:",nmcs
    write(unit=9,fmt=*)"# deltaE  :",del_E
    write(unit=9,fmt=*)"# nbin    :",nbin
    allocate (P(0:nbin))
    ecum  = 0.0
    e2cum = 0.0
    vcum  = 0.0
    P     = 0.0
    accept= 0.0
  end subroutine initial

  subroutine Metropolis() 
    real :: dv,vtrial,de,rnd
    call random_number(rnd)
    dv = (2*rnd - 1) * dvmax              ! trial variation for v
    vtrial = vel + dv                     ! trial velocity v
    de = 0.5 * (vtrial*vtrial - vel*vel)  ! corresponding variation of E
    call random_number(rnd)
    if (de >= 0.0) then
       if ( exp(-beta*de) < rnd ) return ! trial step not accepted
    end if
    vel = vtrial
    accept = accept + 1 
    E = E + de
  end subroutine Metropolis

  subroutine data(vcum,ecum,e2cum)
    real, intent(inout) :: vcum,ecum,e2cum
    Ecum = Ecum + E
    E2cum = E2cum + E*E
    vcum = vcum + vel
    call probability() 
  end subroutine data

  subroutine probability() 
    integer :: ibin
    ibin = int(E/del_E)

    if ( ibin <= nbin )     P(ibin) = P(ibin) + 1
  end subroutine probability

  subroutine averages(nequil,vcum,Ecum,E2cum) 
    integer, intent(in) :: nequil
    real, intent(in) :: vcum,Ecum,E2cum
    real :: znorm, Eave, E2ave, vave, sigma2
    real, parameter :: pigr = acos(-1.)
    integer :: ibin
    znorm  = 1.0/nmcs
    accept = accept / (nmcs+nequil) ! acceptance ratio
    Eave   = Ecum * znorm   ! average energy
    E2ave  = E2cum * znorm  ! 
    vave   = vcum * znorm   ! average velocity
    sigma2 = E2ave - Eave*Eave
    print *,"# <E2>num.:",E2ave
    print *,"# <E> num.:",Eave
    print *,"# <E> th. :",T/2
    print *,"# <v>     :",vave
    print *,"# accept. :",accept
    print *,"# sigma   :",sqrt(sigma2)
    write(unit=9,fmt=*)"# <E2>num:",E2ave
    write(unit=9,fmt=*)"# <E> num.:",Eave
    write(unit=9,fmt=*)"# <E> th. :",T/2
    write(unit=9,fmt=*)"# <v>     :",vave
    write(unit=9,fmt=*)"# accept. :",accept
    write(unit=9,fmt=*)"# sigmaE  :",sqrt(sigma2)
    write(unit=9,fmt=*)"# ibin*del_E, P(E) simul, P(E) th"
    do ibin = 0,nbin
      write(unit=9,fmt=*) ibin*del_E, P(ibin) * znorm / del_E,&
&   exp(-ibin*del_E/T)/sqrt(pigr*T*ibin*del_E)
    end do
    close(unit=9)
  end subroutine averages
end module common

program Boltzmann
  use common
  real :: vcum, ecum, e2cum
  integer :: imcs,nequil
  ! parameters and variable initialization
  call initial(nequil,vcum,ecum,e2cum)
  do  imcs = 1 , nmcs + nequil
     call Metropolis()
     ! data accumulation after each Metropolis step
     if ( imcs > nequil ) call data(vcum,ecum,e2cum) 
  end do
  call averages(nequil,vcum,Ecum,E2cum)
  deallocate(P)
end program Boltzmann

