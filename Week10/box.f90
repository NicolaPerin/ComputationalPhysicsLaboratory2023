!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! box.f90
!
! simulation of the evolution of a physical system towards equilibrium:
! non interacting particles in a box divided into two parts; 
! at each time step,  one and only one particle (randomly choosen) 
! goes from one side to the other one
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module moduli_box
  implicit none
  public :: initial, move
  integer, public :: N,tmax
contains
  subroutine initial()
    integer , dimension(4) :: seed
    print*," total number of particles N >"
    read*,N
    tmax = 10*N  ! we choose the evolution time proportional to N 
    print*," seed (1:4) >"
    read *, seed
    call random_seed(put=seed)
  end subroutine initial

  subroutine move()
    integer :: nl,itime
    real :: r, prob
    nl = N ! we start with all the particles on the left side
    open(unit=2,file="box.out",action="write",status="replace")
    do itime = 1,tmax
       prob = real(nl)/N    ! fraction of particles on the left
       call random_number(r)  
       if (r <= prob) then
          nl = nl - 1
       else
          nl = nl + 1
       end if
       write(unit=2,fmt=*)nl
    end do
    close(unit=2)
  end subroutine move
end module moduli_box

program box
  use moduli_box
  ! compare a random number with the fraction of particles on the left, nl/N:
  ! if r.le.nl/N we move one particle from tyhe left to the right;
  ! elsewhere from the right to the left
  call initial()
  call move()
end program box


