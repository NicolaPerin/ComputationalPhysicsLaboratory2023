module isingf

  implicit none
  public :: init,metropolis,DeltaE
  integer, public, dimension(8) :: seed
  integer, public :: L,N,nmcs,nequil
  real (kind = 8), public, dimension(-8:8) :: w
  integer, public, dimension(:,:), allocatable :: spin
  real (kind = 8), public :: T,E,M
  integer, public :: accept
contains

  subroutine init()
  integer :: dE
   do dE = -8,8,4
       w(dE) = exp(-dE/T)
    end do
  end subroutine init

  subroutine alloc()
   allocate(spin(L,L))
  end subroutine alloc

  subroutine metropolis()
    !  one Monte Carlo step per spin
    integer :: ispin,x,y,dE

    real :: rnd
    do ispin = 1,N
       !     random x and y coordinates for trial spin
       call random_number(rnd)
       x = int(L*rnd) + 1
       call random_number(rnd)
       y = int(L*rnd) + 1
       dE = DeltaE(x,y)
       call random_number(rnd)
       if (rnd <= w(dE)) then
          spin(x,y) = -spin(x,y)
          accept = accept + 1
          M = M + 2*spin(x,y)  ! factor 2 is to account for the variation:
          E = E + dE           ! (-(-)+(+))
       end if
    end do
  end subroutine metropolis

  function DeltaE(x,y) result (DeltaE_result)
    !  periodic boundary conditions
    integer, intent (in) :: x,y
    integer :: DeltaE_result
    integer :: left
    integer :: right
    integer :: up
    integer :: down
    if (x == 1) then
       left = spin(L,y)
       right = spin(2,y)
    else if (x == L) then
       left = spin(L-1,y)
       right = spin(1,y)
    else
       left = spin(x-1,y)
       right = spin(x+1,y)
    end if
    if (y == 1) then
       up = spin(x,2)
       down = spin(x,L)
    else if (y == L) then
       up = spin(x,1)
       down = spin(x,L-1)
    else
       up = spin(x,y+1)
       down = spin(x,y-1)
    end if
    DeltaE_result = 2*spin(x,y)*(left + right + up + down)
! also here the factor 2 is to account for the variation
  end function DeltaE
end module isingf
