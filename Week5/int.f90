!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     int.f90: 
!     integrates f(x)=exp(x) in the interval [vmin,vamx]=[0,1]
!     using trapezoidal and Simpson rule
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module intmod
  public :: f, trapez, simpson
contains

  ! function to be integrated
  !
  function f(x)
    implicit none
    real :: f
    real, intent(in) :: x
    f = exp(x)
    return
  end function f


  ! trapezoidal rule
  !
  function trapez(i, min, max)
    implicit none
    real :: trapez
    integer, intent(in) :: i
    real, intent(in) :: min, max
    integer :: n
    real :: x, interval
    trapez = 0.
    interval= ((max-min) / (i-1))
    ! sum over the internal points (excluding the end-points)
    do n = 2, i-1
       x = interval * (n-1)
       trapez = trapez + f(x) * interval
    end do
    ! add the end-points
    trapez = trapez + 0.5 * (f(min)+f(max)) * interval
    return
  end function trapez


  ! Simpson rule
  !
  function simpson(i, min, max)
    implicit none
    real :: simpson
    integer, intent(in) :: i
    real, intent(in) :: min, max
    integer :: n
    real :: x, interval
    simpson = 0.
    interval = ((max-min) / (i-1))
    ! loop over the internal EVEN points
    do n = 2, i-1, 2
       x = interval * (n-1)
       simpson = simpson + 4*f(x)
    end do
    ! loop over the internal ODD points
    do n = 3, i-1, 2
       x = min + interval * (n-1)
       simpson = simpson + 2*f(x)
    end do
    ! add the end-points
    simpson = simpson + f(min)+f(max)
    simpson = simpson * interval/3
    return
  end function simpson

end module intmod

program int
  use intmod
  !
  !     variable declaration
  !     accuracy limit
  !     min and max in x
  !
  implicit none
  real :: r1, r2, theo, vmin, vmax, t0, t1
  integer :: i, n
  ! exact value 
  vmin = 0.0
  vmax = 1.0
  theo = exp(vmax)-exp(vmin)
  print*,' exact value =',theo
  open(unit=7,file='int-tra-sim.dat',status='unknown')
  !
  write(7,*)"# N,   interval,    exact,   Trap-exact, Simpson-exact"
  call cpu_time(t0)
  do i = 2,8
     n = 2**i
     r1 = trapez(n+1, vmin, vmax)
     r1 = (r1-theo)
     r2 = simpson(n+1, vmin, vmax)
     r2 = (r2-theo)
     write(7,'(i4,4(2x,f10.6))') n, 1./n, theo, r1, r2
  end do
  call cpu_time(t1)
  print*," total time spent:",t1-t0
  close(7)

  print*,' data saved  in int-tra-sim.dat (|diff from exact value|)'
  stop

end program int
