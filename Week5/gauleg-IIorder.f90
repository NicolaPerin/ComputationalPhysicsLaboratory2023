!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  gauleg_test.f90
!  test Gauss-Legendre quadrature II-order formula in a generic interval (a,b);!  Here: integrand function is f(x)=exp(x), but it is easy to change
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module gaussint
  public :: f
contains

  function f(X)
    ! integrand function
    implicit none
    real :: f
    real, intent(in) :: X
    F=exp(x)
  end function f

  subroutine gausquad(f,a,b,estimate)
    real :: f
    real, intent(in) :: a,b
    real, intent(out) :: estimate
    estimate=(b-a)/2*((f(-(b-a)/2/sqrt(3.)+(b+a)/2.)+f((b-a)/2/sqrt(3.)+(b+a)/2.)))
  end subroutine gausquad

end module gaussint

program gauleg_test
  use gaussint
  !
  !     Gauss-Legendre Quadrature
  !
  implicit none
  real :: a,b,estimate, exact
  character :: answ*1
  integer :: J
  write(*,*)' Gauss-Legendre II-order quadrature for exp(x) in [a,b] interval:'
  write(*,*)' Int( F(X) DX  )~  SUM_(k=1)^N w_k F(x_k ) '
10 write(*,*)' insert a, b > '
  read (*,*) a,b
  call GAUSQUAD(f,a,b,estimate)
  write(*,*)' exact value        : ',exp(b)-exp(a)
  write(*,*)' numerical estimate : ', estimate
  write(*,*)' error              : ', estimate-(exp(b)-exp(a))
  write(*,*)' other interval ? <Y/N> '
  read (*,'(A)') answ
  if (answ=='Y' .or. answ=='y') goto 10
  stop
end program gauleg_test
