! gauleg_nr_test.f90
!
! integra f(x)=exp(x) nell'intervallo voluto usando Gauss-Legendre
! e usando routines/modules/interfaces di Numerical Recipes in F90
! (qui semplificate in modo da tenere solo le cose necessarie)
!
program gauleg_nr_test
  use nrtype
  use nr  
  implicit none
  real(sp) :: a,b,quad
  real(sp), dimension(:), allocatable :: x,w
  integer :: npts,i
  print*," insert a, b, npoints>"
  read*, a,b,npts
  allocate(x(npts))    ! l'allocazione di memoria per x e w (ascisse e pesi)
  allocate(w(npts))    ! indica alla subroutine GAULEG quanti punti usare
  call gauleg(a,b,x,w)  
  do i=1,npts
     print*,x(i),w(i)
     quad=quad+exp(x(i))*w(i)
  end do
  print *,' risultato=',quad,' errore=',quad-(exp(b)-exp(a))
end program gauleg_nr_test
