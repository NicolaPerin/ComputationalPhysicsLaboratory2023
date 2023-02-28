program test

! This program tests passing a function name into a subroutine.

implicit none
real, dimension(10) :: x
real :: xbar,s
integer :: i

! This INTERFACE is required in order to identify the function
! names that we're about to hand off to another subroutine which
! will subsequently use those functions.

interface
  function avg(x)
    real :: avg
    real, dimension(:), intent(in) :: x
  end function
  function std(x)
    real :: std
    real, dimension(:), intent(in) :: x
  end function
end interface

! Define the data.  The implied looping variable has to be an
! integer!!

x=(/ ((1.0*i),i=1,10) /)

! Now compute and print out the sample mean and sample
! standard deviation.

call compute(avg,x,xbar)
print*,xbar
call compute(std,x,s)
print*,s

! This statement tells Fortran that *internal* subroutines and 
! functions follow.

contains

!********** The Internal Subroutine **********

  subroutine compute(func,data,ans)

! This is a silly little routine that simply demonstrates how
! you pass a function name which will be used in the subroutine.

  implicit none
  real, dimension(:), intent(in) :: data
  real, intent(out) :: ans

! We have to let Fortran know that FUNC is a dummy name for a
! function to be passed in later.

  interface
    function func(x)
      real :: func
      real, dimension(:), intent(in) :: x
    end function
  end interface

! Use that FUNC.

  ans=func(data)

  end subroutine compute

end program test

!********** These are External Functions **********

function avg(x)

! This function computes an average.

implicit none
real, dimension(:), intent(in) :: x
real :: avg

avg=sum(x)/size(x)

end function

!**********

function std(x)

implicit none
real, dimension(:), intent(in) :: x
real :: std

interface
  function avg(x)
    real :: avg
    real, dimension(:), intent(in) :: x
  end function
end interface

std=( dot_product(x,x)-size(x)*avg(x)*avg(x) )/(size(x)-1.0)
std=sqrt(std)

end function

