!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! a program which calculates the factorial using a recursive function;
! use of module
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module fact

public :: f
contains

recursive function f(n) result (factorial_result)
   integer, intent (in) :: n
   integer :: factorial_result

   if (n <= 0) then
      factorial_result = 1
   else
      factorial_result = n*f(n-1)
   end if
end function f

end module fact

program test_factorial
   use fact

   integer :: n
   print *, "integer n?"
   read *, n
   print "(i4, a, i10)", n, "! = ", f(n)
end program test_factorial

