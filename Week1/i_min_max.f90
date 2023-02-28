      program i_min_max
!c
!c calcola il minimo e massimo intero possibile per una macchina a 32 bit
!c e in una rappresentazione degli interi a 4 bytes 
!c
      implicit none
      integer :: imin, imax, iexp
!c
      print*," massimo intero positivo:"
      iexp = 1
      imax = 1
      do 
         imax = 2**iexp - 1
         print*, " iexp=",iexp,"     2**iexp - 1=",imax
	 if (imax <= 0 ) go to 10
         imax = 2**iexp
         print*, " iexp=",iexp,"     2**iexp    =",imax
	 if (imax <= 0 ) go to 10
         iexp = iexp + 1
         end do
10	continue
!c
      print*, " minimo intero negativo:"
      iexp = 1
      imin = -1
      do 
         imin = - 2**iexp + 1
         print*, " iexp=",iexp,"    -2**iexp + 1=",imin
         if (imin >= 0)exit
         imin = - 2**iexp
         print*, " iexp=",iexp,"    -2**iexp    =",imin
         iexp = iexp + 1
         if (imin >= 0)exit
         end do
!C
      stop
      end program i_min_max
