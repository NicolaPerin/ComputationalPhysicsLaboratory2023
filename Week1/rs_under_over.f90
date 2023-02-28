      program rs_under_over
!c  
!c per un test sui limiti di UNDERFLOW e OVERFLOW della macchina
!c usando 32 bit  (=4 bytes, SINGOLA precisione) 
!c per la rappresentazione dei REALI 
!c

      real :: under, over, factor     ! dichiarazione variabili
      integer :: i, N                   ! e costanti
!c
!c
      under = 1.0                   ! inizializzo  variabili e costanti
      over  = 1.0                   !
      factor = 2.0                  ! fattore con cui moltiplicare e dividere
      N      = 500
      do i=1,N                       !\   
         under = under/factor        ! |  DO LOOP
         over  =  over*factor        !  > moltiplico(divido) iterativamente
         print* ,   i, under, over   ! |  per avere una quantita' sempre
      end do                         !/   piu' grande(piccola)
      stop
      end program rs_under_over
