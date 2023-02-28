      program rd_under_over
!c  
!c per un test sui limiti di UNDERFLOW e OVERFLOW della macchina
!c usando 64 bit  (=8 bytes, DOPPIA precisione) 
!c per la rappresentazione dei REALI 
!c
	integer, parameter :: dp=selected_real_kind(13)
      real(kind=dp) :: under, over, factor     ! dichiarazione variabili
      integer :: i, N                   ! e costanti
!c
!c
      under = 1.0_dp                   ! inizializzo  variabili e costanti
      over  = 1.0_dp                   !
      factor = 2.0_dp                  ! fattore con cui moltiplicare e dividere
      N      = 2000
      do i=1,N                       !\   
         under = under/factor        ! |  DO LOOP
         over  =  over*factor        !  > moltiplico(divido) iterativamente
         print* ,   i, under, over   ! |  per avere una quantita' sempre
      end do                         !/   piu' grande(piccola)
      stop
      end program rd_under_over
