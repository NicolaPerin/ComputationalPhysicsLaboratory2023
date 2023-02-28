      program area
	implicit none
!     
!  calcola l'area del cerchio, dato il raggio
!
!     INPUT : r (raggio)
!     OUTPUT: A (area)
!     costante: pi (pigreco)
!
      real ::  r ,A, pi   ! dichiarazione variabili e costanti.
      pi = 4*atan(1.0)           ! calcola pi.
      print*," raggio ?"          ! sul terminale (*) appare la domanda.
      read *, r                 ! legge r dal terminale 
                                  ! il primo * sta per 'terminale',
                                  ! il secondo * sta per 'formato libero'.
      A = pi * r**2               ! calcola l'area
      print "(a,f10.5,a,f12.7)"," raggio = ",&
     & r, " Area = ", A    ! scrive sul
                                                   ! terminale (*) con un
                                                   ! formato preciso

! 10   format (a10,f10.5,5x,a8,f12.7)               ! formato in f77
                                                   ! programma usa per
                                                   ! scrivere i
                                                   ! risultati sul
                                                   ! terminale
!      stop
      end program area
                                        
      
