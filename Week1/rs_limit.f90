      program rs_limit
!c
!c determina la precisione macchina (eps) entro un certo fattore prestabilito 
!c IN SINGOLA PRECISIONE
!c
!c La precisione macchina e' definita come il MASSIMO NUMERO POSITIVO eps
!c t.c. aggiunto al numero immagazzinato come unita' nella memoria della
!c macchina (one) non varia tal numero (cioe' one = one + eps)
!c
!c allora per determinare "eps" aggiungo un numero sempre piu' piccolo 
!c a "one" finche' questo non succede.
!c
      real :: eps, factor, one     ! dichiarazione variabili
      integer :: i, N                ! e costanti
!c
!c
      eps = 1.0                   ! inizializzo le variabili
      factor = 2.0                ! e le costanti
      N = 100                     !
      do i=1,N                       !\   DO LOOP
         eps = eps/factor            ! |  aggiungo a "one" una quantita' "eps"
         one = 1.0 + eps             !  > sempre piu' piccola e confronto
         print*,  i, one, eps        ! |  "one" e "one + eps"
      end do                         !/
      stop
      end program rs_limit
      
