!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     exp-good.f : un BUON algoritmo per calcolare l'esponenziale e^-x
!                 come somma FINITA della serie
!                 (da confrontarsi con exp-bad.f 
!                  e con la funzione intrinseca della macchina)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program expgood
!c
!c dichiarazioni di variabili:
!c     limite di accuratezza,
!c     minimo e massimo in x, step in x
!c
      implicit none
	integer, parameter :: dp = selected_real_kind(13)
      real(kind=dp) ::  element, min, max, step, sum, x
      integer :: n
      min = 1.e-10_dp
      max = 10.0_dp
      step =0.1_dp
      open(unit=7,file="exp-good.dp.dat",status="replace",action="write")
               write(unit=7,fmt=*) "x, n, sum, exp(-x), abs(sum-exp(-x))/sum" 
!c
!c esecuzione
!c
	x = 0.0
      do 
 	x = x + step
	if (x>max) exit
         sum     = 1
         element = 1
         do  n=1, 10000
            element = element*(-x)/n
            sum = sum + element
            if((abs(element/sum) < min) .and. (sum /= 0)) then
               write(unit=7,fmt=*) x, n, sum, exp(-x), abs(sum-exp(-x))/sum 
               go to 10
            endif
         enddo
 10   continue
       end do
       close(7)
!      stop "dati salvati nel file exp-good.dat"	
	end program expgood
