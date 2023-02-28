!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     exp-bad.f : un cattivo algoritmo per calcolare l'esponenziale e^-x
!                 come somma FINITA della serie
!                 (da confrontarsi con exp-good.f 
!                  e con la funzione intrinseca della macchina)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      program expbad
!
! dichiarazioni di variabili:
!     limite di accuratezza,
!     minimo e massimo in x, step in x, numeratore (up), denominatore (down)
!
      implicit none
      integer, parameter :: dp = selected_real_kind(13)
      real(kind=dp) :: down, min, max, step, sum, up, x
      integer :: i, n
      min = 1.e-10_dp
      max = 10.0_dp
      step =0.1_dp
!
! esecuzione
!     
      x=0.0
      do  
         if(x>max)exit
         x=x+step
         sum = 1
         do  n=1, 10000
            up   = 1
            down = 1
            do  i=1,n
               up   = -up*x
               down = down*i
            enddo
            sum = sum + up/down
            if((abs(up/down) < min) .and. (sum /= 0)) then
               print*,x,sum,exp(-x),up/down
               go to 10
            endif
        enddo
 10    continue
      enddo
      end program expbad
