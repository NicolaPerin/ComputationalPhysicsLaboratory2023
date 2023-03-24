PROGRAM checkloop
  IMPLICIT none
  INTEGER :: i, max
  print *,"initial upper bound loop>"
  read *, max
     print*,' max before loop:',max

     DO  i = 1, max   
           IF (i==3) max = max - 2
           print*,' i, max inside loop:',i, max
     END DO

STOP
END program checkloop
