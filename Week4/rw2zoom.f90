! rw2d.f90
! A simple random walk program  in 2D. 

PROGRAM drunk
  IMPLICIT NONE
  INTEGER :: i, N, stepLength, sizer
  REAL :: phi, rnd
  INTEGER, DIMENSION(:), allocatable :: seed

  REAL :: x=0.0, y=0.0        ! Put drunk initially at the origin
  INTEGER, PARAMETER :: out=1, out1=2 ! Set output unit

  REAL, PARAMETER :: step=1.0, twopi=2.0*3.1415926 ! step size and constants
  CHARACTER(LEN=15), SAVE :: FORMAT1 = "(1i5,1x,2F14.7)"
  
  call random_seed(sizer)
  allocate(seed(sizer))
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  read(*,*)seed
  call random_seed(put=seed)

  PRINT*,"Enter number of steps:"
  READ*, N
  PRINT*,"Enter  stepLenght (integer=>1) :"
  READ*, StepLength
  OPEN(out, FILE='10.dat', STATUS="REPLACE", ACTION="WRITE")
  OPEN(out1, FILE='1.dat', STATUS="REPLACE", ACTION="WRITE")

  i = 0
  WRITE(UNIT=out,FMT=FORMAT1)i,x,y
  WRITE(UNIT=out1,FMT=FORMAT1)i,x,y

  DO  i=1, N*stepLength
     CALL RANDOM_NUMBER(rnd)
     phi=twopi*rnd
     x=x+step*COS(phi)
     y=y+step*SIN(phi)
     if(i<N+1)WRITE(unit=out1,FMT=FORMAT1)i,x,y
     if(mod(i,stepLength)==0)WRITE(UNIT=out,FMT=FORMAT1)i,x,y
  END DO

  CLOSE(out)
   CLOSE(out1)

 deallocate(seed)
 stop
END PROGRAM drunk
