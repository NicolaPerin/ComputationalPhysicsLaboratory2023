! rw2d.f90
! A simple random walk program  in 2D. 

PROGRAM drunk
  IMPLICIT NONE
  INTEGER :: i, N, sizer
  REAL :: phi, rnd
  INTEGER, DIMENSION(:), allocatable :: seed

  REAL :: x=0.0, y=0.0        ! Put drunk initially at the origin
  INTEGER, PARAMETER :: out=1 ! Set output unit

  REAL, PARAMETER :: step=1.0, twopi=2.0*3.1415926 ! step size and constants
  CHARACTER(LEN=15) :: filein
  CHARACTER(LEN=15), SAVE :: FORMAT1 = "(1i5,1x,2F14.7)"
  
  PRINT*,"Enter number of steps:"
  READ*, N
  PRINT*,"Enter file for data"
  READ*,filein
  OPEN(out, FILE=filein, STATUS="REPLACE", ACTION="WRITE")
  call random_seed(sizer)
  allocate(seed(sizer))
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  read(*,*)seed
  CALL RANDOM_SEED(PUT=seed)

  i = 0
  WRITE(UNIT=out,FMT=FORMAT1)i,x,y

  DO  i=1, N
     CALL RANDOM_NUMBER(rnd)
     phi=twopi*rnd
     x=x+step*COS(phi)
     y=y+step*SIN(phi)
     WRITE(UNIT=out,FMT=FORMAT1)i,x,y
  END DO

  CLOSE(out)

  deallocate(seed)
  stop
END PROGRAM drunk
