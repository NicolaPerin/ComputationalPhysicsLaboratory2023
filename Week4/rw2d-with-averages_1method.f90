! rw2d.f90
! A simple random walk program  in 2D. 

PROGRAM drunk
  IMPLICIT NONE
  INTEGER :: istep, irun, N, Nruns, sizer
  REAL :: rnd, x, y, phi
  REAL, DIMENSION(:), allocatable :: x2av, y2av, xav, yav
  INTEGER, DIMENSION(:), allocatable :: seed
  INTEGER, PARAMETER :: out=1, out2=2 ! Set output units
  REAL, PARAMETER :: step=1.0, twopi=2*acos(-1.0) ! constants
  CHARACTER(LEN=15) :: filein, filein2
  CHARACTER(LEN=15), SAVE :: FORMAT1 = "(1i5,1x,2F14.7)"
  
  PRINT*,"Enter number of steps >"
  READ*, N
  PRINT*,"Enter number of walkers >"
  READ*, Nruns
  PRINT*,"Enter filename to store i,x,y of each walker >"
  READ*,filein
  OPEN(out, FILE=filein, STATUS="REPLACE", ACTION="WRITE")
  PRINT*,"Enter filename to store i,<deltar2> averaged over the walkers >"
  READ*,filein2
  OPEN(out2, FILE=filein2, STATUS="REPLACE", ACTION="WRITE")
  call random_seed(sizer)
  allocate(seed(sizer), x2av(N), y2av(N), xav(N), yav(N))
  x2av=0.0; y2av=0.0; xav=0.0; yav=0.0
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  read(*,*)seed
  CALL RANDOM_SEED(PUT=seed)


runs:  do irun = 1, Nruns
   x = 0.0; y=0.0 ! initial position of each walker
  
steps: do istep = 1, N
     call random_number(rnd) ! get a random number
     phi=twopi*rnd
     x=x+step*COS(phi)   ! istantaneous position x of one walker
     y=y+step*SIN(phi)   ! istantaneous position y of one walker
     WRITE(UNIT=out,FMT=FORMAT1)istep,x,y
     xav (istep) = xav (istep) + x     ! accumulate quantities over the runs
     x2av(istep) = x2av(istep) + x**2 
     yav (istep) = yav (istep) + x
     y2av(istep) = y2av(istep) + x**2   
  end do steps

     WRITE(UNIT=out,FMT=*)' '
     WRITE(UNIT=out,FMT=*)' '
     print*,' end of the walker n.',irun
end do runs
     
   print*,' end of the loop over the walkers'
  
   DO istep = 1, N
   WRITE(UNIT=out2,FMT=*)istep,x2av(istep)/nruns+y2av(istep)/nruns-(xav/nruns)**2-(yav/nruns)**2
   END DO
  
   print*,' data saved'
 

  CLOSE(out)
  CLOSE(out2)

  deallocate (seed, x2av, y2av, xav, yav)

  stop
END PROGRAM drunk
