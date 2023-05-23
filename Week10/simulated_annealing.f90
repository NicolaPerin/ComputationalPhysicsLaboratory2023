PROGRAM anneal
! Simulated annealing for minimization of a function
! Metropolis algorithm, J. Chem. Phys. 21 (1953) 1087
! Ulrich Schmitt, 2003-01-15
IMPLICIT NONE
INTEGER :: istep, nsteps
REAL, PARAMETER :: scale=0.5 ! should be chosen for specific function
REAL :: func, fx, fx_min, fx_new, temp, tfactor, x, x_min, x_new
REAL, DIMENSION(2) :: rand ! random numbers

! starting point for search
x = 1.0; fx = func(x); fx_min = fx
PRINT *, 'Starting from x = ', x, ', f(x) = ', fx

! annealing schedule
PRINT *, 'initial (high) temperature (e.g., 10)?'
READ *, temp
PRINT *, 'annealing temperature reduction factor (e.g., 0.9)?'
READ *, tfactor
PRINT *, 'number of steps per block (equilibration, e.g., 1000)?'
READ *, nsteps

DO WHILE (temp > 1E-5) ! anneal cycle

  DO istep = 1, nsteps
    CALL RANDOM_NUMBER(rand) ! 2 random numbers
    x_new  = x + scale*SQRT(temp)*(rand(1) - 0.5) ! stochastic move
    fx_new = func(x_new) ! new object function value
    IF (EXP(-(fx_new - fx)/temp) > rand(2)) THEN ! success, save
      fx = fx_new
      x = x_new
    END IF
      if(temp == 10.)write(1,fmt=*)temp,x,fx 
    IF (fx < fx_min) THEN
      fx_min = fx
      x_min = x
      PRINT '(3ES13.5)', temp, x_min, fx_min
    END IF
  END DO

  temp = temp * tfactor ! decrease temperature
END DO

END PROGRAM anneal

! Function to minimize
REAL FUNCTION func(x)
IMPLICIT NONE
REAL :: x
func = (x + 0.2)*x + cos(14.5*x - 0.3) 
END FUNCTION
