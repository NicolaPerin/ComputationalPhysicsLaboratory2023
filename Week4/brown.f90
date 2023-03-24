PROGRAM Brown
  IMPLICIT NONE
  INTEGER                            		:: npart,it,nit,i,j
  REAL,DIMENSION(:,:),allocatable   	:: pos,pos0,vel,f
  REAL,DIMENSION(:), allocatable              	:: mass
  REAL,DIMENSION(2)                  		:: harvest ! array with 2 random numbers 
  REAL                               		:: dt,gamma,t,w,msq

  WRITE(*,*)"Insert the number of heavy particles  : "
  READ*,npart
  allocate(pos(2,npart))
  allocate(pos0(2,npart))
  allocate(vel(2,npart))
  allocate(f(2,npart))
  allocate(mass(npart))
  WRITE(*,*)"Insert mass of the heavy particles (in kg)  : "
  READ*,mass(1)
  mass(2:npart)=mass(1)
  WRITE(*,*)"Insert time step (in seconds) :"
  READ*,dt
  WRITE(*,*)"Insert number of iterations :"
  READ*,nit 
  WRITE(*,*)"Insert gamma and kT (in J) : "
  READ*,gamma,t

  vel = 0    ! Zero initial positions and velocities
  pos0 = 0
  it = 0       
  pos = pos0

  ! CALL interazione(pos,f)  ! in case of external force to be added
  f = 0                      ! here no external force: only drug and random forces
  WRITE(1,*)"# iteration, time, pos_x, pos_y, vel_x, vel_y of particle 1"
  WRITE(1,*)it,it*dt,pos(1,1),pos(2,1),vel(1,1),vel(2,1)
  WRITE(2,*)"# iteration, time, mean square displacement"

  DO it=1,nit
     DO j=1,npart
        DO i=1,2
           call gasdev(w)
           vel(i,j) = vel(i,j)*( 1 - gamma*dt/mass(j))  + dt  * f(i,j)/mass(j) &
                + w*sqrt(2*gamma*t*dt)/mass(j)
           pos(i,j) = pos(i,j) + vel(i,j) * dt 
        END DO

     END DO

     !CALL interazione(pos,f)
     msq = sum( (pos - pos0)**2 )/npart

     WRITE(1,*)it,it*dt,pos(1,1),pos(2,1),vel(1,1),vel(2,1) 
     WRITE(2,*)it,it*dt,msq

  END DO
  close(1)
  close(2)    
  stop

contains

  SUBROUTINE gasdev(rnd)
    IMPLICIT NONE
    REAL, INTENT(OUT) :: rnd
    REAL :: rsq,v1,v2
    REAL, SAVE :: g
    LOGICAL, SAVE :: gaus_stored=.false.
    if (gaus_stored) then
       rnd=g
       gaus_stored=.false.
    else
       do
          call random_number(v1)
          call random_number(v2)
          v1=2.*v1-1.
          v2=2.*v2-1.
          rsq=v1**2+v2**2
          if (rsq > 0. .and. rsq < 1.) exit
       end do
       rsq=sqrt(-2.*log(rsq)/rsq)
       rnd=v1*rsq
       g=v2*rsq
       gaus_stored=.true.
    end if
  END SUBROUTINE gasdev

  END PROGRAM Brown
