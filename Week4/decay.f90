!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          Simulation of radioactive decay
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
PROGRAM decay
  IMPLICIT none
!  REAL, PARAMETER ::       lambda=0.2
REAL :: lambda
  REAL :: r
  INTEGER :: i, t, nleft, start, sizer
  integer, dimension(:), allocatable :: seed
  !
  call random_seed(sizer)
  allocate(seed(sizer))
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  read(*,*)seed
  call random_seed(put=seed)


  print*,' insert lambda>'
  read *, lambda
  !         initial values
  print *,"initial number of nuclei >"
  read *, start
  t = 1          ! initialize time
  nleft = start  ! at the beginning N(t=0)=start
  ! N(t) nuclei left at time t,
  ! that have a given probability lambda of decay
  ! in the time interval t:t+dt
  !            
  !         open output file
  OPEN(unit=7, FILE="decay.dat", status="replace",action="write")
  WRITE (unit=7,fmt=*) "# t ,       N(t)"
  WRITE (unit=7,fmt=*) "0  ", nleft !REAL(nleft)/start
  !
  !         Execution
  DO                               ! time loop
     DO  i = 1, nleft              ! loop on the nuclei left
        call random_number(r)
        IF (r <= lambda) THEN
           nleft = nleft - 1       ! update the number of nuclei left
        ENDIF
     END DO
     ! 
     WRITE (unit=7,fmt=*) t , nleft ! or REAL(nleft)/start
     if (nleft == 0) exit
     t = t + 1
  END DO
  !
  close(7)
  stop
END program decay
