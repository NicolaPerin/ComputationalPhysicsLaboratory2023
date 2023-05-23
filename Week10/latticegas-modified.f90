! with additional lines to calculate other quantities, 
! including the variance of the <Delta R^2(t)> 
program latticegas

  implicit none

  logical, allocatable :: lattice(:,:) 
  integer, allocatable :: x(:),y(:)
  double precision, allocatable :: dx(:),dy(:)

  integer  :: unit_number
  integer  :: Nsteps,Np,L
  integer  :: istep,isubstep  
  integer  :: dir,i,j,nfail,njumps, sizer
  integer, dimension(:), allocatable :: seed

  integer,  parameter  ::  MAXINT=1000000000  ! Variables  for  counting
  ! allowed    directions   
  integer :: free(4),nfree   
  integer :: dxtrial(4),dytrial(4) 
  integer :: xnew(4),ynew(4)

  real, dimension(2) :: rnd(2)
  real :: rnd1
  double  precision :: dxsum,dysum,dxsqsum,dysqsum,dx4thsum,dy4thsum,vardrsqave
  double  precision :: t,deltat,drsqave,D,a,help,Dave_in_t

  character(len=15)::filename

  ! Set average time  between jumps and jump length Units is  s and cm 
  ! although actually this is not needed for the simulation 
  deltat=1.d0 ! or 1d-9 for  1 ns 
  a=1.d0      ! or 2e-8 for  2 Ang

  print*," Nsteps>"
  read*, Nsteps  
  print*," Np>"
  read*, Np 
  print*," L>"
  read*, L
  call random_seed(sizer)
  allocate(seed(sizer))
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  read(*,*)seed
  call random_seed(put=seed) 
  print *,'Doing lattice gas walk to',Nsteps,'MC steps, initial seed',seed  
  print *,'using',Np,' particles  on a',L,'^2 square lattice'

  if (Np >= L*L) then
     print *,'Number of particles > number of sites !!!' 
     STOP 'Too small lattice' 
  endif

  allocate(lattice(0:L-1,0:L-1))
  allocate(x(Np),y(Np)) 
  allocate(dx(Np),dy(Np))

  ! Mark all positions as empty 
  do i=0,L-1
     do j=0,L-1
        lattice(i,j) = .false. 
     enddo
  enddo

  ! enumeration of directions: 1 left 2 right 3 up 4 down 
  dxtrial(1)=+1; dytrial(1)= 0;   
  dxtrial(2)=-1; dytrial(2)= 0;   
  dxtrial(3)= 0; dytrial(3)=+1; 
  dxtrial(4)= 0; dytrial(4)=-1;

  nfail=0; njumps=0; 
  ! Generate particles on lattice 
  do i=1,Np
     do ! Loop until empty position found
        !  To   be  on  safe  side,   check  that  upper   limit  not  reached
        call random_number(rnd)
        x(i)=int(rnd(1)*L);  if (x(i)>=L) x(i)=L-1;  
        y(i)=int(rnd(2)*L);  if (y(i)>=L) y(i)=L-1; 
        if (lattice(x(i),y(i))) then
           ! Position already filled, loop to find new trial 
           cycle 
        else
           lattice(x(i),y(i))=.true. 
           !  Success, go  to next particle  
           exit 
        endif
     enddo
     dx(i)=0.0d0; dy(i)=0.0d0; 
  enddo

  write(1,fmt=*)'# Density Np/L^2=',real(Np)/L**2,'; Np=',Np,'; L=',L
  write(1,fmt=*)'#  t, drsqave(t), sqrt(var(drsqave(t))/Np), D(t), <D>_T'

  t=0.0; Dave_in_t=0.d0
  do istep=0,Nsteps-1 ! Loop over MC steps
     do isubstep=1,Np ! Do all particles on average once every MC step
        ! Pick one particle at random 
        call random_number(rnd1)
        i=int(rnd1*Np)+1; if (i>Np) i=Np;

        ! Find possible directions, store it in free() 
        nfree=0 
        do j=1,4
           xnew(j)=x(i)+dxtrial(j); 
           if  (xnew(j) >= L) xnew(j)=0; if (xnew(j)<0) xnew(j)=L-1; 
           ynew(j)=y(i)+dytrial(j); 
           if  (ynew(j) >= L) ynew(j)=0; if (ynew(j)<0) ynew(j)=L-1; 
           if (.not. lattice(xnew(j),ynew(j))) then
              ! Success: position free 
              nfree=nfree+1 
              free(nfree)=j 
           endif
        enddo

        ! If no possible directions, get new particle 
        if (nfree == 0) then
           nfail=nfail+1 
           cycle 
        endif
        njumps=njumps+1

        !  Pick  one of  the  possible directions  randomly  
        !  Note that  the dir>nfree  check here  really is  needed!  
        call random_number(rnd1)
        dir=int(rnd1*nfree)+1; if (dir>nfree) dir=nfree 
        j=free(dir)

        ! Now x(i),y(i) is old position and xnew(j),ynew(j) new 
        ! Double check that new site really is free
        if (lattice(xnew(j),ynew(j))) then
           print *,'ERROR:    THIS   SHOULD   BE   IMPOSSIBLE'   
           print *,i,j,dir,nfree  
           print *,free  
           print *,x(i),y(i),xnew(j),ynew(j) 
           STOP 'ERROR  new  site  bug'  
        endif
        !Empty  old  position  and  fill  new
        lattice(x(i),y(i))=.false. 
        lattice(xnew(j),ynew(j))=.true.

        x(i)=xnew(j); y(i)=ynew(j);          
        dx(i)=dx(i)+dxtrial(j); dy(i)=dy(i)+dytrial(j);

!       if you want to follow each particle (but pay attention: you can fill your disk space!)
!        unit_number = 10+i
!        write(unit_number,*)dx(i)**2+dy(i)**2     ! see "part[particle_index]"
     enddo

  !   if (mod(istep*Np,1000000) == 0) then
! Calculate  and print intermediate results - uncomment/modify if you do not want to print every step
        ! Get total displacement from dx,dy 
        dxsum=0.0d0; dysum=0.0d0; 
        dxsqsum=0.0d0; dysqsum=0.0d0;
        dx4thsum=0.0d0;  dy4thsum=0.0d0;
        do i=1,Np
           dxsum=dxsum+dx(i);   dysum=dysum+dy(i);   
            dxsqsum=dxsqsum+dx(i)*dx(i);
            dysqsum=dysqsum+dy(i)*dy(i);
            dx4thsum=dx4thsum+dx(i)**4;
            dy4thsum=dy4thsum+dy(i)**4;
        enddo
        drsqave=(dxsqsum+dysqsum)/Np
        vardrsqave=(dx4thsum+dy4thsum)/Np - ((dxsqsum+dysqsum)/Np)**2

        if (t>0.0) then
           !  Get  diffusion coefficient  by  proper scaling  
           D=drsqave*a*a/(4*t)
            Dave_in_t=Dave_in_t+D
        write(1,fmt='(5(1pe10.2,2x))')t,drsqave,sqrt(vardrsqave/Np),D,Dave_in_t/t  ! see "aveNp.run_number"
        endif

 !    endif

     t=t+deltat  

  enddo

  !  Get   total  displacement   from  dx,dy
  dxsum=0.0d0; dysum=0.0d0; 
  dxsqsum=0.0d0; dysqsum=0.0d0; 
  do i=1,Np
     dxsum=dxsum+dx(i);   dysum=dysum+dy(i);   
     dxsqsum=dxsqsum+dx(i)*dx(i); dysqsum=dysqsum+dy(i)*dy(i);     
  enddo
  print *,'dxsum',dxsum,'   dysum',dysum 
  print *,'dxsqsum',dxsqsum,' dysqsum',dysqsum

  drsqave=(dxsqsum+dysqsum)/(1.0*Np)  
  print *,'drsqave',drsqave 
  print *,'Number of  failed jumps',nfail,' number of  successes',njumps 
  ! Get diffusion  coefficient  by  proper scaling  
  D=drsqave*a*a/(4*t) 
  write(*,fmt='(3(a,1pe10.2))')'Density Np/L^2=',real(Np)/L**2,' :  t=',t,';   D=',D

end  program latticegas
