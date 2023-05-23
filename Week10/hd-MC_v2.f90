module common
implicit none
public :: initial,move,correl,output,save_config

integer, public :: N,nmcs,nbin,accept,nequil
real , dimension(:),allocatable, public :: x,y,gcum
real , public :: Lx,Ly,kx,ky, dr,dxmax,dymax, rho

contains

subroutine initial()
! (NOTE: not very good... to be improved)
      integer :: Nx,Ny,m,i,j
      real  :: rmax,dx,dy,temp
      character(len=1) :: ans
      logical :: expr
      
      print*,' number of particles',&
           '(choose a number that can be factorized in two integers) >'
      read*,N
      rho = float(N)/(Lx*Ly)
      write(*,*)' reduced density (rho): ',rho
      print*,' number of columns Ny  (for which  N can be divided) >'
      read*,Ny
      temp = mod(N,Ny)
      if(mod(N,Ny) .ne. 0)call error('error in initial',' mod(N,Ny)',temp)
      Nx = N / Ny
      allocate(x(N),y(N))

      print*,' number of  MC steps per particle per equilibration>'
      read*,nequil
      print*,' number of  MC steps per particle per data production>'
      read*,nmcs
      rmax=min(Lx/2.0,Ly/2.0)
      print*,' number of bins for g(r) >'
      read*, nbin 
	  allocate(gcum(nbin))
      dr=rmax/nbin           
      print*,' max trial displacement of a particle in x (dxmax) ; equal for y (dymax) >'
      read*,dxmax
      dymax = dxmax
      kx = 2./Lx   ! parameter for PBC
      ky = 2./Ly
      print*,'N =',N,', Ny=',Ny,', Nx=',Nx,'nmcs=',nmcs
      print*,'nbin=',nbin,', rmax=',rmax,', Lx=',Lx,', Ly=',Ly
!      write(*,'(3a,i4,1x,5a,i4,1x,5a,i4,1x,5a,i4)')'# N =',N,', Ny=',Ny,', Nx=',Nx,'nmcs=',nmcs
!      write(*,'(5a,1x,i4,1x,7a,f6.2,1x,a,5x,f6.2,1x,5a,f6.2)')'# nbin=',nbin,', rmax=',rmax,', Lx=',Lx,', Ly=',Ly
	  
      print*,' start from old config. (y/n)? >'
      read(*,'(a)') ans
      if (ans=='y') then
         call read_config()
      else                   ! CREATES a new triangular lattice
         dx = Lx / float(Nx)        
         dy = Ly / float(Ny)
!max. packing 
!dx and dy must be separately .ge. sqrt(3)/2 (in units of sigma)
!and the product dx*dy .ge. sqrt(3)/2
         expr = ((dx<(sqrt(3.)/2)) .or. (dy<(sqrt(3.)/2) ) .or. &
			((dx*dy)<(sqrt(3.)/2)))
		 if (expr) call error('error in initial','  dx*dy = ',dx*dy)
         m = 0
         do j = 1,Ny,2
            do i = 1,Nx
               m = m + 2
               if (m>N) call error('error in initial',' m>N! m = ',float(m))
               ! this should never occur! check!!!
               x(m-1) = (i-0.75)*dx
               y(m-1) = (j-0.5) *dy
               x(m)   = (i-0.25)*dx
               y(m)   = (j+0.5) *dy
           end do
         end do
      end if
      
	  return
end subroutine
      
subroutine move()
! pick up a particle randomly, generates trial cocordinates, 
! apply PBC, check overlap
      integer :: itrial,i
	  real  :: xtrial,ytrial,r,ran1,xold,yold
	  real , dimension(2) :: ran2
      do i=1,N
	call random_number(ran1)
	itrial = int(N*ran1 + 1)    ! choose a particle
	if(itrial>N) then
        print*,'error generating ran1'
        stop     
	end if
                 xold=x(itrial)
		 yold=y(itrial)
		 call random_number(ran2)
		 xtrial = x(itrial) + (2.*ran2(1)-1.)*dxmax  
                 ytrial = y(itrial) + (2.*ran2(2)-1.)*dymax  
         call cell(xtrial,ytrial)
         call overlap(itrial,xtrial,ytrial)
      end do
      return
end subroutine

subroutine cell(xtrial,ytrial)
! input  xtrial e ytrial,
! output : the same  but rescaled using PB (x in [0,Lx[ and y in [0,Ly[)
      real , intent(inout):: xtrial,ytrial
! floor(A) return the largest integer < or = A	        
! What follows is OK whatever xtrial and ytrial are
      xtrial = xtrial - Lx*floor(xtrial/Lx)
      ytrial = ytrial - Ly*floor(ytrial/Ly) 
      if(xtrial<0 .or. xtrial>=Lx .or. ytrial<0 .or. ytrial>=Ly) &
		call error('error   in  cell','  xtrial =',xtrial)
      return
end subroutine

subroutine separation(dx,dy)
! in input  distances dx and dy , which (being already rescaled)
!           are between 0 and L 
! in output the same, but rescaled using minimum image convention
!           and therefore = or < L/2)  
      real , intent(inout) :: dx,dy
      dx = dx - Lx*int(kx*dx)
      dy = dy - Ly*int(ky*dy)
      return
end subroutine

subroutine overlap(itrial,xtrial,ytrial)
! decide to accept or not after cheching overlap;
! in case, update arrays x,y
      
      integer, intent(in) :: itrial
      integer :: j
      real ,intent(in) :: xtrial,ytrial
      real :: r2,dx,dy
	  
! disks overlap if distance < 1 (diameter) 
	  do j = 1,N
         if (itrial /= j) then
            dx = x(j) - xtrial
            dy = y(j) - ytrial
            call separation(dx,dy)    ! calculate true distances
            r2 = dx*dx + dy*dy
      	    if (r2<1.0) return                   ! overlap!
         end if
      end do
         
	accept = accept + 1
        x(itrial) = xtrial
        y(itrial) = ytrial
      return
end subroutine

	subroutine correl()
	! array gcum contains the number of particles at distance
	! between r and r+dr from a given particle.
	! gcum is calculated for N/2 particles (to count once each pair)
      
      integer :: ibin,i,j
      real :: r2,dx,dy
      
      do i=1,N-1
         do j=i+1,N
            dx = x(i) - x(j)
            dy = y(i) - y(j)
            call separation(dx,dy)
            r2 = dx*dx + dy*dy
            ibin = int(sqrt(r2)/dr)+1
            if (ibin<=nbin) then	
		gcum(ibin) = gcum(ibin) + 1 
		!ngcum=ngcum +1
            end if
          end do
      end do
     return
	end subroutine

      subroutine output()
! calculate the normalized pair correlation function  g(r)
! dividing gcum for th enumber of samples, density, area 2\pi r dr
! of the circular strip from the referring particle. Divide by N/2
! because gcum has been calculated by summing over N/2 particles.
! Note: max. separation is max(Lx/2,Ly/2) (in PBC)  
      integer :: ir
      real :: rmax,xnorm,r,g,area,pi,rho_max

      pi = 2*asin(1.)
      rho_max=2/sqrt(3.)
	  
      rmax = nbin * dr 
      write(*,fmt='(a23,1x,f6.2,a1)')' rho/rho_max=',rho/rho_max
      write(*,*)' acceptance ratio:',real(accept)/(N*nmcs)
      open(unit=3,file='g_of_r.dat',status='unknown')
      write(3,*)'# number of particles (even) =',N
      write(3,*)'# dimensions Lx, Ly  of the box =',Lx,Ly
      write(3,*)'# reduced density (rho)=',rho
      write(3,*)'# data production MC steps per particle (nmcs)=',nmcs
      write(3,*)'# bin dr for g(r) =',dr
      write(3,*)'# acceptance ratio =',float(accept)/(N*nmcs)
      write(3,*)'# rmax =',rmax
! g(r) is calculated after each MC step
      xnorm = 2./(rho*nmcs*N)
      do ir = 1,nbin
            r    = ir*dr + 0.5*dr  ! r in the middle of the circular shell
            area = 2.0*pi*r*dr     ! area of the shell
            g    = gcum(ir)*xnorm/area
            write(3,*)r,g
      end do
      close(3)
         return
      endsubroutine


	subroutine save_config()
	!save config in a file
		integer :: i
		character(len=20) :: fileout
      
		print*,' output filename to save config in a file  >'
		read(*,'(a)')fileout
		open(unit=3,file=fileout,status='replace',action='write')
		write(3,*) Lx,Ly
		do i=1,N
			write(3,*)x(i),y(i)
		enddo   
		return
	endsubroutine


      subroutine store_config(fileout)
            !storing to disk each configuration in extended XYZ format https://github.com/libAtoms/extxyz
                  integer :: i
                  character(len=20), intent(in) :: fileout
                  character(len=30) :: LLx
                  write(LLx,*) Lx
                  LLx = 'Lattice="'//TRIM(LLx)
                  open(unit=3,file=fileout,action='write')
                  write(3,*)N
                  write(3,*)LLx,'0.0 ','0.0 ','0.0 ',Ly,'0.0 ', '0.0 ', '0.0 ', '100.0 " Properties=species:S:1:pos:R:3'
                  do i=1,N
                        write(3,*)'C',x(i),y(i),0
                  enddo   
                  return
            endsubroutine
	
	subroutine read_config()
		real :: Lxold,Lyold,xscale,yscale
		character(len=20) :: filein
		integer :: i
		
		print*,' filename old config. >'
         read(*,'(a)') filein
         open(unit=2,file=filein,status='unknown')
         read(2,*)Lxold,Lyold
         xscale = Lx / Lxold
         yscale = Ly / Lyold
         do i = 1,N
            read(2,*) x(i),y(i)
            x(i) = x(i) * xscale  ! uniformly rescale all x positions
            y(i) = y(i) * yscale  ! uniformly rescale all y positions
         enddo
         close(2)
	endsubroutine

	subroutine error(comment1,comment2,dummy)
		character(len =16) :: comment1
		character(len =10) :: comment2
		real :: dummy
    
		write(*,'(1x,a16,a10,1x,f10.4)')comment1,comment2,dummy
		stop
		return
	endsubroutine

endmodule common

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
 program hdMC
!MC simulation of hard disks in 2D with PBC
!unit length is the diameter, radius= 0.5
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
   use common
   implicit none
   integer :: imcs
   character(1) :: ans
   character(len=20) :: fileout

   Lx=4.41*sqrt(0.9503/0.3)
   Ly=0.5*sqrt(3.)*Lx
   print*,'default values for Lx,Ly=',Lx,Ly
   print*,' change? (y/n)'
   read(*,'(a)')ans
   if(ans=='y')then
      print*,' insert Lx >'
      read(*,*)Lx
      print*,' insert Ly >'
      read(*,*)Ly
   end if
   print*,' output filename to save config in a file  >'
   read(*,'(a)') fileout

   call initial()

   do imcs = 1, nequil   !equilibration
       call move()
   enddo
	
   gcum = 0.
   accept = 0
    do imcs = 1, nmcs
       call move()
       call correl()
       call store_config(fileout)
    enddo
	
    call output()
    !call save_config()
    deallocate(x,y)
    deallocate(gcum)
    stop
endprogram hdMC
      
