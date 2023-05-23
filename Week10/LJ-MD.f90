!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  LJ-MD.f90 (from Gould-Tobochnick)
!c
!c  simulate using Molecular Dynamics a 2D system of particles 
!c  interacting via the Lennard-Jones potential;
!c  Use periodic boundary conditions and minimum imagee convention;
!c  Use velocity-Verlet algorithm to integrate the equation of motion
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module periodic
public :: separation,pbc
integer, public, parameter :: double = 8
contains

function separation(ds,L) result (separation_result)
! minimum image convention.
! NOTE that this implementation assumes that the max dist between
! particles is L in each direction, i.e. it is OK if called after  
! the function pbd
   real (kind = double), intent (in) :: ds,L
   real (kind = double) :: separation_result
    if (ds > 0.5*L) then
        separation_result = ds - L
    else if (ds < -0.5*L) then
        separation_result = ds + L
    else
        separation_result = ds
    end if
end function separation

function pbc(pos,L) result (f_pbc)
! use PBC, fold a coordinate from [-L,2L] in [0,L]
   real (kind = double), intent (in) :: pos,L
   real (kind = double) :: f_pbc
    if (pos < 0.0) then
        f_pbc = pos + L
    else if (pos > L) then
        f_pbc = pos - L
    else
        f_pbc = pos
    end if
end function pbc

end module periodic

module common

use periodic
private

public :: initial,allocate_arrays,Verlet,accel,force
public :: check_momentum,save_config,output

integer, public :: N
real (kind = double), public :: Lx,Ly,t,dt,dt2
real (kind = double), public, dimension (:), allocatable :: x,y,vx,vy,ax,ay
integer, dimension(2), public :: seed

contains

subroutine initial(ke,kecum,pecum,vcum,area)
   real (kind = double), intent (out) :: ke,kecum,pecum,vcum,area
   character(len = 20) :: start,file_name
   character(len = 100) :: dum
   integer :: n1,i,row,col
   real :: a_x,a_y,vmax,rnd
   dt = 0.005
   dt2 = dt*dt
!   seed(1) = 1239
!   seed(2) = 1111
!   call random_seed(put=seed)
   print *, "read data (d), read file (f), or lattice start (l) ="
   read *, start
   if (start == "D" .or. start == "d") then
      N = 16  ! for the sake of simplicity, a perfect square
      call allocate_arrays()
      Lx = 6.0
      Ly = 6.0
      x(1:8) = (/ 1.09,3.12,0.08,0.54,2.52,3.03,4.25,0.89 /)
      x(9:16) = (/ 2.76,3.14,0.23,1.91,4.77,5.10,4.97,3.90 /)
      y(1:8) = (/ 0.98,5.25,2.38,4.08,4.39,2.94,3.01,3.11 /)
      y(9:16) = (/ 0.31,1.91,5.71,2.46,0.96,4.63,5.88,0.20 /)
      vx(1:8) = (/ -0.33,0.12,-0.08,-1.94,0.75,1.70,0.84,-1.04 /)
      vx(9:16) = (/ 1.64,0.38,-1.58,-1.55,-0.23,-0.31,1.18,0.46 /)
      vy(1:8) = (/ -0.78,-1.19,-0.10,-0.56,0.34,-1.08,0.47, 0.06 /)
      vy(9:16) = (/ 1.36,-1.24,0.55,-0.16,-0.83,0.65,1.48,-0.51 /)
   else if (start == "l" .or. start == "L") then
      print *, "N (a square...)= "
      read *, N     ! assume that sqr(N) is an integer
      call allocate_arrays()
      print *, "Lx = "
      read *, Lx
      Ly = 0.5*sqrt(3.0)*Lx
      n1 = sqrt(real(N))
      a_x = Lx/n1       ! lattice spacing
      a_y = 0.5*sqrt(3.0)*a_x
      vmax = 1.0
      i = 0
      ! triangular lattice
      do row = 1,n1
         do col = 1,n1
            i = i + 1
            x(i) = (col + 0.5*modulo(row,2) - 1)*a_x
            y(i) = (row - 0.5)*a_y
            ! choose random velocities
            call random_number(rnd)
            vx(i) = (2*rnd - 1)*vmax
            call random_number(rnd)
            vy(i) = (2*rnd - 1)*vmax
         end do
      end do
      do i = 1,N
         x(i) = pbc(x(i),Lx)
         y(i) = pbc(y(i),Ly)
      end do
   else if (start == "f" .or. start == "f") then
      print *, "file name = "
      read *, file_name
      open (unit=7,file=file_name,status="old",action="read")
      read (unit=7,fmt = *) N
      read (unit=7,fmt = *) Lx,Ly
      call allocate_arrays()
      read (unit=7,fmt = *) dum
      do i = 1,N
         read (unit=7,fmt = *) x(i),y(i)
      end do
      read (unit=7,fmt = *) dum
      do i = 1, N
         read (unit=7,fmt = *) vx(i),vy(i)
      end do
      close(unit=7)
   end if
   call check_momentum()
   ke = 0.0                    ! kinetic energy
   do i = 1,N
      ke = ke + vx(i)*vx(i) + vy(i)*vy(i)
   end do
   ke = 0.5*ke
   ! initialize sums
   kecum = 0.0
   pecum = 0.0
   vcum = 0.0
   area = Lx*Ly
   ! print heading for data
   write(unit=8,fmt="(t6,a,t17,a,t27,a,t37,a)")"time","E","T","P"
end subroutine initial

subroutine allocate_arrays()
   allocate(x(N))
   allocate(y(N))
   allocate(vx(N))
   allocate(vy(N))
   allocate(ax(N))
   allocate(ay(N))
end subroutine allocate_arrays

subroutine Verlet(ke,pe,virial)
   real (kind = double), intent (out) :: ke
   real (kind = double), intent (inout) :: pe,virial
   integer :: i
   real (kind = double) :: xnew,ynew

   do i = 1, N
      xnew = x(i) + vx(i)*dt + 0.5*ax(i)*dt2
      ynew = y(i) + vy(i)*dt + 0.5*ay(i)*dt2
      x(i) = pbc(xnew,Lx)
      y(i) = pbc(ynew,Ly)
      ! partially update velocity using old acceleration
      vx(i) = vx(i) + 0.5*ax(i)*dt
      vy(i) = vy(i) + 0.5*ay(i)*dt
   end do
   call accel(pe,virial)         ! new acceleration
   ke = 0.0
   do i = 1, N
      ! complete the update of the velocity using new acceleration
       vx(i) = vx(i) + 0.5*ax(i)*dt
       vy(i) = vy(i) + 0.5*ay(i)*dt
       ke = ke + vx(i)*vx(i) + vy(i)*vy(i)
   end do
   ke = 0.5*ke
   t = t + dt
end subroutine Verlet

subroutine accel(pe,virial)
  real (kind = double), intent (inout) :: pe,virial
  real (kind = double) :: dx,dy,fxij,fyij,pot
  integer :: i,j
  do i = 1, N
     ax(i) = 0.0
     ay(i) = 0.0
  end do
  pe = 0.0
  virial = 0.0
  do i = 1, N - 1            ! compute total force on particle i
     do j = i + 1,N        ! due to particles j > i
        dx = separation(x(i) - x(j),Lx)
        dy = separation(y(i) - y(j),Ly)
        ! acceleration = force because mass = 1 in reduced units
        call force(dx,dy,fxij,fyij,pot)
        ax(i) = ax(i) + fxij
        ay(i) = ay(i) + fyij
        ax(j) = ax(j) - fxij   ! Newton's third law
        ay(j) = ay(j) - fyij
        pe = pe + pot
        virial = virial + dx*ax(i) + dy*ay(i)
     end do
  end do
end subroutine accel

subroutine force(dx,dy,fx,fy,pot)
! Lennard-Jones potential
   real (kind = double), intent (in) :: dx,dy
   real (kind = double), intent (out) :: fx,fy,pot

   real (kind = double) :: r2,rm2,rm6,f_over_r

   r2 = dx*dx + dy*dy
   rm2 = 1.0/r2
   rm6 = rm2*rm2*rm2
   f_over_r = 24*rm6*(2*rm6 - 1)*rm2
   fx = f_over_r*dx
   fy = f_over_r*dy
   pot = 4.0*(rm6*rm6 - rm6)
end subroutine force

subroutine check_momentum()
   real (kind = double) :: vxsum, vysum,vxcm,vycm
   ! compute total center of mass velocity (momentum)
   vxsum = sum(vx)
   vysum = sum(vy)
   vxcm = vxsum/N
   vycm = vysum/N
   vx = vx - vxcm
   vy = vy - vycm
end subroutine check_momentum

subroutine save_config()
   character(len = 32) :: config
   integer :: i
   print *, "file name of configuration?"
   read *, config
   print *, config
   open (unit=1,file=config,status="replace",action="write")
   write (unit=1, fmt="(i4)")N
   write (unit=1, fmt="(2f13.6)")Lx,Ly
   write (unit=1,fmt="(t3,a,t20,a)")"x","y"
   do i=1,N
      write (unit=1, fmt="(2f13.6)") x(i),y(i)
   end do
   write(unit=1, fmt="(t3,a,t20,a)") "vx","vy"
   do i = 1,N
      write(unit=1,fmt="(2f13.6)")vx(i),vy(i)
   end do
   close(unit=1)
end subroutine save_config

subroutine output(ke,pe,virial,kecum,vcum,ncum,area)
   integer, intent (inout) :: ncum
   real (kind = double), intent(in) :: ke,pe,virial,area
   real (kind = double), intent(inout) :: kecum,vcum
   real (kind = double) :: E,mean_ke,P

   ncum = ncum + 1
   E = ke + pe               ! total energy
   kecum = kecum + ke
   vcum = vcum + virial
   mean_ke = kecum/ncum      ! still need to divide by N
   P = mean_ke + (0.5*vcum)/ncum  ! mean pressure * area
   P = P/area
   ! mean_ke/N    ! mean kinetic temperature
   write(unit=8,fmt="(4f10.4)") t,E,mean_ke/N,P
end subroutine output

end module common

program md
!  program adapted from Gould & Tobochnik, Chapter 8
!  simple N(N-1)/2 calculation of forces
   use periodic
   use common
   real (kind = double) :: E
   real (kind = double) :: ke,kecum,pecum,vcum,area,pe,virial
   integer :: ncum
   open(unit=8,file='energy',status='replace')
   call initial(ke,kecum,pecum,vcum,area)
   call accel(pe,virial)
   E = ke + pe                   ! total energy
   ncum = 0                      ! number of times data accumulated
   do
      if (t > 2.0) exit
      call Verlet(ke,pe,virial)
      call output(ke,pe,virial,kecum,vcum,ncum,area)
   end do
   close(8)
   call save_config()
end program md
