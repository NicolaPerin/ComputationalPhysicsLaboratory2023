!
! simulation of 2D hard disks using MD - from Gould-Tobochnick
!
module periodic
public :: separation,pbc
integer, public, parameter  :: double = selected_real_kind(13)
contains

function separation(ds,L) result (separation_result)
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

public :: minimum_collision_time,move,reset_list,check_overlap
public :: uplist,downlist,initial,check_momentum,kinetic_energy
public :: check_collision,contact,save_config,output

integer, public :: N
real (kind = double), public :: Lx,Ly,t,timebig
real(kind = double),public,dimension (:),allocatable :: x,y,vx,vy
real(kind = double),public,dimension (:),allocatable :: collision_time, partner

contains

subroutine initial(vsum,rho,area)
   real (kind = double), intent (out) :: vsum,rho,area
   real :: a_x,a_y,vmax,rnd
   integer :: i,nx,ny,row,col
   character(len = 20) :: start,file_name
   character(len = 100) :: dum
   t = 0.0
   print *, "read file (f) or lattice start (l) = "
   read *, start
   if (start == "f" .or. start == "f") then
      print *, "file name = "
      read *, file_name
      open (unit=7,file=file_name,status="old",action="read")
      read (unit=7,fmt = *) N
      read (unit=7,fmt = *) Lx,Ly
      allocate(x(N))
      allocate(y(N))
      allocate(vx(N))
      allocate(vy(N))
      read (unit=7,fmt = *) dum
      do i = 1,N
         read (unit=7,fmt = *) x(i),y(i)
         print *, x(i),y(i)
      end do
      read (unit=7,fmt = *) dum
      do i = 1, N
         read (unit=7,fmt = *) vx(i),vy(i)
      end do
      close(unit=7)
   else if (start == "l" .or. start == "l") then
      print *, "N (a square...)= "
      read *, N     ! assume that sqrt(N) is an integer
      allocate(x(N))
      allocate(y(N))
      allocate(vx(N))
      allocate(vy(N))
      print *, "Lx = "
      read *, Lx
      print *, "Ly = "
      read *, Ly
      print *, "vmax = "
      read *, vmax
      nx = sqrt(real(N))
      ny = nx
       if (nx >= Lx .or. nx >= Ly) then
          print *, "box too small"
          stop
       end if
       a_x = Lx/nx             ! "lattice" spacing
       a_y = Ly/ny
       i = 0
       do col = 1,nx
           do row = 1,ny
               i = i + 1
               x(i) = (col - 0.5)*a_x
               y(i) = (row - 0.5)*a_y
               ! choose random positions and velocities
               call random_number(rnd)
               vx(i) = (2*rnd - 1)*vmax
               call random_number(rnd)
               vy(i) = (2*rnd - 1)*vmax
           end do
       end do
    end if
    call check_overlap()            ! check if two disks overlap
    call check_momentum()
    allocate(partner(N))
    allocate(collision_time(N))
    area = Lx*Ly
    rho = N/area
    timebig = 1.0e10
    vsum = 0                  ! virial sum
    do i = 1,N
        partner(i) = N
    end do
    collision_time(N) = timebig
    ! set up initial collision lists
    do i = 1,N
        call uplist(i)
    end do
end subroutine initial

subroutine check_momentum()
   real (kind = double) :: vxsum, vysum,vxcm,vycm
   !vxsum = 0.0
   !vysum = 0.0
   ! compute total center of mass velocity (momentum)
   vxsum = sum(vx)
   vysum = sum(vy)
   vxcm = vxsum/N
   vycm = vysum/N
   vx = vx - vxcm
   vy = vy - vycm
end subroutine check_momentum

subroutine minimum_collision_time(i,j,tij)

   real (kind = double), intent (out) :: tij
   integer, intent (out):: i,j
   
   integer :: k
   ! locate minimum collision time
   tij = timebig
   do k = 1, N
      if (collision_time(k) < tij) then
         tij = collision_time(k)
         i = k
      end if
   end do
   j = partner(i)
end subroutine minimum_collision_time

subroutine move(tij)
   real (kind = double), intent (in) :: tij
   integer :: k
   do k = 1,N
      collision_time(k) = collision_time(k) - tij
      x(k) = x(k) + vx(k)*tij
      y(k) = y(k) + vy(k)*tij
      x(k) = pbc(x(k),Lx)
      y(k) = pbc(y(k),Ly)
   end do
end subroutine move

subroutine reset_list(i,j)
   integer, intent (in) :: i,j
   integer :: k,test
    ! reset collision list for relevant particles
    do k = 1,N
        test = partner(k)
        if (k == i .or. test == i .or. k == j .or. test == j) then
           call uplist(k)
        end if
    end do
    call downlist(i)
    call downlist(j)
end subroutine reset_list

subroutine check_overlap()
   real (kind = double) :: tol,dx,dy,r2,r
   integer :: i,j
   tol = 1.0e-4
   do i = 1,N - 1
      do j = i + 1,N
         dx = separation(x(i) - x(j),Lx)
         dy = separation(y(i) - y(j),Ly)
         r2 = dx*dx + dy*dy
         if (r2 < 1.0) then
            r = sqrt(r2)
            if (1.0 - r > tol) then
               print *, "particles ",i," and ",j,"overlap"
               stop
            end if
         end if
      end do
   end do
end subroutine check_overlap

subroutine kinetic_energy(ke)
   real (kind = double), intent (out) :: ke
   integer :: i
   ke = 0.0
   do i = 1,N
      ke = ke + vx(i)*vx(i) + vy(i)*vy(i)
   end do
   ke = 0.5*ke
   print "(a,f13.6)", "kinetic energy =", ke/N
end subroutine kinetic_energy

subroutine uplist(i)
   integer, intent (in) :: i
   integer :: j
   ! look for collisions with particles j > i
   if (i == N) then
      return
   end if
   collision_time(i) = timebig
   do j = i + 1,N
      call check_collision(i,j)
   end do
end subroutine uplist

subroutine downlist(j)
   integer, intent (in) :: j
   integer :: i
   ! look for collisions with particles i < j
   if (j == 1) then
      return
   end if
   do i = 1,j - 1
      call check_collision(i,j)
   end do
end subroutine downlist

subroutine check_collision(i,j)
   integer, intent (in) :: i,j
   real (kind = double) :: dx,dy,dvx,dvy,bij,tij,r2,v2,discr
   integer :: xcell,ycell
   ! consider collisions between i and periodic images of j
   do xcell = -1,1
      do ycell = -1,1
         dx = x(i) - x(j) + xcell*Lx
         dy = y(i) - y(j) + ycell*Ly
         dvx = vx(i) - vx(j)
         dvy = vy(i) - vy(j)
         bij = dx*dvx + dy*dvy
            if (bij < 0) then
               r2 = dx*dx + dy*dy
               v2 = dvx*dvx + dvy*dvy
               discr = bij*bij - v2*(r2 - 1.0)
               if (discr > 0.0) then
                  tij = (-bij - sqrt(discr))/v2
                  if (tij < collision_time(i)) then
                     collision_time(i) = tij
                     partner(i) = j
               end if
            end if
         end if
      end do
   end do
end subroutine check_collision

subroutine contact(i,j,virial)
   real (kind = double), intent (out) :: virial
   integer, intent (in) :: i,j
   real (kind = double) :: dx,dy,dvx,dvy,factor,delvx,delvy
   ! compute collision dynamics for particles i and j at contact
   dx = separation(x(i) - x(j),Lx)
   dy = separation(y(i) - y(j),Ly)
   dvx = vx(i) - vx(j)
   dvy = vy(i) - vy(j)
   factor = dx*dvx + dy*dvy
   delvx = - factor*dx
   delvy = - factor*dy
   vx(i) = vx(i) + delvx
   vx(j) = vx(j) - delvx
   vy(i) = vy(i) + delvy
   vy(j) = vy(j) - delvy
   virial = delvx*dx + delvy*dy
end subroutine contact

subroutine output(collisions,temperature,vsum,rho,area)
   real (kind = double), intent (in) :: temperature,vsum,rho,area
   integer, intent (in) :: collisions
   real :: mean_virial,mean_pressure
   print "(a,i6)", "collisions =",collisions
   print "(a,f10.4)", "t =",t
   print *, ""
   mean_virial = vsum/(2.0*t)
   mean_pressure = rho*temperature + mean_virial/area
   print "(a,f10.4)", "P =", mean_pressure
end subroutine output

subroutine save_config(i,j,tij)
   real (kind = double), intent (inout) :: tij
   integer, intent (inout) :: i,j
   character(len = 32) :: config
   integer :: k
   ! move particles away from collision for final configuration
   call minimum_collision_time(i,j,tij)
   call move(0.5*tij)
   print *, "file name of configuration?"
   read *, config
   print *, config
   open (unit=1,file=config,status="replace",action="write")
   write (unit=1, fmt="(i4)") N
   write (unit=1, fmt="(2f13.6)") Lx,Ly
   write (unit=1,fmt="(t3,a,t20,a)") "x","y"
   do k = 1,N
      write (unit=1, fmt="(2f13.6)") x(k),y(k)
   end do
   write(unit=1, fmt="(t3,a,t20,a)") "vx","vy"
   do k = 1,N
      write(unit=1,fmt="(2f13.6)") vx(k),vy(k)
   end do
   close(unit=1)
end subroutine save_config

end module common

program hd
   ! dynamics of system of hard disks
   ! program based in part on fortran program of Allen and Tildesley
   use periodic
   use common

   real (kind = double) :: virial,temperature,vsum,rho,area,ke,tij
   integer :: collisions, i,j,nshow
   call initial(vsum,rho,area)
   call kinetic_energy(ke)
   temperature = ke/N
   collisions = 0                ! number of collisions
   nshow = 10       ! show output every nshow collisions
   do
      if (collisions > 100) then
         exit
      end if
      call minimum_collision_time(i,j,tij)
      ! move particles forward and reduce collision times by tij
      call move(tij)
      t = t + tij
      collisions = collisions + 1
      call contact(i,j,virial)       ! compute collision dynamics
      vsum = vsum + virial
      if (modulo(collisions,nshow) == 0) then
         call output(collisions,temperature,vsum,rho,area)
      end if
      ! reset collision list for relevant particles
      call reset_list(i,j)
      ! check for overlaps to debug program
      ! call check_overlap()
   end do
   call save_config(i,j,tij)
end program hd

