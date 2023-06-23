module Variables
    implicit none
    
    ! Define directions as 8-bit integers
    integer(1), parameter :: RI     = 1  ! 00000001
    integer(1), parameter :: RD     = 2  ! 00000010
    integer(1), parameter :: LD     = 4  ! 00000100
    integer(1), parameter :: LE     = 8  ! 00001000
    integer(1), parameter :: LU     = 16 ! 00010000
    integer(1), parameter :: RU     = 32 ! 00100000
    integer(1), parameter :: FULL   = 63 ! 00111111

    ! define the collision rules
    integer(1), dimension(6), parameter :: dir = (/RI, RD, LD, LE, LU, RU/)
    integer(1), parameter :: NUM_RULES = 2**6
    integer(1), parameter :: NUM_DIRS = 6
    integer(1), dimension(0: NUM_RULES - 1) :: rules = 0

    real :: start_time, end_time
    real, parameter :: density = 0.2, SQRT3_OVER2 = sqrt(3.0) / 2.0
    real, dimension(6), parameter :: ux = (/1.0, 0.5, -0.5, -1.0, -0.5, 0.5/)
    real, dimension(6), parameter :: uy = (/0.0, -SQRT3_OVER2, -SQRT3_OVER2, 0.0, SQRT3_OVER2, SQRT3_OVER2/)
    real :: vx = 0. ! averaged x velocity for every site configuration
    real :: vy = 0. ! averaged y velocity for every site configuration

    character(len=20) :: filename

    contains

    subroutine InitRules()
        ! two particles cyclic rules
        rules(ior(LE, RI)) = ior(RU, LD) ! rule(9)  = 36
        rules(ior(RU, LD)) = ior(LU, RD) ! rule(36) = 18
        rules(ior(LU, RD)) = ior(LE, RI) ! rule(18) = 9

        ! three particles zero-momentum rules
        rules(ior(ior(LU, LD), RI)) = ior(ior(RU, LE), RD) ! rule(21) = 42
        rules(ior(ior(RU, LE), RD)) = ior(ior(LU, LD), RI) ! rule(42) = 21
    end subroutine InitRules

end module Variables

module printing
    use Variables
    implicit none

    contains

    subroutine stampatutto(Lattice, newLattice, Lx, Ly, avg_vel_site)
        integer :: i, j
        integer, intent(in) :: Lx, Ly
        integer(1), dimension(:,:), allocatable, intent(inout) :: Lattice, newLattice
        real, dimension(:,:,:), allocatable, intent(inout) :: avg_vel_site
        do i = 0, Lx - 1
            do j = 0, Ly - 1
                write(*, '(I3.3, A)', advance='no') Lattice(i, j), ' '
            end do
            write (*,*)
        end do
        write(*,*)
        do i = 0, Lx - 1
            do j = 0, Ly - 1
                write(*, '(I3.3, A)', advance='no') newLattice(i, j), ' '
            end do
            write (*,*)
        end do
        write(*,*)
        do i = 0, Lx - 1
            do j = 0, Ly - 1
                write(*,'(A,F0.2,A,F0.2,A)', advance='no') '(', avg_vel_site(i,j,1), ',', avg_vel_site(i,j,2), ')'
            end do
            write(*,*)
        end do
    end subroutine stampatutto

end module printing

module Dynamic
    use iso_fortran_env
    use Variables
    implicit none

    contains

    subroutine InitLattice1(Lattice, Lx, Ly)
        integer :: i, j
        real :: r
        integer, intent(in) :: Lx, Ly
        integer(1), dimension(:,:), allocatable, intent(inout) :: Lattice

        do i = 0, Lx - 1
            do j = 0, Ly - 1
                call random_number(r)
                if (r < density / 2.0) then
                    Lattice(i, j) = ior(ior(RU, LE), RD)
                else if (r > density) then
                    Lattice(i, j) = 0
                else 
                    Lattice(i, j) = ior(ior(LU, LD), RI)
                end if
            end do
        end do

    end subroutine InitLattice1

    subroutine InitLattice2(Lattice, Lx, Ly)
        integer :: i, j
        real :: r
        integer, intent(in) :: Lx, Ly
        integer(1), dimension(:,:), allocatable, intent(inout) :: Lattice

        do i = 1, Lx - 2
            do j = 0, Ly / 5
                Lattice(i, j) = RD
            end do
            do j = 4 * Ly / 5, Ly - 1
                Lattice(i, j) = LU
            end do
        end do

    end subroutine InitLattice2

    subroutine InitLattice3(Lattice, Lx, Ly)
        integer :: i, j
        real :: r
        integer, intent(in) :: Lx, Ly
        integer(1), dimension(:,:), allocatable, intent(inout) :: Lattice

        do i = 3 * Lx / 8, 5 * Lx / 8
            do j = 3* Ly / 8, 5 * Ly / 8
                Lattice(i, j) = FULL
            end do
        end do

    end subroutine InitLattice3

    subroutine Streaming(Lattice, newLattice, Lx, Ly)

        integer :: i, j
        integer, intent(in) :: Lx, Ly
        integer(1), dimension(:,:), allocatable, intent(inout) :: Lattice, newLattice

        ! Streaming step: the particles in the lattice are moved to their neighboring lattice sites 
        ! based on their current direction of motion. 
        ! Each particle carries a specific direction encoded by a bit pattern

        ! Because the even rows are horizontally displaced 
        ! one half a Lattice spacing from the odd rows, 
        ! we need to treat odd and even rows separately.

        !                       RI        RD         LD         LE         LU         RU
        ! Odd row: (i, j) -> (i, j+1), (i+1, j), (i+1, j-1), (i, j-1), (i-1, j-1), (i-1, j)
        ! Even row:(i, j) -> (i, j+1), (i+1, j+1), (i+1, j), (i, j-1), (i-1, j), (i-1, j+1)
        ! There is a smart way of doing it: we update two consecutive lattice sites at once
        ! we loop j in increments of 2 in order to decrease reads and writes of neighbours
        newLattice = 0
        do j = 0, Ly - 1
            do i = 0, Lx - 1
                if (mod(i,2) == 0) then ! even row
                    newLattice(i,mod(j+1,Ly)) = ior(newLattice(i,mod(j+1,Ly)), iand(lattice(i,j), RI))      ! RI direction
                    newLattice(mod(i+1,Lx),mod(j+1,Ly)) = ior(newLattice(mod(i+1,Lx),mod(j+1,Ly)), iand(lattice(i,j), RD))  ! RD direction
                    newLattice(mod(i+1,Lx),j) = ior(newLattice(mod(i+1,Lx),j), iand(lattice(i,j), LD))      ! LD direction
                    newLattice(i,mod(j-1+Ly,Ly)) = ior(newLattice(i,mod(j-1+Ly,Ly)), iand(lattice(i,j), LE))      ! LE direction
                    newLattice(mod(i-1+Lx,Lx),j) = ior(newLattice(mod(i-1+Lx,Lx),j), iand(lattice(i,j), LU))      ! LU direction
                    newLattice(mod(i-1+Lx,Lx),mod(j+1,Ly)) = ior(newLattice(mod(i-1+Lx,Lx),mod(j+1,Ly)), iand(lattice(i,j), RU))  ! RU direction
                else ! odd row
                    newLattice(i,mod(j+1,Ly)) = ior(newLattice(i,mod(j+1,Ly)), iand(lattice(i,j), RI))      ! RI direction
                    newLattice(mod(i+1,Lx),j) = ior(newLattice(mod(i+1,Lx),j), iand(lattice(i,j), RD))      ! RD direction
                    newLattice(mod(i+1,Lx),mod(j-1+Ly,Ly)) = ior(newLattice(mod(i+1,Lx),mod(j-1+Ly,Ly)), iand(lattice(i,j), LD))  ! LD direction
                    newLattice(i,mod(j-1+Ly,Ly)) = ior(newLattice(i,mod(j-1+Ly,Ly)), iand(lattice(i,j), LE))      ! LE direction
                    newLattice(mod(i-1+Lx,Lx),mod(j-1+Ly,Ly)) = &
                                ior(newLattice(mod(i-1+Lx,Lx),mod(j-1+Ly,Ly)), iand(lattice(i,j), LU))  ! LU direction
                    newLattice(mod(i-1+Lx,Lx),j) = ior(newLattice(mod(i-1+Lx,Lx),j), iand(lattice(i,j), RU))      ! RU direction
                end if
            end do
        end do

        lattice = newLattice

    end subroutine Streaming

    subroutine Collision(Lattice, newLattice, Lx, Ly, avg_vel_site, avg_vel_cell, cell_size)

        integer :: i, j, k, l
        real :: block_sum_x, block_sum_y
        integer, intent(in) :: Lx, Ly, cell_size
        integer(1), dimension(:,:), allocatable, intent(inout) :: Lattice, newLattice
        real, dimension(:,:,:), allocatable, intent(inout) :: avg_vel_site, avg_vel_cell

        ! Collision step: handle interactions between particles.
        ! In this step, collisions between particles are resolved according to a predefined set of collision rules. 
        ! Each collision rule defines the new particle configurations based on the current configurations of neighboring particles. 
        ! The collision rules determine the outcome of particle interactions, including changes in particle directions and velocities.
        do i = 0, Lx - 1
            do j = 0, Ly - 1
                newLattice(i, j) = rules(Lattice(i, j))
                ! Compute avg velocity components at current lattice site
                vx = 0.; vy = 0.
                do k = 0, 5
                    if ( iand(Lattice(i,j), ishft(int(1, int8), k)) /= 0 ) then
                        vx = vx + ux(k+1)
                        vy = vy + uy(k+1)
                    end if
                end do
                avg_vel_site(i,j,1) = vx / 6.
                avg_vel_site(i,j,2) = vy / 6.
            end do
        end do

        do k = 0, Lx / cell_size - 1
            do l = 0, Ly / cell_size - 1
                block_sum_x = 0.0; block_sum_y = 0.0
                do i = k * cell_size, (k + 1) * cell_size - 1
                    do j = l * cell_size, (l + 1) * cell_size - 1
                        block_sum_x = block_sum_x + avg_vel_site(i, j, 1)
                        block_sum_y = block_sum_y + avg_vel_site(i, j, 2)
                    end do
                end do
                avg_vel_cell(k, l, 1) = block_sum_x / (cell_size * cell_size)
                avg_vel_cell(k, l, 2) = block_sum_y / (cell_size * cell_size)
            end do
        end do
        !newLattice = 0
    end subroutine Collision

    subroutine WriteVfield(Lx, Ly, avg_vel, it)
        integer :: i, j
        integer, intent(in) :: Lx, Ly, it
        real, dimension(:,:,:), allocatable, intent(inout) :: avg_vel

        write(filename, '(A, I0, A)') 'avg_vel_', it, '.dat'
        open(10, file=filename, status='replace')
        do j = 0, Ly - 1
            do i = 0, Lx - 1
                write(10, '(2I5, 2F10.6)') i, j, avg_vel(i, j, 1), avg_vel(i, j, 2)
            end do
            write(10, *)
        end do
        close(10)

    end subroutine WriteVfield

end module Dynamic

program FHP1
    use Variables
    use Dynamic
    implicit none

    integer :: Lx, Ly, it, nit, cell_size, frames
    integer(1), dimension(:,:), allocatable :: Lattice, newLattice
    real, dimension(:,:,:), allocatable :: avg_vel_site, avg_vel_cell

    character(len=32) :: arg
    call get_command_argument(1, arg); read(arg, *) Lx        ! Number of rows
    call get_command_argument(2, arg); read(arg, *) Ly        ! Number of columns
    call get_command_argument(3, arg); read(arg, *) nit       ! Number of iterations
    call get_command_argument(4, arg); read(arg, *) cell_size ! Size of the cell for averaging velocities
    call get_command_argument(5, arg); read(arg, *) frames    ! Every how many iterations save the vfield

    ! Check if the lattice dimensions are multiple of the cell size
    if ( mod(Lx, cell_size) /= 0 .or. mod(Ly, cell_size) /= 0) then
        print*, "Bad combination of lattice size and cell size"
        stop
    end if

    allocate(Lattice(0:Lx-1, 0:Ly-1), newLattice(0:Lx-1, 0:Ly-1))
    allocate(avg_vel_site(0 : Lx-1, 0 : Ly-1, 2), avg_vel_cell(0 : Lx/cell_size-1, 0 : Ly/cell_size-1, 2))
    Lattice = 0; newLattice = 0; avg_vel_site = 0.; avg_vel_cell = 0.
    
    call InitRules() ! Create the collision rules
    call InitLattice1(Lattice, Lx, Ly)

    ! Evolution of the system
    call cpu_time(start_time)
    do it = 1, nit
        call Streaming(lattice, newLattice, Lx, Ly)
        call Collision(Lattice, newLattice, Lx, Ly, avg_vel_site, avg_vel_cell, cell_size)
        if ( mod(it, frames) == 0 ) then
            !call WriteVfield(Lx, Ly, avg_vel_site, it / frames)
            call WriteVfield(Lx/cell_size, Ly/cell_size, avg_vel_cell, it / frames)
        end if
    end do
    call cpu_time(end_time)
    print*, "Evolve time:", end_time - start_time

    !call stampatutto(Lattice, newLattice, Lx, Ly, avg_vel_site)

end program FHP1
