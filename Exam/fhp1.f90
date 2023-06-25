module Variables
    implicit none
    
    ! Encode the possible directions using the first 6 bits of an 8 bit integer
    ! 8 bit integers have been use whenever possible
    integer(1), parameter :: RI     = 1  ! 00 000001
    integer(1), parameter :: RD     = 2  ! 00 000010
    integer(1), parameter :: LD     = 4  ! 00 000100
    integer(1), parameter :: LE     = 8  ! 00 001000
    integer(1), parameter :: LU     = 16 ! 00 010000
    integer(1), parameter :: RU     = 32 ! 00 100000
    integer(1), parameter :: B      = 64 ! 01 000000
    integer(1), parameter :: FULL   = 63 ! 00 111111

    ! Useful quantities
    integer(1), dimension(6), parameter :: dir = (/RI, RD, LD, LE, LU, RU/)
    integer(1), parameter :: NUM_RULES = 2**7 -1 ! Size of the lookup table (basically 64 bytes)
    integer(1), dimension(0: NUM_RULES) :: rules1 = 0, rules2 = 0

    double precision :: start_time, end_time ! For profiling
    real, parameter :: density = 0.6, SQRT3_OVER2 = sqrt(3.0) / 2.0
    real, dimension(6), parameter :: ux = (/1.0, 0.5, -0.5, -1.0, -0.5, 0.5/)
    real, dimension(6), parameter :: uy = (/0.0, -SQRT3_OVER2, -SQRT3_OVER2, 0.0, SQRT3_OVER2, SQRT3_OVER2/)
    real :: vx = 0. ! averaged x velocity for every site configuration
    real :: vy = 0. ! averaged y velocity for every site configuration

    character(len=20) :: filename

    contains

    subroutine InitRules() ! Two sets of rules

        integer(1) :: i
        integer(1) :: highBits, lowBits

        ! Set default rules: the particle keeps moving unperturbed if no collision
        do i = 0, NUM_RULES - 1
            rules1(i) = i; rules2(i) = i
        end do

        ! two particles cyclic rules
        rules1(ior(LE, RI)) = ior(RU, LD) ! rule(9)  = 36
        rules1(ior(RU, LD)) = ior(LU, RD) ! rule(36) = 18
        rules1(ior(LU, RD)) = ior(LE, RI) ! rule(18) = 9

        rules2(ior(LE, RI)) = ior(LU, RD) ! rule(9)  = 18
        rules2(ior(LU, RD)) = ior(RU, LD) ! rule(18) = 36
        rules2(ior(RU, LD)) = ior(LE, RI) ! rule(36) = 9

        ! three particles zero-momentum rules
        rules1(ior(ior(LU, LD), RI)) = ior(ior(RU, LE), RD) ! rule(21) = 42
        rules1(ior(ior(RU, LE), RD)) = ior(ior(LU, LD), RI) ! rule(42) = 21

        rules2(ior(ior(LU, LD), RI)) = ior(ior(RU, LE), RD) ! rule(21) = 42
        rules2(ior(ior(RU, LE), RD)) = ior(ior(LU, LD), RI) ! rule(42) = 21

        ! Add rules to bounce back at barriers
        rules1(ior(B,RI)) = ior(B,LE); rules2(ior(B,RI)) = ior(B,LE)
        rules1(ior(B,RD)) = ior(B,RU); rules2(ior(B,RD)) = ior(B,RU)
        rules1(ior(B,LD)) = ior(B,LU); rules2(ior(B,LD)) = ior(B,LU)
        rules1(ior(B,LE)) = ior(B,RI); rules2(ior(B,LE)) = ior(B,RI)
        rules1(ior(B,LU)) = ior(B,LD); rules2(ior(B,LU)) = ior(B,LD)
        rules1(ior(B,RU)) = ior(B,RD); rules2(ior(B,RU)) = ior(B,RD)

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
    use iso_fortran_env ! Needed to perform type casting to int8 at line 195
    use Variables
    implicit none

    contains

    subroutine InitLattice(Lattice, Lx, Ly)
        integer :: i, j
        real :: r
        integer, intent(in) :: Lx, Ly
        integer(1), dimension(:,:), allocatable, intent(inout) :: Lattice

        do i = 1, Lx - 2
            do j = 0, Ly / 2
                call random_number(r)
                if (r < density / 8.0) then
                    Lattice(i, j) = ior(ior(RU, LE), RD)
                else if (density / 8.0 <= r .and. r < density / 4.0) then
                    Lattice(i, j) = ior(ior(LU, LD), RI)
                else if (r > density) then
                    Lattice(i, j) = 0
                else
                    Lattice(i, j) = RI
                end if
            end do
        end do

        do j = 0, Ly - 1
            Lattice(0,j) = B; Lattice(Lx-1, j) = B
        end do

    end subroutine InitLattice

    subroutine InitLattice2(Lattice, Lx, Ly)
        integer :: i, j
        real :: r
        integer, intent(in) :: Lx, Ly
        integer(1), dimension(:,:), allocatable, intent(inout) :: Lattice

        do i = Lx / 2, Lx / 2 + 1
            do j = 0, Ly / 8
                Lattice(i,j) = RU
            end do
        end do

        do j = 0, Ly - 1
            Lattice(0,j) = B; Lattice(Lx-1, j) = B
        end do

    end subroutine InitLattice2

    subroutine Evolve(Lattice, newLattice, Lx, Ly, avg_vel_site, avg_vel_cell, cell_size, do_i_write)

        integer :: i, j, k, l
        real :: block_sum_x, block_sum_y, norm
        integer, intent(in) :: Lx, Ly, cell_size
        logical, intent(in) :: do_i_write
        integer(1), dimension(:,:), allocatable, intent(inout) :: Lattice, newLattice
        real, dimension(:,:,:), allocatable, intent(inout) :: avg_vel_site, avg_vel_cell

        ! Streaming step: the particles in the lattice are moved to their neighboring lattice sites 
        ! based on their current direction of motion. 
        ! Each particle carries a specific direction encoded by a bit pattern

        ! Because the even rows are horizontally displaced 
        ! one half a Lattice spacing from the odd rows, 
        ! we need to treat odd and even rows separately.

        ! These are the six nearest neighbours of site i,j for the two cases:
        !                       RI        RD         LD         LE         LU         RU
        ! Odd row: (i, j) -> (i, j+1), (i+1, j), (i+1, j-1), (i, j-1), (i-1, j-1), (i-1, j)
        ! Even row:(i, j) -> (i, j+1), (i+1, j+1), (i+1, j), (i, j-1), (i-1, j), (i-1, j+1)

        do i = 0, Lx - 1
            do j = 0, Ly - 1
                if (mod(i,2) == 0) then ! even row
                    ! RIGHT
                    newLattice(i,mod(j+1,Ly)) = ior(newLattice(i,mod(j+1,Ly)), iand(lattice(i,j), RI))
                    ! RIGHT DOWN
                    newLattice(mod(i+1,Lx),mod(j+1,Ly)) = ior(newLattice(mod(i+1,Lx),mod(j+1,Ly)), iand(lattice(i,j), RD))
                    ! LEFT DOWN
                    newLattice(mod(i+1,Lx),j) = ior(newLattice(mod(i+1,Lx),j), iand(lattice(i,j), LD))
                    ! LEFT
                    newLattice(i,mod(j-1+Ly,Ly)) = ior(newLattice(i,mod(j-1+Ly,Ly)), iand(lattice(i,j), LE))
                    ! LEFT UP
                    newLattice(mod(i-1+Lx,Lx),j) = ior(newLattice(mod(i-1+Lx,Lx),j), iand(lattice(i,j), LU))
                    ! RIGHT UP
                    newLattice(mod(i-1+Lx,Lx),mod(j+1,Ly)) = ior(newLattice(mod(i-1+Lx,Lx),mod(j+1,Ly)), iand(lattice(i,j), RU))

                else ! odd row (same but with different neighbour indices)
                    newLattice(i,mod(j+1,Ly)) = ior(newLattice(i,mod(j+1,Ly)), iand(lattice(i,j), RI))
                    newLattice(mod(i+1,Lx),j) = ior(newLattice(mod(i+1,Lx),j), iand(lattice(i,j), RD))
                    newLattice(mod(i+1,Lx),mod(j-1+Ly,Ly)) = ior(newLattice(mod(i+1,Lx),mod(j-1+Ly,Ly)), iand(lattice(i,j), LD))
                    newLattice(i,mod(j-1+Ly,Ly)) = ior(newLattice(i,mod(j-1+Ly,Ly)), iand(lattice(i,j), LE))
                    newLattice(mod(i-1+Lx,Lx),mod(j-1+Ly,Ly)) = &
                                ior(newLattice(mod(i-1+Lx,Lx),mod(j-1+Ly,Ly)), iand(lattice(i,j), LU))
                    newLattice(mod(i-1+Lx,Lx),j) = ior(newLattice(mod(i-1+Lx,Lx),j), iand(lattice(i,j), RU))
                end if
            end do
        end do

        ! Collision step: handle interactions between particles.
        ! In this step, collisions between particles are resolved according to the
        ! predefined set of collision rules. Each collision rule defines the new particle configurations 
        ! based on the configurations of the current lattice site.
        ! We switch rule set between even and odd iterations to prevent chirality (particles turning in circles)
        
        do i = 0, Lx - 1
            do j = 0, Ly - 1
                if ( mod(i * (Lx - 1) + j, 2) == 0  ) then
                    Lattice(i, j) = rules1(newLattice(i, j))
                else 
                    Lattice(i, j) = rules2(newLattice(i, j))
                end if
                newLattice(i, j) = 0

                ! Compute avg velocity components at current lattice site (only when we want to save)
                if ( do_i_write ) then
                    vx = 0.; vy = 0.
                    do k = 0, 5
                        ! If there is a particle going in a certain direction
                        if ( iand(Lattice(i, j), ishft(int(1, int8), k)) /= 0 ) then
                            ! Then add its velocity components to the total of the site
                            vx = vx + ux(k+1)
                            vy = vy + uy(k+1)
                        end if
                    end do
                    avg_vel_site(i,j,1) = vx / 6.
                    avg_vel_site(i,j,2) = vy / 6.
                end if
            end do
        end do

        do j = 0, Ly - 1 ! I'm sure there is a better way of doing it
            newLattice(0,j) = B; newLattice(Lx-1, j) = B
        end do

        ! Now compute the 2D block average to reduce noise; this is what will be plotted as vector field
        norm = 1. / real(cell_size * cell_size) ! computed just once here instead of inside the nested loop

        if ( do_i_write ) then
            do k = 0, Lx / cell_size - 1
                do l = 0, Ly / cell_size - 1
                    block_sum_x = 0.0; block_sum_y = 0.0

                    do i = k * cell_size, (k + 1) * cell_size - 1
                        do j = l * cell_size, (l + 1) * cell_size - 1
                            block_sum_x = block_sum_x + avg_vel_site(i, j, 1)
                            block_sum_y = block_sum_y + avg_vel_site(i, j, 2)
                        end do
                    end do
                    avg_vel_cell(k, l, 1) = block_sum_x * norm
                    avg_vel_cell(k, l, 2) = block_sum_y * norm
                end do
            end do
        end if
    end subroutine Evolve

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
    use printing
    use Dynamic
    implicit none

    integer :: Lx, Ly, it, nit, cell_size, frames
    integer(1), dimension(:,:), allocatable :: Lattice, newLattice
    real, dimension(:,:,:), allocatable :: avg_vel_site, avg_vel_cell
    integer :: i
    logical :: do_i_write

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
    call InitLattice(Lattice, Lx, Ly)

    ! Evolution of the system
    call cpu_time(start_time)
    do it = 1, nit
        do_i_write = mod(it, frames) == 0
        call Evolve(Lattice, newLattice, Lx, Ly, avg_vel_site, avg_vel_cell, cell_size, do_i_write)
        if ( do_i_write ) then ! Writing to file is slow!!!
            call WriteVfield(Lx/cell_size, Ly/cell_size, avg_vel_cell, it / frames)
        end if
    end do
    call cpu_time(end_time)
    print*, "Evolve time:", end_time - start_time

    do i = 0, NUM_RULES
        !print*, i, rules1(i), rules2(i)
    end do

    !call stampatutto(Lattice, newLattice, Lx, Ly, avg_vel_site) ! This was for debugging

end program FHP1
