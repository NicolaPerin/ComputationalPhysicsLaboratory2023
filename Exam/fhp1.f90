module utilities
    implicit none
    integer, parameter :: L1 = 200, L2 = 200
    integer, parameter :: N = L1 * L2
    integer, parameter :: EMPTY = 0
    integer, parameter :: FULL = 1
    real :: start_time, end_time

    contains

    function fhp1(s1, s2, s3, s4, s5, s6) result(res)
        integer :: s1, s2, s3, s4, s5, s6
        integer :: res
        integer :: sum
        logical :: tmp1, tmp2, tmp3

        sum = s1 + s2 + s3 + s4 + s5 + s6
        if (sum == 0 .or. sum == 6) then
            res = EMPTY
        else if (sum == 1 .or. sum == 5) then
            res = FULL
        else if (sum == 2) then
            if ((IAND(s1,s4) /= 0) .or. (IAND(s2,s5) /= 0) .or. (IAND(s3,s6) /= 0)) then
                res = EMPTY
            else
                res = FULL
            end if
        else if (sum == 3) then
            if ((IAND(IAND(s1,s2),s3) /= 0) .or. (IAND(IAND(s4,s5),s6) /= 0)) then
                res = FULL
            else
                res = EMPTY
            end if
        else if (sum == 4) then
            tmp1 = (IAND(IAND(IAND(s1,s2),s3),s4) /= 0)
            tmp2 = (IAND(IAND(IAND(s2,s3),s4),s5) /= 0)
            tmp3 = (IAND(IAND(IAND(s3,s4),s5),s6) /= 0)
            if (tmp1 .or. tmp2 .or. tmp3) then
                res = FULL
            else
                res = EMPTY
            end if
        end if

    end function fhp1

    subroutine update_lattice(Lat, newLat)
        integer :: Lat(N), newLat(N)
        integer :: i, j, index, U, D, LE, R, UL, UR

        do i = 1, L1
            do j = 1, L2
                index = (i - 1) * L2 + j 
                U = mod(i - 2 + L1, L1) * L2 + j 
                D = mod(i, L1) * L2 + j 
                LE = (i - 1) * L2 + mod(j - 2 + L2, L2) + 1 
                R = (i - 1) * L2 + mod(j, L2) + 1 
                UL = mod(i - 2 + L1, L1) * L2 + mod(j - 2 + L2, L2) + 1 
                UR = mod(i - 2 + L1, L1) * L2 + mod(j, L2) + 1 
    
                newLat(index) = fhp1(Lat(U), Lat(D), Lat(LE), Lat(R), Lat(UL), Lat(UR))
            end do 
        end do  

    end subroutine update_lattice

end module utilities

program lattice_gas_automata
    use utilities
    implicit none

    integer :: lattice(N) = 0, newlattice(N) = 0
    integer :: i, j, t, sum1, sum2
    real :: r
    logical :: printing

    print*, "Print?"
    read*, printing

    call cpu_time(start_time)
    call random_seed()
    do i = 1, N
        call random_number(r)
        if ( r > 0.5 ) then
            lattice(i) = 1  
        end if
    end do
    sum1 = sum(lattice)
    call cpu_time(end_time)
    print*, "Init time:", end_time - start_time

    call cpu_time(start_time)
    do t = 1, 1000
        if (mod(t, 2) == 1) then
            call update_lattice(lattice, newlattice)
        else
            call update_lattice(newlattice, lattice)
        end if
    end do
    sum2 = sum(lattice)
    call cpu_time(end_time)
    print*, "Evolve time:", end_time - start_time

    if ( printing .eqv. .true. ) then
        do i = 1, L2
            write(*, '(A)', advance='no') '-'
        end do
        write(*,*)

        do i = 1, L1
            do j = 1, L2
                if (lattice((i - 1) * L2 + j) == 1) then
                    write(*, '(i0)', advance='no') 0
                else
                    write(*, '(1x)', advance='no')
                end if
            end do
            write(*,*)
        end do
    
        do i = 1, L2
            write(*, '(A)', advance='no') '-'
        end do
        write(*,*)
    end if

    print*, "N = ", N
    print*, "Number of particles (before):", sum1
    print*, "Number of particles (after):", sum2

end program lattice_gas_automata
