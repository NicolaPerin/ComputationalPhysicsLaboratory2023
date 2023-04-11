module functions
    use iso_fortran_env, only: int64
    implicit none
    integer, parameter :: dp = selected_real_kind(13)
    real(kind=dp), parameter :: pi = 3.14159265358979323846_dp
    contains
    function stirl(z) result(res)
        integer, intent(in) :: z
        real(kind=dp) :: res
        res = sqrt(2.0_dp * pi * real(z,kind=dp)) * (real(z,kind=dp)/exp(1.0_dp))**real(z,kind=dp)
        if (z == 0) res = 1.0_dp
    end function stirl

    recursive function fact(k) result(res)
    integer, intent(in) :: k
    integer(kind=int64) :: k64 = 0
    integer(kind=int64) :: res
    k64 = transfer(k, k64) ! cast to 64 bits integer
    if (k <= 0) then
        res = 1_int64
    else
        res = k * fact(k-1)
    end if
    end function fact
end module functions

module write
    use functions
    use iso_fortran_env, only: int64
    implicit none
    contains
        ! Write the empirical probability distribution P_N (point f of exercise) (O(N))
        subroutine WriteProbabilities(N, nruns, sum_P_N)
            integer, intent(in) :: N, nruns
            integer(kind=int64), dimension(:), allocatable, intent(in) :: sum_P_N ! MUST put allocatable
            integer :: i

            open(unit=13, file='P_N.txt', status='replace', action='write') ! point g of the exercise
            if (N < 21) then ! 20! is the max that 64 bit integers allow
                do i = -N, N
                    write(13,*) i, real(sum_P_N(i)) / nruns, 0.5_dp**N * fact(N) / (fact((N+i)/2) * fact((N-i)/2))
                end do
            else if (N < 100) then ! Use the Stirling approximation
                do i = -N, N
                    write(13,*) i, real(sum_P_N(i))/nruns, & ! empirical distribution
                                0.5_dp**N*stirl(N)/(stirl((N+i)/2)*stirl((N-i)/2))! Using Stirling
                end do
            else ! Big N
                do i = -N, N ! When benchmarking we don't care about the theoretical distribution
                    write(13,*) i, real(sum_P_N(i)) / nruns
                end do
            end if
            close(13)
        end subroutine WriteProbabilities

        subroutine WriteAverageWalkers(N, nruns, sum_x_i, sum_x2_i)
            integer, intent(in) :: N, nruns
            integer(kind=int64), dimension(:), allocatable, intent(in) :: sum_x_i, sum_x2_i
            integer :: i
            real(kind=dp) :: avg, delta

            open(unit=10, file='avg_walkers.txt', status='replace', action='write')
            do i = 1, N ! O(N)
                avg = real(sum_x_i(i),kind=dp) / nruns
                delta = real(sum_x2_i(i),kind=dp) / nruns - avg**2
                write(10,*) i, avg, delta
            end do
            close(10)
        end subroutine WriteAverageWalkers

        subroutine WriteDelta(N, nruns, sum_x_i, sum_x2_i)
            integer, intent(in) :: N, nruns
            integer(kind=int64), dimension(:), allocatable, intent(in) :: sum_x_i, sum_x2_i
            real(kind=dp) :: delta

            ! Calculate the accuracy of the mean square displacement (points c and e of the exercise) (O(1))
            delta = real(sum_x2_i(N),kind=dp) / nruns - (real(sum_x_i(N),kind=dp) / nruns)**2
            print*, "Delta:", abs(delta / N - 1.0_dp) ! N is the theoretical value (4*0.5*0.5*N*1)
            open(unit=12, file='delta.txt', position='append', action='write')
            write(12,*) N, nruns, abs(delta / N - 1.0_dp)
            close(12)
        end subroutine WriteDelta

        subroutine WriteTimings(size, N, nruns, comp, comm, writ)
            integer, intent(in) :: size, N, nruns
            real, intent(in) :: comp, comm, writ

            ! Write the timings in a file (O(1))
            open(unit=11, file='scalability.txt', position='append', action='write')
            write(11,*) size, N, nruns, comp, comm, writ
            close(11)
            print*, "Computation time:", comp, "s"
            print*, "Communication time:", comm, "s"
            print*, "Writing time:", writ, "s"
        end subroutine WriteTimings
end module write
