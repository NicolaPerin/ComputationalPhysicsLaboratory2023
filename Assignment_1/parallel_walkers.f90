program parallel_walkers
    use mpi ! who is this
    use functions
    use write
    use iso_fortran_env, only: int64
    implicit none
    integer :: i, j, N, nruns, clock, ierr, rank, size, sizer, rest = 0
    integer(kind=int64) xN
    integer, dimension(:), allocatable :: seed
    integer(kind=int64), dimension(:), allocatable :: all_seeds, P_N, sum_P_N, x_i, x2_i, sum_x_i, sum_x2_i
    real :: start_time, end_time, comm, comp, writ
    real(kind=dp), dimension(:), allocatable :: p, r

    character(len=32) :: arg ! We like to use command line arguments
    call get_command_argument(1, arg); read(arg, *) N
    call get_command_argument(2, arg); read(arg, *) nruns

    ! Initialize the MPI environment
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    ! If the number of runs is not a multiple of the number of processes, 
    ! the remainder is distributed across the first processes
    if ( rank < mod(nruns, size)) rest = 1

    ! This is a pure MPI implementation -> Distributed memory! (Shared memory is a nightmare)
    ! Allocate a bunch of stuff (each process does it on its own memory)
    call random_seed(sizer)
    allocate(seed(sizer), p(N), r(sizer*nruns), all_seeds(sizer*nruns))
    allocate(x_i(N), x2_i(N), sum_x_i(N), sum_x2_i(N))
    allocate(P_N(-N:N), sum_P_N(-N:N))

    x_i = 0; x2_i = 0; sum_x_i = 0; sum_x2_i = 0; P_N = 0; sum_P_N = 0 ! Init everybody to zero

    call system_clock(count=clock) ! Set the seed using the clock (O(1))
    seed = clock + 37 * (/ (i - 1, i = 1, sizer) /) ! The 37 is still a mistery to this day
    call random_seed(put=seed)

    call random_number(r) ! Generate a random array whose elements will be used as seed components (O(nruns))
    ! Apparently Fortran has list comprehensions but you can't do things like +=
    all_seeds = [(floor(sizer * nruns * r(i)), i = 1, sizer * nruns)]

    start_time = MPI_Wtime() ! Profiling

    ! The computation is carried out in parallel by many processes
    ! This is the heaviest part of the program, having to do nruns * N iterations (easily 10^9-10^10)
    ! Since we are performing O(N) operations on O(nruns) data, if we assume that 
    ! both N and nruns can go to infinity then this is a cpu-bound program.
    do i = 1, nruns/size + rest
        ! Use a different seed for each run
        seed = [(all_seeds((i-1)*sizer + j), j = 1, sizer)]
        call random_seed(put=seed)
        call random_number(p) ! get a sequence of random numbers
        xN = 0 ! reset the initial position for each walker
        do j = 1, N
            if (p(j) < 0.5_dp) then
                xN = xN - 1 ! left
            else
                xN = xN + 1 ! right
            end if
            x_i(j) = x_i(j) + xN
            x2_i(j) = x2_i(j) + xN**2
        end do
        P_N(xN) = P_N(xN) + 1 ! accumulate (only for j = N)
    end do

    end_time = MPI_Wtime()
    comp = end_time - start_time

    ! Reduce values on all processes to a single value (their sum) (O(size * N))
    ! I chose to consider it communication time, even though there is also a computation going on
    ! Usually, computation time refers to the time spent in performing local computations on each process
    start_time = MPI_Wtime()
    call MPI_Reduce(x_i, sum_x_i, N, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(x2_i, sum_x2_i, N, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(P_N, sum_P_N, 2*N+1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    end_time = MPI_Wtime()
    comm = end_time - start_time

    if (rank == 0) then ! Parallel I/O is not yet an option and besides, we only have to write O(N) elements

        ! Write in files, measure the time it takes to do this as well
        start_time = MPI_Wtime()

        call WriteAverageWalkers(N, nruns, sum_x_i, sum_x2_i)
        call WriteProbabilities(N, nruns, sum_P_N)
        call WriteDelta(N, nruns, sum_x_i, sum_x2_i)

        end_time = MPI_Wtime()
        writ = end_time - start_time
        call WriteTimings(size, N, nruns, comp, comm, writ)

    end if

    deallocate(seed, x_i, x2_i, p, sum_x_i, sum_x2_i, r, all_seeds, P_N, sum_P_N)

    call MPI_Finalize(ierr)

end program parallel_walkers
