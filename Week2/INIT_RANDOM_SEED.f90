           SUBROUTINE init_random_seed()
! Credits:
! http://www.ewp.rpi.edu/hartford/~ernesto/S2011/EP/MaterialsforStudents/Zeno/BD_SUITE_1.8_2.5/SOURCE/INIT_RANDOM_SEED.f90
!
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed

            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))

            CALL SYSTEM_CLOCK(COUNT=clock)

            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)

            DEALLOCATE(seed)
          END SUBROUTINE

