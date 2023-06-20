module BitConversionModule
    implicit none
  
    contains
  
    subroutine IntegerToBinaryAndPrint(num)
      implicit none
      integer, intent(in) :: num
      integer :: i, bit

    print*, num
      do i = 32, 1, -1
        bit = BTEST(num, i - 1)
        write(*, '(I0)', advance='no') bit
      end do
      write(*, *)
  
    end subroutine IntegerToBinaryAndPrint
  
end module BitConversionModule

module Variables
    implicit none

    integer, dimension(:), allocatable :: seed
    real :: start_time, end_time
    integer, parameter :: Lx = 20, Ly = 10
    integer, parameter :: RI = 1, RD = 2, LD = 4, LE = 8, LU = 16, RU = 32, S = 64, B = 128
    integer, parameter :: NUM_CHANNELS = 7, NUM_BITS = 8, NUM_RULES = 256
    integer :: i, j, sevenParticles, highBits, lowBits, clock, sizer
    integer, dimension(NUM_RULES) :: rule = 0
    integer, dimension(Lx, Ly) :: Lattice, newLattice
    double precision :: numParticles, density
    real, dimension(Lx, Ly) :: r
    double precision, parameter :: SQRT3_OVER2 = sqrt(3.0d0) / 2.0d0
    double precision, parameter :: SQRT2 = sqrt(2.0d0)
    double precision, dimension(7), parameter :: ux = (/1.0d0, 0.5d0, -0.5d0, -1.0d0, -0.5d0, 0.5d0, 0.0d0/)
    double precision, dimension(7), parameter :: uy = (/0.0d0, -SQRT3_OVER2, -SQRT3_OVER2, 0.0d0, SQRT3_OVER2, SQRT3_OVER2, 0.0d0/)
    double precision, dimension(NUM_RULES) :: vx ! averaged x velocity for every site configuration
    double precision, dimension(NUM_RULES) :: vy ! averaged y velocity for every site configuration

  end module Variables

module RuleModule
    use Variables
    implicit none
  
    contains

    subroutine DefineRules()

      ! three particles zero momentum rules
      rule(IOR(IOR(LU, LD), RI)) = IOR(IOR(RU, LE), RD); rule(IOR(IOR(RU, LE), RD)) = IOR(IOR(LU, LD), RI) ! rule(21) = 42; rule(42) = 21

      ! three particles rules with unperturbed particle
      rule(IOR(IOR(RU, LU), LD)) = IOR(IOR(LU, LE), RI); rule(IOR(IOR(LU, LE), RI)) = IOR(IOR(RU, LU), LD) ! rule(52) = 25; rule(25) = 52
      rule(IOR(IOR(RU, LU), RD)) = IOR(IOR(RU, LE), RI); rule(IOR(IOR(RU, LE), RI)) = IOR(IOR(RU, LU), RD) ! rule(50) = 41; rule(41) = 50
      rule(IOR(IOR(RU, LD), RD)) = IOR(IOR(LE, RD), RI); rule(IOR(IOR(LE, RD), RI)) = IOR(IOR(RU, LD), RD) ! rule(38) = 11; rule(11) = 38
      rule(IOR(IOR(LU, LD), RD)) = IOR(IOR(LE, LD), RI); rule(IOR(IOR(LE, LD), RI)) = IOR(IOR(LU, LD), RD) ! rule(22) = 13; rule(13) = 22
      rule(IOR(IOR(RU, LD), RI)) = IOR(IOR(LU, RD), RI); rule(IOR(IOR(LU, RD), RI)) = IOR(IOR(RU, LD), RI) ! rule(37) = 19; rule(19) = 37
      rule(IOR(IOR(LU, LE), RD)) = IOR(IOR(RU, LE), LD); rule(IOR(IOR(RU, LE), LD)) = IOR(IOR(LU, LE), RD) ! rule(26) = 44; rule(44) = 26

      ! two particles cyclic rules
      rule(IOR(LE, RI)) = IOR(RU, LD) ! rule(9)  = 36
      rule(IOR(RU, LD)) = IOR(LU, RD) ! rule(36) = 18
      rule(IOR(LU, RD)) = IOR(LE, RI) ! rule(18) = 9

      ! four particles cyclic rules
      rule(IOR(IOR(IOR(RU, LU), LD), RD)) = IOR(IOR(IOR(RU, LE), LD), RI) ! rule(54) = 45
      rule(IOR(IOR(IOR(RU, LE), LD), RI)) = IOR(IOR(IOR(LU, LE), RD), RI) ! rule(45) = 27
      rule(IOR(IOR(IOR(LU, LE), RD), RI)) = IOR(IOR(IOR(RU, LU), LD), RD) ! rule(27) = 54

      ! stationary particle creation rules
      rule(IOR(LU, RI)) = IOR(RU, S) ! rule(17) = 96
      rule(IOR(RU, LE)) = IOR(LU, S) ! rule(40) = 80
      rule(IOR(LU, LD)) = IOR(LE, S) ! rule(20) = 72
      rule(IOR(LE, RD)) = IOR(LD, S) ! rule(10) = 68
      rule(IOR(LD, RI)) = IOR(RD, S) ! rule(5)  = 66
      rule(IOR(RD, RU)) = IOR(RI, S) ! rule(34) = 65

      rule(IOR(IOR(IOR(IOR(LU, LE), LD), RD), RI)) = IOR(IOR(IOR(IOR(RU, LE), LD), RD), S) ! rule(31) = 110
      rule(IOR(IOR(IOR(IOR(RU, LE), LD), RD), RI)) = IOR(IOR(IOR(IOR(LU, LD), RD), RI), S) ! rule(47) = 87
      rule(IOR(IOR(IOR(IOR(RU, LU), LD), RD), RI)) = IOR(IOR(IOR(IOR(RU, LE), RD), RI), S) ! rule(55) = 107
      rule(IOR(IOR(IOR(IOR(RU, LU), LE), RD), RI)) = IOR(IOR(IOR(IOR(RU, LU), LD), RI), S) ! rule(59) = 117
      rule(IOR(IOR(IOR(IOR(RU, LU), LE), LD), RI)) = IOR(IOR(IOR(IOR(RU, LU), LE), RD), S) ! rule(61) = 122
      rule(IOR(IOR(IOR(IOR(RU, LU), LE), LD), RD)) = IOR(IOR(IOR(IOR(LU, LE), LD), RI), S) ! rule(62) = 93

      sevenParticles = IOR(IOR(IOR(IOR(IOR(IOR(RU, LU), LE), LD), RD), RI), S) ! 127

      ! add all rules indexed with a stationary particle (dual rules)
      do i = 1, S
        rule(IEOR(i, sevenParticles)) = IEOR(rule(i), sevenParticles)
      end do

      ! add rules to bounce back at the barriers
      do i = B, NUM_RULES
        highBits = IAND(i, IOR(IOR(LE, LU), RU))
        lowBits  = IAND(i, IOR(IOR(RI, RD), LD))
        rule(i) = IOR(IOR(B, RSHIFT(highBits, 3)), LSHIFT(lowBits, 3))
      end do

    end subroutine DefineRules

end module RuleModule

module Utilities
    use Variables
    implicit none

    contains

    subroutine SetAverageSiteVelocities()
    ! For every particle in site configuration i, calculate total net velocity and place it in vx[i], vy[i]

      integer :: j, dir
      
      do j = 1, NUM_RULES
        do dir = 1, NUM_CHANNELS
          if (IAND(j, LSHIFT(1, dir)) /= 0) then
            vx(j) = vx(j) + ux(dir)
            vy(j) = vy(j) + uy(dir)
          end if
        end do
      end do

    end subroutine SetAverageSiteVelocities

    subroutine Initialize(Lx, Ly, density, Lattice, r)
      integer :: i, j
      integer, intent(in) :: Lx, Ly
      double precision, intent(in) :: density
      real, dimension(Lx, Ly) :: r
      integer, dimension(Lx, Ly), intent(inout) :: Lattice

      call random_number(r)

      do j = 1, Ly
        Lattice(1, j) = B; Lattice(Lx, j) = B
      end do

      do i = 1, Lx
        Lattice(i, 1) = B; Lattice(i, Ly) = B
      end do

      do i = 2, Lx - 1
        do j = 2, Ly - 1
          if (r(i, j) < density) then
            Lattice(i, j) = sevenParticles
          else
            Lattice(i, j) = 0
          end if
        end do
      end do

    end subroutine Initialize

    subroutine step(Lx, Ly, Lattice, newLattice)
      integer, intent(in) :: Lx, Ly
      integer :: i, j, site1, site2
      integer, dimension(Lx) :: left, cent, rght
      integer, dimension(Lx, Ly), intent(inout) :: Lattice, newLattice
    
      do i = 1, Lx
        left = Lattice(mod(i-2, Lx) + 1, :)
        cent = Lattice(i, :)
        rght = Lattice(mod(i, Lx) + 1, :)
    
        do j = 2, Ly-2, 2
          
          site1 = Lattice(i, j)
          site2 = Lattice(i, j+1)
    
          rght(j-1) = IAND(site1, RD)
          cent(j-1) = IAND(site1, LD)
          rght(j) = IAND(site1, RI)
          cent(j) = IOR(IAND(site1, IOR(S, B)), IAND(site2, RD))
          left(j) = IOR(IAND(site1, LE), IAND(site2, LD))
          rght(j+1) = IOR(IAND(site1, RI), IAND(site2, RI))
          cent(j+1) = IOR(IAND(site1, LU), IAND(site2, IOR(S, B)))
          left(j+1) = IAND(site2, LE)
          cent(j+2) = IAND(site2, RU)
          left(j+2) = IAND(site2, LU)
        end do
      end do
    
      Lattice = newLattice
    end subroutine step

end module Utilities

program LGAFF
    use BitConversionModule
    use Variables
    use RuleModule
    use Utilities
    implicit none

    call random_seed(sizer)
    allocate(seed(sizer))
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, sizer) /)
    call random_seed(put=seed)

    density = 0.3

    call DefineRules()
    call Initialize(Lx, Ly, density, Lattice, r)

    do i = 1, Lx
      do j = 1, Ly
        write(*, '(I3.3, A)', advance='no') Lattice(i, j), ' '
      end do
      write (*,*)
    end do
  
end program LGAFF

!do i = 1, size(rule)
  !call IntegerToBinaryAndPrint(rule(i))
!end do

!do i = 1, NUM_RULES
!  print*, vx(i), vy(i)
!end do
