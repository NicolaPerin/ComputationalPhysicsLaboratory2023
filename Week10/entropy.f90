!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! entropy.f90
!
! calculates the entropy for each macrostate
! using the "coincidence method" of Ma
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module ma
  implicit none
  public :: start
  integer, public :: nl,nr,nexch,N
  integer, dimension(10), public :: mleft=0, mright=0
  integer, dimension(:), public, allocatable :: micro

contains

  subroutine start()
    !     initialize parameters
    integer :: il,ir
    integer, dimension(4) :: seed
    print*, " total number of particles N (<=10)>"
    read*, N
    print*, " number of particles 0<nl<N initially on the left (MACROstate)>"
    
    read*,nl
    if(nl<=0.or.nl>=N)then
       print*,' not acceptable, wrong nl'
       stop
       end if
    nr = N - nl  !  number of particles on the right
    print*, " number of exchanges >" ! no. of evolution steps of the macrostate
    read*, nexch
    allocate(micro(0:nexch))
    print*, " enter 4 seeds (integer positive)>"
    read*,seed
!    call random_seed(put=seed)
    micro(0) = 0
    write(*,fmt=*)'nleft =',nl
    write(*,fmt=*)'nright=',nr
    do il = 1,nl
       !     list left particles
       mleft(il) = il
       !     quantity characterizing the initial macrostate
       micro(0) = micro(0)*2 + 2
    end do
    do ir = 1,nr
       !     list right particles
       mright(ir) = ir + nl
    end do
!    print*,'microstate(0)=',micro(0)
!    write(*,fmt="(a,i2,a,10(1x,i2))")'nleft =',nl,' labels=',mleft
!    write(*,fmt="(a,i2,a,10(1x,i2))")'nright=',nr,' labels=',mright
  end subroutine start

  subroutine exch()
    !     exchange one particle on the left (ileft) 
    !     with one particle on the right (iright)
    real, dimension(2) :: r
    integer :: iexch,ileft,jleft,iright,jright 
    do  iexch = 1,nexch
       !        choose randomly the labels of the two particles
       call  random_number(r)
       ileft  = int (r(1)*nl + 1)    ! 1 =< ileft  =< nl
       iright = int (r(2)*nr + 1)    ! 1 =< iright =< nr 
       jleft  = mleft (ileft)
       jright = mright(iright)
       mleft (ileft)  = jright  ! new particle on the left
       mright(iright) = jleft   ! new particle on the right
       !     characterizing the microstate:
       micro(iexch) = micro(iexch-1) + 2**jright  - 2**jleft
!    print*,'microstate(',iexch,')=',micro(iexch)
!    write(*,fmt="(a,i2,a,10(1x,i2))")'nleft =',nl,' labels=',mleft
!    write(*,fmt="(a,i2,a,10(1x,i2))")'nright=',nr,' labels=',mright
    end do   
  end subroutine exch

  subroutine output()
    !     calculate the ratio of coincidences with respect to the total number
    !     of possible pairs, and consequently entropy
    !     
    integer :: ncoin, ncomp, iexch, jexch
    real :: rate,S
    ncoin = 0
    ncomp = nexch*(nexch-1)/2
    !     compare microstates: if coincident, count + 1;
    !     upgrade counter
    do  iexch = 1,nexch-1
       do jexch = iexch+1, nexch
          if (micro(iexch) == micro(jexch)) ncoin = ncoin + 1
       end do
    end do
    !     coincidence ratio
    rate = real(ncoin)/real(ncomp)
    if (rate > 0) then 
       S = log(1.0/rate)
       print*, " numerically estimated entropy: S=",S
       else
          print*," no coincidences! estimated entropy infinite!"
          end if
  end subroutine output

end module ma

program entropy
  use ma
  !     N:     total number of particles
  !     nl:    number of left particles (i.e. the MACROstate)
  !     mleft(),mright(): labels of left and right particles
  !     micro: a "global" label for a microstate, here defined through 
  !            mleft() : micro=sum_{il=1,nl} 2**(mleft(il)) 
  !     nexch: total number of exchanges (evolution steps of the macrostate)
  !            microst.)
  call start()
  !     the macrostates evolves (axchanging particles, the microstate changes)
  call exch()
  !     calculate the fraction of coincidence of microstates over all
  !     the possible coincidences with the microstates and the entropy
  call output()
  deallocate(micro)
end program entropy
