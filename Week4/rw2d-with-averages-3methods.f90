program d2rw
 implicit none

	integer, parameter :: rk = selected_real_kind(13)
	real(kind=rk) :: dR, x, y, rnd, rnd2, aux
	integer :: N, nwalks, sizer, i, j
	integer, dimension(:), allocatable :: seed
	real (kind=rk), dimension(:), allocatable :: xn, x2n, yn, y2n

	call random_seed(sizer)
	allocate(seed(sizer))
    print*, "Seed has ", sizer, "components; insert them or type '/' >"
	read*, seed
	call random_seed(put=seed)

	print*, "Number of steps"
	read*, N
	print*, "Number of walks"
	read*, nwalks

	allocate(xn(N))
	allocate(x2n(N))
	allocate(yn(N))
	allocate(y2n(N))
	
	xn = 0.0_rk
	x2n = 0.0_rk
	yn = 0.0_rk
	y2n = 0.0_rk

!----------------------------------------------------------------
! primo modo: step di mod.=1, angolo con distr.unif. tra 0 e pigr
!----------------------------------------------------------------
	do i=1,nwalks
		x = 0.0_rk
		y = 0.0_rk

		do j=1,N

			call random_number(rnd)
			rnd = rnd*2.0_rk*acos(-1.0_rk) ! tra 0 e 2\pi
			x = x + cos(rnd)
			y = y + sin(rnd)
			xn(j) = xn(j) + x
			x2n(j) = x2n(j) + x**2
			yn(j) = yn(j) + y
			y2n(j) = y2n(j) + y**2

		end do

	end do

	do i=1,N
		dR = x2n(i)/float(nwalks)-(xn(i)/float(nwalks))**2 + y2n(i)/float(nwalks)-(yn(i)/float(nwalks))**2
		write(unit=10, fmt=*) i, dR
	end do
!----------------------------------------------------------------
! secondo modo: x,y rnd in[-1;1], poi normalizzare per avere mod.=1
!----------------------------------------------------------------
	xn = 0.0_rk
	x2n = 0.0_rk
	yn = 0.0_rk
	y2n = 0.0_rk

	do i=1,nwalks
		x = 0.0_rk
		y = 0.0_rk

		do j=1,N

            do
			call random_number(rnd)
			call random_number(rnd2)
			if(rnd /= 0 .or. rnd2 /= 0) exit
            end do
            rnd = (rnd - 0.5_rk)*2.0_rk
            rnd2 = (rnd2 - 0.5_rk)*2.0_rk
			x = x + rnd/sqrt(rnd**2 + rnd2**2)
			y = y + rnd2/sqrt(rnd**2 + rnd2**2)
			xn(j) = xn(j) + x
			x2n(j) = x2n(j) + x**2
			yn(j) = yn(j) + y
			y2n(j) = y2n(j) + y**2

		end do

	end do	


	do i=1,N
		dR = x2n(i)/float(nwalks)-(xn(i)/float(nwalks))**2 + y2n(i)/float(nwalks)-(yn(i)/float(nwalks))**2
		write(unit=11, fmt=*) i, dR
	end do
!------------------------------------------------------------
! terzo modo: x,y rnd in[-sqrt(1.5);sqrt(1.5)], è mod.=1 IN MEDIA
!------------------------------------------------------------
	xn = 0.0_rk
	x2n = 0.0_rk
	yn = 0.0_rk
	y2n = 0.0_rk

	do i=1,nwalks
		x = 0.0_rk
		y = 0.0_rk
        aux = 0.0_rk  ! check se IN MEDIA lo step ha mod.=1

		do j=1,N

			call random_number(rnd)
			call random_number(rnd2)
			
			rnd = (rnd - 0.5_rk)*2.0_rk*sqrt(1.5_rk)
			rnd2 = (rnd2 - 0.5_rk)*2.0_rk*sqrt(1.5_rk)

			x = x + rnd
			y = y + rnd2
			xn(j) = xn(j) + x
			x2n(j) = x2n(j) + x**2
			yn(j) = yn(j) + y
			y2n(j) = y2n(j) + y**2
            aux = aux + rnd**2 + rnd2**2

		end do
            aux = aux/N
            print*," metodo n.3: <Δx^2+Δy^2> in walk n.",i," :",aux
	end do	


	do i=1,N
		dR = x2n(i)/float(nwalks)-(xn(i)/float(nwalks))**2 + y2n(i)/float(nwalks)-(yn(i)/float(nwalks))**2
		write(unit=12, fmt=*) i, dR, xn(i)/float(nwalks)
	end do

	
end program d2rw


