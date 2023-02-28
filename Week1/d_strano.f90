       program strano
          implicit none
          integer,parameter:: dp=selected_real_kind(13)
          real(kind=dp) :: dum
          real(kind=dp) :: a=1.0,b=0.0,c=-1.0
          
          print*," 1.0/0.0 = ",a/b
          print*," 0.0/0.0 = ",b/b
          print*,"-1.0/0.0 = ",c/b
          print*,10.0_dp**(-42)
          print*,10.0_dp**(-100)
          print*,epsilon(dum)
	end program strano
        
