       program strano
          implicit none
          real :: a=1.0,b=0.0,c=-1.0
          print*," 1.0/0.0 = ",a/b
          print*," 0.0/0.0 = ",b/b
          print*,"-1.0/0.0 = ",c/b
          print*,10.0**(-42)
          print*,10.0**(-100)
          print*,epsilon(a)
	end program strano
        
