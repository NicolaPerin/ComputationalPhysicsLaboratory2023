program write_files
!  test program to open n files and write data
   integer :: i,n
   character(len = 15) :: file_name
   n = 11
   do i = 1,n
!     assign number.dat to file_name using write statement
      write(unit=file_name,fmt="(i2.2,a)") i,".dat" 
!     // is concatenation operator
      file_name = "config"//file_name
      open (unit=1,file=file_name,action="write",status="replace")
      write (unit=1, fmt=*) i*i,file_name
      close(unit=1)
   end do
end program write_files
