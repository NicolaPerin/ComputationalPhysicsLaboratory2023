First week

In the Makefile, the SRCS variable uses the wildcard function to match all .f90 files in the directory. The OBJS and EXES variables use string substitution to replace the .f90 extension with .o and .x, respectively.

The % sign in the rules for %.o and %.x is a wildcard that matches any filename stem. So, for example, file1.o and file1.x will be generated from file1.f90.

The all target depends on all the executable targets in $ (EXES). The % pattern rules are used to generate the object files (%.o) from the corresponding .f90 files, and then generate the executables (%.x) from the object files. The -c flag is used with the $ (FC) command to generate object files, and -o $@ is used to name the object files the same as the corresponding source files.

Finally, the clean target removes all generated object files ($(OBJS)) and executables ($(EXES)).
