Folder containing the files and the report of the first assignment

  - utils.f90 contains two modules with functions and subroutines used mainly to write to file
  - parallel_walkers.f90 is the file containing the program
  - to compile them run "srun make" or "mpirun make" so that it compiles directly on the node it will be executed on (taking advantage of the -march=native flag)
  - to run the program use mpirun -np #number of processes parallel_walkers N, nruns (it takes two command line arguments, N and nruns, both integers)
  - scalability.sh is a bash script that automates the repeated runs needed to measure the speedup obtained using more than one process
