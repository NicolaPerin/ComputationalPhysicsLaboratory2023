FC = mpif90
FFLAGS = -march=native -O2 -fcheck=all
SRCS = utils.f90 parallel_walkers.f90
EXE = parallel_walkers

$(EXE): $(SRCS)
		$(FC) $(FFLAGS) $(SRCS) -o $(EXE)

clean:
		rm -f $(EXE)
