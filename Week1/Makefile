FC=gfortran
FFLAGS=-Wall -O2
SRCS=$(wildcard *.f90)
OBJS=$(SRCS:.f90=.o)
EXES=$(SRCS:.f90=.x)

.PHONY: all clean

all: $(EXES)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.x: %.o
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(EXES)
