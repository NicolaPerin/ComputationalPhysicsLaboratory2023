.PHONY: all clean flush

FC = gfortran
FFLAGS = -O2 -ftree-vectorize -fcheck=all -Wall

SRC = FHP1.f90
OBJ = $(SRC:.f90=.o)
EXEC = FHP1.x

all: $(EXEC)

$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) $^ -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXEC)

flush:
	rm -f avg_vel*
