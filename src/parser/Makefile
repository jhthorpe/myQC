bindir:= ../../bin

FC= gfortran
FFLAGS= -O2 -fcheck=bounds

all:
	$(FC) $(FFLAGS) -o $(bindir)/parse parser.f90 ../env.f90
