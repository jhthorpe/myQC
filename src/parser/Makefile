bindir:= ../../bin

FC= gfortran
FFLAGS= -O3 -fcheck=bounds
#FC= ifort 
#FFLAGS= -O3

all:
	$(FC) $(FFLAGS) -o $(bindir)/parse parser.f90 ../env.f90
