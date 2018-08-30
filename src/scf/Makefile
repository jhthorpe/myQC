bindir:= ../../bin
FC=gfortran
FFLAGS= -O3 -fcheck=bounds
#FC=ifort
#FFLAGS= -O3 -g
LINAL= -L/Users/jamesthorpe/LAPACK/lapack-3.7.0 -llapack -lblas

all:
	$(FC) -o $(bindir)/scf scf.f90 ../env.f90 ../basis.f90 $(FFLAGS) $(LINAL)

