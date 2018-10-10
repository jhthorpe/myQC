bindir:= ../../bin
#FC=gfortran
#FFLAGS=-O3 -fcheck=bounds
FC=ifort
FFLAGS=-O3

all:
	$(FC) -o $(bindir)/RHFI2G RHFI2G.f90 ../env.f90 $(FFLAGS); \
	$(FC) -o $(bindir)/UHFI2G UHFI2G.f90 ../env.f90 $(FFLAGS)
	
