bindir:= ../../bin
FC=gfortran
FFLAGS=-O2 -fcheck=bounds

all:
	$(FC) -o $(bindir)/RHFI2G RHFI2G.f90 ../env.f90 $(FFLAGS); \
	$(FC) -o $(bindir)/UHFI2G UHFI2G.f90 ../env.f90 $(FFLAGS)
	
