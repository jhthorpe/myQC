bindir:=../../bin
FC=gfortran
FFLAGS= -O3 -fcheck=bounds
#FC=ifort
#FFLAGS= -O3

all:
	$(FC) -o $(bindir)/dens dens.f90 ../env.f90 $(FFLAGS)
