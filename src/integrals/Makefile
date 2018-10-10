bindir:=../../bin
#FC=gfortran
#FFLAGS= -O3 -g -fcheck=bounds
#FFLAGS= -O3 -frecursive
FC= ifort
FFLAGS= -O3

all:
	$(FC) -o $(bindir)/int1e int1e.f90 ../env.f90 ../basis.f90 auxilary.f90 $(FFLAGS); \
	$(FC) -o $(bindir)/int2e int2e.f90 ../env.f90 ../basis.f90 auxilary.f90 $(FFLAGS)
