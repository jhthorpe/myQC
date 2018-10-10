#fortran compiler and flags, change as needed for each subdirectory (sorry)
#FC=gfortran
#FFLAGS=-O3 -fcheck=bounds
FC=ifort
FFLAGS=-O3

export $(FC) $(FFLAGS)

bindir:= ../bin

progdirs:= parser integrals I2G dens scf clean

.PHONY: all prog

all: prog
	$(FC) -o $(bindir)/myQC $(FFLAGS) myQC.f90 

prog:
		@set -e; \
	for i in $(progdirs); do \
		if [ -d $$i ]; then \
		  if [ -f $$i/Makefile ]; then \
		    $(MAKE) -C $$i all;\
		    fi; \
		  fi; \
	done;
