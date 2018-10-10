objects = dens.o
bindir:=../../bin
envdir = ../myQC
envprog = ../myQC/env.o ../myQC/basis.o


all: $(objects)
	cp dens.o $(bindir)/dens

dens.o : dens.f90
	$(FC) $(FFLAGS) -o dens.o dens.f90 $(envprog) -I$(envdir)

clean:
	rm -f *.o *.mod
