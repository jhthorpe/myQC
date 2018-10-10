#makefile for myQC

.PHONY: all

objects = basis.o env.o myQC.o

all: $(objects) 
	cp myQC.o $(bindir)/myQC

basis.o : basis.f90
	$(FC) $(FFLAGS) -c basis.f90
env.o : env.f90
	$(FC) $(FFLAGS) -c env.f90
myQC.o : myQC.f90
	$(FC) $(FFLAGS) myQC.f90 -o myQC.o

clean:
	rm -f *.o *.mod
