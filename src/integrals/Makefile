objects = auxilary.o int1e.o int2e.o 
envdir = ../myQC
envprog = ../myQC/env.o ../myQC/basis.o
bindir = ../../bin

.PHONY: all

all: $(objects)
	cp int1e.o $(bindir)/int1e
	cp int2e.o $(bindir)/int2e
	
auxilary.o : auxilary.f90
	$(FC) $(FFLAGS) -c auxilary.f90 
int1e.o : auxilary.o int1e.f90 
	$(FC) $(FFLAGS) -o int1e.o int1e.f90 auxilary.o $(envprog) -I$(envdir)
int2e.o : auxilary.o int2e.f90 
	$(FC) $(FFLAGS) -o int2e.o int2e.f90 auxilary.o $(envprog) -I$(envdir)

clean:
	rm -r *.o *.mod
