bindir:= ../../bin
envdir = ../myQC
envprog = ../myQC/env.o
objects = RHFI2G.o UHFI2G.o

all: $(objects)
	cp RHFI2G.o $(bindir)/RHFI2G
	cp UHFI2G.o $(bindir)/UHFI2G

RHFI2G.o : RHFI2G.f90
	$(FC) $(FFLAGS) -o RHFI2G.o RHFI2G.f90 $(envprog) -I$(envdir) 
UHFI2G.o : UHFI2G.f90
	$(FC) $(FFLAGS) -o UHFI2G.o UHFI2G.f90 $(envprog) -I$(envdir) 
	
clean:
	rm -f *.o *.mod
