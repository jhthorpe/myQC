objects = ao2mo.o
bindir = ../../bin
envdir = ../myQC
envprog = ../myQC/env.o

all : $(objects)
	cp ao2mo.o $(bindir)/ao2mo

ao2mo.o : ao2mo.f90
	$(FC) $(FFLAGS) -o ao2mo.o ao2mo.f90 $(envprog) -I$(envdir)

clean :
	rm -f *.o *.mod

