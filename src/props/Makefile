objects = cphf.o props.o
bindir = ../../bin
envdir = ../myQC
envprog = ../myQC/env.o

all : $(objects)
	cp cphf.o $(bindir)/cphf
	cp props.o $(bindir)/props

cphf.o : cphf.f90
	$(FC) $(FFLAGS) -o cphf.o cphf.f90 $(envprog) -I$(envdir) $(LINAL)

props.o : props.f90 cphf.f90
	$(FC) $(FFLAGS) -o props.o props.f90 $(envprog) -I$(envdir) 

clean :
	rm -f *.o *.mod

