objects = mp2.o
bindir = ../../bin
envdir = ../myQC
envprog = ../myQC/env.o

all : $(objects)
	cp mp2.o $(bindir)/mp2

mp2.o : mp2.f90
#local
	$(FC) $(FFLAGS) -o mp2.o mp2.f90 -CB $(envprog) -I$(envdir) $(LINAL)

clean :
	rm -f *.o *.mod

