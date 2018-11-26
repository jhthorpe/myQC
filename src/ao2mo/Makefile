objects = ao2mo.o
bindir = ../../bin
envdir = ../myQC
envprog = ../myQC/env.o
MKL=/apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin

all : $(objects)
	cp ao2mo.o $(bindir)/ao2mo

ao2mo.o : ao2mo.f90
#hipergator
#	$(FC) $(FFLAGS) $(LINAL) -o ao2mo.o ao2mo.f90 -CB $(envprog) -I$(envdir) \
	-L$(MKL) \
	-Wl,-R$(MKL) -lmkl_lapack95_lp64 \
	-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm

#local
	$(FC) $(FFLAGS) -o ao2mo.o ao2mo.f90 $(envprog) -I$(envdir) $(LINAL)

clean :
	rm -f *.o *.mod

