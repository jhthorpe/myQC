objects = cis.o
bindir:= ../../bin
envdir = ../myQC
envprog = ../myQC/env.o ../myQC/basis.o
MKL=/apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin

all: $(objects)
#	$(FC) -o $(bindir)/cis cis.f90 ../env.f90 ../basis.f90 $(FFLAGS) $(LINAL) \
#	-L$(MKL) \
#	-Wl,-R$(MKL) -lmkl_lapack95_lp64 \
#	-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm
	cp cis.o $(bindir)/cis

cis.o : cis.f90
#hipergator
#	$(FC) $(FFLAGS) $(LINAL) -o cis.o cis.f90 $(envprog) -I$(envdir) \
	-L$(MKL) \
	-Wl,-R$(MKL) -lmkl_lapack95_lp64 \
	-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm

#local
	$(FC) $(FFLAGS) $(LINAL) -o cis.o cis.f90 $(envprog) -I$(envdir) \
	$(ML)/linal.a -I$(ML)

clean:
	rm -f *.o *.mod


