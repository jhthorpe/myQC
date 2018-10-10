objects = scf.o
bindir:= ../../bin
envdir = ../myQC
envprog = ../myQC/env.o ../myQC/basis.o
#MKL=/apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin

all: $(objects)
#	$(FC) -o $(bindir)/scf scf.f90 ../env.f90 ../basis.f90 $(FFLAGS) $(LINAL) \
#	-L$(MKL) \
#	-Wl,-R$(MKL) -lmkl_lapack95_lp64 \
#	-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm
	cp scf.o $(bindir)/scf

scf.o : scf.f90
	$(FC) $(FFLAGS) $(LINAL) -o scf.o scf.f90 $(envprog) -I$(envdir) \
	-L$(MKL) \
	-Wl,-R$(MKL) -lmkl_lapack95_lp64 \
	-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm

clean:
	rm -f *.o *.mod


