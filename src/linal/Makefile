objects = linal.o linal.mod
bindir:= ../../bin
envdir = ../myQC
envprog = ../myQC/env.o ../myQC/basis.o
MKL=/apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin

#FC=gfortran
#FFLAGS=-O3

all: $(objects)
	$(LC) $(LFLAGS) linal.a $(objects)

linal.o : linal.f90
#mkl on hpg
#	$(FC) $(FFLAGS) $(LINAL) -c linal.f90 $(envprog) -I$(envdir) \
	-L$(MKL) \
	-Wl,-R$(MKL) -lmkl_lapack95_lp64 \
	-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm
#gnu on hpg
#	$(FC) $(FFLAGS) -c linal.f90 $(envprog) -I$(envdir)\
	-L/home/james.thorpe/bin/lapack-gnu/lapack-release/bin
#	$(FC) $(FFLAGS) -c linal.f90 /home/james.thorpe/bin/lapack-gnu/lapack-release/bin/*.a \
	-L/home/james.thorpe/bin/lapack-gnu/lapack-release/bin -llapack -lblas
#local
	$(FC) $(FFLAGS) $(LINAL) -c linal.f90 $(envprog) -I$(envdir) \
	$(LINAL)


linal.mod : linal.f90 linal.o

clean:
	rm -f *.o *.mod *.a


