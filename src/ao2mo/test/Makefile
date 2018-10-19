#hpg
#FC = ifort
#FFLAGS = -O3
#MKL = /apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin

#local
FC = gfortran 
FFLAGS= -O3
LINAL = -L/Users/jamesthorpe/LAPACK/lapack-3.7.0 -llapack -lblas 
tests.o : tests.f90
#hpg
#	$(FC) $(FFLAGS) -o tests.o tests.f90 \
#	-L$(MKL) \
#	-Wl,-R$(MKL) -lmkl_lapack95_lp64 \
#	-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm

#local
	$(FC) $(FFLAGS) -o tests.o tests.f90 $(LINAL)

