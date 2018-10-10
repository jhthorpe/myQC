bindir:= ../../bin
#FC=gfortran
#FFLAGS= -O3 -fcheck=bounds
FC=ifort
FFLAGS= -O3
#LINAL= -L/Users/jamesthorpe/LAPACK/lapack-3.7.0 -llapack -lblas
MKL=/apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin

all:
	$(FC) -o $(bindir)/scf scf.f90 ../env.f90 ../basis.f90 $(FFLAGS) $(LINAL) \
	-L$(MKL) \
	-Wl,-R$(MKL) -lmkl_lapack95_lp64 \
	-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm

