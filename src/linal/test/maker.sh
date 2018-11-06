#!/bin/bash
export MKL=/apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin
ifort -O3 test.f90 ../linal.a -I.. -L$MKL -Wl,-R$MKL -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm
#gfortran test.f90 /home/james.thorpe/bin/lapack-gnu/lapack-release/bin/*.a ../linal.a -I.. -L/home/james.thorpe/bin/lapack-gnu/lapack-release/bin

