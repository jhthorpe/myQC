#!/bin/bash
export Fdir="../bin"
gfortran -o ../bin/parse parser.f90 env.f90
gfortran -o ../bin/int1e -O3 ./integrals/int1e.f90 env.f90 basis.f90 ./integrals/auxilary.f90 -fcheck=bounds
gfortran -o ../bin/int2e -O3 ./integrals/int2e.f90 env.f90 basis.f90 ./integrals/auxilary.f90 -fcheck=bounds
gfortran -o ../bin/dens -O3 ./dens/dens.f90 env.f90 -fcheck=bounds
gfortran -o ../bin/scf -O3 ./scf/scf.f90 env.f90 -fcheck=bounds
gfortran -o ../bin/myQC -O3 myQC.f90 env.f90

