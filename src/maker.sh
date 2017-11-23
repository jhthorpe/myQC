gfortran -o ../bin/parse parser.f90 env.f90
gfortran -o ../bin/int1e ./integrals/int1e.f90 env.f90 basis.f90
gfortran -o ../bin/myQC myQC.f90 env.f90

