#hipergator
#export FC = ifort
#export FFLAGS = -O3
#export LINAL =
#export MKL = /apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin
#export bindir = /home/james.thorpe/bin/myQC/bin
#export LFLAGS=-rcs
#export LC=ar

#local
export FC = gfortran
export FFLAGS = -O3 -fcheck=bounds
export LINAL = -L/Users/jamesthorpe/LAPACK/lapack-3.7.0 -llapack -lblas
export bindir = /Users/jamesthorpe/local-hbar/myQC/myQC/bin
export LFLAGS=-rcs
export LC=ar


progdirs:= myQC parser clean integrals I2G dens scf ao2mo mp2 linal 

.PHONY: all prog clean

all: prog

prog:
	if [ ! -d $(bindir) ]; then mkdir -p $(bindir); fi
	@set -e; \
	for i in $(progdirs); do \
		if [ -d $$i ]; then \
		  if [ -f $$i/Makefile ]; then \
		    $(MAKE) -C $$i all;\
		    fi; \
		  fi; \
	done;

clean:
	@set -e; \
	for i in $(progdirs); do \
		if [ -d $$i ]; then \
		  if [ -f $$i/Makefile ]; then \
		    $(MAKE) -C $$i clean;\
		    fi; \
		  fi; \
	done;
