# myQC
My crummy implementation of Quantum Chemical Techniques

====  RUNNING  ===
1) Compile program with 'make all' in src. If you don't want to use gfortran, then you must edit each subdirectory's Makefile (sorry). The scf subdirectory's Makefile MUST be linked to your blas/lapack libraries!
2) Add bin to your path.
3) Copy 'mybasis' and 'Ftab' (temporary fix) from bin into current directory.
4) In same directory, create or copy 'ZMAT'.
5) Call calculation with the 'myQC' executable.
6) Run 'clean' to remove intermediate files, if desired.

====  NOTES and BUGS  ===
- Current atoms: All 1st and 2nd row elements 
- Current basis sets: STO-3G
- Save considerable time by keeping the intermediate 'XX' file, which contains the two electron integrals.
- Two electron integrals potentially need a more stringent incomplete error function evaluation.

WIKI: https://github.com/jhthorpe/myQC/wiki 
