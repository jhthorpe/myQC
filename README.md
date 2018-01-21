# myQC
My crummy implementation of Hartree Fock

====  RUNNING  ===
1) Compile program with 'make all' in src. If you don't want to use gfortran, then you must edit each subdirectory's Makefile (sorry). The scf subdirectory MUST be linked to your blas/lapack libraries!
2) Add bin to your path.
3) Copy 'mybasis' and 'Ftab' (temporary fix) from bin into current directory.
4) In same directory, create or copy 'input'.
5) Call calculation with the 'myQC' executable.
6) Run 'clean' to remove intermediate files, if desired.

====  NOTES and BUGS  ===
- Current atoms: All 1st and 2nd row elements 
- Current basis sets: STO-3G
- Current systems: Atoms, Diatomics, cartesian molecules, closed shell only.
- Save considerable time by keeping the intermediate 'XX' file, which contains the two electron integrals.
- Two electron integrals need a more stringent incomplete error function evaluation, as there are some numerical problems with MO eigenvalues that should be degenerate.

One electron integrals implimented, two electron integrals implimented, but slow.
RHF and UHF SCF implimented. You can now visualize your molecule with the molden program.
 
