# myQC
My crummy implementation of Hartree Fock

====  RUNNING  ===
1) compile program with 'maker.sh' in src. This MUST be done once, and with the correct linker to BLAS and LAPACK libraries.
2) Add bin to your path.
3) Copy 'mybasis' and 'Ftab' (temporary fix) from bin into current directory.
4) In same directory, create or copy 'input'.
5) Call calculation with the 'myQC' executable.
6) Run 'clean' to remove intermediate files

====  NOTES and BUGS  ===
- Current atoms: All 1st and 2nd row elements 
- Current basis sets: STO-3G
- Current systems: Atoms, Diatomics, cartesian molecules, closed shell only.

One electron integrals implimented, two electron integrals implimented, but slow.
