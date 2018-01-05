# myQC
My crummy implementation of Hartree Fock

====  RUNNING  ===
1) Add bin to your path.
2) Copy 'mybasis' and 'Ftab' (temporary fix) from bin into current directory.
3) In same directory, create or copy 'input'.
4) Call calculation with the 'myQC' executable.
5) Run 'clean' to remove intermediate files

====  NOTES and BUGS  ===
- Current atoms: All 1st and 2nd row elements 
- Current basis sets: STO-3G
- Current systems: Atoms, Diatomics, cartesian molecules, closed shell only.

One electron integrals implimented, two electron integrals implimented, but slow.
