# myQC
My crummy implementation of Hartree Fock

====  RUNNING  ===
1) Add bin to your path.
2) Copy 'mybasis' and 'Ftab' (temporary fix) from bin into current directory.
3) In same directory, create or copy 'input'.
4) Call calculation with the 'myQC' executable.
5) Run 'clean' to remove intermediate files

====  NOTES and BUGS  ===
- Current atoms: H,He,Li,Be,C 
- Current basis sets: STO-3G
- Current systems: Atoms, Diatomics, both closed shell.
- You can get around the diatomic molecule restriction by editing the nucpos and envdat files after running 'parse' but before running 'int1e'.

Currently, all I can do is the one electron integrals 
