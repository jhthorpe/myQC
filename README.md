# myQC
My crummy implementation of Hartree Fock

To run: Add bin to your path.  pull 'mybasis' into current directory. In same directory, create or copy 'input'. Call calculation with the 'myQC' executable.

1) I suggest adding the bin directory into your path
2) edit the input file (in src) to what you need. The only atoms imlimented right now are H, He, Li, Be, and C, and only in STO-3G. 
3) run the 'clean' script (in bin) after each run, or you will get errors. 
4) I can only do calculations on closed shell atoms and molecules: and all molecules must be diatomic. 
5) You can actually get around this by running 'parse', and then editing the nucpos and envdat files before running 'int1e', which calls the basis set constructor

Currently, all I can do is construct the overlap integrals, but more will come soon.
