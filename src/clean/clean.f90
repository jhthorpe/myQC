!///////////////////////////////////////////////////////////////////
!//             Program for running Hartree Fock QC calcualations
!//
!//                     James H Thorpe, in Group of John Stanton
!//                     The University of Florida
!//
!///////////////////////////////////////////////////////////////////

PROGRAM clean
  IMPLICIT NONE

  CALL EXECUTE_COMMAND_LINE('rm fmem envdat nucpos radii error basinfo &
    Suv Huv Sold Hold setinfo XX Da dens.txt Guv Dold eig K_h K_l ijab ijba pqrs_AB 2> /dev/null')

END PROGRAM clean
