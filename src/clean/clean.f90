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
    Suv Huv Sold Hold setinfo XX Da dens.txt Guv Dold eig K_h K_l ijab_AA &
    ijab_BB ijab_AB ajib_AA ajbi_AA ajib_AB ajib_BB ajbi_BB ints_AA ints_BB &
    2> /dev/null')

END PROGRAM clean
