!///////////////////////////////////////////////////////////////////
!//		Program for running Hartree Fock QC calcualations
!//
!//			James H Thorpe, in Group of John Stanton
!//			The University of Florida
!//
!///////////////////////////////////////////////////////////////////

PROGRAM myQC
!  USE parser  !input parser

  IMPLICIT NONE

  ! timeS	: real dp, starting CPU time for program 
  ! timeF	: real dp, ending CPU time for program 
  
 
  REAL(KIND=8) :: timeS, timeF
  LOGICAL :: ex

  CALL CPU_TIME(timeS)
  WRITE(*,*) "Starting myQC"  
  ex = .FALSE.

!~~~~~
! Parse input
  CALL EXECUTE_COMMAND_LINE('parse')
  INQUIRE(file='error', EXIST=ex)
  IF (ex) THEN
    WRITE(*,*) "Error from parse, exiting"
    STOP "Error from parse"
  END IF

!~~~~~~
! 1 electron integrals
  CALL EXECUTE_COMMAND_LINE('int1e')
  INQUIRE(file='error', EXIST=ex)
  IF (ex) THEN
    WRITE(*,*) "Error from int1e, exiting"
    STOP "Error from int1e"
  END IF
!~~~~~
! 2 electron integrals

!~~~~
! Variational Hartree Fock

!~~~
! Output
  CALL CPU_TIME(timeF)

  WRITE(*,*) "myQC ran in (s) : ", (timeF - timeS)

END PROGRAM myQC

