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

999 FORMAT(1x,A15,2x,F8.4)

!~~~~~
! Parse input
  CALL EXECUTE_COMMAND_LINE('parse2')
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
  CALL EXECUTE_COMMAND_LINE('int2e')
  INQUIRE(file='error', EXIST=ex)
  IF (ex) THEN
    WRITE(*,*) "Error from int2e, exiting"
    STOP "Error from int2e"
  END IF

!~~~~
! Variational Hartree Fock
  CALL EXECUTE_COMMAND_LINE('scf')
  INQUIRE (file='error',EXIST=ex)
  IF (ex) THEN
    WRITE(*,*) "Error from scf, exiting"
    STOP "Error from scf"
  END IF

!~~~
! Output
  CALL CPU_TIME(timeF)
  WRITE(*,*)
  WRITE(*,999) "myQC ran in (s)", (timeF - timeS)

END PROGRAM myQC

