!///////////////////////////////////////////////////////////////////
!//		Program for running Hartree Fock QC calcualations
!//
!//			James H Thorpe, in Group of John Stanton
!//			The University of Florida
!//
!///////////////////////////////////////////////////////////////////

PROGRAM myQC
  USE env
  IMPLICIT NONE

  ! timeS	: real dp, starting CPU time for program 
  ! timeF	: real dp, ending CPU time for program 
  
  INTEGER,DIMENSION(:),ALLOCATABLE :: options
  REAL(KIND=8) :: timeS, timeF
  INTEGER :: nopt
  LOGICAL :: ex

  CALL CPU_TIME(timeS)
  WRITE(*,*) "                 Starting myQC"  
  ex = .FALSE.

999 FORMAT(1x,A15,2x,F8.4)

!~~~~~
! Parse input
  CALL EXECUTE_COMMAND_LINE('parse')
  INQUIRE(file='error', EXIST=ex)
  IF (ex) THEN
    WRITE(*,*) "Error from parse, exiting"
    STOP "Error from parse"
  END IF

  OPEN(unit=100,file='envdat',status='old',access='sequential')
  READ(100,*) nopt
  READ(100,*) nopt,nopt
  READ(100,*) nopt
  ALLOCATE(options(0:nopt-1))
  READ(100,*) options
  CLOSE(unit=100)

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

!~~~~
IF (options(1) .GT. 0) THEN
  CALL EXECUTE_COMMAND_LINE('ao2mo')
  INQUIRE (file='error',EXIST=ex)
  IF (ex) THEN
    WRITE(*,*) "Error from ao2mo, exiting"
    STOP "Error from ao2mo"
  END IF
END IF
!~~~
IF (options(1) .EQ. 1) THEN
  CALL EXECUTE_COMMAND_LINE('mp2')
  INQUIRE (file='error',EXIST=ex)
  IF (ex) THEN
    WRITE(*,*) "Error from mp2, exiting"
    STOP "Error from mp2"
  END IF
END IF

!~~~
! Output
  CALL CPU_TIME(timeF)
  WRITE(*,*)
  WRITE(*,999) "myQC ran in (s)", (timeF - timeS)

END PROGRAM myQC

