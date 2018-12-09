!//////////////////////////////////////////////////////////////////
!//		Controls property calculations	
!//
!//		James H. Thorpe, in the Group of John Stanton
!//		The Univeristy of Florida
!//
!//////////////////////////////////////////////////////////////////

!---------------------------------------------------------------------
!	props	
!		James H. Thorpe
!	 	Dec 8, 2018	
!	- control program for property calculations 
!---------------------------------------------------------------------
  ! Variables
  ! noccA, nocc B	:	int, number of alpha,beta occupied orbitals 
  ! options		:	1D int, options array
  ! ntot		:	int, total number of orbitals
  ! nvrtA, nvrtB	:	int, number of alpha,beta virtual orbitals

PROGRAM props
  USE env
  IMPLICIT NONE

  !Internal
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms, options
  INTEGER, DIMENSION(0:1) :: line
  INTEGER, DIMENSION(0:0) :: mem_lvl
  REAL(KIND=8) :: fmem, timeS, timeF
  INTEGER :: nnuc,noccA,noccB,nvrtA,nvrtB,ntot,dummy
  LOGICAL :: flag
999 FORMAT(1x,A24,2x,F8.4)
  CALL CPU_TIME(timeS)
  WRITE(*,*) 
  WRITE(*,*)"               STARTING PROPERTY CALCULATIONS"
  WRITE(*,*)"------------------------------------------------------------"
  WRITE(*,*) "props called"
  WRITE(*,*)

  !Read enviromental data
  CALL getenv(nnuc,noccA,noccB,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  !propetry controls
  IF (options(16) .EQ. 0) THEN
    WRITE(*,*) "Property order   :  none"
  ELSE IF (options(16) .EQ. 1) THEN
    WRITE(*,*) "Property order   :  first"
    IF (options(1) .EQ. 0) THEN
      WRITE(*,*) "Property method  :  cphf"
      CALL EXECUTE_COMMAND_LINE('int1p')
      CALL EXECUTE_COMMAND_LINE('cphf')
    ELSE 
      WRITE(*,*) "Sorry, only cphf properties for RHF available"
    END IF
  ELSE IF (options(16) .EQ. 2) THEN
    WRITE(*,*) "Property order   :  second"
    IF (options(1) .EQ. 0) THEN
      WRITE(*,*) "Property method  :  cphf"
      CALL EXECUTE_COMMAND_LINE('int1p')
      CALL EXECUTE_COMMAND_LINE('cphf')
    ELSE 
      WRITE(*,*) "Sorry, only cphf properties for RHF available"
    END IF
  END IF

  CALL CPU_TIME(timeF)
  WRITE(*,999) "props completed in (s):", timeF-timeS

END PROGRAM props
