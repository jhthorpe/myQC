!//////////////////////////////////////////////////////////////////
!//		Performs CPHF Calculation
!//
!//		James H. Thorpe, in the Group of John Stanton
!//		The Univeristy of Florida
!//
!//////////////////////////////////////////////////////////////////

!---------------------------------------------------------------------
!	cphf	
!		James H. Thorpe
!	 	Dec 8, 2018	
!	- control program for cphf calculations 
!---------------------------------------------------------------------
  ! Variables
  ! noccA, nocc B	:	int, number of alpha,beta occupied orbitals 
  ! options		:	1D int, options array
  ! ntot		:	int, total number of orbitals
  ! nvrtA, nvrtB	:	int, number of alpha,beta virtual orbitals

PROGRAM cphf
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

  CALL CPU_TIME(timeS)
  WRITE(*,*) 
  WRITE(*,*) "cphf called"
  WRITE(*,*)

  !Read enviromental data
  CALL getenv(nnuc,noccA,noccB,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  !get orbital data
  OPEN(unit=100,file='basinfo',status='old',access='sequential')
  READ(100,*) line
  ntot = line(1)
  CLOSE(unit=1)
  nvrtA = ntot - noccA
  nvrtB = ntot - noccB

  CALL CPU_TIME(timeF)
  WRITE(*,*) "cphf completed in (s): ", timeF-timeS

END PROGRAM cphf
