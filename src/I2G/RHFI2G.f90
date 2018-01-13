!//////////////////////////////////////////////////////////////////
!//          Constructs Guv(2e-) Matrix for RHF system 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//             
!//  WORK NOTE - currently this uses a very inefficient algrothim,
!//    at some point the in the future, it should be implimented
!//    using: Raffenetti, Chemical Physics Letters, 1973             
!///////////////////////////////////////////////////////////////////

!=====================================================================
!                       MAIN 

PROGRAM RHFI2G
  USE env
  IMPLICIT NONE

  ! Values
  ! xyz         : 2D dp, array of nuclear positions
  ! atoms       : 1D int, array of which atom is which
  ! fmem        : dp, free memory left in MB
  ! nnuc        : int, number of nuclii
  ! nelc        : int, number of electrons
  ! options     : 1D int, array of options
  ! Guv		: 2D dp, 2e- part of matrix
  ! Da		: 2D dp, density matrix 
  ! II		: 2D dp, 2e- AO integral matrix

  ! Internal
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz,Guv,Da,II
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms, options
  INTEGER, DIMENSION(0:1) :: line
  REAL(KIND=8) :: fmem, timeS, timeF
  INTEGER :: nnuc, nelc, dummy, norb
  LOGICAL :: flag

  CALL getenv(nnuc,nelc,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  !get number of orbitals
  OPEN(unit=1,file='basinfo',status='old',access='sequential')
  READ(1,*) line
  norb = line(1)
  CLOSE(unit=1)

  !Allocate memory


END PROGRAM RHFI2G
