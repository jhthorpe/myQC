!//////////////////////////////////////////////////////////////////
!//            Performs 2 e- integrals for myQC 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//             
!//             Integration scheme a la McMurchie,Davidson 1978 
!//             
!//     WORK NOTE - location of coefficients in bas currently hardcoded 
!///////////////////////////////////////////////////////////////////

!=====================================================================
!                       MAIN 
PROGRAM int2e
  USE env
  USE basis
  USE aux

  IMPLICIT NONE

  ! Values
  ! xyz         : 2D dp, array of nuclear positions
  ! atoms       : 1D int, array of which atom is which
  ! fmem        : dp, free memory left in MB
  ! nnuc        : int, number of nuclii
  ! nelc        : int, number of electrons
  ! norb        : int, number of orbitals in molecule
  ! npri        : int, number of primatives
  ! bas         : 2D dp, basis for each atom: atom : orbital : [d,a]
  ! basinfo     : 2D int, array of basis information
  ! options     : 1D int, array of options
  ! F           : 2D dp, Fock matrix
  ! set		: 2D dp, array of exponential coefficients of sets: atom,orbital numbs
  ! setinfo	: 2D int, array of set information

  ! Variables  
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz,F
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: bas,set
  INTEGER, ALLOCATABLE, DIMENSION(:) :: basinfo,setinfo
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms,options
  REAL(KIND=8) :: timeS, timeF, fmem
  INTEGER :: nnuc,nelc,i,j,k,norb,npri,stat
  LOGICAL :: flag1,flag2,flag

  ! input managment 
  CALL CPU_TIME(timeS)
  WRITE(*,*)
  WRITE(*,*) "int1e called"
  CALL getenv(nnuc,nelc,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  ! build the basis set
  CALL buildBasis(options(2),atoms,bas,basinfo,set,setinfo,.FALSE.)

  ! check that Fock has been created 
  INQUIRE(file='Fuv',EXIST=flag2)
 
  ! check that we need to do these calculations
  INQUIRE(file='Fold',EXIST=flag2)
 
  ! logic gate
  IF (flag2) THEN
    WRITE(*,*) "Reading two electron integrals from Fuv"
  ELSE
    WRITE(*,*) "Calculating two electron integrals"
!    CALL proc2e(F,bas,basinfo,atoms,options,fmem,nnuc,xyz,norb)
  END IF

  CONTAINS
!===================================================================
!                       SUBROUTINES

!---------------------------------------------------------------------
! 		Processes the two electron integrals
!---------------------------------------------------------------------
  SUBROUTINE proc2e(F,bas,basinfo,atoms,options,fmem,nnuc,xyz,norb)
   IMPLICIT NONE
   ! Values
   ! F         : 2D dp, Fock matrix
   ! xyz       : 2D dp, array of nuclear positions
   ! atoms     : 1D int, array of which atom is which
   ! fmem      : dp, free memory left in MB
   ! nnuc      : int, number of nuclii
   ! norb      : int, number of orbitals in molecule
   ! bas       : 2D dp, basis for each atom: atom : orbital : [d,a]
   ! basinfo   : 2D int, array of basis information
   ! options   : 1D int, array of options
 
   !Inout
   REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: F
   REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
   REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
   INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo
   INTEGER, DIMENSION(0:), INTENT(IN) :: options, atoms
   REAL(KIND=8), INTENT(IN) :: fmem
   INTEGER, INTENT(IN) :: norb,nnuc


  END SUBROUTINE proc2e

!===================================================================
!                       FUNCTIONS
!----------



END PROGRAM int2e
