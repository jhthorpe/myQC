 !//////////////////////////////////////////////////////////////////
 !//            Calculates the density matrix for hartree fock 
 !//
 !//              James H Thorpe, in Group of John Stanton
 !//              The University of Florida
 !//             
 !///////////////////////////////////////////////////////////////////

!=====================================================================
!                       MAIN 

PROGRAM dens
  USE env
  IMPLICIT NONE

  ! Values
  ! xyz         : 2D dp, array of nuclear positions
  ! atoms       : 1D int, array of which atom is which
  ! fmem        : dp, free memory left in MB
  ! nnuc        : int, number of nuclii
  ! nelc        : int, number of electrons
  ! options     : 1D int, array of options

  ! Internal
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms, options
  REAL(KIND=8) :: fmem, timeS, timeF
  INTEGER :: nnuc, nelc, dummy
  LOGICAL :: flag

  CALL getenv(nnuc,nelc,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP
  
  !redirect to density subroutines 
  IF (options(2) .EQ. 0) THEN       !RHF
    CALL densRHF(nnuc,nelc,atoms,fmem,options)
  ELSE IF (options(2) .EQ. 1) THEN  !UHF
    CALL densUHF(nnuc,nelc,atoms,fmem,options) 
  ELSE 
    WRITE(*,*) "dens: sorry, that method is not implimented yet"
    STOP "bad method in dens"
  END IF

  CONTAINS

!=====================================================================
!                       SUBROUTINES

!---------------------------------------------------------------------
!		RHF Density Matrix
!---------------------------------------------------------------------
  SUBROUTINE densRHF(nnuc,nelc,atoms,fmem,options)
    IMPLICIT NONE

    ! Values
    ! nnuc	: int, number of nuclii
    ! nelc	: int, number of electrons
    ! atoms	: 1D int array of atoms
    ! options	: 1D int, options of program
    ! fmem	: dp, free memory left in MB
    ! norb	: int, number of orbitals in molecule
    ! Cuv	: 2D dp, molecular orbital coefficients
    
    !Inout
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms,options
    REAL(KIND=8), INTENT(IN) :: fmem
    INTEGER, INTENT(IN) :: nnuc,nelc

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Cuv
    INTEGER, DIMENSION(0:1) :: line
    INTEGER :: norb
    INTEGER :: i,j  

    !get number of orbitals
    OPEN(unit=1,file='basinfo',status='old',access='sequential')
    READ(1,*) line
    norb = line(1)
    CLOSE(unit=1) 

    ALLOCATE(Cuv(0:norb-1,0:norb-1))

    DEALLOCATE(Cuv) 

  END SUBROUTINE densRHF

!---------------------------------------------------------------------
!		UHF Density Matrix
!---------------------------------------------------------------------
  SUBROUTINE densUHF(nnuc,nelc,atoms,fmem,options)
    IMPLICIT NONE

    ! Values
    ! nnuc	: int, number of nuclii
    ! nelc	: int, number of electrons
    ! options	: 1D int, options of program
    ! fmem	: dp, free memory left in MB
    ! atoms	: 1D int, array of atoms
    
    !Inout
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms,options
    REAL(KIND=8), INTENT(IN) :: fmem
    INTEGER, INTENT(IN) :: nnuc,nelc

    WRITE(*,*) "densUHF called"

  END SUBROUTINE densUHF

!---------------------------------------------------------------------
END PROGRAM dens
