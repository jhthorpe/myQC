!///////////////////////////////////////////////////////////////////
!//              Module for dealing with enviromental data 
!//
!//                     James H Thorpe, in Group of John Stanton
!//                     The University of Florida
!//
!///////////////////////////////////////////////////////////////////

MODULE env
  IMPLICIT NONE

  CONTAINS

!~~~~~~
! Get enviromental setup
  SUBROUTINE getenv(nnuc,nelcA,nelcB,xyz,atoms,fmem,options)
    IMPLICIT NONE

    ! Parameters
    ! nnuc	: int, number of nuclii
    ! nelcA,B	: int, number of electrons in orbitals A,B
    ! xyz	: 2D dp, array of nuclear locations (bohr)
    ! atoms	: 1D int, array of atoms
    ! fmem	: dp, MB of free memory left
    ! options	: 1D int, array of options

    ! Inout
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: xyz
    INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms, options
    REAL(KIND=8), INTENT(INOUT) :: fmem
    INTEGER, INTENT(INOUT) :: nnuc, nelcA,nelcB

    ! Internal
    INTEGER :: i,j,ios

    OPEN(unit=1,file='envdat',status='old',access='sequential',iostat=ios)    
    IF (ios .NE. 0) THEN
      WRITE(*,*) "Couldn't open envdat"
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP "Couldn't open envdat"
    END IF
    OPEN(unit=2,file='nucpos',status='old',access='sequential',iostat=ios) 
    IF (ios .NE. 0) THEN
      WRITE(*,*) "Couldn't open nucpos"
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP "Couldn't open nucpos"
    END IF
    OPEN(unit=3,file='fmem',status='old',access='sequential',iostat=ios)
    IF (ios .NE. 0) THEN
      WRITE(*,*) "Couldn't open fmem"
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP "Couldn't open fmem"
    END IF

    READ(3,*) fmem
    READ(1,*) nnuc
    READ(1,*) nelcA, nelcB 
    READ(1,*) j
    ALLOCATE(options(0:j-1))
    READ(1,*) options
    
    ALLOCATE(atoms(0:nnuc-1))
    
    ALLOCATE(xyz(0:nnuc-1,0:2)) 
    DO i=0,nnuc-1
      READ(2,*) atoms(i),  xyz(i,0), xyz(i,1), xyz(i,2)
    END DO

    CLOSE(unit=3)
    CLOSE(unit=2)
    CLOSE(unit=1)

  END SUBROUTINE getenv

!~~~~~~
! clear enviroment
  SUBROUTINE setenv(atoms, xyz, fmem, options)
  !SUBROUTINE setenv(atoms, xyz, radii, fmem)
    IMPLICIT NONE
   
    ! Inout
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: xyz 
    !REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: radii,xyz 
    INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: atoms,options
    REAL(KIND=8), INTENT(IN) :: fmem
    INTEGER :: ios

    OPEN(unit=1,file='fmem',status='replace',access='sequential',iostat=ios)
    IF (ios .NE. 0) THEN
      WRITE(*,*) "Couldn't open fmem"
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP "Couldn't open fmem"
    END IF
    
    WRITE(1,*) fmem
    CLOSE(unit=1)

    DEALLOCATE(xyz)
    DEALLOCATE(atoms)
    DEALLOCATE(options)
    !DEALLOCATE(radii)

  END SUBROUTINE setenv
!~~~~~~
! update memory  mem 
  SUBROUTINE nmem(fmem)
    IMPLICIT NONE
    REAL(KIND=8),INTENT(IN) :: fmem
    OPEN(unit=1,file='fmem',status='old',access='sequential')
    WRITE(1,*) fmem
    CLOSE(unit=1)  
  END SUBROUTINE nmem
!~~~~~~

END MODULE env
