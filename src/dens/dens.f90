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
    ! Cui	: 2D dp, molecular orbital coefficients (u'th AO, i'th MO)
    ! Da	: 2D dp, density matrix for alpha
    
    !Inout
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms,options
    REAL(KIND=8), INTENT(INOUT) :: fmem
    INTEGER, INTENT(IN) :: nnuc,nelc

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Cui, Da
    INTEGER, DIMENSION(0:1) :: line
    REAL(KIND=8) :: temp 
    INTEGER :: norb,stat1,stat2
    INTEGER :: u,v,i,j 

999 FORMAT(1x,A43,F8.5)

    !get number of orbitals
    OPEN(unit=1,file='basinfo',status='old',access='sequential')
    READ(1,*) line
    norb = line(1)
    CLOSE(unit=1) 

    !Allocate memory
    fmem = fmem - 2*norb*norb*8/1.0D6 - norb*4/1.0D6
    IF (fmem .GT. 0) THEN
      ALLOCATE(Cui(0:norb-1,0:norb-1),STAT=stat1)
      ALLOCATE(Da(0:norb-1,0:norb-1),STAT=stat2)
      IF (stat1 + stat2 .NE. 0) THEN
        WRITE(*,*) "densRHF: couldn't allocate space for Cui and density matrix" 
        CALL EXECUTE_COMMAND_LINE('touch error')
        STOP
      END IF
      CALL nmem(fmem)
    ELSE
      WRITE(*,*) "densRHF: max memory reached, exiting"
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP
    END IF

    !get MOs
    OPEN(unit=1,file='Cui',access='sequential',status='old')
    READ(1,*) Cui(:,:) 
    CLOSE(unit=1)

    !zero density matrix
    DO j=0,norb-1
      Da(:,j) = (/ (0.0D0, i=0,norb-1) /)
    END DO

    !generate density matrix
    !loop over orbitals
    DO v=0,norb-1
      DO u=0,norb-1
        !loop over occupied MO
        temp = 0.0D0 
        DO i=0,nelc/2-1
          temp = temp + Cui(u,i)*Cui(v,i) 
        END DO
        Da(u,v) = temp * 2
      END DO
    END DO

    IF (options(7) .GE. 3) THEN
      WRITE(*,*) "Writing human readable density matrix to dens.txt"
      OPEN(unit=1,file='dens.txt',access='sequential',status='replace')
      WRITE(1,*) Da(:,:)
      CLOSE(unit=1,status='keep')
    END IF

    !write density matrix
    OPEN(unit=2,file='Da',access='sequential',status='replace',form='unformatted')
    WRITE(2) Da
    CLOSE(unit=2,status='keep')

    DEALLOCATE(Cui) 
    DEALLOCATE(Da)

    fmem = fmem + 2*norb*norb*8/1.0D6 + norb*4/1.0D6 
    CALL nmem(fmem)

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
