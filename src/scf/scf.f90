!//////////////////////////////////////////////////////////////////
!//           Performs variational Hartree-Fock calculation 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//             
!///////////////////////////////////////////////////////////////////

!=====================================================================
!                       MAIN 

PROGRAM scf
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

  WRITE(*,*) ""
  WRITE(*,*) "          STARTING SCF          "

  CALL getenv(nnuc,nelc,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  !redirect to scf subroutines
  IF (options(2) .EQ. 0) THEN        !RHF
    CALL rhf(nnuc,nelc,atoms,fmem,options)
  ELSE IF (options(2) .EQ. 1) THEN   !UHF
    CALL uhf(nnuc,nelc,atoms,fmem,options)
  ELSE
    WRITE(*,*) "scf: sorry, that method is not implimented yet"
    STOP "bad method in scf"
  END IF
 
  CONTAINS
  
!=====================================================================
!                       SUBROUTINES

!---------------------------------------------------------------------
!               RHF SCF program 
!---------------------------------------------------------------------

  SUBROUTINE rhf(nnuc,nelc,atoms,fmem,options)
    IMPLICIT NONE

    ! Values
    ! nnuc      : int, number of nuclii
    ! nelc      : int, number of electrons
    ! atoms     : 1D int array of atoms
    ! options   : 1D int, options of program
    ! fmem      : dp, free memory left in MB
    ! norb      : int, number of orbitals in molecule
    ! Cui       : 2D dp, molecular orbital coefficients (u'th AO, i'th MO)
    ! Da        : 2D dp, density matrix for alpha
    ! occ       : 1d int, list of occupied orbitals
    ! Mord	: 1d int, list of MO in increasing energy

    !Inout
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms,options
    REAL(KIND=8), INTENT(INOUT) :: fmem
    INTEGER, INTENT(IN) :: nnuc,nelc

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Cui,Huv,Fuv,Guv,Da
    INTEGER, ALLOCATABLE, DIMENSION(:) :: Mord
    INTEGER, DIMENSION(0:1) :: line
    REAL(KIND=8) :: timeS, timeF
    INTEGER :: iter,stat1,stat2,stat3,stat4,stat5,norb
    LOGICAL :: conv, ex  

999 FORMAT(1x,A30,F8.5)
997 FORMAT(1x,A16,F8.5)

    CALL CPU_TIME(timeS)

    iter = 0
    conv = .FALSE.

    !get number of orbitals
    OPEN(unit=1,file='basinfo',status='old',access='sequential')
    READ(1,*) line
    norb = line(1)
    CLOSE(unit=1)
 
    !Allocate memory
    fmem = fmem - (5*norb*norb*8/1.0E6 + norb*4/1.0E6)
    WRITE(*,999) "Allocating memory for rhf (MB) ", (5*norb*norb*8/1.0E6 + norb*4/1.0E6)
    IF (fmem .GT. 0) THEN
      ALLOCATE(Cui(0:norb-1,0:norb-1),STAT=stat1)
      ALLOCATE(Huv(0:norb-1,0:norb-1),STAT=stat2)
      ALLOCATE(Fuv(0:norb-1,0:norb-1),STAT=stat3)
      ALLOCATE(Guv(0:norb-1,0:norb-1),STAT=stat4)
      ALLOCATE(Da(0:norb-1,0:norb-1),STAT=stat5)
      IF (stat1 + stat2 + stat3 + stat4 + stat5 .NE. 0) THEN
        WRITE(*,*) "rhf: couldn't allocate space for matrices"
        CALL EXECUTE_COMMAND_LINE('touch error')
        STOP
      END IF
      CALL nmem(fmem)
    ELSE
      WRITE(*,*) "rhf: max memory reached, exiting"
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP
    END IF

    !check if we need to generate starting point for the calc
    INQUIRE(file='Fuv',exist=ex)
    IF (.NOT. ex) THEN
      CALL initRHF()
    END IF

    ! RHF iterations
    DO WHILE (.NOT. conv)
      iter = iter + 1

      WRITE(*,*) "Starting scf iteration", iter 
!      CALL rhfiter()

      !testing only
      IF (iter .GT. 2) conv = .TRUE.

    END DO

    DEALLOCATE(Cui)
    DEALLOCATE(Fuv)
    DEALLOCATE(Huv)
    DEALLOCATE(Guv)
    DEALLOCATE(Da)

    fmem = fmem - (5*norb*norb*8/1.0E6 + norb*4/1.0E6)
    CALL nmem(fmem)

    CALL CPU_TIME(timeF)
    WRITE(*,997) "rhf ran in (s) :", (timeF - timeS)
    
  END SUBROUTINE rhf

!---------------------------------------------------------------------
!               UHF SCF program 
!---------------------------------------------------------------------
  SUBROUTINE uhf(nnuc,nelc,atoms,fmem,options)
    IMPLICIT NONE

    ! Values
    ! nnuc      : int, number of nuclii
    ! nelc      : int, number of electrons
    ! atoms     : 1D int array of atoms
    ! options   : 1D int, options of program
    ! fmem      : dp, free memory left in MB
    ! norb      : int, number of orbitals in molecule
    ! Cui       : 2D dp, molecular orbital coefficients (u'th AO, i'th MO)
    ! Da        : 2D dp, density matrix for alpha
    ! conv	: bool, false if not converged, true if converged

    !Inout
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms,options
    REAL(KIND=8), INTENT(INOUT) :: fmem
    INTEGER, INTENT(IN) :: nnuc,nelc

    STOP "uhf not implimented yet"

  END SUBROUTINE uhf

!---------------------------------------------------------------------
!		Generate initial RHF Cui and dens 
!---------------------------------------------------------------------
  SUBROUTINE initRHF()
    IMPLICIT NONE

    WRITE(*,*)
    WRITE(*,*) "initRHF called"


  END SUBROUTINE initRHF

!---------------------------------------------------------------------
END PROGRAM scf
