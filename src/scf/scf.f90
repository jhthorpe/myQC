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
    ! LWORK	: int, best size for work arrays, for lapack

    !Inout
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms,options
    REAL(KIND=8), INTENT(INOUT) :: fmem
    INTEGER, INTENT(IN) :: nnuc,nelc

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Cui,Suv,Huv,Fuv,Guv,Da
    INTEGER, ALLOCATABLE, DIMENSION(:) :: Mord
    INTEGER, DIMENSION(0:1) :: line
    REAL(KIND=8) :: timeS, timeF
    INTEGER :: LWORK
    INTEGER :: iter,stat1,stat2,stat3,stat4,stat5,stat6,stat7,norb
    LOGICAL :: conv,ex,flag 

999 FORMAT(1x,A30,F8.5)
997 FORMAT(1x,A16,F8.5)

    CALL CPU_TIME(timeS)

    iter = 0
    LWORK = -1
    conv = .FALSE.

    !get number of orbitals
    OPEN(unit=1,file='basinfo',status='old',access='sequential')
    READ(1,*) line
    norb = line(1)
    CLOSE(unit=1)
 
    !Allocate memory
    fmem = fmem - (6*norb*norb*8/1.0E6 + norb*4/1.0E6)
    WRITE(*,999) "Allocating memory for rhf (MB) ", (6*norb*norb*8/1.0E6 + norb*4/1.0E6)
    IF (fmem .GT. 0) THEN
      ALLOCATE(Cui(0:norb-1,0:norb-1),STAT=stat1)
      ALLOCATE(Suv(0:norb-1,0:norb-1),STAT=stat2)
      ALLOCATE(Huv(0:norb-1,0:norb-1),STAT=stat3)
      ALLOCATE(Fuv(0:norb-1,0:norb-1),STAT=stat4)
      ALLOCATE(Guv(0:norb-1,0:norb-1),STAT=stat5)
      ALLOCATE(Da(0:norb-1,0:norb-1),STAT=stat6)
      ALLOCATE(Mord(0:norb-1),STAT=stat7)
      IF (stat1+stat2+stat3+stat4+stat5+stat6+stat7 .NE. 0) THEN
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

    !get data
    OPEN(unit=2,file='Huv',status='old',access='sequential')
    READ(2,*) Huv(:,:)
    CLOSE(unit=2)
    OPEN(unit=3,file='Suv',status='old',access='sequential')
    READ(3,*) Suv(:,:)
    CLOSE(unit=3)

    !check if we have coefficients 
    INQUIRE(file='Cui',exist=ex)
    IF (.NOT. ex) THEN
      CALL initRHF(norb,Cui,Suv,Huv,LWORK)
!      CALL EXECUTE_COMMAND_LINE('dens') 
    ELSE
      OPEN(unit=1,file='Cui',status='old',access='sequential')
      READ(1,*) Cui(:,:)
      CLOSE(unit=1)
      
      !check if we have density matrix 
      INQUIRE(file='Da',exist=flag)
      IF (.NOT. flag) THEN
!        CALL EXECUTE_COMMAND_LINE('dens')
      END IF
    END IF

    !read in density matrix
!    OPEN(unit=5,file='Da',access='sequential',status='old',form='unformatted')
!    READ(5) Da(:,:)
!    CLOSE(unit=5)

    !check eigenvalues of the overlap matrix
    CALL checkSuv(Suv,norb,fmem,options(7))

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

    fmem = fmem + (6*norb*norb*8/1.0E6 + norb*4/1.0E6)
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
!		Generate initial RHF MO coeffs 
!---------------------------------------------------------------------
  SUBROUTINE initRHF(norb,Cui,Suv,Huv,LWORK)
    IMPLICIT NONE

    !Values
    ! norb	: int, number of orbitals
    ! Cui	: 2D dp, array of MO coefficients
    ! Suv	: 2D dp, array of AO overlaps
    ! Huv	: 2D dp, array of 1e- AO integrals
    ! Da	: 2D dp, array of AO densities
   
    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Cui
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Suv, Huv
    INTEGER, INTENT(INOUT) :: LWORK
    INTEGER, INTENT(IN) :: norb
    
    !Internal
    INTEGER :: i,j,k,u,v 

    WRITE(*,*)
    WRITE(*,*) "initRHF called"
    WRITE(*,*) "Generating Cui from core Hamiltonian." 

    !use core hamiltonian to get initial guess at coefficients

    !

  END SUBROUTINE initRHF

!---------------------------------------------------------------------
!		Check eigenvalues of the overlap matrix	
!---------------------------------------------------------------------
  SUBROUTINE checkSuv(Suv,norb,fmem,verb)
    IMPLICIT NONE
  
    !Values
    !Suv	: 2D dp, array of overlap of orbitals
    !norb	: int, number of orbitals
    !fmem	: dp, free memory in MB
    !verb	: int, verbosity level

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Suv
    REAL(KIND=8), INTENT(INOUT) :: fmem
    INTEGER, INTENT(IN) :: norb,verb

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: A
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: W,WORK
    INTEGER :: INFO,LWORK,stat1,stat2
    INTEGER :: i,j,u,v

999 FORMAT(1x,I3,2x,E15.8)
997 FORMAT(1x,I3,2x,E15.5,2x,A7)
    WRITE(*,*)
    IF (verb .GE. 1) THEN
      WRITE(*,*) "Eigenvalues of the overlap integrals"
    ELSE
      WRITE(*,*) "Checking for dangerous overlap eigenvalues"
    END IF

    ALLOCATE(A(0:norb-1,0:norb-1),STAT=stat1)
    ALLOCATE(W(0:norb-1),STAT=stat2)
    fmem = fmem - (norb*norb*8/1.0E6 + norb*8/1.0E6) 
    IF (fmem .LT. 0.0D0 .OR. stat1+stat2 .NE. 0) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      WRITE(*,*) "scf:checkSuv could not allocate memory"
      STOP
    END IF
    CALL nmem(fmem)

    !put S into A
    DO v=0,norb-1
      DO u=0,v-1
        A(u,v) = 0.0E0
      END DO 
      DO u=v,norb-1
        A(u,v) = Suv(u,v)
      END DO
    END DO

    !get LWORK
    ALLOCATE(WORK(0:1))
    LWORK=-1
    CALL DSYEV('N','L',norb,A,norb,W,WORK,LWORK,INFO)
    LWORK = CEILING(WORK(0))
    DEALLOCATE(WORK)
    ALLOCATE(WORK(0:LWORK-1),STAT=stat1)
    fmem = fmem - LWORK*8/1.0E6 
    IF (fmem .LT. 0.0E6 .OR. stat1 .NE. 0) THEN
      fmem = fmem + LWORK*8/1.0E6
      LWORK=norb
      ALLOCATE(WORK(0:LWORK-1),STAT=stat2)
      fmem = fmem - LWORK*8/1.0E6 
      IF (fmem .LT. 0.0E6 .OR. stat2 .NE. 0) THEN
        CALL EXECUTE_COMMAND_LINE('touch error')
        WRITE(*,*) "scf:checkSuv could not allocate memory"
        STOP
      END IF
    END IF
    
    !get eigenvalues
    CALL DSYEV('N','U',norb,A,norb,W,WORK,LWORK,INFO)
    IF (INFO .GT. 0) THEN
      WRITE(*,*) "Warning: potential contamination of eigenvalues."
      WRITE(*,*) "Number of contaminents: ", INFO
    END IF

    !print eigenvalues
    DO i=0,norb-1
      IF (W(i) .LT. 1.0E-4) THEN
        WRITE(*,997) i, W(i), "WARNING"
      ELSE IF (verb .GE. 1) THEN
        WRITE(*,999) i, W(i)
      END IF
    END DO
    
    DEALLOCATE(A)
    DEALLOCATE(W)
    DEALLOCATE(WORK)

    fmem = fmem + (norb*norb*8/1.0E6 + norb*8/1.0E6) 
    fmem = fmem + (LWORK*8/1.0E6)
    CALL nmem(fmem)

  END SUBROUTINE checkSuv
!---------------------------------------------------------------------
END PROGRAM scf
