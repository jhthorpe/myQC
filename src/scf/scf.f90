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
    CALL RHF(nnuc,nelc,atoms,fmem,options)
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
  SUBROUTINE RHF(nnuc,nelc,atoms,fmem,options)
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
    ! LWORK	: int, best size for work arrays, for lapack
    ! Eig	: 1d dp, orbital eigenvalues
    ! Enr	: dp, nuclear repulsion energy
    ! Etrk	: 1d dp, last couple total energies

    !Inout
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms,options
    REAL(KIND=8), INTENT(INOUT) :: fmem
    INTEGER, INTENT(IN) :: nnuc,nelc

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Cui,Suv,Huv,Fuv,Guv,Da
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Eig
    REAL(KIND=8), DIMENSION(1:3) :: Etrk
    INTEGER, DIMENSION(0:1) :: line
    REAL(KIND=8) :: timeS, timeF, Enr
    INTEGER :: LWORK,i
    INTEGER :: iter,stat1,stat2,stat3,stat4,stat5,stat6,stat7,norb
    LOGICAL :: conv,ex,flag 

999 FORMAT(1x,A30,F8.5)
997 FORMAT(1x,A16,F8.5)

    CALL CPU_TIME(timeS)

    iter = 0
    LWORK = -1
    conv = .FALSE.
    Etrk = [0.0D0, 0.0D0, 0.0D0]

    !get number of orbitals
    OPEN(unit=1,file='basinfo',status='old',access='sequential')
    READ(1,*) line
    norb = line(1)
    CLOSE(unit=1)
 
    !check we will have enough memory
    fmem = fmem - (6*norb*norb*8/1.0D6 + norb*8/1.0D6)
    WRITE(*,999) "Allocating memory for rhf (MB) ", (6*norb*norb*8/1.0D6 + norb*8/1.0D6)
    IF (fmem .GT. 0) THEN
      ALLOCATE(Cui(0:norb-1,0:norb-1),STAT=stat1)
      ALLOCATE(Suv(0:norb-1,0:norb-1),STAT=stat2)
      ALLOCATE(Huv(0:norb-1,0:norb-1),STAT=stat3)
      ALLOCATE(Fuv(0:norb-1,0:norb-1),STAT=stat4)
      ALLOCATE(Guv(0:norb-1,0:norb-1),STAT=stat5)
      ALLOCATE(Da(0:norb-1,0:norb-1),STAT=stat6)
      ALLOCATE(Eig(0:norb-1),STAT=stat7)
      IF (stat1+stat2+stat3+stat4+stat5+stat6 .NE. 0) THEN
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

    !get nuclear repulsion
    CALL getEnr(xyz,atoms,Enr)

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
      CALL initRHF(norb,Cui,Suv,Huv,LWORK,fmem,options(7))
      CALL EXECUTE_COMMAND_LINE('dens') 
    ELSE
      !check if we have density matrix 
      INQUIRE(file='Da',exist=flag)
      IF (.NOT. flag) THEN
        CALL EXECUTE_COMMAND_LINE('dens')
      END IF
    END IF

    !check eigenvalues of the overlap matrix
    CALL checkSuv(Suv,norb,fmem,options(7))
    WRITE(*,*)

    ! RHF iterations
    WRITE(*,*) "RHF calculation..."
    WRITE(*,*) "Interation   Total Energy (hartrees)   Ediff"
    WRITE(*,*) "=============================================="
    DO WHILE (.NOT. conv)
      iter = iter + 1
      CALL RHFiter(Suv,Huv,Guv,Fuv,Cui,Da,Eig,Enr,Etrk,norb,conv,fmem,iter,LWORK,options)
    END DO
    WRITE(*,*) "=============================================="

    !Write output
    WRITE(*,*) "Orbital eigenvalues"
    WRITE(*,*) "=============================================="
    WRITE(*,*) 1, Eig(0)
    DO i=1,norb-1
      IF (Eig(i) .GT. 0.0D0 .AND. Eig(i-1) .LT. 0.0D0) THEN
        WRITE(*,*) "----------------------------------------------"
      END IF 
      WRITE(*,*) i+1, Eig(i)
    END DO
    WRITE(*,*) "=============================================="
    WRITE(*,*)

    DEALLOCATE(Cui)
    DEALLOCATE(Suv)
    DEALLOCATE(Fuv)
    DEALLOCATE(Huv)
    DEALLOCATE(Guv)
    DEALLOCATE(Da)
    DEALLOCATE(Eig)

    fmem = fmem + (6*norb*norb*8/1.0D6 + norb*8.0/1.0D6)
    CALL nmem(fmem)

    CALL CPU_TIME(timeF)
    WRITE(*,997) "rhf ran in (s) :", (timeF - timeS)
    
  END SUBROUTINE RHF

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
!WORK NOTE - should symm the overlap matrix beforehand?
!---------------------------------------------------------------------
  SUBROUTINE initRHF(norb,Cui,Suv,Huv,LWORK,fmem,verb)
    USE env
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
    REAL(KIND=8), INTENT(INOUT) :: fmem
    INTEGER, INTENT(INOUT) :: LWORK
    INTEGER, INTENT(IN) :: norb,verb
    
    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: A,B
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: W,WORK
    REAL(KIND=8) :: temp
    INTEGER :: INFO,stat1,stat2,stat3
    INTEGER :: i,j,k,u,v 

    WRITE(*,*)
    WRITE(*,*) "initRHF called"
    WRITE(*,*) "Generating Cui from core Hamiltonian." 

    ! setup arrays
    ALLOCATE(A(0:norb-1,0:norb-1),STAT=stat1)
    ALLOCATE(B(0:norb-1,0:norb-1),STAT=stat2)
    ALLOCATE(W(0:norb-1),STAT=stat3)
    fmem = fmem - (2*norb*norb*8.0/1.0D6 + norb*8.0/1.0D6)
    IF (fmem .LT. 0.0D6 .OR. stat1+stat2+stat3 .NE. 0) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      WRITE(*,*) "scf:initRHF could not allocate enough memory"
      STOP "scf:initRHF out of memory"
    END IF 
    CALL nmem(fmem)
    DO v=0,norb-1
      A(0:v-1,v) = (/ (0.0D0, u=0,v-1) /)
      A(v:norb-1,v) = (/ (Huv(u,v), u=v,norb-1) /)
      B(0:v-1,v) = (/ (0.0D0, u=0,v-1) /)
      B(v:norb-1,v) = (/ (Suv(u,v), u=v,norb-1) /)
    END DO

    !find LWORK
    ALLOCATE(WORK(0:1))
    LWORK=-1
    CALL DSYGV(1,'V','L',norb,A,norb,B,norb,W,WORK,LWORK,INFO)
    LWORK = CEILING(WORK(0))
    DEALLOCATE(WORK)
    ALLOCATE(WORK(0:LWORK-1),STAT=stat1)
    fmem = fmem - LWORK*8/1.0D6 
    IF (fmem .LT. 0.0D0 .OR. stat1 .NE. 0) THEN
      fmem = fmem + LWORK*8/1.0D6
      LWORK=norb
      ALLOCATE(WORK(0:LWORK-1),STAT=stat2)
      fmem = fmem - LWORK*8/1.0D6 
      IF (fmem .LT. 0.0D0 .OR. stat2 .NE. 0) THEN
        CALL EXECUTE_COMMAND_LINE('touch error')
        WRITE(*,*) "scf:initRHF could not allocate memory"
        STOP
      END IF
    END IF

    !find initial Cui and order of orbitals
    CALL DSYGV(1,'V','L',norb,A,norb,B,norb,W,WORK,LWORK,INFO)
    IF (INFO .GT. 0 .AND. INFO .LT. norb) THEN
      WRITE(*,*) "Warning: potential contamination of eigenvalues."
      WRITE(*,*) "Number of contaminents: ", INFO
    ELSE IF (INFO .GT. norb) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      WRITE(*,*) "Leading minor of order i of B is not positive definite. i=",INFO-norb
      STOP
    END IF

    !print if desired
    IF (verb .GE. 3) THEN
      WRITE(*,*) "Eigenvalues of core hamiltonian"
      DO i=0,norb-1
        WRITE(*,*) i, W(i)
      END DO
      WRITE(*,*)
    END IF

    !write to Cui
    OPEN(unit=7,file='Cui',access='sequential',status='replace')
    WRITE(7,*) A(:,:)
    CLOSE(unit=7,status='keep')
    Cui(:,:) = A(:,:)

    !cleanup memory
    DEALLOCATE(A)
    DEALLOCATE(B)
    DEALLOCATE(W)
    DEALLOCATE(WORK)
    fmem = fmem + (2*norb*norb*8.0/1.0D6 + norb*8.0/1.0D6)
    fmem = fmem + (LWORK*8/1.0D6)
    CALL nmem(fmem)

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
    fmem = fmem - (norb*norb*8/1.0D6 + norb*8/1.0D6) 
    IF (fmem .LT. 0.0D0 .OR. stat1+stat2 .NE. 0) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      WRITE(*,*) "scf:checkSuv could not allocate memory"
      STOP
    END IF
    CALL nmem(fmem)

    !put S into A
    DO v=0,norb-1
      A(0:v-1,v) = (/ (0.0D0, u=0,v-1) /)
      A(v:norb-1,v) = (/ (Suv(u,v), u=v,norb-1) /)
    END DO

    !get LWORK
    ALLOCATE(WORK(0:1))
    LWORK=-1
    CALL DSYEV('N','L',norb,A,norb,W,WORK,LWORK,INFO)
    LWORK = CEILING(WORK(0))
    DEALLOCATE(WORK)
    ALLOCATE(WORK(0:LWORK-1),STAT=stat1)
    fmem = fmem - LWORK*8/1.0D6 
    IF (fmem .LT. 0.0D0 .OR. stat1 .NE. 0) THEN
      fmem = fmem + LWORK*8/1.0D6
      LWORK=norb
      ALLOCATE(WORK(0:LWORK-1),STAT=stat2)
      fmem = fmem - LWORK*8/1.0D6 
      IF (fmem .LT. 0.0D0 .OR. stat2 .NE. 0) THEN
        CALL EXECUTE_COMMAND_LINE('touch error')
        WRITE(*,*) "scf:checkSuv could not allocate memory"
        STOP
      END IF
    END IF
    
    !get eigenvalues
    CALL DSYEV('N','L',norb,A,norb,W,WORK,LWORK,INFO)
    IF (INFO .GT. 0) THEN
      WRITE(*,*) "Warning: potential contamination of eigenvalues."
      WRITE(*,*) "Number of contaminents: ", INFO
    END IF

    !print eigenvalues
    DO i=0,norb-1
      IF (W(i) .LT. 1.0E-4) THEN
        WRITE(*,997) i+1, W(i), "WARNING"
      ELSE IF (verb .GE. 1) THEN
        WRITE(*,999) i+1, W(i)
      END IF
    END DO
    
    DEALLOCATE(A)
    DEALLOCATE(W)
    DEALLOCATE(WORK)

    fmem = fmem + (norb*norb*8/1.0D6 + norb*8/1.0D6) 
    fmem = fmem + (LWORK*8/1.0D6)
    CALL nmem(fmem)

  END SUBROUTINE checkSuv

!---------------------------------------------------------------------
!			RHF iteration	
!---------------------------------------------------------------------
  SUBROUTINE RHFiter(Suv,Huv,Guv,Fuv,Cui,Da,Eig,Enr,Etrk,norb,conv,fmem,iter,LWORK,options)
    IMPLICIT NONE
    
    !Values
    !norb	: int, number of orbitals
    !conv	: bool, have we converged?
    !fmem	: dp, memory remaining, MB
    !options	: 1D int, list of input options
    !Suv	: 2D dp, overlap matrix of AO
    !Huv	: 2D dp, core hamiltonian matrix of AO
    !Guv	: 2D dp, 2e- hamiltonian matrix of AO
    !Fuv	: 2D dp, Fock matrix AO 
    !Cui	: 2D dp, coefficients of MO
    !Da		: 2D dp, density matrix of AO 
    !Eig	: 1D dp, vector of orbital eigenvalues
    !Enr	: dp, nuclear repulsion energy
    !Etrk	: 1D dp, list of latest energies
    !iter	: int, iteration number
    !Eelc	: dp, electronic energy
    !LWORK	: int, length of work array

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Guv,Fuv,Cui,Da
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Suv,Huv
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Eig,Etrk
    INTEGER, DIMENSION(0:), INTENT(IN) :: options
    REAL(KIND=8), INTENT(INOUT) :: fmem
    REAL(KIND=8), INTENT(IN) :: Enr
    INTEGER, INTENT(INOUT) :: LWORK
    LOGICAL, INTENT(INOUT) :: conv 
    INTEGER, INTENT(IN) :: norb,iter

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: A,B
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: WORK
    REAL(KIND=8) :: Eelc,temp
    INTEGER :: INFO,stat1,stat2,stat3
    INTEGER :: i,j,u,v
  
    !get Density matrix  
    CALL EXECUTE_COMMAND_LINE('dens')
    OPEN(unit=8,file='Da',access='sequential',status='old',form='unformatted')
    READ(8) Da(:,:) 
    CLOSE(unit=8)

    !get Guv matrix
    CALL EXECUTE_COMMAND_LINE('RHFI2G')
    OPEN(unit=11,file='Guv',access='sequential',status='old',form='unformatted')
    READ(11) Guv(:,:)
    CLOSE(unit=11) 

    !Create Fock matrix
    Fuv = Huv + Guv

    !Get Electronic energy
    Eelc = 0.0D0
    DO v=0,norb-1
      DO u=0,norb-1
        Eelc = Eelc + Da(u,v)*(Fuv(u,v)+Huv(u,v))
      END DO
    END DO
    Eelc = 0.5D0*Eelc
 
    !write output
    !WORK NOTe - need better output
    WRITE(*,*) iter, Eelc+Enr, Eelc+Enr-Etrk(2)

    !Get our coefficients and eigenvalues

    ! setup arrays
    ALLOCATE(A(0:norb-1,0:norb-1),STAT=stat1)
    ALLOCATE(B(0:norb-1,0:norb-1),STAT=stat2)
    fmem = fmem - (2*norb*norb*8.0/1.0D6) 
    IF (fmem .LT. 0.0D6 .OR. stat1+stat2+stat3 .NE. 0) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      WRITE(*,*) "scf:RHFiter could not allocate enough memory"
      STOP "scf:RHFiter out of memory"
    END IF
    CALL nmem(fmem)
    DO v=0,norb-1
      A(0:v-1,v) = (/ (0.0D0, u=0,v-1) /)
      A(v:norb-1,v) = (/ (Fuv(u,v), u=v,norb-1) /)
      B(0:v-1,v) = (/ (0.0D0, u=0,v-1) /)
      B(v:norb-1,v) = (/ (Suv(u,v), u=v,norb-1) /)
    END DO

    !if first iteration...
    IF (iter .EQ. 0) THEN
      ALLOCATE(WORK(0:1))
      LWORK=-1
      CALL DSYGV(1,'V','L',norb,A,norb,B,norb,Eig,WORK,LWORK,INFO)
      LWORK = CEILING(WORK(0))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(0:LWORK-1),STAT=stat1)
      fmem = fmem - LWORK*8/1.0D6
      IF (fmem .LT. 0.0D0 .OR. stat1 .NE. 0) THEN
        fmem = fmem + LWORK*8/1.0D6
        LWORK=norb
        ALLOCATE(WORK(0:LWORK-1),STAT=stat2)
        fmem = fmem - LWORK*8/1.0D6
        IF (fmem .LT. 0.0D0 .OR. stat2 .NE. 0) THEN
          CALL EXECUTE_COMMAND_LINE('touch error')
          WRITE(*,*) "scf:initRHF could not allocate memory"
          STOP
        END IF
      END IF

    !if not...
    ELSE
      ALLOCATE(WORK(0:LWORK-1),STAT=stat1)
      fmem = fmem - LWORK*8/1.0D6
      IF (fmem .LT. 0.0D0 .OR. stat1 .NE. 0) THEN
        fmem = fmem + LWORK*8/1.0D6
        LWORK=norb
        ALLOCATE(WORK(0:LWORK-1),STAT=stat2)
        fmem = fmem - LWORK*8/1.0D6
        IF (fmem .LT. 0.0D0 .OR. stat2 .NE. 0) THEN
          CALL EXECUTE_COMMAND_LINE('touch error')
          WRITE(*,*) "scf:initRHF could not allocate memory"
          STOP
        END IF
      END IF
    END IF

    !find initial Cui and order of orbitals
    CALL DSYGV(1,'V','L',norb,A,norb,B,norb,Eig,WORK,LWORK,INFO)
    IF (INFO .GT. 0 .AND. INFO .LT. norb) THEN
      WRITE(*,*) "Warning: potential contamination of eigenvalues."
      WRITE(*,*) "Number of contaminents: ", INFO
    ELSE IF (INFO .GT. norb) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      WRITE(*,*) "Leading minor of order i of B is not positive definite. i=",INFO-norb
      STOP
    END IF
 
    !print if desired
    IF (options(7) .GE. 3) THEN
      WRITE(*,*) "Orbital eigenvalues..."
      DO i=0,norb-1
        WRITE(*,*) i, Eig(i)
      END DO
      WRITE(*,*) "---------------------------------------"
    END IF
 
    !write to Cui and Eig
    OPEN(unit=7,file='Cui',access='sequential',status='replace')
    WRITE(7,*) A(:,:)
    CLOSE(unit=7,status='keep')
    Cui(:,:) = A(:,:)
 
    !cleanup memory
    DEALLOCATE(A)
    DEALLOCATE(B)
    DEALLOCATE(WORK)
    fmem = fmem + (2*norb*norb*8.0/1.0D6)
    fmem = fmem + (LWORK*8/1.0D6)
    CALL nmem(fmem)

    !check for convergece
    Etrk(0) = Etrk(1)
    Etrk(1) = Etrk(2)
    Etrk(2) = Eelc + Enr


    IF (iter .GT. 2) THEN
      CALL checkConv(Etrk,options,conv)
    END IF

    IF (iter .GT. 300) THEN
      WRITE(*,*) "SCF failed to converge"
      WRITE(*,*)  
    END IF 

    IF (conv) THEN
      WRITE(*,*) "SCF has converged!"
      WRITE(*,*) 
      WRITE(*,*) "Total SCF energy (a.u)", Eelc+Enr
      WRITE(*,*)
    END  IF


  END SUBROUTINE RHFiter

!---------------------------------------------------------------------
!		Calculate nuclear repulsion energy			
!---------------------------------------------------------------------
  SUBROUTINE getEnr(xyz,atoms,Enr)
    IMPLICIT NONE

    !Values
    ! xyz	: 2D dp, list of atom locations
    ! atoms	: 1D int, list of atom id
    ! Enr	: dp, nuclear repulsion energy
    ! M		: int, number of atoms

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms
    REAL(KIND=8), INTENT(INOUT) :: Enr

    !Internal
    REAL(KIND=8) :: RAB
    INTEGER :: A,B,M
   
    M = SIZE(atoms) 
    Enr = 0.0D0

    DO A=0,M-1
      DO B=0,A-1
        RAB = SQRT((xyz(A,0)-xyz(B,0))**2.0D0+(xyz(A,1)-xyz(B,1))**2.0D0+(xyz(A,2)-xyz(B,2))**2.0D0)
        Enr = Enr + atoms(A)*atoms(B)/RAB
      END DO
    END DO 

    WRITE(*,*) "Nuclear repulsion energy (hartrees) : ", Enr

  END SUBROUTINE

!---------------------------------------------------------------------
!		Check for convergence	
!---------------------------------------------------------------------
  SUBROUTINE checkConv(Etrk,options,conv)
    IMPLICIT NONE

    !Value
    !Etrk	: 1D dp, list of 3 last energies
    !options	: 1D int, options array
    !conv	: bool, converged or not

    !Inout
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Etrk
    INTEGER, DIMENSION(0:), INTENT(IN) :: options
    LOGICAL, INTENT(INOUT) :: conv

    !Internal
    REAL(KIND=8), DIMENSION(0:11) :: lim
    INTEGER :: i 

    !These go up to 11
    lim = [1.0D0,1.0D-1,1.0D-2,1.0D-3,1.0D-4,1.0D-5,1.0D-6,1.0D-7,1.0D-8,1.0D-9,1.0D-10,1.0D-11]
    
    i = options(8)

!    IF (Etrk(0)-Etrk(1) .LT. lim(i) .AND. Etrk(1)-Etrk(2) .LT. lim(i)) THEN
    IF (Etrk(1)-Etrk(2) .LT. lim(i)) THEN
      conv = .TRUE.
    END IF 

  END SUBROUTINE
!---------------------------------------------------------------------
END PROGRAM scf
