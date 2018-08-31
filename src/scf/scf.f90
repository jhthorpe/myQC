!//////////////////////////////////////////////////////////////////
!//           Performs variational Hartree-Fock calculation 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//             
!//  WORK NOTE - it should be noted that the general eigenvalue call
!//    I use to LAPACK produces Cui (and therefor densities) that
!//    are NOT in the orthogonal basis (which is good), but ARE
!//    normalized to Sum(u,v) Cui(u,i)*S(u,v)*Cui(v,i) = 1
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
  ! nelcA,nelcB : int, number of electrons of spin A,B
  ! options     : 1D int, array of options

  ! Internal
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms, options
  REAL(KIND=8) :: fmem, timeS, timeF
  INTEGER :: nnuc, nelcA, nelcB, dummy
  LOGICAL :: flag

  WRITE(*,*) ""
  WRITE(*,*) "                      STARTING SCF"
  !the length is 60
  WRITE(*,*) "------------------------------------------------------------"
  CALL getenv(nnuc,nelcA,nelcB,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  !redirect to scf subroutines
  IF (options(3) .EQ. 0) THEN        !RHF
    CALL RHF(nnuc,nelcA+nelcB,atoms,fmem,options)
  ELSE IF (options(3) .EQ. 1) THEN   !UHF
    CALL UHF(nnuc,nelcA,nelcB,atoms,fmem,options)
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
    REAL(KIND=8) :: timeS,timeF,Enr,mdiff
    INTEGER :: LWORK,i
    INTEGER :: iter,stat1,stat2,stat3,stat4,stat5,stat6,stat7,norb
    LOGICAL :: conv,ex,flag 

999 FORMAT(1x,A30,F8.5)
996 FORMAT(1x,A30,F8.5,F8.5)
997 FORMAT(1x,A16,F8.5)

    CALL CPU_TIME(timeS)

    WRITE(*,*) "rhf called"

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
    ELSE
      WRITE(*,*) "Reading MO coef from Cui"
      OPEN(unit=7,file='Cui',status='old',access='sequential')
      READ(7,*) Cui(:,:)
      CLOSE(unit=7)
    END IF

    CALL EXECUTE_COMMAND_LINE('dens')
    CALL EXECUTE_COMMAND_LINE('RHFI2G')
   
    !check eigenvalues of the overlap matrix
    CALL checkSuv(Suv,norb,fmem,options(7))
    WRITE(*,*)

    ! RHF iterations
    WRITE(*,*) "RHF calculation..."
    WRITE(*,*) "Iteration   Total Energy (hartrees)   Ediff               Ddiff"
    WRITE(*,*) "========================================================================="
    DO WHILE (.NOT. conv)
      CALL RHFiter(Suv,Huv,Guv,Fuv,Cui,Da,Eig,Enr,Etrk,norb,conv,fmem,iter,LWORK,mdiff,options)
      iter = iter + 1
      IF (iter .GE. 500) THEN
        WRITE(*,*) "SCF failed to converge."
        EXIT
      END IF
    END DO
    WRITE(*,*) "========================================================================="

998 FORMAT(15x,I3,4x,F20.15)
    !Write output
    WRITE(*,*) "Orbital eigenvalues (a.u.)"
    WRITE(*,*) "========================================================================="
    WRITE(*,998) 1, Eig(0)
    DO i=1,norb-1
      IF (Eig(i) .GT. 0.0D0 .AND. Eig(i-1) .LT. 0.0D0) THEN
        WRITE(*,*) "-------------------------------------------------------------------------"
      END IF 
      WRITE(*,998) i+1, Eig(i)
    END DO
    WRITE(*,*) "=========================================================================="
    WRITE(*,*)

    !write to molden file
    CALL makeMOLDEN(atoms,xyz,Eig,[0.0D0],0,nnuc,norb,Cui,Cui)

    DEALLOCATE(Cui)
    DEALLOCATE(Suv)
    DEALLOCATE(Fuv)
    DEALLOCATE(Huv)
    DEALLOCATE(Guv)
    DEALLOCATE(Da)
    DEALLOCATE(Eig)
    CALL EXECUTE_COMMAND_LINE('rm Dold')

    fmem = fmem + (6*norb*norb*8/1.0D6 + norb*8.0/1.0D6)
    CALL nmem(fmem)

    CALL CPU_TIME(timeF)
    WRITE(*,997) "rhf ran in (s) :", (timeF - timeS)
    
  END SUBROUTINE RHF

!---------------------------------------------------------------------
!               UHF SCF program 
!---------------------------------------------------------------------
  SUBROUTINE UHF(nnuc,nelcA,nelcB,atoms,fmem,options)
    IMPLICIT NONE

    ! Values
    ! nnuc      : int, number of nuclii
    ! nelcA,B   : int, number of electrons of spin A,B
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
    INTEGER, INTENT(IN) :: nnuc,nelcA,nelcB

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: CuiA,CuiB,Suv,Huv,FuvA,FuvB
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: GuvA,GuvB,Da,Db
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: EigA,EigB
    REAL(KIND=8), DIMENSION(1:3) :: Etrk
    INTEGER, DIMENSION(0:13) :: stat
    INTEGER, DIMENSION(0:1) :: line
    REAL(KIND=8) :: timeS,timeF,Enr,mdiff
    INTEGER :: LWORK,i
    INTEGER :: iter,norb
    LOGICAL :: conv,ex,flag 

999 FORMAT(1x,A30,F8.5)
997 FORMAT(1x,A16,F8.5)

    CALL CPU_TIME(timeS)

    WRITE(*,*) "uhf called"

    iter = 0
    LWORK = -1
    conv = .FALSE.
    Etrk = [0.0D0, 0.0D0, 0.0D0]
    stat = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    !get number of orbitals
    OPEN(unit=1,file='basinfo',status='old',access='sequential')
    READ(1,*) line
    norb = line(1)
    CLOSE(unit=1)
 
    !check we will have enough memory
    fmem = fmem - (10*norb*norb*8/1.0D6 + 2.0D0*norb*8/1.0D6)
    WRITE(*,999) "Allocating memory for uhf (MB) ", (10*norb*norb*8/1.0D6 + 2*norb*8/1.0D6)
    IF (fmem .GT. 0) THEN
      ALLOCATE(CuiA(0:norb-1,0:norb-1),STAT=stat(0))
      ALLOCATE(CuiB(0:norb-1,0:norb-1),STAT=stat(1))
      ALLOCATE(Suv(0:norb-1,0:norb-1),STAT=stat(2))
      ALLOCATE(Huv(0:norb-1,0:norb-1),STAT=stat(3))
      ALLOCATE(FuvA(0:norb-1,0:norb-1),STAT=stat(4))
      ALLOCATE(FuvB(0:norb-1,0:norb-1),STAT=stat(5))
      ALLOCATE(GuvA(0:norb-1,0:norb-1),STAT=stat(6))
      ALLOCATE(GuvB(0:norb-1,0:norb-1),STAT=stat(7))
      ALLOCATE(Da(0:norb-1,0:norb-1),STAT=stat(8))
      ALLOCATE(Db(0:norb-1,0:norb-1),STAT=stat(9))
      ALLOCATE(EigA(0:norb-1),STAT=stat(10))
      ALLOCATE(EigB(0:norb-1),STAT=stat(11))
      IF (SUM(stat(0:11)) .NE. 0) THEN
        WRITE(*,*) "uhf: couldn't allocate space for matrices"
        CALL EXECUTE_COMMAND_LINE('touch error')
        STOP
      END IF
      CALL nmem(fmem)
    ELSE
      WRITE(*,*) "uhf: max memory reached, exiting"
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
      CALL initUHF(norb,CuiA,CuiB,Suv,Huv,LWORK,fmem,options(7))
    ELSE
      WRITE(*,*) "Reading MO coef from Cui"
      OPEN(unit=7,file='Cui',status='old',access='sequential')
      READ(7,*) CuiA(:,:)
      READ(7,*) CuiB(:,:)
      CLOSE(unit=7)
    END IF

    CALL EXECUTE_COMMAND_LINE('dens')
    CALL EXECUTE_COMMAND_LINE('UHFI2G')
   
    !check eigenvalues of the overlap matrix
    CALL checkSuv(Suv,norb,fmem,options(7))
    WRITE(*,*)

    ! UHF iterations
    WRITE(*,*) "UHF calculation..."
    WRITE(*,*) "Iteration   Total Energy (hartrees)   Ediff               Ddiff"
    WRITE(*,*) "========================================================================="
    DO WHILE (.NOT. conv)
      CALL UHFiter(Suv,Huv,GuvA,GuvB,FuvA,FuvB,CuiA,CuiB,Da,Db,EigA,EigB,&
      Enr,Etrk,norb,conv,fmem,iter,LWORK,mdiff,options)
      iter = iter + 1
      IF (iter .GE. 300) THEN
        WRITE(*,*) "SCF failed to converge."
        EXIT
      END IF
    END DO
  
    WRITE(*,*)

    !Write orbitals
998 FORMAT(15x,I3,4x,F20.15)
    WRITE(*,*) "========================================================================="
    WRITE(*,*) "Orbital eigenvalues of Alpha (a.u.)"
    WRITE(*,*) "========================================================================="
    WRITE(*,998) 1, EigA(0)
    DO i=1,norb-1
      IF (EigA(i) .GT. 0.0D0 .AND. EigA(i-1) .LT. 0.0D0) THEN
        WRITE(*,*) "-------------------------------------------------------------------------"
      END IF 
      WRITE(*,998) i+1, EigA(i)
    END DO
    WRITE(*,*) "========================================================================="
    WRITE(*,*)
    WRITE(*,*) "========================================================================="
    WRITE(*,*) "Orbital eigenvalues of Beta (a.u.)"
    WRITE(*,*) "========================================================================="
    WRITE(*,998) 1, EigB(0)
    DO i=1,norb-1
      IF (EigB(i) .GT. 0.0D0 .AND. EigB(i-1) .LT. 0.0D0) THEN
        WRITE(*,*) "-------------------------------------------------------------------------"
      END IF 
      WRITE(*,998) i+1, EigB(i)
    END DO
    WRITE(*,*) "========================================================================="

    !write to molden file
    CALL makeMOLDEN(atoms,xyz,EigA,EigB,1,nnuc,norb,CuiA,CuiB)

    !free the memory we are done with
    DEALLOCATE(FuvA)
    DEALLOCATE(FuvB)
    DEALLOCATE(Huv)
    DEALLOCATE(GuvA)
    DEALLOCATE(GuvB)
    DEALLOCATE(Da)
    DEALLOCATE(Db)
    DEALLOCATE(EigA)
    DEALLOCATE(EigB)

    WRITE(*,*)
    CALL getSpinCont(norb,nelcA,nelcB,CuiA,CuiB,Suv)

    !free last bits of memory
    DEALLOCATE(CuiA)
    DEALLOCATE(CuiB)
    DEALLOCATE(Suv)
    fmem = fmem + (10*norb*norb*8/1.0D6 + 2*norb*8.0/1.0D6)
    CALL nmem(fmem)

    CALL CPU_TIME(timeF)
    WRITE(*,997) "uhf ran in (s) :", (timeF - timeS)
    
  END SUBROUTINE UHF

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
!		Generate initial UHF MO coeffs 
!WORK NOTE - should symm the overlap matrix beforehand?
!---------------------------------------------------------------------
  SUBROUTINE initUHF(norb,CuiA,CuiB,Suv,Huv,LWORK,fmem,verb)
    USE env
    IMPLICIT NONE

    !Values
    ! norb	: int, number of orbitals
    ! Cui	: 2D dp, array of MO coefficients
    ! Suv	: 2D dp, array of AO overlaps
    ! Huv	: 2D dp, array of 1e- AO integrals
    ! Da	: 2D dp, array of AO densities
   
    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: CuiA,CuiB
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
    WRITE(*,*) "initUHF called"
    WRITE(*,*) "Generating Cui from core Hamiltonian." 

    ! setup arrays
    ALLOCATE(A(0:norb-1,0:norb-1),STAT=stat1)
    ALLOCATE(B(0:norb-1,0:norb-1),STAT=stat2)
    ALLOCATE(W(0:norb-1),STAT=stat3)
    fmem = fmem - (2*norb*norb*8.0/1.0D6 + norb*8.0/1.0D6)
    IF (fmem .LT. 0.0D6 .OR. stat1+stat2+stat3 .NE. 0) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      WRITE(*,*) "scf:initUHF could not allocate enough memory"
      STOP "scf:initUHF out of memory"
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
    WRITE(7,*) A(:,:)
    CLOSE(unit=7,status='keep')
    CuiA(:,:) = A(:,:)
    CuiB(:,:) = A(:,:)

    !cleanup memory
    DEALLOCATE(A)
    DEALLOCATE(B)
    DEALLOCATE(W)
    DEALLOCATE(WORK)
    fmem = fmem + (2*norb*norb*8.0/1.0D6 + norb*8.0/1.0D6)
    fmem = fmem + (LWORK*8/1.0D6)
    CALL nmem(fmem)

  END SUBROUTINE initUHF
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

999 FORMAT(1x,I3,2x,ES15.8)
997 FORMAT(1x,I3,2x,ES15.5,2x,A7)
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
  SUBROUTINE RHFiter(Suv,Huv,Guv,Fuv,Cui,Da,Eig,Enr,Etrk,norb,conv,fmem,iter,LWORK,mdiff,options)
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
    !mdiff	: real8, max diff of density matrix

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Guv,Fuv,Cui,Da
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Suv,Huv
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Eig,Etrk
    INTEGER, DIMENSION(0:), INTENT(IN) :: options
    REAL(KIND=8), INTENT(INOUT) :: fmem,mdiff
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
    LOGICAL :: ex
  
    !write Density matrix  
    OPEN(unit=8,file='Da',access='sequential',status='old',form='unformatted')
    READ(8) Da(:,:) 
    CLOSE(unit=8)

    !get Guv matrix
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
 
999 FORMAT(4x,I3,4x,F20.15,8x,ES15.8)
996 FORMAT(4x,I3,4x,F20.15,8x,ES15.8,4x,ES15.8)

    !write output
    !WORK NOTE - need better output
    WRITE(*,996) iter+1, Eelc+Enr, Eelc+Enr-Etrk(2), mdiff

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

    !find Cui 
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
      WRITE(*,*) "Orbital eigenvalues (a.u.)"
      DO i=0,norb-1
        WRITE(*,*) i, Eig(i)
      END DO
      WRITE(*,*) "---------------------------------------"
    END IF
 
    !write MO 
    OPEN(unit=7,file='Cui',access='sequential',status='replace')
    WRITE(7,*) A(:,:)
    CLOSE(unit=7,status='keep')
    Cui(:,:) = A(:,:)

    !write new dens, write old dens, and Guv
    OPEN(unit=9,file='Dold',status='replace',form='unformatted')
    WRITE(9) Da(:,:)
    CLOSE(unit=9)
    CALL EXECUTE_COMMAND_LINE('dens')
    CALL EXECUTE_COMMAND_LINE('RHFI2G')
 
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

    IF (iter .GT. 0) THEN
      CALL checkConv(norb,Etrk,options,conv,mdiff)
    END IF

    IF (conv) THEN
      WRITE(*,*) "SCF has converged!"
      WRITE(*,*) 
      WRITE(*,*) "Total SCF energy (a.u)", Eelc+Enr
      WRITE(*,*)
    END  IF

  END SUBROUTINE RHFiter

!---------------------------------------------------------------------
!			UHF iteration	
!---------------------------------------------------------------------
  SUBROUTINE UHFiter(Suv,Huv,GuvA,GuvB,FuvA,FuvB,CuiA,CuiB,Da,Db,EigA,EigB,&
  Enr,Etrk,norb,conv,fmem,iter,LWORK,mdiff,options)
    IMPLICIT NONE
    
    !Values
    !norb	: int, number of orbitals
    !conv	: bool, have we converged?
    !fmem	: dp, memory remaining, MB
    !options	: 1D int, list of input options
    !Suv	: 2D dp, overlap matrix of AO
    !Huv	: 2D dp, core hamiltonian matrix of AO
    !GuvA,B	: 2D dp, 2e- hamiltonian matrix of AO spin A,B
    !FuvA,B	: 2D dp, Fock matrix AO spin A,B 
    !CuiA,B	: 2D dp, coefficients of MO spin A,B
    !Da,b	: 2D dp, density matrix of AO spin A,B
    !EigA,B	: 1D dp, vector of orbital eigenvalues spin A,B
    !Enr	: dp, nuclear repulsion energy
    !Etrk	: 1D dp, list of latest energies
    !iter	: int, iteration number
    !Eelc	: dp, electronic energy
    !LWORK	: int, length of work array
    !mdiff	: real8, max diff of density matrix

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: GuvA,GuvB,FuvA,FuvB,CuiA,CuiB,Da,Db
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Suv,Huv
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: EigA,EigB,Etrk
    INTEGER, DIMENSION(0:), INTENT(IN) :: options
    REAL(KIND=8), INTENT(INOUT) :: fmem
    REAL(KIND=8), INTENT(IN) :: Enr
    INTEGER, INTENT(INOUT) :: LWORK
    LOGICAL, INTENT(INOUT) :: conv 
    INTEGER, INTENT(IN) :: norb,iter

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: A,B
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: WORK
    REAL(KIND=8) :: Eelc,temp,mdiff
    INTEGER :: INFO,stat1,stat2,stat3
    INTEGER :: i,j,u,v
    LOGICAL :: ex
  

    !write Density matrix  
    OPEN(unit=8,file='Da',access='sequential',status='old',form='unformatted')
    READ(8) Da(:,:) 
    READ(8) Db(:,:)
    CLOSE(unit=8)

    !get Guv matrix
    OPEN(unit=11,file='Guv',access='sequential',status='old',form='unformatted')
    READ(11) GuvA(:,:)
    READ(11) GuvB(:,:)
    CLOSE(unit=11) 

    !Create Fock matrix
    FuvA = Huv + GuvA
    FuvB = Huv + GuvB

    !Get Electronic energy
    Eelc = 0.0D0
    DO v=0,norb-1
      DO u=0,norb-1
        Eelc = Eelc + Da(u,v)*(FuvA(u,v)+Huv(u,v)) + Db(u,v)*(FuvB(u,v)+Huv(u,v))
      END DO
    END DO
    Eelc = 0.5D0*Eelc
 
999 FORMAT(4x,I3,4x,F20.15,8x,ES15.8)
996 FORMAT(4x,I3,4x,F20.15,8x,ES15.8,4x,ES15.8)

    !write energy output
    WRITE(*,996) iter+1, Eelc+Enr, Eelc+Enr-Etrk(2), mdiff

    !Get our new coefficients and eigenvalues

    ! setup arrays
    ALLOCATE(A(0:norb-1,0:norb-1),STAT=stat1)
    ALLOCATE(B(0:norb-1,0:norb-1),STAT=stat2)
    fmem = fmem - (norb*norb*8.0/1.0D6) 
    IF (fmem .LT. 0.0D6 .OR. stat1+stat2 .NE. 0) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      WRITE(*,*) "scf:UHFiter could not allocate enough memory"
      STOP "scf:UHFiter out of memory"
    END IF
    CALL nmem(fmem)

    !---  Spin case ALPHA  ---!
    DO v=0,norb-1
      A(0:v-1,v) = (/ (0.0D0, u=0,v-1) /)
      A(v:norb-1,v) = (/ (FuvA(u,v), u=v,norb-1) /)
      B(0:v-1,v) = (/ (0.0D0, u=0,v-1) /)
      B(v:norb-1,v) = (/ (Suv(u,v), u=v,norb-1) /)
    END DO

    !if first iteration...
    IF (iter .EQ. 0) THEN
      ALLOCATE(WORK(0:1))
      LWORK=-1
      CALL DSYGV(1,'V','L',norb,A,norb,B,norb,EigA,WORK,LWORK,INFO)
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
          WRITE(*,*) "scf:initUHF could not allocate memory"
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
          WRITE(*,*) "scf:initUHF could not allocate memory"
          STOP
        END IF
      END IF
    END IF

    !find Cui for Alpha 
    CALL DSYGV(1,'V','L',norb,A,norb,B,norb,EigA,WORK,LWORK,INFO)
    IF (INFO .GT. 0 .AND. INFO .LT. norb) THEN
      WRITE(*,*) "Warning: potential contamination of eigenvalues."
      WRITE(*,*) "Number of contaminents: ", INFO
    ELSE IF (INFO .GT. norb) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      WRITE(*,*) "Leading minor of order i of B is not positive definite. i=",INFO-norb
      STOP
    END IF

    !assign coefficients for Alpha
    CuiA(:,:) = A(:,:)

    !---  Spin case BETA  ---!
    DO v=0,norb-1
      A(0:v-1,v) = (/ (0.0D0, u=0,v-1) /)
      A(v:norb-1,v) = (/ (FuvB(u,v), u=v,norb-1) /)
      B(0:v-1,v) = (/ (0.0D0, u=0,v-1) /)
      B(v:norb-1,v) = (/ (Suv(u,v), u=v,norb-1) /)
    END DO

    !find Cui for Beta 
    CALL DSYGV(1,'V','L',norb,A,norb,B,norb,EigB,WORK,LWORK,INFO)
    IF (INFO .GT. 0 .AND. INFO .LT. norb) THEN
      WRITE(*,*) "Warning: potential contamination of eigenvalues."
      WRITE(*,*) "Number of contaminents: ", INFO
    ELSE IF (INFO .GT. norb) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      WRITE(*,*) "Leading minor of order i of B is not positive definite. i=",INFO-norb
      STOP
    END IF

    !assign coefficients for Beta 
    CuiB(:,:) = A(:,:)
 
    !---  Both spins  ---!

    !write MO 
    OPEN(unit=7,file='Cui',access='sequential',status='replace')
    WRITE(7,*) CuiA(:,:)
    WRITE(7,*) CuiB(:,:)
    CLOSE(unit=7,status='keep')

    !write new dens, old dens, and Guv
    OPEN(unit=9,file='Dold',status='replace',access='sequential',form='unformatted')
    WRITE(9) Da(:,:)
    WRITE(9) Db(:,:)
    CLOSE(unit=9)
    CALL EXECUTE_COMMAND_LINE('dens')
    CALL EXECUTE_COMMAND_LINE('UHFI2G')
 
    !cleanup memory
    DEALLOCATE(A)
    DEALLOCATE(WORK)
    fmem = fmem + (norb*norb*8.0/1.0D6)
    fmem = fmem + (LWORK*8/1.0D6)
    CALL nmem(fmem)

    !check for convergece
    Etrk(0) = Etrk(1)
    Etrk(1) = Etrk(2)
    Etrk(2) = Eelc + Enr

    IF (iter .GT. 0) THEN
      CALL checkConv(norb,Etrk,options,conv,mdiff)
    END IF

    IF (conv) THEN
      WRITE(*,*) "SCF has converged!"
      WRITE(*,*) 
      WRITE(*,*) "Total SCF energy (a.u)", Eelc+Enr
      WRITE(*,*)
    END  IF

  END SUBROUTINE UHFiter

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
  SUBROUTINE checkConv(norb,Etrk,options,conv,mdiff)
    IMPLICIT NONE

    !Value
    !Etrk	: 1D dp, list of 3 last energies
    !options	: 1D int, options array
    !conv	: bool, converged or not
    !Dna,Dnb	: 2D real8, new density matrix (spin,col,row) 
    !Doa,Dob	: 2D real8, old density matrix
    !norb	: int, number of orbitals
    !mdiff	: real8, max diff of density matrix

    !Inout
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Etrk
    INTEGER, DIMENSION(0:), INTENT(IN) :: options
    LOGICAL, INTENT(INOUT) :: conv
    INTEGER, INTENT(IN) :: norb

    !Internal
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Dna,Dnb,Doa,Dob
    REAL(KIND=8), DIMENSION(0:11) :: lim
    REAL(KIND=8) :: mdiff
    INTEGER :: stat1,stat2
    INTEGER :: i,j,k
    LOGICAL :: flag

    flag = .TRUE.
    conv = .FALSE.

    !These go to 11
    i = options(8)
    lim = [1.0D0,1.0D-1,1.0D-2,1.0D-3,1.0D-4,1.0D-5,1.0D-6,1.0D-7,1.0D-8,1.0D-9,1.0D-10,1.0D-11]

    !check energy convergence
    IF (Etrk(1)-Etrk(2) .LT. lim(i) .AND. Etrk(2) .LE. Etrk(1)) THEN
      flag = .FALSE.
    END IF 

    OPEN(unit=9,file='Dold',status='old',access='sequential',form='unformatted')
    OPEN(unit=10,file='Da',status='old',access='sequential',form='unformatted')
    !get rhf density 
    IF (options(3) .EQ. 0) THEN
      ALLOCATE(Doa(0:norb-1,0:norb-1),STAT=stat1)
      ALLOCATE(Dna(0:norb-1,0:norb-1),STAT=stat2)
      IF (stat1 + stat2 .NE. 0) THEN
        WRITE(*,*) "scf:checkConv could not allocate Dold or Dnew"
        STOP
      END IF 
      READ(9) Doa(:,:)
      READ(10) Dna(:,:)
    !get uhf density
    ELSE
      ALLOCATE(Doa(0:norb-1,0:norb-1),STAT=stat1)
      ALLOCATE(Dob(0:norb-1,0:norb-1),STAT=stat1)
      ALLOCATE(Dna(0:norb-1,0:norb-1),STAT=stat2)
      ALLOCATE(Dnb(0:norb-1,0:norb-1),STAT=stat2)
      IF (stat1 + stat2 .NE. 0) THEN
        WRITE(*,*) "scf:checkConv could not allocate Dold or Dnew"
        STOP
      END IF 
      READ(9) Doa(:,:)
      READ(9) Dob(:,:)
      READ(10) Dna(:,:)
      READ(10) Dnb(:,:)
    END IF

    CLOSE(unit=9)
    CLOSE(unit=10)

    !the old fashioned way
    mdiff = 0.0D0

    !rhf
    IF (options(3) .EQ. 0) THEN
      DO i=0,norb-1
        DO j=0,norb-1
          IF (ABS(Doa(i,j)-Dna(i,j)) .GT. mdiff) THEN
            mdiff = ABS(Doa(i,j)-Dna(i,j))
          END IF
        END DO
      END DO      

    !uhf - currently broken
    ELSE
      DO i=0,norb-1
        DO j=0,norb-1
          IF (ABS(Doa(i,j)-Dna(i,j)) .GT. mdiff) THEN
            mdiff = ABS(Doa(i,j)-Dna(i,j))
          END IF
          IF (ABS(Dob(i,j)-Dnb(i,j)) .GT. mdiff) THEN
            mdiff = ABS(Dob(i,j)-Dnb(i,j))
          END IF
        END DO
      END DO      
    END IF

    !check for convergence
    IF (mdiff .LT. lim(options(8))) THEN
      conv = .TRUE.
    ELSE
      conv = .FALSE.
    END IF

    DEALLOCATE(Doa)
    DEALLOCATE(Dna)
    IF (ALLOCATED(Dob)) THEN
      DEALLOCATE(Dob)
      DEALLOCATE(Dnb)
    END IF 

  END SUBROUTINE

!---------------------------------------------------------------------
!		Get expectation value of S^2 
!---------------------------------------------------------------------
  SUBROUTINE getSpinCont(norb,nelcA,nelcB,CuiA,CuiB,Suv)
    IMPLICIT NONE

    !Values
    !nelcA,B	: int, number of Alpha,Beta electrons
    !CuiA,B	: 2D dp, MO coeff in orthonormal basis
    !Suv	: 2D dp atomic overlap integrals
    !Arr1,Arr2	: 2D dp, dummy arrays

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: CuiA,CuiB,Suv
    INTEGER, INTENT(IN) :: norb,nelcA,nelcB

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Arr1,Arr2
    REAL(KIND=8) ::S,Ssqr,temp
    INTEGER :: i,j,u,v

999 FORMAT(1x,A14,2x,F8.5)
998 FORMAT(1x,A14,2x,F8.5)

    !no need to check memory here, we freed more than enough memory before this
    ALLOCATE(Arr1(0:norb-1,0:norb-1))
    ALLOCATE(Arr2(0:norb-1,0:norb-1))

    WRITE(*,*) "Spin Analysis"

    !exact
    S = 0.5D0*(nelcA - nelcB)
    WRITE(*,999) "<S^2> exact  :", S*(S+1)

    Ssqr = S*(S+1) + nelcB

     !this is how cfour does it, but I've not accounted for occupied orbitals
!    CALL DGEMM('N','N',norb,norb,norb,1.0D0,Suv,norb,CuiB,norb,0.0D0,Arr1,norb) !Suv is symmetric
!    CALL DGEMM('T','N',norb,norb,norb,1.0D0,CuiA,norb,Arr1,norb,0.0D0,Arr2,norb)
    !SHOULD call blas function here, but linker is acting strange...
!    DO i=0,norb-1
!      ssqr = ssqr - DDOT(norb,Arr2(i,:),1,Arr2(i,:),1)
!    END DO
!     DO i=0,nelcB-1
!     temp = 0
!       DO j=0,nelcB-1
!         temp = temp + Arr2(i,j)*Arr2(i,j)
!       END DO
!       Ssqr = Ssqr - temp 
!     END DO

    DO i=0,nelcA-1
      DO j=0,nelcB-1
        temp = 0.0D0
        DO u=0,norb-1
          DO v=0,norb-1
            temp = temp + CuiA(u,i)*Suv(u,v)*CuiB(v,j) 
          END DO
        END DO
        Ssqr = Ssqr - (temp)**2.0
      END DO
    END DO

    WRITE(*,998) "<S^2> actual :", Ssqr
    WRITE(*,*)

    DEALLOCATE(Arr1)
    DEALLOCATE(Arr2)

  END SUBROUTINE getSpinCont

!---------------------------------------------------------------------
!		Make MOLDEN input file	
!---------------------------------------------------------------------
  SUBROUTINE makeMOLDEN(atoms,xyz,EigA,EigB,ref,nnuc,norb,CuiA,CuiB)
    USE basis
    IMPLICIT NONE

    !Values
    !atoms	: 1D int, list of atoms in calculation
    !xyz	: 2D dp, list of nuclei positions
    !EigA,B	: 1D dp, list of MO eigenvalues	
    !ref	: int, reference 0-RHF,1-UHF,2-ROHF
    !nnuc	: int, number of atoms
    !norb	: int, number of orbitals

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz,CuiA,CuiB
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: EigA,EigB
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms
    INTEGER, INTENT(IN) :: ref,nnuc,norb

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: bas,set,val
    INTEGER, ALLOCATABLE, DIMENSION(:) :: basinfo,setinfo,dummy
    CHARACTER(LEN=2), DIMENSION(1:10) :: A
    CHARACTER(LEN=8), DIMENSION(0:0) :: B
    CHARACTER(LEN=1),DIMENSION(0:5) :: C 
    CHARACTER(LEN=8) :: line
    CHARACTER(LEN=2) :: Aname
    INTEGER :: i,j,k,maxN,maxL,bkey,sec,orb,nset
    INTEGER :: func,coef,pri,ang,ori,nelcA,nelcB

999 FORMAT(A2,4x,I3,4x,I3,4x,F15.8,4x,F15.8,4x,F15.8)

    A = ['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne']
    B = ['STO-3G']
    C = ['s','p','d','f','g','h']
    bkey = options(2)
    OPEN(unit=1,file='envdat',status='old',access='sequential')
    READ(1,*)
    READ(1,*) nelcA,nelcB
    CLOSE(unit=1)

    OPEN(unit=12,file='MOLDEN',status='replace',access='sequential')
    OPEN(unit=13,file='mybasis',status='old',access='sequential')

    !write geometry
    WRITE(12,*) "[Molden Format]"
    WRITE(12,*) "[ATOMS] AU"
    DO i=1,nnuc
      WRITE(12,999) A(atoms(i-1)), i, atoms(i-1), xyz(i-1,0:2) 
    END DO

    !write AO 
    WRITE(12,*) '[Molden Format]'
    WRITE(12,*) "[GTO]"

    DO i=0,nnuc-1
      WRITE(12,*) i+1,0
      !get information
      DO WHILE (line .NE. B(bkey))
        READ(13,*) line 
      END DO
      Aname = A(atoms(i))
      DO WHILE(line .NE. Aname)
        READ(13,*) line
      END DO

      READ(13,*) sec,orb,nset
      DO j=0,sec-1
        READ(13,*) func, coef, pri, ang, ori
        ALLOCATE(val(0:coef-1))
        WRITE(12,*) C(ang),func, "1.00"
        DO k=0,func-1 
          ! do stuff 
          READ(13,*) val
          WRITE(12,*) val(1:coef-1), val(0)
        END DO
        DEALLOCATE(val)
      END DO
      WRITE(12,*) 
      REWIND(13)
      READ(13,*) line
    END DO
    WRITE(12,*)
    CLOSE(unit=13) 


    !write MO
998 FORMAT(1x,A4,1x,F15.8)
997 FORMAT(1x,I4,4x,F15.10)
    WRITE(12,*) "[MO]"

    !write for alpha (guarenteed)
    DO i=0,norb-1
      WRITE(12,*) "Sym= A"
      WRITE(12,998) "Ene=",EigA(i)
      WRITE(12,*) " Spin= Alpha"
      IF (EigA(i) .LT. 0.0D0) THEN
        WRITE(12,*) "Occup= 1.0"
      ELSE
        WRITE(12,*) "Occup= 0.0" 
      END IF
      DO j=0,norb-1
        WRITE(12,997) j+1, CuiA(j,i)
      END DO
    END DO

    !write for beta (if needed)
    IF (ref .EQ. 1) THEN
      DO i=0,norb-1
        WRITE(12,*) "Sym= A"
        WRITE(12,998) "Ene=",EigB(i)
        WRITE(12,*) " Spin= Beta"
        IF (EigB(i) .LT. 0.0D0) THEN
          WRITE(12,*) "Occup= 1.0"
        ELSE
          WRITE(12,*) "Occup= 0.0" 
        END IF
        DO j=0,norb-1
          WRITE(12,997) j+1,CuiB(j,i)
        END DO
      END DO
    END IF
    CLOSE(unit=12,status='keep')

  END SUBROUTINE makeMOLDEN
!---------------------------------------------------------------------
END PROGRAM scf
