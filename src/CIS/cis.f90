!//////////////////////////////////////////////////////////////////
!//             Performs CIS Calculation
!//
!//             James H. Thorpe, in the Group of John Stanton
!//             The Univeristy of Florida
!//
!//////////////////////////////////////////////////////////////////

!---------------------------------------------------------------------
!	cis
!		James H. Thorpe
!		Nov 27, 2018
!	-control program for cis calculations
!---------------------------------------------------------------------
  ! Variables
  ! noccA, nocc B       : int, number of alpha,beta occupied orbitals 
  ! options             : 1D int, options array
  ! ntot                : int, total number of orbitals
  ! nvrtA, nvrtB        : int, number of alpha,beta virtual orbitals

PROGRAM cis
  USE env
  USE linal 
  IMPLICIT NONE

  !Internal
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms, options
  INTEGER, DIMENSION(0:1) :: line
  INTEGER, DIMENSION(0:0) :: mem_lvl
  REAL(KIND=8) :: fmem, timeS, timeF
  INTEGER :: nnuc,noccA,noccB,nvrtA,nvrtB,ntot,dummy
  LOGICAL :: flag

999 FORMAT(1x,A22,2x,F8.4)
  CALL CPU_TIME(timeS)
  WRITE(*,*)
  WRITE(*,*)"                      STARTING CIS"
  WRITE(*,*)"-----------------------------------------------------------"
  WRITE(*,*) "cis called"
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

  !control sequence
  IF (options(1) .EQ. 0) THEN !SCF
    IF (options(3) .EQ. 1) THEN !UHF
      CALL cis_scf_uhf(noccA,noccB,nvrtA,nvrtB,ntot,options,fmem)
    ELSE
      WRITE(*,*) "Sorry, only UHF references have been coded."
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP "Bad ref in cis" 
    END IF
  ELSE
    WRITE(*,*) "Sorry, only SCF reference CIS has been coded."
    CALL EXECUTE_COMMAND_LINE('touch error')
    STOP "Bad method in cis"
  END IF

  CALL CPU_TIME(timeF)
  WRITE(*,999) "CIS completed in (s): ", timeF-timeS

  CONTAINS

!---------------------------------------------------------------------
!	cis_scf_uhf
!		James H. Thorpe
!		Nov 27, 2018
!	-cis program for SCF/UHF refernce wavefunction
!---------------------------------------------------------------------
  !Variables
  ! noccA,B		: int, number of alpha,beta occupied orbitals
  ! nvrtA, nvrtB        : int, number of alpha,beta virtual orbitals
  ! ntot                : int, total number of orbitals
  ! options             : 1D int, options array
  SUBROUTINE cis_scf_uhf(noccA,noccB,nvrtA,nvrtB,ntot,options,fmem)
    IMPLICIT NONE
    !Inout
    INTEGER, DIMENSION(0:), INTENT(IN) :: options
    REAL(KIND=8), INTENT(INOUT) :: fmem
    INTEGER, INTENT(IN) :: noccA,noccB,nvrtA,nvrtB,ntot
    !Internal
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Am
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ints,scf_eigs,cis_eigs
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: work
    REAL(KIND=8) :: mem
    INTEGER :: idx,upA,upB,lwork,info,neig
    INTEGER :: a,b,i,j


998 FORMAT(1x,A30,2x,F8.5)

    !Pre-process the needed files
    CALL cis_scf_uhf_proc(noccA,noccB,nvrtA,nvrtB,ntot)

    !Construct A matrix and check memory 
    upA = noccA*nvrtA
    upB = noccB*nvrtB
    mem = ((upA+upB)**2.0D0 + upA)
    mem = mem*(8.0D0/1.0D6)
    WRITE(*,998) "Allocating memory for CIS (MB)", mem
    IF (fmem - mem .LT. 10) THEN
      WRITE(*,*) "Could not allocate enough space, exiting"
      CALL EXECUTE_COMMAND_LINE('touch error')
      RETURN 
    END IF
    ALLOCATE(Am(0:upA+upB-1,0:upA+upB-1)) !and we pray we have memory
    fmem = fmem - mem
    CALL nmem(fmem) 
  
    WRITE(*,*) 
    WRITE(*,*) "Constructing A matrix"
    Am = 0.0D0
    ALLOCATE(scf_eigs(0:ntot-1))

    !-----------------------
    !AA block
    ALLOCATE(ints(0:upA-1))
    OPEN(unit=101,file='ints_AA',status='old',access='sequential',form='unformatted')
    idx = 0
    DO j=0,noccA-1
      DO b=0,nvrtA-1
        READ(101) ints(0:upA-1)
        Am(0:upA-1,idx) = ints(0:upA-1) 
        idx = idx + 1
      END DO
    END DO
    CLOSE(unit=101)
    DEALLOCATE(ints)
    !diagonal terms
    OPEN(unit=104,file='eig',status='old',access='sequential')
    READ(104,*) scf_eigs(0:ntot-1) 
    CLOSE(unit=104)
    idx = 0
    DO i=0,noccA-1
      DO a=0,nvrtA-1
        Am(idx,idx) = Am(idx,idx) + scf_eigs(noccA+a) - scf_eigs(i)
        idx = idx +1
      END DO
    END DO

    !-----------------------
    !AB block
    ALLOCATE(ints(0:upA-1))
    OPEN(unit=102,file='ajib_AB',status='old',access='sequential',form='unformatted')
    idx = upA
    DO j=0,noccB-1
      DO b=0,nvrtB-1
        READ(102) ints(0:upA-1)
        Am(0:upA-1,idx) = ints(0:upA-1) 
        idx = idx + 1
      END DO
    END DO
    CLOSE(unit=102)
    DEALLOCATE(ints)

    !-----------------------
    !BB block
    ALLOCATE(ints(0:upB-1))
    OPEN(unit=103,file='ints_BB',status='old',access='sequential',form='unformatted')
    idx = upA
    DO j=0,noccB-1
      DO b=0,nvrtB-1
        READ(103) ints(0:upB-1)
        Am(upA:upA+upB-1,idx) = ints(0:upB-1) 
        idx = idx + 1
      END DO
    END DO
    CLOSE(unit=103)
    !diagonal terms
    OPEN(unit=104,file='eig',status='old',access='sequential')
    READ(104,*) !dummy read to get past A coefs
    READ(104,*) scf_eigs(0:ntot-1) 
    CLOSE(unit=104)
    idx = 0
    DO j=0,noccB-1
      DO b=0,nvrtB-1
        Am(idx+upA,idx+upA) = Am(idx+upA,idx+upA) + scf_eigs(noccB+b) - scf_eigs(j)
        idx = idx +1
      END DO
    END DO

    !CALL linal_printmat_2Dreal8(Am,upA+upB,upA+upB)

    DEALLOCATE(scf_eigs)
    DEALLOCATE(ints)

    !-----------------------
    !Get eigenvalues - full spectrum for now
    WRITE(*,*) "Diagonalizing A matrix"
    ALLOCATE(cis_eigs(0:upA+upB-1))
    ALLOCATE(work(0:1))
    lwork = -1
    CALL DSYEV('V','U',upA+upB,Am,upA+upB,cis_eigs,work,lwork,info)
    lwork = CEILING(work(0))
    DEALLOCATE(work)
    ALLOCATE(work(0:lwork-1))
    CALL DSYEV('V','U',upA+upB,Am,upA+upB,cis_eigs,work,lwork,info)

    !print eigenvalues and relevent vectors
    neig = options(15)
    IF (neig .GT. upA+upB) neig = upA+upB
    CALL cis_uhf_print(neig,cis_eigs(0:neig-1),Am(0:upA+upB-1,0:neig-1), &
                       noccA,nvrtA,noccB,nvrtB,ntot)
    
    DEALLOCATE(work)
    DEALLOCATE(cis_eigs)
    DEALLOCATE(Am)
    fmem = fmem + mem
    CALL nmem(fmem)

  END SUBROUTINE cis_scf_uhf

!---------------------------------------------------------------------
!	cis_scf_uhf_proc
!		James H. Thorpe
!		Nov 27, 2018
!	-preprocess integral files for CIS
!	-we are going to need the <aj||ib> integrals for AA and BB
!	 spin cases, which are <aj|ib> - <aj|bi>
!	-we do this one block, one vector at a time
!---------------------------------------------------------------------
  !Variables
  ! noccA,B		: int, number of alpha,beta occupied orbitals
  ! nvrtA, nvrtB        : int, number of alpha,beta virtual orbitals
  ! ntot                : int, total number of orbitals
  SUBROUTINE cis_scf_uhf_proc(noccA,noccB,nvrtA,nvrtB,ntot)
    IMPLICIT NONE
    !Inout
    INTEGER, INTENT(IN) :: noccA,noccB,nvrtA,nvrtB,ntot
    !Internal
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ajib,ajbi,proc
    REAL(KIND=8) :: t1,t2
    INTEGER :: up
    INTEGER :: a,b,i,j

998 FORMAT(1x,A25,2x,F8.4)
   
    CALL CPU_TIME(t1)

    WRITE(*,*) 
    WRITE(*,*) "Pre-proccessing 2e- inegrals in MO basis"

    !AA block
    up = noccA*nvrtA
    ALLOCATE(ajib(0:up-1))
    ALLOCATE(ajbi(0:up-1))
    ALLOCATE(proc(0:up-1))
    OPEN(unit=101,file='ajib_AA',status='old',access='sequential',form='unformatted')
    OPEN(unit=102,file='ajbi_AA',status='old',access='sequential',form='unformatted')
    OPEN(unit=103,file='ints_AA',status='replace',access='sequential',form='unformatted')
    DO j=0,noccA-1
      DO b=0,nvrtA-1
        READ(101) ajib(0:up-1)
        READ(102) ajbi(0:up-1) 
        proc(0:up-1) = ajib(0:up-1) - ajbi(0:up-1) 
        WRITE(103) proc(0:up-1) 
      END DO
    END DO
    CLOSE(unit=103)
    CLOSE(unit=102)
    CLOSE(unit=101)
    DEALLOCATE(ajib)
    DEALLOCATE(ajbi)
    DEALLOCATE(proc)

    !BB block
    up = noccB*nvrtB
    ALLOCATE(ajib(0:up-1))
    ALLOCATE(ajbi(0:up-1))
    ALLOCATE(proc(0:up-1))
    OPEN(unit=101,file='ajib_BB',status='old',access='sequential',form='unformatted')
    OPEN(unit=102,file='ajbi_BB',status='old',access='sequential',form='unformatted')
    OPEN(unit=103,file='ints_BB',status='replace',access='sequential',form='unformatted')
    DO j=0,noccB-1
      DO b=0,nvrtB-1
        READ(101) ajib(0:up-1)
        READ(102) ajbi(0:up-1) 
        proc(0:up-1) = ajib(0:up-1) - ajbi(0:up-1) 
        WRITE(103) proc(0:up-1) 
      END DO
    END DO
    CLOSE(unit=103)
    CLOSE(unit=102)
    CLOSE(unit=101)
    DEALLOCATE(ajib)
    DEALLOCATE(ajbi)
    DEALLOCATE(proc)

    CALL CPU_TIME(t2)
    WRITE(*,998) "<aj||ib> processed in (s)", (t2-t1) 
  
  END SUBROUTINE cis_scf_uhf_proc

!---------------------------------------------------------------------
!	cis_uhf_print
!		James H. Thorpe
!		Nov 27, 2018
!	- prints requested nubmer of CIS eigenvalues and eigenvectors
!---------------------------------------------------------------------
  !Variables
  ! neig		: int, number of eigenpairs to print
  ! eval		: 1D real8, list of eigenvalues
  ! evec		: 2D real8, list of eigenvectors (row,col)
  ! noccA,B		: int, number of occupied A,B orbitals		
  ! nvrtA,B		: int, number of virtul A,B orbitals
  ! ntot		: int, number of total orbitals
  SUBROUTINE cis_uhf_print(neig,eval,evec,noccA,nvrtA,noccB,nvrtB,ntot)
    IMPLICIT NONE
    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: evec
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: eval
    INTEGER,INTENT(IN) :: neig,noccA,noccB,nvrtA,nvrtB,ntot 
    !Internal
    INTEGER :: i, up

    WRITE(*,*) "========================================================================="
    WRITE(*,*) "                          CIS Results"
    WRITE(*,*) "========================================================================="

    up = noccA*nvrtA+noccB*nvrtB
    
    DO i=0,neig-1
      WRITE(*,*) 
      WRITE(*,*) "CIS eigenpair #",i+1
      WRITE(*,*) "Eigenvalue (a.u.)", eval(i) 
      CALL cis_uhf_vec(evec(0:up-1,i),noccA,nvrtA,noccB,nvrtB)
      WRITE(*,*) "------------------------------------------------"
    END DO

  END SUBROUTINE cis_uhf_print

!---------------------------------------------------------------------
!	cis_uhf_vec
!		James H. Thorpe
!		Nov 27, 2018
!	-finds relevent contributions to CIS eigenvector
!---------------------------------------------------------------------
  !Variables
  ! evec		: 2D real8, list of eigenvectors (row,col)
  ! noccA,B		: int, number of occupied A,B orbitals		
  ! nvrtA,B		: int, number of virtul A,B orbitals
  SUBROUTINE cis_uhf_vec(evec,noccA,nvrtA,noccB,nvrtB)
    IMPLICIT NONE
    !Inout
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: evec
    INTEGER,INTENT(IN) :: noccA,noccB,nvrtA,nvrtB
    !Internal
    INTEGER, DIMENSION(:), ALLOCATABLE :: sgn1
    CHARACTER(LEN=2), DIMENSION(0:9) :: spin 
    REAL(KIND=8), DIMENSION(0:9) :: Cia
    INTEGER, DIMENSION(0:9,0:3) :: iajb
    INTEGER, DIMENSION(0:9) :: sgn2
    INTEGER :: upA,upB,idx,loc
    INTEGER :: i,j,a,b
999 FORMAT(4x,I3,1x,I3,1x,I3,1x,I3,4x,F15.10,4x,A2)
    spin  = 'NA'
    Cia = 0
    sgn2 = 1
    iajb = 0
    upA = noccA*nvrtA
    upB = noccB*nvrtB

    !Split evec into absolute value and sign for sorting
    ALLOCATE(sgn1(0:upA+upB-1))
    sgn1 = 1
    DO i=0,upA+upB-1
      IF (evec(i) .LT. 0.0D0) sgn1(i) = -1 
    END DO

    !search though and find most significant contributions
    idx = 0
    DO i=0,noccA-1
      DO a=0,nvrtA-1
        IF (ABS(evec(idx)) .GE. MINVAL(Cia)) THEN
          loc = MINLOC(Cia,1)-1
          Cia(loc) = ABS(evec(idx))
          sgn2(loc) = sgn1(idx)
          iajb(loc,0:3) = [i+1,noccA+a+1,0,0]
          spin(loc) = 'AA'    
        END IF  
        idx = idx + 1
      END DO
    END DO

    idx = upA
    DO j=0,noccB-1
      DO b=0,nvrtB-1
        IF (ABS(evec(idx)) .GE. MINVAL(Cia)) THEN
          loc = MINLOC(Cia,1)-1
          Cia(loc) = ABS(evec(idx))
          sgn2(loc) = sgn1(idx)
          iajb(loc,0:3) = [0,0,j+1,noccB+b+1]
          spin(loc) = 'BB'    
        END IF  
        idx = idx + 1
      END DO
    END DO

    !Dirty sort
    

    !write out the results
    WRITE(*,*) "     i   a   j   b      val             spin" 
    WRITE(*,*) "------------------------------------------------"
    
    DO i=0, MIN(upA+upB-1,9)
      IF (Cia(i) .GT. 1.0D-10) WRITE(*,999) iajb(i,0:3),sgn2(i)*Cia(i),spin(i) 
    END DO 

    DEALLOCATE(sgn1)

  END SUBROUTINE cis_uhf_vec
!---------------------------------------------------------------------
END PROGRAM cis

