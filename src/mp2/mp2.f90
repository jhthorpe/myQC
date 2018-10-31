!//////////////////////////////////////////////////////////////////
!//		Performs MP2 Calculation
!//
!//		James H. Thorpe, in the Group of John Stanton
!//		The Univeristy of Florida
!//
!//////////////////////////////////////////////////////////////////

!---------------------------------------------------------------------
!	mp2
!		James H. Thorpe
!		Oct 22, 2018
!	- control program for mp2 calculations
!---------------------------------------------------------------------
  ! Variables
  ! noccA, nocc B	:	int, number of alpha,beta occupied orbitals 
  ! options		:	1D int, options array
  ! ntot		:	int, total number of orbitals
  ! nvrtA, nvrtB	:	int, number of alpha,beta virtual orbitals

PROGRAM mp2
  USE env
  IMPLICIT NONE

  !Internal
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms, options
  INTEGER, DIMENSION(0:1) :: line
  INTEGER, DIMENSION(0:0) :: mem_lvl
  REAL(KIND=8) :: fmem, timeS, timeF
  INTEGER :: nnuc,noccA,noccB,nvrtA,nvrtB,ntot,dummy
  LOGICAL :: flag

  CALL CPU_TIME(timeS)
  WRITE(*,*)
  WRITE(*,*) "mp2 called"
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

  IF (options(3) .EQ. 0) THEN !RHF 
    IF (options(12) .EQ. 0) THEN
      CALL mp2_rhf(noccA,nvrtA,ntot,options)  
    ELSE IF (options(12) .EQ. 1) THEN 
      CALL slow_mp2_rhf(noccA,nvrtA,ntot)
    END IF
  ELSE
    WRITE(*,*) "Sorry, that spin case has not been coded for MP2 yet"
    CALL EXECUTE_COMMAND_LINE('touch error')
    STOP "Bad ref in mp2"
  END IF

  CALL CPU_TIME(timeF)
  WRITE(*,*) "MP2 completed in (s): ", timeF-timeS

  CONTAINS

!---------------------------------------------------------------------
!	mp2_rhf
!		Oct 22, 2018
!		James H. Thorpe
!	- performs mp2 calculation for rhf reference
!---------------------------------------------------------------------
  ! Variables
  ! nocc	:	int, number of occupied orbitals
  ! nvrt	:	int, number of virtual orbitals
  ! ntot	:	int, number ot all orbitals
  ! options	:	int, options array
  ! eig		:	1Dint, SCF eigenvalues
  ! ijab_test,ijba_test	:	2D real8, matrices of <ij|ab> and <ij|ba> integrals
  SUBROUTINE mp2_rhf(nocc,nvrt,ntot,options)
    IMPLICIT NONE

    !inout
    INTEGER, DIMENSION(0:), INTENT(IN) :: options
    INTEGER, INTENT(IN) :: nocc,nvrt,ntot

    !internal
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ijab_test,ijba_test
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: eig
    REAL(KIND=8) :: sum1, sum2, ijab, ijba, sum3
    INTEGER :: i,j,a,b
    LOGICAL :: flag1,flag2
 
    !check that our files are there
    INQUIRE(file='ijab',exist=flag1)
    INQUIRE(file='ijba',exist=flag2)
    IF (.NOT. flag1 .OR. .NOT. flag2) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP "missing ijab or ijba integral files"
    END IF 

    ALLOCATE(ijab_test(0:nvrt-1,0:nvrt-1))
    ALLOCATE(ijba_test(0:nvrt-1,0:nvrt-1))
    ALLOCATE(eig(0:ntot-1))

    sum1 = 0.0D0
    sum2 = 0.0D0
    sum3 = 0.0D0

    !read in eigenvalues
    OPEN(unit=104,file='eig',status='old',access='sequential')
    READ(104,*) eig(:)
    CLOSE(unit=104)

    OPEN(unit=102,file='ijab',status='old',access='sequential',form='unformatted')
    OPEN(unit=103,file='ijba',status='old',access='sequential',form='unformatted')

    WRITE(*,*) "TESTING TESTING TESTING"
    DO i=0,nocc-1
      DO j=0,nocc-1
!    DO i=0,nocc-2
!      DO j=i+1,nocc-1
        ijab_test=0
        ijba_test=0
        DO a=0,nvrt-1
          DO b=0,nvrt-1
!        DO a=0,nvrt-2
!          DO b=a+1,nvrt-1
            !READ(102) ijab_test(a,b)  
            !READ(103) ijba_test(a,b)
            READ(102) ijab
            READ(103) ijba
            ijab_test(a,b) = ijab
            ijba_test(a,b) = ijba
  !          WRITE(*,*) "i,j,a,b | ijab | ijba", i,j,a,b, "|", ijab, "|",  ijba
            sum1 = sum1 + ijab**2.0/(eig(i)+eig(j)-eig(nocc+a)-eig(nocc+b))
            sum2 = sum2 - ijab*ijba/(eig(i)+eig(j)-eig(nocc+a)-eig(nocc+b))
            sum3 = sum3 + ijab*(2*ijab-ijba)/(eig(i)+eig(j)-eig(nocc+a)-eig(nocc+b))
  !          WRITE(*,*) "sum1 | sum2 | denom", sum1 , "|", sum2, "|", eig(i)+eig(j)-eig(nocc+a)-eig(nocc+b)
          END DO
        END DO

        !testing stuff
        IF (i.EQ.0.AND.j.EQ.1)  THEN
          WRITE(*,*) "i=0,j=1,a=0,b=1", ijab_test(0,1)
          WRITE(*,*) "i=0,j=1,b=1,a=0", ijba_test(0,1)
        END IF
        !end testing stuff

      END DO  ! j loop
    END DO ! i loop
    WRITE(*,*) "TESTING TESTING TESTING"

    WRITE(*,*) "TESTING MP2"
      WRITE(*,*) "sum1= ", sum1
      WRITE(*,*) "sum2= ", sum2
      WRITE(*,*) "sum3= ", sum3
    WRITE(*,*) "END TESTING MP2"

     
    CLOSE(unit=102)
    CLOSE(unit=103)

    DEALLOCATE(ijab_test)
    DEALLOCATE(ijba_test)
    
  END SUBROUTINE mp2_rhf

!---------------------------------------------------------------------
!       slow_mp2a_rhf
!               James H. Thorpe
!       - mp2 program for dumb mp2 calcs...
!       - where ao2mo used dumb_ao2mo_rhf
!---------------------------------------------------------------------
  SUBROUTINE slow_mp2_rhf(nocc,nvrt,ntot)
    !inout
    INTEGER, INTENT(IN) :: nocc,nvrt,ntot
    !internal
    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: pqrs
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: eig
    REAL(KIND=8) :: sum1, sum2
    INTEGER :: i,j,a,b,p,q,r,s

    ALLOCATE(pqrs(0:ntot-1,0:ntot-1,0:ntot-1,0:ntot-1))
    ALLOCATE(eig(0:ntot-1))


    OPEN(unit=101,file='moints',status="old",access='sequential')
    DO p=0,ntot-1
      DO q=0,ntot-1
        DO r=0,ntot-1
          DO s=0,ntot-1
            READ(101,*) pqrs(p,q,r,s) 
          END DO
        END DO 
      END DO
    END DO
    CLOSE(unit=101)

    !read in eigenvalues
    OPEN(unit=104,file='eig',status='old',access='sequential')
    READ(104,*) eig(:)
    CLOSE(unit=104)

    sum1 = 0.0D0

    DO i=0,nocc-1
      DO j=0,nocc-1
        DO a=nocc,ntot-1
          DO b=nocc,ntot-1
            sum1 = sum1 + pqrs(i,j,a,b)*(2*pqrs(i,j,a,b)-pqrs(i,j,b,a))/(eig(i)+eig(j)-eig(a)-eig(b))
          END DO
        END DO
      END DO
    END DO

    WRITE(*,*) "MP2 energy : ",sum1

    DEALLOCATE(pqrs)
    DEALLOCATE(eig)

  END SUBROUTINE slow_mp2_rhf
!---------------------------------------------------------------------


END PROGRAM mp2

