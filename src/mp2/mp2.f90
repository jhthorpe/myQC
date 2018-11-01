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
    CALL mp2_rhf(noccA,nvrtA,ntot)
  ELSE
    WRITE(*,*) "Sorry, that spin case has not been coded for MP2 yet"
    CALL EXECUTE_COMMAND_LINE('touch error')
    STOP "Bad ref in mp2"
  END IF

  CALL CPU_TIME(timeF)
  WRITE(*,*) "MP2 completed in (s): ", timeF-timeS

  CONTAINS

!---------------------------------------------------------------------
!       mp2_rhf
!               James H. Thorpe
!       - mp2 program RHF molecules
!       - currently reads in pqrs_AB...
!       - but really only needs ijab_AA, ijba_AA, and ijab_AB
!---------------------------------------------------------------------
  !Variables
  ! nocc        :       int, number of occupied orbitals
  ! nvrt        :       int, number of virtual orbitals
  ! ntot        :       int, total number of orbitals

  SUBROUTINE mp2_rhf(nocc,nvrt,ntot)
    !inout
    INTEGER, INTENT(IN) :: nocc,nvrt,ntot
    !internal
    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: pqrs
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: eig
    REAL(KIND=8) :: sum1, sum2
    INTEGER :: i,j,a,b,p,q,r,s

    ALLOCATE(pqrs(0:ntot-1,0:ntot-1,0:ntot-1,0:ntot-1))
    ALLOCATE(eig(0:ntot-1))

    OPEN(unit=101,file='pqrs_AB',status="old",access='sequential',form='unformatted')
    DO p=0,ntot-1
      DO q=0,ntot-1
        READ(101) pqrs(p,q,:,:)
      END DO
    END DO
    CLOSE(unit=101)

    !read in eigenvalues
    OPEN(unit=104,file='eig',status='old',access='sequential')
    READ(104,*) eig(:)
    CLOSE(unit=104)

    !AA part
    sum1=0.0D0
    DO i=0,nocc-1
      DO j=i+1,nocc-1
        DO a=nocc,ntot-1
          DO b=a+1,ntot-1
            sum1 = sum1 + (pqrs(i,j,a,b)-pqrs(i,j,b,a))**2.0D0/(eig(i)+eig(j)-eig(a)-eig(b))
          END DO
        END DO
      END DO
    END DO
    WRITE(*,*) "E(AA) = ", sum1+0.00

    !AB part
    sum2=0.0D0
    DO i=0,nocc-1
      DO j=0,nocc-1
        DO a=nocc,ntot-1
          DO b=nocc,ntot-1
            sum2 = sum2 + pqrs(i,j,a,b)**2.0D0/(eig(i)+eig(j)-eig(a)-eig(b)) 
          END DO
        END DO
      END DO
    END DO
    WRITE(*,*) "E(AB) = ", sum2+0.00

    WRITE(*,*) "E(MBPT2) = ", 2*sum1+sum2

    DEALLOCATE(pqrs)
    DEALLOCATE(eig)

  END SUBROUTINE mp2_rhf
!---------------------------------------------------------------------


END PROGRAM mp2

