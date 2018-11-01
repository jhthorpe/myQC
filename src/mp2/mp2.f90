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
  WRITE(*,*)"                      STARTING MP2"
  WRITE(*,*)"------------------------------------------------------------"
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
  ELSE IF (options(3) .EQ. 1) THEN
    CALL mp2_uhf(noccA,noccB,nvrtA,nvrtB,ntot)
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
!---------------------------------------------------------------------
  !Variables
  ! nocc        :       int, number of occupied orbitals
  ! nvrt        :       int, number of virtual orbitals
  ! ntot        :       int, total number of orbitals

  SUBROUTINE mp2_rhf(nocc,nvrt,ntot)
    !inout
    INTEGER, INTENT(IN) :: nocc,nvrt,ntot
    !internal
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ijab
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: eig
    REAL(KIND=8) :: sum1, sum2
    INTEGER :: i,j,a,b,p,q,r,s

    ALLOCATE(ijab(0:nvrt-1,0:nvrt-1))
    ALLOCATE(eig(0:ntot-1))

    !read in eigenvalues
    OPEN(unit=104,file='eig',status='old',access='sequential')
    READ(104,*) eig(:)
    CLOSE(unit=104)
 
    !summations
    sum1 = 0.0D0
    sum2 = 0.0D0
    OPEN(unit=102,file='ijab_AB',status='old',access='sequential',form='unformatted')
    DO i=0,nocc-2 !AA/AB i loop
      DO j=0,i    !AB j loop
        READ(102) ijab(:,:)
        DO a=0,nvrt-1
          DO b=0,nvrt-1
            sum2 = sum2 + ijab(a,b)**2.0D0/(eig(i)+eig(j)-eig(a+nocc)-eig(b+nocc))
          END DO
        END DO
      END DO
      DO j=i+1,nocc-1 !AA/AB j loop
        READ(102) ijab(:,:)
        DO a=0,nvrt-2 !AA/AB a loop
          DO b=0,a !AB a loop
            sum2 = sum2 + ijab(a,b)**2.0D0/(eig(i)+eig(j)-eig(a+nocc)-eig(b+nocc))
          END DO
          DO b=a+1,nvrt-1 !AA/AB b loop
            sum2 = sum2 +ijab(a,b)**2.0D0/(eig(i)+eig(j)-eig(a+nocc)-eig(b+nocc))
            sum1 = sum1 +(ijab(a,b)-ijab(b,a))**2.0D0/(eig(i)+eig(j)-eig(a+nocc)-eig(b+nocc))
          END DO
        END DO
        !we missed a=nvrt-1
        DO b=0,nvrt-1
          sum2 = sum2 +ijab(nvrt-1,b)**2.0D0/(eig(i)+eig(j)-eig(ntot-1)-eig(b+nocc))
        END DO
      END DO
    END DO
    !we missed i=nocc-1
    DO j=0,nocc-1
      READ(102) ijab(:,:)
      DO a=0,nvrt-1
        DO b=0,nvrt-1
          sum2 = sum2+ijab(a,b)**2.0D0/(eig(nocc-1)+eig(j)-eig(a+nocc)-eig(b+nocc))
        END DO
      END DO
    END DO
    CLOSE(unit=102)

    WRITE(*,*) "E(AA) = ", sum1*1.0D0
    WRITE(*,*) "E(AB) = ", sum2*1.0D0
    WRITE(*,*) "E(MBPT2) = ", 2*sum1+sum2

    DEALLOCATE(ijab)
    DEALLOCATE(eig)

  END SUBROUTINE mp2_rhf
!---------------------------------------------------------------------
!       mp2_uhf
!               James H. Thorpe
!               Nov. 1, 2018
!       - mp2 program for UHF systems
!---------------------------------------------------------------------
  ! noccA,B        :       int, number of occupied A,B orbitals
  ! nvrtA,B        :       int, number of virtual A,B orbitals
  ! ntot        :       int, total number of orbitals
  SUBROUTINE mp2_uhf(noccA,noccB,nvrtA,nvrtB,ntot)
    IMPLICIT NONE
    !inout
    INTEGER, INTENT(IN) :: noccA,noccB,nvrtA,nvrtB,ntot
    !internal
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ijab
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: eigA,eigB
    REAL(KIND=8) :: sum1,sum2,sum3
    INTEGER :: i,j,a,b,p,q,r,s

    ALLOCATE(eigA(0:ntot-1))
    ALLOCATE(eigB(0:ntot-1))

    !read in eigenvalues
    OPEN(unit=104,file='eig',status='old',access='sequential')
    READ(104,*) eigA(:)
    READ(104,*) eigB(:)
    CLOSE(unit=104)

    !spin case AA
    ALLOCATE(ijab(0:nvrtA-1,0:nvrtA-1))
    sum1 = 0.0D0
    OPEN(unit=102,file='ijab_AA',status='old',access='sequential',form='unformatted')
    DO i=0,noccA-2
      DO j=i+1,noccA-1
        READ(102) ijab(:,:)
        DO a=0,nvrtA-2
          DO b=a+1,nvrtA-1
            sum1 = sum1+(ijab(a,b)-ijab(b,a))**2.0D0/(eigA(i)+eigA(j)-eigA(a+noccA)-eigA(b+noccA))
          END DO
        END DO
      END DO
    END DO
    CLOSE(unit=102)
    DEALLOCATE(ijab)
    
    !spin case BB
    ALLOCATE(ijab(0:nvrtB-1,0:nvrtB-1))
    sum2 = 0.0D0
    OPEN(unit=102,file='ijab_BB',status='old',access='sequential',form='unformatted')
    DO i=0,noccB-2
      DO j=i+1,noccB-1
        READ(102) ijab(:,:)
        DO a=0,nvrtB-2
          DO b=a+1,nvrtB-1
            sum2 = sum2+(ijab(a,b)-ijab(b,a))**2.0D0/(eigB(i)+eigB(j)-eigB(a+noccB)-eigB(b+noccB))
          END DO
        END DO
      END DO
    END DO
    CLOSE(unit=102)
    DEALLOCATE(ijab)

    !spin case AB
    ALLOCATE(ijab(0:nvrtA-1,0:nvrtB-1))
    sum3 = 0.0D0
    OPEN(unit=102,file='ijab_AB',status='old',access='sequential',form='unformatted')
    DO i=0,noccA-1
      DO j=0,noccB-1
        READ(102) ijab(:,:)   
        DO a=0,nvrtA-1
          DO b=0,nvrtB-1
            sum3 = sum3+ijab(a,b)**2.0D0/(eigA(i)+eigB(j)-eigA(a+noccA)-eigB(b+noccB))
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(ijab)
    CLOSE(unit=102)

    WRITE(*,*) "E(AA) = ", sum1
    WRITE(*,*) "E(BB) = ", sum2
    WRITE(*,*) "E(AB) = ", sum3 
    WRITE(*,*) "E(MBPT2) = ", sum1+sum2+sum3

    DEALLOCATE(eigA)
    DEALLOCATE(eigB)

  END SUBROUTINE mp2_uhf

!---------------------------------------------------------------------

END PROGRAM mp2

