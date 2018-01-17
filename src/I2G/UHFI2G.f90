!//////////////////////////////////////////////////////////////////
!//          Constructs Guv(2e-) Matrix for UHF system 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//             
!//  WORK NOTE - currently this uses a very inefficient algrothim,
!//    at some point the in the future, it should be implimented
!//    using: Raffenetti, Chemical Physics Letters, 1973             
!///////////////////////////////////////////////////////////////////

!=====================================================================
!                       MAIN 

PROGRAM UHFI2G
  USE env
  IMPLICIT NONE

  ! Values
  ! xyz         : 2D dp, array of nuclear positions
  ! atoms       : 1D int, array of which atom is which
  ! fmem        : dp, free memory left in MB
  ! nnuc        : int, number of nuclii
  ! nelcA,nelcB : int, number of electrons of spin A,B
  ! options     : 1D int, array of options
  ! GuvA,B	: 2D dp, 2e- part of matrix of spin A,B
  ! Da		: 2D dp, density matrix of spin A,B
  ! XX		: 4D dp, 2e- AO integral matrix

  ! Internal
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: XX 
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz,GuvA,GuvB,Da,Db
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms, options
  INTEGER, DIMENSION(0:1) :: line
  REAL(KIND=8) :: fmem,timeS,timeF,tempA,tempB
  INTEGER :: nnuc,nelcA,nelcB,nelc,dummy,norb,stat1,stat2,stat3,stat4
  INTEGER :: i,j,k,l
  LOGICAL :: flag

  CALL getenv(nnuc,nelcA,nelcB,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  !get number of orbitals
  OPEN(unit=1,file='basinfo',status='old',access='sequential')
  READ(1,*) line
  norb = line(1)
  CLOSE(unit=1)

  !Allocate memory
  fmem = fmem - (norb**4.0D0 + 4*norb*norb)*8.0/1.0D6
  IF (fmem .GT. 0) THEN
    ALLOCATE(GuvA(0:norb-1,0:norb-1),STAT=stat1)
    ALLOCATE(GuvB(0:norb-1,0:norb-1),STAT=stat2)
    ALLOCATE(Da(0:norb-1,0:norb-1),STAT=stat3)
    ALLOCATE(Db(0:norb-1,0:norb-1),STAT=stat4)
    ALLOCATE(XX(0:norb-1,0:norb-1,0:norb-1,0:norb-1),STAT=stat3)
    IF (stat1+stat2+stat3+stat4 .NE. 0) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      WRITE(*,*) "I2G:RHFI2G failed to allocate memory"
      STOP
    END IF
  ELSE
    WRITE(*,*) "I2G:RHFI2G max memory reached"
    CALL EXECUTE_COMMAND_LINE('touch error')
    STOP
  END IF
  CALL nmem(fmem)

  !read in intermediate and density
  OPEN(unit=9,file='XX',status='old',access='sequential',form='unformatted')
  READ(9) XX(:,:,:,:)
  CLOSE(unit=9)
  OPEN(unit=10,file='Da',status='old',access='sequential',form='unformatted')
  READ(10) Da(:,:)
  READ(10) Db(:,:)
  CLOSE(unit=10)

  !my INCREDIBLY dumb implimentation
  DO j=0,norb-1
    DO i=0,norb-1
      tempA = 0.0D0
      tempB = 0.0D0
      DO l=0,norb-1
        DO k=0,norb-1
          tempA = tempA + (Da(k,l)+Db(k,l))*XX(i,j,k,l) - Da(k,l)*XX(i,k,j,l)
          tempB = tempB + (Da(k,l)+Db(k,l))*XX(i,j,k,l) - Db(k,l)*XX(i,k,j,l)
        END DO
      END DO
      GuvA(i,j) = tempA
      GuvB(i,j) = tempB
    END DO
  END DO
  
  !once we have it, write out Guv
  OPEN(unit=11,file='Guv',status='replace',access='sequential',form='unformatted')
  WRITE(11) GuvA(:,:)
  WRITE(11) GuvB(:,:)
  CLOSE(unit=11)

  !Free memory
  DEALLOCATE(GuvA)
  DEALLOCATE(GuvB)
  DEALLOCATE(Da)
  DEALLOCATE(Db)
  DEALLOCATE(XX)

  fmem = fmem - (norb**4.0D0 + 4.0D0*norb*norb)*8.0/1.0D6
  CALL nmem(fmem)

END PROGRAM UHFI2G
