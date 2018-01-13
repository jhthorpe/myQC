!//////////////////////////////////////////////////////////////////
!//          Constructs Guv(2e-) Matrix for RHF system 
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

PROGRAM RHFI2G
  USE env
  IMPLICIT NONE

  ! Values
  ! xyz         : 2D dp, array of nuclear positions
  ! atoms       : 1D int, array of which atom is which
  ! fmem        : dp, free memory left in MB
  ! nnuc        : int, number of nuclii
  ! nelc        : int, number of electrons
  ! options     : 1D int, array of options
  ! Guv		: 2D dp, 2e- part of matrix
  ! Da		: 2D dp, density matrix 
  ! XX		: 4D dp, 2e- AO integral matrix

  ! Internal
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: XX 
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz,Guv,Da
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms, options
  INTEGER, DIMENSION(0:1) :: line
  REAL(KIND=8) :: fmem,timeS,timeF,temp
  INTEGER :: nnuc,nelc,dummy,norb,stat1,stat2,stat3
  INTEGER :: i,j,k,l
  LOGICAL :: flag

  CALL getenv(nnuc,nelc,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  !get number of orbitals
  OPEN(unit=1,file='basinfo',status='old',access='sequential')
  READ(1,*) line
  norb = line(1)
  CLOSE(unit=1)

  !Allocate memory
  fmem = fmem - (norb**4.0D0 + 2*norb*norb)*8.0/1.0D6
  IF (fmem .GT. 0) THEN
    ALLOCATE(Guv(0:norb-1,0:norb-1),STAT=stat1)
    ALLOCATE(Da(0:norb-1,0:norb-1),STAT=stat2)
    ALLOCATE(XX(0:norb-1,0:norb-1,0:norb-1,0:norb-1),STAT=stat3)
    IF (stat1+stat2+stat3 .NE. 0) THEN
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
  CLOSE(unit=10)

  !my INCREDIBLY dumb implimentation
  DO j=0,norb-1
    DO i=0,norb-1
      temp = 0
      DO l=0,norb-1
        DO k=0,norb-1
          temp = temp + Da(k,l)*(XX(i,j,k,l) - 0.25D0*XX(i,k,j,l) - 0.25D0*XX(i,l,j,k)) 
        END DO
      END DO
      Guv(i,j) = temp
    END DO
  END DO
  
  !once we have it, write out Guv
  OPEN(unit=11,file='Guv',status='replace',access='sequential',form='unformatted')
  WRITE(11) Guv(:,:)
  CLOSE(unit=11)

  !Free memory
  DEALLOCATE(Guv)
  DEALLOCATE(Da)
  DEALLOCATE(XX)
  fmem = fmem - (norb**4.0D0 + 2.0D0*norb*norb)*8.0/1.0D6
  CALL nmem(fmem)

END PROGRAM RHFI2G
