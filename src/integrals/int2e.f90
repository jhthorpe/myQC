!//////////////////////////////////////////////////////////////////
!//            Performs 2 e- integrals for myQC 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//             
!//             Integration scheme a la McMurchie,Davidson 1978 
!//             
!//     WORK NOTE - location of coefficients in bas currently hardcoded 
!///////////////////////////////////////////////////////////////////

!=====================================================================
!                       MAIN 
PROGRAM int2e
  USE env
  USE basis
  USE aux

  IMPLICIT NONE

  ! Values
  ! xyz         : 2D dp, array of nuclear positions
  ! atoms       : 1D int, array of which atom is which
  ! fmem        : dp, free memory left in MB
  ! nnuc        : int, number of nuclii
  ! nelc        : int, number of electrons
  ! norb        : int, number of orbitals in molecule
  ! npri        : int, number of primatives
  ! bas         : 2D dp, basis for each atom: atom : orbital : [d,a]
  ! basinfo     : 2D int, array of basis information
  ! options     : 1D int, array of options
  ! set		: 2D dp, array of exponential coefficients of sets: atom,orbital numbs
  ! setinfo	: 2D int, array of set information
  ! maxN	: int, max principle quantum number of basis
  ! maxL	: int, max angular quantum number of basis

  ! Variables  
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: bas,set
  INTEGER, ALLOCATABLE, DIMENSION(:) :: basinfo,setinfo
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms,options
  REAL(KIND=8) :: timeS, timeF, fmem
  INTEGER :: nnuc,nelc,i,j,k,norb,npri,stat,maxN,maxL
  LOGICAL :: flag1,flag2,flag

  ! input managment 
  CALL CPU_TIME(timeS)
  WRITE(*,*)
  WRITE(*,*) "int2e called"
  CALL getenv(nnuc,nelc,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  ! build the basis set
  CALL buildBasis(options(2),atoms,bas,basinfo,set,setinfo,.FALSE.,maxN,maxL)

  ! check that Fock has been created 
  INQUIRE(file='Fuv',EXIST=flag1)
 
  ! check that we need to do these calculations
  INQUIRE(file='Fold',EXIST=flag2)
 
  ! logic gate
  IF (flag1 .AND. flag2) THEN
    WRITE(*,*) "Reading two electron integrals from intermediate"
    !we aren't actually reading anything, but they don't need to know that ;)"
  ELSE
    WRITE(*,*) "Constructing two electron integrals"
    CALL proc2e(bas,basinfo,atoms,options,fmem,nnuc,xyz,norb,set,setinfo,maxL)
  END IF


  CONTAINS
!===================================================================
!                       SUBROUTINES

!---------------------------------------------------------------------
! 		Processes the two electron integrals
!---------------------------------------------------------------------
  SUBROUTINE proc2e(bas,basinfo,atoms,options,fmem,nnuc,xyz,norb,set,setinfo,maxL)
    IMPLICIT NONE
    ! Values
    ! xyz	: 2D dp, array of nuclear positions
    ! atoms	: 1D int, array of which atom is which
    ! fmem	: dp, free memory left in MB
    ! nnuc	: int, number of nuclii
    ! norb	: int, number of orbitals in molecule
    ! bas	: 1D dp, basis for each atom: atom : orbital : [d,a]
    ! basinfo	: 1D int, array of basis information
    ! set	: 1D dp, array of exponential coefficients
    ! setinfo	: 1D int, array of information about sets
    ! options	: 1D int, array of options
    ! II	: 3D dp, intermediate array 
    ! XX	: 4D dp, two electron repulsion integrals
    ! maxL	: int, max principle,angular quantum numbers
 
    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: bas,set
    INTEGER, DIMENSION(0:), INTENT(IN) :: basinfo,setinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: options, atoms
    REAL(KIND=8), INTENT(INOUT) :: fmem
    INTEGER, INTENT(IN) :: norb,nnuc,maxL

    !Internal
    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: XX
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: II
    REAL(KIND=8) :: timeS, timeF, temp
    INTEGER :: nset,setl,OpS,stat1,stat2 
    INTEGER :: a,b,c,d,g,h,i,j

    CALL CPU_TIME(timeS)

998 FORMAT(1x,A47,2x,F8.5)
999 FORMAT(1x,A41,2x,F8.4)
997 FORMAT(4x,A22)
996 FORMAT(4x,I3,1x,I3,1x,I3,1x,I3,4x,F15.8)
    
    nset = setinfo(0)
    setl = setinfo(1)
    OpS = basinfo(0)

    !Iuv will be my intermediate file for the integrals
    OPEN(unit=42,file='Iuv',status='replace',access='sequential',form='unformatted') 
 
    !allocate memory
    temp = (nset**(4.0D0))*8.0/1.0E6 
    temp = temp + (2*maxL)**(3)*8.0D0/1.0E6
    WRITE(*,998) "Allocating space for intermediate matrices (MB)", temp 
    ALLOCATE(XX(0:nset-1,0:nset-1,0:nset-1,0:nset-1),STAT=stat1)
    ALLOCATE(II(0:(2*maxL)**(3),0:nset-1,0:nset-1),STAT=stat2)
    IF (stat1+stat2 .NE. 0) THEN
      WRITE(*,*) "Could not allocate memory in int2e:proc2e"
      CALL EXECUTE_COMMAND_LINE('touch error')
      RETURN
    END IF
    fmem = fmem - temp
    CALL nmem(fmem)

    !zero XX
    DO i=0,nset-1
      DO j=0,nset-1
        DO g=0, nset-1
          XX(i,j,g,:) = (/ (0.0D0, h=0, nset-1) /) 
        END DO
      END DO
    END DO

    !"sum" over a,b
    DO a=0,nset-1
      DO b=0,nset-1      

        !zero II
        DO i=0,nset-1
          DO j=0,nset-1
            II(i,j,:) = (/ (0.0D0, g=0, nset-1) /)
          END DO
        END DO

        !"sum" over c,d
        DO c=0,nset-1
          DO d=0,nset-1
          
          END DO                                   !end loop over d
        END DO                                     !end loop over c

        !"loop" over g,h
        DO g=0,setl-1
          DO h=0,setl-1
            !"loop" over i,j
            DO i=0,setl-1
              DO j=0,setl-1
                XX(i,j,g,h) = XX(i,j,g,h) + 1.0D0
              END DO
            END DO
          END DO
        END DO

      END DO                                       !end loop over b
    END DO                                         !end loop over a

    WRITE(42) XX 
    CLOSE(unit=42,status='keep')

    IF (options(7) .GE. 3) THEN 
      WRITE(*,*) "Two electron Integrals by Set" 
      WRITE(*,997) "I   J   G   H    Value"   
      DO i=0,nset-1
        DO j=0,nset-1
          DO g=0,nset-1
            DO h=0,nset-1
              WRITE(*,996) i,j,g,h,XX(i,j,g,h)
            END DO 
          END DO
        END DO
      END DO
      
    END IF

    DEALLOCATE(XX)
    DEALLOCATE(II)

    !reassign memory
    fmem = fmem + (nset**(4.0D0)+nset**(3.0D0))*8.0/1.0E6
    CALL nmem(fmem)

    CALL CPU_TIME(timeF) 
    WRITE(*,*) 
    WRITE(*,999) "Two electron integrals constructed in (s)", (timeF-timeS)

  END SUBROUTINE proc2e

!===================================================================
!                       FUNCTIONS
!----------



END PROGRAM int2e
