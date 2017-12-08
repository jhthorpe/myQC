!///////////////////////////////////////////////////////////////////
!//            Construct basis set and primatives for myQC 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//
!// - WORK NOTE : Current supported methods/atoms. STO-3G. H, He 
!///////////////////////////////////////////////////////////////////

MODULE basis
  IMPLICIT NONE

  REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931 
  CONTAINS

!~~~~~~~~~~
! alters the leading coefficients of the primatives 
! WORK NOTE:  Currently this constructs the atomic orbitals... I'm not sure if I want to keep it this way
!
!  SUBROUTINE adjBasis(bkey,atoms,B)
!  !SUBROUTINE adjBasis(bkey,atoms,B,xyz,AO)
!
!    !Inout variables
!    !REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT) :: AO
!    !REAL(KIND=8),DIMENSION(0:),INTENT(IN):: xyz
!    REAL(KIND=8),DIMENSION(:,:,:),INTENT(IN) :: B
!    INTEGER,DIMENSION(0:),INTENT(IN) :: atoms
!    INTEGER,INTENT(IN) :: bkey 
!    
!    !Internal variables
!    REAL(KIND=8) :: a,d,val
!    INTEGER :: i,j,k,nnuc
!
!    nnuc=SIZE(atoms)
!     
!    ! One for each basis set
!    IF (bkey .EQ. 1) THEN !STO-3G
!    !  ALLOCATE(AO(0:nnuc-1,0:4))
!      DO i=0,nnuc-1
!        WRITE(*,*) "Atom # ", i
!
!        !1S orbital
!        val = 0.0D0
!        DO j=0,5,2 !loops over 1s
!          d = B(i,0,j)
!          a = B(i,0,j+1) 
!          val = val + d*(2.0D0*a/Pi)**(0.75D0)*EXP(-a*(xyz(0)**2.0D0+xyz(1)**2.0D0+xyz(2)**2.0D0))
!        END DO      !loop over 1s
!        AO(i,0) = val
!        WRITE(*,*) "1S(xyz) = ", val 
!
!        !2S orbital
!        val = 0.0D0
!        DO j=0,5,2 !loops over 2s
!          d = B(i,1,j)
!          a = B(i,1,j+1) 
!          val = val + d*(2.0D0*a/Pi)**(0.75D0)*EXP(-a*(xyz(0)**2.0D0+xyz(1)**2.0D0+xyz(2)**2.0D0))
!        END DO      !loop over 2s
!        AO(i,1) = val
!        WRITE(*,*) "2S(xyz) = ", val 
!
!        !2px,y,z orbital
!        DO k=0,2 
!          val = 0.0D0
!          DO j=0,5,2 !loops over 2p
!            d = B(i,2,j)
!            a = B(i,2,j+1) 
!            val = val + d*(128.0D0*a**5.0D0/Pi)**(0.25D0)*xyz(k)*EXP(-a*(xyz(0)**2.0D0+xyz(1)**2.0D0+xyz(2)**2.0D0))
!          END DO      !loop over 2p
!          WRITE(*,*) "2P?(xyz) ", val
!          AO(i,2+k-1) = val
!        END DO
!        WRITE(*,*)
!      END DO
!
!    ELSE
!      WRITE(*,*) "Sorry, that basis set has not been implimented yet. Exiting..."
!      STOP 'bad basis'
!
!    END IF
!
!    WRITE(*,*) "evalAO completed..."
!  END SUBROUTINE adjBasis
!~~~~~~~~~~
! reads basis set from mybasis 
  SUBROUTINE buildBasis(bkey,atoms,B,basinfo)

    ! WORK NOTE : think about condensing this into a 2D array, and iterating through the vectors using informatino like # of functions and # coefficients per function 

    ! Basis is a 3D array, 1st index is each atom, 2nd index is each orbital section, 3rd index are the individual values within the section in a linear array. They must be iterated through differently for each basis set 

    !Input variables
    ! bkey	: int representation of basis set
    ! atoms	: int 1D array of atoms
    ! B		: dp 3D array of coefficients
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT) :: B
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: basinfo
    INTEGER,DIMENSION(0:),INTENT(IN) :: atoms
    INTEGER,INTENT(IN) :: bkey 

    !Internal Variables
    ! A		: char 1D array of atom names
    ! C		: char 1D array of basis names
    ! Aname	: chr rep of atom
    REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: val
    CHARACTER(LEN=2),DIMENSION(0:9) :: A
    CHARACTER(LEN=8),DIMENSION(0:2) :: C 
    CHARACTER(LEN=8) :: line
    CHARACTER(LEN=2) :: Aname
    INTEGER :: i,j,k,l,m,Anum,Smax,Cmax,func,coef,sec,ang,pri,orb,Omax,ori,temp,nori

    A = ['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne']
    C = ['STO-3G ', 'tester1','tester2']
    Anum = SIZE(atoms)

    WRITE(*,*) "getbasis called..."

    !Open basis file
    OPEN (unit=2,file="mybasis",status="old",access="sequential")
    OPEN (unit=3,file='basinfo',status='replace',access='sequential')
    
    !This is absolutely disgusting code, I will fix later
    
    ! Go through each atom
    DO i=0,Anum-1 
      !Get to the starting point of the basis set
      DO WHILE (line .NE. C(bkey))
        READ(2,*) line
      END DO
      
      !setup basis array from line below title, if first run
      IF (i .EQ. 0) THEN
        READ(2,*) Smax, Cmax, Omax !max number of sections (orbital types), max constants, max orbitals
        ALLOCATE(B(0:Anum-1,0:(Omax-1),0:(Cmax-1)))
        ALLOCATE(basinfo(0:Anum-1,0:4*Omax+3))
        DO j=0,Anum-1
          DO k=0,Omax-1
            DO m=0,Cmax-1
              B(j,k,m) = 0.0D0
            END DO
          END DO
        END DO
      END IF

      !go to atom
      Aname = A(atoms(i)-1)
      DO WHILE (line .NE. Aname)
        READ(2,*) line
      END DO

      !add atom into basinfo
      basinfo(i,0) = atoms(i)

      !We are at atom, get basis stuff      
      READ(2,*) sec, orb !number of sections, number of orbitals

      !add orbitals of atom into basinfo
      basinfo(i,1) = sec
      basinfo(i,2) = orb
      basinfo(i,3) = -1

      temp = 0
      DO j=0,sec-1 !go through each section
        READ(2,*) func, coef, pri, ang, ori !number of functions(primatives), number of coefficients,angular quantum number, number of orientations (-1 spherical, 2 xyz ...)
  
        ! Dealing with x,y,z orbitals in basinfo
        IF (ori .EQ. -1) THEN
          basinfo(i,4*(j+1):4*(j+1)+3) = [pri, ang, ori, func]
          nori = 1
        ELSE IF (ori .EQ. 2) THEN
          DO k=0,ori
            basinfo(i,4*(j+1)+4*k:4*(j+1)+3+4*k) = [pri, ang, k, func] !update basinfo
          END DO
          nori = 3
        ELSE 
          WRITE(*,*) "that angular qunatum number not implimented yet in basis.f90:getbasis."
          WRITE(*,*) "probably also need to re-write coefficient program so that it works with d-orbitals"
          STOP "bad angular quantum number"
        END IF

        ! insert values of primatives info bas
        ALLOCATE(val(0:coef-1)) 
        DO k=0,func-1 !go through each primative
          READ(2,*) val
          DO m=0,coef-1 !assign linear array values
            DO l=0,nori-1
              B(i,j+l,2*k+m) = val(m) !THIS IS WRONG - is it?
              temp = temp + 1
            END DO
          END DO ! m loop (coefficients) 
        END DO ! k loop (primatives)
        DEALLOCATE(val)
      END DO ! j loop (orbitals)

      !cleanup at end of basinfo
      DO j=orb,Omax-1
        basinfo(i,4*(j+1):4*(j+1)+3) = [0, 0, 0, 0]
      END DO

      !write to basinfo
      WRITE(3,*) basinfo(i,:)

      !set up for next atom
      REWIND (2)
      READ(2,*) line
    END DO !i loop (atoms)

    WRITE(3,*)
    WRITE(3,*) "atom, #sections, #orbitals, dummy, {principle quantum number, angular quantum number, orientation, #primatives}"

    WRITE(3,*) 
    WRITE(3,*) 
    DO i=0,Anum-1
      WRITE(3,*) "atom #:", i
      DO j=0,basinfo(i,2)-1
        WRITE(3,*) B(i,j,:) 
      END DO
    END DO
    WRITE(3,*) 
 
    CLOSE (unit=2,status="keep")
    CLOSE (unit=3,status='keep')

    

  END SUBROUTINE buildBasis

!~~~~~~~~~~

END MODULE basis 
