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

!=====================================================================
!                       SUBROUTINES

!---------------------------------------------------------------------
!		Construct Basis, Set, basinfo, and setinfo
!---------------------------------------------------------------------
  SUBROUTINE buildBasis(bkey,atoms,B,basinfo,set,setinfo)

    ! Basis is a 3D array, 1st index is each atom, 2nd index is each orbital section, 3rd index are the individual values within the section in a linear array. They must be iterated through differently for each basis set 

    ! Basis is now ordered by [atom, setnum, {alpha, #orbs, w0,w1,w2,...}]

    ! Values
    ! bkey	: int representation of basis set
    ! atoms	: int 1D array of atoms
    ! B		: dp 3D array of coefficients
    ! set	: 2D dp, array of exponential coefficients
    ! basinfo	: 2D int, array of basis info
    ! setinfo	: 2D int, array of set info
    ! A		: char 1D array of atom names
    ! C		: char 1D array of basis names
    ! Aname	: chr rep of atom
    ! almax	: int, max number of exp coef in set
    ! off	: int, tracks offset of set values

    ! INOUT
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT) :: B
    REAL(KIND=8), DIMENSION(:,:),ALLOCATABLE, INTENT(INOUT) :: set 
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: basinfo,setinfo
    INTEGER,DIMENSION(0:),INTENT(IN) :: atoms
    INTEGER,INTENT(IN) :: bkey 

    !Internal 
    REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: val
    LOGICAL, DIMENSION(:), ALLOCATABLE :: ttab
    CHARACTER(LEN=2),DIMENSION(0:9) :: A
    CHARACTER(LEN=8),DIMENSION(0:3) :: C 
    CHARACTER(LEN=8) :: line
    CHARACTER(LEN=2) :: Aname
    REAL(KIND=8) :: temp
    INTEGER :: i,j,k,l,m,n,Anum,Smax,Cmax,func,coef,sec,ang,pri,orb
    INTEGER :: Omax,ori,nori,almax,nset,setn,setorbs,orbnum,lmax,off

    A = ['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne']
    C = ['STO-3G ', 'tester1','tester2','tester3']
    Anum = SIZE(atoms)

    WRITE(*,*) "getbasis called..."

    !Open basis file
    OPEN (unit=2,file="mybasis",status="old",access="sequential")
    OPEN (unit=3,file='basinfo',status='replace',access='sequential')
    OPEN (unit=4,file='setinfo',status='replace',access='sequential')
    
    !This is absolutely disgusting code, I will fix later

    off = 0 
    
    ! Go through each atom
    DO i=0,Anum-1 
      !Get to the starting point of the basis set
      DO WHILE (line .NE. C(bkey))
        READ(2,*) line
      END DO
      
      !setup basis array from line below title, if first run
      IF (i .EQ. 0) THEN
        READ(2,*) Smax, Cmax, Omax, almax !max number of sections (orbital types), max constants, max orbitals
!        ALLOCATE(B(0:Anum-1,0:(Omax-1),0:(Cmax-1)))
        ALLOCATE(B(0:Anum-1,0:almax-1,0:Omax+1))
        ALLOCATE(basinfo(0:Anum-1,0:4*Omax+3))
        ALLOCATE(set(0:Anum-1,0:almax-1))
        ALLOCATE(setinfo(0:Anum-1,0:almax*(2+Omax)+2))
        ALLOCATE(ttab(0:almax-1))

        DO j=0,Anum-1
          ! zero bas
          DO k=0,almax-1
            DO m=0,Omax-1
              B(j,k,m) = 0.0D0
            END DO
          END DO
          ! zero set,setinfo, and offset
          set(j,:) = (/ (0.0D0, k=0,almax-1) /)
          setinfo(j,:) = (/ (0, k=0,almax*(2+Omax)+2) /)
          basinfo(j,:) = (/ (0, k=0,4*Omax+3) /)
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
      READ(2,*) sec, orb, nset                             !# sections, # orbitals, # sets

      !add orbitals of atom into basinfo
      basinfo(i,1) = sec
      basinfo(i,2) = orb
      basinfo(i,3) = -1

      !update setinfo
      setinfo(i,0) = nset                                  !number of sets
      setinfo(i,1) = orb                                   !number of orbitals
      setinfo(i,2) = Omax+2                                !length of each set

      !tracking variables
      orbnum = 0

      !zero ttab
      DO k=0,almax-1
        ttab(k) = .FALSE. 
      END DO

      !go through each section
      DO j=0,sec-1 
        READ(2,*) func, coef, pri, ang, ori !number of functions(primatives), number of coefficients,angular quantum number, number of orientations (-1 spherical, 2 xyz ...), # of new sets
  
          ! insert values of basis 
        ALLOCATE(val(0:coef-1)) 

        !go through each primative
        DO k=0,func-1 
          READ(2,*) val, temp

          !update set info
          setn = NINT(temp) 
          set(i,setn) = val(coef-1)                        !add exp coef to set 
          setorbs = setinfo(i,2+setn*(Omax+2)+1)           !get #orb in set

          !deal with orientation for set
          ! S-TYPE
          IF (ori .EQ. -1) THEN 
            setorbs = setorbs + 1 
            setinfo(i,2+setn*(Omax+2)+2+setorbs) = orbnum  ! update set orbital count
            B(i,setn,0) = val(coef-1)                      ! add exp coef to bas  
            B(i,setn,1) = setorbs                           ! add number of orbitals to set
            B(i,setn,1+setorbs) = val(0)                    ! WORK NOTE- HARDCODED 

          ! P-TYPE
          ELSE IF (ori .EQ. 2) THEN 
            DO m=0,2   
              setorbs = setorbs + 1
              setinfo(i,2+setn*(Omax+2)+2+setorbs) = orbnum
              B(i,setn,0) = val(coef-1)                    ! add exp coef to bas  
              B(i,setn,1) = setorbs                        ! add number of orbitals to set
              B(i,setn,1+setorbs) = val(0)                  ! WORK NOTE- HARDCODED 
              orbnum = orbnum + 1 
            END DO
            orbnum = orbnum - 3  !reset orbnumber to keep with basinfo
            !check if we need to update lmax
            IF (setinfo(i,2+setn*(Omax+2)+2) .LT. 1) setinfo(i,2+setn*(Omax+2)+2) = 1

          ! D-TPYE
          ELSE 
            WRITE(*,*) "that angular qunatum number not implimented yet in basis.f90:getbasis."
            WRITE(*,*) "probably also need to re-write coefficient program so that it works with d-orbitals"
            STOP "bad angular quantum number"
          END IF

          ! update number of orbs in set
          setinfo(i,setn*(Omax+2)+3) = setorbs

        END DO ! k loop (primatives)
 
        orbnum = orbnum + 1

        DEALLOCATE(val)

        ! Dealing with x,y,z orbitals
        IF (ori .EQ. -1) THEN
          basinfo(i,4*(j+1):4*(j+1)+3) = [pri, ang, ori, func]
          nori = 1
          off = off + 1
        ELSE IF (ori .EQ. 2) THEN
          DO k=0,ori
            basinfo(i,4*(j+1)+4*k:4*(j+1)+3+4*k) = [pri, ang, k, func] !update basinfo
            off = off + 1
          END DO
          nori = 3
        ELSE 
          WRITE(*,*) "that angular qunatum number not implimented yet in basis.f90:getbasis."
          WRITE(*,*) "probably also need to re-write coefficient program so that it works with d-orbitals"
          STOP "bad angular quantum number"
        END IF

      END DO ! j loop (orbitals)

      !cleanup at end
      DO j=orb,Omax-1
        basinfo(i,4*(j+1):4*(j+1)+3) = [0, 0, 0, 0]
      END DO

      !write to basinfo/setinfo/offset
      WRITE(3,*) basinfo(i,:)
      WRITE(4,*) setinfo(i,:)

      !set up for next atom
      REWIND (2)
      READ(2,*) line
    END DO !i loop (atoms)

    WRITE(3,*)
    WRITE(3,*) "atom, #sections, #orbitals, dummy, {principle quantum number, angular quantum number, orientation, #primatives}"

    WRITE(4,*) 
    WRITE(4,*) 
    WRITE(4,*) "#sets, #orbitals, setlength, {#orbs in set, max ang qn, orb1, orb2... maxorb}," 
    !write to bottom of basinfo
!    DO i=0,Anum-1
!      WRITE(3,*) "atom #:", i
!      DO j=0,basinfo(i,2)-1
!        WRITE(3,*) B(i,j,:) 
!      END DO
!    END DO
!    WRITE(3,*) 

 
!    CLOSE (unit=5,status='keep')
    CLOSE (unit=4,status='keep')
    CLOSE (unit=3,status='keep')
    CLOSE (unit=2,status="keep")

  END SUBROUTINE buildBasis

END MODULE basis 
