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
  SUBROUTINE buildBasis(bkey,atoms,bas,basinfo,set,setinfo,verb,maxN,maxL)

    ! Basis is a 3D array, 1st index is each atom, 2nd index is each orbital section, 3rd index are the individual values within the section in a linear array. They must be iterated through differently for each basis set 

    ! Basis is now ordered by [atom, setnum, {alpha, #orbs, w0,w1,w2,...}]

    ! Values
    ! bkey	: int representation of basis set
    ! atoms	: int 1D array of atoms
    ! B		: dp 3D array of coefficients
    ! set	: 2D dp, array of exponential coefficients
    ! basinfo	: 2D int, array of basis info
    ! setinfo	: 2D int, array of set info
    ! verb	: bool, flag for if need to print basis set output
    ! A		: char 1D array of atom names
    ! C		: char 1D array of basis names
    ! Aname	: chr rep of atom
    ! almax	: int, max number of exp coef in set
    ! off	: int, tracks offset of set values
    ! Smax	: int, max sections in basis
    ! Cmax	: int, max coefficients in section of bas
    ! Omax	: int, max orbitals in atom of basis
    ! setl	: int, number of units in one section of set
    ! OpS	: int, max orbs per set
    ! maxN	: int, max principle quantum number
    ! maxL	: int, max angular quantum number

    ! INOUT
    REAL(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT) :: bas,set
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: basinfo,setinfo
    INTEGER, DIMENSION(0:),INTENT(IN) :: atoms
    INTEGER, INTENT(INOUT) :: maxN, maxL
    INTEGER, INTENT(IN) :: bkey 
    LOGICAL, INTENT(IN) :: verb

    !Internal 
    REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: val
    LOGICAL, DIMENSION(:), ALLOCATABLE :: ttab
    CHARACTER(LEN=2),DIMENSION(0:9) :: A
    CHARACTER(LEN=8),DIMENSION(0:3) :: C 
    CHARACTER(LEN=8) :: line
    CHARACTER(LEN=2) :: Aname
    REAL(KIND=8) :: temp
    INTEGER :: i,j,k,l,m,n,Anum,Smax,Cmax,func,coef,sec,ang,pri,orb, OpS
    INTEGER :: Omax,ori,nori,almax,nset,setn,setl,setorbs,orbnum,lmax,setnum

    A = ['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne']
    C = ['STO-3G ', 'tester1','tester2','tester3']
    Anum = SIZE(atoms)

999 FORMAT(2x,A8,2x,I2)
998 FORMAT(2x,A6,2x,A2)
997 FORMAT(4x,A4)
996 FORMAT(4x,I2,1x,I2)
995 FORMAT(4x,F15.8,1x,F15.8)

    IF (verb) WRITE(*,*) "getbasis called..."

    !Open basis file
    OPEN (unit=2,file="mybasis",status="old",access="sequential")
    OPEN (unit=3,file='basinfo',status='replace',access='sequential')
    OPEN (unit=4,file='setinfo',status='replace',access='sequential')
    
    !tracking variables
    setnum = 0                           !number of sets seen so far
    orbnum = 0                           !number of orbitals within set

    IF (verb) WRITE(*,*) "basis:",C(bkey)
    
    ! Go through each atom
    DO i=0,Anum-1 
      !Get to the starting point of the basis set
      DO WHILE (line .NE. C(bkey))
        READ(2,*) line
      END DO
      
      ! If first run, setup the basis set
      IF (i .EQ. 0) THEN
        READ(2,*) Smax, Cmax, Omax, almax, OpS 
        READ(2,*) maxN, maxL
        ALLOCATE(bas(0:(Anum*almax*OpS)-1)) 
        ALLOCATE(basinfo(0:2+(5*Omax*Anum)-1))
        ALLOCATE(set(0:(Anum*almax)-1))
        ALLOCATE(setinfo(0:2+(Anum*almax*(3+OpS))-1))
        ALLOCATE(ttab(0:almax-1))
        ! Zero the arrays
        bas(:) = (/ (0.0D0, j=0, Anum*almax*OpS-1) /)
        set(:) = (/ (0.0D0, j=0, Anum*almax-1) /)
        basinfo(:) = (/ (0, j=0, 2+(5*Omax*Anum)-1) /)
        setinfo(:) = (/ (0, j=0, 2+(Anum*almax*(3+OpS))-1) /)
      END IF

      ! go to atom
      Aname = A(atoms(i)-1)
      DO WHILE (line .NE. Aname)
        READ(2,*) line
      END DO

      IF (verb) WRITE(*,999) "Nuclei #", i+1
      IF (verb) WRITE(*,998) "Atom :", Aname
      IF (verb) WRITE(*,997) "n l"
      
      !We are at atom, adjust basis stuff      
      READ(2,*) sec, orb, nset             !#sections, #orbitals, #sets
      basinfo(0) = OpS                    !number of atoms
      basinfo(1) = basinfo(1) + orb        !number of orbitals
      setinfo(0) = setinfo(0) + nset       !number of sets
      setinfo(1) = 3 + OpS                !length of each set
      setl = 3 + OpS

      !go through each section (orbital)
      DO j=0,sec-1 

        !#primatives, #coefficients, angular qn, #orientations (-1 spherical, 2 xyz ...), # of new sets
        READ(2,*) func, coef, pri, ang, ori 
        ALLOCATE(val(0:coef-1))                        ! structure input array 

        IF (verb) WRITE(*,996) pri, ang
        
        !Primatives - Set updates
        DO k=0,func-1       
          READ(2,*) val, temp                          !coefs, setinfo

          !update set info
          setn = NINT(temp)                            !ID of set 
          set(setnum + setn) = val(coef - 1)           !exp coef of set
          setorbs = setinfo(1+setn*setl+1+setnum*setl) !#orbs in set so far
          setinfo(1+(setnum+setn)*setl+3) = i

          IF (verb) WRITE(*,995) val 

          !S-TYPE
          IF (ori .EQ. -1) THEN
            setinfo(1+setn*setl+4+setorbs+setnum*setl) = orbnum  ! update set orbital count
            bas(setorbs+(setn+setnum)*OpS) = val(0)             !add coefficient to basis
            setorbs = setorbs + 1 

          !P-TYPE
          ELSE IF (ori .EQ. 2) THEN
            DO m=0,2
              setinfo(1+(setn+setnum)*setl+4+setorbs) = orbnum + m
              bas(setorbs+(setn+setnum)*OpS) = val(0)
              setorbs = setorbs + 1
            END DO
            !check if we need to update lmax
            IF (setinfo(1+setn*setl+2+setnum*setl) .LT. 1) setinfo(1+setn*setl+2+setnum*setl) = 1
 
          !D-TYPE
          ELSE
             WRITE(*,*) "that angular qunatum number not implimented yet in basis.f90:getbasis."
             WRITE(*,*) "probably also need to re-write coefficient program so that it works with d-orbitals"
             STOP "bad angular quantum number"
          END IF    

          setinfo(1+(setn+setnum)*setl+1) = setorbs

        END DO                             ! k loop (primatives)

        !Orbitals - basis updates
        ! S-TYPE
        IF (ori .EQ. -1) THEN
          basinfo(2+5*orbnum:2+5*(orbnum+1)-1) = [pri, ang, ori, func,i]
          orbnum = orbnum + 1

        ! P-TYPE
        ELSE IF (ori .EQ. 2) THEN
          DO m=0,2
            basinfo(2+5*(orbnum+m):2+5*(orbnum+m+1)-1) = [pri, ang, m, func,i]
          END DO
          orbnum = orbnum + 3

        ! D-TYPE
        ELSE
          WRITE(*,*) "that angular qunatum number not implimented yet in basis.f90:getbasis."
          WRITE(*,*) "probably also need to re-write coefficient program so that it works with d-orbitals"
          STOP "bad angular quantum number"
        END IF

        DEALLOCATE(val)

      END DO                               ! j loop (orbitals)

      !set up data for next atom
      setnum = setinfo(0) 

      !set up for next atom
      REWIND (2)
      READ(2,*) line
    END DO                                 !end atom loop
     
    ! write output
    WRITE(3,*) basinfo(:)
    WRITE(3,*)
    WRITE(3,*) "#orbitals per set, #orbitals, {principle qn., angular qn., orientation, #primatives, center number},... "
    WRITE(4,*) setinfo(:)
    WRITE(4,*) 
    WRITE(4,*) "#sets, length of each set,{#orbitals, max ang qn., center num, [orbital 0, orbital 1, ...]},..."

    CLOSE (unit=4,status='keep')
    CLOSE (unit=3,status='keep')
    CLOSE (unit=2,status="keep")

    IF (verb) WRITE(*,*) "=================================="

  END SUBROUTINE buildBasis

END MODULE basis 
