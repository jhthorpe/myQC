!//////////////////////////////////////////////////////////////////
!//            Performs 1 e- integrals for myQC 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//		
!//		Integration scheme a la McMurchie,Davidson 1978	
!//		
!//	WORK NOTE - location of coefficients in bas currently hardcoded	
!///////////////////////////////////////////////////////////////////

PROGRAM int1e
  USE env
  USE basis
  IMPLICIT NONE

  ! Values
  ! xyz		: 2D dp, array of nuclear positions
  ! atoms	: 1D int, array of which atom is which
  ! fmem	: dp, free memory left in MB
  ! nnuc	: int, number of nuclii
  ! nelc	: int, number of electrons
  ! norb	: int, number of orbitals in molecule
  ! npri	: int, number of primatives
  ! bas		: 2D dp, basis for each atom: atom : orbital : [d,a]
  ! basinfo	: 2D int, array of basis information
  ! options	: 1D int, array of options
  ! S		: 2D dp, overlap matrix

  ! Variables  
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: bas
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz,S
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: basinfo
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms,options 
  REAL(KIND=8) :: timeS, timeF, fmem
  INTEGER :: nnuc,nelc,i,j,k,norb,npri,stat
  LOGICAL :: flag

! input managment 
  CALL CPU_TIME(timeS)
  WRITE(*,*)
  WRITE(*,*) "int1e called"
  CALL getenv(nnuc,nelc,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

!the actual stuff
! construct the basis
  CALL buildBasis(options(2),atoms,bas,basinfo)
  ! write out basis set for checking, include in verbosity later
  ! WORK NOTE - maybe add in memory of basis set here?
  WRITE(*,*)
  WRITE(*,*) "Basis set"
  DO i=0,nnuc-1
    WRITE(*,*) "Nuclei #", i+1
    WRITE(*,*) "atom : ", atoms(i)
    DO j=1,basinfo(i,1)
      WRITE(*,*) "orbital (n,l) : ", basinfo(i,4*j), basinfo(i,4*j+1)
      DO k=0,basinfo(i,4*j+3)-1
        WRITE(*,*) bas(i,j-1,2*k), bas(i,j-1,2*k+1) 
      END DO
    END DO
  END DO
  WRITE(*,*) "=================================="
  CALL nmem(fmem)

! 1) calculate overlap matrix
  !get number of orbitals
  npri = 0
  norb = 0
  DO i=0,nnuc-1
    norb = norb + basinfo(i,2) 
    DO j=1,basinfo(i,2)
      npri = npri + basinfo(i,4*(j)+3) 
    END DO
  END DO
  WRITE(*,*) 
  WRITE(*,*) "Number of orbitals : ", norb
  WRITE(*,*) "Number of primatives : ", npri
  WRITE(*,*)
  WRITE(*,*) "Allocating space for overlap matrix (MB) : ", norb*norb*8/1.0E6
  ALLOCATE(S(0:norb-1,0:norb-1),STAT=stat)
  IF (stat .NE. 0) STOP "int1e: max memory reached, exiting"
  fmem = fmem - norb*norb*8/1.0E6
  CALL nmem(fmem)
  CALL overlap(S,bas,basinfo,atoms,options,fmem,nnuc,xyz,norb)

! 2) calculate kinetic energy integrals

! 3) calculate coulombic integrals 

! output
  CALL setenv(atoms,xyz,fmem,options)
  CALL CPU_TIME(timeF)
  WRITE(*,*) "int1e ran in (s) : ", (timeF - timeS) 

  CONTAINS 

!~~~~~
  ! calculate overlap matrix
  SUBROUTINE overlap(S,bas,basinfo,atoms,options,fmem,nnuc,xyz,norb)
    IMPLICIT NONE
    ! Values
    ! S		: 2D dp, overlap matrix
    ! xyz	: 2D dp, array of nuclear positions
    ! atoms	: 1D int, array of which atom is which
    ! fmem	: dp, free memory left in MB
    ! nnuc	: int, number of nuclii
    ! norb	: int, number of orbitals in molecule
    ! bas	: 2D dp, basis for each atom: atom : orbital : [d,a]
    ! basinfo	: 2D int, array of basis information
    ! options	: 1D int, array of options

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: S 
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: options, atoms
    REAL(KIND=8), INTENT(IN) :: fmem
    INTEGER, INTENT(IN) :: norb,nnuc

    !Internal
    REAL(KIND=8) :: timeS, timeF
    INTEGER :: u,v,col,row,col0,row0

    CALL CPU_TIME(timeS)
    WRITE(*,*) "constructing overlap matrix"

    !Matrix blocked by atoms
    col0 = 0
    DO u=0,nnuc-1
      row0 = 0
      DO v=0,nnuc-1 
        !determine length, width, and starting points of block
        col = basinfo(u,2) !number of orbitals in a
        row = basinfo(v,2) ! number of orbitals in b
        !construct overlap within block
        CALL Sblock(S(col0:col0+col-1,row0:row0+row-1),bas,basinfo,atoms,xyz,u,v,col,row)
        row0 = row0 + row
      END DO
      col0 = col0 + col
    END DO 

    CALL CPU_TIME(timeF)
    WRITE(*,*) "overlap matrix constructed in (s) :", (timeF-timeS)
  END SUBROUTINE overlap

!~~~~~
  ! calculate block of overlap matrix 
  SUBROUTINE Sblock(Sb,bas,basinfo,atoms,xyz,u,v,na,nb)
    IMPLICIT NONE
    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931

    ! Values
    ! Sb	: 2D dp, subblock of overlap matrix
    ! xyz	: 2D dp, array of nuclear positions
    ! atoms	: 1D int, array of which atom is which
    ! bas	: 3D dp, basis for each atom: atom : orbital : [d,a]
    ! basinfo	: 2D int, array of basis information
    ! u,v	: int, block index of S (col,row) (atom A, atom B)
    ! na, nb	: int, number of orbitals in atom a,b
    ! coef	: 3D dp, array of coefficients for overlap, (xyz, j,k)
    ! PA, PB	: 1D dp, array of molecular locations
    ! EIJ	: dp, gaussian at molecular center
    ! val 	: dp, current evaluation of integral
    ! p		: dp, addition of two coeffiecients
    ! aa,bb	: dp, coefficients of Guassians for atoms a and b

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Sb 
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms
    INTEGER, INTENT(IN) :: u,v,na,nb

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: coef
    REAL(KIND=8), DIMENSION(0:2) :: PA, PB, AB, PP
    REAL(KIND=8) :: EIJ, val, p, m, aa, bb, temp
    INTEGER :: i,j,k,s,t,a,b
    
    WRITE(*,*) "u,v",u,v
    
    DO a=0,na-1 !iterate through orbitals of atom A
      DO b=0,nb-1 !iterate through orbitals of atom B
        WRITE(*,*) "a,b",a,b
        val = 0.0D0
        DO s=0,basinfo(u,4*(a+1)+3)-1 !iterate through primatives of orbital a 
          DO t=0,basinfo(v,4*(b+1)+3)-1 !iteratre thorugh primatives of orbital b

            WRITE(*,*) "s,t", s,t 
            ! 1) get overlap location
            aa = bas(u,a,s*2+1)! WORK NOTE hardcoded to 2nd value in bas of correct u and a 
            bb = bas(v,b,t*2+1)
            WRITE(*,*) "aa,bb", aa, bb
            p = aa + bb 
            m = aa * bb
            DO i=0,2
              AB(i) = xyz(u,i) - xyz(v,i)
              PP(i) = (aa*xyz(u,i) + bb*xyz(v,i))/p
              PA(i) = PP(i) - xyz(u,i)
              PB(i) = PP(i) - xyz(v,i)
            END DO 
            EIJ = EXP(-m*(AB(0)**2.0D0 + AB(1)**2.0D0 + AB(2)**2.0D0)/p) 
 
  ! updated to here

            ! 2) get coefficients for these primatives
            WRITE(*,*) "Sending PA, PB, aa, ab", PA(2),PB(2), aa,bb
            CALL getcoef(coef,PA,PB,aa,bb,basinfo(u,3*(a+1)+1),basinfo(v,3*(b+1)+1),basinfo(u,3*(a+1)+1)+basinfo(v,3*(b+1)+1))

            ! 3) add everything together
            ! precalculated constants
            temp = EIJ*(Pi/p)**(3.0D0/2.0D0)*bas(u,a,s*2)*bas(v,b,t*2) ! WORK NOTE - hardcoded
            ! integral coefficients
            temp = temp*coef(0,0,a,b)*coef(1,0,a,b)*coef(2,0,a,b)
            ! basis set coefficients 
            temp = temp * gtoD(basinfo(u,3*(a+1)+1),aa)
            temp = temp * gtoD(basinfo(v,3*(b+1)+1),bb)
           
            val = val + temp
            STOP "line 221 in integrals.f90. Not sure if px,py,pz are seperate"
            DEALLOCATE(coef)

          END DO !end t
        END DO !end s
        Sb(u,v) = val
        WRITE(*,*) "~~~~~~~~~~"
      END DO !end b
    END DO !end a

   WRITE(*,*) "~~~~~~~~~~~~~~"
   WRITE(*,*) "~~~~~~~~~~~~~~"
  END SUBROUTINE Sblock

!~~~~~
  !Subroutine to calculate coefficient matrix 
  SUBROUTINE getcoef(M,PA,PB,aa,bb,amax,bmax,nmax)
    IMPLICIT NONE
    ! M		: 2D dp, coefficient matrix d_{i}^{j,k} = M(0,i,j,k)
    ! PA	: 1D dp, (xa,ya,za) distance from central gaussian
    ! PB	: 1D dp, (xb,yb,zb) distance from central gaussian
    ! aa	: dp, coefficient of atom A
    ! bb	: dp, coefficient of atom B
    ! amax	: int, max angular quantum number of A needed
    ! amax	: int, max angular quantum number of B needed
    ! nmax	: int, max total quantum number needed
    ! nn	: int, total number of elements
    ! fmat	: 3D dp, boolean array, marks what has been calculated 

    ! Inout
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:), INTENT(INOUT) :: M
    REAL(KIND=8), DIMENSION(0:2), INTENT(IN) :: PA, PB   
    REAL(KIND=8), INTENT(IN) :: aa,bb
    INTEGER, INTENT(IN) :: amax, bmax, nmax
    
    ! Internal
    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: fmat
    REAL(KIND=8) :: pp
    INTEGER :: l,i,j,k,nn

    ! 1) construct initial stuff
    pp = aa+bb

    nn = 0
    DO i=0,nmax
      nn = nn + (i+1) 
    END DO
    WRITE(*,*) "Number of coefficients..., ", nn
    WRITE(*,*) "Coefficients are...", aa, bb
    WRITE(*,*) "nmax, amax, bmax", nmax, amax, bmax
    WRITE(*,*) "PAz, PBz", PA(2), PB(2) 
    ALLOCATE(M(0:2,-1:amax+bmax,-1:amax+1,-1:bmax+1)) 
    ALLOCATE(fmat(-1:amax+bmax,-1:amax+1,-1:bmax+1))
    !initialize fmat
    DO i=-1,amax+bmax
      DO j=-1,amax+1
        fmat(i,j,:) = (/ (.FALSE. , k=-1,bmax+1  ) /)
      END DO
    END DO
    DO l=0,2
      DO i=-1,amax+bmax
        DO j=-1,amax+1
          DO k=-1,bmax+1
            M(l,i,j,k) = 0.0D0
          END DO
        END DO
      END DO
    END DO 

    ! 2) call recursive function on top
    !test call
    CALL lrec(M,0,amax,bmax,PA,PB,pp,fmat)
    CALL rrec(M,0,amax,bmax,PA,PB,pp,fmat)
    
    WRITE(*,*) "===================="
    WRITE(*,*) "0,0,0", M(2,0,0,0) 
    WRITE(*,*) "0,0,1", M(2,0,0,1) 
    WRITE(*,*) "0,1,0", M(2,0,1,0) 
    !WRITE(*,*) "0,0,2", M(2,0,0,2)
    WRITE(*,*) "0,1,1", M(2,0,1,1)
    !WRITE(*,*) "0,2,0", M(2,0,2,0)
    !WRITE(*,*) "0,1,2", M(2,0,1,2)
    !WRITE(*,*) "0,2,1", M(2,0,2,1)
    !WRITE(*,*) "0,2,2", M(2,0,2,2)
    WRITE(*,*) "===================="

    DEALLOCATE (fmat)

  END SUBROUTINE getcoef
!~~~~~
  !left side recursion on coefficients of overlap matrix
  RECURSIVE SUBROUTINE lrec(M,i,j,k,PA,PB,pp,fmat)
    IMPLICIT NONE
   
    ! M		: 4D dp, matrix of coefficients
    ! i,j,k	: int, index of coefficient we need. i=N,j=n,k=nbar
    ! PA, PB	: 1D dp, list of x,y,z distances of atoms A,B from P
    ! fmat	: 3D dp, list of which values we already have
    ! pp	: dp, value of sum of coefficients

    ! INOUT
    REAL(KIND=8), DIMENSION(0:,-1:,-1:,-1:), INTENT(INOUT) :: M
    REAL(KIND=8), DIMENSION(0:2), INTENT(IN) :: PA, PB
    LOGICAL, DIMENSION(-1:,-1:,-1:),INTENT(INOUT) :: fmat
    REAL(KIND=8), INTENT(IN) :: pp
    INTEGER, INTENT(IN) :: i,j,k

    ! Internal
    INTEGER :: l
   
    WRITE(*,*) "lrec called with...", i,j,k
    
    ! 1) base cases
    ! if we've already seen this point
    ! if not, check that it isn't automatically zero
    IF (fmat(i,j,k)) THEN
      WRITE(*,*) "i,j,k  already seen"
      RETURN
    END IF
    IF (i .LT. 0 .OR. j .LT. 0 .OR. k .LT. 0 .OR. i .GT. j+k) THEN
      DO l=0,2
        M(l,i,j,k) = 0.0D0
      END DO
      fmat(i,j,k) = .TRUE.
      WRITE(*,*) "i,j,k is 0"
      RETURN
    ! check for d000 base case
    ELSE IF (i .EQ. 0 .AND. j .EQ. 0 .AND. k .EQ. 0) THEN
      DO l=0,2
        M(l,i,j,k) = 1.0D0
      END DO
      fmat(i,j,k) = .TRUE.
      WRITE(*,*) "hit d000"
      RETURN 
    END IF

    WRITE(*,*) "non-base case" 
    
    ! 2) recursion algorithm 
    ! if i is zero, we will need this in overlap
    IF (j .EQ. k) THEN 
      WRITE(*,*) "going l and r"
      CALL lrec(M,i-1,j-1,k,PA,PB,pp,fmat)
      CALL lrec(M,i,j-1,k,PA,PB,pp,fmat)
      CALL lrec(M,i+1,j-1,k,PA,PB,pp,fmat)
      CALL rrec(M,i-1,j,k-1,PA,PB,pp,fmat)
      CALL rrec(M,i,j,k-1,PA,PB,pp,fmat)
      CALL rrec(M,i+1,j,k-1,PA,PB,pp,fmat)
      DO l=0,2
        M(l,i,j,k) = M(l,i-1,j,k-1)/(2*pp) + PB(l)*M(l,i,j,k-1) + (i+1)*M(l,i+1,j,k-1)   
        fmat(i,j,k) = .TRUE. 
      END DO
      WRITE(*,*) "lvec-here0", i,j,k
      RETURN
    ! get the rest of what we need
    ! if k is bigger, lower it
    ELSE IF (j .LT. k) THEN
      WRITE(*,*) "going r"
      CALL rrec(M,i-1,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,i,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,i+1,j,k-1,PA,PB,pp,fmat) 
      DO l=0,2
        M(l,i,j,k) = M(l,i-1,j,k-1)/(2*pp) + PB(l)*M(l,i,j,k-1) + (i+1)*M(l,i+1,j,k-1)   
        fmat(i,j,k) = .TRUE. 
      END DO
      WRITE(*,*) "lvec-here1"
      RETURN
    ! otherwise, lower j (doesn't really matter which way you do it)
    ELSE IF (j .GE. k) THEN
      WRITE(*,*) "going l"
      CALL lrec(M,i-1,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,i,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,i+1,j-1,k,PA,PB,pp,fmat) 
      DO l=0,2
        M(l,i,j,k) = M(l,i-1,j-1,k)/(2*pp) + PA(l)*M(l,i,j-1,k) + (i+1)*M(l,i+1,j-1,k)   
        fmat(i,j,k) = .TRUE. 
      END DO
      WRITE(*,*) "lvec-here2"
      RETURN
    ELSE 
      WRITE(*,*) "Somehow you broke this at:", i,j,k
      STOP
    END IF

 
  END SUBROUTINE lrec
!~~~~~
  !right side (B) recursion on coefficients of overlap matrix
  RECURSIVE SUBROUTINE rrec(M,i,j,k,PA,PB,pp,fmat)
    IMPLICIT NONE
   
    ! M		: 4D dp, matrix of coefficients
    ! i,j,k	: int, index of coefficient we need. i=N,j=n,k=nbar
    ! PA, PB	: 1D dp, list of x,y,z distances of atoms A,B from P
    ! fmat	: 3D dp, list of which values we already have
    ! pp		: dp, value of sum of coefficients

    ! INOUT
    REAL(KIND=8), DIMENSION(0:,-1:,-1:,-1:), INTENT(INOUT) :: M
    REAL(KIND=8), DIMENSION(0:2), INTENT(IN) :: PA, PB
    LOGICAL, DIMENSION(-1:,-1:,-1:), INTENT(INOUT) :: fmat
    REAL(KIND=8), INTENT(IN) :: pp
    INTEGER, INTENT(IN) :: i,j,k

    ! Internal
    INTEGER :: l
   
    WRITE(*,*) "rrec called with...", i,j,k
    
    ! 1) base cases
    ! if we've already seen this point
    ! if not, check that it isn't automatically zero
    IF (fmat(i,j,k)) THEN
      WRITE(*,*) "i,j,k  already seen"
      RETURN
    END IF
    IF (i .LT. 0 .OR. j .LT. 0 .OR. k .LT. 0 .OR. i .GT. j+k) THEN
      DO l=0,2
        M(l,i,j,k) = 0.0D0
      END DO
      fmat(i,j,k) = .TRUE.
      WRITE(*,*) "i,j,k is 0"
      RETURN
    ! check for d000 base case
    ELSE IF (i .EQ. 0 .AND. j .EQ. 0 .AND. k .EQ. 0) THEN
      DO l=0,2
        M(l,i,j,k) = 1.0D0
      END DO
      fmat(i,j,k) = .TRUE.
      WRITE(*,*) "hit d000"
      RETURN 
    END IF
    
    WRITE(*,*) "non-base case" 

    ! 2) recursion algorithm 
    ! if i is zero, we will need this in overlap
    IF (j .EQ. k) THEN 
      WRITE(*,*) "going l and r"
      CALL lrec(M,i-1,j-1,k,PA,PB,pp,fmat)
      CALL lrec(M,i,j-1,k,PA,PB,pp,fmat)
      CALL lrec(M,i+1,j-1,k,PA,PB,pp,fmat)
      DO l=0,2
        M(l,i,j,k) = M(l,i-1,j-1,k)/(2*pp) + PA(l)*M(l,i,j-1,k) + (i+1)*M(l,i+1,j-1,k)   
        fmat(i,j,k) = .TRUE. 
      END DO
      WRITE(*,*) "rvec-here0", i,j,k
      CALL rrec(M,i-1,j,k-1,PA,PB,pp,fmat)
      CALL rrec(M,i,j,k-1,PA,PB,pp,fmat)
      CALL rrec(M,i+1,j,k-1,PA,PB,pp,fmat)
      RETURN
    ! get the rest of what we need
    ! if k is bigger, lower it
    ELSE IF (j .LT. k) THEN
      WRITE(*,*) "going r"
      CALL rrec(M,i-1,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,i,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,i+1,j,k-1,PA,PB,pp,fmat) 
      DO l=0,2
        M(l,i,j,k) = M(l,i-1,j-1,k)/(2*pp) + PB(l)*M(l,i,j-1,k) + (i+1)*M(l,i+1,j-1,k)   
        fmat(i,j,k) = .TRUE. 
      END DO
      WRITE(*,*) "rvec-here1"
      RETURN
    ! otherwise, lower j (doesn't really matter which way you do it)
    ELSE IF (j .GT. k) THEN
      WRITE(*,*) "going l"
      CALL lrec(M,i-1,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,i,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,i+1,j-1,k,PA,PB,pp,fmat) 
      DO l=0,2
        M(l,i,j,k) = M(l,i-1,j,k-1)/(2*pp) + PA(l)*M(l,i,j,k-1) + (i+1)*M(l,i+1,j,k-1)   
        fmat(i,j,k) = .TRUE. 
      END DO
      WRITE(*,*) "rvec-here2"
      RETURN
    ELSE 
      WRITE(*,*) "Somehow you broke this at:", i,j,k
      STOP
    END IF

  END SUBROUTINE rrec

!~~~~~ 
  !function to generate coefficients of GTO basis sets 
  REAL(KIND=8) FUNCTION gtoD(l,a)
    IMPLICIT NONE
    
    ! l		: int, angular momentum quantum number
    ! a		: dp, coefficient

    ! Inout 
    REAL(KIND=8), INTENT(IN) :: a
    INTEGER, INTENT(IN) :: l

    IF (l .EQ. 0) THEN
      gtoD = (2*a/Pi)**(3.0D0/4.0D0)
    ELSE IF (l .EQ. 1) THEN
      gtoD = (128.0D0*(a**5.0D0)/(Pi**3.0D0))**(1.0D0/4.0D0)
    ELSE IF (l .EQ. 2) THEN
      gtoD = (2048.0D0*(a**7.0D0)/(9.0D0*Pi**(3.0D0)))**(6.0D0/4.0D0)
    ELSE IF (l .GT. 2) THEN
      WRITE(*,*) "this angular momentum not implimented yet"
      STOP
    END IF

  END FUNCTION gtoD

END PROGRAM int1e
