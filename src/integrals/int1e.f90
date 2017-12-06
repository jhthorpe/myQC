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
  USE aux
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
  ! F		: 2D dp, Fock matrix
  ! MOc		: 2D dp, molecular coefficients (u or v, ith MO) 

  ! Variables  
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: bas
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz,S,F,MOc
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: basinfo
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms,options 
  REAL(KIND=8) :: timeS, timeF, fmem
  INTEGER :: nnuc,nelc,i,j,k,norb,npri,stat
  LOGICAL :: flag1,flag2,flag

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

! 1) calculate overlap matrix and kinetic energy integrals
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
  WRITE(*,*) "Allocating space for Fock matrix (MB) : ", norb*norb*8/1.0E6
  ALLOCATE(F(0:norb-1,0:norb-1),STAT=stat)
  IF (stat .NE. 0) STOP "int1e: max memory reached, exiting"
  fmem = fmem - norb*norb*8/1.0E6
  WRITE(*,*) "Allocating space for MO coefficients (MB) : ", norb*norb*8/1.0E6
  ALLOCATE(MOc(0:norb-1,0:norb-1),STAT=stat)
  IF (stat .NE. 0) STOP "int1e: max memory reached, exiting"
  fmem = fmem - norb*norb*8/1.0E6
  WRITE(*,*)
  CALL nmem(fmem)

  !Temporary initial MO
  ! WORK NOTE : should be initialized with noncoulombic Hamiltonian
  WRITE(*,*) "temporarily taking all MO coefficients as one"
  DO i=0,norb-1
    MOc(:,i) = (/ (1.0D0, j=0,norb-1) /)
  END DO

  INQUIRE(file='Suv',EXIST=flag1)
  INQUIRE(file='Fuv',EXIST=flag2)
  flag = flag1 .AND. flag2
  WRITE(*,*)
  IF (.NOT. flag) THEN
    !1) If not there, calculated overlap and kinetic energy integrals
    CALL SaT(S,F,bas,basinfo,atoms,options,fmem,nnuc,xyz,norb)
    WRITE(*,*) "Overlap written to Suv"
    WRITE(*,*) "Kinetic energy written to Fuv"
    WRITE(*,*)
    !2) calculate 1 electron coulomb integrals
    CALL coulomb(F,bas,basinfo,atoms,options,fmem,nnuc,xyz,norb)
    CALL normS(S,MOc,norb,0)
  ELSE
    CALL EXECUTE_COMMAND_LINE('touch Sold')
    CALL EXECUTE_COMMAND_LINE('touch Fold')    
    OPEN(unit=1,file='Suv',status='old',access='sequential')
    OPEN(unit=2,file='Fuv',status='old',access='sequential')
    WRITE(*,*) "Reading overlap matrix from file"
    READ(1,*) S(:,:)
    WRITE(*,*) "Reading Fock matrix from file"
    READ(2,*) F(:,:)
    CLOSE(unit=2)
    CLOSE(unit=1)
  END IF

  DEALLOCATE(F)
  DEALLOCATE(S)

! output
  CALL setenv(atoms,xyz,fmem,options)
  CALL CPU_TIME(timeF)
  WRITE(*,*) "int1e ran in (s) : ", (timeF - timeS) 

  CONTAINS 

!~~~~~
  ! calculate overlap matrix and kinetic energy terms of fock matrix
  SUBROUTINE SaT(S,F,bas,basinfo,atoms,options,fmem,nnuc,xyz,norb)
    IMPLICIT NONE
    ! Values
    ! S		: 2D dp, overlap matrix
    ! F		: 2D dp, Fock matrix
    ! xyz	: 2D dp, array of nuclear positions
    ! atoms	: 1D int, array of which atom is which
    ! fmem	: dp, free memory left in MB
    ! nnuc	: int, number of nuclii
    ! norb	: int, number of orbitals in molecule
    ! bas	: 2D dp, basis for each atom: atom : orbital : [d,a]
    ! basinfo	: 2D int, array of basis information
    ! options	: 1D int, array of options

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: S,F 
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
        !WRITE(*,*) "u,v, col, row", col0, col, row0, row
        CALL block(S(col0:col0+col-1,row0:row0+row-1),F(col0:col0+col-1,row0:row0+row-1),bas,basinfo,atoms,xyz,u,v,col,row)
        row0 = row0 + row
      END DO
      col0 = col0 + col
    END DO 

    OPEN(unit=1,file='Suv',status='replace',access='sequential')
      WRITE(1,*) S(:,:)    
    CLOSE(unit=1)
 
    OPEN(unit=1,file='Fuv',status='replace',access='sequential')
      WRITE(1,*) F(:,:)
    CLOSE(unit=1)

    CALL CPU_TIME(timeF)
    WRITE(*,*) "overlap matrix constructed in (s) :", (timeF-timeS)

  END SUBROUTINE SaT

!~~~~~
  ! calculate block of overlap matrix 
  SUBROUTINE block(Sb,Fb,bas,basinfo,atoms,xyz,u,v,na,nb)
    IMPLICIT NONE
    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931

    ! Values
    ! Sb	: 2D dp, subblock of overlap matrix
    ! Fb	: 2D dp, subblock of Fock matrix
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
    ! na,nb	: 1D int, angular qunatum number of orbital in atom A,B (x,y,z)

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Sb, Fb 
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms
    INTEGER, INTENT(IN) :: u,v,na,nb

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: coef
    REAL(KIND=8), DIMENSION(0:2) :: PA, PB, AB, PP
    INTEGER, DIMENSION(0:2) :: la,lb,amax, bmax
    REAL(KIND=8) :: EIJ, valSb, valFb, p, m, aa, bb, tempSb, tempFb
    INTEGER :: i,j,k,s,t,a,b,ori
    
    WRITE(*,*) "u,v",u,v
    
    DO a=0,na-1 !iterate through orbitals of atom A
      DO b=0,nb-1 !iterate through orbitals of atom B
       ! WRITE(*,*) "a,b",a,b
        valSb = 0.0D0
        valFb = 0.0D0
        DO s=0,basinfo(u,4*(a+1)+3)-1 !iterate through primatives of orbital a 
          DO t=0,basinfo(v,4*(b+1)+3)-1 !iterate thorugh primatives of orbital b

       !     WRITE(*,*) "s,t", s,t 
            ! 1) get overlap location
            aa = bas(u,a,s*2+1)! WORK NOTE hardcoded to 2nd value in bas of correct u and a 
            bb = bas(v,b,t*2+1)
            !WRITE(*,*) "aa,bb", aa, bb
            p = aa + bb 
            m = aa * bb
            DO i=0,2
              AB(i) = xyz(u,i) - xyz(v,i)
              PP(i) = (aa*xyz(u,i) + bb*xyz(v,i))/p
              PA(i) = PP(i) - xyz(u,i)
              PB(i) = PP(i) - xyz(v,i)
            END DO 
            EIJ = EXP(-m*(AB(0)**2.0D0 + AB(1)**2.0D0 + AB(2)**2.0D0)/p) 


! WORK NOTE - at call getcoef 
            ! 2) get coefficients for these primatives
            !WRITE(*,*) "Sending PA, PB, aa, ab", PA(2),PB(2), aa,bb
            
            !get orientation for A 
            ori = basinfo(u,4*(a+1)+2)
       !     WRITE(*,*) basinfo(u,4*(a+1)+1)
            IF (ori .EQ. -1) THEN  ! s type orbital
              la = [basinfo(u,4*(a+1)+1),basinfo(u,4*(a+1)+1),basinfo(u,4*(a+1)+1)]
              amax = la 
            ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN !p type orbital 
              la = [0, 0, 0]
              la(ori) = basinfo(u,4*(a+1)+1)
              amax = la
            END IF

            ! orientation for B
            ori = basinfo(v,4*(b+1)+2)
            IF (ori .EQ. -1) THEN  ! s type orbital
              lb = [basinfo(v,4*(b+1)+1),basinfo(v,4*(b+1)+1),basinfo(v,4*(b+1)+1)]
              bmax = [lb(0)+2, lb(1)+2, lb(2)+2] !+2 due to kinetic term 
            ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN !p type orbital 
              lb = [0,0,0]
              lb(ori) = basinfo(v,4*(b+1)+1)
              bmax = [lb(0)+2, lb(1)+2, lb(2)+2]
            END IF

            ! getcoef(M,PA,PB,aa,bb,amax,bmax)
            CALL getcoef(coef,PA,PB,aa,bb,amax,bmax)

            ! 3) add everything together for overlap/kinetic matrix elements
            
            tempSb = overlap(u,v,a,b,s,t,p,bas,basinfo,coef,lb,la,aa,bb,EIJ) 
            tempFb = kinetic(u,v,a,b,s,t,p,bas,basinfo,coef,lb,la,aa,bb,EIJ)
           
            valSb = valSb + tempSb
            valFb = valFb + tempFb
            DEALLOCATE(coef)

          END DO !end t
        END DO !end s
        Sb(a,b) = valSb
        Fb(a,b) = valFb
        !WRITE(*,*) "~~~~~~~~~~"
      END DO !end b
    END DO !end a

   !WRITE(*,*) "~~~~~~~~~~~~~~"
   !WRITE(*,*) "~~~~~~~~~~~~~~"
  END SUBROUTINE block

!~~~~~
  !Subroutine to calculate coefficient matrix 

  ! WORK NOTE - Rewriting the coefficients section
  ! replace nmax input with x+y+z max
  ! need to treat x,y,z seperately
  ! WORK NOTE : aa - bb , assuming there are all the same for the orbitals, may break with d-type 

  SUBROUTINE getcoef(M,PA,PB,aa,bb,amax,bmax)
    IMPLICIT NONE

    ! M		: 2D dp, coefficient matrix d_{i}^{j,k} = M(0,i,j,k)
    ! PA	: 1D dp, (xa,ya,za) distance from central gaussian
    ! PB	: 1D dp, (xb,yb,zb) distance from central gaussian
    ! aa	: dp, coefficient of atom A -WARNING, currently assuming are same for all orbs
    ! bb	: dp, coefficient of atom B
    ! amax	: 1D int, max angular quantum number of A needed
    ! amax	: 1D int, max angular quantum number of B needed
    ! nmax	: int, max total quantum number needed
    ! nn	: int, total number of elements
    ! fmat	: 3D dp, boolean array, marks what has been calculated 

    ! Inout
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:), INTENT(INOUT) :: M
    REAL(KIND=8), DIMENSION(0:2), INTENT(IN) :: PA, PB   
    INTEGER, DIMENSION(0:2), INTENT(IN) :: amax, bmax
    REAL(KIND=8), INTENT(IN) :: aa,bb
    
    ! Internal
    LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: fmat
    INTEGER, DIMENSION(0:2) :: nmax,nn
    REAL(KIND=8) :: pp
    INTEGER :: l,i,j,k

    ! 1) construct initial stuff
    pp = aa+bb
    nmax = [amax(0)+bmax(0),amax(1)+bmax(1),amax(2)+bmax(2)] 
   ! WRITE(*,*) "Coefficients are...", aa, bb
   ! WRITE(*,*) "x : nmax, amax, bmax", nmax(0), amax(0), bmax(0)
   ! WRITE(*,*) "PAx, PBz", PA(0), PB(0)
   ! WRITE(*,*) "y : nmax, amax, bmax", nmax(1), amax(1), bmax(1)
   ! WRITE(*,*) "PAy, PBy", PA(1), PB(1)
   ! WRITE(*,*) "z : nmax, amax, bmax", nmax(2), amax(2), bmax(2)
   ! WRITE(*,*) "PAz, PBz", PA(2), PB(2) 
    ALLOCATE(M(0:2,-2:MAXVAL(amax)+MAXVAL(bmax),-2:MAXVAL(amax)+2,-2:MAXVAL(bmax)+2)) 
    ALLOCATE(fmat(0:2,-2:MAXVAL(amax)+MAXVAL(bmax),-2:MAXVAL(amax)+2,-2:MAXVAL(bmax)+2))
    !initialize fmat
    DO l=0,2
      DO i=-2,MAXVAL(amax)+MAXVAL(bmax)
        DO j=-2,MAXVAL(amax)+2
          DO k=-2,MAXVAL(bmax)+2
            M(l,i,j,k) = 0.0D0
            fmat(l,i,j,k) = .FALSE.
          END DO
        END DO
      END DO
    END DO 

    ! 2) call recursive function on top
    ! call the recursion for overlap matrix
    DO l=0,2
    !  WRITE(*,*) "~~~~~~~~~~~working on l = ", l
      CALL lrec(M,l,0,amax(l),bmax(l),PA,PB,pp,fmat)
      CALL rrec(M,l,0,amax(l),bmax(l),PA,PB,pp,fmat)
    END DO
        
    !WRITE(*,*) "===================="
    !WRITE(*,*) "z N = 0"
    !WRITE(*,*) "0,0,0", M(2,0,0,0) 
    !WRITE(*,*) "Z N = 1"
    !WRITE(*,*) "0,0,1", M(2,0,0,1) 
    !WRITE(*,*) "0,1,0", M(2,0,1,0) 
    !WRITE(*,*) "Z N = 2"
    !WRITE(*,*) "0,1,1", M(2,0,1,1)
    !WRITE(*,*) "0,0,2", M(2,0,0,2)
    !WRITE(*,*) "0,2,0", M(2,0,2,0)
    !WRITE(*,*) "Z N = 3"
    !WRITE(*,*) "0,1,2", M(2,0,1,2)
    !WRITE(*,*) "0,2,1", M(2,0,2,1)
    !WRITE(*,*) "0,0,3", M(2,0,0,3)
    !WRITE(*,*) "0,3,0", M(2,0,3,0)
    !WRITE(*,*) "Z N = 4"
    !WRITE(*,*) "0,2,2", M(2,0,2,2)
    !WRITE(*,*) "0,1,3", M(2,0,1,3)
    !WRITE(*,*) "0,3,1", M(2,0,3,1)
    !WRITE(*,*) "===================="

    DEALLOCATE (fmat)

  END SUBROUTINE getcoef
!~~~~~
  !left side recursion on coefficients of overlap matrix
  RECURSIVE SUBROUTINE lrec(M,l,i,j,k,PA,PB,pp,fmat)
    IMPLICIT NONE
   
    ! M		: 4D dp, matrix of coefficients
    ! i,j,k	: int, index of coefficient we need. i=N,j=n,k=nbar
    ! l		: int, coordinite of coefficient we need 0=x,1=y,2=z)
    ! PA, PB	: 1D dp, list of x,y,z distances of atoms A,B from P
    ! fmat	: 3D dp, list of which values we already have
    ! pp	: dp, value of sum of coefficients

    ! INOUT
    REAL(KIND=8), DIMENSION(0:,-2:,-2:,-2:), INTENT(INOUT) :: M
    LOGICAL, DIMENSION(0:,-2:,-2:,-2:),INTENT(INOUT) :: fmat
    REAL(KIND=8), DIMENSION(0:2), INTENT(IN) :: PA, PB
    REAL(KIND=8), INTENT(IN) :: pp
    INTEGER, INTENT(IN) :: i,j,k,l

    !WRITE(*,*) "lrec called with (l,i,j,k)", l,i,j,k
    
    ! 1) base cases
    ! if we've already seen this point
    ! if not, check that it isn't automatically zero
    IF (fmat(l,i,j,k)) THEN
    !  WRITE(*,*) "i,j,k  already seen"
      RETURN
    END IF
    IF (i .LT. 0 .OR. j .LT. 0 .OR. k .LT. 0 .OR. i .GT. j+k) THEN
      M(l,i,j,k) = 0.0D0
      fmat(l,i,j,k) = .TRUE.
    !  WRITE(*,*) "i,j,k is 0"
      RETURN
    ! check for d000 base case
    ELSE IF (i .EQ. 0 .AND. j .EQ. 0 .AND. k .EQ. 0) THEN
      M(l,i,j,k) = 1.0D0
      fmat(l,i,j,k) = .TRUE.
     ! WRITE(*,*) "hit d000"
      RETURN 
    END IF

    !WRITE(*,*) "non-base case" 
    
    ! 2) recursion algorithm 
    ! if i is zero, we will need this in overlap
    IF (j .EQ. k) THEN 
    !  WRITE(*,*) "going l and r"
      CALL lrec(M,l,i-1,j-1,k,PA,PB,pp,fmat)
      CALL lrec(M,l,i,j-1,k,PA,PB,pp,fmat)
      CALL lrec(M,l,i+1,j-1,k,PA,PB,pp,fmat)
    !  WRITE(*,*) "checking lrec=rrec, below", i,j,k, M(l,i-1,j-1,k)/(2.0D0*pp) + PA(l)*M(l,i,j-1,k) + (i+1)*M(l,i+1,j-1,k)
      CALL rrec(M,l,i-1,j,k-1,PA,PB,pp,fmat)
      CALL rrec(M,l,i,j,k-1,PA,PB,pp,fmat)
      CALL rrec(M,l,i+1,j,k-1,PA,PB,pp,fmat)
      M(l,i,j,k) = M(l,i-1,j,k-1)/(2.0D0*pp) + PB(l)*M(l,i,j,k-1) + (i+1)*M(l,i+1,j,k-1)   
      fmat(l,i,j,k) = .TRUE. 
    !  WRITE(*,*) "lvec-here0", i,j,k
    !  WRITE(*,*) "i,j,k val + ", i,j,k,M(l,i,j,k)
      RETURN
    ! get the rest of what we need
    ! if k is bigger, lower it
    ELSE IF (j .LT. k) THEN
    !  WRITE(*,*) "going r"
      CALL rrec(M,l,i-1,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i+1,j,k-1,PA,PB,pp,fmat) 
      M(l,i,j,k) = M(l,i-1,j,k-1)/(2.0D0*pp) + PB(l)*M(l,i,j,k-1) + (i+1)*M(l,i+1,j,k-1)   
      fmat(l,i,j,k) = .TRUE. 
    !  WRITE(*,*) "lvec-here1"
    !  WRITE(*,*) "i,j,k val + ", i,j,k,M(l,i,j,k)
      RETURN
    ! otherwise, lower j (doesn't really matter which way you do it)
    ELSE IF (j .GE. k) THEN
    !  WRITE(*,*) "going l"
      CALL lrec(M,l,i-1,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i+1,j-1,k,PA,PB,pp,fmat) 
      M(l,i,j,k) = M(l,i-1,j-1,k)/(2.0D0*pp) + PA(l)*M(l,i,j-1,k) + (i+1)*M(l,i+1,j-1,k)   
      fmat(l,i,j,k) = .TRUE. 
    !  WRITE(*,*) "lvec-here2"
    !  WRITE(*,*) "i,j,k val + ", i,j,k,M(l,i,j,k)
      RETURN
    ELSE 
      WRITE(*,*) "Somehow you broke this at:", i,j,k
      STOP "error in lrec"
    END IF
 
  END SUBROUTINE lrec
!~~~~~
  !right side (B) recursion on coefficients of overlap matrix
  RECURSIVE SUBROUTINE rrec(M,l,i,j,k,PA,PB,pp,fmat)
    IMPLICIT NONE
   
    ! M		: 4D dp, matrix of coefficients
    ! i,j,k	: int, index of coefficient we need. i=N,j=n,k=nbar
    ! l		: int, coordinate we need (0=x,1=y,2=z)
    ! PA, PB	: 1D dp, list of x,y,z distances of atoms A,B from P
    ! fmat	: 3D dp, list of which values we already have
    ! pp		: dp, value of sum of coefficients

    ! INOUT
    REAL(KIND=8), DIMENSION(0:,-2:,-2:,-2:), INTENT(INOUT) :: M
    LOGICAL, DIMENSION(0:,-2:,-2:,-2:), INTENT(INOUT) :: fmat
    REAL(KIND=8), DIMENSION(0:2), INTENT(IN) :: PA, PB
    REAL(KIND=8), INTENT(IN) :: pp
    INTEGER, INTENT(IN) :: i,j,k,l

    !WRITE(*,*) "rrec called with (l,i,j,k)", l,i,j,k
    
    ! 1) base cases
    ! if we've already seen this point
    ! if not, check that it isn't automatically zero
    IF (fmat(l,i,j,k)) THEN
    !  WRITE(*,*) "i,j,k  already seen"
      RETURN
    END IF
    IF (i .LT. 0 .OR. j .LT. 0 .OR. k .LT. 0 .OR. i .GT. j+k) THEN
      M(l,i,j,k) = 0.0D0
      fmat(l,i,j,k) = .TRUE.
    !  WRITE(*,*) "i,j,k is 0"
      RETURN
    ! check for d000 base case
    ELSE IF (i .EQ. 0 .AND. j .EQ. 0 .AND. k .EQ. 0) THEN
      M(l,i,j,k) = 1.0D0
      fmat(l,i,j,k) = .TRUE.
    !  WRITE(*,*) "hit d000"
      RETURN 
    END IF
    
    !WRITE(*,*) "non-base case" 

    ! 2) recursion algorithm 
    ! if i is zero, we will need this in overlap
    IF (j .EQ. k) THEN 
    !  WRITE(*,*) "going l and r"
      CALL lrec(M,l,i-1,j-1,k,PA,PB,pp,fmat)
      CALL lrec(M,l,i,j-1,k,PA,PB,pp,fmat)
      CALL lrec(M,l,i+1,j-1,k,PA,PB,pp,fmat)
      M(l,i,j,k) = M(l,i-1,j-1,k)/(2.0D0*pp) + PA(l)*M(l,i,j-1,k) + (i+1)*M(l,i+1,j-1,k)   
    !  WRITE(*,*) "rvec-here0", i,j,k
    !  WRITE(*,*) "i,j,k val + ", i,j,k,M(l,i,j,k)
      CALL rrec(M,l,i-1,j,k-1,PA,PB,pp,fmat)
      CALL rrec(M,l,i,j,k-1,PA,PB,pp,fmat)
      CALL rrec(M,l,i+1,j,k-1,PA,PB,pp,fmat)
      fmat(l,i,j,k) = .TRUE. 
    !  WRITE(*,*) "checking (i,j,k) above",  M(l,i-1,j,k-1)/(2.0D0*pp) + PB(l)*M(l,i,j,k-1) + (i+1)*M(l,i+1,j,k-1)
      RETURN
    ! get the rest of what we need
    ! if k is bigger, lower it
    ELSE IF (j .LT. k) THEN
    !  WRITE(*,*) "going r"
      CALL rrec(M,l,i-1,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i+1,j,k-1,PA,PB,pp,fmat) 
      M(l,i,j,k) = M(l,i-1,j,k-1)/(2.0D0*pp) + PB(l)*M(l,i,j,k-1) + (i+1)*M(l,i+1,j,k-1)   
      fmat(l,i,j,k) = .TRUE. 
    !  WRITE(*,*) "rvec-here1"
    !  WRITE(*,*) "i,j,k val + ", i,j,k,M(l,i,j,k)
      RETURN
    ! otherwise, lower j (doesn't really matter which way you do it)
    ELSE IF (j .GT. k) THEN
    !  WRITE(*,*) "going l"
      CALL lrec(M,l,i-1,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i+1,j-1,k,PA,PB,pp,fmat) 
      M(l,i,j,k) = M(l,i-1,j-1,k)/(2.0D0*pp) + PA(l)*M(l,i,j-1,k) + (i+1)*M(l,i+1,j-1,k)   
      fmat(l,i,j,k) = .TRUE. 
    !  WRITE(*,*) "rvec-here2"
    !  WRITE(*,*) "i,j,k val + ", i,j,k,M(l,i,j,k)
      RETURN
    ELSE 
      WRITE(*,*) "Somehow you broke this at:", i,j,k
      STOP "error in rrec"
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
      gtoD = (2.0D0*a/Pi)**(3.0D0/4.0D0)
    ELSE IF (l .EQ. 1) THEN
      gtoD = (128.0D0*(a**5.0D0)/(Pi**3.0D0))**(1.0D0/4.0D0)
    ELSE IF (l .EQ. 2) THEN
      gtoD = (2048.0D0*(a**7.0D0)/(9.0D0*Pi**(3.0D0)))**(6.0D0/4.0D0)
    ELSE IF (l .GT. 2) THEN
      WRITE(*,*) "this angular momentum not implimented yet"
      STOP
    END IF

  END FUNCTION gtoD

!~~~~~
  ! Subroutine to normalize overlap matrix for MO #i 
  SUBROUTINE normS(s,c,norb,i)
    IMPLICIT NONE
    ! S		: 2D dp, overlap matrix
    ! c		: 2D dp, molecular coefficients 
    ! norb	: int, number of molecular orbitals
    ! i		: int, the MO that we're looking at

    !inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: S,c
    INTEGER, INTENT(IN) :: norb,i

    !internal
    REAL(KIND=8) :: temp
    INTEGER :: u,v
     
    temp = 0.0D0
 
    ! get sum of elements
    DO u=0,norb-1
      DO v=0,norb-1
        temp = temp + c(u,i)*S(u,v)*c(v,i)
      END DO
    END DO

    !normalize 
    DO u=0,norb-1
      DO v=0,norb-1
        S(u,v) = S(u,v) / temp
      END DO 
    END DO 

  END SUBROUTINE

!~~~~~
  !overlap to calculate matrix element of overlap matrix
  REAL(KIND=8) FUNCTION overlap(u,v,a,b,s,t,p,bas,basinfo,coef,nb,na,aa,bb,EIJ)
    IMPLICIT NONE

    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931

    REAL(KIND=8), DIMENSION(0:,-2:,-2:,-2:), INTENT(IN) :: coef
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: na, nb
    REAL(KIND=8), INTENT(IN) :: aa, bb, EIJ, p
    INTEGER, INTENT(IN) :: u, v, s, t, a, b

    REAL(KIND=8) :: temp

    ! precalculated constants
    temp = EIJ*(Pi/p)**(3.0D0/2.0D0)*bas(u,a,s*2)*bas(v,b,t*2) ! WORK NOTE - hardcoded
    temp = temp*gtoD(basinfo(u,4*(a+1)+1),aa) !basis set coefficients
    temp = temp*gtoD(basinfo(v,4*(b+1)+1),bb) !basis set coefficeints
    ! integral coefficients
    temp = temp*coef(0,0,na(0),nb(0))*coef(1,0,na(1),nb(1))*coef(2,0,na(2),nb(2))

    overlap = temp

  END FUNCTION overlap

!~~~~~
  !function that calculates kinetic energy of a pair of orbitals
  REAL(KIND=8) FUNCTION kinetic(u,v,a,b,s,t,p,bas,basinfo,coef,nb,na,aa,bb,EIJ)
    IMPLICIT NONE

    ! inout
    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931
    REAL(KIND=8), DIMENSION(0:,-2:,-2:,-2:), INTENT(IN) :: coef
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: na, nb
    REAL(KIND=8), INTENT(IN) :: aa, bb, EIJ, p
    INTEGER, INTENT(IN) :: u, v, s, t, a, b
    
    !internal
    REAL(KIND=8) :: temp,val

    temp = 0.0D0
    val = 0.0D0

    !xpart
    temp = nb(0)*(nb(0)-1)*coef(0,0,na(0),nb(0)-2)
    temp = temp - 2.0D0*bb*nb(0)*coef(0,0,na(0),nb(0)) 
    temp = temp - 2.0D0*bb*(nb(0)+1)*coef(0,0,na(0),nb(0))
    temp = temp + 4.0D0*bb**2.0D0*coef(0,0,na(0),nb(0)+2)
    temp = temp * coef(1,0,na(1),nb(1))*coef(2,0,na(2),nb(2)) 
    val = val + temp
    !ypart
    temp = nb(1)*(nb(1)-1)*coef(1,0,na(1),nb(1)-2)
    temp = temp - 2.0D0*bb*nb(1)*coef(1,0,na(1),nb(1)) 
    temp = temp - 2.0D0*bb*(nb(1)+1)*coef(1,0,na(1),nb(1))
    temp = temp + 4.0D0*bb**2.0D0*coef(1,0,na(1),nb(1)+2)
    temp = temp * coef(0,0,na(0),nb(0))*coef(2,0,na(2),nb(2)) 
    val = val + temp
    !zpart
    temp = nb(2)*(nb(2)-1)*coef(2,0,na(2),nb(2)-2)
    temp = temp - 2.0D0*bb*nb(2)*coef(2,0,na(2),nb(2)) 
    temp = temp - 2.0D0*bb*(nb(2)+1)*coef(2,0,na(2),nb(2))
    temp = temp + 4.0D0*bb**2.0D0*coef(2,0,na(2),nb(2)+2)
    temp = temp * coef(0,0,na(0),nb(0))*coef(1,0,na(1),nb(1)) 
    val = val + temp
 !   WRITE(*,*) "zcoef...", coef(2,0,na(2),nb(2)-2), coef(2,0,na(2),nb(2)), coef(2,0,na(2),nb(2)+2)
 !   WRITE(*,*) "x,ycoef...", coef(0,0,na(0),nb(0)), coef(1,0,na(1),nb(1)) 
 !   WRITE(*,*) "zpart is...", temp
    !leading coefficients
    val = val * (-0.5D0)*EIJ*(Pi/p)**(3.0D0/2.0D0) !integration constants
    val = val * bas(u,a,s*2)*bas(v,b,t*2)     !basis set weights
    val = val * gtoD(basinfo(u,4*(a+1)+1),aa) !primative constants 
    val = val * gtoD(basinfo(v,4*(b+1)+1),bb) !primative constants 

 !   WRITE(*,*) "kinetic energy = ", val
   
    kinetic = val
  
  END FUNCTION kinetic

!~~~~~
! subroutine to calculate coulombic integrals between orbitals
  SUBROUTINE coulomb(F,bas,basinfo,atoms,options,fmem,nnuc,xyz,norb)
    IMPLICIT NONE
    ! Values
    ! F		: 2D dp, Fock matrix
    ! xyz	: 2D dp, array of nuclear positions
    ! atoms	: 1D int, array of which atom is which
    ! fmem	: dp, free memory left in MB
    ! nnuc	: int, number of nuclii
    ! norb	: int, number of orbitals in molecule
    ! bas	: 2D dp, basis for each atom: atom : orbital : [d,a]
    ! basinfo	: 2D int, array of basis information
    ! options	: 1D int, array of options

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: F 
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: options, atoms
    REAL(KIND=8), INTENT(IN) :: fmem
    INTEGER, INTENT(IN) :: norb,nnuc

    WRITE(*,*) "starting coulombic integral calculations" 
    !CALL Boys(1,1)

  END SUBROUTINE coulomb
!~~~~~
END PROGRAM int1e
