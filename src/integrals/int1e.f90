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

!=====================================================================
!			MAIN 
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
  ! set		: 2D dp, set of exponent coefficients on each nuclei
  ! setinfo	: 2D int, array of set information

  ! Variables  
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: bas
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz,S,F,MOc,set
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: basinfo, setinfo
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
  CALL buildBasis(options(2),atoms,bas,basinfo,set,setinfo)
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

  INQUIRE(file='Suv',EXIST=flag1)
  INQUIRE(file='Fuv',EXIST=flag2)
  flag = flag1 .AND. flag2
  IF (.NOT. flag) THEN
    !1) If not there, calculated Overlap and Fock 
    CALL proc1e(S,F,bas,basinfo,atoms,options,fmem,nnuc,xyz,norb,set,setinfo)
    WRITE(*,*) "Overlap written to Suv"
    WRITE(*,*) "One electron energy written to Fuv"
    WRITE(*,*)
  ELSE
    CALL EXECUTE_COMMAND_LINE('touch Sold')
    CALL EXECUTE_COMMAND_LINE('touch Fold')    
    OPEN(unit=1,file='Suv',status='old',access='sequential')
    OPEN(unit=2,file='Fuv',status='old',access='sequential')
    WRITE(*,*) "Reading overlap matrix from Suv"
    READ(1,*) S(:,:)
    WRITE(*,*) "Reading Fock matrix from Fuv"
    READ(1,*) F(:,:)
    CLOSE(unit=2)
    CLOSE(unit=1)
  END IF

  ! 2) molecular orbital coefficients
  ! WORK NOTE : should be initialized with 1 electron Hamiltonian?
!    CALL normS(S,MOc,norb,0)
  INQUIRE(file='Cuv',EXIST=flag2)
  IF (flag2) THEN
    WRITE(*,*) "Reading MO coefficients from Cuv"
  ELSE
    WRITE(*,*) "Calculating initial Cuv guess"
    WRITE(*,*) "temporarily taking all MO coefficients as one"
    DO i=0,norb-1
      MOc(:,i) = (/ (1.0D0, j=0,norb-1) /)
    END DO
    OPEN(unit=1,file='Cuv',status='replace',access='sequential')
    WRITE(1,*) MOc(:,:)
    CLOSE(unit=1)
    WRITE(*,*) "MO coefficients written to Cuv"
  END IF
  WRITE(*,*)

  DEALLOCATE(F)
  DEALLOCATE(S)
  DEALLOCATE(MOc)
  DEALLOCATE(bas)
  DEALLOCATE(set)
  DEALLOCATE(basinfo)
  DEALLOCATE(setinfo)

! output
  fmem = fmem + 3*norb*norb*8/1.0E6
  CALL setenv(atoms,xyz,fmem,options)
  CALL CPU_TIME(timeF)
  WRITE(*,*) "int1e ran in (s) : ", (timeF - timeS) 

  CONTAINS 
!=====================================================================
!			SUBROUTINES

!---------------------------------------------------------------------
!		Proccess one electron integrals
!---------------------------------------------------------------------
  SUBROUTINE proc1e(S,F,bas,basinfo,atoms,options,fmem,nnuc,xyz,norb,set,setinfo)
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
    ! set	: 1D dp, array of exponential alphas
    ! setinfo	: 2D int, stores information about set

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: S,F 
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz, set
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo, setinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: options, atoms
    REAL(KIND=8), INTENT(IN) :: fmem
    INTEGER, INTENT(IN) :: norb,nnuc

    !Internal
    REAL(KIND=8) :: timeS, timeF
    INTEGER :: u,v,col,row,col0,row0,i,j

    CALL CPU_TIME(timeS)
    WRITE(*,*) "constructing Overlap and Fock matrix"

    !Matrix blocked by atoms
    col0 = 0
    DO u=0,nnuc-1
      row0 = 0
      DO v=0,nnuc-1 
        !determine length, width, and starting points of block
        col = basinfo(u,2) !number of orbitals in a
        row = basinfo(v,2) ! number of orbitals in b

        !construct overlap within block
        CALL block(S(col0:col0+col-1,row0:row0+row-1),F(col0:col0+col-1,row0:row0+row-1), &
        bas,basinfo,atoms,xyz,u,v,col,row,set,setinfo)

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
    WRITE(*,*) "Overlap and Fock constructed in (s) :", (timeF-timeS)

  END SUBROUTINE proc1e

!---------------------------------------------------------------------
!		Calculate block of Overlap and Fock matrix 
!---------------------------------------------------------------------
  SUBROUTINE block(Sb,Fb,bas,basinfo,atoms,xyz,u,v,na,nb,set,setinfo)
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
    ! coef	: 3D dp, array of coefficients for overlap, (xyz,i,j,k)
    ! PA,PB,PP	: 1D dp, array of distances to overlap, PP = overlap
    ! EIJ	: dp, gaussian at molecular center
    ! val 	: dp, current evaluation of integral
    ! p		: dp, addition of two coeffiecients
    ! aa,bb	: dp, coefficients of Guassians for atoms a and b
    ! na,nb	: 1D int, angular qunatum number of orbital in atom A,B (x,y,z)
    ! set	: 2D dp, set of exponential coefficients
    ! setinfo	: 2D int, stores information about set 

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Sb, Fb 
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz, set
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo, setinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms
    INTEGER, INTENT(IN) :: u,v,na,nb

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: coef
    REAL(KIND=8), DIMENSION(0:2) :: PA, PB, AB, PP
    INTEGER, DIMENSION(0:2) :: la,lb,amax, bmax
    REAL(KIND=8) :: EIJ, valSb, valFb, p, m, aa, bb, tempSb, tempFb
    INTEGER :: i,j,k,s,t,a,b,ori,nset,setlena,setlenb

    WRITE(*,*) "u,v", u, v
    setlena = setinfo(u,2)
    setlenb = setinfo(v,2)

    !zero Sb and Fb
    DO i=0,na-1
      Sb(i,:) = (/ (0.0D0, j=0,nb-1) /)
      Fb(i,:) = (/ (0.0D0, j=0,nb-1) /)
    END DO


    DO a=0, setinfo(u,0)-1 !iterate through set A 
      aa = set(u,a) !alpha a
      amax = setinfo(u,2+a*setlena+2)   !get max ang qn
      DO b=0, setinfo(v,0)-1 !iterate through set B 
!        valSb = 0.0D0
!        valFb = 0.0D0

        bb = set(v,b) !alpha b
        bmax = setinfo(v,2+b*setlenb+2)

        ! 1) get overlap location
        p = aa + bb 
        m = aa * bb

        DO i=0,2
          AB(i) = xyz(u,i) - xyz(v,i)
          PP(i) = (aa*xyz(u,i) + bb*xyz(v,i))/p
          PA(i) = PP(i) - xyz(u,i)
          PB(i) = PP(i) - xyz(v,i)
        END DO 

        ! screen for sufficiently small constant
        EIJ = EXP(-m*(AB(0)**2.0D0 + AB(1)**2.0D0 + AB(2)**2.0D0)/p) 
        IF (EIJ .LT. 1.0D-14) THEN
          CYCLE
        END IF 
        
        CALL getcoef(coef,PA,PB,aa,bb,amax,bmax+2) ! +2 for kinetic term

        ! 2) use setinfo to update Overlap and Fock
        CALL overlap(Sb,u,v,a,b,p,bas,basinfo,coef,setinfo(u,2+a*setlena+1:2+(a+1)*setlena),&
        setinfo(v,2+b*setlenb+1:2+(b+1)*setlenb),aa,bb,EIJ)
 
        DEALLOCATE(coef)

      END DO !end b
    END DO !end a

!        DO s=0,basinfo(u,4*(a+1)+3)-1 !iterate through primatives of set a 
!          DO t=0,basinfo(v,4*(b+1)+3)-1 !iterate thorugh primatives of set b
!
!            ! 3) add everything together for overlap/kinetic matrix elements
!            
!            tempSb = overlap(u,v,a,b,s,t,p,bas,basinfo,coef,la,lb,aa,bb,EIJ) 
!            tempFb = kinetic(u,v,a,b,s,t,p,bas,basinfo,coef,la,lb,aa,bb,EIJ)
!            tempFb = tempFb + coulomb(u,v,a,b,s,t,p,bas,basinfo,PP,la,lb,aa,bb,atoms,EIJ,coef)
!            tempFb = coulomb(u,v,a,b,s,t,p,bas,basinfo,PP,la,lb,aa,bb,atoms,EIJ,coef)
!   
!            valSb = valSb + tempSb
!            valFb = valFb + tempFb
!            DEALLOCATE(coef)
!
!          END DO !end t
!        END DO !end s
!        Sb(a,b) = valSb
!        Fb(a,b) = valFb

  END SUBROUTINE block

!---------------------------------------------------------------------
!		Calculate coefficients of overlap Gaussians 
!---------------------------------------------------------------------

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
    ! fmat	: 3D dp, boolean array, tracks what has been calculated 

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

    pp = aa+bb
    nmax = [amax(0)+bmax(0),amax(1)+bmax(1),amax(2)+bmax(2)] 

    ALLOCATE(M(0:2,-2:MAXVAL(amax)+MAXVAL(bmax),-2:MAXVAL(amax)+2,-2:MAXVAL(bmax)+2)) 
    ALLOCATE(fmat(0:2,-2:MAXVAL(amax)+MAXVAL(bmax),-2:MAXVAL(amax)+2,-2:MAXVAL(bmax)+2))

    !initialize coefficient matrices
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

    ! call recursive function on top
    DO l=0,2
      CALL lrec(M,l,0,amax(l),bmax(l),PA,PB,pp,fmat)
      CALL rrec(M,l,0,amax(l),bmax(l),PA,PB,pp,fmat)
    END DO
        
    DEALLOCATE (fmat)

  END SUBROUTINE getcoef

!---------------------------------------------------------------------
!	Left side (A) recursion on coefficients Hermite Gaussians
!---------------------------------------------------------------------
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

    ! 1) base cases
    ! check if we've already seen this point 
    IF (fmat(l,i,j,k)) THEN
      RETURN
    END IF
    ! check for "bad" points
    IF (i .LT. 0 .OR. j .LT. 0 .OR. k .LT. 0 .OR. i .GT. j+k) THEN
      M(l,i,j,k) = 0.0D0
      fmat(l,i,j,k) = .TRUE.
      RETURN
    ! check for d000 base case
    ELSE IF (i .EQ. 0 .AND. j .EQ. 0 .AND. k .EQ. 0) THEN
      M(l,i,j,k) = 1.0D0
      fmat(l,i,j,k) = .TRUE.
      RETURN 
    END IF

    ! 2) non-base case
    ! if i is zero, we will need this in overlap
    IF (j .EQ. k) THEN 
      CALL lrec(M,l,i-1,j-1,k,PA,PB,pp,fmat)
      CALL lrec(M,l,i,j-1,k,PA,PB,pp,fmat)
      CALL lrec(M,l,i+1,j-1,k,PA,PB,pp,fmat)
      CALL rrec(M,l,i-1,j,k-1,PA,PB,pp,fmat)
      CALL rrec(M,l,i,j,k-1,PA,PB,pp,fmat)
      CALL rrec(M,l,i+1,j,k-1,PA,PB,pp,fmat)
      M(l,i,j,k) = M(l,i-1,j,k-1)/(2.0D0*pp) + PB(l)*M(l,i,j,k-1) + (i+1)*M(l,i+1,j,k-1)   
      fmat(l,i,j,k) = .TRUE. 
      RETURN
    ! get the rest of what we need
    ! if k is bigger, lower it
    ELSE IF (j .LT. k) THEN
      CALL rrec(M,l,i-1,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i+1,j,k-1,PA,PB,pp,fmat) 
      M(l,i,j,k) = M(l,i-1,j,k-1)/(2.0D0*pp) + PB(l)*M(l,i,j,k-1) + (i+1)*M(l,i+1,j,k-1)   
      fmat(l,i,j,k) = .TRUE. 
      RETURN
    ! otherwise, lower j (doesn't really matter which way you do it)
    ELSE IF (j .GT. k) THEN
      CALL lrec(M,l,i-1,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i+1,j-1,k,PA,PB,pp,fmat) 
      M(l,i,j,k) = M(l,i-1,j-1,k)/(2.0D0*pp) + PA(l)*M(l,i,j-1,k) + (i+1)*M(l,i+1,j-1,k)   
      fmat(l,i,j,k) = .TRUE. 
      RETURN
    ELSE 
      WRITE(*,*) "Somehow you broke this at:", i,j,k
      STOP "error in lrec"
    END IF
 
  END SUBROUTINE lrec

!---------------------------------------------------------------------
!	Right side (B) recursion on coefficients of Hermite Gaussians
!---------------------------------------------------------------------
  RECURSIVE SUBROUTINE rrec(M,l,i,j,k,PA,PB,pp,fmat)
    IMPLICIT NONE
   
    ! M		: 4D dp, matrix of coefficients
    ! i,j,k	: int, index of coefficient we need. i=N,j=n,k=nbar
    ! l		: int, coordinate we need (0=x,1=y,2=z)
    ! PA, PB	: 1D dp, list of x,y,z distances of atoms A,B from P
    ! fmat	: 3D dp, list of which values we already have
    ! pp	: dp, value of sum of coefficients

    ! INOUT
    REAL(KIND=8), DIMENSION(0:,-2:,-2:,-2:), INTENT(INOUT) :: M
    LOGICAL, DIMENSION(0:,-2:,-2:,-2:), INTENT(INOUT) :: fmat
    REAL(KIND=8), DIMENSION(0:2), INTENT(IN) :: PA, PB
    REAL(KIND=8), INTENT(IN) :: pp
    INTEGER, INTENT(IN) :: i,j,k,l

    ! 1) base cases
    ! check if we've seen before 
    IF (fmat(l,i,j,k)) THEN
      RETURN
    END IF
    ! check for "bad" points 
    IF (i .LT. 0 .OR. j .LT. 0 .OR. k .LT. 0 .OR. i .GT. j+k) THEN
      M(l,i,j,k) = 0.0D0
      fmat(l,i,j,k) = .TRUE.
      RETURN
    ! check for d000 base case
    ELSE IF (i .EQ. 0 .AND. j .EQ. 0 .AND. k .EQ. 0) THEN
      M(l,i,j,k) = 1.0D0
      fmat(l,i,j,k) = .TRUE.
      RETURN 
    END IF
    
    ! 2) non-base case 
    ! if i is zero, we will need this in overlap
    IF (j .EQ. k) THEN 
      CALL lrec(M,l,i-1,j-1,k,PA,PB,pp,fmat)
      CALL lrec(M,l,i,j-1,k,PA,PB,pp,fmat)
      CALL lrec(M,l,i+1,j-1,k,PA,PB,pp,fmat)
      M(l,i,j,k) = M(l,i-1,j-1,k)/(2.0D0*pp) + PA(l)*M(l,i,j-1,k) + (i+1)*M(l,i+1,j-1,k)   
      CALL rrec(M,l,i-1,j,k-1,PA,PB,pp,fmat)
      CALL rrec(M,l,i,j,k-1,PA,PB,pp,fmat)
      CALL rrec(M,l,i+1,j,k-1,PA,PB,pp,fmat)
      fmat(l,i,j,k) = .TRUE. 
      RETURN
    ! get the rest of what we need
    ! if k is bigger, lower it
    ELSE IF (j .LT. k) THEN
      CALL rrec(M,l,i-1,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i+1,j,k-1,PA,PB,pp,fmat) 
      M(l,i,j,k) = M(l,i-1,j,k-1)/(2.0D0*pp) + PB(l)*M(l,i,j,k-1) + (i+1)*M(l,i+1,j,k-1)   
      fmat(l,i,j,k) = .TRUE. 
      RETURN
    ! otherwise, lower j (doesn't really matter which way you do it)
    ELSE IF (j .GT. k) THEN
      CALL lrec(M,l,i-1,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i+1,j-1,k,PA,PB,pp,fmat) 
      M(l,i,j,k) = M(l,i-1,j-1,k)/(2.0D0*pp) + PA(l)*M(l,i,j-1,k) + (i+1)*M(l,i+1,j-1,k)   
      fmat(l,i,j,k) = .TRUE. 
      RETURN
    ELSE 
      WRITE(*,*) "Somehow you broke this, i,j,k", i,j,k
      STOP "error in rrec"
    END IF

  END SUBROUTINE rrec

!---------------------------------------------------------------------
!      	Normalize the overlap matrix to get MO coefficients 
!---------------------------------------------------------------------
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

!---------------------------------------------------------------------
!	Calculate matrix element of overlap matrix
!---------------------------------------------------------------------
  SUBROUTINE overlap(Sb,u,v,a,b,p,bas,basinfo,coef,seta,setb,aa,bb,EIJ)
    IMPLICIT NONE

    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931

    ! Values
    ! a,b	: int, set number we're on
    ! la,lb	: 1D int, angular quantum numbers 
    ! ta,tb	: int, tracking a and b
    
    ! Inout
    REAL(KIND=8), DIMENSION(0:,-2:,-2:,-2:), INTENT(IN) :: coef
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Sb
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: seta, setb
    REAL(KIND=8), INTENT(IN) :: aa, bb, EIJ, p
    INTEGER, INTENT(IN) :: u, v, a, b

    ! internal
    INTEGER, DIMENSION(0:2) :: la,lb
    REAL(KIND=8) :: temp
    INTEGER :: i,j,k,orba,orbb,ori,ta,tb,prima,primb

    ta = 0
    !update each element in set
    DO i=0,seta(0)-1 !go through set A
      orba = seta(2+i) !id of orbital
      !find where we are in bas weights
!      ta = 0
!      DO k=0,orba-1
!        ta = ta + basinfo(u,3+k*4+4)
!      END DO

      tb = 0
      DO j=0,setb(0)-1 !go through set B
        orbb = setb(2+j) !id of orbital
        primb = basinfo(v,3+orbb*4+4)

!        tb = 0
!        DO k=0,orbb-1
!          tb = tb + basinfo(v,3+k*4+4)
!        END DO
        WRITE(*,*) "a,b,ta,tb,orba,orbb,a-ta,b-tb",a,b,ta,tb,orba,orbb,a-ta,b-tb 
        
        ! get angular quantum numbers for each orbital 
        ori = basinfo(u,4*(orba+1)+2)
        IF (ori .EQ. -1) THEN  ! s type orbital
          la = [basinfo(u,4*(orba+1)+1),basinfo(u,4*(orba+1)+1),basinfo(u,4*(orba+1)+1)]
        ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN !p type orbital 
          la = [0, 0, 0]
          la(ori) = basinfo(u,4*(orba+1)+1)
        END IF

        ori = basinfo(v,4*(orbb+1)+2)
        IF (ori .EQ. -1) THEN  ! s type orbital
          lb = [basinfo(v,4*(orbb+1)+1),basinfo(v,4*(orbb+1)+1),basinfo(v,4*(orbb+1)+1)]
        ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN !p type orbital 
          lb = [0,0,0]
          lb(ori) = basinfo(v,4*(orbb+1)+1)
        END IF

        ! update Suv 
!        temp = EIJ*(Pi/p)**(3.0D0/2.0D0)*bas(u,orba,(a-ta)*2)*bas(v,orbb,(b-tb)*2) ! WORK NOTE - hardcoded in bas
        temp = 1.0D0
        temp = temp * gtoD(basinfo(u,4*(orba+1)+1),aa)                !basis set coefficients
        temp = temp * gtoD(basinfo(v,4*(orbb+1)+1),bb)                !basis set coefficeints
        temp = temp * coef(0,0,la(0),lb(0))*coef(1,0,la(1),lb(1))*coef(2,0,la(2),lb(2))

        Sb(orba,orbb) = Sb(orba,orbb) + temp
        
        tb = primb !number of primatives in previous set
      END DO 
      prima = basinfo(v,3+orba*4+4)
      ta = prima
    END DO

  END SUBROUTINE overlap


!=====================================================================
!                       FUNCTIONS

!---------------------------------------------------------------------
!	Generate coefficients of GTO basis sets 
!---------------------------------------------------------------------
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
      gtoD = (2048.0D0*(a**7.0D0)/(9.0D0*Pi**(3.0D0)))**(1.0D0/4.0D0)
    ELSE IF (l .GT. 2) THEN
      WRITE(*,*) "this angular momentum not implimented yet"
      STOP
    END IF

  END FUNCTION gtoD

!---------------------------------------------------------------------
!	Calculate kinetic energy of a pair of orbitals
!---------------------------------------------------------------------
  REAL(KIND=8) FUNCTION kinetic(u,v,a,b,s,t,p,bas,basinfo,coef,na,nb,aa,bb,EIJ)
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
    !leading coefficients
    val = val * (-0.5D0)*EIJ*(Pi/p)**(3.0D0/2.0D0) !integration constants
    val = val * bas(u,a,s*2)*bas(v,b,t*2)     !basis set weights
    val = val * gtoD(basinfo(u,4*(a+1)+1),aa) !primative constants 
    val = val * gtoD(basinfo(v,4*(b+1)+1),bb) !primative constants 

    kinetic = val
  
  END FUNCTION kinetic

!---------------------------------------------------------------------
!	Calculate the coulomb potential of a gaussian of two orbitals
!---------------------------------------------------------------------
  REAL(KIND=8) FUNCTION coulomb(u,v,a,b,s,t,ap,bas,basinfo,PP,na,nb,aa,bb,atoms,EIJ,coef)
    IMPLICIT NONE
    ! Values
    ! T		: dp, input to RNLMj
    ! ap	: dp, alpha of overlap gausian
    ! Rtab	: 4D dp, table of RNLM values
    ! Rbol	: 4D dp, table input to RNLMj
    ! nb,na	: 1D int, angular momentum values for b and a, {n,l,m} 
    ! nucpos	: 2D dp, list of nuclear positions
    ! PP	: 1D dp, list of overlap x,y,z
    ! CP 	: 1D dp, line segment between nucleus C and overlap PP
    ! Fj	: 1D dp, table of Boys integral 
    ! EIJ	: dp, correction of combining two gaussians
    ! coef	: 4D dp, table of guassian coefficients from overlap 

    ! inout
    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931
    REAL(KIND=8), DIMENSION(0:,-2:,-2:,-2:), INTENT(IN) :: coef
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: PP
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: na, nb, atoms
    REAL(KIND=8), INTENT(IN) :: aa, bb, ap, EIJ
    INTEGER, INTENT(IN) :: u, v, s, t, a, b
    
    !internal
    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: Rtab
    LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: Rbol
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: nucpos
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Fj,CP
    REAL(KIND=8) :: temp,val,TT,hmm
    INTEGER :: c,p,i,j,k,N,L,M,nnuc,dummy

    val = 0.0D0
    nnuc = SIZE(atoms)

    ! combined quantum numbers
    N = na(0) + nb(0)
    L = na(1) + nb(1)
    M = na(2) + nb(2)

    ALLOCATE(Fj(0:N+L+M))
    ALLOCATE(CP(0:2))
    ALLOCATE(Rtab(-2:N,-2:L,-2:M,0:N+L+M))
    ALLOCATE(Rbol(-2:N,-2:L,-2:M,0:N+L+M))
    ALLOCATE(nucpos(0:nnuc-1,0:2))

    ! get nuclear positions
    OPEN(unit=1,file='nucpos',status='old',access='sequential')
    DO c=0,nnuc-1
      READ(1,*) dummy, nucpos(c,0:2) 
    END DO 
    CLOSE(unit=1)

    ! loop through atoms
    DO c=0,nnuc-1

      !construct PC
      DO i=0,2
        CP(i) = nucpos(c,i) - PP(i)
      END DO
      TT = ap*(CP(0)**2.0D0 + CP(1)**2.0D0 + CP(2)**2.0D0)

      !get Boys table
      DO i=0,N+L+M
        Fj(i) = 0.0D0 
      END DO
      CALL Boys(Fj,N+L+M,TT)

      !setup recursion tables
       DO i=-2,N
         DO j=-2,L
           DO k=-2,M
             DO p=0,N+L+M
               Rtab(i,j,k,p) = 0.0D0
               Rbol(i,j,k,p) = .FALSE.
             END DO
           END DO
         END DO
       END DO
      
!      CALL RNLMj(ABS(CP(0)),ABS(CP(1)),ABS(CP(2)),N,L,M,0,ap,Fj,Rtab,Rbol)

      temp = 0.0D0
 
      ! skip cycle if coef is 0
      ! Sum RNLM over all N,L,M 
      DO i=0,N ! loop over N
        !IF (ABS(coef(0,i,na(0),nb(0))) .LT. 1.0D-14) CYCLE
        DO j=0,L !loop over L
          !IF (ABS(coef(1,j,na(1),nb(1))) .LT. 1.0D-14) CYCLE 
          DO k=0,M !loop over M
            !IF (ABS(coef(2,k,na(2),nb(2))) .LT. 1.0D-14) CYCLE
            CALL RNLMj(CP(0),CP(1),CP(2),i,j,k,0,ap,Fj,Rtab,Rbol)
!            CALL RNLMj(ABS(CP(0)),ABS(CP(1)),ABS(CP(2)),i,j,k,0,ap,Fj,Rtab,Rbol)
!            IF (.NOT. Rbol(i,j,k,0)) THEN
!              WRITE(*,*) "here - bad Rijk", i,j,k
!              STOP
!            END IF 
            temp = temp + coef(0,i,na(0),nb(0))*coef(1,j,na(1),nb(1))*coef(2,k,na(2),nb(2))*Rtab(i,j,k,0)      
          END DO
        END DO
      END DO

      ! add to coulombic integral
      temp = temp * (2.0D0*Pi/ap)                  !from Boys
      temp = temp * bas(u,a,s*2)*bas(v,b,t*2)      !basis set weights
      temp = temp * gtoD(basinfo(u,4*(a+1)+1),aa)  !primative constants
      temp = temp * gtoD(basinfo(v,4*(b+1)+1),bb)  !primative constants
      temp = temp * atoms(c)*EIJ                   !proton number and overlap coefficient
      val = val - temp                             !negative sign for coulombic attraction

    END DO
    
    coulomb = val

    DEALLOCATE(Fj)
    DEALLOCATE(CP)
    DEALLOCATE(Rtab)
    DEALLOCATE(Rbol)
    DEALLOCATE(nucpos)

  END FUNCTION coulomb
!---------------------------------------------------------------------
END PROGRAM int1e
