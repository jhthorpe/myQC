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
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz,S,F,MOc
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: bas,set
  INTEGER, ALLOCATABLE, DIMENSION(:) :: basinfo, setinfo
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
  CALL buildBasis(options(2),atoms,bas,basinfo,set,setinfo,.TRUE.)
  ! write out basis set for checking, include in verbosity later
  ! WORK NOTE - maybe add in memory of basis set here?

  CALL nmem(fmem)

! 1) calculate overlap matrix and kinetic energy integrals
  !get number of orbitals
  npri = 0
  norb = basinfo(1) 
  DO i=0,norb-1
      npri = npri + basinfo(1+i*5+4) 
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
    ! bas	: 1D dp, basis weights for each atom
    ! basinfo	: 1D int, array of basis information
    ! options	: 1D int, array of options
    ! set	: 1D dp, array of exponential alphas
    ! setinfo	: 1D int, stores information about set
    ! coef	: 3D dp, array of coefficients for overlap, (xyz,i,j,k)
    ! PA,PB,PP	: 1D dp, array of distances to overlap, PP = overlap
    ! EIJ	: dp, gaussian at molecular center
    ! val 	: dp, current evaluation of integral
    ! p		: dp, addition of two coeffiecients
    ! aa,bb	: dp, coefficients of Guassians for atoms a and b
    ! OpS	: int, max orbitals per set

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: S,F 
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: bas, set
    INTEGER, DIMENSION(0:), INTENT(IN) :: basinfo, setinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: options, atoms
    REAL(KIND=8), INTENT(IN) :: fmem
    INTEGER, INTENT(IN) :: norb,nnuc

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: coef
    REAL(KIND=8), DIMENSION(0:2) :: PA, PB, AB, PP
    INTEGER, DIMENSION(0:2) :: la,lb,amax,bmax
    REAL(KIND=8) :: EIJ, valSb, valFb, p, m, aa, bb, tempSb, tempFb
    REAL(KIND=8) :: timeS, timeF
    INTEGER :: a,b,u,v,col,row,col0,row0,i,j,k,l,setl,nset,ori,OpS

    CALL CPU_TIME(timeS)
    WRITE(*,*) "constructing Overlap and 1e-Fock matrix"

    !zero S and F
    DO i=0,norb-1
     S(i,:) = (/ (0.0D0, j=0,norb-1) /)
     F(i,:) = (/ (0.0D0, j=0,norb-1) /)
    END DO

    nset = setinfo(0)
    setl = setinfo(1)
    OpS = basinfo(0)

    ! left set
    DO a=0,nset-1
      
      !get data
      aa = set(a)                            !alpha a
      u = setinfo(1+a*setl+3)                !center number
      DO l=0,2                               !set amax values
        amax(l) = setinfo(1+a*setl+2)        !get max ang qn
        la(l) = amax(l)
      END DO
     
      ! right set 
      DO b=0,nset-1 

        !get data
        bb = set(b)                           !alpha b
        v = setinfo(1+b*setl+3)               !center number
        DO l=0,2
          bmax(l) = setinfo(1+b*setl+2) + 2  ! + 2 for kinetic terms
          lb(l) = setinfo(1+b*setl+2)
        END DO

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
  
        CALL getcoef(coef,PA,PB,aa,bb,amax,bmax) 

        ! 2) use setinfo to update Overlap and Fock
        CALL overlap(S,u,v,a,b,p,bas(a*OpS:(a+1)*Ops-1),bas(b*Ops:(b+1)*OpS-1),basinfo,coef,&
        setinfo(1+a*setl+1:1+(a+1)*setl),setinfo(1+b*setl+1:1+(b+1)*setl),aa,bb,EIJ)

        DEALLOCATE(coef)

      END DO                                   !loop over left orbital
    END DO                                     !loop over right orbital

    OPEN(unit=1,file='Suv',status='replace',access='sequential')
      WRITE(1,*) S(:,:)    
    CLOSE(unit=1)
 
    OPEN(unit=1,file='Fuv',status='replace',access='sequential')
      WRITE(1,*) F(:,:)
    CLOSE(unit=1)

    CALL CPU_TIME(timeF)
    WRITE(*,*) "Overlap and Fock constructed in (s) :", (timeF-timeS)

  END SUBROUTINE proc1e
! 
!!        CALL kinetic(Fb,u,v,a,b,p,bas,basinfo,coef,setinfo(u,2+a*setlena+1:2+(a+1)*setlena),&
!!        setinfo(v,2+b*setlenb+1:2+(b+1)*setlenb),aa,bb,EIJ)
!
!!        CALL coulomb(Fb,u,v,a,b,p,bas,basinfo,PP,setinfo(u,2+a*setlena+1:2+(a+1)*setlena),&
!!        setinfo(v,2+b*setlenb+1:2+(b+1)*setlenb),aa,bb,atoms,EIJ,coef,setlena,setlenb,la,lb)
!
!        DEALLOCATE(coef)
!
!      END DO !end b
!
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
    INTEGER :: i,j,k,l,s,t,a,b,ori,nset,setlena,setlenb

    WRITE(*,*) "u,v", u, v
    setlena = setinfo(u,2)
    setlenb = setinfo(v,2)

    !zero Sb and Fb
    DO i=0,na-1
      Sb(i,:) = (/ (0.0D0, j=0,nb-1) /)
      Fb(i,:) = (/ (0.0D0, j=0,nb-1) /)
    END DO

    DO a=0, setinfo(u,0)-1                   !iterate through set A 
      aa = set(u,a)                          !alpha a

      DO l=0,2                               !set amax values
        amax(l) = setinfo(u,2+a*setlena+2)   !get max ang qn
        la(l) = amax(l)
      END DO

      DO b=0, setinfo(v,0)-1                 !iterate through set B 
        bb = set(v,b)                        !alpha b

        DO l=0,2
          bmax(l) = setinfo(v,2+b*setlenb+2) + 2
          lb(l) = setinfo(v,2+b*setlenb+2)
        END DO

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
  
        CALL getcoef(coef,PA,PB,aa,bb,amax,bmax) 

        ! 2) use setinfo to update Overlap and Fock
!        CALL overlap(Sb,u,v,a,b,p,bas,basinfo,coef,setinfo(u,2+a*setlena+1:2+(a+1)*setlena),&
!        setinfo(v,2+b*setlenb+1:2+(b+1)*setlenb),aa,bb,EIJ)
 
!        CALL kinetic(Fb,u,v,a,b,p,bas,basinfo,coef,setinfo(u,2+a*setlena+1:2+(a+1)*setlena),&
!        setinfo(v,2+b*setlenb+1:2+(b+1)*setlenb),aa,bb,EIJ)

!        CALL coulomb(Fb,u,v,a,b,p,bas,basinfo,PP,setinfo(u,2+a*setlena+1:2+(a+1)*setlena),&
!        setinfo(v,2+b*setlenb+1:2+(b+1)*setlenb),aa,bb,atoms,EIJ,coef,setlena,setlenb,la,lb)

        DEALLOCATE(coef)

      END DO !end b
    END DO !end a


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
    IF (j .LE. k) THEN                        ! go right
      CALL lrec(M,l,i-1,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i+1,j-1,k,PA,PB,pp,fmat) 
      CALL rrec(M,l,i-1,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i+1,j,k-1,PA,PB,pp,fmat) 
      M(l,i,j,k) = M(l,i-1,j,k-1)/(2.0D0*pp) + PB(l)*M(l,i,j,k-1) + (i+1)*M(l,i+1,j,k-1)   
      fmat(l,i,j,k) = .TRUE. 
      RETURN
    ELSE IF (j .GT. k) THEN                   ! go left
      CALL lrec(M,l,i-1,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i+1,j-1,k,PA,PB,pp,fmat) 
      M(l,i,j,k) = M(l,i-1,j-1,k)/(2.0D0*pp) + PA(l)*M(l,i,j-1,k) + (i+1)*M(l,i+1,j-1,k)   
      CALL rrec(M,l,i-1,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i+1,j,k-1,PA,PB,pp,fmat) 
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
    IF (j .LE. k) THEN
      CALL lrec(M,l,i-1,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i,j-1,k,PA,PB,pp,fmat) 
      CALL lrec(M,l,i+1,j-1,k,PA,PB,pp,fmat) 
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
      CALL rrec(M,l,i-1,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i,j,k-1,PA,PB,pp,fmat) 
      CALL rrec(M,l,i+1,j,k-1,PA,PB,pp,fmat) 
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
  SUBROUTINE overlap(S,u,v,a,b,p,basa,basb,basinfo,coef,seta,setb,aa,bb,EIJ)
    IMPLICIT NONE

    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931

    ! Values
    ! S		: 2D dp, overlap matrix
    ! a,b	: int, set number we're on
    ! la,lb	: 1D int, angular quantum numbers 
    ! ta,tb	: int, tracking a and b
    ! basa,basb	: 1D dp, array of basis set weights for sets a and b
    ! seta,setb	: 1D int, array of setinfo for sets a and b
    
    ! Inout
    REAL(KIND=8), DIMENSION(0:,-2:,-2:,-2:), INTENT(IN) :: coef
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: S
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: basa, basb
    INTEGER, DIMENSION(0:), INTENT(IN) :: basinfo, seta, setb
    REAL(KIND=8), INTENT(IN) :: aa, bb, EIJ, p
    INTEGER, INTENT(IN) :: u, v, a, b

    ! internal
    INTEGER, DIMENSION(0:2) :: la,lb
    REAL(KIND=8) :: temp
    INTEGER :: i,j,k,orba,orbb,ori,prima,primb

    !update each element in set
    DO i=0,seta(0)-1 !go through set A

      orba = seta(3+i) !id of orbital

      DO j=0,setb(0)-1 !go through set B
        orbb = setb(3+j) !id of orbital

        ! get angular quantum numbers for each orbital 
        ori = basinfo(1+5*orba+3)
        !S-TYPE
        IF (ori .EQ. -1) THEN  
          la = [basinfo(1+5*orba+2),basinfo(1+5*orba+2),basinfo(1+5*orba+2)]
        !P-TYPE
        ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN 
          la = [0, 0, 0]
          la(ori) = basinfo(1+5*orba+2)
        END IF

        ori = basinfo(1+5*orbb+3)
        !S-TYPE
        IF (ori .EQ. -1) THEN 
          lb = [basinfo(1+5*orbb+2),basinfo(1+5*orbb+2),basinfo(1+5*orbb+2)]
        !P-TYPE
        ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN 
          lb = [0,0,0]
          lb(ori) = basinfo(1+5*orbb+2)
        END IF

        ! update Suv 
        temp = EIJ*(Pi/p)**(3.0D0/2.0D0)*basa(i)*basb(j)        ! WORK NOTE - hardcoded in bas
        temp = temp * gtoD(basinfo(1+5*orba+2),aa)                !basis set coefficients
        temp = temp * gtoD(basinfo(1+5*orbb+2),bb)                !basis set coefficeints
        temp = temp * coef(0,0,la(0),lb(0))*coef(1,0,la(1),lb(1))*coef(2,0,la(2),lb(2))

        S(orba,orbb) = S(orba,orbb) + temp
        
      END DO 
    END DO

  END SUBROUTINE overlap

!---------------------------------------------------------------------
!	Calculate kinetic energy of a pair of orbitals
!---------------------------------------------------------------------
  SUBROUTINE kinetic(Fb,u,v,a,b,p,bas,basinfo,coef,seta,setb,aa,bb,EIJ)
    IMPLICIT NONE

    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931

    ! na,nb	: 1D int, list of angular quantum number, named n for stupid reasons

    ! inout
    REAL(KIND=8), DIMENSION(0:,-2:,-2:,-2:), INTENT(IN) :: coef
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Fb
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: seta, setb
    REAL(KIND=8), INTENT(IN) :: aa, bb, EIJ, p
    INTEGER, INTENT(IN) :: u, v, a, b

    !internal
    REAL(KIND=8) :: temp,val
    INTEGER, DIMENSION(0:2) :: na,nb
    INTEGER :: i,j,k,orba,orbb,ori,ta,tb,prima,primb

    DO i=0,seta(0)-1
      orba = seta(2+i)

      DO j=0,setb(0)-1
        orbb = setb(2+j)

        temp = 0.0D0
        val = 0.0D0

        ori = basinfo(u,4*(orba+1)+2)
    
        !S-TYPE
        IF (ori .EQ. -1) THEN  
          na = [basinfo(u,4*(orba+1)+1),basinfo(u,4*(orba+1)+1),basinfo(u,4*(orba+1)+1)]
        !P-TYPE
        ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN 
          na = [0, 0, 0]
          na(ori) = basinfo(u,4*(orba+1)+1)
        END IF

        ori = basinfo(v,4*(orbb+1)+2)

        !S-TYPE
        IF (ori .EQ. -1) THEN  
          nb = [basinfo(v,4*(orbb+1)+1),basinfo(v,4*(orbb+1)+1),basinfo(v,4*(orbb+1)+1)]
        !P-TYPE
        ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN 
          nb = [0,0,0]
          nb(ori) = basinfo(v,4*(orbb+1)+1)
        END IF

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
        val = val * bas(u,a,2+i)*bas(v,b,2+j)     !basis set weights
        val = val * gtoD(basinfo(u,4*(orba+1)+1),aa) !primative constants 
        val = val * gtoD(basinfo(v,4*(orbb+1)+1),bb) !primative constants 

        Fb(orba,orbb) = Fb(orba,orbb) + val

      END DO
    END DO

!    WRITE(*,*) "==========="
!    WRITE(*,*) "given", a,b,seta

  END SUBROUTINE kinetic

!---------------------------------------------------------------------
!	Calculate the coulomb potential of a gaussian of two orbitals
!---------------------------------------------------------------------
  SUBROUTINE coulomb(Fb,u,v,a,b,ap,bas,basinfo,PP,seta,setb,aa,bb,atoms,EIJ,coef,lenA,lenB,lmaxA,lmaxB)
    IMPLICIT NONE
    ! Values
    ! Fb	: 2D dp, nuclear block of Fock Matrix
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
    ! lenA,lenB	: int, length of sets in A and B
    ! lmaxA	: 1D int, max ang quantum number of set A, B

    ! inout
    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931
    REAL(KIND=8), DIMENSION(0:,-2:,-2:,-2:), INTENT(IN) :: coef
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Fb
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: PP
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: seta, setb, atoms, lmaxA, lmaxB
    REAL(KIND=8), INTENT(IN) :: aa, bb, ap, EIJ
    INTEGER, INTENT(IN) :: u, v, a, b, lenA, lenB
    
    !internal
    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: Rtab
    LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: Rbol
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: nucpos
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Fj,CP
    INTEGER, DIMENSION(0:2) :: na, nb
    REAL(KIND=8) :: temp,val,TT,hmm
    INTEGER :: c,p,i,j,k,N,L,M,Nmax,Lmax,Mmax,nnuc,dummy,ori,orba,orbb

    val = 0.0D0
    nnuc = SIZE(atoms)

    Nmax = lmaxA(0) + lmaxB(0) 
    Lmax = lmaxA(1) + lmaxB(1) 
    Mmax = lmaxA(2) + lmaxB(2) 

    ALLOCATE(Fj(0:Nmax+Lmax+Mmax))
    ALLOCATE(CP(0:2))
    ALLOCATE(Rtab(-2:Nmax,-2:Lmax,-2:Mmax,0:Nmax+Lmax+Mmax))
    ALLOCATE(Rbol(-2:Nmax,-2:Lmax,-2:Mmax,0:Nmax+Lmax+Mmax))
    ALLOCATE(nucpos(0:nnuc-1,0:2))

    ! get nuclear positions
    OPEN(unit=1,file='nucpos',status='old',access='sequential')
    DO c=0,nnuc-1
      READ(1,*) dummy, nucpos(c,0:2) 
    END DO 
    CLOSE(unit=1)

!    WRITE(*,*) "at a, b",a,b, lmaxA(0), lmaxB(0)
    
    DO c=0,nnuc-1                         ! loop through atoms

      DO i=0,2                            ! construct PC
        CP(i) = nucpos(c,i) - PP(i)
      END DO
      TT = ap*(CP(0)**2.0D0 + CP(1)**2.0D0 + CP(2)**2.0D0)

      DO i=0,Nmax+Lmax+Mmax               ! get Boys table
        Fj(i) = 0.0D0 
      END DO
      CALL Boys(Fj,Nmax+Lmax+Mmax,TT)

       DO i=-2,Nmax
         DO j=-2,Lmax
           DO k=-2,Mmax
             DO p=0,Nmax+Lmax+Mmax
               Rtab(i,j,k,p) = 0.0D0
               Rbol(i,j,k,p) = .FALSE.
             END DO
           END DO
         END DO
       END DO
     
      ! WORK NOTE - check that it is correct to use ABS here with the weird directions of the kinetic integrals
      DO i=0,Nmax
        DO j=0,Lmax
          DO k=0,Mmax
!            CALL RNLMj(ABS(PP(0)),ABS(PP(1)),ABS(PP(2)),i,j,k,0,ap,Fj,Rtab,Rbol)
!            CALL RNLMj(PP(0),PP(1),PP(2),i,j,k,0,ap,Fj,Rtab,Rbol)
!            CALL RNLMj(ABS(CP(0)),ABS(CP(1)),ABS(CP(2)),i,j,k,0,ap,Fj,Rtab,Rbol)
!            CALL RNLMj(CP(0),CP(1),CP(2),i,j,k,0,ap,Fj,Rtab,Rbol)
            CALL RNLMj(-CP(0),-CP(1),-CP(2),i,j,k,0,ap,Fj,Rtab,Rbol)
          END DO
        END DO
      END DO

      !loop over set elements
      DO i=0,seta(0)-1                   ! loop over A
        orba = seta(2+i)

        DO j=0,setb(0)-1                 ! loop over B
          orbb = setb(2+j)

          IF (orba .EQ. 0 .AND. orbb .EQ. 5) THEN
            WRITE(*,*) CP(:)

          END IF

          temp = 0.0D0
          val = 0.0D0
 
          !get ang max for each orbital within set
          ori = basinfo(u,4*(orba+1)+2)
          !S-TYPE
          IF (ori .EQ. -1) THEN  
           na = [basinfo(u,4*(orba+1)+1),basinfo(u,4*(orba+1)+1),basinfo(u,4*(orba+1)+1)]
          !P-TYPE
          ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN 
            na = [0, 0, 0]
            na(ori) = basinfo(u,4*(orba+1)+1)
          END IF

          ori = basinfo(v,4*(orbb+1)+2)
          !S-TYPE
          IF (ori .EQ. -1) THEN  
            nb = [basinfo(v,4*(orbb+1)+1),basinfo(v,4*(orbb+1)+1),basinfo(v,4*(orbb+1)+1)]
          !P-TYPE
          ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN 
            nb = [0,0,0]
            nb(ori) = basinfo(v,4*(orbb+1)+1)
          END IF

          ! loop over all possible solutions
          DO N=0,na(0)+nb(0)
            DO L=0, na(1)+nb(1)
              DO M=0,na(2)+nb(2)
                temp = temp + coef(0,N,na(0),nb(0))*coef(1,L,na(1),nb(1))*coef(2,M,na(2),nb(2))*Rtab(N,L,M,0)
              END DO
            END DO
          END DO

          temp = temp * (2.0D0*Pi/ap)                     !from Boys
          temp = temp * bas(u,a,2+i)*bas(v,b,2+j)         !basis set weights
          temp = temp * gtoD(basinfo(u,4*(orba+1)+1),aa)  !primative constants
          temp = temp * gtoD(basinfo(v,4*(orbb+1)+1),bb)  !primative constants
          temp = temp * atoms(c)*EIJ                      !proton number and overlap coefficient
          
          Fb(orba,orbb) = Fb(orba,orbb) - temp

        END DO                           ! end loop over B
      END DO                             ! end loop over A
    END DO                               ! end nuclei loop

    DEALLOCATE(Fj)
    DEALLOCATE(CP)
    DEALLOCATE(Rtab)
    DEALLOCATE(Rbol)
    DEALLOCATE(nucpos)

  END SUBROUTINE coulomb

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
END PROGRAM int1e
