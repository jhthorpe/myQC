!//////////////////////////////////////////////////////////////////
!//            Auxilary Functions for myQC 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//             
!//              Functions a la McMurchie and Davidson 
!//             
!///////////////////////////////////////////////////////////////////

MODULE aux
  IMPLICIT NONE
 
  CONTAINS

!===================================================================
!			SUBROUTINES

!---------------------------------------------------------------------
!	Generate table of auxilary functions, RNLMj 
!---------------------------------------------------------------------
  RECURSIVE SUBROUTINE RNLMj(a,b,c,N,L,M,j,al,Fj,Rtab,Rbol)
    IMPLICIT NONE
    ! Values
    ! a,b,c	: x,y,z positions
    ! N,L,M	: x,y,z angular momentum numbers
    ! al	: dp, alpha
    ! Fj	: 1d dp table of Fj values
    ! Rtab	: 4d dp table of R matrix values
    ! Rbol	: 4d bool table of if RNLMj has been seen or not

    !Inout
    REAL(KIND=8), DIMENSION(-2:,-2:,-2:,0:), INTENT(INOUT) :: Rtab
    LOGICAL, DIMENSION(-2:,-2:,-2:,0:), INTENT(INOUT) :: Rbol
    REAL(KIND=8), DIMENSION(0:) :: Fj
    REAL(KIND=8), INTENT(IN) :: a,b,c,al
    INTEGER, INTENT(IN) :: N,L,M,j

    ! 1) base cases
    ! If we've seen this before
    IF (Rbol(N,L,M,j)) THEN
      RETURN
    END IF
    ! check for "bad" point
    IF (N .LT. 0 .OR. L .LT. 0 .OR. M .LT. 0) THEN
      Rtab(N,L,M,j) = 0.0D0
      Rbol(N,L,M,j) = .TRUE.
      RETURN
    ! base case
    ELSE IF (N .EQ. 0 .AND. L .EQ. 0 .AND. M .EQ. 0) THEN
      Rtab(N,L,M,j) = (-2.0D0*al)**j*Fj(j)  
      Rbol(N,L,M,j) = .TRUE.
      RETURN
    END IF

    ! 2) non-base case    
    IF (N .NE. 0) THEN
      CALL RNLMj(a,b,c,N-1,L,M,j+1,al,Fj,Rtab,Rbol) 
      CALL RNLMj(a,b,c,N-2,L,M,j+1,al,Fj,Rtab,Rbol)
      Rtab(N,L,M,j) = a*Rtab(N-1,L,M,j+1) + (N-1)*Rtab(N-2,L,M,j+1)
      Rbol(N,L,M,j) = .TRUE.
      RETURN
    ELSE IF (L .NE. 0) THEN
      CALL RNLMj(a,b,c,0,L-1,M,j+1,al,Fj,Rtab,Rbol) 
      CALL RNLMj(a,b,c,0,L-2,M,j+1,al,Fj,Rtab,Rbol)
      Rtab(0,L,M,j) = b*Rtab(0,L-1,M,j+1) + (L-1)*Rtab(0,L-2,M,j+1)
      Rbol(0,L,M,j) = .TRUE.
      RETURN
    ELSE IF (M .NE. 0) THEN
      CALL RNLMj(a,b,c,0,0,M-1,j+1,al,Fj,Rtab,Rbol) 
      CALL RNLMj(a,b,c,0,0,M-2,j+1,al,Fj,Rtab,Rbol)
      Rtab(0,0,M,j) = c*Rtab(0,0,M-1,j+1) + (M-1)*Rtab(0,0,M-2,j+1)
      Rbol(0,0,M,j) = .TRUE.
      RETURN
    ELSE 
      WRITE(*,*) "Somehow you broke RNLMj recursion in auxilary.f90. N,L,M,j =", N,L,M,j
      STOP "Bad value auxilary:RNLMj"
    END IF

  END SUBROUTINE RNLMj

!---------------------------------------------------------------------
!	Control scheme for Boys integrals
!---------------------------------------------------------------------
  SUBROUTINE Boys(Fj,Q,T)
    IMPLICIT NONE
    ! Values
    ! Fj	: 1d dp, table for results of integrals
    ! Q		: int, max j value of Boys integral
    ! T		: dp, T value of Boys integral

    ! Inout
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Fj
    REAL(KIND=8), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: Q

    ! Internal
    REAL(KIND=8) :: val

    ! Catch a potential error
    IF (Q .GT. 16) THEN
      WRITE(*,*) "WARNING WARNING WARNING"
      WRITE(*,*) "J in Boys function > 16"
      WRITE(*,*) "WARNING WARNING WARNING"
      STOP
    END IF

    ! Case 1
    ! if 0 < T < 12 and 0 <= j <= J, use 7 term Taylor expansion
    IF (T .GE. 0 .AND. T .LT. 12) THEN
      CALL Boys1(Fj,Q,T)
    ! for 12 < T < 2*J+36, use  
    ELSE IF (T .GE. 12 .AND. T .LT. 2*Q+36) THEN
      CALL Boys2(Fj,Q,T)
    ! for T < 2*J+36
    ELSE IF (T .GE. 2*Q+36) THEN
      CALL Boys3(Fj,Q,T)
    ELSE 
      WRITE(*,*) "in auxilary.f90, T= ", T
      STOP "Bad logic in auxilary:Boys"
    END IF 
   
  END SUBROUTINE Boys 

!---------------------------------------------------------------------
!	Case 1 of Boys, uses downwards recursion scheme
!---------------------------------------------------------------------
  SUBROUTINE Boys1(Fj,Q,T)
    IMPLICIT NONE
    !Values
    ! Ftab	: 2d dp, table of FJ+k(T*) values, precalculated in Mathmatica, stored T,J

    ! inout
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Fj
    REAL(KIND=8), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: Q

    ! internal
    REAL(KIND=8), DIMENSION(0:120,0:22) :: Ftab
    INTEGER :: j,k,Tk

    !read in Ftab
    OPEN(unit=1,file='Ftab',status='old',access='sequential',form='unformatted')
    READ(1) Ftab
    CLOSE(unit=1)

    Fj = (/ (0.0D0, j=0,SIZE(Fj)-1) /)

    ! Get FJ
    Tk = NINT(T*10)
    DO k=0,6
      Fj(Q) = Fj(Q) + Ftab(Tk,Q+k)*(Tk/10.0D0 - T)**k/factR8(k) 
    END DO

    ! Downwards Recursion
    DO j=Q-1,0,-1
      Fj(j) = (2*T*Fj(j+1) + EXP(-T))/(2.0D0*j+1.0D0)
    END DO

  END SUBROUTINE Boys1

!---------------------------------------------------------------------
!	Case 2 of Boys, uses upwards recursion scheme
!---------------------------------------------------------------------
  SUBROUTINE Boys2(Fj,Q,T)
    IMPLICIT NONE
    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931

    ! inout
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Fj
    REAL(KIND=8), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: Q

    ! internal
    REAL(KIND=8) :: g
    INTEGER :: j

    ! Get g and F0(T)
    g = BoysG(T)
    Fj(0) = 0.5D0*Pi**(0.5D0)*T**(-0.5D0) - EXP(-T)*g/T

    ! Upwards recursion
    DO j=1,Q
      Fj(j) = (2.0D0*T)**(-1.0D0)*((2.0D0*(j-1) + 1)*Fj(j-1) - EXP(-T))
    END DO

  END SUBROUTINE Boys2

!---------------------------------------------------------------------
!	Case 3 of Boys, uses simplified upwards recursion scheme
!---------------------------------------------------------------------
  SUBROUTINE Boys3(Fj,Q,T)
    IMPLICIT NONE
    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931

    ! inout
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Fj
    REAL(KIND=8), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: Q

    ! internal
    INTEGER :: j


    ! no g in this form
    Fj(0) = 0.5D0*Pi**(0.5D0)*T**(-0.5D0)

    !simplified recursive form
    DO j=1,Q
      Fj(j) = (2.0D0*T)**(-1.0D0)*(2.0D0*(j-1)+1.0D0)*Fj(j-1) 
    END DO

  END SUBROUTINE Boys3

!===================================================================
!			FUNCTIONS

!---------------------------------------------------------------------
!	Explicit version of RNLM, currently incorrect and outdated
!---------------------------------------------------------------------
  REAL(KIND=8) FUNCTION RNLM(a,b,c,T,N,L,M)
    IMPLICIT NONE
    ! a,b,c	: dp, x,y,z coordinates of center
    ! N,L,M	: int, x,y,z angular quantum numbers
    
    ! Inout
    REAL(KIND=8), INTENT(IN) :: a,b,c,T
    INTEGER, INTENT(IN) :: N,L,M
    ! Internal
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Fj
    REAL(KIND=8) :: val,temp
    INTEGER :: i,j,k

    STOP "RNLM has not been checked yet, possibly using wrong F"
    ALLOCATE(Fj(0:(N+L+M)))

    val = 0.0D0

    CALL Boys(Fj,N+L+M,T)

    DO i=0,FLOOR(N/2.0D0)
      DO j=0,FLOOR(L/2.0D0)
        DO k=0,FLOOR(M/2.0D0)
          temp = a**(N-2*i)*b**(L-2*j)*c**(M-2*k)
          temp = temp * factR8(N)/(dfactR8(2*i)*factR8(N-2*i)) 
          temp = temp * factR8(L)/(dfactR8(2*j)*factR8(L-2*j))
          temp = temp * factR8(M)/(dfactR8(2*k)*factR8(M-2*k))
          temp = temp * Fj(N+L+M-i-j-k)
          val = val + temp 
        END DO
      END DO
    END DO
    
    RNLM = val 

    DEALLOCATE(Fj)

  END FUNCTION RNLM

!---------------------------------------------------------------------
!		Calculate g value for Boys integrals
!---------------------------------------------------------------------
  REAL(KIND=8) FUNCTION BoysG(T)
    IMPLICIT NONE

    ! Inout
    REAL(KIND=8), INTENT(IN) :: T

    IF (T .GE. 12.0D0 .AND. T .LT. 15.0D0) THEN
      BoysG = 0.4999489092-0.2473631686*T**(-1.0D0)+0.321180909*T**(-2.0D0) &
      -0.3811559346*T**(-3.0D0)
    ELSE IF (T .GE. 15.0D0 .AND. T .LT. 18.0D0) THEN
      BoysG = 0.4998436875-0.24249438*T**(-1.0D0)+0.24642845*T**(-2.0D0)
    ELSE IF (T .GE. 18.0D0 .AND. T .LT. 24.0D0) THEN
      BoysG = 0.499093162-0.2152832*T**(-1.0D0)
    ELSE IF (T .GE. 24.0D0 .AND. T .LT. 30.0D0) THEN
      BoysG = 0.490D0
    ELSE IF (T .LT. 12.0D0) THEN
      WRITE(*,*) "BoysG in auxilary.f90 was called incorrectly. T =", T
      STOP "Bad call to auxilary:BoysG"
    END IF 

  END FUNCTION BoysG

!---------------------------------------------------------------------
!			Real(kind=8) factorial
!---------------------------------------------------------------------
  REAL(KIND=8) FUNCTION factR8(n)
    IMPLICIT NONE
    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931
    ! Inout
    INTEGER, INTENT(IN) :: n
    ! Internal
    REAL(KIND=8) :: val
    INTEGER :: i

    ! Standard factorial
    !IF (n .LT. 21) THEN
    val = 1.0D0
    DO i=n,2,-1
      val = val*i
    END DO
    ! Sterling for larger factorials 
    !ELSE 
    !  val = SQRT(2*Pi*n)*(n/EXP(1.0D0))**n
    !END IF

    factR8 = val

  END FUNCTION factR8

!---------------------------------------------------------------------
!		Real(kind=8) double factorial
!---------------------------------------------------------------------
  REAL(KIND=8) FUNCTION dfactR8(n)
    IMPLICIT NONE
    !inout
    INTEGER, INTENT(IN) :: n
    !internal
    REAL(KIND=8) :: val
    INTEGER :: i
  
    val = 1.0D0
 
    IF (n .EQ. 0 .OR. n .EQ. -1) THEN
      val = 1.0D0
    ELSE IF (MOD(n,2) .EQ. 0) THEN
      DO i=n,2,-2
        val = val * i
      END DO
    ELSE IF (MOD(n,2) .NE. 0) THEN 
      DO i=n,3,-2
        val = val * i
      END DO
    ELSE
      WRITE(*,*) "Somehow you broke the double factorial, n=", n
      STOP "bad n!! in auxilary:dfactR8"
    END IF

    dfactR8 = val

  END FUNCTION dfactR8

!---------------------------------------------------------------------
!               Calculate coefficients of overlap Gaussians 
!---------------------------------------------------------------------
  SUBROUTINE getcoef(M,PA,PB,aa,bb,amax,bmax)
    IMPLICIT NONE

    ! M         : 2D dp, coefficient matrix d_{i}^{j,k} = M(0,i,j,k)
    ! PA        : 1D dp, (xa,ya,za) distance from central gaussian
    ! PB        : 1D dp, (xb,yb,zb) distance from central gaussian
    ! aa        : dp, coefficient of atom A 
    ! bb        : dp, coefficient of atom B
    ! amax      : 1D int, max angular quantum number of A needed
    ! amax      : 1D int, max angular quantum number of B needed
    ! nmax      : int, max total quantum number needed
    ! nn        : int, total number of elements
    ! fmat      : 3D dp, boolean array, tracks what has been calculated 

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
!       Left side (A) recursion on coefficients Hermite Gaussians
!---------------------------------------------------------------------
  RECURSIVE SUBROUTINE lrec(M,l,i,j,k,PA,PB,pp,fmat)
    IMPLICIT NONE

    ! M         : 4D dp, matrix of coefficients
    ! i,j,k     : int, index of coefficient we need. i=N,j=n,k=nbar
    ! l         : int, coordinite of coefficient we need 0=x,1=y,2=z)
    ! PA, PB    : 1D dp, list of x,y,z distances of atoms A,B from P
    ! fmat      : 3D dp, list of which values we already have
    ! pp        : dp, value of sum of coefficients

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
!       Right side (B) recursion on coefficients of Hermite Gaussians
!---------------------------------------------------------------------
  RECURSIVE SUBROUTINE rrec(M,l,i,j,k,PA,PB,pp,fmat)
    IMPLICIT NONE

    ! M         : 4D dp, matrix of coefficients
    ! i,j,k     : int, index of coefficient we need. i=N,j=n,k=nbar
    ! l         : int, coordinate we need (0=x,1=y,2=z)
    ! PA, PB    : 1D dp, list of x,y,z distances of atoms A,B from P
    ! fmat      : 3D dp, list of which values we already have
    ! pp        : dp, value of sum of coefficients

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
!		Compute nonzero combinations of N,L,M in sets 
!---------------------------------------------------------------------
  SUBROUTINE getDk(coef,seta,setb,basa,basb,basinfo,Dk,Ck,Ok,kmax,EIJ,setl,aa,bb)
    IMPLICIT NONE

    !Values
    ! coef	: 4D dp, matrix of coefficients
    ! seta,b	: 1D int, setinfo for seta,setb
    ! basa,b	: 1D dp, basis set weights for sets a,b 
    ! basinfo	: 1D int, basinfo for sets a,b
    ! Dk	: 1D dp, nonzero EIJ*Nk*Lk*Mk values
    ! Ck	: 1D int, tracks Nk,Lk,Mk combinations 
    ! Ok	: 1D int, gives orbs of the combination {left orb, right orb}
    ! kmax	: int, max index of nonzero values
    ! EIJ	: dp, preexponential value
    ! setl	: int, length of set
    ! ll,rr	: 1D int, left,right angular quantum numbers

    !Inout
    REAL(KIND=8), DIMENSION(0:,-2:,-2:,-2:), INTENT(IN) :: coef
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Dk
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: basa,basb
    INTEGER, DIMENSION(0:), INTENT(INOUT) :: Ck,Ok
    INTEGER, DIMENSION(0:), INTENT(IN) :: seta,setb,basinfo
    REAL(KIND=8), INTENT(IN) :: aa,bb,EIJ
    INTEGER, INTENT(INOUT) :: kmax
    INTEGER :: setl

    !Internal
    INTEGER, DIMENSION(0:2) :: ll,rr 
    REAL(KIND=8) :: tol
    INTEGER :: N,L,M,Nmax,Lmax,Mmax,ori
    INTEGER :: i,j,orba,orbb
   
    tol = 0.1D-15
 
    !go through orbitals of set a 
    DO i=0, seta(0)-1
      orba = seta(3+i)

      ! get angular quantum numbers for each orbital 
      ori = basinfo(1+5*orba+3)
      !S-TYPE
      IF (ori .EQ. -1) THEN
        ll = [basinfo(1+5*orba+2),basinfo(1+5*orba+2),basinfo(1+5*orba+2)]
      !P-TYPE
      ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN
        ll = [0, 0, 0]
        ll(ori) = basinfo(1+5*orba+2)
      END IF

      !go through orbitals of set b
      DO j=0, setb(0)-1
        orbb = setb(3+j)

        ori = basinfo(1+5*orbb+3)
        !S-TYPE
        IF (ori .EQ. -1) THEN
          rr = [basinfo(1+5*orbb+2),basinfo(1+5*orbb+2),basinfo(1+5*orbb+2)]
        !P-TYPE
        ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN
          rr = [0,0,0]
          rr(ori) = basinfo(1+5*orbb+2)
        END IF

        Nmax = ll(0) + rr(0)
        Lmax = ll(1) + rr(1)
        Mmax = ll(2) + rr(2)

        DO M=0,Mmax
          DO L=0,Lmax
            IF (ABS(coef(1,L,ll(1),rr(1))) .LT. tol) CYCLE
            DO N=0,Nmax
              IF (ABS(coef(0,N,ll(0),rr(0))) .LT. tol) THEN
                CYCLE
              ELSE
                kmax = kmax + 1
                Ck(kmax) = 300*N + 20*L + 1*M
                Dk(kmax) = coef(0,N,ll(0),rr(0))*coef(1,L,ll(1),rr(1))*coef(2,M,ll(2),rr(2))
                Dk(kmax) = Dk(kmax)*EIJ*gtoD(basinfo(1+5*orba+2),aa)*gtoD(basinfo(1+5*orbb+2),bb) 
                Dk(kmax) = Dk(kmax)*basa(i)*basb(j)
                Ok(2*kmax) = orba
                Ok(2*kmax+1) = orbb
              END IF
            END DO
          END DO
        END DO

      END DO
    END DO

  END SUBROUTINE getDk

!=====================================================================
!                       FUNCTIONS

!---------------------------------------------------------------------
!       Generate coefficients of GTO basis sets 
!---------------------------------------------------------------------
  REAL(KIND=8) FUNCTION gtoD(l,a)
    IMPLICIT NONE
    
    ! l         : int, angular momentum quantum number
    ! a         : dp, coefficient
    
    ! Inout 
    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931
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

END MODULE aux
