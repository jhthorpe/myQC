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

!----------
  ! Recursive subroutine that returns the table of auxilary functions, RNLMj 
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
!----------
! Boys control scheme for integrals
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
    IF (T .GT. 0 .AND. T .LT. 12) THEN
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
!----------
  ! case 1 of boys function, uses downwards recursion scheme
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
    OPEN(unit=1,file='Ftab',status='old',access='sequential')
    READ(1,*) Ftab
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
!----------
  ! case 2 of boys function, uses upwards recursion scheme
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
      Fj(j) = (2.0D0*T)**(-1.0D0)*((2.0D0*j + 1)*Fj(j-1) - EXP(-T))
    END DO

  END SUBROUTINE Boys2
!----------
  ! case 3 of boys function, uses simplified upwards recursion scheme
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
      Fj(j) = (2.0D0*T)**(-1.0D0)*(2.0D0*j+1.0D0)*Fj(j-1) 
    END DO

  END SUBROUTINE Boys3
!----------

!===================================================================
!			FUNCTIONS
!----------
  ! Function that return the explicit version of RNLM
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
!----------
  ! Function that returns the proper g value for Boys integrals
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
!----------
  ! Function that returns real(kind=8) factorial
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
!-----------
  ! Function that return real(kind=8) double factorial
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
!-----------

END MODULE aux
