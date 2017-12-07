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
      WRITE(*,*) "case1"
      CALL Boys1(Fj,Q,T)
    ! for 12 < T < 2*J+36, use  
    ELSE IF (T .GE. 12 .AND. T .LT. 2*Q+36) THEN
      WRITE(*,*) "case2"
      CALL Boys2(Fj,Q,T)
    ! for T < 2*J+36
    ELSE IF (T .GE. 2*Q+36) THEN
      WRITE(*,*) "case3"
      CALL Boys3(Fj,Q,T)
    ELSE 
      WRITE(*,*) "in auxilary.f90, Boys subroutine, bad logic"
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

    WRITE(*,*) Ftab(0,1)
    WRITE(*,*) Ftab(1,1)
    WRITE(*,*) Ftab(2,1)
    WRITE(*,*) Ftab(119,1)
    WRITE(*,*) Ftab(120,1)
    WRITE(*,*) Ftab(0,22)
    WRITE(*,*) Ftab(120,22)
    WRITE(*,*)

    ! Get FJ
!    Tk = NINT(T*10)
!    DO k=0,6
!      Fj(Q) = Fj(Q) + Ftab(Q+k,Tk)*(tk/10.0D0 - T)**k/factR8(k) 
!    END DO

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

    ! Recursive form 
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
    val = n
    DO i=n-1,2,-1
      val = val*i
    END DO
    ! Sterling for larger factorials 
    !ELSE 
    !  val = SQRT(2*Pi*n)*(n/EXP(1.0D0))**n
    !END IF

    factR8 = val

  END FUNCTION factR8
!-----------
END MODULE aux
