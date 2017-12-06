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

    ! Case 1
    ! if 0 < T < 12 and 0 <= j <= J, use 7 term Taylor expansion
    IF (T .GT. 0 .AND. T .LT. 12) THEN
      WRITE(*,*) "case1"
      CALL Boys1(Fj,Q,T)
    ! for 12 < T < 2*J+36, use  
    ELSE IF (T .GE. 12 .AND. T .LT. 2*Q+36) THEN
      WRITE(*,*) "case2"
      CALL Boys2(Fj,Q,T)
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

    ! inout
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Fj
    REAL(KIND=8), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: Q

    ! internal

  END SUBROUTINE Boys1
!----------
  ! case 2 of boys function, uses upwards recursion scheme
  SUBROUTINE Boys2(Fj,Q,T)
    IMPLICIT NONE

    ! inout
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Fj
    REAL(KIND=8), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: Q

    ! internal

  END SUBROUTINE Boys2
!----------
  ! case 3 of boys function, uses simplified upwards recursion scheme
  SUBROUTINE Boys3(Fj,Q,T)
    IMPLICIT NONE

    ! inout
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Fj
    REAL(KIND=8), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: Q

    ! internal

  END SUBROUTINE Boys3
!----------

!===================================================================
!			FUNCTIONS
END MODULE aux
