PROGRAM test
  USE aux
  IMPLICIT NONE 

REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931
  REAL(KIND=8), DIMENSION(0:5) :: Fj
  REAL(KIND=8) :: val
  INTEGER :: i 

  DO i=0,SIZE(Fj)-1
    Fj(i) = 0.0D0
  END DO

!  WRITE(*,*) "Testing factorial of 6, 16, 20"
!  WRITE(*,*) factR8(6), factR8(16), factR8(20), factR8(21)
 
  val = 21.0D0
  DO i=20,2,-1
    val = val*i
  END DO
  WRITE(*,*) val
  val = SQRT(2*Pi*21)*(21/EXP(1.0D0))**21
!  WRITE(*,*) val
  
  WRITE(*,*) "Testing Boys with T = 6"
  CALL Boys(Fj,5,6.0D0)
  WRITE(*,*) Fj
  WRITE(*,*) 

END PROGRAM test
