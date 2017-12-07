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

  DO i=0,120
    CALL Boys(Fj,1,(i+0.3)/10.0D0)
    WRITE(*,*) (i+0.3)/10.0D0, Fj(0)
  END DO

END PROGRAM test
