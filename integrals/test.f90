PROGRAM test
  USE aux
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:5) :: Fj
  REAL(KIND=8), DIMENSION(0:22,0:120) :: Ftab
  
  WRITE(*,*) "Testing Boys with T = 6"
  CALL Boys(Fj,5,6.0)
  WRITE(*,*) Fj
  WRITE(*,*) 

  OPEN(unit=1,file='Ftab',status='old',access='sequential'
  READ(1,*) Ftab
  CLOSE(unit=1)

END PROGRAM test
