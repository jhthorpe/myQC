PROGRAM test
  USE aux
  IMPLICIT NONE 

  REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: Rtab
  LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: Rbol
  REAL(KIND=8), DIMENSION(:),ALLOCATABLE :: Fj
  REAL(KIND=8) :: val,T,al,a,b,c
  INTEGER :: i,j,k,p,N,L,M

  al = 0.75
  N = 2
  L = 1
  M = 1
  a = 0.5
  b = 0.5
  c = 0.5
  T = al*(a**2+b**2+c**2)

  ALLOCATE(Fj(0:N+L+M))
  ALLOCATE(Rtab(-2:N,-2:L,-2:M,0:N+L+M))
  ALLOCATE(Rbol(-2:N,-2:L,-2:M,0:N+L+M))

  DO i=-1,N
    DO j=-1,L
      DO k=-1,M
        DO p=0,N+L+M
          Rtab(i,j,k,p) = 0.0D0
          Rbol(i,j,k,p) = .FALSE.
        END DO
      END DO
    END DO
  END DO

!  WRITE(*,*) RNLM(a,b,c,T,N,L,M)
  
  CALL Boys(Fj,N+L+M,T)

  !WRITE(*,*) (-2*al)**0*Fj(0)
  !WRITE(*,*) c*(-2*al)**1*Fj(1)
  !WRITE(*,*) b*(c*(-2*al)**1*Fj(1))
  
  CALL RNLMj(a,b,c,N,L,M,0,al,Fj,Rtab,Rbol)
  !CALL RNLMj(a,b,c,0,0,0,0,al,Fj,Rtab,Rbol)

  WRITE(*,*) Rtab(2,1,1,0)
  

END PROGRAM test
