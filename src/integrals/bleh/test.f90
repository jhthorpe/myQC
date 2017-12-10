PROGRAM test
  USE aux
  IMPLICIT NONE 

  REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: Rtab
  LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: Rbol
  REAL(KIND=8), DIMENSION(:),ALLOCATABLE :: Fj
  REAL(KIND=8) :: val,T,al,a,b,c
  INTEGER :: i,j,k,p,N,L,M
  LOGICAL :: flag

  al = 74.558086700000004
  N = 0
  L = 0
  M = 1
  a = 0.00D0
  b = 0.00D0
  c = 0.99999433085150613
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

  CALL Boys(Fj,N+M+L,T)
  CALL RNLMj(a,b,c,N,L,M,0,al,Fj,Rtab,Rbol)

  WRITE(*,*) "F1",Fj(1)
  WRITE(*,*) "R", Rtab(0,0,1,0)

!  WRITE(*,*) RNLM(a,b,c,T,N,L,M)
 
!  DO i=0,320
!    T = i*0.1
!    CALL Boys(Fj,N+L+M,T)
!    DO j=0,N+L+M
!      IF (Fj(j) .LE. 0.0 .OR. Fj(j) .GE. 1.0) WRITE(*,*) "T=", T, "j,Fj(j)", j, Fj(j)
!    END DO    
!  END DO 
!  CALL Boys(Fj,N+L+M,T)

!  WRITE(*,*) "j,Fj:", N+L+M,Fj(N+L+M)
!  WRITE(*,*) "j,Fj:", 1,Fj(1)
  


!  WRITE(*,*) "RNLM0:",Rtab(N,L,M,0)
  

END PROGRAM test
