PROGRAM test
  USE linal
  !USE IFPORT
  IMPLICIT NONE
  
  !lanczos
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: A,S,Anew,VL,VR
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: V,W,X,Y,Td,Ts,WR,WI
  INTEGER(KIND=4) :: n,m,i,j
  REAL(KIND=8) :: t1,t2
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: sigma,b

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: WORK,wshit
  INTEGER :: LWORK, INFO

  !gs
  !REAL(KIND=8), DIMENSION(0:2,0:2) :: B
  !REAL(KIND=8), DIMENSION(0:2) :: u,aba,gaba

  !mkl bullshit
  integer iseed /3/

  !testing lanczos
  CALL random_seed(iseed)
 
  ALLOCATE(A(0:1,0:1))
  ALLOCATE(WR(0:1))
  ALLOCATE(WI(0:1))
  ALLOCATE(VL(0:1,0:1))
  ALLOCATE(VR(0:1,0:1))
  ALLOCATE(WORK(0:10))
  A(0,0) = 1
  A(1,0) = 3
  A(0,1) = 3
  A(1,1) = 1
  CALL linal_printmat_2Dreal8(A,2,2)
  CALL DGEEV('N','V',2,A,2,WR,WI,VL,2,VR,2,WORK,11,i)
  WRITE(*,*) "Eigenvalues"
  WRITE(*,*) WR(0)
  WRITE(*,*) WR(1)
  WRITE(*,*) "Eigenvectors"
  !CALL linal_printmat_2Dreal8(VR,2,2)
  WRITE(*,*) VR
  WRITE(*,*) "orthonormal..."
  WRITE(*,*) VR(0,0:1)/SUM(VR(0,0:1)*VR(0,0:1))
  WRITE(*,*) VR(1,0:1)/SUM(VR(1,0:1)*VR(1,0:1))
  WRITE(*,*) VR(0:1,0)/SUM(VR(0:1,0)*VR(0:1,0))
  WRITE(*,*) VR(0:1,1)/SUM(VR(0:1,1)*VR(0:1,1))
  STOP

  n = 3
  m = 1

  ALLOCATE(A(0:n-1,0:n-1))
  ALLOCATE(Anew(0:n-1,0:n-1))
  ALLOCATE(S(0:m-1,0:n-1))
  ALLOCATE(Td(0:m-1))
  ALLOCATE(Ts(0:m-2))
  ALLOCATE(V(0:n-1))
  ALLOCATE(W(0:n-1))
  ALLOCATE(X(0:n-1))
  ALLOCATE(Y(0:n-1))

  LWORK=50000
  ALLOCATE(wshit(0:n-1))
  ALLOCATE(WORK(0:LWORK-1))

  !A = TRANSPOSE(RESHAPE((/ 93,57,93,57,2,12,93,12,19/),(/n,n/)))
  !A = TRANSPOSE(RESHAPE((/ 25,1,0,1,32,0,0,0,93/),(/n,n/)))

  WRITE(*,*) RAN(iseed)
  WRITE(*,*) RAN(iseed)

  
  DO i=0,n-1
    A(i,i) = (RAN(iseed)-0.5)*100.0D0
    DO j=i+1,n-1
      A(i,j) = (RAN(iseed)-0.5)*0.1D0
      A(j,i) = A(i,j) 
    END DO
  END DO
  OPEN(unit=100,file='mat',status='replace',form='unformatted')
  WRITE(100) A
  CLOSE(unit=100)

  WRITE(*,*) 
  WRITE(*,*) "A is:"
  DO i=0,n-1
    WRITE(*,*) A(:,i) 
  END DO

  WRITE(*,*) 
  WRITE(*,*) "-----------------"
  

  !OPEN(unit=100,file='mat',status='old',form='unformatted')
  !READ(100) A
  !CLOSE(unit=100)
 
  !DO i=0,n-1
  !    WRITE(*,*) A(:,i)
  !END DO

  V(0) = 1.0D0
  V(1:n-1) = (/ (0.0D0, i=1,n-1) /)

  call CPU_TIME(t1)
  CALL linal_lanczos_symreal_2Dreal8(A,n,m,V,W,X,S,Td,Ts)
  call CPU_TIME(t2)

  !WRITE(*,*) "S is:"
  !#do i=0,n-1
  !  write(*,*) S(:,i)
  !end do  
  
  write(*,*) 
  !WRITE(*,*) "eigenvalues are"
  !do i=0,m-1
  !  write(*,*) Td(i)
  !end do
  WRITE(*,*) "Smallest eigenvalue..."
  WRITE(*,*) MINVAL(Td(:))
  write(*,*)
  write(*,*) "eigenvectors are"
  
  write(*,*) 
  WRITE(*,*) "my code finished in:", t2-t1 

  call cpu_time(t1)
  Anew = A
  CALL dsyev('N','U',n,Anew,n,wshit,WORK,LWORK,INFO) 
  call cpu_time(t2)
  WRITE(*,*) "The actual values are"
  DO i=0,MIN(m-1,9)
    WRITE(*,*) wshit(i)
  END DO
  WRITE(*,*) 
  WRITE(*,*) "lapack finished in", t2-t1

  !//////////////////////////////////////////////////////////////////
  WRITE(*,*)
  WRITE(*,*) "Testing Davidson"
  ALLOCATE(sigma(0:n-1,0))
  ALLOCATE(b(0:n-1,0))
  b=0
  b(0,0) = 1.0D0
  sigma=0

  DO j=0,n-1
    sigma(j,0) = A(0,j)*b(j,0)
  !  write(*,*) A(0,j)*b(j,0)
  END DO

  !WRITE(*,*) "b is:"
  !WRITE(*,*) b(:,0)

  !WRITE(*,*) "sigma is:"
  !WRITE(*,*) sigma(:,0)
  !WRITE(*,*) sigma

  m = 1
  CALL linal_davidson_2Dreal8(A(0:n-1,0:n-1),sigma(0:n-1,0:m-1),b(0:n-1,0:m-1),m,n,1,10) 

  !//////////////////////////////////////////////////////////////////
  ! testing GS
  !WRITE(*,*) 
  !WRITE(*,*) "////// TESTING GS //////"

  !B = TRANSPOSE(RESHAPE((/1,2,3,0,0,0,0,0,0/),SHAPE(B)))
  
  !WRITE(*,*) "Currect array:"
  !DO i=0,SIZE(B(0,:))-1
  !    WRITE(*,*) B(:,i)
  !END DO

  !WRITE(*,*) "New vector 1"
  !CALL linal_onvec_2Dreal8(B,u,SIZE(w)-1,1)
  !DO i=0,SIZE(u)-1
  !  WRITE(*,*) u(i) 
  !END DO

  !B(1,:) = u

  !WRITE(*,*) "New Vector 2"
  !CALL linal_onvec_2Dreal8(B,u,SIZE(w)-1,2)
  !DO i=0,SIZE(u)-1
  !  WRITE(*,*) u(i) 
  !END DO

  !aba = [1,1,0]
  !gaba = [1,2,3]

  !WRITE(*,*) 
  !WRITE(*,*) SUM(aba*gaba)/SUM(gaba*gaba)*gaba 

  DEALLOCATE(wshit)
  DEALLOCATE(WORK)

END PROGRAM test
