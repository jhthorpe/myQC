! some simple testing
program tests
  implicit none
  
  REAL(KIND=8), DIMENSION(0:2,0:2) :: A
  REAL(KIND=8), DIMENSION(0:2,0:2) :: B
  REAL(KIND=8), DIMENSION(0:2) :: v
  INTEGER :: i,j,Iij

  !blasstuff
  REAL(KIND=8), DIMENSION(0:2) :: Y
  CHARACTER(LEN=1) :: TRANS
  REAL(KIND=8) :: ALPHA,BETA
  INTEGER :: M,N,LDA,INCX,INCY
  
 
  write(*,*) "Matrix order testing"
  
  A = 0
  Iij=0
  DO i=0,2
    DO j=i,2
      A(i,j) = Iij
      Iij = Iij + 1
    END DO
  END DO
  
  write(*,*) "A is..."
  DO i=0,2
!    WRITE(*,*) A(:,j)
     WRITE(*,*) A(i,:)
  END DO
  
  v = [1,2,3]
  WRITE(*,*) 
  WRITE(*,*) "v is..."
  WRITE(*,*) v

  !testing blas matvec
  TRANS='N'
  M=3
  N=3
  ALPHA=1.0D0
  LDA=3
  INCX=1
  BETA=0
  INCY=1 
  CALL DGEMV(TRANS,M,N,ALPHA,A,LDA,v,INCX,BETA,Y,INCY)

  WRITE(*,*)
  WRITE(*,*) "Y is..."
  WRITE(*,*) Y

  !--------------------------
  !unformatted write to intermediate
  B = 0

  WRITE(*,*) 
  WRITE(*,*) "B before is..."
  DO i=0,2
    WRITE(*,*) B(i,:)
  END DO

  OPEN(unit=101,file='inter',status='replace',access='sequential',form='unformatted')
  DO i=0,2
    WRITE(101) A(i,:)
  END DO
  CLOSE(unit=101,status='keep')

  OPEN(unit=101,file='inter',status='old',access='sequential',form='unformatted')
  DO i=0,2
    READ(101) B(i,:)
  END DO
  CLOSE(unit=101)
   
  WRITE(*,*) "B after is..."
  DO i=0,2
    WRITE(*,*) B(i,:)
  END DO

!------
!symmetric matrix stuff
  A(0,:) = [1,2,0]
  A(1,:) = [0,0,3]
  A(2,:) = [0,0,1]
!  A(1,:) = [2,0,3]
!  A(2,:) = [0,3,1]

  WRITE(*,*) 
  WRITE(*,*) "testing sym mat vec"
  WRITE(*,*) "A is"
  DO i=0,2
    WRITE(*,*) A(i,:)
  END DO  

  WRITE(*,*) 
  WRITE(*,*) "v is..."
  DO i=0,2
    WRITE(*,*) v(i)
  END DO

  CALL DSYMV('U',N,ALPHA,A,LDA,v,INCX,BETA,Y,INCY)

  WRITE(*,*) 
  WRITE(*,*) "Y is..."
  WRITE(*,*) Y

  

  

end program tests
