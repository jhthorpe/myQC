!/////////////////////////////////////////////////////////////////////
!//             Module containing linear algebra subroutines 
!//
!//                     James H Thorpe, in Group of John Stanton
!//                     The University of Florida
!//
!/////////////////////////////////////////////////////////////////////

! GENERAL STRUCTURE:
!       linal_operation_matType_dataType
!       example: Lanczos eigenvalues of a Hermitian 2D real*8 matrix 
!               linal_lanczos_hermitian_2Dreal8
!
!	IMPORTANT NOTE : here I am assuming column major arrays

MODULE linal
  CONTAINS

!---------------------------------------------------------------------
!	linal_lanczos_symreal_2Dreal8
!		James H. Thorpe
!		July 21, 2018
!	- Performs Lanczos to find lowest m eigenvalues of ...
!	  ... a (nxn) symmetric, real valued, 2D real8 array 
!	- takes V as an input vector of length n. Program checks if it...
!	  ... its Euclidean Norm is 1
!	- Currently uses my implementation of matrix vector multiplication
!	- DGEM would be prefered, and should perhaps be replaced
!	- T = S*AS
!	- on exit, Td hold the eigenvalues in increasing order
!	- on exit, V holds the eigenvectors in increasing order
!---------------------------------------------------------------------
  ! Values
  ! A		:	2D real8, nxn symetric, real valued array to be tridiagonalized
  ! n		:	int4, size of arrays
  ! m		:	int4, number of lowest eigenvalues to be found
  ! V		:	1D real8, length n input vector 
  ! W		:	1D real8, length n working vector
  ! X		:	1D real8, length n working vector
  ! S		:	2D real8, nxm (row x col) colomns of X's 
  ! Td		:	2D real8, m element vector of diagonals	
  ! Ts		:	1D real8, m-1 element vector of subdiagonals

  SUBROUTINE linal_lanczos_symreal_2Dreal8(A,n,m,V,W,X,S,Td,Ts)
    IMPLICIT NONE
    REAL(KIND=8), PARAMETER :: tol=1.0D-16
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: S
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: A    
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: V,W,X,Td,Ts
    INTEGER(KIND=4), INTENT(IN) :: n,m
    
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: evec
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: diag,work
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: iwork
    INTEGER(KIND=4) :: i,j,lwork,liwork,info
    REAL(KIND=8) :: b,c
    CHARACTER(LEN=1) :: fparam

    fparam='Y'

    Td = 0
    Ts = 0

    !check input vector has Euclidian norm of 1
    b = linal_eunorm_1Dreal8(V,n) 
    IF ( ABS(b - 1.0D0) .GT. tol) THEN
      V(0) = 1.0D0
      V(1:n-1) = (/ (0.0D0, i=1, n-1) /) 
    END IF

    !Initial step
    W = MATMUL(A,V)
    c = SUM(W*V) !this works because we are real valued
    X = W - c*V
    S(0,:) = V
    Td(0) = c

    !change this to while loop
    !start with given m, and climb by 5 iterations each time
    !check convergence every 5 iterations
    !not the best way to do this, but it should work
    DO i=1,m-1
      b = linal_eunorm_1Dreal8(X,n)
      IF (ABS(b - 0.0D0) .GT. tol) THEN
        V = X/b
      ELSE
        ! add a new vector that is orthonormal to the others
        CALL linal_onvec_2Dreal8(S(0:i-1,:),V,n,i)
      END IF
      W = MATMUL(A,V)
      c = SUM(W*V) !this works because we are real valued 
      X = W - c*V-b*S(i-1,:)
      S(i,:) = V
      Td(i) = c
      Ts(i-1) = b 
    END DO 

    !Now we have our tridiagonal matrix, we must get the eigenvalues
    !Currently using LAPACK, though this could probably be self coded

    !allocate memory and find best size for arrays
    ALLOCATE(evec(0:m-1,0:n-1))
    ALLOCATE(work(0:1))
    ALLOCATE(iwork(0:1))
    work = 0
    iwork = 0

    !one at a time?
    CALL DSTEDC('I',m,Td,Ts,evec,m,work,-1,iwork,2,info)
    WRITE(*,*) "work is:", work
    lwork = CEILING(MAX(2.0,work(0)))
    WRITE(*,*) "lwork is", lwork
    DEALLOCATE(work)
    ALLOCATE(work(0:lwork-1))

    CALL DSTEDC('I',m,Td,Ts,evec,m,work,lwork,iwork,-1,info)
    WRITE(*,*) "iwork is:", iwork 
    liwork = MAX(2,iwork(0))
    WRITE(*,*) "liwork is", liwork
    DEALLOCATE(iwork)
    ALLOCATE(iwork(0:liwork-1))
   
    !get eigenvalues
    CALL DSTEDC('I',m,Td,Ts,evec,m,work,lwork,iwork,liwork,info)

    !get eigenvectors
    

    DEALLOCATE(evec)
    DEALLOCATE(work)
    DEALLOCATE(iwork)
   
  END SUBROUTINE linal_lanczos_symreal_2Dreal8

!---------------------------------------------------------------------
!	linal_eunorm_1Dreal8
!		James H. Thorpe
!		July 22, 2018
!	- function that returns the euclidean norm of a real8 vector
!	- might consider changing this for greater numerical stability?
!	- certainly consider changing this for better performance
!---------------------------------------------------------------------
  ! Values
  ! A		:	1D real8, vector who's norm to return
  ! n		:	int4, lenth of vector

  REAL(KIND=8) FUNCTION linal_eunorm_1Dreal8(A,n)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: A
    INTEGER(KIND=4), INTENT(IN) :: n
    
    INTEGER(KIND=4) :: i
    REAL(KIND=8) :: temp

    temp = 0.0D0
    DO i=0, n-1
      temp = temp + A(i)**2.0D0 
    END DO 
   
    temp = SQRT(temp)
    linal_eunorm_1Dreal8 = temp

  END FUNCTION linal_eunorm_1Dreal8

!---------------------------------------------------------------------
!	linal_onvec_2Dreal8
!		James H. Thorpe
!		August 8, 2018
!	- subroutine that finds a vector orthonormal to set of vectors
!	- Array A stores m, n element vectors, and v is the output new vector
!	- uses Gram Schmidt process
!	- assumes that all vectors are themselves already orthonormal
!---------------------------------------------------------------------
  ! Values
  ! A(m,n)	:	2D real8, n-row m-col array of current set 
  ! v		:	1D real8, n element vector to be orthogonalized
  ! n		:	int4, number of rows
  ! m		:	int4, number of vectors we already have

  SUBROUTINE linal_onvec_2Dreal8(A,v,n,m)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: v
    INTEGER(KIND=4), INTENT(IN) :: n,m

    INTEGER(KIND=4) :: i
 
    IF (m .GT. n) THEN
      WRITE(*,*) "linal_onvec_2Dreal8 : WARNING - you are creating a linearly dependent set" 
      STOP 
    END IF

    v = 0.0D0
    v(m) = 1.0D0

    ! go through all currect vectors
    ! note that I *could* assume dot product of A(i) with itself is 1, but
    !   I want to account for numerical problems
    DO i=0,m-1
      v = v - SUM(A(i,:)*v)/SUM(A(i,:)*A(i,:))*A(i,:)
    END DO 

    !normalize
    v = v/SUM(v)

  END SUBROUTINE linal_onvec_2Dreal8

!---------------------------------------------------------------------
!	linal_proj_1Dreal8
!		James H. Thorpe
!		August 8, 2018
!	- subroutine that performs projection operation of vector v onto u
!	- result is stored in vector w
!---------------------------------------------------------------------
  ! Values
  ! u		:	1D real8, vector to be projected upon
  ! v		:	1D real8, vector that is projected
  ! w		:	1D real8, resulting vector of the projection
  ! n		:	int4, number of elements

  SUBROUTINE linal_proj_1Dreal8(u,v,w,n)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: w 
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: u,v
    INTEGER(KIND=4), INTENT(IN) :: n
    
    w = SUM(u*v)/SUM(u*u)*u

  END SUBROUTINE linal_proj_1Dreal8

!---------------------------------------------------------------------

END MODULE linal
