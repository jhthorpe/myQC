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
    !WRITE(*,*) "lwork is", lwork
    DEALLOCATE(work)
    ALLOCATE(work(0:lwork-1))

    CALL DSTEDC('I',m,Td,Ts,evec,m,work,lwork,iwork,-1,info)
    WRITE(*,*) "iwork is:", iwork 
    liwork = MAX(2,iwork(0))
    !WRITE(*,*) "liwork is", liwork
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
!       linal_davidson_2Dreal8
!               James H. Thorpe
!               Nov 18, 2018
!       - subroutine that uses davidson's algorithm to iteratively obtain
!         the k'th eigenvalue from a matrix, A, expanded in subspace of b vectors
!       - uses BLAS and LAPACK for (some) internal linal
!       - label in line wth Davidson, JCP, 1975
!---------------------------------------------------------------------
  ! Values
  ! A           : 2D real8, matrix to be partially diagonalized, size [n,n]
  ! s(_in)      : 2D real8, sigma vectors of the expansion A.s, size [n,m]
  ! b(_in)      : 2D real8, b vectors forming the subspace, size[n,m]
  ! m           : int, size of subspace
  ! n           : int, size of A matrix
  ! k           : int, kth eigenvalue is the one we seek 
  ! texp        : int, exponential of tolerance
  ! An          : 2D real8, pseudo matrix A^~
  ! Ap          : 2D real8, copy of An for lapack 
  ! al          : 1D real8, eigenvector associated with kth eigenvalue
  ! la          : real8, k'th eigenvalue
  ! qm          : 1D real8, intermediate vector, size [n]
  ! qnrm        : real8, eucleian norm of qm
  ! xi          : 1D real9, intermediate vector, size [n]
 
  SUBROUTINE linal_davidson_2Dreal8(A,s_in,b_in,m,n,k,texp)
    IMPLICIT NONE
    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: A,s_in,b_in
    INTEGER(KIND=4), INTENT(INOUT) :: m
    INTEGER(KIND=4), INTENT(IN) :: k,n,texp
    !Internal
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: An,Ap,temp1,temp2,temp3,one
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: s,b,stemp,btemp
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: w,work,al,qm,xi,dm,bm
    REAL(KIND=8), DIMENSION(0:1) :: eigs 
    INTEGER(KIND=4) :: i,j
    REAL(KIND=8) :: la, tol, qnrm
    INTEGER :: lwork, info
    LOGICAL :: conv

    !Initial step
    ALLOCATE(s(0:n-1,0:m-1))
    ALLOCATE(b(0:n-1,0:m-1))
    ALLOCATE(An(0:m-1,0:n-1))
    ALLOCATE(Ap(0:m-1,0:m-1))
    ALLOCATE(qm(0:n-1))
    ALLOCATE(xi(0:n-1))
    ALLOCATE(dm(0:n-1))
    ALLOCATE(bm(0:n-1))
    ALLOCATE(one(0:n-1,0:n-1))
    ALLOCATE(temp1(0:n-1,0:n-1))
    ALLOCATE(temp2(0:n-1,0:n-1))
    ALLOCATE(temp3(0:n-1,0:n-1))

    !FIRST ITERATION STUFF
    s(0:n-1,0:m-1) = s_in(0:n-1,0:m-1)
    b(0:n-1,0:m-1) = b_in(0:n-1,0:m-1)
    An = 0
    one = 1.0D0
    temp2 = 0.0
    DO i=0,n-1
      temp2(i,i) = 1.0D0
    END DO
    conv = .FALSE.

    CALL set_tol(texp,tol) 

    WRITE(*,*) "Starting Davidson"
    WRITE(*,*) "Iteration       eigenvalue      norm" 

    !form An
    DO j=0,m-1
      DO i=0,j
        An(i,j) = linal_vTv_1Dreal8(b(0:n-1,i),s(0:n-1,j),n)
      END DO
    END DO
    Ap = An

    !WRITE(*,*) "An is..."
    !CALL linal_printmat_2Dreal8(An(0:m-1,0:m-1),m,m) 
   
    !diagonalize An, using LAPACK 
    ALLOCATE(w(0:m-1))
    ALLOCATE(work(0:1))
    lwork = -1
    CALL DSYEV('V','U',m,Ap,m,w,work,lwork,info)
    lwork = CEILING(MAX(2.0,work(0)))
    DEALLOCATE(work)
    ALLOCATE(work(0:lwork-1))
    CALL DSYEV('V','U',m,Ap,m,w,work,lwork,info)

    !grab kth eigenvalue and eigenvector
    ALLOCATE(al(0:m-1))
    al = Ap(0:m-1,k-1)
    la = w(k-1)

    !WRITE(*,*) "Eigenvector is..."
    !WRITE(*,*) al(0:m-1)
    !WRITE(*,*) "Eigenvalue is"
    !WRITE(*,*) la
    
    !Initial eigs, force and extra iteration
    eigs(0) = la
    eigs(1) = la + 2*tol

    !BEGIN ITERATIONS
    DO  WHILE (.NOT. conv) 

      !check we aren't over the limit
      IF (m .GT. 100) THEN
        WRITE(*,*) "Davidson failed to converge"
        EXIT
      END IF
   
      !construct qm
      qm = 0
      DO i=0,m-1
        qm = qm + al(i)*(s(0:n-1,i) - la*b(0:n-1,i)) 
      END DO 
      qnrm = norm2(qm)

      !WRITE(*,*) "qm is:"
      !WRITE(*,*) qm(0:n-1)
  
      !check convergence
      WRITE(*,*) m,"  ", la,"  ",qnrm
      IF (ABS(eigs(0) - eigs(1)) .LT. tol .OR. qnrm .LT. tol) THEN
        WRITE(*,*) "Davidson has converged!"     
        conv = .TRUE.
      END IF
    
      !form xi
      xi = 0.0
      DO i=0,n-1
        xi(i) = qm(i)/(la - A(i,i))
      END DO

      !WRITE(*,*) "xi is"
      !WRITE(*,*) xi(0:n-1)

      !form dm, this is updated each time
      temp1 = 0.0D0
      i = 0
      CALL linal_vvT_1Dreal8(b(0:n-1,i),b(0:n-1,i),n,temp1)
      temp2 = MATMUL(one(0:n-1,0:n-1) - temp1(0:n-1,0:n-1),temp2(0:n-1,0:n-1))  
      CALL linal_smv_2Dreal8(temp2(0:n-1,0:n-1),xi(0:n-1),n)
      dm(0:n-1) = xi(0:n-1)
      !WRITE(*,*) "dm is:"
      !WRITE(*,*) dm(0:n-1)

      !form bm
      bm = dm/norm2(dm)
      !WRITE(*,*) "bm is:"
      !WRITE(*,*) bm(0:n-1)

      !add new vector to subspace 
      m = m+1
      !WRITE(*,*) "B before is..."
      !CALL linal_printmat_2Dreal8(b(0:n-1,0:m-1),n,m)  
      ALLOCATE(btemp(0:n-1,0:m-1))   
      btemp(0:n-1,0:m-2) = b(0:n-1,0:m-2)
      btemp(0:n-1,m-1) = bm(0:n-1)
      DEALLOCATE(b)
      ALLOCATE(b(0:n-1,0:m-1))
      b = btemp
      DEALLOCATE(btemp)
      !WRITE(*,*) "B after is"
      !CALL linal_printmat_2Dreal8(b(0:n-1,0:m-1),n,m)  

      !generate new sigma vectors
      !WRITE(*,*)
      !WRITE(*,*) "S before is..."
      !CALL linal_printmat_2Dreal8(s(0:n-1,0:m-1),n,m)  
      ALLOCATE(stemp(0:n-1,0:m-1))
      stemp(0:n-1,0:m-2) = s(0:n-1,0:m-2)
      CALL linal_smv_2Dreal8(A(0:n-1,0:n-1),bm(0:n-1),n)
      stemp(0:n-1,m-1) = bm(0:n-1)
      DEALLOCATE(s)
      ALLOCATE(s(0:n-1,0:m-1))
      s(0:n-1,0:m-1) = stemp(0:n-1,0:m-1)
      DEALLOCATE(stemp)
      !WRITE(*,*) "S after is..."
      !CALL linal_printmat_2Dreal8(s(0:n-1,0:m-1),n,m)  
    
      !generate new A~ matrix
      !WRITE(*,*) 
      !WRITE(*,*) "Old A~"
      !CALL linal_printmat_2Dreal8(An(0:m-2,0:m-2),m-1,m-1)
      DEALLOCATE(Ap)
      ALLOCATE(Ap(0:m-1,0:m-1))
      Ap = 0.0D0
      Ap(0:m-2,0:m-2) = An(0:m-2,0:m-2)
      DEALLOCATE(An)
      DO i=0,m-1
        Ap(i,m-1) = linal_vTv_1Dreal8(b(0:n-1,i),s(0:n-1,m-1),n)
      END DO
      ALLOCATE(An(0:m-1,0:m-1))
      An = Ap
      !WRITE(*,*) "New A~"
      !CALL linal_printmat_2Dreal8(An(0:m-1,0:m-1),m,m)

      !diagonalize An, using LAPACK 
      DEALLOCATE(w)
      DEALLOCATE(work)
      ALLOCATE(w(0:m-1))
      ALLOCATE(work(0:1))
      lwork = -1
      CALL DSYEV('V','U',m,Ap,m,w,work,lwork,info)
      lwork = CEILING(MAX(2.0,work(0)))
      DEALLOCATE(work)
      ALLOCATE(work(0:lwork-1))
      CALL DSYEV('V','U',m,Ap,m,w,work,lwork,info)

      !grab kth eigenvalue and eigenvector
      DEALLOCATE(al)
      ALLOCATE(al(0:m-1))
      al = Ap(0:m-1,k-1)
      la = w(k-1)

      !WRITE(*,*) "new eigenvalue:"
      !WRITE(*,*) la
      !WRITE(*,*) "new eigenvector:"
      !WRITE(*,*) al(0:m-1)
      !WRITE(*,*) "----------------------"

    END DO

    !CLEANUP
    WRITE(*,*) 

    DEALLOCATE(al)
    DEALLOCATE(w)
    DEALLOCATE(work)
    DEALLOCATE(bm)
    DEALLOCATE(dm)
    DEALLOCATE(xi)
    DEALLOCATE(qm)
    DEALLOCATE(temp1)
    DEALLOCATE(temp2)
    DEALLOCATE(temp3)
    DEALLOCATE(An)
    DEALLOCATE(Ap)
    DEALLOCATE(s)
    DEALLOCATE(b)

  END SUBROUTINE linal_davidson_2Dreal8

!---------------------------------------------------------------------
!       set_tol
!               James H. Thorpe
!               Nov 19, 2018
!       - sets the tolerance of variable tol
!---------------------------------------------------------------------
  SUBROUTINE set_tol(texp,tol)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(INOUT) :: tol
    INTEGER, INTENT(IN) :: texp
    IF (texp .EQ. 1) THEN
      tol = 1.0D-1
    ELSE IF (texp .EQ. 2) THEN
      tol = 1.0D-2
    ELSE IF (texp .EQ. 3) THEN
      tol = 1.0D-3
    ELSE IF (texp .EQ. 4) THEN
      tol = 1.0D-4
    ELSE IF (texp .EQ. 5) THEN
      tol = 1.0D-5
    ELSE IF (texp .EQ. 6) THEN
      tol = 1.0D-6
    ELSE IF (texp .EQ. 7) THEN
      tol = 1.0D-7  
    ELSE IF (texp .EQ. 8) THEN
      tol = 1.0D-8  
    ELSE IF (texp .EQ. 9) THEN
      tol = 1.0D-9
    ELSE IF (texp .EQ. 10) THEN
      tol = 1.0D-10
    ELSE IF (texp .EQ. 11) THEN
      tol = 1.0D-11
    ELSE IF (texp .EQ. 12) THEN
      tol = 1.0D-12
    ELSE IF (texp .EQ. 13) THEN
      tol = 1.0D-13
    ELSE 
      WRITE(*,*) "Sorry, that convergence is not coded."
      WRITE(*,*) "Nor should it be, shame on you for thinking that"
      WRITE(*,*) "I coded this well enough for E-13 convergence"
    END IF
  END SUBROUTINE set_tol

!---------------------------------------------------------------------
!       linal_vvT_1Dreal8
!               James H. Thorpe
!               Nov 23, 2018
!       -Subroutine that performs v1*v2^T, to form matrix A
!---------------------------------------------------------------------
  !values
  !v1           : 1D real8, vector1 [N]
  !v2           : 1D real8, vector2 [N]
  !N            : 1D real8, size of v1,2
  !A            : 2D real8, array created, [NxN]
  SUBROUTINE linal_vvT_1Dreal8(v1,v2,N,A)
    IMPLICIT NONE
    !Inout
    REAL(KIND=8),DIMENSION(0:,0:),INTENT(INOUT) :: A
    REAL(KIND=8),DIMENSION(0:),INTENT(IN) :: v1,v2
    INTEGER(KIND=4),INTENT(IN) :: N
    !Internal
    INTEGER(KIND=4) :: i
    DO i=0,N-1
      A(0:N-1,i) = v1(0:N-1)*v2(i)
    END DO
  END SUBROUTINE linal_vvT_1Dreal8

!---------------------------------------------------------------------
!       linal_smv_2Dreal8
!               James H. Thorpe
!               Nov 23, 2018
!       -Subroutine that calculates A*v, where A is an [N,N] matrix,
!         b [N] vector
!       -results are stored in the v vector
!---------------------------------------------------------------------
  !Values
  !A            : 2Dreal8, square matrix to multiply, [N,N]
  !v            : 1Dreal8, vector to multiply [N]
  !N            : int, size of matrix
  SUBROUTINE linal_smv_2Dreal8(A,v,N)
    IMPLICIT NONE
    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: v
    INTEGER(KIND=4), INTENT(IN) :: N
    !Internal
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: temp
    INTEGER(KIND=4) :: i
    ALLOCATE(temp(0:N-1))
    temp = 0.0D0
    DO i=0,N-1
      temp(i) = linal_vTv_1Dreal8(A(i,0:N-1),v(0:N-1),N) 
    END DO
    v(0:N-1) = temp(0:N-1)
    DEALLOCATE(temp)
  END SUBROUTINE linal_smv_2Dreal8

!---------------------------------------------------------------------
!       linal_vTv_1Dreal8
!               James H. Thorpe
!               Nov 23, 2018
!       -Function that caluclates v1^T*v2
!---------------------------------------------------------------------
  !Values
  !v1           : 1D real8, vector1 [N]
  !v2           : 1D real8, vector2 [N]
  !N            : int, size of v1,2 
  REAL(KIND=8) FUNCTION linal_vTv_1Dreal8(v1,v2,N)
    IMPLICIT NONE
    !Inout
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: v1,v2
    INTEGER(KIND=4), INTENT(IN) :: N
    !Internal
    REAL(KIND=8) :: temp
    INTEGER(KIND=4) :: i
    temp = 0.0D0
    DO i=0,N-1
      temp = temp + v1(i)*v2(i)
    END DO
    linal_vTv_1Dreal8 = temp
  END FUNCTION linal_vTv_1Dreal8

!---------------------------------------------------------------------
!       linal_printmat_2Dreal8
!               James H. Thorpe
!               Nov 25, 2018
!       -prints matrix in row major form
!---------------------------------------------------------------------
  SUBROUTINE linal_printmat_2Dreal8(A,N,M)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: A
    INTEGER(KIND=4), INTENT(IN) :: N,M
    INTEGER(KIND=4) :: i 
    DO i=0,N-1
      WRITE(*,*) A(i,0:M-1)
    END DO
  END SUBROUTINE 
!---------------------------------------------------------------------
END MODULE linal
