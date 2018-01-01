program test
  use aux
  implicit none

  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: coef
  REAL(KIND=8), DIMENSION(0:2) :: PA, PB
  INTEGER, DIMENSION(:), ALLOCATABLE  :: Dk
  INTEGER, DIMENSION(0:2)  :: left, right
  REAL(KIND=8) :: aa,bb
  INTEGER :: kmax,nl,nr,ll,lr,ml,mr,k
  INTEGER :: i,j,N,L,M

  nl = 1
  nr = 1
  ll = 1
  lr = 1
  ml = 1
  mr = 1

  left = [nl,ll,ml]
  right = [nr,lr,mr] 

  kmax = 0

  ALLOCATE(Dk(0:3**3))
  
  PA = [0.0D0, 0.0D0, 0.0D0]
  PB = [0.0D0, 0.0D0, 0.0D0] 

  aa = 1.0D0
  bb = 1.0D0

  CALL getcoef(coef,PA, PB, aa, bb, [nl,ll,ml],[nr,lr,mr]) 
  CALL getDk(coef,left,right,Dk,kmax)

  WRITE(*,*) "================"
  !print coefficients
  DO N=0,nl+nr
    DO L=0,ll+lr
      DO M=0,ml+mr
        WRITE(*,*) "N,L,M", N,L,M, coef(0,N,nl,nr), coef(1,L,ll,lr), coef(2,M,ml,mr)
      END DO
    END DO
  END DO
  WRITE(*,*) "================"


  ! print nonzero combinations
  DO k=0,kmax
    WRITE(*,*) "----------"
    WRITE(*,*) "k = ", k
    WRITE(*,*) "Dk = ", Dk(k) 
  END DO

  DEALLOCATE(Dk)
 

end program test
