program test
  use aux
  implicit none

  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: coef
  REAL(KIND=8), DIMENSION(0:2) :: PA, PB
  INTEGER, DIMENSION(:), ALLOCATABLE  :: Dk
  INTEGER, DIMENSION(0:2)  :: left, right
  REAL(KIND=8) :: aa,bb
  INTEGER :: kmax,nl,nr,ll,lr,ml,mr,k,temp1,temp2
  INTEGER :: i,j,N,L,M,Np,Lp,Mp

  nl = 1
  nr = 2
  ll = 1
  lr = 2
  ml = 1
  mr = 2

  left = [nl,ll,ml]
  right = [nr,lr,mr] 

  kmax = 0

  ALLOCATE(Dk(0:3**3))
  
  PA = [1.0D0, 5.0D0, 0.4D0]
  PB = [-0.8D0, -4.5D0, -0.6D0] 

  aa = 1.0D0
  bb = 1.0D0

  CALL getcoef(coef, PA, PB, aa, bb, [nl,ll,ml],[nr,lr,mr]) 
  coef(0,1,nl,nr) = 42.0D0
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

  kmax = kmax + 1
  Dk(kmax) = 322 
  kmax = kmax + 1
  Dk(kmax) = 646 

  ! print nonzero combinations
  DO k=0,kmax
    WRITE(*,*) "----------"
    WRITE(*,*) "k = ", k
    WRITE(*,*) "Dk = ", Dk(k) 
    temp1 = Dk(k)/300
    Np = temp1
    temp2 = Dk(k) - temp1*300
    temp1 = temp2/20
    Lp = temp1
    temp2 = temp2 - temp1*20
    Mp = temp2
    WRITE(*,*) "N,L,M :", Np, Lp, Mp 
  END DO

  DEALLOCATE(Dk)
 

end program test
