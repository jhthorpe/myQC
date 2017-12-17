program test
  use basis
  use env
  implicit none

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: bas
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz,S,F,MOc,set
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: basinfo, setinfo
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms,options
  REAL(KIND=8) :: timeS, timeF, fmem
  INTEGER :: nnuc,nelc,i,j,k,norb,npri,stat,setlen
  LOGICAL :: flag1,flag2,flag

  CALL getenv(nnuc,nelc,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  CALL buildBasis(options(2),atoms,bas,basinfo,set,setinfo)
  DO i=0,nnuc-1
    DO j=0,setinfo(i,0)-1
      WRITE(*,*) "set:", j
      WRITE(*,*) "alpha: ", set(i,j)
      WRITE(*,*) "orbs in set"
      setlen = setinfo(i,2)
      DO k=1,setinfo(i,3+j*setlen)
        WRITE(*,*) setinfo(i,4+k + j*setlen) 
      END DO
      WRITE(*,*) "setinfo: ",setinfo(i,2+j*setlen+1:2+(j+1)*setlen)
      WRITE(*,*) 
    END DO 
    WRITE(*,*) "~~~~~~~~~~"
    WRITE(*,*) "~~~~~~~~~~"
  END DO


  CALL setenv(atoms,xyz,fmem,options)

end program test
