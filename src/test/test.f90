program test
  use basis
  implicit none

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: bas,set
  INTEGER, DIMENSION(:), ALLOCATABLE :: basinfo,setinfo
  INTEGER,  DIMENSION(0:1) :: atoms
  INTEGER :: bkey

  atoms = [3,1]
  bkey = 0

  CALL buildBasis(bkey,atoms,bas,basinfo,set,setinfo)

  WRITE(*,*) "Set"
  WRITE(*,*) set(:)
  WRITE(*,*) 
  WRITE(*,*) "Basis"
  WRITE(*,*) bas(:)

  contains

     SUBROUTINE printBasis(nnuc,atoms,setinfo,basinfo,bas,set)
     IMPLICIT NONE
 
     ! inout
     REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: bas
     REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: set
     INTEGER, DIMENSION(0:,0:), INTENT(IN) :: basinfo, setinfo
     INTEGER, DIMENSION(0:), INTENT(IN) :: atoms
     INTEGER, INTENT(IN) :: nnuc
 
     !internal
     INTEGER :: i,j,k,orb,setlen,n,l
 
 999 FORMAT(2x,A8,1x,I2)
 998 FORMAT(4x,A8,1x,I2,4x,A8,1x,F15.8)
 997 FORMAT(3x,A6,1x,I2)
 996 FORMAT(4x,F15.8,2x,I2,2x,I2,2x,I2)
 
     WRITE(*,*)
     WRITE(*,*) "Basis set"
 
     ! go through atoms
     DO i=0,nnuc-1
       WRITE(*,999) "Nuclei #", i+1
       WRITE(*,997) "atom :", atoms(i)
 
       ! go through sets
       DO j=0,setinfo(i,0)-1
         WRITE(*,998) "set #:", j+1, "alpha :", set(i,j)
         ! go print elements of each set
         setlen = setinfo(i,2)
         DO k=0,NINT(bas(i,j,1))-1
           orb = setinfo(i,2+k+j*setlen+3)
           WRITE(*,996) bas(i,j,2+k), basinfo(i,4*(orb+1)+0),basinfo(i,4*(orb+1)+1),basinfo(i,4*(orb+1)+2)
         END DO
 
       END DO
     END DO
     WRITE(*,*) "=================================="
 
   END SUBROUTINE printBasis


end program test
