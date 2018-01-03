program test
  implicit none
  
  integer, dimension (0:2) :: ll

  ll = [0,1,2]

  WRITE(*,*) ll
  WRITE(*,*) ll+2

end program test
