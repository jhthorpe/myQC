program test
  implicit none
  
  integer :: i,j,k
  real(kind=8), dimension(0:3) :: set, coef 
  integer, dimension(0:3) :: setinfo
  integer, dimension(0:1) :: orb

  set = [1.0,2.0,3.0,4.0]
  coef = [0.25, 0.33, 0.42, 0.56]
  setinfo = [0,0,1,1]
  orb = [0,1]

  do i=0,3
    do j=0,3
      write(*,*)i,j, orb(setinfo(i)),orb(setinfo(j)), coef(i), coef(j)

    end do
  end do
end program test
