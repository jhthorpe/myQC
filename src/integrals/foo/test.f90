program test
  implicit none

  integer :: i,j,g,h,n

  n = 6

  do i=0,n-1
    do j=i,n-1
      do h=j,n-1
        write(*,*) i,j,i,h
      end do
      do g=i+1,n-1
        do h=g,n-1
          write(*,*) i,j,g,h
        end do
      end do
    end do
  end do

end program test
