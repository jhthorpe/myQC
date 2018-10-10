program R2B
  implicit none

  REAL(KIND=8), DIMENSION(0:120,0:22) :: Ftab

  OPEN(unit=1,file='old-Ftab',status='old',access='sequential')
  READ(1,*) Ftab
  CLOSE(unit=1)

  OPEN(unit=2,file='Ftab',status='replace',access='sequential',form='unformatted')
  WRITE(2) Ftab
  CLOSE(unit=2)

end program R2B
