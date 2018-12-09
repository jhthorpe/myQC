!//////////////////////////////////////////////////////////////////
!//		Performs CPHF Calculation
!//
!//		James H. Thorpe, in the Group of John Stanton
!//		The Univeristy of Florida
!//
!//////////////////////////////////////////////////////////////////

!---------------------------------------------------------------------
!	cphf	
!		James H. Thorpe
!	 	Dec 8, 2018	
!	- control program for cphf calculations 
!---------------------------------------------------------------------
  ! Variables
  ! noccA, nocc B	:	int, number of alpha,beta occupied orbitals 
  ! options		:	1D int, options array
  ! ntot		:	int, total number of orbitals
  ! nvrtA, nvrtB	:	int, number of alpha,beta virtual orbitals

PROGRAM cphf
  USE env
  IMPLICIT NONE

  !Internal
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms, options
  INTEGER, DIMENSION(0:1) :: line
  INTEGER, DIMENSION(0:0) :: mem_lvl
  REAL(KIND=8) :: fmem, timeS, timeF
  INTEGER :: nnuc,noccA,noccB,nvrtA,nvrtB,ntot,dummy
  LOGICAL :: flag
999 FORMAT(1x,A24,2x,F8.4)

  CALL CPU_TIME(timeS)
  WRITE(*,*) 
  WRITE(*,*) "cphf called"
  WRITE(*,*)

  !Read enviromental data
  CALL getenv(nnuc,noccA,noccB,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  !get orbital data
  OPEN(unit=100,file='basinfo',status='old',access='sequential')
  READ(100,*) line
  ntot = line(1)
  CLOSE(unit=1)
  nvrtA = ntot - noccA
  nvrtB = ntot - noccB

  IF (options(3) .EQ. 0) THEN !RHF
    CALL cphf_dipole_rhf(nnuc,noccA,nvrtA,ntot,xyz,atoms,fmem,options) 
  ELSE
    WRITE(*,*) "Sorry, only RHF CPHF equations are coded"
  END IF

  CALL CPU_TIME(timeF)
  WRITE(*,999) "cphf completed in (s) :", timeF-timeS


  CONTAINS
!---------------------------------------------------------------------
!	cphf_dipole_rhf
!		James H. Thorpe
!		Dec 9, 2018
!	-calculates cphf for dipole related properties, rhf
!---------------------------------------------------------------------
  ! nnuc		: int, number of nuclei
  ! noccA,nvrtA,ntot	: int, number of occupied, virtual, total orb
  ! xyz			: 2D real*8, xyz coords of nuclei
  ! atoms		: 1D int, ids (Z) of atoms 
  ! fmem		: real*8, free memory in MB
  ! options		: 1D int, list of options 
  ! Px			: 3D real*8, dipole integrals
  ! Ux,y,z		: real*8, x,y,z components of dipole moment 
  ! Duv			: 2D real*8, SCF density matrix, AO basis

  SUBROUTINE cphf_dipole_rhf(nnuc,noccA,nvrtA,ntot,xyz,atoms,fmem,options)
    IMPLICIT NONE
    !inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
    INTEGER, DIMENSION(0:), INTENT(IN) :: atoms,options
    REAL(KIND=8), INTENT(INOUT) :: fmem
    INTEGER, INTENT(IN) :: nnuc,noccA,nvrtA,ntot
    !internal
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: Px
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Duv
    REAL(KIND=8) :: Ux,Uy,Uz
    INTEGER :: i,j,k,l,a,b,c,d,p,q,r,s,u,v
999 FORMAT(F15.6,F15.6,F15.6)

    ALLOCATE(Px(0:2,0:ntot-1,0:ntot-1))
    ALLOCATE(Duv(0:ntot-1,0:ntot-1))

    OPEN(unit=100,file='Pxuv',status='old',access='sequential')
    READ(100,*) Px(0,:,:)
    READ(100,*) Px(1,:,:)
    READ(100,*) Px(2,:,:)
    CLOSE(unit=100)
    OPEN(unit=101,file='Da',status='old',access='sequential',form='unformatted')
    READ(101) Duv
    CLOSE(unit=101)

    !Dipole moment
    Ux = 0
    Uy = 0
    Uz = 0
    DO u=0,ntot-1
      DO v = 0,ntot-1
        Ux = Ux - Duv(u,v)*Px(0,v,u) 
        Uy = Uy - Duv(u,v)*Px(1,v,u) 
        Uz = Uz - Duv(u,v)*Px(2,v,u) 
      END DO
    END DO

   DO i=0,nnuc-1
     Ux = Ux + atoms(i)*xyz(i,0)
     Uy = Uy + atoms(i)*xyz(i,1) 
     Uz = Uz + atoms(i)*xyz(i,2)
   END DO

   WRITE(*,*) "                Dipole Moments                "
   WRITE(*,*) "----------------------------------------------" 
   WRITE(*,*) "       X              Y              Z        "
   WRITE(*,999) Ux,Uy,Uz
   WRITE(*,*) "----------------------------------------------" 


  END SUBROUTINE cphf_dipole_rhf
!---------------------------------------------------------------------

  
END PROGRAM cphf
