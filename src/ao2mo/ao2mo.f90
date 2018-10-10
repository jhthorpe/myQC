!//////////////////////////////////////////////////////////////////
!//           Performs atomic to molecular integral transformations 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//             
!///////////////////////////////////////////////////////////////////

!---------------------------------------------------------------------
!       ao2mo
!               James H. Thorpe
!               Oct 10, 2018
!       - reads in jobdata to determine which integral transforms...
!               ...need to be done
!---------------------------------------------------------------------
  !  Variables
  ! nelcA       :       int, number of electrons in alpha orbitals
  ! nelcB       :       int, number of electrons in beta orbitals
  ! options     :       1D int, options array

PROGRAM ao2mo
  USE env
  IMPLICIT NONE

  !Internal
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms, options
  REAL(KIND=8) :: fmem, timeS, timeF
  INTEGER :: nnuc, nelcA, nelcB, dummy
  LOGICAL :: flag

  CALL CPU_TIME(timeS)
  WRITE(*,*) 
  WRITE(*,*) "Starting AO to MO integral transform"
  WRITE(*,*) 
  CALL getenv(nnuc,nelcA,nelcB,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  IF (options(1) .EQ. 1) THEN
    IF (options(3) .EQ. 0) THEN
      CALL uvld_ijab_RHF(nelcA,options)
    ELSE
      WRITE(*,*) "Sorry, that spin case for MP2 has not been coded"      
      STOP "Bad ref in ao2mo"
    END IF
  ELSE
    WRITE(*,*) "Sorry, that transform type has not been coded yet"
    STOP "Bad calc in ao2mo"
  END IF
  
  CALL CPU_TIME(timeF)
  WRITE(*,*) "ao2mo completed in (s): ", timeF-timeS

  CONTAINS
  
!---------------------------------------------------------------------
!       uvld_ijab_RHF
!               James H. Thorpe
!               Oct 10, 2018
!       - subroutine that performs the partial integral transform
!         from AO (uv|ld) to the partial MO integrals <ij|ab>
!       - only for RHF references
!---------------------------------------------------------------------
  !Variables
  ! noccA       :       int, number of occupied orbitals of spin alpha
  ! nvrtA       :       int, number of virtual orbitals of spin alpha
  ! ntot        :       int, number of total orbitals
  ! options     :       1D int, options array
  ! line        :       1D int, dummy readline

  SUBROUTINE uvld_ijab_RHF(noccA,options)
    IMPLICIT NONE
    
    !Inout
    INTEGER, DIMENSION(0:), INTENT(IN) :: options
    INTEGER, INTENT(IN) :: noccA

    !Internal
    INTEGER, DIMENSION(0:1) :: line
    INTEGER :: nvrtA,ntot
    INTEGER :: u,v,l,d,i,j,a,b

    !get orbital data
    OPEN(unit=100,file='basinfo',status='old',access='sequential')
    READ(100,*) line
    ntot = line(1)
    CLOSE(unit=1)

    nvrtA = ntot - noccA 
    WRITE(*,*) "Number of AO itegrals (uv|ld)   : ", CEILING(((noccA+1)*noccA)**2/8.0+((noccA+1)*noccA)/4.0) 
    WRITE(*,*) "Number of MO integrals <ij|ab>  : ", CEILING((noccA+1)*noccA/2.0*((nvrtA+1)*nvrtA/2.0)) 
    WRITE(*,*) "Number of MO integrals <ij|ba>  : ", CEILING((noccA+1)*noccA/2.0*((nvrtA+1)*nvrtA/2.0)) 

  END SUBROUTINE uvld_ijab_RHF

END PROGRAM ao2mo

