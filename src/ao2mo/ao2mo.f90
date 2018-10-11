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
  ! noccA       :       int, number of electrons in alpha orbitals
  ! noccB       :       int, number of electrons in beta orbitals
  ! options     :       1D int, options array
  ! ntot        :       int, total number of orbitals
  ! 

PROGRAM ao2mo
  USE env
  IMPLICIT NONE

  !Internal
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms, options
  INTEGER, DIMENSION(0:1) :: line
  INTEGER, DIMENSION(0) :: mem_lvl
  REAL(KIND=8) :: fmem, timeS, timeF
  INTEGER :: nnuc,noccA,noccB,nvrtA,nvrtB,ntot,dummy
  LOGICAL :: flag

  CALL CPU_TIME(timeS)
  WRITE(*,*) 
  WRITE(*,*) "ao2mo called"
    WRITE(*,*)
  WRITE(*,*) "Starting AO to MO integral transform"

  !Read enviromental data
  CALL getenv(nnuc,noccA,noccB,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  !Determine memory needs
  !get orbital data
  OPEN(unit=100,file='basinfo',status='old',access='sequential')
  READ(100,*) line
  ntot = line(1)
  CLOSE(unit=1)
  nvrtA = ntot - noccA
  nvrtB = ntot - noccB
  
  CALL mem_analysis(mem_lvl,noccA,noccB,nvrtA,nvrtB,ntot,fmem,options)  

  IF (options(1) .EQ. 1) THEN
    IF (options(3) .EQ. 0) THEN
!      CALL uvld_ijab_AAAA(noccA,options)
    ELSE
      WRITE(*,*) "Sorry, that spin case for MP2 has not been coded"      
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP "Bad ref in ao2mo"
    END IF
  ELSE
    WRITE(*,*) "Sorry, that transform type has not been coded yet"
    CALL EXECUTE_COMMAND_LINE('touch error')
    STOP "Bad calc in ao2mo"
  END IF
  
  CALL CPU_TIME(timeF)
  WRITE(*,*) "ao2mo completed in (s): ", timeF-timeS

  CONTAINS
!---------------------------------------------------------------------
!       mem_analysis
!               James H. Thorpe 
!               Oct 11, 2018
!       - analyis the memory requirements of various ao2mo integral
!         transformations
!       - outputs: mem_lvl [ <ij|ab> ]
!       - mem_levels:   1 - lowmem, 2 - medmem, 3 - highmem 
!---------------------------------------------------------------------
  ! Variables
  ! mem_lvl     :       1D int, list of memory levels
  ! noccA,B     :       int, number of alpha,beta occupied orbitals
  ! nvrtA,B     :       int, number of alpha,beta virtual orbitals
  ! ntot        :       int, number of total orbitals
  ! fmem        :       real8, free memory in MB
  ! options     :       1D int, list of options
  ! hmem        :       real8, high memory requirement
  ! lmem        :       real8, low memory requirement

  SUBROUTINE mem_analysis(mem_lvl,noccA,noccB,nvrtA,nvrtB,ntot,fmem,options)
    IMPLICIT NONE
  
    INTEGER, DIMENSION(0:), INTENT(INOUT) :: mem_lvl
    INTEGER, DIMENSION(0:), INTENT(IN) :: options
    REAL(KIND=8), INTENT(IN) :: fmem
    INTEGER, INTENT(IN) :: noccA,noccB,nvrtA,nvrtB,ntot

    REAL(KIND=8) :: hmem, lmem

    WRITE(*,*) "--------------------------------------------------------------"
    WRITE(*,*) "Starting memory analysis"
    WRITE(*,*) "Your free memory (MB)               :", fmem 
    

    !MP2 analysis
    IF (options(1) .EQ. 1) THEN
      !RHF
      IF (options(3) .EQ. 0) THEN
        lmem = (3*ntot**2+nvrtA*ntot+nvrtA*nvrtA)*8/1.D6
        hmem = (2*ntot**2+nvrtA*ntot+nvrtA*nvrtA+(noccA-1)*ntot**2)*8/1.D6
        WRITE(*,*) "<ij|ab> miminum memory (MB)         :", lmem 
        WRITE(*,*) "<ij|ab> high memory (MB)            :", hmem 
        IF (hmem .LT. fmem-1) THEN
          WRITE(*,*) "L matrix                            : held in memory"
          mem_lvl(0) = 3
        ELSE
          WRITE(*,*) "L matrix                            : R/W to disk"
          mem_lvl(1) = 1
        END IF
      !UHF
      ELSE
        WRITE(*,*) "Sorry, that spin case not coded yet" 
        CALL EXECUTE_COMMAND_LINE('touch error')
        STOP "Bad spin case in ao2mo"
      END IF   
    !other analysis
    ELSE
      WRITE(*,*) "Sorry, beyond MP2 not coded yet" 
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP "Bad calc in ao2mo"
    END IF

    WRITE(*,*) "--------------------------------------------------------------"

  END SUBROUTINE mem_analysis

!---------------------------------------------------------------------
  
!---------------------------------------------------------------------
!       uvld_ijab_AAAA_high
!               James H. Thorpe
!               Oct 10, 2018
!       - subroutine that performs the partial integral transform
!         from AO (uv|ld) to the partial MO integrals <ij|ab>
!         and <ij|ba>
!       - only for (AA|AA) and (BB|BB) spin cases 
!       - high memory case, hold everything in core memory
!---------------------------------------------------------------------
  !Variables
  ! noccA       :       int, number of occupied orbitals of spin alpha
  ! nvrtA       :       int, number of virtual orbitals of spin alpha
  ! ntot        :       int, number of total orbitals
  ! options     :       1D int, options array
  ! line        :       1D int, dummy readline

  SUBROUTINE uvld_ijab_AAAA_high(noccA,nvrtA,ntot)
    IMPLICIT NONE

    !Inout
    INTEGER, INTENT(IN) :: noccA,nvrtA,ntot
    
    !Internal
    INTEGER :: u,v,l,d,i,j,a,b

    WRITE(*,*) 
    WRITE(*,*) "Number of AO itegrals (uv|ld)   : ", CEILING(((noccA+1)*noccA)**2/8.0+((noccA+1)*noccA)/4.0) 
    WRITE(*,*) "Number of MO integrals <ij|ab>  : ", CEILING((noccA-1)*noccA/2.0*((nvrtA-1)*nvrtA/2.0)) 
    WRITE(*,*) 
    WRITE(*,*) "Spin Case <AA|AA>" 

  END SUBROUTINE uvld_ijab_AAAA_high

END PROGRAM ao2mo

