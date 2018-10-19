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
  INTEGER, DIMENSION(0:0) :: mem_lvl
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
      IF (mem_lvl(0) .EQ. 3) THEN
        CALL uvld_ijab_AAAA_high(noccA,nvrtA,ntot,options)
      ELSE
        WRITE(*,*) "Sorry, that memory case not coded yet"
        CALL EXECUTE_COMMAND_LINE('touch error')
        STOP "Bad mem case in ao2mo"
      END IF
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

    REAL(KIND=8) :: hmem,mmem,lmem

    WRITE(*,*) "--------------------------------------------------------------"
    WRITE(*,*) "Starting memory analysis"
    WRITE(*,*) "Your free memory (MB)               :", fmem 

    !MP2 analysis
    IF (options(1) .EQ. 1) THEN
      !RHF
      IF (options(3) .EQ. 0) THEN
        lmem = (3*ntot**2+nvrtA*ntot+nvrtA*nvrtA)*8/1.D6
        mmem = (2*ntot**2+nvrtA*ntot+nvrtA*nvrtA+(noccA-1)*ntot**2)*8/1.D6
        hmem = mmem + (ntot**4)*8/1.D6 
        WRITE(*,*) "<ij|ab> minimum memory (MB)         :", lmem 
        WRITE(*,*) "<ij|ab> medium memory (MB)          :", mmem
        WRITE(*,*) "<ij|ab> maximum memory (MB)         :", hmem 
        IF (hmem .LT. fmem-1) THEN
          WRITE(*,*) "K matrix                            : held in memory"
          WRITE(*,*) "L matrix                            : held in memory"
          mem_lvl(0) = 3
        ELSE IF (mmem .LT. fmem-1) THEN
          WRITE(*,*) "K matrix                            : R/W to disk"
          WRITE(*,*) "L matrix                            : held in memory"
          mem_lvl(0) = 2
        ELSE
          WRITE(*,*) "K matrix                            : R/W to disk"
          WRITE(*,*) "L matrix                            : R/W to disk"
          mem_lvl(0) = 1
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
  ! Km          :       4D real8, K matrix of (l,d,u,v). ntot x ntot x ntot x ntot
  ! Lm          :       3D real8, L matrix of (d,v,l). ntot x ntot x ntot
  ! Cm          :       2D real8, Cui matrix of coefficients ntot x ntot
  ! Mm          :       2D real8, M matrix of ij(v,d). ntot x ntot
  ! Nm          :       2D real8, N matrix of ij(a,d). nvrt x ntot
  ! Om          :       2D real8, M matrix of ij(a,b). nvrt x nvrt
  ! f[1,4-6]	:	1D real8, dummy matrix
  ! f[2-3]	:	1D real8, dummy vector

  SUBROUTINE uvld_ijab_AAAA_high(noccA,nvrtA,ntot,options)
    IMPLICIT NONE

    !Inout
    INTEGER, DIMENSION(0:), INTENT(IN) :: options
    INTEGER, INTENT(IN) :: noccA,nvrtA,ntot
    
    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: Km
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: Lm
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Cm,Mm,Nm,Om,f1,f4,f5
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: f2,f3
    INTEGER :: u,v,l,d,i,j,a,b

    WRITE(*,*) 
    WRITE(*,*) "Number of AO itegrals (uv|ld)   : ", CEILING(((noccA+1)*noccA)**2/8.0+((noccA+1)*noccA)/4.0) 
    WRITE(*,*) "Number of MO integrals <ij|ab>  : ", CEILING((noccA-1)*noccA/2.0*((nvrtA-1)*nvrtA/2.0)) 
    WRITE(*,*) 
    WRITE(*,*) "Spin Case <AA|AA>" 

    !Sort the 2e- integrals into Ild,u,v order
    CALL build_K(ntot,3)
    ALLOCATE(Km(0:ntot-1,0:ntot-1,0:ntot-1,0:ntot-1))
    ALLOCATE(Lm(0:ntot-1,0:ntot-1,0:ntot-1))
    ALLOCATE(Om(0:nvrtA-1,0:nvrtA-1))
    ALLOCATE(Nm(0:nvrtA-1,0:ntot-1))
    ALLOCATE(Mm(0:ntot-1,0:ntot-1))
    ALLOCATE(Cm(0:ntot-1,0:ntot-1))
    ALLOCATE(f1(0:ntot-1,0:ntot-1))
    ALLOCATE(f4(0:ntot-1,0:nvrtA-1))
    ALLOCATE(f5(0:nvrtA-1,0:ntot-1))
    ALLOCATE(f2(0:ntot-1))
    ALLOCATE(f3(0:ntot-1))

    Km = 0
    Lm = 0
    Mm = 0
    Cm = 0

    !get K matrix
    OPEN(unit=100,file='K_h',status='old',access='sequential',form='unformatted')
    READ(100) Km(:,:,:,:)
    CLOSE(unit=100)

    !get coef matrix
    OPEN(unit=101,file='Cui',status='old',access='sequential')
    READ(101,*) Cm(:,:)
    CLOSE(unit=101)

    WRITE(*,*) "TESTING TESTING"
    WRITE(*,*) "noccA", noccA
    WRITE(*,*) "Cm(:,1)"
    WRITE(*,*) Cm(:,1)

    OPEN(unit=102,file='ijab',status='replace',access='sequential',form='unformatted')
    OPEN(unit=103,file='ijba',status='replace',access='sequential',form='unformatted')

    DO i=0,noccA-2
      DO d=0,ntot-1
        DO l=0,ntot-1
          !construct L(d,:,l) 
          IF (d .EQ. -1) THEN
            WRITE(*,*) "TESTING TESTING TESTING"
            WRITE(*,*) "Km is:"
            WRITE(*,*) Km(:,:,l,d) 
            WRITE(*,*)
            WRITE(*,*) "C(:,i) is"
            WRITE(*,*) Cm(:,i)
          END IF
         
          f1 = Km(:,:,l,d)
          f2 = Cm(:,i) 
          f3 = Lm(d,:,l)
          CALL DSYMV('U',ntot,1.0D0,f1,ntot,f2,1,0.0D0,f3,1)
          Lm(d,:,l) = f3(:)

          IF (d .EQ. -1) THEN
            WRITE(*,*) 
            WRITE(*,*) "Li(d,:,:) is..."
            WRITE(*,*) Lm(d,:,:)
            WRITE(*,*) "TESTING TESTING TESTING"
          END IF

        END DO !lambda loop
      END DO !delta loop
      DO j=i+1,noccA-1
        DO d=0,ntot-1
          !construct M(v,d)
          f1 = Lm(d,:,:)
          f2 = Cm(:,j)
          f3 = Mm(:,d)
!          CALL DGEMV('N',ntot,ntot,1.0D0,Lm(d,:,:),ntot,Cm(:,j),1,0.0D0,Mm(:,d),1)
          CALL DGEMV('N',ntot,ntot,1.0D0,f1,ntot,f2,1,0.0D0,f3,1)
          Mm(:,d) = f3(:)

          IF (d .EQ. -1) THEN
            WRITE(*,*) "j,d=",j,d
            WRITE(*,*) "Lm(d,:,:) is:"
            WRITE(*,*) Lm(d,:,:)
            WRITE(*,*) ""
            WRITE(*,*) "Cm(:,j) is"
            WRITE(*,*) Cm(:,j)
            WRITE(*,*) ""
            WRITE(*,*) "Mij(:,:) is"
            WRITE(*,*) Mm(:,:)
          END IF
        END DO !delta loop
      DO a=noccA+1,ntot-2
        !construct N(a,d)
        f1 = Mm(:,:) 
        f4 = Cm(0:ntot-1,noccA+1:ntot-1)
        !check that this behaves properly
        CALL DGEMM('T','N',nvrtA,ntot,ntot,1.0D0,f4,ntot,f1,ntot,0.0D0,Nm,nvrtA)
        DO b=a+1,ntot-1
          !construct O(a,b) 
          f5 = Nm
          f4 = Cm(0:ntot-1,noccA+1:ntot-1)
          CALL DGEMM('N','N',nvrtA,nvrtA,ntot,1.0D0,f5,nvrtA,f4,ntot,0.0,Om,nvrtA)
          !Write integrals to disk
          CALL write_Oij(Om,nvrtA)
        END DO !b loop
      END DO !a loop
      END DO !j loop
    END DO !i loop

    CLOSE(unit=102,status='keep') 
    CLOSE(unit=103,status='keep')

    DEALLOCATE(Km)
    DEALLOCATE(Lm)
    DEALLOCATE(Om)
    DEALLOCATE(Nm)
    DEALLOCATE(Mm)
    DEALLOCATE(Cm)
    DEALLOCATE(f1)
    DEALLOCATE(f2)
    DEALLOCATE(f3)
    DEALLOCATE(f4)
    DEALLOCATE(f5)
    
  END SUBROUTINE uvld_ijab_AAAA_high

!---------------------------------------------------------------------
!       build_K
!               James H. Thorpe
!               Oct 11, 2018
!       - read in (uv|ld) ao integrals and sort into K matrices
!       - NOTE : this algorithm is based on my crappy 4 index storage
!       - it WILL need to be fixed when I fix that
!---------------------------------------------------------------------
  ! Variables
  ! ntot        :       int, number of orbitals
  ! mem_lvl     :       int, memory level
  ! K_h         :       2d real8, K matrix storage in high memory case

  SUBROUTINE build_K(ntot,mem_lvl)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ntot, mem_lvl
  
    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: XX, K_h
    INTEGER :: ios1,ios2
    INTEGER :: Ild,l,d

    ALLOCATE(XX(0:ntot-1,0:ntot-1,0:ntot-1,0:ntot-1),STAT=ios1)
    IF (ios1 .NE. 0) THEN
      WRITE(*,*) "ao2mo:build_K failed to allocate memory"
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP "memory problem in ao2mo:build_K"
    END IF

    OPEN(unit=100,file='XX',status='old',access='sequential',form='unformatted')
    READ(100) XX(:,:,:,:)
    CLOSE(unit=100)

    IF (mem_lvl .EQ. 3) THEN

      !high memory case
      ALLOCATE(K_h(0:ntot-1,0:ntot-1,0:ntot-1,0:ntot-1),STAT=ios2)
      IF (ios2 .NE. 0) THEN
        WRITE(*,*) "ao2mo:build_K failed to allocate memory"
        CALL EXECUTE_COMMAND_LINE('touch error')
        STOP "memory problem in ao2mo:build_K"
      END IF

      DO l=0,ntot-1
        DO d=l,ntot-1 !using symmetry
          K_h(l,d,:,:) = XX(:,:,l,d)
        END DO
      END DO  

      OPEN(unit=101,file='K_h',status='replace',access='sequential',form='unformatted')
      WRITE(101) K_h(:,:,:,:)
      CLOSE(unit=101,status='keep')

    ELSE
      !low/medium memory case
      IF (ios2 .NE. 0) THEN
        WRITE(*,*) "ao2mo:build_K failed to allocate memory"
        CALL EXECUTE_COMMAND_LINE('touch error')
        STOP "memory problem in ao2mo:build_K"
      END IF

      OPEN(unit=101,file='K_l',status='replace',access='sequential',form='unformatted')
      DO l=0,ntot-1
        DO d=0,ntot-1 !not using symmetry, we don't care about disk space
          WRITE(101) XX(:,:,l,d) 
        END DO 
      END DO
      CLOSE(unit=101,status='keep')
    END IF
   
    DEALLOCATE(XX)
    IF (ALLOCATED(K_h)) DEALLOCATE(K_h)

  END SUBROUTINE build_K

!---------------------------------------------------------------------
!	write_Oij
!		James H. Thorpe
!		Oct 19, 2018
!	- write Oij matrix of <ij|ab> and <ij|ba> integrals to disk	
!	- assumes that we are writing to files 102 and 103
!---------------------------------------------------------------------
  ! Variables
  ! Oij		:	2d real8, matrix of <ij|ab> integrals
  ! n		:	int, dimension of Oij matrix

  SUBROUTINE write_Oij(Oij,n)
    IMPLICIT NONE
    !inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Oij
    INTEGER, INTENT(IN) :: n
    !internal
    INTEGER :: a,b

    WRITE(*,*)
    WRITE(*,*) 
    WRITE(*,*) "Oij is:"
    WRITE(*,*) Oij

    DO a=0,n-2 
      DO b=a+1,n-1 
        WRITE(102) Oij(a,b) 
        WRITE(103) Oij(b,a)
      END DO
    END DO

  END SUBROUTINE write_Oij
!---------------------------------------------------------------------

END PROGRAM ao2mo

