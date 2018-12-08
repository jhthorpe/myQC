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

999 FORMAT(1x,A24,2x,F8.4)


  CALL CPU_TIME(timeS)
  WRITE(*,*) 
  WRITE(*,*)"                STARTING AO TRANSFORM"
  WRITE(*,*)"------------------------------------------------------------"
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

!  CALL mem_analysis(mem_lvl,noccA,noccB,nvrtA,nvrtB,ntot,fmem,options)  

  IF (options(1) .EQ. 1) THEN !MP2
    IF (options(3) .EQ. 0) THEN !RHF
      CALL slow_ao2mo_MP2_RHF(noccA,nvrtA,ntot,options)
      !IF (mem_lvl(0) .EQ. 3) THEN !memory
      !  IF (options(12) .EQ. 0) THEN !fast
      !    !CALL uvld_ijab_AAAA_high(noccA,nvrtA,ntot,options)
      !  ELSE IF (options(12) .EQ. 1) THEN !slow
      !    CALL slow_ao2mo_RHF(noccA,nvrtA,ntot)
      !  END IF
      !ELSE
      !  WRITE(*,*) "Sorry, that memory case not coded yet"
      !  CALL EXECUTE_COMMAND_LINE('touch error')
      !  STOP "Bad mem case in ao2mo"
      !END IF
    ELSE IF (options(3) .EQ. 1) THEN !UHF
      CALL slow_ao2mo_MP2_UHF(noccA,noccB,nvrtA,nvrtB,ntot)
    ELSE !not RHF or UHF
      WRITE(*,*) "Sorry, that reference not coded yet"
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP "Bad ref in ao2mo"
    END IF
  ELSE IF (options(13) .EQ. 1 .AND. options(1) .EQ. 0) THEN !SCF CIS
    IF (options(3) .EQ. 1 ) THEN !UHF
      CALL slow_ao2mo_CIS_UHF(noccA,noccB,nvrtA,nvrtB,ntot)
    ELSE
      WRITE(*,*) "Sorry, only UHF CIS is coded" 
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP "Bad ref in ao2mo"
    END IF
  ELSE
    WRITE(*,*) "Sorry, that transform type has not been coded yet"
    CALL EXECUTE_COMMAND_LINE('touch error')
    STOP "Bad calc in ao2mo"
  END IF
  
  CALL CPU_TIME(timeF)
  WRITE(*,999) "ao2mo completed in (s): ", timeF-timeS

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
    WRITE(*,*) "Starting memory analysis - currently unused and incorrect"
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
    !Cm = TRANSPOSE(Cm)

    WRITE(*,*) "TESTING TESTING"
    WRITE(*,*) "noccA", noccA
    WRITE(*,*) "Cm(:,1)"
    WRITE(*,*) Cm(:,1)

    OPEN(unit=102,file='ijab',status='replace',access='sequential',form='unformatted')
    OPEN(unit=103,file='ijba',status='replace',access='sequential',form='unformatted')

!    DO i=0,noccA-2
    DO i=0,noccA-1
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
!      DO j=i+1,noccA-1
      DO j=0,noccA-1
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

        !construct N(a,d)
        f1 = Mm(:,:) 
        f4 = Cm(0:ntot-1,noccA+1:ntot-1)
        !check that this behaves properly
        CALL DGEMM('T','N',nvrtA,ntot,ntot,1.0D0,f4,ntot,f1,ntot,0.0D0,Nm,nvrtA)

        !construct O(a,b) 
        f5 = Nm
        f4 = Cm(0:ntot-1,noccA+1:ntot-1)
        CALL DGEMM('N','N',nvrtA,nvrtA,ntot,1.0D0,f5,nvrtA,f4,ntot,0.0,Om,nvrtA)

        !Write integrals to disk
        CALL write_Oij(Om,nvrtA)
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
        !DO d=0,ntot-1
          K_h(l,d,:,:) = XX(:,:,l,d)  !why did I do this?
          !K_h(:,:,l,d) = XX(:,:,l,d)  !why did I do this?
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

    CALL print_moints(noccA,noccB,nvrtA,nvrtB,ntot,options)
   
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

!    WRITE(*,*)
!    WRITE(*,*) 
!    WRITE(*,*) "Oij is:"
!    WRITE(*,*) Oij

!    DO a=0,n-2 
!      DO b=a+1,n-1 
    DO a=0,n-1
     DO b=0,n-1 
        WRITE(102) Oij(a,b) 
        WRITE(103) Oij(b,a)
      END DO
    END DO

  END SUBROUTINE write_Oij
!---------------------------------------------------------------------
!       slow_ao2mo
!               James H. Thorpe
!               Oct 31, 2018
!       -slow version, explicit do loops, O(N^5) ao -> mo transform
!       -useful for checking code
!---------------------------------------------------------------------
  SUBROUTINE slow_ao2mo_MP2_RHF(noccA,nvrtA,ntot,options)
    ! ij        :       int, occupied indicies
    ! ab        :       int, virtual indicies
    ! uvld      :       int, ao indicies
    ! Km        :       4D real8, matrix of (uv|ld) ao integrals
    ! Cm        :       2D real8, matrix of Cui's from SCF
    ! options	: 	1D int, array of options

    !Inout
    INTEGER, INTENT(IN) :: noccA,nvrtA,ntot

    !Internal
    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: Km,Lm,Mm,Nm,Om
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Cm
    INTEGER, DIMENSION(0:) :: options
    REAL(KIND=8) :: sum1
    INTEGER :: u,v,l,d,p,q,r,s,i,j,a,b
997 FORMAT(4x,A22)
996 FORMAT(4x,I3,1x,I3,1x,I3,1x,I3,4x,F15.8)
    !spin case AB, the only one we need now
    WRITE(*,*)
    WRITE(*,*) "Spin Case AB"
  
    ALLOCATE(Km(0:ntot-1,0:ntot-1,0:ntot-1,0:ntot-1))
    ALLOCATE(Cm(0:ntot-1,0:ntot-1))
    Km = 0
    Cm = 0

    !get K matrix
    OPEN(unit=100,file='XX',status='old',access='sequential',form='unformatted')
    READ(100) Km(:,:,:,:)
    CLOSE(unit=100)

    !get coef matrix
    OPEN(unit=101,file='Cui',status='old',access='sequential')
    READ(101,*) Cm(:,:)
    CLOSE(unit=101)
     
    !first index (uv|ld) -> (iv|ld)
    WRITE(*,*) "Transforming (uv|ld) -> <ij|ab>"
    ALLOCATE(Lm(0:noccA-1,0:ntot-1,0:ntot-1,0:ntot-1))
    CALL idx1_trans(noccA,ntot,ntot,ntot,ntot,Km,Cm(0:ntot-1,0:noccA-1),Lm)
    DEALLOCATE(Km)

    !second index (iv|ld) -> (ia|ld) 
    ALLOCATE(Mm(0:noccA-1,0:nvrtA-1,0:ntot-1,0:ntot-1))
    !  DO a=0,nvrtA-1
    !    DO l=0,ntot-1
    !      DO d=0,ntot-1
    !        sum1 = 0.0D0
    !        DO v=0,ntot-1
    !          sum1 = sum1 + (Lm(i,v,l,d) &
    !                 *Cm(v,a+noccA))
    !        END DO
    !        Mm(i,a,l,d) = sum1
    !      END DO
    !    END DO
    !  END DO 
    !END DO
    CALL idx2_trans(noccA,nvrtA,ntot,ntot,ntot,Lm,Cm(0:ntot-1,noccA:ntot-1),Mm)
    DEALLOCATE(Lm)

    !third index (ia|ld) -> (ia|jd) 
    ALLOCATE(Nm(0:noccA-1,0:nvrtA-1,0:noccA-1,0:ntot-1))
    !DO i=0,noccA-1
    !  DO a=0,nvrtA-1
    !    DO j=0,noccA-1
    !      DO d=0,ntot-1
    !        sum1 = 0.0D0
    !        DO l=0,ntot-1
    !          sum1 = sum1 + (Mm(i,a,l,d) &
    !                 *Cm(l,j))
    !        END DO
    !        Nm(i,a,j,d) = sum1
    !      END DO
    !    END DO
    !  END DO
    !END DO
    CALL idx3_trans(noccA,nvrtA,noccA,ntot,ntot,Mm,Cm(0:ntot-1,0:noccA-1),Nm)
    DEALLOCATE(Mm)

    !fourth index (ia|jd) -> (ia|jb) 
    ALLOCATE(Om(0:noccA-1,0:nvrtA-1,0:noccA-1,0:nvrtA-1))
    !DO i=0,noccA-1
    !  DO a=0,nvrtA-1
    !    DO j=0,noccA-1
    !      DO b=0,nvrtA-1
    !        sum1 = 0.0D0
    !        DO d=0,ntot-1
    !          sum1 = sum1 + (Nm(i,a,j,d) &
    !                 *Cm(d,b+noccA))
    !        END DO
    !        Om(i,a,j,b) = sum1 
    !      END DO
    !    END DO
    !  END DO
    !END DO
    CALL idx4_trans(noccA,nvrtA,noccA,nvrtA,ntot,Nm,Cm(0:ntot-1,noccA:ntot-1),Om)
    DEALLOCATE(Nm)
   
    !for now only
    WRITE(*,*) "Writing to ijab_AA" 
    OPEN(unit=106,file="ijab_AA",status="replace",form="unformatted")
    DO i=0,noccA-2
      DO j=i+1,noccA-1
        WRITE(106) Om(i,:,j,:)
      END DO
    END DO
    CLOSE(unit=106)

    WRITE(*,*) "Writing to ijab_AB" 
    OPEN(unit=107,file="ijab_AB",status="replace",form="unformatted")
    DO i=0,noccA-1
      DO j=0,noccA-1
        WRITE(107) Om(i,:,j,:)
      END DO
    END DO
    CLOSE(unit=107)

    IF (options(7) .GE. 3) THEN
    WRITE(*,*) "Two electron Integrals by Orbital (<ij|ab>)"
    WRITE(*,997) "i   j   a   b    Value"
    DO i=0,noccA-1
      DO j=i,noccA-1
        DO a=0,nvrtA-1
          DO b=a,nvrtA-1
            IF (Om(i,a,j,b) .NE. 0.0D0) WRITE(*,996) i,j,a,b,Om(i,a,j,b)
          END DO
        END DO
      END DO
    END DO

    END IF

    WRITE(*,*)
    DEALLOCATE(Om)

  END SUBROUTINE slow_ao2mo_MP2_RHF
!---------------------------------------------------------------------
!       slow_ao2mo_MP2_UHF
!               James H. Thorpe
!               Nov. 1, 2018
!       -uses do loops to create ijab_AA, ijab_BB, and ijab_AB
!---------------------------------------------------------------------
  !Variables
  ! noccA,B     :       int, number of occupied A,B orbitals
  ! nvrtA,B     :       int, number of virtual A,B orbitals
  ! ntot        :       int, total number of orbitals
  
  SUBROUTINE slow_ao2mo_MP2_UHF(noccA,noccB,nvrtA,nvrtB,ntot)
    IMPLICIT NONE
    !Inout
    INTEGER, INTENT(IN) :: noccA,noccB,nvrtA,nvrtB,ntot
    !Internal
    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: Km,Lm,Mm,Nm,Om
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: CmA,CmB
    REAL(KIND=8) :: sum1
    INTEGER :: i,j,a,b,u,v,l,d

    ALLOCATE(CmA(0:ntot-1,0:ntot-1))
    ALLOCATE(CmB(0:ntot-1,0:ntot-1))
    ALLOCATE(Km(0:ntot-1,0:ntot-1,0:ntot-1,0:ntot-1))
    Km = 0
    CmA = 0
    CmB = 0

    !get K matrix
    OPEN(unit=100,file='XX',status='old',access='sequential',form='unformatted')
    READ(100) Km(:,:,:,:)
    CLOSE(unit=100)

    !get coef matrix
    OPEN(unit=101,file='Cui',status='old',access='sequential')
    READ(101,*) CmA(:,:)
    READ(101,*) CmB(:,:)
    CLOSE(unit=101)

    !Spin Case AA
    WRITE(*,*) "Spin Case AA"
    !first index (uv|ld) -> (iv|ld)
    WRITE(*,*) "Transforming (uv|ld) -> (iv|ld)"
    ALLOCATE(Lm(0:noccA-1,0:ntot-1,0:ntot-1,0:ntot-1))
    Lm = 0.0D0
    DO i=0,noccA-1
      DO v=0,ntot-1
        DO l=0,ntot-1
          DO d=0,ntot-1
            sum1 = 0.0D0
            DO u=0,ntot-1
              sum1 = sum1 + (Km(u,v,l,d) & 
                            *CmA(u,i))
            END DO
            Lm(i,v,l,d) = sum1
          END DO
        END DO
      END DO
    END DO
    
    !second index (iv|ld) -> (ia|ld) 
    WRITE(*,*) "Transforming (iv|ld) -> (ia|ld)"
    ALLOCATE(Mm(0:noccA-1,0:nvrtA-1,0:ntot-1,0:ntot-1))
    DO i=0,noccA-1
      DO a=0,nvrtA-1
        DO l=0,ntot-1
          DO d=0,ntot-1
            sum1 = 0.0D0
            DO v=0,ntot-1
              sum1 = sum1 + (Lm(i,v,l,d) &
                     *CmA(v,a+noccA))
            END DO
            Mm(i,a,l,d) = sum1
          END DO
        END DO
      END DO 
    END DO
    DEALLOCATE(Lm)
    !third index (ia|ld) -> (ia|jd) 
    WRITE(*,*) "Transforming (ia|ld) -> (ia|jd)"
    ALLOCATE(Nm(0:noccA-1,0:nvrtA-1,0:noccA-1,0:ntot-1))
    DO i=0,noccA-1
      DO a=0,nvrtA-1
        DO j=0,noccA-1
          DO d=0,ntot-1
            sum1 = 0.0D0
            DO l=0,ntot-1
              sum1 = sum1 + (Mm(i,a,l,d) &
                     *CmA(l,j))
            END DO
            Nm(i,a,j,d) = sum1
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(Mm)
    !fourth index (ia|jd) -> (ia|jb) 
    WRITE(*,*) "Transforming (ia|jd) -> (ia|jb)"
    ALLOCATE(Om(0:noccA-1,0:nvrtA-1,0:noccA-1,0:nvrtA-1))
    DO i=0,noccA-1
      DO a=0,nvrtA-1
        DO j=0,noccA-1
          DO b=0,nvrtA-1
            sum1 = 0.0D0
            DO d=0,ntot-1
              sum1 = sum1 + (Nm(i,a,j,d) &
                     *CmA(d,b+noccA))
            END DO
            Om(i,a,j,b) = sum1 
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(Nm)
    WRITE(*,*) "Transforming (ia|jdb) -> <ij|ab>"
    WRITE(*,*) "Writing to ijab_AA" 
    OPEN(unit=106,file="ijab_AA",status="replace",form="unformatted")
    DO i=0,noccA-2
      DO j=i+1,noccA-1
        WRITE(106) Om(i,:,j,:)
      END DO
    END DO
    CLOSE(unit=106)
    DEALLOCATE(Om)

    !Spin Case BB
    WRITE(*,*) 
    WRITE(*,*) "Spin Case BB"
    !first index (uv|ld) -> (iv|ld)
    WRITE(*,*) "Transforming (uv|ld) -> (iv|ld)"
    ALLOCATE(Lm(0:noccB-1,0:ntot-1,0:ntot-1,0:ntot-1))
    Lm = 0.0D0
    DO i=0,noccB-1
      DO v=0,ntot-1
        DO l=0,ntot-1
          DO d=0,ntot-1
            sum1 = 0.0D0
            DO u=0,ntot-1
              sum1 = sum1 + (Km(u,v,l,d) & 
                            *CmB(u,i))
            END DO
            Lm(i,v,l,d) = sum1
          END DO
        END DO
      END DO
    END DO
    !second index (iv|ld) -> (ia|ld) 
    WRITE(*,*) "Transforming (iv|ld) -> (ia|ld)"
    ALLOCATE(Mm(0:noccB-1,0:nvrtB-1,0:ntot-1,0:ntot-1))
    DO i=0,noccB-1
      DO a=0,nvrtB-1
        DO l=0,ntot-1
          DO d=0,ntot-1
            sum1 = 0.0D0
            DO v=0,ntot-1
              sum1 = sum1 + (Lm(i,v,l,d) &
                     *CmB(v,a+noccB))
            END DO
            Mm(i,a,l,d) = sum1
          END DO
        END DO
      END DO 
    END DO
    DEALLOCATE(Lm)
    !third index (ia|ld) -> (ia|jd) 
    WRITE(*,*) "Transforming (ia|ld) -> (ia|jd)"
    ALLOCATE(Nm(0:noccB-1,0:nvrtB-1,0:noccB-1,0:ntot-1))
    DO i=0,noccB-1
      DO a=0,nvrtB-1
        DO j=0,noccB-1
          DO d=0,ntot-1
            sum1 = 0.0D0
            DO l=0,ntot-1
              sum1 = sum1 + (Mm(i,a,l,d) &
                     *CmB(l,j))
            END DO
            Nm(i,a,j,d) = sum1
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(Mm)
    !fourth index (ia|jd) -> (ia|jb) 
    WRITE(*,*) "Transforming (ia|jd) -> (ia|jb)"
    ALLOCATE(Om(0:noccB-1,0:nvrtB-1,0:noccB-1,0:nvrtB-1))
    DO i=0,noccB-1
      DO a=0,nvrtB-1
        DO j=0,noccB-1
          DO b=0,nvrtB-1
            sum1 = 0.0D0
            DO d=0,ntot-1
              sum1 = sum1 + (Nm(i,a,j,d) &
                     *CmB(d,b+noccB))
            END DO
            Om(i,a,j,b) = sum1 
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(Nm)
    WRITE(*,*) "Transforming (ia|jb) -> <ij|ab>"
    WRITE(*,*) "Writing to ijab_BB" 
    OPEN(unit=106,file="ijab_BB",status="replace",form="unformatted")
    DO i=0,noccB-2
      DO j=i+1,noccB-1
        WRITE(106) Om(i,:,j,:)
      END DO
    END DO
    CLOSE(unit=106)
    DEALLOCATE(Om)

    !Spin Case AB: i=A,a=A,j=B,b=B
    WRITE(*,*) 
    WRITE(*,*) "Spin Case AB"
    !first index (uv|ld) -> (iv|ld)
    WRITE(*,*) "Transforming (uv|ld) -> (iv|ld)"
    ALLOCATE(Lm(0:noccA-1,0:ntot-1,0:ntot-1,0:ntot-1))
    Lm = 0.0D0
    DO i=0,noccA-1
      DO v=0,ntot-1
        DO l=0,ntot-1
          DO d=0,ntot-1
            sum1 = 0.0D0
            DO u=0,ntot-1
              sum1 = sum1 + (Km(u,v,l,d) & 
                            *CmA(u,i))
            END DO
            Lm(i,v,l,d) = sum1
          END DO
        END DO
      END DO
    END DO
    !second index (iv|ld) -> (ia|ld) 
    WRITE(*,*) "Transforming (iv|ld) -> (ia|ld)"
    ALLOCATE(Mm(0:noccA-1,0:nvrtA-1,0:ntot-1,0:ntot-1))
    DO i=0,noccA-1
      DO a=0,nvrtA-1
        DO l=0,ntot-1
          DO d=0,ntot-1
            sum1 = 0.0D0
            DO v=0,ntot-1
              sum1 = sum1 + (Lm(i,v,l,d) &
                     *CmA(v,a+noccA))
            END DO
            Mm(i,a,l,d) = sum1
          END DO
        END DO
      END DO 
    END DO
    DEALLOCATE(Lm)
    !third index (ia|ld) -> (ia|jd) 
    WRITE(*,*) "Transforming (ia|ld) -> (ia|jd)"
    ALLOCATE(Nm(0:noccA-1,0:nvrtA-1,0:noccB-1,0:ntot-1))
    DO i=0,noccA-1
      DO a=0,nvrtA-1
        DO j=0,noccB-1
          DO d=0,ntot-1
            sum1 = 0.0D0
            DO l=0,ntot-1
              sum1 = sum1 + (Mm(i,a,l,d) &
                     *CmB(l,j))
            END DO
            Nm(i,a,j,d) = sum1
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(Mm)
    !fourth index (ia|jd) -> (ia|jb) 
    WRITE(*,*) "Transforming (ia|jd) -> (ia|jb)"
    ALLOCATE(Om(0:noccA-1,0:nvrtA-1,0:noccB-1,0:nvrtB-1))
    DO i=0,noccA-1
      DO a=0,nvrtA-1
        DO j=0,noccB-1
          DO b=0,nvrtB-1
            sum1 = 0.0D0
            DO d=0,ntot-1
              sum1 = sum1 + (Nm(i,a,j,d) &
                     *CmB(d,b+noccB))
            END DO
            Om(i,a,j,b) = sum1 
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(Nm)
    WRITE(*,*) "Transforming (ia|jdb) -> <ij|ab>"
    WRITE(*,*) "Writing to ijab_AB" 
    OPEN(unit=106,file="ijab_AB",status="replace",form="unformatted")
    DO i=0,noccA-1
      DO j=0,noccB-1
        WRITE(106) Om(i,:,j,:)
      END DO
    END DO
    CLOSE(unit=106)
    DEALLOCATE(Om)

    DEALLOCATE(Km)
    DEALLOCATE(CmA)
    DEALLOCATE(CmB) 

  END SUBROUTINE slow_ao2mo_MP2_UHF

!---------------------------------------------------------------------
!       slow_ao2mo_CIS_UHF
!               James H. Thorpe
!               Nov. 1, 2018
!       -uses do loops to create ajib_AA, ajbi_AA, ajib_AB, ajib_BB,
!          and ajbi_BB 
!	-this is terrible, aweful code of which I am not proud
!---------------------------------------------------------------------
  !Variables
  ! noccA,B     :       int, number of occupied A,B orbitals
  ! nvrtA,B     :       int, number of virtual A,B orbitals
  ! ntot        :       int, total number of orbitals
  
  SUBROUTINE slow_ao2mo_CIS_UHF(noccA,noccB,nvrtA,nvrtB,ntot)
    IMPLICIT NONE
    !Inout
    INTEGER, INTENT(IN) :: noccA,noccB,nvrtA,nvrtB,ntot
    !Internal
    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: Km,Lm,Mm,Nm,Om
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: CmA,CmB,testM
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: vec
    REAL(KIND=8) :: sum1
    INTEGER :: i,j,a,b,u,v,l,d,idx,idxA,idxB

    ALLOCATE(CmA(0:ntot-1,0:ntot-1))
    ALLOCATE(CmB(0:ntot-1,0:ntot-1))
    ALLOCATE(Km(0:ntot-1,0:ntot-1,0:ntot-1,0:ntot-1))
    Km = 0
    CmA = 0
    CmB = 0

    !get K matrix
    OPEN(unit=100,file='XX',status='old',access='sequential',form='unformatted')
    READ(100) Km(:,:,:,:)
    CLOSE(unit=100)

    !get coef matrix
    OPEN(unit=101,file='Cui',status='old',access='sequential')
    READ(101,*) CmA(0:ntot-1,0:ntot-1)
    READ(101,*) CmB(0:ntot-1,0:ntot-1)
    CLOSE(unit=101)

    !------------------------------
    !Spin Case AA : (uv|ld) -> (ai|jb) -> <aj|ib>
    WRITE(*,*) "Spin Case AA"
    WRITE(*,*) "Transforming (uv|ld) -> <aj|ib>"

    !first index (uv|ld) -> (av|ld)
    ALLOCATE(Lm(0:nvrtA-1,0:ntot-1,0:ntot-1,0:ntot-1))
    CALL idx1_trans(nvrtA,ntot,ntot,ntot,ntot,Km,CmA(0:ntot-1,noccA:ntot-1),Lm)

    !second index (av|ld) -> (ai|ld) 
    ALLOCATE(Mm(0:nvrtA-1,0:noccA-1,0:ntot-1,0:ntot-1))
    CALL idx2_trans(nvrtA,noccA,ntot,ntot,ntot,Lm,CmA(0:ntot-1,0:noccA-1),Mm)
    DEALLOCATE(Lm)

    !third index (ai|ld) -> (ai|jd) 
    ALLOCATE(Nm(0:nvrtA-1,0:noccA-1,0:noccA-1,0:ntot-1))
    CALL idx3_trans(nvrtA,noccA,noccA,ntot,ntot,Mm,CmA(0:ntot-1,0:noccA-1),Nm)
    DEALLOCATE(Mm)

    !fourth index (ai|jd) -> (ai|jb) 
    ALLOCATE(Om(0:nvrtA-1,0:noccA-1,0:noccA-1,0:nvrtA-1))
    CALL idx4_trans(nvrtA,noccA,noccA,nvrtA,ntot,Nm,CmA(0:ntot-1,noccA:ntot-1),Om)
    DEALLOCATE(Nm)
    
    !We are writing vector by vector
    WRITE(*,*) "Writing to ajib_AA" 
    ALLOCATE(vec(0:noccA*nvrtA-1))
    OPEN(unit=106,file="ajib_AA",status="replace",form="unformatted")
    DO j=0,noccA-1
      DO b=0,nvrtA-1
        vec = 0.0D0
        idx = 0
        DO i=0,noccA-1
          DO a=0,nvrtA-1
            vec(idx) = Om(a,i,j,b)
            idx = idx + 1  
          END DO
        END DO
        WRITE(106) vec(0:noccA*nvrtA-1)
      END DO
    END DO
    CLOSE(unit=106)
    !testing testing testgin
    !WRITE(*,*) "TESTING TESTING TESTING"
    !ALLOCATE(testM(0:noccA*nvrtA-1,0:noccA*nvrtA-1))
    !idxB = 0
    !DO j=0,noccA-1
    !  DO b=0,nvrtA-1
    !    idxA = 0
    !    DO i=0,noccA-1
    !      DO a=0,nvrtA-1
    !        testM(idxA,idxB) = Om(a,i,j,b)
    !        idxA = idxA + 1
    !      END DO
    !    END DO 
    !    idxB = idxB + 1
    !  END DO
    !END DO
    !CALL linal_printmat_2Dreal8(testM,noccA*nvrtA,noccA*nvrtA)
    !WRITE(*,*) "bleh"
    !CALL linal_printmat_2Dreal8(CmA,ntot,ntot)
    !WRITE(*,*) "CmA occupied" 
    !WRITE(*,*) CmA(0:ntot-1,0:noccB-1)
    !WRITE(*,*) "CmB virtual"
    !WRITE(*,*) CmA(0:ntot-1,noccB:ntot-1)
    !DEALLOCATE(testM)
    !WRITE(*,*) "TESTING TESTING TESTING"
    DEALLOCATE(Om)
    DEALLOCATE(vec)
   
    !------------------------------
    !Spin Case AA : (uv|ld) -> (ab|ji) -> <aj|bi>
    WRITE(*,*) "Transforming (uv|ld) -> <aj|bi>"

    !first index (uv|ld) -> (av|ld)
    ALLOCATE(Lm(0:nvrtA-1,0:ntot-1,0:ntot-1,0:ntot-1))
    CALL idx1_trans(nvrtA,ntot,ntot,ntot,ntot,Km,CmA(0:ntot-1,noccA:ntot-1),Lm)

    !second index (av|ld) -> (ab|ld) 
    ALLOCATE(Mm(0:nvrtA-1,0:nvrtA-1,0:ntot-1,0:ntot-1))
    CALL idx2_trans(nvrtA,nvrtA,ntot,ntot,ntot,Lm,CmA(0:ntot-1,noccA:ntot-1),Mm)
    DEALLOCATE(Lm)

    !third index (ab|ld) -> (ab|jd) 
    ALLOCATE(Nm(0:nvrtA-1,0:nvrtA-1,0:noccA-1,0:ntot-1))
    CALL idx3_trans(nvrtA,nvrtA,noccA,ntot,ntot,Mm,CmA(0:ntot-1,0:noccA-1),Nm)
    DEALLOCATE(Mm)

    !fourth index (ab|jd) -> (ab|ji) 
    ALLOCATE(Om(0:nvrtA-1,0:nvrtA-1,0:noccA-1,0:noccA-1))
    CALL idx4_trans(nvrtA,nvrtA,noccA,noccA,ntot,Nm,CmA(0:ntot-1,0:noccA-1),Om)
    DEALLOCATE(Nm)

    !We are writing vector by vector
    WRITE(*,*) "Writing to ajbi_AA" 
    ALLOCATE(vec(0:noccA*nvrtA-1))
    OPEN(unit=107,file="ajbi_AA",status="replace",form="unformatted")
    DO j=0,noccA-1
      DO b=0,nvrtA-1
        vec = 0.0D0
        idx = 0
        DO i=0,noccA-1
          DO a=0,nvrtA-1
            vec(idx) = Om(a,b,j,i)
            idx = idx + 1  
          END DO
        END DO
        WRITE(107) vec(0:noccA*nvrtA-1)
      END DO
    END DO
    CLOSE(unit=107)
    DEALLOCATE(Om)
    DEALLOCATE(vec)

    !------------------------------
    !Spin Case AB : (uv|ld) -> (ai|jb) -> <aj|ib>
    WRITE(*,*) "Spin Case AB"
    WRITE(*,*) "Transforming (uv|ld) -> <aj|ib>"
  
    !first index (uv|ld) -> (av|ld)
    ALLOCATE(Lm(0:nvrtA-1,0:ntot-1,0:ntot-1,0:ntot-1))
    CALL idx1_trans(nvrtA,ntot,ntot,ntot,ntot,Km,CmA(0:ntot-1,noccA:ntot-1),Lm)

    !second index (av|ld) -> (ai|ld) 
    ALLOCATE(Mm(0:nvrtA-1,0:noccA-1,0:ntot-1,0:ntot-1))
    CALL idx2_trans(nvrtA,noccA,ntot,ntot,ntot,Lm,CmA(0:ntot-1,0:noccA-1),Mm)
    DEALLOCATE(Lm)

    !third index (ai|ld) -> (ai|jd) 
    ALLOCATE(Nm(0:nvrtA-1,0:noccA-1,0:noccB-1,0:ntot-1))
    CALL idx3_trans(nvrtA,noccA,noccB,ntot,ntot,Mm,CmB(0:ntot-1,0:noccB-1),Nm)
    DEALLOCATE(Mm)

    !fourth index (ai|jd) -> (ai|jb) 
    ALLOCATE(Om(0:nvrtA-1,0:noccA-1,0:noccB-1,0:nvrtB-1))
    CALL idx4_trans(nvrtA,noccA,noccB,nvrtB,ntot,Nm,CmB(0:ntot-1,noccB:ntot-1),Om)
    DEALLOCATE(Nm)
    
    !We are writing vector by vector
    WRITE(*,*) "Writing to ajib_AB" 
    ALLOCATE(vec(0:noccA*nvrtA-1))
    OPEN(unit=108,file="ajib_AB",status="replace",form="unformatted")
    DO j=0,noccB-1
      DO b=0,nvrtB-1
        vec = 0.0D0
        idx = 0
        DO i=0,noccA-1
          DO a=0,nvrtA-1
            vec(idx) = Om(a,i,j,b)
            idx = idx + 1  
          END DO
        END DO
        WRITE(108) vec(0:noccA*nvrtA-1)
      END DO
    END DO
    CLOSE(unit=108)
    DEALLOCATE(Om)
    DEALLOCATE(vec)

    !------------------------------
    !Spin Case BB : (uv|ld) -> (ai|jb) -> <aj|ib>
    WRITE(*,*) "Spin case BB"
    WRITE(*,*) "Transforming (uv|ld) -> <aj|ib>"

    !first index (uv|ld) -> (av|ld)
    ALLOCATE(Lm(0:nvrtB-1,0:ntot-1,0:ntot-1,0:ntot-1))
    CALL idx1_trans(nvrtB,ntot,ntot,ntot,ntot,Km,CmB(0:ntot-1,noccB:ntot-1),Lm)

    !second index (av|ld) -> (ai|ld) 
    ALLOCATE(Mm(0:nvrtB-1,0:noccB-1,0:ntot-1,0:ntot-1))
    CALL idx2_trans(nvrtB,noccB,ntot,ntot,ntot,Lm,CmB(0:ntot-1,0:noccB-1),Mm)
    DEALLOCATE(Lm)

    !third index (ai|ld) -> (ai|jd) 
    ALLOCATE(Nm(0:nvrtB-1,0:noccB-1,0:noccB-1,0:ntot-1))
    CALL idx3_trans(nvrtB,noccB,noccB,ntot,ntot,Mm,CmB(0:ntot-1,0:noccB-1),Nm)
    DEALLOCATE(Mm)

    !fourth index (ai|jd) -> (ai|jb) 
    ALLOCATE(Om(0:nvrtB-1,0:noccB-1,0:noccB-1,0:nvrtB-1))
    CALL idx4_trans(nvrtB,noccB,noccB,nvrtB,ntot,Nm,CmB(0:ntot-1,noccB:ntot-1),Om)
    DEALLOCATE(Nm)
    
    !We are writing vector by vector
    WRITE(*,*) "Writing to ajib_BB" 
    ALLOCATE(vec(0:noccB*nvrtB-1))
    OPEN(unit=109,file="ajib_BB",status="replace",form="unformatted")
    DO j=0,noccB-1
      DO b=0,nvrtB-1
        vec = 0.0D0
        idx = 0
        DO i=0,noccB-1
          DO a=0,nvrtB-1
            vec(idx) = Om(a,i,j,b)
            idx = idx + 1  
          END DO
        END DO
        WRITE(109) vec(0:noccB*nvrtB-1)
      END DO
    END DO
    CLOSE(unit=109)
    !testing testing testgin
    !WRITE(*,*) "TESTING TESTING TESTING"
    !ALLOCATE(testM(0:noccB*nvrtB-1,0:noccB*nvrtB-1))
    !idxB = 0
    !DO j=0,noccB-1
    !  DO b=0,nvrtB-1
    !    idxA = 0
    !    DO i=0,noccB-1
    !      DO a=0,nvrtB-1
    !        testM(idxA,idxB) = Om(a,i,j,b)
    !        idxA = idxA + 1
    !      END DO
    !    END DO 
    !    idxB = idxB + 1
    !  END DO
    !END DO
    !CALL linal_printmat_2Dreal8(testM,noccB*nvrtB,noccB*nvrtB)
    !WRITE(*,*) "bleh"
    !CALL linal_printmat_2Dreal8(CmB,ntot,ntot)
    !WRITE(*,*) "CmB occupied" 
    !WRITE(*,*) CmB(0:ntot-1,0:noccB-1)
    !WRITE(*,*) "CmB virtual"
    !WRITE(*,*) CmB(0:ntot-1,noccB:ntot-1)
    !DEALLOCATE(testM)
    !WRITE(*,*) "TESTING TESTING TESTING"
    DEALLOCATE(Om)
    DEALLOCATE(vec)
   
    !------------------------------
    !Spin Case BB : (uv|ld) -> (ab|ji) -> <aj|bi>
    WRITE(*,*) "Transforming (uv|ld) -> <aj|bi>"

    !first index (uv|ld) -> (av|ld)
    ALLOCATE(Lm(0:nvrtB-1,0:ntot-1,0:ntot-1,0:ntot-1))
    CALL idx1_trans(nvrtB,ntot,ntot,ntot,ntot,Km,CmB(0:ntot-1,noccB:ntot-1),Lm)

    !second index (av|ld) -> (ab|ld) 
    ALLOCATE(Mm(0:nvrtB-1,0:nvrtB-1,0:ntot-1,0:ntot-1))
    CALL idx2_trans(nvrtB,nvrtB,ntot,ntot,ntot,Lm,CmB(0:ntot-1,noccB:ntot-1),Mm)
    DEALLOCATE(Lm)

    !third index (ab|ld) -> (ab|jd) 
    ALLOCATE(Nm(0:nvrtB-1,0:nvrtB-1,0:noccB-1,0:ntot-1))
    CALL idx3_trans(nvrtB,nvrtB,noccB,ntot,ntot,Mm,CmB(0:ntot-1,0:noccB-1),Nm)
    DEALLOCATE(Mm)

    !fourth index (ab|jd) -> (ab|ji) 
    ALLOCATE(Om(0:nvrtB-1,0:nvrtB-1,0:noccB-1,0:noccB-1))
    CALL idx4_trans(nvrtB,nvrtB,noccB,noccB,ntot,Nm,CmB(0:ntot-1,0:noccB-1),Om)
    DEALLOCATE(Nm)

    !We are writing vector by vector
    WRITE(*,*) "Writing to ajbi_BB" 
    ALLOCATE(vec(0:noccB*nvrtB-1))
    OPEN(unit=110,file="ajbi_BB",status="replace",form="unformatted")
    DO j=0,noccB-1
      DO b=0,nvrtB-1
        vec = 0.0D0
        idx = 0
        DO i=0,noccB-1
          DO a=0,nvrtB-1
            vec(idx) = Om(a,b,j,i)
            idx = idx + 1  
          END DO
        END DO
        WRITE(110) vec(0:noccB*nvrtB-1)
      END DO
    END DO
    CLOSE(unit=110)
    DEALLOCATE(Om)
    DEALLOCATE(vec)

    !------------------------------
    !Cleanup
    DEALLOCATE(CmA)
    DEALLOCATE(CmB)
    DEALLOCATE(Km)

  END SUBROUTINE slow_ao2mo_CIS_UHF

!---------------------------------------------------------------------
!       slow_ao2mo_CIS_RHF
!               James H. Thorpe
!               Nov 3., 2018
!       -calculates AO -> MO transforms for RHF CIS calculations
!---------------------------------------------------------------------
  SUBROUTINE slow_ao2mo_CIS_RHF(noccA,nvrtA,ntot)
  ! Variables
  ! noccA       :       int, number of A occupied orbitals
  ! nvrtA       :       int, number of A virtual orbitals
  ! ntot        :       int, total nubmer of orbitals
    IMPLICIT NONE
    ! Inout
    INTEGER, INTENT(IN) :: noccA,nvrtA,ntot

  END SUBROUTINE
!---------------------------------------------------------------------
!       print_moints
!               James H. Thorpe
!               Nov. 1, 2018
!       -prints the numbers of MO integrals needed for each spin case
!---------------------------------------------------------------------
  ! noccA,noccB         :       int, number of occupied A,B orbitals
  ! nvrtA,nvrtB         :       int, number of virtual A,B orbitals
  ! ntot                :       int, total number of orbitals
  ! options             :       1D int, options array

  SUBROUTINE print_moints(noccA,noccB,nvrtA,nvrtB,ntot,options)
    IMPLICIT NONE
    !inout
    INTEGER, DIMENSION(0:), INTENT(IN) :: options
    INTEGER, INTENT(IN) :: noccA,noccB,nvrtA,nvrtB,ntot
    !internal
    INTEGER :: i,j,k,l,a,b,c,d
    INTEGER ::no,nv

    WRITE(*,*) "MO Integrals"
    WRITE(*,*) "-------------------------------"
    !RHF print
    IF (options(3) .EQ. 0) THEN
      !AA
      WRITE(*,*) "Spin Case AA"
      !the incredibly lazy way of doing things
      no = 0
      nv = 0
      DO i=0,noccA-2
        DO j=i+1,noccA-1
          no = no + 1
        END DO 
      END DO
      DO a=0,nvrtA-2
        DO b=a+1,nvrtA-1
          nv = nv + 1
        END DO 
      END DO
      WRITE(*,*) "<ij|ab>         ", no*nv
      !AB
      WRITE(*,*) "Spin Case AB"
      WRITE(*,*) "<ij|ab>         ",noccA*noccA*nvrtA*nvrtA
      WRITE(*,*)
    END IF

  END SUBROUTINE

!---------------------------------------------------------------------
!	idx1_trans
!		James H. Thorpe
!		Nov 26,2018
!	-transforms index 1 of array A into array B
!	- NOTE: be careful that you pass in the x vector with the
!         correct bounds!
!---------------------------------------------------------------------
  !Values
  ! n1..5	: int, upper bounds on index 1-5
  ! A		: 2Dreal8, matrix to be transformed
  ! x		: 1Dreal8, vector that defines coefs of transform
  ! B		: 2Dreal8, new array that has been transformed
  SUBROUTINE idx1_trans(n1,n2,n3,n4,n5,A,x,B)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:,0:,0:,0:),INTENT(INOUT) :: B
    REAL(KIND=8), DIMENSION(0:,0:,0:,0:),INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: n1,n2,n3,n4,n5
    REAL(KIND=8) :: temp
    INTEGER :: p,q,r,s,t
    B = 0.0D0
    DO p=0,n1-1
      DO q=0,n2-1
        DO r=0,n3-1
          DO s=0,n4-1
            temp = 0.0D0
            DO t=0,n5-1
              temp = temp + A(t,q,r,s)*x(t,p)
            END DO
            B(p,q,r,s) = temp
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE idx1_trans

!---------------------------------------------------------------------
!	idx2_trans
!		James H. Thorpe
!		Nov 26, 2018
!	-transforms index 2 of array A into array B
!	- NOTE: be careful that you pass in the x vector with the
!         correct bounds!
!---------------------------------------------------------------------
  !Values
  ! n1..5	: int, upper bounds on index 1-5
  ! A		: 2Dreal8, matrix to be transformed
  ! x		: 1Dreal8, vector that defines coefs of transform
  ! B		: 2Dreal8, new array that has been transformed
  SUBROUTINE idx2_trans(n1,n2,n3,n4,n5,A,x,B)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:,0:,0:,0:),INTENT(INOUT) :: B
    REAL(KIND=8), DIMENSION(0:,0:,0:,0:),INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: n1,n2,n3,n4,n5
    REAL(KIND=8) :: temp
    INTEGER :: p,q,r,s,t
    B = 0.0D0
    DO p=0,n1-1
      DO q=0,n2-1
        DO r=0,n3-1
          DO s=0,n4-1
            temp = 0.0D0
            DO t=0,n5-1
              temp = temp + A(p,t,r,s)*x(t,q)
            END DO
            B(p,q,r,s) = temp
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE idx2_trans

!---------------------------------------------------------------------
!	idx3_trans
!		James H. Thorpe
!		Nov 26, 2018
!	-transforms index 3 of array A into array B
!	- NOTE: be careful that you pass in the x vector with the
!         correct bounds!
!---------------------------------------------------------------------
  !Values
  ! n1..5	: int, upper bounds on index 1-5
  ! A		: 2Dreal8, matrix to be transformed
  ! x		: 1Dreal8, vector that defines coefs of transform
  ! B		: 2Dreal8, new array that has been transformed
  SUBROUTINE idx3_trans(n1,n2,n3,n4,n5,A,x,B)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:,0:,0:,0:),INTENT(INOUT) :: B
    REAL(KIND=8), DIMENSION(0:,0:,0:,0:),INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: n1,n2,n3,n4,n5
    REAL(KIND=8) :: temp
    INTEGER :: p,q,r,s,t
    B = 0.0D0
    DO p=0,n1-1
      DO q=0,n2-1
        DO r=0,n3-1
          DO s=0,n4-1
            temp = 0.0D0
            DO t=0,n5-1
              temp = temp + A(p,q,t,s)*x(t,r)
            END DO
            B(p,q,r,s) = temp
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE idx3_trans

!---------------------------------------------------------------------
!	idx4_trans
!		James H. Thorpe
!		Nov 26, 2018
!	-transforms index 4 array A into array B
!	- NOTE: be careful that you pass in the x vector with the
!         correct bounds!
!---------------------------------------------------------------------
  !Values
  ! n1..5	: int, upper bounds on index 1-5
  ! A		: 2Dreal8, matrix to be transformed
  ! x		: 1Dreal8, vector that defines coefs of transform
  ! B		: 2Dreal8, new array that has been transformed
  SUBROUTINE idx4_trans(n1,n2,n3,n4,n5,A,x,B)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:,0:,0:,0:),INTENT(INOUT) :: B
    REAL(KIND=8), DIMENSION(0:,0:,0:,0:),INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: n1,n2,n3,n4,n5
    REAL(KIND=8) :: temp
    INTEGER :: p,q,r,s,t
    B = 0.0D0
    DO p=0,n1-1
      DO q=0,n2-1
        DO r=0,n3-1
          DO s=0,n4-1
            temp = 0.0D0
            DO t=0,n5-1
              temp = temp + A(p,q,r,t)*x(t,s)
            END DO
            B(p,q,r,s) = temp
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE idx4_trans
!---------------------------------------------------------------------
!       linal_printmat_2Dreal8
!               James H. Thorpe
!               Nov 25, 2018
!       -prints matrix in row major form
!---------------------------------------------------------------------
  SUBROUTINE linal_printmat_2Dreal8(A,N,M)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:N-1,0:M-1), INTENT(IN) :: A
    INTEGER(KIND=4), INTENT(IN) :: N,M
    INTEGER(KIND=4) :: i
    DO i=0,M-1
      WRITE(*,*) A(i,0:N-1)
    END DO
  END SUBROUTINE
!---------------------------------------------------------------------

END PROGRAM ao2mo

