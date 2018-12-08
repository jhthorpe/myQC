!//////////////////////////////////////////////////////////////////
!//            Performs 1 e- property integrals for myQC 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//		
!//		Integration scheme a la McMurchie,Davidson 1978	
!//		
!//	WORK NOTE - location of coefficients in bas currently hardcoded	
!///////////////////////////////////////////////////////////////////

!=====================================================================
!			MAIN 
PROGRAM int1p
  USE env
  USE basis
  USE auxilary

  IMPLICIT NONE

  ! Values
  ! xyz		: 2D dp, array of nuclear positions
  ! atoms	: 1D int, array of which atom is which
  ! fmem	: dp, free memory left in MB
  ! nnuc	: int, number of nuclii
  ! nelcA,B	: int, number of electrons of spin A,B
  ! norb	: int, number of orbitals in molecule
  ! npri	: int, number of primatives
  ! bas		: 2D dp, basis for each atom: atom : orbital : [d,a]
  ! basinfo	: 2D int, array of basis information
  ! options	: 1D int, array of options
  ! Px		: 3D real*8, dipole integrals 
  ! set		: 2D dp, set of exponent coefficients on each nuclei
  ! setinfo	: 2D int, array of set information

  ! Variables  
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: Px
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: bas,set
  INTEGER, ALLOCATABLE, DIMENSION(:) :: basinfo, setinfo
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms,options 
  REAL(KIND=8) :: timeS, timeF, fmem
  INTEGER :: nnuc,nelcA,nelcB,i,j,k,norb,npri,stat1,stat2,maxN,maxL
  LOGICAL :: flag1,flag2,flag

! input managment 
  CALL CPU_TIME(timeS)
  WRITE(*,*) "int1p called"
  CALL getenv(nnuc,nelcA,nelcB,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

!format stuff
999 FORMAT(1x,A31,F8.5)
998 FORMAT(1x,A21,2x,I4)
997 FORMAT(1x,A18,F8.5)

!the actual stuff
! construct the basis
  CALL buildBasis(options(2),atoms,bas,basinfo,set,setinfo,.FALSE.,maxN,maxL)
  CALL nmem(fmem)

! 1) calculate property integrals 
  !get number of orbitals
  npri = 0
  norb = basinfo(1) 
  DO i=0,norb-1
      npri = npri + basinfo(1+i*5+4) 
  END DO
  WRITE(*,*) "               STARTING PROPERTY INTEGRALS" 
  WRITE(*,*) "------------------------------------------------------------"

  WRITE(*,999) "Allocating space for int1p (MB)", 2*norb*norb*8/1.0D6
  fmem = fmem - 2*norb*norb*8.0/1.D6
  IF (fmem .LT. 0.0D0) THEN
    CALL EXECUTE_COMMAND_LINE('touch error')
    WRITE(*,*) "int1p: max memory reached"
    STOP
  ELSE
    ALLOCATE(Px(0:2,0:norb-1,0:norb-1),STAT=stat1)
    IF(stat1+stat2 .NE. 0) THEN
      CALL EXECUTE_COMMAND_LINE('touch error')
      WRITE(*,*) "int1p: could not allocate memory"
      STOP
    END IF
  END IF
  WRITE(*,*)
  CALL nmem(fmem)

  INQUIRE(file='Pxuv',EXIST=flag)
  IF (.NOT. flag) THEN
    !1) If not there, calculated Overlap and Fock 
    CALL proc1p(Px,bas,basinfo,atoms,options,fmem,nnuc,xyz,norb,set,setinfo)
    WRITE(*,*) "Dipole moment integrals written to Pxuv"
  ELSE
    WRITE(*,*) "Reading Dipole integrals from Pxuv"
  END IF

  DEALLOCATE(Px)
  DEALLOCATE(bas)
  DEALLOCATE(set)
  DEALLOCATE(basinfo)
  DEALLOCATE(setinfo)

! output
  fmem = fmem + 2*norb*norb*8/1.0D6
  CALL setenv(atoms,xyz,fmem,options)
  CALL CPU_TIME(timeF)

  CONTAINS 

!---------------------------------------------------------------------
!		proc1p
!			James H. Thorpe
!			Dec 8, 2018	
!		-process 1e- property integrals
!---------------------------------------------------------------------
    ! Values
    ! Px	: 3D real*8, x,y,z P(u,v) property integrals
    ! xyz	: 2D dp, array of nuclear positions
    ! atoms	: 1D int, array of which atom is which
    ! fmem	: dp, free memory left in MB
    ! nnuc	: int, number of nuclii
    ! norb	: int, number of orbitals in molecule
    ! bas	: 1D dp, basis weights for each atom
    ! basinfo	: 1D int, array of basis information
    ! options	: 1D int, array of options
    ! set	: 1D dp, array of exponential alphas
    ! setinfo	: 1D int, stores information about set
    ! coef	: 3D dp, array of coefficients for overlap, (xyz,i,j,k)
    ! PA,PB,PP	: 1D dp, array of distances to overlap, PP = overlap
    ! EIJ	: dp, gaussian at molecular center
    ! val 	: dp, current evaluation of integral
    ! p		: dp, addition of two coeffiecients
    ! aa,bb	: dp, coefficients of Guassians for atoms a and b
    ! OpS	: int, max orbitals per set
    ! Dk	: 1D dp, tracks nonzero EIJ*N*L*M
    ! Ck	: 1D int, nonzero combinations of Nk,Lk,Mk
    ! Ok	: 1D int, lists orbitals of Dk,Ck
  SUBROUTINE proc1p(Px,bas,basinfo,atoms,options,fmem,nnuc,xyz,norb,set,setinfo)
    IMPLICIT NONE

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: Px
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: bas, set
    INTEGER, DIMENSION(0:), INTENT(IN) :: basinfo, setinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: options, atoms
    REAL(KIND=8), INTENT(IN) :: fmem
    INTEGER, INTENT(IN) :: norb,nnuc

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: coef
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: nucpos
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Dk
    INTEGER, DIMENSION(:), ALLOCATABLE :: Ck,Ok
    REAL(KIND=8), DIMENSION(0:2) :: PA,PB,AB,PP
    INTEGER, DIMENSION(0:2) :: la,lb,amax,bmax
    REAL(KIND=8) :: EIJ, valSb, valFb, p, m, aa, bb, tempSb, tempFb
    REAL(KIND=8) :: timeS, timeF
    INTEGER :: setl,nset,OpS,kmax
    INTEGER ::  a,b,u,v,i,j,k,l

    CALL CPU_TIME(timeS)
    WRITE(*,*) "constructing Dipole Matrix" 

996 FORMAT(1x,A37,2x,F8.5)

    !zero Px 
    Px = 0.0D0

    nset = setinfo(0)
    setl = setinfo(1)
    OpS = basinfo(0)

    kmax = Ops**2 * (2*maxL*(2*maxL+1)/2)**3

    ALLOCATE(Dk(0:kmax)) 
    ALLOCATE(Ck(0:kmax))
    ALLOCATE(Ok(0:2*kmax+1))

    ! left set
    DO a=0,nset-1
      
      !get data
      aa = set(a)                            !alpha a
      u = setinfo(1+a*setl+3)                !center number
      DO l=0,2                               !set amax values
        amax(l) = setinfo(1+a*setl+2)        !get max ang qn
        la(l) = amax(l)
      END DO
     
      ! right set 
      DO b=0,nset-1 

        !get data
        bb = set(b)                           !alpha b
        v = setinfo(1+b*setl+3)               !center number
        DO l=0,2
          bmax(l) = setinfo(1+b*setl+2)
          lb(l) = setinfo(1+b*setl+2)
        END DO

        ! 1) get overlap location
        p = aa + bb 
        m = aa * bb
        DO i=0,2
          AB(i) = xyz(u,i) - xyz(v,i)
          PP(i) = (aa*xyz(u,i) + bb*xyz(v,i))/p
          PA(i) = PP(i) - xyz(u,i)
          PB(i) = PP(i) - xyz(v,i)
        END DO 

        ! screen for sufficiently small constant
        EIJ = EXP(-m*(AB(0)**2.0D0 + AB(1)**2.0D0 + AB(2)**2.0D0)/p) 
        IF (EIJ .LT. 1.0D-14) THEN
          CYCLE
        END IF 

        kmax = -1

        CALL getcoef_dipole(coef,PA,PB,aa,bb,amax,bmax) 
        CALL getDk(coef,setinfo(1+a*setl+1:1+(a+1)*setl),setinfo(1+b*setl+1:1+(b+1)*setl), &
        bas(a*OpS:(a+1)*Ops-1),bas(b*Ops:(b+1)*OpS-1),basinfo,Dk,Ck,Ok,kmax,EIJ,setl,aa,bb)

        !dipole
        CALL dipole(Px,u,v,a,b,p,bas(a*OpS:(a+1)*Ops-1),bas(b*Ops:(b+1)*OpS-1),basinfo,coef,&
        setinfo(1+a*setl+1:1+(a+1)*setl),setinfo(1+b*setl+1:1+(b+1)*setl),aa,bb,EIJ,nnuc,PP,AB)

        DEALLOCATE(coef)

      END DO                                   !loop over left orbital
    END DO                                     !loop over right orbital

    OPEN(unit=100,file='Pxuv',status='replace',access='sequential')
      WRITE(100,*) Px(0,:,:)    
      WRITE(100,*) Px(1,:,:)
      WRITE(100,*) Px(2,:,:)
    CLOSE(unit=1)

    IF (options(7) .GE. 3) CALL print_Px(Px,norb)
 
    DEALLOCATE(Dk)
    DEALLOCATE(Ck)
    DEALLOCATE(Ok)

    CALL CPU_TIME(timeF)
    WRITE(*,996) "Dipole Integrals constructed in (s) :", (timeF-timeS)

  END SUBROUTINE proc1p

!---------------------------------------------------------------------
!      dipole
!		James H. Thorpe
!		Dec 8, 2018 
!	-calculates primative dipole integrals
!---------------------------------------------------------------------
  ! Values
  ! Px          : 3D real*8, dipole matrix (coord,u,v)
  ! a,b         : int, set number we're on
  ! la,lb       : 1D int, angular quantum numbers 
  ! ta,tb       : int, tracking a and b
  ! basa,basb   : 1D dp, array of basis set weights for sets a and b
  ! seta,setb   : 1D int, array of setinfo for sets a and b
  ! nncu	: int, number of nuclei
  ! PP		: 1D real*8, list of overlap locations, [x,y,z]
  ! CC		: 1D real*8, radius between center A and center B (x,y,z)

  SUBROUTINE dipole(Px,u,v,a,b,p,basa,basb,basinfo,coef,seta,setb,aa,bb,EIJ,nnuc,PP,CC)
    IMPLICIT NONE
    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931
    ! Inout
    REAL(KIND=8), DIMENSION(0:,-2:,-2:,-2:), INTENT(IN) :: coef
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: Px
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: basa, basb,PP,CC
    INTEGER, DIMENSION(0:), INTENT(IN) :: basinfo, seta, setb
    REAL(KIND=8), INTENT(IN) :: aa, bb, EIJ, p
    INTEGER, INTENT(IN) :: u,v,a,b, nnuc

    ! internal
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: nucpos
    REAL(KIND=8), DIMENSION(0:2) :: CP
    INTEGER, DIMENSION(0:2) :: la,lb
    REAL(KIND=8) :: temp
    INTEGER :: dummy
    INTEGER :: c,i,j,k,orba,orbb,ori,prima,primb
    
    ALLOCATE(nucpos(0:nnuc-1,0:2))

    ! get nuclear positions
    OPEN(unit=1,file='nucpos',status='old',access='sequential')
    DO c=0,nnuc-1
      READ(1,*) dummy, nucpos(c,0:2)
    END DO
    CLOSE(unit=1)

    !update each element in set
    DO i=0,seta(0)-1 !go through set A
      orba = seta(3+i) !id of orbital

      DO j=0,setb(0)-1 !go through set B
        orbb = setb(3+j) !id of orbital

        ! get angular quantum numbers for each orbital 
        ori = basinfo(1+5*orba+3)
        !S-TYPE
        IF (ori .EQ. -1) THEN
          la = [basinfo(1+5*orba+2),basinfo(1+5*orba+2),basinfo(1+5*orba+2)]
        !P-TYPE
        ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN
          la = [0, 0, 0]
          la(ori) = basinfo(1+5*orba+2)
        END IF

        ori = basinfo(1+5*orbb+3)
        !S-TYPE
        IF (ori .EQ. -1) THEN
          lb = [basinfo(1+5*orbb+2),basinfo(1+5*orbb+2),basinfo(1+5*orbb+2)]
        !P-TYPE
        ELSE IF (ori .GE. 0 .AND. ori .LE. 2) THEN
          lb = [0,0,0]
          lb(ori) = basinfo(1+5*orbb+2)
        END IF

        !get distances to center
        DO k=0,2
          CP(k) = CC(k) - PP(k)
        END DO

        !xpart
        temp = EIJ*(Pi/p)**(3.0D0/2.0D0)*basa(i)*basb(j)          !pre-exponential and basis weights
        temp = temp * gtoD(basinfo(1+5*orba+2),aa)                !basis set coefficients
        temp = temp * gtoD(basinfo(1+5*orbb+2),bb)                !basis set coefficeints
        temp = temp * coef(1,0,la(1),lb(1))*coef(2,0,la(2),lb(2))
        temp = temp * (coef(0,1,la(0),lb(0)) + CP(0)*coef(0,0,la(0),lb(0)))
        Px(0,u,v) =  Px(0,u,v) + temp

        !ypart
        temp = EIJ*(Pi/p)**(3.0D0/2.0D0)*basa(i)*basb(j)          !pre-exponential and basis weights
        temp = temp * gtoD(basinfo(1+5*orba+2),aa)                !basis set coefficients
        temp = temp * gtoD(basinfo(1+5*orbb+2),bb)                !basis set coefficeints
        temp = temp * coef(0,0,la(0),lb(0))*coef(2,0,la(2),lb(2))
        temp = temp * (coef(1,1,la(1),lb(1)) + CP(1)*coef(1,0,la(1),lb(1)))
        Px(1,u,v) =  Px(1,u,v) + temp

        !zpart
        temp = EIJ*(Pi/p)**(3.0D0/2.0D0)*basa(i)*basb(j)          !pre-exponential and basis weights
        temp = temp * gtoD(basinfo(1+5*orba+2),aa)                !basis set coefficients
        temp = temp * gtoD(basinfo(1+5*orbb+2),bb)                !basis set coefficeints
        temp = temp * coef(0,0,la(0),lb(0))*coef(1,0,la(1),lb(1))
        temp = temp * (coef(2,1,la(2),lb(2)) + CP(2)*coef(2,0,la(2),lb(2)))
        Px(2,u,v) =  Px(2,u,v) + temp

      END DO
    END DO


  END SUBROUTINE dipole
!---------------------------------------------------------------------
!	print_Px
!		James H. Thorpe
!		Dec 8, 2018
!	-print Px integrals
!---------------------------------------------------------------------
  ! Px		:  3D real*8, dipole integrals (x,u,v)
  ! norb	:  int, number of orbitals
  SUBROUTINE print_Px(Px,norb)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: Px
    INTEGER, INTENT(IN) :: norb
    INTEGER :: u,v
999 FORMAT(4x,I3,1x,I3,1x,F15.8,1x,F15.8,4x,F15.8)
998 FORMAT(4x,A50)
    WRITE(*,*) "One electron dipole moments"
    WRITE(*,998) "  u   v      Dx              Dy                 Dz"
    DO u=0,norb-1
      DO v=0,norb-1
        WRITE(*,999) u,v,Px(0,u,v),Px(1,u,v),Px(2,u,v)
      END DO
    END DO 

  END SUBROUTINE print_Px
!---------------------------------------------------------------------


END PROGRAM int1p
