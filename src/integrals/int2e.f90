!//////////////////////////////////////////////////////////////////
!//            Performs 2 e- integrals for myQC 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//             
!//             Integration scheme a la McMurchie,Davidson 1978 
!//             
!//     WORK NOTE - location of coefficients in bas currently hardcoded 
!///////////////////////////////////////////////////////////////////

!=====================================================================
!                       MAIN 
PROGRAM int2e
  USE env
  USE basis
  USE aux

  IMPLICIT NONE

  ! Values
  ! xyz         : 2D dp, array of nuclear positions
  ! atoms       : 1D int, array of which atom is which
  ! fmem        : dp, free memory left in MB
  ! nnuc        : int, number of nuclii
  ! nelc        : int, number of electrons
  ! norb        : int, number of orbitals in molecule
  ! npri        : int, number of primatives
  ! bas         : 2D dp, basis for each atom: atom : orbital : [d,a]
  ! basinfo     : 2D int, array of basis information
  ! options     : 1D int, array of options
  ! set		: 2D dp, array of exponential coefficients of sets: atom,orbital numbs
  ! setinfo	: 2D int, array of set information
  ! maxN	: int, max principle quantum number of basis
  ! maxL	: int, max angular quantum number of basis

  ! Variables  
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xyz
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: bas,set
  INTEGER, ALLOCATABLE, DIMENSION(:) :: basinfo,setinfo
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms,options
  REAL(KIND=8) :: timeS, timeF, fmem
  INTEGER :: nnuc,nelc,i,j,k,norb,npri,stat,maxN,maxL
  LOGICAL :: flag1,flag2,flag

  ! input managment 
  CALL CPU_TIME(timeS)
  WRITE(*,*)
  WRITE(*,*) "int2e called"
  CALL getenv(nnuc,nelc,xyz,atoms,fmem,options)
  INQUIRE(file='error',EXIST=flag)
  IF (flag) STOP

  ! build the basis set
  CALL buildBasis(options(2),atoms,bas,basinfo,set,setinfo,.FALSE.,maxN,maxL)

  ! check that Fock has been created 
  INQUIRE(file='Fuv',EXIST=flag1)
 
  ! check that we need to do these calculations
  INQUIRE(file='Fold',EXIST=flag2)
 
  ! logic gate
  IF (flag1 .AND. flag2) THEN
    WRITE(*,*) "Reading two electron integrals from intermediate"
    !we aren't actually reading anything, but they don't need to know that ;)"
  ELSE
    WRITE(*,*) "Constructing two electron integrals"
    CALL proc2e(bas,basinfo,atoms,options,fmem,nnuc,xyz,set,setinfo,maxL)
  END IF


  CONTAINS
!===================================================================
!                       SUBROUTINES

!---------------------------------------------------------------------
! 		Processes the two electron integrals
!---------------------------------------------------------------------
  SUBROUTINE proc2e(bas,basinfo,atoms,options,fmem,nnuc,xyz,set,setinfo,maxL)
    IMPLICIT NONE
    ! Values
    ! xyz	: 2D dp, array of nuclear positions
    ! atoms	: 1D int, array of which atom is which
    ! fmem	: dp, free memory left in MB
    ! nnuc	: int, number of nuclii
    ! norb	: int, number of orbitals in molecule
    ! bas	: 1D dp, basis for each atom: atom : orbital : [d,a]
    ! basinfo	: 1D int, array of basis information
    ! set	: 1D dp, array of exponential coefficients
    ! setinfo	: 1D int, array of information about sets
    ! options	: 1D int, array of options
    ! II	: 3D dp, intermediate array 
    ! XX	: 4D dp, two electron repulsion integrals, on orbitals i,j,g,h
    ! maxL	: int, max principle,angular quantum numbers
    ! PP,QQ,PQ	: 1D dp, arrays of overlap integral distances
    ! la,lb,...	: 1D int, arrays that contain the max angular quantum numbers of sets 
    ! ABk,CDk	: 1D int, contains nonzero indices of coefficients
    ! kmax	: int, max k value that is nonzero
    ! Dk,Dkp	: 1D dp, array of nonzero coefficients for sets AB and CD
    ! Ok,Okp	: 1D int, array of orbitals of sets AB, CD, stored {left, right}
 
    !Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: bas,set
    INTEGER, DIMENSION(0:), INTENT(IN) :: basinfo,setinfo
    INTEGER, DIMENSION(0:), INTENT(IN) :: options, atoms
    REAL(KIND=8), INTENT(INOUT) :: fmem
    INTEGER, INTENT(IN) :: nnuc,maxL

    !Internal
    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: XX,coefAB,coefCD
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: II
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Dk, Dkp
    REAL(KIND=8), DIMENSION(0:2) :: PA, PB, AB, QC, QD, CD
    REAL(KIND=8), DIMENSION(0:2) :: PP,QQ,PQ
    INTEGER, DIMENSION(:), ALLOCATABLE :: Ck,Ckp,Ok,Okp
    INTEGER, DIMENSION(0:2) :: la,lb,lc,ld,na,nb,nc,nd
    REAL(KIND=8) :: timeS,timeF,temp,EIJ,EGH
    REAL(KIND=8) :: aa,bb,cc,dd,p,q
    INTEGER :: nset,setl,OpS,stat1,stat2,setK,kmaxAB,kmaxCD 
    INTEGER :: orba,orbb,orbc,orbd,norb
    INTEGER :: a,b,c,d,g,h,i,j,k,m,n,s,t,u,v

    CALL CPU_TIME(timeS)

998 FORMAT(1x,A47,2x,F8.5)
999 FORMAT(1x,A41,2x,F8.4)
997 FORMAT(4x,A22)
996 FORMAT(4x,I3,1x,I3,1x,I3,1x,I3,4x,F15.8)
    
    nset = setinfo(0)
    setl = setinfo(1)
    OpS = basinfo(0)
    norb = basinfo(1)
    setK = Ops**2*(2*maxL*(2*maxL+1)/2)**3

    ALLOCATE(Ck(0:setK))
    ALLOCATE(Ckp(0:setK))
    ALLOCATE(Dk(0:setK)) 
    ALLOCATE(Dkp(0:setK))
    ALLOCATE(Ok(0:2*setK+1))
    ALLOCATE(Okp(0:2*setK+1))

    !Iuv will be my intermediate file for the integrals
    OPEN(unit=42,file='Iuv',status='replace',access='sequential',form='unformatted') 
 
    !allocate memory
    temp = (norb**(4))*8.0/1.0E6 
    temp = temp + setK*norb*norb*8.0D0/1.0E6
    WRITE(*,998) "Allocating space for intermediate matrices (MB)", temp 
    ALLOCATE(XX(0:norb-1,0:norb-1,0:norb-1,0:norb-1),STAT=stat1)
    ALLOCATE(II(0:setK,0:norb-1,0:norb-1),STAT=stat2)
    IF (stat1+stat2 .NE. 0) THEN
      WRITE(*,*) "Could not allocate memory in int2e:proc2e"
      CALL EXECUTE_COMMAND_LINE('touch error')
      RETURN
    END IF
    fmem = fmem - temp
    CALL nmem(fmem)

    !zero XX
    DO i=0,norb-1
      DO j=0,norb-1
        DO g=0,norb-1
          XX(i,j,g,:) = (/ (0.0D0, h=0, norb-1) /) 
        END DO
      END DO
    END DO

    !"sum" over a,b sets
    DO a=0,nset-1
      aa = set(a)                            !alpha a
      u = setinfo(1+a*setl+3)                !center number
      !max angular momentum
      la(:) = [setinfo(1+a*setl+2),setinfo(1+a*setl+2),setinfo(1+a*setl+2)]

      DO b=0,nset-1      
        bb = set(b)                           !alpha b
        v = setinfo(1+b*setl+3)
        lb = [setinfo(1+b*setl+2),setinfo(1+b*setl+2),setinfo(1+b*setl+2)] 

        kmaxAB = -1

        !zero II
        DO k=0,setK
          DO g=0,norb-1
            II(k,g,:) = (/ (0.0D0, h=0, norb-1) /)
          END DO
        END DO

        !left overlap location
        p = aa + bb
        m = aa * bb
        DO i=0,2
          AB(i) = xyz(u,i) - xyz(v,i)
          PP(i) = (aa*xyz(u,i) + bb*xyz(v,i))/p
          PA(i) = PP(i) - xyz(u,i)
          PB(i) = PP(i) - xyz(v,i)
        END DO
        EIJ = EXP(-m*(AB(0)**2.0D0 + AB(1)**2.0D0 + AB(2)**2.0D0)/p)

        CALL getcoef(coefAB,PA,PB,aa,bb,la,lb)
        CALL getDk(coefAB,setinfo(1+a*setl+1:1+(a+1)*setl),setinfo(1+b*setl+1:1+(b+1)*setl), &
        bas(a*OpS:(a+1)*Ops-1),bas(b*Ops:(b+1)*OpS-1),basinfo,Dk,Ck,Ok,kmaxAB,EIJ,setl,aa,bb)

        !"sum" over c,d sets
        DO c=0,nset-1
          cc = set(c)                            !alpha c
          s = setinfo(1+c*setl+3)
          lc(:) = [setinfo(1+c*setl+2),setinfo(1+c*setl+2),setinfo(1+c*setl+2)]

          DO d=0,nset-1
            dd = set(d)                          !alpha d
            t = setinfo(1+d*setl+3)
            ld(:) = [setinfo(1+d*setl+2),setinfo(1+d*setl+2),setinfo(1+d*setl+2)]

            !Right overlap location
            q = cc + dd
            n = cc * dd
            DO i=0,2
              CD(i) = xyz(s,i) - xyz(t,i)
              QQ(i) = (cc*xyz(s,i) + dd*xyz(t,i))/q
              QC(i) = QQ(i) - xyz(s,i)
              QD(i) = QQ(i) - xyz(t,i)
            END DO
            EGH = EXP(-n*(CD(0)**2.0D0 + CD(1)**2.0D0 + CD(2)**2.0D0)/q)
            IF (EGH * EIJ .LT. 1.0D-14) CYCLE    
     
            kmaxCD = -1

            CALL getcoef(coefCD,QC,QD,cc,dd,lc,ld)
            CALL getDk(coefCD,setinfo(1+c*setl+1:1+(c+1)*setl),setinfo(1+d*setl+1:1+(d+1)*setl), &
            bas(c*OpS:(c+1)*Ops-1),bas(d*Ops:(d+1)*OpS-1),basinfo,Dkp,Ckp,Okp,kmaxCD,EGH,setl,cc,dd)

            ! get overlap of overlaps 
            DO i=0,2
              PQ(i) = PP(i) - QQ(i)
            END DO
   
            CALL clmII(II,p,q,la,lb,lc,ld,PQ,Dk,Dkp,Ck,Ckp,kmaxAB,kmaxCD,&
            setinfo(1+c*setl+1:1+(c+1)*setl),setinfo(1+d*setl+1:1+(d+1)*setl),setl)

!            CALL clmII(II,a,b,c,d,p,q,PQ,setinfo(1+a*setl+1:1+(a+1)*setl),&
!            setinfo(1+b*setl+1:1+(b+1)*setl),setinfo(1+c*setl+1:1+(c+1)*setl),&
!            setinfo(1+d*setl+1:1+(d+1)*setl),setl,basinfo,bas,EIJ,EGH,&
!            coefAB,coefCD,la,lb,lc,ld)
             DEALLOCATE(coefCD)

          END DO                                   !end loop over d
        END DO                                     !end loop over c

        !"loop" over g,h orbitals
        DO h=0,norb-1
          DO g=0,norb-1
            !"loop" over i,j orbitals
            DO j=0,norb-1
              DO i=0,norb-1
!               CALL inner loop on k,i,j
              END DO
            END DO
          END DO
        END DO

        DEALLOCATE(coefAB)

      END DO                                       !end loop over b
    END DO                                         !end loop over a

    WRITE(42) XX 
    CLOSE(unit=42,status='keep')

    !Optional Printing
    IF (options(7) .GE. 3) THEN 
      WRITE(*,*) "Two electron Integrals by Orbital" 
      WRITE(*,997) "I   J   G   H    Value"   
      DO i=0,norb-1
        DO j=0,norb-1
          DO g=0,norb-1
            DO h=0,norb-1
              WRITE(*,996) i,j,g,h,XX(i,j,g,h)
            END DO
          END DO
        END DO
      END DO
    END IF

    DEALLOCATE(XX)
    DEALLOCATE(II)
    DEALLOCATE(Ck)
    DEALLOCATE(Ckp)
    DEALLOCATE(Dk)
    DEALLOCATE(Dkp)
    DEALLOCATE(Ok)
    DEALLOCATE(Okp)

    !reassign memory
    fmem = fmem + (norb**(4))*8.0/1.0E6
    fmem = fmem + (2*setK)*norb*norb*8.0D0/1.0E6
    CALL nmem(fmem)

    CALL CPU_TIME(timeF) 
    WRITE(*,*) 
    WRITE(*,999) "Two electron integrals constructed in (s)", (timeF-timeS)

  END SUBROUTINE proc2e

!---------------------------------------------------------------------
!		Generates intermediate II for electron repulsion	
!---------------------------------------------------------------------
  SUBROUTINE clmII(II,p,q,la,lb,lc,ld,PQ,Dk,Dkp,Ck,Ckp,kmaxAB,kmaxCD,setc,setd,setl)
             
    IMPLICIT NONE
    REAL(KIND=8),PARAMETER :: Pi = 3.1415926535897931
    !Values
    ! II	: 3D dp, values of intermediate
    ! p,q	: dp, alpha p and q, exponential coefficients
    ! PQ	: 1D dp, P-Q array
    ! la,b,c,d	: 1D int, max angular QN of set a,b,c,d
    ! Dk,Dkp	: 1D dp, nonzero coefficients of setAB,CD
    ! Ck,Ckp	: 1D int, nonzero N,L,M of setAB,CD 
    ! setc,setd	: 1D int, setinfo c,d 
    ! setl	: int, length of set
    ! Fj	: 1D dp, Boys function table
    ! Rtab	: 4D dp, recursive auxil. table
    ! Nmax	: int, n + nbar + n' + n'bar
    ! foo	: int, dummy

    !Inout
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: II
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: PQ, Dk,Dkp
    INTEGER, DIMENSION(0:), INTENT(IN) :: Ck,Ckp,setc,setd 
    INTEGER, DIMENSION(0:), INTENT(IN) :: la,lb,lc,ld
    REAL(KIND=8), INTENT(IN) :: p,q
    INTEGER, INTENT(IN) :: kmaxAB,kmaxCD,setl

    !Internal
    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: Rtab
    LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: Rbol
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Fj
    REAL(KIND=8), DIMENSION(0:2) :: nc,nd
    REAL(KIND=8) :: ll,TT
    INTEGER :: N,Np,L,Lp,M,Mp,Nmax,Lmax,Mmax,orbc,orbd,foo
    INTEGER :: g,h,i,j,k,kp,z

    Nmax = la(0) + lb(0) + lc(0) + ld(0)
    Lmax = la(1) + lb(1) + lc(1) + ld(1)
    Mmax = la(2) + lb(2) + lc(2) + ld(2)

    ALLOCATE(Fj(0:Nmax+Lmax+Mmax))
    ALLOCATE(Rtab(-2:Nmax,-2:Lmax,-2:Mmax,0:Nmax+Lmax+Mmax))
    ALLOCATE(Rbol(-2:Nmax,-2:Lmax,-2:Mmax,0:Nmax+Lmax+Mmax))
 
    !new overlap values
    ll = 2*Pi**(5.0D0/2.0D0)/(p*q*(p+q)**(0.5D0)) 
    TT = p*q*(PQ(0)**2.0D0 + PQ(1)**2.0D0 + PQ(2)**2.0D0)/(p+q) 

    !get Boys Table
    Fj(:) = (/ (0.0D0, i=0, Nmax+Lmax+Mmax) /) 
    CALL Boys(Fj,Nmax+Lmax+Mmax,TT)

    !setup recursive table
    DO N=-2,Nmax
      DO L=-2,Lmax
        DO M=-2,Mmax
          DO z=0,Nmax+Lmax+Mmax
            Rtab(N,L,M,z) = 0.0D0
            Rbol(N,L,M,z) = .FALSE.
          END DO
        END DO
      END DO
    END DO

    !setup the intermediate array
    DO g=0,setc(0)-1
      orbc = setc(3+g)

      DO h=0,setd(0)-1
        orbd = setd(3+h)

        !loop over k'
        DO kp=0,kmaxCD
          foo = Ckp(kp)
          Np = foo/300
          foo = foo - Np*300
          Lp = foo/20
          Mp = foo - Lp*20
         
          DO k=0,kmaxAB
            foo = Ck(k)
            N = foo/300
            foo = foo - N*300
            L = foo/20
            M = foo - L*20
            
           CALL RNLMj(PQ(0),PQ(1),PQ(2),N+Np,L+Lp,M+Mp,0,p+q,Fj,Rtab,Rbol)

           II(k,orbc,orbd) = II(k,orbc,orbd) + (-1)**(Np+Lp+Mp)*&
           ll*Dkp(kp)*Rtab(N+Np,L+Lp,M+Mp,0)

          END DO                                        !end loop over k
        END DO                                          !end loop over kp
      END DO                                            !end loop over h
    END DO                                              !end loop over g

    DEALLOCATE(Fj)
    DEALLOCATE(Rtab)
    DEALLOCATE(Rbol)

  END SUBROUTINE clmII
!===================================================================
!                       FUNCTIONS

END PROGRAM int2e
