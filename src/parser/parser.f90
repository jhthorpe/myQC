!///////////////////////////////////////////////////////////////////
!//            Parses input for myQC 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//
!///////////////////////////////////////////////////////////////////

PROGRAM parser

  IMPLICIT NONE

  !angstrom to borh conversion
  REAL(KIND=8), PARAMETER :: A2B=(1.8897161646320724D0)

  !WORK NOTE : hardcoded, needs fixing

  ! atoms 	: 1D int, list of atoms
  ! radii	: 1D dp, list of radii (A)
  ! options	: 1D in, list of calculation options

  !options currently:
  ! 0) geometry type    : 0- internal, 1- cartesian 
  ! 1) calculation type : 0-SCF, 1-MP2
  ! 2) basis set        : 0 - STO-3G
  ! 3) referecnce       : 0 - RHF, 1 - UHF, 2 - ROHF
  ! 4) parallel alg     : 0 - None, 1 - OMP, 2 - MPI
  ! 5) number of procs  : number
  ! 6) memory           : number in MB
  ! 7) verbosity        : 0 - none, 1 - some, 2 - all, 3 - wtf
  ! 8) SCF_Conv		: 10^-x, x:[0-12]
  ! 9) charge		: int
  !10) multiplicity	: int
  !11) units		: 0-angstrom,1-bohr
  !12) ao2mo alg        : 0-fast,1-slow
  !13) excitation       : 0-none,1-cis
  !14) root algorithm   : 0-lanczos
  !15) number ex. states: number
  !16) prop level	: 0-none, 1-first order, 2-second order

  ! Variables
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: radii
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xyz 
  INTEGER,DIMENSION(:),ALLOCATABLE :: atoms
  INTEGER(KIND=8),DIMENSION(0:16) :: options 

  !internal variables
  CHARACTER(LEN=20),DIMENSION(0:1) :: line
  CHARACTER(LEN=20) :: str
  CHARACTER(LEN=2) :: id
  REAL(KIND=8) :: val
  INTEGER :: i,nnuc,fline,nline
  LOGICAL :: flag

999 FORMAT(1x,A19,F8.5)

  i = 0

  !defaults
  options = [0,0,0,0,0,1,1000,1,7,0,1,0,1,0,0,0,0]

  CALL EXECUTE_COMMAND_LINE('cat ZMAT')

  WRITE(*,*) ""
  WRITE(*,*) "Input parameters"

  !get number of lines in file
  CALL getfline(flag,fline) 
  IF (flag) STOP "parser encountered an error"

  OPEN (unit=1,file="ZMAT",status="old",access="sequential")
  
  !WORK NOTE : hardcoded, needs fixing
  READ(1,*) str       !read what type of molecule we have
  options(0) = getsys(str)

  !read in geometry
  IF (options(0) .EQ. 0) THEN
    CALL cartesian(fline,nline,nnuc,atoms,xyz)
  ELSE
    WRITE(*,*) "Sorry, that input style not supported yet"
    CALL EXECUTE_COMMAND_LINE('touch error')
    STOP
  END IF

  !read in options
  IF (nline .LT. fline) CALL read_options(fline-nline-1,options)
  CLOSE(unit=1)

  !build molecule
  CALL build(atoms,radii,options,xyz)

  !check no contradicting options
  flag = check_options(options)
  IF (flag) THEN
    WRITE(*,*) "==========================================="
    CALL EXECUTE_COMMAND_LINE('touch error')
    STOP "Bad options in ZMAT"
  END IF

  WRITE(*,*) "==========================================="

  !if requested, print options
  CALL print_options(options)

  DEALLOCATE(atoms)
  DEALLOCATE(xyz)

  CONTAINS

!~~~~~~~~~~
  ! Function to get the elemental value of text input
  INTEGER FUNCTION getelem(chr)
    IMPLICIT NONE
    CHARACTER(LEN=2),INTENT(IN) :: chr
    ! WORK NOTES : ugly - impliment dictionary?
    getelem = -1
    IF (chr .EQ. 'H') THEN
      getelem = 1
    ELSE IF (chr .EQ. 'He') THEN
      getelem = 2
    ELSE IF (chr .EQ. 'Li') THEN
      getelem = 3
    ELSE IF (chr .EQ. 'Be') THEN
      getelem = 4
    ELSE IF (chr .EQ. 'B') THEN
      getelem = 5
    ELSE IF (chr .EQ. 'C') THEN
      getelem = 6
    ELSE IF (chr .EQ. 'N') THEN
      getelem = 7
    ELSE IF (chr .EQ. 'O') THEN
      getelem = 8
    ELSE IF (chr .EQ. 'F') THEN
      getelem = 9
    ELSE IF (chr .EQ. 'Ne') THEN
      getelem = 10
    ELSE
      WRITE(*,*) "Sorry, that element is not implemented yet. Exiting ..."
      getelem = -1
      STOP 'bad element'
    END IF
  END FUNCTION getelem

!~~~~~~~~~~
  ! Function to return the system option value
  INTEGER FUNCTION getsys(chr)
    IMPLICIT NONE
    CHARACTER(LEN=9),INTENT(IN) :: chr
    ! WORK NOTE : ugly - impliment a dictionary?
    getsys = 0
    IF (chr .EQ. 'CARTESIAN') THEN
      getsys = 0
  !    WRITE(*,*) "System Type : CARTESIAN"
    ELSE IF (chr .EQ. 'INTERNAL') THEN
      getsys = 1
   !   WRITE(*,*) "System Type : INTERNAL"
    ELSE
      WRITE(*,*) "Bad system type input. Exiting..."
      STOP 'bad system'
    END IF    
  END FUNCTION getsys

!~~~~~~~~~~
  ! Function to return the calculation option value
  INTEGER FUNCTION getcalc(chr)
    IMPLICIT NONE
    CHARACTER(LEN=10),INTENT(IN) :: chr
    ! WORK NOTE : ugly - impliment dictionary?
    IF (chr .EQ. 'SCF' .OR. chr .EQ. 'HF') THEN
      getcalc = 0
   !   WRITE(*,*) "Method : SCF"
    ELSE IF (chr .EQ. 'MP2') THEN
      getcalc = 1
    ELSE IF (chr .EQ. 'CIS') THEN
      getcalc=2
    ELSE
      WRITE(*,*) "Sorry, that method has not been implimented. Exiting..."
      STOP 'bad method'
    END IF
  END FUNCTION getcalc

!~~~~~~~~~~
  ! Function to return the basis set option value
  INTEGER FUNCTION getbasis(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20),INTENT(IN) :: chr
    ! WORK NOTE : ugly - impliment dictionary?
    getbasis = 0
    IF (chr .EQ. 'STO-3G') THEN
      getbasis = 0
  !    WRITE(*,*) "Basis : STO-3G"
    ELSE IF (chr .EQ. 'tester1') THEN
      getbasis = 1
  !    WRITE(*,*) "Basis : tester1"
    ELSE IF (chr .EQ. 'tester2') THEN
      getbasis = 2
  !    WRITE(*,*) "basis : tester2"
    ELSE IF (chr .EQ. 'tester3') THEN
      getbasis = 3
  !    WRITE(*,*) "basis : tester3"
    ELSE
      WRITE(*,*) "Sorry, that basis has not been implimented. Exiting..."
      STOP 'bad basis'
    END IF
  END FUNCTION getbasis

!~~~~~~~~~~
  ! Function to return the basis set option value
  INTEGER FUNCTION getref(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20),INTENT(IN) :: chr
    ! WORK NOTE : ugly - impliment dictionary?
    getref = 0
    IF (chr .EQ. 'RHF') THEN
      getref = 0
  !    WRITE(*,*) "Ref : RHF"
    ELSE IF (chr .EQ. 'UHF') THEN
      getref = 1
  !    WRITE(*,*) "Ref : UHF"
    ELSE IF (chr .EQ. 'ROHF') THEN
      getref = 2
  !    WRITE(*,*) "Ref : ROHF"
      WRITE(*,*) "Sorry, this reference is not implimented yet. Exiting..."
      STOP 'bad ref'
    ELSE
      WRITE(*,*) "Sorry, that reference has not been implimented. Exiting..."
      STOP 'bad ref'
    END IF
  END FUNCTION getref

!~~~~~~~~~~
  ! Function to return the parallel option value
  INTEGER FUNCTION getpar(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20),INTENT(IN) :: chr
    getpar=1
    IF (chr .EQ. 'OMP') THEN
      getpar = 2
  !    WRITE(*,*) "Parallel : OMP"
    ELSE IF (chr .EQ. 'MPI') THEN
      getpar = 3
  !    WRITE(*,*) "Parallel : MPI"
    ELSE
      getpar = 1
  !    WRITE(*,*) "Parallel : none"
    END IF
  END FUNCTION getpar

!~~~~~~~~~~
  ! Function to set the number of nodes
  INTEGER FUNCTION getnode(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20),INTENT(IN) :: chr
    WRITE(*,*) "Nodes : not implimented yet"
    getnode = 1
  END FUNCTION getnode

!~~~~~~~~~~
  ! Function to get how much memory to use
  INTEGER FUNCTION getmem(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20), INTENT(IN) :: chr
    INTEGER :: val
    READ (chr,'(I8)') val 
    IF (val .LT. 0) THEN
      WRITE(*,*) "You specificed a memory less than zero."
      getmem = 1000
    ELSE
 !     WRITE(*,*) "Memory per CPU (MB) : ", val
      getmem = val
    END IF
  END FUNCTION getmem

!~~~~~~~~~~
  ! Function to get SCF_CONV criteria 
  INTEGER FUNCTION getSCF_Conv(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20), INTENT(IN) :: chr
    INTEGER :: val
    READ (chr,'(I8)') val 
    IF (val .GT. 11 .OR. val .LT. 0) THEN
  !    WRITE(*,*) "SCF Convergence (default) :", 7
    ELSE IF( val .EQ. 11) THEN
  !    WRITE(*,*) "SCF Convergence (these go to...) : ", val
    ELSE
  !    WRITE(*,*) "SCF Convergence : ", val
    END IF
    getSCF_CONV = val
  END FUNCTION getSCF_Conv

!~~~~~~~~~~
  ! Function to get charge 
  INTEGER FUNCTION getcharge(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20), INTENT(IN) :: chr
    INTEGER :: val
    READ (chr,'(I8)') val 
 !   WRITE(*,*) "Charge : ", val
    getcharge = val
  END FUNCTION getcharge

!~~~~~~~~~~
  ! Function to get charge 
  INTEGER FUNCTION getmulti(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20), INTENT(IN) :: chr
    INTEGER :: val
    READ (chr,'(I8)') val 
    IF (val .GT. 0) THEN
 !     WRITE(*,*) "Multiplicity : ", val
    ELSE
      WRITE(*,*) "bad value for multiplicity, exiting."
      CALL EXECUTE_COMMAND_LINE('touch error')
      STOP
    END IF
    getmulti = val
  END FUNCTION getmulti

!---------------------------------------------------------------------
!			Get verbosity level	
!---------------------------------------------------------------------
  INTEGER FUNCTION getverb(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20),INTENT(IN) :: chr
    ! WORK NOTE : ugly - impliment dictionary?
    getverb = 0
    IF (chr .EQ. '1') THEN
      getverb = 1
  !    WRITE(*,*) "Verbosity : 1 (some)"
    ELSE IF (chr .EQ. '2') THEN
      getverb = 2
  !    WRITE(*,*) "Verbosity : 2 (detailed)"
    ELSE IF (chr .EQ. '3') THEN
      getverb = 3
  !    WRITE(*,*) "Verbosity : 3 (wtf)"
    ELSE
  !    WRITE(*,*) "Verbosity : 0 (none)"
    END IF
    END FUNCTION getverb

!---------------------------------------------------------------------
!			Get units
!---------------------------------------------------------------------
  INTEGER FUNCTION getunits(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20),INTENT(IN) :: chr
    IF (chr .EQ. 'Bohr') THEN
      getunits = 1
  !    WRITE(*,*) "Units : Bohr"
    ELSE
      getunits = 0
  !    WRITE(*,*) "Units : Angstrom"
    END IF
  END FUNCTION getunits
!---------------------------------------------------------------------
!			Get ao2mo algorithm	
!---------------------------------------------------------------------
  INTEGER FUNCTION getao2mo(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20),INTENT(IN) :: chr
    IF (chr .EQ. '1') THEN
 !     WRITE(*,*) "ao2mo alg : 1 (slow)"
      getao2mo=1
    ELSE
      WRITE(*,*) "Sorry, only slow ao2mo supported"
      getao2mo=1
      !WRITE(*,*) "ao2mo alg : 0 (fast)"
      !getao2mo=0
    END IF
  END FUNCTION getao2mo
!---------------------------------------------------------------------
!                       Get excite
!---------------------------------------------------------------------
  INTEGER FUNCTION getexcite(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20),INTENT(IN) :: chr
    IF (chr .EQ. 'CIS' .OR. chr .EQ. "1") THEN
      getexcite=1
    ELSE IF (chr .EQ. 'NONE' .OR. chr .EQ. "0") THEN
      getexcite=0
    ELSE
      WRITE(*,*) "Sorry, only CIS supported"
      getexcite=0
    END IF
  END FUNCTION getexcite
!---------------------------------------------------------------------
!                       Get root_alg
!---------------------------------------------------------------------
  INTEGER FUNCTION getroot_alg(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20),INTENT(IN) :: chr
    IF (chr .EQ. 'lanczos' .OR. chr .EQ. "0") THEN
      getroot_alg=0
    ELSE
      WRITE(*,*) "Sorry, only Lanczos root_alg supported"
      getroot_alg=0
    END IF
  END FUNCTION getroot_alg
!---------------------------------------------------------------------
!                       Get E_NUM
!---------------------------------------------------------------------
  INTEGER FUNCTION gete_num(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20), INTENT(IN) :: chr
    INTEGER :: val
    READ (chr,'(I8)') val
    IF (val .LT. 0) THEN
      WRITE(*,*) "You specificed a number of states less than zero."
      gete_num = 1
    ELSE
      gete_num = val
    END IF
  END FUNCTION gete_num
!---------------------------------------------------------------------
!                       Get PROP
!---------------------------------------------------------------------
  INTEGER FUNCTION get_prop(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20), INTENT(IN) :: chr
    IF (chr .EQ. 'NONE' .OR. chr .EQ. "ZERO" .OR. chr .EQ. "0") THEN 
      get_prop = 0
    ELSE IF (chr .EQ. 'FIRST' .OR. chr .EQ. "1") THEN  
      get_prop = 1
    ELSE IF (chr .EQ. 'SECOND' .OR. chr .EQ. "2") THEN
      get_prop = 2
    ELSE
      WRITE(*,*) "Sorry, only 0,1,2 order properties are coded"
      get_prop = 0 
    END IF
  END FUNCTION get_prop
!---------------------------------------------------------------------
!		Build the molecule (nucpos and envdat files)	
!---------------------------------------------------------------------
  SUBROUTINE build(atoms, radii, options, xyz)
    IMPLICIT NONE

    ! Inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: xyz 
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: radii
    INTEGER(KIND=8), DIMENSION(0:) :: options
    INTEGER, DIMENSION(0:) :: atoms

    ! Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: r 
    REAL(KIND=8), DIMENSION(0:2) :: COM
    REAL(KIND=8), DIMENSION(1:10) :: mass
    REAL(KIND=8) :: temp
    INTEGER :: i,j,nnuc,nelc,nelcA,nelcB,charge,unpr

    nnuc = SIZE(atoms)

    IF (.NOT. checkgeom(xyz,nnuc,options(11))) THEN 
      CALL EXECUTE_COMMAND_LINE('touch error')
    END IF

    ALLOCATE(r(0:nnuc-1))
    COM = [0.0D0, 0.0D0, 0.0D0] 
    mass = [1.000, 4.000, 7.000, 9.000, 11.000, 12.000, 14.000, 16.000, 19.000, 20.000]

    !write xyz file
    OPEN(unit=1,file='nucpos',status='replace',access='sequential') 
!    OPEN(unit=2,file='radii',status='replace',access='sequential')
    OPEN(unit=3,file='envdat',status='replace',access='sequential')
    OPEN(unit=4,file='fmem',status='replace',access='sequential')

    !get COM 
    DO i=0,nnuc-1
     COM(0) = COM(0) + mass(atoms(i))*xyz(i,0) 
     COM(1) = COM(1) + mass(atoms(i))*xyz(i,1)
     COM(2) = COM(2) + mass(atoms(i))*xyz(i,2) 
     temp = temp + mass(atoms(i)) 
    END DO

    COM(:) = COM(:) / temp

    !recenter
    DO i=0,nnuc-1
      xyz(i,0) = (xyz(i,0) - COM(0))
      xyz(i,1) = (xyz(i,1) - COM(1)) 
      xyz(i,2) = (xyz(i,2) - COM(2)) 
    END DO

    !adjust for input units
    IF (options(11) .EQ. 0) THEN
      xyz(:,:) = xyz(:,:)*A2B
    END IF

    !write to nucpos file
    DO i=0,nnuc-1
      WRITE(1,*) atoms(i), xyz(i,:)
    END DO

    !write radii file
!    DO i=0,nnuc-1 
!      r = (/ (SQRT((xyz(i,0) - xyz(j,0))**2.0D0 &
!       + (xyz(i,1) - xyz(j,1))**2.0D0 &
!       + (xyz(i,2) - xyz(j,2))**2.0D0 ), j=0,nnuc-1) /)  
!      WRITE(2,*) r*A2b
!    END DO

    !charge
    charge = options(9)
    unpr = options(10)-1
    nelc = SUM(atoms)-charge
    nelcA = (nelc-unpr)/2
    nelcB = (nelc-unpr)/2
    
    !multiplicity
    nelcA = nelcA + unpr

    !checking
    IF (nelcA + nelcB .NE. nelc) THEN
      WRITE(*,*) "That charge and multiplicity is not allowed."
      CALL EXECUTE_COMMAND_LINE('touch error')
    END IF
    IF (nelc .LT. 0 .OR. nelcA .LT. 0 .OR. nelcB .LT. 0) THEN
      WRITE(*,*) "You have less than 0 electrons. :)"
      CALL EXECUTE_COMMAND_LINE('touch error')
    END IF 
    IF (nelcA .NE. nelcB .AND. options(3) .EQ. 0) THEN
      WRITE(*,*) "You cannot use RHF for open shell molecules!"
      CALL EXECUTE_COMMAND_LINE('touch error')
    END IF
    
    !write to enviroment file
    WRITE(3,*) nnuc        !number of nuclei
    WRITE(3,*) nelcA,nelcB
    WRITE(3,*) SIZE(options)
    WRITE(3,*) options
    WRITE(3,*) 
    WRITE(3,*) "#number of nuclei" 
    WRITE(3,*) "#number of electrons"
    WRITE(3,*) "#length of options array"
    WRITE(3,*) "options array"

    !write to free memory file
    WRITE(4,*) options(6) 

    CLOSE(unit=1)
!    CLOSE(unit=2)
    CLOSE(unit=3)
    CLOSE(unit=4)

    DEALLOCATE(r)

  END SUBROUTINE build
  
!---------------------------------------------------------------------
!			Check if atoms are too close	
! 			Returns false if the geom is good
!---------------------------------------------------------------------
  ! function that returns false if the geom is good
  LOGICAL FUNCTION checkgeom(xyz,nnuc,units)
    IMPLICIT NONE
    REAL(KIND=8), PARAMETER :: A2B=(1.8897161646320724D0)
 
    ! inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
    INTEGER(KIND=8), INTENT(IN) :: units
    INTEGER, INTENT(IN) :: nnuc

    ! internal
    REAL(KIND=8) :: r, tol
    INTEGER :: i,j
    LOGICAL :: flag

    flag = .TRUE.

    DO i=0, nnuc-1
      DO j=i+1,nnuc-1
        r = (xyz(i,0)-xyz(j,0))**2.0D0 
        r = r + (xyz(i,1)-xyz(j,1))**2.0D0
        r = r + (xyz(i,2)-xyz(j,2))**2.0D0
        IF (units == 0 .AND. SQRT(r) .LT. 0.20D0) THEN
          WRITE(*,*) "These atoms are too close (r < 0.2 A) :", i,j
          flag = .FALSE.
        ELSE IF (units == 1 .AND. SQRT(r) .LT. 0.2D0*A2B) THEN
          WRITE(*,*) "These atroms are too close (r < 0.2A) :", i,j 
          flag = .FALSE. 
        END IF
      END DO
    END DO    

    checkgeom = flag

  END FUNCTION checkgeom

!---------------------------------------------------------------------
!		Check if file exists and get number of lines	
!---------------------------------------------------------------------
  SUBROUTINE getfline(flag,fline)
    IMPLICIT NONE
    !Inout
    INTEGER, INTENT(INOUT) :: fline
    LOGICAL, INTENT(INOUT) :: flag
    !Internal
    INTEGER :: io
    LOGICAL :: ex
    flag = .FALSE. 
    INQUIRE(file='ZMAT',EXIST=ex)
    IF (.NOT. ex) THEN
      WRITE(*,*) "You need to create the input file : 'ZMAT'"
      CALL EXECUTE_COMMAND_LINE('touch error')
      flag = .TRUE. 
      fline = -1
      RETURN 
    END IF
    fline = 0
    io = 0
    OPEN(unit=1,file='ZMAT',status='old',access='sequential')
    DO WHILE (io .EQ. 0)
      READ(1,*,iostat=io)
      IF (io .EQ. 0) fline = fline + 1
    END DO
    CLOSE(unit=1)
  END SUBROUTINE getfline

!---------------------------------------------------------------------
!			Read in cartesian style input		
!---------------------------------------------------------------------
  SUBROUTINE cartesian(fline,nline,nnuc,atoms,xyz)
    IMPLICIT NONE

    !Inout
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: xyz
    INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: atoms
    INTEGER, INTENT(INOUT) :: nline, nnuc
    INTEGER, INTENT(IN) :: fline

    !Internal
    CHARACTER(LEN=2) :: id
    INTEGER :: i
    LOGICAL :: flag
    
    flag = .TRUE.
    nnuc = 0
    i = 0
    nline = 0

    !get number of nuclei
    DO WHILE (flag) 
      i = i + 1
      IF (i .GT. fline) THEN
        WRITE(*,*) "You need to put 'END' marker in ZMAT"
        CALL EXECUTE_COMMAND_LINE('touch error')
        RETURN 
      END IF
      READ(1,*) str
      IF (str .EQ. 'END') THEN
        flag = .FALSE.
      ELSE
        nnuc = nnuc + 1
      END IF
    END DO
   
    !allocate memory 
    IF (nnuc .LE. 0 ) STOP "No atoms in system"
    ALLOCATE(atoms(0:nnuc-1))
    ALLOCATE(xyz(0:nnuc-1,0:2))
    REWIND(1)
    READ(1,*)
    nline = nline + 1

    !read in the molecule
    DO i=0,nnuc-1
      READ(1,*) id, xyz(i,:)
  !    WRITE(*,*) id, xyz(i,:)
      atoms(i) = getelem(id)
      nline = nline + 1
    END DO
 
    !extra lines
    nline = nline + 1
 !   WRITE(*,*)
    READ(1,*) 
 
  END SUBROUTINE cartesian

!---------------------------------------------------------------------
!			Read remaining options	
!---------------------------------------------------------------------
  SUBROUTINE read_options(rline,options)
    IMPLICIT NONE

    !Inout
    INTEGER(KIND=8),DIMENSION(0:), INTENT(INOUT) :: options
    INTEGER, INTENT(IN) :: rline

    !Internal
    CHARACTER(LEN=20),DIMENSION(0:1) :: line
    INTEGER :: i

    READ(1,*)

    DO i=0, rline-1
      READ(1,*) line
      IF (line(0) == 'CALC=') THEN
        options(1) = getcalc(line(1))
      ELSE IF (line(0) == 'BASIS=') THEN
        options(2) = getbasis(line(1))
      ELSE IF (line(0) == 'CHARGE=') THEN
        options(9) = getcharge(line(1))
      ELSE IF (line(0) == 'MULTI=') THEN
        options(10) = getmulti(line(1))
      ELSE IF (line(0) == 'REF=') THEN
        options(3) = getref(line(1))
      ELSE IF (line(0) == 'PAR=') THEN
        options(4) = getpar(line(1))
      ELSE IF (line(0) == 'NODES=') THEN
        options(5) = getnode(line(1))
      ELSE IF (line(0) == 'MEMORY=') THEN
        options(6) = getmem(line(1))
      ELSE IF (line(0) == 'VERB=') THEN
        options(7) = getverb(line(1))
      ELSE IF (line(0) == 'SCF_Conv=') THEN
        options(8) = getSCF_Conv(line(1))
      ELSE IF (line(0) == 'UNITS=') THEN
        options(11) = getunits(line(1))
      ELSE IF (line(0) == 'AO2MO=') THEN
        options(12) = getao2mo(line(1))
      ELSE IF (line(0) == 'EXCITE=') THEN
        options(13) = getexcite(line(1))
      ELSE IF (line(0) == 'ROOT_ALG=') THEN
        options(14) = getroot_alg(line(1))
      ELSE IF (line(0) == 'E_NUM=') THEN
        options(15) = gete_num(line(1))
      ELSE IF (line(0) == 'PROP=') THEN
        options(16) = get_prop(line(1))
      ELSE
        WRITE(*,*) "parser could not understand options line ", i
      END IF
    END DO

  END SUBROUTINE read_options

!---------------------------------------------------------------------
!			Print options	
!---------------------------------------------------------------------
  SUBROUTINE print_options(options)
    IMPLICIT NONE
    INTEGER(KIND=8), DIMENSION(0:), INTENT(IN) ::  options
999 FORMAT(1x,A22,I8)
    WRITE(*,*) "Options"
    WRITE(*,*) "Read type           : ",options(0)
    WRITE(*,*) "Calculation         : ",options(1)
    WRITE(*,*) "Basis               : ",options(2)
    WRITE(*,*) "Reference           : ",options(3)
    WRITE(*,*) "Parallel Algorithm  : ",options(4)
    WRITE(*,*) "Nodes               : ",options(5)
    WRITE(*,*) "Memory (MB)         : ",options(6)
    WRITE(*,*) "SCF Convergence     : ",options(8)
    WRITE(*,*) "Charge              : ",options(9)
    WRITE(*,*) "Multiplicity        : ",options(10)
    WRITE(*,*) "Units               : ",options(11)
    WRITE(*,*) "ao2mo               : ",options(12)
    WRITE(*,*) "excite              : ",options(13)
    WRITE(*,*) "root algorithm      : ",options(14)
    WRITE(*,*) "num excite          : ",options(15)
    WRITE(*,*) "property order      : ",options(16)
    WRITE(*,*) "==========================================="
  END SUBROUTINE print_options
!---------------------------------------------------------------------
!			Check options	
!---------------------------------------------------------------------
  LOGICAL FUNCTION check_options(options)
    IMPLICIT NONE
    INTEGER(KIND=8), DIMENSION(0:), INTENT(IN) :: options
    LOGICAL :: flag

    flag = .FALSE.

    !Excitation checks
    IF (options(13) .NE. 0) THEN
      IF (options(1) .NE. 0) THEN
        WRITE(*,*) "Sorry, only SCF CIS coded"
        flag = .TRUE.
      END IF
      IF (options(3) .NE. 1) THEN
        WRITE(*,*) "Sorry, only UHF CIS references are coded"
        flag = .TRUE.
      END IF
    END IF

    !property checks
    IF (options(16) .NE. 0) THEN
      IF (options(1) .NE. 0) THEN
        WRITE(*,*) "Sorry, only SCF properties are coded"
        flag = .TRUE.
      END IF
      IF (options(3) .NE. 0) THEN
        WRITE(*,*) "Sorry, only RHF reference properites are coded"
        flag = .TRUE.
      END IF
    END IF
    
    check_options = flag

  END FUNCTION check_options

!---------------------------------------------------------------------
END PROGRAM parser
