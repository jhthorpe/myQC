!///////////////////////////////////////////////////////////////////
!//            Parses input for myQC 
!//
!//              James H Thorpe, in Group of John Stanton
!//              The University of Florida
!//
!///////////////////////////////////////////////////////////////////

PROGRAM parser
  USE env

  IMPLICIT NONE

  REAL(KIND=8), PARAMETER :: A2b=(1.8897161646320724D0)

  !WORK NOTE : hardcoded, needs fixing

  ! atoms 	: 1D int, list of atoms
  ! radii	: 1D dp, list of radii (A)
  ! options	: 1D in, list of calculation options

  !options currently:
  ! 0) geometry type    : 1 - atom, 2 - diatomic, 3 - linear, 4- internal, 5- cartesian 
  ! 1) calculation type : 0 - scf 
  ! 2) basis set        : 0 - STO-3G
  ! 3) referecnce       : 0 - RHF, 1 - UHF, 2 - ROHF
  ! 4) parallel alg     : 0 - None, 1 - OMP, 2 - MPI
  ! 5) number of procs  : number
  ! 6) memory           : number in MB
  ! 7) verbosity        : 0 - none, 1 - some, 2 - all, 3 - wtf

  ! Variables
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE:: radii
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xyz 
  INTEGER,DIMENSION(:),ALLOCATABLE :: atoms
  INTEGER(KIND=8),DIMENSION(0:7) :: options 

  !internal variables
  CHARACTER(LEN=20),DIMENSION(0:1) :: line
  CHARACTER(LEN=20) :: str
  CHARACTER(LEN=2) :: id
  REAL(KIND=8) :: val, timeS, timeF
  INTEGER :: i,nnuc
  LOGICAL :: flag

  CALL CPU_TIME(timeS)

  i = 0

  WRITE(*,*) "parse called"
  WRITE(*,*) "=================================="
  WRITE(*,*) "Input parameters"

  OPEN (unit=1,file="input",status="old",access="sequential")
  
  !WORK NOTE : hardcoded, needs fixing
  READ(1,*) str       !read what type of molecule we have
  options(0) = getsys(str)

  !~~~~~~~~~~~~~~
  ! Read the atoms/geometry
  ! atom read
  IF (options(0) .EQ. 1) THEN
    ALLOCATE(atoms(0:0))
    ALLOCATE(radii(0:0))
    ALLOCATE(xyz(0:0,0:0))
    READ(1,*) str
    atoms(0) = getelem(str) 
    READ(1,*) 

  !Diatomic read
  ELSE IF (options(0) .EQ. 2) THEN
    !get atoms
    ALLOCATE(atoms(0:1))
    ALLOCATE(radii(0:0))
    ALLOCATE(xyz(0:1,0:2))
    DO i=0,1
      READ(1,*) str
      atoms(i) = getelem(str)
    END DO
    READ(1,*)
    !get radii
    READ(1,*) radii(0)
    READ(1,*)

  ! Cartesian Read
  ELSE IF (options(0) .EQ. 5) THEN
    flag = .TRUE.
    nnuc = 0

    ! get number of atoms
    DO WHILE (flag) 
      READ(1,*) str
      IF (str .EQ. 'END') THEN
        flag = .FALSE.
      ELSE
        nnuc = nnuc + 1
      END IF
    END DO
    
    IF (nnuc .LE. 0 ) STOP "No atoms in system"
    ALLOCATE(atoms(0:nnuc-1))
    ALLOCATE(xyz(0:nnuc-1,0:2))
    ALLOCATE(radii(0:0))

    REWIND(1)
    READ(1,*)

    !construct molecule
    DO i=0,nnuc-1
      READ(1,*) id, xyz(i,:)
      WRITE(*,*) id, xyz(i,:)
      atoms(i) = getelem(id)
    END DO
 
    !check it isn't crap
    IF (.NOT. checkgeom(xyz,nnuc)) STOP 'atoms too close'
    READ(1,*)
   
  !Not implimented yet 
  ELSE
    WRITE(*,*) "Sorry, that system type not implimented yet. Exiting parser."
    STOP 'bad system'

  END IF 

  ! WORK NOTE : this is currently hardcoded for dev. improve later
  !~~~~~~~~~~~~~~
  ! Read the options 
  READ(1,*) line
  options(1) = getcalc(line(1)) 
  READ(1,*) line
  options(2) = getbasis(line(1))
  READ(1,*) line
  options(3) = getref(line(1))
  READ(1,*) line
  options(4) = getpar(line(1))
  READ(1,*) line
  options(5) = getnode(line(1)) 
  READ(1,*) line
  options(6) = getmem(line(1))
  READ(1,*) line
  options(7) = getverb(line(1)) 
  !~~~~~~~~~~~~~~
    
  IF (options(0) .NE. 5) THEN
    WRITE(*,*) "Atom array" 
    WRITE(*,*) atoms
    WRITE(*,*) "Radii array"
    WRITE(*,*) radii
  END IF


  CALL build(atoms,radii,options,xyz)

  CLOSE (unit=1, status="keep")
  WRITE(*,*) "=================================="

  IF (options(0) .NE. 1) DEALLOCATE(radii)
  DEALLOCATE(atoms)
  DEALLOCATE(xyz)

  CALL CPU_TIME(timeF)

  WRITE(*,*) "parse ran in (s) : ", (timeF - timeS)

  CONTAINS

!~~~~~~~~~~
  ! Function to get the elemental value of text input
  INTEGER FUNCTION getelem(chr)
    IMPLICIT NONE
    CHARACTER(LEN=2),INTENT(IN) :: chr
    ! WORK NOTES : ugly - impliment dictionary?
    IF (options(0) .NE. 5) WRITE(*,*) chr
    getelem = 1
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
    getsys = 3
    IF (chr .EQ. 'ATOM') THEN
      getsys = 1
      WRITE(*,*) "System Type :  ATOM" 
    ELSE IF (chr .EQ. 'DIATOMIC') THEN
      getsys = 2
      WRITE(*,*) "System Type :  DIATOMIC" 
    ELSE IF (chr .EQ. 'LINEAR') THEN
      getsys = 3
      WRITE(*,*) "System Type :  LINEAR" 
    ELSE IF (chr .EQ. 'MOLECULE') THEN
      getsys = 4
      WRITE(*,*) "System Type :  INTERNAL" 
    ELSE IF (chr .EQ. 'CARTESIAN') THEN
      getsys = 5
      WRITE(*,*) "System Type : CARTESIAN"
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
    getcalc = 1
    IF (chr .EQ. 'SCF' .OR. chr .EQ. 'HF') THEN
      getcalc = 2
      WRITE(*,*) "Method : SCF"
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
      WRITE(*,*) "Basis : STO-3G"
    ELSE IF (chr .EQ. 'tester1') THEN
      getbasis = 1
      WRITE(*,*) "Basis : tester1"
    ELSE IF (chr .EQ. 'tester2') THEN
      getbasis = 2
      WRITE(*,*) "basis : tester2"
    ELSE IF (chr .EQ. 'tester3') THEN
      getbasis = 3
      WRITE(*,*) "basis : tester3"
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
    getref = 2
    IF (chr .EQ. 'RHF') THEN
      getref = 1
      WRITE(*,*) "Ref : RHF"
    ELSE IF (chr .EQ. 'UHF') THEN
      getref = 2
      WRITE(*,*) "Ref : UHF"
      WRITE(*,*) "Sorry, this reference is not implimented yet. Exiting..."
      STOP 'bad ref'
    ELSE IF (chr .EQ. 'ROHF') THEN
      getref = 3
      WRITE(*,*) "Ref : ROHF"
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
      WRITE(*,*) "Parallel : OMP"
    ELSE IF (chr .EQ. 'MPI') THEN
      getpar = 3
      WRITE(*,*) "Parallel : MPI"
    ELSE
      getpar = 1
      WRITE(*,*) "Parallel : none"
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
    !WRITE(*,*) val 
    !WRITE(*,*) "Memory per CPU (MB) : ", FLOOR(val*1.0D0*1000.0D0)
    !getmem = FLOOR(val*1.0D0*1000.0D0)
    WRITE(*,*) "Memory per CPU (MB) : ", val
    getmem = val
  END FUNCTION getmem

!~~~~~~~~~~
  ! Function to return the calculation option value
  INTEGER FUNCTION getverb(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20),INTENT(IN) :: chr
    ! WORK NOTE : ugly - impliment dictionary?
    getverb = 1
    IF (chr .EQ. '1') THEN
      getverb = 1
      WRITE(*,*) "Verbosity : 1"
    ELSE IF (chr .EQ. '0') THEN
      getverb = 0
      WRITE(*,*) "Verbosity : 0"
    END IF
    END FUNCTION getverb
!~~~~~~~~~~
! Subroutine that actually build the xyz and radii files 
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
    INTEGER :: i,j,nnuc

    nnuc = SIZE(atoms)

    ALLOCATE(r(0:nnuc-1))
    COM = [0.0D0, 0.0D0, 0.0D0] 
    mass = [1.000, 4.000, 7.000, 9.000, 11.000, 12.000, 14.000, 16.000, 19.000, 20.000]

    !write xyz file
    OPEN(unit=1,file='nucpos',status='replace',access='sequential') 
    OPEN(unit=2,file='radii',status='replace',access='sequential')
    OPEN(unit=3,file='envdat',status='replace',access='sequential')
    OPEN(unit=4,file='fmem',status='replace',access='sequential')

    !WORK NOTE - currently hardcoded, only works for linear molecules
    IF (options(0) .EQ. 1) THEN
      DO i=1,nnuc-1 
        xyz(i,:) = [0.0D0, 0.0D0, radii(i-1)*A2b]
      END DO
    END IF

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
      xyz(i,0) = (xyz(i,0) - COM(0))*A2B 
      xyz(i,1) = (xyz(i,1) - COM(1))*A2B 
      xyz(i,2) = (xyz(i,2) - COM(2))*A2B 
      WRITE(1,*) atoms(i), xyz(i,:)
    END DO

    !write radii file
    DO i=0,nnuc-1 
      r = (/ (SQRT((xyz(i,0) - xyz(j,0))**2.0D0 &
       + (xyz(i,1) - xyz(j,1))**2.0D0 &
       + (xyz(i,2) - xyz(j,2))**2.0D0 ), j=0,nnuc-1) /)  
      WRITE(2,*) r*A2b
    END DO

    !WORK NOTE - hardcoded
    !write to enviroment file
    WRITE(3,*) nnuc        !number of nuclei
    WRITE(3,*) SUM(atoms)  !number of electrons, currently assuming closed shell 
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
    CLOSE(unit=2)
    CLOSE(unit=3)
    CLOSE(unit=4)

    DEALLOCATE(r)

  END SUBROUTINE build
  
!~~~~~~~~~~
  ! function that returns false if the geom is good
  LOGICAL FUNCTION checkgeom(xyz,nnuc)
    IMPLICIT NONE
 
    ! inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: xyz
    INTEGER, INTENT(IN) :: nnuc

    ! internal
    REAL(KIND=8) :: r, tol
    INTEGER :: i,j

    checkgeom = .TRUE.

    DO i=0, nnuc-1
      DO j=i+1,nnuc-1
        r = (xyz(i,0)-xyz(j,0))**2.0D0 
        r = r + (xyz(i,1)-xyz(j,1))**2.0D0
        r = r + (xyz(i,2)-xyz(j,2))**2.0D0
        IF (SQRT(r) .LT. 0.20D0) THEN
          WRITE(*,*) "These atoms are too close (r < 0.2 A) :", i,j
          checkgeom = .FALSE.
        END IF
      END DO
    END DO    

  END FUNCTION checkgeom
!~~~~~~~~~~

END PROGRAM parser
