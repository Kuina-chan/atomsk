MODULE options
!
!**********************************************************************************
!*  OPTIONS                                                                       *
!**********************************************************************************
!* This module applies one or several transformations to the system.              *
!*                                                                                *
!* WHAT MUST BE PROVIDED:                                                         *
!* options_array                                                                  *
!*           if allocated, string array of unknown dimension containing the       *
!*           names and parameters of options to be applied to the system.         *
!*           If not allocated, then no option will be applied.                    *
!* H         3 x 3 real array containing vectors of the supercell. H(1,:) is      *
!*           the first vector, H(2,:) the second, H(3,:) the third.               *
!* P         N x 4 real array containing atom positions. N=number of atoms,       *
!*           each line contains the coordinates (X,Y,Z) and atomic number.        *
!* S         if allocated, N x 4 real array containing positions of shells        *
!*           (for ionic core/shell model), and number of shells N must be         *
!*           identical to N in P. If S is not allocated, then no shell exist.     *
!* AUXNAMES  if allocated, string array of dimension M, M being the number of     *
!*           existing auxiliary properties (must be equal to M in AUX). If        *
!*           AUXNAMES is not allocated, then AUX must also be non-allocated.      *
!* AUX       if allocated, N x M real array containing values of auxiliary        *
!*           properties of atoms, e.g. their velocity, forces, etc. The first     *
!*           dimension N must be equal to the number of atoms N in P, the second  *
!*           dimension M must be equal to the number of auxiliary properties.     *
!*           If AUX is not allocated then no auxiliary property exists, and       *
!*           AUXNAMES must also be non-allocated.                                 *
!*                                                                                *
!* WHAT IS RETURNED BY THIS ROUTINE:                                              *
!* modified arrays H, P, S, AUXNAMES and AUX.                                     *
!*                                                                                *
!**********************************************************************************
!* (C) October 2010 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 Nov. 2026                                     *
!**********************************************************************************
!* This program is free software: you can redistribute it and/or modify           *
!* it under the terms of the GNU General Public License as published by           *
!* the Free Software Foundation, either version 3 of the License, or              *
!* (at your option) any later version.                                            *
!*                                                                                *
!* This program is distributed in the hope that it will be useful,                *
!* but WITHOUT ANY WARRANTY; without even the implied warranty of                 *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  *
!* GNU General Public License for more details.                                   *
!*                                                                                *
!* You should have received a copy of the GNU General Public License              *
!* along with this program.  If not, see <http://www.gnu.org/licenses/>.          *
!**********************************************************************************
!
!Load modules
USE atoms
USE comv
USE constants
USE strings
USE crystallography
USE elasticity
USE exprev
USE messages
USE files
USE subroutines
USE guess_form
USE resize
!
!Modules managing options
USE addatom
USE addshells
USE alignx
USE bindshells
USE carttofrac
USE center
USE cell
USE crack
USE cut_cell
USE deform
!USE denoise
USE dislocation
USE disturb
USE duplicate
USE freeze
USE mirror
USE orient
USE orthocell
USE properties
USE reduce_cell
USE remdoubles
USE rmatom
USE rmprop
USE rmshells
USE roll
USE rotate
USE roundoff
USE select
USE separate
USE shift
USE sort
USE spacegroup
USE stress
USE substitute
USE swap
USE torsion
USE unit
USE unskew
USE velocity
USE wrap
! -- add other options in alphabetical order --
!
!
!
CONTAINS
!
SUBROUTINE OPTIONS_AFF(options_array,Huc,H,P,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
!
!Declare variables
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=32):: optionname
CHARACTER(LEN=4096):: temp, msg
CHARACTER(LEN=4096),DIMENSION(10):: treal !text containing a real number and maybe a word
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=32),DIMENSION(:),ALLOCATABLE::  substr
LOGICAL:: fileexists !does a file exist?
INTEGER:: i, ioptions, j, k
INTEGER:: status
INTEGER:: strlength
REAL(dp):: tempreal !temporary real number
REAL(dp),DIMENSION(3):: vector
REAL(dp),DIMENSION(3,3),INTENT(IN):: Huc !Base vectors of unit cell (MUST NOT be modified by options)
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: HS  !Copy of H for shells
REAL(dp),DIMENSION(3,3):: ORIENT    !crystalographic orientation
REAL(dp),DIMENSION(9,9):: C_tensor  !stiffness tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P  !atomic positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S  !shell positions (if any)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX, AUXdummy !auxiliary properties of atoms/shells
!
!Variables relative to Option: add-atom
CHARACTER(LEN=2):: addatom_species !species of atom(s) to add
CHARACTER(LEN=8):: addatom_type    !"at" or "near" or "random"
REAL(dp),DIMENSION(4):: addatom_prop  !properties of atom(s) to add
                                      !if addatom_type=="at", position x,y,z of new atom
                                      !if addatom_type=="relative", index of atom and x,y,z
                                      !if addatom_type=="near", index of atom
                                      !if addatom-type=="random", number of atoms to add
!
!Variables relative to Option: crack
CHARACTER(LEN=1):: crackline !direction of crack tip line, must be x, y or z
CHARACTER(LEN=1):: crackplane !normal to the plane of cut, must be x, y or z
CHARACTER(LEN=3):: crackmode !crack mode, must be I, II or III
CHARACTER(LEN=11):: cracktype !planestress or planestrain
REAL(dp):: crackK  !stress intensity factor
REAL(dp):: mu   !shear modulus of the material
!
!Variables for option: center
INTEGER:: center_atom  !index of atom that must be centered; if 0, center the center of mass
!
!Variables for option: cell
CHARACTER(LEN=5):: celldir  !direction in which the cell is modified: x,y,z,xy,yx,xz,zx,yz,zy,xyz,all
CHARACTER(LEN=5):: cellop   !operation to perform on the cell: add,rm,set,auto (or "rebox")
REAL(dp):: celllength  !value by which the cell will be modified
!
!Variables for option: create shells
CHARACTER(LEN=2):: cs_species
!
!Variables relative to Option: cut
CHARACTER(LEN=5):: cut_dir   !above or below
CHARACTER(LEN=16):: cutdir   !x, y, z, or crystallographic direction
REAL(dp):: cutdistance
!
!Variables relative to Option: deform
CHARACTER(LEN=2):: def_dir       !direction of applied strain (X, Y or Z)
REAL(dp):: def_strain, def_poisson  !applied strain and Poisson's ratio
!
!Variables relative to Option: denoise
!INTEGER:: denoise_Niter  !max. number of iterations
!REAL(dp):: denoise_dmax  !max. displacement under which denoise is stopped
!
!Variables relative to Option: dislocation
CHARACTER(LEN=16):: dislocplane   !(x, y, z, or Miller vector)
CHARACTER(LEN=4096):: dislocline    !(x, y, z, or Miller vector, or file name)
CHARACTER(LEN=4096):: disloctype     !edge, edge_add, edge_rm, screw, mixed, loop
REAL(dp):: nu
REAL(dp),DIMENSION(5):: pos !pos(1:3) = position of dislocation; pos(4) = radius of disloc.loop
REAL(dp),DIMENSION(3):: b !Burgers vector
!
!Variables relative to Option: disturb
REAL(dp),DIMENSION(3):: dist_dmax  !maximum translation vector along each Cartesian direction
!
!Variables relative to Option: duplicate
INTEGER, DIMENSION(3):: dupmatrix  !number of times the system must be
                                   !repeated in each direction of space
!
!Variables relative to Option: freeze
CHARACTER(LEN=5):: freeze_coord    !coordinate(s) to freeze (x,y,z,xy,xz,yz,xyz)
CHARACTER(LEN=5):: freeze_side     !'above' or 'below'
CHARACTER(LEN=16):: freeze_normal  !direction normal to plane (x,y,z, or [hkl])
REAL(dp):: freeze_dist             !distance to plane
!
!Variables relative to Option: mirror
CHARACTER(LEN=16):: mirror_dir   !x, y, z, or crystallographic direction
REAL(dp):: mirror_d  !distance between mirror plane and origin of coordinates
!
!Variables relative to Option: orient
CHARACTER(LEN=32),DIMENSION(3):: oldvec, newvec !old and new Miller indices
!
!Variables relative to Option: properties
CHARACTER(LEN=128):: propfile  !file containing the properties
REAL(dp),DIMENSION(3):: lat_a0   !lattice constants a, b, c
REAL(dp),DIMENSION(9,9):: C_tensor_dummy !elastic tensor (dummy)
!
!Variables relative to Option: reduce-cell
CHARACTER(LEN=1):: reduce_dir  !direction along which the cell will be reduced (may be empty)
!
!Variables relative to Option: rmatom
CHARACTER(LEN=32):: rmatom_prop  !what atom(s) to erase
!
!Variables relative to Option: remove doubles
REAL(dp):: rmd_radius            !in angstroms
!
!Variables relative to Option: rmprop
CHARACTER(LEN=32):: rmprop_prop  !property that must be erased
!
!Variables relative to Option: rmshells
CHARACTER(LEN=6):: rmshells_prop  !species on which shells are removed
!
!Variables relative to Option: rotate
CHARACTER(LEN=4096):: rot_axis    !Cartesian x, y or z, Miller vector, or 3 real coordinates
CHARACTER(LEN=4096):: rot_file    !Cartesian x, y or z, Miller vector, or 3 real coordinates
REAL(dp):: rot_angle              !in degrees
!
!Variables relative to Option: select
CHARACTER(LEN=16):: select_multiple !'union' or 'subtract' or 'intersect'
CHARACTER(LEN=32):: region_dir     !x, y, z
CHARACTER(LEN=128):: region_side   !'in' or 'out' or 'all' or 'property'
CHARACTER(LEN=4096):: region_geom  !geometry of the region: "box" or "sphere"
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
REAL(dp),DIMENSION(3):: region_1 !First corner for'box', or center of sphere
REAL(dp),DIMENSION(3):: region_2 !Last corner for'box', or radius of sphere
!
!Variables relative to Option: separate
REAL(dp):: sep_radius  !atoms closer than this distance will be pulled apart
REAL(dp):: sep_dist    !atoms will be pulled apart by this distance
!
!Variables relative to Option: shear
CHARACTER(LEN=1):: shear_dir  !direction of applied shear strain (X, Y or Z)
CHARACTER(LEN=1):: shear_surf !normal to the surface to be sheared (X, Y or Z)
REAL(dp):: shear_strain  !applied shear strain in %
!
!Variables relative to Option: shift
CHARACTER(LEN=5):: shift_dir       !'above' or 'below'
CHARACTER(LEN=16):: shift_axis     !x, y, z or crystallographic direction
REAL(dp):: shift_dist       !distance of the "cut plane" to the origin
REAL(dp):: shift_tau1, shift_tau2, shift_tau3  !shift vector
!
!Variables relative to Option: sort
CHARACTER(LEN=8):: sortorder
CHARACTER(LEN=16):: sortcol
!
!Variables relative to Option: stress
CHARACTER(LEN=4096):: stress_in  !Voigt stress component, or name of file
REAL(dp):: stress                !value of applied stress (GPa)
!
!Variables relative to Option: substitute
CHARACTER(LEN=3):: sp1, sp2
!
!Variables relative to Option: swap
CHARACTER(LEN=16),DIMENSION(2):: swap_id  !Cartesian axes or indices of atoms to swap
!
!Variables relative to Option: unit
CHARACTER(LEN=128):: unit1, unit2
!
!Variables relative to Option: velocity
REAL(dp):: vel_T  !target temperature for Maxwell-Boltzmann distribution
!
! -- add other options variables in alphabetical order --
!
j=0
!
!
WRITE(msg,*) "===================================================="
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
WRITE(msg,*) "===   NOW APPLYING OPTIONS IN SEQUENTIAL ORDER   ==="
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Initialize variables
lat_a0(:) = 1.d0
!If orientation of the system is unknown, default orientation is assumed to be [100] [010] [001]
IF( NINT(VECLENGTH(ORIENT(1,:))) == 0 .OR.          &
  & NINT(VECLENGTH(ORIENT(2,:))) == 0 .OR.          &
  & NINT(VECLENGTH(ORIENT(3,:))) == 0       ) THEN
  ORIENT(:,:) = Id_Matrix(:,:)
ENDIF
!By default the elastic tensor is not set
 !C_tensor(:,:) = 0.d0
 C_tensor_dummy(:,:) = 0.d0
!
!If there is no option then exit this module
IF(.NOT.ALLOCATED(options_array)) GOTO 1000
IF( verbosity==4 ) THEN
  WRITE(temp,*) SIZE(options_array)
  WRITE(msg,*) "Number of options: "//TRIM(ADJUSTL(temp))
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  IF( SIZE(options_array)>0 ) THEN
    DO i=1,SIZE(options_array)
      WRITE(temp,*) i
      WRITE(msg,*) "options_array("//TRIM(ADJUSTL(temp))//") : "//TRIM(ADJUSTL(options_array(i)))
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
  ENDIF
ENDIF
IF(SIZE(options_array)<=0) THEN
  WRITE(msg,*) "options_array has zero size, skipping"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  GOTO 1000
ENDIF
!
IF( .NOT.ALLOCATED(P) .OR. SIZE(P)<=0 ) THEN
  !The user wants to apply an option, but no system exists in memory
  !=> error
  nerr=nerr+1
  WRITE(msg,*) "ERROR array P is unallocated or has zero size"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  CALL ATOMSK_MSG(2814,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Verify the order of options
DO ioptions=1,SIZE(options_array)-1
  READ(options_array(ioptions),*,END=100,ERR=100) temp
  READ(options_array(ioptions+1),*,END=100,ERR=100) msg
  IF( msg=="-orthogonal-cell" .OR. msg=="-orthorhombic-cell" .OR. &
    & msg=="-orthogonal-box" .OR. msg=="-orthorhombic-box" .OR. msg=="-orthobox" ) THEN
    msg="-orthocell"
  ENDIF
  IF( temp=="-duplicate" .AND. (msg=="-orthocell" .OR. msg=="-wrap")  ) THEN
    !Advise to use other options BEFORE "-duplicate"
    CALL ATOMSK_MSG(2600,(/temp,options_array(ioptions+1)/),(/0.d0/))
  ELSEIF( temp=="-unskew" .AND. msg=="-alignx"  ) THEN
    !Advise to use other option "-alignx"  BEFORE "-unskew"
    CALL ATOMSK_MSG(2600,(/temp,msg/),(/0.d0/))
  ENDIF
ENDDO
!
!
!
100 CONTINUE
!Loop on all options
!NOTE: each option must be applied to shells and auxiliary properties at the same time
!     as they are applied to atoms! All modules that change the number
!     of atoms or change their order in the array P must also add, remove or
!     re-order the shells in S and/or auxiliary properties in AUX accordingly.
!
!
DO ioptions=1,SIZE(options_array)
  !Initialisations
  fileexists = .FALSE.
  treal(:) = ""
  disloctype = ""
  dislocline = ""
  dislocplane = ""
  select_multiple = ""
  region_dir = ""
  region_geom = ""
  region_side = ""
  status = 0
  b(:) = 0.d0
  pos(:) = 0.d0
  region_1(:) = 0.d0
  region_2(:) = 0.d0
  !
  IF( verbosity==4 ) THEN
    WRITE(msg,*) "===================================================="
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(temp,*) ioptions
    msg = " ===   OPTION # "//TRIM(ADJUSTL(temp))//" : "//TRIM(options_array(ioptions))
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  !
  !READ(options_array(ioptions),*,END=1000,ERR=1000) optionname
  optionname = options_array(ioptions)
  j = SCAN(optionname," ")
  optionname = TRIM(optionname(1:j))
  optionname = TRIM(ADJUSTL(optionname))
  IF(optionname=='') EXIT
  HS=H
  !Make sure that the dummy arrays are empty or deallocated
  IF(ALLOCATED(AUXdummy)) DEALLOCATE(AUXdummy)
  C_tensor_dummy(:,:) = 0.d0
  !
  SELECT CASE(StrDnCase(optionname))
  !
  CASE('')
    !empty entry => just ignore it
    CONTINUE
    !
  CASE('-add-atom','-add-atoms','-addatom','-addatoms')
    !READ(options_array(ioptions),*,END=800,ERR=800) optionname, addatom_species, addatom_type
    k=0
    CALL STRCUT(options_array(ioptions),' ',substr)
    addatom_species = TRIM(substr(2))
    addatom_type = TRIM(substr(3))
    IF( addatom_type=="at" .OR. addatom_type=="AT" .OR. addatom_type=="@" ) THEN
      !"-add-atom species at x y z"
      !Read coordinates x, y, z of atom to add
      !Can be real numbers or fractional coordinates with the keyword "box"
      !READ(options_array(ioptions),*,END=800,ERR=800) optionname, addatom_species, &
      !  & addatom_type, treal(1), treal(2), treal(3)
      IF( SIZE(substr)<6 ) GOTO 800
      treal(1) = TRIM(substr(4))
      treal(2) = TRIM(substr(5))
      treal(3) = TRIM(substr(6))
      CALL BOX2DBLE( H(:,1) , treal(1) , addatom_prop(1) , status )
      IF(status>0) THEN
        temp = treal(1)
        GOTO 810
      ENDIF
      CALL BOX2DBLE( H(:,2) , treal(2) , addatom_prop(2) , status )
      IF(status>0) THEN
        temp = treal(2)
        GOTO 810
      ENDIF
      CALL BOX2DBLE( H(:,3) , treal(3) , addatom_prop(3) , status )
      IF(status>0) THEN
        temp = treal(3)
        GOTO 810
      ENDIF
    ELSEIF( addatom_type=="relative" .OR. addatom_type=="rel" ) THEN
      !"-add-atom species relative id"
      !Read index of atom near which the new atom must be added
      !and relative coordinates x, y, z of atom to add
      !Can be real numbers or fractional coordinates with the keyword "box"
      !READ(options_array(ioptions),*,END=800,ERR=800) optionname, addatom_species, &
      !  & addatom_type, j, treal(1), treal(2), treal(3)
      IF( SIZE(substr)<7 ) GOTO 800
      READ(substr(4),*,END=800,ERR=800) addatom_prop(1)
      treal(1) = TRIM(substr(5))
      treal(2) = TRIM(substr(6))
      treal(3) = TRIM(substr(7))
      CALL BOX2DBLE( H(:,1) , treal(1) , addatom_prop(2) , status )
      IF(status>0) THEN
        temp = treal(1)
        GOTO 810
      ENDIF
      CALL BOX2DBLE( H(:,2) , treal(2) , addatom_prop(3) , status )
      IF(status>0) THEN
        temp = treal(2)
        GOTO 810
      ENDIF
      CALL BOX2DBLE( H(:,3) , treal(3) , addatom_prop(4) , status )
      IF(status>0) THEN
        temp = treal(3)
        GOTO 810
      ENDIF
    ELSEIF( addatom_type=="near" .OR. addatom_type=="NEAR" .OR.      &
          & addatom_type=="random" .OR. addatom_type=="RANDOM" .OR.  &
          & addatom_type=="rand" .OR. addatom_type=="RAND"           ) THEN
      !"-add-atom species near id"
      !"-add-atom species random N"
      !Read index of atom near which the new atom must be added,
      !or the number of atoms to insert randomly
      !READ(options_array(ioptions),*,END=800,ERR=800) optionname, addatom_species, &
      !  & addatom_type, addatom_prop(1)
      READ(substr(4),*,END=800,ERR=800) addatom_prop(1)
    ELSE
      nerr=nerr+1
      GOTO 1000
    ENDIF
    CALL ADDATOM_XYZ(H,P,S,AUXNAMES,AUX,addatom_species,addatom_type,addatom_prop,SELECT)
  !
  CASE('-add-shells','-addshells','-as')
    !"-add-shells species"
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp
    IF(temp=='all') THEN
      cs_species = 'XX'
    ELSE
      READ(temp,*,ERR=800,END=800) cs_species
    ENDIF
    CALL ADDSHELLS_XYZ(P,S,cs_species)
  !
  CASE('-alignx')
    CALL ALIGN_X(H,P,S,C_tensor)
  !
  CASE('-bind-shells','-bs')
    CALL BSHELLS_XYZ(H,P,S,AUXNAMES,AUX,SELECT)
  !
  CASE('-cell','-box','change_cell','change_box')
    !"-cell <X|Y|Z> <add|rm|set> <d>"
    !READ(options_array(ioptions),*,END=800,ERR=800) optionname, celldir, cellop, treal(1)
    CALL STRCUT(options_array(ioptions),' ',substr)
    IF( SIZE(substr)<4 ) GOTO 800
    celldir = TRIM(substr(2))
    cellop = TRIM(substr(3))
    treal(1) = TRIM(substr(4))
    SELECT CASE(StrDnCase(celldir))
    CASE('a','x',"xx","xy","xz")
      CALL BOX2DBLE( H(:,1) , treal(1) , celllength , status )
    CASE('b','y','yy',"yx","yz")
      CALL BOX2DBLE( H(:,2) , treal(1) , celllength , status )
    CASE('c','z',"zz","zx","zy")
      CALL BOX2DBLE( H(:,3) , treal(1) , celllength , status )
    CASE DEFAULT
      j = SCAN(treal(1),'°')
      IF(j>0) treal(1)(j:j) = ' '
      CALL BOX2DBLE( (/1.d0,0.d0,0.d0/) , treal(1) , celllength , status )
    END SELECT
    IF(status>0) THEN
      temp = treal(1)
      GOTO 810
    ENDIF
    CALL CELL_XYZ(H,celldir,cellop,celllength)
  !
  CASE('-center')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, center_atom
    CALL CENTER_XYZ(H,P,S,center_atom,SELECT)
  !
  CASE('-crack')
    !"-crack <I|II|III> <stress|strain> <Ks> <p1> <p2> <line> <normal> <mu> <nu>"
    !READ(options_array(ioptions),*,END=800,ERR=800) optionname, crackmode, &
    !    & cracktype, crackK, treal(1), treal(2), crackline, crackplane, mu, nu
    CALL STRCUT(options_array(ioptions),' ',substr)
    IF( SIZE(substr)<10 ) GOTO 800
    crackmode = TRIM(substr(2))
    cracktype = TRIM(substr(3))
    READ(substr(4),*,END=800,ERR=800) crackK
    treal(1) = TRIM(substr(5))
    treal(2) = TRIM(substr(6))
    crackline = TRIM(substr(7))
    crackplane = TRIM(substr(8))
    READ(substr(9),*,END=800,ERR=800) mu
    READ(substr(10),*,END=800,ERR=800) nu
    CALL BOX2DBLE( H(:,1) , treal(1) , pos(1) , status )
    IF(status>0) THEN
      temp = treal(1)
      GOTO 810
    ENDIF
    CALL BOX2DBLE( H(:,2) , treal(2) , pos(2) , status )
    IF(status>0) THEN
      temp = treal(2)
      GOTO 810
    ENDIF
    CALL CRACK_XYZ(H,P,S,crackmode,cracktype,crackK,crackline,crackplane,mu,nu,pos(1),pos(2),SELECT,AUXNAMES,AUX)
  !
  CASE('-cut')
    !"-cut <above|below> <d> <normal>"
    cut_dir = TRIM(ADJUSTL(options_array(ioptions)(5:)))
    IF( cut_dir(1:5)=='above' .OR. cut_dir(1:5)=='below' ) THEN
      !READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
      !    & cut_dir, treal(1), cutdir
      CALL STRCUT(options_array(ioptions),' ',substr)
      IF( SIZE(substr)<4 ) GOTO 800
      cut_dir = TRIM(substr(2))
      treal(1) = TRIM(substr(3))
      cutdir = TRIM(substr(4))
      !Check if numbers contain a keyword like "BOX" or "INF"
      i=1
      SELECT CASE(cutdir)
      CASE('x','X')
        CALL BOX2DBLE( H(:,1) , treal(1) , cutdistance , status )
      CASE('y','Y')
        CALL BOX2DBLE( H(:,2) , treal(1) , cutdistance , status )
      CASE('z','Z')
        CALL BOX2DBLE( H(:,3) , treal(1) , cutdistance , status )
      CASE DEFAULT
        !Should be Miller indices
        CALL MILLER2VEC(H,cutdir,ORIENT,addatom_prop(:),j)
        IF(j>0) THEN
          !Other error, unable to convert this string into a proper vector
          CALL ATOMSK_MSG(817,(/TRIM(cutdir)/),(/0.d0/))
          GOTO 1000
        ENDIF
        !Check that vector is not [000]
        IF( VECLENGTH(addatom_prop)<1.d-12 ) THEN
          CALL ATOMSK_MSG(814,(/""/),(/0.d0/))
          nerr=nerr+1
          GOTO 1000
        ENDIF
        vector(:) = addatom_prop(1)*H(1,:) + addatom_prop(2)*H(2,:) + addatom_prop(2)*H(3,:)
        CALL BOX2DBLE( vector(:) , treal(1) , cutdistance , status )
      END SELECT
      IF(status>0) THEN
        temp = treal(1)
        GOTO 810
      ENDIF
    ELSE
      cut_dir=""
    ENDIF
    CALL CUTCELL(H,P,S,AUX,cut_dir,cutdistance,cutdir,ORIENT,SELECT)
  !
  CASE('-def', '-deform')
    !"-deform <x|y|z> <def> [<nu>]"
    !"-deform <xy|xz|yz|yx|zx|zy> <def>"
    j=0
    def_poisson = 0.d0
    !READ(options_array(ioptions),*,END=800,ERR=800) optionname, def_dir, temp
    CALL STRCUT(options_array(ioptions),' ',substr)
    IF( SIZE(substr)<3 ) GOTO 800
    def_dir = TRIM(substr(2))
    temp = TRIM(substr(3))
    IF( SIZE(substr)>=4 ) THEN
      !User provided a Poisson ratio
      READ(substr(4),*,END=800,ERR=800) def_poisson
    ENDIF
    IF( temp=="untilt" .OR. temp=="unshear" ) THEN
      !Set the shear strain to untilt the box
      SELECT CASE(StrDnCase(def_dir))
      CASE("yx")
        def_strain = -1.d0*H(1,2)/H(1,1)
      CASE("zx")
        def_strain = -1.d0*H(1,3)/H(1,1)
      CASE("xy")
        def_strain = -1.d0*H(2,1)/H(2,2)
      CASE("zy")
        def_strain = -1.d0*H(2,3)/H(2,2)
      CASE("xz")
        def_strain = -1.d0*H(3,1)/H(3,3)
      CASE("yz")
        def_strain = -1.d0*H(3,2)/H(3,3)
      END SELECT
      !
    ELSE
      !"temp" should contain a number: read it
      treal(1) = TRIM(temp)
      !treal may contain the strain expressed in percent, e.g. "3%"
      SELECT CASE(StrDnCase(def_dir))
      CASE('x')
        i=1
      CASE('y')
        i=2
      CASE('z')
        i=3
      CASE("yx","zx")
        i=1
        j=1
      CASE("xy","zy")
        i=2
        j=1
      CASE("xz","yz")
        i=3
        j=1
      END SELECT
      CALL BOX2DBLE( H(:,i) , treal(1) , def_strain , status )
      IF(status>0) THEN
        temp = treal(1)
        GOTO 810
      ENDIF
    ENDIF
    IF( def_poisson.NE.0.d0 .AND. j==1 ) THEN
      !user specified a Poisson ratio, but asked for shear strain
      !=> display warning message that value of Poisson ratio will be ignored
      CALL ATOMSK_MSG(2766,(/""/),(/0.d0/))
    ENDIF
    102 CONTINUE
    CALL DEFORM_XYZ(H,P,S,def_dir,def_strain,def_poisson,SELECT)
  !
!  CASE('-denoise')
!     denoise_Niter = 0
!     denoise_dmax = 0.d0
!     !Read one value
!     READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp
!     !Detect if an integer or
!     temp = TRIM(ADJUSTL(temp))
!     IF( SCAN(temp,'.') > 0 ) THEN
!       !It is a real value
!       READ(temp,*,END=800,ERR=800) denoise_dmax
!     ELSE
!       !It is an integer value
!       READ(temp,*,END=800,ERR=800) denoise_Niter
!     ENDIF
!    CALL DENOISE_XYZ(H,P,S,SELECT)
  !
  CASE("-dislocation","-dislo","-disloc")
    !"-dislocation <p1> <p2> screw <line> <normal> <b>"
    !"-dislocation <p1> <p2> <edge|edge_add|edge_rm> <line> <normal> <b> <nu>"
    !"-dislocation <p1> <p2> mixed <line> <normal> <b1> <b2> <b3>"
    !"-dislocation loop <x> <y> <z> <normal> <radius> <bx> <by> <bz> <nu>"
    !"-dislocation file <file> <nu>"
    nu = 0.d0
    status=0
    !Read the first two keywords
    !READ(options_array(ioptions),*,END=800,ERR=800) optionname, treal(9), treal(10)
    CALL STRCUT(options_array(ioptions),' ',substr)
    !Depending on type, read further parameters
    IF( substr(2)=="loop" ) THEN
      disloctype = "loop"
      !Read coordinates of loop center, loop radius, components of Burgers vector, and Poisson ratio
      !READ(options_array(ioptions),*,END=800,ERR=800) optionname, treal(9), &
      !    & treal(1), treal(2), treal(3), dislocline, pos(4), b(1), b(2), b(3), nu
      IF( SIZE(substr)<11 ) GOTO 800
      treal(1) = TRIM(substr(3))
      treal(2) = TRIM(substr(4))
      treal(3) = TRIM(substr(5))
      dislocline = TRIM(substr(6))
      READ(substr(7),*,END=800,ERR=800) pos(4)
      READ(substr(8),*,END=800,ERR=800) b(1)
      READ(substr(9),*,END=800,ERR=800) b(2)
      READ(substr(10),*,END=800,ERR=800) b(3)
      READ(substr(11),*,END=800,ERR=800) nu
      !Convert the treal(:) into real numbers, if they contain the keyword "box"
      DO i=1,3
        CALL BOX2DBLE( H(:,i) , treal(i) , pos(i) , status )
        IF(status>0) THEN
          temp = treal(i)
          nerr=nerr+1
          GOTO 810
        ENDIF
      ENDDO
      !
    ELSEIF( substr(2)=="array" .OR. substr(2)=="file" ) THEN
      !Save file name into dislocline
      !READ(options_array(ioptions),*,END=800,ERR=800) optionname, treal(9), treal(10), nu
      IF( SIZE(substr)<4 ) GOTO 800
      dislocline = TRIM(substr(3))
      READ(substr(4),*,END=800,ERR=800) nu
    ELSE
      !Screw, edge, or mixed dislocation: the first 2 fields must be the dislocation position
      IF( SIZE(substr)<7 ) GOTO 800
      treal(9) = TRIM(substr(2))
      treal(10) = TRIM(substr(3))
      disloctype = TRIM(substr(4))
      !READ(options_array(ioptions),*,END=800,ERR=800) optionname, treal(9), treal(10), disloctype
      IF( disloctype=="mixed" ) THEN
        !Read the three components of Burgers vector
        !READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        !    & treal(9), treal(10), disloctype, dislocline, dislocplane, b(1), b(2), b(3)
        IF( SIZE(substr)<9 ) GOTO 800
        dislocline = TRIM(substr(5))
        dislocplane = TRIM(substr(6))
        READ(substr(7),*,END=800,ERR=800) b(1)
        READ(substr(8),*,END=800,ERR=800) b(2)
        READ(substr(9),*,END=800,ERR=800) b(3)
      ELSE
        !Only the norm of the Burgers vector is given
        !Construct the complete Burgers vector b(:)
        b(:) = 0.d0
        IF(disloctype=="screw") THEN
          !don't attempt to read nu
          !READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
          !    & treal(9), treal(10), disloctype, dislocline, dislocplane, tempreal
          IF( SIZE(substr)<7 ) GOTO 800
          dislocline = TRIM(substr(5))
          dislocplane = TRIM(substr(6))
          READ(substr(7),*,END=800,ERR=800) tempreal
          !only one component is given, the one along dislocline
          SELECT CASE(dislocline)
          CASE("x","X")
            b(1) = tempreal
          CASE("y","Y")
            b(2) = tempreal
          CASE("z","Z")
            b(3) = tempreal
          END SELECT
        ELSEIF(disloctype(1:4)=="edge") THEN
          !(it may be edge, edge_add or edge_rm)
          !Edge character => read nu
          !READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
          !    & treal(9), treal(10), disloctype, dislocline, dislocplane, tempreal, nu
          IF( SIZE(substr)<8 ) GOTO 800
          dislocline = TRIM(substr(5))
          dislocplane = TRIM(substr(6))
          READ(substr(7),*,END=800,ERR=800) tempreal
          READ(substr(8),*,END=800,ERR=800) nu
          !only one component is given, normal to dislocline and plane of cut
          SELECT CASE(dislocline)
          CASE("x","X")
            SELECT CASE(dislocplane)
            CASE('z','Z')
              b(2) = tempreal
            CASE('y','Y')
              b(3) = tempreal
            END SELECT
          CASE("y","Y")
            SELECT CASE(dislocplane)
            CASE('x','X')
              b(3) = tempreal
            CASE('z','Z')
              b(1) = tempreal
            END SELECT
          CASE("z","Z")
            SELECT CASE(dislocplane)
            CASE('x','X')
              b(2) = tempreal
            CASE('y','Y')
              b(1) = tempreal
            END SELECT
          END SELECT
        ELSE
          !Error: unknown dislocation type
          GOTO 800
        ENDIF
      ENDIF  !en if "mixed, "screw"
      !Check if numbers contain a keyword like "BOX" or "INF"
      SELECT CASE(dislocplane)
      CASE('x','X')
        j=1
      CASE('y','Y')
        j=2
      CASE('z','Z')
        j=3
      END SELECT
      SELECT CASE(dislocline)
      CASE('x','X')
        treal(2) = treal(9)
        treal(3) = treal(10)
        IF(j==2) THEN
          i=3
        ELSEIF(j==3) THEN
          i=2
        ENDIF
      CASE('y','Y')
        treal(3) = treal(9)
        treal(1) = treal(10)
        IF(j==1) THEN
          i=3
        ELSEIF(j==3) THEN
          i=1
        ENDIF
      CASE('z','Z')
        treal(1) = treal(9)
        treal(2) = treal(10)
        IF(j==1) THEN
          i=2
        ELSEIF(j==2) THEN
          i=1
        ENDIF
      END SELECT
      !Convert strings containing "box" with appropriate box dimension
      CALL BOX2DBLE( H(:,i) , treal(i) , pos(1) , status )
      CALL BOX2DBLE( H(:,j) , treal(j) , pos(2) , status )
    ENDIF
    CALL DISLOC_XYZ(H,P,S,disloctype,dislocline,dislocplane,b,nu,pos,SELECT,ORIENT,AUXNAMES,AUX,C_tensor)
  !
  CASE('-disturb')
    !Read one value
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp
    !Detect if more values follow
    strlength = LEN_TRIM(temp) + 10
    temp = TRIM(ADJUSTL(options_array(ioptions)(strlength:)))
    IF( LEN_TRIM(temp) > 0 ) THEN
      !There seem to be more values: try to read them
      j=0
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, dist_dmax(1), dist_dmax(2), dist_dmax(3)
    ELSE
      !There seems to be only one value : read it
      j=1
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, dist_dmax(1)
    ENDIF
    CALL DISTURB_XYZ(P,S,SELECT,j,dist_dmax)
  !
  CASE('-duplicate', '-dup', '-e', '-expand','-replicate')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        & dupmatrix(1), dupmatrix(2), dupmatrix(3)
    CALL DUPLICATECELL(H,P,S,dupmatrix,SELECT,AUX)
  !
  CASE('-fix','-freeze')
    !"-freeze <x|y|z|xy|xz|yz|xyz> [<above|below> <d> <x|y|z>]"
    !READ(options_array(ioptions),*,END=800,ERR=800) optionname, freeze_coord
    !strlength = LEN_TRIM(optionname) + LEN_TRIM(freeze_coord)
    !IF( LEN_TRIM(options_array(ioptions)) > strlength+2 ) THEN
    !  READ(options_array(ioptions),*,END=800,ERR=800) optionname, freeze_coord, freeze_side
    !ENDIF
    CALL STRCUT(options_array(ioptions),' ',substr)
    IF(SIZE(substr)<2) GOTO 800
    freeze_coord = TRIM(substr(2))
    IF( SIZE(substr)>=5 ) THEN
      !READ(options_array(ioptions),*,END=800,ERR=800) optionname, freeze_coord, freeze_side, treal(1), freeze_normal
      freeze_side = TRIM(substr(3))
      treal(1) = TRIM(substr(4))
      freeze_normal = TRIM(substr(5))
      !Check if treal(1) contains a keyword like "BOX" or "INF"
      SELECT CASE(StrDnCase(freeze_normal))
      CASE('x')
        CALL BOX2DBLE( H(:,1) , treal(1) , freeze_dist , status )
      CASE('y')
        CALL BOX2DBLE( H(:,2) , treal(1) , freeze_dist , status )
      CASE('z')
        CALL BOX2DBLE( H(:,3) , treal(1) , freeze_dist , status )
      CASE DEFAULT
        !Should be Miller indices
        CALL MILLER2VEC(H,freeze_normal,ORIENT,addatom_prop(:),j)
        IF(j>0) THEN
          !Other error, unable to convert this string into a proper vector
          CALL ATOMSK_MSG(817,(/TRIM(freeze_normal)/),(/0.d0/))
          GOTO 1000
        ENDIF
        !Check that vector is not [000]
        IF( VECLENGTH(addatom_prop)<1.d-12 ) THEN
          CALL ATOMSK_MSG(814,(/""/),(/0.d0/))
          nerr=nerr+1
          GOTO 1000
        ENDIF
        vector(:) = addatom_prop(1)*H(1,:) + addatom_prop(2)*H(2,:) + addatom_prop(2)*H(3,:)
        CALL BOX2DBLE( vector(:) , treal(1) , freeze_dist , status )
      END SELECT
      IF(status>0) THEN
        temp = treal(1)
        GOTO 810
      ENDIF
    ELSE
      freeze_side = ""
      freeze_normal = ""
    ENDIF
    CALL FREEZE_XYZ(H,P,AUXNAMES,AUX,freeze_coord,freeze_side,freeze_dist,freeze_normal,ORIENT,SELECT)
  !
  CASE('-frac','-fractional')
    CALL CART2FRAC_XYZ(H,P,S)
  !
  CASE('-mirror')
    !"-mirror <d> <normal>"
    !READ(options_array(ioptions),*,END=800,ERR=800) optionname, treal(1), mirror_dir
    CALL STRCUT(options_array(ioptions),' ',substr)
    IF(SIZE(substr)<3) GOTO 800
    treal(1) = TRIM(substr(2))
    mirror_dir = TRIM(substr(3))
    !Check if numbers contain a keyword like "BOX" or "INF"
    i=1
    SELECT CASE(mirror_dir)
    CASE('x','X')
      CALL BOX2DBLE( H(:,1) , treal(1) , mirror_d , status )
    CASE('y','Y')
      CALL BOX2DBLE( H(:,1) , treal(1) , mirror_d , status )
    CASE('z','Z')
      CALL BOX2DBLE( H(:,1) , treal(1) , mirror_d , status )
    CASE DEFAULT
      !Should be Miller indices
      CALL MILLER2VEC(H,mirror_dir,ORIENT,addatom_prop(:),j)
      IF(j>0) THEN
        !Other error, unable to convert this string into a proper vector
        CALL ATOMSK_MSG(817,(/TRIM(mirror_dir)/),(/0.d0/))
        GOTO 1000
      ENDIF
      !Check that vector is not [000]
      IF( VECLENGTH(addatom_prop)<1.d-12 ) THEN
        CALL ATOMSK_MSG(814,(/""/),(/0.d0/))
        nerr=nerr+1
        GOTO 1000
      ENDIF
      vector(:) = addatom_prop(1)*H(1,:) + addatom_prop(2)*H(2,:) + addatom_prop(2)*H(3,:)
      CALL BOX2DBLE( vector(:) , treal(1) , mirror_d , status )
    END SELECT
    IF(status>0) THEN
      temp = treal(1)
      GOTO 810
    ENDIF
    CALL MIRROR_XYZ(P,S,AUXNAMES,AUX,mirror_dir,mirror_d,ORIENT,SELECT)
  !
  CASE('-orient')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, oldvec(1), oldvec(2),  &
        & oldvec(3), newvec(1), newvec(2), newvec(3)
    CALL ORIENT_XYZ(H,P,S,oldvec,newvec,SELECT,C_tensor)
  !
  CASE('-orthogonal-cell','-orthorhombic-cell','-orthogonal-box','-orthorhombic-box','-orthocell','-orthobox')
    CALL ORTHOCELL_XYZ(H,P,S,AUX,SELECT)
  !
  CASE('-prop','-properties','-property')
    READ(options_array(ioptions),'(a128)',END=800,ERR=800) temp
    READ(temp(6:),'(a128)') propfile
    propfile = TRIM(ADJUSTL(propfile))
    CALL READ_PROPERTIES(propfile,H,P,S,ORIENT,C_tensor,AUXNAMES,AUX,SELECT)
  !
  CASE('-rebox','-shrink-wrap')
    CALL ATOMSK_MSG(2153,(/""/),(/0.d0/))
    CALL DETERMINE_H(H,P)
  !
  CASE('-reduce-cell','-reducecell','-reduce-box','-reducebox')
    i = LEN_TRIM(optionname)
    IF( LEN_TRIM(options_array(ioptions)(i+1:)) > 0 ) THEN
      READ(options_array(ioptions)(i+1:),*,END=800,ERR=800) reduce_dir
    ELSE
      reduce_dir = ""
    ENDIF
    CALL REDUCECELL(H,P,S,AUX,reduce_dir,SELECT)
  !
  CASE('-remove-atom','-remove-atoms','-rmatom','-rmatoms','-delete-atoms','-delete_atoms')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, rmatom_prop
    CALL RMATOM_XYZ(P,S,AUX,rmatom_prop,SELECT)
  !
  CASE('-remove-doubles', '-remove-double', '-rmd')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        & rmd_radius
    CALL REMDOUBLES_XYZ(H,P,S,AUX,rmd_radius,SELECT)
  !
  CASE('-remove-property', '-rmprop')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, rmprop_prop
    CALL RMPROP_XYZ(AUXNAMES,AUX,rmprop_prop,SELECT)
  !
  CASE('-remove-shells','-remove-shell','-rmshells','-rmshell')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, rmshells_prop
    rmshells_prop = ADJUSTL(rmshells_prop)
    CALL RMSHELLS_XYZ(P,S,rmshells_prop,SELECT)
  !
  CASE('-roll')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, shift_dir, rot_angle, rot_axis
    CALL ROLL_XYZ(H,P,S,AUXNAMES,AUX,shift_dir,rot_angle,rot_axis,SELECT)
  !
  CASE('-rot', '-rotate')
    rot_file=""
    rot_axis=""
    rot_angle=0.d0
    j=SCAN(options_array(ioptions),' ')
    temp =  TRIM(ADJUSTL( options_array(ioptions)(j:) ))
    IF( temp(1:3)=="com" ) THEN
      j=INDEX(options_array(ioptions),'com') + 3
      temp = TRIM(ADJUSTL( options_array(ioptions)(j:) ))
      j=1
    ELSE
      j=0
    ENDIF
    READ(temp,*,END=800,ERR=800) msg
    msg = TRIM(ADJUSTL(msg))
    INQUIRE(FILE=msg,EXIST=fileexists)
    IF( fileexists ) THEN
      rot_file = TRIM(ADJUSTL(msg))
    ELSE
      IF( SCAN(msg,'xXyYzZ[]')>0 ) THEN
        READ(temp,*,END=800,ERR=800) rot_axis, rot_angle
      ELSEIF( SCAN(msg,'1234567890')>0 ) THEN
        READ(temp,*,END=800,ERR=800) treal(1), treal(2), treal(3), rot_angle
        rot_axis = TRIM(treal(1))//' '//TRIM(treal(2))//' '//TRIM(treal(3))
      ELSE
        !Inconsistent parameters
        GOTO 800
      ENDIF
    ENDIF
    CALL ROTATE_XYZ(H,P,S,AUXNAMES,AUX,j,rot_axis,rot_angle,rot_file,ORIENT,SELECT,C_tensor)
  !
  CASE('-roundoff','-round-off')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, tempreal
    CALL ROUNDOFF_XYZ(H,P,S,AUXNAMES,AUX,temp,tempreal,SELECT)
  !
  CASE('-select')
    CALL STRCUT(options_array(ioptions),' ',substr)
    select_multiple = TRIM(substr(2))
    !READ(options_array(ioptions),*,END=800,ERR=800) optionname, select_multiple
    !Check if current selection must be combined with previous selection
    k=1
    IF( select_multiple=="union" .OR. select_multiple=="UNION" .OR.           &
      & select_multiple=="add" .OR. select_multiple=="ADD" .OR.               &
      & select_multiple=="or" .OR. select_multiple=="OR"          ) THEN
      select_multiple = "add"
      region_side = TRIM(substr(3))
    ELSEIF( select_multiple=="subtract" .OR. select_multiple=="SUBTRACT" .OR.            &
          & select_multiple=="substract" .OR. select_multiple=="SUBSTRACT" .OR.          &
          & select_multiple=="delete" .OR. select_multiple=="DELETE" .OR.                &
          & select_multiple=="remove"   .OR. select_multiple=="REMOVE"   .OR.            &
          & select_multiple=="rm" .OR. select_multiple=="RM" .OR. select_multiple=="del" ) THEN
      select_multiple = "rm"
      region_side = TRIM(substr(3))
    ELSEIF( select_multiple=="intersect" .OR. select_multiple=="INTERSECT" .OR. &
          & select_multiple=="and" .OR. select_multiple=="AND"                 ) THEN
      select_multiple = "intersect"
      region_side = TRIM(substr(3))
    ELSEIF( select_multiple=="xor"  .OR. select_multiple=="XOR"  ) THEN
      select_multiple = "xor"
      region_side = TRIM(substr(3))
    ELSEIF( select_multiple=="among"  .OR. select_multiple=="AMONG"  ) THEN
      select_multiple = "among"
      region_side = TRIM(substr(3))
    ELSE
      !No keyword after "-select" => save following string in region_side
      select_multiple = ""
      region_side = TRIM(substr(2))
      k=0
    ENDIF
    !
    IF( region_side=="above" .OR. region_side=="below" ) THEN
      !"-select <above|below> <d> <dir>"
      !IF( LEN_TRIM(select_multiple)>0 ) THEN
      !  READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, region_side, treal(1), region_dir
      !ELSE
      !  READ(options_array(ioptions),*,END=800,ERR=800) optionname, region_side, treal(1), region_dir
      !ENDIF
      treal(1) = TRIM(substr(3+k))
      region_dir = TRIM(substr(4+k))
      i=1
      SELECT CASE(region_dir)
      CASE('x','X')
        CALL BOX2DBLE( H(:,1) , treal(1) , region_1(1) , status )
      CASE('y','Y')
        CALL BOX2DBLE( H(:,2) , treal(1) , region_1(1) , status )
      CASE('z','Z')
        CALL BOX2DBLE( H(:,2) , treal(1) , region_1(1) , status )
      CASE DEFAULT
        !Should be Miller indices
        CALL MILLER2VEC(H,region_dir,ORIENT,addatom_prop(:),j)
        IF(j>0) THEN
          !Other error, unable to convert this string into a proper vector
          CALL ATOMSK_MSG(817,(/TRIM(region_dir)/),(/0.d0/))
          GOTO 1000
        ENDIF
        !Check that vector is not [000]
        IF( VECLENGTH(addatom_prop)<1.d-12 ) THEN
          CALL ATOMSK_MSG(814,(/""/),(/0.d0/))
          nerr=nerr+1
          GOTO 1000
        ENDIF
        vector(:) = addatom_prop(1)*H(1,:) + addatom_prop(2)*H(2,:) + addatom_prop(2)*H(3,:)
        CALL BOX2DBLE( vector(:) , treal(1) , region_1(1) , status )
      END SELECT
      IF(status>0) THEN
        temp = treal(1)
        GOTO 810
      ENDIF
    ELSEIF( region_side=="in" .OR. region_side=="out" ) THEN
      !"-select <in|out> cell"
      !"-select <in|out> box <x> <y> <z> <x'> <y'> <z'>"
      !"-select <in|out> sphere <x> <y> <z> <R>"
      !"-select <in|out> cylinder <axis> <x1> <x2> <R>"
      !"-select <in|out> cone <axis> <x> <y> <z> <alpha>"
      !"-select <in|out> torus <axis> <x> <y> <z> <R> <r>"
      region_dir = ""
      region_1(:) = 0.d0
      region_2(:) = 0.d0
      !IF( LEN_TRIM(select_multiple)>0 ) THEN
      !  READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, region_side, region_geom
      !ELSE
      !  READ(options_array(ioptions),*,END=800,ERR=800) optionname, region_side, region_geom
      !ENDIF
      region_geom = TRIM(substr(3+k))
      IF(region_geom=='cell') THEN
        !No other parameter
      ELSEIF(region_geom=='box' .OR. region_geom=='block') THEN
        !IF( LEN_TRIM(select_multiple)>0 ) THEN
          !READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp,&
          !    & region_side, region_geom, treal(1), treal(2), treal(3), &
          !    & treal(4), treal(5), treal(6)
!         ELSE
!           READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
!               & region_side, region_geom, treal(1), treal(2), treal(3), &
!               & treal(4), treal(5), treal(6)
!         ENDIF
        !Check if numbers contain a keyword like "BOX" or "INF"
        IF(SIZE(substr)<9+k) GOTO 800
        treal(1) = TRIM(substr(4+k))  !xmin
        treal(2) = TRIM(substr(5+k))  !ymin
        treal(3) = TRIM(substr(6+k))  !zmin
        treal(4) = TRIM(substr(7+k))  !xmax
        treal(5) = TRIM(substr(8+k))  !ymax
        treal(6) = TRIM(substr(9+k))  !zmax
        DO i=1,3
          CALL BOX2DBLE( H(:,i) , treal(i) , region_1(i) , status )
          IF(status>0) THEN
            temp = treal(i)
            GOTO 810
          ENDIF
        ENDDO
        DO i=1,3
          CALL BOX2DBLE( H(:,i) , treal(i+3) , region_2(i) , status )
          IF(status>0) THEN
            temp = treal(i+3)
            GOTO 810
          ENDIF
        ENDDO
      ELSEIF(region_geom=='sphere') THEN
!         IF( LEN_TRIM(select_multiple)>0 ) THEN
!           READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, &
!               & region_side, region_geom, treal(1), treal(2), treal(3), treal(4)
!         ELSE
!           READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
!               & region_side, region_geom, treal(1), treal(2), treal(3), treal(4)
!         ENDIF
        IF(SIZE(substr)<7+k) GOTO 800
        treal(1) = TRIM(substr(4+k))  !x
        treal(2) = TRIM(substr(5+k))  !y
        treal(3) = TRIM(substr(6+k))  !z
        READ(substr(7+k),*,END=800,ERR=800) region_2(1)  !radius (in angströms)
        !Check if numbers contain a keyword like "BOX" or "INF"
        DO i=1,3
          CALL BOX2DBLE( H(:,i) , treal(i) , region_1(i) , status )
          IF(status>0) THEN
            temp = treal(i)
            GOTO 810
          ENDIF
        ENDDO
      ELSEIF(region_geom=='cylinder') THEN
!         IF( LEN_TRIM(select_multiple)>0 ) THEN
!           READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, &
!               & region_side, region_geom, region_dir, treal(1), treal(2), treal(4)
!         ELSE
!           READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
!               & region_side, region_geom, region_dir, treal(1), treal(2), treal(4)
!         ENDIF
        IF(SIZE(substr)<7+k) GOTO 800
        region_dir = TRIM(substr(4+k)) !axis
        treal(1) = TRIM(substr(5+k))   !x
        treal(2) = TRIM(substr(6+k))   !y
        READ(substr(7+k),*,END=800,ERR=800) region_2(1)  !radius (in angströms)
        !Check if numbers contain a keyword like "BOX" or "INF"
        SELECT CASE(region_dir)
        CASE('x','X')
          CALL BOX2DBLE( H(:,2) , treal(1) , region_1(1) , status )
          IF(status>0) THEN
            temp = treal(1)
            GOTO 810
          ENDIF
          CALL BOX2DBLE( H(:,3) , treal(2) , region_1(2) , status )
          IF(status>0) THEN
            temp = treal(2)
            GOTO 810
          ENDIF
        CASE('y','Y')
          CALL BOX2DBLE( H(:,3) , treal(1) , region_1(1) , status )
          IF(status>0) THEN
            temp = treal(1)
            GOTO 810
          ENDIF
          CALL BOX2DBLE( H(:,1) , treal(2) , region_1(2) , status )
          IF(status>0) THEN
            temp = treal(2)
            GOTO 810
          ENDIF
        CASE('z','Z')
          CALL BOX2DBLE( H(:,1) , treal(1) , region_1(1) , status )
          IF(status>0) THEN
            temp = treal(1)
            GOTO 810
          ENDIF
          CALL BOX2DBLE( H(:,2) , treal(2) , region_1(2) , status )
          IF(status>0) THEN
            temp = treal(2)
            GOTO 810
          ENDIF
        CASE DEFAULT
          nerr=nerr+1
          GOTO 1000
        END SELECT
      ELSEIF(region_geom=='cone') THEN
!         IF( LEN_TRIM(select_multiple)>0 ) THEN
!           READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, region_side, &
!               & region_geom, region_dir, treal(1), treal(2), treal(3), region_2(1)
!         ELSE
!           READ(options_array(ioptions),*,END=800,ERR=800) optionname, region_side, &
!               & region_geom, region_dir, treal(1), treal(2), treal(3), region_2(1)
!         ENDIF
        IF(SIZE(substr)<8+k) GOTO 800
        region_dir = TRIM(substr(4+k)) !axis
        treal(1) = TRIM(substr(5+k))   !x
        treal(2) = TRIM(substr(6+k))   !y
        treal(3) = TRIM(substr(7+k))   !z
        treal(4) = TRIM(substr(8+k))   !z
        READ(substr(8+k),*,END=800,ERR=800) region_2(1)  !angle alpha
        !Check if numbers contain a keyword like "BOX" or "INF"
        DO i=1,3
          CALL BOX2DBLE( H(:,i) , treal(i) , region_1(i) , status )
          IF(status>0) THEN
            temp = treal(i)
            GOTO 810
          ENDIF
        ENDDO
      ELSEIF(region_geom=='torus' .OR. region_geom=='pyramid') THEN
!         IF( LEN_TRIM(select_multiple)>0 ) THEN
!           READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, region_side, &
!               & region_geom, region_dir, treal(1), treal(2), treal(3), region_2(1), region_2(2)
!         ELSE
!           READ(options_array(ioptions),*,END=800,ERR=800) optionname, region_side, &
!               & region_geom, region_dir, treal(1), treal(2), treal(3), region_2(1), region_2(2)
!         ENDIF
        IF(SIZE(substr)<9+k) GOTO 800
        region_dir = TRIM(substr(4+k)) !axis
        treal(1) = TRIM(substr(5+k))   !x
        treal(2) = TRIM(substr(6+k))   !y
        treal(3) = TRIM(substr(7+k))   !z
        READ(substr(8+k),*,END=800,ERR=800) region_2(1)  !primary radius (in angströms)
        READ(substr(9+k),*,END=800,ERR=800) region_2(2)  !secondary radius (in angströms)
        !Check if numbers contain a keyword like "BOX" or "INF"
        DO i=1,3
          CALL BOX2DBLE( H(:,i) , treal(i) , region_1(i) , status )
          IF(status>0) THEN
            temp = treal(i)
            GOTO 810
          ENDIF
        ENDDO
      ENDIF
    ELSEIF( region_side=="prop" .OR. region_side=="property" ) THEN
      !"-select prop <property> <value1[:value2]>"
      !store property name in "region_geom"
!       IF( LEN_TRIM(select_multiple)>0 ) THEN
!         READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, region_side, region_geom, temp
!       ELSE
!         READ(options_array(ioptions),*,END=800,ERR=800) optionname, region_side, region_geom, temp
!       ENDIF
      !store value(s) of the property into
      IF(SIZE(substr)<4+k) GOTO 800
      region_geom = TRIM(substr(3+k))
      temp = TRIM(substr(4+k))
      IF( temp(1:3)=="min" .OR. temp(1:3)=="max" ) THEN
        region_dir = TRIM(ADJUSTL(temp(1:3)))
      ELSE
        strlength = SCAN(temp,":")
        IF( strlength>0 ) THEN
          !user gives a range => read min and max values (can be -INF or +INF)
          treal(1) = temp(1:strlength-1)
          CALL BOX2DBLE( H(:,1) , treal(1) , region_1(1) , status )
          treal(2) = temp(strlength+1:)
          CALL BOX2DBLE( H(:,1) , treal(2) , region_1(2) , status )
          region_2(1) = 10.d0
          !make sure region_1(1) is smaller than region_1(2)
          IF( region_1(1) > region_1(2) ) THEN
            tempreal = region_1(2)
            region_1(2) = region_1(1)
            region_1(1) = tempreal
          ENDIF
        ELSE
          !user gives only one value => save it to region_1(1)
          READ(temp,*,END=800,ERR=800) region_1(1)
          region_2(1) = 0.d0
        ENDIF
      ENDIF
    ELSEIF( region_side=="random" ) THEN
      !"-select random <N> <species>"
      !store number of atoms N to region_1(1) and species to region_geom
!       IF( LEN_TRIM(select_multiple)>0 ) THEN
!         READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, &
!             & region_side, treal(1), region_geom
!       ELSE
!         READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
!             & region_side, treal(1), region_geom
!       ENDIF
      IF(SIZE(substr)<4+k) GOTO 800
      treal(1) = TRIM(substr(3+k))
      region_geom = TRIM(substr(4+k))
      j = INDEX(treal(1),"%")
      IF( j>0 ) THEN
        !treal(1) contains a percentage
        !=> append the symbol "%" to region_side
        region_side = "random%"
        !Read the number from treal(1)
        treal(1)(j:j) = " "
        READ(treal(1),*,END=800,ERR=800) region_1(1)
        region_1(1) = region_1(1) / 100.d0
      ELSE
        READ(treal(1),*,END=800,ERR=800) region_1(1)
      ENDIF
      !Read the number in treal(1) and save it into region_1(1)
      !CALL BOX2DBLE( (/1.d0,0.d0,0.d0/) , treal(1) , region_1(1) , status )
      !IF(status>0) THEN
      !  temp = treal(1)
      !  GOTO 810
      !ENDIF
    ELSEIF( region_side=="stl" ) THEN
      !"-select stl [center] [rescale|fill] <file>"
!       IF( LEN_TRIM(select_multiple)>0 ) THEN
!         READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, region_side, region_geom
!       ELSE
!         READ(options_array(ioptions),*,END=800,ERR=800) optionname, region_side, region_geom
!       ENDIF
      region_geom = TRIM(substr(3+k))
      IF( INDEX(region_geom,"center")>0 .OR. INDEX(region_geom,"scale")>0  &
        & .OR. INDEX(region_geom,"fill")>0  ) THEN
        IF( LEN_TRIM(select_multiple)>0 ) THEN
          READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, region_side, region_dir, region_geom
        ELSE
          READ(options_array(ioptions),*,END=800,ERR=800) optionname, region_side, region_dir, region_geom
        ENDIF
        region_dir = TRIM(substr(3+k))
        region_geom = TRIM(substr(4+k))
      ENDIF
    ELSEIF( region_side=="neigh" ) THEN
      !"-select <NNN> <species> neighbors <index>"
      !Store number of neighbors, or cutoff radius for neighbor search, into region_1(1)
      !Store species of neighbors (can be "all" or "any") in region_geom
      !Store index of atoms whose neighbors must be searched in region_1(2)
!       IF( LEN_TRIM(select_multiple)>0 ) THEN
!         READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, &
!             & region_side, region_1(1), region_geom, region_1(2)
!       ELSE
!         READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
!             & region_side, region_1(1), region_geom, region_1(2)
!       ENDIF
      IF(SIZE(substr)<5+k) GOTO 800
      READ(substr(3+k),*,END=800,ERR=800) region_1(1)
      region_geom = TRIM(substr(4+k))
      READ(substr(5+k),*,END=800,ERR=800) region_1(2)
    ELSEIF( region_side=="modulo" .OR. region_side=="mod" ) THEN
      !"-select <i> modulo <j>"
      !Store index into region_1(1) and modulo into region_1(2)
!       IF( LEN_TRIM(select_multiple)>0 ) THEN
!         READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, &
!             & region_side, region_1(1), region_1(2)
!       ELSE
!         READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
!             & region_side, region_1(1), region_1(2)
!       ENDIF
      IF(SIZE(substr)<4+k) GOTO 800
      READ(substr(3+k),*,END=800,ERR=800) region_1(1)
      READ(substr(4+k),*,END=800,ERR=800) region_1(2)
    ELSEIF( region_side=="list" .OR. region_side=="grid" ) THEN
      !store the name of file containing list of atoms in "region_geom"
!       IF( LEN_TRIM(select_multiple)>0 ) THEN
!         READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp, region_side, region_geom
!       ELSE
!         READ(options_array(ioptions),*,END=800,ERR=800) optionname, region_side, region_geom
!       ENDIF
      IF(SIZE(substr)<3+k) GOTO 800
      READ(substr(3+k),*,END=800,ERR=800) region_geom
    ENDIF
    IF(nerr.NE.0) GOTO 1000
    CALL SELECT_XYZ(H,P,AUXNAMES,AUX,select_multiple,region_side,region_geom,region_dir,region_1,region_2,ORIENT,SELECT)
  !
  CASE('-separate','sep')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, sep_radius, sep_dist
    CALL SEPARATE_XYZ(H,P,S,sep_radius,sep_dist,SELECT)
  !
  CASE('-shift')
    !"-shift [above d <normal>] dx dy dz"
    !shift_dir = TRIM(ADJUSTL(options_array(ioptions)(7:)))
    !j = SCAN(shift_dir," ")
    CALL STRCUT(options_array(ioptions),' ',substr)
    shift_dir = TRIM(substr(2))
    IF( shift_dir=='above' .OR. shift_dir=='below' ) THEN
      IF( SIZE(substr)<7 ) GOTO 800
      treal(4) = TRIM(substr(3))
      shift_axis = TRIM(substr(4))
      treal(1) = TRIM(substr(5))
      treal(2) = TRIM(substr(6))
      treal(3) = TRIM(substr(7))
      !Check if numbers contain a keyword like "BOX" or "INF"
      i=1
      SELECT CASE(shift_axis)
      CASE('x','X')
        CALL BOX2DBLE( H(:,1) , treal(4) , shift_dist , status )
      CASE('y','Y')
        CALL BOX2DBLE( H(:,2) , treal(4) , shift_dist , status )
      CASE('z','Z')
        CALL BOX2DBLE( H(:,3) , treal(4) , shift_dist , status )
      CASE DEFAULT
        !Should be Miller indices
        CALL MILLER2VEC(H,shift_axis,ORIENT,addatom_prop(:),j)
        IF(j>0) THEN
          !Other error, unable to convert this string into a proper vector
          CALL ATOMSK_MSG(817,(/TRIM(shift_axis)/),(/0.d0/))
          GOTO 1000
        ENDIF
        !Check that vector is not [000]
        IF( VECLENGTH(addatom_prop)<1.d-12 ) THEN
          CALL ATOMSK_MSG(814,(/""/),(/0.d0/))
          nerr=nerr+1
          GOTO 1000
        ENDIF
        vector(:) = addatom_prop(1)*H(1,:) + addatom_prop(2)*H(2,:) + addatom_prop(2)*H(3,:)
        CALL BOX2DBLE( vector(:) , treal(4) , shift_dist , status )
      END SELECT
      IF(status>0) THEN
        temp = treal(4)
        GOTO 810
      ENDIF
    ELSE
      shift_dir = ''
      shift_axis = ''
      IF( SIZE(substr)<4 ) GOTO 800
      treal(1) = TRIM(substr(2))
      treal(2) = TRIM(substr(3))
      treal(3) = TRIM(substr(4))
    ENDIF
    CALL BOX2DBLE( H(:,1) , treal(1) , shift_tau1 , status )
    IF(status>0) THEN
      temp = treal(1)
      GOTO 810
    ENDIF
    CALL BOX2DBLE( H(:,2) , treal(2) , shift_tau2 , status )
    IF(status>0) THEN
      temp = treal(2)
      GOTO 810
    ENDIF
    CALL BOX2DBLE( H(:,3) , treal(3) , shift_tau3 , status )
    IF(status>0) THEN
      temp = treal(3)
      GOTO 810
    ENDIF
    CALL SHIFT_XYZ(H,P,S,shift_dir,shift_dist,shift_axis,shift_tau1,shift_tau2,shift_tau3,ORIENT,SELECT)
  !
  CASE('-sort')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp
    IF( temp=="random" ) THEN
      sortorder="random"
      sortcol=""
    ELSEIF( temp=="reverse" .OR. temp=="r" ) THEN
      sortorder="reverse"
      sortcol=""
    ELSE
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, sortcol, sortorder
    ENDIF
    CALL SORT_XYZ(P,S,AUXNAMES,AUX,SELECT,sortcol,sortorder)
  !
  CASE('-spacegroup')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp
    CALL SPACEGROUP_XYZ(H,P,S,AUXNAMES,AUX,temp)
  !
  CASE('-stress')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, stress_in
    stress_in = ADJUSTL(stress_in)
    SELECT CASE(stress_in)
    CASE('x','X','xx','XX','y','Y','yy','YY','z','Z','zz','ZZ','xy','XY','yx','YX','zx','ZX','xz','XZ','zy','ZY','yz','YZ','p','P')
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, stress_in, stress
    CASE DEFAULT
      !stress_in contains the name of a file => file will be read within the option
    END SELECT
    CALL STRESS_XYZ(H,P,S,stress_in,stress,C_tensor)
  !
  CASE('-substitute', '-sub')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, sp1, sp2
    CALL SUBSTITUTE_XYZ(P,S,sp1,sp2,SELECT)
  !
  CASE('-swap')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, swap_id(1), swap_id(2)
    CALL SWAP_XYZ(H,P,S,AUXNAMES,AUX,swap_id,SELECT)
  !
  CASE('-torsion')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, rot_axis, rot_angle
    CALL TORSION_XYZ(H,P,S,AUXNAMES,AUX,rot_axis,rot_angle,SELECT)
  !
  CASE('-unit', '-u')
    READ(options_array(ioptions),'(a128)',END=800,ERR=800) temp  !optionname, unit1, unit2
    temp = ADJUSTL(temp)
    strlength = SCAN(temp,' ')
    temp = temp(strlength:)
    temp = ADJUSTL(temp)
    strlength = SCAN(temp,' ')
    unit1 = temp(1:strlength)
    temp = temp(strlength:)
    unit2 = ADJUSTL(temp)
    CALL UNIT_XYZ(H,P,S,AUXNAMES,AUX,unit1,unit2,SELECT)
  !
  CASE('-unskew')
    CALL UNSKEW_XYZ(H)
  !
  CASE('-velocity')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, vel_T
    CALL VELOCITY_XYZ(P,AUXNAMES,AUX,vel_T,SELECT)
  !
  CASE('-wrap')
    CALL WRAP_XYZ(H,P,S,SELECT)
  !
  ! -- please add other options in alphabetical order --
  !
  CASE DEFAULT
    !Could not understand option => display a warning
    nwarn = nwarn+1
    CALL ATOMSK_MSG(2700,(/options_array(ioptions)/),(/0.d0/))
  !
  END SELECT
  !
  !
  IF(verbosity==4) THEN
    !DEBUGGING
    !Write atoms positions to "atomsk.xyz"
    temp = 'atomsk.xyz'
    OPEN(UNIT=36,FILE=temp,STATUS="UNKNOWN",FORM="FORMATTED")
    IF( ALLOCATED(P) ) THEN
      WRITE(36,*) SIZE(P,1)
      msg = '#Debug file for Atomsk, after option: '//TRIM(optionname)
      WRITE(36,*) TRIM(msg)
      DO i=1, SIZE(P,1)
        CALL ATOMSPECIES(P(i,4),species)
        WRITE(36,'(a2,2X,3(f16.8,1X))') species, P(i,1:3)
      ENDDO
      IF( ALLOCATED(S) ) THEN
        DO i=1, SIZE(S,1)
          CALL ATOMSPECIES(S(i,4),species)
          species = TRIM(ADJUSTL(species))//"_s"
          WRITE(36,'(a2,2X,3(f16.8,1X))') species, S(i,1:3)
        ENDDO
      ENDIF
      WRITE(36,'(a4)') 'alat'
      WRITE(36,'(a3)') '1.0'
      WRITE(36,'(a9)') 'supercell'
      WRITE(36,'(3f16.6)') H(1,:)
      WRITE(36,'(3f16.6)') H(2,:)
      WRITE(36,'(3f16.6)') H(3,:)
    ELSE
      msg = '# ERROR: no atom remaining after option: '//TRIM(optionname)
      WRITE(36,*) TRIM(msg)
    ENDIF
    CLOSE(36)
    !
    IF( ANY(C_tensor(:,:)>0.1d0) ) THEN
      !Write elastic tensor to logfile
      msg = 'Elastic tensor (GPa):'
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO i=1,9
        WRITE(msg,'(9(e10.3,2X))') (C_tensor(i,j), j=1,9)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
    ENDIF
    !
    IF( ALLOCATED(AUXNAMES) ) THEN
      WRITE(msg,*) 'Size AUXNAMES, AUX:', SIZE(AUXNAMES), ',', &
                 & SIZE(AUX(:,1)), SIZE(AUX(1,:))
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !Write auxiliary properties to logfile
      msg = 'Auxiliary properties:'
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      WRITE(msg,*) (TRIM(ADJUSTL(AUXNAMES(i)))//'  ', i=1,MIN(SIZE(AUXNAMES),6) )
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO i=1,MIN( 20,SIZE(AUX(:,1)) )
        msg = ""
        DO j=1,MIN(SIZE(AUX,2),6)
          IF( DABS(AUX(i,j))<10000.d0 ) THEN
            WRITE(temp,'(f9.3)') AUX(i,j)
          ELSE
            WRITE(temp,'(e9.3)') AUX(i,j)
          ENDIF
          msg = TRIM(ADJUSTL(msg))//"  "//temp(1:9)
        ENDDO
        IF(j>=6) msg = TRIM(ADJUSTL(msg))//"  (...)"
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
      IF(i>=20) THEN
        WRITE(msg,*) '  (...discontinued...)'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDIF
    ENDIF
    !
    IF( ALLOCATED(SELECT) .AND. SIZE(SELECT)>0 ) THEN
      !Write positions of selected atoms to "atomsk_select.xyz"
      temp = 'atomsk.xyz'
      OPEN(UNIT=36,FILE=temp,STATUS="UNKNOWN",FORM="FORMATTED")
      WRITE(36,*) SIZE(P,1)
      msg = '#Debug file for Atomsk: SELECTed atoms only'
      j=0
      DO i=1,SIZE(SELECT)
        IF( SELECT(i) ) THEN
          j = j+1
          CALL ATOMSPECIES(P(i,4),species)
          WRITE(36,'(a2,2X,3(f16.8,1X))') species, P(i,1:3)
        ENDIF
      ENDDO
      WRITE(36,'(a4)') 'alat'
      WRITE(36,'(a3)') '1.0'
      WRITE(36,'(a9)') 'supercell'
      WRITE(36,'(3f16.6)') H(1,:)
      WRITE(36,'(3f16.6)') H(2,:)
      WRITE(36,'(3f16.6)') H(3,:)
      CLOSE(36)
      WRITE(msg,*) 'SELECTed atoms:', j
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
    !
  ENDIF !end v==4
  !
  !
  !In case of error break loop and exit
  IF(nerr.NE.0) EXIT
  !
  !
  IF( .NOT.ALLOCATED(P) .OR. SIZE(P,1)<=0 ) THEN
    !No atom left in system: warn user
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2768,(/""/),(/0.d0/))
    !
  ELSE
    !If shells exist, check that array S has a smaller size than P
    IF( ALLOCATED(S) ) THEN
      IF( SIZE(S,1) .NE. SIZE(P,1) ) THEN
        nerr = nerr+1
        CALL ATOMSK_MSG(2810,(/""/),(/0.d0/))
      ENDIF
    ENDIF
    !
    !If a selection exists, check for size consistency
    IF( ALLOCATED(SELECT) ) THEN
      IF( SIZE(SELECT)<=0 .OR. .NOT.ANY(SELECT(:)) ) THEN
        !No atom is selected => clear selection by de-allocating SELECT
        nwarn=nwarn+1
        IF( ALLOCATED(SELECT)) DEALLOCATE(SELECT)
        CALL ATOMSK_MSG(2750,(/""/),(/0.d0/))
      ELSEIF( SIZE(SELECT).NE.SIZE(P,1) ) THEN
        !Resize SELECT array so that its size matches the number of atoms
        CALL RESIZE_LOGICAL1(SELECT,SIZE(P,1),i)
        IF( i>0 ) THEN
          nerr=nerr+1
          CALL ATOMSK_MSG(818,(/"SELECT"/),(/0.d0/))
          GOTO 1000
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  !
  !If an elastic tensor was defined, check it
  IF( ANY( C_tensor(1:3,1:3).NE.0.d0 ) ) THEN
    !Check that there are no "NaN" values
    CALL CHECKNAN(C_tensor,i)
    IF(i.NE.0) THEN
      nerr=nerr+1
      CALL ATOMSK_MSG(2809,(/""/),(/0.d0/))
      GOTO 1000
    ENDIF
    !Check if tensor is symmetric, i.e. Cij = Cji
    CALL CTENSOR_SYMMETRIC(C_tensor,i)
    IF(i.NE.0) THEN
      nwarn=nwarn+1
      CALL ATOMSK_MSG(2740,(/""/),(/0.d0/))
    ENDIF
  ENDIF
  !
  !
  IF(ALLOCATED(substr)) DEALLOCATE(substr)
!
ENDDO  !End Loop on all Options
!
!
!End of options => check sizes of the arrays
IF(ALLOCATED(substr)) DEALLOCATE(substr)
IF( .NOT.ALLOCATED(P) .OR. SIZE(P,1)==0 ) THEN
  
  nerr = nerr+1
  GOTO 1000
ENDIF
CALL CHECK_ARRAY_CONSISTENCY(P,S,AUX,AUXNAMES,i)
IF( i.NE.0 ) THEN
  IF( i==1 ) THEN
    msg = 'S'
  ELSEIF( i==2 ) THEN
    msg = 'AUX'
  ELSEIF( i==3 ) THEN
    msg = 'AUXNAMES'
  ENDIF
  nerr=nerr+1
  CALL ATOMSK_MSG(1802,(/msg/),(/0.d0/))
ENDIF
!
GOTO 1000
!
!
!
800 CONTINUE
nerr=nerr+1
temp = options_array(ioptions)
CALL ATOMSK_MSG(2806,(/TRIM(temp)/),(/0.d0/))
CALL DISPLAY_HELP(temp)
GOTO 1000
!
810 CONTINUE
CALL ATOMSK_MSG(2813,(/TRIM(temp)/),(/0.d0/))
GOTO 1000
!
!
!
1000 CONTINUE
WRITE(msg,*) "================   END OF OPTIONS   ================"
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
WRITE(msg,*) "===================================================="
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
END SUBROUTINE OPTIONS_AFF
!
!
END MODULE
