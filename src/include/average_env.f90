MODULE avgenv
!
!**********************************************************************************
!*  AVGENV                                                                        *
!**********************************************************************************
!* This module parses all atoms in the system, detects the different types        *
!* of sites (based on central atom, number of neighbours...)                      *
!* and construct "reference" sites (Pref) by averaging all atoms that have        *
!* the same or similar environments.                                              *
!* If different lattices exist in the system (e.g. fcc and bcc), this module      *
!* should be able to recognize each as a different site/environment.              *
!* Atoms that have an odd number of neighbours (due to free surfaces,             *
!* vacancies, etc.) are not included in the averaging.                            *
!* Averaged environments are returned in array Pref(:,:,:)                        *
!* - first index is the site index [1-Nsites]                                     *
!* - second index is the atom index [1-Natoms];                                   *
!*   #1 = central atom, next ones are neighbours                                  *
!* - third index [1-6] is for atom position ([1-3]=x,y,z), atomic number [4],     *
!*   number of atoms in this type of site [5], and number of neighbours [6]       *
!* So, data for the first type of site is in Pref(1,:,:):                         *
!* - central atom position is Pref(1,1,1:3), atomic number Pref(1,1,4)            *
!* - number of atoms matching this type of site is in Pref(1,1,5)                 *
!* - number of neighbours for this type of site is Nneigh=Pref(1,1,6)             *
!* - rel.positions of neighbours in Pref(1,n,1:4), with n=2 up to Nneigh          *
!* Similarly, data for the second type of site is in Pref(2,:,:), and so on.      *
!* Values in Pref(:,2:N,5:6) are left empty.                                      *
!**********************************************************************************
!* (C) February 2023 - Pierre Hirel                                               *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 04 Nov. 2026                                     *
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
USE atoms
USE comv
USE constants
USE messages
USE files
USE neighbors
USE sorting
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE AVG_ENV(H,P,NeighList,Pref,NeighFactor,siteindex)
!
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=4096):: conffile
LOGICAL:: fileexists  !does the file exist?
INTEGER:: NeighIndex  !index of a neighbour
INTEGER:: iat, i, j, k, m, n
INTEGER:: Nneighbors
INTEGER:: NNmin=3     !minimum number of neighbours to keep
INTEGER:: Nsites, NneighMax
INTEGER:: progress      !To show calculation progress
INTEGER,DIMENSION(20):: Ctable !correspondence table
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex  !list of index after sorting
INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT):: siteindex !for each atom in P, index of its type of site in Pref
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList !neighbour list can be provided as input, or computed below
REAL(dp):: alpha, theta  !angle between two vectors
REAL(dp),INTENT(INOUT):: NeighFactor !%of tolerance in the radius for neighbor search
REAL(dp):: radius=6.d0  !R for neighbor search: 6 A should be enough to find some neighbors in any system
REAL(dp):: tempreal, tempreal2
REAL(dp),DIMENSION(9),PARAMETER:: known_angle= &
       & (/30.d0,45.d0,60.d0,90.d0,109.28d0,109.47122d0,120.d0,150.d0,180.d0/)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(3):: vector                !just a vector
REAL(dp),DIMENSION(3,3),INTENT(IN):: H        !cell vectors
REAL(dp),DIMENSION(:,:),INTENT(IN):: P        !atom positions
REAL(dp),DIMENSION(20,21,6):: Tref            !references for atoms environments (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList !list of positions of neighbors of one atom
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT):: Pref  !references for atoms environments (final)
!
WRITE(msg,*) "  ENTERING  AVG_ENV"
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
Nsites=0
Nneighmax=0
Tref(:,:,:) = 0.d0
IF(ALLOCATED(Pref)) DEALLOCATE(Pref)
!
IF(ALLOCATED(siteindex)) DEALLOCATE(siteindex)
ALLOCATE(siteindex(SIZE(P,1)))
siteindex(:) = 0
!
!If provided array NeighList already contains the neighbour list, use it
!Otherwise, construct the neighbour list now
IF( .NOT.ALLOCATED(NeighList) .OR. SIZE(NeighList,1).NE.SIZE(P,1) ) THEN
  msg = "  NeighList empty: constructing neighbour list ..."
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  CALL NEIGHBOR_LIST(H,P,radius,NeighList)
ENDIF
!
!

!!!!  Set up NeighFactor  !!!!
!Work with neighbours that are closer than neighbour #NNmin
!NOTE: first neighbours shell count 8 atoms in bcc, 12 atoms in fcc
!     However it counts only 6 atoms in NaCl structure. This will be adjusted later on
NNmin = 3
!After the neighbor search, neighbors will be sorted by increasing distance
!Neighbor #NNmin will be at distance d
!Only neighbors within a distance d*NeighFactor will be kept and used to compute G
!We want to captude 1st neighbors, but leave out the 2nd neighbors
!In fcc lattice the ratio between 2nd and 1st neighbors is sqrt(2)=~1.41
!In bcc lattice it is about 1.155, which is much smaller
NeighFactor = 1.1d0
!
!Check if user wrote specific values for these parameters in a configuration file
!(UNIX/Linux: "atomsk.conf", Windows: "atomsk.ini"), with the format:
! nye nnmin <value>
! nye neighfactor <value>
#if defined(WINDOWS)
conffile = "atomsk.ini"
#else
conffile = "atomsk.conf"
#endif
INQUIRE(FILE=conffile,EXIST=fileexists)
IF( fileexists ) THEN
  CALL ATOMSK_MSG(16,(/conffile/),(/0.d0/))
  OPEN(UNIT=31,FILE=conffile,STATUS="OLD",FORM="FORMATTED")
  DO
    n=1
    READ(31,'(a)',ERR=110,END=110) temp
    temp = TRIM(ADJUSTL(temp))
    IF( temp(1:1).NE.'#' ) THEN
      IF( StrDnCase(temp(1:3))=="nye" ) THEN
        temp = TRIM(ADJUSTL(temp(4:)))
        IF( StrDnCase(temp(1:5))=="nnmin" ) THEN
          READ(temp(6:),*,ERR=105,END=105) NNmin
        ELSEIF( StrDnCase(temp(1:11))=="neighfactor" ) THEN
          READ(temp(12:),*,ERR=105,END=105) NeighFactor
        ENDIF
      ENDIF
    ENDIF
    n=0
    105 CONTINUE
    IF( n>0 ) CALL ATOMSK_MSG(808,(/TRIM(temp)/),(/0.d0/))
  ENDDO
  110 CONTINUE
  CLOSE(31)
  !
ELSE
  !User did not define values in a file
  !Default values above fail for complex or distorted systems
  !=> Try to adjust them if system contains more than one atom type
  !Search how many atom types exist in the system
  CALL FIND_NSP(P(:,4),aentries)
  !Get number of different elements
  n = SIZE(aentries,1)
  IF( n>1 ) THEN
    !There is more than one element in this system
    !Get max. number of atoms for one element, multiply by 10%
    tempreal = 0.1d0*MAXVAL(aentries(:,2))
    IF( n==2 ) THEN
      !Compare relative concentrations of the two elements
      IF( aentries(1,2)>tempreal .AND. aentries(2,2)>tempreal ) THEN
        !All elements are present in large number => binary material
        !=> the default NeighFactor will probably be too small, use larger value
        IF( DABS(aentries(1,2)-aentries(2,2)) < tempreal ) THEN
          !Both elements are in comparable concentration
          !Assume rock-salt (NaCl, MgO...), where angle between first neighbours is 90°
          NNmin = 3
          NeighFactor = 1.3d0
        ELSE
          !One element is more concentrated than the other, e.g. SiO2, Al2O3...
          NNmin = 3
          NeighFactor = 1.3d0
        ENDIF
      ELSE
        !Otherwise, one element is present only in small concentration (<10%)
        !=> probably a unitary compound with few impurities, don't change NeighFactor
      ENDIF
      !
    ELSE IF( n>=3 ) THEN
      !Good chances that it is a complex material
      !=> boost NeighFactor
      NNmin = 3
      NeighFactor = 1.35d0   !1.33
    ENDIF
  ENDIF
ENDIF
!
!
100 CONTINUE
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "  AVG_ENV , NeighFactor = ", NeighFactor
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "  AVG_ENV ,       NNmin = ", NNmin
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  msg = "  AVG_ENV Parsing atoms to build reference averaged environments ..."
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!Loop on atoms in the system
!!!!!$OMP !PARALLEL DO DEFAULT(SHARED) &
!!!!!$OMP& !PRIVATE(i,iat,j,k,n,alpha,tempreal,PosList,newindex,Nneighbors) &
!!!!!$OMP& !REDUCTION(+:Tref(:,:,5))
DO iat=1,SIZE(P,1)
  !
  progress = progress+1
  !
  IF( SIZE(P,1) > 100000 ) THEN
    !If there are many atoms, display a fancy progress bar
    CALL ATOMSK_MSG(10,(/""/),(/DBLE(progress),DBLE(SIZE(P,1))/))
  ENDIF
  !
  !Using neighbor list, save Cartesian neighbours positions in PosList
  CALL NEIGHBOR_POS(H,P,iat,P(iat,1:3),NeighList(iat,:),ALLOCATED(NeighList),radius,PosList)
  !Now PosList(:,:) contains for each neighbor, its position (x,y,z), distance, and index in P
  !
  !Sort neighbours by increasing distance
  CALL BUBBLESORT(PosList,4,'up  ',newindex)
  IF(ALLOCATED(newindex)) DEALLOCATE(newindex)
  !
  !Get the relative positions of neighbors with respect to atom #iat
  !Keep only neighbors that are within NeighFactor times the distance of neighbor #NNmin
  Nneighbors = 0
  DO j=1,MIN(15,SIZE(PosList,1))
    !Shift all neighbours positions so the central atom is at (0,0,0)
    PosList(j,1:3) = PosList(j,1:3) - P(iat,1:3)
    IF( PosList(j,4) <= NeighFactor*PosList(NNmin,4) ) THEN
      !IF( NINT(Psecond(NINT(PosList(j,5)),4))==NINT(Psecond(NINT(PosList(1,5)),4)) ) THEN
        !Keep this neighbor
        Nneighbors = Nneighbors+1
      !ELSE
        !Ignore this neighbor (not same species as first neighbor)
      !ENDIF
    ELSE
      !Atom #j is too far, so will be the next ones: exit loop
      EXIT
    ENDIF
  ENDDO
  !
  !Wipe out the rest of the list
  PosList(Nneighbors+1:,:) = 0.d0
  !
  !
  !Consider atom #iat only if its number of neighbors is 4, 6, 8, 12
  !(only atoms in near-perfect environment are used to compute average, other atoms are discarded)
  !IF( Nneighbors>3 .AND. Nneighbors<14 .AND. MOD(Nneighbors,2)==0 .AND. MOD(Nneighbors,5)>0 ) THEN
  SELECT CASE(Nneighbors)
  CASE(4,6,8,12)
    !Search in Tref if central atom #iat already matches an environment in Tref
    !If not, create a new type of site in Tref
    DO i=1,SIZE(Tref,1)
!     IF(i>SIZE(Tref,1)) THEN
!       nerr=nerr+1
!       CALL ATOMSK_MSG(4070,(/""/),(/DBLE(i)/))
!       EXIT
!     ENDIF
      IF( NINT(P(iat,4))==NINT(Tref(i,1,4)) .AND. NINT(Tref(i,1,6))==Nneighbors ) THEN
        !Atom species and number of neighbors match: atom #iat occupies a site of the type #i
        !Save the index of the site for atom #iat
        siteindex(iat) = i
        !Increment counter of atoms for this type of site
        Tref(i,1,5) = Tref(i,1,5) + 1.d0
        EXIT
      ENDIF
    ENDDO
    !
    !If siteindex is still zero, no suitable site was found for this atom: create a new one
    IF( siteindex(iat)==0 ) THEN
      !Parse Tref(:,:,:) for an empty slot
      DO i=1,SIZE(Tref,1)
        IF( NINT(Tref(i,1,5))==0 ) THEN
          !This slot is empty: save central atom + neighbors
          !Position of central atom = (0,0,0)
          !Save species of central atom for this site
          Tref(i,1,4) = P(iat,4)
          !Increment counter of atoms for this type of site
          Tref(i,1,5) = Tref(i,1,5) + 1.d0
          !Save number of neighbors for this site
          Tref(i,1,6) = DBLE(Nneighbors)
          !Save the index of the site for atom #iat
          siteindex(iat) = i
          EXIT
        ENDIF
      ENDDO
    ENDIF
    !
    IF( i>0 ) THEN
      !Now we know atom #iat is in a site of type #i
      !Add neighbors positions to already known positions for averaging
      IF( NINT(Tref(i,1,5))<=1 ) THEN
        !No atom was found for this site yet, this is the first one
        !Simply store rel. pos. of neighbors of atom #iat in Tref(i,:,:)
        DO k=1,Nneighbors
          NeighIndex = NINT(PosList(k,5))
          IF(NeighIndex<=0) NeighIndex=iat
          Tref(i,k+1,1:3) = PosList(k,1:3)
          Tref(i,k+1,4) = P(NeighIndex,4)
        ENDDO
      ELSE
        !Some atoms were already found in this site
        !For each neighbor in PosList, find its best match in Tref
        !(i.e. the one that maximizes the dot product)
        !and save it into Ctable(:)
        Ctable(:) = 0
        DO k=1,Nneighbors  !Loop on neighbors in PosList
          alpha = 0.d0
          j=2
          DO WHILE( j<=NINT(Tref(i,1,6))+1 )  !Loop on sites in Tref(i,j,:)
            tempreal = DOT_PRODUCT( PosList(k,1:3) , Tref(i,j,1:3) )
            IF( tempreal > 1.01d0*alpha ) THEN
              !Site #j is a better match than previous sites
              Ctable(k) = j
              alpha = tempreal
            ENDIF
            j=j+1
          ENDDO
        ENDDO
        !Each value in Ctable(:) should be unique, otherwise it means that
        !two different atoms were assigned to the same reference atom
        n=0
        DO k=1,Nneighbors-1
          DO j=k+1,Nneighbors
            IF( Ctable(k)==Ctable(j) ) THEN
              !This is a duplicate: do not consider this atom for averaging
              n=1
            ENDIF
            IF (Ctable(j)==0) n=1
          ENDDO
          IF (Ctable(k)==0) n=1
        ENDDO
        !
        IF( n==0 ) THEN
          !Each neighbor was paired with a reference: add this site to the average
          DO k=1,Nneighbors
            !Ctable(k) is the index of the neighbor in Tref that is the best match to atom #k in PosList
            !Save this index in Ctable
            n = Ctable(k)
            !Add position of atom #n and perform averaging
            Tref(i,n,1:3) = ( Tref(i,n,5)*Tref(i,n,1:3) + PosList(k,1:3) ) / (Tref(i,n,5)+1.d0)
            !Save atomic number of this neighbor (only if it's empty)
            IF( Tref(i,n,4)<1.d0 ) THEN
              NeighIndex = NINT(PosList(k,5))
              IF(NeighIndex<=0) NeighIndex=iat
              Tref(i,n,4) = P(NeighIndex,4)
            ENDIF
            !Increment number of neighbors at this position
            Tref(i,n,5) = Tref(i,n,5) + 1.d0
          ENDDO !k
        ENDIF !end if n==0
        !
      ENDIF  !end if (NINT(Tref<=1))
    ENDIF
    !
  END SELECT  !CASE(Nneighbors)
  !ENDIF !end if (Nneighbors...)
  !
ENDDO !loop on iat
!!!!!$OMP !END PARALLEL DO
!
!
Nsites = MAXVAL(siteindex)
!
IF(verbosity>=4) THEN
  WRITE(msg,*) "  AVG_ENV Sites detected BEFORE CLEANUP: ", Nsites
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(Tref,1)
    WRITE(msg,*) "  AVG_ENV  Tref(", i, ") :", Tref(i,1,5), " atoms found in this env. BEFORE CLEANUP"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    IF( NINT(Tref(i,1,5))>0 ) THEN
      DO j=1,NINT(Tref(i,1,6))+1
        WRITE(msg,'(6X,6f9.3)') Tref(i,j,1:6)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
    ENDIF
  ENDDO
ENDIF
!
!Check that some sites are not duplicate, and discard sites with low occurrence
n=Nsites
IF( Nsites>1 ) THEN
  !Double loop on all types of sites
  DO i=1,n-1
    DO j=i+1,n
      !Compare the atomic number of central atoms in sites #i and #j
      IF( NINT(Tref(i,1,4)) == NINT(Tref(j,1,4)) .AND.   &
        & NINT(Tref(i,1,5))>0 .AND. NINT(Tref(j,1,5))>0  ) THEN
        !Sites i and j have the same central atom
        !Verify if one of them has a very low occurrence, or has very few neighbours
        IF( NINT(Tref(i,1,5)) <= MAX(10,NINT(Tref(j,1,5)/10.d0)) .OR. Tref(i,1,6)<3 &
          & .OR. .NOT.ANY(NINT(Tref(i,2:,4))>0)  ) THEN
          !Site of type #i has a very low occurrence compared to site #j
          !Add its occurrence to the site type #j
          !Tref(j,1,5) = Tref(j,1,5)+Tref(i,1,5)
          !Delete site #i from the list
          Tref(i,:,:) = 0.d0
          !Update total number of types of sites
          Nsites = Nsites-1
          !Update siteindex: for all atoms in site #i, change site to type #j
          DO k=1,SIZE(siteindex)
            IF(siteindex(k)==i) siteindex(k) = j
          ENDDO
        ELSEIF( NINT(Tref(j,1,5)) <= MAX(10,NINT(Tref(i,1,5)/10.d0)) .OR. Tref(j,1,6)<3 &
          & .OR. .NOT.ANY(NINT(Tref(j,2:,4))>0)  ) THEN
          !Site of type #j has a very low occurrence compared to site #i
          !Add its occurrence to the site type #i
          !Tref(i,1,5) = Tref(i,1,5)+Tref(j,1,5)
          !Delete site #j from the list
          Tref(j,:,:) = 0.d0
          !Update total number of types of sites
          Nsites = Nsites-1
          !Update siteindex: for all atoms in site #j, change site to type #i
          DO k=1,SIZE(siteindex)
            IF(siteindex(k)==j) siteindex(k) = i
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDIF
!
!
!Get number of different sites and max. number of neighbours
!Nsites = MAXVAL(siteindex)
Nneighmax = MAXVAL(Tref(:,1,6))
IF(verbosity>=4) THEN
  WRITE(msg,*) "  AVG_ENV Sites detected: ", Nsites
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "  AVG_ENV  Max. N.neigh.:  ", Nneighmax
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(Tref,1)
    WRITE(msg,*) "  AVG_ENV  Tref(", i, ") :", Tref(i,1,5), " atoms found in this environment"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    IF( NINT(Tref(i,1,5))>0 ) THEN
      DO j=1,NINT(Tref(i,1,6))+1
        WRITE(msg,'(6X,6f9.3)') Tref(i,j,1:6)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
    ENDIF
  ENDDO
ENDIF
!
!
!Save averaged sites in final array Pref(:,:,:)
IF( Nsites>0 ) THEN
  ALLOCATE( Pref(Nsites,Nneighmax+1,6) )
  Pref(:,:,:) = 0.d0
  n=0
  DO i=1,SIZE(Tref,1)  !Loop on all sites in Tref
    IF( NINT(Tref(i,1,5)) > 0 .AND. NINT(Tref(i,1,6)) > 0 ) THEN
      !This is a valid site, save it in Pref
      n=n+1
      DO j=1,NINT(Tref(i,1,6))+1
        Pref(n,j,:) = Tref(i,j,:)
      ENDDO
      !Update siteindex: for all atoms in site #i, change site to type #n
      DO k=1,SIZE(siteindex)
        IF(siteindex(k)==i) siteindex(k) = n
      ENDDO
    ENDIF
  ENDDO
ELSE
  !No atomic environment found: make sure Pref and siteindex are returned unallocated
  IF(ALLOCATED(Pref)) DEALLOCATE(Pref)
  IF(ALLOCATED(siteindex)) DEALLOCATE(siteindex)
ENDIF
!
!
200 CONTINUE
GOTO 1000  !this is not fully implemented yet, skip for now
IF( ALLOCATED(Pref) ) THEN
  !
  !Now Pref contains the averaged relative positions of neighbors for each type of atom site
  !For each central atom, compute angles between neighbors, and if they are close to
  !some known/expected angles, correct them
  !
  DO i=1,SIZE(Pref,1)
    !Loop on all neighbors
    DO j=2,NINT(Pref(i,1,6))-1
      !Loop on all other neighbors
      DO k=j,NINT(Pref(i,1,6))
        n=0
        !Compute angle between neighbors #j and #k
        alpha = ANGVEC( Pref(i,j,1:3)-Pref(i,1,1:3) , Pref(i,k,1:3)-Pref(i,1,1:3) )
        !Check if angle alpha is close to a known value
        DO m=1,SIZE(known_angle)
          IF( DABS(RAD2DEG(alpha)-known_angle(m)) < 0.1d0 ) THEN
            n=m
            EXIT
          ENDIF
        ENDDO
        IF( n>0 ) THEN
          !Correct position of neighbor #k
          vector(:) = Pref(i,j,1:3)-Pref(i,1,1:3)
          tempreal = VECLENGTH( Pref(i,k,1:3)-Pref(i,1,1:3) )

          Pref(i,k,1:3) = Pref(i,1,1:3) + vector(1:3)
        ENDIF
      ENDDO !k
    ENDDO !j
  ENDDO !i
  !
ENDIF !end if allocated(Pref)
!
!
1000 CONTINUE
!Write some debugging information
IF( verbosity==4 ) THEN
  IF( ALLOCATED(Pref) .AND. SIZE(Pref,1)>0 ) THEN
    !For each site, write position of central atom (that should always be (0,0,0))
    !and positions of neighbors into a XYZ file named "atomsk_site_i.xyz" for visualization
    WRITE(msg,*) " ATOMIC  ENVIRONMENTS  DETECTED:", Nsites, SIZE(Pref,1)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO i=1,MIN(Nsites,SIZE(Pref,1))
      !Check if number of neighbours for site #i is positive
      WRITE(msg,*) i
      OPEN(UNIT=23,FILE="atomsk_site_"//TRIM(ADJUSTL(msg))//".xyz",FORM="FORMATTED")
      WRITE(23,*) NINT(Pref(i,1,6))+1
      WRITE(23,*) "# Atomsk AVG_ENV: Averaged environment for site #", i
      CALL ATOMSPECIES(Pref(i,1,4),species)
      WRITE(23,'(a2,3X,3f16.8)') species, Pref(i,1,1:3)  !central atom
      IF( NINT(Pref(i,1,6))>0 ) THEN
        WRITE(msg,*) i
        WRITE(temp,*) NINT(Pref(i,1,6))
        WRITE(msg,*) "  Site # "//TRIM(ADJUSTL(msg))//" occupied by "//species// &
                    & " atom, has "//TRIM(ADJUSTL(temp))//" neighbors:"
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        DO j=2,NINT(Pref(i,1,6))+1  !loop on neighbors
          CALL ATOMSPECIES(Pref(i,j,4),species)
          WRITE(23,'(a2,3X,3f16.8)') species, Pref(i,j,1:3)
          WRITE(msg,'(8X,a2,3X,3f12.3)') species, Pref(i,j,1:3)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDDO
      ENDIF
      CLOSE(23)
    ENDDO
    OPEN(UNIT=23,FILE="atomsk_siteindex.txt",FORM="FORMATTED")
    WRITE(23,*) SIZE(P,1)
    WRITE(23,*) "# Atomsk AVG_ENV: atom id and type of site they occupy"
    DO i=1,SIZE(siteindex,1)
      WRITE(23,*) i, siteindex(i)
    ENDDO
    CLOSE(23)
  ELSE
    msg = " ERROR:  NO ATOM ENVIRONMENT DETECTED!"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
ENDIF
!
msg = "  EXITING  AVG_ENV"
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!
END SUBROUTINE AVG_ENV
!
END MODULE avgenv
