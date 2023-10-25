! ===========================================================================================
! A module that contains functions to move from reduced repesentation of hydrogen bond
! topology to an all atom structure
! 
! To use in a fortran program insert
! use output_potential
! 
! Then call setup_potential and setup_struct to allocate all the various arrays
! The exact form of outut can then be chosen (current options: DCD/LAMMPS/XYZ/gromacs)
! Call destroy_struct when finished
!
! ===========================================================================================

MODULE OUTPUT_POTENTIAL
IMPLICIT NONE

INTEGER,PARAMETER        :: dp = selected_real_kind(p=15,r=300)
REAL(KIND=dp),PARAMETER  :: pi = 4.0_dp*ATAN(1.0_dp)


! Parameters of potential model

TYPE POTENTIAL
  PRIVATE
  REAL*8, DIMENSION(0:4) :: Hangle, Hlength
  LOGICAL                :: OH, H3O, NH4, Msite
  REAL*8                 :: OM, dphi
END TYPE POTENTIAL

TYPE STRUC
   PRIVATE
   INTEGER, DIMENSION(:,:),   ALLOCATABLE :: conns
   REAL,    DIMENSION(:,:),   ALLOCATABLE :: fixedPos
   INTEGER, DIMENSION(:),     ALLOCATABLE :: val
   REAL*8,  DIMENSION(:),     ALLOCATABLE :: msite
   INTEGER                                :: nmols, natoms
   REAL                                   :: x, y, maxz, minz
END TYPE STRUC

CONTAINS

! ===========================================================================================
SUBROUTINE SETUP_POTENTIAL(H2Oang, H2Olen, pot, Nang, HlenN, H3Oang, HlenH3O, HlenOH, dOM)
IMPLICIT NONE

TYPE(POTENTIAL),  INTENT(OUT) :: pot
REAL*8,           INTENT(IN)  :: H2Oang, H2Olen
REAL*8, OPTIONAL, INTENT(IN)  :: Nang, HlenN, H3Oang, HlenH3O, HlenOH, dOM

pot%Hangle  = (H2Oang/180.0_dp)*pi   ! Angle between two hydrogens
pot%Hlength = H2Olen                 ! Length of TIP4P bond 

pot%OH = .FALSE. ; pot%H3O = .FALSE. ; pot%NH4 = .FALSE. ; pot%Msite = .FALSE.

IF ((PRESENT(Nang)) .AND. (PRESENT(HlenN))) THEN
   pot%Hangle(4)  = (Nang/180.0_dp)*pi     ! Angle between two hydrogens
   pot%Hlength(4) = HlenN                  ! Length of NH4+ bond 
   pot%NH4        = .TRUE.
END IF


IF ((PRESENT(H3Oang)) .AND. (PRESENT(HlenH3O))) THEN
   pot%Hangle(3)  = (H3Oang/180.0_dp)*pi      ! Half the angle between two hydrogens
   pot%Hlength(3) = HlenH3O                   ! Length of H3O bond 
   pot%H3O        = .TRUE.
   ! Want to work out the dphi for angluar separation of H3Oang for two Hs with the same theta
   ! Used whenever we refine an H-O-H angle in H3O so caluculate only once for efficiency
   ! From the great angle formula
   pot%dphi = ACOS(((COS(pot%Hangle(3)) - COS(pot%Hangle(3))**2)/(SIN(pot%Hangle(3))**2)))
END IF

IF (PRESENT(HlenOH)) THEN
   pot%Hangle(1)  = (H2Oang/180.0_dp)*pi ! Angle between two hydrogens of water molecule (OH has no angle)
   pot%Hlength(1) = HlenOH               ! Length of OH bond 
   pot%OH         = .TRUE.
END IF

IF (PRESENT(dOM)) THEN
   pot%OM    = dOM
   pot%Msite = .TRUE.
END IF

END SUBROUTINE SETUP_POTENTIAL
! =========================================================================
SUBROUTINE SETUP_STRUCT(pos, pot, nmols, ions, x, y, z, zoffset, struct, incMs)
IMPLICIT NONE

TYPE(POTENTIAL),         INTENT(INOUT) :: pot
TYPE(STRUC),             INTENT(OUT)   :: struct
REAL, DIMENSION(:,:),    INTENT(INOUT) :: pos
INTEGER, DIMENSION(:,:), INTENT(IN)    :: ions
INTEGER,                 INTENT(IN)    :: nmols
REAL,                    INTENT(IN)    :: x, y, z, zoffset
LOGICAL, OPTIONAL,       INTENT(IN)    :: incMs

INTEGER               :: i, j, k, O, Hpos
REAL                  :: mindist, invx, invy, invz, dist, zp
INTEGER, DIMENSION(2) :: numIons
REAL,    DIMENSION(3) :: tmpvect
LOGICAL               :: Ms

struct%x = x ; struct%y = y ; struct%maxz = z ; struct%minz = zoffset

zp = z - zoffset

invx = 1/x ; invy = 1/y ; invz = 1/zp

Ms = .FALSE.
IF (PRESENT(incMs)) THEN
   Ms = incMs
   IF (Ms) WRITE(*, *) '#INFO: Using M-sites.'
END IF

struct%nmols = nmols

ALLOCATE(struct%val(nmols))
ALLOCATE(struct%msite(nmols))

struct%val = 2       ! Default is H2O with valence of 2

struct%msite = MERGE(pot%OM, 0.0_8, Ms)    ! To avoid having to input potential later on

numIons = SHAPE(ions)

DO i = 1, numIons(1)
   struct%val(ions(i, 1))   = ions(i, 2)
   struct%msite(ions(i, 1)) = 0
END DO

IF ((COUNT(struct%val(:) .EQ. 1) .NE. 0) .AND. (.NOT. pot%OH)) THEN
   WRITE(*, *) 'ERROR: Structure contains OH- ion but potential has not been specified.'
   STOP
END IF

IF ((COUNT(struct%val(:) .EQ. 3) .NE. 0) .AND. (.NOT. pot%H3O)) THEN
   WRITE(*, *) 'ERROR: Structure contains H3O+ ion but potential has not been specified.'
   STOP
END IF

IF ((COUNT(struct%val(:) .EQ. 4) .NE. 0) .AND. (.NOT. pot%NH4)) THEN
   WRITE(*, *) 'ERROR: Structure contains NH4+ ion but potential has not been specified.'
   STOP
END IF

IF ((Ms) .AND. (.NOT. pot%Msite)) THEN
   WRITE(*, *) 'ERROR: Resquested output of TIP4P oxygen M-site but potential has not been specified.'
   STOP
END IF

struct%natoms = SUM(struct%val) + struct%nmols
IF (Ms) struct%natoms = struct%natoms + COUNT(struct%val(:) .EQ. 2) ! Account for M-sites

ALLOCATE(struct%conns(1:struct%nmols, 1:4))
ALLOCATE(struct%fixedPos(1:3, 1:struct%nmols))

struct%conns = -1

DO i = 1, nmols
   ! Find connections matrix
   struct%fixedPos(:, i) = pos(:, (3*i - 2))

   DO k = 1, 2
      O = i  ! The O this hydrogen is covalently bonded to

      ! Now need to determine the next closest O (i.e. this H is H bonded to)

      mindist = 500.0         ! Set a high minimum distance for finding next closest oxygen  
      Hpos    = 3*O - 2 + k   ! Position of H in pos array
      DO j = 1, nmols
         IF (j .EQ. O) CYCLE ! Want H bond not covalent bond

         tmpvect = pos(:, Hpos) - pos(:, (3*j - 2)) ! This is the position of O in the pos array
         ! Check images

         tmpvect(1) = tmpvect(1) -  x*ANINT(tmpvect(1)*invx)
         tmpvect(2) = tmpvect(2) -  y*ANINT(tmpvect(2)*invy)
         tmpvect(3) = tmpvect(3) - zp*ANINT(tmpvect(3)*invz)

         dist       = SQRT(tmpvect(1)**2 + tmpvect(2)**2 + tmpvect(3)**2)
         
         IF (dist .LT. mindist) THEN
            struct%conns(O, k) = j
            mindist = dist
         END IF
      END DO

      ! Now str%connection(O, 1) and str%connection(O, 2) are ones where O is cov bonded to H
      ! where H is then H bonded to j

     IF (struct%conns(struct%conns(O, k), 3) .NE. -1) THEN
         struct%conns(struct%conns(O, k), 4) = O
      ELSE
         struct%conns(struct%conns(O, k), 3) = O
      END IF
   END DO
END DO

! Don't need to worry about basal planes or anything like that as just want connections matrix

END SUBROUTINE SETUP_STRUCT
! =========================================================================
SUBROUTINE DESTROY_STRUCT(struct)
IMPLICIT NONE

TYPE(STRUC), INTENT(INOUT) :: struct

DEALLOCATE(struct%conns, struct%val, struct%fixedPos, struct%msite)

END SUBROUTINE DESTROY_STRUCT
! =========================================================================
!SUBROUTINE EXPAND_HASH(str, hash, full_hash)
!IMPLICIT NONE
!TYPE(STRUC),               INTENT(INOUT) ::
!CHARACTER(LEN=*),          INTENT(IN)    :: hash
!CHARACTER(LEN=4*str%nmol), INTENT(OUT)   :: full_hash

!CHARACTER(LEN=*), PARAMETER :: mapping = 'abcdefghijklmnop'

! Need this temporarily but turns out this won't work for noPBCs (as well
! as removing an independent check) so need to write the bitmask thing anyway



! =========================================================================
SUBROUTINE OUTPUT_POSITIONS(hash, str, pot, pos, midO)
IMPLICIT NONE

TYPE(STRUC),                        INTENT(INOUT) :: str
TYPE(POTENTIAL),                    INTENT(INOUT) :: pot
CHARACTER(LEN=*),                   INTENT(IN)    :: hash  ! Hash string detailing bonding configuration
REAL, DIMENSION(1:3, 1:str%natoms), INTENT(OUT)   :: pos
LOGICAL, OPTIONAL,                  INTENT(IN)    :: midO

INTEGER                    :: i, n, j, k
CHARACTER                  :: curr_hash
CHARACTER(LEN=4*str%nmols) :: full_hash
REAL, DIMENSION(3)         :: tmpvec, cross, H1vec, H2vec
REAL, DIMENSION(3, 4)      :: all_vec
REAL                       :: invx, invy, invz, zp, theta, phi, thetap, phip, phia, phib, AA, BB, a, b, c
LOGICAL                    :: pmz, mdO

zp = str%maxz - str%minz

invx = 1/str%x ; invy = 1/str%y ; invz = 1/zp
mdO = .FALSE.
IF (PRESENT(midO)) THEN
   mdo = midO
   IF (mdO) WRITE(*, *) '#INFO: Distribulting Hs between O-O bonds.'
END IF


!IF (LEN(hash) .LT. str%nmols) THEN
!   ! Using a reduced hash for memory requirements 
!   full_hash = EXPAND_HASH(hash)
!ELSE
!   full_hash = hash
!END IF

full_hash = hash

n = 1
DO i = 1, str%nmols
   pos(:, n) = str%fixedPos(:, i)
   n = n + 1
   IF (str%val(i) .EQ. 1) THEN
      ! F- has no attached hydrogens, H2O/H3O/NH4 require finessing
      DO j = 1, 4
         curr_hash = full_hash(4*i - 4 + j:4*i - 4 + j)
         IF (curr_hash .EQ. '1') THEN
            tmpvec    = str%fixedPos(:, str%conns(i, j)) - str%fixedPos(:, i) 
            tmpvec(1) = tmpvec(1) - str%x*ANINT(tmpvec(1)*invx)
            tmpvec(2) = tmpvec(2) - str%y*ANINT(tmpvec(2)*invy)
            tmpvec(3) = tmpvec(3) -    zp*ANINT(tmpvec(3)*invz)
            tmpvec(:) = tmpvec/SQRT(DOT_PRODUCT(tmpvec, tmpvec))

            pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*tmpvec(:)         
            n         = n + 1
         END IF
      END DO
   ELSE IF (str%val(i) .EQ. 2) THEN
      ! Now have all of the H positions, need to refine to ensure consistent lengths and angles
      k = 0
      DO j = 1, 4
         curr_hash = full_hash(4*i - 4 + j:4*i - 4 + j)
         IF ((curr_hash .EQ. '1') .AND. (k .EQ. 0)) THEN
            H1vec    = str%fixedPos(:, str%conns(i, j)) - str%fixedPos(:, i) 
            H1vec(1) = H1vec(1) - str%x*ANINT(H1vec(1)*invx)
            H1vec(2) = H1vec(2) - str%y*ANINT(H1vec(2)*invy)
            H1vec(3) = H1vec(3) -    zp*ANINT(H1vec(3)*invz)
            H1vec(:) = H1vec/SQRT(DOT_PRODUCT(H1vec, H1vec))
            k = k + 1
         ELSE IF ((curr_hash .EQ. '1') .AND. (k .EQ. 1)) THEN
            H2vec    = str%fixedPos(:, str%conns(i, j)) - str%fixedPos(:, i) 
            H2vec(1) = H2vec(1) - str%x*ANINT(H2vec(1)*invx)
            H2vec(2) = H2vec(2) - str%y*ANINT(H2vec(2)*invy)
            H2vec(3) = H2vec(3) -    zp*ANINT(H2vec(3)*invz)
            H2vec(:) = H2vec/SQRT(DOT_PRODUCT(H2vec, H2vec))
            EXIT
         END IF
      END DO

      IF (mdO) THEN
         ! Want to find the positions that ensures the angular distance of each H from its ideal position
         ! (along O-O bond) is the same
         
         ! The two ideal points are H1vec and H2vec - want to move first H along this great circle by
         ! half of the extra bit of the angle and second back by this (AND THEN CHECK ANGLE)

         theta = ACOS(DOT_PRODUCT(H1vec, H2vec))  ! The OOO angle
         phi   = theta - pot%Hangle(2)            ! How much bigger the OOO angle is than the HOH angle

         cross    = CROSS_PRODUCT(H2vec, H1vec)
         tmpvec   = CROSS_PRODUCT(cross, H2vec)

         H1vec    = tmpvec*SIN(pot%Hangle(str%val(i))+0.5*phi) + H2vec*COS(pot%Hangle(str%val(i))+0.5*phi)  ! Unit vector

         cross    = CROSS_PRODUCT(H1vec, H2vec)
         tmpvec   = CROSS_PRODUCT(cross, H1vec)
         
         H2vec    = tmpvec*SIN(pot%Hangle(str%val(i))) + H1vec*COS(pot%Hangle(str%val(i)))  ! Unit vector
      
      ELSE
         IF (ABS(H2vec(3)) .GE. 0.95) THEN
            ! Don't move if in +/- z for dipole consistency (-z is always 1)
            cross    = CROSS_PRODUCT(H2vec, H1vec)
            tmpvec   = CROSS_PRODUCT(cross, H2vec)
            
            H1vec    = tmpvec*SIN(pot%Hangle(str%val(i))) + H2vec*COS(pot%Hangle(str%val(i)))  ! Unit vector
         ELSE
            cross    = CROSS_PRODUCT(H1vec, H2vec)
            tmpvec   = CROSS_PRODUCT(cross, H1vec)
         
            H2vec    = tmpvec*SIN(pot%Hangle(str%val(i))) + H1vec*COS(pot%Hangle(str%val(i)))  ! Unit vector         
         END IF
      END IF

         pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*H1vec
         n         = n + 1
         pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*H2vec
         n         = n + 1

         
   ELSE IF (str%val(i) .EQ. 4) THEN
      
      DO j = 1, 4
         all_vec(:, j) = str%fixedPos(:, str%conns(i, j)) - str%fixedPos(:, i) 
         all_vec(1, j) = all_vec(1, j) - str%x*ANINT(all_vec(1, j)*invx)
         all_vec(2, j) = all_vec(2, j) - str%y*ANINT(all_vec(2, j)*invy)
         all_vec(3, j) = all_vec(3, j) -    zp*ANINT(all_vec(3, j)*invz)
         all_vec(:, j) = all_vec(:, j)/SQRT(DOT_PRODUCT(all_vec(:, j), all_vec(:, j)))
         
         IF (ABS(all_vec(3, j)) .GE. 0.95) k = j
      END DO

      H1vec = all_vec(:, k)
      
      IF (k .EQ. 1) THEN
         H2vec = all_vec(:, 2)
      ELSE
         H2vec = all_vec(:, 1)
      END IF

      cross    = CROSS_PRODUCT(H1vec, H2vec)
      tmpvec   = CROSS_PRODUCT(cross, H1vec)
      
      H2vec    = tmpvec*SIN(pot%Hangle(str%val(i))) + H1vec*COS(pot%Hangle(str%val(i)))  ! Unit vector

      pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*H1vec
      n         = n + 1
      pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*H2vec
      n         = n + 1

      DO j = 3, str%val(i)
         
         H1vec = H2vec 

         IF (j .EQ. k) THEN
            H2vec = all_vec(:, 2)
         ELSE
            H2vec    = all_vec(:, j)
         END IF

         cross    = CROSS_PRODUCT(H1vec, H2vec)
         tmpvec   = CROSS_PRODUCT(cross, H1vec)
         
         H2vec    = tmpvec*SIN(pot%Hangle(str%val(i))) + H1vec*COS(pot%Hangle(str%val(i)))
                  
         pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*H2vec
         n         = n + 1

      END DO

   ! Assume correct definition of ions
   
   ELSE 

      k   = 0
      pmz = .FALSE.

      DO j = 1, 4
         IF (full_hash(4*i - 4 + j:4*i - 4 + j) .NE. '1') k = j

         all_vec(:, j) = str%fixedPos(:, str%conns(i, j)) - str%fixedPos(:, i) 
         all_vec(1, j) = all_vec(1, j) - str%x*ANINT(all_vec(1, j)*invx)
         all_vec(2, j) = all_vec(2, j) - str%y*ANINT(all_vec(2, j)*invy)
         all_vec(3, j) = all_vec(3, j) -    zp*ANINT(all_vec(3, j)*invz)
         all_vec(:, j) = all_vec(:, j)/SQRT(DOT_PRODUCT(all_vec(:, j), all_vec(:, j)))
         
         IF (ABS(all_vec(3, j)) .GE. 0.95) THEN
            ! We care if it's +/- z for finding the other values
            theta = MERGE(pot%Hangle(str%val(i)), pi-pot%Hangle(str%val(i)), (all_vec(3, j) .GT. 0))

            H1vec(1) = 0  ; H1vec(2) = 0  ; H1vec(3) = MERGE(1, -1, all_vec(3, j) .GT. 0)

            IF (k .NE. j) THEN
               ! For ease of calculation, want to ensure that this always stays at +/- z
               pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*H1vec
               n         = n + 1
               pmz       = .TRUE.
            END IF
         END IF
         IF (k .NE. j) WRITE(*,*) str%fixedPos(:, i)+pot%Hlength(str%val(i))*all_vec(:, j), pot%Hlength(str%val(i))*all_vec(:, j)
      END DO

      IF (pmz) THEN
         ! If pmz is true:
         ! 1. Find another connected H, refine its theta to alpha/pi-alpha, and determine phi
         ! 2. Use current position of other connected atom to determine if this is +dphi or -dphi from other atom
         
         DO j = 1, 4
            IF (k .EQ. j) CYCLE            
            IF (ABS(all_vec(3, j)) .GT. 0.95) CYCLE
            
            EXIT
         END DO

         ! This is an H lying on the theta line
         ! Refine tmpvec(3) to be COS(theta), then determine phi
            
         all_vec(3, j) = COS(theta)
         all_vec(:, j) = all_vec(:, j)/SQRT(DOT_PRODUCT(all_vec(:, j), all_vec(:, j)))
         
         pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*all_vec(:, j)
         n         = n + 1
         
         phi = ATAN2(all_vec(2, j), all_vec(1, j))


         DO k = j+1, 4
            IF (full_hash(4*i - 4 + k:4*i - 4 + k) .NE. '1') CYCLE         
            IF (ABS(tmpvec(3)) .GT. 0.95) CYCLE
            
            EXIT
         END DO
         
         ! Ok, there is a potential issue here due to signs maybe?

         ! This is the third H - we know its theta and one of two options for phi
         phip = phi-pot%dphi         
         
         H1vec(1) = SIN(theta)*COS(phip) 
         H1vec(2) = SIN(theta)*SIN(phip)
         H1vec(3) = COS(theta)

         phip = phi+pot%dphi
         
         H2vec(1) = SIN(theta)*COS(phip) 
         H2vec(2) = SIN(theta)*SIN(phip)
         H2vec(3) = COS(theta)

         IF ((ABS(H1vec(1) - all_vec(1, k)) .LE. 0.25) .AND. (ABS(H1vec(2) - all_vec(2, k)) .LE. 0.25)) THEN
            ! At phi-dphi
            pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*H1vec
         ELSE IF ((ABS(H2vec(1) - all_vec(1, k)) .LE. 0.25) .AND. (ABS(H2vec(2) - all_vec(2, k)) .LE. 0.25)) THEN
            ! At phi+dphi
            pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*H2vec
         ELSE
            WRITE(*,*) "ERROR DETERMINING H3O POSITIONS" 
            STOP
         END IF
         n = n + 1


      ELSE
         ! If pmz is false:
         ! 1. Two connected Hs are ~ in positions (x, y, z) (-x, y, z). Allow for (x, y, z) (x, -y, z) as well.
         ! 2. Refine z to cos(theta) and angle between them
         ! 2. Use derived formula to obtain the position of the other atom
         
         ! Then create new positions, and ensure that the offset from old positions is negligible (CAN REMOVE ONCE TESTED)
         
         H1vec(1) = 0 

         DO j = 1, 4
            IF (k .EQ. j) CYCLE            
            IF ((ABS(all_vec(1, j)) .LT. 0.05) .OR. (ABS(all_vec(2, j)) .LT. 0.05)) CYCLE
            
            IF (ABS(H1vec(1)) .LT. 0.05) THEN
               H1vec = all_vec(:, j)
            ELSE
               H2vec = all_vec(:, j)
            END IF
         END DO
         
         ! We have the positions of the +/- bond and the x=0 bond
         
         phi  = ATAN2(H1vec(2), H1vec(1))
         phip = ATAN2(H2vec(2), H2vec(1))

         ! Difference in angle between current angle and angle it's supposed to be.

         phia = ABS(phi - phip) - pot%dphi

         IF (phi .GT. phip) THEN
           phi  = phi  - 0.5*phia
           phip = phip + 0.5*phia
         ELSE
            phi  = phi  + 0.5*phia
            phip = phip - 0.5*phia
         END IF
         
         H1vec(1) = COS(phi)*SIN(theta)
         H1vec(2) = SIN(phi)*SIN(theta)
         H1vec(3) = COS(theta)

         H2vec(1) = COS(phip)*SIN(theta)
         H2vec(2) = SIN(phip)*SIN(theta)
         H2vec(3) = COS(theta)
         
         pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*H1vec
         n         = n + 1
         pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*H2vec
         n         = n + 1

         ! Confident in these two positions

         phia = 0.5*(phi + phip)
         phib = MODULO(phia + 3*pi, 2.0*pi)
       
         a = COS(theta)
         b = SIN(theta)*COS(ABS(phia-phi))
         c = COS(pot%Hangle(str%val(i)))
       
         AA = 2*b*c
         BB = (a*a - b*b)

         BB = BB*(MERGE(1, -1, (a*c .GT. 0)))
         
         phip   = phia
         thetap = ATAN2(AA, BB)

         IF (ABS(COS(thetap) - COS(theta)) .GT. 0.25) THEN
            ! This thetap is the pole
            b = SIN(theta)*COS(ABS(phib-phi))
            
            AA = 2*b*c
            BB = (a*a - b*b)
            BB = BB*(MERGE(1, -1, (a*c .GT. 0)))
         
            phip   = phib
            thetap = ATAN2(AA, BB)

         END IF

         H2vec(1) = COS(phip)*SIN(thetap)
         H2vec(2) = SIN(phip)*SIN(thetap)
         H2vec(3) = COS(thetap)

      
         pos(:, n) = str%fixedPos(:, i) + pot%Hlength(str%val(i))*H2vec
         n         = n + 1
      END IF
      
   END IF

   ! For outputting positions, this distinguishes between H3O+ and H2O with an M-site  - now need to deal with those
   IF (str%msite(i) .GE. 0.05) THEN     ! Avoid FPEs
      tmpvec = (pos(:, n - 2) - str%fixedPos(:, i)) + (pos(:, n - 1) - str%fixedPos(:, i))
      tmpvec = tmpvec/SQRT(DOT_PRODUCT(tmpvec, tmpvec))
      
      pos(:, n) = str%fixedPos(:, i) + tmpvec*str%msite(i)
      n         = n + 1
   END IF
END DO

IF ((n - 1) .NE. str%natoms) THEN
   WRITE(*,*) 'ERROR: Number of atoms is not as expected. Check ions have been correctly defined.'
   STOP
END IF

END SUBROUTINE OUTPUT_POSITIONS
! =========================================================================
SUBROUTINE WRITE_XYZ(hash, str, pot, outfile)
IMPLICIT NONE

TYPE(STRUC),      INTENT(INOUT) :: str
TYPE(POTENTIAL),  INTENT(INOUT) :: pot
CHARACTER(LEN=*), INTENT(IN)    :: hash
CHARACTER(LEN=*), INTENT(IN)    :: outfile

REAL, DIMENSION(1:3, 1:str%natoms) :: pos
INTEGER                            :: fnum = 93, i, k, n

! Get positions
CALL OUTPUT_POSITIONS(hash, str, pot, pos)

! Write an xyz file of atomic positions
OPEN(fnum, file = TRIM(outfile), status='new')

WRITE(fnum, *) str%natoms
WRITE(fnum, *) '# xyz file produced by refine_bonds.f90 from atomic hash'

n = 1

DO i = 1, str%nmols
   IF (str%val(i) .EQ. 0) THEN
      ! This is an F- 
      WRITE(fnum, '(A5, F12.8, A5, F12.8, A5, F12.8)') 'F ', pos(1, n), ' ', pos(2, n), ' ', pos(3, n)
   ELSE IF (str%val(i) .EQ. 4) THEN
      ! This is an NH4+
      WRITE(fnum, '(A5, F12.8, A5, F12.8, A5, F12.8)') 'N ', pos(1, n), ' ', pos(2, n), ' ',  pos(3, n)
   ELSE
      ! OH-/H2O/H3O+ - central atom is O
      WRITE(fnum, '(A5, F12.8, A5, F12.8, A5, F12.8)') 'O ', pos(1, n), ' ', pos(2, n), ' ', pos(3, n)
   END IF
   n = n + 1
   DO k = 1, str%val(i)
      ! Should write right number of H (i.e. skip for F)
      WRITE(fnum, '(A5, F12.8, A5, F12.8, A5, F12.8)') 'H ', pos(1, n), ' ', pos(2, n), ' ',  pos(3, n)
      n = n + 1
   END DO
   ! Deal with M-sites if appropriate
   IF (str%msite(i) .GE. 0.05) THEN     ! Avoid FPEs
      WRITE(fnum, '(A5, F12.8, A5, F12.8, A5, F12.8)') 'M ', pos(1, n), ' ', pos(2, n), ' ',  pos(3, n)
      n = n + 1
   END IF
END DO

WRITE(fnum, *)
WRITE(fnum, '(A8, F12.8, A1, F12.8, A1, F12.8)') 'Vector1 ', str%x, ' ', 0.0, ' ', 0.0
WRITE(fnum, '(A8, F12.8, A1, F12.8, A1, F12.8)') 'Vector2 ', 0.0, ' ', str%y, ' ', 0.0
WRITE(fnum, '(A8, F12.8, A1, F12.8, A1, F12.8)') 'Vector3 ', 0.0, ' ', 0.0, ' ', str%maxz-str%minz
WRITE(fnum, '(A8, F12.8, A1, F12.8, A1, F12.8)') 'Offset  ', 0.0, ' ', 0.0, ' ', str%minz

CLOSE(fnum)


END SUBROUTINE WRITE_XYZ
! =========================================================================
SUBROUTINE WRITE_DCD_HEADER(dcd_file, str, snaps)
IMPLICIT NONE

! From D. Quigley
TYPE(STRUC),      INTENT(INOUT) :: str
INTEGER,          INTENT(IN)    :: snaps ! Number of snapshots to save
CHARACTER(LEN=*), INTENT(IN)    :: dcd_file

CHARACTER(LEN=4)           :: hdr='CORD'
INTEGER, DIMENSION(20)     :: dcd_header_info     ! called icntrl in Q's code
INTEGER                    :: ierr, dcd=208, headerlines=20, nlines, i
CHARACTER*80,DIMENSION(32) :: dcdtitle

WRITE(*,*) '#WARNING: Unit cell information for dcd files currently assumes cell with orthogonal axes.'

OPEN(unit = dcd, file = dcd_file, form = 'unformatted', status = 'new')

dcd_header_info(1)     = snaps                ! Number of snapshots in history file
dcd_header_info(2)     = 0                    ! Starting timestep
dcd_header_info(3)     = 1                    ! Gap in steps between snapshots (doesn't matter)
dcd_header_info(4)     = snaps                ! Total numberof steps (VMD ignores this)
dcd_header_info(5:7)   = 0
dcd_header_info(8)     = str%natoms           ! Ndeg 
dcd_header_info(9)     = 0                    ! Number of fixed atoms
dcd_header_info(10)    = 0                    ! Timestep
dcd_header_info(11)    = 1                    ! 1/0 for unit cell presence
dcd_header_info(12:19) = 0
dcd_header_info(20)    = 24                   ! CHARMM version number (fixes dcd format)


WRITE(dcd) hdr, dcd_header_info

WRITE(dcd) 1, (dcdtitle(i), i = 1, 1)
WRITE(dcd) str%natoms

CLOSE(dcd)

END SUBROUTINE WRITE_DCD_HEADER
! =========================================================================
SUBROUTINE WRITE_DCD_SNAPSHOT(hash, str, pot, dcd_file)
IMPLICIT NONE

TYPE(STRUC),      INTENT(INOUT) :: str
TYPE(POTENTIAL),  INTENT(INOUT) :: pot
CHARACTER(LEN=*), INTENT(IN)    :: hash
CHARACTER(LEN=*), INTENT(IN)    :: dcd_file

REAL(KIND=dp), DIMENSION(6)                 :: hmatrix
REAL,          DIMENSION(1:3, 1:str%natoms) :: pos
INTEGER                                     :: ierr, dcd=109, headerlines=20, nlines, i

! Get positions
CALL OUTPUT_POSITIONS(hash, str, pot, pos)

! Create Hmatrix
hmatrix    = 1.0_dp*90
hmatrix(1) = str%x*1.0_dp ; hmatrix(3) = str%y*1.0_dp ; hmatrix(6) = (str%maxz - str%minz)*1.0_dp

OPEN(unit=dcd, file=dcd_file, status='old', position='append', form='unformatted')
WRITE(dcd)hmatrix
WRITE(dcd)(pos(1, i), i = 1, str%natoms)
WRITE(dcd)(pos(2, i), i = 1, str%natoms)
WRITE(dcd)(pos(3, i), i = 1, str%natoms)

CLOSE(dcd)

END SUBROUTINE WRITE_DCD_SNAPSHOT
! =========================================================================
SUBROUTINE WRITE_LAMMPS(hash, str, pot, outfile, mH, mO, qH2O_H, qH2O_O, mN, mF, qF, qN, qNH, qOH_O, qOH_H, &
                        qH3O_O, qH3O_H, minE, ionO, flipz, repx, repy)
IMPLICIT NONE

TYPE(STRUC),       INTENT(INOUT) :: str
TYPE(POTENTIAL),   INTENT(INOUT) :: pot
CHARACTER(LEN=*),  INTENT(IN)    :: hash
CHARACTER(LEN=*),  INTENT(IN)    :: outfile
REAL,              INTENT(IN)    :: mH, qH2O_O, mO, qH2O_H
REAL, OPTIONAL,    INTENT(IN)    :: mN, mF, qF, qN, qNH, qOH_O, qOH_H, qH3O_O, qH3O_H
LOGICAL, OPTIONAL, INTENT(IN)    :: minE, ionO, flipz
INTEGER, OPTIONAL, INTENT(IN)    :: repx, repy  ! Create this number of explicit replications in x/y for non-periodic system 

REAL,    DIMENSION(1:3,1:str%natoms) :: pos
REAL,    DIMENSION(1:3)              :: dipole
REAL*8,  DIMENSION(0:4)              :: charges, Hcharges
INTEGER, DIMENSION(0:4)              :: indices, Hindices
INTEGER, DIMENSION(1:4)              :: bonds
INTEGER, DIMENSION(2:4)              :: angles
INTEGER                              :: fnum = 74, i, j, k, n, atn, atk, ry, rx, xreps, yreps
CHARACTER(LEN=63)                    :: form = '(I9, A1, I9, A2, I1, A2, F9.6, A1, F12.8, A1, F12.8, A1, F12.8)'
LOGICAL                              :: Ms, IOs, flip, mnE


charges = qH2O_O ; Hcharges = qH2O_H ; indices = 1 ; Hindices = 2
bonds = 1 ; angles = 1

Ms = ((COUNT(str%msite .GE. 0.05) .NE. 0))

IOs = .FALSE.
IF (PRESENT(ionO)) IOs = ionO

flip = .FALSE.
IF (PRESENT(flipz)) flip = flipz

mnE = .FALSE.
IF (PRESENT(minE)) mnE = minE


! Get positions
CALL OUTPUT_POSITIONS(hash, str, pot, pos, midO=.TRUE.)

IF (flip) pos(3, :) = pos(3, :)*(-1) + str%maxz + str%minz

dipole = GET_DIPOLE(str, pos, Ms, qH2O_H, qH2O_O, qF, qN, qNH, qOH_O, qOH_H, qH3O_O, qH3O_H)
WRITE(*,*) dipole

! Write a LAMMPS file of atomic positions
OPEN(fnum, file = TRIM(outfile), status = 'new')
WRITE(fnum, *) '# input file produced by refine_bonds.f90 from atomic hash'

! UC parameters always go from 0 but due to lack of z periodicity account for bonds pointing down
WRITE(fnum, *)
IF (PRESENT(repx)) THEN
   WRITE(fnum, '(F12.8, A1, F12.8, A8)') str%minz, ' ', repx*str%x+ABS(str%minz), ' xlo xhi'   ! Take into account maximum extent of bonds in non-periodic direction
ELSE
   WRITE(fnum, '(F12.8, A1, F12.8, A8)') 0.0, ' ', str%x, ' xlo xhi' 
END IF

IF (PRESENT(repy)) THEN
   WRITE(fnum, '(F12.8, A1, F12.8, A8)') str%minz, ' ', repy*str%y+ABS(str%minz), ' ylo yhi'   ! Take into account maximum extent of bonds in non-periodic direction
ELSE
   WRITE(fnum, '(F12.8, A1, F12.8, A8)') 0.0, ' ', str%y, ' ylo yhi'  
END IF

IF (mnE) THEN
   WRITE(*, *) '#INFO: Adding vacuum gap to LAMMPS input for energy minimisation.'
   WRITE(fnum, '(F12.8, A1, F12.8, A8)') str%minz, ' ', 4*str%maxz, ' zlo zhi'  

ELSE
   WRITE(fnum, '(F12.8, A1, F12.8, A8)') str%minz, ' ', str%maxz, ' zlo zhi'  
END IF

WRITE(fnum, *)

xreps = MERGE(repx, 1, PRESENT(repx))
yreps = MERGE(repy, 1, PRESENT(repy))

! Number of atoms, bonds, angles is self-explanatory
WRITE(fnum, *)
WRITE(fnum, '(I4, A6)') xreps*yreps*str%natoms, ' atoms'

IF (Ms) THEN
   ! Number of bonds is total valence plus number of m-sites
   WRITE(fnum, '(I4, A6)') xreps*yreps*(SUM(str%val(:)) + COUNT(str%val .EQ. 2)), ' bonds'
   ! Number of angles is nmol (H2O) + 3*H3O + 5*NH4 (6 angles per NH4 mol, 4C2) plus two per m-site
   WRITE(fnum, '(I4, A7)') xreps*yreps*(3*(COUNT(str%val(:).EQ.2) + (COUNT(str%val(:) .EQ. 3)) + &
                                        2*(COUNT(str%val(:) .EQ. 4)))), ' angles' 
   charges(2) = 0.0
ELSE
   ! Number of bonds is total valence
   WRITE(fnum, '(I4, A6)') xreps*yreps*(SUM(str%val(:))), ' bonds'
   ! Number of angles is nmol (H2O) + 3*H3O + 5*NH4 (6 angles per NH4 mol, 4C2)
   WRITE(fnum, '(I4, A7)') xreps*yreps*(COUNT(str%val(:).EQ.2) + 3*(COUNT(str%val(:) .EQ. 3)) + &
                                     6*(COUNT(str%val(:) .EQ. 4))), ' angles'  
END IF

WRITE(fnum, *)

! Number of types of things
WRITE(fnum, *)
! For atom types only need to check NH4/F- and one of H3O+/OH- (can't use TIP4P pairstyle for ions)
i = 2 ! Atom types
j = 1 ! Bond types
k = 1 ! Angle types

IF (Ms) THEN
   i = i + 1  ; j = j + 1  ; k = k + 1
END IF

indices    = 2 ! To deal with several ions
indices(2) = 1

IF (COUNT(str%val(:) .EQ. 0) .NE. 0) THEN
   IF((.NOT. PRESENT(mF)) .OR. (.NOT. PRESENT(qF))) THEN
      WRITE(*,*) 'ERROR: Structure contains F- ion but details have not been specified.'
      STOP
   END IF
   charges(0) = qF
   i = i + 1  ; indices(0) = i ! Some F-
END IF

IF (COUNT(str%val(:) .EQ. 1) .NE. 0) THEN
   IF ((.NOT. PRESENT(qOH_O)) .OR. (.NOT. PRESENT(qOH_H))) THEN
      WRITE(*,*) 'ERROR: Structure contains OH- ion but details have not been specified.'
      STOP
   END IF
   charges(1) = qOH_O ; Hcharges(1) = qOH_H
   IF ((.NOT. Ms) .AND. (.NOT. IOs)) WRITE(*, *) 'WARNING: Using OH- or H3O+ ions in LAMMPS TIP4P input ', &
                                                  'requires either M-site or ionic oxygen creation.'
   IF (IOs) i = i + 1 
   indices(1) = MERGE(i, 2, IOs)
   indices(3) = MERGE(i, 2, IOs)
   j = j + 1  ; bonds(1)   = j ! New bond, OH-
END IF

IF (COUNT(str%val(:) .EQ. 3) .NE. 0) THEN ! H3O+
   IF ((.NOT. PRESENT(qH3O_O)) .OR. (.NOT. PRESENT(qH3O_H))) THEN
      WRITE(*,*) 'ERROR: Structure contains H3O+ ion but details have not been specified.'
      STOP
   END IF
   charges(3) = qH3O_O ; Hcharges(3) = qH3O_H
   IF (COUNT(str%val(:) .EQ. 1) .EQ. 0) THEN
      ! No OH-, need to determine if we need to count for a new ionic O
      IF ((.NOT. Ms) .AND. (.NOT. IOs)) WRITE(*, *) 'WARNING: Using OH- or H3O+ ions in LAMMPS TIP4P input ', &
                                                    'requires either M-site or ionic oxygen creation.'
      IF (IOs) i = i + 1 
   END IF
   indices(3) = MERGE(i, 2, IOs)
   j = j + 1     ; k = k + 1 ! New bond, new angle
   bonds(3) = j  ; angles(3) = k
END IF

IF (COUNT(str%val(:) .EQ. 4) .NE. 0) THEN ! NH4+
   IF (((.NOT. PRESENT(qN)) .OR. (.NOT. PRESENT(qNH))) .OR. (.NOT. PRESENT(mN))) THEN
      WRITE(*,*) 'ERROR: Structure contains NH4+ ion but details have not been specified.'
      STOP
   END IF
   charges(4) = qN  ; Hcharges(4) = qNH
   i = i + 1        ; j = j + 1     ; k = k + 1 ! New atom, new bond, new angle 
   indices(4) = i   ; bonds(4) = j  ; angles(4) = k
END IF

WRITE(fnum, '(I1, A11)') i, ' atom types'
WRITE(fnum, '(I1, A11)') j, ' bond types' 
WRITE(fnum, '(I1, A12)') k, ' angle types'
WRITE(fnum, *) 

! Define atoms with their masses
WRITE(fnum, *) 
WRITE(fnum, *) 'Masses'
WRITE(fnum, *)
WRITE(fnum, '(I1, A1, F12.8, A10)') indices(2),     ' ', mO, ' # Water O'
WRITE(fnum, '(I1, A1, F12.8, A4)')  indices(2) + 1, ' ', mH, ' # H'

IF (Ms)                                          WRITE(fnum, '(I1, A1, F12.8, A4)')  indices(2) + 2, ' ', 1e-8, ' # M'    ! LAMMPS will not accept a truly massless M-site
IF (indices(0)  .NE. 2)                          WRITE(fnum, '(I1, A1, F12.8, A4)')  indices(0),     ' ',  mF, ' # F'
IF (indices(1)  .NE. 2)                          WRITE(fnum, '(I1, A1, F12.8, A10)') indices(1),     ' ',  mO, ' # Ionic O'
IF ((indices(3) .NE. 2) .AND. indices(1) .EQ. 2) WRITE(fnum, '(I1, A1, F12.8, A10)') indices(3),     ' ',  mO, ' # Ionic O'
IF (indices(4)  .NE. 2)                          WRITE(fnum, '(I1, A1, F12.8, A4)')  indices(4),     ' ',  mN, ' # N'

WRITE(fnum, *)

! Write atomic positions
WRITE(fnum, *)
WRITE(fnum, *) 'Atoms'
WRITE(fnum, *) '# id mol type  charge         x            y            z ' 

n = 1
DO rx = 0, xreps - 1  
   DO ry = 0, yreps - 1 
      atn = 1
      DO i = 1, str%nmols
         IF ((pos(3, atn) .LT. str%minz) .OR. (pos(3, atn) .GT. str%maxz)) &
              WRITE(*,*) 'WARNING: z outside box with coordinate ', pos(3, atn)
         WRITE(fnum, form) n, ' ' , i + str%nmols*(rx*yreps + ry), '  ', indices(str%val(i)), '  ', charges(str%val(i)), &
                              ' ', pos(1, atn) + rx*str%x, ' ', pos(2, atn) + ry*str%y, ' ', pos(3, atn)
         n   = n   + 1
         atn = atn + 1 
         DO k = 1, str%val(i)
            IF ((pos(3, atn) .LT. str%minz) .OR. (pos(3, atn) .GT. str%maxz)) &
                 WRITE(*,*) 'WARNING: z outside box with coordinate ', pos(3, atn)
            WRITE(fnum, form) n, ' ', i + str%nmols*(rx*yreps + ry), '  ', Hindices(str%val(i)), '  ', Hcharges(str%val(i)), &
                                 ' ', pos(1, atn) + rx*str%x, ' ', pos(2, atn) + ry*str%y, ' ', pos(3, atn)
            n   = n   + 1
            atn = atn + 1
         END DO
         IF ((Ms) .AND. (str%val(i) .EQ. 2)) THEN
            WRITE(fnum, form) n, ' ', i + str%nmols*(rx*yreps + ry),  ' ', Hindices(str%val(i)) + 1, '  ', qH2O_O, & 
                                 ' ', pos(1, atn) + rx*str%x, ' ', pos(2, atn) + ry*str%y, ' ',pos(3, atn) ! M-site carries all of the charge on O
            n   = n   + 1 
            atn = atn + 1
         END IF
      END DO
   END DO
END DO

WRITE(fnum, *)
WRITE(fnum, *)

! Write bonds information
WRITE(fnum, *)
WRITE(fnum, *) 'Bonds'
WRITE(fnum, *) '# bondid type atm1 atm2'

n   = 1
atn = 1
DO rx = 1, xreps
   DO ry = 1, yreps
      DO i = 1, str%nmols
         IF (str%val(i) .EQ. 0) THEN
            ! F- has no bonds, but still need to skip over to ensure correct numbering of atoms
            atn = atn + 1
            CYCLE
         END IF
         DO k = 1, str%val(i)
            WRITE(fnum, '(I8, A2, I1, A2, I4, A1, I4)') n, '  ', bonds(str%val(i)), '  ', atn, ' ', atn + k
            n = n + 1
         END DO
         IF ((Ms) .AND. (str%val(i) .EQ. 2)) THEN
            WRITE(fnum, '(I8, A2, I1, A2, I4, A1, I4)') n, '  ', bonds(str%val(i)) + 1, '  ', atn, ' ', atn + 3
            n   = n + 1 
            atn = atn + 2 + str%val(i) 
         ELSE
            atn = atn + 1 + str%val(i)
         END IF
      END DO
   END DO
END DO

WRITE(fnum, *)
WRITE(fnum, *)

! Write angles information
WRITE(fnum, *)
WRITE(fnum, *) 'Angles'
WRITE(fnum, *) '# anglid type atm1 atm2 atm3'

n   = 1
atn = 1
atk = 2 
DO rx = 1, xreps
   DO ry = 1, yreps
      DO i = 1, str%nmols
         IF (str%val(i) .LT. 2) THEN
            ! OH- and F- have no angles, but still need to ensure atoms skipped
            atn = atn + 1 + str%val(i)
            atk = atn + 1
            CYCLE
         END IF
         DO k=str%val(i) - 1, 1, -1
            DO j = 1, k
               WRITE(fnum, '(I8, A2, I1, A2, I4, A1, I4, A1, I4)') n, '  ', angles(str%val(i)), '  ', atk, ' ', atn, ' ', atk + j 
               n = n + 1
            END DO
            atk = atk + 1 ! To ensure correct progression through atoms
         END DO
         IF ((Ms) .AND. (str%val(i) .EQ. 2)) THEN
            WRITE(fnum, '(I8, A2, I1, A2, I4, A1, I4, A1, I4)') n, '  ', angles(str%val(i)) + 1, '  ', atn + 1, &
                                                                    ' ', atn, ' ', atn + 3 
            n = n + 1 
            WRITE(fnum, '(I8, A2, I1, A2, I4, A1, I4, A1, I4)') n, '  ', angles(str%val(i)) + 1, '  ', atn + 2, &
                                                                    ' ', atn, ' ', atn + 3
            n   = n + 1 
            atn = atn + 2 + str%val(i) 
         ELSE
            atn = atn + 1 + str%val(i)
         END IF
         atk = atn + 1
      END DO
   END DO
END DO

CLOSE(fnum)


END SUBROUTINE WRITE_LAMMPS
! ===========================================================================================
SUBROUTINE WRITE_GROMACS(hash, str, pot, outfile, mH, mO, qH2O_H, qH2O_O, mN, mF, qF, qN, qNH, qOH_O, qOH_H, &
                         qH3O_O, qH3O_H, repx, repy, minE, ionO, wall)
IMPLICIT NONE

TYPE(STRUC),       INTENT(INOUT) :: str
TYPE(POTENTIAL),   INTENT(INOUT) :: pot
CHARACTER(LEN=*),  INTENT(IN)    :: hash
CHARACTER(LEN=*),  INTENT(IN)    :: outfile
REAL,              INTENT(IN)    :: mH, qH2O_O, mO, qH2O_H
REAL, OPTIONAL,    INTENT(IN)    :: mN, mF, qF, qN, qNH, qOH_O, qOH_H, qH3O_O, qH3O_H
INTEGER, OPTIONAL, INTENT(IN)    :: repx, repy
LOGICAL, OPTIONAL, INTENT(IN)    :: minE, ionO, wall

REAL,             DIMENSION(1:3,1:str%natoms) :: pos
REAL*8,           DIMENSION(0:4)              :: charges, Hcharges, masses
CHARACTER(LEN=5), DIMENSION(0:4)              :: res, atom, Hatom, number, Msite
CHARACTER(LEN=2), DIMENSION(0:4)              :: name
INTEGER                                       :: fnum = 74, i, j, k,  n, atn, atk, xreps, yreps, rx, ry
CHARACTER(LEN=20)                             :: formg = '(I5, 2A5, I5, 3F8.3)'
CHARACTER(LEN=66)                             :: formt = '(I4, A1, A2, A1, I4, A1, A5, A1, A5, A1, I4, A1, F12.8, A1, F12.8)'
LOGICAL                                       :: Ms, IOs, mnE, wll
REAL, PARAMETER                               :: cp = 0.1 ! Conversion parameter - output position in nm, not Angstrom as in LAMMPS

charges = qH2O_O   ; Hcharges = qH2O_H
res     = 'SOL'  ; atom     = 'OW'     ; Hatom = 'HW'  ; name = 'O'  ; Msite = 'MW'

charges = qH2O_O ; Hcharges = qH2O_H ; masses = mO


Ms = ((COUNT(str%msite .GE. 0.05) .NE. 0))

IOs = .FALSE.
IF (PRESENT(ionO)) IOs = ionO

mnE = .FALSE.
IF (PRESENT(minE)) mnE = minE

wll = .TRUE.
IF (PRESENT(wall)) wll = wall

xreps = 3  ; yreps = 4
IF (PRESENT(repx)) xreps = repx
IF (PRESENT(repy)) yreps = repy

IF ((xreps .LT. 2) .OR. (yreps .LT. 4)) WRITE(*,*) "#INFO: Lack of repeats in short directions may lead to errors in GROMACS"

IF (Ms) charges(2) = 0.0

! Gromacs output gives a GRO file - atom coordinates - and a TOP file
! Define residue and atom names

IF (COUNT(str%val(:) .EQ. 0) .NE. 0) THEN
   IF((.NOT. PRESENT(mF)) .OR. (.NOT. PRESENT(qF))) THEN
      WRITE(*,*) 'ERROR: Structure contains F- ion but details have not been specified.'
      STOP
   END IF
   charges(0) = qF      ; masses(0) = mF 
   res(0)     = 'FION'  ; atom(0)   = 'F'  ; name(0) = 'F'
END IF

IF (COUNT(str%val(:) .EQ. 1) .NE. 0) THEN
   IF ((.NOT. PRESENT(qOH_O)) .OR. (.NOT. PRESENT(qOH_H))) THEN
      WRITE(*,*) 'ERROR: Structure contains OH- ion but details have not been specified.'
      STOP
   END IF
   charges(1) = qOH_O ; Hcharges(1) = qOH_H
   IF ((.NOT. Ms) .AND. (.NOT. IOs)) WRITE(*, *) 'WARNING: Using OH- or H3O+ ions in GROMACS without M-site or ionic oxygen.'
   IF (IOs) THEN
      atom(1) = 'IONO'  ; atom(3) = 'IONO'
   END IF
   res(1) = 'OH'
END IF

IF (COUNT(str%val(:) .EQ. 3) .NE. 0) THEN ! H3O+
   IF ((.NOT. PRESENT(qH3O_O)) .OR. (.NOT. PRESENT(qH3O_H))) THEN
      WRITE(*,*) 'ERROR: Structure contains H3O+ ion but details have not been specified.'
      STOP
   END IF
   charges(3) = qH3O_O ; Hcharges(3) = qH3O_H
   IF (COUNT(str%val(:) .EQ. 1) .EQ. 0) THEN
      ! No OH-, need to determine if we need to count for a new ionic O
      IF ((.NOT. Ms) .AND. (.NOT. IOs)) WRITE(*, *) 'WARNING: Using OH- or H3O+ ions in GROMACS without M-site or ionic oxygen.'
      IF (IOs) atom(3) = 'IONO'
   END IF
   res(3) = 'HHHO'
END IF

IF (COUNT(str%val(:) .EQ. 4) .NE. 0) THEN ! NH4+
   IF (((.NOT. PRESENT(qN)) .OR. (.NOT. PRESENT(qNH))) .OR. (.NOT. PRESENT(mN))) THEN
      WRITE(*,*) 'ERROR: Structure contains NH4+ ion but details have not been specified.'
      STOP
   END IF
   charges(4) = qN   ; Hcharges(4) = qNH  ; masses(4) = mN
   res(4) = 'AMMON'  ; atom(4) = 'N'      ; name(4) = 'N'   ; Hatom(4) = 'HA'
END IF

! ================ !
! WRITING GRO FILE !
! ================ !

! The GRO file gives the position of atoms

number(1) = '1'  ; number(2) = '2'  ; number(3) = '3'   ; number(4) = '4'

! Get positions
CALL OUTPUT_POSITIONS(hash, str, pot, pos, midO=.TRUE.)


WRITE(*,*) GET_DIPOLE(str, pos, Ms, qH2O_H, qH2O_O, qF, qN, qNH, qOH_O, qOH_H, qH3O_O, qH3O_H)

OPEN(fnum, file = TRIM(outfile)//'.gro', status = 'new')
WRITE(fnum, *) '; Produced by refine_bonds.f90 from atomic hash'
IF (wll) THEN
   WRITE(fnum, *) (str%natoms+8)*xreps*yreps
ELSE
   WRITE(fnum, *) str%natoms*xreps*yreps
END IF
! Contains residue [molecule] number, residue name, atom name, atom number, position *in nm not A* 

WRITE(*,*) '#INFO: writing ion information first - note that this is untested for multiple ions'
WRITE(*,*) '#INFO: explicitly replicating box by ', xreps, ' in x and ', yreps, 'in y direction'
IF (wll) WRITE(*,*) '#INFO: explicitly placing atomic sites for a charged wall under O positions at z = -3 A. ', &
                    'Note this assumes 8 atoms in basal plane'

n   = 1
atn = 1
atk = 1
DO i = 1, str%nmols
   IF (str%val(i) .NE. 2) THEN
      IF ((pos(3, atk) .LT. str%minz) .OR. (pos(3, atk) .GT. str%maxz)) &
           WRITE(*,*) 'WARNING: z outside box with coordinate ', pos(3, atk)
      DO rx = 0, xreps-1
         DO ry = 0, yreps-1
            WRITE(fnum, formg) atn, res(str%val(i)), TRIM(atom(str%val(i))), n, (pos(1, atk)+rx*str%x)*cp, &
                             (pos(2, atk)+ry*str%y)*cp, pos(3, atk)*cp
            n   = n + 1
            atk = atk + 1
            DO k = 1, str%val(i)
               IF ((pos(3, atk) .LT. str%minz) .OR. (pos(3, atk) .GT. str%maxz)) &
                    WRITE(*,*) 'WARNING: z outside box with coordinate ', pos(3, atk)            
               WRITE(fnum, formg) atn, res(str%val(i)), TRIM(Hatom(str%val(i)))//TRIM(number(k)), &
                               n, (pos(1, atk)+rx*str%x)*cp, (pos(2, atk)+ry*str%y)*cp, pos(3, atk)*cp
               n   = n + 1
               atk = atk + 1
            END DO
            atn = atn + 1
            atk = atk - 1 - str%val(i)
         END DO
      END DO
   END IF
   atk = atk + 3
   IF (Ms) atk = atk + 1
END DO
         
DO rx = 0, xreps - 1
   DO ry = 0, yreps - 1
      atk = 1
      DO i = 1, str%nmols
         IF (str%val(i) .EQ. 2) THEN
            IF ((pos(3, atk) .LT. str%minz) .OR. (pos(3, atk) .GT. str%maxz)) &
                 WRITE(*,*) 'WARNING: z outside box with coordinate ', pos(3, atk)
            WRITE(fnum, formg) atn, res(str%val(i)), TRIM(atom(str%val(i))), n, (pos(1, atk)+rx*str%x)*cp, &
                             (pos(2, atk)+ry*str%y)*cp, pos(3, atk)*cp
            n   = n   + 1
            atk = atk + 1 
            DO k = 1, str%val(i)
               IF ((pos(3, atk) .LT. str%minz) .OR. (pos(3, atk) .GT. str%maxz)) &
                    WRITE(*,*) 'WARNING: z outside box with coordinate ', pos(3, atk)
               WRITE(fnum, formg) atn, res(str%val(i)), TRIM(Hatom(str%val(i)))//TRIM(number(k)), &
                               n, (pos(1, atk)+rx*str%x)*cp, (pos(2, atk)+ry*str%y)*cp, pos(3, atk)*cp
               n   = n   + 1
               atk = atk + 1
            END DO
            IF (Ms) THEN
               WRITE(fnum, formg) atn, res(str%val(i)), TRIM(Msite(str%val(i))), n, &
                                  (pos(1, atk)+rx*str%x)*cp, (pos(2, atk)+ry*str%y)*cp, pos(3, atk)*cp
               n   = n   + 1 
               atk = atk + 1
            END IF
            atn = atn + 1            
         ELSE
            atk = atk + 1 + str%val(i)
         END IF
      END DO
   END DO
END DO


IF (wll) THEN
   DO rx = 0, xreps - 1
      DO ry = 0, yreps - 1
         atk = 1
         DO i = 1, str%nmols
            IF (pos(3, atk) .EQ. 0) THEN
               WRITE(fnum, formg) atn, 'WW   ', 'WA', n, (pos(1, atk)+rx*str%x)*cp, &
                    (pos(2, atk)+ry*str%y)*cp, -3*cp
               n   = n   + 1
               atn = atn + 1            
            END IF
            IF (str%val(i) .EQ. 2) THEN
               atk = atk + 4
            ELSE
               atk = atk + 1 + str%val(i)
            END IF
         END DO
      END DO
   END DO
END IF




! Write UC parameters - DO NOT KNOW HOW GROMACS WILL DEAL WITH Z OFFSET

IF (mnE) THEN
   WRITE(*, *) '#INFO: Adding vacuum gap to GROMACS input for energy minimisation.'
   WRITE(fnum, '(F12.8, A1, F12.8, A1, F12.8)') xreps*str%x*cp, ' ', yreps*str%y*cp, ' ',  255*cp
ELSE
   WRITE(fnum, '(F12.8, A1, F12.8, A1, F12.8)') xreps*str%x*cp, ' ', yreps*str%y*cp, ' ', (str%maxz - str%minz)*cp
END IF


CLOSE(fnum)

!fnum = fnum + 5
!
!! ================ !
!! WRITING TOP FILE ! 
!! ================ !
!
!OPEN(fnum, file = TRIM(outfile)//'.top', status = 'new')
!
!WRITE(fnum, *) ';'
!WRITE(fnum, *) '; Topology file written by refine_bonds.f90 from atomic hash'
!WRITE(fnum, *) ';'
!WRITE(fnum, *) '; Force-field files to be included'
!WRITE(fnum, *) '#include '
!WRITE(fnum, *)
!
!WRITE(fnum, *) '[ moleculetype ]'
!WRITE(fnum, *) '; name nrexcl'
!WRITE(fnum, *) 'SOL 1'
!IF (COUNT(str%val(:) .EQ. 0) .NE. 0) WRITE(fnum, *) 'FION 1'
!IF (COUNT(str%val(:) .EQ. 1) .NE. 0) WRITE(fnum, *) 'OH 1'
!IF (COUNT(str%val(:) .EQ. 3) .NE. 0) WRITE(fnum, *) 'HHHO 1'
!IF (COUNT(str%val(:) .EQ. 4) .NE. 0) WRITE(fnum, *) 'AMMON 1'
!WRITE(fnum, *) 
!
!
!WRITE(fnum, *) '[ atoms ]'
!WRITE(fnum, *) '; counter atomtype residuenum residuename atomname chargegroup charge mass'
!n = 1
!DO i = 1, str%nmols
!   WRITE(fnum, formt) n, ' ', TRIM(name(str%val(i))), ' ', i, ' ', res(str%val(i)), ' ', TRIM(atom(str%val(i))), ' ', n, ' ', & 
!                     charges(str%val(i)), ' ', masses(str%val(i))
!   n = n + 1
!   DO k = 1, str%val(i)
!      WRITE(fnum, formt) n, ' ', 'H', ' ', i, ' ', res(str%val(i)), ' ', TRIM(Hatom(str%val(i)))//TRIM(number(k)), ' ', &
!                         n, ' ', Hcharges(str%val(i)), ' ', mH
!      n = n + 1
!   END DO
!   IF ((Ms) .AND. (str%val(i) .EQ. 2)) THEN
!      WRITE(fnum, formt) n, ' ', 'M', ' ', i, ' ', res(str%val(i)), ' ', TRIM(Msite(str%val(i))), ' ', n, ' ', qH2O_O, ' ', 0.0
!      n = n + 1 
!   END IF
!END DO
!WRITE(fnum, *)
!
!! Not entirely sure about the next bit
!
!WRITE(fnum, *) '[ bonds ]'
!WRITE(fnum, *) '; atom1 atom2'
!atn = 1
!DO i = 1, str%nmols
!   IF (str%val(i) .EQ. 0) THEN
!      ! F- has no bonds, but still need to skip over to ensure correct numbering of atoms
!      atn = atn + 1
!      CYCLE
!   END IF
!   DO k = 1, str%val(i)
!      WRITE(fnum, '(I4, A1, I4)') atn, ' ', atn + k
!   END DO
!   IF ((Ms) .AND. (str%val(i) .EQ. 2)) THEN
!      WRITE(fnum, '(I4, A1, I4)') atn, ' ', atn + 3
!      atn = atn + 2 + str%val(i) 
!   ELSE
!      atn = atn + 1 + str%val(i)
!   END IF
!END DO
!WRITE(fnum, *)
!
!! Write angles information
!WRITE(fnum, *) '[ Angles ]'
!WRITE(fnum, *) '; atom1 atom2 atom3'
!atn = 1
!atk = 2 
!DO i = 1, str%nmols
!   IF (str%val(i) .LT. 2) THEN
!      ! OH- and F- have no angles, but still need to ensure atoms skipped
!      atn = atn + 1 + str%val(i)
!      atk = atn + 1
!      CYCLE
!   END IF
!   DO k=str%val(i) - 1, 1, -1
!      DO j = 1, k
!         WRITE(fnum, '(I4, A1, I4, A1, I4)') atk, ' ', atn, ' ', atk + j 
!         n = n + 1
!      END DO
!      atk = atk + 1 ! To ensure correct progression through atoms
!   END DO
!   IF ((Ms) .AND. (str%val(i) .EQ. 2)) THEN
!      WRITE(fnum, '(I4, A1, I4, A1, I4)') atn + 1, ' ', atn, ' ', atn + 3 
!      n = n + 1 
!      WRITE(fnum, '(I4, A1, I4, A1, I4)') atn + 2, ' ', atn, ' ', atn + 3
!      n   = n + 1 
!      atn = atn + 2 + str%val(i) 
!   ELSE
!      atn = atn + 1 + str%val(i)
!   END IF
!   atk = atn + 1
!END DO
!
!CLOSE(fnum)

END SUBROUTINE WRITE_GROMACS
! ===========================================================================================
FUNCTION CROSS_PRODUCT(a, b) RESULT(cross)
  IMPLICIT NONE
  REAL, DIMENSION(3), INTENT(IN) :: a, b
  REAL, DIMENSION(3)             :: cross

  cross(1) =   a(2)*b(3) - a(3)*b(2)
  cross(2) = -(a(1)*b(3) - a(3)*b(1))
  cross(3) =   a(1)*b(2) - a(2)*b(1)
  cross(:) =   cross/SQRT(DOT_PRODUCT(cross, cross))
  
END FUNCTION CROSS_PRODUCT
! ===========================================================================================
FUNCTION GET_DIPOLE(str, pos, Ms, qH2O_H, qH2O_O, qF, qN, qNH, qOH_O, qOH_H, qH3O_O, qH3O_H) RESULT(dipole)
  IMPLICIT NONE

  TYPE(STRUC),                       INTENT(INOUT) :: str
  REAL, DIMENSION(1:3,1:str%natoms), INTENT(IN)    :: pos
  LOGICAL,                           INTENT(IN)    :: Ms
  REAL,                              INTENT(IN)    :: qH2O_O, qH2O_H
  REAL, OPTIONAL,                    INTENT(IN)    :: qF, qN, qNH, qOH_O, qOH_H, qH3O_O, qH3O_H
  REAL, DIMENSION(1:3)                             :: dipole

  INTEGER              :: i, j, n
  REAL, DIMENSION(0:4) :: charges, Hcharges

  charges = qH2O_O  ; Hcharges = qH2O_H

  IF (COUNT(str%val(:) .EQ. 0) .NE. 0) THEN ! F-
     IF (.NOT. PRESENT(qF)) THEN
        WRITE(*,*) 'ERROR: Structure contains F- ion but details have not been specified.'
        STOP
     END IF
     charges(0) = qF
  END IF

  IF (COUNT(str%val(:) .EQ. 1) .NE. 0) THEN ! OH-
     IF ((.NOT. PRESENT(qOH_O)) .OR. (.NOT. PRESENT(qOH_H))) THEN
        WRITE(*,*) 'ERROR: Structure contains OH- ion but details have not been specified.'
        STOP
     END IF
     charges(1) = qOH_O  ; Hcharges(1) = qOH_H
  END IF

  IF (COUNT(str%val(:) .EQ. 3) .NE. 0) THEN ! OH-
     IF ((.NOT. PRESENT(qH3O_O)) .OR. (.NOT. PRESENT(qH3O_H))) THEN
        WRITE(*,*) 'ERROR: Structure contains H3O+ ion but details have not been specified.'
        STOP
     END IF
     charges(3) = qH3O_O  ; Hcharges(3) = qH3O_H
  END IF

  IF (COUNT(str%val(:) .EQ. 4) .NE. 0) THEN ! NH4+
     IF ((.NOT. PRESENT(qN)) .OR. (.NOT. PRESENT(qNH))) THEN
        WRITE(*,*) 'ERROR: Structure contains NH4+ ion but details have not been specified.'
        STOP
     END IF
     charges(4) = qN  ; Hcharges(4) = qNH
  END IF

  dipole = 0.0
  n = 1

  DO i = 1, str%nmols
     IF (.NOT. ((Ms) .AND. (str%val(i) .EQ. 2))) dipole = dipole + charges(str%val(i))*pos(:, n)
     n      = n + 1 
     DO j = 1, str%val(i)
        dipole = dipole + Hcharges(str%val(i))*pos(:, n)
        n      = n + 1
     END DO
     IF ((Ms) .AND. (str%val(i) .EQ. 2)) THEN
        dipole = dipole + charges(str%val(i))*pos(:, n)
        n      = n + 1
     END IF
  END DO

END FUNCTION GET_DIPOLE
! ===========================================================================================
END MODULE OUTPUT_POTENTIAL
