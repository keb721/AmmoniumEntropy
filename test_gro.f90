PROGRAM MAIN
  
  USE OUTPUT_POTENTIAL
  IMPLICIT NONE

  INTEGER, PARAMETER         :: recordlen = 130, hashes = 1 ! record length is +2 for formatting?
  CHARACTER(LEN=6)           :: recform ='(A130)'
  INTEGER, DIMENSION(hashes) :: hashpos = (/ 2102 /)!0, 150, 500, 1000, 2500, 5000, 10000, 25000, 50000 /)
  INTEGER, DIMENSION(1, 2)   :: ions
  TYPE(POTENTIAL)            :: pot
  TYPE(STRUC)                :: struct
  REAL*8                     :: H2Oang, H3Oang, NH4ang, NH4len, H3Olen, H2Olen, OHlen, dOM
  REAL                       :: mH, mO, mN, mF 
  REAL                       :: qF, qOH_O, qOH_H, qH2O_O, qH2O_H, qH3O_O, qH3O_H, qN, qHN
  CHARACTER(LEN=3)           :: outform  ! 'lmp' or 'dcd' or 'xyz' or 'gro' - Whether to output as a GROMACS/LAMMPS input or an xyz/dcd file
  CHARACTER(LEN=500)         :: icePos, hashesFile, outputPos, hash, dump
  CHARACTER(LEN=7)           :: outnum
 
  REAL, DIMENSION(:,:), ALLOCATABLE :: pos
  INTEGER                           :: fnum, atom, j, i, k, xreps, yreps
  REAL                              :: x, y, z, zoffset
  LOGICAL                           :: incMs = .FALSE., ionO = .FALSE., minE = .FALSE., reps=.FALSE.

  ! H2O data comes from TIP4P/Ice (LAMMPS webpage)
  ! NH4+ and OH- data comes from CHARMM 36
  ! F- data comes from DOI:10.1021/acs.jctc.1c00550
  ! H3O+ data comes from DOI:10.3389/fchem.2019.00439

  H2Oang = 104.52
  H3Oang = 111.4
  NH4ang = 109.5

  OHlen  = 0.97     !FFbonded L 1773
  H2Olen = 0.9572
  H3Olen = 0.98
  NH4len = 1.04
  dOM    = 0.1577
  
  mH     = 1.008
  mO     = 15.9994
  mN     = 14.007
  mF     = 18.998

  qH2O_O = -1.1794
  qN     = -0.24
  qH2O_H = 0.5897
  qHN    = 0.31
  qF     = -1.0  ! Table 2
  qH3O_H = 0.8
  qH3O_O = -1.4
  qOH_O  = -1.320
  qOH_H  = 0.320


  ions(1, 1) = 3  ; ions(1, 2) = 4

  icePos     = "ordered_o211.xyz"
  hashesFile = "all_strucs_c211_N3_b8.txt" 
  
  !CALL GET_COMMAND_ARGUMENT(1, icePos) ! File with original positions of the atoms
  !CALL GET_COMMAND_ARGUMENT(2, hashesFile) ! File with the hashes 
  CALL GET_COMMAND_ARGUMENT(1, outputPos) ! Prefix of the file to save positions
  !CALL GET_COMMAND_ARGUMENT(4, outform) ! Whether to save as .xyz or lammps data file

  fnum=77
  
  READ(outputPos,*)j

  hashpos(1) = j


  OPEN(fnum, file = TRIM(icePos), status = 'old')  
  READ(fnum,*) ; READ(fnum,*) ; READ(fnum, *) atom ; READ(fnum, *) ! Two comment lines, then number of atoms, then periodicity

  ! Allocate positions array
  ALLOCATE(pos(1:3, 1:atom))
  
  DO j = 1, atom
     READ(fnum,*) dump, pos(1, j), pos(2, j), pos(3, j)
  END DO
  READ(fnum, *) dump, x ; READ(fnum, *) dump, dump, y ; READ(fnum, *) dump, dump, dump, z ; READ(fnum, *) dump, dump, dump, zoffset
  CLOSE(fnum)
  

  CALL SETUP_POTENTIAL(H2Oang, H2Olen, pot, NH4ang, NH4len, H3Oang, H3Olen, OHlen, dOM)

  CALL SETUP_STRUCT(pos, pot, INT(atom/3), ions, x, y, z, zoffset, struct, .TRUE.)

  fnum = fnum + 10
  OPEN(fnum, file = TRIM(hashesFile), status = 'old', form = 'FORMATTED', access = 'DIRECT', recl = recordlen)

  DO j = 1, hashes
     READ(fnum, fmt = recform, rec = hashpos(j)) hash
     hash = hash(2:)         ! For some reason, the hash starts with a leading space

     write(*,*) hash
     WRITE(outnum, '(I0)') hashpos(j)

     CALL WRITE_GROMACS(hash, struct, pot, 'input_'//TRIM(outnum)//'_', mH, mO, qH2O_H, qH2O_O, &
                       mN, mF, qF, qN, qHN, qOH_O, qOH_H, qH3O_O, qH3O_H, minE=.TRUE.)

  END DO

  CLOSE(fnum)



  CALL DESTROY_STRUCT(struct)

  DEALLOCATE(pos)

END PROGRAM MAIN
