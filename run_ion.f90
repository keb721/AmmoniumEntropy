program main

  use hbond_module
  
  implicit none

  integer, parameter :: atoms=96
  type(system) :: str
  real, dimension(:,:), allocatable :: pos
  real :: x, y, z, z_off
  integer :: maxt, unchanged, i, strucs, k, stored, maxsame
  integer :: atom, fnum, ifail, j, innum, strnum
  character(len=100) :: dump, ice_file, inputt, output_fl, basenum
  character(len=4*atoms/3) :: hashout, hashin
  character(len=4*atoms/3), dimension(:), allocatable :: hashes

  integer, dimension(1,2) :: ions  ! Nx2 array where first column gives atomic location of ions and column 2 gives number of donated H bonds
  integer :: base != 7 ! Number of connections of base to surface
  logical :: layer = .TRUE. ! Bulk ice or a layer
  
  ions(1, 1) = 2 
  ions(1, 2) = 4 ! Ammonium (4 donated hydrogen bonds) in position of atom 2 


  maxt = 50000000 ! Maximum number of moves to try - will depend on system
  maxsame = 250 ! If can't generate a new structure after this many moves, exit

  CALL GET_COMMAND_ARGUMENT(1, ice_file)
  CALL GET_COMMAND_ARGUMENT(2, output_fl) ! Suffix of output files
  CALL GET_COMMAND_ARGUMENT(3, basenum)
  
  READ(basenum, *) base

  write(*,*) "#INFO: Reading atomic positions from ", TRIM(ice_file), " for ", atoms/3, " molecules"
  write(*,*) "#INFO: Writing to file strucs",output_fl
  write(*,*) "#INFO: Running for ", maxt, " moves or until structure unchanged after ", maxsame, "moves"

  fnum=307
  innum=97
  strnum=427
  
  open(fnum, file=TRIM(ice_file), status='old')  
  read(fnum,*) ; read(fnum,*); read(fnum, *) atom ; read(fnum, *) ! Two comment lines, then number of atoms, then periodicity

  if (atom .NE. atoms) then
     write(*,*) "Incorrect number of atoms - check fortran parameter"
     stop
  end if

  ! Allocate positions array
  allocate(pos(1:3, 1:atoms))
  
  do j=1, atoms
     read(fnum,*) dump, pos(1, j), pos(2, j), pos(3, j)
  end do
  read(fnum, *) dump, x ; read(fnum, *) dump, dump, y ; read(fnum, *) dump, dump, dump, z  ; read(fnum, *) dump, dump, dump, z_off
  close(fnum)
  
  
  CALL createHbonds(pos, x, y, z, z_off, atoms/3, str, layer, base, ions)
  deallocate(pos)
  call output(str, hashout)

  if (layer) then
     allocate(hashes(1:maxt+1))
  else
     allocate(hashes(1:2*(maxt+1)))
  end if
  strucs    = 1
  hashes = hashout 
  open(innum,  file="increases"//TRIM(output_fl), status='new')
  open(strnum, file="strucs"//TRIM(output_fl), status='new')
  write(strnum, *) hashout
  
  if (layer) then
     stored = 2
  else
     call equivalent(str, hashout)
     hashes(2) = hashout
     stored = 3
  end if
  
  unchanged = 0


  do i=1, maxt
     call changeHbonds(str, ifail, hashout)
     if(COUNT(hashes(1:stored) .EQ. hashout) .EQ. 0) then
        ! Not currently in the list of hashes - add it and its inverse, then check for isomorphisms
        hashes(stored) = hashout
        write(strnum, *) hashout
        write(innum, *) i
        stored = stored + 1
        unchanged = 0
        if (.NOT. layer) then
           call equivalent(str,  hashout)
           hashes(stored) = hashout
           stored = stored + 1
        end if
     else
        unchanged = unchanged + 1
     end if
     
     if (unchanged .EQ. maxsame) then
        ! Assume structure saturated
        write(strnum, *) "#INFO: After ", i, " moves, no change in struture within timeframe given"
        exit
     end if
  end do


  write(strnum, *) "#INFO: Total number of structures found is ", stored - 1
  close(innum)
  close(strnum)
  call destroyHbonds(str)
 
  deallocate(hashes)

end program main
