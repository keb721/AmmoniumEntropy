! ===========================================================================================
! A module that contains functions that change the hydrogen bonding topology
! using the algorithm due to Rick:
! 
! To use in a fortran program insert
! use hbond_module, only setupHbonds,changeHbonds,destroyHbonds
! 
! Then call setupHbonds to allocate all the various arrays
! call changeHbonds whenever you want to change the hydrogen bond topology 
! call destoryHbonds at the end of your program 
! 
! WARNING - The program assumes that you are using a TIP4P potential
! 
! ===========================================================================================
module hbond_module 
implicit none

type system
  private
  integer :: nmol, plane_sz ! Plane_sz is 0 if unfixed, else gives the number of atoms in basal plane
  integer,allocatable,dimension(:,:) :: connection, base, ghost  ! Base gives the atoms in basal plane (i.e. the base). Ghost deals with the ghost atoms and their connections to real atoms
  integer,allocatable,dimension(:) :: bonds, reverses, stoich
  logical :: move_base ! If any of the atoms in the base can move (prevent for loop if all of base is fixed)
end type system

type node
   private
   integer :: data = 0
   type(node), pointer :: next=>NULL()
end type node

type FIFO_list
   private
   type(node), pointer :: head=>NULL(), tail=>NULL()
end type FIFO_list

! ==========================================================================================
contains
subroutine setupHbonds( nmol, structure, ions )
use random
implicit none
integer,intent(in) :: nmol
integer, dimension(:,:), intent(in) :: ions
!logical,intent(in) :: base ! If there is a fixed base
type(system),intent(out) :: structure
real*8 :: rand1
integer :: num_conns, i
integer, dimension(2) :: num_ions


rand1=random_uniform_random()

! Initialize jumble array - used in move
jumble=(/1,2,3,4/)

structure%nmol= nmol
num_conns = SUM(ions(:, 2))

allocate( structure%connection(1:structure%nmol, 1:4) )
allocate( structure%bonds(1:4*structure%nmol), structure%reverses(1:4*structure%nmol), structure%stoich(1:structure%nmol))

if (ABS(num_conns) .GT. (0.5*nmol)) then
   ! Maximum/minimum deviation of bonds to create a valid structure
   write(*,*) "ERROR: Too many ions to place into the structure."
   stop
end if

num_ions = SHAPE(ions)

structure%stoich(:) = 2

do i=1, num_ions(1)
   structure%stoich(ions(i, 1)) = ions(i, 2) 
end do

end subroutine setupHbonds

! Bookeeping stuff for changing hydrogen bonding topology
! ==========================================================================================
subroutine destroyHbonds(structure)
implicit none
type(system),intent(inout) :: structure

deallocate( structure%bonds,structure%connection,structure%ghost  )
deallocate( structure%reverses,structure%stoich,structure%base )

end subroutine destroyHbonds

! ==========================================================================================
subroutine createHbonds(pos, x, y, z, z_offset, mols, str, layer, base, ions)
use random, only : random_uniform_random
implicit none
type(system),intent(out) :: str
real, intent(in) :: x, y, z, z_offset
integer, intent(in) :: mols, base
integer, dimension(:,:), intent(in) :: ions
real,dimension(:,:),intent(inout) :: pos
logical, intent(in)      :: layer

! Str is the structure to create
! pos is the GenIce structure (no Ns, no constraints)
! x, y, z are the lattice dimensions 
! ions is the index of oxygens to replace with ions with specified stoichiometry - job of user to check constraints etc.
! mol is total number of molecules
! Layer is if we want a layer of ice (as opposed to bulk ice)
! Base is the fixed number of bonds in the basal plane (THIS CANNOT CHANGE WITHIN THE SIMULATION)

integer :: i, j, k, O, Hpos, coord, allowed, Ni, ci
real, dimension(3) :: tmpvect
real :: invx, invy, invz, mindist, dist, minz, diff
integer, dimension(:, :), allocatable :: temp
integer, dimension(1) :: lenbas, location
integer, dimension(2) :: num_ions

integer, dimension(:), allocatable :: basatoms, basatomstmp, tmp

invx = 1/x ; invy = 1/y ; invz = 1/(z-z_offset)

CALL setupHbonds(mols, str,  ions)

diff = 0.15 ! How much leeway to allow for differences in y and for x/z displacements

! Ok, first thing I want to do is set up an array detailing how all of the Os are connected together
! (using genice H information - then this can be discarded) - can toggle periodicity later
! I think because we discard H positions, we can ignore knowing which image things are in

str%connection = -1
str%reverses   = -1

do i=1, str%nmol
   do k=1, 2
   
      O = i  ! The O this hydrogen is covalently bonded to
      ! Nw need to determine the next closest O (i.e. this H is H bonded to)

      mindist = 500.0 ! Set a high minimum distance for finding next closest oxygen
      Hpos = 3*O-2+k   ! Position of H in pos array
      
      do j=1, str%nmol

         if (j .EQ. O) CYCLE ! Want H bond not covalent bond
      
         ! ASSUMING POS IS (3, NMOL*3) - Be careful here!

         tmpvect = pos(:, Hpos) - pos(:,(3*j-2)) ! This is the position of O in the pos array
   
         ! Check images
         tmpvect(1) = tmpvect(1) - x*anint(tmpvect(1)*invx)
         tmpvect(2) = tmpvect(2) - y*anint(tmpvect(2)*invy)
         tmpvect(3) = tmpvect(3) - (z-z_offset)*anint(tmpvect(3)*invz)
      
         dist   = sqrt(tmpvect(1)**2 + tmpvect(2)**2 + tmpvect(3)**2)
         if (dist.lt.mindist) then
            str%connection(O, k) = j         
            mindist = dist
         end if
      end do

      ! Now str%connection(O, 1) and str%connection(O, 2) are ones where O is cov bonded to H
      ! where H is then H bonded to j

      if (str%connection(str%connection(O, k), 3) .NE. -1) then
         str%connection(str%connection(O, k), 4) = O
         str%reverses(4*str%connection(O,k)) = 4*O-4+k
         str%reverses(4*O-4+K) = 4*str%connection(O,k)
      else
         str%connection(str%connection(O, k), 3) = O
         str%reverses(4*str%connection(O,k)-1) = 4*O-4+k
         str%reverses(4*o-4+K) = 4*str%connection(O,k)-1
      end if
   end do
end do

! Now we have an array giving the connections between all the Os, and the direction of these connections

if (COUNT(str%reverses .EQ. -1) .NE. 0) then
   write(*,*) 'Error creating bonds'
   stop
end if

do i=1, str%nmol*4
   if (COUNT(str%reverses .EQ. i) .NE. 1) then
      write(*,*) 'Error creating bonds'
      stop
   end if
end do

str%bonds(1::4) = 1 ! 4O-3 and 3O-2 correspond to covalent bonds (IN)  wrt THIS O
str%bonds(2::4) = 1 ! 4O-3 and 3O-2 correspond to covalent bonds (IN)  wrt THIS O
str%bonds(3::4) = 0 ! 4O-1 and 4O correspond to hydrogen bonds (OUT) wrt THIS O
str%bonds(4::4) = 0 ! 4O-1 and 4O correspond to hydrogen bonds (OUT)  wrt THIS O


if (layer) then

   write(*,*) "#INFO: Layer structure selected. Input is required to be ordered ice structure"

   ! Now we have all the connections, we want to consider the position of the structure layer
   ! In genice, the "straight" connection is in the y direction
   ! Want the base layer to be the lowest (z) atoms where one of the connections is -z

   ! Incorporate the stoichiometry of ions

   str%stoich(:) = 2 ! Default is H20 (2 in 2 out)
   num_ions = SHAPE(ions)
   write(*,*) "#INFO: Incorporating ", num_ions(1), " ions"
   do k=1, num_ions(1)
      str%stoich(ions(k, 1)) = ions(k, 2)
   end do

   ! If there are any anti-ammonium in the base (i.e. 0 covalent bonds to H) need to make sure that it is disconnected at the start

   
   allocate(temp(1:int(0.5*str%nmol), 2))   ! We know half up half down
   j = 1

   minz = 50.0
   
   do i=1, str%nmol
      O = 3*i - 2   ! Oxygen position
      do k = 1, 4
         tmpvect =  pos(:,3*str%connection(i, k)-2) - pos(:, O) ! This is the position of O in the pos array
      
         ! Check images
         tmpvect(1) = tmpvect(1) - x*anint(tmpvect(1)*invx)
         tmpvect(2) = tmpvect(2) - y*anint(tmpvect(2)*invy)
         tmpvect(3) = tmpvect(3) - (z-z_offset)*anint(tmpvect(3)*invz)
      
         if ((ABS(tmpvect(1)) .LT. diff) .AND. (ABS(tmpvect(2)) .LT. diff)) then
            ! This connection is in +/- z
            if (tmpvect(3) .LT. 0) then
               ! This is a connection in -z - guaranteed to be a hydrogen bond
               temp(j, 1) = i
               temp(j, 2) = k
               j = j + 1
               minz = merge(pos(3, O), minz, (pos(3, O) .LT. minz))
            end if
         end if
      end do
   end do

!   write(*,*) minz

   ! Temp now contains all atoms and connections pointing down - now find lowest

   do i=1, int(0.5*str%nmol)
      if (pos(3, 3*temp(i, 1)-2) - minz .GT. diff) temp(i, 1) = 0
   end do
   
   str%plane_sz = COUNT(temp(:, 1)  .NE. 0) ! Number of atoms in the basal plane
   
   allocate(str%base(1:str%plane_sz, 2)) ! IDs of atoms in basal plane and bonds to "surface"
   allocate(str%ghost(1:str%plane_sz, 2)) ! Allocate the ghost atom arry - know we need one ghost atom per real atom in the base

   str%move_base = .TRUE. 

   write(*,*) "#INFO: The basal plane contains ", str%plane_sz, "atoms. ", &
        base, " of these are connected to the surface at all time."

   if (str%plane_sz .EQ. base) then
      str%move_base = .FALSE. ! Can't change basal layer
   else if (base .GT. str%plane_sz) then
      write(*,*) "ERROR: Too many atoms to be fixed to the surface."
      stop
   end if

   j = 1
   do i=1, int(0.5*str%nmol)
      if (temp(i, 1) .NE. 0) then
         str%base(j, :) = temp(i, :)
         ! Now need to get rid of y periodicity. Luckily, we know all of the affected atoms (base layer and connections)
         ! Introduce ficticious atoms, either fixed (must be covalent) (-1) or free (0) (can be the position of an H atom, or a fake H bond)
         str%connection(str%base(j, 1), str%base(j, 2)) = -1
         str%bonds(4*str%base(j, 1)-4+str%base(j, 2)) = 1
         k = str%reverses(4*str%base(j, 1)-4+str%base(j, 2))
         ! Don't need to take account of reverses which have already been deallocated (as only ever straight down)
         if (MOD(k, 4) .EQ. 0) then
            ! Add these connections to the ghost atom array (then add ghost atoms to real array later on)
            str%connection(FLOOR(0.25*k), 4) = str%nmol+j
            str%ghost(j, 1) = FLOOR(0.25*k)
            str%ghost(j, 2) = 4
         else
            str%connection(FLOOR(0.25*k)+1, MOD(k, 4)) = str%nmol+j
            str%ghost(j, 1) = FLOOR(0.25*k)+1
            str%ghost(j, 2) = MOD(k, 4)
         end if
         str%reverses(k) = 0
         str%reverses(4*str%base(j, 1)-4+str%base(j, 2)) = 0
         j = j+1
      end if
   end do

   deallocate(temp)

   allocate(tmp(1:str%plane_sz)) ! Give the layer 2 atoms
   
   tmp(:) = 0
   k = 1
   
   ! Loop over all base layer connections
   do i=1, str%plane_sz
      do j=1, 4
         if ((str%connection(str%base(i, 1), j) .NE. -1) .AND. & ! Not the base
            (COUNT(tmp(:) .EQ. str%connection(str%base(i, 1), j)) .EQ. 0)) then   ! Not already in tmp
            tmp(k) = str%connection(str%base(i, 1), j)
            k = k + 1
         end if
      end do
   end do

   Ni    = 0 ! number of ammonia
   ci    = 0 ! number of anti-ammonia
   coord = 0 ! stoichiometry of basal plane

   do i=1, str%plane_sz
      coord = coord + str%stoich(str%base(i, 1))
      if (str%stoich(str%base(i, 1)) .EQ. 4) Ni = Ni + 1
      if (str%stoich(str%base(i, 1)) .EQ. 0) ci = ci + 1
   end do

   if ((Ni .GT. base) .OR. (ci .GT. (str%plane_sz - base))) then
      write(*,*) "ERROR: Too many ions in basal plane"
      stop
   end if
   
   
   do i=1, str%plane_sz
      coord = coord +str%stoich(tmp(i))
   end do
   

   if (coord .GT. (base + 4*str%plane_sz)) then
      write(*,*) "ERROR: Too many ions in basal plane"
      stop
   end if

   deallocate(tmp)


   if (str%move_base) then
      ! There are fewer connections to the base than there are base atoms

      ! Now make sure we have the right number of connections to the surface
      j = 0
      ! First, any ammonium *MUST* be connected
      
      allocate(basatoms(1:str%plane_sz)) ! These are the possible atoms we can disconnect from the base
      do i=1, str%plane_sz
         if (str%stoich(str%base(i, 1)) .NE. 4) then
            ! This isn't an ammonium  atom, we can have it as disconnecting atoms
            j = j+1 
            basatoms(j) = str%base(i, 1)
         end if
      end do
   
     if (j .NE. str%plane_sz) then
         allocate(basatomstmp(1:j))
         basatomstmp(1:j) = basatoms(1:j)
         deallocate(basatoms)
         allocate(basatoms(1:j))
         basatoms(1:j) = basatomstmp(1:j)
         deallocate(basatomstmp)
      end if

      j = str%plane_sz
      ! Now, any anti-ammonium *MUST* be disconnected
      do k=1, str%plane_sz
         if (str%stoich(str%base(k, 1)) .EQ. 0) then
            str%bonds(4*str%base(k, 1)-4+str%base(k, 2)) = 0
            j = j-1
            lenbas = SHAPE(basatoms)
            allocate(basatomstmp(1:lenbas(1)-1))
            allowed = 1
            do i=1, lenbas(1)
               if (basatoms(i) .NE. str%base(k, 1)) then
                  basatomstmp(allowed) = basatoms(i)
                  allowed = allowed +1
               end if
            end do
            deallocate(basatoms)
            allocate(basatoms(1:lenbas(1)))
            basatoms = basatomstmp
            deallocate(basatomstmp)
         end if
      end do

      ! Can then proceed as normal

      do while (j .GT. base)
         lenbas = SHAPE(basatoms)
         k = (lenbas(1)*random_uniform_random())+1
         i = basatoms(k)
         location = MAXLOC(str%base(:,1), mask=(str%base(:,1) .EQ. i))
         i = location(1)
         ! Know this is an allowed atom
         str%bonds(4*str%base(i, 1)-4+str%base(i, 2)) = 0
         j = j - 1
         allocate(basatomstmp(1:lenbas(1)-1))
         allowed = 1
         do i=1, lenbas(1)
            if (i .NE. k) then
               basatomstmp(allowed) = basatoms(i)
               allowed = allowed +1
            end if
         end do
         deallocate(basatoms)
         allocate(basatoms(1:lenbas(1)))
         basatoms = basatomstmp
         deallocate(basatomstmp)
      end do
      
      deallocate(basatoms)
   end if


   ! In this new formalism, the ion positions are given - we just need to make sure they are allowed. 
   ! Probably want to give the connections too (i.e. information about where they are)

   num_ions = SHAPE(ions)

   do k=1, num_ions(1)
      Ni = str%stoich(ions(k, 1))
      ! Position of ion now chosen by user input - only need to return an error message if too close together
      ! First, any (anti) ammonium connected to another (anti) ammonium
      do i=1, 4
         if (str%connection(ions(k, 1), i) .EQ. -1) cycle
         if (str%connection(ions(k, 1), i) .GT. str%nmol) cycle
         if (((Ni .EQ. 0) .OR. (Ni .EQ. 4)) .AND.(Ni .EQ. str%stoich(str%connection(ions(k, 1), i)))) then
            write(*,*) "ERROR: Two connected NH4s"
            stop
         end if
         ! Then ensure that there aren't too many bonds for ions

         ci = Ni + str%stoich(str%connection(ions(k, 1), i))
         do j=1,4
            if (str%connection(str%connection(ions(k, 1), i), j) .EQ. ions(k, 1)) cycle
            if (str%connection(str%connection(ions(k, 1), i), j) .EQ. -1) cycle
            if (str%connection(str%connection(ions(k, 1), i), j) .GT. str%nmol) cycle
            ci = Ni
            ci = ci + str%stoich(str%connection(ions(k, 1), i))
            ci = ci + str%stoich(str%connection(str%connection(ions(k, 1), i), j))
         end do
         if (ci .GT. 10) then
            write(*,*) "ERROR: Too many connected ions"
            stop
         end if
      end do
   end do


   ! NOTE: This checks limited compatability ONLY. Other issues may arise but not be flagged.
   !       However this will simply cause the program to hang rather than give incorrect results
   
   CALL make_valid(str)
   write(*,*) "#INFO: Structure initialised"

   coord=test(str)
   if (coord.ne.0) then
      write(0,*)"Configuration input was not a valid hydrogen bonding configuration"
      stop
   end if
   
else
   allocate(str%base(0,0))
   allocate(str%ghost(0,0))
   str%base = 0
   num_ions = SHAPE(ions)
   if (num_ions(1) .NE. 0) then
      write(*,*) "ERROR: Code currently not written to allow ions in a bulk structure."
      stop
   end if
   write(*,*) "#INFO: Bulk structure modelled"
end if

end subroutine createHbonds

! ==========================================================================================

subroutine changeHbonds( structure,ifail, hash)
implicit none
type(system),intent(inout) :: structure
integer,intent(out) :: ifail 
integer :: coord
integer, dimension(1) :: ntg
character(len=4*structure%nmol), intent(out) :: hash

ifail = 0

coord=test(structure)
if (coord.ne.0) then
   write(0,*)"Configuration input was not a valid hydrogen bonding configuration"
   ifail = -2
   return
end if

call move(structure, ifail)

do while (ifail .NE. 0)
   call move(structure, ifail)
end do

coord=test(structure)
if (coord.ne.0) then
   write(0,*)"Configuration output was not a valid hydrogen bonding configuration"
   ifail = -2
   return
end if

!write(*,*) "success"

call output(structure, hash)

end subroutine changeHbonds

! ==========================================================================================

function test(str)
  implicit none
  type(system), intent(inout) :: str
  integer :: i, test

  test = 0
  
  do i=1, str%nmol
     if (COUNT(str%bonds(4*i-3:4*i) .EQ. 1 ) .NE. str%stoich(i)) then
        test = -2
        write(*,*) "Atom ", i, str%bonds(4*i-3:4*i)
        exit
     end if
  end do
end function test


! ==========================================================================================
subroutine move(str, ifail)
! Subroutine to move to a new hydrogen bonded configuration.
! Works by randomly picking a molecule then moving around a 
! closed loop swapping the orientations of bonds.
use random, only : random_uniform_random,jumble,random_seq
implicit none
type(system) :: str
integer, intent(out) :: ifail
integer,dimension(1:str%nmol*10) :: mol_l,bond_l, extra_bond  ! Allow long loops
integer :: rand1, rand2, i, j,k,coord,l_start,l_end, bondout, nextmol
logical :: found_l, notstart

call random_seq(10)
notstart = .TRUE.
extra_bond = 0 ! Need to know where ghost layer traverses happen as discard if not part of the loop
nextmol = 0

do while (notstart)
   ! Select a molecule to start on.
   rand1=int(str%nmol*random_uniform_random())+1 
 
   ! Select one of the hydrogen bonds from that water.
   rand2=int(4*random_uniform_random())+1  ! 2 in, 2 out

   if (str%connection(rand1, rand2) .GT. str%nmol)  then
      ! Goes into ghost layer - check that we can exit before proceeding!
      CALL ESCAPE_LAYER(str, str%bonds(4*rand1 - 4 + rand2), rand1, bondout, nextmol)
      if (nextmol .NE. 0) then
         extra_bond(1) = bondout
         notstart = .FALSE.
      end if
   else if (str%connection(rand1, rand2) .EQ. -1) then
      if (str%move_base) then
         CALL ESCAPE_LAYER(str, str%bonds(4*rand1 - 4 + rand2), rand1, bondout, nextmol)
         if (nextmol .NE. 0) then
            extra_bond(1) = bondout
            notstart = .FALSE.
         end if
      end if
   else if ((str%stoich(rand1) .NE. 4) .OR. (str%stoich(rand1) .NE. 0)) then
      notstart = .FALSE.
   end if
end do
      
mol_l = 0
bond_l = 0
ifail = 0

mol_l(1)=rand1                          ! Choose a molecule
bond_l(1)=4*mol_l(1) - 4 + rand2           ! The correct H bond for the molecule
if (nextmol .NE. 0) then
   ! Goes through ghost layer
   mol_l(2) = nextmol
else
   mol_l(2)=str%connection(mol_l(1),rand2) ! The right adjacent O for the molecule 
end  if
found_l=.false.


! Step one of move algorithm - we attempt to find a loop. 
if (str%bonds(bond_l(1)) .eq. 1) then
! Bond is covalent from starting "oxygen" - want to go along covalent bonds along string
! to replace like with like
do j=2,str%nmol*10  ! A long loop, dont need weird backwards indexing
   do k=1,4          ! Over all possible bonds
    if (str%connection(mol_l(j),jumble(k)).eq.mol_l(j-1)) cycle   ! Don't go backwards
    if ((str%stoich(str%connection(mol_l(j), jumble(k))) .EQ. 0) .OR. &
        (str%stoich(str%connection(mol_l(j), jumble(k))) .EQ. 4)) cycle ! Ignore (anti) ammonium

    if (str%bonds(4*mol_l(j)-4+jumble(k)).eq.1) then !  We must select a hydrogen bond donor - replace like with like
       bond_l(j)=4*mol_l(j)-4+jumble(k)
       !write(*,*) "bond number", bond_l(j), "to", str%connection(mol_l(j), jumble(k)), "from", mol_l(j)
       if (str%connection(mol_l(j), jumble(k)) .GT. str%nmol) then
          !write(*,*) "ghost to", str%connection(mol_l(j), jumble(k)), "from", mol_l(j)
          call escape_layer(str, 1, mol_l(j), bondout, nextmol)
          if (nextmol .NE. 0) then
             extra_bond(j) = bondout
             mol_l(j+1) = nextmol
             exit
          end if
       else if (str%connection(mol_l(j), jumble(k)) .EQ. -1) then
          if (str%move_base) then
             call escape_layer(str, 1, mol_l(j), bondout, nextmol)
             if (nextmol .NE. 0) then
                extra_bond(j) = bondout
                mol_l(j+1) = nextmol
                exit
             end if
          end if
       else
          mol_l(j+1)=str%connection(mol_l(j),jumble(k))
          exit
       end if
    end if
 end do


 if (mol_l(j+1) .EQ. 0) then
    ! Haven't found a link
    ifail = 2
    exit
 end if
 ! Test whether or not we have gone around a complete loop.

 do k=1,j
    if (mol_l(k).eq.mol_l(j+1)) then
       l_start=k
       l_end=j
       found_l=.true.
       exit
    end if
 end do

  if (found_l) exit

  call random_seq(10)

end do

else if(str%bonds(bond_l(1)).eq.0) then

!Hydrogen bond wrt starting oxygen
do j=2,str%nmol*10
  do k=1,4
    if (str%connection(mol_l(j),jumble(k)).eq.mol_l(j-1)) cycle   ! Don't go backwards

    if ((str%stoich(str%connection(mol_l(j), jumble(k))) .EQ. 0) .OR. &
       (str%stoich(str%connection(mol_l(j), jumble(k))) .EQ. 4)) cycle ! Ignore (anti) ammonium
    
    if (str%bonds(4*mol_l(j)-4+jumble(k)).eq.0) then !  We must select a hydrogen bond donor
       bond_l(j)=4*mol_l(j)-4+jumble(k)
       if (str%connection(mol_l(j), jumble(k)) .GT. str%nmol) then
          call escape_layer(str, 0, mol_l(j), bondout, nextmol)
          if (nextmol .NE. 0) then
             extra_bond(j) = bondout
             mol_l(j+1) = nextmol
             exit
          end if
       else if (str%connection(mol_l(j), jumble(k)) .EQ. -1) then
          if (str%move_base) then
             call escape_layer(str, 0, mol_l(j), bondout, nextmol)
             if (nextmol .NE. 0) then
                extra_bond(j) = bondout
                mol_l(j+1) = nextmol
                exit
             end if
          end if
       else
          mol_l(j+1)=str%connection(mol_l(j),jumble(k))
          exit
       end if
    end if
 end do
 if (mol_l(j+1) .EQ. 0) then
    ! Haven't found a link
    ifail = 2
    exit
 end if
 ! Test whether or not we have gone around a complete loop.
 do k=1,j
    if (mol_l(k).eq.mol_l(j+1)) then
       l_start=k
       l_end=j
       found_l=.true.
       exit
    end if
 end do

 if (found_l) exit

 call random_seq(10)

end do

else
   write(*,*) 'bond', str%bonds(bond_l(j))
   write(*,*) 'molecule', mol_l(1), rand2
   write(*,*) str%bonds(4*mol_l(1)-3:4*mol_l(1)), str%connection(mol_l(1), :)
   write(6,*)'Bond unassigned'
   stop
end if

! Move algorithm step 2 - change orientations of all bonds in loop

if (ifail .EQ. 0) then
   do j=l_start,l_end
      ! We have found a complete loop, now reverse if they are "forwards" or "backwards" - can't move as they are on a lattice
      ! Proton disorder thing
      str%bonds(bond_l(j))=MOD(str%bonds(bond_l(j))+1, 2)
      if (str%reverses(bond_l(j)) .NE. 0) then
         ! Can go into ghost layer with reverses of 0      
         str%bonds(str%reverses(bond_l(j)))=MOD(str%bonds(str%reverses(bond_l(j)))+1, 2)  ! Have to change the reverse of each bond too!
      else
         ! Went into the ghost layer (as fixed bonds can't be changed) => there is an extra bond to deal with (also won't have a reverse)
         str%bonds(extra_bond(j)) = MOD(str%bonds(extra_bond(j))+1, 2)
      end if
   end do

   coord=test(str)
   if (coord.ne.0) then
      write(*,*) bond_l
      write(6,*)'Structure generated by move routine is not a valid structure'
      stop
   end if
   
end if

end subroutine move

!=============================================================

subroutine output(str, hash)
  implicit none
  type(system), intent(inout) :: str
  character(len=4*str%nmol), intent(out) :: hash   ! For ease of regenerating structures
  integer :: i
  character(len=1) :: tmp
  
  ! Create a has string for comparison and output to file 
  hash = ''

  do i=1,4*str%nmol
     ! Including all possible connections due to partially open PBCs
     write(tmp, '(I1)')str%bonds(i)
     hash = TRIM(hash)//tmp
  end do
     
end subroutine output

!=============================================================
! Take a new structure and output the equivalent structures by changing 0/1
! Insist on sensbile values of shift outside this subroutine!

subroutine equivalent(str, hashout)
  implicit none
  type(system), intent(inout) :: str
  character(len=4*str%nmol), intent(out) :: hashout
  integer :: i
  character(len=1) :: tmp

  hashout = ''
  do i=1,str%nmol*4
     write(tmp, '(I1)')MOD(str%bonds(i)+1, 2)  !MERGE(1, 0, (str%bonds(i).EQ.0))
     hashout = TRIM(hashout)//tmp
  end do
  
end subroutine equivalent

!============================================================
subroutine escape_layer(str, bonddir, molin, bondout, nextmol)
use random, only : random_uniform_random
  implicit none
  type(system), intent(inout) :: str
  integer, intent(in) :: bonddir ! If it's a 1 or 0
  integer, intent(in) :: molin ! Where the atom enters the ghost/surface layer
  integer, dimension(:), allocatable :: tmp, untried_ghosts
  integer :: i, k, l, ii, ghosts, gs, randn
  integer, dimension(1) :: location
  integer, intent(out) :: bondout ! Extra bond from exiting layer
  integer, intent(out) :: nextmol ! Next molecule in the chain

  integer, dimension(str%plane_sz, 2) :: test_set

  ! Move enters a ficticious layer. We want it to be able to move around in this layer as much as it wants really
  ! We just need to make sure it comes out at some point other than where it went in
  ! Try to exit the layer (at random), once per atom in the layer
  
  if (COUNT(str%ghost(:, 1) .EQ. molin) .EQ. 1) then
     ! This is in top layer
     test_set(:, :) = str%ghost(:, :)
  else
     ! In surface layer
     test_set(:, :) = str%base(:, :)
  end if

  allocate(untried_ghosts(1:str%plane_sz-1))
  k = 1
  ghosts = str%plane_sz-1
  do i=1, str%plane_sz
     if (test_set(i, 1) .NE. molin) then
        ! Test to make sure atom isn't a(n) (anti) ammonium (can't be changed)
        if ((str%stoich(test_set(i, 1)) .EQ. 0) .OR. (str%stoich(test_set(i, 1)) .EQ. 4)) then
           ghosts = ghosts - 1
        else
           untried_ghosts(k) = test_set(i, 1)
           k = k+1
        end if
     end if
  end do
  ! Stop (anti) ammonium from being changed!
  if (ghosts .NE. str%plane_sz-1) then
     ! (Anti) ammonium near ghost layer - can't change!
     allocate(tmp(1:ghosts))
     tmp(1:ghosts) = untried_ghosts(1:ghosts)
     deallocate(untried_ghosts)
     allocate(untried_ghosts(1:ghosts))
     untried_ghosts=tmp
     deallocate(tmp)
  end if

  !write(*,*) "coming from", molin, "untried", untried_ghosts

  gs = ghosts
  bondout = 0 
  nextmol = 0
  do i=1, gs
     randn = int(ghosts*random_uniform_random())+1
     ! Don't need to check if it's the same atom, as already taken out
     location = MAXLOC(test_set(:,1), mask=(test_set(:,1) .EQ. untried_ghosts(randn)))
     l = location(1)
     ! Now check the bond
     if (str%bonds(4*test_set(l, 1)-4+test_set(l, 2)) .EQ. bonddir) then  ! Bond goes the wrong way, as into ghost layer, not from
        allocate(tmp(1:ghosts))
        tmp = untried_ghosts
        deallocate(untried_ghosts)
        ghosts=ghosts-1
        allocate(untried_ghosts(1:ghosts))
        ii = 1
        do k=1, ghosts+1
           if (k .NE. randn) then
              untried_ghosts(ii) = tmp(k)
              ii = ii+1
           end if
        end do
        deallocate(tmp)
     else
        bondout = 4*test_set(l, 1)-4+test_set(l, 2)
        ! Need to take account of the unmatched bonds as travelling through ghost layer
        nextmol = test_set(l, 1)
        exit 
     end if
  end do
  deallocate(untried_ghosts)
  
end subroutine escape_layer

!=============================================================

subroutine make_valid(str)
  use random, only : random_uniform_random,jumble,random_seq
  implicit none
  type(system), intent(inout) :: str
  
  ! Take in a system which is not valid (i.e. base is now constrained) and make it valid
  ! NOTE - CONNECTIONS AND REVERSES MUST BE SET CORRECTLY
  ! currently ignoring polarisation

  type(FIFO_list) :: defects
  type(node), pointer :: current, dat
  integer :: i, atom, connected
  logical :: invalid, more 

  invalid = .FALSE.
  
  do atom=1, str%nmol
     if (COUNT(str%bonds(4*atom-3:4*atom) .EQ. 1) .EQ. str%stoich(atom)) then
        cycle
     else
        invalid = .TRUE.
        ALLOCATE(current)
        current%data = atom
        CALL add_item(current, defects)
     end if
  end do

  do while (invalid) 
     call random_seq(10)
     ! This is an invalid structure, pop the first atoms which is a defect (as in genice)
     call pop(defects, dat, more)
     atom = dat%data
     NULLIFY(dat%next)
     DEALLOCATE(dat)
     invalid = more
     if (str%stoich(atom) .EQ. 4) then
        ! deal with ammonium
        str%bonds(4*atom-3:4*atom) = 1
        do i=1, 4
           if (str%reverses(4*atom-4+jumble(i)) .NE. 0 )  str%bonds(str%reverses(4*atom-4+jumble(i))) = 0  ! Change the opposite of this bond - take into account lack of y periodicity. Only called if base
           ! Now need to check if this is now a defect
           connected = str%connection(atom, jumble(i))
           ! Discount edges, where conneted is 0 or -1 (i.e. .LE. 0), and know there are no 4-4 connections
           if ((str%reverses(4*atom-4+jumble(i)) .NE. 0) .AND. &
              (COUNT(str%bonds(4*connected-3:4*connected) .EQ. 1 ) .NE. str%stoich(connected))) then
              ! Add this to the queue - if it gets sorted later that will sort itself out
              allocate(current)
              current%data = connected
              call add_item(current, defects)
              invalid = .TRUE.
           end if
        end do
     else if (str%stoich(atom) .EQ. 0) then
        ! deal with anti ammonium
        str%bonds(4*atom-3:4*atom) = 0
        do i=1, 4
           if (str%reverses(4*atom-4+jumble(i)) .NE. 0 )  str%bonds(str%reverses(4*atom-4+jumble(i))) = 1  ! Change the opposite of this bond - take into account lack of y periodicity. Only called if base
           ! Now need to check if this is now a defect
           connected = str%connection(atom, jumble(i))
           ! Discount edges, where conneted is 0 or -1 (i.e. .LE. 0), and know there are no 0-0 connections
           if ((str%reverses(4*atom-4+jumble(i)) .NE. 0) .AND. &
              (COUNT(str%bonds(4*connected-3:4*connected) .EQ. 1 ) .NE. str%stoich(connected))) then
              ! Add this to the queue - if it gets sorted later that will sort itself out
              allocate(current)
              current%data = connected
              call add_item(current, defects)
              invalid = .TRUE.
           end if
        end do
     else
        ! Not an (anti) ammonium
        if (COUNT(str%bonds(4*atom-3:4*atom) .EQ. 1) .LT. str%stoich(atom)) then
           ! We have too many H bonds (not enough covalent bonds - don't need to worry about fixed base)
           do i=1,4
              if (str%connection(atom, jumble(i)) .EQ. -1) cycle ! Can have H bonds to surface, but want to ensure the right number
              if (str%bonds(4*atom-4+jumble(i)) .EQ. 1) cycle  ! Don't change covalent bond - should be to anti-ammonium so don't explicitly check
              if (str%stoich(str%connection(atom, i)) .EQ. 4) cycle   ! Don't change ammoniums bonds - cannot be reversed
              
              str%bonds(4*atom-4+jumble(i)) = 1 ! Make an H bond a covalent bond
              if (str%reverses(4*atom-4+jumble(i)) .NE. 0 )   then
                 ! This is NOT an edge bond (must be a 0 as -1 can only have a bond value of 1)
                 str%bonds(str%reverses(4*atom-4+jumble(i))) = 0  ! Change the opposite of this bond - take into account lack of y periodicity. Only called if base
                 ! Now need to check if this is now a defect
                 connected = str%connection(atom, jumble(i))
                 if (COUNT(str%bonds(4*connected-3:4*connected) .EQ. 1 ) .NE. str%stoich(connected)) then
                    ! Add this to the queue - if it gets sorted later that will sort itself out
                    allocate(current)
                    current%data = connected
                    !write(*,*) "connected", connected
                    call add_item(current, defects)
                    invalid = .TRUE.
                 end if
              end if
              ! Need to remember to check this atom and, if needed, add back in
              if (COUNT(str%bonds(4*atom-3:4*atom) .EQ. 1 ) .NE. str%stoich(atom)) then
                 ! Add this to the queue - if it gets sorted later that will sort itself out
                 allocate(current)
                 current%data = atom
                 !write(*,*) "connected", connected
                 call add_item(current, defects)
                 invalid = .TRUE.
              end if
              exit
           end do
        else if (COUNT(str%bonds(4*atom-3:4*atom) .EQ. 1) .GT. str%stoich(atom)) then
           ! We have too many covalent bonds (lacking H bonds). Need to take into account fixed bonds
           do i=1,4
              if (str%bonds(4*atom-4+jumble(i)) .EQ. 0) cycle ! Don't change H bond - should mean don't change ammonim
              if (str%connection(atom, jumble(i)) .EQ. -1) cycle ! This is a fixed base bond, don't change (ensure right number)
              if (str%stoich(str%connection(atom, i)) .EQ. 0) cycle ! Ignore anti-ammonium
              str%bonds(4*atom-4+jumble(i)) = 0 ! Make covalent bond an H bond
              if (str%reverses(4*atom-4+jumble(i)) .NE. 0 ) then
                 ! This is NOT an edge bond (must be a 0 as -1 already ignored)
                 str%bonds(str%reverses(4*atom-4+jumble(i))) = 1  ! Change the opposite of this bond - take into account lack of y periodicity. Only called if base
                 ! Now need to check if this is now a defect
                 connected = str%connection(atom, jumble(i))
                 if (COUNT(str%bonds(4*connected-3:4*connected) .EQ. 1) .NE. str%stoich(connected)) then
                    ! Add this to the queue - if it gets sorted later that will sort itself out
                    allocate(current)
                    current%data = connected
                    call add_item(current, defects)
                    invalid = .TRUE.
                 end if
              end if
              ! Need to remember to check this atom and, if needed, add back in
              if (COUNT(str%bonds(4*atom-3:4*atom) .EQ. 1 ) .NE. str%stoich(atom)) then
                 ! Add this to the queue - if it gets sorted later that will sort itself out
                 allocate(current)
                 current%data = atom
                 !write(*,*) "connected", connected
                 call add_item(current, defects)
                 invalid = .TRUE.
              end if
              exit
           end do
        end if

     end if
  end do
!  write(*,*) str%bonds
!  write(*,*) str%nitrogen
!  write(*,*) str%base
end subroutine make_valid

!=============================================================

subroutine add_item(item, list)
  implicit none
  type(node), pointer :: item
  type(FIFO_list) :: list
  
  IF(ASSOCIATED(list%tail)) THEN
     !List tail is not null, i.e. there is at least one item

     ! Next link of old tail is item
     list%tail%next => item
     ! item becomes new tail
     list%tail => item
     ! Just in case item%next was not null
     NULLIFY(item%next)

    ELSE
      ! New item is both head and tail!
      list%tail => item
      list%head => item
      ! Just in case item pointers are not null
      NULLIFY(item%next)
    END IF
end subroutine add_item

!=============================================================

subroutine pop(list, dat, more)
  implicit none
!  integer, intent(out) :: dat
  logical, intent(out) :: more
  ! Pop the head of the list
  type(FIFO_list) :: list
  type(node), pointer :: dat

  more = .TRUE.
  dat => list%head
 
  IF (ASSOCIATED(list%head, TARGET=list%tail)) then
     ! Only item in the list
     NULLIFY(list%head%next)
     NULLIFY(list%head, list%tail)
     more = .FALSE.
  else
     list%head => list%head%next
  end IF
  
end subroutine pop

! ==========================================================================================
    
end module hbond_module 

