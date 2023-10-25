! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                          R A N D O M                                        !
!=============================================================================!
!                                                                             !
! $Id: random.f90,v 1.2 2007/11/06 16:30:58 ccseac Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! This is the serial version of the random number generator written by Matt   !
! Probert and used in CASTEP. This module is basically a copy of the relevant !
! routines from the CASTEP Algor module.                                      !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log$
!
!
!=============================================================================!

Module Random

!  Use Constants , Only : dp
  Implicit None                                 !Impose strong typing

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)  
  
  Private                                       !Everything is private ...

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  Public :: random_set_random_seed               !... unless exposed here.
  Public :: random_uniform_random
  Public :: random_seq
  Public :: random_normal_random
  Public :: random_test_uniform

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  Public :: jumble 
  integer,dimension(1:4) :: jumble

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  !For the uniform random number generator, we presume a 32-bit integer in order
  !to do bit-shuffling, so define and store here.
  Integer, Parameter        :: bit_32=Selected_int_kind(9)
  Integer(kind=bit_32),Save :: ix,iy                     !intermediate randoms

  Logical,Save              :: random_initialised=.False.
   
Contains

  Subroutine Random_set_random_seed(seed,dumrank,dumsize)
    !=========================================================================!
    ! Set module internal random number generator seed.                       !
    ! If seed=0, then use system clock to get a unique seed every call.       !
    ! If seed<>0 then use as deterministic seed.                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   seed intent(in)=> =0 random seed, else deterministic seed.            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   iseed, ix, iy and random_initialised are set here.                    !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   comms for different iseed on different nodes iff seed=0               !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   iseed as starting 32-bit integer seed for all generators              !
    !   ix as seed 32-bit integer in Marsaglia generator                      !
    !   iy as seed 32-bit integer in Park-Miller generator                    !
    !-------------------------------------------------------------------------!
    ! Architecture dependencies:                                              !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   comms must be initialised first.                                      !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 01/07/2000                               !
    !=========================================================================!
    Implicit None
    Integer, Intent(in) :: seed                !may or may not be 32 bit...

    !32 bit integer version of seed
    Integer(kind=bit_32) :: iseed              !random number seed
    integer,intent(in),optional   :: dumrank 
    integer,intent(in),optional   :: dumsize   

    !f90 intrinsic time call
    Character(len=10) :: system_time           !length is crucial ...
    Real (kind=dp)    :: rtime

    !comms variables
    integer :: myrank=0,size=1

    if (present(dumrank)) myrank = dumrank
    if (present(dumsize)) size   = dumsize


    !We only want to set the seed if it has not already been set
    !and the flag is set to .false. by the compiler
    If (.Not.random_initialised) Then

       If (seed == 0 ) Then
          Call Date_and_time(time=system_time)  !character string hhmmss.xxx
          Read (system_time,*) rtime            !convert to real
          rtime = rtime * 1000.0_dp             !0<rtime<235959999.0 which fits within huge(1)
       Else
          rtime = Real(Abs(seed),kind=dp)          !convert seed to real
       End If
       
       !and then convert to bit_32 size integer
       iseed = Int(rtime,kind=bit_32)              !must fit within huge(1)
       !and make sure it is different on each node
       iseed = iseed + Int(MyRank*1000,kind=bit_32) -1_bit_32

       ix=Ieor(777755555_bit_32,iseed)                   !Marsaglia generator
       iy=Ior(Ieor(888889999_bit_32,iseed),1_bit_32)     !Parks-Miller generator

       random_initialised = .true.

       ! Initialize jumble array - used in move
       jumble=(/1,2,3,4/)

    End If

    Return

  End Subroutine random_set_random_seed

  Function random_uniform_random()
    !=========================================================================!
    ! Return a single random deviate ~ uniform [0,1].                         !
    ! Based on Park-Miller "minimal standard" generator with Schrage's method !
    !  to do 32-bit multiplication without requiring higher precision, plus   !
    !  additional Marsaglia shift to suppress any weaknesses & correlations.  !
    ! Using two independent methods greatly increases the period of the       !
    !  generator, s.t. resulting period ~2*10^18                              !
    ! NB Routine is only set to work with 32 bit integers!                    !
    !-------------------------------------------------------------------------!
    ! References:                                                             !
    !   S.K. Park and K.W. Miller, Commun. ACM, 31, p1192-1201 (1988)         !
    !   L. Schrage, ACM Trans. Math. Soft., 5, p132-138 (1979)                !
    !   G. Marsaglia, Linear Algebra and its Applications, 67, p147-156 (1985)!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Return value:                                                           !
    !   algor_uniform_random => required random deviate                       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ix as next 32-bit integer in Marsaglia generator (updated)            !
    !   iy as next 32-bit integer in Park-Miller generator (updated)          !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 01/07/2000                               !
    !=========================================================================!
    Implicit None
    Real(kind=dp)                 :: Random_uniform_random

    !NB We use 3 logical XOR operations to do Marsaglia shift
    !=> hard-wire 3 shift values (not all triplets any good)
    !=> entire routine preset to only work with 32 bit integers.

    !local variables ...
    Integer(kind=bit_32)            :: iy_tmp       !working value to force integer division
    Integer(kind=bit_32), Parameter :: iy_max=2147483647 !2^31-1
    Real(kind=dp), Parameter        :: inv_iy_max=1.0_dp/2147483647.0_dp

    !Catch uninitialised random number generator and set to random seed
    If (.Not.random_initialised) Call random_set_random_seed(0)

    !do Marsaglia shift sequence, period 2^32-1, to get new ix
    ix=Ieor(ix,Ishft(ix, 13_bit_32))
    ix=Ieor(ix,Ishft(ix,-17_bit_32))
    ix=Ieor(ix,Ishft(ix,  5_bit_32))

    !Park-Miller sequence, period iy_max-1, to get new iy
    iy_tmp=iy/127773_bit_32                         !NB integer division
    iy=16807_bit_32*(iy-iy_tmp*127773_bit_32)-iy_tmp*2836_bit_32  !next value of iy
    If (iy < 0_bit_32) iy=iy+iy_max                 !integer modulo iy_max

    !Combine ix and iy to get new random number, rescale onto range [0,1]
    !with masking to ensure non-zero
    random_uniform_random=inv_iy_max*Ior(Iand(iy_max,Ieor(ix,iy)),1_bit_32)
    
    Return
  End Function random_uniform_random

! Subroutine randomizes a sequence of four numbers
  subroutine random_seq(num)
    implicit none
    integer,intent(in) :: num
    integer :: rand1,rand2,j,tmp

    do j=1,num
       rand1=min(int(4*random_uniform_random())+1,4)
       rand2=min(int(4*random_uniform_random())+1,4)
       tmp=jumble(rand1)
       jumble(rand1)=jumble(rand2)
       jumble(rand2)=tmp
    end do

    return

  end subroutine random_seq

  function random_normal_random()
    !=========================================================================!
    ! Return a single random deviate ~ Normal(0,1)                            !
    ! Based on the Box-Muller method, using the uniform_random generator.     !
    ! Actually generates a pair of uniform deviates, checks that they lie     !
    !      within the unit circle, s.t. can map onto cos/sin without being    !
    !      explicit, and then transforms into a pair of Gaussian deviates.    !
    ! NB As second deviate is stored in routine we only need to do the        !
    !    calculation every other call!                                        !
    !-------------------------------------------------------------------------!
    ! References:                                                             !
    !   G.E.P. Box & M.E. Muller, Ann. math. stat. 29, p610-611 (1958)        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Return value:                                                           !
    !   algor_normal_random => required random deviate                        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   v1,v2  = cos/sin theta if radius=|v1^2+v2^2|<1                        !
    !   normal = sqrt(-2.0_dp*log(radius)/radius)*cos(theta) (or sin(theta))  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 01/07/2000                               !
    !=========================================================================!

    implicit none
    real(kind=dp) :: random_normal_random

    !local variables ...
    logical, save       :: done2=.false.
    real(kind=dp), save :: next_random
    real(kind=dp)       :: fac, v1, v2, radius
    integer             :: seed

    seed=0
    if (.not.random_initialised) call random_set_random_seed(seed)

    if (done2) then
       random_normal_random=next_random
       done2=.false.
    else
       do
          v1 = 2.0_dp*random_uniform_random() - 1.0_dp
          v2 = 2.0_dp*random_uniform_random() - 1.0_dp
          radius = v1**2 + v2**2
          if (radius<1.0_dp) exit
       end do
       fac = sqrt(-2.0_dp*log(radius)/radius)
       random_normal_random = v1*fac
       next_random          = v2*fac
       done2=.true.
    end if

    return

  end function random_normal_random


  subroutine random_test_uniform(log)
    !=========================================================================!
    ! Tests the distribution of random numbers generated by the above         !
    ! random_uniform_random routine.                                          !
    !-------------------------------------------------------------------------!
    ! D.Quigley - April 2011                                                  !
    !=========================================================================!
    implicit none
    integer,parameter :: ntest = 1000000
    integer,parameter :: nbins = 100
    
    integer,intent(in) :: log  ! unit number to report results to
    
    integer,dimension(nbins) :: histogram
    integer :: itest,k
    real(kind=dp) :: xi

    write(log,'("#                                                              #")')   
    write(log,'("# Testing random number generator with ",I12,"   samples  #")')ntest
    write(log,'("# -------------------------------------------------------------#")')
    write(log,'("#                                                              #")') 

    histogram = 0

    do itest = 1,ntest
       
       xi = random_uniform_random()
       k  = int(xi*real(nbins,kind=dp)) + 1
       histogram(k) = histogram(k) + 1

    end do

    write(log,'("# Ntrials / Nbins            = ",I12,10x" samples  #")')int(real(ntest,kind=dp)/real(nbins,kind=dp))
    write(log,'("# Min hits per histogram bin = ",I12,10x" samples  #")')minval(histogram,1)
    write(log,'("# Max hits per histogram bin = ",I12,10x" samples  #")')maxval(histogram,1)
    write(log,'("#                                                              #")') 

    return

  end subroutine random_test_uniform

End Module random
