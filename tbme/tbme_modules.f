CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  modules for code BASHO
C  
C  CWJ @ SDSU  started 3/2007
C  based upon REDSTICK by W.E.Ormand 
C  & ELDORADO by CWJ
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module sporbit
C
C  single-particle ORBIT information
C
      implicit none

      integer numorb          ! # of s.p. orbits

C------------ CREATE A DEFINED TYPE----------------
      type orb
        integer :: nr       ! radial quantum number
        integer :: j        ! 2 x j
        integer :: l        ! L
        integer :: w        ! excitation 

      end type orb

      type (orb),allocatable :: orbqn(:)
      logical spinless  ! flag when fermions do not have intrinsic spin

      end module sporbit

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module spstate
C
C  single-particle STATE information
C
      implicit none

      integer nsps          ! # of s.p. states

C------------ CREATE A DEFINED TYPE----------------
      type spst
        integer :: nr       ! radial quantum number
        integer :: j        ! 2 x j
        integer :: m        ! 2 x jz
        integer :: l        ! L
        integer :: w        ! excitation 
        integer :: par      ! parity
        integer :: orb      ! orbit label
      end type spst

      type (spst),allocatable :: spsqn(:)

      end module spstate

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module tbmes

      implicit none
      integer nme   ! # of tbmes

      real, allocatable :: spe(:)  ! single-particle energies

      type tbme 
         integer j
         integer t
         integer n(4)  
         real v
      end type tbme

      type (tbme),allocatable :: vtbme(:)
		  
		  real,allocatable :: obme(:,:)
      end module tbmes
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc

      module Vrelative

      integer lmax    ! max l, starting at 0

      integer nmax    ! max n, starting at 0

      type cousin
         real, allocatable :: v(:,:)
      end type cousin
      type (cousin), allocatable,target :: vrel(:,:)  ! l,s
      integer moshphase

      end module Vrelative
!=======================================================