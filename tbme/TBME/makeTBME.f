CCCCCCCCCCCCCCCCCALVINCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program makeTBME

C
C  creates TBME (+ spe's) for input into Oxbash, Glasgow codes 
C

C -- creates or reads in list of single-particle orbits 
C -- choose from a variety of allowed interactions
C -- output as either Oxbash or Glasgow
C
C  Calvin Johnson
C
C  version 1.0  Feb 1999, LSU
C  version 1.1  March 2007, SDSU : updated:
C                   uses cleb, threej,sixj routines from libra.f
C                   improved input
C                   uses memory allocation & modules
!  version 2.0 Dec 2014: add in "auto" sps option
C
C  SUBROUTINE CALLS:
C	MAKESPLIST: compiles list of single-particle states 
C	GetOxSPS: reads oxbash .sps file of orbits
C	TBMElist: from s.p. orbits, creates all independent tbme list
C	TBME: actually computes the tbme's	
C

      use sporbit
      use tbmes
      implicit none

C---- SINGLE-PARTICLE LIST------

C----TBME LIST
      integer nme1			! # of TBME
C..................PROGRAM CONTROL.................

     				
C..................FILE CONTROL...................

      character*1 ychar
      character filename*70             ! name of file
      integer ilast

C---------- DUMMY COUNTERS ETC -------

      integer I,K
      
C BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN 
      
      write(6,*)' '
      write(6,*)' * * * * * * * * * * * * * * * * * * *  '
      write(6,*)' *                                   * '
      write(6,*)' *  Welcome !                        * '
      write(6,*)' *                                   * '
      write(6,*)' *  Version 2.0 (Dec 2014)           * '
      write(6,*)' *                                   * '
      write(6,*)' * * * * * * * * * * * * * * * * * * *  '
      write(6,*)' '

C.... 
      call get_sp_info
      
C---- CREATE INTERNAL LIST OF ALL INDEPENDENT MATRIX ELEMENTS & Q#s

      call TBMElist
      
C      do i = 1,nme
C      	write(6,*)jtot(i),ttot(i),(map(i,k),k=1,4)
C      enddo

C---- CREATE TBME's + ASSOCIATED SPE's----------------------------

      call TBMEcalc


C------------------ WRITE OUTPUT  ------------

      call writeout
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine get_sp_info
C
C  reads in s.p. space information;
C  fills in q#s into orbqn and spsqn
C
C  reads in either .spo or.sps files
C  it first looks for .spo file; if that fails, 
C  automatically looks for .sps file of same name
C
C----------------------------------------------

      use sporbit
!      use spstate
      implicit none

C------ FILE CONTROL ---------------------------

      character*15 filename
      character*1 achar
       
      character*70 title
      integer ilast


C------------ TEMP -------------

      real :: xn,xl,xj,xw    ! orbital q#s
C-----------------------------------------------
C      dummy counters
C--------------------------------------------------
      integer i,j,k,m
      logical success
C---------------BEGIN -------------------------

C-------------- OPEN A FILE ----------------------------------
      success = .false.

      do while(.not.success)
          print*,' Enter file with s.p. orbit information (.spo/.sps)'
          read(5,'(a)')filename
          ilast = index(filename,' ')-1
          if ( filename(1:ilast)=='auto') then
            call autofillsps
            return
          end if
C........... ATTEMPT TO OPEN .spo FILE.................
          open(unit=1,file=filename(1:ilast)//'.spo',status='old', 
     &     err=101)
          success = .true.
          cycle
101       continue
C...........ATTEMPT TO OPEN .sps FILE..................
          open(unit=1,file=filename(1:ilast)//'.sps',status='old',
     &      err=102)
          success = .true.
          cycle
102       continue
          print*,filename(1:ilast),'.spo/.sps file does not exist '

      enddo
C-------------- READ PAST TITLE CARDS ---------------------------

      success = .false.
      do while(.not.success)
        read(1,'(a)')achar
        backspace(1)
        if(achar /= '#' )then
           success = .true.
        else
           read(1,'(a)')title
           write(6,*)title
        endif
      enddo

C-------------- READ PAST POSSIBLE LABEL OF ISO/PN

      read(1,'(a)')achar
      if(achar == 'p' .or. achar=='P')then
          print*,' .sps file in pn formalism, cannot handle '
          stop
      elseif(achar /= 'i' .and. achar/='I')then
          backspace(1)
      endif

C............ READ # OF ORBITS----------

      read(1,*)numorb

C----------------ALLOCATE MEMORY ------------

      allocate(orbqn(numorb))

C---------------READ IN---------------------

      do i = 1,numorb
        read(1,*,end=2001)xn,xl,xj,xw
        orbqn(i)%nr = nint(xn)
        orbqn(i)%l = nint(xl)
        orbqn(i)%j = nint(2*xj)
        orbqn(i)%w = nint(xw)
      enddo

      close(unit=1)

!--------------- CHECK IF "SPINLESS"-------------
!
! sometimes fermions do not carry intrinsic spin; this is used
! for atomic cases or cold atoms
!
      spinless = .false.
      if( 2*orbqn(1)%l == orbqn(1)%j) spinless = .true.

      if(spinless)then
         print*,' WARNING -- for many cases antisymmeterization',
     &' not yet corrected '
      endif
      if( spinless)then      ! double check
        if( (orbqn(1)%j/2)*2 /= orbqn(1)%j)then
           print*,' l and j do not make sense ',spinless
           print*,orbqn(1)%l, orbqn(1)%j, (orbqn(1)%j/2)*2
           print*,orbqn(1)%nr
           stop
        endif
      else
        if( (orbqn(1)%j/2)*2 == orbqn(1)%j)then
           print*,' l and j do not make sense ',spinless
           print*,orbqn(1)%j, orbqn(1)%l
           stop
        endif

      endif

      return
C------------- ERROR TRAP FOR END OF FILE -----------
2001  continue
      print*,' sudden end of file in ',filename(1:ilast)
      stop

      end subroutine get_sp_info
!=====================================================================

      subroutine autofillsps

      use sporbit
      implicit none
      integer :: ierr
      integer :: nprinc
      integer :: i,j,n,l,lparity
      integer :: it

      print*,' Enter maximum principle quantum number N '
      print*,' (starting with 0s = 0, 0p = 1, 1s0d = 2, etc. ) '
      read*,nprinc


      numorb = (nprinc+1)*(nprinc+2)/2
      allocate( orbqn(numorb ))
      n = 0
      do i = 0,nprinc
         lparity = mod(i,2)

         do j = 0,i
            n = n + 1
               orbqn(n)%j  = 2*j+1
               l = j + 1
               if(mod(l,2)/=lparity)l = l-1
               orbqn(n)%l = l
               orbqn(n)%w  = i
               orbqn(n)%nr = (i - l)/2

         end do   ! j 
      end do  ! i

      return
      end subroutine autofillsps
!==================================================================

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine TBMElist     
C  creates list of all possible TBME matrix elements 
C

      use sporbit
      use tbmes
      implicit none

      
      type pairmaster
         integer j
         integer t
         integer ia
         integer ib
      end type pairmaster
      type (pairmaster),allocatable :: pair(:)

      integer npairs,jmin,jmax,kk

      integer i,k
      integer li,lf

      logical compare			! function for checking



      integer ia,ib
      integer ja,jb
 
C------------------

      allocate(spe(numorb))
      spe = 0.
C--------------------- sort through and count all possible pairs 

      npairs = 0
      do ia = 1,numorb
        ja = orbqn(ia)%j
	do ib = ia,numorb
            jb = orbqn(ib)%j
	    jmin = abs(ja-jb)
	    jmax = ja+jb
            do k = jmin,jmax,2
		if(ia.eq.ib)then
		   npairs = npairs +1
		else
		    npairs = npairs + 2
		endif
	    enddo
	enddo
      enddo

      write(6,*)npairs, ' two-particle pairs possible '
C-------------------------- now create those pairs ---------

      allocate(pair(npairs))
      npairs = 0
      do ia = 1,numorb
        ja = orbqn(ia)%j
	do ib = ia,numorb
            jb = orbqn(ib)%j
	    jmin = abs(ja-jb)
	    jmax = ja+jb
            do k = jmin,jmax,2
		if(ia.eq.ib)then
		   npairs = npairs +1
		   pair(npairs)%ia=ia
		   pair(npairs)%ib=ib
		   kk = k/2
		   pair(npairs)%j = kk
                   if(spinless)then
		   if((kk/2)*2.eq.kk)then  
			pair(npairs)%t = 0
		   else
			pair(npairs)%t = 1
		   endif

                   else
		   if((kk/2)*2.eq.kk)then  
			pair(npairs)%t = 1
		   else
			pair(npairs)%t = 0
		   endif
                   endif
		else
		   npairs = npairs + 1
		   pair(npairs)%ia=ia
		   pair(npairs)%ib=ib
                   pair(npairs)%j = k/2
		   pair(npairs)%t = 0
   	           npairs = npairs + 1
		   pair(npairs)%ia=ia
		   pair(npairs)%ib=ib
                   pair(npairs)%j = k/2
		   pair(npairs)%t = 1
		endif
	    enddo
	enddo
      enddo

C--------------------- now count all possible matrix elements
C----------- include cuts on l
      nme = 0
      do i = 1,npairs
        li = orbqn( pair(i)%ia )%l + orbqn( pair(i)%ib)%l
	do k = i,npairs
            lf = orbqn( pair(k)%ia )%l + orbqn( pair(k)%ib)%l

	    if(pair(i)%j == pair(k)%j .and. 
     &         pair(i)%t == pair(k)%t .and. 
     &         mod(li+lf,2) == 0)then
		nme= nme+1	
	    endif
	enddo
      enddo		

      write(6,*)nme,' two-body matrix elements '

C--------------------- now list all possible matrix elements
      allocate(vtbme(nme))
      nme = 0
      do i = 1,npairs
        li = orbqn( pair(i)%ia )%l + orbqn( pair(i)%ib)%l
	do k = i,npairs
            lf = orbqn( pair(k)%ia )%l + orbqn( pair(k)%ib)%l
	    if(pair(i)%j == pair(k)%j .and. 
     &         pair(i)%t == pair(k)%t .and. 
     &         mod(li+lf,2) == 0)then
		nme= nme+1	
                vtbme(nme)%j = pair(i)%j
                vtbme(nme)%t = pair(i)%t
                vtbme(nme)%n(1) = pair(i)%ia
                vtbme(nme)%n(2) = pair(i)%ib
                vtbme(nme)%n(3) = pair(k)%ia
                vtbme(nme)%n(4) = pair(k)%ib
                vtbme(nme)%v = 0.0
	    endif
	enddo
      enddo		
      print*,' check nme = ',nme
      return     
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine writeout

      use sporbit
      use tbmes
      implicit none
C..................FILE CONTROL...................

      character*1 ychar
      character filename*15             ! name of file
      integer ilast
      integer i,k
      integer nme1
      logical success

      success = .false.

      do while(.not.success)

      write(6,*)' Name for interaction file? (.int) '
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

C      check number of nonzero matrix elements


C      nme1=nme
      open(unit=13,file=filename(1:ilast)//'.int',status='new',
     & err=1007)
      success = .true.
      cycle

1007  continue
      write(6,*)' That files exists; overwrite (y/n)? '
      read(5,'(a)')ychar
      if(ychar == 'y' .or. ychar == 'Y')then
      open(unit=13,file=filename(1:ilast)//'.int',status='unknown')
      success = .true.
      endif

      enddo

      nme1 = 0
      do i =1,nme
	if(vtbme(i)%v.ne.0.0)nme1=nme1+1
      enddo	


      write(13,1002) nme1, (spe(i),i=1,min(10,numorb))
      if(numorb > 10)then
        write(13,1004)(spe(i),i=11,numorb)
      endif

      do i = 1,nme
      	if(vtbme(i)%v.ne.0.0)then
  	write(13,1003)(vtbme(i)%n(k),k=1,4),vtbme(i)%j,vtbme(i)%t,
     &                                      vtbme(i)%v
        endif
      enddo
      close(unit=13)
 1002 format(i9,10f10.5)
 1003 format(4i4,5x,2i4,4x,f14.7)    	
 1004 format(10f10.5)


      
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
