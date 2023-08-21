
    implicit none
    integer Npts
    real dr

    call get_sp_info

    print*,' enter # pts, dr '
    read*,npts, dr
    call test_contact_integrals(Npts,dr)

    end


      subroutine get_sp_info
!C
!C  reads in s.p. space information;
!C  fills in q#s into orbqn and spsqn
!C
!C  reads in either .spo or.sps files
!C  it first looks for .spo file; if that fails, 
!C  automatically looks for .sps file of same name
!
!----------------------------------------------

      use sporbit
      use spstate
      implicit none

!------ FILE CONTROL ---------------------------

      character*15 filename
      character*1 achar
       
      character*70 title
      integer ilast


!------------ TEMP -------------

      real :: xn,xl,xj,xw    ! orbital q#s
!-----------------------------------------------
!      dummy counters
!--------------------------------------------------
      integer i,j,k,m
      logical success
!---------------BEGIN -------------------------

!-------------- OPEN A FILE ----------------------------------
      success = .false.

      do while(.not.success)
          print*,' Enter file with s.p. orbit information (.spo/.sps)'
          read(5,'(a)')filename
          ilast = index(filename,' ')-1
!........... ATTEMPT TO OPEN .spo FILE.................
          open(unit=1,file=filename(1:ilast)//'.spo',status='old',  err=101)
          success = .true.
          cycle
101       continue
!...........ATTEMPT TO OPEN .sps FILE..................
          open(unit=1,file=filename(1:ilast)//'.sps',status='old', err=102)
          success = .true.
          cycle
102       continue
          print*,filename(1:ilast),'.spo/.sps file does not exist '

      enddo
!-------------- READ PAST TITLE CARDS ---------------------------

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

!-------------- READ PAST POSSIBLE LABEL OF ISO/PN

      read(1,'(a)')achar
      if(achar == 'p' .or. achar=='P')then
          print*,' .sps file in pn formalism, cannot handle '
          stop
      elseif(achar /= 'i' .and. achar/='I')then
          backspace(1)
      endif

!............ READ # OF ORBITS----------

      read(1,*)numorb

!----------------ALLOCATE MEMORY ------------

      allocate(orbqn(numorb))

!---------------READ IN---------------------

      do i = 1,numorb
        read(1,*,end=2001)xn,xl,xj,xw
        orbqn(i)%nr = int(xn)
        orbqn(i)%l = int(xl)
        orbqn(i)%j = int(2*xj)
        orbqn(i)%w = int(xw)
      enddo

      close(unit=1)

!------------ SET UP S.P. STATE INFO -------------

      nsps = 0
      do i = 1,numorb
          nsps = nsps+orbqn(i)%j+1
      enddo

      allocate(spsqn(nsps))
      k = 0
      do i = 1,numorb
        j = orbqn(i)%j
        do m = -j,j,2
           k = k +1
           spsqn(k)%nr = orbqn(i)%nr
           spsqn(k)%l = orbqn(i)%l
           spsqn(k)%par = (-1)**(orbqn(i)%l)
           spsqn(k)%j = orbqn(i)%j
           spsqn(k)%w = orbqn(i)%w
           spsqn(k)%m = m
           spsqn(k)%orb = i
        enddo
      enddo

      return
!------------- ERROR TRAP FOR END OF FILE -----------
2001  continue
      print*,' sudden end of file in ',filename(1:ilast)
      stop

      end subroutine get_sp_info