CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C

      subroutine transformRelToLab
C
C transforms from relative (s-wave only) to Lab
C

      implicit none
C-------------- BEGIN INTERFACES ----------

      interface

      subroutine countcreateorbpairs(countflag,npairs,pairs,jmin,jmax)
         implicit none
         logical :: countflag
         integer :: npairs
         integer, pointer :: pairs(:,:),jmin(:),jmax(:)

      end subroutine countcreateorbpairs
      
      subroutine countcreatetbmes(countflag,npairs,pairs,jmin,jmax)
         implicit none
         logical :: countflag
         integer :: npairs
         integer, pointer :: pairs(:,:),jmin(:),jmax(:)

      end subroutine countcreatetbmes
      end interface

C---------------END INTERFACES


      integer :: numorbpairs
      integer, pointer :: orbpairs(:,:), pairjmin(:),pairjmax(:)

      call countcreateorbpairs(.true.,numorbpairs,orbpairs,
     &  pairjmin,pairjmax)
      call countcreateorbpairs(.false.,numorbpairs,orbpairs,
     &  pairjmin,pairjmax)

      call countcreatetbmes(.true.,numorbpairs,orbpairs,pairjmin,
     & pairjmax)

      call countcreatetbmes(.false.,numorbpairs,orbpairs,pairjmin,
     & pairjmax)
      return

      end subroutine transformRelToLab

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine countcreateorbpairs(countflag,npairs,pairs,jmin,jmax)

      use sporbit

      implicit none

      logical :: countflag
      integer :: npairs
      integer, pointer :: pairs(:,:),jmin(:),jmax(:)

C---------- intermediate

      integer a,b
      integer j, minj, maxj

      npairs = 0

      do a = 1,numorb
        do b = a, numorb
             npairs = npairs + 1
             minj = abs(orbqn(a)%j - orbqn(b)%j)/2
             maxj = (orbqn(a)%j + orbqn(b)%j)/2
             if(.not.countflag)then
                pairs(npairs,1) = a
                pairs(npairs,2) = b
                jmin(npairs) = minj
                jmax(npairs) = maxj
             endif
        enddo

      enddo
      if(countflag)then
         print*,' # pairs = ',npairs
         allocate(pairs(npairs,2),jmin(npairs),jmax(npairs))
      endif

      return
      end subroutine countcreateorbpairs

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine countcreatetbmes(countflag,npairs,pairs,jmin,jmax)

      use sporbit
      use tbmes
      use hamrel
      use parameters
      implicit none

      logical :: countflag
      integer :: npairs
      integer, pointer :: pairs(:,:),jmin(:),jmax(:)

C---------- intermediate -------------

      integer i1,i2
      integer L1,L2
      integer L
      integer ntmp
      logical ident
      integer a,b,c,d
      integer ea,eb,ec,ed
      integer la,lb,lc,ld
      integer n1,n2
      integer  E1, E2,Ecm  ! the principle quantum numbers, really
      real fact1,fact2
      real vtmp

      real*8 tmb, tmb1,tmb2

      ntmp = 0
      do i1 = 1,npairs

        a = pairs(i1,1)
        b = pairs(i1,2)
        ea = orbqn(a)%w
        la = orbqn(a)%l
        eb = orbqn(b)%w
        lb = orbqn(b)%l
        if(a /= b)then
            fact1 = sqrt(2.0)
        else
            fact1 = 1.0
        endif
        do i2 = i1,npairs
C           print*,i1,i2, max(jmin(i1),jmin(i2)),min(jmax(i1),jmax(i2))
           if( max(jmin(i1),jmin(i2)) > min(jmax(i1),jmax(i2))) cycle  ! can't get ang mom right 
           c = pairs(i2,1)
           d = pairs(i2,2)
           ec = orbqn(c)%w
           lc = orbqn(c)%l
           ed = orbqn(d)%w
           ld = orbqn(d)%l

           if(c /= d)then
              fact2 = sqrt(2.0)
           else
              fact2 = 1.0
           endif
           E2 = orbqn(c)%w + orbqn(d)%w
           if( a == b .or. c == d)then
               ident = .true.
           else
               ident = .false.
           endif
           do L = max(jmin(i1),jmin(i2)), min(jmax(i1),jmax(i2))
             if(ident .and. (L/2)*2 /=L)cycle  ! make sure even for identical particles

             vtmp = 0.0

             if(countflag)then
              ntmp = ntmp +1 
             endif
             if(.not.countflag)then
             do E1 = 0, ea + eb -L,2
                Ecm = ea + eb - E1
                E2 =  ec + ed -Ecm
                if(E2 < 0) cycle
                n1 = e1/2 +1
                n2 = e2/2 +1
                if(n1 > q+1 .or. n2 > q+1)cycle

                tmb1 = tmb(Ecm, L, E1, 0, ea,la,eb,lb,L,1.d0)
                tmb2 = tmb(Ecm, L, E2, 0, ec,lc,ed,ld,L,1.d0)
                vtmp = vtmp +  tmb1*tmb2*fact1*fact2*vnew(n1,n2)

             enddo

             if(abs(vtmp) > 0.00001)then
                ntmp = ntmp + 1

                   vtbme(ntmp)%v = vtmp
                   vtbme(ntmp)%j = l
                   vtbme(ntmp)%orb(1) = a
                   vtbme(ntmp)%orb(2) = b
                   vtbme(ntmp)%orb(3) = c
                   vtbme(ntmp)%orb(4) = d
                endif

             endif         
           enddo ! loop over J
        enddo ! loop over i2
      enddo ! loop over i1

      if(countflag)then
        print*,' there are ',ntmp,' matrix elements '
       
        allocate(vtbme(ntmp))
      else
        ntbme = ntmp
      endif

      return
      end subroutine countcreatetbmes