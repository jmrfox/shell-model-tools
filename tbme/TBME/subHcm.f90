 subroutine Hcm_atom

 use sporbit
 use tbmes
 implicit none

 real b
 integer A
 logical speflag

 print*,' '
 print*,' Enter total number of particles '
 read*,A

 print*,' Enter osc. parameter b of c.m. relative to b for trap '
 read*,b

 if(b <= 0.0)then
   print*,' Problem with b ',b
   stop
 endif

 call R2atom(b**2/A)
 call P2atom(1./A/b**2)

 speflag = .false. ! .true.

 if(b /= 1.0)speflag = .false.
 if(speflag)then
      call speCMatom_diag(A)
 else
	call homaster(A,b)
 endif

 return

 end subroutine Hcm_atom

!==================================================
!
! alternative version, starting from relative frame
!
 subroutine Hcm_atom2

 use sporbit
 use Vrelative
 implicit none

 real b
 integer A

 print*,' '
 print*,' Enter total number of particles '
 read*,A

 print*,' Enter osc. parameter b of c.m. relative to b for trap '
 read*,b

 if(b <= 0.0)then
   print*,' Problem with b ',b
   stop
 endif
 print*,' Which Moshinsky phase do you want? (0= GB, 1 = std) '
 read*,moshphase
 call homaster2(A,b)
 call setVrelcm(b,A)
  call binom2
  call trinom
 call Vlabfromrel

 return

 end subroutine Hcm_atom2

!==================================================



!==================================================
!
! alternative version, starting from relative frame
!
 subroutine Rintr_atom

 use sporbit
 use Vrelative
 implicit none

 real b
 integer A


 print*,' Enter osc. parameter b of c.m. relative to b for trap '
 read*,b

 if(b <= 0.0)then
   print*,' Problem with b ',b
   stop
 endif
 print*,' Which Moshinsky phase do you want? (0= GB, 1 = std) '
 read*,moshphase
 call setVrelR2(b,A)
  call binom2
  call trinom
 call Vlabfromrel

 return

 end subroutine Rintr_atom
!=================================================


!=================================================
 subroutine Rrel_atom_lab

 use sporbit
 use tbmes
 implicit none

 real, allocatable::  rawme(:,:)
 integer i,j
 integer ni,li,nj,lj
 integer a,b,c,d,L,S
 integer la,lb,lc,ld 
 real vdir,vex
 real sj ! 6-j symbol from libra
 real xl,xla,xlb,xlc,xld
 real fact
 real r2me   ! function

 allocate( rawme(numorb,numorb) )
! print*,scaling
 rawme(:,:) = 0.0
 do i = 1,numorb
     ni = orbqn(i)%nr
     li = orbqn(i)%l

     do j = 1, numorb
         nj = orbqn(j)%nr
         lj = orbqn(j)%l

         if(li == lj)cycle
         if( abs(li-lj) /= 1)cycle
         if( abs(ni-nj) > 1)cycle

         if( ni == nj )then
             if(li == lj+1)then
                rawme(i,j) = sqrt( float(lj+1)*(nj+lj+1.5) )
             else
                rawme(i,j) = -sqrt( float(lj)*(nj+lj+0.5) )
             endif
         else
              if(ni == nj -1 .and.li == lj+1)then
                rawme(i,j) = -sqrt(float((lj+1)*nj))
              endif
              if(ni == nj +1 .and.li == lj-1)then
                rawme(i,j) = sqrt(float(lj*(nj+1)))
              endif


         endif

     end do ! j
 enddo ! i
!------------ CYCLE OVER ALLOWED TBMEs------------

 do i = 1,nme
     a = vtbme(i)%n(1)
     b = vtbme(i)%n(2)
     c = vtbme(i)%n(3)
     d = vtbme(i)%n(4)
     L = vtbme(i)%J
     S = vtbme(i)%T
     la = orbqn(a)%l
     lb = orbqn(b)%l
     lc = orbqn(c)%l
     ld = orbqn(d)%l
     xL = float(L)
     xla =float(la)
     xlb =float(lb)
     xlc =float(lc)
     xld =float(ld)
     fact = 1.0
     if( a == b)fact = fact/sqrt(2.)
     if(c == d) fact = fact/sqrt(2.)
     vdir = (-1)**(L+lb+lc)*sj(xl,xlb,xla,1.,xlc,xld)*rawme(a,c)*rawme(b,d)
     vex  = (-1)**(L+lb+ld)*sj(xl,xlb,xla,1.,xld,xlc)*rawme(a,d)*rawme(b,c)
     vtbme(i)%v = vtbme(i)%v - fact*(vdir-(-1)**(l+s+lc+ld+1)*vex)*2  ! don't know
                                                                      ! about this x2

 enddo

!...... RESET RAWME to matrix elements of r2

 rawme(:,:) = 0.0
 do i = 1,numorb
     ni = orbqn(i)%nr
     li = orbqn(i)%l

     do j = 1, numorb
         nj = orbqn(j)%nr
         lj = orbqn(j)%l
         if(li /= lj)then
              rawme(i,j) = 0.0
              rawme(j,i) = 0.0
              cycle
         end if
         rawme(i,j) = r2me(nj,lj,ni,li)
         rawme(j,i) = rawme(i,j)
     end do
 end do
 call cnvt1to2(rawme,2)

 deallocate(rawme)

 return
 end subroutine Rrel_atom_lab
!==================================================

 subroutine R2atom(scaling)

 use sporbit
 use tbmes
 implicit none

 real scaling !

 real, allocatable::  rawme(:,:)
 integer i,j
 integer ni,li,nj,lj
 integer a,b,c,d,L,S
 integer la,lb,lc,ld 
 real vdir,vex
 real sj ! 6-j symbol from libra
 real xl,xla,xlb,xlc,xld
 real fact

 allocate( rawme(numorb,numorb) )
 print*,scaling
 rawme(:,:) = 0.0
 do i = 1,numorb
     ni = orbqn(i)%nr
     li = orbqn(i)%l

     do j = 1, numorb
         nj = orbqn(j)%nr
         lj = orbqn(j)%l

         if(li == lj)cycle
         if( abs(li-lj) /= 1)cycle
         if( abs(ni-nj) > 1)cycle

         if( ni == nj )then
             if(li == lj+1)then
                rawme(i,j) = sqrt( float(lj+1)*(nj+lj+1.5) )
             else
                rawme(i,j) = -sqrt( float(lj)*(nj+lj+0.5) )
             endif
         else
              if(ni == nj -1 .and.li == lj+1)then
                rawme(i,j) = -sqrt(float((lj+1)*nj))
              endif
              if(ni == nj +1 .and.li == lj-1)then
                rawme(i,j) = sqrt(float(lj*(nj+1)))
              endif


         endif

     end do ! j
 enddo ! i
!------------ CYCLE OVER ALLOWED TBMEs------------

 do i = 1,nme
     a = vtbme(i)%n(1)
     b = vtbme(i)%n(2)
     c = vtbme(i)%n(3)
     d = vtbme(i)%n(4)
     L = vtbme(i)%J
     S = vtbme(i)%T
     la = orbqn(a)%l
     lb = orbqn(b)%l
     lc = orbqn(c)%l
     ld = orbqn(d)%l
     xL = float(L)
     xla =float(la)
     xlb =float(lb)
     xlc =float(lc)
     xld =float(ld)
     fact = 1.0
     if( a == b)fact = fact/sqrt(2.)
     if(c == d) fact = fact/sqrt(2.)
     vdir = (-1)**(L+lb+lc)*sj(xl,xlb,xla,1.,xlc,xld)*rawme(a,c)*rawme(b,d)
     vex  = (-1)**(L+lb+ld)*sj(xl,xlb,xla,1.,xld,xlc)*rawme(a,d)*rawme(b,c)
     vtbme(i)%v = vtbme(i)%v + scaling*fact*(vdir-(-1)**(l+s+lc+ld+1)*vex)

 enddo


 deallocate(rawme)
 return
 end subroutine R2atom
!==================================================

 subroutine P2atom(scaling)

 use sporbit
 use tbmes
 implicit none

 real scaling !

 real, allocatable::  rawme(:,:)
 integer i,j
 integer ni,li,nj,lj
 integer a,b,c,d,L,S
 integer la,lb,lc,ld 
 real vdir,vex
 real sj ! 6-j symbol from libra
 real xl,xla,xlb,xlc,xld
 real fact

 allocate( rawme(numorb,numorb) )

 do i = 1,numorb
     ni = orbqn(i)%nr
     li = orbqn(i)%l

     do j = 1, numorb
         rawme(i,j) = 0.0
         nj = orbqn(j)%nr
         lj = orbqn(j)%l

         if(li == lj)cycle
         if( abs(li-lj) /= 1)cycle
         if( abs(ni-nj) > 1)cycle

         if( ni == nj )then
             if(li == lj+1)then
                rawme(i,j) = -sqrt( float(lj+1)*(nj+lj+1.5) )
             else
                rawme(i,j) = -sqrt(float(lj)*(nj+lj+0.5))
             endif
         else
              if(ni == nj -1 .and.li == lj+1)then
                rawme(i,j) = -sqrt(float((lj+1)*nj))
              endif
              if(ni == nj +1 .and.li == lj-1)then
                rawme(i,j) = -sqrt(float(lj*(nj+1)))
              endif


         endif!'

     end do ! j
 enddo ! i

!------------ CYCLE OVER ALLOWED TBMEs------------

 do i = 1,nme
     a = vtbme(i)%n(1)
     b = vtbme(i)%n(2)
     c = vtbme(i)%n(3)
     d = vtbme(i)%n(4)
     L = vtbme(i)%J
     S = vtbme(i)%T
     la = orbqn(a)%l  
     lb = orbqn(b)%l
     lc = orbqn(c)%l
     ld = orbqn(d)%l
     xL = float(L)
     xla =float(la)
     xlb =float(lb)
     xlc =float(lc)
     xld =float(ld)
     fact = 1.0
     if( a == b)fact = fact/sqrt(2.)
     if(c == d) fact = fact/sqrt(2.)

     vdir = (-1)**(L+lb+lc)*sj(xl,xlb,xla,1.,xlc,xld)*rawme(a,c)*rawme(b,d)
     vex  = (-1)**(L+lb+ld)*sj(xl,xlb,xla,1.,xld,xlc)*rawme(a,d)*rawme(b,c)
!     print*,vdir,vex,vtbme(i)%v

     vtbme(i)%v = vtbme(i)%v - scaling*fact*(vdir-(-1)**(l+s+lc+ld+1)*vex)

 enddo


 deallocate(rawme)
 return
 end subroutine P2atom

!========================================================

 subroutine speCMatom_diag(A)

 use sporbit
 use tbmes
 implicit none
 integer i
 integer n,l
 integer A

 do i = 1,numorb
     n = orbqn(i)%nr
     l = orbqn(i)%l

     spe(i) = spe(i)+(2*n+l+1.5)/A

 enddo ! i

 return
 end subroutine speCMatom_diag

!========================================================
 subroutine Hcm_nuke

 use sporbit
 use tbmes
 implicit none

 logical speflag
 real b
 integer A

 print*,' '
 print*,' Enter total number of particles '
 read*,A

 print*,' Enter osc. parameter b of c.m. relative to b for trap '
 read*,b

 if(b <= 0.0)then
   print*,' Problem with b ',b
   stop
 endif

 call R2nuke(b**2/A)
 call P2nuke(1./A/b**2)

 speflag = .true.

 if(b /= 1.0)speflag = .false.

 if(speflag)then
      call speCMnuke_diag(A)
 else
      print*,' not ready yet ' 
 endif
 
 print*,' '
 print*,' ATTENTION!  If you want to create (Hcm - 3/2 hw)'
 print*,' create the number operator (option 21) '
 print*,' and scale by - ',A*1.5

 return

 end subroutine Hcm_nuke

!========================================================
 subroutine Hcm_nuke2

 use sporbit
 use Vrelative
 implicit none

 real b
 integer A

 print*,' '
 print*,' Enter total number of particles '
 read*,A

 print*,' Enter osc. parameter b of c.m. relative to b for trap '
 read*,b

 if(b <= 0.0)then
   print*,' Problem with b ',b
   stop
 endif
 print*,' Which Moshinsky phase do you want? (0= GB, 1 = std) '
 read*,moshphase

 call homaster2(A,b)
 call setVrelcm(b,A)
  call binom2
  call trinom
 call Vlabfromrel
 return

 end subroutine Hcm_nuke2

!==================================================
 subroutine Trel_nuke2

 use sporbit
 use Vrelative
 implicit none

 real b
 integer A

 print*,' '
 print*,' Enter total number of particles '
 read*,A

 print*,' Enter osc. parameter b of c.m. relative to b for trap '
 read*,b

 if(b <= 0.0)then
   print*,' Problem with b ',b
   stop
 endif
 print*,' Which Moshinsky phase do you want? (0= GB, 1 = std) '
 read*,moshphase

 call setTrelcm(b,A)
  call binom2
  call trinom
 call Vlabfromrel
 return

 end subroutine Trel_nuke2

!==================================================
!========================================================
 subroutine Trel_nuke

 use sporbit
 use tbmes
 implicit none

 logical speflag
 real b
 integer A

 print*,' '
 print*,' Enter total number of particles '
 read*,A

 print*,' Enter osc. parameter b of c.m. relative to b for trap '
 read*,b

 if(b <= 0.0)then
   print*,' Problem with b ',b
   stop
 endif

 call P2nuke(1./A/b**2)



 return

 end subroutine Trel_nuke

!========================================================
!==================================================

 subroutine R2nuke(scaling)

 use sporbit
 use tbmes
 implicit none

 real scaling !

 real, allocatable::  rawme(:,:)
 integer i
 integer ni,li,nj,lj
 integer a,b,c,d,J,T
 integer la,lb,lc,ld 
 integer ja,jb,jc,jd
 real xja,xjb,xjc,xjd
 real vdir,vex
 real sj ! 6-j symbol from libra
 real xJ,xla,xlb,xlc,xld
 real fact

 allocate( rawme(numorb,numorb) )
 print*,scaling
 rawme(:,:) = 0.0
 do i = 1,numorb
     ni = orbqn(i)%nr
     li = orbqn(i)%l

     do j = 1, numorb
         nj = orbqn(j)%nr
         lj = orbqn(j)%l

         if(li == lj)cycle
         if( abs(li-lj) /= 1)cycle
         if( abs(ni-nj) > 1)cycle

         if( ni == nj )then
             if(li == lj+1)then
                rawme(i,j) = sqrt( float(lj+1)*(nj+lj+1.5) )
             else
                rawme(i,j) = -sqrt( float(lj)*(nj+lj+0.5) )
             endif
         else
              if(ni == nj -1 .and.li == lj+1)then
                rawme(i,j) = -sqrt(float((lj+1)*nj))
              endif
              if(ni == nj +1 .and.li == lj-1)then
                rawme(i,j) = sqrt(float(lj*(nj+1)))
              endif


         endif

     end do ! j
 enddo ! i
!------------ CYCLE OVER ALLOWED TBMEs------------

 do i = 1,nme
     a = vtbme(i)%n(1)
     b = vtbme(i)%n(2)
     c = vtbme(i)%n(3)
     d = vtbme(i)%n(4)
     J = vtbme(i)%J
     T = vtbme(i)%T
     la = orbqn(a)%l
     lb = orbqn(b)%l
     lc = orbqn(c)%l
     ld = orbqn(d)%l
     xJ = float(J)
     xla =float(la)
     xlb =float(lb)
     xlc =float(lc)
     xld =float(ld)
     ja = orbqn(a)%j
     jb = orbqn(b)%j
     jc = orbqn(c)%j
     jd = orbqn(d)%j
     xja =float(ja)/2.
     xjb =float(jb)/2.
     xjc =float(jc)/2.
     xjd =float(jd)/2.

     fact = 1.0
     if( a == b)fact = fact/sqrt(2.)
     if(c == d) fact = fact/sqrt(2.)
     vdir = (-1)**(J+la+lb+(jb+jd)/2)*sj(xJ,xjb,xja,1.,xjc,xjd) & 
 * sj(xla,xja,0.5,xjc,xlc,1.)*sj(xlb,xjb,0.5,xjd,xld,1.0)   &
 *  sqrt(float( (ja+1)*(jb+1)*(jc+1)*(jd+1))) * rawme(a,c)*rawme(b,d)

     vex = (-1)**(J+la+lb+(jb+jc)/2)*sj(xJ,xjb,xja,1.,xjd,xjc) & 
 * sj(xla,xja,0.5,xjd,xld,1.)*sj(xlb,xjb,0.5,xjc,xlc,1.0)   &
 *  sqrt(float( (ja+1)*(jb+1)*(jc+1)*(jd+1))) * rawme(a,d)*rawme(b,c)
     vtbme(i)%v = vtbme(i)%v + scaling*fact*(vdir-(-1)**(J+T+(jc+jd)/2+1)*vex)

!     write(6,'(6i3,4f8.4)')a,b,c,d,j,t,scaling*fact*vdir, scaling*fact*vex &
! , rawme(a,d), rawme(b,c)

 enddo


 deallocate(rawme)
 return
 end subroutine R2nuke
!==================================================

 subroutine P2nuke(scaling)

 use sporbit
 use tbmes
 implicit none

 real scaling !

 real, allocatable::  rawme(:,:)
 integer i,j
 integer ni,li,nj,lj
 integer a,b,c,d,T
 integer la,lb,lc,ld 
 real vdir,vex
 real sj ! 6-j symbol from libra
 real xJ,xla,xlb,xlc,xld
 integer ja,jb,jc,jd
 real xja,xjb,xjc,xjd
 real fact

 allocate( rawme(numorb,numorb) )

 do i = 1,numorb
     ni = orbqn(i)%nr
     li = orbqn(i)%l

     do j = 1, numorb
         rawme(i,j) = 0.0
         nj = orbqn(j)%nr
         lj = orbqn(j)%l

         if(li == lj)cycle
         if( abs(li-lj) /= 1)cycle
         if( abs(ni-nj) > 1)cycle

         if( ni == nj )then
             if(li == lj+1)then
                rawme(i,j) = -sqrt( float(lj+1)*(nj+lj+1.5) )
             else
                rawme(i,j) = -sqrt(float(lj)*(nj+lj+0.5))
             endif
         else
              if(ni == nj -1 .and.li == lj+1)then
                rawme(i,j) = -sqrt(float((lj+1)*nj))
              endif
              if(ni == nj +1 .and.li == lj-1)then
                rawme(i,j) = -sqrt(float(lj*(nj+1)))
              endif


         endif!'

     end do ! j
 enddo ! i

!------------ CYCLE OVER ALLOWED TBMEs------------

 do i = 1,nme
     a = vtbme(i)%n(1)
     b = vtbme(i)%n(2)
     c = vtbme(i)%n(3)
     d = vtbme(i)%n(4)
     J = vtbme(i)%J
     T = vtbme(i)%T
     la = orbqn(a)%l
     lb = orbqn(b)%l
     lc = orbqn(c)%l
     ld = orbqn(d)%l
     xJ = float(J)
     xla =float(la)
     xlb =float(lb)
     xlc =float(lc)
     xld =float(ld)
     ja = orbqn(a)%j
     jb = orbqn(b)%j
     jc = orbqn(c)%j
     jd = orbqn(d)%j
     xja =float(ja)/2.
     xjb =float(jb)/2.
     xjc =float(jc)/2.
     xjd =float(jd)/2.

     fact = 1.0
     if( a == b)fact = fact/sqrt(2.)
     if(c == d) fact = fact/sqrt(2.)
     vdir = (-1)**(J+la+lb+(jb+jd)/2)*sj(xJ,xjb,xja,1.,xjc,xjd) & 
 * sj(xla,xja,0.5,xjc,xlc,1.)*sj(xlb,xjb,0.5,xjd,xld,1.0)   &
 *  sqrt(float( (ja+1)*(jb+1)*(jc+1)*(jd+1))) * rawme(a,c)*rawme(b,d)

     vex = (-1)**(J+la+lb+(jb+jc)/2)*sj(xJ,xjb,xja,1.,xjd,xjc) & 
 * sj(xla,xja,0.5,xjd,xld,1.)*sj(xlb,xjb,0.5,xjc,xlc,1.0)   &
 *  sqrt(float( (ja+1)*(jb+1)*(jc+1)*(jd+1))) * rawme(a,d)*rawme(b,c)
     vtbme(i)%v = vtbme(i)%v - scaling*fact*(vdir-(-1)**(J+T+(jc+jd)/2+1)*vex)

 enddo


 deallocate(rawme)
 return
 end subroutine P2nuke

!========================================================

 subroutine speCMnuke_diag(A)

 use sporbit
 use tbmes
 implicit none
 integer i
 integer n,l
 integer A

 do i = 1,numorb
     n = orbqn(i)%nr
     l = orbqn(i)%l

     spe(i) = spe(i)+(2*n+l+1.5)/A

 enddo ! i

 return
 end subroutine speCMnuke_diag

!========================================================

      subroutine homaster(A,b)
      use sporbit
      implicit none

      integer A
      real b
      real spemat(numorb,numorb)



      spemat = 0.0
      call hoonebody(b,spemat,1./float(A))
      call cnvt1to2(spemat,A)
      return

      end subroutine  homaster

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine hoonebody(b,spemat,scale)

      use sporbit
      implicit none
      real spemat(numorb,numorb)
      real b
      real r2me ,p2me
      integer i1,i2
      integer n1,l1,n2,l2,j1,j2
      real scale
      integer factor  ! = 2 - delta(a,b)
      

      do i1 = 1,numorb
        n1 = orbqn(i1)%nr
        l1 = orbqn(i1)%l
        j1 = orbqn(i1)%j
        do i2 = 1,numorb
           n2 = orbqn(i2)%nr
           l2 = orbqn(i2)%l
           j2 = orbqn(i2)%j
           if(l1 == l2 .and. j1 == j2 )then
           spemat(i1,i2) = spemat(i1,i2)+(0.5*(b*b*r2me(n1,l1,n2,l2)+ &
                                  p2me(n1,l1,n2,l2)/(b*b)))*scale
           endif
        enddo
      enddo

      return
      end subroutine hoonebody

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function r2me(n1,l1,n2,l2)
!
!  matrix elements of r2 between h.o. wfns
!   < n1, l1 | r^2 | n2, l2 >
!  INPUT:
!    n1,l1,n2,l2  : n = nodal quantum number (0, 1, 2,...)
!
      implicit none
      integer n1,l1,n2,l2

      r2me = 0.0

      if(n1 == n2 .and. l1 == l2)then
        r2me =(2*n1 +l1 + 1.5)
      endif

      if(n1 == n2 .and. l2 == l1 +2)then  ! ?? not confirmed
        r2me = sqrt((n1+l2 + 1.5)*(n1+l2+2.5))
      endif
      if(n1 == n2 .and. l2 == l1 -2)then
        r2me = sqrt((n1+l1 + 1.5)*(n1+l1+2.5))
      endif

      if(l1 == l2 .and. n2 == n1+1)then
        r2me = -sqrt(n2*(n2+l1+0.5))
      endif
      if(l1 == l2 .and. n2 == n1-1)then
        r2me = -sqrt(n1*(n1+l1+0.5))
      endif

      if(l2 == l1+2 .and. n2 == n1 -1)then  ! ?? not confirmed
        r2me = -sqrt(2*n1+(2*n1+2*l1+3.0))
      endif
      if(l2 == l1-2 .and. n2 == n1 +1)then
        r2me = -sqrt(2*n2+(2*n2+2*l2+3.0))
      endif
      write(17,*)n1,n2,r2me
      return
      end function r2me

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      real function p2me(n1,l1,n2,l2)
!
!  matrix elements of r2 between h.o. wfns
!   < n1, l1 | p^2 | n2, l2 >
!  INPUT:
!    n1,l1,n2,l2  : n = nodal quantum number (0, 1, 2,...)
!
      implicit none
      integer n1,l1,n2,l2

      p2me = 0.0
      if(n1 == n2 .and. l1 == l2)then
        p2me =(2*n1 +l1 + 1.5)
      endif
      if(n1 == n2 .and. l2 == l1 +2)then ! ?? not confirmed
        p2me = -sqrt((n1+l2 + 1.5)*(n1+l2+2.5))
      endif
      if(n1 == n2 .and. l2 == l1 -2)then
        p2me = -sqrt((n1+l1 + 1.5)*(n1+l1+2.5))
      endif

      if(l1 == l2 .and. n2 == n1+1)then
        p2me = sqrt(n2*(n2+l1+0.5))
      endif
      if(l1 == l2 .and. n2 == n1-1)then
        p2me = sqrt(n1*(n1+l1+0.5))
      endif

      if(l2 == l1+2 .and. n2 == n1 -1)then ! ?? not confirmed
        p2me = sqrt(2*n1+(2*n1+2*l1+3.0))
      endif
      if(l2 == l1-2 .and. n2 == n1 +1)then
        p2me = sqrt(2*n2+(2*n2+2*l2+3.0))
      endif

      return
      end function p2me


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cnvt1to2(spemat,Np)
!
!  basic routine for converting from 1-body to 2-body
!

      use sporbit
      use tbmes
      implicit none
      real spemat(numorb,numorb)

      integer Np
      integer n
      integer a,b,c,d
      integer j,t
      logical isocheck
      integer icount
      integer iphase
      real zetaab,zetacd
      integer la,lb,lc,ld
      icount = 0
      do n = 1,nme
        a = vtbme(n)%n(1)
        b = vtbme(n)%n(2)
        c = vtbme(n)%n(3)
        d = vtbme(n)%n(4)

        J = vtbme(n)%j
        T = vtbme(n)%t
        
        if( a /= c .and. a /= d .and. b /= c .and. b/=d)cycle  ! no repeated indices
!        print*,n,a,b,c,d
!        print*,' made it '
        la = orbqn(a)%l
        lb = orbqn(b)%l
        lc = orbqn(c)%l
        ld = orbqn(d)%l

        zetaab = 1.
        zetacd = 1.
        if(a==b)zetaab = 1./sqrt(2.)
        if(c==d)zetacd = 1./sqrt(2.)

        if( b == d)then
           vtbme(n)%v = vtbme(n)%v+ zetaab*zetacd*spemat(a,c)/float(Np -1)
        endif   ! b == d 

        if( a == d)then
                  iphase = (-1)**(j+la + lb + t)

           vtbme(n)%v = vtbme(n)%v+ &
                      iphase*zetaab*zetacd*spemat(b,c)/float(Np -1)

        endif   ! a == d 

        if( b == c)then

                  iphase = (-1)**(j+lc + ld + t)

           vtbme(n)%v = vtbme(n)%v+ iphase*zetaab*zetacd*spemat(a,d)/float(Np -1)

        endif   ! b == c

        if( a == c)then

           vtbme(n)%v = vtbme(n)%v+ zetaab*zetacd*spemat(b,d)/float(Np -1)

        endif   ! a == c
      enddo
      return
      end


!========================================================
!
! used when converting from rel
!
      subroutine homaster2(A,b)
      use sporbit
      implicit none

      integer A
      real b
      real spemat(numorb,numorb)



      spemat = 0.0
      call hoonebody(b,spemat,1.)
      call cnvt1to2(spemat,A)
      return

      end subroutine  homaster2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  used when converting from rel
! 
      subroutine hoonebody2(b,spemat,A)

      use sporbit
      implicit none
      real spemat(numorb,numorb)
      real b
      real r2me ,p2me
      integer i1,i2
      integer n1,l1,n2,l2,j1,j2
      integer A
      integer factor  ! = 2 - delta(a,b)
      

      do i1 = 1,numorb
        n1 = orbqn(i1)%nr
        l1 = orbqn(i1)%l
        j1 = orbqn(i1)%j
        do i2 = 1,numorb
           n2 = orbqn(i2)%nr
           l2 = orbqn(i2)%l
           j2 = orbqn(i2)%j
           if(l1 == l2 .and. j1 ==  j2)then
           spemat(i1,i2) = spemat(i1,i2)+(0.5*(b*b*r2me(n1,l1,n2,l2)+ &
                                  p2me(n1,l1,n2,l2)/(b*b)))

           endif
        enddo
      enddo
      return
      end subroutine hoonebody2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 subroutine Vlabfromrel
 use sporbit
 use tbmes
 implicit none
 integer i
 integer a,b,c,d
 integer L,S
 real ans
!$OMP PARALLEL DO PRIVATE(a,b,c,d,L,S,ans)
 do i = 1,nme
     a = vtbme(i)%n(1)
     b = vtbme(i)%n(2)
     c = vtbme(i)%n(3)
     d = vtbme(i)%n(4)
     L = vtbme(i)%J
     S = vtbme(i)%T
     if(spinless)then
        call vALSfromrel(a,b,c,d,L,S,ans)
     else
        call vAjjfromrel(a,b,c,d,L,S,ans)

     endif
     vtbme(i)%v = vtbme(i)%v+ans
 enddo

 return
 end subroutine Vlabfromrel

!========================================================
!
!  sets Vrel matrix elements
!

      subroutine setVrelcm(b,A)

      use sporbit
      use Vrelative
      implicit none
      integer A
      real b
      integer lrel,s
      real p2me,r2me
      integer nf,ni
      integer i
      integer emax

      if(allocated(vrel))deallocate(vrel)

!.............. ESTIMATE MAX L, MAX N, BY GOING THROUGH S.P. ORBITS

      lmax = 0
      Emax = 0
      do i = 1,numorb
          lmax = max(lmax,orbqn(i)%l)
          Emax = max(Emax,2*orbqn(i)%nr+orbqn(i)%l)
      end do ! i
      lmax = lmax*2 +2 
      nmax = Emax +2 

!      print*,' Enter max rel L '
!      read*,lmax

      allocate(vrel(0:lmax,0:1))

!      print*,' Enter max n in relative '
!      read*,nmax

      do lrel = 0,lmax
          do s = 0,1
             allocate(vrel(lrel,s)%v(0:nmax,0:nmax))
!................ NOW FILLL
             do nf = 0,nmax
                 do ni = nf,nmax
                     vrel(lrel,s)%v(nf,ni) = -((b*b*r2me(nf,lrel,ni,lrel)+ &
                                  p2me(nf,lrel,ni,lrel)/(b*b) ))/float(A)
                     vrel(lrel,s)%v(ni,nf) = vrel(lrel,s)%v(nf,ni)
                 end do ! ni
             end do ! nf
          enddo
      enddo  !lrel
!      print*,vrel(0,0)%v(:,:)
      return

      end subroutine setVrelcm


!========================================================
!========================================================
!
!  sets Vrel matrix elements
!

      subroutine setTrelcm(b,A)

      use sporbit
      use Vrelative
      implicit none
      integer A
      real b
      integer lrel,s
      real p2me,r2me
      integer nf,ni
      integer i
      integer emax

      if(allocated(vrel))deallocate(vrel)

!.............. ESTIMATE MAX L, MAX N, BY GOING THROUGH S.P. ORBITS

      lmax = 0
      Emax = 0
      do i = 1,numorb
          lmax = max(lmax,orbqn(i)%l)
          Emax = max(Emax,2*orbqn(i)%nr+orbqn(i)%l)
      end do ! i
      lmax = lmax*2 +2 
      nmax = Emax +2 

!      print*,' Enter max rel L '
!      read*,lmax

      allocate(vrel(0:lmax,0:1))

!      print*,' Enter max n in relative '
!      read*,nmax

      do lrel = 0,lmax
          do s = 0,1
             allocate(vrel(lrel,s)%v(0:nmax,0:nmax))
!................ NOW FILLL
             do nf = 0,nmax
                 do ni = nf,nmax
                     vrel(lrel,s)%v(nf,ni) = (( &
                                  p2me(nf,lrel,ni,lrel)/(b*b) ))*.5 !/float(A)
                     vrel(lrel,s)%v(ni,nf) = vrel(lrel,s)%v(nf,ni)
                 end do ! ni
             end do ! nf
          enddo
      enddo  !lrel
!      print*,vrel(0,0)%v(:,:)
      return

      end subroutine setTrelcm


!========================================================
!
!  sets Vrel matrix elements
!

      subroutine setVrelR2(b,A)

      use sporbit
      use Vrelative
      implicit none
      integer A
      real b
      integer lrel,s
      real p2me,r2me
      integer nf,ni
      integer i
      integer emax

      if(allocated(vrel))deallocate(vrel)

!.............. ESTIMATE MAX L, MAX N, BY GOING THROUGH S.P. ORBITS

      lmax = 0
      Emax = 0
      do i = 1,numorb
          lmax = max(lmax,orbqn(i)%l)
          Emax = max(Emax,2*orbqn(i)%nr+orbqn(i)%l)
      end do ! i
      lmax = lmax*2 +2 
      nmax = Emax +2 

!      print*,' Enter max rel L '
!      read*,lmax

      allocate(vrel(0:lmax,0:1))

!      print*,' Enter max n in relative '
!      read*,nmax

      do lrel = 0,lmax
          do s = 0,1
             allocate(vrel(lrel,s)%v(0:nmax,0:nmax))
!................ NOW FILLL
             do nf = 0,nmax
                 do ni = nf,nmax
                     vrel(lrel,s)%v(nf,ni) = 2.*b*b*r2me(nf,lrel,ni,lrel)
                     vrel(lrel,s)%v(ni,nf) = vrel(lrel,s)%v(nf,ni)
                 end do ! ni
             end do ! nf
          enddo
      enddo  !lrel
!      do nf = 0,nmax
!         print*,vrel(0,0)%v(nf,:)
!      end do
      return

      end subroutine setVrelR2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine vAjjfromrel(a,b,c,d,J,T,ans)
      use sporbit

      implicit none
      integer a,b,c,d
      integer J,T
      real ans
      real vdir, vex
      real fact
      integer jc,jd,iphase

      fact = 1.0
      if(a == b) fact = fact/sqrt(2.)
      if(c == d) fact = fact/sqrt(2.)

      jc = orbqn(c)%j
      jd = orbqn(d)%j

      iphase = (-1)**(J+T+1+(jc+jd)/2)
      call vNSjjfromrel(a,b,c,d,J,T,vdir)
      call vNSjjfromrel(a,b,d,c,J,T,vex)

      ans = fact*(vdir -iphase*vex)
      return

      end subroutine vAjjfromrel

!========================================================


      subroutine vALSfromrel(a,b,c,d,L,S,ans)
      use sporbit
      use Vrelative
      implicit none
      integer a,b,c,d
      integer L,S
      real ans
      real vdir, vex, vtmp
      real fact
      integer iphase
      integer lrel

      integer na,nb,nc,nd
      integer la,lb,lc,ld

      la = orbqn(a)%l
      lb = orbqn(b)%l
      lc = orbqn(c)%l
      ld = orbqn(d)%l


      na = orbqn(a)%nr
      nb = orbqn(b)%nr
      nc = orbqn(c)%nr
      nd = orbqn(d)%nr

      fact = 1.0
      if(a == b) fact = fact/sqrt(2.)
      if(c == d) fact = fact/sqrt(2.)

      iphase = (-1)**(L+S+1+(lc+ld))
      vdir = 0.0
      vex = 0.0
      do lrel = 0,lmax
           if( (-1)**(S+lrel) == -1)cycle  ! must be antisym in relative
           call gen_rel_to_lab(na,la,nb,lb,nc,lc,nd,ld,L,lrel,S,vtmp)
           vdir = vdir + vtmp
           call gen_rel_to_lab(na,la,nb,lb,nd,ld,nc,lc,L,lrel,S,vtmp)
           vex = vex + vtmp
      end do ! lrel
      ans = fact*(vdir -iphase*vex)
!      print*,a,b,c,d,l,s,vdir,vex,ans

      return

      end subroutine vAlsfromrel

!========================================================
!
!  compute non-antisymmeterized jj coupled ME from relative

      subroutine vNSjjfromrel(a,b,c,d,J,T,ans)

      use sporbit
      use Vrelative
      implicit none

      real ans,tmp,tmpans

      integer a,b,c,d
      integer na,nb,nc,nd
      integer la,lb,lc,ld
      integer ja,jb,jc,jd
      integer J,S,T
      integer L,lrel
      real*8 :: c9j

      ja = orbqn(a)%j
      jb = orbqn(b)%j
      jc = orbqn(c)%j
      jd = orbqn(d)%j

      la = orbqn(a)%l
      lb = orbqn(b)%l
      lc = orbqn(c)%l
      ld = orbqn(d)%l

      na = orbqn(a)%nr
      nb = orbqn(b)%nr
      nc = orbqn(c)%nr
      nd = orbqn(d)%nr

      ans = 0.0
      do S = 0,1
      do L = abs(J-S),abs(J+S)
        tmpans = 0.0
        do lrel = 0,lmax
           if( (-1)**(S+T+lrel) == 1)cycle  !must be antisym in rel frame
           call gen_rel_to_lab(na,la,nb,lb,nc,lc,nd,ld,L,lrel,S,tmp)
           tmpans = tmpans + tmp
        end do  ! lrel
        ans = ans + tmpans*sqrt( float((ja+1)*(jb+1)*(jc+1)*(jd+1)))  & 
              * (2*L+1)*(2*J+1)*    c9j(la*2,1,ja,lb*2,1,jb,2*L,2*S,2*J)* &
         c9j(lc*2,1,jc,ld*2,1,jd,2*L,2*S,2*J)
      end do ! l
      end do ! s
      return
      end subroutine vNSjjfromrel
!====================================================

      subroutine gen_rel_to_lab(na,la,nb,lb,nc,lc,nd,ld,L,lrel,S,ans)

      use sporbit 
      use vrelative
      implicit none
      integer na,nb,nc,nd
      integer la,lb,lc,ld,L,lrel,s
      real ans
!................ INTERMEDIATE

      integer Ea,Eb,Ec,Ed
      integer Eab, Ecd
      integer Ecm,minEcm,maxEcm
      integer Ncm, minNcm,maxNcm
      integer Ei,Ef,ni,nf
      integer Lcm, minLcm, maxLcm
      real*8 tmb,tmb0  ! function Talmi-Moshinsky bracket
      real,pointer::v(:,:)

      real tmbi,tmbf

      ans = 0.0
      if(lrel > lmax)return
      v => vrel(lrel,s)%v

!................ error trap--only for contact interactions

!      if( ((la + lb +L)/2)*2 /= la + lb +L)return
!      if( ((lc + ld +L)/2)*2 /= lc + ld +L)return

!..........................
      Ea = 2*na + la
      Eb = 2*nb + lb         
      Ec = 2*nc + lc
      Ed = 2*nd + ld

      
      minEcm = 0
      maxEcm = max(Ea+Eb,Ec+Ed)

      ans = 0.0
      maxLcm = max(la+lb,lc+ld)
      maxLcm = max(maxLcm,lrel+L)
      minLcm = 0
!      minLcm = max(abs(la-lb),abs(lc-ld))
!      minLcm = max(maxLcm,abs(lrel-L))
!      write(37,*)L,S,lrel,minLcm,maxLcm
      if(minLcm > maxLcm) return
      
      do Lcm = minLcm, maxLcm
!      write(37,*)' lcm ',lcm,maxEcm
      if(Lcm > maxEcm)cycle
      do Ecm = minEcm,maxEcm
         if(Ecm -Lcm < 0)cycle
         ef = ea+eb -Ecm 
         ei = ec+ed-Ecm 
         if(ef < 0 .or. ei < 0) cycle! then   
!             print*,' woopps ',ef,ei
!             print*,ea,eb,ec,ed,Ecm
!             stop
!         endif
         if(ef < lrel .or. ei < lrel)cycle
         nf = (ef-lrel) /2
         ni = (ei-lrel) /2
         if(nf > nmax .or. ni > nmax)cycle
!         print*,Ecm,Lcm, Ef,lrel,ea,la,eb,lb,L
         if(moshphase == 0)then
         tmbf = tmb0(Ecm, Lcm, Ef, lrel, ea,la,eb,lb,L,1.d0)
         tmbi = tmb0(Ecm, Lcm, Ei, lrel, ec,lc,ed,ld,L,1.d0)
         else
         tmbf = tmb(Ecm, Lcm, Ef, lrel, ea,la,eb,lb,L,1.d0)
         tmbi = tmb(Ecm, Lcm, Ei, lrel, ec,lc,ed,ld,L,1.d0)
         endif
         if(s==0 .and.lrel==0)write(6,'(6i3,f8.3)')ea,eb,ec,ed,ni,nf,tmbf*tmbi
         ans = ans + tmbf*tmbi*v(ni,nf)
!         print*,tmbf,tmbi,v(ni,nf)
!         write(37,'(9i4)')ecm,lcm,ef,lrel,ea,la,eb,lb,l
!         write(37,*)' sum ',tmbf,tmbi,v(ni,nf)
      enddo ! Ecm
      end do ! lcm
      return
      end subroutine gen_rel_to_lab
!==============================================================
!
!  computes for harmonic oscilator basis
!  non-antisymmeterized < na la, nb lb, L | r1 * r2 | nc lc, nd lc, L >
!
!
      subroutine vr1dotr2NA(na,la,nb,lb,nc,lc,nd,ld,ltot,vme)

      implicit none
!......... INPUTS
      integer :: na, la, nb, lb, nc, lc, nd, ld, ltot
!..........OUTPUT
      real    :: vme
!........ INTERNAL

!......... FUNCTIONS CALLED
      real SJ2I   ! called from LIBRA.F
      real rvecredme

      vme = (-1)**(lc+lb+Ltot)*SJ2I(2*Ltot,2*lb,2*la, 2, 2*lc, 2*ld)  & 
       * rvecredme(na,la,nc,lc)*rvecredme(nb,lb,nd,ld)

      return
      end subroutine vr1dotr2NA

!==============================================================
!
! function to generate the reduced matrix element
! < nf, lf || vec(r) || ni, li >
! in harmonic oscillator basis

      function rvecredme(nf,lf,ni,li)
      implicit none

      real rvecredme
      integer nf, lf, ni,li

      rvecredme = 0.0
      if (ni == nf .and. lf == li+1)then
         rvecredme =sqrt( float(li+1)*(ni+li+1.5) )
      end if
      if (ni == nf .and. lf == li-1)then
         rvecredme =-sqrt( float(li)*(ni+li+.5) )
      end if
      if (ni+1 == nf .and. lf == li-1)then
         rvecredme =sqrt( float(li*(ni+1)) )
      end if
      if (ni-1 == nf .and. lf == li+1)then
         rvecredme =-sqrt( float(li+1)*float(ni) )
      end if

      return
      end function rvecredme

!==============================================================
