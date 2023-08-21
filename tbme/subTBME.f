CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALVINCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine TBMEcalc
C
C  computes TBME's for a variety of interactions
C
C
C  CALLS:
C	GAUSS: generates gaussian random numbers for noise 
C
C

      use sporbit
      use tbmes
      implicit none

C---- SINGLE-PARTICLE LIST ------


      real xj1,xj2,xj3,xj4
      real xl1,xl2,xl3,xl4
      real n1,n2,n3,n4
C      logical active(100) 		! which s.p. states are active

C----TBME LIST -------


      real V 				! strength of any given interactoin 
      
C      real mspe(100),mme(5000)
      
C----MULTIPOLE INFORMATION 

C      real b
C      integer L,T
C      logical phase_convention 
      
C---- RANDOM NOISE -----
C      integer iseed
C      real gauss

C.......................MULTIPOLE INFORMATION..........

C      integer sizeme
C      parameter (sizeme = 12)
C      real opme(sizeme,sizeme)

C..................FILE CONTROL...................

      character*1 ychar

      character filename*15  		! name of  file 
      integer ilast 
      logical describe 		!flag to write out parameters
      character name*60 

C----- DUMMY COUNTERS -----

      integer i,k,n,iact
      integer ichoice
      real temp

C BEGIN PROGRAM  BEGIN PROGRAM  BEGIN PROGRAM  BEGIN PROGRAM  BEGIN PROGRAM 

      if(.not.allocated(spe))allocate(spe(numorb))

      write(6,*)' '
      write(6,*)' -------------------------------------------------' 
      write(6,*)' |                                               |'
      write(6,*)' |   Now to describe the two-body interactions   |' 
      write(6,*)' |                                               |'
      write(6,*)' -------------------------------------------------' 
      write(6,*)' '

      
      write(6,*)' Do you want a record of your choices written out? '
      read(5,'(a)')ychar 
      
      
      if(ychar.eq.'y' .or. ychar .eq. 'Y')then
      	describe = .true.
      	

      write(6,*)' Enter record name (.descript) '
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

      open(unit=7,file=filename(1:ilast)//'.descript',status='new')
C      write(6,'('' Shell-Model space file name'',1x,a15)')spsfil
      write(6,*)' Enter a brief title or description of interaction '
      read(5,'(a)')name
      write(7,*)name
      else
      	describe = .false.
      endif
   
      
    1 continue

      write(6,*)' '
      write(6,*)' Menu of currently available interactions: '
      write(6,*)' (0) random noise (TBRE, two-body random ensemble) '
      write(6,*)' (1) Pairing '
      write(6,*)' (2) Multipole-multipole (you choose L) '
      write(6,*)' (3) S^2 (total spin)   '
      write(6,*)' (4) L^2 (total orbital angular momentum ) '
      write(6,*)' (5) J^2 (total angular momentum ) '
      write(6,*)' (6) L*S (spin-orbit) 1+2 body '
      write(6,*)' (7) "monopole" (n_a x n_b) '
      write(6,*)' (8) constant displacement '
      write(6,*)' (9) random multipole '
      write(6,*)' (10)  T^2 (isospin) '
      write(6,*)' (11) contact interaction '
      write(6,*)' (13) special 2-level interaction '
      if(spinless)write(6,*)' (14) Hcm for cold atoms '
      if(spinless)write(6,*)' (16) R^2 (intrinsic) for cold atoms '
      if(.not.spinless)write(6,*)' (15) Hcm for nucleons '
      if(.not.spinless)write(6,*)' (17) Prel**2 for nucleons '
      write(6,*)' (18) S^2 T^2 '
      write(6,*)' (19) isoscalar pairing '
	  write(6,*)' (20) Read in one-body operator'
	  write(6,*)' (21) number operator N '
      write(6,*)' (99) no further interactions '
      write(6,*)' '

      write(6,*)' Enter choice from menu '
      read(5,*)ichoice
      if(ichoice.gt.21.or. ichoice.lt.0) then
      	if(describe) close (unit=7)
        return
      endif
      	
  
      write(6,*)' NOTE: For ALL interactions, V < 0 is attractive, ',
     &'V > 0 repulsive '
      write(6,*)' So be sure to include correct sign! '
      
C ---------- RANDOM NOISE ------------------------------------

      if(ichoice.eq.0)then 
         call randomnoise(describe)
      endif 
           
C----------- PAIRING -------------------------------------------

      if(ichoice.eq.1)then 
        call pairing(describe)
      endif      

C---------- MULTIPOLE --------------------------------------

      if(ichoice.eq.2)then
         call mpole(describe)
      endif

C--------- SIGMA*SIGMA ------------------------------------

      if(ichoice.eq.3)then 
        call spin2(describe)

      endif   

C--------- SIGMA*SIGMA ------------------------------------

      if(ichoice.eq.18)then 
        call gt(describe)

      endif 

C--------- L^2 ------------------------------------

      if(ichoice.eq.4)then 
        call l2(describe)

      endif               
C--------- J^2 ------------------------------------

      if(ichoice.eq.5)then 
        call j2(describe)

      endif      

C--------- L*S ------------------------------------

      if(ichoice.eq.6)then 
        call spinorbit(describe)

      endif               
      
C----------MONOPOLE -----------------------------------

      if(ichoice.eq.7)then
        call monopole(describe)
      endif
C---------------- CONSTANT DISPLACEMENT ------------------
      if(ichoice.eq.8)then
         print*,' Not currently working '

      endif

C------------------ GENERAL MULTIPOLE ------------------------


      if(ichoice.eq.9)then
            write(6,*)'  Random multipole '

  876 continue
            write(6,*)' Enter ang mom, Isospin rank of multipole '
C            read(5,*)L,T
C            if(t.ne.0 .and.t.ne.1)then
C              write(6,*)' Isospin must be 0 or 1 '
C              goto 876
C            endif
C            write(6,*)' Enter random seed '
C            read(5,*)iseed
C            if(iseed.gt.0)iseed=-iseed
C            write(6,*)' Enter strength of interaction ',
C     & ' ( < 0 attractive, > 0 repulsive '
C             read(5,*)v
C            call make_random_multipole(iseed,L,jt,sizej,nsp,
C     & sizeme,opme)

C      write(23,*)' random multipole '
C      write(23,*)nsp
C      do i = 1,nsp
C         write(23,*)nt(i),jt(i)
C      enddo
C      write(23,*)L,0
C      do i = 1,nsp
C        do k = 1,nsp
CC           write(23,*)k,i,opme(k,i)
C        enddo
C      enddo

C            call genMultipole(L,T,sizej, jt, nsp,
C     &  sizeme,opme,sizetbme,nme,jtot,ttot,map,mme,mspe)
C        do i = 1,nsp        
C          spe(i) = spe(i) +v*mspe(i)
C          mspe(i)=0.0	! for future use
C        enddo
C        do i = 1,nme
C          xme(i) = xme(i) +v*mme(i)
          
C          mme(i) = 0.0   ! for future use
C        enddo	

      endif

C--------- J^2 ------------------------------------

      if(ichoice.eq.10)then 
        call T2(describe)

      endif  

C.................... CONTACT

      if( ichoice .eq. 11)then
        call contact_master(describe)

      endif

C---------------- special interaction

      if(ichoice.eq.13)then




      endif
C----------------- Hcm for cold atoms-------------

      if(ichoice == 14 )then
         if(.not.spinless)then
            print*,' INVALID '   ! eventually put in Hcm for regular fermions
         else
            print*,' Calculate directly in lab frame (y/n)?'
            read(5,'(a)')ychar
            if(ychar=='n' .or. ychar=='N')then
            call Hcm_atom2
            else
            call Hcm_atom
            endif

         endif

      endif
C----------------- Hcm for nucleons-------------

      if(ichoice == 15 )then
         if(spinless)then
            print*,' INVALID '   ! eventually put in Hcm for regular fermions
         else
            print*,' Calculate directly in lab frame (y/n)?'
            read(5,'(a)')ychar
            if(ychar=='n' .or. ychar=='N')then
            call Hcm_nuke2
            else
            call Hcm_nuke
            endif

         endif

      endif
C----------------- R^2 intrinsic for cold atoms-------------

      if(ichoice == 16 )then
         if(.not.spinless)then
            print*,' INVALID '   ! eventually put in Hcm for regular fermions
         else
            print*,' Calculate directly in lab frame (y/n)?'
            read(5,'(a)')ychar
            if(ychar=='n' .or. ychar=='N')then
            call Rintr_atom

            else
            call Rrel_atom_lab
            endif
         endif

      endif
C---------------------------------------------------------
C----------------- Prel^2 for nucleons-------------

      if(ichoice == 17 )then
         if(spinless)then
            print*,' INVALID '   ! eventually put in Hcm for regular fermions
         else
            print*,' Calculate directly in lab frame (y/n)?'
            read(5,'(a)')ychar
            if(ychar=='n' .or. ychar=='N')then
            call Trel_nuke2
            else
            call Trel_nuke
            endif


         endif

      endif

!.............. "isoscalar pairing "---------
      

      if(ichoice == 19 )then
         if(spinless)then
            print*,' INVALID '   ! eventually put in Hcm for regular fermions
         else

            call scalarpairing


         endif

      endif
!....... READ IN ONE-BODY OPERATOR
      if(ichoice== 20)then
		  call readphonon
		  
	  end if	  
!....... Number operator n
      if(ichoice== 21)then
		  call number1
		  
	  end if	  


      goto 1
      
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine randomnoise(describe)

      use sporbit
      use tbmes
      implicit none
      real v
      real gauss
      integer iseed
      integer i,j,t
      logical describe
      integer ichoice
      print*,nme
!      write(6,*)' Note this version does not add noise to pairing'

      write(6,*)' Enter one of the following weightings '
      write(6,*)' '
      write(6,*)' (1) No weighting '
      write(6,*)' (2) Weight ~ sqrt(2J+1)(2T+1) '
      write(6,*)' (3) Weight ~ 1/sqrt(2J+1)(2T+1) '
      write(6,*)' (4) No J = 0, T = 1 matrix elements '
      write(6,*)' (5) All matrix elements < 0 '
      write(6,*)' '
      read(5,*)ichoice
      if(ichoice < 1 .or. ichoice > 5) then
         print*,' Not acceptable choice', ichoice
         return
      endif

      write(6,*)' Enter strength of random noise '
      read(5,*)V
      if(v.eq.0.0)return
        
        write(6,*)' Enter random seed '
        read(5,*) iseed
        iseed = -iseed
        print*,iseed
!        if(describe)write(7,*)' noise: ',V, ' (seed: ',iseed,')'
        print*,nme
        do i = 1,nme
!            if(vtbme(i)%t==0)then
            j = vtbme(i)%j
            t = vtbme(i)%t

            select case (ichoice)

            case(1)
            vtbme(i)%v = vtbme(i)%v+v*gauss(iseed)

            case(2)
            vtbme(i)%v = vtbme(i)%v+ v*gauss(iseed)*
     &  sqrt(2.*j+1.)*sqrt(2.*t+1.)

            case(3)
            vtbme(i)%v = vtbme(i)%v+v*gauss(iseed)/
     & sqrt(2.*j+1.)*sqrt(2.*t+1.)
            case(4)
            if(j/= 0 .or. t /=1)then
            vtbme(i)%v = vtbme(i)%v+v*gauss(iseed)
            endif 

            case(5)
            vtbme(i)%v = vtbme(i)%v- v*abs(gauss(iseed))

            end select

!            endif
        enddo

      return
      end subroutine randomnoise
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine pairing(describe)

      use sporbit
      use tbmes
      implicit none
      real v
      integer i
      real xj1,xj3
      logical describe

      	write(6,*)' Enter pairing strength G'
      	read(5,*)V
        if(v.eq.0.0)return
      	if(describe)write(7,*)' pairing: ',V 
        do i = 1,nme
          if( vtbme(i)%j == 0
     &   .and. vtbme(i)%n(1) == vtbme(i)%n(2)
     &   .and.  vtbme(i)%n(3) == vtbme(i)%n(4) )then
             xj1 = float( orbqn( vtbme(i)%n(1) )%j )/2.
             xj3 = float( orbqn( vtbme(i)%n(3) )%j )/2.
 
             vtbme(i)%v=vtbme(i)%v+v*sqrt((xj1+0.5)*(xj3+0.5))
          endif
        
        enddo

      return
      end subroutine pairing
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine scalarpairing

      use sporbit
      use tbmes
      implicit none
      real v
      integer i
      real xj1,xj3
      integer a, b, c, d,ja,jb,jc,jd,li,lf
      logical describe
      real xme
      real sj2i

        print*,' Here couple up pairs to L = 0, S = 1, T = 0 '
      	write(6,*)' Enter pairing strength G (make < 0 for attractive)'
      	read(5,*)V
        if(v.eq.0.0)return
      	if(describe)write(7,*)' isoscalar pairing: ',V 
        do i = 1,nme
          if(vtbme(i)%j/=1 .or. vtbme(i)%t/=0)cycle

          a = vtbme(i)%n(1)
          b = vtbme(i)%n(2)
          c = vtbme(i)%n(3)
          d = vtbme(i)%n(4)
          if( orbqn(a)%nr /= orbqn(b)%nr )cycle
          if(  orbqn(c)%nr /= orbqn(d)%nr)cycle

          if( orbqn(a)%l /= orbqn(b)%l)cycle
          if( orbqn(c)%l /= orbqn(d)%l)cycle

          ja = orbqn(a)%j
          jb = orbqn(b)%j
          jc = orbqn(c)%j
          jd = orbqn(d)%j

          lf = orbqn(a)%l
          li = orbqn(c)%l

          xme = (-1)**(lf+li+1+(ja+jc)/2)*2.
          xme = xme*sqrt( float( (ja+1)*(jb+1)*(jc+1)*(jd+1) ))
          xme = xme* sj2i(jb,ja,2, 1, 1, 2*lf)*sj2i(jd,jc,2,1,1,2*li)
          if(ja ==jb)xme =xme/sqrt(2.)
          if(jc==jd)xme =xme/sqrt(2.)
             vtbme(i)%v=vtbme(i)%v+v*xme
        
        enddo

      return
      end subroutine scalarpairing
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mpole(describe)

      use sporbit
      use tbmes
      implicit none
      real v
      real gauss
      integer iseed
      integer i
      integer L
      real b
      character ychar*1
      logical phase_convention
      logical describe
      integer J,T			!  for 2-particle states
      integer i1,i2,i3,i4
      integer n1p,n2p,n1,n2		! radial q# of s.p. states
      integer l1p,l2p,l1,l2 		! orbital ang mom
      real xj1p,xj2p,xj1,xj2 		! j's 

      real vx				! 2-body matrix element 


      	write(6,*)' Multipole interaction is of the form ',
     &' V(r1,r2)= 1/2 r1^L r2^L sum_M (-1)^M  Y_LM(r1) Y_L-M (r2) ',
     &' x (2L+1)/4 pi '
        write(6,*)' (using h.o. wfns) '
        write(6,*)' Enter multipolarity L '
        read(5,*)L
        if(L.gt.8)then
          write(6,*)' Multipolarity too high '
          return
        endif
C        write(6,*)' Choice of "isoscalar" (T=0) P+N densities ',
C     &'or "isovector" (T=1) P-N densities '
    2   write(6,*)' Enter T=0,1'
        read(5,*)T
        T = 0 			! for now only allow isoscalar densities 
        if(t.ne.0 .and. t.ne.1)then
            write(6,*)' T must be 0 or 1 '
            goto 2
        endif
C        if(t.eq.1)then
C            write(6,*)' WARNING so-called isovector not consistently '
C            write(6,*)' defined, so use as your own risk. '
C        endif

        write(6,*)' You must choose a phase convention of ordering '
        write(6,*)' (l,1/2) = j or (1/2,l)= j; former is the default .'
        write(6,*)' Do you wish default phase? (y/n) '
        read(5,'(a)')ychar
        if(ychar.eq.'n' .or.ychar.eq.'N')then
           phase_convention=.true.
           write(6,*)' Choosing alternate phase convention '
        else
           phase_convention=.false.
        endif

        write(6,*)' Enter oscillator length parameter b, in fm '
        write(6,*)' (b = hbar c /sqrt( m c^2 hbar w), where ',
     &  'typical hbar w for nuclei = 41 A^1/3 MeV; '
        write(6,*)' Then numerically b ~ 1.0 A^-1/6 fm ) '
        if(L.eq.2)then 
          write(6,*)' Note: to get SU(3) Casimir_2, set b = 2^(3/4) ',
     & '= 1.681792831 and set strength = 2 '
        
        endif
        read(5,*)b 
        if(b.eq.0.0)return
        write(6,*)' Enter strength '
        read(5,*)V
        if(v.eq.0.0)return        
        if(describe)write(7,*)' MM: L = ',L,', T = ',T,',b = ',b,
     &  ' fm; ',' V: ',V
C        write(6,*)active(1)
C        call Multipole(L,T,nsp,nt,lt,jt,nme,jtot,ttot,map,mme,mspe,
C     &  active,phase_convention)

C         call  make_std_multipole(L,jt,lt,nt,sizej,nsp,
C     & sizeme,opme)

C............. write to file............

      write(23,*)' multipole from make tbme '

C      write(23,*)nsp
C      do i = 1,nsp
C         write(23,*)nt(i),jt(i)
C      enddo
C      write(23,*)L,0
C      do i = 1,nsp
C        do k = 1,nsp
C           write(23,*)k,i,opme(k,i)*b**L
C        enddo
C      enddo


C            call genMultipole(L,T,sizej, jt, nsp,
C     &  sizeme,opme,sizetbme,nme,jtot,ttot,map,mme,mspe)
        
C        do i = 1,nsp        
C          spe(i) = spe(i) +0.5* V*b**(2*L)*mspe(i)
C          mspe(i) = 0.0  ! for future use
C        enddo
        do i = 1,nme
      J = vtbme(i)%j
      T = vtbme(i)%t
      i1 = vtbme(i)%n(1)
      i2 = vtbme(i)%n(2)
      i3 = vtbme(i)%n(3)
      i4 = vtbme(i)%n(4)
      n1 = orbqn(i1)%nr
      n2 = orbqn(i2)%nr
      n1p = orbqn(i3)%nr
      n2p= orbqn(i4)%nr

      xj1 = float(orbqn(i1)%j)/2.
      xj2 = float(orbqn(i2)%j)/2.
      xj1p = float(orbqn(i3)%j)/2.
      xj2p = float(orbqn(i4)%j)/2.

      l1 = orbqn(i1)%l
      l2 = orbqn(i2)%l
      l1p = orbqn(i3)%l
      l2p = orbqn(i4)%l
            vx = 0.0
            call multipoleME(L,J,n1p,l1p,xj1p,n2p,l2p,xj2p,
     & n1,l1,xj1,n2,l2,xj2,t,vx) 

C 
          vtbme(i)%v = vtbme(i)%v+ 0.5*V*b**(2*L)*vx
          
        enddo


      return
      end subroutine mpole

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine spin2(describe)
      
      use sporbit
      use tbmes
      implicit none
      integer i
      real v
      logical describe

      	write(6,*)' Enter strength for S^2 '
      	read(5,*)V
        if(v == 0.0)return
      	if(describe)write(7,*)' S^2: ',V 
      	do i = 1,numorb
      	    spe(i) = spe(i)+.75*v
      	enddo

        do i = 1,nme
	    call ss(i,v)
     
        enddo

      return
      end subroutine spin2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gt(describe)
      
      use sporbit
      use tbmes
      implicit none
      integer i
      real v
      logical describe

      	write(6,*)' Enter strength for S^2 T^2 '
      	read(5,*)V
        if(v == 0.0)return
      	if(describe)write(7,*)' S^2 T^2 : ',V 
      	do i = 1,numorb
      	    spe(i) = spe(i)+0.75**2*v
      	enddo

        do i = 1,nme
	    call sstt(i,v)
     
        enddo

      return
      end subroutine gt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine l2(describe)
      
      use sporbit
      use tbmes
      implicit none
      integer i
      real v
      logical describe
      integer l
 
      	write(6,*)' Enter strength for L^2 '
      	read(5,*)V
        if(v.eq.0.0)return
      	if(describe)write(7,*)' L^2: ',V 
      	do i = 1,numorb
            l  = orbqn(i)%l
      	    spe(i) = spe(i)+v*float(l*(l+1))
      	enddo
        do i = 1,nme
	    call ll(i,V)
        
        enddo

      return
      end subroutine l2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine j2(describe)
      
      use sporbit
      use tbmes
      implicit none
      integer i
      real v
      logical describe
      real xj1

      	write(6,*)' Enter strength for J^2 '
      	read(5,*)V
        if(v.eq.0.0)return
      	if(describe)write(7,*)' J^2: ',V 
      	do i = 1,numorb
	    xj1 = float(orbqn(i)%j)/2.
      	    spe(i) = spe(i)+v*(xj1*(xj1+1.))
      	enddo
        do i = 1,nme
	    call jsq(i,V)
C	    xme(i) =xme(i)+2*v*temp
        
        enddo 

      return
      end subroutine j2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine t2(describe)
      
      use sporbit
      use tbmes
      implicit none
      integer i
      real v
      logical describe
      real xj1

      	write(6,*)' Enter strength for T^2 '
      	read(5,*)V
        if(v.eq.0.0)return
      	if(describe)write(7,*)' T^2: ',V 
      	do i = 1,numorb
	    xj1 = float(orbqn(i)%j)/2.
      	    spe(i) = spe(i)+v*0.75
      	enddo
        do i = 1,nme
	    call tsq(i,V)
C	    xme(i) =xme(i)+2*v*temp
        
        enddo 

      return
      end subroutine t2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine spinorbit(describe)

      use sporbit
      use tbmes
      implicit none
      integer i
      real v
      real xl1,xj1
      logical describe
	  integer :: it,tfactor

	    write(6,*)' Enter 0 or 1 for isospin dependence '
		read(5,*)it
		
      	write(6,*)' Enter strength for L^2 '
      	write(6,*)' (Note, this is total L * total S) '
      	read(5,*)V
        if(v.eq.0.0)return
      	if(describe)write(7,*)' L*S: ',V 
      	do i = 1,numorb
	    xl1 = float( orbqn(i)%l )     	   
	    xj1 = float( orbqn(i)%j )/2.   	   
      	    spe(i) = spe(i)-0.5*v*(xl1*(xl1+1.)+0.75-xj1*(xj1+1))
      	enddo
        do i = 1,nme
			if(it==1 .and. vtbme(i)%t==1)call ls(i,2*V)
	        if(it==0) call ls(i,V)
		  
        enddo

      return
      end subroutine spinorbit

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine monopole(describe)


      use sporbit
      use tbmes
      implicit none
      logical describe

      integer i,k,n
      integer a,b,c,d
      real v

      	write(6,*)' This option is for a force of the type ',
     & '  W_ab  n_a * (n_b-delta_ab), where n_a = number in shell a '
       

        write(6,*)' Enter W_ab for the following a,b single '
     & ,'particle states: '
        write(6,*)'N(a) J(a) ; N(b) J(b) '
        do i= 1,numorb
           do k=i,numorb
               write(6,*)orbqn(i)%j,'/2  ; ',
     &                   orbqn(k)%j,'/2   ?'
               read(5,*)v
C               if(describe .and. v.ne.0.0)then
C               	 write(7,*)' MONOPOLE (',nt(i),', ',
C     &  jt(i),'/2; ',nt(k),', ',jt(k),'/2): ',v
C               endif
               
	if(v.ne.0.0)then

      do n = 1,nme 
        a = vtbme(n)%n(1)
        b = vtbme(n)%n(2)
        c = vtbme(n)%n(3)
        d = vtbme(n)%n(4)
        if(a /=c)cycle
        if(b /=d)cycle
        if( ( i == a .and. k == b) .or. (i ==b .and. k ==a))then
           vtbme(n)%v = v
C        if((i.eq. map(n,1) .and. i.eq.map(n,3) 
C     &  .and. k .eq. map(n,2) .and. k .eq. map(n,4))
C     & .or.
C     & (i.eq. map(n,2) .and. i.eq.map(n,4) 
C     &  .and. k .eq. map(n,1) .and. k .eq. map(n,3)) )then
C         xme(n) = V
         endif
        
        enddo	! loop over n
        
C          if(i.eq.k)spe(i) = spe(i)+0.5*V
          endif
           
           
           enddo ! loop over k	
        enddo	!  loop over i
        
      
      return
      end subroutine monopole

	  
C==========================================================================
	  subroutine readphonon
	      use sporbit
	      use tbmes
          implicit none
		  
		  character*25 obmefile  
		  integer  ilast
		  integer iorb
		  character*4 title
		  integer i,n,l
		  real    xj,x
		  integer jtran,ttran
		  integer a,b
		  logical success
		  real xscale
		  
		  print*,' '
		  print*,' Need to enter a 1-body operator generated by ',
     & ' program tropic1b.x '
		  print*,' '
1         continue
		  print*,' Enter name of .opme file (do not enter suffix)'
		  read(5,'(a)')obmefile
		  ilast = index(obmefile,' ')-1
		  open(unit=22,file=obmefile(1:ilast)//'.opme',status='old',
     & err=3)
		  goto 4
3         continue
          print*,' File ',obmefile(1:ilast),'.opme does not exist '
		  stop
4         continue
	      
!........ READ BEGINNING......................
          success = .false.
          do i = 1,1000
			  read(22,'(a)')title
			  if(title(1:4)==' iso')then
				  success=.true.
				  exit
			  end if
		  end do ! i
		  if(.not.success)then
			  print*,' File does not appear correctly formatted'
			  stop
		  end if
		  read(22,*)n
		  if(n/=numorb)then
		      print*,' Mismatch in # of orbits ',n,numorb
			  stop	  
	  	  end if	        		  
!........ FIRST CHECK SINGLE-PARTICLE ORBITS AGREE

          do iorb = 1,numorb
              read(22,*)i,n,l,xj
			  if(i/=iorb)then
				  print*,' unordered s.p. orbits ',i,iorb
				  stop
			  end if
			  if(n/=orbqn(iorb)%nr)then
				  print*,' bad value of n ',iorb,n,orbqn(iorb)%nr
				  stop
			  end if
			  if(l/=orbqn(iorb)%l)then
				  print*,' bad value of l ',iorb,l,orbqn(iorb)%l
				  stop
			  end if
			  if(int(2*xj)/=orbqn(iorb)%j)then
				  print*,' bad value of j ',iorb,xj,orbqn(iorb)%j
				  stop
			  end if
          end do  ! iorb
!........READ IN J,T OF TRANSITION.....
          read(22,*)jtran,ttran
		  print*,' '
		  print*,' J; T of transition = ',jtran,ttran

!........READ IN one-body matrix elements........
          allocate(obme(numorb,numorb))
          obme(:,:)=0.0
		  do i = 1,numorb**2
			  read(22,*,end=45)a,b,x
			  obme(a,b)=x
		  end do
45        continue
          close(22)	  
!............ READ IN OVERALL SCALE....
          print*,' '
		  print*,' Enter overall scale '
		  read*,xscale
!........... COMPUTE SINGLE PARTICLE ENERGIES.... assuming only diagonal
          do a = 1,numorb
			  call spephonon(a,xscale,spe(a))
          end do

!........... COMPUTE 2-body matrix elements...........
          do i = 1,nme
			  call v2phonon(i,jtran,ttran,xscale)
		  end do
		  return
      
	  end subroutine readphonon
!===============================================================
!
!  some multiple of the number operator N
!

      subroutine number1
	      use sporbit
	      use tbmes
          implicit none
		  real xscale
		  integer :: i
!............ READ IN OVERALL SCALE....
		  print*,' '
		  print*,' Enter overall scale '
		  read*,xscale
		  do i = 1,numorb
			  spe(i)=spe(i)+xscale
			  
		  end do ! i
		  		  

		  return
	  end subroutine number1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION RAN3(IDUM)
C         IMPLICIT REAL*4(M)
C         PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=2.5E-7)
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      common/randumb/inext,inextp,ma     ! needed on HP
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE
        DO 13 K=1,4
          DO 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
        INEXT=0
        INEXTP=31
        IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      RETURN
      END

      real function gauss(iseed)
C************************************************************
c      returns Gaussian random variable.                      *
C      Algorithm is that on p. 209-210 Computational Physics  *
C      by Koonin and Meridith (fortran version).              *
c      Note that the algorithm generates two variables,       *
C      here we randomly select between the two.               *
c**************************************************************
      implicit real (a-h,o-z)
      pi=2.00*asin(1.0)
      two=-2.0*log(1.0-ran3(iseed))
      radius=sqrt(two)
      theta=2.0*pi*ran3(iseed)
      gauss1=radius*cos(theta)
      gauss2=radius*sin(theta)
      xchoose=ran3(iseed)
      if(xchoose.le.0.50)then
         gauss=gauss1
      else
         gauss=gauss2
      end if
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

