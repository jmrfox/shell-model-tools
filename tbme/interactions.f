

C==========================================================================

      subroutine multipoleME(L,J,n1p,l1p,xj1p,n2p,l2p,xj2p,
     & n1,l1,xj1,n2,l2,xj2,t,tbme) 
C
C  computes matrix element between h.o. wfns 
C  < n1p (1/2 l1p) j1p, n2p (1/2 l2p) jp2; JT ||| [O(L,iso) x O(L,iso)]_00 
C   ||| n1 (1/2 l1) j1, n2 (1/2 l2) j2; JT >
C    
C  where O(L,iso) is a multipole operator of the form r^L Y_Lm   
C
C  NOTE -- to ease consistency with SU(3) and other operators,
C  multiplied by (4pi/2L+1)

      implicit none
   
      integer L 			! ang mom of multipole 
      integer J,T			!  for 2-particle states

      integer n1p,n2p,n1,n2		! radial q# of s.p. states
      integer l1p,l2p,l1,l2 		! orbital ang mom
      real xj1p,xj2p,xj1,xj2 		! j's 

      real tbme				! 2-body matrix element 
      real sj,tj,HORadInt
      real v1,v2

      tbme =  1./float(2*L+1) !   (4.*3.1415926) 
     &  *sqrt((2*xj1+1)*(2*xj2+1)*(2*xj1p+1)*(2*xj2p+1))
     &  *sqrt((2.*l1p+1.)*(2.*l2p+1.)*(2.*l1+1)*(2*l2+1))
     &  *sj(float(l1p),xj1p,0.5,xj1,float(l1),float(L))
     &  *sj(float(l2p),xj2p,0.5,xj2,float(l2),float(L))
     &  *tj(float(l1p),float(L),float(l1),0.,0.0,0.)
     &  *tj(float(l2p),float(L),float(l2),0.,0.0,0.)
      v1 = HORadInt(L,n1,l1,n1p,l1p)
      v2 = HORadInt(L,n2,l2,n2p,l2p)
C      write(6,*)' to date ',tbme,v1,v2
      tbme = tbme*v1*v2
     & *(-1)**(1+int(xj1+xj2))
      if(t==1)tbme = tbme * *0.75*(-1)**(1+T)/float(2*T+1)

C      write(6,*)' tbme ',tbme
      return
      end

CCCCCCCCCCCCCCCCALVINCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real function HORadInt(L,n1,l1,n2,l2)
C
C   computes radial integral R(n1,l1, alpha*r) R(n2,l2, alpha*r) r^(L+2) 
C   where R's are h.o. radial functions 
C
C   uses Lawson 1.11a; alpha = 1  must scale by 1/alpha^2 to get final 
C   answer  (alpha = 1/b, b = h.o. scale parameter)
C
      implicit none 
    
      integer L			! weighting of r ( = L in Ylm)
      integer n1,n2,l1,l2       ! q#'s of states 

      real sum			! result 
      integer q			! dummy for summation 
      integer qmin,qmax		! limits of summation 

      real lnprefact,lnsum,lnfact,ln2fact		

C............let's get going..................

C      write(6,*)' in HOradint ',l1,l2,l
      if( mod(l1+l2+L+2,2) .ne. 0)then 	  ! must be even overall
	HORadInt = 0.0    
	return
      endif

C............. l1,l2,L must satisfy triangle relation....
      if( ( l1+l2 .lt. L) .or. (abs(l1-l2) .gt. L) )then
	HORadInt = 0.0    
	return
      endif

      lnprefact = (lnfact(n1)+lnfact(n2) + log(2.)*(n1+n2-L) 
     & - ln2fact(2*n1+2*l1+1) - ln2fact(2*n2+2*l2+1))/2. 
     & + lnfact( (l2-l1+L)/2) + lnfact( (l1-l2+L)/2)
      
      qmax = min(n1,n2)
      qmin = max(0,max( n1-(l2-l1+L)/2,n2-(l1-l2+L)/2))
      sum = 0.0
C      write(6,*)' prefact ',lnprefact,qmin,qmax
      lnsum = 0.0
      do q = qmin,qmax
	lnsum =  ln2fact(l1+l2+L+2*q+1) 
     & -q*log(2.) - lnfact(q) - lnfact(n1-q) -lnfact(n2-q) 
     & - lnfact(q+(l1-l2+L)/2-n2) - lnfact(q+(l2-l1+L)/2-n1)
C        write(6,*)q,lnsum
	sum = sum + exp( lnprefact+lnsum) 
      enddo

      if( mod( abs(n1-n2),2) .ne.0) sum = -sum
     
      HORadInt =   sum
C      write(6,*)' radial integral ',sum
      return
      end

C==========================================================================
      real function lnfact(n)

      implicit none 

      integer n
      integer i

      lnfact = 0.0

      do i = 2,n
	lnfact = lnfact + log(float(i))
      enddo
      return
      end

C==========================================================================
      real function ln2fact(n)
C
C   log of double factorial ln n!!
C
      implicit none 

      integer n
      integer i

      ln2fact = 0.0
      if(n .eq. 0 .or. n .eq. 1) return 

      do i = n,1,-2
	ln2fact = ln2fact + log(float(i))
      enddo
      return
      end

C==========================================================================

      integer function delta(i,j)
      implicit none
      integer i,j
      
      if(i.eq.j)then
      	delta = 1
      else
      	delta = 0
      endif
      return
      end
C==========================================================================
      real function zeta(i,j)
      
      implicit none
      integer i,j,delta
      zeta=sqrt(1.0+float(delta(i,j)) )
      return
      end
C==========================================================================

      subroutine ss(i,V)
C
C  computes antisymmeterized matrix element of sigma(1)*sigma(2)
C  IMPORTANT- NOTE: assumes  coupling  l x s = j not s x l 
C
C  INPUT:
C       I:    which matrix element in list is it
C
C  OUTPUT:
C	Vss:  matrix element
C
C  CALLS:
C	delta(i,j)	delta function
C	sixj
C  
      use sporbit
      use tbmes
      implicit none

      integer i
      real v,vss
      integer J,T

C---- SINGLE-PARTICLE LIST ------
      integer i1,i2,i3,i4
      real xj1,xj2,xj3,xj4
      real xl1,xl2,xl3,xl4  
      integer l1,l2,l3,l4
      integer n1,n2,n3,n4     
      
C---------functions called

      integer delta
      real sj  ! found in libra
 
      J = vtbme(i)%j
      T = vtbme(i)%t
      i1 = vtbme(i)%n(1)
      i2 = vtbme(i)%n(2)
      i3 = vtbme(i)%n(3)
      i4 = vtbme(i)%n(4)

      n1 = orbqn(i1)%nr
      n2 = orbqn(i2)%nr
      n3 = orbqn(i3)%nr
      n4 = orbqn(i4)%nr

      xj1 = float(orbqn(i1)%j)/2.
      xj2 = float(orbqn(i2)%j)/2.
      xj3 = float(orbqn(i3)%j)/2.
      xj4 = float(orbqn(i4)%j)/2.

      l1 = orbqn(i1)%l
      l2 = orbqn(i2)%l
      l3 = orbqn(i3)%l
      l4 = orbqn(i4)%l

      xl1 = float(orbqn(i1)%l)
      xl2 = float(orbqn(i2)%l)
      xl3 = float(orbqn(i3)%l)
      xl4 = float(orbqn(i4)%l)

C -- direct piece first
      vss =   ( delta(n1,n3)* delta(n2,n4) * delta(l1,l3)*delta(l2,l4)
     &  * (-1)**( int(l1+l2+ xj1+xj3)+2 )*
     &   sj(0.5, xj1,xl1,xj3,.5,1.)*sj(.5,xj2,xl2,xj4,.5,1.)
     &   *sj(float(J),xj2,xj1,1.0,xj3,xj4)
C      write(6,*)vss
C      vss = vss
     &   -(-1)**(2+ int( J+T +xj3-xj4))*
C---indirect piece      
     &   delta(n1,n4)*delta(n2,n3)*delta(l1,l4)*delta(l2,l3)
     &   *(-1)**( int(l1+l2+ xj1+xj4) +2)*
     &   sj(float(J),xj2,xj1,1.0,xj4,xj3)*
     &   sj(0.5, xj1,xl1,xj4,.5,1.)*sj(.5,xj2,xl2,xj3,.5,1.) )
       
       vss = vss *3. /sqrt( float ( (1+ delta(i1,i2 ) )
     & * (1+delta(i3,i4) ) ))
     &  *sqrt( (2*xj1+1)* (2*xj2+1)* (2*xj3+1)* (2*xj4+1))  * 
     & (-1)**( J)
C      write(6,*)vss    ,sixj(float(jtot(i))/2.,xj2,xj1,1.0,xj4,xj3)
      vtbme(i)%v = vtbme(i)%v + v*vss
      return
      end          
C==========================================================================

      subroutine sstt(i,V)
C
C  computes antisymmeterized matrix element of sigma(1)*sigma(2)
C  IMPORTANT- NOTE: assumes  coupling  l x s = j not s x l 
C
C  INPUT:
C       I:    which matrix element in list is it
C
C  OUTPUT:
C	Vss:  matrix element
C
C  CALLS:
C	delta(i,j)	delta function
C	sixj
C  
      use sporbit
      use tbmes
      implicit none

      integer i
      real v,vsstt
      integer J,T

C---- SINGLE-PARTICLE LIST ------
      integer i1,i2,i3,i4
      real xj1,xj2,xj3,xj4
      real xl1,xl2,xl3,xl4  
      integer l1,l2,l3,l4
      integer n1,n2,n3,n4     
      
C---------functions called

      integer delta
      real sj  ! found in libra
 
      J = vtbme(i)%j
      T = vtbme(i)%t
      i1 = vtbme(i)%n(1)
      i2 = vtbme(i)%n(2)
      i3 = vtbme(i)%n(3)
      i4 = vtbme(i)%n(4)

      n1 = orbqn(i1)%nr
      n2 = orbqn(i2)%nr
      n3 = orbqn(i3)%nr
      n4 = orbqn(i4)%nr

      xj1 = float(orbqn(i1)%j)/2.
      xj2 = float(orbqn(i2)%j)/2.
      xj3 = float(orbqn(i3)%j)/2.
      xj4 = float(orbqn(i4)%j)/2.

      l1 = orbqn(i1)%l
      l2 = orbqn(i2)%l
      l3 = orbqn(i3)%l
      l4 = orbqn(i4)%l

      xl1 = float(orbqn(i1)%l)
      xl2 = float(orbqn(i2)%l)
      xl3 = float(orbqn(i3)%l)
      xl4 = float(orbqn(i4)%l)

C -- direct piece first
      vsstt =   ( delta(n1,n3)* delta(n2,n4) * delta(l1,l3)*delta(l2,l4)
     &  * (-1)**( int(l1+l2+ xj1+xj3)+2 )*
     &   sj(0.5, xj1,xl1,xj3,.5,1.)*sj(.5,xj2,xl2,xj4,.5,1.)
     &   *sj(float(J),xj2,xj1,1.0,xj3,xj4)
C      write(6,*)vss
C      vss = vss
     &   -(-1)**(2+ int( J+T +xj3-xj4))*
C---indirect piece      
     &   delta(n1,n4)*delta(n2,n3)*delta(l1,l4)*delta(l2,l3)
     &   *(-1)**( int(l1+l2+ xj1+xj4) +2)*
     &   sj(float(J),xj2,xj1,1.0,xj4,xj3)*
     &   sj(0.5, xj1,xl1,xj4,.5,1.)*sj(.5,xj2,xl2,xj3,.5,1.) )
       
       vsstt = vsstt *3. /sqrt( float ( (1+ delta(i1,i2 ) )
     & * (1+delta(i3,i4) ) ))
     &  *sqrt( (2*xj1+1)* (2*xj2+1)* (2*xj3+1)* (2*xj4+1))  * 
     & (-1)**( J)
C      write(6,*)vss    ,sixj(float(jtot(i))/2.,xj2,xj1,1.0,xj4,xj3)
      vsstt = vsstt*0.75*(-1)**(1+T)/float(2*T+1)
      vtbme(i)%v = vtbme(i)%v + v*vsstt
      return
      end          

C==========================================================================

      subroutine ll(i,V)
C
C  computes antisymmeterized matrix element of L(1)*L(2)
C  IMPORTANT- NOTE: assumes  coupling  l x s = j not s x l 
C  INPUT:
C       I:    which matrix element in list is it
C	NT(i):  PRINCIPLE radial Q# (not NODAL #) for i'th s.p. orbit
C 			N = 2n+l, where n is nodal q#
C	LT(i):  orbital ang mom for i'th s.p. orbit
C	JT(i):  2xJ for i'th s.p. orbit
C	JTOT(k):  total J for k'th tbme
C	TTOT(k):  total T for k'th TBME
C	MAP(k,i): k'th matrix element < 1 2 | V| 3 4 > i = 1,4 
C				  maps s.p. states
C
C  OUTPUT:
C	V:  matrix element
C
C  CALLS:
C	delta(i,j)	delta function
C	sixj
C  

      use sporbit
      use tbmes
      implicit none

      integer i
      real v,vll
      integer J,T

C---- SINGLE-PARTICLE LIST ------
      integer i1,i2,i3,i4
      real xj1,xj2,xj3,xj4
      real xl1,xl2,xl3,xl4  
      integer l1,l2,l3,l4
      integer n1,n2,n3,n4     
     
C---------functions called

      integer delta
      real sj  ! found in libra
      real zeta

      J = vtbme(i)%j
      T = vtbme(i)%t
      i1 = vtbme(i)%n(1)
      i2 = vtbme(i)%n(2)
      i3 = vtbme(i)%n(3)
      i4 = vtbme(i)%n(4)
      n1 = orbqn(i1)%nr
      n2 = orbqn(i2)%nr
      n3 = orbqn(i3)%nr
      n4 = orbqn(i4)%nr

      xj1 = float(orbqn(i1)%j)/2.
      xj2 = float(orbqn(i2)%j)/2.
      xj3 = float(orbqn(i3)%j)/2.
      xj4 = float(orbqn(i4)%j)/2.

      l1 = orbqn(i1)%l
      l2 = orbqn(i2)%l
      l3 = orbqn(i3)%l
      l4 = orbqn(i4)%l

      xl1 = float(orbqn(i1)%l)
      xl2 = float(orbqn(i2)%l)
      xl3 = float(orbqn(i3)%l)
      xl4 = float(orbqn(i4)%l)

C -- direct piece first
      vll =  (delta(n1,n3)*delta(n2,n4)*delta(l1,l3)*delta(l2,l4)
     &   *sj(float(j),xj2,xj1,1.0,xj3,xj4)
     &   *(-1)**( int( xj2+xj4)+2 )*
     &   sj(xl3, xj1,0.5,xj3,xl1,1.)*sj(xl2,xj2,0.5,xj4,xl2,1.)
C
     &   -(-1)**(2+ int(j+t+ int(xj3-xj4)))*
C---indirect piece      
     &   delta(n1,n4)*delta(n2,n3)*delta(l1,l4)*delta(l2,l3)*
     &   (-1)**( int( xj2+xj3) +2)*
     &   sj(float(j),xj2,xj1,1.0,xj4,xj3)*
     &   sj(xl3, xj1,0.5,xj4,xl1,1.)*sj(xl2,xj2,0.5,xj3,xl2,1.) )
       vll = vll /zeta(i1,i2) /zeta(i3,i4 ) 
     &  *sqrt( (2*xj1+1)* (2*xj2+1)* (2*xj3+1)* (2*xj4+1))  
     &  * sqrt( xl1*xl2*(xl1+1.)*(xl2+1)*(2*xl1+1)*(2*xl2+1) )
     &  *(-1)**( j+2 +l1+l2  )
     
      vtbme(i)%v = vtbme(i)%v + v*vll*2
      return
      end          
C==========================================================================

      subroutine ls(i,V)
C
C  computes antisymmeterized matrix element of L(1)*S(2)+L(2)*S(1)
C  IMPORTANT- NOTE: assumes  coupling  l x s = j not s x l 
C  INPUT:
C       I:    which matrix element in list is it
C	NT(i):  PRINCIPLE radial Q# (not NODAL #) for i'th s.p. orbit
C 			N = 2n+l, where n is nodal q#
C	LT(i):  orbital ang mom for i'th s.p. orbit
C	JT(i):  2xJ for i'th s.p. orbit
C	JTOT(k):  total J for k'th tbme
C	TTOT(k):  total T for k'th TBME
C	MAP(k,i): k'th matrix element < 1 2 | V| 3 4 > i = 1,4 
C				  maps s.p. states
C
C  OUTPUT:
C	V:  matrix element
C
C  CALLS:
C	delta(i,j)	delta function
C	sixj
C  

      use sporbit
      use tbmes
      implicit none

      integer i
      real v,vls
      integer J,T

C---- SINGLE-PARTICLE LIST ------
      integer i1,i2,i3,i4
      real xj1,xj2,xj3,xj4
      real xl1,xl2,xl3,xl4  
      integer l1,l2,l3,l4
      integer n1,n2,n3,n4     
     
C---------functions called

      integer delta
      real sj  ! found in libra
      real zeta

      J = vtbme(i)%j
      T = vtbme(i)%t
      i1 = vtbme(i)%n(1)
      i2 = vtbme(i)%n(2)
      i3 = vtbme(i)%n(3)
      i4 = vtbme(i)%n(4)
      n1 = orbqn(i1)%nr
      n2 = orbqn(i2)%nr
      n3 = orbqn(i3)%nr
      n4 = orbqn(i4)%nr

      xj1 = float(orbqn(i1)%j)/2.
      xj2 = float(orbqn(i2)%j)/2.
      xj3 = float(orbqn(i3)%j)/2.
      xj4 = float(orbqn(i4)%j)/2.

      l1 = orbqn(i1)%l
      l2 = orbqn(i2)%l
      l3 = orbqn(i3)%l
      l4 = orbqn(i4)%l

      xl1 = float(orbqn(i1)%l)
      xl2 = float(orbqn(i2)%l)
      xl3 = float(orbqn(i3)%l)
      xl4 = float(orbqn(i4)%l)
 
C -- direct piece first
      vls =   ( delta(n1,n3)*delta(n2,n4)*delta(l1,l3)*delta(l2,l4)
     &  * sj(float(j),xj2,xj1,1.0,xj3,xj4)
     & *(  sqrt( (2*xl1+1)*xl1*(xl1+1) )
     & *  sj(xl3, xj1,0.5,xj3,xl1,1.)*sj(0.5,xj2,xl2,xj4,0.5,1.)
     & + (-1)**(int( xj1+xj2+xj3+xj4+2))*sqrt((2*xl2+1)*xl2*(xl2+1)) 
     & *  sj(xl2, xj2,0.5,xj4,xl2,1.)*sj(0.5,xj1,xl1,xj3,0.5,1.)
     & ) 
     &   -(-1)**(2+ int(j+t+xj3-xj4))
C---indirect piece      
     &   *delta(n1,n4)*delta(n2,n3)*delta(l1,l4)*delta(l2,l3)
     &    *sj(float(j),xj2,xj1,1.0,xj4,xj3)
     & *(  sqrt( (2*xl1+1)*xl1*(xl1+1) )
     & *  sj(xl1, xj1,0.5,xj4,xl1,1.)*sj(0.5,xj2,xl2,xj3,0.5,1.)
     & + (-1)**(int( xj1+xj2+xj3+xj4+2))*sqrt((2*xl2+1)*xl2*(xl2+1)) 
     & *  sj(xl2, xj2,0.5,xj3,xl2,1.)*sj(0.5,xj1,xl1,xj4,0.5,1.)
     & )
     & )

      vtbme(i)%v =vtbme(i)%v+ v*vls*(-1)**(j +int(xl1+xl2)+1)*sqrt(1.5)  
     &   /zeta(i1,i2 ) /zeta(i3,i4 ) 
     &  *sqrt( (2*xj1+1)* (2*xj2+1)* (2*xj3+1)* (2*xj4+1))  
     
      return
      end          
C==========================================================================

      subroutine jsq(i,V)
C
C  computes antisymmeterized matrix element of J(1)*J(2)
C
C  INPUT:
C       I:    which matrix element in list is it
C	NT(i):  PRINCIPLE radial Q# (not NODAL #) for i'th s.p. orbit
C 			N = 2n+l, where n is nodal q#
C	LT(i):  orbital ang mom for i'th s.p. orbit
C	JT(i):  2xJ for i'th s.p. orbit
C	JTOT(k):  total J for k'th tbme
C	TTOT(k):  total T for k'th TBME
C	MAP(k,i): k'th matrix element < 1 2 | V| 3 4 > i = 1,4 
C				  maps s.p. states
C
C  OUTPUT:
C	V:  matrix element
C
C  CALLS:
C	delta(i,j)	delta function
C       zeta(i,j) = sqrt(1+delta)
C	sixj
C  

            use sporbit
      use tbmes
      implicit none

      integer i
      real v,vjj
      integer J,T

C---- SINGLE-PARTICLE LIST ------
      integer i1,i2,i3,i4
      real xj1,xj2,xj3,xj4
      real xl1,xl2,xl3,xl4  
      integer l1,l2,l3,l4
      integer n1,n2,n3,n4     
     
C---------functions called

      integer delta
      real sj  ! found in libra
      real zeta

      J = vtbme(i)%j
      T = vtbme(i)%t
      i1 = vtbme(i)%n(1)
      i2 = vtbme(i)%n(2)
      i3 = vtbme(i)%n(3)
      i4 = vtbme(i)%n(4)
      n1 = orbqn(i1)%nr
      n2 = orbqn(i2)%nr
      n3 = orbqn(i3)%nr
      n4 = orbqn(i4)%nr

      xj1 = float(orbqn(i1)%j)/2.
      xj2 = float(orbqn(i2)%j)/2.
      xj3 = float(orbqn(i3)%j)/2.
      xj4 = float(orbqn(i4)%j)/2.

      l1 = orbqn(i1)%l
      l2 = orbqn(i2)%l
      l3 = orbqn(i3)%l
      l4 = orbqn(i4)%l

      xl1 = float(orbqn(i1)%l)
      xl2 = float(orbqn(i2)%l)
      xl3 = float(orbqn(i3)%l)
      xl4 = float(orbqn(i4)%l)

      if(i1 /= i3 .or. i2 /= i4)return

      vjj = float(j*(j+1)) - xj1*(xj1+1.) - xj2*(xj2+1)
!    vjj = sqrt( (2*xj1+1)*(2*xj2+1))* sqrt(xj1*(xj1+1)*xj2*(xj2+1)) 
!     & * (-1)**( j + int(xj1+xj2) )
C -- direct piece first
!     &  *(delta(i1,i3)*delta(i2,i4)
!     &   *sj(float(j),xj2,xj1,1.0,xj3,xj4)
!     &   -(-1)**(2+ int( j+t+ xj3-xj4))*
C---indirect piece 
!     & delta(i1,i4)*delta(i2,i3)*
!     &   sj(float(j),xj2,xj1,1.0,xj4,xj3) )
!       vjj=vjj /( zeta(i1,i2) * zeta(i3,i4) )
      vtbme(i)%v = vtbme(i)%v+ v*vjj  !*2  
      return
      end          
C==========================================================================
      
      subroutine tsq(i,V)
C
C  computes antisymmeterized matrix element of T(1)*T(2)
C
C  INPUT:
C       I:    which matrix element in list is it
C	NT(i):  PRINCIPLE radial Q# (not NODAL #) for i'th s.p. orbit
C 			N = 2n+l, where n is nodal q#
C	LT(i):  orbital ang mom for i'th s.p. orbit
C	JT(i):  2xJ for i'th s.p. orbit
C	JTOT(k):  total J for k'th tbme
C	TTOT(k):  total T for k'th TBME
C	MAP(k,i): k'th matrix element < 1 2 | V| 3 4 > i = 1,4 
C				  maps s.p. states
C
C  OUTPUT:
C	V:  matrix element
C
C  CALLS:
C	delta(i,j)	delta function
C       zeta(i,j) = sqrt(1+delta)
C	sixj
C  

            use sporbit
      use tbmes
      implicit none

      integer i
      real v,vjj
      integer J,T

C---- SINGLE-PARTICLE LIST ------
      integer i1,i2,i3,i4
      real xj1,xj2,xj3,xj4
      real xl1,xl2,xl3,xl4  
      integer l1,l2,l3,l4
      integer n1,n2,n3,n4     
     
C---------functions called

      integer delta
      real sj  ! found in libra
      real zeta

      J = vtbme(i)%j
      T = vtbme(i)%t
      i1 = vtbme(i)%n(1)
      i2 = vtbme(i)%n(2)
      i3 = vtbme(i)%n(3)
      i4 = vtbme(i)%n(4)


      if(i1 /= i3 .or. i2 /= i4)return

      vjj = float(t*(t+1)) - 1.5

      vtbme(i)%v = vtbme(i)%v+ v*vjj  !*2  
      return
      end          
C==========================================================================
!
!  computes single-particle energies from "phonon" read from file
!
      subroutine spephonon(a,xscale,xspe)
		  use sporbit
		  use tbmes
		  implicit none
		  integer a
		  real xscale,xspe
		  integer b
		  real xj
		  
		  xj = orbqn(a)%j*0.5
		  
		  do b = 1,numorb
			  xspe = xspe + 0.5*xscale*obme(a,b)**2/(2*xj+1.)
		  end do
		  return
	  end subroutine spephonon
!=========================================================
C==========================================================================
      
      subroutine v2phonon(i,jtran,ttran,V)
C
C  computes antisymmeterized matrix element of a phonon
C
C  INPUT:
C       I:    which matrix element in list is it
C	NT(i):  PRINCIPLE radial Q# (not NODAL #) for i'th s.p. orbit
C 			N = 2n+l, where n is nodal q#
C	LT(i):  orbital ang mom for i'th s.p. orbit
C	JT(i):  2xJ for i'th s.p. orbit
C	JTOT(k):  total J for k'th tbme
C	TTOT(k):  total T for k'th TBME
C	MAP(k,i): k'th matrix element < 1 2 | V| 3 4 > i = 1,4 
C				  maps s.p. states
C
C  OUTPUT:
C	V:  matrix element
C
C  CALLS:
C	delta(i,j)	delta function
C       zeta(i,j) = sqrt(1+delta)
C	sixj
C  

      use sporbit
      use tbmes
      implicit none

      integer i
	  integer jtran,ttran
	  real xscale
      real v,vjj
      integer J,T

C---- SINGLE-PARTICLE LIST ------
      integer i1,i2,i3,i4
      real xj1,xj2,xj3,xj4
      real xl1,xl2,xl3,xl4  
      integer l1,l2,l3,l4
      integer n1,n2,n3,n4     
     
C---------functions called

      integer delta
      real sj  ! found in libra
      real zeta

      J = vtbme(i)%j
      T = vtbme(i)%t
      i1 = vtbme(i)%n(1)
      i2 = vtbme(i)%n(2)
      i3 = vtbme(i)%n(3)
      i4 = vtbme(i)%n(4)
	  xj1 = 0.5*orbqn(i1)%j
	  xj2 = 0.5*orbqn(i2)%j
	  xj3 = 0.5*orbqn(i3)%j
	  xj4 = 0.5*orbqn(i4)%j

      vjj = sj(xj1,xj2,float(j),xj4,xj3,float(jtran))*obme(i1,i3)
     &  *obme(i2,i4)
      vjj =vjj-(-1)**(j+T)*sj(xj1,xj2,float(j),xj3,xj4,float(jtran))*
     & obme(i1,i4)*obme(i2,i3)
	   vjj = (-1)**(j +( orbqn(i2)%j+orbqn(i3)%j)/2 )*vjj
	   vjj = vjj /zeta(i1,i2)/zeta(i3,i4)
	   if(ttran==1)vjj = vjj*(-1)**(1-T)/(2*T+1)

      vtbme(i)%v = vtbme(i)%v+ v*vjj  !*2  
      return
      end subroutine v2phonon 
C==========================================================================

