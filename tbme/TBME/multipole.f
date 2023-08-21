
      subroutine Multipole(L,T,nsp,nt,lt,jt,nme,jtot,ttot,map,mme,mspe,
     & active,convention)
C
C  computes antisymmeterized two-body matrix elements, and 
C  associated s.p. energies, from multipole-multipole interaction 
C
C  NOTE: THIS ASSUMES coupling (l,1/2)=j, NOT (1/2,l)=j
C
C  INPUT:
C	NSP: 	# of single particle states
C	T: = "isospin" scalar or isoscalar densities (note: not TRUE isospin!)
C	NT(i):  PRINCIPLE radial Q# (not NODAL #) for i'th s.p. orbit
C 			N = 2n+l, where n is nodal q#
C	LT(i):  orbital ang mom for i'th s.p. orbit
C	JT(i):  2xJ for i'th s.p. orbit
C       nme  : # of two body-matrix elements
C	JTOT(k):  total J for k'th tbme
C	TTOT(k):  total T for k'th TBME
C	MAP(k,i): k'th matrix element < 1 2 | V| 3 4 > i = 1,4 
C				  maps s.p. states
C
C  OUTPUT:
C	MME(k):  k'th tbme
C	MSPE(i):  i'th s.p. energy
C
C  CALLS:
C	multipoleME: computes in density-density format
C	sixj
C	delta

      implicit none
C----- INTERACTION

      integer L 
      integer T
      logical active(100) 		! which s.p. states are active
      logical convention		! if true then alternate phase on ls
C---- SINGLE-PARTICLE LIST ------

      integer nsp			! # of s.p. states
      integer jt(100),nt(100),lt(100)  ! 2*J, N, L of s.p. state
      real xj1,xj2,xj3,xj4
      integer l1,l2,l3,l4
      integer n1,n2,n3,n4
      
      real sixj
      integer delta
      real zeta
      real fact
      
      integer maxj,minj 
C----TBME LIST -------

      integer nme			! # of TBME
      integer map(5000,4)		! maps TBME to 4 s.p. states
      integer jtot(5000),ttot(5000) 	! JT of i'th TBME
      real mme(5000) 			! two-body matrix elements
      real mspe(100)			! single-particle energies
      
            
C----dummy integers 

      integer i,j,k
      real vv
      real temp1,temp2
      
C      write(6,*)active(1)
      do i = 1,nme 
        
      	vv = 0.0
      	J = jtot(i)/2
        xj1 = float(jt(map(i,1)))/2.
        xj2 = float(jt(map(i,2)))/2.
        xj3 = float(jt(map(i,3)))/2.
        xj4 = float(jt(map(i,4)))/2.

        l1 = (lt(map(i,1)))
        l2 = (lt(map(i,2)))
        l3 = (lt(map(i,3)))
        l4 = (lt(map(i,4)))
 	
 	n1 = (nt(map(i,1))-lt(map(i,1)))/2
 	n2 = (nt(map(i,2))-lt(map(i,2)))/2
 	n3 = (nt(map(i,3))-lt(map(i,3)))/2
 	n4 = (nt(map(i,4))-lt(map(i,4)))/2
     
        fact = 1./zeta(map(i,1),map(i,2)) /zeta(map(i,3),map(i,4)) 
     &  * (-1)**( int(J +xj2+xj3 )+2  )*(2*L+1)

        call multipoleME(L,J,n1,l1,xj1,n2,l2,xj2,
     & n3,l3,xj3,n4,l4,xj4,temp1) 
        call multipoleME(L,J,n1,l1,xj1,n2,l2,xj2,
     & n4,l4,xj4,n3,l3,xj3,temp2)          
        temp2 = temp2*(-1)**(J+2 + ttot(i)/2)
        
            vv = temp1*sixj(xj1,xj2,float(j),xj4,xj3,float(L))
     &    -temp2*sixj(xj1,xj2,float(j),xj3,xj4,float(L))
C        if(ttot(i).eq.0)vv=vv*(-1)**(2+t)  ! additional phase 
        vv =vv*(-1)**(2+t)
      	mme(i) = vv*fact
      	if(convention)mme(i)=mme(i)*(-1)**(2+int(xj1+xj2+xj3+xj4))
C        if(ttot(i).eq.0)then
C             write(6,*)i,j,ttot(i)/2,temp1,temp2,vv*fact
C     & ,sixj(xj1,xj2,float(j),xj4,xj3,float(L))
C     &  ,sixj(xj1,xj2,float(j),xj3,xj4,float(L))
C         endif
      enddo
      
      do i = 1,nsp
	vv = 0.0

        do k = 1,nsp
        if(active(k))then
        xj1 = float(jt(i))/2.
        xj2 = float(jt(k))/2.

        l1 = (lt(i))
        l2 = (lt(k))
 	
 	n1 = (nt(i)-lt(i))/2
 	n2 = (nt(k)-lt(k))/2
        maxj=int(xj1+xj2)
        minj = int(abs(xj1-xj2))
        call multipoleME(L,J,n1,l1,xj1,n2,l2,xj2,
     & n2,l2,xj2,n1,l1,xj1,temp1) 
        do J = minj,maxj 
        

	  vv = vv+0.5*(2*j+1)*(2*L+1)*(-1)**(int(xj1+xj2+1))
     &  / (2.*xj1+1.) *sixj(xj1,xj2,float(J),xj1,xj2,float(l))
     & *temp1
	
	enddo
	endif
	enddo
	mspe(i) = vv
      
      enddo
      return
      end          

C==========================================================================

      subroutine multipoleME(L,J,n1p,l1p,xj1p,n2p,l2p,xj2p,
     & n1,l1,xj1,n2,l2,xj2,tbme) 
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
      real sixj,threej,HORadInt
      real v1,v2

      tbme =  1./float(2*L+1) !   (4.*3.1415926) 
     &  *sqrt((2*xj1+1)*(2*xj2+1)*(2*xj1p+1)*(2*xj2p+1))
     &  *sqrt((2.*l1p+1.)*(2.*l2p+1.)*(2.*l1+1)*(2*l2+1))
     &  *sixj(float(l1p),xj1p,0.5,xj1,float(l1),float(L))
     &  *sixj(float(l2p),xj2p,0.5,xj2,float(l2),float(L))
     &  *threej(float(l1p),float(L),float(l1),0.,0.0,0.)
     &  *threej(float(l2p),float(L),float(l2),0.,0.0,0.)
      v1 = HORadInt(L,n1,l1,n1p,l1p)
      v2 = HORadInt(L,n2,l2,n2p,l2p)
C      write(6,*)' to date ',tbme,v1,v2
      tbme = tbme*v1*v2
     & *(-1)**(1+int(xj1+xj2))

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

      subroutine ss(i,Vss)
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
      real vss
      integer J,T

C---- SINGLE-PARTICLE LIST ------
      integer i1,i2,i3,i4
      real xj1,xj2,xj3,xj4
      real xl1,xl2,xl3,xl4  
      integer l1,l2,l3,l4

      
      
C---------functions called

      integer delta
      real sj  ! found in libra
 
      J = vtbme(i)%j
      T = vtbme(i)%t
      i1 = vtbme(i)%n(1)
      i2 = vtbme(i)%n(2)
      i3 = vtbme(i)%n(3)
      i4 = vtbme(i)%n(4)

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
      vss =   ( delta(i1,i3)* delta(i2,i4) * delta(l1,l3)*delta(l2,l4)
     &   (-1)**( int(l1+l2+ xj1+xj3)+2 )*
     &   sj(0.5, xj1,xl1,xj3,.5,1.)*sixj(.5,xj2,xl2,xj4,.5,1.)
     &   *sj(float(J).,xj2,xj1,1.0,xj3,xj4)
C      write(6,*)vss
C      vss = vss
     &   -(-1)**(2+ int( J+T +xj3-xj4))*
C---indirect piece      
     &   delta(i1,i4)*delta(i2,i3)*delta(l1,l4)*delta(l2,l3)
     &   (-1)**( int(l1+l2+ xj1+xj4) +2)*
     &   sj(float(J),xj2,xj1,1.0,xj4,xj3)*
     &   sj(0.5, xj1,xl1,xj4,.5,1.)*sixj(.5,xj2,xl2,xj3,.5,1.) )
       
       vss = vss *3. /sqrt( float ( (1+ delta(i1,i2 ) )
     & * (1+delta(i3,i4) ) ))
     &  *sqrt( (2*xj1+1)* (2*xj2+1)* (2*xj3+1)* (2*xj4+1))  * 
     & (-1)**( J)
C      write(6,*)vss    ,sixj(float(jtot(i))/2.,xj2,xj1,1.0,xj4,xj3)
     
      return
      end          

C==========================================================================

      subroutine ll(i,nt,lt,jt,jtot,ttot,map,V)
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
      implicit none

      integer i
      real v
      
C---- SINGLE-PARTICLE LIST ------

      integer jt(100),nt(100),lt(100)  ! 2*J, N, L of s.p. state
      real xj1,xj2,xj3,xj4
      real xl1,xl2,xl3,xl4  

C----TBME LIST -------

      integer map(5000,4)		! maps TBME to 4 s.p. states
      integer jtot(5000),ttot(5000) 	! JT of i'th TBME

C---------functions called

      integer delta
      real zeta
      real sixj 
      
      xj1 = float(jt(map(i,1)))/2.
      xj2 = float(jt(map(i,2)))/2.
      xj3 = float(jt(map(i,3)))/2.
      xj4 = float(jt(map(i,4)))/2.

      xl1 = float(lt(map(i,1)))
      xl2 = float(lt(map(i,2)))
      xl3 = float(lt(map(i,3)))
      xl4 = float(lt(map(i,4)))
 
C -- direct piece first
      v =   ( delta(nt(map(i,1)),nt(map(i,3)))* 
     &    delta(nt(map(i,2)),nt(map(i,4)))*
     &    delta(lt(map(i,1)),lt(map(i,3)))*
     &    delta(lt(map(i,2)),lt(map(i,4)))
     &   *sixj(float(jtot(i))/2.,xj2,xj1,1.0,xj3,xj4)
     &   *(-1)**( int( xj2+xj4)+2 )*
     &   sixj(xl3, xj1,0.5,xj3,xl1,1.)*sixj(xl2,xj2,0.5,xj4,xl2,1.)
C
     &   -(-1)**(2+ int( jtot(i)+ttot(i)+jt(map(i,3))-jt(map(i,4)))/2)*
C---indirect piece      
     &   delta(nt(map(i,1)),nt(map(i,4)))* 
     &    delta(nt(map(i,2)),nt(map(i,3)))*
     &    delta(lt(map(i,1)),lt(map(i,4)))*
     &    delta(lt(map(i,2)),lt(map(i,3)))*
     &   (-1)**( int( xj2+xj3) +2)*
     &   sixj(float(jtot(i))/2.,xj2,xj1,1.0,xj4,xj3)*
     &   sixj(xl3, xj1,0.5,xj4,xl1,1.)*sixj(xl2,xj2,0.5,xj3,xl2,1.) )
       v = v /zeta(map(i,1),map(i,2) ) /zeta(map(i,3),map(i,4) ) 
     &  *sqrt( (2*xj1+1)* (2*xj2+1)* (2*xj3+1)* (2*xj4+1))  
     &  * sqrt( xl1*xl2*(xl1+1.)*(xl2+1)*(2*xl1+1)*(2*xl2+1) )
     &  *(-1)**( (jtot(i)+4 )/2+lt(map(i,1))+lt(map(i,2)) )
     
      return
      end          
C==========================================================================
      

      subroutine ls(i,nt,lt,jt,jtot,ttot,map,V)
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
      implicit none

      integer i
      real v
      
C---- SINGLE-PARTICLE LIST ------

      integer jt(100),nt(100),lt(100)  ! 2*J, N, L of s.p. state
      real xj1,xj2,xj3,xj4
      real xl1,xl2,xl3,xl4  

C----TBME LIST -------

      integer map(5000,4)		! maps TBME to 4 s.p. states
      integer jtot(5000),ttot(5000) 	! JT of i'th TBME

C---------functions called

      integer delta
      real zeta
      real sixj 
      
      xj1 = float(jt(map(i,1)))/2.
      xj2 = float(jt(map(i,2)))/2.
      xj3 = float(jt(map(i,3)))/2.
      xj4 = float(jt(map(i,4)))/2.

      xl1 = float(lt(map(i,1)))
      xl2 = float(lt(map(i,2)))
      xl3 = float(lt(map(i,3)))
      xl4 = float(lt(map(i,4)))
 
C -- direct piece first
      v =   ( delta(nt(map(i,1)),nt(map(i,3)))* 
     &    delta(nt(map(i,2)),nt(map(i,4)))*
     &    delta(lt(map(i,1)),lt(map(i,3)))*
     &    delta(lt(map(i,2)),lt(map(i,4)))
     &   *sixj(float(jtot(i))/2.,xj2,xj1,1.0,xj3,xj4)
     & *(  sqrt( (2*xl1+1)*xl1*(xl1+1) )
     & *  sixj(xl3, xj1,0.5,xj3,xl1,1.)*sixj(0.5,xj2,xl2,xj4,0.5,1.)
     & + (-1)**(int( xj1+xj2+xj3+xj4+2))*sqrt((2*xl2+1)*xl2*(xl2+1)) 
     & *  sixj(xl2, xj2,0.5,xj4,xl2,1.)*sixj(0.5,xj1,xl1,xj3,0.5,1.)
     & ) 
     &   -(-1)**(2+ int( jtot(i)+ttot(i)+jt(map(i,3))-jt(map(i,4)))/2)
C---indirect piece      
     &   *delta(nt(map(i,1)),nt(map(i,4)))* 
     &    delta(nt(map(i,2)),nt(map(i,3)))*
     &    delta(lt(map(i,1)),lt(map(i,4)))*
     &    delta(lt(map(i,2)),lt(map(i,3)))
     &    *sixj(float(jtot(i))/2.,xj2,xj1,1.0,xj4,xj3)
     & *(  sqrt( (2*xl1+1)*xl1*(xl1+1) )
     & *  sixj(xl1, xj1,0.5,xj4,xl1,1.)*sixj(0.5,xj2,xl2,xj3,0.5,1.)
     & + (-1)**(int( xj1+xj2+xj3+xj4+2))*sqrt((2*xl2+1)*xl2*(xl2+1)) 
     & *  sixj(xl2, xj2,0.5,xj3,xl2,1.)*sixj(0.5,xj1,xl1,xj4,0.5,1.)
     & )
     & )

       v = v*(-1)**(jtot(i)/2 +int(xl1+xl2)+1)*sqrt(1.5)  
     &   /zeta(map(i,1),map(i,2) ) /zeta(map(i,3),map(i,4) ) 
     &  *sqrt( (2*xj1+1)* (2*xj2+1)* (2*xj3+1)* (2*xj4+1))  
     
      return
      end          
C==========================================================================
      

      subroutine jsq(i,nt,lt,jt,jtot,ttot,map,V)
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
      implicit none

      integer i
      real v
      
C---- SINGLE-PARTICLE LIST ------

      integer jt(100),nt(100),lt(100)  ! 2*J, N, L of s.p. state
      real xj1,xj2,xj3,xj4

C----TBME LIST -------

      integer map(5000,4)		! maps TBME to 4 s.p. states
      integer jtot(5000),ttot(5000) 	! JT of i'th TBME

C---------functions called

      integer delta
      real zeta
      real sixj 
      
      xj1 = float(jt(map(i,1)))/2.
      xj2 = float(jt(map(i,2)))/2.
      xj3 = float(jt(map(i,3)))/2.
      xj4 = float(jt(map(i,4)))/2.
      v = sqrt( (2*xj1+1)*(2*xj2+1))* sqrt(xj1*(xj1+1)*xj2*(xj2+1)) 
     & * (-1)**( (jtot(i)+4)/2 + int(xj1+xj2) )
C -- direct piece first
     &  *( delta(map(i,1),map(i,3))* delta(map(i,2),map(i,4)) 
     &   *sixj(float(jtot(i))/2.,xj2,xj1,1.0,xj3,xj4)
     &   -(-1)**(2+ int( jtot(i)+ttot(i)+jt(map(i,3))-jt(map(i,4)))/2)*
C---indirect piece 
     & delta(map(i,1),map(i,4))* delta(map(i,2),map(i,3)) *
     &   sixj(float(jtot(i))/2.,xj2,xj1,1.0,xj4,xj3) )
       v = v /( zeta(map(i,1),map(i,2)) * zeta(map(i,3),map(i,4)))
     
      return
      end          



      

