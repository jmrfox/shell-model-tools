      subroutine genMultipole(K,iso,sizej, j1body, nsp,
     &  sizeme,opme,sizetbme,nme,jtot,ttot,map,mme,mspe)
    
C
C  computes antisymmeterized two-body matrix elements, and 
C  associated s.p. energies, from multipole-multipole interaction 
C
C  INPUT:
C       K : "multipolarity"
C	iso: = "isospin" of "multipole"
C       sizej    : dim of j1body
C	J1body(i):  2xJ for i'th s.p. orbit
C	NSP: 	# of single particle states
C       sizeme:	dimension of opme(i,j)
C       opme(i,j):  "multipole" 1-body matrix elements
C       nme  : # of two body-matrix elements
C       sizetbme:  dimension of jtot,ttot,map
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
C    genMMtbme
C    genMMspe
C

      implicit none

C................ INPUT...................

      integer k,iso
      integer sizej
      integer j1body(sizej)
      integer nsp
      integer sizeme  ! dimension of array opme
      real opme(sizeme,sizeme)

      
      integer sizetbme
      integer nme
      integer jtot(sizetbme)
      integer ttot(sizetbme)
      integer map(sizetbme,4)

C.................OUTPUT..................

      real mme(sizetbme)
      real mspe(sizej)

C................. INTERMEDIATE...........

      integer ia,ib,ic,id
      integer n
      real xme

C.........................................
      write(6,*)' Multipolarity = ',K
      do n =1,nme
        ia = map(n,1)
        ib = map(n,2)
        ic = map(n,3)
        id = map(n,4)
	call genMMtbme(K,iso,opme,ia,ib,ic,id,j1body,
     & sizeme,sizej,Jtot(n),Ttot(n),xme)
           mme(n)=xme

      enddo

      do ia =1,nsp
	call genMMspe(opme,ia,j1body,sizeme,sizej,nsp,xme)
 	mspe(ia) =xme
      enddo

      return 
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine genMMtbme(K,iso,opme,ia,ib,ic,id,j1body,
     & sizeme,sizej,J,T,xme)
C
C  this program accepts a "multipole" operator
C and computes a "multipole-multipole" 2-body interaction
C
C  let O(K,iso) be a "mutipole" operator; read in 
C  the doubly-reduced matrix elements < a ||| O(K,iso) ||| b >
C
C  then compute the antisymmeterized matrix elements
C  < a b; JT ||| 1/2 0(k,iso)^dagger * O(k,iso) ||| c d ; JT >
C
C   INPUT:  
C        K,iso:	rank, isospin of "multipole"
C       opme(i,j):  "multipole" 1-body matrix elements
C       ia,ib,ic,id: labels for two-body matrix elements
C       j1body(i) : ang moment for one-body states
C       J,T	: for 2 body matrix element
C       sizeme:	dimension of opme(i,j)
C        sizej: dimension of j1body(i)
C
C   OUTPUT:
C        xme:    two-body matrix element
C
C   functions called: 
C       sixj

      implicit none

C....................INPUT ................................

      integer K, iso    ! ang mom and isospin rank of "multipole"
      integer sizeme  ! dimension of array opme
      real opme(sizeme,sizeme)
    
      integer ia,ib,ic,id
      integer J,T
      integer sizej
      integer j1body(sizej)

C...................OUTPUT...........................

      real xme	                                                                                     

C..................functions called................

      real sixj
      real xja,xjb,xjc,xjd
      real xJ,xK

C..................intermediate.....................

      real zetaab,zetacd

C..................................................

      zetaab=1.0
      zetacd=1.0
      if(ia.eq.ib)zetaab=sqrt(2.)
      if(ic.eq.id)zetacd=sqrt(2.)

      xj = float(j)/2.
      xk = float(k)

      xja = 0.5*float(j1body(ia))
      xjb = 0.5*float(j1body(ib))
      xjc = 0.5*float(j1body(ic))
      xjd = 0.5*float(j1body(id))

      xme = opme(ia,ic)*opme(ib,id)*sixj(xja,xjb,xj,xjd,xjc,xk)
     &-(-1)**((j+t)/2+2)*opme(ia,id)*opme(ib,ic)*
     & sixj(xja,xjb,xj,xjc,xjd,xk)
 
      xme = 0.5*xme/zetaab/zetacd*(-1)**( (J+j1body(ib)+j1body(ic))/2)
C      write(6,*)xme
C      write(6,*)j,t
C      write(6,*)xja,xjb,xj,xjd,xjc,xme
      if(iso.eq.1)then
        xme = -xme/(2.*t+1.)
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  
      subroutine genMMspe(opme,ia,j1body,sizeme,sizej,nsp,xspe)
C
C  this program accepts a "multipole" operator
C and computes a "multipole-multipole" single-particle energy
C
C  let O(K,iso) be a "mutipole" operator; read in 
C  the doubly-reduced matrix elements < a ||| O(K,iso) ||| b >
C
C
C   INPUT:  
C       opme(i,j):  "multipole" 1-body matrix elements
C       ia: labels for spe matrix elements
C       j1body(i) : ang moment for one-body states
C       sizeme:	dimension of opme(i,j)
C        sizej: dimension of j1body(i)
C       nsp: # of single-particle states
C
C   OUTPUT:
C        xspe:    single-particle energy
C

      implicit none

C....................INPUT ................................

      integer sizeme  ! dimension of array opme
      real opme(sizeme,sizeme)
    
      integer ia,ib
      integer J,T
      integer sizej
      integer j1body(sizej)
      integer nsp

C....................OUTPUT................................


      real xspe

C.......................................................

   
      xspe = 0.0
      do ib = 1,nsp
	xspe = xspe+opme(ia,ib)**2
      enddo
      xspe = 0.25*xspe/float(j1body(ia)+1)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc

      subroutine make_random_multipole(iseed,K,j1body,sizej,nsp,
     & sizeme,opme)
C
C   INPUT:  
C      iseed:  random seed
C      K: ang momentum rank of interaction
C       j1body(i) : ang moment for one-body states
C       sizeme:	dimension of opme(i,j)
C        sizej: dimension of j1body(i)
C       nsp: # of single-particle states
C
C   OUTPUT:
C       opme(i,j):  "multipole" 1-body matrix elements
C 
C   SUBROUTINES CALLED
C	GAUSS: generates gaussian random numbers for noise 
C

C....................INPUT ................................

      integer iseed
      integer K
      integer sizej
      integer j1body(sizej)
      integer nsp

C....................OUTPUT................................

      integer sizeme  ! dimension of array opme
      real opme(sizeme,sizeme)

C.................functions called
      real gauss

C ............. INTERMEDIATE.......................

      integer ia,ib
      integer phase

      do ia =1,nsp
          do ib =ia,nsp
            if(j1body(ia)+j1body(ib).ge.2*K .and.
     &    abs(j1body(ia)-j1body(ib)).le.2*K)then
                opme(ia,ib)=gauss(iseed)
                opme(ib,ia)=opme(ia,ib)*
     &  (-1)**((j1body(ia)-j1body(ib))/2+2)
            else
                opme(ia,ib)=0.0
                opme(ib,ia)=0.0
            endif
           enddo
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine make_spin_multipole(j1body,l1body,n1body,sizej,nsp,
     & sizeme,opme)
C
C  notes this assumes (l,1/2)j ordering
C
C   INPUT:  
C       j1body(i) :2x ang moment for one-body states
C       l1body(i):  org ang momentum for one-body states
C       n1body(i):  PRINCIPAL quantum number
C       sizeme:	dimension of opme(i,j)
C        sizej: dimension of j1body(i)
C       nsp: # of single-particle states
C
C   OUTPUT:
C       opme(i,j):  "multipole" 1-body matrix elements
C 
C   SUBROUTINES CALLED
C	sixj
C

C....................INPUT ................................

      integer sizej
      integer j1body(sizej)
      integer n1body(sizej)
      integer l1body(sizej)
      integer nsp

C....................OUTPUT................................

      integer sizeme  ! dimension of array opme
      real opme(sizeme,sizeme)

C----------------- functions called--------------

      real sixj

C-----------------intermediate

      integer ia,ib
      real xja,xjb

      do ia =1,nsp
          xja=float(j1body(ia))/2.
          do ib =ia,nsp
            if( ( j1body(ia)+j1body(ib).ge.2) .and.
     &    (abs(j1body(ia)-j1body(ib)).le.2) .and.
     &     (l1body(ia).eq.l1body(ib)) .and.
     &    (n1body(ia).eq.n1body(ib)) )then
                xjb=float(j1body(ib))/2.
                opme(ia,ib)=sixj(0.5,xja,float(l1body(ia)),
     &  xjb,0.5,1.0) * sqrt(6.)*  sqrt( (2*xja+1.)*(2*xjb+1.) )*
     & (-1)**(l1body(ia)+1+(j1body(ia)+2)/2)
                opme(ib,ia)=opme(ia,ib)*
     &  (-1)**((j1body(ia)-j1body(ib))/2+2)
            else
                opme(ia,ib)=0.0
                opme(ib,ia)=0.0
            endif
           enddo
      enddo


      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine make_std_multipole(L,j1body,l1body,n1body,sizej,nsp,
     & sizeme,opme)
C
C  notes this assumes (l,1/2)j ordering
C
C   INPUT:  
C      L : order of multipole
C       j1body(i) :2x ang moment for one-body states
C       l1body(i):  org ang momentum for one-body states
C       n1body(i):  PRINCIPAL quantum number
C       sizeme:	dimension of opme(i,j)
C        sizej: dimension of j1body(i)
C       nsp: # of single-particle states
C
C   OUTPUT:
C       opme(i,j):  "multipole" 1-body matrix elements
C 
C   SUBROUTINES CALLED
C	sixj
C	horadint

C....................INPUT ................................

      integer L
      integer sizej
      integer j1body(sizej)
      integer n1body(sizej)
      integer l1body(sizej)
      integer nsp

C....................OUTPUT................................

      integer sizeme  ! dimension of array opme
      real opme(sizeme,sizeme)

C----------------- functions called--------------

      real sixj,threej

C-----------------intermediate

      integer ia,ib

      real xja,xjb,xla,xlb
      integer na,nb	! radial quantum number

      do ia =1,nsp
          xja=float(j1body(ia))/2.
          na = (n1body(ia)-l1body(ia))/2
          xla = float(l1body(ia))
          do ib =ia,nsp
              
            if( ( j1body(ia)+j1body(ib).ge.2*L) .and.
     &    (abs(j1body(ia)-j1body(ib)).le.2*L) .and.
     &     (l1body(ia)+l1body(ib).ge.L) .and.
     &     (abs(l1body(ia)-l1body(ib)).le.L) )then
                xjb=float(j1body(ib))/2
                 xlb=float(l1body(ib))
                nb =(n1body(ib)-l1body(ib))/2
C               write(6,*)xla,float(l),xlb
                opme(ia,ib)=sixj(xla,xja,0.5,
     &  xjb,xlb,float(L))* sqrt( 2.*L+1.) /sqrt(4.*3.1415926) *
C     &  xjb,xlb,float(L))* ! sqrt( 2.*L+1.) /sqrt(4.*3.1415926) *  ! to agree with other qq
     & sqrt( (2*xja+1.)*(2*xjb+1.)*(2*xla+1.)*(2*xlb+1.) )*
     & threej(xla,float(L),xlb,0.0, 0.0, 0.0) *
     & (-1)**(L+(j1body(ib)+2)/2)
     & * horadint(L,na,l1body(ia),nb,l1body(ib))
C                write(6,*)' me' ,ia,ib,opme(ia,ib) 
                opme(ib,ia)=opme(ia,ib)*
     &  (-1)**((j1body(ia)-j1body(ib))/2+2)
            else
                opme(ia,ib)=0.0
                opme(ib,ia)=0.0
            endif
           enddo
      enddo


      return
      end

