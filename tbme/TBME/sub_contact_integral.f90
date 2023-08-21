!=====================================================
!
! file sub_contract_integral.f90
!
! routines for computing and integrals for contact interactions
!
! initiated 11/09 CWJ @ SDSU 
!
!========================================================

  subroutine contact_master(describe)
  use sporbit
  use tbmes

  implicit none

  logical describe
  integer Npts
  real dr
  integer ncme
  real, allocatable :: cints(:) 
  real ::  bosc
  real :: V0, V1   ! isoscalar, isovector
  integer n
  

  print*,' Enter oscillator length b (in fm) '
  read*,bosc
  print*,' Enter isoscalar and isovector strengths (in MeV -fm^3 ) '
  read*,V0, V1
!  print*,' Enter # of points, dr '
!  read*,Npts,dr

  V0 = V0/bosc**3
  V1 = V1/bosc**3
  call binom2
  call trinom
!  call check_integrals(Npts,dr)
  ncme = numorb*(numorb+1)*(numorb+2)*(numorb+3)/24
  allocate(cints(nme))
!  call get_contact_integrals(Npts,dr,ncme,cints)

  do n = 1,nme
     call make_contact_tbme(n,ncme,cints,V0,V1)
  end do ! n
  return
  end subroutine contact_master

!==========================================================
  subroutine make_contact_tbme(ime,ncme,cints,V0,v1)
  use sporbit
  use tbmes

  implicit none

  integer ncme
  real :: cints(ncme) 
  real :: V0, V1   ! isoscalar, isovector
  integer ime
  integer indx
  integer a,b,c,d
  integer J, T, S
  real vdir,vexc
  real vme
  real zeta
  integer phase

  a = vtbme(ime)%n(1)
  b = vtbme(ime)%n(2)
  c = vtbme(ime)%n(3)
  d = vtbme(ime)%n(4)

  zeta = 1.
  if (a == b)zeta = zeta/sqrt(2.)
  if( c == d) zeta = zeta/sqrt(2.)

  J = vtbme(ime)%j
  T = vtbme(ime)%t
  if(T == 0)then
       S = 1
       vme = v0
  else
       S = 0
       vme = v1
  endif

  vdir = 0.0
  vexc = 0.0
  call contact_loops2(a,b,c,d,S,J,vdir)
  call contact_loops2(a,b,d,c,S,J,vexc)

!.......... GET PHASE FOR EXCHANGE TERM..........
  phase = (-1)**( T + J + (orbqn(c)%j + orbqn(d)%j)/2 )
!........ FIND INDEX BY SORTING a,b,c,d


!  if (a == 3 .and. b ==3 .and. c == 3 .and. d == 8)then
!     write(22,*)zeta,vdir,phase,vexc
!  endif

  call swap1(a,b)
  call swap1(b,c)
  call swap1(c,d)
  call swap1(b,c)
  call swap1(a,b)
  call swap1(b,c)
  call swap1(c,d)
  if( b < c .or. a < b .or. c < d) then
    print*,' huh ',a,b,c,d
    stop
  endif
  indx = d + c*(c-1)/2 + (b+1)*b*(b-1)/6 + (a+2)*(a+1)*a*(a-1)/24

  if(indx < 1 .or. indx > ncme)then
      print*,' woops ',indx,a,b,c,d
      stop
  endif
!-------------- NOW PUT ALL TOGETHER

!  vtbme(ime)%v = vme*zeta*cints(indx)*( vdir  + phase*vexc)
  vtbme(ime)%v = vme*zeta*( vdir  + phase*vexc)
!  write(95,*)zeta,vdir,phase,vexc
!  if(s == 0)write(96,*)a,b,c,d,vdir*cints(indx)

  end subroutine make_contact_tbme
!==========================================================
  subroutine swap1(i,j)
  integer i,j,k
  if(j > i)then
     k = j
     j = i
     i = k
  endif
  return
  end subroutine swap1
!==========================================================

  subroutine contact_loops(a,b,c,d,S,J,ans)

  use sporbit

  implicit none
  integer a,b,c,d
  real :: ans
  integer la,lb,lc,ld
  real :: xla,xlb,xlc,xld
  real :: xja,xjb,xjc,xjd
  integer :: lambda,lambdamin,lambdamax
  integer J, S
  real :: xj, xs
  integer :: L, Lmin,Lmax
  real :: xl, xlambda
  real :: fact1, fact2
  real :: tj,coefr9j,sj  ! functions called

  la = orbqn(a)%l
  xla = float(la)
  xja = float( orbqn(a)%j) /2.

  lb = orbqn(b)%l
  xlb = float(lb)
  xjb = float( orbqn(b)%j) /2.

  lc = orbqn(c)%l
  xlc = float(lc)
  xjc = float( orbqn(c)%j) /2.

  ld = orbqn(d)%l
  xld = float(ld)
  xjd = float( orbqn(d)%j) /2.

  xj = float(j)
  xs = float(s)

  Lmax = min(la + lb, lc + ld)
  Lmin = max( abs(la-lb),abs(lc-ld))
  Lmin = max( Lmin, abs(J -S ))
  Lmax = min( Lmax, J + S)

  if( Lmin > Lmax) return

  lambdamax = min( la + lc, lb + ld)
  lambdamin = max (abs(la-lc), abs(lb-ld))

  if(lambdamin > lambdamax) return

  ans = 0.0
  do L = Lmin,Lmax,2
     xl = float(l)
     fact1 = (2*L+1)*(-1)**L * &  
        coefr9j( xla, 0.5, xlb, 0.5, xja,xjb, xL, xS, xJ) * & 
        coefr9j( xlc, 0.5, xld, 0.5, xjc, xjd, xL, xS, xJ) 
     fact2 = 0.0
     do lambda = lambdamin, lambdamax,2
        xlambda = float(lambda)
        fact2 = fact2 + (2*lambda+1)  *(-1)**lambda   & 
            * sj(xla,xlc,xlambda,xld,xlb,xl)  & 
            * tj(xla,xlambda,xlc,0.,0.,0.)  & 
            * tj(xlb,xlambda,xld,0.,0.,0. )
     enddo ! lambda
     ans = ans + fact1*fact2
  enddo  ! L
!  check
  if( la == 0 .and. lb == 0 .and. lc == 0 .and. ld == 0)then
     if( abs(ans -0.25/(2*xs+1)) > 0.0001)then
        print*,' problem ',ans, s
        print*,a,b,c,d, J,S
        stop
     endif
  endif

  ans = ans*(2*s+1)/(4.*3.1415926) * & 
     sqrt(  (2*xja+1)*(2*xjb+1)*(2*xjc+1)*(2*xjd +1) * & 
            (2*xla+1)*(2*xlb+1)*(2*xlc+1)*(2*xld +1) )

  return

  end subroutine contact_loops


!==========================================================
!
!  this version uses conversion from relative coordinates
!
  subroutine contact_loops2(a,b,c,d,S,J,ans)

  use sporbit

  implicit none
  integer a,b,c,d
  real :: ans
  integer na,nb,nc,nd
  integer la,lb,lc,ld
  real :: xla,xlb,xlc,xld
  real :: xja,xjb,xjc,xjd
  integer :: ja,jb,jc,jd
  integer J, S
  real :: xj, xs
  integer :: L, Lmin,Lmax
  real :: xl
  real :: fact1, fact2
  real :: tj,coefr9j,sj  ! functions called
  real*8 :: c6j,c9j

  na = orbqn(a)%nr
  la = orbqn(a)%l
  ja = orbqn(a)%j
  xla = float(la)
  xja = float( orbqn(a)%j) /2.

  nb = orbqn(b)%nr
  lb = orbqn(b)%l
  jb = orbqn(b)%j
  xlb = float(lb)
  xjb = float( orbqn(b)%j) /2.

  nc = orbqn(c)%nr
  lc = orbqn(c)%l
  jc = orbqn(c)%j
  xlc = float(lc)
  xjc = float( orbqn(c)%j) /2.

  nd = orbqn(d)%nr
  ld = orbqn(d)%l
  jd = orbqn(d)%j
  xld = float(ld)
  xjd = float( orbqn(d)%j) /2.

  xj = float(j)
  xs = float(s)

  Lmax = min(la + lb, lc + ld)
  Lmin = max( abs(la-lb),abs(lc-ld))
  Lmin = max( Lmin, abs(J -S ))
  Lmax = min( Lmax, J + S)

  if( Lmin > Lmax) return

  ans = 0.0
  do L = Lmin,Lmax,2
     xl = float(l)
     fact1 = (2*L+1)* &  
         c9j(la*2,1,ja,lb*2,1,jb,2*L,2*S,2*J)* &
         c9j(lc*2,1,jc,ld*2,1,jd,2*L,2*S,2*J)
!        coefr9j( xla, 0.5, xlb, 0.5, xja,xjb, xL, xS, xJ) * & 
!        coefr9j( xlc, 0.5, xld, 0.5, xjc, xjd, xL, xS, xJ) 
     call rel_to_lab(na,la,nb,lb,nc,lc,nd,ld,L,fact2)
!     if(fact1 == 0.0)then
!       print*,la,ja,lb,jb,L,S,J
!       print*,lc,jc,ld,jd
!       print*,         c9j(la*2,1,ja,lb*2,1,jb,2*L,2*S,2*J), &
!        coefr9j( xla, 0.5, xlb, 0.5, xja,xjb, xL, xS, xJ), & 
!         c9j(lc*2,1,jc,ld*2,1,jd,2*L,2*S,2*J), & 
!        coefr9j( xlc, 0.5, xld, 0.5, xjc, xjd, xL, xS, xJ) 
!     endif
!     print*,fact1,fact2
     ans = ans + fact1*fact2
  enddo  ! L


  ans = ans*(2*s+1) * & 
     sqrt(  (2*xja+1)*(2*xjb+1)*(2*xjc+1)*(2*xjd +1) )

  return

  end subroutine contact_loops2


!==========================================================

  subroutine test_contact_integrals(Npts,dr)

  use sporbit
  implicit none
  integer Npts
  real dr
  integer nme
  real, allocatable :: cints(:)

  call check_integrals(Npts,dr)
  nme = numorb*(numorb+1)*(numorb+2)*(numorb+3)/24
  allocate(cints(nme))
  call get_contact_integrals(Npts,dr,nme,cints)

  return
  end subroutine test_contact_integrals
!========================================================

  subroutine get_contact_integrals(Npts,dr,nme,cints)
  use sporbit
  implicit none
  integer Npts
  real(kind = 4) :: dr
  integer :: nme 
  real(kind = 4) :: cints(nme)
  integer a,b,c,d
  integer na,nb,nc,nd  ! radial nodal quantum numbers
  integer la,lb,lc,ld  ! orbital ang momenta  
  integer l
  real(kind = 4) :: ans
  integer indx
  
  cints(:) = 0.0
  do a = 1,numorb
     na = orbqn(a)%nr
     la = orbqn(a)%l 
     do b = 1,a
        nb = orbqn(b)%nr
        lb = orbqn(b)%l 
        do c = 1,b
           nc = orbqn(c)%nr
           lc = orbqn(c)%l 
           do d = 1,c
              nd = orbqn(d)%nr
              ld = orbqn(d)%l 
!.................... CHECK PARITY CONSERVED
              L = la+lb+lc+ld
              if( (L/2)*2 /= L)cycle

              call basic_contact_integral(na,la,nb,lb,nc,lc,nd,ld,Npts,dr,ans)
              indx = d + c*(c-1)/2 + (b+1)*b*(b-1)/6 + (a+2)*(a+1)*a*(a-1)/24
              if(cints(indx) /= 0.0)then
                 print*,' whoopsie ',a,b,c,d,indx
                 stop
              endif
              cints(indx) = ans
           enddo ! d

        enddo ! c
     enddo ! b

  end do ! a

  return
  end subroutine get_contact_integrals
!========================================================

  subroutine check_integrals(Npts,dr)
!
! master subroutine to see if integrals are good enough
!
   use sporbit
   implicit none
   integer Npts
   real(kind = 4) :: dr

   integer a,b
   integer na,nb,la,lb

   do a = 1,numorb
     na = orbqn(a)%nr
     la = orbqn(a)%l
     do b = 1,a
        nb = orbqn(b)%nr
        lb = orbqn(b)%l
        call test_integral(na,la,nb,lb,Npts,dr)

     enddo ! b
   enddo ! a

   print*,' All integrals look okay '
   return
   end subroutine check_integrals

!========================================================

  subroutine basic_contact_integral(na,la,nb,lb,nc,lc,nd,ld,Npts,dr,ans)
!
! computes basis integral 
!
!  = int R_a(r)R_b(r) R_c(r) R_d(r) r^2 dr 
! INPUT
!  na,nb,nc,nd  = nodal quantum numbers
!  la,lb,lc,ld  = angular momenta
!  Npts = # of points in integral
!  dr = step size in dr
!
!  OUTPUT 
!    ans = ans 
!
!  SUBROUTINES CALLED
!      dquart_int    ! boole's rule for integration, double-precision
!
!  FUNCTIONS CALLED
!      MyRadHOwfn2(n,l,rd,bd,vald) double precision for radial wfn
!
!=============================================

   implicit none

   integer na,nb,nc,nd  ! radial nodal quantum numbers
   integer la,lb,lc,ld  ! orbital ang momenta
   integer Npts
   real(kind = 4) :: dr
   real (kind = 4) :: ans

!........... INTERMEDIATE

   real(kind =8), allocatable :: array(:)
   real(kind =8) :: dans
   real(kind = 8) :: ddr
   real(kind = 8) :: r
   real(kind =8)  :: vala,valb,valc,vald

   integer i

   if(.not.allocated(array)) allocate(array(0:Npts))
   ddr = dble(dr)
   array(0) = 0.d0
   do i = 1,npts
      r = i*ddr
      call MyRadHOwfn2(na,la,r,1.d0,vala)
      call MyRadHOwfn2(nb,lb,r,1.d0,valb)
      call MyRadHOwfn2(nc,lc,r,1.d0,valc)
      call MyRadHOwfn2(nd,ld,r,1.d0,vald)
      array(i) = vala*valb*valc*vald*r*r
   enddo


   call dquart_int(array,npts,ddr,dans)

   ans = sngl(dans)
   return
   end subroutine basic_contact_integral
!=======================================================================
   subroutine test_integral(na,la,nb,lb,Npts,dr)
!
!  subroutine to check accuracy of integrals by checking orthonormality
!
!  INPUT:
!    na,nb  = nodal quantum  numbers
!    la,lb  = orb ang mom
!    Npts = # of points 
   integer na,nb  ! radial nodal quantum numbers
   integer la,lb  ! orbital ang momenta
   integer Npts
   real(kind = 4) :: dr

!........... INTERMEDIATE

   real(kind =8), allocatable :: array(:)
   real(kind =8) :: dans
   real(kind = 8) :: ddr
   real(kind = 8) :: r
   real(kind =8)  :: vala,valb
   real(kind = 8), parameter :: dtol = 0.00001
   integer i

   if(la /= lb) return
   if(.not.allocated(array)) allocate(array(0:Npts))
   ddr = dble(dr)
   array(0) = 0.d0
   do i = 1,npts
      r = i*ddr
      call MyRadHOwfn2(na,la,r,1.d0,vala)
      call MyRadHOwfn2(nb,lb,r,1.d0,valb)

      array(i) = vala*valb*r*r
   enddo


   call dquart_int(array,npts,ddr,dans)

   if( na == nb)then
      if( abs(  dans - 1.d0) > dtol)then
         print*,' problem with normalization ',la,na,dans
         stop
      endif
   else
      if( abs(  dans) > dtol)then
         print*,' problem with orthogonality ',la,na,nb,dans
         stop
      endif

   endif

   return

   end subroutine test_integral

!=======================================================================


      subroutine dquart_int(f,n,h,ans)

! quartic integration -- Boole's rule (cf. Koonin and Meredith eq. 1.13b0

      implicit none

      integer n,i,steps,k
      real(kind=8) :: f(0:n),ans,h

      ans=0.

      if (mod(n,4).ne.0)then
          stop 999
      endif
      steps=n/4
      do i=1,steps
         k=i*4
         ans=ans+7.*(f(k-4)+f(k))+32.*(f(k-3)+f(k-1))   +12*f(k-2)
      enddo
      ans=2.*h*ans/45.
      return
      end 

!========================================================================
!  ALTERNATE FORMULATION via relative coordinates and Talmi-Moshinsky brackets
!========================================================================

!  subroutine rel_to_lab  converts contact interaction from relative frame to lab
!
!  call TMB
!
      subroutine rel_to_lab(na,la,nb,lb,nc,lc,nd,ld,L,ans)

      use sporbit 
      implicit none
      integer na,nb,nc,nd
      integer la,lb,lc,ld,L
      real ans
!................ INTERMEDIATE

      integer Ea,Eb,Ec,Ed
      integer Eab, Ecd
      integer Ecm
      integer Ncm, minNcm,maxNcm
      integer Ei,Ef,ni,nf
      real*8 tmb,tmb0  ! function Talmi-Moshinsky bracket
      real hopsi0  ! function to get psi(0)
      real tmbi,tmbf

      ans = 0.0
!................ error trap--only for contact interactions

      if( ((la + lb +L)/2)*2 /= la + lb +L)return
      if( ((lc + ld +L)/2)*2 /= lc + ld +L)return

!..........................
      Ea = 2*na + la
      Eb = 2*nb + lb         
      Ec = 2*nc + lc
      Ed = 2*nd + ld

      minNcm = 0
      maxNcm = min( na + nb +(la + lb -L)/2, nc + nd +(lc + ld -L)/2)
      ans = 0.0
      do Ncm = minNcm, maxNcm
         Ecm = 2* Ncm + L
         ef = ea+eb -Ecm 
         ei = ec+ed-Ecm 
         if(ef < 0 .or. ei < 0) then   
             print*,' woopps ',ef,ei
             print*,ea,eb,ec,ed,Ecm
             stop
         endif
         nf = ef /2
         ni = ei /2
!         if(ni > 3 .or. nf > 3) cycle
         tmbf = tmb0(Ecm, L, Ef, 0, ea,la,eb,lb,L,1.d0)
         tmbi = tmb0(Ecm, L, Ei, 0, ec,lc,ed,ld,L,1.d0)
!         write(87,*)nf,ni,tmbf,tmbi

         ans = ans + tmbf*tmbi*hopsi0(nf)*hopsi0(ni)
      enddo ! Ecm
      ans = ans *(0.5)**1.5

      return
      end subroutine rel_to_lab
!========================================================================

      function hopsi0(n)

      implicit none
      integer n
      real hopsi0
      real lndfact,lnfact

      hopsi0 = exp( 0.5 * ( lndfact(2*n+1)- log(2.)*n -lnfact(n)))
      hopsi0 = hopsi0* 0.42377721

      end function hopsi0
!===================================================================================
! ln factorial
!      function lnfact(n)
!      implicit none
!      integer n,i
!      real lnfact
!      lnfact = 0.0
!      if(n <2)return
!      do i = 2,n
!        lnfact = lnfact + log(float(i))
!      enddo
!      return
!      end function lnfact
!===================================================================================
! log double factorial
      function lndfact(n)
      implicit none
      integer n,i
      real lndfact
      lndfact = 0.0
      if(n <2)return
      do i = n,1, -2
        lndfact = lndfact + log(float(i))
      enddo
      return
      end function lndfact
!======================================================
