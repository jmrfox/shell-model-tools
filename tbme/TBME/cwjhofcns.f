CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function dlnfactor(n)
C
C  COMPUTES ln n! in double precision
C
      implicit none
      integer :: n
      real(kind=8) :: dlnfactor
      integer i

      dlnfactor = 0.d0

      if(n == 0 .or. n==1)return
      do i = 2,n
        dlnfactor = dlnfactor + dlog(dfloat(i))
      enddo ! i
      return
      end function dlnfactor

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dlndfactor(n)
C
C  COMPUTES ln n!! in double precision
C
      implicit none
      integer :: n
      real(kind=8) :: dlndfactor
      integer i

      dlndfactor = 0.d0

      if(n == 0 .or. n==1)return
      do i = n,1,-2
        dlndfactor = dlndfactor + dlog(dfloat(i))
      enddo ! i
      return
      end function dlndfactor

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine dhalfasslaguerre(n,l,xd,ald)
C
C  computes special associated laguerre polynomial
C  ald = L_n^{l+1/2}(xd)
C  input xld and output ald are double precision
C

      implicit none
      integer n,l
      integer m
      real(kind = 8) :: xd,ald
      real(kind = 8) :: xdm
      real(kind = 8) :: dln2, darg
      real(kind = 8) :: dlnfactor,dlndfactor
      integer iphase

      dln2 = dlog(2.d0)
      xdm = 1.d0
      if(n == 0)then
        ald = xdm
        return
      endif
      ald = 0.d0
      iphase = 1
      do m = 0,n
        darg = dlndfactor(2*n+2*l+1) - dlndfactor(2*m+2*l+1)
     &  - dlnfactor(m) -dlnfactor(n-m) - dln2*(n-m)
        ald = ald + iphase*xdm*dexp(darg)
        xdm = xdm*xd
        iphase = - iphase
      enddo
      return
      end subroutine dhalfasslaguerre

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function dhonorm(n,l,bd)  
C
C  required normalization for spherical radial wfns for 
C  h.o.
C  here n = nodal quantum number
C  l      = orbital angular momentum
C  bd     = oscillator length parameter (in double precision)
C
C  output:
C    dhonorm  = norm (in double precision)
C
      implicit none
      integer n,l
      real(kind = 8) :: bd
      real(kind = 8) :: dhonorm
      real(kind = 8)  :: dln2, dlnpi
      real(kind = 8) :: dlnfactor,dlndfactor

      dln2 = dlog(2.d0)
      dlnpi = dlog(2.d0*dasin(1.d0))

      dhonorm = dlnfactor(n) + (n+l+2)*dln2
     & -0.5d0*dlnpi - dlndfactor(2*n+2*l+1) - 3*dlog(bd)
      dhonorm = dexp(0.5d0 * dhonorm)
      return
      end function dhonorm
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine MyRadHOwfn1(n,l,r,b,val)
C
C  single-precision  computes R_nl(r)  (not u_nl !)
C  n = nodal quantum #
C  l = orbital ang mom
C  b = h.o. length parameters
C
C  OUTPUT:
C    val = value of radial wavefunction R_nl (single-precision)
C
C  SUBROUTINES CALLED:
C     dhalfasslaguerre  :: computes value of special associated laguerre
C                          specifically L_n^{l+1/2}
C  FUNCTIONS CALLED:
C     dhonorm :: (double-precision) normalization of h.o. wfn
C


      real(kind = 4 ) :: b,val,r
      real(kind = 8 ) :: bd,vald,x

      integer :: n,l

      real(kind = 8) :: dhonorm

      bd = dble(b)
      x = dble(r)/bd
      call dhalfasslaguerre(n,l,x*x,vald)
      vald = vald*x**l * dexp(-0.5d0*x*x)
      vald = vald*dhonorm(n,l,bd)
      val = sngl(vald)
      return
      end subroutine MyRadHOwfn1


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine MyRadHOwfn2(n,l,rd,bd,vald)
C
C  double-precision  computes R_nl(r)  (not u_nl !)
C  n = nodal quantum #
C  l = orbital ang mom
C  b = h.o. length parameters
C
C  OUTPUT:
C    vald = value of radial wavefunction R_nl (double-precision)
C
C  SUBROUTINES CALLED:
C     dhalfasslaguerre  :: computes value of special associated laguerre
C                          specifically L_n^{l+1/2}
C  FUNCTIONS CALLED:
C     dhonorm :: (double-precision) normalization of h.o. wfn
C


C      real(kind =  ) :: b,val,r
      real(kind = 8 ) :: rd,bd,vald,x

      integer :: n,l

      real(kind = 8) :: dhonorm

C      bd = dble(b)
      x = rd/bd
      call dhalfasslaguerre(n,l,x*x,vald)
      vald = vald*x**l * dexp(-0.5d0*x*x)
      vald = vald*dhonorm(n,l,bd)
      return
      end subroutine MyRadHOwfn2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc

      subroutine RARadHOwf(n,l,bd,npts,xdl,a0,a1,a2,wf)
C
C  computes an array of radial wfns recursively
C  everything in double precision
C
C  INPUT:
C
C  n:   : current target
C  l:     orb ang moment
C bd:     h.o. length parameter
C npts:   # of points in array, starting at 0; declared dimension of arrays
C xdl:    array of x-values
C
C  OUTPUT:
C  wf(0:npts): = radial h.o. wfn for n,l
C  a0(0:npts), a1(0:npts), a2(0:npts): dummy arrays for recursively computing Laguerre
C

      implicit none
      integer :: nsize
      integer :: n,l
      real (kind = 8) :: bd
      integer npts
      real (kind = 8) :: xdl(0:npts)
      real(kind = 8) :: a0(0:npts),a1(0:npts),a2(0:npts)
      real(kind =8)  :: wf(0:npts)

      integer i
      real(kind = 8) :: dhonorm, dfact
      real(kind =8) :: x

      select case (n)

      case (0)
         a0(:) = 1.d0
      case (1)
         a1(:) = 1.d0
         do i = 0,npts
           a0(i) = - xdl(i)*xdl(i) + l + 1.5d0
         enddo

      case default
         a2(:) = a1(:)
         a1(:) = a0(:)
         dfact = 1.d0/dfloat(n)
         do i = 0,npts
           x = xdl(i)
           a0(i)=dfact*((2*n+l-0.5d0 -x*x)*a1(i)-(n+l-0.5d0)*a2(i))
         enddo
      end select
C------------ now convert to radial wfn

      dfact = dhonorm(n,l,bd)

      do i = 0,npts
         x = xdl(i)
         
         wf(i) = dfact*x**l*dexp(-0.5d0*x*x)*a0(i)
      enddo ! i

      return
      end subroutine RARadHOwf

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc