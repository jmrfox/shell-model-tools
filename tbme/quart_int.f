
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine quart_int(f,n,h,ans)

! quartic integration -- Bode's rule (cf. Koonin and Meredith eq. 1.13b0

      implicit none

      integer n,i,steps,k
      real f(0:n),ans,h

      ans=0.

      if (mod(n,4).ne.0)then
          stop 999
      endif
      steps=n/4
      do i=1,steps
         k=i*4
         ans=ans+7.*(f(k-4)+f(k))+32.*(f(k-3)+f(k-1))
     1                  +12*f(k-2)
      enddo
      ans=2.*h*ans/45.
      return
      end 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
