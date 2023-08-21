CCC here are a set of programs to calculate 3j,6j,9j, functions,
C also gamma


      function xninej(ej1,ej2,ej3,ej7,ej8,ej9)
c     calculates the special 9j symbol (j1,j2,j3,j4=.5,j5=.5,j6=1.,j7,
c     j8,j9) as defined by edmonds
      xninej=0.
      t1=tria(ej1,ej2,ej3)
      t2=tria(ej3,1.,ej9)
      t3=tria(ej1,.5,ej7)
      t4=tria(ej2,.5,ej8)
      t5=tria(ej7,ej8,ej9)
      j1=ej1*2.
      j2=ej2*2.
      j3=ej3*2.
      j7=ej7*2.
      j8=ej8*2.
      j9=ej9*2.
      if ((t1.eq.1.).and.(t2.eq.1.).and.(t3.eq.1.).and.(t4.eq.1.).and.
     1 (t5.eq.1.).and.(2*((j1+j2+j3)/2).eq.(j1+j2+j3)).and.(2*((j3+j9)
     2 /2).eq.(j3+j9)).and.(2*((j1+1+j7)/2).eq.(j1+1+j7)).and.(2*((j2+
     3 1+j8)/2).eq.(j2+1+j8)).and.(2*((j7+j8+j9)/2).eq.(j7+j8+j9)))
     4 go to 1
      return
    1 if (ej3.eq.ej9) go to 2
      el=(ej3+ej9)/2.
      xninej=(-1.)**((j3+j9)/2)*sixj(ej7,ej8,ej9,.5,el,ej2)*sixj(ej2,
     1 ej1,ej3,.5,el,ej7)/3./sixj(ej9,ej3,1.,.5,.5,el)
      return
    2 el=ej3+.5
      xninej=((-1.)**(j3+1)*sixj(ej7,ej2,el,.5,ej9,ej8)*sixj(ej2,ej7,
     1 el,.5,ej3,ej1)-(-1.)**((j1+j8-j3-1)/2)/2./(2.*ej9+1.)*sixj(ej9,
     2 ej2,ej1,.5,ej7,ej8))/3./sixj(ej9,ej3,1.,.5,.5,el)
      return
      end
 
c     ***************************************************
 
      function tria(a,b,c)
      tria=0.
      if (((a+b-c).lt.0.).or.((a+c-b).lt.0.).or.((b+c-a).lt.0.)) return
      tria=1.
      return
      end
 
c     ***************************************************
 
      function sixj(ej1,ej2,ej3,ej4,ej5,ej6)
      dimension fi(111),sf(111)
c     calculates 6j (j.lt.10) symbols defined by edmonds
      data jmin/2/
      sixj=0.
      if ((ej1.lt.0.).or.(ej2.lt.0.).or.(ej3.lt.0.).or.(ej4.lt.0.).or.
     1 (ej5.lt.0.).or.(ej6.lt.0.)) go to 5
      t1=tria(ej1,ej2,ej3)
      t2=tria(ej3,ej4,ej5)
      t3=tria(ej2,ej4,ej6)
      t4=tria(ej1,ej5,ej6)
      j1=ej1*2.01
      j2=ej2*2.01
      j3=ej3*2.01
      j4=ej4*2.01
      j5=ej5*2.01
      j6=ej6*2.01
      if ((t1.eq.1.).and.(t2.eq.1.).and.(t3.eq.1.).and.(t4.eq.1.).and.
     1 ((2*((j1+j2+j3)/2)).eq.(j1+j2+j3)).and.((2*((j3+j4+j5)/2)).eq.
     2 (j3+j4+j5)).and.((2*((j2+j4+j6)/2)).eq.(j2+j4+j6)).and.
     3 ((2*((j1+j5+j6)/2)).eq.(j1+j5+j6))) go to 1
      return
    1 jmax=ej1+ej2+ej3+ej4+ej5+ej6+1
c     the following defines fi(2j+1)=1./fact(j) and sf(2j+1)=
c     sqrt(fact(j)) for j=2,...,jmax
      fi(1)=1.
      fi(3)=1.
      sf(1)=1.
      sf(3)=1.
      if (jmin.gt.jmax) go to 3
      do 2 i=jmin,jmax
      f=fi(2*i-1)/i
      fi(2*i+1)=f
    2 sf(2*i+1)= sqrt(1./f)
      jmin=jmax+1
 
    3 if (ej3.eq.(ej4+ej5)) go to 7
c     we remember to cast our 6j symbols in the main program in this
c     form whenever possible
      jzmax=min0(j1+j2+j4+j5,j2+j3+j5+j6,j1+j3+j4+j6) +1
      jzmin=max0(j1+j2+j3,j1+j5+j6,j2+j4+j6,j3+j4+j5) +1
      k11=j1+j2+1-j3
      k12=j1+j3-j2+1
      k13=j2+j3+1-j1
      k14=j1+j2+j3+3
      k21=j1+j5-j6+1
      k22=j1-j5+j6+1
      k23=j5+j6-j1+1
      k24=j1+j5+j6+3
      k31=j4+j2-j6+1
      k32=j4-j2+j6+1
      k33=j2+j6-j4+1
      k34=j4+j2+j6+3
      k41=j4+j5-j3+1
      k42=j4-j5+j3+1
      k43=j5+j3-j4+1
      k44=j4+j5+j3+3
      a=sf(k11)*sf(k12)*(sf(k13)/sf(k14))
      b=sf(k21)*sf(k22)*(sf(k23)/sf(k24))
      c=sf(k31)*sf(k32)*(sf(k33)/sf(k34))
      d=sf(k41)*sf(k42)*(sf(k43)/sf(k44))
      i=jzmin-1
      if(i.ne.(4*(i/4))) a=-a
      do 6 jz=jzmin,jzmax,2
      k11=jz-j1-j2-j3
      k12=jz-j1-j5-j6
      k13=jz-j4-j2-j6
      k14=jz-j4-j5-j3
      k15=j1+j2+j4+j5-jz+2
      k16=j2+j3+j5+j6-jz+2
      k17=j1+j3+j4+j6-jz+2
      k18=jz+2
      sixj=sixj+(a*fi(k11)*fi(k12))*(b*fi(k13)*fi(k14))*
     1 (c*fi(k15)*d*fi(k16))*(fi(k17)/fi(k18))
    6 a=-a
c     we have used equation 6.3.7 in edmonds
      return
    7 k11=j3+j1-j2+1
      k12=j3+j2-j1+1
      k13=j1+j2+j3+3
      k14=j1+j2-j3+1
      k21=j5+j1-j6+1
      k22=j5+j6-j1+1
      k23=j1+j5+j6+3
      k24=j1+j6-j5+1
      k31=j4+j2-j6+1
      k32=j4+j6-j2+1
      k33=j4+j2+j6+3
      k34=j2+j6-j4+1
      k41=j3+j4-j5+1
      k42=j3+j5-j4+1
      k43=j3+j4+j5+3
      e=sf(k11)*sf(k12)*(sf(k13)/sf(k14))
      f=sf(k21)*sf(k22)*(sf(k23)/sf(k24))
      g=sf(k31)*sf(k32)*(sf(k33)/sf(k34))
      h=sf(k41)*(sf(k42)/sf(k43))
      sixj=(e/f)*(h/g)*(-1.)**((j1+j2+j3)/2)
c     we have used equation 6.3.1 in edmonds
      return
    5 write (61,10) ej1,ej2,ej3,ej4,ej5,ej6
   10 format (5x,39hthese quantum numbers not allowed: j1= ,f5.1,3x,
     1 4hj2= ,f5.1,3x,4hj3= ,f5.1,3x,4hj4= ,f5.1,3x,4hj5= ,f5.1,3x,
     2 4hj6= ,f5.1)
      stop
      end
 
c     ***************************************************
 
      function tria1(a,b,c,d,e,f)
      tria1=0.
      if (((d*d).gt.(a*a)).or.((e*e).gt.(b*b)).or.((f*f).gt.(c*c)))
     1 return
      tria1=1.
      return
      end
 
c     ***************************************************
 
      function threej(ej1,ej2,ej3,em1,em2,em3)
c     calculates 3j (j.lt.19) symbols defined by edmonds
      dimension fi(111),sf(111)
      data jmin/2/
      threej=0.
      em=em1+em2+em3
      t=tria(ej1,ej2,ej3)
      t1=tria1(ej1,ej2,ej3,em1,em2,em3)
      if ((ej1.lt.0.).or.(ej2.lt.0.).or.(ej3.lt.0.)) go to 5
      j1=ej1*2.01
      j2=ej2*2.01
      j3=ej3*2.01
      m1=em1*2.01
      m2=em2*2.01
      m3=em3*2.01
      if (((2*((j1-m1)/2)).ne.(j1-m1)).or.((2*((j2-m2)/2)).ne.(j2-m2))
     1 .or.((2*((j3-m3)/2)).ne.(j3-m3))) go to 5
      if ((em.eq.0.).and.(t.eq.1.).and.(t1.eq.1.)) go to 1
      return
    1 jmax=ej1+ej2+ej3+1
c     the following defines fi(2j+1)=1./fact(j) and sf(2j+1)=
c     sqrt(fact(j)) for j=2,...,jmax
      fi(1)=1.
      fi(3)=1.
      sf(1)=1.
      sf(3)=1.
      if (jmin.gt.jmax) go to 3
      do 2 i=jmin,jmax
      f=fi(2*i-1)/i
      fi(2*i+1)=f
    2 sf(2*i+1)= sqrt(1./f)
      jmin=jmax+1
    3 if ((em1.eq.0.).and.(em2.eq.0.)) go to 6
      jzmin=max0(0,j2-j3-m1,j1-j3+m2) +1
      jzmax=min0(j1+j2-j3,j1-m1,j2+m2) +1
      k11=j1+m1+1
      k12=j1-m1+1
      k21=j2+m2+1
      k22=j2-m2+1
      k31=j3+m3+1
      k32=j3-m3+1
      k41=j1-j2+j3+1
      k42=j2-j1+j3+1
      k43=j1+j2-j3+1
      k44=j1+j2+j3+3
      a=sf(k11)*sf(k12)
      b=sf(k21)*sf(k22)
      c=sf(k31)*sf(k32)
      d=sf(k41)*sf(k42)*(sf(k43)/sf(k44))
      i=jzmin+j1-j2-m3-1
      if (i.ne.4*(i/4)) a=-a
      do 4 jz=jzmin,jzmax,2
      k11=j1-m1-jz+2
      k12=j2+m2-jz+2
      k13=j1+j2-j3-jz+2
      k14=j3-j2+m1+jz
      k15=j3-j1-m2+jz
      threej=threej+(a*fi(k11)*b*fi(k12))*(c*fi(k13))*((d*fi(jz)*
     1 fi(k14))*fi(k15))
    4 a=-a
c     we have used equations 3.6.11 and 3.7.3 in edmonds
      return
    6 threej=0.
      j=j1+j2+j3
      if (((j/4)*4).eq.j) go to 7
      return
    7 k11=j-2*j1+1
      k12=j/2-j1+1
      k13=j-2*j2+1
      k14=j/2-j2+1
      k15=j-2*j3+1
      k16=j/2-j3+1
      k17=j+3
      k18=j/2+1
      threej=(-1.)**(j/4)*(sf(k11)*fi(k12))*(sf(k13)*fi(k14))*
     1 (sf(k15)*fi(k16))*(1./sf(k17)/fi(k18))
c      we have used equation 3.7.17 in edmonds
      return
    5 write (61,10) ej1,ej2,ej3,em1,em2,em3
   10 format (5x,39hthese quantum numbers not allowed: j1= ,f5.1,3x,
     1 4hj2= ,f5.1,3x,4hj3= ,f5.1,3x,4hm1= ,f5.1,3x,4hm2= ,f5.1,3x,
     2 4hm3= ,f5.1)
      stop
      end
 
c     ***************************************************
 
 
c     ***************************************************
 
      function x9j(ej11,ej12,ej13,ej21,ej22,ej23,ej31,ej32,ej33)
c     we use formula (6.4.3) of edmonds
      if ((ej11.lt.0.).or.(ej12.lt.0.).or.(ej13.lt.0.).or.(ej21.lt.0.)
     1 .or.(ej22.lt.0.).or.(ej23.lt.0.).or.(ej31.lt.0.).or.
     2 (ej32.lt.0.).or.(ej33.lt.0.)) go to 30
      if (tria(ej11,ej12,ej13)+tria(ej21,ej22,ej23)+tria(ej31,ej32,ej33)
     1 +tria(ej11,ej21,ej31)+tria(ej12,ej22,ej32)+tria(ej13,ej23,ej33)
     2 -6.) 1,2,1
    2 j11=ej11*2.
      j12=ej12*2.
      j13=ej13*2.
      j21=ej21*2.
      j22=ej22*2.
      j23=ej23*2.
      j31=ej31*2.
      j32=ej32*2.
      j33=ej33*2.
      if ((2*((j11+j12+j13)/2).eq.j11+j12+j13).and.
     1 (2*((j21+j22+j23)/2).eq.j21+j22+j23).and.(2*((j31+j32+j33)/2)
     2 .eq.j31+j32+j33).and.(2*((j11+j21+j31)/2).eq.j11+j21+j31).and.
     3 (2*((j12+j22+j32)/2).eq.j12+j22+j32).and.(2*((j13+j23+j33)/2)
     4 .eq.j13+j23+j33)) go to 3
    1 x9j=0.
      return
    3 xkmin=amax1(abs(ej11-ej33),abs(ej32-ej21),abs(ej12-ej23))
      xkmax=amin1(ej11+ej33,ej32+ej21,ej12+ej23)
      val=0.
      z=xkmin
   10 iph=2.*z
      val=val+sixj(ej11,ej21,ej31,ej32,ej33,z)*sixj(ej12,ej22,ej32,
     1 ej21,z,ej23)*sixj(ej13,ej23,ej33,z,ej11,ej12)*(2.*z+1.)*(-1.)**
     2 iph
      z=z+1.
      if (z-xkmax) 10,10,20
   20 x9j=val
      return
   30 write (61,31)
   31 format(2x,14hstopped in x9j)
      stop
      end
 
c     ***************************************************
 
      function gamma(y)
      double precision c
      dimension c(14)
      data c/                          0.99999 99999 99990 44d0,
     1   0.42278 43351 02334 79d0,     0.41184 03301 66781 29d0,
     2   0.08157 69261 24155 46d0,     0.07424 89154 19444 74d0,
     3  -0.00026 61865 94953 06d0,     0.01114 97143 35778 93d0,
     4  -0.00283 64625 30372 82d0,     0.00206 10918 50225 54d0,
     5  -0.00083 75646 85135 17d0,     0.00037 53650 52263 07d0,
     6  -0.00012 14173 48706 32d0,     0.00002 79832 88993 83d0,
     7  -0.00000 30301 90810 28d0/
      n=0
      f=1.
    1 if (y+n-2.) 2,5,3
    2 f=f*(y+n)
      n=n+1
      go to 1
    3 if (y+n-3.) 6,4,4
    4 n=n-1
      f=f*(y+n)
      go to 3
    5 gamma=1.
      go to 8
    6 gamma=c(14)
      z=y+n-2
      do 7 ii=1,13
      i=13-ii+1
    7 gamma=gamma*z+c(i)
    8 if (n) 11,12,9
    9 if (abs(f).lt.1.e-30) go to 10
      gamma=gamma/f
      go to 12
   10 gamma=1.e30
      go to 12
   11 gamma=gamma*f
   12 return
      end
 
c     ***************************************************
 
