      function jpi(jt,par)
      implicit none
      character*5 jpi
      integer par
      real*8:: jt
      CHARACTER LNAME(0:14),PARITY(3)
      DATA LNAME / 's','p','d','f','g','h','i','j','k','l','m',
     &                'n','o','p','q' /
      DATA PARITY / '-','?','+' / 
      Write(jpi,'(1f3.1,a1)')jt,parity(par+2)
      end function


	subroutine writemat(mat,ndim,kout)
c writes square real matrix ndim x ndim to file
	implicit none
	integer ndim,kout,n,m	
	real*8:: mat(ndim,ndim)
	do n=1,ndim
	  write(kout,'(100f10.5)') (mat(n,m),m=1,ndim)
        enddo
	end subroutine
	


c *** writes j/parity in the form:
c 1/2+  -> 12p
c 1+    -> 1p 
c 1/2-  -> 12n
       function jname(j,par)
       implicit none
       integer l, par
       real*8 j       
       character*3 jname
       character*1 parchar
       character LNAME(0:14),PARITY(3)
       DATA LNAME / 's','p','d','f','g','h','i','j','k','l','m',
     &                'n','o','p','q' /
       DATA PARITY / '-','?','+' / 
       if (par.gt.0) then
         write(parchar,'(a1)') 'p'
       else
         write(parchar,'(a1)') 'n'
       endif

       if ((-1)**(2*j).lt.0) then ! j integer
         write(jname,'(i1,a1)') nint(j),parchar
       else
         write(jname,'(i1,a1,a1)') nint(2*j),'2',parchar
       endif
	end function


	function chname(jc,pcore,l,j)
	implicit none
       character*10 chname
	integer l,pcore
	real*8:: jc,j
       CHARACTER LNAME(0:14),PARITY(3)
       DATA LNAME / 's','p','d','f','g','h','i','j','k','l','m',
     &                'n','o','p','q' /
       DATA PARITY / '-','?','+' / 
       Write(chname,'(i1, a1,"x", a1,1f3.1)')
     &  jc,parity(pcore+1),trim(lname(l)),j
	end function


!       subroutine realloc_real(rv,nin,nout)
!       implicit none
!       real*8 rv(nin)
!       real*8 temp(nout)
!       temp(1:nout)=rv(1:nout)
!       deallocate(rv)
!       allocate(rv(nout))
!       rv(1:nout)=temp(1:nout)
!       return
!       end subroutine
       





c *** ---------------------------------------------------------------------------
       function fail3(x,y,z)
        implicit none
        logical fail3 
        real*8 x,y,z
        logical frac ! (x) = c
	fail3= frac(x+y+z) .or. x.gt.y+z .or. x.lt.abs(y-z)
       end function 

        function frac(x)
         implicit none
         real*8 x
         logical:: frac
         frac=abs(x-nint(x)).gt.1e-5 
        end function 



*************************************************************************
* Integral of y() (complex) 
* -----------------------------------------------------------------------
      subroutine simc(y,res,m,n,h,nramax)
      implicit real*8(a-h,o-z)
      complex*16:: y(nramax),res
      real*8 :: h,fr(nramax),fi(nramax)
      real*8 :: r1, r2
      fr(:)=real(y(:))
      fi(:)=dimag(y(:))
      r1=0d0; r2=0d0
      call sim(fr,r1,m,n,h,nramax)
      call sim(fi,r2,m,n,h,nramax)
      res=cmplx(r1,r2)
      return
      end subroutine


*************************************************************************
      subroutine sim(fa,res,m,n,h,nramax) !! REAL
*------------------------------------------------------------------------
*     subroutine does the integral of fa stored
*     h=step length
*     res=result
*     fa=array with "y" values
*     m,n= initial and final values of fa to consider on integration
*     nramax=dimension of fa
*
*     in the array of the same name using simpsons rule. the step length
*     is h and the integral is between the elements m and n of the arrays
*     only. resulting integral is placed in res.
*------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      integer :: m,n,nramax
      real*8 ::res,fa(nramax),dq(nramax),h

      do 90 i=m,n
      dq(i)=fa(i)
   90 continue
      rq1=dq(m+1)
      rq2=dq(m+2)
      i=m+3
   98 continue
      if(i.ge.n) go to 99
      rq1=rq1+dq(i)
      rq2=rq2+dq(i+1)
      i=i+2
      go to 98
   99 continue
      res=(1d0/3d0)*h*(dq(m)+4.d0*rq1+2.d0*rq2-dq(n))
      return
      end


c *** ----------------------------------------------------
      function cmp(r1,r2)
c *** ----------------------------------------------------
        real*8:: r1,r2,eps=1e-5
        logical:: cmp
        cmp=.false.
        if (abs(r1-r2).lt.eps) then
           cmp=.true.
        endif
      end function cmp


! Gamma function in double precision
!
      function dgamma(x)
      implicit real*8 (a - h, o - z)
      parameter (
     &    p0 = 0.999999999999999990d+00, 
     &    p1 = -0.422784335098466784d+00, 
     &    p2 = -0.233093736421782878d+00, 
     &    p3 = 0.191091101387638410d+00, 
     &    p4 = -0.024552490005641278d+00, 
     &    p5 = -0.017645244547851414d+00, 
     &    p6 = 0.008023273027855346d+00)
      parameter (
     &    p7 = -0.000804329819255744d+00, 
     &    p8 = -0.000360837876648255d+00, 
     &    p9 = 0.000145596568617526d+00, 
     &    p10 = -0.000017545539395205d+00, 
     &    p11 = -0.000002591225267689d+00, 
     &    p12 = 0.000001337767384067d+00, 
     &    p13 = -0.000000199542863674d+00)
      n = nint(x - 2)
      w = x - (n + 2)
      y = ((((((((((((p13 * w + p12) * w + p11) * w + p10) * 
     &    w + p9) * w + p8) * w + p7) * w + p6) * w + p5) * 
     &    w + p4) * w + p3) * w + p2) * w + p1) * w + p0
      if (n .gt. 0) then
          w = x - 1
          do k = 2, n
              w = w * (x - k)
          end do
      else
          w = 1
          do k = 0, -n - 1
              y = y * (x + k)
          end do
      end if
      dgamma = w / y
      end



************************************************************************
*     REAL 4-point lagrange interpolation routine.
*     interpolates thr FUNCTION value fival at point r from an
*     array of points stored in fdis(ndm). this array is assumed
*     to be defined such that the first element fdis(1) CONTAINS
*     the FUNCTION value at r=xv(1) and xv(2 .. ndm) are monotonically
*     increasing.
************************************************************************
      FUNCTION fival(r,xv,fdis,ndm,alpha)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 fdis(ndm),y1,y2,y3,y4
      DIMENSION xv(ndm)
      IF(r.GT.xv(ndm)) go to 9
      DO 5 k=1,ndm-2
 5    IF(r.LT.xv(k)) go to 6
      k=ndm-2
 6    nst=MAX(k-1,1)
      x1=xv(nst)
      x2=xv(nst+1)
      x3=xv(nst+2)
      x4=xv(nst+3)
      y1=fdis(nst+0)
      y2=fdis(nst+1)
      y3=fdis(nst+2)
      y4=fdis(nst+3)
      pii1=(x1-x2)*(x1-x3)*(x1-x4)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)
      xd1=r-x1
      xd2=r-x2
      xd3=r-x3
      xd4=r-x4
      pi1=xd2*xd3*xd4
      pi2=xd1*xd3*xd4
      pi3=xd1*xd2*xd4
      pi4=xd1*xd2*xd3
      fival=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4
      RETURN
 9    fival=fdis(ndm) * EXP(alpha*(xv(ndm)-r))
      RETURN
      END


!************************************************************************
!*     REAL 4-point lagrange interpolation routine.
!*     interpolates thr FUNCTION value fival at point r from an
!*     array of points stored in fdis(ndm). this array is assumed
!*     to be defined such that the first element fdis(1) CONTAINS
!*     the FUNCTION value at r=xv(1) and xv(2 .. ndm) are monotonically
!*     increasing.
!************************************************************************
      FUNCTION cfival(r,xv,fdis,ndm,alpha)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 cfival,fdis(ndm),y1,y2,y3,y4
      DIMENSION xv(ndm)
      IF(r.GT.xv(ndm)) go to 9
      DO 5 k=1,ndm-2
 5    IF(r.LT.xv(k)) go to 6
      k=ndm-2
 6    nst=MAX(k-1,1)
      x1=xv(nst)
      x2=xv(nst+1)
      x3=xv(nst+2)
      x4=xv(nst+3)
      y1=fdis(nst+0)
      y2=fdis(nst+1)
      y3=fdis(nst+2)
      y4=fdis(nst+3)
      pii1=(x1-x2)*(x1-x3)*(x1-x4)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)
      xd1=r-x1
      xd2=r-x2
      xd3=r-x3
      xd4=r-x4
      pi1=xd2*xd3*xd4
      pi2=xd1*xd3*xd4
      pi3=xd1*xd2*xd4
      pi4=xd1*xd2*xd3
      cfival=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4
      RETURN
 9    cfival=fdis(ndm) * EXP(alpha*(xv(ndm)-r))
      RETURN
      END




      real*8 function threej(ja,jb,jc,ma,mb,mc)
      implicit none
      real*8::ja,jb,jc,ma,mb,mc,cleb
      threej=(-1)**(jb+mc-ja)*
     &   cleb(ja,ma,jb,mb,jc,-mc)/sqrt(2d0*jc+1d0)
      end function

c *** --------------------------------------------
c     Clebsch-Gordan coefficient  <l',m',l,m|JM>
c *** --------------------------------------------
      real*8 function cleb(ria,rid,rib,rie,ric,rif)
      implicit real*8(a-h,o-z)
      real*8:: flog,ria,rid,rib,rie,ric,rif
!      common/clebma/faclog(500)
!      COMMON /PRAHA/ FLOG(100), GM(100), DG(25)
      ia=2.d0*(ria+.0001d0)
      ib=2.d0*(rib+.0001d0)
      ic=2.d0*(ric+.0001d0)
      id=int(sign(1.d0,rid)*2.d0*(abs(rid)+.0001d0))
      ie=int(sign(1.d0,rie)*2.d0*(abs(rie)+.0001d0))
      if=int(sign(1.d0,rif)*2.d0*(abs(rif)+.0001d0))
      wwww=-1.0d0
      cleb=0.0d0
      if(id+ie-if) 7000,105,7000
  105 k1=ia+ib+ic
      if((-1)**k1) 7000,107,107
  107 if(.not.((id.eq.0).and.(ie.eq.0))) go to 110
      k1=k1/2
      if((-1)**k1) 7000,110,110
  110 k1=ia+ib-ic
      k2=ic-iabs(ia-ib)
      k3=min0(k1,k2)
      if(k3) 7000,130,130
  130 if((-1)**(ib+ie)) 7000,7000,140
  140 if((-1)**(ic+if)) 7000,7000,150
  150 if(ia-iabs (id)) 7000,152,152
  152 if(ib-iabs (ie)) 7000,154,154
  154 if(ic-iabs (if)) 7000,160,160
  160 if(ia) 7000,175,165
  165 if(ib) 7000,175,170
  170 if(ic) 7000,180,250
  175 cleb=1.0d0
      go to 7000
  180 fb=float(ib+1)
      cleb=((wwww)**((ia-id)/2))/sqrt(fb)
      go to 7000
  250 fc2=ic+1
      iabcp=(ia+ib+ic)/2+1
      iabc=iabcp-ic
      icab=iabcp-ib
      ibca=iabcp-ia
      iapd=(ia+id)/2+1
      iamd=iapd-id
      ibpe=(ib+ie)/2+1
      ibme=ibpe-ie
      icpf=(ic+if)/2+1
      icmf=icpf-if
      vvv=0.5d0
      sqfclg=vvv*(log(fc2)-flog(iabcp+1)
     1      +flog(iabc)+flog(icab)+flog(ibca)
     2      +flog(iapd)+flog(iamd)+flog(ibpe)
     3      +flog(ibme)+flog(icpf)+flog(icmf))
      nzmic2=(ib-ic-id)/2
      nzmic3=(ia-ic+ie)/2
      nzmi= max0(0,nzmic2,nzmic3)+1
      nzmx= min0(iabc,iamd,ibpe)
      if(nzmx.lt.nzmi) go to 7000
      s1=(wwww)**(nzmi-1)
      do 400 nz=nzmi,nzmx
      nzm1=nz-1
      nzt1=iabc-nzm1
      nzt2=iamd-nzm1
      nzt3=ibpe-nzm1
      nzt4=nz-nzmic2
      nzt5=nz-nzmic3
      termlg=sqfclg-flog(nz)-flog(nzt1)-flog(nzt2)
     1           -flog(nzt3)-flog(nzt4)-flog(nzt5)
      ssterm=s1*exp (termlg)
      cleb=cleb+ssterm
  400 s1=-s1
 7000 return
      end


c     ------------------------------------------
c     Six-J coefficient {{j1,j2,j3},{j4,j5,j6}}
c     ------------------------------------------
      real*8 function sixj(r1,r2,r3,r4,r5,r6)
	implicit real*8(a-h,o-z)
        real*8:: phase,rac
	phase=(-1)**(r1+r2+r5+r4)
        sixj=phase*rac(r1,r2,r5,r4,r3,r6)
      end function	      

      real*8 function rac(ria,rib,ric,rid,rie,rif)
      implicit real*8(a-h,o-z)
!      COMMON /PRAHA/ FLOG(100), GM(100), DG(25)                    
c     -------------------------------------------------------------------
c     subroutine calculates the racah coefficient w(abcd;ef) defined
c     according to the convention of brink and satchler.
c     the arguments are real and are the actual values of the angular
c     momenta, ( i.e. they can take half integer values )
c     -------------------------------------------------------------------
!      common/clebma/faclog(500)
      dimension lt(6)

!      write(*,*)'In rac:'
!      write(*,*) "2",flog(2),exp(flog(2))
!      write(*,*)"10",flog(10),exp(flog(10))
! 	write(*,*)flog(10),exp(flog(10))

      rac=0.0d0
      ia=2.d0*(ria+.0001d0)
      ib=2.d0*(rib+.0001d0)
      ic=2.d0*(ric+.0001d0)
      id=2.d0*(rid+.0001d0)
      ie=2.d0*(rie+.0001d0)
      if=2.d0*(rif+.0001d0)
      k1=ia+ib-ie
      k2=ie-iabs (ia-ib)
      k3=ic+id-ie
      k4=ie-iabs (ic-id)
      k5=ia+ic-if
      k6=if-iabs (ia-ic)
      k7=ib+id-if
      k8=if-iabs(ib-id)
      k9= min0 (k1,k2,k3,k4,k5,k6,k7,k8)
      if(k9) 7000,20,20
   20 k2=k1-2*(k1/2)
      k4=k3-2*(k3/2)
      k6=k5-2*(k5/2)
      k8=k7-2*(k7/2)
      if(max0(k2,k4,k6,k8)) 7000,25,7000
   25 ltmin=min0(ia,ib,ic,id,ie,if)
      if(ltmin) 7000,30,150
   30 lt(1)=ia
      lt(2)=ib
      lt(3)=ic
      lt(4)=id
      lt(5)=ie
      lt(6)=if
      ltmin=lt(1)
      kmin=1
      do 40 n=2,6
      if(lt(n)-ltmin) 35,40,40
   35 ltmin=lt(n)
      kmin=n
   40 continue
      s1=1.0d0
      f1=ie
      f2=if
      go to (55,55,55,55,45,50),kmin
   45 f1=ia
      f2=ic
      s1=(-1.d0)**(k5/2)
      go to 55
   50 f1=ia
      f2=ib
      s1=(-1.d0)**(k1/2)
   55 rac=s1/dsqrt((f1+1.d0)*(f2+1.d0))
      go to 7000
  150 iabep=(ia+ib+ie)/2+1
      icdep=(ic+id+ie)/2+1
      iacfp=(ia+ic+if)/2+1
      ibdfp=(ib+id+if)/2+1
      iabe=iabep-ie
      ieab=iabep-ib
      ibea=iabep-ia
      icde=icdep-ie
      iecd=icdep-id
      idec=icdep-ic
      iacf=iacfp-if
      ifac=iacfp-ic
      icfa=iacfp-ia
      ibdf=ibdfp-if
      ifbd=ibdfp-id
      idfb=ibdfp-ib
      iabcd1=(ia+ib+ic+id+4)/2
      iefmad=(ie+if-ia-id)/2
      iefmbc=(ie+if-ib-ic)/2
      nzmax=min0(iabe,icde,iacf,ibdf)
      nzmi1=-iefmad
      nzmi2=-iefmbc
      nzmin=max0(0,nzmi1,nzmi2)+1
      if(nzmax.lt.nzmin) go to 7000
      sqlog=flog(iabe)+flog(ieab)+flog(ibea)+flog(icde)+
     &      flog(iecd)+flog(idec)+flog(iacf)+flog(ifac)+
     &      flog(icfa)+flog(ibdf)+flog(ifbd)+flog(idfb)-
     &      flog(iabep+1)-flog(icdep+1)-flog(iacfp+1)-flog(ibdfp+1)
      sqlog=0.5d0*sqlog
      do 200 nz=nzmin,nzmax
      nzm1=nz-1
      k1=iabcd1-nzm1
      k2=iabe-nzm1
      k3=icde-nzm1
      k4=iacf-nzm1
      k5=ibdf-nzm1
      k6=nz
      k7=iefmad+nz
      k8=iefmbc+nz
      sslog=sqlog+flog(k1)-flog(k2)-flog(k3)-flog(k4)
     &           -flog(k5)-flog(k6)-flog(k7)-flog(k8)
      ssterm=((-1.d0)**nzm1)*dexp(sslog)
      rac=rac+ssterm
!      write(96,'(8i3,20g12.6)') k1,k2,k3,k4,k5,k6,k7,k8,
!     &                      flog(k1),flog(k2),flog(k3),flog(k4),
!     &                      flog(k5),flog(k6),flog(k7),flog(k8) 
  200 continue
 7000 return
      end


      function flog(i)
        use factorials
        implicit none
        integer::i 
        real*8:: flog,logfac
c changed in v2.2 (use stored factorials !) 
!        flog=logfac(i-1)
        flog=dlfac(i-1)   ! already stored with factorialgen
!        if (abs(logfac(i-1)-dlfac(i-1)).gt.1e-5) then
!        write(*,*)'flog:',i,logfac(i-1),dlfac(i-1)
!        endif
      end function flog



c *** ---------------------------------------------
c Factorial LOG in recursive form 
c *** --------------------------------------------
      RECURSIVE FUNCTION logfac2(n) RESULT(lfact)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
!        INTEGER :: nfact
        REAL*8 :: lfact,xn
        
        xn=n
        IF(n > 0) THEN
           lfact = log(xn) + logfac2(n-1)
        ELSE
           lfact = 0d0
        END IF
!	write(*,*)'logfac: n,logfac=',n,lFACT
      END FUNCTION logfac2
           

c *** ---------------------------------------------
c Factorial LOG 
c *** --------------------------------------------
      real*8 function logfac(n) ! FL(N)                                       
      implicit real*8(a-h,o-z),integer*4(i-n) 
       fl=0
       if(n>1) then
       FN = 1.                                                           
       DO 10 I = 2,N                                                     
       FN = FN + 1.                                                      
   10  FL = FL +  LOG(FN)    
      endif 
      logfac=fl                                 
      END FUNCTION


c===================================================================
      subroutine sbesjh(x,lmax,xj,xjp,xh1,xh1p,ifail)
c ***                                                       i.j.thompson
c ***                                                       31 may 1985.
c ***  complex spherical bessel functions from l=0 to l=lmax
c ***    for x in the upper half plane ( im(x) > -3)
c ***
c ***    xj(l)   = j/l(x)          regular solution: xj(0)=sin(x)/x
c ***    xjp(l)  = d/dx j/l(x)
c ***    xh1(l)  = h(1)/l(x)       irregular hankel function:
c ***    xh1p(l) = d/dx h(1)/l(x)            xh1(0) = j0(x) + i. y0(x)
c ***                                               =(sin(x)-i.cos(x))/x
c ***                                               = -i.exp(i.x)/x
c ***  using complex cf1, and trigonometric forms for l=0 solutions.
c ***
      implicit complex*16 (a-h,o-z)
      parameter (limit=20000)
      dimension xj(0:lmax),xjp(0:lmax),xh1(0:lmax),xh1p(0:lmax)
      real*8 zero,one,accur,tm30,absc
      data zero,one/ 0.0d0,1.0d0 /, accur /1.0d-12/, tm30 / 1d-30 /,
     #     ci / (0d0,1d0) /
      absc(w) = abs(real(w)) + abs(dimag(w))
      ifail= -1
      if(absc(x).lt.accur .or. dimag(x).lt.-3.0) go to 5
      xi = one/x
      w  = xi + xi
      pl = lmax*xi
      f = pl + xi
      b  = f + f + xi
      d  = zero
      c  = f
      do 1 l=1,limit
      d  = b - d
      c  = b - one/c
         if(absc(d).lt.tm30) d = tm30
         if(absc(c).lt.tm30) c = tm30
      d = one / d
      del= d * c
      f = f * del
      b = b + w
    1 if(absc(del-one).lt.accur) go to 2
        ifail = -2
        go to 5
c
    2 xj(lmax)   = tm30
      xjp(lmax)  = f * xj(lmax)
c
c *** downward recursion to l=0 (n.b.  coulomb functions)
c
      do 3 l = lmax-1,0,-1
      xj(l) = pl*xj(l+1) + xjp(l+1)
      xjp(l)= pl*xj(l)   - xj(l+1)
    3 pl = pl - xi
c *** calculate the l=0 bessel functions
      xj0  = xi * sin(x)
      xh1(0) = exp(ci*x) * xi * (-ci)
      xh1p(0)= xh1(0) * (ci - xi)
c
c *** rescale xj, xjp,  converting to spherical bessels.
c *** recur   xh1,xh1p             as spherical bessels.
c
        w = one/xj(0)
         pl = xi
      do 4  l = 0,lmax
      xj(l)  =  xj0*(w*xj(l))
      xjp(l) =  xj0*(w*xjp(l)) - xi*xj(l)
         if(l.eq.0) go to 4
      xh1(l) = (pl-xi) * xh1(l-1) - xh1p(l-1)
         pl = pl + xi
      xh1p(l)=- pl     * xh1(l)   + xh1(l-1)
   4  continue
      ifail = 0
      return
    5      write(6,10) ifail
   10      format( 'sbesjh : ifail = ',i4)
      return
      end

c*****************************************************************
      subroutine gauss(a,b,npoint,xri,wri)
c*****************************************************************
      implicit real*8(a-h,o-z)
      dimension xg(npoint),wg(npoint),xri(npoint),wri(npoint)
      if(npoint.le.0) return
      call setmgl(npoint,xg,wg)
      do 20 j=1,npoint
      xri(j) = (a+b)/2.0 + (b-a)/2.0*xg(j)
      wri(j) = (b-a)/2.0*wg(j)
   20 continue
      return
      end

c********************************************************************
      subroutine setmgl(n,points,weight)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     setmgl returns points and weights for n-point gauss-legendre
c     quadrature.
c     on the interval [-1,1] rather than the usual [0,infinity).
c     see abramowitz and stegun for the grubby details.
c     n      -- the number of points to use in the numerical quadrature.
c               it is also the degree of the needed legendre polynomial.
c               this is passed to the routine.
c     points -- the points of integration are returned here.
c     weight -- the weights for integration are returned here.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      implicit none
!	use globals, only: pi
      implicit real*8(a-h,o-z)
      dimension poin16(n)
      real*8, parameter :: pi=acos(-1d0)
c      parameter  (pi = 3.14159 26535 89793 23846 26433 83279 50288)
      dimension points(n), weight(n)
c
      if (n.le.0)  return

      m = ( n + 1 ) / 2
      e1 = n * ( n + 1 )
      do 1 i = 1, m
        t = ( 4*i - 1 ) * pi / ( 4*n + 2 )
        x0 = ( 1.0 - ( 1.0 - 1.0/n ) / ( 8.0*n*n ) ) * cos(t)
        pkm1 = 1.0
        pk = x0
        do 3 k = 2, n
          t1 = x0 * pk
          pkp1 = t1 - pkm1 - ( t1-pkm1 )/k + t1
          pkm1 = pk
          pk = pkp1
3       continue
        den = 1.0 - x0*x0
        d1 = n * ( pkm1 - x0*pk )
        dpn = d1 / den
        d2pn = ( 2.0*x0*dpn - e1*pk ) / den
        d3pn = ( 4.0*x0*d2pn + (2.0-e1)*dpn ) / den
        d4pn = ( 6.0*x0*d3pn + (6.0-e1)*d2pn ) / den
        u = pk / dpn
        v = d2pn / dpn
        h = -u * ( 1.0 + 0.5*u*(v+u*(v*v-u*d3pn/(3.0*dpn))))
        p = pk + h*(dpn+0.5*h*(d2pn+h/3.0*(d3pn+0.25*h*d4pn)))
        dp = dpn + h*(d2pn+0.5*h*(d3pn+h*d4pn/3.0))
        h = h - p / dp
        poin16(i) = x0 + h
        fx = d1 - h*e1*(pk+0.5*h*(dpn+h/3.0*
     1                   (d2pn+0.25*h*(d3pn+0.2*h*d4pn))))
        weight(i) = 2.0 * ( 1.0 - poin16(i)*poin16(i)) / (fx*fx)
1     continue
      if ( m + m .gt. n ) poin16(m) = 0.0
      do 10 i = n/2 + 1, n
        poin16(i) = poin16( n + 1 - i )
        weight(i) = weight( n + 1 - i )
        poin16( n + 1 - i ) = -poin16( n + 1 - i )
10    continue
c      do 20 i = 1, n
c        print 1000, i, points(i), weight(i)
c1000  format( 1x, i4, 2(5x, d22.15) )
c        temp = pi4 * ( points(i) + 1.0 )
c        points(i) = dtan(temp)
c        weight(i) = pi4 * weight(i) / ( dcos(temp) ** 2 )
c20    continue
        do 30 i=1,n
 30     points(i)=poin16(i)
      return
      end
c*
c* ==============================================================================
c*
!----------------------------------------------GAUSS-LAGUERRE HYPERRADIAL BASIS
	subroutine qualag(x,w,n,al)
! gauss-laguerre quadrature formulas
! al is an integer number
	use factorials
	implicit real*8 (a-h,o-z)
	dimension x(n),w(n)
	y=1.
	eps=1.e-13
! initial values for nodes
	h=.01d0/n
	f=2.d0*n+al+1.d0
	g=f+sqrt(f*f+.25d0-al*al)
	m=g/h+1
	j=0
	do 1 i=1,m
	z=h*i
	v=polag(n,al,z)
	c=v/y
	if(c.lt.0) then
	j=j+1
	x(j)=z
	end if
1	y=v
! exact nodes
	do 2 i=1,n
	xs=x(i)
	call rootlag(xs,y,eps,n,al)
2	x(i)=y
! the weights
	lal=al
	s=exp(dlfac(n+lal)-dlfac(n+1))/(n+1)
	do 3 i=1,n
	a=x(i)
	b=polag(n+1,al,a)
3	w(i)=s*x(i)/(b*b)
	return
	end
	real*8 function lagnorm(n,al)	
	use factorials
	implicit real*8 (a-h,o-z)
c norm. coeff. for laguerre polynomials
	lal=al
	lagnorm=sqrt(exp(dlfac(n)-dlfac(n+lal)))
	return
	end

      subroutine rootlag(xs,x,eps,n,al)
c solution of eq. f(x)=0 by the newton method. eps shows the accuracy
c of the solution (see the program). the program calls the subr. fun which
c calculates f(y) and its derivative f1(y). xs is the starting value.
      implicit real*8 (a-h,o-z)
      y=xs
 1    continue
      call fun(f,f1,y,al,n)
      x=y-f/f1
c      print 2, x
 2    format(d24.14)
      a=dabs(x-y)
      if(a.lt.eps) go to 3
      y=x
      go to 1
 3    continue
      return
      end	
	
	subroutine fun(f,f1,y,al,n)
	implicit real*8 (a-h,o-z)
	f=polag(n,al,y)
	f1=-polag(n-1,al+1.d0,y)
	return
	end
	
	
	real*8 function polag(n,al,x)
	implicit real*8 (a-h,o-z)
calculates laguerre polynomials	
	a=1.d0
	if(n.eq.0) then
	   polag=a
	   return
	end if
	b=al+1.d0-x
	if(n.eq.1) then
	   polag=b
	   return
	end if
	do 1 i=2,n
	   c=((2.d0*i+al-1.d0-x)*b-(i+al-1.d0)*a)/i
	   a=b
1	   b=c
	polag=c
	return
	end
	real*8 function ekin(n1,n2,k)
c without h^2/(2m*b^2) factor
	implicit real*8 (a-h,o-z)
	real*8 lagnorm
	ma=max0(n1,n2)
	mi=min0(n1,n2)
	de=0.
	if(ma.eq.mi) de=1.
	f=lagnorm(ma,5.d0)/lagnorm(mi,5.d0)
	ekin=f*(.5d0-.25d0*de+mi/6.d0+
     &k*(k+4.d0)*(5.d0*(ma-mi+1.d0)+ma+mi+1.d0)/120.d0)
	return
	end      	


!------------------------------------------------------- MISC ROUTINES
      subroutine factorialgen(n)
! n>0
! ln(i!) and ln((2i+1)!)
      use factorials
      implicit real*8 (a-h,o-z)
      dlfac(0)=0.
      fact(0)=1.
      dlfac2(0)=0
      dlfac2(1)=0
      if (n.lt.1) return
      do 1 i=1,n 
	a=i
	dlfac(i)=dlfac(i-1)+log(a)
        fact(i)=dexp(dlfac(i))
1     if (i.ge.2) dlfac2(i)=dlfac2(i-2)+log(a)

	dl2fac(0)=0.d0
	do 2 i=1,n
	a=2.d0*i+1.d0
2	dl2fac(i)=dl2fac(i-1)+log(a)
	continue
	return
	end

      function ddexp(x)
      implicit doubleprecision(a-h,o-z)
      if(abs(x).gt.80.d0) x = sign(80.d0,x)
      ddexp=exp(x)
      return
      end


      SUBROUTINE PLM(X,N,M,NA,PL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 PL(NA,M+1),L,X
      N1 = N+1
      M1 = M+1
      DO 10 J=1,M1
      DO 10 I=1,N1
10    PL(I,J)=0.
      PL(1,1) = 1.
      PL(2,1) = X
      SX = SQRT(1.-X*X)
      PL(2,2) = SX
	  FACT=1.
	  PMM=1.
	  DO 15 J=2,min(M1,N1)
	    mm = J-1
	    PMM = PMM*FACT*SX
		FACT=FACT+2.
		PL(J,J) = PMM
		if(J+1.le.N1) PL(J+1,J) = X*(2*mm+1.) * PL(J,J)
15      CONTINUE
	  
	  DO 20 J=1,M1
	   mm = J-1
	  DO 20 I=J+2,N1
	   ll = I-1
      PL(I,J)=((2.*ll-1.)*X*PL(I-1,J) - (ll+mm-1.)*PL(I-2,J))/(ll-mm)	  
20    CONTINUE
      RETURN
      END


*----------------------------------------------------------------------
      real*8 function plm_nr(il,imm,c,c2,s,s2)
      implicit real*8(a-h,o-z)
*-----------------------------------------------------------------------
* Associated Legendre functions for l<=8
*-----------------------------------------------------------------------
      plm_nr=1.d0
      if(il.eq.0) return
      im=abs(imm)
*-----------------------------------------------------------------------
      if(il.eq.1) then
       if(im.eq.0) then
        plm=c
       else if(im.eq.1)then
        plm=s
       endif
*-----------------------------------------------------------------------
      else if(il.eq.2) then
       if(im.eq.0) then
        plm=(3.d0*c2-1.d0)/2.d0
       else if(im.eq.1) then
        plm=3.d0*s*c
       else if(im.eq.2) then
        plm=3.d0*s2
       endif
*-----------------------------------------------------------------------
      else if(il.eq.3) then
       if(im.eq.0) then
        plm=(5.d0*c2-3.d0)*c/2.d0
       else if(im.eq.1) then
        plm=s*(15.d0*c2-3.d0)/2.d0
       else if(im.eq.2) then
        plm=15.d0*s2*c
       else if(im.eq.3) then
        plm=15.d0*s2*s
       endif
*-----------------------------------------------------------------------
      else if(il.eq.4) then
       if(im.eq.0) then
        plm=(35.d0*c2*c2-30.d0*c2+3.d0)/8.d0
       else if(im.eq.1) then
        plm=s*(35.d0*c2-15.d0)*c/2.d0
       else if(im.eq.2) then
        plm=s2*(105.d0*c2-15.d0)/2.d0
       else if(im.eq.3) then
        plm=105.d0*s*s2*c
       else if(im.eq.4) then
        plm=105.d0*s2*s2
       endif
*-----------------------------------------------------------------------
      else if(il.eq.5) then
       if(im.eq.0) then
        plm=(63.d0*c2*c2-70.d0*c2+15.d0)*c/8.d0
       else if(im.eq.1) then
        plm=s*(315.d0*c2*c2-210.d0*c2+15.d0)/8.d0
       else if(im.eq.2) then
        plm=s2*(315.d0*c2-105.d0)*c/2.d0
       else if(im.eq.3) then
        plm=s*s2*(945.d0*c2-105.d0)/2.d0
       else if(im.eq.4) then
        plm=945.d0*s2*s2*c
       else if(im.eq.5) then
        plm=945.d0*s2*s2*s
       endif
*-----------------------------------------------------------------------
      else if(il.eq.6) then
       if(im.eq.0) then
        plm=(231.d0*c2*c2*c2-315.d0*c2*c2+105.d0*c2-5.d0)/16.d0
       else if(im.eq.1) then
        plm=21.d0*s*c*(33.d0*c2*c2-30.d0*c2+5.d0)/8.d0
       else if(im.eq.2) then
        plm=(33.d0*c2*c2*c2-51.d0*c2*c2+19.d0*c2-1.d0)
        plm=-105.d0*plm/8.d0
       else if(im.eq.3) then
        plm=315.d0*s*s2*c*(11.d0*c2-3.d0)/2.d0
       else if(im.eq.4) then
        plm=(11.d0*c2*c2*c2-23.d0*c2*c2+13.d0*c2-1.d0)
        plm=945.d0*plm/2.d0
       else if(im.eq.5) then
        plm=10395.d0*s2*s2*s*c
       else if(im.eq.6) then
        plm=10395.d0*s2*s2*s2
       endif
*-----------------------------------------------------------------------
      else if(il.eq.7) then 
         if(im.eq.0) then 
            plm=(-35*z + 315*z**3 - 693*z**5 + 429*z**7)/16.
         else if(im.eq.1) then
            plm=(7*Sqrt(1 - c**2)*
     -    (-5 + 135*c**2 - 495*c**4 + 429*c**6))/16.
         else if(im.eq.2) then
            plm=(-63*(-1 + c**2)*(15*c - 110*c**3 + 143*c**5))/8.
         else if(im.eq.3) then
            plm= -(315*Sqrt(1 - c**2)*(-1 + c**2)*
     -    (3 - 66*c**2 + 143*c**4))/8.
         else if(im.eq.4) then
            plm=(3465*(-1 + c**2)**2*(-3*c + 13*c**3))/2.
         else if(im.eq.5) then
            plm= (10395*Sqrt(1 - c**2)*(-1 + c**2)**2*
     -    (-1 + 13*c**2))/2.
         else if(im.eq.6) then
            plm=-135135*c*(-1 + c**2)**3
         else if(im.eq.7) then
            plm=135135*(1 - c**2)**3.5
         endif
*-----------------------------------------------------------------------
      else if(il.eq.8) then
         if(im.eq.0) then 
            plm= (35 - 1260*c**2 + 6930*c**4 - 12012*c**6 + 
     -           6435*c**8)/128.
         else if(im.eq.1) then 
!            plm=(9*Sqrt(1 - c**2)*
!     -           (-35*c + 385*c**3 - 1001*c**5 + 715*c**7))/16.
            plm=(9*s*
     -           (-35*c + 385*c**3 - 1001*c**5 + 715*c**7))/16.

         else if(im.eq.2) then 
            plm=(-315*(-1 + c**2)*
     -           (-1 + 33*c**2 - 143*c**4 + 143*c**6))/16.
         else if(im.eq.3) then 
!            plm=-(3465*Sqrt(1 - c**2)*(-1 + c**2)*
!     -                 (3*c - 26*c**3 + 39*c**5))/8.
            plm=-(3465*s*(-1 + c**2)*
     -                 (3*c - 26*c**3 + 39*c**5))/8.
         else if(im.eq.4) then 
            plm=(10395*(-1 + c**2)**2*(1 - 26*c**2 + 65*c**4))/8.
         else if(im.eq.5) then
!            plm=(135135*Sqrt(1 - c**2)*(-1 + c**2)**2*
!     -           (-c + 5*c**3))/2.
            plm=(135135*s*(-1 + c**2)**2*
     -           (-c + 5*c**3))/2.
         else if(im.eq.6) then 
!            plm=(-135135*(-1 + c**2)**3*(-1 + 15*c**2))/2.
            plm=(135135*s2**3*(-1 + 15*c**2))/2.
         else if(im.eq.7) then
!            plm=2027025*c*(1 - c**2)**3.5
             plm=2027025*c*s2**3.5
         else if(im.eq.8) then
!            plm=2027025*(-1 + c**2)**4
            plm=-2027025*(s2)**4
         endif
      endif
      plm_nr=plm
      return
      end



c Another subroutine for PLM and its derivative
        SUBROUTINE LPMN(MM,M,N,X,PM,PD)
C
C       =====================================================
C       Purpose: Compute the associated Legendre functions 
C                Pmn(x) and their derivatives Pmn'(x)
C       Input :  x  --- Argument of Pmn(x)
C                m  --- Order of Pmn(x),  m = 0,1,2,...,n
C                n  --- Degree of Pmn(x), n = 0,1,2,...,N
C                mm --- Physical dimension of PM and PD
C       Output:  PM(m,n) --- Pmn(x)
C                PD(m,n) --- Pmn'(x)
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PM(0:MM,0:N),PD(0:MM,0:N)
        DO 10 I=0,N
        DO 10 J=0,M
           PM(J,I)=0.0D0
10         PD(J,I)=0.0D0
        PM(0,0)=1.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 I=1,N
              PM(0,I)=X**I
15            PD(0,I)=0.5D0*I*(I+1.0D0)*X**(I+1)
           DO 20 J=1,N
           DO 20 I=1,M
              IF (I.EQ.1) THEN
                 PD(I,J)=1.0D+300
              ELSE IF (I.EQ.2) THEN
                 PD(I,J)=-0.25D0*(J+2)*(J+1)*J*(J-1)*X**(J+1)
              ENDIF
20         CONTINUE
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XQ=DSQRT(LS*(1.0D0-X*X))
        XS=LS*(1.0D0-X*X)
        DO 30 I=1,M
30         PM(I,I)=-LS*(2.0D0*I-1.0D0)*XQ*PM(I-1,I-1)
        DO 35 I=0,M
35         PM(I,I+1)=(2.0D0*I+1.0D0)*X*PM(I,I)
        DO 40 I=0,M
        DO 40 J=I+2,N
           PM(I,J)=((2.0D0*J-1.0D0)*X*PM(I,J-1)-
     &             (I+J-1.0D0)*PM(I,J-2))/(J-I)
40      CONTINUE
        PD(0,0)=0.0D0
        DO 45 J=1,N
45         PD(0,J)=LS*J*(PM(0,J-1)-X*PM(0,J))/XS
        DO 50 I=1,M
        DO 50 J=I,N
           PD(I,J)=LS*I*X*PM(I,J)/XS+(J+I)
     &             *(J-I+1.0D0)/XQ*PM(I-1,J)
50      CONTINUE
        RETURN
        END


!
!  CALCULATE THE COEFFICIENT OF P(L,M)*E(I*M*PHI) IN Y(L,M)
!  AMoro: Uses |M|  
!
      FUNCTION YLMC2(L,M)
	use factorials
      IMPLICIT REAL*8(A-H,O-Z)
      integer,intent(IN):: L,M
      PHASE(I) = (-1)**I
      pi=acos(-1d0)
      LF1 = L + 1
      MA = ABS(M)
! changed by AMoro
!      R =  FACT(LF1-MA)-FACT(LF1+MA)
      R =  FLOG(LF1-MA)-FLOG(LF1+MA)
      R = SQRT((2*L+1)/(4*PI)*EXP(R))* PHASE(M)
      IF(M.LT.0) R = R * PHASE(MA)
      YLMC2 = R
      RETURN
      END




c-----------------------------------------------------------------------
c    Al parecer a L y M le corresponde Y(L,M) y PL(L+1,abs(M)+1)
c    en virtud de que P(L,-M)=factor(L,M)*P(L,M)
c------------------------------------------------------------------------

      FUNCTION YLMC(L,M)
	use factorials
        use globals
      IMPLICIT REAL*8(A-H,O-Z)
      integer,intent(IN):: L,M
      PHASE(I) = (-1)**I
C
C      CALCULATE THE COEFFICIENT OF P(L,M)*E(I*M*PHI) IN Y(L,M)
C
         LF1 = L + 1
         MA = ABS(M)
                             R =  FACT(LF1-MA)-FACT(LF1+MA)
      R = SQRT((2*L+1)/(4*PI)*EXP(R))
     &      * PHASE(M)
      IF(M.LT.0) R = R * PHASE(MA) 
      YLMC = R
      RETURN
      END

c----------------------------------Hay que acordarse de llamar a logfac3(L)

c----------------------------------¿que pasa con FAC(0)?

      SUBROUTINE LOGFAC3(L)
	use factorials
C  FACT(I)= LN((I-1)!)
C  THE VALUE OF L CORRESPONDS TO THE DIMENSION OF (FACT(
      FACT(1)=0.
      DO 100 J=2,L
  100 FACT(J)=FACT(J-1)+LOG(J-1D0)
      RETURN
      END




!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!	Subroutine to calculate rl wich is the integrals of r^l*f(r)g(r)     !
!       it also gives the norm: f(r)g(r) and the root mean square radii rms  !
!       h is the step, r0 the first point and nr the number of points        !
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
        subroutine normfun(f,g,h,r0,nr,l,norm,rms,rl)
        implicit none
        integer :: nr,l,i
        real*8:: f(nr),g(nr),r0,h
        real*8:: res1,res2,norm,rms,rl
        real*8:: faux(nr),r(nr)
!        write(*,*)h,r0,nr,l
        do i=1,nr
           r(i)=r0+h*dble(i-1)
        enddo
        faux(:)=r(:)**2*f(:)*g(:)
        call sim(faux,norm,1,nr,h,nr)
        faux(:)=r(:)**4*f(:)*g(:)
        call sim(faux,res1,1,nr,h,nr)
        rms=sqrt(res1)
        faux(:)=r(:)**(2+l)*f(:)*g(:)
        call sim(faux,rl,1,nr,h,nr)
        end




      SUBROUTINE GAUSS5(N,NR,A,M,SING,DET,EPS,SHOW)
C
C    SOLVE BY GAUSSIAN ELIMINATION SUM(J): A(I,J).P(J) = A(I,N+1)
C             WHERE P(J) IS LEFT IN MAT(J,N+1)
C
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(KO=99)
      COMPLEX*16 A(NR,N+M),DET,RA
      LOGICAL SING,SHOW
      SING  = .FALSE.
      NPM = N + M
      DO 201 I=1,N
201   IF(SHOW) WRITE(KO,402) (A(I,J)  ,J=1,NPM)
      DET = 0.
      DO 9 K = 1, N
      IF (abs(A(K,K)) .NE. 0.0 ) GO TO 5
         SING = .TRUE.
         WRITE(KO,3) K,DET/LOG(10.)
    3  FORMAT(//' THE MATRIX IS SINGULAR AT',I3,', Log10 determinant is
     &  ',2E16.8/)
      DO 401 I=1,N
401   IF(SHOW) WRITE(KO,402) (            A(I,J)  ,J=1,NPM)
C402   FORMAT( 1X,20F6.1/(1X,20F6.3))
402   FORMAT( 1X,1P,14E9.1/(1X,24E9.1))
         RETURN
    5 KP1 = K + 1
      DET = DET + LOG(A(K,K))
         RA = 1.0/A(K,K)
      DO 6 J = KP1, NPM
    6 A(K,J) = A(K,J) * RA
      A(K,K) = 1
      DO 9 I = 1, N
      IF (I .EQ. K  .OR. ABS (A(I,K)) .EQ. 0) GO TO 9
CDIR$ IVDEP
         DO 8 J = KP1, NPM
    8    A(I,J) = A(I,J) - A(I,K)*A(K,J)
         A(I,K) = 0
    9 CONTINUE
      IF(SHOW) WRITE(KO,15 ) DET/LOG(10.)
15    FORMAT(/' Log10 determinant is ',2F10.5)
      IF(SHOW) WRITE(KO,402) (A(I,N+1),I=1,N)
      RETURN
      END

c Real Kronecker function
      function kron(i,j)
      integer i,j
      real*8 kron
      kron=0d0
      if (i.eq.j) kron=1d0
      return
      end function



c *** From DCE
      double precision function dkron(n,m)
      implicit real*8 (a-h,o-z)
      if(n.eq.m) then
      dkron=1.d0
      else
      dkron=0.d0
      endif
      return
      end

       double precision function fleg(x,nn)
       implicit real*8 (a-h,o-z)
       p0=1.d0
       if(nn.eq.0) then
       fleg=p0
       return
       endif
       p1=x
       if(nn.eq.1) then
       fleg=p1
       return
       endif
       do n=1,nn-1
       pn1=(x*(2.d0*n+1.d0)*p1-n*p0)/dble(n+1)
       p0=p1
       p1=pn1
       enddo
       fleg=pn1
       return
       end



c *** ------------------------------------------------------------------
c *** Splines 
c ***
      subroutine spline2(n,x,y,b,c,d)
      implicit real*8(a-h,o-z)
      dimension x(n),y(n),b(n),c(n),d(n)
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
c
      nm1 = n-1
      if (n.lt.2) return
      if (n.lt.3) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1)=x(2)-x(1)
      c(2)=(y(2)-y(1))/d(1)
      do 10 i =2,nm1
      d(i)=x(i+1)-x(i)
      b(i)=2.d0*(d(i-1)+d(i))
      c(i+1)=(y(i+1)-y(i))/d(i)
      c(i)=c(i+1)-c(i)
10    continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1)=-d(1)
      b(n)=-d(n-1)
      c(1)=0.d0
      c(n)=0.d0
      if(n.eq.3) go to 15
      c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
      c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
      c(1)=c(1)*d(1)**2/(x(4)-x(1))
      c(n)=-c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
15    do 20 i=2,n
      t=d(i-1)/b(i-1)
      b(i)=b(i)-t*d(i-1)
      c(i)=c(i)-t*c(i-1)
20    continue
c
c  back substitution
c
      c(n)=c(n)/b(n)
      do 30 ib=1,nm1
      i=n-ib
      c(i)=(c(i)-d(i)*c(i+1))/b(i)
30    continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n)=(y(n)-y(nm1))/d(nm1)+d(nm1)*(c(nm1)+2.d0*c(n))
      do 40 i=1,nm1
      b(i)=(y(i+1)-y(i))/d(i)-d(i)*(c(i+1)+2.d0*c(i))
      d(i)=(c(i+1)-c(i))/d(i)
      c(i)=3.d0*c(i)
40    continue
      c(n)=3.d0*c(n)
      d(n)=d(n-1)
      return
c
50    b(1)=(y(2)-y(1))/(x(2)-x(1))
      c(1)=0.d0
      d(1)=0.d0
      b(2)=b(1)
      c(2)=0.d0
      d(2)=0.d0
      return
      end 

      SUBROUTINE gauleg(x1,x2,x,w,n)
      implicit real*8 (a-h,o-z)
      PARAMETER (EPS=3.d-14)
      dimension x(n),w(n)
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      pi=4.d0*datan(1.d0)
      do 12 i=1,m
      z=cos(pi*(i-0.25d0)/(n+0.5d0))
1     continue
      p1=1.d0
      p2=0.d0
      do 11 j=1,n
      p3=p2
      p2=p1
      p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11    continue
      pp=n*(z*p1-p2)/(z*z-1.d0)
      z1=z
      z=z1-p1/pp
      if(dabs(z-z1).gt.EPS)goto 1
      x(i)=xm-xl*z
      x(n+1-i)=xm+xl*z
      w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
      w(n+1-i)=w(i)
12    continue
      return
      end

c ... Interpolation routines from FRESCO 
c ... Complex interpolation function for complex
c ... (for uniform grids only!) 
c .... (credits to by Prof IJ Thompson)
      FUNCTION FFC4(Y,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 F(N),FFC4
      REAL*8 P,Y,P1,P2,Q,X
      DATA X/.16666666666667E0/
      P=Y
      I=P
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFC4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFC4=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFC4=F(N)
      RETURN
      END
c ... Interpolation of real function
      FUNCTION FFR4(Y,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 F(N),P,P1,P2,Q,X,FFR4
      REAL*8 Y
      DATA X/.16666666666667E0/
      P=Y
      I=P
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFR4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFR4=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFR4=F(N)
      RETURN
      END
      FUNCTION FIC(PP,F,N)
      COMPLEX*16 FIC,F(N)
      REAL*8 PP
      PARAMETER(X=.16666666666667)
      I=PP
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=PP-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FIC=(-P2*F(I)+Q*F(I+3))*(P*P1*X)+(P1*F(I+1)-P*F(I+2))*(Q*P2*.5)
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FIC=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FIC=F(N)
      RETURN
      END
c ... Another interpolation routine for real functions??? 
      FUNCTION FFR(PP,F,N)
      COMPLEX*16 F(N)
      REAL*8 PP
      REAL*8 FFR
      PARAMETER(X=.16666666666667)
      I=PP
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=PP-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFR = ( - P2*DBLE(F(I)) + Q*DBLE(F(I + 3)))*(P*P1*X)
     &        + (P1*DBLE(F(I+1)) - P*DBLE(F(I+2)))*(Q*P2*.5)
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFR=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFR=F(N)
      RETURN
      END
c ... Interpolation for complex functions
      COMPLEX*16 FUNCTION FFCI(PP,F,N,L)
      INTEGER, INTENT(IN):: N,L
      COMPLEX*16, INTENT(IN):: F(N)
      REAL*8, INTENT(IN):: PP
      PARAMETER(X=.16666666666667)
      I=PP
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=PP-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFCI=(-P2*F(I)+Q*F(I+3))*(P*P1*X)+(P1*F(I+1)-P*F(I+2))*(Q*P2*.5)
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFCI=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFCI=F(N) * ((N-1.)/PP)**(L+1)
      RETURN
      END



