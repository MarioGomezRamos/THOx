!      SUBROUTINE rmatel(nfil1,mset,xlambdac,nch1,nch2,K,nQ,lambda,
      SUBROUTINE rmatel(nset,mset,xlambdac,nch,mch,K,nQ,lambda,
     .xjp1,xjp2,xKc,rmatn,rmatc) 
!      use xcdcc, only: chann!,rmatn,rmatc
      use factoriales
      use channels,only:spchan,jpiset,sn
      
      implicit real*8(a-h,o-z) 
      real*8 rmatn,rmatc,xjp1,xjp2
      integer nch,mch,nset
      pi=dacos(-1.d0)
      rmatn=0d0; rmatc=0d0
!      xI1=chann(nch,1,nset)
!      xl1=chann(nch,2,nset)
!      s1=chann(nch,3,nset)
!      xj1=chann(nch,4,nset)
!      xI2=chann(mch,1,mset)
!      xl2=chann(mch,2,mset)
!      xj2=chann(mch,4,mset)

!       s1=spchan(nset,1)%sn
       s1=sn 
c AMM: initial???
       xl1=jpiset(nset)%lsp(nch)
       xj1=jpiset(nset)%jsp(nch)
       xI1=jpiset(nset)%jc(nch)
       
c AMM: final??
       xl2=jpiset(mset)%lsp(mch)
       xj2=jpiset(mset)%jsp(mch)
       xI2=jpiset(mset)%jc(mch)

      lambdac=nint(xlambdac)

!      if (lambdac.eq.2) write(197,'(2f8.4,i3,2f8.4,4i3)') 
!     &  xi1,xl1,nch,xi2,xl2,mch,k,nq,lambda

      xKg=dsqrt(2.d0*dble(K)+1.d0)
      Qg=dsqrt(2.d0*dble(nQ)+1.d0)

      xl1g=dsqrt(2.d0*xl1+1.d0)
      xl2g=dsqrt(2.d0*xl2+1.d0)

      xj1g=dsqrt(2.d0*xj1+1.d0)
      xj2g=dsqrt(2.d0*xj2+1.d0)

      xI1g=dsqrt(2.d0*xI1+1.d0)
      xI2g=dsqrt(2.d0*xI2+1.d0)

 
      if(lambdac.gt.(K+lambda).or.lambdac.lt.abs(K-lambda).or.
     .(-1)**nint(dble(K)+dble(lambda)+xlambdac).ne.1) return
      if(xKc.eq.0.d0.and.(-1)**nint(dble(nQ)+xI1+xI2).ne.1) return
      if((nq.eq.0).and.(xI1.ne.xI2)) return

      xfac=(-1.d0)**nint(xj2+xl1+xl2+s1+dble(nQ))
     .*Qg**2*xKg*xj1g*xj2g*xl1g*xl2g
     .*dsqrt(fact(2*nQ)/(fact(2*lambda)*fact(2*(nQ-lambda))))
     .*WIGN3J(dble(K),dble(lambda),xlambdac,0.d0,0.d0,0.d0)
      n1=nint(dabs(xl2-xl1))
      n2=nint(xl1+xl2)
      red=0.d0
      do i=n1,n2
      if((-1)**nint(dble(K)+dble(nQ-lambda)+dble(i)).ne.1) cycle
      if((-1)**nint(dble(i)+xl1+xl2).ne.1) cycle
      if(WIGN6J(dble(i),xlambdac,dble(nQ),dble(lambda),
     .dble(nQ-lambda),dble(K)).eq.0.d0) cycle
      if(WIGN6J(xj1,xj2,dble(i),xl2,xl1,s1).eq.0.d0) cycle
      if(WIGN9J(xjp1,xjp2,xlambdac,xj1,xj2,dble(i),xI1,xI2,dble(nQ))
     ..eq.0.d0) cycle
      red=red+(2.d0*dble(i)+1.d0)
     .*WIGN3J(dble(K),dble(nQ-lambda),dble(i),0.d0,0.d0,0.d0)
     .*WIGN3J(dble(i),xl1,xl2,0.d0,0.d0,0.d0)
     .*WIGN6J(dble(i),xlambdac,dble(nQ),dble(lambda),dble(nQ-lambda),
     .dble(K))
     .*WIGN6J(xj1,xj2,dble(i),xl2,xl1,s1)
     .*WIGN9J(xjp1,xjp2,xlambdac,xj1,xj2,dble(i),xI1,xI2,dble(nQ))
      enddo
!!! CHECK IF THIS IS VALID ALSO FOR NQ=0!!!!!!! AMORO!!!!!!!!!!!!!!!1
!      strmat=WIGN3J(dble(nQ),xI1,xI2,0.d0,-xKc,xKc)
!     .*(-1.d0)**nint(xI2+xKc+dble(nQ))*xI1g*xI2g        ! nuclear
      strmat=rotor(xI1,xKc,nQ,xI2,1)
      rmatn=xfac*red*strmat
      strmat=rotor(xI1,xKc,nQ,xI2,2)                     ! coulomb
      rmatc=xfac*red*strmat
123   return
      end

!Reduced matrix element for target excitation NOT INCLUDING <It|Q|It'>
      SUBROUTINE rmatel_tdef(nset,mset,ilam,nch,mch,K,nQ,lambda,
     .xjp1,xjp2,rmatn,ignorespin) 
!      use xcdcc, only: chann!,rmatn,rmatc
      use factoriales
      use channels,only:spchan,jpiset,sn
      use xcdcc, only: lambdahpr
      
      implicit none
      logical ignorespin
      real*8 rmatn,xjp1,xjp2,pi,xKg,Qg,xl1g,xl2g,xj1g,xj2g,sn1
      integer nch,mch,nset,mset,ilam,K,nQ,lambda
      real*8 xlamb,xlambp,xl1,xl2,xj1,xj2,xI1,xI2,xfac
      real*8,external:: threej,sixj
      real*8 xj1aux,xj2aux,xjp1aux
      
      pi=dacos(-1.d0)
      rmatn=0d0
      xfac=0d0


       xl1=jpiset(nset)%lsp(nch)
       xj1=jpiset(nset)%jsp(nch)
       xI1=jpiset(nset)%jc(nch)
       
       xl2=jpiset(mset)%lsp(mch)
       xj2=jpiset(mset)%jsp(mch)
       xI2=jpiset(mset)%jc(mch)
       
       sn1=sn
       
       if (ignorespin) then
       xj1=xl1
       xI1=0d0
       xj2=xl2
       xI2=0d0
       sn1=0d0
       endif
       
       
       if(xI2.ne.xI1) then
         write(*,*) 'Something wrong. Core i',xI1,'not equal to core f',
     &   xI2 
         return
       endif


!      if (lambdac.eq.2) write(197,'(2f8.4,i3,2f8.4,4i3)') 
!     &  xi1,xl1,nch,xi2,xl2,mch,k,nq,lambda
      xlamb=lambdahpr(ilam)%lambda+0d0
      xlambp=lambdahpr(ilam)%lambdap+0d0
      if (lambdahpr(ilam)%q.ne.nQ) then
      write(*,*) 'Something wrong: Q in hipervector',lambdahpr(ilam)%q,
     & 'not equal to Q',nq
      endif
      
      if((nint(xlamb).gt.(K+lambda)).or.(nint(xlamb).lt.abs(K-lambda))
     &.or.(mod(nint(xlamb)+K+lambda,2).ne.0)) return
      if(nint(xlambp).gt.(K+nq-lambda).or.(nint(xlambp)
     & .lt.abs(K-nq+lambda)).or.(mod(nint(xlambp)+K+nq-lambda,2).ne.0)) 
     & return
      
      xKg=dsqrt(2.d0*dble(K)+1.d0)
      Qg=dsqrt(2.d0*dble(nQ)+1.d0)
      xl1g=dsqrt(2.d0*xl1+1.d0)
      xl2g=dsqrt(2.d0*xl2+1.d0)
      xj1g=dsqrt(2.d0*xj1+1.d0)
      xj2g=dsqrt(2.d0*xj2+1.d0)



      xfac=(-1d0)**(xj2+xj1+xl1+xl2+sn1+xjp1+xI1)*
     & dsqrt(2.d0*xjp1+1.d0)*dsqrt(2.d0*xjp2+1.d0)*2d0*sqrt(pi)*
     & dsqrt(2.d0*xlamb+1.d0)*xKg*Qg**2*xl1g*xl2g*xj1g*xj2g*
     & dsqrt(fact(2*nQ)/(fact(2*lambda)*fact(2*(nQ-lambda))))*
     & threej(dble(K),dble(lambda),xlamb,0.d0,0.d0,0.d0)*
     & (2.d0*xlambp+1.d0)*threej(xl1,xl2,xlambp,0.d0,0.d0,0.d0)*
     & threej(dble(K),dble(nq-lambda),xlambp,0.d0,0.d0,0.d0)*
     & sixj(xjp2,xjp1,xlambp,xj1,xj2,xI1)*
     & sixj(xj1,xj2,xlambp,xl2,xl1,sn1)*
     & sixj(xlambp,xlamb,dble(nQ),dble(lambda),dble(nQ-lambda),dble(K))
     

!      write(330,*)'li',xl1,'lf',xl2,'ji',xj1,'jf',xj2,'s',sn,'K',k
!      write(330,*)'lamp',xlambp,'Q',nq,'lambda',lambda,'lam',xlamb
!      write(330,*) 'xfac',xfac
 
!!! CHECK IF THIS IS VALID ALSO FOR NQ=0!!!!!!! AMORO!!!!!!!!!!!!!!!1
!      strmat=WIGN3J(dble(nQ),xI1,xI2,0.d0,-xKc,xKc)
!     .*(-1.d0)**nint(xI2+xKc+dble(nQ))*xI1g*xI2g        ! nuclear
!      strmat=rotor(xI1,xKc,nQ,xI2,1)
      rmatn=xfac
!      if (abs(rmatn).gt.1e-6) then
!      write(330,*) 'Coupling states'
!      write(330,*) '[(',nint(xl1),',',sn,')',xj1,',',xI1,']',xjp1
!      write(330,*) '[(',nint(xl2),',',sn,')',xj2,',',xI2,']',xjp2
!      write(330,*) 'K',k,'Q',nq,'lamb',lambda,'Lambda',xlamb,'Labmda_p',
!     & xlambp
!      endif
!      strmat=rotor(xI1,xKc,nQ,xI2,2)                     ! coulomb
123   return
      end

c ------------------------------------------------------------
c Reduced matrix element in ROTOR model
C nc=1: coulomb
C nc=2: nuclear
c
c ROTOR= sqrt(2I+1) < I  K lambda 0 | I' K >
c with some extra phase???
c --------------------------------------------------------------
      function rotor(xI1,xK,Q,xI2,nc)
!      use sistema, only: bel
      implicit none
!      implicit real*8(a-h,o-z)
      real*8:: rotor,pi,zero,mnq
      real*8:: clebsg,WIGN3J
      real*8:: xI1,xI2,xI1g,xI2g,xK,xQ,Qg
      integer:: nc,Q
      zero=0d0
      pi=acos(-1d0)
      rotor=0d0
      xI1g=dsqrt(2.d0*xI1+1.d0)
      xI2g=dsqrt(2.d0*xI2+1.d0)
      xQ=Q
      Qg=dsqrt(2.d0*xQ+1.d0)

      if ((q.eq.0).and.(xI1.ne.xI2)) return
      if (nc.eq.1) then      ! ------ nuclear
       rotor=WIGN3J(xq,xI1,xI2,zero,-xK,xK)
     .      *(-1.d0)**nint(xI2+xK+xq)*xI1g*xI2g 
      else if (nc.eq.2) then ! ------ Coulomb
      if (Q.eq.0) then
!  AMoro (CHECK FACTORS HERE!!!!!) 
!       rotor=1d0/sqrt(4d0*pi)
       rotor=WIGN3J(xq,xI1,xI2,zero,-xK,xK)
     .*(-1.d0)**nint(xI2+xK+xq)*xI1g*xI2g 
      else 
!       Mnq=bel(xq)
!       rotor=Mnq*dsqrt(2.d0*xI1+1.d0)*
       rotor=dsqrt(2.d0*xI1+1.d0)*
     . clebsg(xI1,xq,xI2,xK,zero,xK)*
     . (-1.d0)**(0.5d0*(xI1-xI2+dabs(xI1-xI2)))
!!!! AMoro: these factors already included in radial formfactor
!       rotor=rotor*dsqrt(4.d0*pi)/Qg
!!!!!
!cc      rotor=WIGN3J(xq,xI1,xI2,zero,-xK,xK)
!cc     .      *(-1.d0)**nint(xI2+xK+xq)*xI1g*xI2g 
!cc      rotor=Mnq*dsqrt(4.d0*pi)/Qg*rotor

      endif
      endif
      return
      end

      
