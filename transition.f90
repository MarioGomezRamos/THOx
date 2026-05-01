c  DCE code 
c  Calculates coupling potentials for 3-body CDCC calculations with core excitation
c  according to the formalism first proposed in PRC74, 014606 (2006) 
c
c  R. de Diego, A.M. Moro (2013-2014)
c  ----------------------------------------------------------------------------------




c *** Use projectile WFs from THOx *****************************
c     and interpolate them at quadrature points
c *** ---------------------------------------------------------
      subroutine wf2quad
      use sistema
      use channels, only: jpiset,jpsets,nchmax
      use wfs, only: wfr,wfc,energ,idx
      use xcdcc,only:nquad,rquad,frad,exch,parch,jpch,nex,elab
      use wfs, only:  wfeig,ebin,nr,rvec,dr,rmin
      use parameters, only: maxeset
      use trace, only: cdccwf
      implicit none
      integer:: jset,j,nchan,ne,nst,iq,ie,ich,iparity
c changed in v2.3
!      real*8,pointer:: yvp(:)
!      real*8 :: faux
      complex*16,pointer:: yvp(:)
      real*8      :: raux,xjtot
      complex*16  :: faux,caux,ffc4
      character*5 :: jpi
      integer,parameter:: alpha=0d0

!      write(*,*)'jpsets=',jpsets,' nquad=',nquad,' nchmax=',nchmax
!      allocate(frad(jpsets,nener,nquad,nchmax),chann(nchmax,4,jpsets))
      allocate(frad(jpsets,maxeset,nquad,nchmax))
      write(*,'(//,2x,"***  PROJECTILE WAVEFUNCTIONS *** ")')
      
      nst=0
      do jset=1,jpsets
      ne=jpiset(jset)%nex
      nst=nst+ne
      nchan=jpiset(jset)%nchan
      xjtot=jpiset(jset)%jtot
      iparity=jpiset(jset)%partot
      write(*,*)
      write(*,'(4x,"Set",i3,"  Jpi=",1x,a5,
     &    "=>",i3," state(s)" )')
     &    jset,jpi(xjtot,iparity),ne


      do ie=1,ne
      write(*,'(8x,"#",i3,2x,"Ex=",1f8.4)')
     & ie,energ(jset,ie)
      do ich=1,nchan
c changed in v2.3
!     yvp=>wfr(jset,ie,ich,:)
      yvp=>wfc(jset,ie,ich,:)

!      if(jset.eq.2) then
!      do iq=1,nr
!      write(95,'(1f8.3,2x,100g14.6)')rvec(iq),
!     &   rvec(iq)*wfc(jset,ie,ich,iq)
!      enddo
!      write(95,*)'&'
!      endif


      do iq=1,nquad
      raux=rquad(iq)
!      faux=cfival(raux,rvec,yvp,nr,alpha) ! R(r)
      faux=FFC4((raux-rmin)/dr,yvp,nr)
      frad(jset,ie,iq,ich)=faux*raux     ! u(r)=r*R(r)
!      write(95,'(1f8.3,2x,100g14.6)')raux,frad(jset,ie,iq,ich)
      enddo ! iq=quadrature point
      enddo ! ich=channel
      enddo ! ie = energy
!      write(*,'(5x,"Set:"i2," => ",i3," state(s) considered")')
!     &  jset,ne
      enddo ! jset=j/pi set
      write(*,'(/, 5x,"=>", i3,1x,"state(s) in ",i2," set(s)")')
     & nst,jpsets


c Write energies and internal WFS for full CDCC WF
      if (cdccwf) then
         write(85,'(5a,1f8.2,a)') '** CDCC WF for: ', 
     &   trim(namep),'+',trim(namet),'  at E=',elab,' MeV'
         write(85,'(a,i3)') 
     &   '# Nb. states: ', nst
         write(85,'(a,i5,2f8.3)') '# Radial grid:', 
     &   nr,dr,rmin
      endif

c     store excitation energies, parities for coupled-channels calculations
      nex=nst
      allocate(exch(nex))
      allocate(parch(nex))
      allocate(jpch(nex))
      allocate(idx(jpsets,maxval(jpiset(:)%nex)))
      
      nst=0
      do jset=1,jpsets
      ne=jpiset(jset)%nex
      do ie=1,ne
      nst=nst+1
      exch(nst)   = energ(jset,ie)
      parch(nst)  = jpiset(jset)%partot
      jpch(nst)   = jpiset(jset)%jtot
      idx(jset,ie)= nst

      if (cdccwf) then
       write(85,'(a,i3,1f10.5,1f5.1)') '# n, Ex, Jp:   ',
     &      nst,exch(nst),jpch(nst)
      do ich=1,jpiset(jset)%nchan
      write(85,'(6g14.6)')
     &   (rvec(iq)*wfc(jset,ie,ich,iq),iq=1,nr)
      enddo !nch
      endif !cdccwf?

      enddo ! ie
      enddo ! jset        

!      deallocate(wfr)
      return
      end subroutine

c nc=1: nuclear part
c nc=2: coulomb part
      complex*16 function pot(r1,r2,iq,k,coef,nc)
      use xcdcc,only:nquad,xquad,wquad
      implicit real*8 (a-h,o-z)
      complex*16 ss,func
      integer nc
      common/vari/nq,kk
      common/varr/xcoef,xr1,xr2
!      write(*,*)'pot: nc=',nc
      kk=k
      xcoef=coef
      xr1=r1
      xr2=r2
      nq=iq
      ind=mod(nquad,2)
      if(ind.eq.1) nquad=nquad-1
      ss=(0.d0,0.d0)
      do j=1,nquad/2
      dx=xquad(j)
      ss=ss+wquad(j)*(func(dx,nc)+func(-dx,nc))
      enddo
      pot=ss
      return
      end

c This is the integrand in eq.(22) of Summers et al
c Vq is the Q-multipole of the core-target or valence-target potential
c according to :
c
c Vct(\vec rc,\vec xi)=Sum_{Q,q} \hat{Q} V^Q(rc) C_{Qq}(\hat rc) C*_{Qq}(\hat xi)
c
c and likewise for the valence-target
c
c nc=1: nuclear part
c nc=2: coulomb part
      complex*16 function func(xx,nc)
      use ptpots, only:vcore,vcorei,vval,vvali,vcorec,vvalc,rfrag,nr
!      use wfs, only:rfrag,nr
      use constants, only:e2
      use sistema, only:zc,zv,zt,rcc,zp
      implicit real*8 (a-h,o-z)
      real*8 faux,zfrag,vcoul,raux,xpos
      complex*16 xi,vq
      common/vari/iq,k
      common/varr/coef,r1,r2 ! r,R ???
c     --------------------------------------------------------------
      xi=(0.d0,1.d0)
      alpha=0.d0
      func =0 

      xrvec=dsqrt((coef*r1)**2+r2**2-2.d0*r1*r2*coef*xx)
      if (coef.lt.0) then
         zfrag=zv
      else
         zfrag=zc
      endif
 
      if ((nc.eq.2).and.(zfrag.lt.1e-4)) return

      if(nc.eq.1) then                    ! nuclear part
      xpos=(xrvec-rfrag(1))/(rfrag(2)-rfrag(1)) 
      if(coef.lt.0.d0) then
!      vvalp=fival(xrvec,rfrag,vval(:,iq),nr,alpha)
        vvalp=FFR4(xpos,vval(:,iq),nr)
!       vvalip=fival(xrvec,rfrag,vvali(:,iq),nr,alpha)
       vvalip=ffr4(xpos,vvali(:,iq),nr)
       vq=vvalp + xi*vvalip
      else             
!        vcorep=fival(xrvec,rfrag,vcore(:,iq),nr,alpha)
!        vcoreip=fival(xrvec,rfrag,vcorei(:,iq),nr,alpha)
        vcorep=FFR4(xpos,vcore(:,iq),nr)
        vcoreip=FFR4(xpos,vcorei(:,iq),nr)
        vq=vcorep+xi*vcoreip
      endif
      func=Vq*fleg(xx,k)/(2.d0*xrvec**iq) !c.f. Eq. (22) of Summers 
      endif
 
      if((nc.eq.2).and.(zfrag.gt.0)) then      ! Coulomb part
      xpos=(xrvec-rfrag(1))/(rfrag(2)-rfrag(1)) 
      if (coef.lt.0) then ! valence-target
!        Vq=fival(xrvec,rfrag,vvalc(:,iq),nr,alpha)
        Vq=FFR4(xpos,vvalc(:,iq),nr)
       else                ! core-target
!        Vq=fival(xrvec,rfrag,vcorec(:,iq),nr,alpha)
        Vq=FFR4(xpos,vcorec(:,iq),nr)
!        if (iq.eq.2) write(199,*)xx, vq
       endif
       func=Vq*fleg(xx,k)/(2.d0*xrvec**iq) !c.f. Eq. (22) of Summers
      endif     
      RETURN
      end
