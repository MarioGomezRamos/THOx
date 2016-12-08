c  DCE code 
c  Calculates coupling potentials for 3-body CDCC calculations with core excitation
c  according to the formalism first proposed in PRC74, 014606 (2006) 
c
c  R. de Diego, A.M. Moro (2013-2014)
c  ----------------------------------------------------------------------------------

      subroutine transition_targdef
      use xcdcc
      use sistema
      use wfs, only:energ,idx
      use ptpots
      use globals
      use channels
      use memory
      implicit none
      logical skip,writeff
      integer qcmin,qcmax,kcmax,coups,ilam,ilamp,iq,idtf,idti!,lambmax
      logical,external :: fail3
c     -------------------------------------------------------------------------------
      character*5:: jpi
      character*40 filename,couplings,comment,fileamps
c     -------------------------------------------------------------------------------
      real*8 factor,ptr,ttr
      real*8,allocatable:: xrad2(:),xrad3(:)
c changed to complex in v2.3
      complex*16,pointer:: ui(:),uf(:)
      complex*16 :: fprod
      complex*16 xi,pot,xsum, xsumn,xsumc !,gunc
!      parameter(Kmax=6,xKrot=0.d0,nkmax=300)
      integer nkmax
      parameter(nkmax=300)
c     ------------------------------------------------------------------------------
       logical ignorespin,noclosed
       real*8 a,b,bel,coef,coefc,coefv,e1,e2,emaxf,eminf,ep,et,finish
       real*8 fmem,fscale,drvec,jt,pi,r1,r2,radmax,rextrap,rfirst,rjp
       real*8 rmat,rmatc,rmatn,rmax,rvecin,rvecmax,start,t1,t2,vmon
       real*8 ximax,xjp1,xjp2,xjpmax,xl1,xl2,xlc,xlcp,xj1,xj2,sn1
       real*8 rmelc,rmeln,jti,jtf
       real*8,external ::vcoul,threej,sixj,rotor
       integer i1,i2,i,ibandp,ibandt,icopyt,icpot,id1,id2,ie,ie1,ie2,if1
       integer if2,Ifi,Iin,ik,ipar1,ipar2,iquad,ir,irad,irad2,irvec,k
       integer ken,kfr,kmax,kptype,kread,l,j,lc,ltr,m,m1,m2,n,nc,nchann1
       integer nchann2,nff,nIcmax,nk,npa1,npa2,npe,nq,nqmax,nqmin,nrin
       integer itex,icopyp,mproj1,mproj2,iexclude,exclude(100)
       integer itex1,itex2,nonzero
       integer initstates,isc1,isc2,jsc1,jsc2,ksc1,ksc2
       real*8 krotor,ecmi,ech
! v0.6c
!      dimension energ(nener,numfilmax),!,xrad2(:),elapsed(2),
!     .np(numfilmax)
c      real*8, allocatable:: energ(:,:)
c      integer, allocatable:: np(:)
      TYPE GFAC  ! derived variable to store the geometrical factors
       INTEGER i,f,k,l,nq
       REAL*8:: rmatc,rmatn
      END TYPE
      TYPE(GFAC) PK(nkmax)
      interface
       subroutine potdef(vc,qmax,rnc,vdef,betar)
       integer::qmax
       real*8:: vc(:),rnc,vdef(:),betar
      end subroutine
      end interface

      interface
        subroutine extpot(filename,vr,vi,nr)
        integer::nr
        character*40 filename
        real*8 :: vr(:),vi(:)
       end subroutine
      end interface

! Modified by AMM in v2.2
!      namelist /trans/ mp,mt,zt,zp,rcc,namep,namet
      namelist /trans/ skip,rcc,writeff, noclosed
      namelist /grid/ nquad,radmax,rmax,rextrap,rstep,rvecin,drvec,hin
      namelist /wf/ filename,eminf,emaxf
      namelist /coupling/ qcmin,qcmax,kcmax,lamax,coups,ncoul,ignorespin
      namelist /BEQ/ bel
!      interface
!      subroutine gauleg(x1,x2,x,w,n)
!      integer::n
!      real*8:: x1,x2,x(:),w(:)
!      end subroutine
!      end interface
      call cpu_time(start)
      debug=.false.
      xi=(0.d0,1.d0)
      kin=5
      kread=84
      kfr=4
      ken=8
      pi=acos(-1d0)
      writeff=.false.
      noclosed=.false.
c     ----------------------------------------------------------------------------      
      write(*,'(/,5x,"*********** COUPLING POTENTIALS  *********** ")')
      
c ---------------------------------------------------------------------------------
c AMoro: partial preread to set dimensions
      rextrap=0
      skip=.false. 
      read(kin,nml=trans)

      if(jpsets.eq.0) then
      write(*,'(/,5x,"** No basis states have been defined",
     &    "=> skipping this section **" )') 
      return
      endif

      if (skip) then
      write(*,'(5x,"[No coupling potentials requested ]")')
      return
      endif

c *** Quantum number limits for transition potentials
      lamax=-1
      qcmin=0
      coups=0
      qcmax=-1
      iftrans=.true.
      kcmax=-1
      ignorespin=.false.
      read(kin,nml=coupling)
      if (qcmax.lt.0) then
       write(*,*)'[=> No coupling potentials requested]'
       iftrans=.false.
       return
      endif
!      write(*,'(/,2x,"o Couplings considered for formfactors:")')
      select case(coups)
      case (0)
       write(*,*) ' [ coups=0=> All transitions will be considered]'
      case (1)
       write(*,*) ' [ coups=1=> ONLY gs->cont couplings considered]'
      end select
      
      select case(ncoul)
      case (0)
        write(*,*)' [ ncoul=0 => C+N couplings ]'
      case (1)
        write(*,*)' [ ncoul=1 => NUCLEAR couplings only ]'
      case (2)
        write(*,*)' [ ncoul=2 => COULOMB couplings only]'
      end select



c *** Radial grids ------------------------------------------------------------
      read(kin,nml=grid)
      coefc=av/mp
      coefv=-ac/mp ! = coefc-1
!      write(0,*) 'MESSING WITH SIGNS CAREFUL!!!!!!!'
!      write(*,*) 'MESSING WITH SIGNS CAREFUL!!!!!!!'
      coef=max(abs(coefc),dabs(coefv))
!      coefc=-coefc
!      coefv=-coefv
!      coef=-coef
      if (rstep.gt.1e-5) then 
      nrad2=nint((rmax-rstep)/rstep)+1
      if (rextrap.lt.rmax) rextrap=rmax
      nrad3=nint((rextrap-rstep)/rstep)+1
      allocate(xrad2(nrad2))
      allocate(xrad3(nrad3))
      else
      write(*,*)'radial step (rstep) too small!'; stop
      endif
      allocate(xquad(nquad),wquad(nquad),rquad(nquad))
      allocate(xintgn(nquad),xintgc(nquad)) ! AMoro
      a=-1.d0
      b=1.d0
      call gauleg(a,b,xquad,wquad,nquad)

      rquad=radmax*0.5d0*(1.d0+xquad)
      rvecmax=coef*radmax+rmax
      nr=nint((rvecmax-rvecin)/drvec)
      nrin=nint(radmax/hin)
      allocate(rfrag(nr))
      allocate(rvin(nrin))
      do irvec=1,nr
      rfrag(irvec)=rvecin+drvec*irvec
      enddo
      do irvec=1,nrin
      rvin(irvec)=hin*irvec
      enddo

      if(.not.allocated(vcore))  allocate(vcore(1:nr,0:nmult))
      if(.not.allocated(vcorei)) allocate(vcorei(1:nr,0:nmult))
      if(.not.allocated(vcorec)) allocate(vcorec(1:nr,0:nmult))
      if(.not.allocated(vval))   allocate(vval(1:nr,0:nmult))
      if(.not.allocated(vvali))  allocate(vvali(1:nr,0:nmult))
      if(.not.allocated(vvalc))  allocate(vvalc(1:nr,0:nmult))
      kptype=0
      vfrag => vcore
      vfragi=> vcorei
      vcoup => vcorec
      write(*,'(/,2x,"** CORE-TARGET potential:")')
      call read_fragpot(kin,kptype,zc)
      melc(1:nmult)=mel(1:nmult)

      nullify(vcoup)
      kptype=1
      vfrag =>vval
      vfragi=>vvali
      vcoup =>vvalc
      write(*,'(/,2x,"** VALENCE-TARGET potential:")')
      call read_fragpot(kin,kptype,zv)
      melv(1:nmult)=mel(1:nmult)
      do irad2=1,nrad2
      xrad2(irad2)=rmax/dble(nrad2)*irad2
      enddo

c *** Interpolate projectile wfs at quadrature points
      call wf2quad()


c *** Set limits for quantum numbers & multipoles ----------------------------------
      npe    =maxval(jpiset(:)%nex)*(ntex+1)
      xImax  =maxval(tset(:)%jtarg) !Maximum I
      xJpmax =maxval(jpiset(:)%jtot) !Maximum Jp
      nIcmax =nint(2.d0*xImax)       
!      nQmax=min(nmult,nIcmax) ! RDD

      if (ignorespin) then
      xJpmax=0
      do m=1,jpsets
      xJpmax =max(maxval(jpiset(m)%lsp(:))+0d0,xJpmax)
      enddo
      endif

!!!! AMM: Iin, Iff, need to be checked
      nQmin=max(0,qcmin)
      nQmin=qcmin
      nQmax=min(qcmax,nIcmax)
!MGR lambda prima
      if (lamax.ge.0) then
          lambpmax=min(lamax,nint(2.d0*xJpmax))
      else
          lambpmax=nint(2.d0*xJpmax)
      endif

!MGR lambda
      if (lamax.ge.0) then
          lambmax=min(lamax,lambpmax+qcmax)
      else
          lambmax=lambpmax+qcmax
      endif

!!! NEW IN v2.2g (CHECK!!)
      kmax=nqmax+lambpmax
      if (kcmax.ge.0) kmax=min(kmax,kcmax)  

!      write(*,*)'qcmin,nqmin,Iin,Ifi=',qcmin,nqmin,Iin,Ifi
!      write(*,400) nqmin,nqmax,0,min(kmax,kcmax),0,lambmax
      write(*,400) nqmin,nqmax,0,kmax,0,lambpmax,0,lambmax

400   format(/,2x,"o Quantum numbers limits:",/,5x,
     &  "- Core multipoles:",i2,"-",i2,/,5x,
     &  "- K:",i2,"-",i2,/,5x,
     &  "- Proj.-proj multipoles:", i2,"-", i2,/,5x,
     &  "- Proj.-target multipoles:", i2,"-", i2)
     
c -------------------------------------------------------------------------------------

!MGR Build lambda hipervector
       nlambhipr=0
       krotor=0
       do ilamp=0,lambpmax
         do ilam=0,lambmax
           do iq=qcmin,qcmax
             if (fail3(iq+0d0,ilam+0d0,ilamp+0d0)) cycle
             if (mod(ilamp+ilam+iq,2).ne.0) cycle
             
             nonzero=0
             loop1:do itex1=1,ntex+1
               do itex2=itex1,ntex+1
             if((abs(tset(itex1)%Krot-tset(itex2)%Krot)).gt. 1e-1) cycle
              krotor=tset(itex1)%Krot 
              if (rotor(tset(itex2)%jtarg,krotor,iq,tset(itex1)%jtarg,1)
     &              .ne.0d0) then
                   nonzero=1
                   exit loop1
                 endif
               enddo
             enddo loop1
             
             if (nonzero.eq.0) cycle
             
             nlambhipr=nlambhipr+1
             lambdahpr(nlambhipr)%lambda=ilam
             lambdahpr(nlambhipr)%lambdap=ilamp
             lambdahpr(nlambhipr)%q=iq
             write(*,*) 'Hipervector',nlambhipr
             write(*,*) 'lam',ilam,'lamp',ilamp,'Q',iq
           enddo
         enddo
       enddo
!--------------------------------------------------------------------------
       

c *** START CALCULATION OF COUPLING POTENTIALS ---------------------------------------
      write(*,'(/,2x,"** FORMFACTORS **")') 
      fmem=nex*nex*(nlambhipr+1)*nrad3*lc16/1e6
      write(*,'(5x," [ formfactors need",1f7.1," Mbytes ]")') fmem

      allocate(potQKn_targ(nquad,nrad2,0:nQmax,0:Kmax,2))
      allocate(potQKc_targ(nquad,nrad2,0:nQmax,0:Kmax,2))
      allocate(Fc(nex,nex,nlambhipr,nrad3),Fn(nex,nex,nlambhipr,nrad3))
c commented by AMoro, to save memory
!      allocate(Fc(nex,nex,0:lambmax,nrad3))
!      allocate(Fn(nex,nex,0:lambmax,nrad3)) 
!      Fc=0.d0;  Fn=0.d0
      Fc=0d0; Fn=0d0
      potQKn_targ=0d0; potQKc_targ=0d0
      rmat=0d0

       
      write(*,'(/,2x,"o V^{K,Q}(r,R) ",$)')
      call cpu_time(t1)
      call factorial(2*nQmax)

!      write(*,*)'nqmin,nqmax=',nqmin,nqmax,' kmax=',kmax,kcmax
!      write(*,*)'nquad,nrad2=',nquad,nrad2
      do iquad=1,nquad
      if (debug) write(*,*) 'nquad,iquad=',nquad,iquad
      do irad2=1,nrad2
!      write(*,*) 'irad2=',nrad2,irad2
!      do nq=0,min(nQmax,qcmax) !AMM
      do nq=nqmin,nqmax
      do k=0,kmax
! AMoro: 
! ncoul=0: coulomb + nuclear
! ncoul=1: nuclear only
! ncoul=2: Coulomb only  
      potnv=0d0; potcv=0d0; potncor=0d0;  potccor=0d0
      r1=rquad(iquad); ! internal coordinate (r)
      r2=xrad2(irad2)  ! p-t coordinate (R)
      
      if ((ncoul.eq.0).or.(ncoul.eq.1)) then  !nuclear part
        nc=1 ! nuclear part
        potncor= pot(r1,r2,nq,k,coefc,nc)
        potnv=   pot(r1,r2,nq,k,coefv,nc)
      endif
      if ((ncoul.eq.0).or.(ncoul.eq.2)) then ! coulomb part
       nc=2 ! coulomb part
       potccor= pot(r1,r2,nq,k,coefc,nc)
       potcv  = pot(r1,r2,nq,k,coefv,nc)
!      write(*,*)'coulomb:',pot(r1,r2,nq,k,coefc,nc),
!     & pot(r1,r2,nq,k,coefv,nc)
      endif       
!core
      potQKc_targ(iquad,irad2,nq,k,1)=potccor  
      potQKn_targ(iquad,irad2,nq,k,1)=potncor
!valence
      potQKc_targ(iquad,irad2,nq,k,2)=potcv  
      potQKn_targ(iquad,irad2,nq,k,2)=potnv
      enddo !k 
      enddo !nq
      
      enddo ! nrad2 (R)
      enddo ! nquad (r)
      call cpu_time(t2)
      write(*,'(4x,"[",1f6.2," secs ]")') t2-t1
     
      rfirst=rstep 
      write(*,300) rfirst,rmax,rstep
300   format(/,5x,"[Radial grid: Rmin=",1f6.3," fm,", 
     &      " Rmax=",1f6.1," fm," 
     &      " Step=",1f6.3," fm]",/)

      if (rextrap.gt.rmax) then
       write(*,301) rextrap
301   format(5x,"[ Coulomb couplings extrapolated to",1f6.1," fm]",/)
      endif
      nff=0
      exclude=1
      do m=1,jpsets
        nchann1=jpiset(m)%nchan ! nch(m)
        xjp1=jpiset(m)%jtot! xjp(m)
!        if (m.eq.2) write(0,*)'Going through 2'
      do n=m,jpsets
        nchann2=jpiset(n)%nchan ! nch(n)
        xjp2=jpiset(n)%jtot     ! xjp(n)
!        write(0,*) 'Nchan1',nchann1,'Nchan2',nchann2
!      xnlc=nint(xjp1+xjp2-abs(xjp1-xjp2))
!      do ilc=0,nlc
!        write(*,*)'nchans:',nchann1, nchann2
!     nint(dabs(xjp1-xjp2)),nint(xjp1+xjp2)
!      if(m.eq.2)write(0,*) 'Going through 2'
!      write(0,*) 'nlambhipr',nlambhipr
      do lc=1,nlambhipr
        xlc=dble(lambdahpr(lc)%lambda)
        nq=lambdahpr(lc)%q
        xlcp=dble(lambdahpr(lc)%lambdap)
        
!        if(m.eq.2)write(0,*) 'Going through 2'
        nk=0  !Nb. of P(IK) for each n,m,Lambda
        call cpu_time(t1)
        do i=1,nchann1
!      xl1=chann(i,2,n)
       xl1=jpiset(m)%lsp(i) 
       xj1=jpiset(m)%jsp(i) 
!      Iin=nint(chann(i,1,m))
! AMM: check this!!!, jc could be half-integer
       Iin=jpiset(m)%jc(i) 

       if (ignorespin) then
       if (nchann1.gt.1) 
     &   write(0,*) 'WARNING!!!, ignorespin not properly implemented for
     & multichannel states'
       xj1=jpiset(m)%lsp(i)
       Iin=0d0
       xjp1=jpiset(m)%lsp(i)
       endif


       do j=1,nchann2      
!      xl2=chann(j,2,n)
      xl2=jpiset(n)%lsp(j)
      xj2=jpiset(n)%jsp(j)
      sn1=sn
!      Ifi=nint(chann(j,1,n))
! AMM: check this!!!, jc could be half-integer
      Ifi=jpiset(n)%jc(j) 
!      WRITE(*,*)'xl1,Iin=',xl1,Iin
!      WRITE(*,*)'xl1,Iin=',xl2,Ifn
 
      if (ignorespin) then
        if (nchann2.gt.1) 
     &   write(0,*) 'WARNING!!!, ignorespin not properly implemented for
     & multichannel states'
       xj2=jpiset(n)%lsp(j)
       Ifi=0d0
       sn1=0d0
       xjp2=jpiset(n)%lsp(j)
       endif
 
      if (abs(Ifi-Iin).gt.1e-3) then
       write(*,*) 'Different cores Cycle'
       cycle
      endif 
      do k=0,kmax
!      do nq=abs(Iin-Ifi),min(nQmax,Iin+Ifi) ! AMM
!      do nq=nqmin,nqmax
      do l=0,nq !small \lambda
!      if ((qcmax.eq.0).and. .false.) then
!      write(*,*) 'Using Q=0 expression'
!      rmatn=0d0
!      rmatn=sixj(xj2,xj1,dble(K),xl1,xl2,sn)*sqrt(2d0*K+1d0)*
!     & (-1d0)**(sn1+xj2+xjp1+xj1+dble(k)+Iin)*sqrt(2d0*xj1+1d0)*
!     & sqrt(2d0*xj2+1d0)*sqrt(2d0*xl1+1d0)*sqrt(2d0*xl2+1d0)*
!     & sqrt(2d0*xjp1+1d0)*sqrt(2d0*xjp2+1d0)*
!     & sixj(xj2,xj1,dble(K),xjp1,xjp2,Iin)
!      else
!      write(*,*) 'Jpiset',m,'chan',i,'l',
!     & jpiset(m)%lsp(i)
!      write(*,*) 'To Jpiset',n,'chan',j,'l',
!     & jpiset(n)%lsp(j)
      call rmatel_tdef(m,n,lc,i,j,k,nq,l,xjp1,xjp2,rmatn,ignorespin) !I have excluded <It|Q|It'> from the matrix element, be careful
!      endif
      if (rmatn.ne.rmatn) then
        write(*,*) 'Rmat NaN for Jpi',xjp1,'Jpf',xjp2
        write(*,*) 'K',k,'Q',nq,'lambda',xlc,'lambdap',xlcp
      endif
      if (abs(rmatn).lt.1e-6) cycle
      nk=nk+1
      
      if (nk.gt.nkmax) then
        write(*,*)'NK > NKMAX=',NKMAX,' SO INCREASE NKMAX'
        stop
      endif
      pk(nk)%i     =i
      pk(nk)%f     =j
      pk(nk)%k     =k
      pk(nk)%nq    =nq  
      pk(nk)%l     =l    ! \lambda
      pk(nk)%rmatn =rmatn
      if ((m.eq.1).and.(n.eq.1).and.(verb.gt.0)) then
        if (nk.eq.1) then
        write(198,'("j=",1f4.2, " Jf=",1f4.2)')xjp1,xjp2
        write(198,'(7A7,2A10)')
     & "NK","Ch1","Ch2","LAM","K","QC","lambda","P(N)", "P(C)"
        endif
      endif
      if (verb.gt.0) write(198,1100) nk,i,j,lc,k,nq,l,
     & rmatn,xl1,xl2
1100  format(7i7,2x,2f10.3,2x,2f4.1,2i5)
      enddo ! l (lambda)
      enddo ! K
      enddo ! j (nchan)
      enddo ! i (nchan)


      if (nk.gt.0) then
      write(*,'(4x,a5,"-> ",a5," with LAM=",1i2,
     &   2x,"NK=",i3,2x, "non-zero geom. P(K)  ",$)') 
     &   jpi(xjp1,jpiset(m)%partot),
     &   jpi(xjp2,jpiset(n)%partot),nint(xlc),nk
      else        
        if (debug) write(*,'(4x,"m=",i2,"-> n=",i2," with LAM=",1i2,
     &  2x,"(NO allowed transitions)")') m,n,lc
      endif

      if (nk.gt.0) exclude(lc)=0
      



c ---------------------------------------------------------
c    Calculate formfactors F(R) for each |n,ie1> ->|m,ie2> 
c ---------------------------------------------------------
      if (nk.eq.0) cycle     
      if (debug) write(*,*)'radial integrals'
      do ie1=1,jpiset(m)%nex ! np(m)
      id1=idx(m,ie1)
      do ie2=1,jpiset(n)%nex ! np(n)
      id2=idx(n,ie2)
      nff=nff+1 !number of FF
      do irad2=1,nrad2
	do ik=1,NK
        i=pk(ik)%i 
  	j=pk(ik)%f
        k=pk(ik)%k
	nq=pk(ik)%nq
	l=pk(ik)%l
	rmatn=pk(ik)%rmatn
!  radial integral: R(ik)=f*V(K,Q)*f'
        xsumn=(0.d0,0.d0)
        xsumc=(0.d0,0.d0)
        if (abs(rmatn).lt.1e-5) goto 1200
      
!        ui=>frad(m,ie1,:,i)
!        uf=>frad(n,ie2,:,j)
!      call radquad(ui,uf,potQKn(:,irad2,nq,k))
      
      do iquad=1,nquad
      vmon=0d0
      fprod=frad(m,ie1,iquad,i)*frad(n,ie2,iquad,j)
!      write(*,*) 'CAREFUL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!      write(*,*) 'Cutting lambda'
!      if(l.ne.nq) cycle !CAREFUL
      if(l.eq.nq) then
      xintgn(iquad)=fprod*(potQKn_targ(iquad,irad2,nq,k,1)+
     &              potQKn_targ(iquad,irad2,nq,k,2))*xrad2(irad2)**l
c subtract projectile-target monopole Coulomb
      if ((nq.eq.0).and.zp*zt.gt.1e-3.and.k.eq.0.and.ncoul.ne.1)
     & then
!     &  .and.(m.eq.n).and.(ie1.eq.ie2)) then 
       vmon=VCOUL(xrad2(irad2),zp,zt,Rcc)     
      endif
       xintgc(iquad)=fprod*
     &((potQKc_targ(iquad,irad2,nq,k,1)+potQKc_targ(iquad,irad2,nq,k,2))
     & *xrad2(irad2)**l-vmon)
      
      else 
      xintgn(iquad)=fprod*
     .potQKn_targ(iquad,irad2,nq,k,1)*xrad2(irad2)**l*
     .(coefc*rquad(iquad))**(nq-l)+fprod*
     .potQKn_targ(iquad,irad2,nq,k,2)*xrad2(irad2)**l*
     .(coefv*rquad(iquad))**(nq-l)

      xintgc(iquad)=fprod*
     .potQKc_targ(iquad,irad2,nq,k,1)*xrad2(irad2)**l*
     .(coefc*rquad(iquad))**(nq-l)+fprod*
     .potQKc_targ(iquad,irad2,nq,k,2)*xrad2(irad2)**l*
     .(coefv*rquad(iquad))**(nq-l)
      endif
      xsumn=xsumn + xintgn(iquad)*wquad(iquad)
      xsumc=xsumc + xintgc(iquad)*wquad(iquad)
      enddo ! iquad

1200  fauxn=rmatn*xsumn*dsqrt(2.d0*dble(k)+1.d0)*0.5d0*radmax
      fauxc=rmatn*xsumc*dsqrt(2.d0*dble(k)+1.d0)*0.5d0*radmax
!       if (fauxn.ne.fauxn) write(*,*) 'Fauxn NaN'
!      Fc(id1,id2,lc,irad2)=Fc(id1,id2,lc,irad2)+fauxc    ! coulomb
!      Fn(id1,id2,lc,irad2)=Fn(id1,id2,lc,irad2)+fauxn    ! nuclear
      Fn(id1,id2,lc,irad2)=Fn(id1,id2,lc,irad2)+fauxn ! total
      Fc(id1,id2,lc,irad2)=Fc(id1,id2,lc,irad2)+fauxc
      enddo ! ik
      enddo ! irad2 (R)
      enddo ! ie2
      enddo ! ie1
      call cpu_time(t2)
      write(*,'("(",1f6.2," secs )")')t2-t1
      call flush(6)
      enddo ! lc (\Lamda)
      enddo ! n (jpsets)
      enddo ! m (jpsets)

      call cpu_time(finish)
!      print*, 'Time = ',finish-start,'seconds'
      write(*,1220) nff,finish-start
1220  format(/,5x,"-> [",i7, " formfactors generated in ",
     & 1x,1f6.1," secs ]",/)


!      if (.false.) then !Skipping input for fresco

c    -----------------------------------------------------------
c    Write formfactors
c    ----------------------------------------------------------- 
      if (writeff) open(kfr,file="ff.fr",status='unknown')
      ttr=0.d0
      fscale=1.d0
      npa1=0
      
      write(*,*)' Output files:'
      select case(ncoul)
      case(0) ! nuc + coul
      write(*,*) ' ff.fr:  COUL+NUC couplings for Fresco'
!      ffr=>Ff(:,:,:,:,:,:)
      case(1) ! nuc
      write(*,*) ' ff.fr: NUCLEAR couplings for Fresco'
!      ffr=>Fn(:,:,:,:,:,:)
      case(2) ! coul
      write(*,*) ' ff.fr: COULOMB couplings for Fresco'
!      ffr=>Fc(:,:,:,:,:,:)
      end select 
      
      ecmi=elab*mt/(mp+mt)
      isc1=0;isc2=0;jsc1=0;jsc2=0;ksc1=0;ksc2=0
      do idti=1,1+ntex
      jti=tset(idti)%jtarg
      do idtf=1,1+ntex
      jtf=tset(idtf)%jtarg
      npa1=0
      do if1=1,jpsets
      npa1=npa1+jpiset(if1)%nex ! np(if1)
      npa2=0
      xjp1 =jpiset(if1)%jtot
      ipar1=jpiset(if1)%partot
      

      
      if (ignorespin) then
      xjp1 =jpiset(if1)%lsp(1)
      ipar1=(-1)**xjp1
      endif
      do if2=1,jpsets
      npa2  =npa2+jpiset(if2)%nex ! np(if2)
      xjp2  =jpiset(if2)%jtot
      ipar2 =jpiset(if2)%partot
      
      if (ignorespin) then
      xjp2 =jpiset(if2)%lsp(1)
      ipar2=(-1)**xjp2
      endif
      
      do ie1=1,jpiset(if1)%nex  ! np(if1)
      e1=energ(if1,ie1)
            
     
      do ie2=1,jpiset(if2)%nex  ! np(if2)

      e2=energ(if2,ie2)

!      m1=if1+ie1-1
!      m2=if2+ie2-1
      mproj1=ie1+npa1-jpiset(if1)%nex
      mproj2=ie2+npa2-jpiset(if2)%nex
      if (mproj1.eq.1) then
      m1=idti
      else
      m1=ie1+npa1-jpiset(if1)%nex+(idti-1)*sum(jpiset(1:jpsets)%nex)
     & +ntex+1-idti !np(if1)
      endif
      if (mproj2.eq.1) then
      m2=idtf
      else
      m2=ie2+npa2-jpiset(if2)%nex+(idtf-1)*sum(jpiset(1:jpsets)%nex)
     & +ntex+1-idtf !np(if2)
      endif
      
      
      select case(coups)
      case (0) ! all couplings (do nothing)
      case (1) ! gs->cont + diagonal
        if ((if1.ne.1).or.(ie1.ne.1)) then
          if ((if2.ne.if1).or.(ie2.ne.ie1)) cycle
            if ((idti.ne.idtf)) cycle
        endif
      end select
!!!!!!!!!!!!!!!!
      if (m1.gt.m2) cycle !!!! AMORO
!      write(*,*) 'm1,m2:',m1,m2
!      write(*,*) 'mproj1,mproj2:',mproj1,mproj2
      call flush(6)
c for fresco
      do lc=1,nlambhipr
!        if (exclude(lc).eq.1) then
!          do iexclude=lc,nlambhipr-1
!             lambdahpr(iexclude)%lambda=lambdahpr(iexclude+1)%lambda
!             lambdahpr(iexclude)%lambdap=lambdahpr(iexclude+1)%lambdap
!             lambdahpr(iexclude)%q=lambdahpr(iexclude+1)%q
!             exclude(iexclude)=exclude(iexclude+1)
!          enddo     
!           nlambhipr=nlambhipr-1
!        endif
        xlc=dble(lambdahpr(lc)%lambda)
        nq=lambdahpr(lc)%q
        xlcp=dble(lambdahpr(lc)%lambdap)
      
      if((abs(tset(idti)%Krot-tset(idtf)%Krot)).gt. 1e-1) cycle
      krotor=tset(idti)%Krot
      rmelc=rotor(jtf,krotor,nq,jti,2)*(-1d0)**(jti-jtf)
      rmeln=rotor(jtf,krotor,nq,jti,1)*(-1d0)**(jti-jtf)
      
      if (noclosed) then
      ech=ecmi+energ(1,1)-e1-tset(idti)%extarg
      if (ech.lt.0) then
      rmelc=0d0
      rmeln=0d0
      endif
      ech=ecmi+energ(1,1)-e2-tset(idtf)%extarg
      if (ech.lt.0) then
      rmelc=0d0
      rmeln=0d0
      endif   
      endif
      
      if ((abs(rmeln)+abs(rmelc)).lt.1e-10) cycle  
      
!c \hat{J_p}*hat{J'_p}*(2*\Lambda+1)*(-1)^\Lambda
!      factor=(2d0*ptr+1)*sqrt(2*xjp1+1)*sqrt(2*xjp2+1)*(-1)**ptr
!      factor=factor*(-1)**ptr

      if (debug) then
      write(*,*)' Conversion factor in F(r)[fres]=factor*F(r)',factor
      endif
!      write(330,*) 'mproj1',mproj1,'mproj2',mproj2,'jti',jti,'jtf',jtf,
!     & 'm1',m1,'m2',m2
!      write(330,*) 'lam',xlc,'Q',nq,'lamp',xlcp
!      write(330,*)

!      if ((m1.eq.3).and.((m2.eq.4).or.(m2.eq.5))) 
!     & write(330,*)'Rmel',rmeln,rmelc
      if (mproj1.le.mproj2) then
      if ((sum(abs(Fc(mproj1,mproj2,lc,:)))+
     & sum(abs(Fn(mproj1,mproj2,lc,:)))).lt.1e-10) cycle
      else
      if ((sum(abs(Fc(mproj2,mproj1,lc,:)))+
     & sum(abs(Fn(mproj2,mproj1,lc,:)))).lt.1e-10) cycle
      endif
!      write(330,*) 'Non-zero'
!      write(330,*)
      write(comment,'(2x,"<",i3,"|",i2,"|",i3,">")') m1,nint(xlc),m2
      if (writeff) 
     & write(kfr,500) nrad3,rstep,rfirst,fscale,nint(xlc),xlcp,nq+0d0
     & ,m2,m1,comment
      if (verb.ge.4) write(120,'("# <",i3,"|",i2,"|",i3,">")') m1,lc,m2
      
      do irad=1,nrad2
      r2=xrad2(irad)
      if (mproj1.le.mproj2) then
      fauxc=Fc(mproj1,mproj2,lc,irad)*rmelc
      fauxn=Fn(mproj1,mproj2,lc,irad)*rmeln  
      else
      fauxc=Fc(mproj2,mproj1,lc,irad)*rmelc
      fauxn=Fn(mproj2,mproj1,lc,irad)*rmeln 
      endif
 
!      fauxn=0d0
!      if ((ncoul.eq.0).or.(ncoul.eq.1)) then
!       fauxn=Fn(m1,m2,lc,irad)
!      endif      
!      if ((ncoul.eq.0).or.(ncoul.eq.2)) then
!       fauxc=Fc(m1,m2,lc,irad)   
!      endif          
!      write(kfr,'(2x,1g16.10,2x,1g16.10)') factor*(fauxc+fauxn)
!      write(120,'(1x,1f8.3,2x,2g16.8)') r2,factor*(fauxc+fauxn)
       if (writeff) 
     &      write(kfr,'(2x,1g16.10,2x,1g16.10)') fauxc+fauxn
      if (verb.ge.4) 
     & write(120,'(1x,1f8.3,2x,2g16.8)') r2,fauxc+fauxn
      enddo !irad
!      write(440,*) 'm1',m1,'m2',m2,'xlc',nint(xlc),'xlcp',xlcp,'Q',nq,
!     & 'rmeln',rmeln,'rmelc',rmelc
c Extrapolate Coulomb formactors from R=Rmax to Rextrap
      fauxc=fauxc/rmelc
      if (rextrap.gt.rmax) then
      caux=real(fauxc)*rmax**(xlc+1) 
      do ir=nrad2+1,nrad3
!      r2=rmax+dble(rstep)*ir
      r2=rstep+rstep*(ir-1)
      if ((nint(xlc).gt.0).and.(ncoul.ne.1)) then
      fauxc=caux/r2**(xlc+1d0)
      endif
     
      if (writeff) write(kfr,'(2x,1g16.10,2x,1g16.10)')fauxc*rmelc
      if (verb.ge.4) write(120,'(1x,1f8.3,2x,2g16.8)')r2,
     & fauxc*rmelc
      Fc(mproj1,mproj2,lc,ir)=fauxc
      Fn(mproj1,mproj2,lc,ir)=0
      if ((ir.eq.nrad3).and.(abs(fauxc).gt.0.01d0)) 
     & write(330,500) nrad3,rstep,rfirst,fscale,nint(xlc),xlcp,nq+0d0
     & ,m2,m1,comment
      enddo !nrad3
      endif ! rextrap> rmax

      enddo !lc (next multipole)
     
      if (verb.ge.4) then
c ---------- Nuclear
!      write(10,800) '#','Jp1=',xjp1,'parity1=',ipar1,'Jp2=',xjp2,
!     & 'parity2=',ipar2,'ener1=',e1,'ener2=',e2
c ---------- Coulomb
!      write(11,800) '#','Jp1=',xjp1,'parity1=',ipar1,'Jp2=',xjp2,
!     & 'parity2=',ipar2,'ener1=',e1,'ener2=',e2
c ---------- nuclear + coulomb
      write(12,800) '#','Jp1=',xjp1,'parity1=',ipar1,'Jp2=',xjp2,
     & 'parity2=',ipar2,'ener1=',e1,'ener2=',e2
      endif


      i1=nint(dabs(xjp1-xjp2))
!      i2=nint(xjp1+xjp2)
      i2=min(nint(xjp1+xjp2),lambmax)
      if (verb.ge.4) then
      do irad=1,nrad3
!      r2=xrad2(irad)
      r2=rstep+rstep*dble(irad-1)
!      write(12,900) r2, (Fn(m1,m2,i,irad),i=i1,i2)
      write(10,900) r2, (Fn(m1,m2,i,irad),i=i1,i2)
      write(11,900) r2, (Fc(m1,m2,i,irad),i=i1,i2)
500   format(i4,3f8.4,i4,2f4.0,2i4,a35)
800   format (a,2(a,(f8.4),2x,a,i3,2x),2(a,(f8.4),2x),/)
820   format (a,a5,a,a5,2(a,(f8.4),2x))
900   format (1f8.3,2x,100(e12.4,2x,e12.4,4x))
      enddo !irad
!      write(10,*) '&'
!      write(11,*) '&'
      write(12,*) '&'
      if (verb.ge.4)  write(120,*)'&'
      endif
      enddo !ie2
      enddo !ie1
      enddo !if2
      enddo !if1
      enddo !idtf
      enddo !idti
      write(*,*)
      if (verb.ge.2) then
!        write(*,*)'fort.10: nuclear formfactors'
!        write(*,*)'fort.11: Coulomb formfactors'
        write(*,*)'fort.12: nuclear + Coulomb formfactors'      
      endif ! verb
      call flush(330)
c     -----------------------------------------------------------
c     Write energies in Fresco fromat
c     ----------------------------------------------------------- 
      open(ken,file="states.fr",status='unknown')     
      i=0
      icopyt=0
      icpot=1
      icopyp=0
!      ibandt=0
!      et=0.d0
!      jt=0d0
      write(ken,*) m1
!Let's move inelastic to the front   
      do itex=1,ntex+1
      jt=tset(itex)%jtarg
      ibandt=tset(itex)%partarg
      et=tset(itex)%extarg
      i=0
      icopyt=0
      rjp   =jpiset(1)%jtot    ! xjp(n)
      ibandp=jpiset(1)%partot  ! ipar(n)
      if (ignorespin) then
      rjp =jpiset(1)%lsp(1)
      ibandp=(-1)**rjp
      endif
      i=i+1
      if (itex.ne.1) icopyp=1
      ep=0
      write(ken,11511) rjp,ibandp,ep,icpot
      if(icopyp.ne.0) write(ken,1156) icopyp
      if(icopyt==0) write(ken,1153) jt,ibandt,et
      enddo !itex
!-----------------------------------------------------------------   


!And now breakup
      do itex=1,ntex+1
      jt=tset(itex)%jtarg
      ibandt=tset(itex)%partarg
      et=tset(itex)%extarg
      i=ntex+1
      icopyt=itex
      do n=1,jpsets
      rjp   =jpiset(n)%jtot    ! xjp(n)
      ibandp=jpiset(n)%partot  ! ipar(n)
      if (ignorespin) then
      rjp =jpiset(n)%lsp(1)
      ibandp=(-1)**rjp
      endif
      initstates=1
      if (n.eq.1) initstates=2
      do ie=initstates,jpiset(n)%nex    ! np(n)
      i=i+1
      if (itex.ne.1) icopyp=i
      if (itex.eq.1) icopyp=0
      ep=energ(n,ie)-energ(1,1)
      write(ken,11511) rjp,ibandp,ep,icpot
      if(icopyp.ne.0) write(ken,1156) icopyp
      if(icopyt==0) write(ken,1153) jt,ibandt,et
      if(icopyt.ne.0) write(ken,1154) icopyt
      enddo !ie
      enddo !n
      enddo !itex

11511 format(' &States jp=',f4.1,' ptyp=',i2,' ep=',f8.4,'  cpot=',i3)
 1153 format('         jt=',f4.1,' ptyt=',i2,' et=',f8.4,' /')
 1154 format('         copyt=',i2,' /')
 1155 format(' Skipping state #',i3,' J/pi',f4.1,i2,' and Ex=',f8.4)
 1156 format('         copyp=',i2)
c     ----------------------------------------------------------------
!      endif  !End skipping

c      total = etime(elapsed)
c      print*,'total=',total,'user=',elapsed(1),
c     .'system=',elapsed(2)
      call flush(ken)
      call flush(kfr)
      close(ken); close(kfr)
!      deallocate(fn,fc)
      end


