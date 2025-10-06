c *** --------------------------------------------------------------
c **  Scattering amplitudes & cross sections
c *** --------------------------------------------------------------
      subroutine xsecs(kin,ncc,iexgs)
      use xcdcc,    only: smats,elab,jpch,jtmin,jtmax,nex,exch,parch,
     &                    famps0,jbord,jump,binset,rel
      use channels, only: jptset,jpiset,jpsets
      use factorials
      use sistema
      use constants
      use parameters, only:  maxchan,maxsets
      use globals,  only: debug,written,verb,mu12
      use nmrv,     only: hort
      use wfs,      only: nr,dr,energ,rvec,wfr
      implicit none
      logical :: doublexs,triplexs,phixs,jsets(maxsets),kinset(maxsets)
      logical :: alphaxs
c     ---------------------------------------------------------
      integer, parameter :: klog=99
      integer   , parameter:: kamp=136, kxs=300,kfam=137,
     &                        ksmat=170, !ksj=157,
     &                        nearfa=1
      integer :: istat,kin,iset,n
      integer :: inc,iexgs,icc,ncc,nch,ni,nf,lmax,iex,lfmax,ith,icore
      integer :: partot,nth,njt,parold,maxleg,ner
      integer :: iam,im,imp,mlp,li,lf,nmi,nmf,mi,mf,ijt,l,nmfmax
c     ------------------------------------------------------------
      real*8  :: kron,cleb,ylmc,ylmc2 !external functions
      real*8    , parameter:: zero=0
      real*8    , parameter:: theps=1e-3,alpha=0d0
      real*8    , allocatable:: lfv(:),xv(:)
      real*8  :: jtot,jpgs,jpi,jpf,jpold,jtarg
      real*8  :: kcmi,kcmf,exc,vi,vf,tkch
      real*8  :: ecmi,ecmf,etai,etaf,mupt,ermin,ermax
      real*8  :: excore,ei,ef,ki,kf
      real*8  :: thmin,thmax,thcut,dth,cth,sth,s2,th,thrad
      real*8, allocatable :: pl(:,:),delci(:),delcf(:)
      real*8  :: rm,rmp,rmlp,lri,lrf
      real*8  :: factor, r1,r2,delif
      real*8  :: xsruth,xs,xstot,xrtot,xsj,xr,xrj,sigex(nex),xsaux
c     -----------------------------------------------------------
      complex*16, parameter:: zz=cmplx(0.,1.)
      complex*16,allocatable:: ampl(:,:,:),fam(:,:,:)
      complex*16, allocatable:: ampaux(:,:,:),gaux(:)
      complex*16 smat,c1,faux,fc,phci,phcf,caux,cfival,resc                     
c     -----------------------------------------------------------
      CHARACTER  PARITY(3)
      character*40 fileamp
      DATA PARITY / '-','?','+' /
c ---------------------------------------------------------------
c PLM
!      REAL*8 PM(0:1000,0:1000),PD(0:1000,0:1000)

      complex*16,allocatable:: ampauxJ(:,:,:,:) !MGR
      integer :: njtot,nlp,ijtot,ilp,njtotcalc,ijtotcalc,nspline,ispline
      integer :: jspl
      integer,allocatable :: interpJ(:)
      real*8 ljmax


      namelist/xsections/ thmin,thmax,dth,fileamp,doublexs,phixs,jsets,
     &                    ermin,ermax,ner,icore,thcut,
     &                    triplexs,rel,alphaxs!GR

c initialize -------------------------------------------------
      written(kamp)=.false.
      written(kxs)=.true.
      written(kfam)=.true.
      
!      written(ksj)=.false.
      written(kxs+1:kxs+min(nex,9))=.true.
      pi=acos(-1d0)
      thmin=0; thmax=0; dth=0; thcut=0
      doublexs=.false. ; triplexs=.false. ; phixs=.false.
      rel=.false.
      ermin=-1; ermax=-1
      jsets(:)=.true.        
      fileamp=""
      ner=0
      
c     -----------------------------------------------------------
      rewind(kin); 
      read(kin,nml=xsections)
! Commented nov/19     
!      if ((.not.doublexs).and.(.not.triplexs)) return

      if (dth.lt.1e-6) then
       write(*,*)' Angular step (dth) too small!'; stop
      endif
       
      if (thmax.gt.thmin) then
        nth=nint((thmax-thmin)/dth)+1
      else
        write(*,*)'thmax<thmin!; Aborting'; stop
      endif
      if (thcut.lt.1e-6) thcut=thmax 

      if (ner.eq.0) then
         write(*,*)'NER=0!!';
      endif

! AMM: July '19
      if (phixs) then
        write(*,*)'Triple diff. xsection with NO phi integration!'
      endif 
      
      if (ncc.eq.0) goto 1000

      jpgs =jpch(iexgs)
      jtarg=0.
      mupt =mp*mt/(mp+mt)*amu 
      ecmi =elab*mt/(mp+mt)
      kcmi =sqrt(2*mupt*ecmi)/hc
      etai =zp*zt*mupt/kcmi/hc/finec 
      lmax =nint(jtmax+maxval(jpch))



      write(*,'(/,5x,"*** SCATTERING AMPLITUDES AND X-SECTIONS ***",/)')
      write(*,320) elab, ecmi,kcmi, etai
320   format(3x,"Elab=",1f7.2, 5x,'Ecm=',1f8.3,5x,
     &  'Kcm=',1f6.3, 4x,'Eta param. =',1f7.4)
      write(*,'(/,3x, "o Jtmin=",1f4.1,3x,"Jtmax=",1f5.1,/ )')
     &  jtmin,jtmax


c *** Factorials (CHECK IF THIS IS NEEDED!!!!!!!!!!!!!!)
!      lmax=nint(jtmax+maxval(jpch))     
!      write(*,*)'famps0: lmax=',lmax
!      call factorialgen(2*lmax)


c *** Compute Coulomb phases for l< lmax
       allocate(delci(0:lmax),stat=istat)
       if (istat>0) then
         write(*,*) 'readsmat: allocating memory for DELC failed!';
         stop
       endif
       call coulph (etai,delci,lmax) !Coulomb phase-shifts for inc channel



c *** Compute & store amplitudes (Fresco, CPC Eq. (3.30))   
c *** (zero target spin assumed here!) 
      nmi    =2*jpgs+1
      nmfmax =2*maxval(jpiset(:)%jtot)+1
      allocate(famps0(nex,nth,nmi*nmfmax))

      if (nex.lt.1) then
         write(*,*)'solvecc: nex<1!'; stop
      endif
      if (iexgs.lt.1) then
         write(*,*)'solvecc: inc=',iexgs; stop
      endif

      do iex=1,nex  ! loop in projectile states
      jpf  = jpch(iex) 
      exc  = exch(iex)
      tkch=ecmi+exch(iexgs)-exc
      if (tkch.lt.0) cycle    !closed channel
      kcmf = sqrt(2*mupt*tkch)/hc
      etaf = zp*zt*mupt/kcmf/hc/finec 
      nmf  = 2*jpf +1
      lfmax= nint(jtmax+jpf)
      maxleg=max(1,lfmax)
      allocate(delcf(0:lfmax),stat=istat)
      call coulph(etaf,delcf,lfmax) !Coulomb phase-shifts for final channel

!MGR let's find the maximum number of L' that there are for any Jtot
      nlp=0
      njtot=0
      ljmax=0d0
      do icc=1,ncc
        jtot  =jptset(icc)%jtot
        njtot =max(njtot,floor(jtot))
        nch   =jptset(icc)%nchan
        do nf=1,nch
          if (jptset(icc)%nchan.eq.0) cycle
          if (jptset(icc)%idx(nf).ne.iex) cycle 
          lf  =jptset(icc)%l(nf);
          ljmax=max(ljmax,abs(jtot-lf-0d0))           
        enddo  
      enddo
      njtot=njtot-floor(jbord(1))
      nlp=nint(2d0*ljmax+1d0)
!      write(0,*) 'nlp',nlp,'njtot',njtot
!-----------------------------------------------------------------------


      write(kfam,'(4f6.1,i5,i2,1f10.3)') jpgs,jtarg,jpf,
     & jtarg,nth,nearfa,elab
!      write(1370,'(a,i3,a,1f8.4)') 'iex=',iex,' Ex=',exc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST AMPS
      write(kxs,302) iex,exc,jpf,parity(parch(iex)+2)
      if (iex.le.10) write(kxs+iex,302) 
     &      iex,exc,jpf,parity(parch(iex)+2)
300   format(3x,'o Angle-indep amplitudes for state #',i3,
     &   ' Ex=',1f7.3,' MeV',' J/pi=',1f5.1,a1)     
301   format(1x,'# f(theta) amplitudes for state #',i3,
     &   ' Ex=',1f7.3,' MeV',' J/pi=',1f5.1,a1) 
302   format(1x,'# Cross sections for state #',i3,
     &   ' Ex=',1f7.3,' MeV',' J/pi=',1f5.1,a1)

      allocate(ampl(0:lfmax,nmi,nmf)) ! angle-independent amplitudes
c1     allocate(ampaux(0:lfmax,nmi,nmf),lfv(0:lfmax) )
             allocate(ampauxJ(0:njtot,nlp,nmi,nmf)) !MGR
       allocate(interpJ(0:njtot)) !MGR
      ampl(:,:,:)=0; 
c1      ampaux(:,:,:)=0; lfv(:)=0
      ampauxJ(:,:,:,:)=0d0;interpJ(:)=1 !MGR
      ijtot=-1 !MGR
      do icc=1,ncc
       partot=jptset(icc)%partot
       jtot  =jptset(icc)%jtot
       nch   =jptset(icc)%nchan
       if (jptset(icc)%interp) then
!          write(*,*)'skipping jtot=',jtot
          cycle 
       endif
       if (jtot.ge.jbord(1)) then
       ijtot=floor(jtot-jbord(1))!MGR
       interpJ(ijtot)=0!MGR
       endif
       do ni=1,nch
       if (jptset(icc)%idx(ni).ne.iexgs) cycle ! no inc waves 
       li   =jptset(icc)%l(ni)
       lri  =li
       jpi  =jptset(icc)%jp(ni)             ! should be just jpgs
       nmi  =2*jpi+1

       do nf=1,nch
       if (jptset(icc)%nchan.eq.0) cycle
       if (jptset(icc)%idx(nf).ne.iex) cycle 
       lf  =jptset(icc)%l(nf);   
!       nlf=nlf+1;     
!       lfv(icc,nlf)=lf 
       lrf =lf
       jpf =jptset(icc)%jp(nf)
       nmf =2*jpf+1
       smat=smats(icc,ni,nf) 
       
       if (lf.gt.lfmax) stop 'internal error; lf >lfmax in famps0!' 
325    format(' Chan:',2i5,2x,2f6.1,2x,"S=",2f11.7)

       do im=1,nmi
        rm= -jpi + (im -1)
        do imp=1,nmf
        rmp=-jpf + (imp-1)
        mlp=nint(rm-rmp)
        rmlp=rm-rmp
        if (abs(mlp).gt.lf) cycle        
        r1=cleb(lri,zero,  jpi, rm, jtot,rm)
     &    *cleb(lrf,rm-rmp,jpf, rmp,jtot,rm)
        r2=sqrt(kcmf/kcmi)  ! this needs to be generalized for different partitions!!
        c1=sqrt(pi)/zz/kcmi
        phci =exp(zz*(delci(li)-delci(0)))
        phcf =exp(zz*(delcf(lf)-delcf(0)))
c Satchler (4.58) for zero-spin target:
        ampl(lf,im,imp)= ampl(lf,im,imp)  
     &                 + phci*phcf*r1*r2*c1*sqrt(2*lri+1)  
     &                 * (smat-kron(ni,nf))*YLMC2(lf,mlp) 
c1        ampaux(lf,im,imp)= ampaux(lf,im,imp)  
c1     &                   + phci*exp(-zz*(delcf(lf)-delcf(0)))
c1     &                   * r1*r2*c1*sqrt(2*lri+1)           
c1    &                   * (smat-kron(ni,nf))*YLMC2(lf,mlp) 
c1        lfv(lf)=lf

!        write(777,*)'lf,mlp,ylm=', lf,mlp,YLMC(lf,mlp),YLMC2(lf,mlp)

        if (debug) then
        write(klog,'(6x,"Contrib to AMP: ni,li,iex=",
     &  3i3," => nf,lf,iex=",3i3,6x,"S=",2f12.6," r1,r2=",2f12.5)')
     &  ni,li,1,nf,lf,iex,smat,r1,r2
       endif
       
!MGR--------------------------------------------------------------------
       if(ijtot.ge.0) then
       ilp=nint(ljmax+jtot-lrf)+1 
!       write(0,*) 'ijtot',ijtot,'ilp',ilp,'im',im,'imp',imp, 'jtot',jtot
       ampauxJ(ijtot,ilp,im,imp)= ampauxJ(ijtot,ilp,im,imp)  
     &                 + phci*conjg(phcf)*r1*r2*c1*sqrt(2*lri+1)  
     &                 * (smat-kron(ni,nf))*YLMC2(lf,mlp)
       endif 
!-----------------------------------------------------------------------        
       enddo ! imp
       enddo ! imp
       enddo ! nchf
       enddo ! nchi
       enddo ! icc (j/pi sets)
    

!!!!!!!!!!!!!!!!!!!!!!!!! INTERPOLATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if (.true.) then
      if (jump(1).gt.0.5d0) then ! fixed 21/nov/16
!First we need to know the points that are already calculated
      njtotcalc=0
      do ijtot=0,njtot
      if (interpJ(ijtot).eq.0) njtotcalc=njtotcalc+1
      enddo
      
!------------------------------------------------------------------------      
      allocate(gaux(njtotcalc),xv(njtotcalc))
      njtotcalc=0
      do ijtot=0,njtot
      if (interpJ(ijtot).eq.0) then
       njtotcalc=njtotcalc+1
       xv(njtotcalc)=ijtot
      endif
      enddo
!      write(0,*) 'Passed njtotcalc,', njtotcalc,'njtot',njtot
      do ijtot=0,njtot
       jtot  =ijtot+jbord(1)
       if (interpJ(ijtot).eq.0) cycle

!MGR Let's try with splines
!      pspline=(jtot-jbord(1))/jump(1)
!      call SPLINT(pspline,njtotcalc,ispline,nspline,spline)
      
!      write(0,*) 'Spline N',njtotcalc,' NP',nspline,'I',ispline,'P',
!     & pspline,'SP',spline(1:nspline)
       
       do im=1,nmi
        rm= -jpi + (im -1)
        do imp=1,nmf
        rmp= -jpf + (imp -1)
        do ilp=1,nlp

        lf=nint(ljmax+jtot-ilp)+1
        if (abs(nint(rm-rmp)).gt.lf) cycle
             
        phcf =exp(2d0*zz*(delcf(lf)-delcf(0)))
        gaux(:)=0d0
!        write(0,*) 'Interpolating Jtot:', jtot,'L',lf 
        do ijtotcalc=1,njtotcalc
          gaux(ijtotcalc)=ampauxJ(nint(xv(ijtotcalc)),ilp,im,imp)
        enddo
!MGR splines        
!        ampauxJ(ijtot,ilp,im,imp)=0d0
!        do jspl=1,nspline
!        if(interpJ(ijtot).eq.1) then
!        ampl(lf,im,imp)=ampl(lf,im,imp)+spline(jspl)*
!     &   gaux(jspl+ispline-1)*phcf
!        endif
!        ampauxJ(ijtot,ilp,im,imp)=ampauxJ(ijtot,ilp,im,imp)+spline(jspl)
!     &  *gaux(jspl+ispline-1)
!        enddo
                   
!        write(*,'(10(i3,2f10.5))') (nint(xv(l)),gaux(l),l=1,lfmax)
        caux=cfival(ijtot+0d0,xv,gaux,njtotcalc,0d0)
        ampl(lf,im,imp)= ampl(lf,im,imp)+caux*phcf
        ampauxJ(ijtot,ilp,im,imp)=ampauxJ(ijtot,ilp,im,imp)+caux
         
        enddo !ilp
       enddo ! imp
       enddo ! im
       enddo ! ijtot (j/pi sets)
       deallocate(gaux,xv)
       endif 
!       write(0,*) 'Exit interpolation'

!       do ijtot=0,njtot
!       write(888,*) ijtot, dble(ampauxJ(ijtot,1,1,1)),
!     &  aimag(ampauxJ(ijtot,1,1,1))
!       enddo
!       write(888,*) '&'    
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c write angle-indep amplitudes ----------------------------------------
       if (written(kamp)) then
       write(kamp,300) iex,exc,jpf,parity(parch(iex)+2)
       do lf=0,lfmax
!225    WRITE(kamp,226) LF,((AMPL(lf,im,imp)*phcf,imp=1,nmf), im=1,nmi)
225    WRITE(kamp,226) lf,((ampl(lf,im,imp),imp=1,nmf), im=1,nmi)
226    FORMAT(' LP=',I5,1X,1P,10E12.4 / (10X,10E12.4),/)
       do im=1,nmi
        rm= -jpi + (im -1)
        do imp=1,nmf
        rmp=-jpf + (imp-1)
        enddo !imp
       enddo !im
       enddo ! lf
       endif !write?


c Angle-dependent amplitudes: f(theta,m,mp) (Fresco CPC eq. (3.31)
       allocate(fam(nth,nmi,nmf))
!       write(*,*)'allocate pl with lfmax+1=',lfmax+1
       allocate(pl(maxleg+1,maxleg+1))
       fam(:,:,:)=0; pl=0;
       do ith=1,nth
          th = thmin + dth*(ith-1)
          if (th<theps) th=theps
          thrad = th*pi/180.
          cth=cos(thrad)
          sth=sin(thrad)
          s2= sin(thrad/2)**2
          call PLM(cth,maxleg,maxleg,maxleg+1,PL)

c Point Coulomb (Satchler 4.11a), but with different phase convention 
c (no e^{2 i delc(0)}  factor as in Fresco !) 
          r1=-etai*log(s2)  !  2*delc(0)
          fc=-etai/(2*kcmi)*exp(zz*r1)/s2
c 
       iam=0
       do im=1,nmi
       rm= -jpi + (im -1)
       do imp=1,nmf
       iam=iam+1 ! global index for spin projections
       rmp=-jpf + (imp-1)
       faux=0
       rmlp=rm-rmp
       mlp=nint(rm-rmp)
       if (abs(mlp).gt.lfmax) stop'mlp>lfmax'

       do lf=abs(mlp),lfmax
         if(mlp.lt.0) then
          r1=pl(lf+1,-mlp+1)
         else 
          r1=pl(lf+1,mlp+1)          
         endif 
        faux=faux + ampl(lf,im,imp)*r1
       enddo !lf
       fam(ith,im,imp)=fc*kron(im,imp)*kron(iexgs,iex) + faux 
       famps0(iex,ith,iam)= fam(ith,im,imp)
       enddo !im
       enddo !imp
       write(kfam,'(2x,1f8.4)') th 
!       if (ith.le.2) write(1370,*)'theta=',th
       write(kfam,228) ((fam(ith,im,imp),imp=1,nmf),im=1,nmi)
!       if (ith.le.2) write(1370,'(20g12.4)')
!     &   ((fam(ith,im,imp),imp=1,nmf),im=1,nmi)
228    format(1P,6E12.4)
       enddo !angles
!       write(1370,*)
777    format(f8.3,1e12.3,1f10.4,1e12.3)



c *** Differential cross sections for each state
       do ith=1,nth
        th = thmin + dth*(ith-1)
        if (th<theps) th=theps
        factor=1./nmi
        thrad=th*pi/180        
        xs=0d0
        do im=1,nmi
        do imp=1,nmf
         xs=xs+10*factor*abs(fam(ith,im,imp))**2
        enddo !im
        enddo !imp
        if ((iex.eq.1).and.(zp*zt.gt.0)) then 
            s2=sin(thrad/2.)**2
            xsruth=10*(etai/(2.*kcmi*s2))**2  ! point Coulomb(mb/sr)
            write(kxs,800) th,xs/xsruth,xs!,xsruth*10
            write(kxs+1,800) th,xs/xsruth,xs!,xsruth*10
        else
            write(kxs,800) th,xs
            if (iex.le.10) write(kxs+iex,800) th,xs
        endif
        enddo !ith
800     format(1x,1f8.3,2g15.5)
        write(kxs,*)'&'
        deallocate(ampl,fam,pl,delcf)
        deallocate(ampauxJ)
         deallocate(interpJ)
!        deallocate(ampaux,xv,lfv)
        if(allocated(gaux)) deallocate(gaux)
      enddo ! iex 


c *** Angle-integrated cross sections from S-matrices
      sigex(:)=0
      xstot = 0.
      xrtot = 0.
      jpi   = jpch(iexgs)
      nmi   = 2*jpi+1
      icc=0
      njt  =jtmax-jtmin+1
      do ijt=1,njt
      jtot=jtmin+dble(ijt-1)
      xsj=0
      xrj=0
      do partot=1,-1,-2
      xs=0
      xr=0
       icc=icc+1 ! CC set (total J/pi)
       if (jptset(icc)%interp) cycle
!       partot=jptset(icc)%partot
!       jtot  =jptset(icc)%jtot
       nch   =jptset(icc)%nchan
       do ni=1,nch
       if (jptset(icc)%idx(ni).ne.iexgs) cycle ! no inc waves 
       li= jptset(icc)%l(ni)
       do nf=1,nch
       iex  = jptset(icc)%idx(nf)
       exc  = exch(iex)
       tkch=ecmi+exch(iexgs)-exc
       if (tkch.lt.0) cycle    !closed channel
       kcmf = sqrt(2*mupt*tkch)/hc
       lf   = jptset(icc)%l(nf)
       jpf  = jptset(icc)%jp(nf) 
       
! Fixed July 1st 2020.        
!       if (iex.eq.iexgs) then ! elastic channel 
       if ((iex.eq.iexgs).and.(ni.eq.nf)) then ! elastic channel 
         written(ksmat)=.true.
         smat=smats(icc,ni,nf)
!         write(0,*)'jtot, ni, nf, smat=',jtot,partot,ni,nf,smat,
!     &  (pi/kcmi**2)*(2.*jtot+1)/nmi*(1-abs(smat)**2)
         xr=xr+10*(pi/kcmi**2)*(2.*jtot+1)/nmi*(1-abs(smat)**2)
         write(ksmat,'(1f6.1,2x,2i3,1f6.1," ->",2i3,1f6.1,2x,3f14.8)')
     &   jtot,iexgs,li,jpi,iex,lf,jpf,smat,abs(smat)
       else      
         smat=smats(icc,ni,nf)! *sqrt(kcmf/kcmi)
         xsaux=10*(pi/kcmi**2)*(2.*jtot+1)/nmi*abs(smat)**2*kcmf/kcmi
         xs=xs + xsaux 
         sigex(iex)=sigex(iex) + xsaux
         if (verb.gt.1) then
         write(ksmat,'(1f6.1,2x,2i3,1f6.1," ->",2i3,1f6.1,2x,3f14.8)')
     &   jtot,iexgs,li,jpi,iex,lf,jpf,smat,abs(smat)
         endif
       endif
       enddo !nf
       enddo !ni
       xstot= xstot + xs
       xrtot= xrtot + xr
       xsj  = xsj   + xs
       xrj  = xrj   + xr
       write(*,'(5x," J/pi=",1f5.1,a1,3x,
     &   " xs(inel)=",1e12.3, " mb"," reac=",1e12.3)') 
     &   jtot,parity(partot+2),xs,xr
      enddo !partot
! COMMENTED, becase we have printed the x-sections after each J/pi set
!      write(ksj, '(1x,1f5.1,3f10.4,2x,i2)') jtot,xrj-xsj,xrj,xsj
      enddo !ijt (JT)
      
      write(*,'(//,3x," INTEGRATED CROSS SECTIONS:")')
      write(*,782) xrtot
      write(*,780) xstot
      write(*,784) xrtot-xstot      
      
780   format(5x," o Total inel+bu x-section =",1f8.3," mb")
782   format(5X," o Reaction cross section  =",1f8.3," mb")
784   format(5X," o Absorption cross section=",1f8.3," mb")


c *** Angle-integrated cross section for each projectile state
      write(*,'(//,2x,"o Angle-integrated xs for each state:")')
      jpold=-1 ; parold=-2
      do iex=1,nex
        partot=parch(iex)
        exc   =exch(iex)
        jpf   =jpch(iex)
        if (iex.eq.1) jpold=jpf
        if ((iex.gt.1).and.(jpf.ne.jpold).and.(partot.ne.parold))
     & write(130,*) '&'
        jpold=jpf
        write(130,795) exc, sigex(iex), iex,jpf,partot
        write(*,790) iex,jpf,parity(partot+2),exc, sigex(iex)
790     format(8x,"#",i4,3x,1f5.1, a1,f10.3,"->",f12.6,' mb')
795     format(1x,1f10.3,2x,1f12.6,1x,i3,1f5.1,2x,i2)
      enddo

c Generate an approximate distribution dsigma/dE
c need to fix the calculation of widths!!!!!
!      call budist(exch,nex,sigex)
!      write(*,*)'call d2sigma with ner',ner


c Simple ds/de for bins 
      open(131,file='dsde_bin.xs')
      iex=0
      xstot=0
      do iset=1,jpsets
      xs=0; 
      nex=jpiset(iset)%nex
      inc=jpiset(iset)%inc
      excore =jpiset(iset)%exc(inc) ! core energy 
      do n=1,nex
      if(n.eq.1)write(131,*)'# Set',iset,'j=',jpiset(iset)%jtot
      iex=iex+1
      if (jpiset(iset)%nho.ne.0) cycle
      ki=binset(iset)%klow(n)
      kf=binset(iset)%kup(n)
      ei=hc**2/2./mu12*ki**2
      ef=hc**2/2./mu12 *kf**2
      exc=excore+ binset(iset)%ebin(n) 
      write(131,'(1x,1f10.3,2x,2f12.6,1x,i3)') 
     & exc,sigex(iex)/(ef-ei),sigex(iex),iex 
      xs=xs + sigex(iex)
      enddo ! n
      xstot= xstot + xs 
      write(131,'(5x,a,1f10.4,a)') '# -> in this set xsec=',xs,' mb'
      write(131,*)'&'
      enddo ! iset
      write(131,*)'   # Total Bins xs=',xstot,' mb'
      close(131)

      
c Double x-sections 
1000  if ((doublexs).or.(triplexs).or.(alphaxs)) then
!       write(0,*)'calling d2sigma:'
        if(.not. alphaxs) then
        call d2sigma(nth,dth,thmin,thcut,ermin,ermax,ner,
     &                 icore,jsets,fileamp,doublexs,triplexs,phixs)
        else
        write(*,*) 'Computing angle alpha'
        call d2sigma_alpha_1d(nth,dth,thmin,thcut,ermin,ermax,ner,
     &                 icore,jsets,fileamp,doublexs,triplexs,phixs)
        endif
      endif
      if (allocated(famps0)) deallocate(famps0)
      end subroutine





c ----------------------------------------------------------------
c Double differential x-sections dsigma/dEx dW  (NEW VERSION) 
c ----------------------------------------------------------------
      subroutine d2sigma(nth,dth,thmin,thcut,emin,emax,ncont,
     &                   icore,jsets,fileamp,doublexs,triplexs,phixs)
      use xcdcc,    only: elab,jpch,nex,exch,parch,famps0,binset,rel
      use channels, only: jpiset,jpsets,nchmax,sn,qnc
      use sistema
      use constants
      use parameters, only:  maxchan,maxsets,maxeset
      use globals,  only: debug,written,verb,kin
      use wfs,      only: nr,dr,energ,rvec,wfc !,wbin
      use memory,   only: t3d,lr8
      use scattering, only: method,nbas,ns
      implicit none

      include "omp_lib.h"

c     ----------------------------------------------------------------------- 
      character*40:: line,fileamp
      logical :: writewf=.false.
      logical, parameter:: energy=.true.   ! NEEDED FOR FFC4 INTERPOLATION!! 
      logical :: jsets(maxsets),doublexs,triplexs,phixs
      integer , parameter:: kfam=137  ! file containing amplitudes
      integer , parameter:: eps=1e-3
      integer , parameter:: maxstat=400
      integer :: ic,icore,ncont,idet,ier
      integer ::nsetmax,ir,iset,n,nchan,iam,ich,iecv,nth,iex,ith,mth
      integer ::nmi,nmf,nmfmax,nm,li,lf,im,imp,inc,nearfa
      integer,allocatable::ips(:,:)
c     -------------------------------------------------------------------------
      real*8:: kpt,kcv,krel,dkrdk,wbin
      real*8:: ji,jf,jci,jcf,raux,xsaux,jac,jpi,jpf,jtarg
      real*8:: ecmi,ecmf,kcmi,kcmf,excore,th,thmin,dth,thcut
      real*8:: dec,ebind,ecv,ethr,emin,emax,exc
      real*8:: facK,mv,mc,mucv,mupt,f2t,xstot,sumr
      real*8,allocatable::dsdew(:,:),dsdew_b(:,:,:)
      real*8:: delta(ncont)
      real*8:: ti,tf
c     -------------------------------------------------------------------------
      complex*16,allocatable:: wfcont(:,:,:),ampaux(:,:,:) !,gaux(:)
      complex*16,allocatable,target:: gsolap(:,:,:,:)
      complex*16 :: smate(ncont,maxchan),gaux(nr),haux(ncont)
c     -------------------------------------------------------------------------
c     Triple x-sections
c     -------------------------------------------------------------------------
      logical zert,fail3,wrt,kinset(maxsets)
      character*3 :: sys ! 'com' or 'lab'
c    ......................................................................
      integer imu,imo,ival,j,iii,ii,iflag,iroot,iroot1n,iroots,iroot2n
      integer iten,isig,ip,inu,iang,ilval,itphi,iv,itv,itphim,itc,lmax
      integer ind,ien,nchann,niag,nial,ibin,nieg,niel,ik,maxne,nk
      integer, parameter:: ncleb=2000, nord=2
c     .....................................................................
      real*8,allocatable:: angsr(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      real*8,allocatable:: xyt(:,:)
!      real*8 xyt(2,ncont)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8:: Kthr,dk,dkb,kmax,a,b
      real*8:: dot,m,plmrec,ylmc2,plm_nr
      real*8:: dphir,dtc,dtv,phiu,tcl,tvl,tcu,tvu,dtcr,dtvr
      real*8:: acc,degrad,mtot,totp
      real*8:: erel(ncont),ki,kf,ebin
      real*8:: erelmax,tcd,tc,costc,sintc,tvd,tv,costv,sintv
      real*8:: En,Enlow,Enup,dEn,p1L,p2L,phid,phil,dphi,phi,cospv,sinpv
      real*8:: cospc,sinpc,aq,bq,cq,dq,Ec,Ev,eKb,eks
      real*8:: qp(3),kcL(3),kvL(3),pcL(3),pvL(3),tem(3)
      real*8:: ptot(3),ktot(3),bkp(3),kp(3),sigphi(4000),xsig
      real*8:: mkp,co,phKb,phks,co2,si2,si,aleg,mult,xb,xbd,ybo,yb
      real*8:: tmatsq,sc,mu,sig,mo,rnu,cleb,rli
      real*8:: cplm(0:10,21),cgc(ncleb)
c     ..........................................................................
      complex*16 :: resc,caux,tmat_aux,faux,cfival,ffc4
      complex*16:: ampt,xi,f2c,sumps,phbin(maxsets),wkfac(maxsets)
      complex*16:: phask1,phask,phas,ylm(0:10,21)
      complex*16,allocatable:: tmat(:,:,:)
      complex*16,pointer:: overlap(:)
      complex*16,allocatable:: fxyc(:,:)
!      complex*16:: fxyc(10,ncont)
      complex*16,allocatable:: fv(:)
c  *** Relativistic kinematics
! included now in xcdcc module
!      logical  :: rel=.true.
      real*8   :: sinv,etcm,pcmi,pcmf,mpx
c     ..........................................................................

!!! TESTING
      real*8 tmatsq2,ampt2
      integer icount,if2c
      complex*16,allocatable:: fin(:,:,:)
      complex*16:: bindk(maxsets,maxeset)
      real*8:: xtab(6),ytab(6),kcmf_cont
      complex*16:: fxytab(6,6),fint2d,fint2dd,fint2db_jin

!MGR---------------------------------------------------------------------
      real*8,allocatable :: xs3body(:,:,:), energ3(:)
! AMM July 19
      real*8,allocatable :: xs3body_phi(:,:,:,:), energ4(:,:,:,:)
!------------------------------------------------------------------------
c ------------------------------------------------------------------------------
      namelist /framework/ sys,idet ! ,atarget not needed
      namelist /gridener/ Enlow,Enup,dEn
      namelist /gridthetac/ tcl,tcu,dtc
      namelist /gridthetav/ tvl,tvu,dtv
      namelist /gridphi/ phil,phiu,dphi
c ------------------------------------------------------------------------------
      nsetmax=maxval(jpiset(:)%nex) !max number of states within j/pi set
      allocate(fv(nsetmax))
      nmfmax =2*nint(maxval(jpiset(:)%jtot))+1
      mc=ac*amu
      mv=av*amu
      mt=mt*amu
      mp=mp*amu  
      mupt=mp*mt/(mp+mt)
      mucv=mc*mv/(mc+mv)
      nmi=nint(2*jpiset(1)%jtot+1)

      write(*,'(//,3x, "** DOUBLE DIFFERENTIAL CROSS SECTIONS **",/ )')

c     If no amplitudes have been calculated before, we try to read them
c     from an external file, specified by the user 
      if (.not.allocated(famps0).and.(fileamp.ne."")) then
      if (fileamp.eq."") then
          write(*,*)'No amplitudes calculated and no file specified';
          stop
      endif
      write(*,'(3x," Reading f-amps from file ",a20)')fileamp
      open(kfam,file=fileamp,status='old')
      read(kfam,'(a)',err=900) line 
!      write(*,'(a)') line
      read(line,*) jpi,jtarg,jpf,jtarg,mth,nearfa,elab     
!      read(kfam,*,err=900) jpi,jtarg,jpf,jtarg,mth,nearfa,elab
!     & ,a,b
!        write(*,*) 'ie,jpf,mth,elab,a=',iex,jpf,mth,elab,a
      rewind(kfam)

      allocate(famps0(maxstat,mth,nmi*nmfmax))
      iex=0
      famps0(:,:,:)=0
      do iset=1,jpsets  
      nex=jpiset(iset)%nex
      jpf=jpiset(iset)%jtot
      nmf=nint(2.*jpf+1.)
      nm=nmi*nmf
      write(*,'(a,i3,a,i3,a,1f4.1)') ' Expected: iset=',iset, 
     & '  nex=', nex, ' jpf=',jpf
      do n=1,nex
      iex=iex+1
      if (iex.gt.maxstat) stop'increase dimension of famps0' ! temporary solution
!      write(*,*)'expected jpf=',jpf
      read(kfam,'(a)',err=900) line
      read(line,*) jpi,jtarg,jpf,jtarg,mth,nearfa,elab    
!      read(kfam,*) jpi,jtarg,jpf,jtarg,mth,nearfa,elab
!     & ,a,b
      write(*,'(a,i3,a,i3,a,1f4.1)') ' Found: ie=',iex, 
     & ' n=',n,' jpf=',jpf

!      if (iex.eq.1) write(*,*) ' amplitudes are for elab=',elab
      if(nth.ne.mth) then
         write(*,*) 'mismatch in number of angles: specified:',nth,
     &              ' found:',mth
      endif 
      do ith=1,nth
      read(kfam,*) th  
      read(kfam,228) (famps0(iex,ith,iam),iam=1,nm) 
228   format(1P,6E12.4)
      enddo !states
      enddo !ith
      write(*,*)' ' 
      enddo !iset
650   write(*,'(5x,"[f amplitudes read for ",i3," states and",i5,
     & " angles ]")')iex,nth
      endif
c     -------------------------------------------------------------------------
 

      if (.not.rel) then       ! non-relativistic
        ecmi=elab*mt/(mp+mt)
        mupt=mp*mt/(mp+mt)!*amu
        kcmi=sqrt(2*mupt*ecmi)/hc
        write(*,*)' ** USING NON-RELATIVISTIC KINEMATICS **'
      else                     ! relativistic
        write(*,*)' ** USING RELATIVISTIC KINEMATICS **'
        sinv=((mp+mt))**2+2*mt*elab
        etcm=sqrt(sinv)          ! total relativistic energy
        ecmi=etcm-mp-mt       ! kinetic energy in initial channel
        Pcmi=sqrt((sinv-mp**2-mt**2)**2-4.*(mp*mt)**2)/2/sqrt(sinv)
        Kcmi=pcmi/hc
        write(*,'(3x,a,1f10.3,a)')
     & 'Relativistic invariant mass: sqrt(s)=',
     & sqrt(sinv), ' MeV'
        write(*,'(3x,"Ecmi=",1f8.3, "MeV   Kcmi=",1f8.3," fm-1" )') 
     &  Ecmi, kcmi
      endif

      ebind=energ(1,1) ! exch(1)          
!      Kcm=sqrt(2.d0*mupt*ecm)/hc        
      Ethr=Ecmi+ebind
      Kthr=sqrt(2.d0*mupt*Ethr)/hc

!      write(99,*) 'Ethr=',Ethr, ' Kthr=',Kthr
!      write(99,'(a,1f8.3,a,3f10.5)')'Elab=',elab,
!     & ' mp,mc,mv=',mp/amu,mc/amu,mv/amu
!      write(99,*)'Ecmi=',Ecmi, "Kcm=",kcmi, 'Ethr=',Ethr


c *** Assign global index to states of all j/pi sets iPS=iPS(iset,n)
      iex=0
      allocate(ips(jpsets,nsetmax))
      do iset=1,jpsets
      nex=jpiset(iset)%nex
      do n=1,nex
      iex=iex+1
      ips(iset,n)=iex
      if ((verb.ge.3).and.(iset.gt.1).and.(n.eq.1)) write(*,*)
      if (verb.ge.3) write(*,*)'Set:',iset,' state=',n,"-> iex=",iex
      enddo
      enddo


c *** Sum discrete breakup angular distributions    
      if (any(jpiset(1:jpsets)%nho>0)) then  
         allocate(dsdew_b(jpsets,nsetmax,nth))
         dsdew_b(:,:,:)=0
      endif
      
      open(90,file='dsdw_bu.xs',status='unknown')
      do ith=1,nth
      xstot=0     
      do iset=1,jpsets
!      if (jsets(iset)) write(*,*)'include j/pi set', iset
      if (.not.jsets(iset)) cycle
!      if(jpiset(iset)%nho>0) cycle ! not bins (PS)
!      if (.not.allocated(dsdew_b)) allocate(dsdew_b(jpsets,nsetmax,nth))
      nex=jpiset(iset)%nex
      jpf=jpiset(iset)%jtot
      nmf=nint(2*jpf+1)
      do n=1,nex  
      xsaux=0
      iex=ips(iset,n)
      
      if (energ(iset,n).lt.0)    cycle  ! bound state 
      if (energ(iset,n).gt.Ethr) cycle  ! closed channel
      iam=0
      raux=0
      do im=1,nmi
      do imp=1,nmf
      iam=iam+1  
      raux=10.*abs(famps0(iex,ith,iam))**2/nmi
      xsaux = xsaux  + raux 
      xstot = xstot  + raux
      enddo !imp 
      enddo !im
      
      if(jpiset(iset)%nho.eq.0)then
          wbin=binset(iset)%wbin(n)
          dsdew_b(iset,n,ith)=  dsdew_b(iset,n,ith) + xsaux/wbin
      endif
      enddo ! states within j/pi set
      enddo !j/pi sets
      write(90,*) thmin + dth*(ith-1),xstot
      enddo !ith
      call flush(90)
      close(90)
c ----------------------------------------------------------------------


c Print ds/dE for bins !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (any(jpiset(1:jpsets)%nho>0)) then
      open(92,file='dsde_bin.xs')
      write(*,*)' ds/dE for bins sets: thcut=',thcut
      do iset=1,jpsets
      if(jpiset(iset)%nho>0) cycle ! not bins (PS)
      nex=jpiset(iset)%nex
      do n=1,nex    
      iex=ips(iset,n)
      exc  = energ(iset,n)
      xsaux=0
	  do ith=1,nth
      raux=0.
      th = thmin + dth*(ith-1)
      if (th.gt.thcut) cycle      
      xsaux=xsaux + dsdew_b(iset,n,ith)*2*pi*sin(th*pi/180.)*dth*pi/180.   
      enddo !ith
      write(92,*) exc,xsaux,iex
      enddo !n within set
      write(92,*)'&'
	  enddo ! iset
      endif ! are there bin sets?



      
c *** Overlaps between PS's and Scattering states  -----------------------------
      allocate(gsolap(jpsets,nsetmax,maxchan,ncont))
      gsolap(:,:,:,:)=0
      if (emax<0)     emax=maxval(energ(:,:))
      if (emax>Ethr)  emax=Ethr-0.001
      if (emin.le.0.) emin=0.01
      dec=(emax-emin)/dble(ncont-1)
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      if (jpiset(iset)%nho.eq.0) cycle
      nchan=jpiset(iset)%nchan
      nex =jpiset(iset)%nex
      do inc=1,nchan 
      if (inc.gt.maxchan) stop'inc too large in d2sigma!'
      jci    =jpiset(iset)%jc(inc)
      li     =jpiset(iset)%lsp(inc)
      excore =jpiset(iset)%exc(inc) ! core energy
      ic     =jpiset(iset)%cindex(inc)
      if (ic.ne.icore) then
         if (ic.eq.abs(icore)) then
         call solapread(gsolap(iset,:,inc,:),jpiset(iset)%partot,
     & inc,nsetmax,ncont,jpiset(iset)%jtot,dec,emin+excore)
         else if (verb.ge.2) then
         write(*,'(4x,"-> skipping core state",i3)')ic 
         endif
      cycle
      endif
      
!          write(0,*)'method=',method
      method=4;
      nbas  =50
      ns    =1; 

      if (allocated(wfcont)) deallocate(wfcont)
      if (ncont.gt.0) then
        allocate(wfcont(ncont,nchan,nr)); wfcont=0.
        call wfrange(iset,nchan,inc,emin,emax,ncont,energy,wfcont,
     & smate,delta,writewf)
      else
       write(*,*)'*ERROR* ncont,nchan=',ncont,nchan
      endif

!!!! CHECK Continuum WFS
!      write(505,*)'#WF for Ecm;J=',emin,jpiset(iset)%jtot
!      do ir=1,nr
!         write(505,111)rvec(ir),wfcont(1,1,ir)
!      enddo
!!!!!!!!!!!!!!!!!!!!!!!!

      do n=1,nex
      if (energ(iset,n).lt.0)    cycle ! bound state
      if (energ(iset,n).gt.Ethr) cycle ! closed channel
         
      erel(:)=0
      do iecv=1,ncont
      ecv=emin+(iecv-1)*dec
      erel(iecv)=ecv
      if (ecv.lt.eps) ecv=eps
!      kcv=sqrt(2.d0*mucv*ecv)/hc    
      sumr=0
      caux=0.
      do ich=1,nchan ! sum in final channels
      jcf =jpiset(iset)%jc(ich)
      lf  =jpiset(iset)%lsp(ich)
      jf  =jpiset(iset)%jc(ich)
      do ir=1,nr   
      gaux(ir)=wfc(iset,n,ich,ir)*wfcont(iecv,ich,ir)*rvec(ir)*
     &         (-1)**(li+jci+lf+jcf) ! AMM: I think this phase is always 1! 
      enddo !ir

      call simc(gaux,resc,1,nr,dr,nr)
      caux=caux+resc
      sumr = sumr+abs(caux)**2
      enddo ! ich 
      gsolap(iset,n,inc,iecv)=caux
      enddo ! iecv (c-v relative energy)
      enddo ! n  (PS within this j/pi set)
      enddo ! inc  (incoming channel)
      enddo ! iset (j/pi set)
      icore=abs(icore)
     
c *** -------------- PRINT OVERLAPS FOR TESTING ----------------------------------
      if (verb.ge.3) then
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nex  =jpiset(iset)%nex
      nchan=jpiset(iset)%nchan
      if (jpiset(iset)%nho.eq.0) cycle
      do n=1,nex
      raux=0
      write(97,'(a,i3,a,i3,a,1f10.4)')'# Set:',iset,
     & ' n=',n,' Ex=',energ(iset,n)
      do iecv=1,ncont
       ecv=emin+(iecv-1)*dec
       if (energ(iset,n).lt.0)    cycle ! bound state
       if (energ(iset,n).gt.Ethr) cycle ! closed channel
       if (ecv.lt.eps) ecv=eps
       kcv=sqrt(2.d0*mucv*ecv)/hc
       jac=mucv/(hc**2)/kcv
       write(97,111)ecv,
     &  (jac*abs(gsolap(iset,n,inc,iecv))**2,inc=1,nchan)
111    format(2x,1f12.6,2x,10g14.6)
       do inc=1,nchan
       raux=raux+jac*abs(gsolap(iset,n,inc,iecv))**2*dec
       enddo ! inc
      enddo !iecv
      if (verb.ge.3) then
      write(97,'(3x,a,i3,a,i3,a,1f8.5)')'# -> Set:',iset, '  PS=',n,
     & ' Norm solap=',raux*2/pi
      endif
      write(97,*)'&'
      enddo ! n
      enddo !iset
      endif ! verb
c-------------------------------------------------------------------------------


c *** Compute DOUBLE x-sections dsigma/(dEdW) FOR SELECTED CORE STATE ICORE
      if (ncont.gt.0) then
        allocate(dsdew(ncont,nth))
      else  
         write(*,*)'*ERROR* ncont,nchan=',ncont,nchan
      endif
        
      dsdew=0
      raux=0.5/(pi**3)/hc**6*mupt**2*mucv/nmi
      raux=raux*10 ! fm^2 -> mb

      excore=qnc(icore)%exc
      write(*,'(3x,a,i3,a,1f7.3,a)')
     & '=> Double diff. x-sections for final core state',
     &  icore,' with energy=',excore, ' MeV'

      do iecv=1,ncont
      ecv=emin+(iecv-1)*dec
      if (ecv.lt.1e-3) ecv=1e-3
      kcv=sqrt(2.d0*mucv*ecv)/hc  ! above threshold!!!
      do iset=1,jpsets
      if (jpiset(iset)%nho.eq.0) cycle ! bins
!      if (jsets(iset)) write(*,*) 'include j/pi set',iset
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      nex  =jpiset(iset)%nex
      jpf  =jpiset(iset)%jtot
      nmf  =nint(2.*jpf+1.)
      iam=0
      do im =1,nmi
      do imp=1,nmf
      iam=iam+1
      nchan=jpiset(iset)%nchan
      do inc=1,nchan  ! sum over open channels
      ic  =jpiset(iset)%cindex(inc)
      if (ic.ne.icore) cycle
      li     =jpiset(iset)%lsp(inc)
      ji     =jpiset(iset)%jsp(inc)
      jci    =jpiset(iset)%jc(inc)
      excore =jpiset(iset)%exc(inc) ! core energy 
      ic     =jpiset(iset)%cindex(inc)

      ecmf   =ecmi-dabs(ebind)-ecv-excore ! p-t c.m. energy in final channel !! CHEEEEEEEEEEEEEEEEECK
      if (ecmf.lt.0) then
          write(*,*)'Ecm=',Ecmf, 'Ecv=',Ecv
          stop
      endif
      if (ecmf.lt.1e-4) ecmf=1e-4

      if (.not.rel) then ! non-relativistic
        Kcmf=sqrt(2.d0*mupt*ecmf)/hc   
        if ((ith.eq.1).and.(iecv.eq.1))
     &     write(*,*)'Kcmf(nr)=',kcmf    
      else               ! relativistic
        mpx=mp+abs(ebind)+ecv+excore
        Pcmf=sqrt((sinv-mpx**2-mt**2)**2-4.*(mpx*mt)**2)/2/sqrt(sinv)
        Kcmf=pcmf/hc
!        if ((ith.eq.1).and.(iecv.eq.1))
!         write(*,*)'Kcmf(rel)=',kcmf,sqrt(2.d0*mupt*ecmf)/hc
      endif

     

      Kcmf=sqrt(2.d0*mupt*ecmf)/hc 
      facK=kcmf/kcmi/kcv
! commented v26
!      f2t=-2.*pi*hc**2/mupt*sqrt(Kcmi/Kcmf) ! conversion factor f(theta) -> T-mat
!
      do ith=1,nth  
      sumps=0d0
      do n=1,nex ! PS states within this j/pi set
      if (energ(iset,n).lt.0)    cycle ! bound states
      if (energ(iset,n).gt.Ethr) cycle ! closed channels
      iex=ips(iset,n)
!      if((iecv.eq.1).and.(imp.eq.1)) 
!     &   write(*,*)'inc,iset,n,iex=',inc,iset,n,iex
c v26 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ecmf=ecmi-dabs(ebind)-energ(iset,n)-excore
      Kcmf=sqrt(2.d0*mupt*ecmf)/hc   
      f2t=-2.*pi*hc**2/mupt*sqrt(Kcmi/Kcmf) ! conversion factor f(theta) -> T-mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      tmat_aux=f2t*famps0(iex,ith,iam)
      sumps=sumps+gsolap(iset,n,inc,iecv)*tmat_aux
!      if (ith.eq.1) write(*,*) iset,n,inc,iecv,
!     & abs(gsolap(iset,n,inc,iecv))
      enddo ! PS's
      dsdew(iecv,ith)= dsdew(iecv,ith)
     &               + raux*facK*abs(sumps)**2
      enddo ! ith
      enddo ! inc
      enddo ! imf
      enddo ! imi
      enddo ! j/pi 
      enddo ! iecv (valence-core relative energy)  


c Integrate in ENERGY, to get dsigma/dOmega for core state ICORE
      open(90,file='dsdw_conv.xs',status='unknown')
      open(91,file='dsdwe_conv.xs',status='unknown')
      open(92,file='dsdwe_conv.gnu',status='unknown')
      write(91,'("nel ", i5,2x, "ang ", i5)') ncont,nth 
      do ith=1,nth
      raux=0.
      th = thmin + dth*(ith-1)
      write(91,'(a,1f8.4)') '#thcm=',th
      do iecv=1,ncont
      ecv=emin+(iecv-1)*dec
      raux=raux + dsdew(iecv,ith)*dec
      write(91,'(1f10.4,1g14.6)') ecv, dsdew(iecv,ith)      
      write(92,'(2f10.4,1g14.6)') th,ecv, dsdew(iecv,ith)      
      enddo !iecv
      write(92,*)' ' 
      write(90,*) th, raux
      enddo !ith
      call flush(90)
      close(90); close(91); close(92)



c Integrate in ANGLE, to get dsigma/dE  for core state ICORE
      open(93,file='dsde_conv.xs',status='unknown')
      xstot=0
      do iecv=1,ncont
      ecv=emin+(iecv-1)*dec
      raux=0.
      do ith=1,nth
      th = thmin + dth*(ith-1)
      if (th.gt.thcut) cycle
      raux=raux + 2*pi*dsdew(iecv,ith)*
     &     sin(th*pi/180.)*dth*pi/180.   
      enddo !iecv
      write(93,*) ecv +excore, raux
      xstot=xstot + raux*dec
      enddo !ith
      write(*,*)' Integrated dsigma/dEx=',xstot
      call flush(93)
      close(93)
      
      

c---------------------------------------------------------------------------------
c     Start calculation of triple differential cross sections 
c     (inspired in cdcn code by JAT)
c---------------------------------------------------------------------------------
      if (.not.triplexs) return
      write(*,'(//,3x, "** TRIPLE DIFFERENTIAL CROSS SECTIONS **",/ )')

!      open(99,file='thox.log')
      
      read(kin,nml=framework)
      write(77,'(a)') sys
      write(77,316) mt/amu,zt,elab/(ac+av)
      write(77,316) ac,zc
      write(77,316) av,abs(ebind)
      write(77,316) elab
      write(77,*) idet
      
      if (phixs) then
      write(777,'(a)') sys
      write(777,316) mt/amu,zt,elab/(ac+av)
      write(777,316) ac,zc
      write(777,316) av,abs(ebind)
      write(777,316) elab
      write(777,*) idet
      end if      
      
*     ---------------------------------------------------------------
      call cpu_time(ti)
      xi=(0.d0,1.d0)
      acc=1.d-6
      degrad=pi/180.d0
      mtot=mc+mv+mt
!!!! REDEFINE Mp as IN TESTN 
      mp=mc+mv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      excore=qnc(icore)%exc
      write(*,'(5x,a,i3,a,1f7.3,a)')
     & '=> Triple x-sections for final core state',
     &  icore,' with energy=',excore, ' MeV'


      allocate(angsr(nth))
      nchann=maxval(jpiset(:)%nchan)
      if (allocated(tmat)) deallocate(tmat)
      allocate(tmat(jpsets,nchann,nmi*nmfmax))
      allocate(fin(jpsets,nchann,nmi*nmfmax)) ! TESTING

*     --------------------------------------------------------------- 
      do iang=1,nth
      angsr(iang)=thmin+(iang-1)*dth
      angsr(iang)=angsr(iang)*degrad
      enddo
*     ---------------------------------------------------------------
      read(kin,nml=gridener)
      read(kin,nml=gridthetac)
      read(kin,nml=gridthetav)
      read(kin,nml=gridphi)
*     ---------------------------------------------------------------
      iten=nint((Enup-Enlow)/dEn)+1
      dEn=(Enup-Enlow)/(iten-1)
*     ---------------------------------------------------------------
*     particle with specified energy defines phi=0 degrees
*     ---------------------------------------------------------------
      if(idet.eq.1) then
      cospc=1.d0
      sinpc=0.d0
      else
      cospv=1.d0
      sinpv=0.d0
      endif
*     ---------------------------------------------------------------           
      itc=nint((tcu-tcl)/dtc)+1
      dtc=(tcu-tcl)/(itc-1)
      dtcr=dtc*degrad
*     ---------------------------------------------------------------           
      itv=nint((tvu-tvl)/dtv)+1
      dtv=(tvu-tvl)/(itv-1)
      dtvr=dtv*degrad
*     ---------------------------------------------------------------           
      itphi=nint((phiu-phil)/dphi)+1
      itphim=itphi/2+1 ! AMM
      dphi=(phiu-phil)/(itphi-1)
      dphir=dphi*degrad
*     ---------------------------------------------------------------
      write(77,307) itc,dtcr,tcl,tcu,dtc
      write(77,307) itv,dtvr,tvl,tvu,dtv
      write(77,307) iten,dEn,Enlow,Enup
      write(77,307) itphi,dphir,phil,phiu,dphi
      
      if (phixs) then
      write(777,307) itc,dtcr,tcl,tcu,dtc
      write(777,307) itv,dtvr,tvl,tvu,dtv
      write(777,307) iten,dEn,Enlow,Enup
      write(777,307) itphim,dphir,phil,phil+(itphim-1)*dphi,dphi
      end if
      
*     ---------------------------------------------------------------
*     specify total incident momentum (momentum of c.m.) in MeV/c
      if(sys.eq.'lab') totp=sqrt(2.d0*mp*elab)
*     if detector positions refer to c.m. frame 
      if(sys.eq.'com') totp=0.d0
*     ---------------------------------------------------------------
*     set up projectile cm K vector (qp) (in z-direction)
*     set up ptot vector and ktot (also in z-direction)             
      do j=1,2      
      qp(j)=0.d0
      ptot(j)=0.d0       
      ktot(j)=0.d0
      enddo
      qp(3)=Kcmi
      ptot(3)=totp
      ktot(3)=totp/hc
      erelmax=0.d0
!      write(99,*)'mp,elab,ptot(3)=',mp/amu,elab,totp
*     -------------------------------------------------------------- 
*     compute and store CG coefficients
*     -------------------------------------------------------------- 
      sc=qnc(icore)%jc
      ind=0
      do imu=1,nint(2*sc+1)
      mu=imu-sc-1
      do isig=1,nint(2*sn+1)
      sig=isig-sn-1
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      jpf  =jpiset(iset)%jtot
      do inc=1,nchan
      li     =jpiset(iset)%lsp(inc)
      ji     =jpiset(iset)%jsp(inc)
      jci    =jpiset(iset)%jc(inc)
      ic  =jpiset(iset)%cindex(inc)
      if (ic.ne.icore) cycle
      rli=li
      if (fail3(rli,sn,ji)) cycle
      if (fail3(ji,jci,jpf)) cycle
      do inu=1,2*li+1
      rnu=inu-li-1
      if (abs(rnu+sig).gt.ji.or.abs(rnu+sig+mu).gt.jpf) cycle
      ind=ind+1
      cgc(ind)=cleb(rli,rnu,sn,sig,ji,rnu+sig)*
     .cleb(ji,rnu+sig,jci,mu,jpf,rnu+sig+mu)
      enddo
      enddo
      enddo
      enddo
      enddo
      if(ind.gt.ncleb) then
      write(*,*)'increase ncleb up to',ind
      stop 
      endif
*     --------------------------------------------------------------- 
*     compute largest orbital angular momentum (lmax), maximum relative 
*     momentum (kmax), and maximum nb of states per set (maxne)
*     ---------------------------------------------------------------
!      write(99,'(/,4(a,1f9.5,4x),/)')'Ecmi=',Ecmi," Kcmi=",Kcmi, 
!     &"Kthr=",Kthr, 'Ethr=',Ethr

      lmax=0
      kmax=0
      maxne=0
!      write(*,*)'energ(1,1)=',energ(1,1)
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      nex=jpiset(iset)%nex
      jpf=jpiset(iset)%jtot
      nmf=nint(2.*jpf+1.)
      iam=nmi*nmf 
      if (nex.gt.maxne) maxne=nex
      do n=1,nex     
      ebin=energ(iset,n) 
      if (ebin.lt.0) cycle ! bound state
      iex=ips(iset,n)
!      write(99,*)'iset,n=',iset,n
!      write(99,308) (famps0(iex,1,iii),iii=1,iam)
308   format(1p,6e12.4)
      if (jpiset(iset)%nho.eq.0) then ! bin set
!        ebin=hc**2*(binset(iset)%khat(n))**2/(2.d0*mucv)
!        energ(iset,n)=ebin
        if (binset(iset)%kmax.gt.kmax) 
     &  kmax=binset(iset)%kmax
      else                            ! PS set
        kmax=sqrt(2*mucv*emax)/hc
      endif
      ecmf=ecmi-abs(ebind)-ebin
      Kcmf=sqrt(2.d0*mupt*ecmf)/hc    
      f2t=-2.*pi*hc**2/mupt*sqrt(Kcmi/Kcmf) ! conversion factor f(theta) -> T-mat
      famps0(iex,:,:)=f2t*famps0(iex,:,:)      ! convert f(theta) -> T(theta)
      if (ecmf.lt.1e-4) ecmf=1e-4
 3130 format(a,i3,a,i3,a,1f11.6,a,2f11.5,a,1f11.5,a,1f11.5)
      enddo ! n
      do inc=1,nchan 
      if(jpiset(iset)%lsp(inc).gt.lmax) then
      lmax=jpiset(iset)%lsp(inc)
      endif
      enddo !inc
      enddo !jpsets

      write(*,'(3x,a,i3,a)') 'There are',jpsets,' sets with:'
      write(*,'(5x,a,1f7.3, a,1f8.3,a)') 'o Maximum relative momentum=',
     & kmax, ' fm-1 and relative energy:',hc**2*kmax**2/2/mucv,' MeV'
      write(*,'(5x,a,i3)') 'o Maximum orbital ang. momentum=',lmax
      write(*,'(5x,a,i3)') 'o Max. number of states=',maxne

      if (allocated(xyt)) deallocate(xyt)
      allocate(xyt(2,maxne))
      if (allocated(fxyc)) deallocate(fxyc)
      allocate(fxyc(10,maxne))

*     ---------------------------------------------------------------
*     calculate and store coefficients for spherical harmonics
*     ---------------------------------------------------------------
      do ilval=0,lmax
      do inu=1,2*ilval+1
      cplm(ilval,inu)=ylmc2(ilval,inu-ilval-1)
      enddo
      enddo
*     ---------------------------------------------------------------
*     loop over core particle detection angle
      wrt=.false. ; icount=0; if2c=0 !MGR
      
      write(0,*) '- Triple diff xsections need ',
     &           itc*itv*iten*lr8/1e6,'MB' !MGR  
      allocate(xs3body(iten,itv,itc),energ3(iten))
!      allocate(energ4(itphim*itc*itv*iten))

      xs3body=0d0
!AMM 
      if (phixs) then
       itphim=itphi/2+1
       write(0,*) '- Triple diff xsections WITHOUT phi int. need ',
     &           itc*itv*iten*itphim*lr8/1e6,'MB' !AMM
      allocate(xs3body_phi(iten,itv,itc,itphim))
      allocate(energ4(iten,itv,itc,itphim))
      xs3body_phi=0.0
      end if
      
      do ien=1,iten
        energ3(ien)=Enlow+(ien-1)*dEn
      enddo
!      write(0,*) 'Ener', energ3(:)
    
!$OMP PARALLEL DO FIRSTPRIVATE(tcd,tc,costc,sintc,iv,tvd,tv,costv,
!$OMP& sintv,itphim,iflag,zert,ien,En,p1L,p2L,ip,sigphi,phid,phi,cospv,
!$OMP& sinpv,cospc,sinpc,pcL,pvL,tem,aq,bq,cq,dq,iroots,iroot,iroot1n,
!$OMP& iroot2n,Ec,Ev,kvL,kcL,bkp,kp,ekb,eks,
!$OMP& tmatsq,tmatsq2,ybo,xb,xbd,yb,phkb,xyt,iang,nial,mkp,co,phks,co2,
!$OMP& si2,si,phask1,phask,phas,ilval,j,inu,aleg,ylm,fxyc,tmat,
!$OMP& fin,iset,nchan,jpf,nex,nmf,iam,inc,fv,n,iii,ic,excore,ier,nieg,
!$OMP& niel,iecv,ecv,ecmf,kcmf,th,sumps,iex,ibin,ki,kf,dk,dkb,ik,phbin,
!$OMP& if2c,icount,imo,mo,ind,imu,mu,isig,sig,ampt,ampt2,li,ji,jci,rli,
!$OMP& rnu,mult,xsig,nk,kinset,wkfac,haux,fxytab) PRIVATE(ii) 

!      write(*,*)'itc,itv,iten=',itc,itv,iten
      
      do 10 ii=1,itc
      tcd=tcl+(ii-1)*dtc
      tc=tcd*degrad
      costc=cos(tc)
      sintc=sin(tc)
*     ---------------------------------------------------------------
*     loop over valence particle detection angle
      do 20 iv=1,itv
      tvd=tvl+(iv-1)*dtv
      tv=tvd*degrad
      costv=cos(tv)
      sintv=sin(tv)
*     ---------------------------------------------------------------
*     check if either theta is zero since no need to do multiple phi
*     calculations (if wanted, itphi>1) in this geometry
*     ---------------------------------------------------------------
      itphim=itphi/2+1
      iflag=0
      zert=(abs(tcd).le.acc.or.abs(tvd).le.acc)
      if(itphi.gt.1.and.zert) then
      itphim=1
      iflag=1
      endif
*     ---------------------------------------------------------------
*     loop over detected energy of chosen particle
      do 30 ien=1,iten
      En=Enlow+(ien-1)*dEn
      if(idet.eq.1) then
      p1L=sqrt(2.d0*mc*En)
      else
      p2L=sqrt(2.d0*mv*En)       
      endif
*     ---------------------------------------------------------------
*     the detected particle defines phi=0 degrees
*     ---------------------------------------------------------------
*     loop over azimuth detection angle          
      do 40 ip=1,itphim
      icount= icount+1
      sigphi(ip)=0.d0
      phid=phil+(ip-1)*dphi
      phi=phid*degrad  
      if(idet.eq.1) then
      cospv=cos(phi)
      sinpv=sin(phi)
      else
      cospc=cos(phi)
      sinpc=sin(phi)
      endif
*     -----------------------------------------------------------------    
*     core particle momentum UNIT vector in chosen frame             
      pcL(1)=sintc*cospc
      pcL(2)=sintc*sinpc
      pcL(3)=costc
*     valence particle momentum UNIT vector in chosen frame
      pvL(1)=sintv*cospv
      pvL(2)=sintv*sinpv
      pvL(3)=costv
      do j=1,3
      if(idet.eq.1) then
      pcL(j)=p1L*pcL(j)
      else
      pvL(j)=p2L*pvL(j)
      endif
      enddo
!      if (wrt)write(99,'("tcd=",20g12.5)') tcd,tvd,en,phid
!      if (wrt)write(99,'(/,"iv,ii,ien=",4i3," pcL=",3g14.6,/)') 
!     & iv,ii,ien,ip,pcL(1:3)! , pvL(1:3)
*     -----------------------------------------------------------------
*     calculate second particle energy and momentum depending on idet
      do j=1,3
      if(idet.eq.1) then              
      tem(j)=ptot(j)-pcL(j)
      else
      tem(j)=ptot(j)-pvL(j)        
      endif
      enddo
!      if (wrt)write(99,'(/,"tem=",3g16.6)') tem(1:3)
*     coefficients in quadratic for second particle momentum
      if(idet.eq.1) then
      aq=(1.d0/mv+1.d0/mt)/2.d0
      bq=-dot(tem,pvL)/mt
      cq=(dot(tem,tem)/mt+dot(pcL,pcL)/mc)/2.d0
      else
      aq=(1.d0/mc+1.d0/mt)/2.d0
      bq=-dot(tem,pcL)/mt
      cq=(dot(tem,tem)/mt+dot(pvL,pvL)/mv)/2.d0
      endif
!      if (wrt)write(99,*)'a,b,c=',aq,bq,cq
      if(sys.eq.'lab') cq=cq-(elab+ebind-excore) !!! CHECK CHANGE WITH JAT CODE!!!!!!!!!!!!!!!!!!!!
      if(sys.eq.'com') cq=cq-(ecmi+ebind-excore)  !!! CHECK CHANGE WITH JAT CODE!!!!!!!!!!!!!!!!!!!!
*     -----------------------------------------------------------------
*     discriminant of the quadratic equation
*     -----------------------------------------------------------------
      dq=bq*bq-4.d0*aq*cq
*     -----------------------------------------------------------------
*     possible there are no roots for light targets and some angles
*     and energies in the lab frame
*     -----------------------------------------------------------------
      iroots=0
      if(dq.lt.0.d0) goto 52
*     -----------------------------------------------------------------
*     otherwise proceeed to the two roots and check for positive values
*     -----------------------------------------------------------------
      iroots=2
      iroot=1
*     -----------------------------------------------------------------
*     finalise momentum vector of second particle
*     -----------------------------------------------------------------
      iroot1n=0
      if(idet.eq.1) then
      p2L=(-bq+sqrt(dq))/(2.d0*aq)
      if(p2L.le.0.d0) iroot1n=1
      else
      p1L=(-bq+sqrt(dq))/(2.d0*aq)          
      if(p1L.le.0.d0) iroot1n=1
      endif
      if(iroot1n.gt.0) iroots=1
*     -----------------------------------------------------------------
53    if(iroot.eq.2.or.iroot1n.gt.0) then
      iroot2n=0
      if(idet.eq.1) then
      p2L=(-bq-sqrt(dq))/(2.d0*aq)
      if(p2L.le.0.d0) iroot2n=1
      else
      p1L=(-bq-sqrt(dq))/(2.d0*aq)          
      if(p1L.le.0.d0) iroot2n=1
      endif
      if(iroot2n.gt.0) then
      iroots=iroots-1
      goto 52
      endif
      endif
*     -----------------------------------------------------------------
*     core and valence particle energies and momenta in chosen frame
*     -----------------------------------------------------------------
      Ec=(p1L*p1L)/(2.d0*mc)
      Ev=(p2L*p2L)/(2.d0*mv)
      if ((phixs).and.(idet.eq.1)) energ4(ien,iv,ii,ip)=Ev
      if ((phixs).and.(idet.eq.2)) energ4(ien,iv,ii,ip)=Ec

      
      if(idet.eq.1) then
      do j=1,3
      kvL(j)=p2L*pvL(j)/hc                 
      kcL(j)=pcL(j)/hc
      enddo                  
      else
      do j=1,3
      kcL(j)=p1L*pcL(j)/hc         
      kvL(j)=pvL(j)/hc              
      enddo
      endif

!      if (phixs)
!     & write(99,'(i8,a,i2,a,i2,a,i2,a,i2,a,4g14.7,a,2g14.7)') 
!     & icount,' iv=',iv,' ic=',ii,' ien=',ien,' ip=',ip,
!     & 'ec,ev,ener3,ener4=',ec,ev,energ3(ien),energ4(ien,iv,ii,ip),
!     & 'bq=',bq
*     -----------------------------------------------------------------
*     construct remaining vectors
*     -----------------------------------------------------------------
      do j=1,3
*     wavevector of cm of core and valence particles
      bkp(j)=kcL(j)+kvL(j)-(mc+mv)/mtot*ktot(j)
*     wavevector of relative motion of core and valence particles
      kp(j)=mc/(mc+mv)*kvL(j)-mv/(mc+mv)*kcL(j)
      enddo
      eKb=hc**2*dot(bkp,bkp)/(2.d0*mupt)
      eks=hc**2*dot(kp,kp)/(2.d0*mucv)
!     write(99,'(5i3,20g14.5)') iv,ii,ien,ip,iroot,ekb,eks,bkp(1:3)
      phbin(:)=1
*     -----------------------------------------------------------------
*     increment maximum relative energy 
*     -----------------------------------------------------------------
      if(eks.gt.erelmax) erelmax=eks

!      if (wrt)
!     & write(99,'(5i3,20g12.5)') iv,ic,ien,ip,iroot,ekb,eks,bkp(1:3)
      tmatsq=0; tmatsq2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ((eks.gt.emax).or.(eks.lt.emin)) goto 500  ! NEEDS CHECKING (AMM)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*     -----------------------------------------------------------------
*     CALCULATION OF T-MATRIX <bkp,kp|T|qp> STARTS HERE
*     -----------------------------------------------------------------
*     ybo = Kout_cm, xb =theta_cm (rad), xbd (degrees)
*     -----------------------------------------------------------------
*     compute polar angles (th,ph) of cm wave vector
*     -----------------------------------------------------------------
      ybo=m(bkp)
      xb=min(1d0,bkp(3)/ybo)
      xb=max(-1d0,xb)
      xb=acos(xb)
      xbd=xb/degrad
      yb=abs(ybo-Kthr)
      if(abs(bkp(1)).lt.1.d-6.and.abs(bkp(2)).lt.1.d-6) then
      phKb=0.d0
      else 
      phKb=atan2(bkp(2),bkp(1))
      endif
      if(wrt) 
     & write(99,'(4i3," phikb,xbd=",20g12.5)') iv,ii,ien,ip,phKb,xbd  
*     -----------------------------------------------------------------
*     pick out nearest angles indices in array for interpolations
*     -----------------------------------------------------------------
      xyt(:,:)=0
!      iang=nint((xbd-thmin)/dth)+1
!      if (wrt) write(99,*)'thmin,thinc=',thmin,dth
      iang=int((xbd-thmin)/dth)+1 ! CHECK
      nial=max0(iang-2,1);  
!      if (wrt)write(99,*) 'iang,nial=',iang,nial
      niag=min0(nial+5,nth)
      nial=niag-5
      do iang=1,6
      xyt(1,iang)=angsr(nial+iang-1)
c v26
      xtab(iang)=angsr(nial+iang-1)
      enddo
      if (wrt)write(99,*) 'iang,nial,niag,angsr=',iang,nial,niag
      if (wrt)write(99,*)'angrs=',angsr(1),angsr(nth),angsr(nial+iang-1)
      
*     -----------------------------------------------------------------
*     compute polar angles (theta,phi) of relative wave vector
*     -----------------------------------------------------------------
      mkp=m(kp)
      co=kp(3)/mkp
      if(abs(kp(1)).lt.1.d-6.and.abs(kp(2)).lt.1.d-6) then
      phks=0.d0
      else 
      phks=atan2(kp(2),kp(1))
      endif
      co2=co*co
      si2=1.d0-co2
      si=sqrt(abs(si2))
      if (wrt)write(99,'(4i3," phks=",20g12.5)') iv,ii,ien,ip,phKs
*     --------------------------------------------------------------- 
*     calculate the spherical harmonics for this relative wave vector
*     ---------------------------------------------------------------
      phask1=(1.d0,0.d0)
      phask=exp(xi*phks)
      do ilval=0,lmax
      phask1=phask1/phask
      phas=phask1
      do inu=1,2*ilval+1
      phas=phas*phask
c AMM: Use explicit forms for PLM (l<9); otherwise, original recursive formula
       if (ilval.le.10) then
       aleg=plm_nr(ilval,inu-ilval-1,co,co2,si,si2)
       else
       aleg=plmrec(ilval,inu-ilval-1,co,co2,si,si2)
       endif
!      ylm(ilval,inu)=ylmc2(ilval,inu-ilval-1)*phas*aleg
      ylm(ilval,inu)=cplm(ilval,inu)*phas*aleg
!     write(99,'(6i3,20g12.5)') iv,ii,ien,ip,ilval,inu,
!     & ylmc2(ilval,inu-ilval-1),cplm(ilval,inu)
      enddo
      enddo
      fxyc(:,:)=0; fxytab=0
      tmat(:,:,:)=0; fin=0;
      kinset(:)=.false.
      wkfac(:)=0
*     -------------------------------------------- 
*     interpolate the sum \sum_i gsolap(i)*tmat(i)
*     --------------------------------------------
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      jpf=jpiset(iset)%jtot
      nex=jpiset(iset)%nex
      nmf=nint(2.*jpf+1.)
      iam=nmi*nmf 
      if (iam.lt.0) stop 'iam<0! '
      do inc=1,nchan
      ic=jpiset(iset)%cindex(inc)
      if (ic.ne.icore) cycle

      IF ((jpiset(iset)%nho).gt.0) then ! PS's -> need overlaps
        if (eks.gt.emax) then
        fv(1:nex)=0
        else
        do n=1,nex
        haux(1:ncont)=gsolap(iset,n,inc,1:ncont)
        fv(n)=ffc4((eks-emin)/dec,haux,ncont)
        enddo
        endif
      ELSE                              ! BINs
        if (inc.ne.jpiset(iset)%inc) cycle  ! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ik=0
! determine and store bin w(k) factor for this j/pi set
      do n=1,nex
      ki=binset(iset)%klow(n)
      kf=binset(iset)%kup(n)
      dk=kf-ki
      nk=binset(iset)%nk(n)
      dkb=(kf-ki)/(nk - 1)

      if ((mkp.ge.ki).and.(mkp.le.kf)) then
        ik=nint((mkp-ki)/dkb+1)     
        if (ik.gt.nk) stop 'ik too large'
        if (ik.le.0) stop 'ik <=0!!'
        kinset(iset)=.true.
!        phbin(iset)=exp(xi*binset(iset)%pshel(n,ik))
        wkfac(iset)=binset(iset)%wk(n,ik)
       if (wrt)write(99,'(a,1f6.2,i3,4f12.6)')'eks,iset,phbin=',
     &   eks,iset,phbin(iset)
       if (wrt)write(99,*)'iset,ik,wkfac=',mkp,iset,ik,wkfac(iset)
      endif !mkp within this bin
      enddo !n 

      endif  ! PS or Bins?

c Interpolate m-dependent amplitudes
      do iii=1,iam

c *** PS's .......................................................
      if ((jpiset(iset)%nho).gt.0) then 
      ier=nint((eks-emin)/dec)+1
      niel=max0(ier-2,1)
      nieg=min0(niel+5,ncont)
      niel=nieg-5
      if (wrt) write(99,*)'eks=',eks, 'ier=',ier,'niel=',niel

      do ier=1,6
      iecv=niel+ier-1 
      ecv=emin+(iecv-1)*dec 
      xyt(2,ier)=emin+(iecv-1)*dec 
c v26
      ytab(ier)=emin+(iecv-1)*dec
!      if (wrt) write(99,*)'ier,ecv=',ier,ecv
       ecmf=ecmi-abs(ebind)-ecv-excore
!      if (excore.gt.0) stop' excore!'
!      if (ecmf.lt.0) then
!          write(*,*) 'Ecm=',Ecmf,'Ecv=',Ecv; stop
!      endif
!      if (ecmf.lt.1e-4) ecmf=1e-4
       Kcmf_cont=sqrt(2.d0*mupt*ecmf)/hc     
c1      f2t=-2.*pi*hc**2/mupt*sqrt(Kcmi/Kcmf) ! conversion factor f(theta) -> T-mat
c1 ------------------------------------
      do ith=0,5
      sumps=0.d0
      do n=1,nex
      if (energ(iset,n).lt.0)    cycle ! bound states
      if (energ(iset,n).gt.Ethr) cycle ! closed channels
      iex=ips(iset,n)
! famps0() already transformed to T-matrix
!      sumps=sumps+fv(n)*f2t*famps0(iex,nial+ith,iii)

!       ecmf=ecmi-abs(ebind)-energ(iset,n)
!       Kcmf=sqrt(2.d0*mupt*ecmf)/hc  
!       raux=sqrt(kcmf/kcmf_cont)
!       sumps=sumps+fv(n)*raux*famps0(iex,nial+ith,iii)
      sumps=sumps+fv(n)*famps0(iex,nial+ith,iii)


!!!!

      enddo !  n (PS's)
      fxyc(ith+1,ier)=sumps
c v26
      fxytab(ith+1,ier)=sumps
!      if (wrt) write(99,*)'ith,ier,fxy=',ith,ier,sumps
      enddo !  ith (theta)
      enddo ! ier
c v26
!      tmat(iset,inc,iii)=f2c(xb,eks,xyt,fxytab,6,6,nord,6,6)
       tmat(iset,inc,iii)=f2c(xb,eks,xyt,fxyc,6,6,nord,10,maxne)
!      tmat(iset,inc,iii)=fint2dd(xtab,ytab,fxytab,xb,eks,6,6,6)
!      tmat(iset,inc,iii)=fint2d(xtab,ytab,fxytab,xb,eks,6,6,nord,6)

!             ijcmsec=fint2db2(ecmb,thcm,xsecmb,iecm,jthcm,necm+1,
!     +       nthcm+1,1.0d0)
!       tmat(iset,inc,iii)=fint2db_jin(xtab,ytab,fxytab,xb,eks,6,6,1.0d0)
     

      if (wrt) then
          write(99,*)'iii=',iii,'xb,yb=',xb,eks
          write(99,'("xyt(1,:)=",20g12.5)') (xyt(1,ith),ith=1,6)
          write(99,'("xyt(2,:)=",20g12.5)') (xyt(2,ier),ier=1,6)
!          write(99,'("fxyc=",20g12.5)') (fxyc(1:,ier),ier=1,6)
!          write(99,*)' -> tmat=',tmat(iset,inc,iii)
      endif

      ELSE 
c *** BINs  .......................................................................
      ibin=0
      if (kinset(iset)) then 
!      if (abs(wkfac(iset)).lt.1e-6) 
!     &    write(99,*)'iset,wkfac=',iset,wkfac(iset)
!       if (eks < binset(iset)%emin-5) cycle
!      if (eks > binset(iset)%emax+5) cycle
!      if (mkp < binset(iset)%klow(1))   cycle
!      if (mkp > binset(iset)%kup(nex)) cycle

      do n=1,nex
!      ki=binset(iset)%klow(n)
!      kf=binset(iset)%kup(n)
!      dk=kf-ki
!      nk=binset(iset)%nk(n)
      ecmf=ecmi-abs(ebind)-energ(iset,n)
      if(wrt) write(99,'(a,i3,6f12.6)')'n,k,bindk=',
     & n,binset(iset)%khat(n),wkfac(iset) !phbin(iset)/sqrt(dk)
      if (ecmf.lt.0) then
          write(*,*) 'Ecm=',Ecmf,'Ecv=',Ecv; stop
      endif
      if (ecmf.lt.1e-4) ecmf=1e-4
      Kcmf=sqrt(2.d0*mupt*ecmf)/hc     
c1      f2t=-2.*pi*hc**2/mupt*sqrt(Kcmi/Kcmf) ! conversion factor f(theta) -> T-mat
      xyt(2,n)=abs(Kcmf-Kthr)
!      call flush(99)
      do ith=0,5
        iex=ips(iset,n)
        fxyc(ith+1,n)=famps0(iex,nial+ith,iii)
     &               *sqrt(pi/2.)*wkfac(iset) ! includes phase e^(-i delta_el)
!     &             *sqrt(pi/2.)/sqrt(dk)  * phbin(iset)
      enddo !ith

      if(wrt)write(99,'(a,4i4,14g12.6)') 'iset,n,nial,ik,',
     &    iset,n,nial,ik,fxyc(1:6,n)

      if(wrt)write(99,'(a,i3,20e12.4)') '3 famp,iii=',iii,
     &  (famps0(iex,nial+ith,iii),ith=0,5)


      enddo !n 
      if2c=if2c+1
      tmat(iset,inc,iii)=f2c(xb,yb,xyt,fxyc,6,nex,nord,10,maxne) ! CHECK DIMS HERE!!!!!!!!!!!!!!!!
c1      fin(iset,inc,iii)=f2c(xb,yb,xyt,fxyc2,6,nex,nord,10,50) ! TESTING
      if (wrt) then
          write(99,*)'iii=',iii,'xb,yb=',xb,yb
          write(99,'("xyt(1,:)=",20g12.5)') (xyt(1,ith),ith=1,6)
          write(99,'("xyt(2,:)=",20g12.5)') (xyt(2,ier),ier=1,nex)
          write(99,'("fxyc(1,:)=",12g12.5)') fxyc(1,1:nex)
          write(99,'("fxyc(2,:)=",12g12.5)') fxyc(2,1:nex)
          write(99,*)
          endif
!          if(abs(tmat(iset,inc,iii)).gt.1e5) then
!           write(*,*)'Big tmat for iset,inc,iii=',iset,inc,iii,
!     &     tmat(iset,inc,iii)
!          endif
      else 
         tmat(iset,inc,iii)=0
      endif ! k in this set?
         
      endif ! Bins or PS's     
c    ...............................................................................
      enddo ! iii
      enddo ! inc channel
      enddo ! jpset
*     ----------------------------------------------------------------- 
*     |t-matrix|**2 summed on spin projections / (2sp+1)
*     ----------------------------------------------------------------- 
      tmatsq=0.d0
      tmatsq2=0. !!TESTING
      do imo=1,nmi
      mo=imo-jpiset(1)%jtot-1
      ind=0
      do imu=1,nint(2*sc+1)
      mu=imu-sc-1
      do isig=1,nint(2*sn+1)
      sig=isig-sn-1
      ampt=(0.d0,0.d0)
      ampt2=0. !!!! TESTING
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      nex  =jpiset(iset)%nex
      jpf  =jpiset(iset)%jtot
      do inc=1,nchan  ! sum over open channels
      ic  =jpiset(iset)%cindex(inc)
      if (ic.ne.icore) cycle
      li     =jpiset(iset)%lsp(inc)
      ji     =jpiset(iset)%jsp(inc)
      jci    =jpiset(iset)%jc(inc)
      rli=li
      do inu=1,2*li+1
      rnu=inu-li-1
      if (fail3(rli,sn,ji)) cycle
      if (fail3(ji,jci,jpf)) cycle
      if (abs(rnu+sig).gt.ji.or.abs(rnu+sig+mu).gt.jpf) cycle
      ind=ind+1
!      if (abs(cgc(ind).lt.eps)) cycle
      iam=nint((imo-1)*(2*jpf+1)+rnu+sig+mu+jpf+1)
      if (iam.lt.0) then
          write(*,*)'solvecc error: iam<0!';
          write(*,*)'iam,imo,jpf,rnu,sig,mu=',iam,imo,jpf,rnu,sig,mu
          write(*,*)'shape(tmat)=',shape(tmat)
          stop
      endif

!!! CHECK
!      raux=cleb(rli,rnu,sn,sig,ji,rnu+sig)*
!     .cleb(ji,rnu+sig,jci,mu,jpf,rnu+sig+mu)
!      if (abs(raux-cgc(ind)).gt.1e-3) then
!      write(99,*)'cgc mismatch!: rli,rnu,sig,mu,jpf=',rli,rnu,sig,mu,jpf     
!      endif
!AMMoro Sept 16
!      ampt=ampt+(-xi)**li*cgc(ind)*ylm(li,inu)*  ! NOTE phbin() inserted here!!!!!
      ampt=ampt+(-xi)**li*cgc(ind)*ylm(li,inu)* ! phbin(iset)*  ! NOTE phbin() inserted here!!!!!
     .tmat(iset,inc,iam)*exp(xi*nint(mo-rnu-sig-mu)*phKb)

!TESTING
c1      ampt2=ampt2+(-xi)**li*cgc(ind)*ylm(li,inu)*phbin(iset)*
c1     & fin(iset,inc,iam)*exp(xi*nint(mo-rnu-sig-mu)*phKb)
!!
      enddo ! inu
      enddo ! inc
      enddo ! j/pi
      tmatsq=tmatsq+abs(ampt)**2
c2      tmatsq2=tmatsq2+abs(ampt2)**2   !!! TESTING
      enddo ! isigma
      enddo ! imu
      enddo ! imo
c2      tmatsq2=tmatsq2*((2.d0*pi)**3)/dble(nmi)
      tmatsq=(4.d0*pi/mkp)**2*tmatsq/dble(nmi)
c1      if (icount.lt.100)
c1     & write(*,*)'tmatsq=',tmatsq,tmatsq2,tmatsq2/tmatsq
      

500   continue
*     -----------------------------------------------------------------
*     Phase Space Factor
*     -----------------------------------------------------------------
      mult=(mupt*mv/(2.d0*pi*hc**2)**2)*m(kvL)/m(qp)
      mult=mult*mc*m(kcL)*hc/(2.d0*pi*hc)**3
*     -----------------------------------------------------------------
*     vector (k[detected]-ktot) for phase space factor
*     -----------------------------------------------------------------
      do j=1,3
      if(idet.eq.1) then
      tem(j)=kcL(j)-ktot(j)
      else
      tem(j)=kvL(j)-ktot(j)  
      endif
      enddo
*     -----------------------------------------------------------------
*     detected particle dependent parts
*     -----------------------------------------------------------------
      if(idet.eq.1) then
      mult=mult*abs((mt/((mv+mt)+mv*dot(tem,kvL)/dot(kvL,kvL))))
      else
      mult=mult*abs((mt/((mc+mt)+mc*dot(tem,kcL)/dot(kcL,kcL))))       
      endif
*     -----------------------------------------------------------------
*     triple differential cross section (in mb)
*     -----------------------------------------------------------------
      sigphi(ip)=sigphi(ip)+10.d0*mult*tmatsq
52    if(iroots.eq.0) sigphi(ip)=0.d0
*     -----------------------------------------------------------------
*     prepare to go around again at this phi if two roots required
*     -----------------------------------------------------------------
      if(iroots.eq.2.and.iroot.eq.1) then
      iroot=2
      goto 53
      endif
*     -----------------------------------------------------------------
*     close the phi loop
*     -----------------------------------------------------------------
40    continue
*     -----------------------------------------------------------------
*     integrate over phi angles NOW if more than one phi angle
*     -----------------------------------------------------------------
      if(itphi.eq.1) then
      xsig=sigphi(1)
      else if(iflag.eq.0) then
      call sim(sigphi,xsig,1,itphim,dphir,itphim)
      xsig=2.d0*xsig
      else
      xsig=sigphi(1)*(phiu-phil)*degrad
      endif
      if(ien.eq.iten) then
!      print 315,ien,iv,ii,xsig,En,iflag
      endif
*     -----------------------------------------------------------------
       xs3body(ien,iv,ii)=xsig

! AMM July 2019
       if (phixs) then
       do ip=1,itphim
       xs3body_phi(ien,iv,ii,ip)=sigphi(ip)
       enddo !ip
       endif
       
!      write(77,315) ien,iv,ii,xsig,En
!      STOP
!      call flush(77)
*     -----------------------------------------------------------------
*     close the angle (core and valence thetas) and energy loops
*     -----------------------------------------------------------------
!      energ4(ien)=ev
!      write(0,*)ec,energ3(ien),ev
30    continue      
20    continue
10    continue
!$OMP END PARALLEL DO
      print 510
      print*,'  erelmax =',erelmax
      write(*,*)'icount=',icount,' if2c=',if2c
      write(*,*)'Dims=',itc,itv,iten,itphim
      print 510
!MGR--------------------------------------------------------------------
      icount=0
      do ii=1,itc
        do iv=1,itv
          do ien=1,iten
            write(77,315) ien,iv,ii,xs3body(ien,iv,ii),energ3(ien)
            do ip=1,itphim
            icount=icount+1
            if (phixs) write(777,318) ien,iv,ii,ip,
     &       xs3body_phi(ien,iv,ii,ip),energ3(ien),energ4(ien,iv,ii,ip)
            enddo !ip
          enddo
        enddo
      enddo    
      deallocate(xs3body,energ3)
      
      if (phixs) deallocate(xs3body_phi,energ4)
      write(*,*)'xs3body_phi has',icount, 'elements'
!-----------------------------------------------------------------------      
      call cpu_time(tf)
      t3d=tf-ti
      return 
307   format(i5,f15.8,3f10.3)
315   format(3i5,d19.8,1x,f12.6,i5)
316   format(3f10.3)
318   format(4i5,d19.8,1x,2f12.6,i5)
510   format('  ------------------------------------------------------',
     +'-----')
900   write(*,*)'Error reading f-amps file!'; stop
      end subroutine
*----------------------------------------------------------------------
      real*8 function dot(a,b)
*     scalar product of two vectors
      real*8 a(3),b(3)
      dot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end
*----------------------------------------------------------------------
      real*8 function m(a)
*     modulus of vector
      real*8 a(3)
      m=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      return
      end
*----------------------------------------------------------------------
      recursive real*8 function plmrec(il,imm,c,c2,s,s2) result(res)
      implicit real*8(a-h,o-z)
*-----------------------------------------------------------------------
*     Associated Legendre functions
*     ECS 24th Nov 2008
*     l<=2 calculated explicitly
*     l>2  calculated recursively
*-----------------------------------------------------------------------
      res=1.d0
      if(il.eq.0) return    !il=0
      im=abs(imm)
*-----------------------------------------------------------------------
      if(il.eq.1) then      !il=1
      if(im.eq.0) then
      res=c
      else if(im.eq.1)then
      res=s
      endif
      else if(il.eq.2) then !il=2
      if(im.eq.0) then
      res=(3.d0*c2-1.d0)/2.d0
      else if(im.eq.1) then
      res=3.d0*s*c
      else if(im.eq.2) then
      res=3.d0*s2
      endif
      elseif(im.gt.1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ((1.-c2).lt.1e-3) c2=0.99 ! AMoro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      res=(1.d0/sqrt(1.d0-c2))*
     &(2*(im-1)*c*plmrec(il,im-1,c,c2,s,s2))-
     &((il+im-1)*(il-im+2))*plmrec(il,im-2,c,c2,s,s2)
      elseif(il.gt.2) then
      res=(1.d0/real(il-im))*((2*il-1)*c*plmrec(il-1,im,c,c2,s,s2)
     &-(il+im-1)*plmrec(il-2,im,c,c2,s,s2))
*-----------------------------------------------------------------------
      endif
      return
      end
*----------------------------------------------------------------------
      complex*16 function f2c(xb,yb,xyt,fxyt,nnx,nny,nord,mmx,mmy)
      implicit real*8 (a-h,o-z)
      complex*16 fxyt(mmx,mmy),fxy(33,33),cm1,cj3
      real*8 x(33),y(33),xyb(2),xyt(2,mmy)
      integer nxy(2),nbg(2)
      nnn=nord+1
      nxy(1)=nnx
      nxy(2)=nny
      xyb(1)=xb
      xyb(2)=yb
*     ---------------------------------------------------------------
*     loop on x and y grid finding sub-matrix needed
*     ---------------------------------------------------------------
      do i=1,2
      if(xyb(i).lt.xyt(i,1)) then
      num=1
      else if(xyb(i).gt.xyt(i,nxy(i))) then
      num=nxy(i)-nnn+1
      else
      min=1
      max=nxy(i)
      num=nxy(i)/2
      if(num .lt.1) num=1
55    if(max-min.lt.4) goto 70
      if(xyt(i,num)-xyb(i)) 60,82,61
60    min=num
      num=num+(max-min+1)/2
      goto 55
61    max=num
      num=num-(max-min+1)/2
      goto 55
70    num=max
71    if(xyt(i,num)-xyb(i)) 82,82,81
81    num=num-1
      goto 71
*     ---------------------------------------------------------------
*     change/correction July 2008
*     ---------------------------------------------------------------
82    num=max0(1,num-nord/2)
      num=min0(num,nxy(i)-nord)
      endif
      nbg(i)=num-1
      enddo
*     ---------------------------------------------------------------
      nx=nbg(1)
      ny=nbg(2)
      do i=1,nnn
      x(i)=xyt(1,nx+i)-xb
      y(i)=xyt(2,ny+i)-yb
      do j=1,nnn
      fxy(i,j)=fxyt(nx+i,ny+j)
      enddo
      enddo
*     ---------------------------------------------------------------
      do ii=2,nnn
      iii=ii-1
      do jj=ii,nnn
      cj2=x(jj)-x(iii)
      cj3=fxy(iii,iii)*x(jj)-fxy(jj,iii)*x(iii)
      do mm=ii,nnn
      cm1=fxy(iii, mm)*x(jj)-fxy(jj, mm)*x(iii)
      fxy(jj,mm)=(cj3*y(mm)-cm1*y(iii))/(y(mm)-y(iii))/cj2 ! ****
      enddo
      enddo
      enddo
      f2c=fxy(nnn,nnn)
      return
      end




c *** Calculates approximate energy distribution from PS's
      subroutine budist(enex,nset,sigset)
      implicit real*8(a-h,o-z)
      parameter (nmax=500)
      integer nb
      real*8  enex(1:nset),sigset(1:nset)
      real*8 d1,di,dnex,be
      be=0.0
c     -----------------------------------------------
      nb=0
      do ia=1,nset
!      write(*,*)'ia=',ia
      if(enex(ia).lt.0) nb=nb+1 ! bound states
      enddo

!      print*,'budist: nb=',nb

      if (nset-nb.lt.3) return
      open(15,file='dsde_ps.xs') 
!      d1=2.d0*ENEX(2) - 1.5d0*ENEX(1) - 0.5d0*ENEX(3)
      d1=enex(nb+1)
!      print*,'d1=',d1
      write(15,50) ENEX(nb+1)-be,sigset(nb+1)/d1,enex(1)+d1/2.d0 
      write(18,50) ENEX(nb+1)-be,sigset(nb+1),
     X     enex(nb+1)+d1/2.d0 

c     Interpolation for 1<=IA<NEX-1
      do IA=nb+2,nset-1
         if (enex(ia).lt.enex(ia-1)) write(15,*)'&'
         di=(ENEX(IA+1)-ENEX(IA-1))/2.d0
         write(15,50) ENEX(IA)-be,sigset(IA)/di,
     X           enex(ia)+di/2.d0,di
         write(18,50) ENEX(IA)-be,sigset(IA),
     X           enex(ia)+di/2.d0
      enddo
c     Interpolation for IA=NSET
      dnex=-2.d0*ENEX(NSET-1)+1.5d0*ENEX(NSET)+0.5d0*ENEX(NSET-2)
      write(15,50) ENEX(NSET)-be,sigset(NSET)/dnex,
     X        enex(nset)+dnex/2.d0,dnex
      write(18,50) ENEX(NSET)-be,sigset(NSET),
     X        enex(nset)+dnex/2.d0

 50   format(6f12.4)
      sig=0
      do i=1,nset
         sig=sig+sigset(i)
      enddo
      write(*,*) ' => Total inel+bu x-section=',sig
      write(*,*)  '------------------------------------'
      write(15,*) '&'
      write(18,*) '&'
      close(15)
      return
      end


!MGR amplitudes with target excitation or spin of target

c *** --------------------------------------------------------------
c **  Scattering amplitudes & cross sections
c *** --------------------------------------------------------------
      subroutine xsecs_tdef(kin,ncc,iexgs)
      use xcdcc,    only: smats,elab,jpch,jtmin,jtmax,nex,exch,parch,
     &                    famps,jbord,jump,rel
      use channels, only: jptset,jpiset,jpsets,tset
      use factorials
      use sistema
      use constants
      use parameters, only:  maxchan,maxsets
      use globals,  only: debug,written,verb
      use nmrv,     only: hort
      use wfs,      only: nr,dr,energ,rvec,wfr
      implicit none
      logical :: doublexs,triplexs,phixs,jsets(maxsets)
c     ---------------------------------------------------------
      integer, parameter :: klog=99
      integer   , parameter:: kamp=136, kxs=300,kfam=137,
     &                        ksmat=170, !ksj=157,
     &                        nearfa=1
      integer :: istat,kin
      integer :: inc,iexgs,icc,ncc,nch,ni,nf,lmax,iex,lfmax,ith,icore
      integer :: partot,nth,njt,parold,maxleg,ner,itex,itarg
      integer :: iam,im,imp,mlp,li,lf,nmi,nmti,nmtf,nmf,mi,mf,ijt,l
      integer :: nmtfmax,nmfmax,imt,imtp,idt
c     ------------------------------------------------------------
      real*8  :: kron,cleb,ylmc,ylmc2 !external functions
      real*8    , parameter:: zero=0
      real*8    , parameter:: theps=1e-3,alpha=0d0
      real*8    , allocatable:: lfv(:),xv(:)
      real*8  :: jtot,jpgs,jpi,jpf,jpold,jtarg,jti,jtf,jlpi,jlpf
      real*8  :: kcmi,kcmf,exc,vi,vf,tkch,etarg,partarg,extarg
      real*8  :: ecmi,ecmf,etai,etaf,mupt,ermin,ermax,jwtarg
      real*8  :: thmin,thmax,dth,cth,sth,s2,th,thrad
      real*8, allocatable :: pl(:,:),delci(:),delcf(:)
      real*8  :: rm,rmp,rmlp,lri,lrf,rmt,rmtp
      real*8  :: factor, r1,r2,delif,yffc
      real*8  :: xsruth,xs,xstot,xsj,xr,xrj,sigex(nex),xsaux
      real*8  :: sigtex(ntex+1),sigextex(nex,ntex+1)
c     -----------------------------------------------------------
      complex*16, parameter:: zz=cmplx(0.,1.)
      complex*16,allocatable:: ampl(:,:,:,:,:),fam(:,:,:,:,:)
      complex*16, allocatable:: ampaux(:,:,:,:,:),gaux(:)
      complex*16 smat,c1,faux,fc,phci,phcf,caux,cfival,resc,ffc4                     
c     -----------------------------------------------------------
      CHARACTER  PARITY(3)
      character*40 fileamp
      DATA PARITY / '-','?','+' /
      
      complex*16,allocatable:: ampauxJ(:,:,:,:,:,:) !MGR
      integer :: njtot,nlp,ijtot,ilp,njtotcalc,ijtotcalc,nspline,ispline
      integer :: jspl
      integer,allocatable :: interpJ(:)
      real*8 ljmax      
c ---------------------------------------------------------------
c PLM
!      REAL*8 PM(0:1000,0:1000),PD(0:1000,0:1000)

      namelist/xsections/ thmin,thmax,dth,fileamp,doublexs,jsets,
     &                    ermin,ermax,ner,icore,itarg,
     &                    triplexs,rel

c initialize -------------------------------------------------
      written(kamp)=.false.
      written(kxs)=.true.
      written(kfam)=.true.
      
!      written(ksj)=.false.
      written(kxs+1:kxs+min(nex,9))=.true.
      pi=acos(-1d0)
      thmin=0; thmax=0; dth=0; 
      doublexs=.false. ; triplexs=.false. ; phixs=.false.
      rel=.false.
      ermin=-1; ermax=-1
      jsets(:)=.true.        
      fileamp=""
      ner=0
      
c     -----------------------------------------------------------
      rewind(kin)
      itarg=1
      read(kin,nml=xsections)
      write(0,*)'reading xsections namelist for target def'

      if (thmax.gt.thmin) then
        nth=nint((thmax-thmin)/dth)+1
      else
        write(*,*)'thmax<thmin!; Aborting'; stop
      endif

      if (ner.eq.0) then
         write(*,*)'NER=0!!';
      endif

      if (ncc.eq.0) goto 1000

      jpgs =jpch(iexgs)
      mupt =mp*mt/(mp+mt)*amu 
      ecmi =elab*mt/(mp+mt)
      kcmi =sqrt(2*mupt*ecmi)/hc
      etai =zp*zt*mupt/kcmi/hc/finec 
      lmax =nint(jtmax+maxval(jpch)+maxval(tset(:)%jtarg))



      write(*,'(/,5x,"*** SCATTERING AMPLITUDES AND X-SECTIONS ***",/)')
      write(*,320) elab, ecmi,kcmi, etai
320   format(3x,"Elab=",1f7.2, 5x,'Ecm=',1f8.3,5x,
     &  'Kcm=',1f6.3, 4x,'Eta param. =',1f7.4)
      write(*,'(/,3x, "o Jtmin=",1f4.1,3x,"Jtmax=",1f5.1,/ )')
     &  jtmin,jtmax
      call flush(6)

c *** Factorials (CHECK IF THIS IS NEEDED!!!!!!!!!!!!!!)
!      lmax=nint(jtmax+maxval(jpch))     
!      write(*,*)'famps: lmax=',lmax
!      call factorialgen(2*lmax)


c *** Compute Coulomb phases for l< lmax
       allocate(delci(0:lmax),stat=istat)
       if (istat>0) then
         write(*,*) 'readsmat: allocating memory for DELC failed!';
         stop
       endif
       call coulph (etai,delci,lmax) !Coulomb phase-shifts for inc channel



c *** Compute & store amplitudes (Fresco, CPC Eq. (3.30))   
c *** (zero target spin assumed here!) !MGR let's try and change it
      nmi    =nint(2*jpgs+1)
      nmti   =nint(2*jtgs+1)
      nmfmax =nint(2*maxval(jpiset(:)%jtot)+1)
      nmtfmax=nint(2*maxval(tset(:)%jtarg)+1)
      allocate(famps(nex,ntex+1,nth,nmi*nmfmax*nmti*nmtfmax))

      
      if (nex.lt.1) then
         write(*,*)'solvecc: nex<1!'; stop
      endif
      if (iexgs.lt.1) then
         write(*,*)'solvecc: inc=',iexgs; stop
      endif

      do iex=1,nex  ! loop in projectile states
      jpf  = jpch(iex) 
      exc  = exch(iex)
      do itex=1,ntex+1
      jtarg =tset(itex)%jtarg
      etarg =tset(itex)%extarg
      tkch=ecmi+exch(iexgs)-exc-etarg
      
      
      
      if (tkch.lt.0) cycle    !closed channel
      kcmf = sqrt(2*mupt*tkch)/hc
      etaf = zp*zt*mupt/kcmf/hc/finec 
      nmf  = nint(2*jpf +1)
      nmtf = nint(2*jtarg+1)
      lfmax= nint(jtmax+jpf+jtarg)
      maxleg=max(1,lfmax)
      write(*,*) 'Projectile state',iex,'Target state',itex,'Ecm',tkch,
     & 'k_f:',kcmf
!      write(0,*) 'lfmax',lfmax,'nmi',nmi,'nmti',nmti,'nmf',nmf,'nmtf',
!     & nmtf
      allocate(ampl(0:lfmax,nmi,nmti,nmf,nmtf)) ! angle-independent amplitudes
      allocate(lfv(0:lfmax))
      allocate(delcf(0:lfmax),stat=istat)
      call coulph(etaf,delcf,lfmax) !Coulomb phase-shifts for final channel

!MGR let's find the maximum number of L' that there are for any Jtot
      nlp=0
      njtot=0
      ljmax=0d0
      do icc=1,ncc
        jtot  =jptset(icc)%jtot
        njtot =max(njtot,floor(jtot))
        nch   =jptset(icc)%nchan
        do nf=1,nch
          if (jptset(icc)%nchan.eq.0) cycle
          if (jptset(icc)%idx(nf).ne.iex) cycle 
          if (jptset(icc)%idt(nf).ne.itex) cycle
          lf  =jptset(icc)%l(nf);
          ljmax=max(ljmax,abs(jtot-lf-0d0))           
        enddo  
      enddo
      njtot=njtot-floor(jbord(1))
      nlp=nint(2d0*ljmax+1d0)
!      write(0,*) 'nlp',nlp,'njtot',njtot
!-----------------------------------------------------------------------


!      write(kfam,301) iex,exc,jpf,parity(parch(iex)+2)
!      write(kfam,*) elab,nth
      write(kfam,'(4f6.1,i5,i2,1f10.3)') jpgs,jtgs,jpf,
     & jtarg,nth,nearfa,elab
      write(kxs,302) iex,exc,jpf,parity(parch(iex)+2)
      if (iex+itex.le.10) write(kxs+(iex-1)*(ntex+1)+itex,302) 
     &      iex,exc,jpf,parity(parch(iex)+2),itex,tset(itex)%extarg
300   format(3x,'o Angle-indep amplitudes for state #',i3,
     &   ' Ex=',1f7.3,' MeV',' J/pi=',1f5.1,a1)     
301   format(1x,'# f(theta) amplitudes for state #',i3,
     &   ' Ex=',1f7.3,' MeV',' J/pi=',1f5.1,a1) 
302   format(1x,'# Cross sections for state #',i3,
     &   ' Ex=',1f7.3,' MeV',' J/pi=',1f5.1,a1,'Targ. state',i3,
     &   ' E_targ=',1f7.3)


      allocate(ampauxJ(0:njtot,nlp,nmi,nmti,nmf,nmtf),interpJ(0:njtot)) !MGR
      ampauxJ(:,:,:,:,:,:)=0d0;interpJ(:)=1 !MGR
      ijtot=-1 !MGR
      ampl(:,:,:,:,:)=0; 
      lfv(:)=0

      do icc=1,ncc
       partot=jptset(icc)%partot
       jtot  =jptset(icc)%jtot
       nch   =jptset(icc)%nchan
       if (jptset(icc)%interp) then
!          write(*,*)'skipping jtot=',jtot
          cycle 
       endif
       if (jtot.ge.jbord(1)) then
       ijtot=floor(jtot-jbord(1))!MGR
       interpJ(ijtot)=0!MGR
       endif

       do ni=1,nch
       if (jptset(icc)%idx(ni).ne.iexgs) cycle ! no inc waves 
       if (jptset(icc)%idt(ni).ne.1) cycle !no inc waves NEEDS REPHRASING MGR
       li   =jptset(icc)%l(ni)
       lri  =li
       jpi  =jptset(icc)%jp(ni)             ! should be just jpgs
       jlpi =jptset(icc)%jlp(ni)
       jti  =jptset(icc)%jt(ni)
       nmi  =nint(2*jpi+1)
       nmti =nint(2*jti+1)

       do nf=1,nch
       if (jptset(icc)%nchan.eq.0) cycle
       if (jptset(icc)%idx(nf).ne.iex) cycle 
       if (jptset(icc)%idt(nf).ne.itex) cycle 
       lf  =jptset(icc)%l(nf);   
!       nlf=nlf+1;     
!       lfv(icc,nlf)=lf 
       lrf =lf
       jpf =jptset(icc)%jp(nf)
       jlpf=jptset(icc)%jlp(nf)
       jtf= jptset(icc)%jt(nf)
       nmf =nint(2*jpf+1)
       nmtf=nint(2*jtf+1)
       smat=smats(icc,ni,nf) 
       
       if (lf.gt.lfmax) stop 'internal error; lf >lfmax in famps!' 
325    format(' Chan:',2i5,2x,2f6.1,2x,"S=",2f11.7)

       do im=1,nmi
        rm= -jpi + (im -1)
        do imp=1,nmf
        rmp=-jpf + (imp-1)
        do imt=1,nmti
         rmt= -jti + (imt -1)
        do imtp=1,nmtf
         rmtp=-jtf + (imtp-1)
        mlp=nint(rm+rmt-rmp-rmtp)
        rmlp=rm+rmt-rmp-rmtp
        if (abs(mlp).gt.lf) cycle        
        r1=cleb(lri,zero,  jpi, rm, jlpi,rm)*
     &     cleb(jlpi,rm,jti,rmt,jtot,rm+rmt)
     &    *cleb(lrf,rmlp,jpf, rmp,jlpf,rmlp+rmp)*
     &     cleb(jlpf,rmlp+rmp,jtf,rmtp,jtot,rm+rmt)
        r2=sqrt(kcmf/kcmi)  ! this needs to be generalized for different partitions!!
        c1=sqrt(pi)/zz/kcmi
        phci =exp(zz*(delci(li)-delci(0)))
        phcf =exp(zz*(delcf(lf)-delcf(0)))
c Satchler (4.58) for spin zero target:
        ampl(lf,im,imt,imp,imtp)= ampl(lf,im,imt,imp,imtp)  
     &                 + phci*phcf*r1*r2*c1*sqrt(2*lri+1)  
     &                 * (smat-kron(ni,nf))*YLMC2(lf,mlp) 
!        ampaux(lf,im,imt,imp,imtp)= ampaux(lf,im,imt,imp,imtp)  
!     &                   + phci*exp(-zz*(delcf(lf)-delcf(0)))
!     &                   * r1*r2*c1*sqrt(2*lri+1)           
!     &                   * (smat-kron(ni,nf))*YLMC2(lf,mlp) 
        lfv(lf)=lf

!        write(777,*)'lf,mlp,ylm=', lf,mlp,YLMC(lf,mlp),YLMC2(lf,mlp)

        if (debug) then
        write(klog,'(6x,"Contrib to AMP: ni,li,iex=",
     &  3i3," => nf,lf,iex=",3i3,6x,"S=",2f12.6," r1,r2=",2f12.5)')
     &  ni,li,1,nf,lf,iex,smat,r1,r2
       endif
       
!MGR--------------------------------------------------------------------
       if(ijtot.ge.0) then
       ilp=nint(ljmax+jtot-lrf)+1 
!       write(0,*) 'ijtot',ijtot,'ilp',ilp,'im',im,'imp',imp, 'jtot',jtot
!     & ,'imt',imt,'imtp',imtp
       ampauxJ(ijtot,ilp,im,imt,imp,imtp)= 
     &  ampauxJ(ijtot,ilp,im,imt,imp,imtp)  
     &                 + phci*conjg(phcf)*r1*r2*c1*sqrt(2*lri+1)  
     &                 * (smat-kron(ni,nf))*YLMC2(lf,mlp)
       endif 
!-----------------------------------------------------------------------               
       enddo ! imtp
       enddo ! imt
       enddo ! imp
       enddo ! imp
       enddo ! nchf
       enddo ! nchi
       enddo ! icc (j/pi sets)
    

!!!!!!!!!!!!!!!!!!!!!!!!! INTERPOLATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (.true.) then
!First we need to know the points that are already calculated
      njtotcalc=0
      do ijtot=0,njtot
      if (interpJ(ijtot).eq.0) njtotcalc=njtotcalc+1
      enddo
!------------------------------------------------------------------------      
      allocate(gaux(njtotcalc),xv(njtotcalc))
      njtotcalc=0
      do ijtot=0,njtot
      if (interpJ(ijtot).eq.0) then
       njtotcalc=njtotcalc+1
       xv(njtotcalc)=ijtot
      endif
      enddo
!      write(0,*) 'Passed njtotcalc,', njtotcalc,'njtot',njtot
      do ijtot=0,njtot
       jtot  =ijtot+jbord(1)
       if (interpJ(ijtot).eq.0) cycle

!MGR Let's try with splines
!      pspline=(jtot-jbord(1))/jump(1)
!      call SPLINT(pspline,njtotcalc,ispline,nspline,spline)
      
!      write(0,*) 'Spline N',njtotcalc,' NP',nspline,'I',ispline,'P',
!     & pspline,'SP',spline(1:nspline)
       do im=1,nmi
        rm= -jpi + (im -1)
        do imp=1,nmf
        rmp=-jpf + (imp-1)
        do imt=1,nmti
         rmt= -jti + (imt -1)
        do imtp=1,nmtf
         rmtp=-jtf + (imtp-1)
        mlp=nint(rm+rmt-rmp-rmtp)
        rmlp=rm+rmt-rmp-rmtp
        do ilp=1,nlp

        lf=nint(ljmax+jtot-ilp)+1
        if(lf.gt.lfmax)write(0,*) 'lf',lf,'> lfmax',lfmax,'!!'
!       write(0,*)'jtot',jtot,'m',rm,'mp',rmp,'mt',rmt,'mtp',rmtp,'l',lf
        if (abs(nint(rm+rmt-rmp-rmtp)).gt.lf) cycle
             
        phcf =exp(2d0*zz*(delcf(lf)-delcf(0)))
        gaux(:)=0d0
!        write(0,*) 'Interpolating Jtot:', jtot,'L',lf 
        do ijtotcalc=1,njtotcalc
        gaux(ijtotcalc)=ampauxJ(nint(xv(ijtotcalc)),ilp,im,imt,imp,imtp)
        enddo
!MGR splines        
!        ampauxJ(ijtot,ilp,im,imp)=0d0
!        do jspl=1,nspline
!        if(interpJ(ijtot).eq.1) then
!        ampl(lf,im,imp)=ampl(lf,im,imp)+spline(jspl)*
!     &   gaux(jspl+ispline-1)*phcf
!        endif
!        ampauxJ(ijtot,ilp,im,imp)=ampauxJ(ijtot,ilp,im,imp)+spline(jspl)
!     &  *gaux(jspl+ispline-1)
!        enddo
                   
!        write(*,'(10(i3,2f10.5))') (nint(xv(l)),gaux(l),l=1,lfmax)
        caux=cfival(ijtot+0d0,xv,gaux,njtotcalc,0d0)
        ampl(lf,im,imt,imp,imtp)= ampl(lf,im,imt,imp,imtp)+caux*phcf
        ampauxJ(ijtot,ilp,im,imt,imp,imtp)=
     &  ampauxJ(ijtot,ilp,im,imt,imp,imtp)+caux
         
        enddo !ilp
        enddo !imtp
        enddo !imt
       enddo ! imp
       enddo ! im
       enddo ! ijtot (j/pi sets)
       deallocate(gaux,xv,ampauxJ,interpJ)
       endif 
!       write(0,*) 'Exit interpolation'

!       do ijtot=0,njtot
!       write(888,*) ijtot, dble(ampauxJ(ijtot,1,1,1)),
!     &  aimag(ampauxJ(ijtot,1,1,1))
!       enddo
!       write(888,*) '&'    
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c write angle-indep amplitudes ----------------------------------------
       if (written(kamp)) then
       write(kamp,300) iex,exc,jpf,parity(parch(iex)+2)
       do lf=0,lfmax
!225    WRITE(kamp,226) LF,((AMPL(lf,im,imp)*phcf,imp=1,nmf), im=1,nmi)
225    WRITE(kamp,226) lf,((((ampl(lf,im,imt,imp,imtp),imtp=1,nmtf),
     &  imp=1,nmf),imt=1,nmti), im=1,nmi)
226    FORMAT(' LP=',I5,1X,1P,10E12.4 / (10X,10E12.4),/)
       enddo ! lf
       endif


c Angle-dependent amplitudes: f(theta,m,mp) (Fresco CPC eq. (3.31)
       allocate(fam(nth,nmi,nmti,nmf,nmtf))
!       write(*,*)'allocate pl with lfmax+1=',lfmax+1
       allocate(pl(maxleg+1,maxleg+1))
       fam(:,:,:,:,:)=0; pl=0;
       do ith=1,nth
          th = thmin + dth*(ith-1)
          if (th<theps) th=theps
          thrad = th*pi/180.
          cth=cos(thrad)
          sth=sin(thrad)
          s2= sin(thrad/2)**2
          call PLM(cth,maxleg,maxleg,maxleg+1,PL)

c Point Coulomb (Satchler 4.11a), but with different phase convention 
c (no e^{2 i delc(0)}  factor as in Fresco !) 
          r1=-etai*log(s2)  !  2*delc(0)
          fc=-etai/(2*kcmi)*exp(zz*r1)/s2
c 
       iam=0
       do im=1,nmi
       rm= -jpi + (im -1)
       do imt=1,nmti
         rmt= -jti + (imt -1)
       do imp=1,nmf
       ! global index for spin projections
       rmp=-jpf + (imp-1)
       do imtp=1,nmtf
         rmtp=-jtf + (imtp-1)
       iam=iam+1 
       faux=0
       rmlp=rm+rmt-rmp-rmtp
       mlp=nint(rm+rmt-rmp-rmtp)
       if (abs(mlp).gt.lfmax) stop'mlp>lfmax'

       do lf=abs(mlp),lfmax
         if(mlp.lt.0) then
          r1=pl(lf+1,-mlp+1)
         else 
          r1=pl(lf+1,mlp+1)          
         endif 
        faux=faux + ampl(lf,im,imt,imp,imtp)*r1
       enddo !lf
       fam(ith,im,imt,imp,imtp)=fc*kron(im,imp)*kron(imt,imtp)*
     &  kron(iexgs,iex)*kron(itex,1) + faux !Needs rephrasing MGR
       famps(iex,itex,ith,iam)= fam(ith,im,imt,imp,imtp)
       enddo !imtp
       enddo !im
       enddo !imt
       enddo !imp
!       write(kfam,'(2x,"theta=",1f6.2)') th 
       write(kfam,'(2x,1f6.2)') th ! changed in v2.2e
       write(kfam,228) ((((fam(ith,im,imt,imp,imtp),imtp=1,nmtf),
     &  imp=1,nmf),imt=1,nmti),im=1,nmi)
228    format(1P,6E12.4)
!       write(kfam,*) ' '
       enddo !angles
777    format(f8.3,1e12.3,1f10.4,1e12.3)
!       write(0,*) 'Exited angle-dependent'


c *** Differential cross sections for each state
       do ith=1,nth
        th = thmin + dth*(ith-1)
        if (th<theps) th=theps
        factor=1./nmi/nmti
        thrad=th*pi/180        
        xs=0d0
        do im=1,nmi
        do imt=1,nmti
        do imp=1,nmf
        do imtp=1,nmtf
         xs=xs+10*factor*abs(fam(ith,im,imt,imp,imtp))**2
        enddo !imtp
        enddo !imp
        enddo !imt
        enddo !imp
        if ((iex.eq.1).and.(itex.eq.1).and.(zp*zt.gt.0)) then 
            s2=sin(thrad/2.)**2
            xsruth=10*(etai/(2.*kcmi*s2))**2  ! point Coulomb(mb/sr)
            write(kxs,800) th,xs/xsruth,xs!,xsruth*10
            write(kxs+1,800) th,xs/xsruth,xs!,xsruth*10
        else
            write(kxs,800) th,xs
            if (iex.le.10) write(kxs+(iex-1)*(ntex+1)+itex,800) th,xs
        endif
        enddo !ith
800     format(1x,1f8.3,2g15.5)
        write(kxs,*)'&'
        
        if(allocated(fam))deallocate(fam)
        if(allocated(pl))deallocate(pl)
        if(allocated(delcf))deallocate(delcf)
        if(allocated(lfv))deallocate(lfv)
        if(allocated(ampl))deallocate(ampl)
        
      enddo !itex
      enddo ! iex 


c *** Angle-integrated cross sections 
      sigex(:)=0
      sigtex(:)=0d0
      sigextex(:,:)=0d0
      xstot = 0.
      jpi   = jpch(iexgs)
      nmi   = 2*jpi+1
      jti   = jtgs
      nmti  =nint(2*jtgs+1)
      icc=0
      njt  =jtmax-jtmin+1
      do ijt=1,njt
      jtot=jtmin+dble(ijt-1)
      xsj=0
      xrj=0
      do partot=1,-1,-2
      xs=0
      xr=0
       icc=icc+1 ! CC set (total J/pi)
       if (jptset(icc)%interp) cycle
!       partot=jptset(icc)%partot
!       jtot  =jptset(icc)%jtot
       nch   =jptset(icc)%nchan
       do ni=1,nch
       if (jptset(icc)%idx(ni).ne.iexgs) cycle ! no inc waves 
       if (jptset(icc)%idt(ni).ne.1) cycle !Needs rephrasing MGR 
       li= jptset(icc)%l(ni)
       jlpi=jptset(icc)%jlp(ni)
       jti=jptset(icc)%jt(ni)
       do nf=1,nch
       iex  = jptset(icc)%idx(nf)
       idt  = jptset(icc)%idt(nf)
       exc  = exch(iex)
       etarg= jptset(icc)%ext(nf)
       tkch=ecmi+exch(iexgs)-exc-etarg
       if (tkch.lt.0) cycle    !closed channel
       kcmf = sqrt(2*mupt*tkch)/hc
       lf   = jptset(icc)%l(nf)
       jpf  = jptset(icc)%jp(nf) 
       jlpf=jptset(icc)%jlp(nf)
       jtf=jptset(icc)%jt(nf)

       if ((iex.eq.iexgs).and.(idt.eq.1)) then ! elastic channel 
         written(ksmat)=.true.
         smat=smats(icc,ni,nf)
         xr=xr+10*(pi/kcmi**2)*(2.*jtot+1)/nmi/nmti*(1-abs(smat)**2)
         write(ksmat,'(1f6.1,2x,3i3,3f6.1," ->",3i3,3f6.1,2x,3f14.8)')
     &   jtot,iexgs,1,li,jpi,jlpi,jti,iex,idt,lf,jpf,jlpf,jtf,smat,
     &   abs(smat)
       else      
         smat=smats(icc,ni,nf)! *sqrt(kcmf/kcmi)
         xsaux=10*(pi/kcmi**2)*(2.*jtot+1)/nmi/nmti
     &    *abs(smat)**2*kcmf/kcmi
         xs=xs + xsaux 
         sigex(iex)=sigex(iex) + xsaux
         sigtex(idt)=sigtex(idt)+ xsaux
         sigextex(iex,idt)=sigextex(iex,idt)+ xsaux
         if (verb.gt.1) then
         write(ksmat,'(1f6.1,2x,3i3,3f6.1," ->",3i3,3f6.1,2x,3f14.8)')
     &   jtot,iexgs,1,li,jpi,jlpi,jti,iex,itex,lf,jpf,jlpf,jtf,smat,
     &   abs(smat)
         endif
       endif
       enddo !nf
       enddo !ni
       xstot= xstot + xs
       xsj  = xsj   + xs
       xrj  = xrj   + xr
       write(*,'(5x," J/pi=",1f5.1,a1,3x,
     &   " xs(inel)=",1e12.3, " mb"," reac=",1e12.3)') 
     &   jtot,parity(partot+2),xs,xr
      enddo !partot
! COMMENTED, becase we have printed the x-sections after each J/pi set
!      write(ksj, '(1x,1f5.1,3f10.4,2x,i2)') jtot,xrj-xsj,xrj,xsj
      enddo !ijt (JT)
      write(*,780) xstot
780   format(/,3x," => total inel+bu x-section=",1f8.3," mb")


c *** Angle-integrated cross section for each projectile state
      write(*,'(/,2x,"o Angle-integrated xs for each state:")')
      jpold=-1 ; parold=-2
      do iex=1,nex
        partot=parch(iex)
        exc   =exch(iex)
        jpf   =jpch(iex)
        if (iex.eq.1) jpold=jpf
        if ((iex.gt.1).and.(jpf.ne.jpold).and.(partot.ne.parold))
     & write(130,*) '&'
        jpold=jpf
        write(130,795) exc, sigex(iex), iex,jpf,partot
        write(*,790) iex,jpf,parity(partot+2),exc, sigex(iex)
789     format(8x,"#",i4,3x,1f5.1, a1,f10.3,"->")
790     format(8x,"#",i4,3x,1f5.1, a1,f10.3,"->",f12.6,' mb')
794     format(1x,1f10.3,2x,12x,1x,i3,1f5.1,2x,i2)
795     format(1x,1f10.3,2x,1f12.6,1x,i3,1f5.1,2x,i2)
      enddo

!   Angle-integrated cross section for each target state
      write(*,'(/,2x,"o Angle-integrated xs for each target state:")')
      jpold=-1 ; parold=-2
      do idt=1,ntex+1
        partot= tset(idt)%partarg
        exc   = tset(idt)%extarg
        jtf   =tset(idt)%jtarg
        jpold=jpf
        write(1130,795) exc, sigtex(idt), idt,jtf,partot
        write(*,790) idt,jtf,partot,exc, sigtex(idt)
      enddo

!   Angle-integrated cross section for each target and projectile state
      write(*,'(/,2x,"o Angle-integrated xs for each target and 
     & projectile state:")')
      jpold=-1 ; parold=-2
      do iex=1,nex
        partot=parch(iex)
        exc   =exch(iex)
        jpf   =jpch(iex)
        if (iex.eq.1) jpold=jpf
        if ((iex.gt.1).and.(jpf.ne.jpold).and.(partot.ne.parold))
     & write(130,*) '&'
        jpold=jpf
        
        do idt=1,ntex+1
          partarg= tset(idt)%partarg
          extarg   = tset(idt)%extarg
          jwtarg   =tset(idt)%jtarg
          write(1230,794) exc, iex,jpf,partot
        write(*,789) iex,jpf,parity(partot+2),exc
          write(1230,795) extarg, sigextex(iex,idt), idt,jtf,partot
          write(*,790) idt,jtf,partot,extarg, sigextex(iex,idt)
        enddo
      enddo

1000  continue

      if (allocated(famps)) deallocate(famps)



      end subroutine





! two dimension interpolation function
! based on the method of fival function
! f(xbar,ybar) from f(x,y)
      function fint2db_jin(xtab,ytab,fxytab,xbar,ybar,nnx,nny,alpha)
        implicit none
        integer :: nnx,nny ! dimesnion of x and y
        real*8 ::  xtab(nnx),ytab(nny)  ! vector of x and y
!        real*8 :: fxytab(nnx,nny) ! original function
        complex*16 :: fxytab(nnx,nny) ! original function
        real*8 :: xbar,ybar,f
        complex*16:: fint2db_jin ! interpolted points and value
        real*8 :: alpha
!        real*8 :: fival2d(nny)
        complex*16 :: fival2d(nny)
        complex*16 :: cfival

        fint2db_jin=0d0


        if(xbar>xtab(nnx) .or. xbar<xtab(1) .or. ybar>ytab(nny)
     +                  .or. ybar<ytab(1) ) return


        call fival2c(xbar,xtab,fxytab,nnx,nny,alpha,fival2d)

        fint2db_jin = cfival(ybar,ytab,fival2d,nny,alpha)


      end function



************************************************************************
*     subroutine modified for the two dimension interpolation
*     REAL 4-point lagrange interpolation routine.
*     interpolates thr FUNCTION value fival at point r from an
*     array of points stored in fdis(ndm). this array is assumed
*     to be defined such that the first element fdis(1) CONTAINS
*     the FUNCTION value at r=xv(1) and xv(2 .. ndm) are monotonically
*     increasing.
************************************************************************
      subroutine fival2c(r,xv,fdis,ndm,ndm2,alpha,fival2d)
      IMPLICIT REAL*8(A-H,O-Z)
!      REAL*8 fdis(ndm,ndm2),y1(ndm2),y2(ndm2),y3(ndm2),y4(ndm2)
      complex*16 fdis(ndm,ndm2),y1(ndm2),y2(ndm2),y3(ndm2),y4(ndm2)
      DIMENSION xv(ndm)
!     real*8 fival2d(ndm2)  !ndm2 stands for the other dimension of 2-dimension interpolation
      complex*16 fival2d(ndm2)  !ndm2 stands for the other dimension of 2-dimension interpolation
      IF(r.GT.xv(ndm)) go to 9
      DO 5 k=1,ndm-2
 5    IF(r.LT.xv(k)) go to 6
      k=ndm-2
 6    nst=MAX(k-1,1)
      x1=xv(nst)
      x2=xv(nst+1)
      x3=xv(nst+2)
      x4=xv(nst+3)
      y1=fdis(nst+0,:)
      y2=fdis(nst+1,:)
      y3=fdis(nst+2,:)
      y4=fdis(nst+3,:)
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
      fival2d=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4
      RETURN
 9    fival2d=fdis(ndm,:) * EXP(alpha*(xv(ndm)-r))
      RETURN
      END subroutine





c **** -------------------------------------------------------------------
c     2 dimension aitkin interpolation routine. 
c     Note that mesh points (xtab,ytab) must be in increasing order.
c *** ---------------------------------------------------------------------
      function fint2d(xtab,ytab,fxytab,xbar,ybar,nnx,nny,nord,mmx)
      implicit none
      integer ix,iy,nnx,nny,mmx,ixref,iyref,nord,i,j,ii,jj
      integer min,max,mm,num,nx,ny,nnn
      integer nxy(2),nbg(2)
!      implicit real*8 (a-h,o-z),integer*4(i-n)
      real*8 xbar,ybar
      real*8 xytab(2,100),x(10),y(10),
     1 xybar(2),xtab(mmx),ytab(mmx)
      complex*16:: fxytab(mmx,mmx),fxy(10,10)
      complex*16 p1,p2,p3,p4,u,t,fint2d

      write(95,'(50f12.5)') xtab(1:5)
      write(95,'(50f12.5)') ybar,ytab(1:nny)
      write(95,'(50f12.5)') fxytab(1:5,1:1)
      
      
      nnn=nord+1
c          set up arrays for loop
      do 40 j=1,nnx
 40   xytab(1,j)=xtab(j)
      do 45 j=1,nny
 45   xytab(2,j)=ytab(j)
      nxy(1)=nnx
      nxy(2)=nny
      xybar(1)=xbar
c          begin loop to determine the (order+1) x and y points over
c          which interpolation is made; 1=x, 2=y.
      xybar(2)=ybar
      do 10 i=1,2
30    num=1
      if(xybar(i).lt.xytab(i,1))go to 85
      num=nxy(i)-nnn+1
      if(xybar(i).gt.xytab(i,nxy(i)))go to 85
50     min=1
      max=nxy(i)
      num=nxy(i)/2
55    if(max-min.lt.4)goto70
      if(xytab(i,num)-xybar(i))60,82,61
60     min=num
          num=num+(max-min+1)/2
      goto55 
 61   max=num
      num=(max-min+1)/2
      goto55
  70   num=max
  71  if(xytab(i,num)-xybar(i))82,82,81
 81   num=num-1
      goto71
 82   num=max0(1,num-nnn/2)
      num=min0(num,nxy(i)-nnn+1)
 85   nbg(i)=num
 10   continue
c          end loop.  set up x, y, function arrays of required
c          (order+1)**2 points
      nx=nbg(1)
      ny=nbg(2)
      do 20 i=1,nnn
      x(i)=xytab(1,nx+i-1)
      y(i)=xytab(2,ny+i-1)
      do 20 j=1,nnn
      fxy(i,j)=fxytab(nx+i-1,ny+j-1)
 20   continue
c          do interpolation
 90   do 95 ii=2,nnn
      do 95 jj=ii,nnn
      do 95 mm=ii,nnn
         if (abs(y(mm)-y(ii-1)).lt.1e-6) then
            write(0,*)'error 1 in fint2d';stop
         endif
         if (abs(y(ii-1)-xybar(2)).lt.1e-6) then
            write(0,*) y(ii-1),xybar(2)
            write(0,*)'error 2 in fint2d';stop
         endif
         if (abs(x(jj)-xybar(1)).lt.1e-6) then
            write(0,*)x(jj),xybar(1)
            write(0,*)'error 3 in fint2d';stop
         endif
         
      fxy(jj,mm)=((fxy(ii-1,ii-1)*(x(jj)-xybar(1))-fxy(jj,ii-1)*(x(ii-1)
     1 -xybar(1)))*(y(mm)-xybar(2))-(fxy(ii-1,mm)*(x(jj)-xybar(1))
     2 -fxy(jj,mm)*(x(ii-1)-xybar(1)))*(y(ii-1)-xybar(2)))
     3 /(x(jj)-x(ii-1))/(y(mm)-y(ii-1))
  95   continue
      fint2d=fxy(nnn,nnn)
      return
      end




c *** ------------------------------------------------------
c ***  Simple 2-dim interpolation
c ***  ------------------------------------------------------
      function fint2dd(xtab,ytab,fxytab,xbar,ybar,nnx,nny,mmx)
      implicit none
!      implicit real*8 (a-h,o-z)
      integer ix,iy,nnx,nny,mmx,ixref,iyref
      real*8 xtab(nnx),ytab(nny)
      real*8 xbar,ybar
      complex*16 p1,p2,p3,p4,u,t,fint2dd,fxytab(mmx,nny)

!      write(95,'(100e12.4)') xtab(1:neset)
!      write(96,'(100e12.4)') fxytab(1:nnx,2)
!      return
      fint2dd=0d0
     
      
      do iy=1,nny-1
         iyref=iy
         if((ybar.ge.ytab(iy)).and.(ybar.le.ytab(iy+1)))goto 120
      enddo
      

 120  do ix=1,nnx-1
         ixref=ix
         if ((xbar.ge.xtab(ix)).and.(xbar.le.xtab(ix+1)))goto 100
      enddo
!      write(*,*)'fint2db error'
      return
 100  p1=fxytab(ixref,iyref)
      p2=fxytab(ixref+1,iyref)
      p3=fxytab(ixref,iyref+1)
      p4=fxytab(ixref+1,iyref+1)

      

      t=(xbar-xtab(ixref))/(xtab(ixref+1)-xtab(ixref))
      u=(ybar-ytab(iyref))/(ytab(iyref+1)-ytab(iyref))
      
!      fint2d=(p1+p2+p3+p4)/4d0
c bilinear interpolation
       fint2dd=(1-t)*(1-u)*p1+t*(1-u)*p2+t*u*p3+(1-t)*u*p4
!      write(96,'(10f12.4)') ytab(iyref),ybar,ytab(iyref+1)
!      write(96,'(3f12.4)') xtab(ixref),xbar,xtab(ixref+1)
!      write(96,'(2i4,5e12.4/)')ixref,iyref,p1,p2,p3,p4,fint2db
      end

c *** ------------------------------------------------------
c ***  Read overlaps between PS's and scattering states
c      Read from file g+.n with n=jtot-0.5 for partot>0
c      and from g-.n for partot<0.
c      This subroutine is called if icore<0
c ***  ------------------------------------------------------
      subroutine solapread(gsolap,partot,inc,nsetmax,ncont,jtot,dec,e0)
      implicit none
      integer:: partot,inc,nsetmax,ncont
      real*8:: jtot, dec, e0
      complex*16:: gsolap(nsetmax,ncont)
      integer:: i, ie, ich
      integer:: nk, nset
      real*8:: eminr, emaxr, alpha, greal(inc), gimag(inc)
      real*8, allocatable:: ein(:)
      complex*16:: cfival
      complex*16, allocatable:: gvec(:)
      character*40 filename, line


      alpha=0

      if (partot.eq.1) then
      write(filename,'("g+.",i0)') int(jtot-0.5) 
      else 
      write(filename,'("g-.",i0)') int(jtot-0.5) 
      endif
      write(*,*) 'Reading from ', filename
      open(52,file=trim(filename))
      read(52,'(a)') line
      read(line,*) eminr, emaxr, nk, nset
      
      if (nset>nsetmax) then
      write(*,*)'*ERROR* nset>nsetmax:',nset,'>',nsetmax
      stop
      endif
      
      allocate(ein(nk),gvec(nk))
      
      do i=1,nset
      read(52,*) line
      do ie=1,nk
      read(52,*) ein(ie),(greal(ich),gimag(ich),ich=1,inc)
      gvec(ie)=cmplx(greal(inc),gimag(inc))
      end do !ie
 
!    interpolate in THOx grid    
      do ie=1,ncont
      gsolap(i,ie)=cfival(e0+(ie-1)*dec,ein(:),gvec(:),nk,alpha)
      enddo !ie
      enddo !i

      close(52)
      return
      end subroutine


