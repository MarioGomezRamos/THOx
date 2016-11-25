      module nmrv
      logical:: debug    ! if true,print debug information
      integer:: verb     ! determines the amount of trace output
      integer:: nch      ! number of channels
      integer:: nopen    ! number of open channels ! added in v2.4
      integer:: nr       ! number of radial points
      real*8 :: h        ! radial step for integration
      real*8 :: conv     ! conversion factor 2*mu/hbar^2
      real*8 :: mu       ! reduced mass
      real*8 :: rmin     ! min radius for integration
      real*8 :: hort     ! step for orhogonalization
      real*8 :: rmort    ! max radius for orhogonalization
      real*8 :: cutr=-20 ! cutoff
      real*8, allocatable :: rvec(:)
      complex*16,allocatable,target :: vcoup(:,:,:)   ! coupling matrix formfactors
      real*8,allocatable       :: ech(:)   ! channel energies
      real*8,allocatable       :: kch(:)   ! channel wave-numbers      
      complex*16,allocatable,target:: y(:,:,:) !
      integer,allocatable      :: ql(:)    ! l-values for channels
      end module 



!-------------------------------------------------------------------------------
c***Auxiliary routine for calculation of scat. wf. for arbitrary jpi set and a 
!   range of continuum energies (NOTE normalization below!)
!
!   o nset     =jpiset
!   o ecmi,ecmf=energy range (with respect to corresponding threshod )
!   o ne       =nb of energies
!   o nset     =j/pi set of projectile (c+v system)
!   o nchan    =number of channels of this j/pi set
!   o inc      =incoming channel
!   o energy   =energies evenly spaced in energy (T) or momentum (F)
!   o wfc      =continuum wfs with norm:
!               <k|k'> =(pi/2) delta(k-k')  
!   o phase    =phase-shifts (complex)
!   o wftype   =1 : Scattering states (complex in general)
!              =2 : Real CC wfs
!-------------------------------------------------------------------------------
      subroutine wfrange(nset,nchan,inc,ecmi,ecmf,ne,energy,wfc,
     & smate,delta)
      use nmrv,only:nch,ech,vcoup,hort,cutr,nopen
      use constants
      use channels, only:jpiset,jtot,ql,qspl,qj,qspj,qjc,exc,
     &                   cindex,spindex,lmax
      use parameters, only:  maxchan
      use sistema
      use potentials, only: ccmat
      use wfs, only: nr,dr,rvec,wfsp,rint,rmin
      use globals, only: written
      use forbidden ! added in v2.2g (TESTING)
      implicit none
c     ------------------------------------------------------------
      character*5 jpi
      logical :: info,energy,realcc,debug
c ... 
      integer ir,n,nchan,method,ne,ie,nrint
      integer :: nset,inc,partot,ifail,wftype
      real*8  :: ecm,ecmi,ecmf,de,excore,econt
      real*8  :: dk,kf,ki,kcm,kcont,krm,tkch
      real*8  :: z12,rm,factor,r0
      real*8  :: deladd,deltap,deltai,delta(ne),cph(0:lmax)
      complex*16 :: phase(nchan),smat(nchan),smate(ne,maxchan)
      complex*16 :: wfc(ne,nchan,nr),wf(nchan,nr),wfop(nchan,nchan,nr)
      complex*16 :: ci
c *** R-matrix solutions -----------------------------------------
      integer*4, parameter:: ninc=1,nr2=0,nmax=100
      logical :: iftr,twf,ifrmat
      integer*4 :: nbas,ich,ichp,nvc(inc)
      real*8 :: r,zrma(nmax),eta(nchan),kch(nchan),jci,jcf,rmatch
      real*8,parameter :: alpha=0.
      real*8 :: aux,auxp,fival
      complex*16 :: cfival,caux,phc,ph2,cauxp
      complex*16 :: cpot(nmax,nchan,nchan),cu(nchan,nchan),
     &              cf(nmax,nchan,ninc),faux(nr)
      real*8 fc(500),dfc(500),gc(500),dgc(500),xfac(nchan,nchan)
      complex*16 cfx(nchan)
c *** Pauli blocking (TESTING) -----------------------------------------
      integer ispi,ip,ncp2,ii,ncorei,i,j
      integer:: lp,np,li,linc
      real*8:: ji,jp,rp
      real*8,allocatable:: wfpau(:)  
      integer :: ns
      complex*16, allocatable:: cpnl(:,:,:)     
c *** Extrapolation of WFS
      real*8,parameter :: zero=0.0
      real*8 :: kron
      real*8, dimension(0:11):: f,g,gp,fp
      real*8, dimension(1:12+1):: f1,fp1
      complex*16:: ch
c *** Eigenphases
      integer:: iwrite=90
      real*8 :: eph(nchan),ener(ne),phz(nchan,ne)
      real*8 QUAD(NCHAN),QUAD1(NCHAN)
c     ------------------------------------------------------------------
      
      if (ne.lt.2) stop 'wfconts: NE must be >1'

c *** Initialize some variables 
      rm=av*ac/(av+ac)
      factor=(2*amu/hc**2)*rm
      wfc   =cmplx(0.,0.)
      info  =.false. ! silent output
      deladd=0.
      ci=cmplx(0.,1.)
      z12=zc*zv
      wftype=1  !(1=scat. states; 2=real CC states)
      debug=.false.
c *** -------------------------

c     ------------------------------------------------------------
      if (npauli.gt.0) then
        ifrmat=.true.  ! use R-matrix routine by P.Descouvemnt
      else 
        ncp2=0
        ifrmat=.false. ! use Numerov
      endif 
c     ------------------------------------------------------------

      if (info) then
        write(*,100)ne,jpi(jtot,partot),inc,ecmi,ecmf
100     format(3x,"Generating", i4, " continuum wfs for J/pi=",
     &  a5,", inc. chan.",i3,
     &  " in [ECMmin=",1f6.3," ECMmax=",1f6.2," MeV]")
      endif 


      jtot    =jpiset(nset)%jtot
      partot  =jpiset(nset)%partot
      excore  =jpiset(nset)%exc(inc)
!      bastype =jpiset(nset)%bastype

      if (nchan.ne.jpiset(nset)%nchan) then
         write(*,*)'Wrong number of channels passed to wfrange!'
      endif
      nch=nchan

      if (allocated(vcoup)) deallocate(vcoup)
      if (allocated(ech))   deallocate(ech)
      allocate(vcoup(nch,nch,1:nr))

      vcoup(:,:,:)=0d0 
      allocate(ech(nchan))
      ql(1:nch)    =jpiset(nset)%lsp(1:nch)  
      ech(1:nch)   =jpiset(nset)%exc(1:nch)
      qj(1:nch)    =jpiset(nset)%jsp(1:nch)
      qjc(1:nch)   =jpiset(nset)%jc(1:nch)
      cindex(1:nch)=jpiset(nset)%cindex(1:nch)


      call coefmat(nch)
      vcoup(:,:,1:nr)=ccmat(:,:,1:nr)/factor

 
      if (energy) then
        de=(ecmf-ecmi)/(ne-1)
      else
        kf=sqrt(factor*ecmf)
        ki=sqrt(factor*ecmi)
        dk=(kf-ki)/(ne-1)
      endif


c Initialize variables for Numerov integration
      r0=rvec(1)
      method=4       ! enhanced Numerov as used in Fresco
      hort  =0       ! no stabilization
      info  =.false. ! silent output
      cutr  =-50 
      if ((rint.gt.0).and.(rint.lt.rvec(nr))) then
        nrint=ceiling((rint-rmin)/dr)+1
        rmatch=rint
      else
        rmatch=rvec(nr)
        nrint=nr
      endif
c ------------------------------------------------------


c *** Initialize R-matrix ----------------------------------
      if (ifrmat) then
      nvc(1) =inc       ! incoming channel
      nbas   =100       ! Lagrange functions per radial interval
      rmatch =rvec(nr)  ! rvec(nr)
      twf    =.true.    ! generate also WFS
      ncp2   =0         ! non-local couplings
      ns     =1         ! nb. of intervals
      if (npauli.gt.0) ns=1

      if(nbas*ns.gt.nmax)then
      print*,'wfrange: nmax too small'
      stop
      end if
      call rmat_ini(nbas,ns,rmatch,zrma)

      do ich=1,nch
      do ichp=ich,nch
      phc=(0d0,1d0)**(ql(ichp)-ql(ich))
      jci= qjc(ich)
      jcf= qjc(ichp) 
      ph2=(-1)**NINT((JCI-JCF - ABS(JCI-JCF))/2.)
      faux(1:nr)=ccmat(ich,ichp,1:nr)
      do 4 ir=1,nbas*ns
      r=zrma(ir)
      caux=cfival(r,rvec,faux,nr,alpha)
      cpot(ir,ich,ichp)=caux
      if (ich.ne.ichp) cpot(ir,ichp,ich)=caux
    4 continue
      enddo !ich
      enddo !ichp
c      do ir=1,nbas*ns
c        write(1,'(1f8.3,2x,50f12.8)')zrma(ir),cpot(ir,1,1)
c      enddo
c      do ir=1,nr
c        write(2,'(1f8.3,2x,50f12.8)')rvec(ir),ccmat(1,1,ir)
c      enddo      
      endif ! rmatrix initialization


c *** Pauli blocking operator (used only with R-matrix solutions!) --------
      if (npauli.gt.0) then 
        ncp2=nbas**2
        allocate(cpnl(ncp2,nch,nch)) ! Non-local coupling matrix
        cpnl=(0.,0.)
        allocate(wfpau(nr))
	do ip=1,npauli
	lp=paul(ip) 
        jp=pauj(ip)
 	np=paun(ip)

        do ich=1,nch
        li=ql(ich)	
        ji=qj(ich)
        jci=qjc(ich)
        ncorei=cindex(ich)
        ispi=spindex(ich)
	if (li.ne.lp) cycle
	if (ji.ne.jp) cycle
        write(*,*)' -> remove pauli ip',ip,'chan=',ich,
     &   'core=',ncorei,' shift=',pshift(ip,ncorei) 
       
        if (wfpread(ip)) then
          wfpau(:)=wfpauli(ip,:)
        else
          wfpau(:)=wfsp(ispi,paun(ip),:) !Pauli forbidden 
        endif
        write(97,'(f6.3,2x,g12.6)') (rvec(ir),wfpau(ir),ir=1,nr)
        write(97,*)'&'

        ii=0
        do 5 i=1,nbas
        do 5 j=1,nbas
        ii=ii+1
        r=zrma(i)
        rp=zrma(j)
        aux =fival(r, rvec,wfpau,nr,alpha)
        auxp=fival(rp,rvec,wfpau,nr,alpha)
        cpnl(ii,ich,ich)=  cpnl(ii,ich,ich)
     &                  + aux*auxp*pshift(ip,ncorei)*r*rp*factor
    5   continue
        enddo !ich
        enddo !ip

       endif ! Pauli
c ----------------------------------------------------------------------------


c ---------------------------------------------------------------------------------
c Energy/momentum grid
      do ie=1,ne
         wf=cmplx(0.,0.)
         if (energy) then 
         ecm=ecmi + de*(ie-1.)
         else
         kcm=ki+dk*(ie-1)
         ecm=kcm**2/factor
         endif
         if (ecm.lt.0) cycle
         ener(ie)=ecm

         select case(wftype) 
         case (2) 
c ... Real multichannel states from Numerov solutions ...........
          if (ie.eq.1) print*,'calling erinwrc'
          call schcc_erwinrc(nch,ecm,zc*zv,inc,ql,factor,dr,r0,
     &                    nrint,wfop,method,info)
          wfc(ie,1:nch,1:nrint)=wfop(1:nch,inc,1:nrint)
          if (rint.gt.0.and.rint.lt.rvec(nr)) then
             write(*,*)' ** warning: continuum wfs set to 0 from rint=',
     &       rint, ' to rmax=', rvec(nr)
          endif
          cycle ! go to next energy
 
         case (1)  ! CC Bins
c *** Scattering states with numerov solution ................
         if (.not.ifrmat) then
!          call schcc_erwin(nch,ecm,zc*zv,inc,ql,factor,dr,r0,
!     &                    nr,wf,phase,smat,method,info)

          call schcc_erwin(nch,ecm,zc*zv,inc,ql,factor,dr,r0,
     &                    nrint,wf,phase,smat,method,info)

!          call schcc(nch,ecm,zc*zv,inc,ql,factor,dr,r0,
!     &              nrint,wf,phase,smat,info)

           wfc(ie,1:nch,1:nrint)=wf(1:nch,1:nrint)
!           delta(ie,:)=phase(:)*pi/180.

         else    
c *** Scattering states with R-matrix solution ................                                              
          nopen=0
          do ich=1,nch
           aux=ecm+ech(inc)-ech(ich)
           kch(ich)=sqrt(factor*abs(aux))
           eta(ich)=factor*zc*zv*e2/kch(ich)/2.
           if (aux.lt.0) kch(ich)=-kch(ich)
          enddo !ich
          call rmatrix(nch,ql,kch,eta,rmatch,nbas,ns,cpot,cu,nmax,nch,
     &         nopen,twf,cf,nmax,nch,ninc,nvc,ncp2,cpnl)

c ... interpolate wfs from Lagrange to uniform grid
          do ich=1,nch
          if (kch(ich)>0.) then ! open channel -> get phase-shifts and S-matrix
          phase(ich)=(0.,-0.5)*LOG(cu(ich,inc)
     &               *sqrt(kch(inc)/kch(ich)))*180/pi        
          smat(ich)=cu(ich,inc)*sqrt(kch(inc)/kch(ich))
          endif
c ... store wfs
          do ir=1,nrint
          r=rvec(ir)
          if (r.lt.zrma(nbas)) then
             wfc(ie,ich,ir)=cfival(r,zrma,cf(:,ich,inc),nbas,alpha)
     &                      *0.5*(0.,1.)
          endif                         
          enddo !nrint
          enddo !nch

c ...................................................................
          if (ie.eq.1) then
!          write(*,*)'Nopen=',nopen,' Elastic S-matrix=',cu(1,inc)
!          do ich=1,nch
!             write(*,*)'Chan=',ich,' S-mat=',
!     &      cu(ich,inc)*sqrt(kch(inc)/kch(ich))
!          enddo
          do  ir=1,nbas*ns
          write(501,1002)zrma(ir),
     &      (cf(ir,ich,1)*(0.,-1.),ich=1,nch)
1002      format(f8.3,10(2es12.4,2x))
          enddo
          write(501,*)'&'
          endif

         ENDIF ! Choose integration method (Numerov / R-matrix)

        case default
         write(*,*)'wfrange: wftype=',wftype,' not valid'
         stop
          
        end select 
c ....................................................................


c ... Extrapolate wfs from rint to rmax 
      if ((rint.gt.0).and.(rint.lt.rvec(nr)).and.wftype.ne.4) then
!      write(0,*)'extrapolating from ',rint,' to=',rvec(nr)
      do ich=1,nch
      li=ql(ich) 
      if (li.gt.10) stop'increase lmax in wfrange!'
      tkch=ecm+ech(inc)-ech(ich) 
      kch(ich)=sqrt(factor*abs(tkch))
      eta(ich)=factor*z12*e2/kch(ich)/2.      
      do ir=nrint,nr
      krm=kch(ich)*rvec(ir)
      if (tkch.gt.0d0) then   ! open channel
        call coul90(krm,eta(ich),zero,li,f,g,fp,gp,0,IFAIL)
        ch=cmplx(g(li),f(li))  ! H^(+)
        wfc(ie,ich,ir)=(0.,0.5)*(conjg(ch)*kron(ich,inc)-smat(ich)*ch)
        if (ifail.ne.0) then 
        write(*,*) 'coul90: ifail=',ifail; stop
        endif
      else                    ! closed channel
       call whit(eta(ich),rvec(ir),kch(ich),tkch,li,f1,fp1,0)
       ch = smat(ich)*f1(li+1)* (0.,.5) ! (i/2)*C*Whittaker
       wfc(ie,ich,ir)=ch
      endif

      enddo !ir
      enddo !chans
      endif ! rlast > rmax
c ...............................................


c ... Add Coulomb phase 
c    (could be done within ERWIN, but I do it here in order not to affect other parts of the code using ERWIN
      if (wftype.ne.4) then
      linc=ql(inc) 
      kch(inc)=sqrt(factor*abs(ecm))
      eta(inc)=factor*z12*e2/kch(inc)/2.    
!      print*,'z12,kch(inc),eta,linc=',z12,kch(inc),eta,linc
      call coulph(eta(inc),cph,linc)
!      print*,'cph(linc)=',cph(linc)
      phc=exp(ci*cph(linc))
! AMM: should we restrit this to open channels???
      wfc(ie,:,:)= phc*wfc(ie,:,:)
      endif ! wftype

c ... If wfs are to be 'normalized' as <k|k'> = delta(k-k'), require additional factor sqrt(2/pi)
c      wf(:,:)=wf(:,:)*sqrt(2./pi)

c Store S-matrix and phase-shifts 
      smate(ie,1:nch)=smat(1:nch)

c Avoid pi jumps in phase-shifts
      deltai=phase(inc)
      if (ie.eq.1) then
          deltap=deltai
      else
         if(deltai<deltap-90) deladd=deladd+180
         if(deltai>deltap+90) deladd=deladd-180
      endif
      deltap=deltai
      deltai = deltai + deladd
      delta(ie)=deltai !*pi/180.

!!!! Diagnostic TEST PHASE SHIFT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (ie.eq.1) then
        written(45)=.true.
        write(45,205)jpi(jtot,partot),inc,ql(1:nchan)
205     format("#Phase-shifts for J/pi=",a5,"  Inc. chan.=",
     &         i3, " ql=",5i3)
        endif
        write(45,'(1f6.3,2x,10f12.5)') ecm, (phase(ich),ich=1,nchan)
        

        
!!!! Diagnostic TEST WF
        if ((ie.eq.0).and.(inc.eq.1)) then
!        write(*,*)'Ecm=',Ecm,' Elastic S-matrix=',smat(1)
!        write(*,*)'kch=',kch(1:nch)
!        do ich=1,nch
!        write(*,*)'Chan=',ich,' S-mat=',
!     &    smat(ich)
!        enddo

        write(500,200)jpi(jtot,partot),ecm,inc,ql(1:nchan)
200     format("#WF for J/pi=",a5," Ecm=",1f7.3," MeV  Inc. chan.=",i3,
     &          " ql=",5i3)

        
        do ir=1,nr
        write(500,'(1f8.3,2x,50f12.8)') rvec(ir),(wfc(ie,n,ir),n=1,nch)
        enddo
        write(500,*)'&'
      endif
!        write(45,'(1f10.3,3x,10g14.6)') ecm,(phase(n),n=1,nch) 
!!!! TEST
      enddo ! energy loop (ie)

c Sort eigenphases
!      write(*,*)'Try making ephases continuous'

!      call ESMOOV(NCH,QUAD,QUAD1,PHZ,ENER,NE,IWRITE)
!      do ie=1,ne
!      write(94,'(1f8.4,20f12.6)') ener(ie),
!     &   (phz(n,ie),n = 1, nch)
!      enddo ! ie 

      write(450,*)'&'
      deallocate(ccmat,vcoup)
      end subroutine



!-------------------------------------------------------------------------------
!  Auxiliary routine for calculation of scat. wf. for arbitrary jpi set and ecm
!-------------------------------------------------------------------------------
      subroutine test_cont(nset,nchan,inc,ecm,wf)
      use nmrv,only:nch,ech,vcoup,hort
      use constants
      use channels, only: jpiset, ql,qspl, qj, qspj,qjc,exc,cindex
      use sistema
      use potentials, only: ccmat
      use wfs, only: nr,dr,rvec
      implicit none
c     ------------------------------------------------------------
      logical :: info
      integer ir,n,nchan,method
      integer nset,inc,partot
      real*8  ecm,z12,jtot,rm,factor,r0
      complex*16 :: phase(nchan),smat(nchan)
      complex*16 :: wf(nchan,nr)
c     ------------------------------------------------------------
      rm=av*ac/(av+ac)
      factor=(2*amu/hc**2)*rm

      jtot    =jpiset(nset)%jtot
      partot  =jpiset(nset)%partot
      if (nchan.ne.jpiset(nset)%nchan) then
         write(*,*)'Wrong number of channels passed to test_cont!'
      endif
      nch=nchan

      if (allocated(vcoup)) deallocate(vcoup)
      if (allocated(ech))   deallocate(ech)
      allocate(vcoup(nch,nch,1:nr)) ! changed in v2.2

      vcoup(:,:,:)=0d0 
      allocate(ech(nch))
      ql(1:nchan)    =jpiset(nset)%lsp(1:nchan)  
      ech(1:nch)     =jpiset(nset)%exc(1:nch)
      qj(1:nchan)    =jpiset(nset)%jsp(1:nch)
      qjc(1:nchan)   =jpiset(nset)%jc(1:nch)
      cindex(1:nchan)=jpiset(nset)%cindex(1:nch)

      call coefmat(nch)
      vcoup(:,:,1:nr)=ccmat(:,:,1:nr)/factor

      r0=rvec(1)
      method=4       ! enhanced Numerov as used in Fresco
      hort  =0       ! no stabilization
      info  =.false. ! silent output
!     call schcc(nch,ecm,zv*zc,inc,ql,factor,dr,r0,nr-1,wf,phase,smat)
!      call schcc(nch,ecm,zv*zc,inc,ql,factor,dr,r0,
!     & nr,wf,phase,smat,info)
      call schcc_erwin(nch,ecm,zc*zv,inc,ql,factor,dr,r0,
     & nr,wf,phase,smat,method,info)
c to get continuum wfs 'normalized' as <k|k'> = delta(k-k')
!      wf(:,:)=wf(:,:)*sqrt(2./pi)

!      if ((ecm.gt.0.5).and.(ecm.lt.1)) then 
!      write(500,*)'# WF wf for ecm,ql(:)=',ecm,ql(1:nchan)
!      do ir=1,nr
!         write(500,'(1f8.3,2x,50f12.8)') rvec(ir),
!     &       (real(wf(n,ir)),n=1,nch)
!      enddo
!      write(500,*)'&'
!      endif
      deallocate(ccmat,vcoup)
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c ------------------------------------------------------------------
c   Solve multichannel Schrodinger equation callint R-matrix routines
c   by P. Descouvemont CPC XXX XXX
c
c   NLAG = nb. of Laguerre bases per radial interval
c   NS   = nb. of radial intervals
c   -----------------------------------------------------------------
      subroutine schcc_rmat(nch,ecm,z12,inc,ql,conv,dr,
     & r0,nr,wf,phase,smat,info,nlag,ns)
      use nmrv,only: vcoup,ech
      use xcdcc, only: rvcc
      use constants , only: e2
      use memory
      implicit none
      logical :: iftr,twf,info
c     -------------------------------------------------------------
      integer :: ir,nma,nch,inc,nr,nmax
      integer*4, parameter:: ninc=1,nr2=0
      integer*4 :: nlag,ns,ich,ichp,nopen,nvc(ninc)
c     -----------------------------------------------------------
      real*8 :: ecm,tkch,z12,conv
      real*8  :: r,rmax,dr,r0
      real*8  :: zrma(nlag*ns),eta(nch),kch(nch),ql(nch)
      real*8  :: aux,jci,jcf
      real*8,parameter :: alpha=0.
!      real*8 fc(500),dfc(500),gc(500),dgc(500),xfac(nchan,nchan)
c     -----------------------------------------------------------
      complex*16 :: cfival,caux,phc,ph2
      complex*16 :: cpot(nlag*ns,nch,nch),cu(nch,nch),
     &              cf(nlag*ns,nch,ninc),cpnl,faux(nr)
      complex*16 cfx(nch)
      complex*16, intent(out):: wf(nch,nr)
      complex*16 :: phase(nch),smat(nch) !,svel(nch)
c     ----------------------------------------------------------
!      call cpu_time(start)
c *** Initialize some variables & R-matrix
      nmax=nlag*ns
      rmax=rvcc(nr)
      twf=.false.
      nvc(1)=inc 
      call rmat_ini(nlag,ns,rmax,zrma)
c *** -----------------------------------


c *** Interpolation of coupling matrix in Lagrange mesh
      do ich=1,nch
      aux=ecm+ech(inc)-ech(ich)
      kch(ich)=sqrt(conv*abs(aux))
      eta(ich)=conv*z12*e2/kch(ich)/2.
      if (aux.lt.0) kch(ich)=-kch(ich)

      do ichp=ich,nch
      phc=(0d0,1d0)**(ql(ichp)-ql(ich))
      faux(1:nr)=conv*vcoup(ich,ichp,1:nr)
      do 4 ir=1,nmax
      r=zrma(ir)
      caux=cfival(r,rvcc,faux,nr,alpha)
      cpot(ir,ich,ichp)=caux
      if (ich.ne.ichp) cpot(ir,ichp,ich)=caux
!    &   *conjg(phc)
    4 continue
      enddo !ich
      enddo !ichp
      deallocate(vcoup) ! not needed any more ?????????????

! TEST
!      do ir=1,nmax
!        write(1,'(1f8.3,2x,50f12.8)')zrma(ir),cpot(ir,1,1)
!      enddo
!      do ir=1,nr
!        write(2,'(1f8.3,2x,50f12.8)')rvec(ir),ccmat(1,1,ir)
!      enddo 
      call rmatrix(nch,ql,kch,eta,rmax,nlag,ns,cpot,cu,nmax,nch,
     &             nopen,twf,cf,nmax,nch,ninc,nvc,0,cpnl)

       smat(1:nopen)=cu(1:nopen,inc)*sqrt(kch(inc)/kch(1:nopen))

       do ich=1,nch     
       write(*,500)ich,smat(ich),cu(ich,inc),abs(smat(ich))
500   format(10x,"Chan. #",i3,5x,
     &    "S=(",1f8.5,",",1f8.5,")", 5x, 
     &    "S=(",1f8.5,",",1f8.5,") ->  |S|=",1f9.6)     
       enddo

       if (twf) then
       do  ir=1,nmax
       write(501,1002)zrma(ir),
     &      (cf(ir,ich,1)*(0.,-1.),ich=1,nch)
!     &      (cf(ir,ich,1)*(0.,-1.)*sqrt(kch(ich)/kch(inc)),ich=1,nch)
1002   format(f8.3,10(2es12.4,2x))
       enddo
       write(501,*)'&'
       endif
          
     
      return
      end subroutine


c -------------------------------------------------------
c   Solve multichannel Schrodinger equation
c   ecm     = kinetic energy
c           (E= Ecm + core energy)
c   nch     = number of channels
c   inc     = incoming channel
c   ql(nch) = array with l'values for nch channels
c             (NEEDS TO BE ALLOCATED BY USER)
c   nr      = number of radial points
c   dr      = radial step 
c   z1*z2   = product of charges
c   wf      = array to store wfs(nch,nr) (out)
c   phase   = array to store phase shifts (out)
c   factor  = conversion factor 2*mu/hbar^2
c   
c   The following variables must be passed through the nmrv module:
c    ech(nch) = channel energies (NEEDS TO BE ALLOCATED BY USER)
C    vcoup   = matrix with couplings (nchan,nchan,nr)
c   ----------------------------------------------------
c   ym is the component of the solution at r-h (similarly for Z)
c   y0 is the component of the solution at r
c   yp                                     r+h
c ------------------------------------------------------
      subroutine schcc(nch,ecm,z12,inc,ql,factor,dr,r0,
     & npt,wf,phase,smat,info)
      use nmrv,only: nr,h,vcoup,conv,ech,debug,rmin,hort,
     &               rvec,rmort
      use globals, only: verb
      use factorials
      use constants , only: e2
      use memory
      implicit none
      logical:: copen(nch),orto,info
      integer:: nch,klog,npt,ql(nch)
      integer:: ir,irmin,is,isp,ich,inc,l,lmax
      integer:: i,j,k,im,ip,i0,norto
c     ----------------------------------------------
      real*8 :: ymem,ti,tf
      real*8 :: ecm,tkch,eta,z12
      real*8 :: dr,r,r0,rm,h2,factor
      real*8 :: rturn,rmorto,rmatch
      real*8 :: kch(nch),kch2(nch),rturnv(nch),etav(nch)
      real*8 :: start,finish,aux,big,small
      real*8 , parameter:: pi=acos(-1d0)
c     ----------------------------------------------
      complex*16, intent(out):: wf(nch,nr)
      complex*16, dimension(nch,nch):: yp,yc,y0,ym
      complex*16, dimension(nch,nch):: zp,zc,z0,zm
      complex*16 :: y(nch,nch,npt)
      complex*16 :: phase(nch),smat(nch)
      complex*16 :: a(nch,nch),c1,c2
c     ---------------------------------------------
      y=0d0      
      debug=.false.

      big=huge(big)
      small=epsilon(small)
      h=dr
      h2=h*h
      conv=factor
      copen(:)=.true.
      nr=npt
      klog=99
      rmin=0
      lmax=maxval(ql)
      orto=.false.

      if (debug) then
        write(*,*)'schcc: ql=',ql(1:nch)
        write(*,*)'schcc: ecm=',ecm
      endif

      if (allocated(rvec)) deallocate(rvec)
      allocate(rvec(nr))
      do ir=1,nr
        rvec(ir)=r0+(ir-1)*dr
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius for derivative

      norto=nint(hort/h)
      if (abs(hort).gt.0.) orto=.true.

c     ------------------------------------------------
      if (verb.ge.3) then
        ymem=nch*nch*(npt+2)*lc16/1e6
        write(*,190) ymem
190     format(5x,"[ WFS require", 1f8.2," Mbytes ]")
        ymem=nch*nch*lc16/1e6
!        write(*,*)'Auxiliary arrays require 8x',
!     &  ymem,'=',ymem*8, 'Mbyes'
!        write(*,*)'Max real=',big
      endif 
c     ---------------------------------------------------
      if (ecm.lt.0) then 
      write(*,200) ecm
200   format(4x,'Skipping negative energy:',1f7.3, ' MeV')
      wf(:,:)=0d0
      return
      endif

      if (rmort.gt.0) rmorto=rmort
         
      do ich=1,nch
      aux=ecm+ech(inc)-ech(ich)
      kch2(ich)=conv*aux
      if (aux.gt.0) then
        if (debug) write(99,300) ich,aux,"open"
        kch(ich)=sqrt(conv*aux)
        etav(:)=conv*z12*e2/kch(:)/2.
        rturnv(ich) =(etav(ich)+SQRT(etav(ich)**2 + 
     &               ql(ich)*(ql(ich)+1d0)))/kch(ich)
        if (rturnv(ich).gt.rmorto) rmorto=rturnv(ich) 
        copen(ich)=.true.
      else
        if (debug) write(99,300) ich,aux,"closed"
        kch(ich)=sqrt(-conv*aux)
        rmorto=rm
        copen(ich)=.false.
      endif
300   format(3x,"Channel",i2," Final Ecm=",1f8.3," MeV ->",a6)
      enddo
      if (rmort.lt.0) rmorto=rmorto+abs(rmort) 



c classical turning point for incoming channel
      eta=conv*z12*e2/kch(inc)/2.
      l=ql(inc)
      RTURN =(eta+SQRT(eta**2 + L*(L+1d0)))/kch(inc)  
c turning point for all channels
!      rmorto=maxval(rturnv(:))  ! outermost turning point

      if (debug) write(*,'(6(5x,i3,2x,1f8.3))')
     &   (ich,rturnv(ich),ich=1,nch) 
 

      call factorialgen(lmax+1)
      call cpu_time(start)
c 
c Starting values and initial integration point
c 
!      rmin=h ;       irmin=1
      rmin=max(h,rturn-10.) ! CHECK h or r0 ?
      irmin=nint(rmin/h)+1

!      write(*,*) 'PC-Numerov  from ir=',irmin,' r=',rmin,' fm',
!    & ' to  R=',Rm,'in steps of',h,' and matching at Rmatch=',rmatch


      if (debug) write(95,*)'Starting value ir=,',irmin,'r=',rvec(irmin)
      ym(:,:)=0d0; zm(:,:)=0d0; 
      y0(:,:)=0d0; z0(:,:)=0d0;
      do ich=1,nch
        l=ql(ich)
!        if (debug) write(klog,*)'ich,l,dlfac(l)',ich,l,dlfac2(l)
 
50      aux= (kch(ich)*h)**(l+1)/exp(0.5*dlfac2(l))
!        if (aux.lt.1e-30) then
!          rmin=rmin+5*h
!          goto 50        
!        endif

!        irmin=nint(rmin/h) !+1
!        y0(ich,ich)=(kch(ich)*h)**(l+1)/exp(0.5*dlfac2(l))
!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        y0(ich,ich)=1e-10
        y0(ich,ich)=1e-10
!!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
        if (debug) then
        write(95,320) (y0(ich,is),is=1,min(10,nch)) 
320     format(5x,50("(",2g16.6,")",5x))
        endif
      enddo

!      if (debug) write(95,*)'Initial dot product(ABS):'
!      write(95,'(10(2x,"is=",i2,3x,1g10.5))') (is,
!     &     abs(dot_product(y0(:,is),y0(:,is-1)))
!     & , is=2,nch)


      if (info) then
      write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
     &  5x,"outer=",1f8.2," fm",
     &  5x,"R-min=",1f6.2," fm")') 
     &  minval(rturnv(:)), maxval(rturnv(:)), rmin      
      endif
c
c Zm=F[r=0,ym] and Z0=F[r=h,y0]
c
      if (irmin.lt.2) write(*,*)'wrong irmin=',irmin
      call zfun(nch,ql,kch2,irmin-1,ym,zm)
      call zfun(nch,ql,kch2,irmin  ,y0,z0)
      
      do ir=irmin,nr-1  
      r=rvec(ir) ! (ir-1)*h    
      im=ir-1  ! r-h
      i0=ir    ! r
      ip=ir+1  ! r+h
      call zfun(nch,ql,kch2,im,ym,zm)
!      call zfun(nch,ql,kch2,i0,y0,z0)

      if (debug.and.(mod(ir,10).eq.0)) then
         write(95,'(/,"y0(x) ir=",i4,3x,"r=",1f7.2)') ir,r
         do ich=1,min(nch,10)
      write(95,'(5x,50g14.5)')(y0(ich,is),is=1,min(nch,10))
 !        write(*,320) (y0(ich,is),is=1,min(nch,10))
         enddo
      endif !debug

c
c predictor
c                         
      YP(:,:)=2*Y0(:,:)-YM(:,:)+h2*Z0(:,:)  
 
      if (maxval(abs(yp)).gt.big) then
         write(*,*)' ** WARNING ** Solution too big at ir,r=',ir,r
      endif
      if (debug.and.(mod(ir,10).eq.0)) then
        write(95,*)'Predictor for ir,r=',ir,r
        do ich=1,min(nch,10)
!        write(95,'(5x,50g14.5)') (yp(ich,is),is=1,nch)
        write(95,320) (yp(ich,is),is=1,min(nch,10))
        enddo
        call flush(95)
      endif !debug
                         
c
c corrector
c         
      call zfun(nch,ql,kch2,ip,yp,zp)
      yc(:,:)=yp(:,:)+(h2/12.)*(zp(:,:)+zm(:,:)-2*z0(:,:))

      if (debug.and.(mod(ir,10).eq.0)) then
        write(95,*)'Corrector:'
        do ich=1,min(nch,10)
        write(95,'(5x,50g14.5)') (yc(ich,is),is=1,nch)
        enddo
      endif ! debug

c 
c modify final solution 
c
      call zfun(nch,ql,kch2,ip,yc,zc)
      yp(:,:)=yc(:,:)+(h2/30.)*(zc(:,:)-zp(:,:))
      call zfun(nch,ql,kch2,ip,yp,zp)

      if (debug.and.(mod(ir,10).eq.0)) then
        write(95,*)'Modified corrector:'
        do ich=1,min(nch,10)
        write(95,'(5x,50g14.5)') (yp(ich,is),is=1,nch)
        enddo
      
        write(95,*)'Max|Yc-Yp|:',maxval(abs(yc-yp))
!        do ich=1,nch
!        write(94,'(5x,50g14.5)') 
!     &  (abs(yc(ich,is)-yp(ich,is)),is=1,nch)
!        enddo
        call flush(95)

      endif ! debug

c store solution
      y(1:nch,1:nch,ip)=yp(1:nch,1:nch) 


c Re-orthogonalize solution vectors
!      orto =.true.
!      if (orto) write(*,*)ir,ir-irmin,norto,r,rmorto,mod(ir-irmin,norto)
      if (orto.and.(r.lt.rmorto).and.(ir.gt.irmin)
     & .and.(mod(ir-irmin,norto).eq.0)) then 
!      if (orto.and.(r.le.rmorto).and.(mod(ir-irmin,norto).eq.0)) then 
        write(*,'(5x,"-> orthogonalizing at r=",1f7.2)')r

        if (debug) then
        write(94,*)'|yp x yp| (before G-S) for ir,r',ir,r
        do is=2,nch
        write(94,'(5x,50g14.5)') 
     &    ( abs(dot_product(yp(:,is),yp(:,ich)))
     &      / abs(dot_product(yp(:,is),yp(:,is)))
     &      / abs(dot_product(yp(:,ich),yp(:,ich))), ich=1,nch) 
!         write(94,'(5x,i2,"x",i2,3x,50g14.5)') is,is-1,
!     &    abs(dot_product(yp(:,is),yp(:,is-1)))
        enddo
        endif !debug

  
!      call gsorto(yp,a,nch)  ! standard Gram-Schmidt
!      call mgsorto(yp,a,nch)  ! modified Gram-Schmidt
!      do i=irmin,ir+1
!      y(:,:,i)=matmul(y(:,:,i),a(:,:))
!      enddo ! i

!       write(*,*)'calling orto for ir,r=',ir,r
       call cpu_time(ti)
       call qrfac(yp,y,nch,ip,irmin,nr)      ! QR factorization
       call cpu_time(tf)
       torto=torto+tf-ti

      if (debug) then
      write(94,*)'Ortog. matrix at ir=',ir
      do ich=1,nch
        write(94,'(5x,50g14.5)') 
     &  (a(ich,is),is=1,nch)
           write(94,'(5x,50g14.5)') 
     &  (a(ich,is),is=1,nch)
      enddo
      call flush(94)
      endif 

!      y(:,:,ip)=yp

      if (debug) then
      write(90,*)'y at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(90,'(5x,i3,50g14.5)') ich,
     &  (abs(y(ich,is,ip)), is=1,min(15,nch))
      enddo
      write(90,*)' '

      write(90,*)'yp at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(90,'(5x,i3,50g14.5)') ich,
     &  (abs(yp(ich,is)), is=1,min(15,nch))
      enddo
      write(*,*)' '

      write(90,*)'|y-yp|/y at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(90,'(5x,i3,50g14.5)') ich,
     &  (abs(y(ich,is,ip)-yp(ich,is))
     &   /abs(y(ich,is,ip)), is=1,min(15,nch))
      enddo
      write(*,*)' '
      endif

      yp=y(:,:,ip);  call zfun(nch,ql,kch2,ip,yp,zp)
      y0=y(:,:,i0);  call zfun(nch,ql,kch2,ip,y0,z0)


c dot products after orthogonalization
      if (debug) then
        write(94,*)'|yp x yp| (after G-S)'
        do is=1,nch
!        write(95,'(5x,i2,"x",i2,3x,50g14.5)') 
        write(94,'(5x,50g14.5)') 
     &    ( abs(dot_product(yp(:,is),yp(:,ich)))
     &      / abs(dot_product(yp(:,is),yp(:,is)))
     &      / abs(dot_product(yp(:,ich),yp(:,ich))), ich=1,nch) 

        enddo
       endif !debug

!       y(1:nch,1:nch,ir+1)=yp(1:nch,1:nch)  

      if (debug) then
        write(94,*)'Gram-Schmidt ortogonalized:'
        do ich=1,nch
        write(94,'(5x,50g14.5)') (yp(ich,is),is=1,nch)
        enddo
        write(94,*)' '
        call flush(94)
      endif !debug


      endif !orto
! ------------------------------------------------------------


c Alternative Gramm-Schmitt orthogonalization of solution vectors
!      orto=.false.
!      rmorto=maxval(rturnv(:))
!      if (orto.and.(r.lt.rmorto).and.(mod(ir,5).eq.0)) then 
!         call gsorto(yp,a,nch)
!      end if 

!!!!!!!!!!!!!!!!!!!!!!!!! TEST
      if (debug) then
!       write(91,*) r,maxval(abs(y(:,:,ip)))
       write(91,'(5x,1f6.3,3x,50g14.5)') r,(abs(y(1,ich,ip)),ich=1,nch)
!       y(:,:,:)=y(:,:,:)/1.e5
!       yp=y(:,:,ip)
!       y0=y(:,:,i0)
!       call zfun(nch,ql,kch2,ip,yp,zp)
!       call zfun(nch,ql,kch2,ip,y0,z0)
      endif       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c
c ym <- y0 <- yp 
c
!      ym(:,:)=  y(1:nch,1:nch,ir)
!      zm(:,:)
      ym(:,:)=y0(:,:) 
      zm(:,:)=z0(:,:)
      y0(:,:)=yp(:,:) 
      z0(:,:)=zp(:,:)


      enddo !ir
      call cpu_time(finish)

!!! TEST
!      if (maxval(abs(y(:,:,nr))).gt.1e20) then
!        y(:,:,:)=y(:,:,:)/1.e20
!      endif
!!! TEST

      
c 
c Match with asymptotic solution (derivative)
c Using derivative and gauss5
       call cpu_time(start)
!      call matching(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! gauss5
c Adjacent points
      call matching3(ecm,z12,nch,ql,lmax,inc,nr,y,wf,phase,smat,info) ! gauss5

!      call matching2(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! LAPACK 

       call cpu_time(finish)
       tmatch=tmatch+finish-start
      if (debug) write(klog,'(2f10.3,3x,10g14.6)') ecm,z12,ql(1),
     &  (phase(ich),ich=1,nch) 


      call flush(6)
      end subroutine schcc



c -------------------------------------------------------
c   Solve multichannel Schrodinger equation
c   ecm     = kinetic energy
c           (E= Ecm + core energy)
c   nch     = number of channels
c   inc     = incoming channel
c   ql(nch) = array with l'values for nch channels
c             (NEEDS TO BE ALLOCATED BY USER)
c   nr      = number of radial points
c   dr      = radial step 
c   z1*z2   = product of charges
c   wf      = array to store wfs(nch,nr) (out)
c   phase   = array to store phase shifts (out)
c   factor  = conversion factor 2*mu/hbar^2
c   
c   The following variables must be passed through the nmrv module:
c    ech(nch) = channel energies (NEEDS TO BE ALLOCATED BY USER)
C    vcoup   = matrix with couplings (nchan,nchan,nr)
c   ----------------------------------------------------
c   ym is the component of the solution at r-h (similarly for Z)
c   y0 is the component of the solution at r
c   yp                                     r+h
c ------------------------------------------------------
      subroutine schcc_ena(nch,ecm,z12,inc,ql,factor,dr,r0,
     & npt,wf,phase,smat,method,info)
      use nmrv,only: nr,h,vcoup,conv,ech,debug,rmin,hort,rvec,rmort
      use globals, only: verb
      use factorials
      use constants , only: e2
      use memory
      implicit none
      integer:: method
      logical:: copen(nch),orto,raynal,info
      integer:: nch,klog,npt,ql(nch)
      integer:: ir,irmin,is,isp,ich,inc,l,lmax
      integer:: i,j,k,im,ip,i0,norto

      real*8 :: ecm,tkch,eta,z12
      real*8 :: dr,r,r0,rm,h2,factor
      real*8 :: rturn,rmorto,rmatch
      real*8 :: kch(nch),kch2(nch),rturnv(nch),etav(nch)
      real*8 :: start,finish,aux

      complex*16, intent(out):: wf(nch,nr)
      complex*16, dimension(nch,nch):: yp,yc,y0,ym
      complex*16, dimension(nch,nch):: zp,zc,z0,zm
      complex*16, dimension(nch,nch):: wp,w0,wm

      complex*16 :: y(nch,nch,npt),w(nch,nch,npt)
      complex*16 :: phase(nch),smat(nch)
      real*8 , parameter:: pi=acos(-1d0)
      complex*16 :: a(nch,nch),c1,c2
      real*8     :: ymem,ti,tf
c     -------------------------------------------------------
      y=0d0      
      debug=.false.
      if (method.eq.3) raynal=.true.
      h=dr
      h2=h*h
      conv=factor
      copen(:)=.true.
      nr=npt
      klog=99
      rmin=0
      lmax=maxval(ql)
      orto=.false.

      
      if (allocated(rvec)) deallocate(rvec)
      allocate(rvec(nr))
      do ir=1,nr
        rvec(ir)=r0+(ir-1)*dr
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius

      norto=nint(hort/h)
      if (abs(hort).gt.0.) orto=.true.

      if (verb.ge.3) then
        ymem=nch*nch*(npt+2)*lc16/1e6
        write(*,190) ymem
190     format(5x,"[ WFS require", 1f8.2," Mbytes ]")
        ymem=nch*nch*lc16/1e6
!        write(*,*)'Auxiliary arrays require 8x',
!     &  ymem,'=',ymem*8, 'Mbyes'
      endif 

      if (ecm.lt.0) then 
      write(*,200) ecm
200   format(4x,'Skipping negative energy:',1f7.3, ' MeV')
      wf(:,:)=0d0
      return
      endif

      rmorto=rmort
      do ich=1,nch
      aux=ecm+ech(inc)-ech(ich)
      kch2(ich)=conv*aux
      if (aux.gt.0) then
        if (debug) write(99,300) ich,aux,"open"
        kch(ich)=sqrt(conv*aux)
        etav(ich)=conv*z12*e2/kch(ich)/2.
        rturnv(ich) =(etav(ich)+SQRT(etav(ich)**2 + 
     &              ql(ich)*(ql(ich)+1d0)))/kch(ich) 
        if(rturnv(ich).gt.rmorto) rmorto=rturnv(ich) 
        copen(ich)=.true.
      else
        if (debug) write(99,300) ich,aux,"closed"
        kch(ich)=sqrt(-conv*aux)
        rmorto=rm
        copen(ich)=.false.
      endif
300   format(3x,"Channel",i2," Final Ecm=",1f8.3," MeV ->",a6)

      enddo


c classical turning point for incoming channel
      eta=conv*z12*e2/kch(inc)/2.
      l=ql(inc)
      RTURN =(eta+SQRT(eta**2 + L*(L+1d0)))/kch(inc)  

c turning point for each channel

      if (debug) write(*,'(6(5x,i3,2x,1f8.3))')
     &   (ich,rturnv(ich),ich=1,nch) 
 

      call factorialgen(lmax+1)
      call cpu_time(start)
c 
c Starting values and initial integration point
c 
!      rmin=h ;       irmin=1
      rmin=max(h,rturn-10.) ! CHECK h or r0 ?
      irmin=nint(rmin/h)+1

!      write(*,*) 'Enhanced Numerov from ir=',irmin,' r=',rmin,' fm',
!     & ' to  R=',Rm,'in steps of',h,' and matching at Rmatch=',rmatch


      if (debug) write(95,*)'Starting value ir=,',irmin,'r=',rvec(irmin)
      ym(:,:)=0d0; zm(:,:)=0d0; 
      y0(:,:)=0d0; z0(:,:)=0d0;
      do ich=1,nch
        l=ql(ich)
!        if (debug) write(klog,*)'ich,l,dlfac(l)',ich,l,dlfac2(l)
 
!50      aux= (kch(ich)*h)**(l+1)/exp(0.5*dlfac2(l))
!        if (aux.lt.1e-30) then
!          rmin=rmin+5*h
!          goto 50        
!        endif

!        irmin=nint(rmin/h) !+1
!        y0(ich,ich)=(kch(ich)*h)**(l+1)/exp(0.5*dlfac2(l))
!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        y0(ich,ich)=1e-10
!!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
        if (debug) then
!        write(95,'("Starting values:")')
        write(95,320) (y0(ich,is),is=1,min(10,nch)) 
320     format(5x,50("(",2g16.6,")",5x))
        endif
      enddo

!      if (debug) write(95,*)'Initial dot product(ABS):'
!      write(95,'(10(2x,"is=",i2,3x,1g10.5))') (is,
!     &     abs(dot_product(y0(:,is),y0(:,is-1)))
!     & , is=2,nch)

      write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
     &  5x,"outer=",1f8.2," fm",
     &  5x,"R-min=",1f6.2," fm")') 
     &  minval(rturnv(:)), maxval(rturnv(:)), rmin      

c 
c Initial values: w(rmin) and w(rmin+h) 
c 
      if (irmin.lt.2) write(*,*)'wrong irmin=',irmin
      call zfun(nch,ql,kch2,irmin-1,ym,zm)
      wm=ym-(h2/12.)*zm
      call zfun(nch,ql,kch2,irmin  ,y0,z0)
      w0=y0-(h2/12.)*z0
 
      do ir=irmin,nr-1  
      r=rvec(ir) ! (ir-1)*h    
      im=ir-1  ! r-h
      i0=ir    ! r
      ip=ir+1  ! r+h
      call csh(nch,ql,kch2,i0,w0,z0,method) ! Cosh[Sqrt(sqrt(h2 F(x))]

      if (debug.and.(mod(ir,10).eq.0)) then
         write(95,'(/,"y0(x) ir=",i4,3x,"r=",1f7.2)') ir,r
         do ich=1,min(nch,10)
      write(95,'(5x,50g14.5)')(y0(ich,is),is=1,min(nch,10))
         enddo
      endif !debug

      WP(:,:)=2*Z0(:,:)-WM(:,:)  

      if (debug.and.(mod(ir,10).eq.0)) then
        write(95,*)'wp at ir,r=',ir,r
        do ich=1,min(nch,10)
        write(95,320) (wp(ich,is),is=1,min(nch,10))
        enddo
        call flush(95)
      endif !debug
                         

c calculate and store Y=(1-T)W  
      if ((ir.ge.nr-8).or.(orto)) then      
        call w2y(wp,yp,ql,kch2,ir,nch)
        y(1:nch,1:nch,ip)=yp(1:nch,1:nch) 
      endif

c Re-orthogonalize solution vectors
      if (orto.and.(r.lt.rmorto).and.(ir.gt.irmin)
     & .and.(mod(ir-irmin,norto).eq.0)) then 
        write(*,'(5x,"-> orthogonalizing at r=",1f7.2)')r
        if (debug) then
        write(94,*)'|yp x yp| (before G-S) for ir,r',ir,r
        do is=2,nch
        write(94,'(5x,50g14.5)') 
     &    ( abs(dot_product(yp(:,is),yp(:,ich)))
     &      / abs(dot_product(yp(:,is),yp(:,is)))
     &      / abs(dot_product(yp(:,ich),yp(:,ich))), ich=1,nch) 
        enddo
        endif !debug

  
!      call gsorto(yp,a,nch)  ! standard Gram-Schmidt
!      call mgsorto(yp,a,nch)  ! modified Gram-Schmidt
!      do i=irmin,ir+1
!      y(:,:,i)=matmul(y(:,:,i),a(:,:))
!      enddo ! i

       call cpu_time(ti)
       call qrfac(yp,y,nch,ip,irmin,nr)      ! QR factorization
!       call qrfac(wp,w,nch,ip,irmin,nr)      ! QR factorization



       call cpu_time(tf)
       torto=torto+tf-ti

      if (debug) then
      write(94,*)'Ortog. matrix at ir=',ir
      do ich=1,nch
        write(94,'(5x,50g14.5)') 
     &  (a(ich,is),is=1,nch)
           write(94,'(5x,50g14.5)') 
     &  (a(ich,is),is=1,nch)
      enddo
      call flush(94)
      endif 

      if (debug) then
      write(90,*)'y at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(90,'(5x,i3,50g14.5)') ich,
     &  (abs(y(ich,is,ip)), is=1,min(15,nch))
      enddo
      write(90,*)' '

      write(90,*)'yp at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(90,'(5x,i3,50g14.5)') ich,
     &  (abs(yp(ich,is)), is=1,min(15,nch))
      enddo
      write(*,*)' '

      write(90,*)'|y-yp|/y at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(90,'(5x,i3,50g14.5)') ich,
     &  (abs(y(ich,is,ip)-yp(ich,is))
     &   /abs(y(ich,is,ip)), is=1,min(15,nch))
      enddo
      write(*,*)' '
      endif

      yp=y(:,:,ip);  !call zfun(nch,ql,kch2,ip,yp,zp)
      y0=y(:,:,i0);  !call zfun(nch,ql,kch2,ip,y0,z0)
      call zfun(nch,ql,kch2,ip,yp,zp)
      wp=yp-(h2/12.)*zp
      call zfun(nch,ql,kch2,i0,y0,z0)
      w0=y0-(h2/12.)*z0


!      wp=w(:,:,ip);  !call zfun(nch,ql,kch2,ip,yp,zp)
!      w0=w(:,:,i0);  !call zfun(nch,ql,kch2,ip,y0,z0)


c dot products after orthogonalization
      if (debug) then
        write(94,*)'|yp x yp| (after G-S)'
        do is=1,nch
!        write(95,'(5x,i2,"x",i2,3x,50g14.5)') 
        write(94,'(5x,50g14.5)') 
     &    ( abs(dot_product(yp(:,is),yp(:,ich)))
     &      / abs(dot_product(yp(:,is),yp(:,is)))
     &      / abs(dot_product(yp(:,ich),yp(:,ich))), ich=1,nch) 

        enddo
       endif !debug

!       y(1:nch,1:nch,ir+1)=yp(1:nch,1:nch)  

!      if (debug) then
!        write(94,*)'Gram-Schmidt ortogonalized:'
!        do ich=1,nch
!        write(94,'(5x,50g14.5)') (yp(ich,is),is=1,nch)
!        enddo
!        write(94,*)' '
!        call flush(94)
!      endif !debug

      endif !orto
! ------------------------------------------------------------


c Alternative Gramm-Schmitt orthogonalization of solution vectors
!      orto=.false.
!      rmorto=maxval(rturnv(:))
!      if (orto.and.(r.lt.rmorto).and.(mod(ir,5).eq.0)) then 
!         call gsorto(yp,a,nch)
!      end if 

!!!!!!!!!!!!!!!!!!!!!!!!! TEST
      if (debug) then
       write(91,'(5x,1f6.3,3x,50g14.5)') r,(abs(y(1,ich,ip)),ich=1,nch)
      endif       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c
c wm <- w0 <- wp 
c
      wm(:,:)=w0(:,:) 
      w0(:,:)=wp(:,:) 

      enddo !ir
      call cpu_time(finish)


c 
c Match y with asymptotic solution -> S-matrix
c            
      call cpu_time(start)
!      call matching(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! gauss5 + derivative
!      call matching2(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! LAPACK 
      call matching3(ecm,z12,nch,ql,lmax,inc,nr,y,wf,phase,smat,info) ! gauss5 + 2 points

       call cpu_time(finish)
       tmatch=tmatch+finish-start

      if (debug) write(klog,'(2f10.3,3x,10g14.6)') ecm,z12,ql(1),
     &  (phase(ich),ich=1,nch) 

      call flush(6)
      end subroutine schcc_ena


c ... Calculate y from w(r)=[1-T(r)] y(r)
c     with T(r)= (h^2/12) F(x)
c     and y''= F(x) y(x)     
c 
c     Uses Lapack to solve the equation
c     A X =B
c     with B(LDB,NRHS) = W(r)
c          A(LDA,LDA)  = 1-T(r)
c          X(LDA,LDA)  = Y(r)
      subroutine w2y(w,y,ql,kch2,ir,n)
      use nmrv,only   : conv,vcoup,h,debug,rvec
      implicit none
      real*8    :: r,kch2(n),h2
      complex*16             :: w(n,n)
      complex*16             :: y(n,n)
      complex*16             :: a(n,n),ainv(n,n)
!      real*8                 :: c(n,n),cinv(n,n)
      integer   :: ich,l,ir,n,is
      integer   :: ql(n)
      real,parameter:: eps=1e-6
      integer   :: nrhs, ldb, ipiv(n), info    
      character*1  trans
      EXTERNAL         ZGETRF, ZGETRS
      debug=.false.
      a=0
      h2=h*h
      r=rvec(ir) 
      if (r.lt.eps) r=eps
c ... For zgetrs
      trans='N'
      ldb=n
      nrhs=n
c ... 
      

c A= I - T(x) = 1 -(h**2/12) F(x)
c 
      do ich=1,n ! diagonal part of A
      a(ich,ich)=1.-(h2/12.)*(l*(l+1)/r**2 - kch2(ich))  
      enddo
      a(:,:)=a(:,:)-(h2/12.)*conv*vcoup(:,:,ir)
       

!!!! TEST lapack routines
!      a=0
!      a(1,1)=4.; a(1,2)=3.
!      a(2,1)=6.; a(2,2)=3.     
!      w(1,1)=1.; w(1,2)=0.
!      w(2,1)=0.; w(2,2)=1.

!      c(1,1)=4; c(1,2)=3
!      c(2,1)=6; c(2,2)=3
      
!!! 

!      Cinv = C
!      write(*,*)'call dgetrf with n=',n
!      call DGETRF(n, n, Cinv, n, ipiv, info)
!      if (info /= 0) then
!       write(*,*)'dgetrf returned error',info; stop
!      else
!       write(*,*)'dgetrf ok: cinv=',cinv
!       write(*,*)'ipiv=',ipiv
!      endif
      
      if (debug) then
      write(*,*)'  W at r=',r
      do ich=1,n
      write(*,'(5x,50g16.5)') (W(ich,is),is=1,n)
      enddo

      write(*,*)'  A at r=',r
      do ich=1,n
      write(*,'(5x,50g16.5)') (A(ich,is),is=1,n)
      enddo
      endif !debug

  ! Compute an LU factorization of matrix A
  ! using partial pivoting with row interchanges.
  ! Store A in Y to prevent it from being overwritten by LAPACK
!      Y = A
!      n = size(A,1)
!      write(*,*)'call zgetrf with n=',n
      call ZGETRF(n, n, a, n, ipiv, info)
      if (info /= 0) then
       write(*,*)'zgetrf returned error',info; stop
!Matrix is numerically singular!'
      else
!       write(*,*)'zgetrf info=',info,' ipiv=',ipiv
!       write(*,*)'Ainv=',Ainv
      endif


!      write(*,*)'trans,n,nrhs,ldb=',trans,n,nrhs,ldb
      y(:,:)=w(:,:)
      call zgetrs(trans,n,nrhs,a,n,ipiv,y,ldb,info)

      if (info /= 0) then
          stop 'Matrix inversion failed!'
      else
!       write(*,*)'cgetrs info=',info
!       write(*,*)'ipiv=',ipiv
      endif


      if (debug) then
      write(*,*)'  Y at r=',r
      do ich=1,n
      write(*,'(5x,50g16.5)') (Y(ich,is),is=1,n)
      enddo
      write(*,*)' '
      endif !debug

      
      end subroutine 



c ----------------------------------------
c  Calculates function Cosh(Sqrt[h^2/12 F(x)]) 
c  and performs multiplication
c  Cosh[] Wm(:,:) = ZM(:,:)
c ----------------------------------------
      subroutine csh(nch,ql,kch2,ir,wm,zm,method)
      use nmrv,only   : conv,vcoup,h,debug,rvec
      use memory, only: tcsh
      implicit none
      integer:: method
      integer:: l,ir,ich,is,nch,ql(nch)
      real*8 :: r,kch2(nch),h2
      real*8 :: ti,tf
      real,parameter:: eps=1e-6
      complex*16,intent(in)  :: wm(nch,nch)
      complex*16,intent(out) :: zm(nch,nch)
      complex*16             :: t(nch,nch) !,id(nch,nch)
      complex*16,dimension(nch,nch):: s,s2,s3,s4
      complex*16             :: saux, waux
!      raynal=.true.
      h2=h*h
      S(:,:)=0
!      id(:,:)=0
      debug=.false.
      call cpu_time(ti)
!      write(*,*)'zfun: nch,ir',nch,ir
!      write(*,*)'zfun: h=',h
!      write(*,*)'zfun: kch2=',kch2(1:nch)
!      write(*,*)'zfun: conv=',conv
      r=rvec(ir) 
      if (r.lt.eps) r=eps

c kinetic energy part
      do ich=1,nch
      l=ql(ich)
      s(ich,ich)= (l*(l+1)/r**2 - kch2(ich))*h2
      enddo
      
c diagonal part of potential matrix
      s(:,:)=s(:,:) + h2*conv*vcoup(1:nch,1:nch,ir)

      select case(method)
      case(1)   ! Thorlacius with 5 terms in Cos[]
       s2=matmul(s,s)
       s3=matmul(s,s2)
       s4=matmul(s,s3)        
       zm=wm +      0.5*matmul(s,wm)
     &      + (1./24.)*matmul(s2,wm)
     &      + (1./720.)*matmul(s3,wm)
     &      + (1./40320.)*matmul(s4,wm)

      case(2)  ! Thorlacius with 5 terms in Cos[], only diag V()
       zm(:,:)=wm(:,:)  + 0.5*matmul(s,wm)   
       do is=1,nch
       do ich=1,nch 
       saux=s(ich,ich)
       waux=wm(ich,is)
       zm(ich,is)= zm(ich,is)  
     &      + (1./24.)*saux**2*waux
     &      + (1./720.)*saux**3*waux
     &      + (1./40320.)*saux**4*waux
       enddo
       enddo

      case(3)   !Raynal (eg. Thorlacius eq. (11))
       t(:,:)=s(:,:)/12.
       zm=wm + matmul(6.*t + 6.*matmul(t,t),wm) 
      end select
 

      call cpu_time(tf)
      tcsh=tcsh+(tf-ti)
      return

      if (debug.and.((ir.eq.1).or.mod(ir,5).eq.0)) then
      write(99,*)'r=',r

      write(99,*)'  V(r):'
      do ich=1,nch
      write(99,'(5x,50f12.5)') (vcoup(ich,is,ir),is=1,nch)
      enddo

      write(99,*)'  y(r):'
      do ich=1,nch
      write(99,'(5x,50g14.5)') (wm(ich,is),is=1,nch)
      enddo

      write(99,*)'  S(r):'
      do ich=1,nch
      write(99,'(5x,50g14.5)') (s(ich,is),is=1,nch)
      enddo

      write(99,*)'  Z(r):'
      do ich=1,nch
      write(99,'(5x,50g14.5)') (zm(ich,is),is=1,nch)
      enddo
      endif ! debug

      end 




c ----------------------------------------
c  Calculates function y''(ri)=F[ri,y(ri)]=G(ri,yi)*y(ri)
c  and stores result in zm(:,:) [= y''(ri)]
c ----------------------------------------
      subroutine zfun_clean(nch,ql,kch2,ir,ym,zm)
      use nmrv,only   : conv,vcoup,h,debug,rvec
      use memory, only: tzfun
      implicit none
      integer:: l,ir,ich,is,nch,ql(nch)
      real*8 :: r,kch2(nch)
      real*8 :: ti,tf
      real,parameter:: eps=1e-6
      complex*16,intent(in)  :: ym(nch,nch)
      complex*16,intent(out) :: zm(nch,nch)
      complex*16             :: gmat(nch,nch)
      debug=.false.
      gmat(:,:)=0d0
      call cpu_time(ti)
!      r=h*ir
      r=rvec(ir)

      if (r.lt.eps) r=eps
      do ich=1,nch
      l=ql(ich)
      gmat(ich,ich)= kch2(ich)-l*(l+1)/r**2  
      enddo
      
      gmat(1:nch,1:nch)=gmat(1:nch,1:nch)
     &                 -conv*vcoup(1:nch,1:nch,ir)
      zm=-matmul(gmat,ym)

      call cpu_time(tf)
      tzfun=tzfun+(tf-ti)
      return
      end



c ----------------------------------------
c  Calculates function F[ri,y(ri)]  [= y''(ri)]
c  and stores result in zm(:,:)
c ----------------------------------------
      subroutine zfun(nch,ql,kch2,ir,ym,zm)
      use nmrv,only   : conv,vcoup,h,debug,rvec
      use memory, only: tzfun
      implicit none
      integer:: l,ir,ich,is,nch,ql(nch)
      real*8 :: r,kch2(nch)
      real*8 :: ti,tf
      real,parameter:: eps=1e-6
      complex*16,intent(in)  :: ym(nch,nch)
      complex *16,intent(out) :: zm(nch,nch)
      complex*16             :: gmat(nch,nch)

      gmat(:,:)=0d0
      debug=.false.
      call cpu_time(ti)
!      write(*,*)'zfun: nch,ir',nch,ir
!      write(*,*)'zfun: h=',h
!      write(*,*)'zfun: kch2=',kch2(1:nch)
!      write(*,*)'zfun: conv=',conv
      r=rvec(ir) 
      if (r.lt.eps) r=eps

      do ich=1,nch
      l=ql(ich)
      gmat(ich,ich)= kch2(ich)-l*(l+1)/r**2  
!      write(*,*) 'zfun:r,ich,zdiag=',r,ich,gmat(ich,ich),
!     &  l*(l+1)/r**2  
      enddo
      
      gmat(1:nch,1:nch)=gmat(1:nch,1:nch)
     &                 -conv*vcoup(1:nch,1:nch,ir)
      zm=-matmul(gmat,ym)

      if (debug.and.((ir.eq.1).or.mod(ir,5).eq.0)) then
      write(99,*)'r=',r

      write(99,*)'  V(r):'
      do ich=1,nch
      write(99,'(5x,50f12.5)') (vcoup(ich,is,ir),is=1,nch)
      enddo

      write(99,*)'  y(r):'
      do ich=1,nch
      write(99,'(5x,50g14.5)') (ym(ich,is),is=1,nch)
      enddo

      write(99,*)'  G(r):'
      do ich=1,nch
      write(99,'(5x,50g14.5)') (gmat(ich,is),is=1,nch)
      enddo

      write(99,*)'  Z(r):'
      do ich=1,nch
      write(99,'(5x,50g14.5)') (zm(ich,is),is=1,nch)
      enddo
      endif
      call cpu_time(tf)
      tzfun=tzfun+(tf-ti)
      return
      end



c ...
c ... Enhanced Numerov (version of IJ Thompson used in Fresco) 
c     
c     (simplified version for real couplings, e.g. projectile wfs)
c ... 
      subroutine schcc_erwinrc(nch,ecm,z12,inc,ql,factor,dr,r0,
     & npt,wf,method,info)
      use nmrv,only: nr,h,vcoup,conv,ech,debug,rmin,hort,rvec,cutr,nopen
      use globals, only: verb
      use factorials
      use constants , only: e2
      use memory
      implicit none
      integer:: method
      logical:: copen(nch),orto,raynal,info
      integer:: nch,klog,npt,ql(nch)
      integer:: ir,irmin,is,isp,ich,inc,l,lmax
      integer:: i,j,k,im,ip,i0,norto
c     .........................................................
      real*8     :: ymem,ti,tf
      real*8 , parameter:: pi=acos(-1d0)
      real*8 :: ecm,tkch,eta,z12,r12,ena2,ena3
      real*8 :: dr,r,r0,rm,h2,factor
      real*8 :: rturn,rmorto,rmatch
      real*8 :: kch(nch),kch2(nch),rturnv(nch),etav(nch)
      real*8 :: start,finish,aux,big,small,tmin,tmax
      real*8, dimension(nch,nch):: y0
      real*8, dimension(nch,nch):: z0,zm
      real*8:: y(nch,nch,npt),coupl(nch,nch),v(nch,nch)
c ... for QR factorization
      real*8:: w0(nch,nch),wm(nch,nch)
c ... TEST (delete after debugging)
      real*8:: yaux(nch,nch)
c     .........................................................
      complex*16, intent(out):: wf(nch,nch,nr)
      complex*16 :: phase(nch),smat(nch)
      real*8 :: s(nch),ziv,zpv
      complex*16 :: a(nch,nch)
      complex*16 :: c1,c2,c
      complex*16,parameter:: onec=(1D0,0D0),zero=(0d0,0d0)

c     ---------------------------------------------------------
      R12 = 1D0/12D0
      ENA2 = 2D0/5D0 * R12**2
      ENA3 = - 4D0/35D0 * R12**3
!      ONEC = (1D0,0D0)
      TMAX = 20.
      TMIN = -125.
      big=huge(big)
      small=epsilon(small)
     
      h=dr
      h2=h*h
      conv=factor
      copen(:)=.true.
      nr=npt
      klog=99
      rmin=0
      lmax=maxval(ql)
      orto=.false.
c ... Initialize some variables 
      y=0d0      
      coupl=0
c ... Set TRUE for debugging
      debug=.false.
c ..............................................................................

!      do ich=1,nch
!       write(*,*)'erwinRC: Ecm,ich,inc=',Ecm,ich,inc
!      enddo

!      if (nch.lt.1) then 
!        write(*,*)'Scatcc_erwin: nch=',nch; stop
!      endif

      if (ecm.lt.1e-6) then
        write(*,*)'Scatcc_erwin: Ecm too small! =',ecm
        stop
      endif 
 
      if (allocated(rvec)) deallocate(rvec)
      allocate(rvec(nr))
      do ir=1,nr
        rvec(ir)=r0+(ir-1)*dr
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius if derivative is used

      norto=nint(hort/h)
      if (abs(hort).gt.0.) orto=.true.

      if (info) then
        ymem=nch*nch*(npt+2)*lc16/1e6
        if (verb.ge.3) write(*,190) ymem
190     format(5x,"[ WFS require", 1f8.2," Mbytes ]")
        ymem=nch*nch*lc16/1e6
!        write(*,*)'Auxiliary arrays require 8x',
!     &  ymem,'=',ymem*8, 'Mbyes'
!        write(*,*)' Biggest real=',big
      endif 

      if (ecm.lt.0) then 
      write(*,200) ecm
200   format(4x,'Skipping negative energy:',1f7.3, ' MeV')
      wf(:,:,:)=0d0
      return
      endif

      rmorto=0
      nopen=0
      do ich=1,nch
      aux=ecm+ech(inc)-ech(ich)
      kch2(ich)=conv*aux
      if (aux.gt.0) then
        copen(ich)=.true.
        nopen=nopen+1
        if (debug) write(99,300) ich,aux,"open"
        write(99,*)'ecm,inc,aux=',ecm,inc,aux
        kch(ich)=sqrt(conv*aux)
        etav(ich)=conv*z12*e2/kch(ich)/2.
        rturnv(ich) =(etav(ich)+SQRT(etav(ich)**2 
     &              + ql(ich)*(ql(ich)+1d0)))/kch(ich) 
        if(rturnv(ich).gt.rmorto) rmorto=rturnv(ich)
      else
        copen(ich)=.false.
        if (debug) write(99,300) ich,aux,"closed"
        kch(ich)  =sqrt(-conv*aux)
        rturnv(ich)=rm
        rmorto     =rm
      endif
300   format(3x,"Channel",i2," Final Ecm=",1f8.3," MeV ->",a6)
      enddo !ich


c classical turning point for incoming channel
      eta=conv*z12*e2/kch(inc)/2.
      l=ql(inc)
      RTURN =(eta+SQRT(eta**2 + L*(L+1d0)))/kch(inc)  

      if (1>2) write(*,'(6(5x,i3,2x,1f8.3))')
     &   (ich,rturnv(ich),ich=1,nch) 

 
      call factorialgen(lmax+1)
      call cpu_time(start)
c 
c Starting values and initial integration point
c 
      rmin  =max(h,minval(rturnv(:)-abs(cutr)))
      irmin=nint(rmin/h)+1

!      write(*,*) 'Enhanced Numerov from ir=',irmin,' r=',rmin,' fm',
!     & ' to  R=',Rm,'in steps of',h,' and Rmax==',rm


      if (debug) write(95,*)'Starting value ir=,',irmin,'r=',rvec(irmin)
      z0(:,:)=zero; zm(:,:)=zero; 
      w0(:,:)=zero; wm(:,:)=zero
      y0(:,:)=zero; v(:,:)=zero 
      do ich=1,nch
        l=ql(ich)
        z0(ich,ich)= H**(l+1.) / EXP(0.5 * FACT(l+1))  
      enddo
      
      if (info) then
      write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
     &  5x,"outer=",1f8.2," fm",
     &  5x,"R-min=",1f6.2," fm")') 
     &  minval(rturnv(:)), maxval(rturnv(:)), rmin      
      endif

      if (irmin.lt.2) write(*,*)'** ERWIN: wrong irmin=',irmin
 
c Start radial integration
      do ir=irmin,nr
      r=rvec(ir) ! (ir-1)*h    
      im=ir-1  ! r-h
      i0=ir    ! r
      ip=ir+1  ! r+h
      if (r.lt.1e-6) r=1e-6

c kinetic energy part
c Fresco: SMAT(ich) = -LL1(K)*RI2 + (-ECM(ich) + C) * H
      do ich=1,nch    
      l=ql(ich)
      s(ich)= h2*(-l*(l+1)/r**2         ! centrifugal
     &          + kch2(ich)             ! energy
     &        - conv*vcoup(ich,ich,ir)) ! potential
      enddo

    
      do is=1,nch
      do ich=1,nch      
      y0(ich,is)= z0(ich,is)*
     X      (ONEC - S(ich) * (R12 - S(ich)*(ENA2 + S(ich)*ENA3)))
      enddo
      enddo

      w0(:,:)=zero
      do is=1,nch
      w0(is,is)=(ONEC - S(is) * (R12 - S(is)*(ENA2 + S(is)*ENA3)))
      enddo 
!      y0=matmul(w0,z0)

c non-diagonal part of V matrix 
      coupl(:,:)=-h2*conv*vcoup(:,:,ir)  
      do ich=1,nch
         coupl(ich,ich)=0d0
      enddo 

      w0=w0-r12*coupl
      
      if (1>2) then
      DO 24 ICH=1,NCH
	 DO 24 IS=1,NCH
        C = COUPL(is,ich) * R12
	 if(C/=ZERO) y0(ich,1:nch) = y0(ich,1:nch) - C * z0(is,1:nch)
24    CONTINUE
        DO 34 ich=1,nch ! k
        V(ich,:)  = ZERO 
        DO 34 is=1,nch   ! j
	 C = COUPL(is,ich)
 	 if(C/=ZERO) V(ich,1:nch) = V(ich,1:nch) + C * y0(is,1:nch)
34    CONTINUE
      else  !..........Idem using builtin matrix multiplication 
        V(1:nch,1:nch)  = ZERO 
        y0=y0  - r12*matmul(coupl,z0)
        v =   matmul(coupl,y0)
      endif


c store solution 
      y(:,:,ir)=y0(:,:)       

c zm <- z0, wm <-w0
      do is=1,nch
      do ich=1,nch
      ziv=z0(ich,is)
      zpv=2.0*ziv - zm(ich,is) - v(ich,is) - s(ich)*y0(ich,is)
      zm(ich,is)=ziv
      z0(ich,is)=zpv
      wm(ich,is)=w0(ich,is)
      enddo  !ich
      enddo  !is


c Re-orthogonalize solution vectors ....................................
!      if (orto.and.(r.lt.rmorto).and.(ir.gt.irmin).and.(hort.gt.0)
!     & .and.(mod(ir-irmin,norto).eq.0).and.(ir.lt.nr-10)) then 
!        write(*,'(5x,"->orthogonalizing at ir,r=",i4,1f7.2)')ir,r  
!       call qrerwinz(y0,z0,zm,w0,nch,ir)      ! QR factorization
!       y(:,:,ir)=y0(:,:) ! new (orthogonolized) solution at ir
!      endif !orto
!........................................................................

      if (debug) then
       write(91,'(5x,1f6.3,3x,50g14.5)') r,(y(1,ich,ir),ich=1,nch)
      endif       
      enddo !ir

c 
c Match y with asymptotic solution 
c            
!      call matching4(ecm,z12,nch,ql,lmax,inc,nr,y,wf,phase,smat,
!     &               info,eph) ! gauss5 + 2 points
! Under development!
      call match_real(ecm,z12,nch,nopen,ql,lmax,inc,nr,y,wf,info)

      if (debug) write(klog,'(2f10.3,3x,10g14.6)') ecm,z12,ql(1),
     &  (phase(ich),ich=1,nch) 

      call flush(6)
      end subroutine schcc_erwinrc




c ...
c ... Enhanced Numerov (version of IJ Thompson used in Fresco) 
c ... 
      subroutine schcc_erwin(nch,ecm,z12,inc,ql,factor,dr,r0,
     & npt,wf,phase,smat,method,info)
      use nmrv,only: nr,h,vcoup,conv,ech,debug,rmin,hort,rvec,cutr
      use globals, only: verb
      use factorials
      use constants , only: e2
      use memory
      implicit none
      integer:: method
      logical:: copen(nch),orto,raynal,info
      integer:: nch,klog,npt,ql(nch)
      integer:: ir,irmin,is,isp,ich,inc,l,lmax
      integer:: i,j,k,im,ip,i0,norto
c     .........................................................
      real*8     :: ymem,ti,tf
      real*8 , parameter:: pi=acos(-1d0)
      real*8 :: ecm,tkch,eta,z12,r12,ena2,ena3
      real*8 :: dr,r,r0,rm,h2,factor
      real*8 :: rturn,rmorto,rmatch
      real*8 :: kch(nch),kch2(nch),rturnv(nch),etav(nch)
      real*8 :: start,finish,aux,big,small,tmin,tmax
c     .........................................................
      complex*16, intent(out):: wf(nch,nr)
      complex*16, dimension(nch,nch):: y0
      complex*16, dimension(nch,nch):: z0,zm
      complex*16 :: y(nch,nch,npt),coupl(nch,nch),v(nch,nch)
      complex*16 :: phase(nch),smat(nch)
      complex*16 :: s(nch),ziv,zpv
      complex*16 :: a(nch,nch)
      complex*16 :: c1,c2,c
      complex*16,parameter:: onec=(1D0,0D0),zero=(0d0,0d0)
c ... for QR factorization
      complex*16 :: HDIAG(nch),w0(nch,nch),wm(nch,nch)
c ... TEST (delete after debugging)
      complex*16:: yaux(nch,nch)
c     ---------------------------------------------------------
      R12 = 1D0/12D0
      ENA2 = 2D0/5D0 * R12**2
      ENA3 = - 4D0/35D0 * R12**3
!      ONEC = (1D0,0D0)
      TMAX = 20.
      TMIN = -125.
      big=huge(big)
      small=epsilon(small)
     
      h=dr
      h2=h*h
      conv=factor
      copen(:)=.true.
      nr=npt
      klog=99
      rmin=0
      lmax=maxval(ql)
      orto=.false.

c ... Initialize some variables 
      y=0d0      
      coupl=0

c ... Set TRUE for debugging
      debug=.false.

c ..............................................................................
      if (nch.lt.1) then 
        write(*,*)'Scatcc_erwin: nch=',nch; stop
      endif

      if (ecm.lt.1e-6) then
        write(*,*)'Scatcc_erwin: Ecm too small! =',ecm
        stop
      endif 
 
      if (allocated(rvec)) deallocate(rvec)
      allocate(rvec(nr))

      do ir=1,nr
        rvec(ir)=r0+(ir-1)*dr
!        write(222,*)rvec(ir),real(vcoup(1,1,ir)),ir
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius if derivative is used
     
   
      norto=nint(hort/h)
      if (abs(hort).gt.0.) orto=.true.

      if (info) then
        ymem=nch*nch*(npt+2)*lc16/1e6
        if (verb.ge.3) write(*,190) ymem
190     format(5x,"[ WFS require", 1f8.2," Mbytes ]")
        ymem=nch*nch*lc16/1e6
!        write(*,*)'Auxiliary arrays require 8x',
!     &  ymem,'=',ymem*8, 'Mbyes'
!        write(*,*)' Biggest real=',big
      endif 

      if (ecm.lt.0) then 
      write(*,200) ecm
200   format(4x,'Skipping negative energy:',1f7.3, ' MeV')
      wf(:,:)=0d0
      return
      endif

      rmorto=0
      do ich=1,nch
      aux=ecm+ech(inc)-ech(ich)
      kch2(ich)=conv*aux
!      write(223,*)'ich,E,kch2=',ich,aux,kch2(ich)

      if (aux.gt.0) then
        copen(ich)=.true.
        if (debug) write(99,300) ich,aux,"open"
        kch(ich)=sqrt(conv*aux)
        etav(ich)=conv*z12*e2/kch(ich)/2.
        rturnv(ich) =(etav(ich)+SQRT(etav(ich)**2 
     &              + ql(ich)*(ql(ich)+1d0)))/kch(ich) 
        if(rturnv(ich).gt.rmorto) rmorto=rturnv(ich)
      else
        copen(ich)=.false.
        if (debug) write(99,300) ich,aux,"closed"
        kch(ich)  =sqrt(-conv*aux)
        rturnv(ich)=rm
        rmorto     =rm
      endif
300   format(3x,"Channel",i2," Final Ecm=",1f8.3," MeV ->",a6)
      enddo


c classical turning point for incoming channel
      eta=conv*z12*e2/kch(inc)/2.
      l=ql(inc)
      RTURN =(eta+SQRT(eta**2 + L*(L+1d0)))/kch(inc)  

      if (1>2) write(*,'(6(5x,i3,2x,1f8.3))')
     &   (ich,rturnv(ich),ich=1,nch) 

 
      call factorialgen(lmax+1)
      call cpu_time(start)
c 
c Starting values and initial integration point
c 
!      rmin=h ;       irmin=1
!      rmin=max(h,rturn-10.) ! CHECK h or r0 ?
      rmin  =max(h,minval(rturnv(:)-abs(cutr)))
      irmin=nint(rmin/h)+1

!      write(*,*) 'Enhanced Numerov from ir=',irmin,' r=',rmin,' fm',
!     & ' to  R=',Rm,'in steps of',h,' and Rmax==',rm


      if (debug) write(95,*)'Starting value ir=,',irmin,'r=',rvec(irmin)
      z0(:,:)=zero; zm(:,:)=zero; 
      w0(:,:)=zero; wm(:,:)=zero
      y0(:,:)=zero; v(:,:)=zero 
      do ich=1,nch
        l=ql(ich)
!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ZI(L,L) = H**(QNF(9,L)+1) / EXP(0.5 * FACT(QNF(9,L)+1))
        z0(ich,ich)= 1e-10  !TEEEEEEEEEEEEEEEEEEEEEESSSSST!        
!        z0(ich,ich)= H**(l+1.) / EXP(0.5 * FACT(l+1.))  
!         write(*,*)'ich,z0=:',  ich,z0(ich,ich)   
!!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
      enddo
      
      if (info) then
      write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
     &  5x,"outer=",1f8.2," fm",
     &  5x,"R-min=",1f6.2," fm")') 
     &  minval(rturnv(:)), maxval(rturnv(:)), rmin      
      endif

      if (irmin.lt.2) write(*,*)'** ERWIN: wrong irmin=',irmin
 
c Start radial integration
      do ir=irmin,nr
      r=rvec(ir) ! (ir-1)*h    
      im=ir-1  ! r-h
      i0=ir    ! r
      ip=ir+1  ! r+h
      if (r.lt.1e-6) r=1e-6

c kinetic energy part
c Fresco: SMAT(ich) = -LL1(K)*RI2 + (-ECM(ich) + C) * H
      do ich=1,nch    
      l=ql(ich)
      s(ich)= h2*(-l*(l+1)/r**2         ! centrifugal
     &          + kch2(ich)             ! energy
     &        - conv*vcoup(ich,ich,ir)) ! potential
      enddo
    
      do is=1,nch
      do ich=1,nch      
      y0(ich,is)= z0(ich,is)*
     X      (ONEC - S(ich) * (R12 - S(ich)*(ENA2 + S(ich)*ENA3)))
      enddo
      enddo

      w0(:,:)=zero
      do is=1,nch
      w0(is,is)=(ONEC - S(is) * (R12 - S(is)*(ENA2 + S(is)*ENA3)))
      enddo 
!      y0=matmul(w0,z0)

c non-diagonal part of V matrix 
      coupl(:,:)=-h2*conv*vcoup(:,:,ir)  
      do ich=1,nch
         coupl(ich,ich)=0d0
      enddo 

      w0=w0-r12*coupl
      
      if (1>2) then
      DO 24 ICH=1,NCH
	 DO 24 IS=1,NCH
        C = COUPL(is,ich) * R12
	 if(C/=ZERO) y0(ich,1:nch) = y0(ich,1:nch) - C * z0(is,1:nch)
24    CONTINUE
        DO 34 ich=1,nch ! k
        V(ich,:)  = ZERO 
        DO 34 is=1,nch   ! j
	 C = COUPL(is,ich)
 	 if(C/=ZERO) V(ich,1:nch) = V(ich,1:nch) + C * y0(is,1:nch)
34    CONTINUE
      else  !..........Idem using builtin matrix multiplication 
        V(1:nch,1:nch)  = ZERO 
        y0=y0  - r12*matmul(coupl,z0)
        v =   matmul(coupl,y0)
      endif

!! TEST 
      if (debug) then
      yaux=matmul(w0,z0)
      write(*,*)'y0 at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (y0(ich,is), is=1,min(15,nch))
      enddo

!      write(*,*)'w0*z0 at ir,r=',ir,r
!      do ich=1,min(15,nch)
!        write(*,'(5x,i3,50g14.5)') ich,
!     &  (yaux(ich,is), is=1,min(15,nch))
!      enddo
      endif

!!!!!!!!!!!!!!!!!!!!!!1


      if (debug) then
!      write(*,*)'coupl at ir,r=',ir,r
!      do ich=1,min(15,nch)
!        write(*,'(5x,i3,50g14.5)') ich,
!     &  (coupl(ich,is), is=1,min(15,nch))
!      enddo

!      write(*,*)'V at ir,r=',ir,r
!      do ich=1,min(15,nch)
!        write(*,'(5x,i3,50g14.5)') ich,
!     &  (V(ich,is), is=1,min(15,nch))
!      enddo

      write(*,*)'z0 at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (z0(ich,is), is=1,min(15,nch))
      enddo

      endif ! debug


c store solution 
      y(:,:,ir)=y0(:,:)       


c zm <- z0, wm <-w0
      do is=1,nch
      do ich=1,nch
      ziv=z0(ich,is)
      zpv=2.0*ziv - zm(ich,is) - v(ich,is) - s(ich)*y0(ich,is)
      zm(ich,is)=ziv
      z0(ich,is)=zpv
      wm(ich,is)=w0(ich,is)
      enddo  !ich
      enddo  !is


c Re-orthogonalize solution vectors ....................................
      if (orto.and.(r.lt.rmorto).and.(ir.gt.irmin).and.(hort.gt.0)
     & .and.(mod(ir-irmin,norto).eq.0).and.(ir.lt.nr-10)) then 
        if (verb.ge.2) 
     &   write(*,'(5x,"->orthogonalizing at ir,r=",i4,1f7.2)')ir,r
  
      if(1>2) then
c check orthogonality
      a=matmul(conjg(transpose(y0)),y0)
      write(*,*)'Unitarity (before) at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (a(ich,is)/a(1,1), is=1,min(15,nch))
      enddo
      endif

      if (1>2) then 
      write(*,*)'y0 (before) at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (y0(ich,is), is=1,min(15,nch))
      enddo
      endif

       call cpu_time(ti)
       call qrerwinz(y0,z0,zm,w0,nch,ir)      ! QR factorization

       y(:,:,ir)=y0(:,:) ! new (orthogonolized) solution at ir

      if (debug) then
      write(*,*)'y0 (after) at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (y0(ich,is), is=1,min(15,nch))
      enddo

        write(*,*)'|yp x yp| (after QR)'
        do is=1,nch
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(y0(:,is),y0(:,ich)))
     &      / abs(dot_product(y0(:,is),y0(:,is)))
     &      / abs(dot_product(y0(:,ich),y0(:,ich))), ich=1,nch) 
        enddo
      endif

       call cpu_time(tf)
       torto=torto+tf-ti

      endif !orto
!........................................................................


!!!!!!!!!!!!!!!!!!!!!!!!! TEST
      if (debug) then
       write(91,'(5x,1f6.3,3x,50g14.5)') r,(abs(y(1,ich,ir)),ich=1,nch)
      endif       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         DO 46 IT=1,NEQS
!         DO 43 K=1,NEQS
!            ZIV = ZI(IT,K)
!            ZPV = ZIV + ZIV - ZM(IT,K) - V(IT,K) - SMAT(K) * FI(IT,K)
!            ZM(IT,K) = ZIV
!            ZI(IT,K) = ZPV
!43            CONTINUE
!46    CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      enddo !ir
      call cpu_time(finish)


c 
c Match y with asymptotic solution -> S-matrix
c            
      call cpu_time(start)
!      call matching(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! gauss5 + derivative
!      call matching2(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! LAPACK 
      call matching3(ecm,z12,nch,ql,lmax,inc,nr,y,wf,phase,smat,info) ! gauss5 + 2 points
       call cpu_time(finish)
       tmatch=tmatch+finish-start

      if (debug) write(klog,'(2f10.3,3x,10g14.6)') ecm,z12,ql(1),
     &  (phase(ich),ich=1,nch) 


      call flush(6)
      end subroutine schcc_erwin


c
c ... Orthogonalization by QR factorization
c
      subroutine qrerwin(yi,zi,zm,wi,wm,y,n,ir,irmin,nr)
      use nmrv, only: h,rvec
      implicit none
      logical    debug
      integer    nr
      complex*16 y(n,n,nr)
c ... for zgeqrf
      integer info,lda,lwork,n,m,i,is,ich,ir,irmin
      complex*16 a(n,n),tau(n),work(2*n),rt(n,n)
      complex*16,dimension(n,n):: yi,ym,zi,zm,wi,wm
c ... for ztrsm
      complex*16 alpha
      integer ldb
      character diag,side,transa,uplo
c ...
      integer   :: nrhs, ipiv(n)
      character*1  trans
      EXTERNAL  ZGETRF, ZGETRS

c ---- TEST      
      integer ndim,lworkx
      complex*16 aux(3,3),taux(3),raux(3,3),workx(6)
      complex*16 zaux(n,n),waux(n,n)
c .................................................
      lwork=2*n
      debug=.false.
     

c ... QR factorization YI=Q.R 
c     (Q=orthogonal matrix; R= triangular UPPER matrix)
      rt=yi
      call ZGEQRF( N, N, RT, n, TAU, WORK, LWORK, INFO )
      if (info.ne.0)  then
         write(*,*)'zgeqrf failed; info=',info
      endif


c ... Apply transformation backwards
c       Y(old)= Y(new) R 
c     Solve:
c        X    *  RT  = Y(old)
c     and make: Y(new)= X 
      
c ... X*op( A ) = alpha*B
      side  ='R' ! A acts on the right
      uplo  ='U' ! R is upper triangle matrix
      transa='N' ! no transpose
      diag  ='N' !??????????
      m=n
      alpha =1d0
      lda   =n   ! leading dim of A=RT
      ldb   =n   ! leading dim of B=Y(old)

c Apply transformation to YI,YM
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,YI,LDB)
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,YM,LDB)
!      y(:,:,ir  )=yi(:,:) ! new (orthogonolized) solution at ir
!      y(:,:,ir-1)=ym(:,:) ! idem at ir-1
      
 
c Apply the same transformation to ZI,ZM
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,ZI,LDB)
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,ZM,LDB)

      if (debug) then
      write(*,*)'yi(new) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(yi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

        write(*,*)'|yp x yp| (after QR)'
        do is=1,n
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(yi(:,is),yi(:,ich)))
     &      / abs(dot_product(yi(:,is),yi(:,is)))
     &      / abs(dot_product(yi(:,ich),yi(:,ich))), ich=1,n) 

        enddo
       endif !debug



      RETURN  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      do i=irmin,ir
      yi(:,:)=y(:,:,i)
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,YI,LDB)
      y(:,:,i)=yi(:,:) ! override 
      enddo ! ir
      yi(:,:)=y(:,:,ir  ) ! new (orthogonolized) solution at ir
      ym(:,:)=y(:,:,ir-1) ! idem at ir-1






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c ... Recalculate ZI functions for new Y solutions applying the inverse transformation
c       YI(new) = W ZI(new)   -> ZI(new)
c       YM(new) = W ZM(new)   -> ZM(new) 
c     Use Lapack to solve:
c       A * X  =  alpha*B
c       W * Z  =  1    *Y

      write(*,*)'Zi (old) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

c ... First, compute LU factorization of matrix W to get IPIV
      call ZGETRF(n, n, wi, n, ipiv, info)
      if (info /= 0) then 
       print*, 'zgetrf returned error',info
      endif

c ... Second, use IPIV & zgetrs to solve system  W*Z= Y  -> Z 
      trans='N';  ldb=n;  nrhs=n
      zi(:,:)= yi(:,:) 
      call zgetrs(trans,n,nrhs,wi,n,ipiv,zi,ldb,info)

      if (info /= 0) then
        print*,'zgetrs failed at qrerwin!'
      endif

      write(*,*)'Zi (new) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

c ... Repeat for  YM(new) = W ZM(new) 
      call ZGETRF(n, n, wm, n, ipiv, info)
      if (info /= 0) stop 'zgetrf returned error'
      zm(:,:)= ym(:,:) 

      call zgetrs(trans,n,nrhs,wm,n,ipiv,zm,ldb,info)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      Return

c ......................... TEST ........................

      if (1>2) then
      write(*,*)'R (triangular) matrix:'
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(RT(ich,is)), is=1,min(15,n))
      enddo

      write(*,*)'yi(old) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(y(ich,is,ir)), is=1,min(15,n))
      enddo
      write(*,*)' '

      endif


! TEST THAT Y=W*Z
      if (1>2) then 
      write(*,*)'y(old) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(y(ich,is,ir)), is=1,min(15,n))
      enddo
      write(*,*)' '

      write(*,*)'yi(old) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(yi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

!! TEST that Z=W^-1*Y
      write(*,*)'Zi (arg) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

      waux=wi
      call ZGETRF(n, n, waux, n, ipiv, info)
      if (info /= 0) then
       write(*,*)'zgetrf returned error',info; stop
      endif

c ... Use IPIV & zgetrs to solve system  W*X= Y   (X=Z)
      trans='N';  ldb=n;  nrhs=n
      zaux(:,:)= yi(:,:) 
      call zgetrs(trans,n,nrhs,waux,n,ipiv,zaux,ldb,info)

      write(*,*)'Zi = W^-1*Y at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zaux(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

      if (info /= 0) then
          stop 'Matrix inversion failed!'
      else
       write(*,*)'zgetrs info=',info
       write(*,*)'ipiv=',ipiv
      endif


!!! TEST  Y= W*Z
      YM=MATMUL(WM,ZM)
      write(*,*)'W*ZM at ir-1=',ir-1
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(ym(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

      ym=y(:,:,ir-1)
      write(*,*)'y(old) at ir-1',ir-1
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(ym(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '
      endif

      if (1>2) then
        write(*,*)'|yi x yi| (before QR)'
        do is=1,n
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(y(:,is,ir),y(:,ich,ir)))
     &      / abs(dot_product(y(:,is,ir),y(:,is,ir)))
     &      / abs(dot_product(y(:,ich,ir),y(:,ich,ir))), 
     &     ich=1,n) 

        enddo
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (0>1) then
      aux(1,1)=12; aux(1,2)=-51; aux(1,3)=4
      aux(2,1)=6;  aux(2,2)=167; aux(2,3)=-68
      aux(3,1)=-4; aux(3,2)=24; aux(3,3)=-41
      lda=3
      lworkx=6
      raux=aux
      call ZGEQRF( lda, lda, raux, lda, TAUX, WORKX, LWORKX, INFO )
      if (info.ne.0)  then
         write(*,*)'zgeqrf failed; info=',info
      endif

      write(*,*)'R matrix:'
      do ich=1,3
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(raux(ich,is)), is=1,lda)
      enddo

c ... X*op( A ) = alpha*B
      side  ='R'
      uplo  ='U'
      transa='N'
      diag  ='N' !??????????
      m=n
      alpha =1d0
      lda   =3
      ldb   =3
      m     =lda
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,m,ALPHA,RAUX,LDA,AUX,LDB)

      write(*,*)'Q matrix:'
      do ich=1,3
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(aux(ich,is)), is=1,lda)
      enddo


        write(*,*)'|yp x yp| (before QR)'
        do is=1,3
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(aux(:,is),aux(:,ich)))
     &      / abs(dot_product(aux(:,is),aux(:,is)))
     &      / abs(dot_product(aux(:,ich),aux(:,ich))), ich=1,3) 

        enddo

      stop
      endif

c .......................................................
      end subroutine


c *** ----------------------------------------
C ... Re-orthogonalize Z's by QR factorization
c *** ----------------------------------------
      subroutine qrerwinz(yi,zi,zm,wi,n,ir)
      use nmrv, only: h,rvec
      implicit none
      logical    debug
c ... for zgeqrf
      integer info,lda,lwork,n,m,i,is,ich,ir
      complex*16 a(n,n),tau(n),work(2*n),rt(n,n)
      complex*16,dimension(n,n):: yi,zi,zm,wi
c ... for ztrsm
      complex*16 alpha
      integer ldb
      character diag,side,transa,uplo
c ...
      integer   :: nrhs, ipiv(n)
      character*1  trans
      EXTERNAL  ZGETRF, ZGETRS

c ---- TEST      
      complex*16 zaux(n,n),waux(n,n)
c .................................................
      lwork=2*n
      debug=.false.
     

c ... QR factorization ZI=Q.RT 
c     (Q=orthogonal matrix; R= triangular UPPER matrix)
      rt=zi
      call ZGEQRF( N, N, RT, n, TAU, WORK, LWORK, INFO )
      if (info.ne.0)  then
         write(*,*)'zgeqrf failed with exit code=',info
      endif
 
c ... Apply same transformation to ZM
c       Z(old)= Z(new) RT 
c     Solve:
c       X*RT  = Z(old)
c     and make: Y(new)= X 
      
c ... X*op( A ) = alpha*B
      side  ='R' ! A acts on the right
      uplo  ='U' ! R is upper triangle matrix
      transa='N' ! no transpose
      diag  ='N' !??????????
      m=n
      alpha =1d0
      lda   =n   ! leading dim of A=RT
      ldb   =n   ! leading dim of B=Y(old)

      if ((debug)) then ! .and.(mod(ir,5).eq.0)) then

      write(*,*)'Re[zi(old)] at ir',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '


        write(*,*)'|zi x zi| (before QR)'
        do is=1,n
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(zi(:,is),zi(:,ich)))
     &      / sqrt(abs(dot_product(zi(:,is),zi(:,is))))
     &      / sqrt(abs(dot_product(zi(:,ich),zi(:,ich)))), ich=1,n) 
        enddo
      endif
      
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,ZI,LDB)
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,ZM,LDB)


      if ((debug)) then ! .and.(mod(ir,5).eq.0)) then

      write(*,*)'Re[zi(new)] at ir',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '


        write(*,*)'|zi x zi| (after QR)'
        do is=1,n    
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(zi(:,is),zi(:,ich)))
     &      / abs(dot_product(zi(:,is),zi(:,is)))
     &      / abs(dot_product(zi(:,ich),zi(:,ich))), ich=1,n) 
        enddo
        write(*,*)'' 
      endif
     
c ... Recalculate yi,ym  
! MGR nov/16
!      yi=matmul(wi,zi) ! new (orthogonolized) solution at ir
      yi=matmul(transpose(wi),zi) ! new (orthogonolized) solution at ir

!      ym=matmul(wm,zm) ! idem at ir-1 (NOT NEEDED)

      RETURN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c ... Use IPIV & zgetrs to solve system  W*X= Y   (X=Z)
      trans='N';  ldb=n;  nrhs=n
      zaux(:,:)= yi(:,:) 
      call zgetrs(trans,n,nrhs,waux,n,ipiv,zaux,ldb,info)

      if (debug) then 
      write(*,*)'Zi = W^-1*Y at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zaux(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '
      endif

      if (info /= 0) then
          stop 'Matrix inversion failed!'
      else
       write(*,*)'zgetrs info=',info
       write(*,*)'ipiv=',ipiv
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (1>2) then
      write(*,*)'R (triangular) matrix:'
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(RT(ich,is)), is=1,min(15,n))
      enddo
      endif


      if ((debug).and.(mod(ir,5).eq.0)) then
      write(*,*)'yi(new) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(yi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

        write(*,*)'|yp x yp| (after QR)'
        do is=1,n
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(yi(:,is),yi(:,ich)))
     &      / abs(dot_product(yi(:,is),yi(:,is)))
     &      / abs(dot_product(yi(:,ich),yi(:,ich))), ich=1,n) 

        enddo
       endif !debug

c ... Recalculate ZI functions for new Y solutions applying the inverse transformation
c       YI(new) = W ZI(new)   -> ZI(new)
c       YM(new) = W ZM(new)   -> ZM(new) 
c     Use Lapack to solve:
c       A * X  =  alpha*B
c       W * Z  =  1    *Y

      if (debug) then
      write(*,*)'Zi (old) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '
      endif

c ... First, compute LU factorization of matrix A (=W) -> IPIV
      call ZGETRF(n, n, wi, n, ipiv, info)
      if (info /= 0) then
       write(*,*)'zgetrf returned error',info; stop
      endif

c ... Use IPIV & zgetrs to solve system  W*X= Y   (X=Z)
      trans='N';  ldb=n;  nrhs=n
      zi(:,:)= yi(:,:) 
      call zgetrs(trans,n,nrhs,wi,n,ipiv,zi,ldb,info)

      if (info /= 0) stop 'zgetrs failed at qrerwin!'
     
      if (debug) then
      write(*,*)'Zi (new) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '
      endif

c ... Repeat for  YM(new) = W ZM(new) 
!      call ZGETRF(n, n, wm, n, ipiv, info)
!      if (info /= 0) then
!       write(*,*)'zgetrf returned error',info; stop
!      endif

!      zm(:,:)= ym(:,:) 
!      call zgetrs(trans,n,nrhs,wm,n,ipiv,zm,ldb,info)

      end subroutine






	SUBROUTINE POTWF2(COUPL,ZI,FI,V,NEQS,MAXB,MAXCH,R12)
	IMPLICIT NONE
	REAL*8 R12
	INTEGER NEQS,MAXB,K,J,IT,MAXCH
	COMPLEX*16 COUPL(MAXB,MAXB),ZI(MAXB,MAXCH),FI(MAXB,MAXCH),
     X		   V(MAXB,MAXCH),C,ZERO
	PARAMETER(ZERO = (0d0,0d0))
C					VECTORISED MATRIX MULTIPLIES
C					ALLOWS SKIPS IF ZERO COUPLING

         DO 24 K=1,NEQS
	     DO 24 J=1,NEQS
             C = COUPL(K,J) * R12
	      if(C/=ZERO) FI(1:NEQS,K) = FI(1:NEQS,K) - C * ZI(1:NEQS,J)
 24       CONTINUE

         DO 34 K=1,NEQS
	    V(:,K)  = ZERO
            DO 34 J=1,NEQS
	      C = COUPL(K,J)
	      if(C/=ZERO) V(1:NEQS,K) = V(1:NEQS,K) + C * FI(1:NEQS,J)
34       CONTINUE
	return
	end




c -------------------------------------------------------------
c Match with asymptotics solution to get wfs and S-matrix
c The wfs are obtained as a linear combination of nchan independent
c solutions, verifying: 
c
c f(r)_{n,inc} -> 0.5 i[H(-)_{li} delta(li,ln) + S_{n,i}H(+)_{l}]
c 
c Note that some authors include additional factor: e^{sigma_c }
c in the definition of f(r). 
c
c ----------------------------------------------- --------------
      subroutine matching(ecm,z12,nch,ql,lmax,inc,
     & rmatch,y,wf,phase,smat,show)
      use nmrv, only: nr,mu,ech,h,conv
      use constants, only : e2,pi
      use globals, only: verb,debug
      implicit none
      logical :: sing,show
      integer ich,is,inc,ir,klog
      integer nch,ql(nch),l,linc,lmax,ifail,ie,nmatch

      real*8, dimension(0:lmax):: f,g,gp,fp
      real*8, dimension(1:lmax+1):: f1,fp1
      real*8 :: kron,rmatch
      real*8 :: ecm,tkch,eta,krm,z12,kch(nch)
      real*8 :: phr,cph(0:lmax) !!!! changed to lmax check !!!!
      real*8 :: fpmax,small,acc8,flux,magn(nr)

      complex*16:: wf(nch,nr)
      complex*16:: yd(nch,nch),yaux(5)
      complex*16:: y(nch,nch,nr)
      complex*16:: ch,chd,mat(2*nch,2*nch+1),det
      complex*16:: ci,test,yfac,tmat
      complex*16:: phase(nch),phc
      complex*16:: smat(nch),delta,scoul,anc,svel

      acc8 = epsilon(acc8);
      fpmax = huge(acc8)**0.8d0 
      SMALL=1D0/FPMAX
      ci=(0d0,1d0)
      pi=acos(-1d0) 
      debug=.false.
!      show=.false.
      klog=99
      magn(:)=0      

      if (nch.eq.0) then
       write(*,*)'Matching: nch=0! Abort'; stop     
      endif
      mat(:,:)=0d0
!      mat(1:nch,1:nch)=y(1:nch,1:nch,nr)
      mat(1:nch,1:nch)=y(1:nch,1:nch,nr-2) ! matrix solution at r=rmatch
      nmatch=nr-2                        ! matching done at ir=nr-2 

c
c derivatives at matching point
c 
      if (debug) write(99,*)'Derivatives at  rmatch=',rmatch
      do ich=1,nch
      do is=1,nch
!      yaux(1:5)= y(ich,is,nr-2:nr+2)
      yaux(1:5)= y(ich,is,nr-4:nr) ! use last 5 points for derivative
!      yd(ich,is)=deriv1(yaux,h,nr+2,nr)
      yd(ich,is)=(yaux(1)
     &   -8*yaux(2)+8*yaux(4)-yaux(5))/(h*12)
      mat(nch+ich,is)= yd(ich,is)
c3      mat(nch+ich,is)= y(ich,is,nr-1)
      enddo
      if (debug) write(99,'(5x,50f12.6)') (yd(ich,is),is=1,nch) 
      enddo
      
      linc=ql(inc)
      do ich=1,nch
      l=ql(ich) 
      tkch=ecm+ech(inc)-ech(ich) 

      kch(ich)=sqrt(conv*abs(tkch))
      krm=kch(ich)*rmatch
      eta=conv*z12*e2/kch(ich)/2.

      if (tkch.gt.0d0) then   ! open channel
       call coulph(eta,cph,linc)
       phc=exp(ci*cph(linc))

!       write(190,'(6f10.4)') tkch,kch(ich),eta,krm,cph(linc)
       call coul90(krm,eta,0d0,l,f,g,fp,gp,0,IFAIL)


       if (ifail.ne.0) then 
       write(*,*) 'coul90: ifail=',ifail; stop
       endif
       ch  = dcmplx(g(l),f(l))  * (0.,.5)   ! outgoing (i/2) H(+)
       chd = dcmplx(gp(l),fp(l)) * (0.,.5) ! (i/2) H'(+)
c3      call coul90(kch(ich)*(rm-h),eta,0d0,l,f,g,fp,gp,0,IFAIL)
c3      chd = dcmplx(g(l),f(l)) * (0.,.5) 
c lhs
      mat(ich,nch+ich)    = ch
      mat(nch+ich,nch+ich)= chd*kch(ich)
c rhs
      mat(ich,2*nch+1)    =-kron(ich,inc)*conjg(ch) 
      mat(nch+ich,2*nch+1)=-kron(ich,inc)*
     &                    conjg(chd)*kch(ich) 
      else                      ! closed channel
       IE = 0  
       call whit(eta,rmatch,kch(ich),tkch,l,f1,fp1,ie)
       ch = -f1(l+1)  ! Whittaker
       chd= -fp1(l+1) ! derivative 
!       write(0,*) 'l,k,f',l,kch(ich),f1(l+1)
c lhs
      mat(ich,nch+ich)    = -ch
      mat(nch+ich,nch+ich)= -chd*kch(ich)
c rhs
      mat(ich,2*nch+1)    =0. 
      mat(nch+ich,2*nch+1)=0. 
      endif
c2       mat(ich,nch+ich)    = -(g(l)+ci*f(l))
c2       mat(nch+ich,nch+ich)= -(gp(l)+ci*fp(l))*kch(ich)
c3        mat(ich,nch+ich)    = ch
c3        mat(nch+ich,nch+ich)    = chd
c2      mat(ich,2*nch+1)    =kron(ich,inc)*f(l)
c2      mat(nch+ich,2*nch+1)=kron(ich,inc)*fp(l)*kch(ich)
c3      mat(ich,2*nch+1)    =-kron(ich,inc)*conjg(ch)
c3      mat(nch+ich,2*nch+1)=-kron(ich,inc)*conjg(chd)

      enddo
      call GAUSS5(2*nch,2*nch,mat,1,sing,det,small,debug)
      if (debug) write(klog,'(5x,"Coefficients:",100f8.3)') 
     & (mat(ich,2*nch+1),ich=1,nch)

      IF(SING) then
      DO 605 Ir=1,nr
605   IF(MAGN(Ir).NE.0.0) MAGN(Ir) = LOG10(MAGN(Ir))
      WRITE(*,610) (MAGN(Ir),Ir=1,nr)
610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
      endif

c
c S-matrix and phase-shifts
c 
      flux=0d0
      phase(1:nch)=0.
      if (show) 
!     &  write(*,'(8x,"S-matrix (wo/with velocity factors):")') 
     &   write(*,'(a40,6x,a45)')
     &   "S-matrix (WITHOUT velocity factors)", 
     &   "S-matrix (WITH velocity factors)"
      do ich=1,nch
      tkch=ecm+ech(inc)-ech(ich) 
c2      tmat=mat(nch+ich,2*nch+1)
c2      smat(ich)=(tmat-1)*(0.,-0.5)
      if (tkch.gt.0) then                                   ! open channel
!commented v2.1b
!      smat(ich)=mat(nch+ich,2*nch+1)*ci**(QL(INC)-QL(ICH))
       smat(ich)=mat(nch+ich,2*nch+1)
       svel=smat(ich)*sqrt(kch(ich)/kch(inc)) ! S-matrix with velocity factors
       phase(ich)=(0.,-0.5) * LOG(mat(nch+ich,2*nch+1))
       test=(0.,-0.5) * LOG(mat(nch+ich,2*nch+1)) 
       if (((show).and.(ich.eq.inc)).or.(verb.ge.1)) then
!       write(*,500) ich,smat(ich),svel,abs(svel),test*180/pi
       write(*,500) ich,smat(ich),svel,abs(svel) !,test*180/pi
       endif
       flux=flux + abs(svel)**2
500    format(10x,"Chan. #",i3,5x,
     &    "S=(",1f10.6,",",1f10.6,")", 5x, 
     &    "S=(",1f10.6,",",1f10.6,") ->  |S|=",1f9.6)
!     &    ,"delta=(",1f8.3,",",1f8.3,")")

      else                                                 ! closed channel
        smat(ich)=mat(nch+ich,2*nch+1)
        write(*,500) ich,smat(ich),svel,abs(svel),test*180/pi
        if (verb.gt.3) then
        anc=mat(nch+ich,2*nch+1)
        if (debug)write(*,510) ich,anc,abs(anc)
510     format(5x,"Channel",i3,3x,
     &    "ANC=(",1g12.6,",",1g12.6,") -> |ANC|=",1g12.6)
      endif
      endif
      enddo !ich

      phase(1:nch)=phase(1:nch)*180/pi
!      write(45,'(1f10.3,3x,10g14.6)') ecm,
!    &  (phase(ich),ich=1,nch) 
       
      write(*,'(10x,"=> Unitarity=",1f12.8)')flux


c
c scattering wavefunctions
c
      wf(1:nch,1:nr)=(0d0,0d0)
      if(abs(smat(inc)).gt.1e-20) then
        phr = (0.,-0.5) * LOG(smat(inc))
        yfac=exp(-ci*phr)
      endif

      do ich=1,nch
      do is=1,nch
      wf(ich,1:nr)=wf(ich,1:nr)+mat(is,2*nch+1)*y(ich,is,1:nr) 
!      wf(ich,1:nr)=wf(ich,1:nr)+xmat(is,1)*y(ich,is,1:nr)
      enddo
      enddo

!      wf(:,:)=dsqrt(2d0/pi)*phc*yfac*wf(:,:)
ccc      wf(:,:)=dsqrt(2d0/pi)*yfac*wf(:,:)
! TEST
!      wf(:,:)=dsqrt(2d0/pi)*wf(:,:)
! TO MATCH WITH FRESCO FORT.17 , NO FACTORS ARE NEEDED HERE!


c      do ir=1,nr
c      write(51,'(2x,1f6.3,2x,10f14.8)')ir*h,(wf(ich,ir),ich=1,nch)
c      enddo

      end subroutine



c same as matching, but using LAPACK to solve system of eqns
      subroutine matching2(ecm,z12,nch,ql,lmax,inc,
     & rmatch,y,wf,phase,smat)
      use nmrv, only: nr,mu,ech,h,conv
      use constants, only : e2,pi
      use globals, only: verb
      implicit none
!!!!!!!!!!!!!!!!!!!  LAPACK    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer  nch
      integer  iter,info
      integer,  parameter:: nrhs=1
!      integer          ::  n=2*nch
      integer    ipiv(2*nch)
      double precision :: rwork(2*nch) 
      complex    swork(2*nch*(2*nch+nrhs))
      complex*16 mat(2*nch,2*nch),bmat(2*nch,nrhs),xmat(2*nch,nrhs),
     &           work(2*nch)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      logical :: debug !sing,show
      integer ich,is,inc,ir,klog
      integer ql(nch),l,linc,lmax,ifail,ie

      real*8, dimension(0:lmax):: f,g,gp,fp
      real*8, dimension(1:lmax+1):: f1,fp1
      real*8 :: kron,rmatch
      real*8 :: ecm,tkch,eta,krm,z12,kch(nch)
      real*8 :: phr,cph(0:lmax) !!!! changed to lmax check !!!!
      real*8 :: fpmax,small,acc8,flux,magn(nr)

      complex*16:: wf(nch,nr)
      complex*16:: yd(nch,nch),yaux(5)
      complex*16:: y(nch,nch,nr+2)
      complex*16:: ch,chd,det
      complex*16:: ci,test,yfac,tmat
      complex*16:: phase(nch),phc
      complex*16:: smat(nch),delta,scoul,anc,svel

      acc8 = epsilon(acc8);
      fpmax = huge(acc8)**0.8d0 
      SMALL=1D0/FPMAX
      ci=(0d0,1d0)
      pi=acos(-1d0) 
      debug=.false.
      klog=99
      magn(:)=0      

      if (nch.eq.0) then
       write(*,*)'Matching: nch=0! Abort'; stop     
      endif
      mat(:,:)=0d0
      mat(1:nch,1:nch)=y(1:nch,1:nch,nr-2)

c
c derivatives at end point
c 
      if (debug) write(99,*)'Derivatives at end point rm=',rmatch
      do ich=1,nch
      do is=1,nch
!      yaux(1:5)= y(ich,is,nr-2:nr+2)
      yaux(1:5)= y(ich,is,nr-4:nr) ! use last 5 points for derivative
!      yd(ich,is)=deriv1(yaux,h,nr+2,nr)
      yd(ich,is)=(yaux(1)-
     &   8*yaux(2)+8*yaux(4)-yaux(5))/(h*12)
      mat(nch+ich,is)= yd(ich,is)
c3      mat(nch+ich,is)= y(ich,is,nr-1)
      enddo
      if (debug) write(99,'(5x,50f12.6)') (yd(ich,is),is=1,nch) 
      enddo
      
      linc=ql(inc)
      do ich=1,nch
      l=ql(ich) 
      tkch=ecm+ech(inc)-ech(ich) 

      kch(ich)=sqrt(conv*abs(tkch))
      krm=kch(ich)*rmatch
      eta=conv*z12*e2/kch(ich)/2.

      if (tkch.gt.0d0) then   ! open channel
       call coulph(eta,cph,linc)
       phc=exp(ci*cph(linc))
       call coul90(krm,eta,0d0,l,f,g,fp,gp,0,IFAIL)


       if (ifail.ne.0) then 
       write(*,*) 'coul90: ifail=',ifail; stop
       endif
       ch  = dcmplx(g(l),f(l))  * (0.,.5)   ! outgoing (i/2) H(+)
       chd = dcmplx(gp(l),fp(l)) * (0.,.5) ! (i/2) H'(+)
c3      call coul90(kch(ich)*(rm-h),eta,0d0,l,f,g,fp,gp,0,IFAIL)
c3      chd = dcmplx(g(l),f(l)) * (0.,.5) 
c lhs
      mat(ich,nch+ich)    = ch
      mat(nch+ich,nch+ich)= chd*kch(ich)
c rhs
      bmat(ich,1)      =-kron(ich,inc)*conjg(ch) 
      bmat(nch+ich, 1) =-kron(ich,inc)*conjg(chd)*kch(ich) 
      else                      ! closed channel
       IE = 0  
       call whit(eta,rmatch,kch(ich),tkch,l,f1,fp1,ie)
       ch = -f1(l+1) !* (0.,.5) ! Whittaker*i/2
       chd= -fp1(l+1)!* (0.,.5) ! derivative 
c lhs
      mat(ich,nch+ich)    = -ch
      mat(nch+ich,nch+ich)= -chd*kch(ich)
c rhs
      bmat(ich,    1) =0. 
      bmat(nch+ich,1) =0. 
      endif
      enddo


c LAPACK SOLVER
      call zcgesv(2*nch,nrhs,mat,2*nch,ipiv,bmat,2*nch,xmat,2*nch,
     &           work,swork,rwork,iter,info)

c
c S-matrix and phase-shifts
c 
      flux=0d0
      phase(1:nch)=0.
      if (verb.ge.0) 
!     &  write(*,'(8x,"S-matrix (wo/with velocity factors):")') 
     &   write(*,'(a40,6x,a45)')
     &   "S-matrix (WITHOUT velocity factors)", 
     &   "S-matrix (WITH velocity factors)"
      do ich=1,nch
      tkch=ecm+ech(inc)-ech(ich) 
c2      tmat=mat(nch+ich,2*nch+1)
c2      smat(ich)=(tmat-1)*(0.,-0.5)
      if (tkch.gt.0) then                                   ! open channel
      smat(ich)=xmat(nch+ich,1)
      svel=smat(ich)*sqrt(kch(ich)/kch(inc)) ! S-matrix with velocity factors
      phase(ich)=(0.,-0.5) * LOG(smat(ich))
      test=(0.,-0.5) * LOG(smat(ich)) 
      if (((verb.eq.0).and.(ich.eq.inc)).or.(verb.ge.1)) then
       write(*,500) ich,smat(ich),svel,abs(svel),test*180/pi
      endif
      flux=flux + abs(svel)**2
500   format(10x,"Chan. #",i3,5x,
     &    "S=(",1f10.6,",",1f10.6,")", 5x, 
     &    "S=(",1f10.6,",",1f10.6,") ->  |S|=",1f9.6,3x, 
     &    "delta=(",1f8.4,",",1f8.4,")")

      else                                                  ! closed channel
        if (verb.gt.3) then
        anc=xmat(nch+ich,1)
        smat(ich)=anc
        if (debug)write(*,510) ich,anc,abs(anc)
510     format(5x,"Channel",i3,3x,
     &    "ANC=(",1g12.6,",",1g12.6,") -> |ANC|=",1g12.6)
      endif
      endif
      enddo !ich
      phase(1:nch)=phase(1:nch)*180/pi
       
      write(*,'(10x,"=> Unitarity=",1f12.8)')flux
      
c
c scattering wavefunctions
c
      wf(1:nch,1:nr)=(0d0,0d0)
      if(abs(smat(inc)).gt.1e-20) then
        phr = (0.,-0.5) * LOG(smat(inc))
        yfac=exp(-ci*phr)
      endif

      do ich=1,nch
      do is=1,nch
      wf(ich,1:nr)=wf(ich,1:nr)+xmat(is,1)*y(ich,is,1:nr)
      enddo
      enddo

!      wf(:,:)=dsqrt(2d0/pi)*phc*yfac*wf(:,:)
ccc      wf(:,:)=dsqrt(2d0/pi)*yfac*wf(:,:)
! TEST
!      wf(:,:)=dsqrt(2d0/pi)*wf(:,:)
! TO MATCH WITH FRESCO FORT.17 , NO FACTORS ARE NEEDED HERE!

!      do ir=1,nr
!      write(51,'(2x,1f6.3,2x,10f14.8)')ir*h,(wf(ich,ir),ich=1,nch)
!      enddo

      end subroutine





c -------------------------------------------------------------
c Match with asymptotics solution to get wfs and S-matrix
c The wfs are obtained as a linear combination of nchan independent
c solutions, verifying: 
c
c f(r)_{n,inc} -> 0.5 i[H(-)_{li} delta(li,ln) + S_{n,i}H(+)_{l}]
c              =
c 
c Note that some authors include additional factor: e^{sigma_c }
c in the definition of f(r). 
c
c USE 2  POINTS FOR MATCHING INSTEAD OF DERIVATIVE 
c ----------------------------------------------- --------------
      subroutine matching3(ecm,z12,nch,ql,lmax,inc,
     & n1,y,wf,phase,smat,show)
      use nmrv, only: nr,mu,ech,h,conv,rvec
      use constants, only : e2,pi
      use globals, only: verb,debug
      implicit none
      logical :: sing,show
      integer ich,is,inc,ir,klog,nd,n1,n2
      integer nch,ql(nch),l,linc,lmax,ifail,ie
c     --------------------------------------------------------------
      real*8, dimension(0:lmax):: f,g,gp,fp
      real*8, dimension(1:lmax+1):: f1,fp1
      real*8 :: kron,r1,r2
      real*8 :: ecm,tkch,eta,krm1,krm2,z12,kch(nch)
      real*8 :: phr,cph(0:lmax) !!!! changed to lmax check !!!!
      real*8 :: fpmax,small,acc8,flux,magn(nr)
c     ----------------------------------------------------------------
      complex*16:: wf(nch,nr)
      complex*16:: yd(nch,nch),yaux(5)
      complex*16:: y(nch,nch,nr)
      complex*16:: ch1,ch2,mat(2*nch,2*nch+1),det
      complex*16:: ci,test,yfac,tmat
      complex*16:: phase(nch),phc
      complex*16:: smat(nch),delta,scoul,anc,svel
c     ---------------------------------------------------------------
      acc8 = epsilon(acc8);
      fpmax = huge(acc8)**0.8d0 
      SMALL=1D0/FPMAX
      ci=(0d0,1d0)
      pi=acos(-1d0) 
      debug=.false.
!      show=.false.
      klog=99
      magn(:)=0      
      nd=6
      n2=n1-nd
 

      if (nch.eq.0) then
       write(*,*)'Matching: nch=0! Abort'; stop     
      endif
      mat(:,:)=0d0
!      mat(1:nch,1:nch)=y(1:nch,1:nch,nr)
      mat(1:nch,1:nch)=y(1:nch,1:nch,n1) 

      do ich=1,nch
      do is=1,nch
      mat(nch+ich,is)= y(ich,is,n1-nd)
      enddo
      enddo
      
      linc=ql(inc)
      do ich=1,nch
      l=ql(ich) 
      tkch=ecm+ech(inc)-ech(ich) 

      kch(ich)=sqrt(conv*abs(tkch))
      eta=conv*z12*e2/kch(ich)/2.

      r1=rvec(n1)
      r2=rvec(n2)
      krm1=kch(ich)*r1
      krm2=kch(ich)*r2

      if (tkch.gt.0d0) then   ! open channel
       call coulph(eta,cph,linc)
       phc=exp(ci*cph(linc))

!       write(190,'(6f10.4)') tkch,kch(ich),eta,krm1,cph(linc)
       call coul90(krm1,eta,0d0,l,f,g,fp,gp,0,IFAIL)
       if (ifail.ne.0) then 
       write(*,*) 'coul90: ifail=',ifail; stop
       endif
       ch1  = dcmplx(g(l),f(l))  * (0.,.5)   ! outgoing (i/2) H(+)
       call coul90(krm2,eta,0d0,l,f,g,fp,gp,0,IFAIL)
       ch2 = dcmplx(g(l),f(l)) * (0.,.5) 
c lhs
       mat(ich,nch+ich)    = ch1
       mat(nch+ich,nch+ich)= ch2
c rhs
       mat(ich,2*nch+1)    =-kron(ich,inc)*conjg(ch1) 
       mat(nch+ich,2*nch+1)=-kron(ich,inc)*conjg(ch2) 
      else !......................................... closed channel
       IE = 0  
       call whit(eta,r1,kch(ich),tkch,l,f1,fp1,ie)
       ch1 = f1(l+1)* (0.,.5) ! Whittaker*i/2
       call whit(eta,r2,kch(ich),tkch,l,f1,fp1,ie)
       ch2=  f1(l+1)* (0.,.5)  
c lhs
        mat(ich,nch+ich)    = -ch1
        mat(nch+ich,nch+ich)= -ch2
c rhs
        mat(ich,2*nch+1)    =0. 
        mat(nch+ich,2*nch+1)=0. 
      endif ! OPEN/CLOSED CHANNELS

c Commented in v2.3
c1      mat(ich,nch+ich)    = ch1
c1      mat(nch+ich,nch+ich)    = ch2
c1      mat(ich,2*nch+1)    =-kron(ich,inc)*conjg(ch1) ! repetido?
c1      mat(nch+ich,2*nch+1)=-kron(ich,inc)*conjg(ch2) ! repetido?? (ver arriba)

      enddo
      if (debug) show=.true.
      call GAUSS5(2*nch,2*nch,mat,1,sing,det,small,debug)
      if (debug) write(klog,'(5x,"Coefficients:",100f8.3)') 
     & (mat(ich,2*nch+1),ich=1,nch)

      IF(SING) then
      write(*,*)'Singular matrix for Ecm=',ecm
      DO 605 Ir=1,nr
605   IF(MAGN(Ir).NE.0.0) MAGN(Ir) = LOG10(MAGN(Ir))
      WRITE(*,610) (MAGN(Ir),Ir=1,nr)
610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
      endif

c
c S-matrix and phase-shifts
c 
      flux=0d0
      phase(1:nch)=0.
      if (show) 
!     &  write(*,'(8x,"S-matrix (wo/with velocity factors):")') 
     &   write(*,'(a40,6x,a45)')
     &   "S-matrix (WITHOUT velocity factors)", 
     &   "S-matrix (WITH velocity factors)"
      do ich=1,nch
      tkch=ecm+ech(inc)-ech(ich) 
      if (tkch.gt.0) then                    ! open channel
!commented v2.1b
!      smat(ich)=mat(nch+ich,2*nch+1)*ci**(QL(INC)-QL(ICH))
      smat(ich)=mat(nch+ich,2*nch+1)
      svel=smat(ich)*sqrt(kch(ich)/kch(inc)) ! S-matrix with velocity factors
      if (abs(mat(nch+ich,2*nch+1)).gt.small) then
        phase(ich)=(0.,-0.5) * LOG(mat(nch+ich,2*nch+1))
        test=(0.,-0.5) * LOG(mat(nch+ich,2*nch+1)) 
      endif

      if (show) then
       if (((ich.eq.inc).and.(verb.ge.1)).or.(verb.ge.3)) then
       write(*,500) ich,smat(ich),svel,abs(svel)!!!,test*180/pi
       endif
      endif
      flux=flux + abs(svel)**2
500   format(10x,"Chan. #",i3,5x,
     &    "S=(",1f8.5,",",1f8.5,")", 5x, 
     &    "S=(",1f8.5,",",1f8.5,") ->  |S|=",1f9.6) 
!     &    "delta=(",1f8.3,",",1f8.3,")")

      else                                   ! closed channel
        smat(ich)=mat(nch+ich,2*nch+1) ! Coefficient of Whittaker function
        svel=0; test=0
        if (show) write(*,501) ich,smat(ich)
501     format(10x,"Chan. #",i3,5x,
     &    "S=(",1f10.6,",",1f10.6,")", 5x, ' (closed channel)') 
        if (verb.gt.3) then
        anc=mat(nch+ich,2*nch+1)
        write(*,510) ich,anc,abs(anc)
510     format(5x,"Channel",i3,3x,
     &    "ANC=(",1g12.6,",",1g12.6,") -> |ANC|=",1g12.6)
      endif
      endif
      enddo !ich
      phase(1:nch)=phase(1:nch)*180/pi
!      write(45,'(1f10.3,3x,10g14.6)') ecm,
!    &  (phase(ich),ich=1,nch) 
       
      if (show) write(*,'(10x,"=> Unitarity=",1f12.8)')flux


c
c scattering wavefunctions
c
      wf(1:nch,1:nr)=(0d0,0d0)
      if(abs(smat(inc)).gt.1e-20) then
        phr = (0.,-0.5) * LOG(smat(inc))
        yfac=exp(-ci*phr)
      endif

      do ich=1,nch
      do is=1,nch
      wf(ich,1:nr)=wf(ich,1:nr)+mat(is,2*nch+1)*y(ich,is,1:nr) ! gauss5
      enddo
      enddo

!      wf(:,:)=dsqrt(2d0/pi)*phc*yfac*wf(:,:)
c To make the wfs real in the single-channel case:
!       wf(:,:)=dsqrt(2d0/pi)*yfac*wf(:,:)
! TEST
!      wf(:,:)=dsqrt(2d0/pi)*wf(:,:)
! TO MATCH WITH FRESCO FORT.17 , NO FACTORS ARE NEEDED HERE!

c TEST
!      do ir=1,nr
!      write(500,'(2x,1f6.3,2x,10f14.8)')ir*h,(wf(ich,ir),ich=1,nch)
!      enddo
!      write(500,*)'&'
      end subroutine




c -----------------------------------------------------------------
c Match with asymptotics solution to get wfs. Assume real indep. solutions
c
c USE 2  POINTS FOR MATCHING INSTEAD OF DERIVATIVE 
c ----------------------------------------------- --------------
      subroutine matching4(ecm,z12,nch,ql,lmax,iref,
     & n1,y,wf,phase,smat,show,eph)
      use nmrv, only: nr,mu,ech,h,conv,rvec
      use constants, only : e2,pi
      use globals, only: verb,debug
      implicit none
      logical :: sing,show,ifop(nch)
      integer ich,is,inc,iref,ir,klog,nd,n1,n2,iop
      integer nch,ql(nch),l,linc,lmax,ifail,ie,nopen
c     --------------------------------------------------------------
      real*8, dimension(0:lmax):: f,g,gp,fp
      real*8, dimension(1:lmax+1):: f1,fp1
      real*8 :: kron,r1,r2
      real*8 :: ecm,tkch,eta,krm1,krm2,z12,kch(nch)
      real*8 :: phr,cph(0:lmax) !!!! changed to lmax check !!!!
      real*8 :: fpmax,small,acc8,flux,magn(nr)
      real*8:: yd(nch,nch),yaux(5)
      real*8:: y(nch,nch,nr)
c     ----------------------------------------------------------------
      complex*16:: wf(nch,nr)
      complex*16:: ch1,ch2,mat(2*nch,2*nch+1),det
      complex*16:: ci,test,yfac,tmat
      complex*16:: phase(nch),phc
      complex*16:: smati(nch),delta,scoul,anc,svel,smat(nch,nch)
c     ---------------------------------------------------------------
c     For S-matrix diagonalization
      integer sort,n
      parameter (sort = 1)
      complex*16:: umat(nch,nch),evec(nch)
      real*8 :: eph(nch)
      real*8 :: WORK(NCH,NCH),
     &          WORK1(NCH),WORK2(3*NCH)
c     ...............................................................


c ... Initialize some variables
      acc8 = epsilon(acc8);
      fpmax = huge(acc8)**0.8d0 
      SMALL=1D0/FPMAX
      ci=(0d0,1d0)
      pi=acos(-1d0) 
      debug=.false.
!      show=.false.
      klog=99
      magn(:)=0      
      nd=6
      n2=n1-nd
      nopen=0; 
      ifop(:)=.false.
c .............................
 
      do ich=1,nch
       write(*,*)'matching4: inc,nch=',ich,ech(ich)
      enddo

      if (nch.eq.0) then
       write(*,*)'Matching: nch=0! Abort'; stop     
      endif
      
! NEW: calculated S-matrix for all open channels
      iop=0
      do inc=1,nch
      linc=ql(inc)
      tkch=ecm+ech(iref)-ech(inc) 
      if (tkch.lt.0) then
        smat(inc,:)=0
        write(*,*)' No inc waves for chan=',inc
        cycle
      endif

      mat(:,:)=0d0
      mat(1:nch,1:nch)=y(1:nch,1:nch,n1) 

      do ich=1,nch
      do is=1,nch
      mat(nch+ich,is)= y(ich,is,n1-nd)
      enddo
      enddo

      print*,'inc=',inc,'ech=',ech(inc),'iref=',iref,' Ecm=',tkch

      do ich=1,nch
      l=ql(ich) 
      tkch=ecm+ech(iref)-ech(ich) 

      kch(ich)=sqrt(conv*abs(tkch))
      eta=conv*z12*e2/kch(ich)/2.

      r1=rvec(n1)
      r2=rvec(n2)

      krm1=kch(ich)*r1
      krm2=kch(ich)*r2

      IF (tkch.gt.0d0) THEN   ! open channel
!       call coulph(eta,cph,linc)
!       phc=exp(ci*cph(linc))
!       write(190,'(6f10.4)') tkch,kch(ich),eta,krm1,cph(linc)
       call coul90(krm1,eta,0d0,l,f,g,fp,gp,0,IFAIL)
       if (ifail.ne.0) then 
       write(*,*) 'coul90: ifail=',ifail; stop
       endif
       ch1  = dcmplx(g(l),f(l))  * (0.,.5)   ! outgoing (i/2) H(+)
       call coul90(krm2,eta,0d0,l,f,g,fp,gp,0,IFAIL)
       ch2 = dcmplx(g(l),f(l)) * (0.,.5) 
c lhs
       mat(ich,nch+ich)    = ch1
       mat(nch+ich,nch+ich)= ch2
c rhs
       mat(ich,2*nch+1)    =-kron(ich,inc)*conjg(ch1) 
       mat(nch+ich,2*nch+1)=-kron(ich,inc)*conjg(ch2) 

      ELSE !......................................... closed channel
       IE = 0  
       call whit(eta,r1,kch(ich),tkch,l,f1,fp1,ie)
       ch1 = f1(l+1)* (0.,.5) ! Whittaker*i/2
       call whit(eta,r2,kch(ich),tkch,l,f1,fp1,ie)
       ch2=  f1(l+1)* (0.,.5)  
c lhs
        mat(ich,nch+ich)    = -ch1
        mat(nch+ich,nch+ich)= -ch2
c rhs
        mat(ich,2*nch+1)    =0. 
        mat(nch+ich,2*nch+1)=0. 
      ENDIF ! OPEN/CLOSED CHANNELS

      enddo ! ich

      if (debug) show=.true.
      call GAUSS5(2*nch,2*nch,mat,1,sing,det,small,debug)
      if (debug) write(klog,'(5x,"Coefficients:",100f8.3)') 
     & (mat(ich,2*nch+1),ich=1,nch)

      IF(SING) then
      write(*,*)'Singular matrix for Ecm=',ecm
      DO 605 Ir=1,nr
605   IF(MAGN(Ir).NE.0.0) MAGN(Ir) = LOG10(MAGN(Ir))
      WRITE(*,610) (MAGN(Ir),Ir=1,nr)
610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
      endif

c
c S-matrix and phase-shifts
c 
      flux=0d0
      phase(1:nch)=0.
      if (show) 
!     &  write(*,'(8x,"S-matrix (wo/with velocity factors):")') 
     &   write(*,'(a40,6x,a45)')
     &   "S-matrix (WITHOUT velocity factors)", 
     &   "S-matrix (WITH velocity factors)"
      do ich=1,nch
      tkch=ecm+ech(iref)-ech(ich) 
      if (tkch.gt.0) then                    ! open channel
!commented v2.1b
!      smat(ich)=mat(nch+ich,2*nch+1)*ci**(QL(INC)-QL(ICH))
       smati(ich)=mat(nch+ich,2*nch+1)
       svel=smati(ich)*sqrt(kch(ich)/kch(inc)) ! S-matrix with velocity factors
       print*,'ich=',ich,'Ech=',ech(ich),' Open-> Svel=',svel

      if (abs(mat(nch+ich,2*nch+1)).gt.small) then
        phase(ich)=(0.,-0.5) * LOG(mat(nch+ich,2*nch+1))
        test=(0.,-0.5) * LOG(mat(nch+ich,2*nch+1)) 
      endif
      if (show) then
       if (((ich.eq.inc).and.(verb.ge.1)).or.(verb.ge.3)) then
       write(*,500) ich,smati(ich),svel,abs(svel)
       endif
      endif
      flux=flux + abs(svel)**2
500   format(10x,"Chan. #",i3,5x,
     &    "S=(",1f8.5,",",1f8.5,")", 5x, 
     &    "S=(",1f8.5,",",1f8.5,") ->  |S|=",1f9.6) 
!     &    "delta=(",1f8.3,",",1f8.3,")")

      iop=iop+1 
      ifop(inc)=.true.
      smat(inc,ich)=svel

      else                                   ! closed channel
        smati(ich)=mat(nch+ich,2*nch+1) ! Coefficient of Whittaker function
        svel=0; test=0
        if (show) write(*,501) ich,smati(ich)
501     format(10x,"Chan. #",i3,5x,
     &    "S=(",1f10.6,",",1f10.6,")", 5x, ' (closed channel)') 
        if (verb.gt.3) then
        anc=mat(nch+ich,2*nch+1)
        write(*,510) ich,anc,abs(anc)
510     format(5x,"Channel",i3,3x,
     &    "ANC=(",1g12.6,",",1g12.6,") -> |ANC|=",1g12.6)
        endif
        print*,'ich=',ich,' Closed: ANC=',anc
      smat(inc,ich)=1 ! No physical meaning; just for convenience! 
      endif !open/closed channel
      enddo !ich
      phase(1:nch)=phase(1:nch)*180/pi
!      write(45,'(1f10.3,3x,10g14.6)') ecm,
!    &  (phase(ich),ich=1,nch) 
       
      if (show) write(*,'(10x,"=> Unitarity=",1f12.8)')flux
      
!      write(96,'(50F12.5)') (smat(inc,ich),ich=1,nch)

      enddo !inc


c
c Diagonalize S-matrix to get eigenphases.
c
      eph(:)=0
!      call SEigensystem(nch,smat,nch,evec, Umat,nch, sort)

!       call dsyev('v','l',nch,smat,nch,work1,work2,3*nch,
!     &              ifail)
!*          write (6,*)
!*          write (6,*)'optimal LWORK is',work2(1)

!      DO 2 N=1,NCH
!      EPH(n) = ATAN(WORK1(NCH+1-N))
! 2    continue


      do ich=1,nch
         if(abs(evec(ich)).gt.small) eph(ich)=(0.,-0.5)*LOG(evec(ich))
         print*,'Eph=',eph(ich)
      enddo
      print 101, "Eigen:",(n,evec(n),abs(evec(n)),eph(n),n = 1, nch)
101   format(/A, 10(/"d(", I1, ") = ", F10.6, SP, F10.6, SS, " I", 
     & "|d|=",f10.6, " Phas=",f10.6))
      write(95,'(1f8.4,20f10.6)') ecm,(eph(n),n = 1, nch)


c
c scattering wavefunctions
c
      wf(1:nch,1:nr)=(0d0,0d0)
      if(abs(smati(iref)).gt.1e-20) then
        phr = (0.,-0.5) * LOG(smati(iref))
        yfac=exp(-ci*phr)
      endif

      do ich=1,nch
      do is=1,nch
      wf(ich,1:nr)=wf(ich,1:nr)+mat(is,2*nch+1)*y(ich,is,1:nr) ! gauss5
      enddo
      enddo

      end subroutine



c -----------------------------------------------------------------
c Matching to real WFS (real couplings!) 
c
c ----------------------------------------------- --------------
      subroutine match_real(ecm,z12,nch,nop,ql,lmax,inc,n1,y,wf,show)
      use nmrv, only: nr,mu,ech,h,conv,rvec
      use constants, only : e2,pi
      use globals, only: verb,debug
      implicit none
      logical :: sing,show
      integer ich,icl,is,inc,ir,klog,nd,n1,n2
      integer nch,nop,ql(nch),l,linc,lmax,ifail,ie
      integer ncl,iop,sgn,i,j
c     --------------------------------------------------------------
      real*8, dimension(0:lmax):: f,g,gp,fp
      real*8, dimension(1:lmax+1):: f1,fp1
      real*8 :: kron,r1,r2,r,whitp,raux,d1,d2,recdet,det
      real*8 :: ecm,tkch,eta,krm1,krm2,z12,kch(nch)
      real*8 :: phr,cph(0:lmax) !!!! changed to lmax check !!!!
      real*8 :: fpmax,small,acc8,flux,magn(nr)
      real*8:: yd(nch,nch),yaux(5),y1,y2
      real*8:: y(nch,nch,nr),yop(nch,nop,nr),yor(nch,nop,nr)
      real*8:: ccoef,dcoef,gaux(nch,nr),amod,rcoef,eph(nch)
      real*8:: cmat(nch-nop,nch),dmat(nch-nop,nch)
      real*8:: caux(nch-nop,nch-nop),cinv(nch-nop,nch-nop)
      real*8:: kmat(nch-nop,nop)
      real*8:: subnmat(nop,nop),nmat(nop,nop)
c     ----------------------------------------------------------------
      complex*16, parameter::iu=cmplx(0.,1.)
      complex*16:: wf(nch,nch,nr) 
      complex*16:: ch1,ch2,chp1,chp2,zaux
      complex*16:: ci,test,cj
      complex*16:: afas,cnorm,acoef
      complex*16:: anc
      complex*16:: amat(nop,nch),bmat(nop,nop)
c     ---------------------------------------------------------------
      debug=.true.
!      show=.false.

      acc8 = epsilon(acc8);
      fpmax = huge(acc8)**0.8d0 
      SMALL=1D0/FPMAX
      pi=acos(-1d0) 
      klog=99
      magn(:)=0      
      nd=6
      n2=n1-nd
      ncl=nch-nop
 
!      print*,'Match_real for Ecm=',ecm,' nop=',nop

      if (nch.eq.0) then
       write(*,*)'Matching: nch=0! Abort'; stop     
      endif

      amat(:,:)=0d0; cmat(:,:)=0; dmat(:,:)=0
      linc=ql(inc)

      do ich=1,nch
        tkch=ecm+ech(inc)-ech(ich) 
        kch(ich)=sqrt(conv*abs(tkch))
      enddo

      do is=1,nch
      do ich=1,nch
      l=ql(ich) 
      tkch=ecm+ech(inc)-ech(ich) 
!      kch(ich)=sqrt(conv*abs(tkch))
      eta=conv*z12*e2/kch(ich)/2.

      r1=rvec(n1)
      r2=rvec(n2)

      krm1=kch(ich)*r1
      krm2=kch(ich)*r2

      y1=y(ich,is,n1)
      y2=y(ich,is,n2)
   
      IF (tkch.gt.0d0) THEN   ! open channel
       call coulph(eta,cph,linc)
!       phc=exp(ci*cph(linc))
!       write(190,'(6f10.4)') tkch,kch(ich),eta,krm1,cph(linc)
       call coul90(krm1,eta,0d0,l,f,g,fp,gp,0,IFAIL)
       if (ifail.ne.0) then 
       write(*,*) 'coul90: ifail=',ifail; stop
       endif
       ch1  = dcmplx(g(l),f(l))     !  H(+)(r2)
       call coul90(krm2,eta,0d0,l,f,g,fp,gp,0,IFAIL)
       ch2 = dcmplx(g(l),f(l))      !  H(+)(r2)

!       afas=(y1*ch2-y2*ch1)/(y1*conjg(ch2)-y2*conjg(ch1))
!       amod=(y1*ch2-y2*ch1)/(y1*conjg(ch2)-y2*conjg(ch1))

!       acoef=amod*sqrt(afas)

! Alt method
       acoef=(y1*ch2-y2*ch1)/(conjg(ch2)*ch1-conjg(ch1)*ch2)
       acoef=acoef*2*iu*sqrt(kch(ich)*kch(inc))
       amat(ich,is)=acoef

!       print*,'Sol=',is,' Chan=',ich,'(open) a=',amod*sqrt(afas),acoef
!       acoef=(y2*conjg(ch1)-y1*conjg(ch2))/
!     &  (ch2*conjg(ch1)-ch1*conjg(ch2))
!     &  *2*iu*sqrt(kch(ich)*kch(inc))
!       print*,'Sol=',is,' Chan=',ich,'(open) a=',amod*sqrt(afas),acoef

      ELSE !......................................... closed channel
       IE = 0  
       call whit(eta,r1,kch(ich),tkch,l,f1,fp1,ie)
       ch1 = f1(l+1) ! Whittaker at r1
       call whit(eta,r2,kch(ich),tkch,l,f1,fp1,ie)
       ch2=  f1(l+1) ! Whittaker at r2  
 
       chp1=whitp(eta,r1,kch(ich)) ! exp[+k r1]
       chp2=whitp(eta,r2,kch(ich)) ! exp[+k r2]

       ccoef=(y1*ch2-y2*ch1)/(chp1*ch2-chp2*ch1)*kch(ich)
       dcoef=(y1*chp2-y2*chp1)/(ch1*chp2-ch2*chp1)*kch(ich)

       icl=ich-nop
       cmat(icl,is)=ccoef
       dmat(icl,is)=dcoef
       
!      print*,'Sol=',is,' Chan=',ich,' (closed) c,d=',ccoef,dcoef
       ENDIF ! OPEN/CLOSED CHANNELS      
      enddo ! ich
      enddo ! is

c ... Open channels solutions
      yop(:,1:nop,1:nr)=y(:,1:nop,1:nr)
!      do iop=1,nop
!       write(96,*)'#Open regular solution ',iop,nop,' Ecm=',ecm
!       do ir=1,nr
!       r=rvec(ir)
!       write(96,'(1f8.3,3x,50g14.6)')r,(yop(ich,iop,ir),ich=1,nch)
!       enddo !ir 
!       write(96,*)'&'
!      enddo !iop


c ... Regularize in case of closed channels
      if (ncl.gt.0) then
        cinv(1:ncl,1:ncl)=cmat(1:ncl,nop+1:nch)
        caux(:,:)=0.0
        call cmatin(cinv,caux,ncl)
! TEST      
!        caux=matmul(cinv,cmat(1:ncl,nop+1:nch))
!        write(*,*)'cinv.c=',caux 

        kmat=-matmul(cinv,cmat(1:ncl,1:nop))

c     Regular solutions
       do iop=1,nop
       gaux(:,1:nr)=0.0
       do icl=1,ncl
        gaux(:,1:nr)= gaux(:,1:nr)
     &              + kmat(icl,iop)*y(:,nop+icl,1:nr) ! CHECH!!!
       enddo ! closed channels
       yop(:,iop,1:nr)=yop(:,iop,1:nr)+gaux(:,1:nr)
      enddo ! open channels 
      endif ! open/closed channels?

c B-matrix
!      print*,'B-matrix: ncl=',ncl
      do iop=1,nop
      do ich=1,nop
       zaux=0
       do icl=1,ncl
       zaux=zaux+kmat(icl,iop)*amat(ich,nop+icl)
       enddo !icl
       bmat(ich,iop)=amat(ich,iop)+zaux
!       write(*,*)'iop,ich,b=',iop,ich,bmat(ich,iop)
      enddo !ich
      enddo !iop

c Overlap matrix
      nmat(:,:)=0
      do iop=1,nop
      do ich=1,nop
    
      nmat(iop,ich)=(sum(conjg(bmat(iop,1:nop))*bmat(ich,1:nop))
     +             +sum(bmat(iop,1:nop)*conjg(bmat(ich,1:nop))))/2.0
!      do   j=1,nop  
!      nmat(iop,ich)=nmat(iop,ich)
!     &             +0.5*(conjg(bmat(iop,j))*bmat(ich,j)
!     &             +bmat(iop,j)*conjg(bmat(ich,j)))
!      enddo ! j
!      write(*,*)'iop,ich,N=',iop,ich,nmat(iop,ich)
      enddo !ich
      enddo !iop

c orthogonalize and normalize solutions using Gram-Schmidt 
      if (1>2) then
! GS in determinant form
      rcoef=1./sqrt(nmat(1,1))
      yor(:,:,:)=0.
      yor(1:nch,1,:)=yop(1:nch,1,:)*rcoef
      sgn=+1
      do iop=2,nop
        d1=recdet(nmat(1:iop,1:iop),iop)
        d2=recdet(nmat(1:iop-1,1:iop-1),iop-1)
        do is=1,iop
        subnmat(1:iop-1,1:is-1)   =nmat(1:iop-1,1:is-1)
        subnmat(1:iop-1,is:iop-1) =nmat(1:iop-1,is+1:iop) 
        det=recdet(subnmat(1:iop-1,1:iop-1),iop-1)
        yor(1:nch,iop,:)=yor(1:nch,iop,:)+sgn*det*yop(1:nch,is,:)
        sgn=-sgn
        enddo !is
        yor(1:nch,iop,:)=yor(1:nch,iop,:)/sqrt(d1*d2)
       
      if (debug) then
      write(97,*)'# Regular solution ',iop,'of',nop,' Ecm=',ecm
      do ir=1,nr
      r=rvec(ir)
      write(97,'(1f8.3,3x,50g14.6)')r,(yor(ich,iop,ir),ich=1,nch)
      enddo !ir 
      write(97,*)'&'
      endif
      enddo !iop
 
      ELSE 
! Conventional GS form 
      rcoef=1./sqrt(nmat(1,1))
!      yor(1:nch,1:nop,1:nr)=yop(1:nch,1:nop,1:nr)
      wf(1:nch,1,:)=yop(1:nch,1,:)*rcoef
      if (debug) then
      write(97,'(a,i2,a,i2,a,1f8.3)')'# Regular solution ',1,
     & 'of',nop,' Ecm=',ecm
      do ir=1,nr
      r=rvec(ir)
      write(97,'(1f8.3,3x,50g14.6)')r,(wf(ich,1,ir),ich=1,nch)
      enddo !ir 
      write(97,*)'&'
      endif
      do iop=2,nop
        wf(1:nch,iop,:)=yop(1:nch,iop,:)
        do is=1,iop-1
        wf(:,iop,:)=wf(:,iop,:)-nmat(iop,is)*yop(:,is,:)
        enddo ! is
c norm
        raux=0
        do i=1,iop
        do j=1,iop
        ci=-nmat(iop,i)/nmat(i,i)
        cj=-nmat(iop,j)/nmat(j,j)
        if (i.eq.iop) ci=1.
        if (j.eq.iop) cj=1.
        raux=raux + ci*cj*nmat(i,j)
        enddo !j
        enddo !i
        wf(:,iop,:)=wf(:,iop,:)/sqrt(raux)
!!!!!!!!
      if (debug) then
      write(97,'(a,i2,a,i2,a,1f8.3)')'# Regular solution ',
     &  iop,'of',nop,' Ecm=',ecm
      do ir=1,nr
      r=rvec(ir)
      write(97,'(1f8.3,3x,50g14.6)')r,(wf(ich,iop,ir),ich=1,nch)
      enddo !ir 
      write(97,*)'&'
      endif
!!!!!!!!
      enddo ! iop
      endif
      end subroutine

c ... Recursive determinanat
      recursive function recdet( mat, n ) result( accum )
      integer :: n
      real*8    :: mat(n, n)
      real*8    :: submat(n-1, n-1), accum
      integer :: i, sgn

      if ( n == 1 ) then
        accum = mat(1,1)
      else
        accum = 0.0
        sgn = 1
        do i = 1, n
            submat( 1:n-1, 1:i-1 ) = mat( 2:n, 1:i-1 )
            submat( 1:n-1, i:n-1 ) = mat( 2:n, i+1:n )

            accum = accum + sgn * mat(1, i) * recdet( submat, n-1 )
            sgn = - sgn
        enddo
      endif
      end function

c
c ... Orthogonalization by QR factorization
c
      subroutine qrfac(yp,y,n,ip,irmin,nr)
      use nmrv, only: h,rvec
      implicit none
      logical    debug
      integer    nr
      complex*16 y(n,n,nr)
c ... for zgeqrf
      integer info,lda,lwork,n,m,is,ich,ip,ir,irmin
      complex*16 a(n,n),tau(n),work(2*n),rt(n,n),yp(n,n)
c ... for ztrsm
      complex*16 alpha
      integer ldb
      character diag,side,transa,uplo
!      complex*16 a(lda,*),b(ldb, *)
c ---- TEST
      integer ndim,lworkx
      complex*16 aux(3,3),taux(3),raux(3,3),workx(6)

      debug=.false.
      

      lwork=2*n

c ... QR factorization YP=Q.R 
c     (Q=orthogonal matrix; R= triangular upper matrix)
      rt=yp
      call ZGEQRF( N, N, RT, n, TAU, WORK, LWORK, INFO )
      if (info.ne.0)  then
         write(*,*)'zgeqrf failed; info=',info
      endif

      if (debug) then 
      write(*,*)'|y(old)| at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (abs(y(ich,is,ip)), is=1,min(15,n))
      enddo
      write(*,*)' '


        write(*,*)'|yp x yp| (before QR)'
        do is=1,n
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(y(:,is,ip),y(:,ich,ip)))
     &      / abs(dot_product(y(:,is,ip),y(:,is,ip)))
     &      / abs(dot_product(y(:,ich,ip),y(:,ich,ip))), 
     &     ich=1,n) 

        enddo


      write(*,*)'|R| matrix:'
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (abs(RT(ich,is)), is=1,min(15,n))
      enddo
      endif
      
c ... X*op( A ) = alpha*B
      side  ='R'
      uplo  ='U'
      transa='N'
      diag  ='N' !??????????
      m=n
      alpha =1d0
      lda   =n
      ldb   =n
      
c ... Apply transformation backwards
      do ir=irmin,ip
      yp(:,:)=y(:,:,ir)
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,YP,LDB)
      y(:,:,ir)=yp(:,:) ! override 
      enddo ! ir
      
      if (.not.debug) return

      if (debug) then
      write(*,*)'yp(new) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (abs(yp(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '


      if (debug) then
        write(*,*)'|yp x yp| (after QR)'
        do is=1,n
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(yp(:,is),yp(:,ich)))
     &      / abs(dot_product(yp(:,is),yp(:,is)))
     &      / abs(dot_product(yp(:,ich),yp(:,ich))), ich=1,n) 

        enddo
       endif !debug



 
      write(*,*)'|y-yp|/y at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (abs(y(ich,is,ip)-yp(ich,is))
     &   /abs(y(ich,is,ip)), is=1,min(15,n))
      enddo
      write(*,*)' '
      endif

c ......................... TEST ........................
      if (0>1) then
      aux(1,1)=12; aux(1,2)=-51; aux(1,3)=4
      aux(2,1)=6;  aux(2,2)=167; aux(2,3)=-68
      aux(3,1)=-4; aux(3,2)=24; aux(3,3)=-41
      lda=3
      lworkx=6
      raux=aux
      call ZGEQRF( lda, lda, raux, lda, TAUX, WORKX, LWORKX, INFO )
      if (info.ne.0)  then
         write(*,*)'zgeqrf failed; info=',info
      endif

      write(*,*)'R matrix:'
      do ich=1,3
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(raux(ich,is)), is=1,lda)
      enddo

c ... X*op( A ) = alpha*B
      side  ='R'
      uplo  ='U'
      transa='N'
      diag  ='N' !??????????
      m=n
      alpha =1d0
      lda   =3
      ldb   =3
      m     =lda
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,m,ALPHA,RAUX,LDA,AUX,LDB)

      write(*,*)'Q matrix:'
      do ich=1,3
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(aux(ich,is)), is=1,lda)
      enddo


        write(*,*)'|yp x yp| (before QR)'
        do is=1,3
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(aux(:,is),aux(:,ich)))
     &      / abs(dot_product(aux(:,is),aux(:,is)))
     &      / abs(dot_product(aux(:,ich),aux(:,ich))), ich=1,3) 

        enddo

      stop
      endif

c .......................................................
      end subroutine


c
c Gram-Smidt orthogonalization 
c
      subroutine gsorto(yp,a,nch)
      implicit none
      integer nch,ic,icp,n,i,j
      complex*16 yp(nch,nch),a(nch,nch),yaux(nch)

      a(:,:)=0;
      a(1,1)=1
      do j=2,nch      ! columns
         a(j,j)=1
         yaux(:)=0
         do i=1,j-1   ! rows
         a(i,j)=-dot_product(yp(:,i),yp(:,j))
     &          /dot_product(yp(:,i),yp(:,i))
         yaux(:)=yaux(:) +  yp(:,i)*a(i,j)
         enddo ! i
         yp(:,j)=yp(:,j)+ yaux(:)
      enddo !j
      return
      end subroutine



c
c Gram-Smidt orthogonalization 
c
      subroutine mgsorto(yp,a,nch)
      implicit none
      integer nch,ic,icp,n,i,j
      complex*16 yp(nch,nch),s(nch,nch),a(nch,nch)
!!! TEST
      character*1 jobvs,sort
      logical select
      integer lda,ldvs,sdim,lwork,info
      complex*16 w(nch),vs(nch,nch),work(2*nch)
      double precision rwork(nch)
      logical bwork(nch)
!!!!!
      a(:,:)=0; a(1,1)=1
      do j=2,nch      ! columns
         a(j,j)=1
!         yaux(:)=0
         do i=1,j-1   ! rows
         a(i,j)=-dot_product(yp(:,i),yp(:,j))
     &          /dot_product(yp(:,i),yp(:,i))
!         yaux(:)=yaux(:) +  yp(:,i)*a(i,j)
         yp(:,j)=yp(:,j)-yp(:,i)*a(i,j)
         enddo ! i
!         yp(:,j)=yp(:,j)+ yaux(:)
      enddo !j

      RETURN
c overlap matrix 
      do ic=1,nch
      do icp=1,nch
         s(ic,icp)=dot_product(yp(:,ic),yp(:,icp))
      enddo ! icp
      enddo ! ic


      lwork=2*nch
      lda=nch
      ldvs=nch
      jobvs='N'
      sort='N'
      call ZGEES(JOBVS, SORT, SELECT, NCH, A, LDA, SDIM, W, VS,
     $                  LDVS, WORK, LWORK, RWORK, BWORK, INFO )
*
      if (info.eq.0) then
      write(*,'(5(i3,2x,1g12.6,3x))') (n,w(n),n=1,nch)
      else
         write(*,*)'error: info=',info
      endif
            
      end subroutine
 
      function select(a)
      complex*16 a
      logical select
      end function 




c----------------------------------------------------------------------
c calculates Coulomb phase-shifts (by C. Bertulani)
c----------------------------------------------------------------------
      subroutine coulph(eta,cph,lmax)
      implicit none
      real*8 eta, sto,fi
      integer lmax,ii
!      implicit real*8(a-h,o-z)
      real*8 cph(0:lmax)

!      include 'crossec.dim'
!      dimension cph(0:lmax)
c
      sto=16.d0+eta*eta
      cph(0)=-eta+(eta/2.d0)*dlog(sto)+3.5d0*datan(eta/4.d0)-(
     1 datan(eta)+datan(eta/2.d0)+datan(eta/3.d0))-(eta/(12.d0*
     2 sto))*(1.d0+(1.d0/30.d0)*(eta**2-48.d0)/sto**2+(1.d0/105.d0)
     3 *(eta**4-160.d0*eta**2+1280.d0)/sto**4)
c
      do 1 ii=1,lmax
         fi=dble(ii)
         cph(ii)=cph(ii-1)+datan(eta/fi)
1     continue
      return
      end

c -------------------------------------------------------
c   Solve multichannel Schrodinger equation
c   ecm     = kinetic energy
c           (E= Ecm + core energy)
c   nch     = number of channels
c   inc     = incoming channel
c   ql(nch) = array with l'values for nch channels
c             (NEEDS TO BE ALLOCATED BY USER)
c   nr      = number of radial points
c   dr      = radial step 
c   z1*z2   = product of charges
c   wf      = array to store wfs(nch,nr) (out)
c   phase   = array to store phase shifts (out)
c   factor  = conversion factor 2*mu/hbar^2
c   
c   The following variables must be passed through the nmrv module:
c    ech(nch) = channel energies (NEEDS TO BE ALLOCATED BY USER)
C    vcoup   = matrix with couplings (nchan,nchan,nr)
c   ----------------------------------------------------
c   ym is the component of the solution at r-h (similarly for Z)
c   y0 is the component of the solution at r
c   yp                                     r+h
c ------------------------------------------------------
      subroutine schcc_MGR(nch,ecm,z12,incvec,ql,factor,dr,r0,
     & npt,wf,phase,smat,info,einc,icc)
      use xcdcc, only: smats
      use nmrv,only: nr,h,vcoup,conv,ech,debug,rmin,hort,
     &               rvec,rmort
      use globals, only: verb
      use factorials
      use constants , only: e2
      use memory
      implicit none
      logical:: copen(nch),orto,info
      integer:: nch,klog,npt,ql(nch)
      integer:: ir,irmin,is,isp,ich,inc,l,lmax
      integer:: i,j,k,im,ip,i0,norto,icc
c     ----------------------------------------------
      real*8 :: ymem,ti,tf
      real*8 :: ecm,tkch,eta,z12
      real*8 :: dr,r,r0,rm,h2,factor
      real*8 :: rturn,rmorto,rmatch
      real*8 :: kch(nch),kch2(nch),rturnv(nch),etav(nch)
      real*8 :: start,finish,aux,big,small
      real*8 , parameter:: pi=acos(-1d0)
      real*8 :: einc,kinc
c     ----------------------------------------------
      complex*16, intent(out):: wf(nch,nr)
      complex*16, dimension(nch,nch):: yp,yc,y0,ym
      complex*16, dimension(nch,nch):: zp,zc,z0,zm
      complex*16 :: y(nch,nch,npt)
      complex*16 :: phase(nch),smat(nch)
      complex*16 :: a(nch,nch),c1,c2
      logical :: incvec(nch)
c     ---------------------------------------------
      y=0d0      
      debug=.false.

      big=huge(big)
      small=epsilon(small)
      h=dr
      h2=h*h
      conv=factor
      copen(:)=.true.
      nr=npt
      klog=99
      rmin=0
      lmax=maxval(ql)
      orto=.false.

      if (debug) then
        write(*,*)'schcc: ql=',ql(1:nch)
        write(*,*)'schcc: ecm=',ecm
      endif

      if (allocated(rvec)) deallocate(rvec)
      allocate(rvec(nr))
      do ir=1,nr
        rvec(ir)=r0+(ir-1)*dr
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius for derivative

      norto=nint(hort/h)
      if (abs(hort).gt.0.) orto=.true.

c     ------------------------------------------------
      if (verb.ge.3) then
        ymem=nch*nch*(npt+2)*lc16/1e6
        write(*,190) ymem
190     format(5x,"[ WFS require", 1f8.2," Mbytes ]")
        ymem=nch*nch*lc16/1e6
!        write(*,*)'Auxiliary arrays require 8x',
!     &  ymem,'=',ymem*8, 'Mbyes'
!        write(*,*)'Max real=',big
      endif 
c     ---------------------------------------------------
      if (ecm.lt.0) then 
      write(*,200) ecm
200   format(4x,'Skipping negative energy:',1f7.3, ' MeV')
      wf(:,:)=0d0
      return
      endif

      if (rmort.gt.0) rmorto=rmort
      
      kinc=sqrt(conv*ecm)   
      do ich=1,nch
      aux=ecm+einc-ech(ich)
      kch2(ich)=conv*aux
      if (aux.gt.0) then
        if (debug) write(99,300) ich,aux,"open"
        kch(ich)=sqrt(conv*aux)
        etav(:)=conv*z12*e2/kch(:)/2.
        rturnv(ich) =(etav(ich)+SQRT(etav(ich)**2 + 
     &               ql(ich)*(ql(ich)+1d0)))/kch(ich)
        if (rturnv(ich).gt.rmorto) rmorto=rturnv(ich) 
        copen(ich)=.true.
      else
        if (debug) write(99,300) ich,aux,"closed"
        kch(ich)=sqrt(-conv*aux)
        rmorto=rm
        copen(ich)=.false.
      endif
300   format(3x,"Channel",i2," Final Ecm=",1f8.3," MeV ->",a6)
      enddo
      if (rmort.lt.0) rmorto=rmorto+abs(rmort) 



c classical turning point for incoming channel
      eta=conv*z12*e2/kinc/2.
      l=1000
      do ich=1,nch
      if (incvec(ich)) then
      l=min(ql(ich),l)
      endif
      enddo
      RTURN =(eta+SQRT(eta**2 + L*(L+1d0)))/kinc  
c turning point for all channels
!      rmorto=maxval(rturnv(:))  ! outermost turning point

      if (debug) write(*,'(6(5x,i3,2x,1f8.3))')
     &   (ich,rturnv(ich),ich=1,nch) 
 

      call factorialgen(lmax+1)
      call cpu_time(start)
c 
c Starting values and initial integration point
c 
!      rmin=h ;       irmin=1
      rmin=max(h,rturn-10.) ! CHECK h or r0 ?
      irmin=nint(rmin/h)+1

!      write(*,*) 'PC-Numerov  from ir=',irmin,' r=',rmin,' fm',
!    & ' to  R=',Rm,'in steps of',h,' and matching at Rmatch=',rmatch


      if (debug) write(95,*)'Starting value ir=,',irmin,'r=',rvec(irmin)
      ym(:,:)=0d0; zm(:,:)=0d0; 
      y0(:,:)=0d0; z0(:,:)=0d0;
      do ich=1,nch
        l=ql(ich)
!        if (debug) write(klog,*)'ich,l,dlfac(l)',ich,l,dlfac2(l)
 
50      aux= (kch(ich)*h)**(l+1)/exp(0.5*dlfac2(l))
!        if (aux.lt.1e-30) then
!          rmin=rmin+5*h
!          goto 50        
!        endif

!        irmin=nint(rmin/h) !+1
!        y0(ich,ich)=(kch(ich)*h)**(l+1)/exp(0.5*dlfac2(l))
!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        y0(ich,ich)=1e-10
        y0(ich,ich)=1e-10
!!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
        if (debug) then
        write(95,320) (y0(ich,is),is=1,min(10,nch)) 
320     format(5x,50("(",2g16.6,")",5x))
        endif
      enddo

!      if (debug) write(95,*)'Initial dot product(ABS):'
!      write(95,'(10(2x,"is=",i2,3x,1g10.5))') (is,
!     &     abs(dot_product(y0(:,is),y0(:,is-1)))
!     & , is=2,nch)


      if (info) then
      write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
     &  5x,"outer=",1f8.2," fm",
     &  5x,"R-min=",1f6.2," fm")') 
     &  minval(rturnv(:)), maxval(rturnv(:)), rmin      
      endif
c
c Zm=F[r=0,ym] and Z0=F[r=h,y0]
c
      if (irmin.lt.2) write(*,*)'wrong irmin=',irmin
      call zfun(nch,ql,kch2,irmin-1,ym,zm)
      call zfun(nch,ql,kch2,irmin  ,y0,z0)
      
      do ir=irmin,nr-1  
      r=rvec(ir) ! (ir-1)*h    
      im=ir-1  ! r-h
      i0=ir    ! r
      ip=ir+1  ! r+h
      call zfun(nch,ql,kch2,im,ym,zm)
!      call zfun(nch,ql,kch2,i0,y0,z0)

      if (debug.and.(mod(ir,10).eq.0)) then
         write(95,'(/,"y0(x) ir=",i4,3x,"r=",1f7.2)') ir,r
         do ich=1,min(nch,10)
      write(95,'(5x,50g14.5)')(y0(ich,is),is=1,min(nch,10))
 !        write(*,320) (y0(ich,is),is=1,min(nch,10))
         enddo
      endif !debug

c
c predictor
c                         
      YP(:,:)=2*Y0(:,:)-YM(:,:)+h2*Z0(:,:)  
 
      if (maxval(abs(yp)).gt.big) then
         write(*,*)' ** WARNING ** Solution too big at ir,r=',ir,r
      endif
      if (debug.and.(mod(ir,10).eq.0)) then
        write(95,*)'Predictor for ir,r=',ir,r
        do ich=1,min(nch,10)
!        write(95,'(5x,50g14.5)') (yp(ich,is),is=1,nch)
        write(95,320) (yp(ich,is),is=1,min(nch,10))
        enddo
        call flush(95)
      endif !debug
                         
c
c corrector
c         
      call zfun(nch,ql,kch2,ip,yp,zp)
      yc(:,:)=yp(:,:)+(h2/12.)*(zp(:,:)+zm(:,:)-2*z0(:,:))

      if (debug.and.(mod(ir,10).eq.0)) then
        write(95,*)'Corrector:'
        do ich=1,min(nch,10)
        write(95,'(5x,50g14.5)') (yc(ich,is),is=1,nch)
        enddo
      endif ! debug

c 
c modify final solution 
c
      call zfun(nch,ql,kch2,ip,yc,zc)
      yp(:,:)=yc(:,:)+(h2/30.)*(zc(:,:)-zp(:,:))
      call zfun(nch,ql,kch2,ip,yp,zp)

      if (debug.and.(mod(ir,10).eq.0)) then
        write(95,*)'Modified corrector:'
        do ich=1,min(nch,10)
        write(95,'(5x,50g14.5)') (yp(ich,is),is=1,nch)
        enddo
      
        write(95,*)'Max|Yc-Yp|:',maxval(abs(yc-yp))
!        do ich=1,nch
!        write(94,'(5x,50g14.5)') 
!     &  (abs(yc(ich,is)-yp(ich,is)),is=1,nch)
!        enddo
        call flush(95)

      endif ! debug

c store solution
      y(1:nch,1:nch,ip)=yp(1:nch,1:nch) 


c Re-orthogonalize solution vectors
!      orto =.true.
!      if (orto) write(*,*)ir,ir-irmin,norto,r,rmorto,mod(ir-irmin,norto)
      if (orto.and.(r.lt.rmorto).and.(ir.gt.irmin)
     & .and.(mod(ir-irmin,norto).eq.0)) then 
!      if (orto.and.(r.le.rmorto).and.(mod(ir-irmin,norto).eq.0)) then 
        write(*,'(5x,"-> orthogonalizing at r=",1f7.2)')r

        if (debug) then
        write(94,*)'|yp x yp| (before G-S) for ir,r',ir,r
        do is=2,nch
        write(94,'(5x,50g14.5)') 
     &    ( abs(dot_product(yp(:,is),yp(:,ich)))
     &      / abs(dot_product(yp(:,is),yp(:,is)))
     &      / abs(dot_product(yp(:,ich),yp(:,ich))), ich=1,nch) 
!         write(94,'(5x,i2,"x",i2,3x,50g14.5)') is,is-1,
!     &    abs(dot_product(yp(:,is),yp(:,is-1)))
        enddo
        endif !debug

  
!      call gsorto(yp,a,nch)  ! standard Gram-Schmidt
!      call mgsorto(yp,a,nch)  ! modified Gram-Schmidt
!      do i=irmin,ir+1
!      y(:,:,i)=matmul(y(:,:,i),a(:,:))
!      enddo ! i

!       write(*,*)'calling orto for ir,r=',ir,r
       call cpu_time(ti)
       call qrfac(yp,y,nch,ip,irmin,nr)      ! QR factorization
       call cpu_time(tf)
       torto=torto+tf-ti

      if (debug) then
      write(94,*)'Ortog. matrix at ir=',ir
      do ich=1,nch
        write(94,'(5x,50g14.5)') 
     &  (a(ich,is),is=1,nch)
           write(94,'(5x,50g14.5)') 
     &  (a(ich,is),is=1,nch)
      enddo
      call flush(94)
      endif 

!      y(:,:,ip)=yp

      if (debug) then
      write(90,*)'y at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(90,'(5x,i3,50g14.5)') ich,
     &  (abs(y(ich,is,ip)), is=1,min(15,nch))
      enddo
      write(90,*)' '

      write(90,*)'yp at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(90,'(5x,i3,50g14.5)') ich,
     &  (abs(yp(ich,is)), is=1,min(15,nch))
      enddo
      write(*,*)' '

      write(90,*)'|y-yp|/y at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(90,'(5x,i3,50g14.5)') ich,
     &  (abs(y(ich,is,ip)-yp(ich,is))
     &   /abs(y(ich,is,ip)), is=1,min(15,nch))
      enddo
      write(*,*)' '
      endif

      yp=y(:,:,ip);  call zfun(nch,ql,kch2,ip,yp,zp)
      y0=y(:,:,i0);  call zfun(nch,ql,kch2,ip,y0,z0)


c dot products after orthogonalization
      if (debug) then
        write(94,*)'|yp x yp| (after G-S)'
        do is=1,nch
!        write(95,'(5x,i2,"x",i2,3x,50g14.5)') 
        write(94,'(5x,50g14.5)') 
     &    ( abs(dot_product(yp(:,is),yp(:,ich)))
     &      / abs(dot_product(yp(:,is),yp(:,is)))
     &      / abs(dot_product(yp(:,ich),yp(:,ich))), ich=1,nch) 

        enddo
       endif !debug

!       y(1:nch,1:nch,ir+1)=yp(1:nch,1:nch)  

      if (debug) then
        write(94,*)'Gram-Schmidt ortogonalized:'
        do ich=1,nch
        write(94,'(5x,50g14.5)') (yp(ich,is),is=1,nch)
        enddo
        write(94,*)' '
        call flush(94)
      endif !debug


      endif !orto
! ------------------------------------------------------------


c Alternative Gramm-Schmitt orthogonalization of solution vectors
!      orto=.false.
!      rmorto=maxval(rturnv(:))
!      if (orto.and.(r.lt.rmorto).and.(mod(ir,5).eq.0)) then 
!         call gsorto(yp,a,nch)
!      end if 

!!!!!!!!!!!!!!!!!!!!!!!!! TEST
      if (debug) then
!       write(91,*) r,maxval(abs(y(:,:,ip)))
       write(91,'(5x,1f6.3,3x,50g14.5)') r,(abs(y(1,ich,ip)),ich=1,nch)
!       y(:,:,:)=y(:,:,:)/1.e5
!       yp=y(:,:,ip)
!       y0=y(:,:,i0)
!       call zfun(nch,ql,kch2,ip,yp,zp)
!       call zfun(nch,ql,kch2,ip,y0,z0)
      endif       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c
c ym <- y0 <- yp 
c
!      ym(:,:)=  y(1:nch,1:nch,ir)
!      zm(:,:)
      ym(:,:)=y0(:,:) 
      zm(:,:)=z0(:,:)
      y0(:,:)=yp(:,:) 
      z0(:,:)=zp(:,:)


      enddo !ir
      call cpu_time(finish)

!!! TEST
!      if (maxval(abs(y(:,:,nr))).gt.1e20) then
!        y(:,:,:)=y(:,:,:)/1.e20
!      endif
!!! TEST

      
c 
c Match with asymptotic solution (derivative)
c Using derivative and gauss5
       call cpu_time(start)
!      call matching(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! gauss5
c Adjacent points
      do ich=1,nch
      if (incvec(ich)) then
      inc=ich
      call matching3(ecm,z12,nch,ql,lmax,inc,nr,y,wf,phase,smat,info) ! gauss5
      endif
            smats(icc,inc,1:nch)=smat(1:nch) 
      enddo
!      call matching2(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! LAPACK 

       call cpu_time(finish)
       tmatch=tmatch+finish-start
      if (debug) write(klog,'(2f10.3,3x,10g14.6)') ecm,z12,ql(1),
     &  (phase(ich),ich=1,nch) 


      call flush(6)
      end subroutine schcc_MGR

!-----------------------------------------------------------------------

c ...
c ... Enhanced Numerov (version of IJ Thompson used in Fresco) 
c ... 
      subroutine schcc_erwin_MGR(nch,ecm,z12,incvec,ql,factor,dr,r0,
     & npt,wf,phase,smat,method,info,einc,icc)
      use xcdcc, only: smats
      use nmrv,only: nr,h,vcoup,conv,ech,debug,rmin,hort,rvec,cutr
      use globals, only: verb
      use factorials
      use constants , only: e2
      use memory
      implicit none
      integer:: method
      logical:: copen(nch),orto,raynal,info
      integer:: nch,klog,npt,ql(nch)
      integer:: ir,irmin,is,isp,ich,inc,l,lmax
      integer:: i,j,k,im,ip,i0,norto,icc
c     .........................................................
      real*8     :: ymem,ti,tf
      real*8 , parameter:: pi=acos(-1d0)
      real*8 :: ecm,tkch,eta,z12,r12,ena2,ena3
      real*8 :: dr,r,r0,rm,h2,factor,einc,kinc
      real*8 :: rturn,rmorto,rmatch
      real*8 :: kch(nch),kch2(nch),rturnv(nch),etav(nch)
      real*8 :: start,finish,aux,big,small,tmin,tmax
c     .........................................................
      complex*16, intent(out):: wf(nch,nr)
      complex*16, dimension(nch,nch):: y0
      complex*16, dimension(nch,nch):: z0,zm
      complex*16 :: y(nch,nch,npt),coupl(nch,nch),v(nch,nch)
      complex*16 :: phase(nch),smat(nch)
      complex*16 :: s(nch),ziv,zpv
      complex*16 :: a(nch,nch)
      complex*16 :: c1,c2,c
      complex*16,parameter:: onec=(1D0,0D0),zero=(0d0,0d0)
c ... for QR factorization
      complex*16 :: HDIAG(nch),w0(nch,nch),wm(nch,nch)
c ... TEST (delete after debugging)
      complex*16:: yaux(nch,nch)
      logical incvec(nch)
c     ---------------------------------------------------------
      R12 = 1D0/12D0
      ENA2 = 2D0/5D0 * R12**2
      ENA3 = - 4D0/35D0 * R12**3
!      ONEC = (1D0,0D0)
      TMAX = 20.
      TMIN = -125.
      big=huge(big)
      small=epsilon(small)
     
      h=dr
      h2=h*h
      conv=factor
      copen(:)=.true.
      nr=npt
      klog=99
      rmin=0
      lmax=maxval(ql)
      orto=.false.
c ... Initialize some variables 
      y=0d0      
      coupl=0
c ... Set TRUE for debugging
      debug=.false.
c ..............................................................................
      if (nch.lt.1) then 
        write(*,*)'Scatcc_erwin: nch=',nch; stop
      endif

      if (ecm.lt.1e-6) then
        write(*,*)'Scatcc_erwin: Ecm too small! =',ecm
        stop
      endif 
 
      if (allocated(rvec)) deallocate(rvec)
      allocate(rvec(nr))

      do ir=1,nr
        rvec(ir)=r0+(ir-1)*dr
!        write(222,*)rvec(ir),real(vcoup(1,1,ir)),ir
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius if derivative is used
     
   
      norto=nint(hort/h)
      if (abs(hort).gt.0.) orto=.true.

      if (info) then
        ymem=nch*nch*(npt+2)*lc16/1e6
        write(*,190) ymem
190     format(5x,"[ WFS require", 1f8.2," Mbytes ]")
        ymem=nch*nch*lc16/1e6
!        write(*,*)'Auxiliary arrays require 8x',
!     &  ymem,'=',ymem*8, 'Mbyes'
!        write(*,*)' Biggest real=',big
      endif 

      if (ecm.lt.0) then 
      write(*,200) ecm
200   format(4x,'Skipping negative energy:',1f7.3, ' MeV')
      wf(:,:)=0d0
      return
      endif

      rmorto=0
      kinc=sqrt(conv*ecm)
!      write(0,*) 'einc',einc,'kinc',kinc,'incvec',incvec(:)
      do ich=1,nch
      aux=ecm+einc-ech(ich)
      kch2(ich)=conv*aux
!      write(223,*)'ich,E,kch2=',ich,aux,kch2(ich)

      if (aux.gt.0) then
        copen(ich)=.true.
        if (debug) write(99,300) ich,aux,"open"
        kch(ich)=sqrt(conv*aux)
        etav(ich)=conv*z12*e2/kch(ich)/2.
        rturnv(ich) =(etav(ich)+SQRT(etav(ich)**2 
     &              + ql(ich)*(ql(ich)+1d0)))/kch(ich) 
        if(rturnv(ich).gt.rmorto) rmorto=rturnv(ich)
      else
        copen(ich)=.false.
        if (debug) write(99,300) ich,aux,"closed"
        kch(ich)  =sqrt(-conv*aux)
        rturnv(ich)=rm
        rmorto     =rm
      endif
300   format(3x,"Channel",i2," Final Ecm=",1f8.3," MeV ->",a6)
      enddo


c classical turning point for incoming channel
      eta=conv*z12*e2/kinc/2.
      l=1000
      do ich=1,nch
      if (incvec(ich)) then
      l=min(ql(ich),l)
      endif
      enddo
      RTURN =(eta+SQRT(eta**2 + L*(L+1d0)))/kinc  
!      write(0,*) 'eta',eta,'l',l,'rturn',rturn

      if (1>2) write(*,'(6(5x,i3,2x,1f8.3))')
     &   (ich,rturnv(ich),ich=1,nch) 

 
      call factorialgen(lmax+1)
      call cpu_time(start)
c 
c Starting values and initial integration point
c 
!      rmin=h ;       irmin=1
!      rmin=max(h,rturn-10.) ! CHECK h or r0 ?
      rmin  =max(h,minval(rturnv(:)-abs(cutr)))
      irmin=nint(rmin/h)+1

!      write(*,*) 'Enhanced Numerov from ir=',irmin,' r=',rmin,' fm',
!     & ' to  R=',Rm,'in steps of',h,' and Rmax==',rm


      if (debug) write(95,*)'Starting value ir=,',irmin,'r=',rvec(irmin)
      z0(:,:)=zero; zm(:,:)=zero; 
      w0(:,:)=zero; wm(:,:)=zero
      y0(:,:)=zero; v(:,:)=zero 
      do ich=1,nch
        l=ql(ich)
!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ZI(L,L) = H**(QNF(9,L)+1) / EXP(0.5 * FACT(QNF(9,L)+1))
        z0(ich,ich)= 1e-10  !TEEEEEEEEEEEEEEEEEEEEEESSSSST!        
!        z0(ich,ich)= H**(l+1.) / EXP(0.5 * FACT(l+1.))  
!         write(*,*)'ich,z0=:',  ich,z0(ich,ich)   
!!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
      enddo
      
      if (info) then
      write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
     &  5x,"outer=",1f8.2," fm",
     &  5x,"R-min=",1f6.2," fm")') 
     &  minval(rturnv(:)), maxval(rturnv(:)), rmin      
      endif

      if (irmin.lt.2) write(*,*)'** ERWIN: wrong irmin=',irmin
 
c Start radial integration
      do ir=irmin,nr
      r=rvec(ir) ! (ir-1)*h    
      im=ir-1  ! r-h
      i0=ir    ! r
      ip=ir+1  ! r+h
      if (r.lt.1e-6) r=1e-6

c kinetic energy part
c Fresco: SMAT(ich) = -LL1(K)*RI2 + (-ECM(ich) + C) * H
      do ich=1,nch    
      l=ql(ich)
      s(ich)= h2*(-l*(l+1)/r**2         ! centrifugal
     &          + kch2(ich)             ! energy
     &        - conv*vcoup(ich,ich,ir)) ! potential
      enddo
    
      do is=1,nch
      do ich=1,nch      
      y0(ich,is)= z0(ich,is)*
     X      (ONEC - S(ich) * (R12 - S(ich)*(ENA2 + S(ich)*ENA3)))
      enddo
      enddo

      w0(:,:)=zero
      do is=1,nch
      w0(is,is)=(ONEC - S(is) * (R12 - S(is)*(ENA2 + S(is)*ENA3)))
      enddo 
!      y0=matmul(w0,z0)

c non-diagonal part of V matrix 
      coupl(:,:)=-h2*conv*vcoup(:,:,ir)  
      do ich=1,nch
         coupl(ich,ich)=0d0
      enddo 

      w0=w0-r12*coupl
      
      if (2>1) then
      DO 24 ICH=1,NCH
	 DO 24 IS=1,NCH
        C = COUPL(is,ich) * R12
	 if(C/=ZERO) y0(ich,1:nch) = y0(ich,1:nch) - C * z0(is,1:nch)
24    CONTINUE
        DO 34 ich=1,nch ! k
        V(ich,:)  = ZERO 
        DO 34 is=1,nch   ! j
	 C = COUPL(is,ich)
 	 if(C/=ZERO) V(ich,1:nch) = V(ich,1:nch) + C * y0(is,1:nch)
34    CONTINUE
      else  !..........Idem using builtin matrix multiplication 
        V(1:nch,1:nch)  = ZERO 
        y0=y0  - r12*matmul(coupl,z0)
        v =   matmul(coupl,y0)
      endif

!! TEST 
      if (debug) then
      yaux=matmul(w0,z0)
      write(*,*)'y0 at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (y0(ich,is), is=1,min(15,nch))
      enddo

!      write(*,*)'w0*z0 at ir,r=',ir,r
!      do ich=1,min(15,nch)
!        write(*,'(5x,i3,50g14.5)') ich,
!     &  (yaux(ich,is), is=1,min(15,nch))
!      enddo
      endif

!!!!!!!!!!!!!!!!!!!!!!1


      if (debug) then
!      write(*,*)'coupl at ir,r=',ir,r
!      do ich=1,min(15,nch)
!        write(*,'(5x,i3,50g14.5)') ich,
!     &  (coupl(ich,is), is=1,min(15,nch))
!      enddo

!      write(*,*)'V at ir,r=',ir,r
!      do ich=1,min(15,nch)
!        write(*,'(5x,i3,50g14.5)') ich,
!     &  (V(ich,is), is=1,min(15,nch))
!      enddo

      write(*,*)'z0 at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (z0(ich,is), is=1,min(15,nch))
      enddo

      endif ! debug


c store solution 
      y(:,:,ir)=y0(:,:)       


c zm <- z0, wm <-w0
      do is=1,nch
      do ich=1,nch
      ziv=z0(ich,is)
      zpv=2.0*ziv - zm(ich,is) - v(ich,is) - s(ich)*y0(ich,is)
      zm(ich,is)=ziv
      z0(ich,is)=zpv
      wm(ich,is)=w0(ich,is)
      enddo  !ich
      enddo  !is


c Re-orthogonalize solution vectors ....................................
      if (orto.and.(r.lt.rmorto).and.(ir.gt.irmin).and.(hort.gt.0)
     & .and.(mod(ir-irmin,norto).eq.0).and.(ir.lt.nr-10)) then 
        write(*,'(5x,"->orthogonalizing at ir,r=",i4,1f7.2)')ir,r
  
      if(1>2) then
c check orthogonality
      a=matmul(conjg(transpose(y0)),y0)
      write(*,*)'Unitarity (before) at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (a(ich,is)/a(1,1), is=1,min(15,nch))
      enddo
      endif

      if (1>2) then 
      write(*,*)'y0 (before) at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (y0(ich,is), is=1,min(15,nch))
      enddo
      endif

       call cpu_time(ti)
       call qrerwinz(y0,z0,zm,w0,nch,ir)      ! QR factorization

       y(:,:,ir)=y0(:,:) ! new (orthogonolized) solution at ir

      if (debug) then
      write(*,*)'y0 (after) at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (y0(ich,is), is=1,min(15,nch))
      enddo

        write(*,*)'|yp x yp| (after QR)'
        do is=1,nch
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(y0(:,is),y0(:,ich)))
     &      / abs(dot_product(y0(:,is),y0(:,is)))
     &      / abs(dot_product(y0(:,ich),y0(:,ich))), ich=1,nch) 
        enddo
      endif

       call cpu_time(tf)
       torto=torto+tf-ti

      endif !orto
!........................................................................


!!!!!!!!!!!!!!!!!!!!!!!!! TEST
      if (debug) then
       write(91,'(5x,1f6.3,3x,50g14.5)') r,(abs(y(1,ich,ir)),ich=1,nch)
      endif       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         DO 46 IT=1,NEQS
!         DO 43 K=1,NEQS
!            ZIV = ZI(IT,K)
!            ZPV = ZIV + ZIV - ZM(IT,K) - V(IT,K) - SMAT(K) * FI(IT,K)
!            ZM(IT,K) = ZIV
!            ZI(IT,K) = ZPV
!43            CONTINUE
!46    CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      enddo !ir
      call cpu_time(finish)


c 
c Match y with asymptotic solution -> S-matrix
c            
      call cpu_time(start)
!      call matching(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! gauss5 + derivative
!      call matching2(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! LAPACK 
      do ich=1,nch
      if (incvec(ich)) then
      inc=ich
      call matching3(ecm,z12,nch,ql,lmax,inc,nr,y,wf,phase,smat,info) ! gauss5 + 2 points
      endif
            smats(icc,inc,1:nch)=smat(1:nch) 
      enddo
       call cpu_time(finish)
       tmatch=tmatch+finish-start

      if (debug) write(klog,'(2f10.3,3x,10g14.6)') ecm,z12,ql(1),
     &  (phase(ich),ich=1,nch) 


      call flush(6)
      end subroutine schcc_erwin_MGR


!-----------------------------------------------------------------------

c -------------------------------------------------------
c   Solve multichannel Schrodinger equation
c   ecm     = kinetic energy
c           (E= Ecm + core energy)
c   nch     = number of channels
c   inc     = incoming channel
c   ql(nch) = array with l'values for nch channels
c             (NEEDS TO BE ALLOCATED BY USER)
c   nr      = number of radial points
c   dr      = radial step 
c   z1*z2   = product of charges
c   wf      = array to store wfs(nch,nr) (out)
c   phase   = array to store phase shifts (out)
c   factor  = conversion factor 2*mu/hbar^2
c   
c   The following variables must be passed through the nmrv module:
c    ech(nch) = channel energies (NEEDS TO BE ALLOCATED BY USER)
C    vcoup   = matrix with couplings (nchan,nchan,nr)
c   ----------------------------------------------------
c   ym is the component of the solution at r-h (similarly for Z)
c   y0 is the component of the solution at r
c   yp                                     r+h
c ------------------------------------------------------
      subroutine schcc_ena_MGR(nch,ecm,z12,incvec,ql,factor,dr,r0,
     & npt,wf,phase,smat,method,info,einc,icc)
      use xcdcc, only: smats
      use nmrv,only: nr,h,vcoup,conv,ech,debug,rmin,hort,rvec,rmort
      use globals, only: verb
      use factorials
      use constants , only: e2
      use memory
      implicit none
      integer:: method
      logical:: copen(nch),orto,raynal,info
      integer:: nch,klog,npt,ql(nch)
      integer:: ir,irmin,is,isp,ich,inc,l,lmax
      integer:: i,j,k,im,ip,i0,norto,icc

      real*8 :: ecm,tkch,eta,z12,einc,kinc
      real*8 :: dr,r,r0,rm,h2,factor
      real*8 :: rturn,rmorto,rmatch
      real*8 :: kch(nch),kch2(nch),rturnv(nch),etav(nch)
      real*8 :: start,finish,aux

      complex*16, intent(out):: wf(nch,nr)
      complex*16, dimension(nch,nch):: yp,yc,y0,ym
      complex*16, dimension(nch,nch):: zp,zc,z0,zm
      complex*16, dimension(nch,nch):: wp,w0,wm

      complex*16 :: y(nch,nch,npt),w(nch,nch,npt)
      complex*16 :: phase(nch),smat(nch)
      real*8 , parameter:: pi=acos(-1d0)
      complex*16 :: a(nch,nch),c1,c2
      real*8     :: ymem,ti,tf
      logical incvec(nch)
c     -------------------------------------------------------
      y=0d0      
      debug=.false.
      if (method.eq.3) raynal=.true.
      h=dr
      h2=h*h
      conv=factor
      copen(:)=.true.
      nr=npt
      klog=99
      rmin=0
      lmax=maxval(ql)
      orto=.false.

      
      if (allocated(rvec)) deallocate(rvec)
      allocate(rvec(nr))
      do ir=1,nr
        rvec(ir)=r0+(ir-1)*dr
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius

      norto=nint(hort/h)
      if (abs(hort).gt.0.) orto=.true.

      if (verb.ge.3) then
        ymem=nch*nch*(npt+2)*lc16/1e6
        write(*,190) ymem
190     format(5x,"[ WFS require", 1f8.2," Mbytes ]")
        ymem=nch*nch*lc16/1e6
!        write(*,*)'Auxiliary arrays require 8x',
!     &  ymem,'=',ymem*8, 'Mbyes'
      endif 

      if (ecm.lt.0) then 
      write(*,200) ecm
200   format(4x,'Skipping negative energy:',1f7.3, ' MeV')
      wf(:,:)=0d0
      return
      endif

      rmorto=rmort
      kinc=sqrt(conv*ecm)
      do ich=1,nch
      aux=ecm+einc-ech(ich)
      kch2(ich)=conv*aux
      if (aux.gt.0) then
        if (debug) write(99,300) ich,aux,"open"
        kch(ich)=sqrt(conv*aux)
        etav(ich)=conv*z12*e2/kch(ich)/2.
        rturnv(ich) =(etav(ich)+SQRT(etav(ich)**2 + 
     &              ql(ich)*(ql(ich)+1d0)))/kch(ich) 
        if(rturnv(ich).gt.rmorto) rmorto=rturnv(ich) 
        copen(ich)=.true.
      else
        if (debug) write(99,300) ich,aux,"closed"
        kch(ich)=sqrt(-conv*aux)
        rmorto=rm
        copen(ich)=.false.
      endif
300   format(3x,"Channel",i2," Final Ecm=",1f8.3," MeV ->",a6)

      enddo


c classical turning point for incoming channel
      eta=conv*z12*e2/kinc/2.
      l=1000
      do ich=1,nch
      if (incvec(ich)) then
      l=min(ql(ich),l)
      endif
      enddo
      RTURN =(eta+SQRT(eta**2 + L*(L+1d0)))/kinc  

c turning point for each channel

      if (debug) write(*,'(6(5x,i3,2x,1f8.3))')
     &   (ich,rturnv(ich),ich=1,nch) 
 

      call factorialgen(lmax+1)
      call cpu_time(start)
c 
c Starting values and initial integration point
c 
!      rmin=h ;       irmin=1
      rmin=max(h,rturn-10.) ! CHECK h or r0 ?
      irmin=nint(rmin/h)+1

!      write(*,*) 'Enhanced Numerov from ir=',irmin,' r=',rmin,' fm',
!     & ' to  R=',Rm,'in steps of',h,' and matching at Rmatch=',rmatch


      if (debug) write(95,*)'Starting value ir=,',irmin,'r=',rvec(irmin)
      ym(:,:)=0d0; zm(:,:)=0d0; 
      y0(:,:)=0d0; z0(:,:)=0d0;
      do ich=1,nch
        l=ql(ich)
!        if (debug) write(klog,*)'ich,l,dlfac(l)',ich,l,dlfac2(l)
 
!50      aux= (kch(ich)*h)**(l+1)/exp(0.5*dlfac2(l))
!        if (aux.lt.1e-30) then
!          rmin=rmin+5*h
!          goto 50        
!        endif

!        irmin=nint(rmin/h) !+1
!        y0(ich,ich)=(kch(ich)*h)**(l+1)/exp(0.5*dlfac2(l))
!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        y0(ich,ich)=1e-10
!!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
        if (debug) then
!        write(95,'("Starting values:")')
        write(95,320) (y0(ich,is),is=1,min(10,nch)) 
320     format(5x,50("(",2g16.6,")",5x))
        endif
      enddo

!      if (debug) write(95,*)'Initial dot product(ABS):'
!      write(95,'(10(2x,"is=",i2,3x,1g10.5))') (is,
!     &     abs(dot_product(y0(:,is),y0(:,is-1)))
!     & , is=2,nch)

      write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
     &  5x,"outer=",1f8.2," fm",
     &  5x,"R-min=",1f6.2," fm")') 
     &  minval(rturnv(:)), maxval(rturnv(:)), rmin      

c 
c Initial values: w(rmin) and w(rmin+h) 
c 
      if (irmin.lt.2) write(*,*)'wrong irmin=',irmin
      call zfun(nch,ql,kch2,irmin-1,ym,zm)
      wm=ym-(h2/12.)*zm
      call zfun(nch,ql,kch2,irmin  ,y0,z0)
      w0=y0-(h2/12.)*z0
 
      do ir=irmin,nr-1  
      r=rvec(ir) ! (ir-1)*h    
      im=ir-1  ! r-h
      i0=ir    ! r
      ip=ir+1  ! r+h
      call csh(nch,ql,kch2,i0,w0,z0,method) ! Cosh[Sqrt(sqrt(h2 F(x))]

      if (debug.and.(mod(ir,10).eq.0)) then
         write(95,'(/,"y0(x) ir=",i4,3x,"r=",1f7.2)') ir,r
         do ich=1,min(nch,10)
      write(95,'(5x,50g14.5)')(y0(ich,is),is=1,min(nch,10))
         enddo
      endif !debug

      WP(:,:)=2*Z0(:,:)-WM(:,:)  

      if (debug.and.(mod(ir,10).eq.0)) then
        write(95,*)'wp at ir,r=',ir,r
        do ich=1,min(nch,10)
        write(95,320) (wp(ich,is),is=1,min(nch,10))
        enddo
        call flush(95)
      endif !debug
                         

c calculate and store Y=(1-T)W  
      if ((ir.ge.nr-8).or.(orto)) then      
        call w2y(wp,yp,ql,kch2,ir,nch)
        y(1:nch,1:nch,ip)=yp(1:nch,1:nch) 
      endif

c Re-orthogonalize solution vectors
      if (orto.and.(r.lt.rmorto).and.(ir.gt.irmin)
     & .and.(mod(ir-irmin,norto).eq.0)) then 
        write(*,*) 'orthogonalizing at ir,r=',ir,r

        if (debug) then
        write(94,*)'|yp x yp| (before G-S) for ir,r',ir,r
        do is=2,nch
        write(94,'(5x,50g14.5)') 
     &    ( abs(dot_product(yp(:,is),yp(:,ich)))
     &      / abs(dot_product(yp(:,is),yp(:,is)))
     &      / abs(dot_product(yp(:,ich),yp(:,ich))), ich=1,nch) 
        enddo
        endif !debug

  
!      call gsorto(yp,a,nch)  ! standard Gram-Schmidt
!      call mgsorto(yp,a,nch)  ! modified Gram-Schmidt
!      do i=irmin,ir+1
!      y(:,:,i)=matmul(y(:,:,i),a(:,:))
!      enddo ! i

       call cpu_time(ti)
       call qrfac(yp,y,nch,ip,irmin,nr)      ! QR factorization
!       call qrfac(wp,w,nch,ip,irmin,nr)      ! QR factorization



       call cpu_time(tf)
       torto=torto+tf-ti

      if (debug) then
      write(94,*)'Ortog. matrix at ir=',ir
      do ich=1,nch
        write(94,'(5x,50g14.5)') 
     &  (a(ich,is),is=1,nch)
           write(94,'(5x,50g14.5)') 
     &  (a(ich,is),is=1,nch)
      enddo
      call flush(94)
      endif 

      if (debug) then
      write(90,*)'y at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(90,'(5x,i3,50g14.5)') ich,
     &  (abs(y(ich,is,ip)), is=1,min(15,nch))
      enddo
      write(90,*)' '

      write(90,*)'yp at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(90,'(5x,i3,50g14.5)') ich,
     &  (abs(yp(ich,is)), is=1,min(15,nch))
      enddo
      write(*,*)' '

      write(90,*)'|y-yp|/y at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(90,'(5x,i3,50g14.5)') ich,
     &  (abs(y(ich,is,ip)-yp(ich,is))
     &   /abs(y(ich,is,ip)), is=1,min(15,nch))
      enddo
      write(*,*)' '
      endif

      yp=y(:,:,ip);  !call zfun(nch,ql,kch2,ip,yp,zp)
      y0=y(:,:,i0);  !call zfun(nch,ql,kch2,ip,y0,z0)
      call zfun(nch,ql,kch2,ip,yp,zp)
      wp=yp-(h2/12.)*zp
      call zfun(nch,ql,kch2,i0,y0,z0)
      w0=y0-(h2/12.)*z0


!      wp=w(:,:,ip);  !call zfun(nch,ql,kch2,ip,yp,zp)
!      w0=w(:,:,i0);  !call zfun(nch,ql,kch2,ip,y0,z0)


c dot products after orthogonalization
      if (debug) then
        write(94,*)'|yp x yp| (after G-S)'
        do is=1,nch
!        write(95,'(5x,i2,"x",i2,3x,50g14.5)') 
        write(94,'(5x,50g14.5)') 
     &    ( abs(dot_product(yp(:,is),yp(:,ich)))
     &      / abs(dot_product(yp(:,is),yp(:,is)))
     &      / abs(dot_product(yp(:,ich),yp(:,ich))), ich=1,nch) 

        enddo
       endif !debug

!       y(1:nch,1:nch,ir+1)=yp(1:nch,1:nch)  

!      if (debug) then
!        write(94,*)'Gram-Schmidt ortogonalized:'
!        do ich=1,nch
!        write(94,'(5x,50g14.5)') (yp(ich,is),is=1,nch)
!        enddo
!        write(94,*)' '
!        call flush(94)
!      endif !debug

      endif !orto
! ------------------------------------------------------------


c Alternative Gramm-Schmitt orthogonalization of solution vectors
!      orto=.false.
!      rmorto=maxval(rturnv(:))
!      if (orto.and.(r.lt.rmorto).and.(mod(ir,5).eq.0)) then 
!         call gsorto(yp,a,nch)
!      end if 

!!!!!!!!!!!!!!!!!!!!!!!!! TEST
      if (debug) then
       write(91,'(5x,1f6.3,3x,50g14.5)') r,(abs(y(1,ich,ip)),ich=1,nch)
      endif       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c
c wm <- w0 <- wp 
c
      wm(:,:)=w0(:,:) 
      w0(:,:)=wp(:,:) 

      enddo !ir
      call cpu_time(finish)


c 
c Match y with asymptotic solution -> S-matrix
c            
      call cpu_time(start)
!      call matching(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! gauss5 + derivative
!      call matching2(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! LAPACK 
      do ich=1,nch
      if (incvec(ich)) then
      inc=ich
      call matching3(ecm,z12,nch,ql,lmax,inc,nr,y,wf,phase,smat,info) ! gauss5 + 2 points
            smats(icc,inc,1:nch)=smat(1:nch) 
      endif
      enddo
       call cpu_time(finish)
       tmatch=tmatch+finish-start
       


      if (debug) write(klog,'(2f10.3,3x,10g14.6)') ecm,z12,ql(1),
     &  (phase(ich),ich=1,nch) 

      call flush(6)
      end subroutine schcc_ena_MGR

!-----------------------------------------------------------------------


c ------------------------------------------------------------------
c   Solve multichannel Schrodinger equation callint R-matrix routines
c   by P. Descouvemont CPC XXX XXX
c
c   NLAG = nb. of Laguerre bases per radial interval
c   NS   = nb. of radial intervals
c   -----------------------------------------------------------------
      subroutine schcc_rmat_MGR(nch,ecm,z12,incvec,ql,conv,dr,
     & r0,nr,wf,phase,smat,info,nlag,ns,einc,icc)
      use nmrv,only: vcoup,ech
      use xcdcc, only: rvcc,smats
      use constants , only: e2
      use memory
      implicit none
      logical :: iftr,twf,info
c     -------------------------------------------------------------
      integer :: ir,nma,nch,inc,nr,nmax,icc,ninc,i
      integer*4, parameter:: nr2=0
      integer*4 :: nlag,ns,ich,ichp,nopen
      integer*4,allocatable:: nvc(:)
c     -----------------------------------------------------------
      real*8 :: ecm,tkch,z12,conv
      real*8  :: r,rmax,dr,r0,einc,kinc
      real*8  :: zrma(nlag*ns),eta(nch),kch(nch),ql(nch)
      real*8  :: aux,jci,jcf
      real*8,parameter :: alpha=0.
!      real*8 fc(500),dfc(500),gc(500),dgc(500),xfac(nchan,nchan)
c     -----------------------------------------------------------
      complex*16 :: cfival,caux,phc,ph2
      complex*16 :: cpot(nlag*ns,nch,nch),cu(nch,nch),
     &              cpnl,faux(nr)
      complex*16 cfx(nch)
      complex*16, intent(out):: wf(nch,nr)
      complex*16 :: phase(nch),smat(nch) !,svel(nch)
      complex*16,allocatable :: cf(:,:,:)
      logical :: incvec(nch)
c     ----------------------------------------------------------
!      call cpu_time(start)
c *** Initialize some variables & R-matrix
      nmax=nlag*ns
      rmax=rvcc(nr)
      twf=.false.
      
      ninc=0
      do i=1,nch
      if (incvec(i)) ninc=ninc+1
      enddo
      allocate(nvc(ninc))
      i=1
      do ich=1,nch
      if (incvec(ich)) then
      nvc(i)=ich
      i=i+1
      endif
      enddo
      allocate(cf(nlag*ns,nch,ninc))
      
      
      call rmat_ini(nlag,ns,rmax,zrma)
c *** -----------------------------------


c *** Interpolation of coupling matrix in Lagrange mesh
      kinc=sqrt(conv*ecm)
      do ich=1,nch
      aux=ecm+einc-ech(ich)
      kch(ich)=sqrt(conv*abs(aux))
      eta(ich)=conv*z12*e2/kch(ich)/2.
      if (aux.lt.0) kch(ich)=-kch(ich)

      do ichp=ich,nch
      phc=(0d0,1d0)**(ql(ichp)-ql(ich))
      faux(1:nr)=conv*vcoup(ich,ichp,1:nr)
      do 4 ir=1,nmax
      r=zrma(ir)
      caux=cfival(r,rvcc,faux,nr,alpha)
      cpot(ir,ich,ichp)=caux
      if (ich.ne.ichp) cpot(ir,ichp,ich)=caux
!    &   *conjg(phc)
    4 continue
      enddo !ich
      enddo !ichp
      deallocate(vcoup) ! not needed any more ?????????????

! TEST
!      do ir=1,nmax
!        write(1,'(1f8.3,2x,50f12.8)')zrma(ir),cpot(ir,1,1)
!      enddo
!      do ir=1,nr
!        write(2,'(1f8.3,2x,50f12.8)')rvec(ir),ccmat(1,1,ir)
!      enddo 
      call rmatrix(nch,ql,kch,eta,rmax,nlag,ns,cpot,cu,nmax,nch,
     &             nopen,twf,cf,nmax,nch,ninc,nvc,0,cpnl)

       do ich=1,ninc
       i=nvc(ich)
       smat(1:nopen)=cu(1:nopen,i)*sqrt(kinc/kch(1:nopen))
       smats(icc,inc,1:nopen)=smat(1:nopen) 
       enddo

       do i=1,ninc
       inc=nvc(i)
       do ich=1,nch     
       write(*,500)ich,smat(ich),cu(ich,inc),abs(smat(ich))
500   format(10x,"Chan. #",i3,5x,
     &    "S=(",1f8.5,",",1f8.5,")", 5x, 
     &    "S=(",1f8.5,",",1f8.5,") ->  |S|=",1f9.6)     
       enddo
       enddo

       if (twf) then
       do  ir=1,nmax
       write(501,1002)zrma(ir),
     &      (cf(ir,ich,1)*(0.,-1.),ich=1,nch)
!     &      (cf(ir,ich,1)*(0.,-1.)*sqrt(kch(ich)/kch(inc)),ich=1,nch)
1002   format(f8.3,10(2es12.4,2x))
       enddo
       write(501,*)'&'
       endif
          
     
      return
      end subroutine


      SUBROUTINE ESMOOV(NC,QUAD,QUAD1,PHZ,EN,NUME,NFTA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C***********************************************************************
C     Downloaded from:  http://www.tampa.phys.ucl.ac.uk/rmat/old/eigenp.f
C
C     ESMOOV was formerly part of ETABLE and carries out smoothing
C     of eigenphases using second order finite differences
C
C***********************************************************************
C
      DIMENSION QUAD(NC),QUAD1(NC),PHZ(NC,NUME),EN(NUME)
      DATA PI/3.1415926535897932D+00/,ZERO/0.D0/,ONE/1.0D0/
C
C     INITIALIZE ARRAY QUAD1 WITH THE MULTIPLES OF PI WHICH ARE
C     TO BE ADDED TO THE EIGENPHASES AT THE FIRST ENERGY
C
      scale = one
      ieav = 0
      DO 3 I=1,NC
      QUAD1(I)=QUAD(I)
    3 CONTINUE
C
C     CHECK CONTINUITY OF EIGENPHASE PARAMETERS AS A FUNCTION OF
C     SCATTERING ENERGY BY THE USE OF SECOND-ORDER DIVIDED DIFFERENCES
C
      DO 190 NO=1,NC
C
C     INITIATE SMOOTHING ; CONTINUITY OF FIRST THREE ENERGY POINTS
C     IS ASSUMED
C
      Q1=QUAD1(NO)*PI
C
      II=3
      E1=EN(1)
      E2=EN(2)
      E3=EN(3)
      P1=PHZ(NO,1)+Q1
      P2=PHZ(NO,2)+Q1
      P3=PHZ(NO,3)+Q1
      PHZ(NO,1)=P1
      PHZ(NO,2)=P2
      PHZ(NO,3)=P3
C
C     FORM DIVIDED DIFFERENCES
C
      FD1=(P2-P1)/(E2-E1)
      FD2=(P3-P2)/(E3-E2)
      SD2=(FD2-FD1)/(E3-E1)
C
C     BEGIN ITERATION OVER ENERGY POINTS
C
    5 II=II+1
      IF(II .GT. NUME) GO TO 190
      Q1=QUAD1(NO)*PI
      P2=P3
      P3=PHZ(NO,II)+Q1
C
      E1=E2
      E2=E3
      E3=EN(II)
      FD1=FD2
      FD2=(P3-P2)/(E3-E2)
      SD1=SD2
      SD2=(FD2-FD1)/(E3-E1)
      EPS1=SCALE*DABS(SD1)
C
C     CHECK CONTINUITY
C
      IF(DABS(SD2-SD1) .GT. EPS1) GO TO 7
C
C     CONTINUITY CHECKS CORRECTLY
C
    6 PHZ(NO,II)=P3
      GO TO 5
C
C                     **************************
C
C     CORRECTION SEQUENCE :
C
C     (1) ADD QUADRANT CORRECTION OF PI
C
    7 P3=P3+PI
      FD2=(P3-P2)/(E3-E2)
      SD2=(FD2-FD1)/(E3-E1)
      IF(DABS(SD2-SD1) .GT. EPS1) GO TO 9
C
      QUAD1(NO)=QUAD1(NO)+ONE
      GO TO 6
C
C     (2) SUBTRACT QUADRANT CORRECTION OF PI
C
    9 P3=P3-PI-PI
      FD2=(P3-P2)/(E3-E2)
      SD2=(FD2-FD1)/(E3-E1)
      IF(DABS(SD2-SD1) .GT. EPS1) GO TO 11
C
      QUAD1(NO)=QUAD1(NO)-ONE
      GO TO 6
C
C     (3) CHECK WHETHER CHANNELS HAVE BEEN SWAPPED
C
   11 NO1=NO+1
      IF(NO1 .GT. NC) GO TO 20
C
      DO 18 NOX=NO1,NC
      IPT=1
      NOP=NOX
C
      P3=PHZ(NOP,II)+Q1
   12 FD2=(P3-P2)/(E3-E2)
      SD2=(FD2-FD1)/(E3-E1)
C
C     JUMP OUT IF SEARCH SUCCEEDS
C
      IF(DABS(SD2-SD1) .LE. EPS1) GO TO 25
C
      if(ipt.eq.1) then
C
C     TRY ADDING PI
C
        P3=P3+PI
        IPT=2
        GO TO 12
C
C     TRY SUBTRACTING PI
C
      else if(ipt.eq.2) then
        P3=P3-PI-PI
        IPT=3
        GO TO 12
      endif
C
   18 CONTINUE
C
C     SMOOTHING FAILURE
C
   20 WRITE(NFTA,22)NO,II,E3
   22 FORMAT(' CONTINUITY ERROR ENCOUNTERED IN CHANNEL ',I3,
     1       ' EIGENPHASE AT ENERGY E(',I3,') =',D16.8)
      P3=PHZ(NO,II)+Q1
      GO TO 5
C
C     SWAP EIGENPHASE COLUMNS NO AND NOP BEGINNING AT ENERGY POINT II
C
   25 DO 27 IEP=II,NUME
      SAVE=PHZ(NO,IEP)
      PHZ(NO,IEP)=PHZ(NOP,IEP)
      PHZ(NOP,IEP)=SAVE
   27 CONTINUE
C
C     REVERSE QUAD1 VALUES
C
      SAVE=QUAD1(NOP)
      QUAD1(NOP)=QUAD1(NO)
      QUAD1(NO)=SAVE
      if(ipt.eq.2) then
        QUAD1(NO)=QUAD1(NO)+ONE
      else if(ipt.eq.3) then
        QUAD1(NO)=QUAD1(NO)-ONE
      endif
C
      PHZ(NO,II)=P3
C
      WRITE(NFTA,32)NO,NOP,II
   32 FORMAT(' CHANNELS NO =',I3,' AND NOP =',I3,' SWAPPED BEGINNING',
     1       ' AT ENERGY POINT IEN =',I4)
      GO TO 5
C
  190 CONTINUE
      RETURN
      END

