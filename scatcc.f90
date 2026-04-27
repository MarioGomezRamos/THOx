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
!              =2 : Real CC wfs (NOT FUNCTIONAL!) 
!-------------------------------------------------------------------------------
      subroutine wfrange(nset,nchan,inc,ecmi,ecmf,ne,energy,wfc,
     & smate,delta,writewf)
      use nmrv,only:nch,ech,vcoup,hort,cutr,nopen
      use constants
      use channels, only:jpiset,jtot,ql,qspl,qj,qspj,qjc,exc,
     &                   cindex,spindex,partot,lscoup,qspstot
      use parameters, only:  maxchan
      use sistema
      use potentials, only: ccmat
      use scattering, only: ifcont,method,nbas,ns
      use wfs, only: nr,dr,rvec,wfsp,rint,rmin
      use globals, only: written,verb
      use forbidden ! added in v2.2g (TESTING)
      implicit none
c     ------------------------------------------------------------
      character*5 jpi
      logical :: info,energy,realcc,debug,writewf
c ... 
      integer ir,n,nchan,ne,ie,nrint,lmax !method
      integer :: nset,inc,ifail,wftype
      real*8  :: ecm,ecmi,ecmf,de,excore,econt
      real*8  :: dk,kf,ki,kcm,kcont,krm,tkch
      real*8  :: z12,rm,factor,r0,vscale
      real*8  :: deladd,deltap,deltai,delta(ne) !,cph(0:500)
      real*8 ,   allocatable  :: cph(:)
      complex*16 :: phase(nchan),smat(nchan),smate(ne,maxchan)
      complex*16 :: wfc(ne,nchan,nr),wf(nchan,nr),wfop(nchan,nchan,nr)
      complex*16 :: ci
      
c *** R-matrix solutions -----------------------------------------
      integer*4, parameter:: ninc=1,nr2=0,nmax=200
      logical :: iftr,twf,ifrmat
      integer*4 :: ich,ichp,nvc(inc)
      
      real*8 :: r,zrma(nmax),eta(nchan),kch(nchan),jci,jcf,rmatch
      real*8,parameter :: alpha=0.
      real*8 :: aux,auxp,fival
      complex*16 :: cfival,caux,phc,ph2,cauxp
      complex*16 :: cpot(nmax,nchan,nchan),cu(nchan,nchan),
     &              cf(nmax,nchan,ninc),faux(nr)
!      real*8 fc(500),dfc(500),gc(500),dgc(500),xfac(nchan,nchan)
      complex*16 cfx(nchan)
c *** Pauli blocking (TESTING) -----------------------------------------
      integer ispi,ip,ncp2,ii,ncorei,i,j
      integer:: lp,np,li,linc
      real*8:: ji,jp,rp
      real*8,allocatable:: wfpau(:)  
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
      
c     For S-matrix diagonalization
      integer sort
      parameter (sort = 1)
      complex*16:: umat(nchan,nchan),evec(nchan)
!      complex*16, allocatable::umat(:,:),evec(:),saux(:,:)
      real*8 :: WORK(NCHAN,NCHAN),
     &          WORK1(NCHAN),WORK2(3*NCHAN)

      real*8 :: big,small
      big=huge(big)
      small=epsilon(small)      
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
      wftype=1  !(1=scat. states; 2=real CC states 3=F(kr))
      debug=.false.
     
c *** -------------------------

c     ------------------------------------------------------------
      if ((npauli.gt.0).or.(method.eq.5)) then
        ifrmat=.true.  ! use R-matrix routine by P.Descouvemnt
      else 
        ncp2=0
        ifrmat=.false. ! use Numerov
      endif 
c     ------------------------------------------------------------

!      write(0,*)'method=',method; stop

      if (info) then
        write(*,100)ne,jpi(jtot,partot),inc,ecmi,ecmf
100     format(3x,"Generating", i4, " continuum wfs for J/pi=",
     &  a5,", inc. chan.",i3,
     &  " in [ECMmin=",1f6.3," ECMmax=",1f6.2," MeV]")
      endif 

!      bastype =jpiset(nset)%bastype
      jtot    =jpiset(nset)%jtot
      partot  =jpiset(nset)%partot
      excore  =jpiset(nset)%exc(inc)
      lmax    =maxval(jpiset(nset)%lsp(1:nchan))
      vscale  =jpiset(nset)%vscale
      
      allocate(cph(0:lmax))
!      write(0,*)'lsp()',jpiset(nset)%lsp(1:5)
!      write(0,*)'wfrange: lmax=',lmax
      
     
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
      qspstot(1:nch)=jpiset(nset)%stot(1:nch)
      cindex(1:nch)=jpiset(nset)%cindex(1:nch)


      call coefmat(nset,nch)
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
      cutr  =-100 
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
!      nbas   =60       ! Lagrange functions per radial interval
      rmatch =rvec(nr)  ! rvec(nr)
      twf    =.true.    ! generate also WFS
      ncp2   =0         ! non-local couplings
!      ns     =1         ! nb. of intervals

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
!      ph2=(-1)**NINT((JCI-JCF - ABS(JCI-JCF))/2.)
      faux(1:nr)=ccmat(ich,ichp,1:nr)
      do 4 ir=1,nbas*ns
      r=zrma(ir)
      caux=cfival(r,rvec,faux,nr,alpha)
      cpot(ir,ich,ichp)=caux
      if (ich.ne.ichp) cpot(ir,ichp,ich)=caux
    4 continue
      enddo !ich
      enddo !ichp
      do ir=1,nbas*ns
        write(1,'(1f8.3,2x,50f12.8)')zrma(ir),cpot(ir,1,1)
      enddo
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
 
         case (1)  ! CC continuum wfs
c *** Scattering states with numerov solution ................
         if (.not.ifrmat) then
!                  if (.not.ifrmat.and.(ecm>0.2)) then
!          call schcc_erwin(nch,ecm,zc*zv,inc,ql,factor,dr,r0,
!     &                    nr,wf,phase,smat,method,info)

          call schcc_erwin(nch,ecm,zc*zv,inc,ql,factor,dr,r0,
     &                    nrint,wf,phase,smat,method,info,eph)

!          call schcc(nch,ecm,zc*zv,inc,ql,factor,dr,r0,
!     &              nrint,wf,phase,smat,info)

           wfc(ie,1:nch,1:nrint)=wf(1:nch,1:nrint)
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
           li=ql(ich) 
          if (kch(ich)>0.) then ! open channel -> get phase-shifts and S-matrix
          phase(ich)=(0.,-0.5)*LOG(cu(ich,inc)
     &               *sqrt(kch(inc)/kch(ich)))*180/pi        
          smat(ich)=cu(ich,inc)*sqrt(kch(inc)/kch(ich))
          endif
c ... store wfs
!          write(0,*)'nrint,zrma(nbas)=',nrint,zrma(nbas)
          do ir=1,nrint
          r=rvec(ir)
          if (r.lt.zrma(nbas)) then
             wfc(ie,ich,ir)=cfival(r,zrma,cf(:,ich,inc),nbas,alpha)
     &                      *0.5*(0.,1.)
          else     
          
          krm=kch(ich)*rvec(ir)
          if (kch(ich).gt.0d0) then   ! open channel
!          write(0,*)'krm,ich,li,eta,smat,inc',krm,ich,li,eta(ich),
!     & smat(ich),inc
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
         endif ! r< zrma
                            
          enddo !nrint
          enddo !nch

c ...................................................................
!          if (ie.eq.1) then
!!          write(*,*)'Nopen=',nopen,' Elastic S-matrix=',cu(1,inc)
!!          do ich=1,nch
!!             write(*,*)'Chan=',ich,' S-mat=',
!!     &      cu(ich,inc)*sqrt(kch(inc)/kch(ich))
!!          enddo
!          do  ir=1,nbas*ns
!          write(501,1002)zrma(ir),
!     &      (cf(ir,ich,1)*(0.,-1.),ich=1,nch)
!1002      format(f8.3,10(2es12.4,2x))
!          enddo
!          write(501,*)'&'
!          endif

         ENDIF ! Choose integration method (Numerov / R-matrix)

! Added July '24
      case(3) ! Use free/Coulomb functions
      do ich=1,nch
      li=ql(ich) 
      if (li.gt.10) stop'increase lmax in wfrange!'
      tkch=ecm+ech(inc)-ech(ich) 
      kch(ich)=sqrt(factor*abs(tkch))
      eta(ich)=factor*z12*e2/kch(ich)/2.      
      do ir=1,nr
      krm=kch(ich)*rvec(ir)
      if (krm.lt.1e-6) krm=1e-6
      if (tkch.gt.0d0) then   ! open channel
         write(*,*)'krm=',krm, 'rvec=',rvec(ir)
        call coul90(krm,eta(ich),zero,li,f,g,fp,gp,0,IFAIL)
!        ch=cmplx(g(li),f(li))  ! H^(+)
!        wfc(ie,ich,ir)=(0.,0.5)*(conjg(ch)*kron(ich,inc)-smat(ich)*ch)
        wfc(ie,ich,ir)=f(li)
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


        case default
         write(*,*)'wfrange: wftype=',wftype,' not valid'
         stop
          
        end select 
c ....................................................................


c ... Extrapolate wfs from rint to rmax 
      if ((rint.gt.0).and.(rint.lt.rvec(nr)).and.wftype.ne.4) then
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


!!! Commented in June 2020
!!! Re-added in July 2020 
c ... Add Coulomb phase 
c    (could be done within ERWIN, but I do it here in order not to affect other parts of the code using ERWIN
      if (wftype.ne.4) then
      linc=ql(inc) 
      kch(inc)=sqrt(factor*abs(ecm))
      eta(inc)=factor*z12*e2/kch(inc)/2.    
      call coulph(eta(inc),cph,linc)
      phc=exp(ci*cph(linc))    
!      write(*,*)'phic=',phc
! AMM: should we restrit this to open channels???
      wfc(ie,:,:)= phc*wfc(ie,:,:)
      endif ! wftype

c ... If wfs are to be 'normalized' as <k|k'> = delta(k-k'), require additional factor sqrt(2/pi)
c      wf(:,:)=wf(:,:)*sqrt(2./pi)

c ... Store S-matrix and phase-shifts 
      smate(ie,1:nch)=smat(1:nch)

c ... Avoid pi jumps in phase-shifts
      deltai=phase(inc)  ! elastic phase shift
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
        write(45,205)jpi(jtot,partot),inc,min(ql(1:nchan),10)
205     format("#Phase-shifts for J/pi=",a5,"  Inc. chan.=",
     &         i3, " ql=",10i3)
        endif
        write(45,'(1f6.3,2x,10g14.5)') ecm, (phase(ich),ich=1,nchan)
!        write(46,'(1f6.3,2x,10f12.5)') ecm, 
!    & (sin(phase(ich)*pi/180),ich=1,nchan)
        


!	write(0,*)'writewf=',writewf; stop
        
!!!! Diagnostic TEST WF
        if (writewf) then
        written(500)=.true.
        if (ie.eq.1) then
        write(500,400) jpi(jtot,partot), nchan,inc,energy !, ql(1:nchan)
400     format(5x,"# Continuum WFS for Jp=",a5,"; Channels=",i2,
     & " ;  Inc. chan.=",i2, ' Uniform Ecm=',l1)
        write(500,'("#",1x,i5,2x,2f8.4)') nr,dr,r0
        if (energy) then
          write(500,'("#",1x,i5,2x,2f8.4)') ne,de,ecmi
        else
          write(500,'("#",1x,i5,2x,2f8.4)') ne,dk,ki
        endif
        endif

!        write(*,*)'Ecm=',Ecm,' Elastic S-matrix=',smat(1)
!        write(*,*)'kch=',kch(1:nch)
!        do ich=1,nch
!        write(*,*)'Chan=',ich,' S-mat=',
!     &    smat(ich)
!        enddo

        write(500,200)jpi(jtot,partot),ecm,inc,ql(1:min(nchan,10))
200     format("#WF for J/pi=",a5," Ecm=",1f7.3," MeV  Inc. chan.=",i3,
     &          " ql=",10i3)

        phc=exp(ci*cph(linc))  
        do ir=1,nr
        write(500,'(1f8.3,2x,50f12.8)') rvec(ir),
     & (exp(-ci*phase(inc)*pi/180.)*conjg(phc)*wfc(ie,n,ir),n=1,nch)
        enddo
        write(500,*)'&'
      endif
!        write(45,'(1f10.3,3x,10g14.6)') ecm,(phase(n),n=1,nch) 
!!!! TEST

      
      enddo ! energy loop (ie)

      write(45,*)'&'
      deallocate(ccmat,vcoup)
      RETURN
               
c     **********************************************************************
c Diagonalize S-matrix to get eigenphases.
c
      write(0,*)'ie,diagonalize S-matrix=, ifcont',ie,ifcont

      if (ifcont) then
!      write(0,*)'diagonalize S-matrix'
!      do ich=1,nopen
!      write(*,*)ich,'smat=',smat(ich)
!      enddo
!      write(0,*) 'nopen,sort=',nopen,sort
!      allocate (saux(nopen,nopen),Umat(nopen,nopen),evec(nopen))
!      saux(1:nopen,1:nopen)=smat(1:nopen,1:nopen)
      eph(:)=0
      call SEigensystem(nopen,smat,nch,evec, Umat,nopen, sort)

      do ich=1,nopen
         if(abs(evec(ich)).gt.small) eph(ich)=(0.,-0.5)*LOG(evec(ich))
!         print*,ich,' -> eigenphase=',eph(ich)
      enddo
!      print 101, "Eigen:",(n,evec(n),abs(evec(n)),eph(n),n = 1, nopen)
101   format(/A, 10(/"d(", I1, ") = ", F10.6, SP, F10.6, SS, " I", 
     & "|d|=",f10.6, " Phas=",f10.6))
      written(46)=.true.
      write(46,'(1f8.4,50f12.6)') ecm,(eph(n)*180./pi,n = 1, nchan)
      endif
!      deallocate(saux,umat,evec)
c Sort eigenphases
!      write(*,*)'Try making ephases continuous'

!      call ESMOOV(NCH,QUAD,QUAD1,PHZ,ENER,NE,IWRITE)
!      do ie=1,ne
!      write(94,'(1f8.4,20f12.6)') ener(ie),
!     &   (phz(n,ie),n = 1, nch)
!      enddo ! ie 


!      enddo ! energy loop (ie)

!      write(45,*)'&'
!      deallocate(ccmat,vcoup)
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
      real*8  ecm,z12,jtot,rm,factor,r0,vscale
      complex*16 :: phase(nchan),smat(nchan)
      complex*16 :: wf(nchan,nr)
c     ------------------------------------------------------------
      real*8 :: eph(nchan)
c     ------------------------------------------------------------
      rm=av*ac/(av+ac)
      factor=(2*amu/hc**2)*rm

      jtot    =jpiset(nset)%jtot
      partot  =jpiset(nset)%partot
      vscale  =jpiset(nset)%vscale
      
      
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

      call coefmat(nset,nch)
      vcoup(:,:,1:nr)=ccmat(:,:,1:nr)/factor

      r0=rvec(1)
      method=4       ! enhanced Numerov as used in Fresco
      hort  =0       ! no stabilization
      info  =.false. ! silent output
!     call schcc(nch,ecm,zv*zc,inc,ql,factor,dr,r0,nr-1,wf,phase,smat)
!      call schcc(nch,ecm,zv*zc,inc,ql,factor,dr,r0,
!     & nr,wf,phase,smat,info)
      call schcc_erwin(nch,ecm,zc*zv,inc,ql,factor,dr,r0,
     & nr,wf,phase,smat,method,info,eph)
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
      integer :: ir,nma,nch,inc,nr,nmax,ql(nch)
      integer*4, parameter:: ninc=1,nr2=0
      integer*4 :: nlag,ns,ich,ichp,nopen,nvc(ninc)
c     -----------------------------------------------------------
      real*8 :: ecm,tkch,z12,conv
      real*8  :: r,rmax,dr,r0
      real*8  :: zrma(nlag*ns),eta(nch),kch(nch)
      real*8  :: aux,jci,jcf
      real*8,parameter :: alpha=0.
!      real*8 fc(500),dfc(500),gc(500),dgc(500),xfac(nchan,nchan)
c     -----------------------------------------------------------
      complex*16 :: cfival,caux,phc,ph2
      complex*16 :: cpot(nlag*ns,nch,nch),cu(nch,nch),
     &              cf(nlag*ns,nch,ninc),faux(nr)
      complex*16, allocatable:: cpnl(:,:,:)    
      complex*16 cfx(nch)
      complex*16, intent(out):: wf(nch,nr)
      complex*16 :: phase(nch),smat(nch) !,svel(nch)
c     ----------------------------------------------------------
!      call cpu_time(start)
c *** Initialize some variables & R-matrix
      write(*,*) '=== schcc_rmat received: nlag=',nlag,' ns=',ns
      nmax=nlag*ns
      write(*,*) '=== nmax=nlag*ns=',nmax
      rmax=rvcc(nr)
      write(*,*) '=== rmax=',rmax,' nr=',nr
!      write(0,*)'rmat: rmax=',rmax
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

! TEST - print cpot for R-matrix
      write(*,*) '=== R-matrix cpot (first 10 points, channel 1,1) ==='
      do ir=1,min(10,nmax)
        write(*,'(a,i4,a,f8.3,a,2f12.6)')' ir=',ir,' r=',zrma(ir),
     &   ' cpot=',real(cpot(ir,1,1)),aimag(cpot(ir,1,1))
      enddo
      write(*,*) '=================================================='

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
       if (verb.ge.4)write(*,'(5x,"-> orthogonalizing at r=",1f7.2)')r

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
      integer:: nch,klog,npt,ql(nch),nopen
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
      nopen=0
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
        nopen=nopen + 1
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
      if (verb.ge.4)write(*,'(5x,"-> orthogonalizing at r=",1f7.2)')r
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
      l=ql(ich)
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
      use nmrv,only: nr,h,vcoup,conv,ech,debug,rmin,hort,rvec,cutr
      use globals, only: verb
      use factorials
      use constants , only: e2
      use memory
      implicit none
      integer:: method
      logical:: copen(nch),orto,raynal,info
      integer:: nch,klog,npt,ql(nch),nopen
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
c ... Eigenphases
      real*8 eph(nch)

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
      nopen =0
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

! CHanged by AMoro to prevent undefined mod(ir-irmin,norto) below
!      norto=nint(hort/h)
      norto=max(nint(hort/h),1)
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
!      write(0,*)'schccrc: calling matching4'
      call matching4(ecm,z12,nch,ql,lmax,inc,nr,y,wf,phase,smat,
     &               info,eph) ! gauss5 + 2 points
! Under development!
!      call match_real(ecm,z12,nch,nopen,ql,lmax,inc,nr,y,wf,info)

      if (debug) write(klog,'(2f10.3,3x,10g14.6)') ecm,z12,ql(1),
     &  (phase(ich),ich=1,nch) 

      call flush(6)
      end subroutine schcc_erwinrc




c ...
c ... Enhanced Numerov (version of IJ Thompson used in Fresco) 
c ... 
      subroutine schcc_erwin(nch,ecm,z12,inc,ql,factor,dr,r0,
     & npt,wf,phase,smat,method,info,eph)
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
c ... Eigenphases
      real*8 :: eph(nch)
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
      nopen=0
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
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius if derivative is used
     
! CHanged by AMoro to prevent undefined mod(ir-irmin,norto) below
!      norto=nint(hort/h)
      norto=max(nint(hort/h),1)   
      
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

      if (aux.gt.0) then
        copen(ich)=.true.
        nopen=nopen+1
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
      endif

!!!!!!!!!!!!!!!!!!!!!!1


      if (debug) then
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
      if (norto.gt.0) then
      if (orto.and.(r.lt.rmorto).and.(ir.gt.irmin).and.(hort.gt.0)
     & .and.(mod(ir-irmin,norto).eq.0).and.(ir.lt.nr-10)) then 
        if (verb.ge.4) 
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
      endif ! norto>0
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
!      call matching3(ecm,z12,nch,ql,lmax,inc,nr,y,wf,phase,smat,info) ! gauss5 + 2 points
      
!      write(0,*)'schcc_erwin: calling matching3_eph: ecm=',ecm
      call matching3_eph(ecm,z12,nch,ql,lmax,inc,nr,y,wf,phase,smat,
     &     info,eph) ! gauss5 + 2 points; eigenphase calculation

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

!      IF(SING) then
!      DO 605 Ir=1,nr
!605   IF(MAGN(Ir).NE.0.0) MAGN(Ir) = LOG10(MAGN(Ir))
!      WRITE(*,610) (MAGN(Ir),Ir=1,nr)
!610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
!      endif

c
c S-matrix and phase-shifts
c 
      flux=0d0
      phase(1:nch)=0.
      if (show.and.(verb.ge.1)) 
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
      if (verb.gt.0)
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
c in the definition of f(r). We do not include this factor here, 
c because is multiplied later on in wfrange() 
c
c USE 2  POINTS FOR MATCHING INSTEAD OF DERIVATIVE 
c ----------------------------------------------- --------------
      subroutine matching3(ecm,z12,nch,ql,lmax,inc,
     & n1,y,wf,phase,smat,show)
      use nmrv, only: nr,mu,ech,h,conv,rvec,nopen
      use constants, only : e2,pi
      use globals, only: verb,debug
      use trace , only: cdccwf
      use xcdcc, only : wfcdcc
      implicit none
      logical :: sing,show
      integer ich,is,inc,ir,klog,nd,n1,n2
      integer nch,ql(nch),l,linc,lmax,ifail,ie
      integer m1 ! coufg
c     --------------------------------------------------------------
      real*8, dimension(0:lmax):: f,g,gp,fp
      real*8, dimension(1:lmax+1):: f1,fp1
      real*8 :: kron,r1,r2
      real*8 :: ecm,tkch,eta,krm1,krm2,z12,kch(nch)
      real*8 :: phr,cph(0:lmax) !!!! changed to lmax check !!!!
      real*8 :: fpmax,small,acc8,flux,magn(nr)
      real*8 :: whit_rescale,asympt_whit
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
!!!!TEST
!        write(*,*)'krm1,eta,l=',krm1,eta,l
        call coul90(krm1,eta,0d0,l,f,g,fp,gp,0,IFAIL)
!        call coulfg(krm1,eta,0d0,1d0*l,f,g,fp,gp,1,0,ifail,m1)  !This was giving a problem

       if (ifail.ne.0) then 
       write(*,*) 'coul90: ifail=',ifail; ! stop
       endif
       ch1  = dcmplx(g(l),f(l))  * (0.,.5)   ! outgoing (i/2) H(+)
!!!! TEST
        call coul90(krm2,eta,0d0,l,f,g,fp,gp,0,IFAIL)
!        call coulfg(krm2,eta,0d0,1d0*l,f,g,fp,gp,1,0,ifail,m1)  !This was giving a problem

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

!MGR  Since the ANC is not relevant, let's try to avoid small values by rescaling both Whittaker values setting the one for r1 to the incoming channel
!      asympt_whit = 5d0*(eta+sqrt(eta**2 + l*(l+1)))/kch(ich)
!      if (r2.gt.asympt_whit) then
!      whit_rescale = exp(-kch(ich)*(r2-r1))*(r2/r1)**eta
!      mat(ich,nch+ich)    = -mat(inc,nch+inc)
!      mat(nch+ich,nch+ich)= -whit_rescale*mat(inc,nch+inc) 
!      else
!
!      mat(ich,nch+ich)    = -mat(inc,nch+inc)
!      mat(nch+ich,nch+ich)= -ch2/ch1*mat(inc,nch+inc)  
!      endif      
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
      write(*,*)'** WARNING ** Singular matrix for Ecm=',ecm
!      DO 605 Ir=1,nr
!605   IF(MAGN(Ir).NE.0.0) MAGN(Ir) = LOG10(MAGN(Ir))
!      WRITE(*,610) (MAGN(Ir),Ir=1,nr)
!610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
      endif

c
c S-matrix and phase-shifts
c 
      flux=0d0
      phase(1:nch)=0.
      if (show.and.(verb.ge.1)) 
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
       
      if (show) write(*,520)inc,flux
520   format(10x,"For incoming chan.",i2," => unitarity=",1f12.8) 

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

      if (cdccwf)then
          wfcdcc(inc,1:nch,1:nr)=wf(1:nch,1:nr)
          do is=1,nch
          write(85,*)'# Initial chan =>  Final chan'
          write(85,'(2i7)') inc,is
          write(85,'(6g16.6)') (wf(is,ir),ir=1,nr)
          enddo !is (final channel)
      endif

!!!! TEST
      if (0>1)then
          do ir=1,nr
!          write(91,*)'# Initial chan =>  Final chan'
!          write(91,'(2i7)') inc,is
          write(91,'(1f8.3,10g16.6)')rvec(ir), (wf(ich,ir),ich=1,5)
          enddo !is (final channel)
      endif


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
      subroutine matching3_eph(ecm,z12,nch,ql,lmax,iref,
     & n1,y,wf,phref,smati,show,eph)
      use nmrv, only: nr,mu,ech,h,conv,rvec
      use constants, only : e2,pi
      use globals, only: verb,debug,written
      use trace , only: cdccwf
      use scattering, only: ifcont
      implicit none
      logical :: sing,show,ifop(nch)
      integer ich,is,inc,iref,ir,klog,nd,n1,n2,iop,nopen
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
      complex*16:: phref(nch),phase(nch,nch),phc
      complex*16:: smati(nch),delta,scoul,anc,svel,smat(nch,nch)

c     ---------------------------------------------------------------
c     For S-matrix diagonalization
      integer sort,n
      parameter (sort = 1)
      complex*16:: umat(nch,nch),evec(nch)
      real*8 :: eph(nch)
!      real*8 :: WORK(NCH,NCH),
!    &          WORK1(NCH),WORK2(3*NCH)
c     ...............................................................

      acc8 = epsilon(acc8);
      fpmax = huge(acc8)**0.8d0 
      SMALL=1D0/FPMAX
      ci=(0d0,1d0)
      pi=acos(-1d0) 
      debug=.false.
      show=.false. !!!!!!!!! TEST 
      klog=99
      magn(:)=0      
      nd=6
      n2=n1-nd
      nopen=0; 
      ifop(:)=.true.
      phase(1:nch,1:nch)=0.

      if (nch.eq.0) then
       write(*,*)'In matching3_eph: nch=0! Abort'; stop     
      endif
      
! NEW: calculated S-matrix for all open channels in order to 
!      later compute the eigenphases 
      iop=0
      do inc=1,nch
      if (.not.ifcont.and.(iref.ne.inc)) cycle
      linc=ql(inc)
      tkch=ecm+ech(iref)-ech(inc) 
      if (tkch.lt.0) then
        smat(inc,:)=0
!        write(*,*)' No inc waves for chan=',inc
        ifop(inc)=.false.
        cycle
      endif      
      nopen=nopen+1
      
      mat(:,:)=0d0
      mat(1:nch,1:nch)=y(1:nch,1:nch,n1) 

      do ich=1,nch
      do is=1,nch
      mat(nch+ich,is)= y(ich,is,n1-nd)
      enddo
      enddo
      
      linc=ql(inc)
      do ich=1,nch
      l=ql(ich) 
      tkch=ecm+ech(iref)-ech(ich) 

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
      write(*,*)'** WARNING ** Singular matrix for Ecm=',ecm
!      DO 605 Ir=1,nr
!605   IF(MAGN(Ir).NE.0.0) MAGN(Ir) = LOG10(MAGN(Ir))
!      WRITE(*,610) (MAGN(Ir),Ir=1,nr)
!610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
      endif

c
c S-matrix and phase-shifts
c 
      flux=0d0

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
      if (abs(smati(ich)).gt.small) then
        phase(inc,ich)=(0.,-0.5) * LOG(smati(ich))
        test=(0.,-0.5) * LOG(smati(ich)) 
      endif

      if (show) then
       if (((ich.eq.inc).and.(verb.ge.1)).or.(verb.ge.3)) then
       write(*,500) ich,smati(ich),svel,abs(svel)!!!,test*180/pi
       endif
      endif
      flux=flux + abs(svel)**2
500   format(10x,"Chan. #",i3,5x,
     &    "S=(",1f8.5,",",1f8.5,")", 5x, 
     &    "S=(",1f8.5,",",1f8.5,") ->  |S|=",1f9.6) 

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
      endif
      enddo !ich
!      phase(1:nch)=phase(1:nch)*180/pi
      phref(1:nch)=phase(iref,1:nch)*180/pi

!      write(45,'(1f10.3,3x,10g14.6)') ecm,
!    &  (phase(ich),ich=1,nch) 
       
      if (show) write(*,'(10x,"=> Unitarity=",1f12.8)')flux


c
c scattering wavefunctions
c
      if (inc.eq.iref) then
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

!      wf(:,:)=dsqrt(2d0/pi)*phc*yfac*wf(:,:)
c To make the wfs real in the single-channel case:
!       wf(:,:)=dsqrt(2d0/pi)*yfac*wf(:,:)
! TEST
!      wf(:,:)=dsqrt(2d0/pi)*wf(:,:)
! TO MATCH WITH FRESCO FORT.17 , NO FACTORS ARE NEEDED HERE!

      if (cdccwf)then
          do is=1,nch
          write(85,*)'# Initial chan =>  Final chan'
          write(85,'(2i7)') inc,is
          write(85,'(6g16.6)') (wf(1,ir),ir=1,nr)
          enddo !is (final channel)
      endif
      endif ! inc=iref?


      enddo !inc


c
c Diagonalize S-matrix to get eigenphases.
c
      if (ifcont) then
!      do iop=1,nopen
!      write(*,*)iop,'smat=',smat(iop,1:nopen)
!      enddo

      eph(:)=0
!      write(0,*)'matching3_eph: calling seigensystem'
! change 3rd argument from open to nch
!      call SEigensystem(nopen,smat,nopen,evec, Umat,nopen, sort)
      call SEigensystem(nopen,smat,nch,evec, Umat,nopen, sort)

!       call dsyev('v','l',nopen,smat,nopen,work1,work2,3*nopen,
!     &              ifail)
!*          write (6,*)
!*          write (6,*)'optimal LWORK is',work2(1)

!      DO 2 N=1,NCH
!      EPH(n) = ATAN(WORK1(NCH+1-N))
! 2    continue


      do ich=1,nopen
         if(abs(evec(ich)).gt.small) eph(ich)=(0.,-0.5)*LOG(evec(ich))
!         print*,ich,' -> eigenphase=',eph(ich)
      enddo
!      print 101, "Eigen:",(n,evec(n),abs(evec(n)),eph(n),n = 1, nopen)
101   format(/A, 10(/"d(", I1, ") = ", F10.6, SP, F10.6, SS, " I", 
     & "|d|=",f10.6, " Phas=",f10.6))
     
      written(46)=.true.
      write(46,'(1f8.4,50f12.6)') ecm,(eph(n)*180./pi,n = 1, nch)
      endif



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

      if (nch.eq.0) then
       write(*,*)'Matching: nch=0! Abort'; stop     
      endif
      
      write(0,*)' Entering matching4:'
      
! NEW: calculated S-matrix for all open channels
      iop=0
      do inc=1,nch
      linc=ql(inc)
      tkch=ecm+ech(iref)-ech(inc) 
      if (tkch.lt.0) then
        smat(inc,:)=0
!        write(*,*)' No inc waves for chan=',inc
        cycle
      endif

      mat(:,:)=0d0
      mat(1:nch,1:nch)=y(1:nch,1:nch,n1) 

      do ich=1,nch
      do is=1,nch
      mat(nch+ich,is)= y(ich,is,n1-nd)
      enddo
      enddo

!      write(0,*)'inc=',inc,'ech=',ech(inc),'iref=',iref,' Ecm=',tkch

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
        nopen=nopen+1
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
      write(*,*)'** WARNING ** Singular matrix for Ecm=',ecm
!      DO 605 Ir=1,nr
!605   IF(MAGN(Ir).NE.0.0) MAGN(Ir) = LOG10(MAGN(Ir))
!      WRITE(*,610) (MAGN(Ir),Ir=1,nr)
!610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
      endif

c
c S-matrix and phase-shifts
c 
      flux=0d0
      phase(1:nch)=0.
      if (show.and.(verb.ge.1)) 
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
       
      if (show) write(*,'(10x,"=> unitarity=",1f12.8)')flux
      
!      write(96,'(50F12.5)') (smat(inc,ich),ich=1,nch)

      enddo !inc


c
c Diagonalize S-matrix to get eigenphases.
c
      nopen=iop
      write(0,*)'matching4: There are', nopen,' channels'

      if (nopen>0) then
!      do iop=1,nopen
!      write(*,*)iop,smat(iop,:)
!      enddo

      eph(:)=0
      call SEigensystem(nopen,smat,nch,evec, Umat,nch, sort)

!       call dsyev('v','l',nch,smat,nch,work1,work2,3*nch,
!     &              ifail)
!*          write (6,*)
!*          write (6,*)'optimal LWORK is',work2(1)
!      DO 2 N=1,NCH
!      EPH(n) = ATAN(WORK1(NCH+1-N))
! 2    continue


      do ich=1,nopen
         if(abs(evec(ich)).gt.small) eph(ich)=(0.,-0.5)*LOG(evec(ich))
!         print*,ich,' -> eigenphase=',eph(ich)
      enddo
!      print 101, "Eigen:",(n,evec(n),abs(evec(n)),eph(n),n = 1, nopen)
101   format(/A, 10(/"d(", I1, ") = ", F10.6, SP, F10.6, SS, " I", 
     & "|d|=",f10.6, " Phas=",f10.6))
      write(95,'(1f8.4,20f10.6)') ecm,(eph(n),n = 1, nopen)
      endif

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
!      debug=.true.
      show=.false.

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

! CHanged by AMoro to prevent undefined mod(ir-irmin,norto) below
!      norto=nint(hort/h)
      norto=max(nint(hort/h),1)
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


      if (info.and.(verb.ge.2)) then
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
       if (verb.ge.4) 
     &   write(*,'(5x,"-> orthogonalizing at r=",1f7.2)')r

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
c Using two adjacent points
      do ich=1,nch
      if (incvec(ich)) then
      inc=ich
!     call matching2(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! LAPACK 
      call matching3(ecm,z12,nch,ql,lmax,inc,nr,y,wf,phase,smat,info) ! gauss5
      endif
      smats(icc,inc,1:nch)=smat(1:nch) 
      enddo

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
      use trace , only: cdccwf
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
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius if derivative is used
     
! CHanged by AMoro to prevent undefined mod(ir-irmin,norto) below
!      norto=nint(hort/h)
      norto=max(nint(hort/h),1)   

      if (abs(hort).gt.0.) orto=.true.

      if (info.and.(verb.ge.3)) then
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
      do ich=1,nch
      aux=ecm+einc-ech(ich)
      kch2(ich)=conv*aux

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
 
      call factorialgen(lmax+1)
      call cpu_time(start)
c 
c Starting values and initial integration point
c 
!      rmin=h ;       irmin=1
!      rmin=max(h,rturn-10.) ! CHECK h or r0 ?
      rmin  =max(h,minval(rturnv(:)-abs(cutr)))
      irmin=nint(rmin/h)+1

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
      
      if (info.and.(verb.ge.2)) then
      write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
     &  5x,"outer=",1f8.2," fm",
     &  5x,"R-min=",1f6.2," fm ",i4)') 
     &  minval(rturnv(:)), maxval(rturnv(:)), rmin,irmin      
      endif

      if (irmin.lt.2) write(*,*)'** ERWIN: wrong irmin=',irmin

! TEST - print S-matrix for Numerov (MGR version)
      write(*,*) '=== Numerov will calculate S-matrix ==='

c Start radial integration
      do ir=irmin,nr
      r=rvec(ir) ! (ir-1)*h
      im=ir-1  ! r-h
      i0=ir    ! r
      ip=ir+1  ! r+h
      if (r.lt.1e-10) r=1e-10

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
      endif

!!!!!!!!!!!!!!!!!!!!!!1


      if (debug) then
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
      if (verb.ge.4) 
     &   write(*,'(5x,"->orthogonalizing at ir,r=",i4,1f7.2)')ir,r
  
       call cpu_time(ti)
       call qrerwinz(y0,z0,zm,w0,nch,ir)      ! QR factorization

       y(:,:,ir)=y0(:,:) ! new (orthogonolized) solution at ir


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

       call cpu_time(tf)
       torto=torto+tf-ti

      endif !orto
!........................................................................


!!!!!!!!!!!!!!!!!!!!!!!!! TEST
      if (0>1) then
       write(91,'(5x,1f7.3,3x,50g14.5)') r,(y(1,ich,ir),ich=1,nch)
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
! DEBUG - print S-matrix from Numerov
      write(*,'(a,i2,a,2f10.5,a,f10.6)')
     &  ' Numerov ch',1,': S=(',real(smat(1)),aimag(smat(1)),
     &  ') |S|=',abs(smat(1))
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

! CHanged by AMoro to prevent undefined mod(ir-irmin,norto) below
!      norto=nint(hort/h)
      norto=max(nint(hort/h),1)

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

      if (verb.ge.2) 
     & write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
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
        if (verb.ge.4) 
     &   write(*,*) 'orthogonalizing at ir,r=',ir,r

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

      call flush(6)
      end subroutine schcc_ena_MGR

!-----------------------------------------------------------------------


c ------------------------------------------------------------------
c   Solve multichannel Schrodinger equation calling R-matrix routines
c   by P. Descouvemont CPC XXX XXX
c
c   NLAG = nb. of Laguerre bases per radial interval
c   NS   = nb. of radial intervals
c   -----------------------------------------------------------------
      subroutine schcc_rmat_MGR(nch,ecm,z12,incvec,ql,conv,dr,
     & r0,nr,wf,phase,smat,info,nlag,ns,einc,icc)
      use globals, only: written
      use nmrv,only: vcoup,ech
      use xcdcc, only: rvcc,smats
      use constants , only: e2,pi
      use scattering, only: ifcont
      use memory
      implicit none
      logical :: iftr,twf,info
c     -------------------------------------------------------------
      integer :: ir,nma,nch,inc,nr,nmax,icc,ninc,i,ql(nch)
      integer, parameter  :: ncp2=0  ! nb of non-local couplings
      integer*4, parameter:: nr2=0
      integer*4 :: nlag,ns,ich,ichp,nopen
      integer*4,allocatable:: nvc(:)
c     -----------------------------------------------------------
      real*8 :: ecm,tkch,z12,conv
      real*8  :: r,rmax,dr,r0,einc,kinc
      real*8  :: zrma(nlag*ns),eta(nch),kch(nch)
      real*8  :: aux,jci,jcf
      real*8,parameter :: alpha=0.
!      real*8 fc(500),dfc(500),gc(500),dgc(500),xfac(nchan,nchan)
c     -----------------------------------------------------------
      complex*16 :: cfival,caux,phc,ph2
      complex*16 :: cpot(nlag*ns,nch,nch),cu(nch,nch),faux(nr)
      complex*16, allocatable:: cpnl(:,:,:)    
!      complex*16:: cpnl(0,nch,nch)   
      complex*16 cfx(nch)
      complex*16, intent(out):: wf(nch,nr)
      complex*16 :: phase(nch),smat(nch) !,svel(nch)
      complex*16,allocatable :: cf(:,:,:)
      logical :: incvec(nch)
c     ----------------
c     For S-matrix diagonalization
      integer sort,n
      parameter (sort = 1)
      complex*16:: umat(nch,nch),evec(nch)
      real*8 :: eph(nch)
!      real*8 :: WORK(NCH,NCH),
!     &          WORK1(NCH),WORK2(3*NCH)
      real*8 :: big,small
      big=huge(big)
      small=epsilon(small)
c     ----------------------------------------------------------
!      call cpu_time(start)
c *** Initialize some variables & R-matrix
      write(*,*) '=== schcc_rmat_MGR: nlag=',nlag,' ns=',ns
      nmax=nlag*ns
      write(*,*) '=== nmax=nlag*ns=',nmax
      rmax=rvcc(nr)
      write(0,*)'Rmat: rmax=',rmax
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
!      deallocate(vcoup) ! not needed any more ?????????????

! TEST - print parameters passed to rmatrix
      write(*,*) '=== R-matrix parameters ==='
      write(*,*) 'nch=',nch,' nlag=',nlag,' ns=',ns,' rmax=',rmax
      write(*,*) 'conv=',conv,' z12=',z12
      do ich=1,min(3,nch)
        write(*,'(a,i2,a,f10.5,a,f10.5,a,i2)')
     &   ' ch',ich,': kch=',kch(ich),' eta=',eta(ich),' l=',ql(ich)
      enddo
      write(*,*) '=== zrma (first/last 5 pts) ==='
      do ir=1,min(5,nmax)
        write(*,'(a,i4,a,f12.6)') ' zrma(',ir,')=',zrma(ir)
      enddo
      do ir=max(1,nmax-4),nmax
        write(*,'(a,i4,a,f12.6)') ' zrma(',ir,')=',zrma(ir)
      enddo
      write(*,*) '=== R-matrix cpot (first 10 pts, ch 1,1) ==='
      do ir=1,min(10,nmax)
        write(*,'(a,i4,a,f8.3,a,2f12.6)')' ir=',ir,' r=',zrma(ir),
     &   ' cpot=',real(cpot(ir,1,1)),aimag(cpot(ir,1,1))
      enddo
      write(*,*) '=== vcoup raw (first 10 pts, ch 1,1) ==='
      do ir=1,min(10,nr)
        write(*,'(a,i4,a,f8.3,a,2f12.6)')' ir=',ir,' r=',rvcc(ir),
     &   ' vcoup*conv=',conv*vcoup(1,1,ir),0.d0
      enddo

      call rmatrix(nch,ql,kch,eta,rmax,nlag,ns,cpot,cu,nmax,nch,
     &             nopen,twf,cf,nmax,nch,ninc,nvc,ncp2,cpnl)

  
       do ich=1,ninc
       inc=nvc(ich)
       smat(1:nopen)=cu(1:nopen,inc)*sqrt(kinc/kch(1:nopen))
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
       
       

c
c Diagonalize S-matrix to get eigenphases.
c
!      write(0,*)'diagonalize S-matrix ifcont=',ifcont

      if (ifcont) then
!      write(0,*)'diagonalize S-matrix'
!      do ich=1,nopen
!      write(*,*)ich,'smat=',smat(ich,1:nopen)
!      enddo

      eph(:)=0
      call SEigensystem(nopen,smat,nch,evec, Umat,nopen, sort)

      do ich=1,nopen
         if(abs(evec(ich)).gt.small) eph(ich)=(0.,-0.5)*LOG(evec(ich))
         print*,ich,' -> eigenphase=',eph(ich)
      enddo
      print 101, "Eigen:",(n,evec(n),abs(evec(n)),eph(n),n = 1, nopen)
101   format(/A, 10(/"d(", I1, ") = ", F10.6, SP, F10.6, SS, " I", 
     & "|d|=",f10.6, " Phas=",f10.6))
      written(46)=.true.
      write(46,'(1f8.4,50f12.6)') ecm,(eph(n)*180./pi,n = 1, nch)
      endif
          
     
      return
      end subroutine


c   -----------------------------------------------------------------
c   Solve multichannel Schrodinger equation calling HPRMAT routines
c   (High-Performance R-matrix)
c
c   NLAG = nb. of Lagrange bases per radial interval
c   NS   = nb. of radial intervals
c   -----------------------------------------------------------------
      subroutine schcc_rmat_hp_MGR(nch,ecm,z12,incvec,ql,conv,dr,
     & r0,nr,wf,phase,smat,info,nlag,ns,einc,icc)
      use globals, only: written
      use nmrv,only: vcoup,ech
      use xcdcc, only: rvcc,smats
      use constants , only: e2,pi
      use scattering, only: ifcont
      use memory
      use rmat_hp_mod  ! HPRMAT module
      implicit none
      logical :: iftr,twf,info
c     -------------------------------------------------------------
      integer :: ir,nma,nch,inc,nr,nmax,icc,ninc,i,ql(nch)
      integer, parameter  :: ncp2=0  ! nb of non-local couplings
      integer*4, parameter:: nr2=0
      integer*4 :: nlag,ns,ich,ichp,nopen
      integer*4,allocatable:: nvc(:)
c     -----------------------------------------------------------
      real*8 :: ecm,tkch,z12,conv
      real*8  :: r,rmax,dr,r0,einc,kinc
      real*8  :: zrma(nlag*ns),eta(nch),kch(nch)
      real*8  :: aux,jci,jcf
      real*8,parameter :: alpha=0.
c     -----------------------------------------------------------
      complex*16 :: cfival,caux,phc,ph2
      complex*16 :: cpot(nlag*ns,nch,nch),cu(nch,nch),faux(nr)
      complex*16, allocatable:: cpnl(:,:,:)
      complex*16 cfx(nch)
      complex*16, intent(out):: wf(nch,nr)
      complex*16 :: phase(nch),smat(nch)
      complex*16,allocatable :: cf(:,:,:)
      logical :: incvec(nch)
c     ----------------
c     For S-matrix diagonalization
      integer sort,n
      parameter (sort = 1)
      complex*16:: umat(nch,nch),evec(nch)
      real*8 :: eph(nch)
      real*8 :: big,small
      big=huge(big)
      small=epsilon(small)
c     ----------------------------------------------------------
c *** Initialize some variables & R-matrix
      nmax=nlag*ns
      rmax=rvcc(nr)
      write(0,*)'HPRMAT Rmat: rmax=',rmax,' solver_type=',solver_type
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

      call rmat_ini_hp(nlag,ns,rmax,zrma)
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
      do ir=1,nmax
      r=zrma(ir)
      caux=cfival(r,rvcc,faux,nr,alpha)
      cpot(ir,ich,ichp)=caux
      if (ich.ne.ichp) cpot(ir,ichp,ich)=caux
      enddo
      enddo !ich
      enddo !ichp

! DEBUG - print parameters for HPRMAT
      write(*,*) '=== HPRMAT parameters ==='
      write(*,*) 'nr=',nr,' rvcc(nr)=',rvcc(nr)
      write(*,*) 'nch=',nch,' nlag=',nlag,' ns=',ns,' rmax=',rmax
      write(*,*) 'conv=',conv,' z12=',z12
      do ich=1,min(3,nch)
        write(*,'(a,i2,a,f10.5,a,f10.5,a,i2)')
     &   ' ch',ich,': kch=',kch(ich),' eta=',eta(ich),' l=',ql(ich)
      enddo
      write(*,*) 'zrma(1)=',zrma(1),' zrma(nmax)=',zrma(nmax)
      write(*,*) '=== cpot(ir,1,1) imaginary part ==='
      do ir=1,nmax,nmax/10
        write(*,'(a,i4,a,f8.2,a,2f12.4)') ' ir=',ir,' r=',zrma(ir),
     &   ' cpot=',real(cpot(ir,1,1)),aimag(cpot(ir,1,1))
      enddo

      call flush(6)
      call rmatrix_hp(nch,ql,kch,eta,rmax,nlag,ns,cpot,cu,nmax,nch,
     &             nopen,twf,cf,nmax,nch,ninc,nvc,ncp2,cpnl)


       do ich=1,ninc
       inc=nvc(ich)
       smat(1:nopen)=cu(1:nopen,inc)*sqrt(kinc/kch(1:nopen))
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
1002   format(f8.3,10(2es12.4,2x))
       enddo
       write(501,*)'&'
       endif


c
c Diagonalize S-matrix to get eigenphases.
c
      if (ifcont) then
      eph(:)=0
      call SEigensystem(nopen,smat,nch,evec, Umat,nopen, sort)

      do ich=1,nopen
         if(abs(evec(ich)).gt.small) eph(ich)=(0.,-0.5)*LOG(evec(ich))
         print*,ich,' -> eigenphase=',eph(ich)
      enddo
      print 101, "Eigen:",(n,evec(n),abs(evec(n)),eph(n),n = 1, nopen)
101   format(/A, 10(/"d(", I1, ") = ", F10.6, SP, F10.6, SS, " I",
     & "|d|=",f10.6, " Phas=",f10.6))
      written(46)=.true.
      write(46,'(1f8.4,50f12.6)') ecm,(eph(n)*180./pi,n = 1, nch)
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


c ...
c ... Enhanced Numerov (version of IJ Thompson used in Fresco) 
c ... 
      subroutine schcc_erwin_cuts(nch,ecm,z12,incvec,ql,factor,dr,r0,
     & npt,wf,phase,smat,method,info,einc,icc)
      use xcdcc, only: smats
      use nmrv,only: nr,h,vcoup,conv,ech,debug,rmin,hort,rvec,cutr
      use globals, only: verb
      use factorials
      use constants , only: e2
      use memory
      use trace , only: cdccwf
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
      real*8 :: kch(nch),kch2(nch),rturnv(nch),etav(nch),barrier(nch)
      real*8 :: start,finish,aux,big,small,tmin,tmax
      real*8 :: tol_mag,tol_mag_inc,tol_closed
      real*8, dimension(nch):: wkb
      integer:: rt_points(nch)
      integer:: irminaux
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
      complex*16:: yaux(nch,nch),ytolsma(nch,nch,5)
      logical incvec(nch)
      real*8:: renorm,tolsma
      integer jr,count_tolsma
      real*8:: maxdiff,maxz0,minz0
      complex*16:: qry0(nch,nch),qrz0(nch,nch),qrzm(nch,nch)
      complex*16:: qrw0(nch,nch),sqr(nch)
      integer::  permute(nch)
c     ---------------------------------------------------------
      R12 = 1D0/12D0
      ENA2 = 2D0/5D0 * R12**2
      ENA3 = - 4D0/35D0 * R12**3
!      ONEC = (1D0,0D0)
      TMAX = 20.
      TMIN = -125.
      big=huge(big)**0.8
      small=epsilon(small)**0.8
     
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
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius if derivative is used
     
! CHanged by AMoro to prevent undefined mod(ir-irmin,norto) below
!      norto=nint(hort/h)
      norto=max(nint(hort/h),1)   

      if (abs(hort).gt.0.) orto=.true.

      if (info.and.(verb.ge.3)) then
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
      do ich=1,nch
      aux=ecm+einc-ech(ich)
      kch2(ich)=conv*aux

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

! In order to not fully kill closed channels set turning point to rturn
      do ich=1,nch
      if (.not. copen(ich)) then
            rturnv(ich)=rturn
      endif
      enddo


      call factorialgen(lmax+1)
      call cpu_time(start)
c 
c Starting values and initial integration point
c 
!      rmin=h ;       irmin=1
!      rmin=max(h,rturn-10.) ! CHECK h or r0 ?
      rmin  =max(h,minval(rturnv(:)-abs(cutr)),rturn-abs(cutr)) 
      irmin=nint(rmin/h)+1

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
      
      if (info.and.(verb.ge.2)) then
      write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
     &  5x,"outer=",1f8.2," fm",
     &  5x,"R-min=",1f6.2," fm")') 
     &  minval(rturnv(:)), maxval(rturnv(:)), rmin      
      endif

!Displace rturnv by cutr MGR
      do ich=1,nch
            rturnv(ich)=rturnv(ich)-abs(cutr)
      enddo


      if (irmin.lt.2) write(*,*)'** ERWIN: wrong irmin=',irmin
 
! Use WKB approximation to find cut point. Set a tolerance, which is the number of orders of magnitude we allow our wf to
! decrease within the classically forbidden region and add all decreases until we reach that tolerance
      tol_mag=20d0
      tol_mag_inc=20d0
      tol_closed= 20d0    


      rt_points(:)=0
      wkb(:)=0d0
      do ir=nr-1,2,-1
      r=rvec(ir) ! (ir-1)*h 
      do ich=1,nch
      if (kch2(ich).lt.0d0)cycle
        barrier(ich)= (-l*(l+1)/r**2         ! centrifugal
     &          + kch2(ich)             ! energy
     &        - conv*dble(vcoup(ich,ich,ir))) ! potential
        if (barrier(ich).lt.0d0.and.rt_points(ich).eq.0) then !Compute WKB within forbidden region
             wkb(ich)=wkb(ich)+
     &        sqrt(abs(barrier(ich)))*(rvec(3)-rvec(2))*log10(exp(1d0))
            if(.not. incvec(ich).and.wkb(ich).ge.tol_mag .and.
     &        rt_points(ich).eq.0) then
               rt_points(ich)=ir
            endif
            if(incvec(ich).and.wkb(ich).ge.tol_mag_inc .and.
     &        rt_points(ich).eq.0) then
               rt_points(ich)=ir
            endif
        endif
      enddo
      enddo


!For closed channels find the first point where the kinetic energy is positive
      do ich=1,nch
      if (kch2(ich).lt.0d0) then
      do ir=2,nr-1
      r=rvec(ir)
      if (ir.eq.nr-1) then
        rt_points(ich)=nr !Channel is always closed, do not bother with it 
        exit
      endif
      barrier(ich)= (-l*(l+1)/r**2         ! centrifugal
     &          + kch2(ich)             ! energy
     &        - conv*dble(vcoup(ich,ich,ir))) ! potential
        if (barrier(ich).ge.(-tol_closed*conv).and.rt_points(ich).eq.0) 
     &   then
            rt_points(ich)=max(2,ir-nint(abs(cutr)/h)) !set cut point CUTR before
            exit  
      endif
      enddo
      endif
      enddo


      irminaux=1000000
      do ich=1,nch
      if (rt_points(ich).eq.0) rt_points(ich)=2
      rturnv(ich)=max(rvec(rt_points(ich)),rturnv(ich))!Set turning point to WKB cut point if larger
      !if (incvec(ich)) irminaux=min(irminaux,rt_points(ich))! Fix irminaux to smallest turning point of incoming channels
      irminaux=min(irminaux,rt_points(ich))
      enddo
      
      write(444,*) 'irmin',irmin

!Set all turning points to the minimum value, it seems different channels should have the same :(
      do ich=1,nch
            rturnv(ich)=minval(rturnv(:))
      enddo

      irmin=max(irmin,irminaux)
      !irmin=irminaux  !ignore cutr

      if (2>1) then
           write(444,*) 'Turning points Jtot'
      do ich=1,nch
            irminaux=min(rt_points(ich),nr)
      write(444,*) ich,rt_points(ich),wkb(ich),rturnv(ich),irmin,
     & copen(ich),kch2(ich)
      enddo      
      endif

      write(*,*) 'BIG',big,' SMALL',small
      count_tolsma=0
c Start radial integration
      do ir=irmin,nr
      r=rvec(ir) ! (ir-1)*h    
      im=ir-1  ! r-h
      i0=ir    ! r
      ip=ir+1  ! r+h
      if (r.lt.1e-6) r=1e-6

!Check if renormalization is needed
      renorm=1d0
      do ich=1,nch
        do is=1,nch
          renorm=max(renorm,abs(y0(ich,is)),abs(z0(ich,is)),
     &      abs(zm(ich,is)))
        enddo
      enddo
!      if (mod(ir,10) .eq. 0) write(*,*)'Maximum at',ir,':',renorm
      if (renorm.gt.big) then
      renorm=1d0/renorm
      write(*,*)'Renormalizing y0 at ir,r=',ir,r,' renorm=',renorm
      do is=1,nch
      do ich=1,nch
      z0(ich,is)=z0(ich,is)*renorm
      zm(ich,is)=zm(ich,is)*renorm
      y0(ich,is)=y0(ich,is)*renorm
        do jr=irmin,ir-1
            y(ich,is,jr)=y(ich,is,jr)*renorm
        enddo
      enddo
      enddo
      endif

! Check maximum difference between points
      maxz0=1e-30
      minz0=1e30
      do is=1,nch
      do ich=1,nch
      maxz0=max(maxz0,abs(z0(ich,is)))
      if (abs(z0(ich,is)).gt.1d0/big) then
       minz0=min(minz0,abs(z0(ich,is)))
      endif
      enddo
      enddo
      maxdiff=maxz0/minz0
      if (maxdiff.gt.1e20.and.maxz0.gt.1e10)then
      write(*,*)'Ratio between components at ir,r=',ir,r,
     & ' max=',maxz0,' min=',minz0,' ratio=',maxdiff

      write(*,*)'Renormalizing y0 at ir,r=',ir,r,' renorm=',1d0/maxz0
      renorm=1d0/maxz0
      do is=1,nch
      do ich=1,nch
      z0(ich,is)=z0(ich,is)*renorm
      zm(ich,is)=zm(ich,is)*renorm
      y0(ich,is)=y0(ich,is)*renorm
        do jr=irmin,ir-1
            y(ich,is,jr)=y(ich,is,jr)*renorm
        enddo
      enddo
      enddo
      endif

      




c kinetic energy part
c Fresco: SMAT(ich) = -LL1(K)*RI2 + (-ECM(ich) + C) * H
      do ich=1,nch    
      l=ql(ich)
      s(ich)= h2*(-l*(l+1)/r**2         ! centrifugal
     &          + kch2(ich)             ! energy
     &        - conv*vcoup(ich,ich,ir)) ! potential
      enddo

!Correct for rturnv MGR
      do ich=1,nch
      if (r.lt.rturnv(ich)) s(ich)=0d0
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

!Correct for rturnv MGR
      do ich=1,nch
      if (r.lt.rturnv(ich))w0(ich,ich)=ONEC
      enddo

!      y0=matmul(w0,z0)

c non-diagonal part of V matrix 
      coupl(:,:)=-h2*conv*vcoup(:,:,ir)  
      do ich=1,nch
         coupl(ich,ich)=0d0
      enddo 

!  Correct for rturnv MGR
      do ich=1,nch
      do is=1,nch
         if (r.lt.rturnv(ich).or. r.lt.rturnv(is)) then
         coupl(ich,is)=0d0
         if (ich.eq.is) then
            y0(ich,ich)=1e-10
         else
            y0(ich,is)=0d0
         endif
         endif
      enddo
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
      endif

!!!!!!!!!!!!!!!!!!!!!!1


      if (debug) then
      write(*,*)'z0 at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (z0(ich,is), is=1,min(15,nch))
      enddo

      endif ! debug


c store solution 
      y(:,:,ir)=y0(:,:)       

!     Compute nu factor from Tolsma CPC 40 233 (1986)
      if (ir.eq. irmin) write(447,*) '#R        Tolsma'
      count_tolsma=count_tolsma+1!Avoid using for the derivative points before and after orthogonalization
      tolsma=1d0
      if (count_tolsma.ge.5) then
      ytolsma(:,:,:) =y(:,:,ir-4:ir)
      call nutolsma(nch,s,ytolsma,tolsma,H,ir)
      write(447,*) r,tolsma
      endif


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
      if (ir.eq.nr) write(447,*)
      if (ir.eq.nr) write(447,*)

c Re-orthogonalize solution vectors ....................................
!      if (orto.and.(r.lt.rmorto).and.(ir.gt.irmin).and.(hort.gt.0)
!     & .and.(mod(ir-irmin,norto).eq.0).and.(ir.lt.nr-10)) then 
      if (orto .and. tolsma .lt.5e-4 .and. (ir.gt.irmin).and.
     & (ir.lt.nr-10)) then
      if (verb.ge.4) 
     &   write(*,'(5x,"->orthogonalizing at ir,r=",i4,1f7.2)')ir,r
  
         write(447,'(5x,"#->orthogonalizing at ir,r=",i4,1f7.2)')ir,r
       call cpu_time(ti)

!Permute columns and rows according to energy (lowest energy at last channel)

      sqr=s
      permute=0
      do ich=1,nch
      is=minloc(dble(sqr),1)
        permute(nch+1-ich)=is
        sqr(is)=onec*1e30
      enddo

!Now perform the permutation

      do ich=1,nch
      do is=1,nch
      qry0(ich,is)=y0(permute(ich),permute(is))
      qrz0(ich,is)=z0(permute(ich),permute(is))
      qrzm(ich,is)=zm(permute(ich),permute(is))
      qrw0(ich,is)=w0(permute(ich),permute(is))
      enddo
      enddo

       call qrerwinz(qry0,qrz0,qrzm,qrw0,nch,ir)      ! QR factorization
       count_tolsma=0!Reset count for derivative

!Undo permutation

      do ich=1,nch
      do is=1,nch
      y0(permute(ich),permute(is))=qry0(ich,is)
      z0(permute(ich),permute(is))=qrz0(ich,is)
      zm(permute(ich),permute(is))=qrzm(ich,is)
      w0(permute(ich),permute(is))=qrw0(ich,is)
      enddo
      enddo


       y(:,:,ir)=y0(:,:) ! new (orthogonolized) solution at ir


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

       call cpu_time(tf)
       torto=torto+tf-ti

      endif !orto
!........................................................................


!!!!!!!!!!!!!!!!!!!!!!!!! TEST
      if (0>1) then
       write(91,'(5x,1f7.3,3x,50g14.5)') r,(y(1,ich,ir),ich=1,nch)
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
      end subroutine schcc_erwin_cuts

      subroutine schcc_erwin_cuts_Q(nch,ecm,z12,incvec,ql,factor,dr,r0,
     & npt,wf,phase,smat,method,info,einc,icc)
      use xcdcc, only: smats
      use nmrv,only: nr,h,vcoup,conv,ech,debug,rmin,hort,rvec,cutr
      use globals, only: verb
      use factorials
      use constants , only: e2
      use memory
      use trace , only: cdccwf
      implicit none
      
      integer:: method
      logical:: copen(nch),orto,raynal,info
      integer:: nch,klog,npt,ql(nch)
      integer:: ir,irmin,is,isp,ich,inc,l,lmax
      integer:: i,j,k,im,ip,i0,norto,icc
c     .........................................................
      real*8     :: ymem,ti,tf
      real*8 , parameter:: pi=acos(-1d0)
      real*8 :: ecm,tkch,eta,z12
      real*16 ::r12,ena2,ena3,h2
      real*8 :: dr,r,r0,rm,factor,einc,kinc
      real*8 :: rturn,rmorto,rmatch
      real*8 :: kch(nch),kch2(nch),rturnv(nch),etav(nch),barrier(nch)
      real*8 :: start,finish,aux,tmin,tmax
      real*8 :: tol_mag,tol_mag_inc,tol_closed
      real*8, dimension(nch):: wkb
      integer:: rt_points(nch)
      integer:: irminaux
c     .........................................................
      complex*16, intent(out):: wf(nch,nr)
      
! MODIFIED: y remains complex*16 as requested
      complex*16 :: y(nch,nch,npt)

! MODIFIED: y0, z0, zm promoted to complex*32
      complex*32, dimension(nch,nch):: y0
      complex*32, dimension(nch,nch):: z0,zm
      complex*16 :: y0qr(nch,nch),z0qr(nch,nch),zmqr(nch,nch)
      complex*16 :: w0qr(nch,nch)

! MODIFIED: v moved from complex*16 to complex*32
      complex*16, dimension(nch,nch):: coupl
      complex*32, dimension(nch,nch):: v
      
      complex*16 :: phase(nch),smat(nch)
      complex*16 :: a(nch,nch)

! MODIFIED: s, ziv, zpv promoted to complex*32
      complex*32 :: s(nch),ziv,zpv
! MODIFIED: c1, c2, c promoted to complex*32
      complex*32 :: c1,c2,c
      
! MODIFIED: onec and zero are now complex*32 parameters
      complex*32,parameter:: onec=(1.0Q0, 0.0Q0),zero=(0.0Q0, 0.0Q0)
c ... for QR factorization
! MODIFIED: w0, wm promoted to complex*32
      complex*16 :: HDIAG(nch)
      complex*32 :: w0(nch,nch),wm(nch,nch)
c ... TEST (delete after debugging)
! MODIFIED: yaux promoted to complex*32
      complex*32:: yaux(nch,nch)
      
      logical incvec(nch)
      real*16:: renorm
      integer jr
      complex*16:: ytolsma(nch,nch,5),stolsma(nch)
      real*8:: tolsma
      integer count_tolsma
      real*16:: maxdiff,maxz0,minz0
      complex*32:: qry0(nch,nch),qrz0(nch,nch),qrzm(nch,nch)
      complex*32:: qrw0(nch,nch),sqr(nch)
      integer::  permute(nch)
      real*16 big,small
c     ---------------------------------------------------------
      R12 = 1Q0/12Q0
      ENA2 = 2Q0/5Q0 * R12**2
      ENA3 = - 4Q0/35Q0 * R12**3
!
! MODIFIED: This line is redundant with the parameter, but kept for logic
!      ONEC = (1.0Q0, 0.0Q0) 
      TMAX = 20.
      TMIN = -125.
      big=huge(big)**(0.8Q0)
      small=epsilon(small)**0.8Q0
     
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
      
! MODIFIED: y remains complex*16 and is initialized with 0d0 (double)
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
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius if derivative is used
     
! CHanged by AMoro to prevent undefined mod(ir-irmin,norto) below
!      norto=nint(hort/h)
      norto=max(nint(hort/h),1)   

      if (abs(hort).gt.0.) orto=.true.
      if (info.and.(verb.ge.3)) then
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
      do ich=1,nch
      aux=ecm+einc-ech(ich)
      kch2(ich)=conv*aux

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

! In order to not fully kill closed channels set turning point to rturn
      do ich=1,nch
      if (.not. copen(ich)) then
            rturnv(ich)=rturn
      endif
      enddo

 
      call factorialgen(lmax+1)
      call cpu_time(start)
c 
c Starting values and initial integration point
c 
!      rmin=h ;       irmin=1
!      rmin=max(h,rturn-10.) ! CHECK h or r0 ?
      rmin  =max(h,minval(rturnv(:)-abs(cutr)),rturn-abs(cutr)) 
      irmin=nint(rmin/h)+1

      if (debug) write(95,*)'Starting value ir=,',irmin,'r=',rvec(irmin)
      
! MODIFIED: Initializing complex*32 variables with quad 'zero'
      z0(:,:)=zero;
      zm(:,:)=zero; 
      w0(:,:)=zero; wm(:,:)=zero
      y0(:,:)=zero; v(:,:)=zero 
      do ich=1,nch
        l=ql(ich)
!!!!!      CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ZI(L,L) = H**(QNF(9,L)+1) / EXP(0.5 * FACT(QNF(9,L)+1))
! MODIFIED: Use quad precision for constant
        z0(ich,ich)= 1.0Q-20  !TEEEEEEEEEEEEEEEEEEEEEESSSSST!
!        z0(ich,ich)= H**(l+1.) / EXP(0.5 * FACT(l+1.))  
!         write(*,*)'ich,z0=:',  ich,z0(ich,ich)   
!!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
      
      if (info.and.(verb.ge.2)) then
      write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
     &  5x,"outer=",1f8.2," fm",
     &  5x,"R-min=",1f6.2," fm")') 
     &  minval(rturnv(:)), maxval(rturnv(:)), rmin      
      endif

!Displace rturnv by cutr MGR
      do ich=1,nch
            rturnv(ich)=rturnv(ich)-abs(cutr)
      enddo


      if (irmin.lt.2) write(*,*)'** ERWIN: wrong irmin=',irmin
 
! Use WKB approximation to find cut point for open channels.
! Set a tolerance, which is the number of orders of magnitude we allow our wf to
!
! decrease within the classically forbidden region and add all decreases until we reach that tolerance
      tol_mag=20d0
      tol_mag_inc=20d0
      tol_closed= 20d0  

      rt_points(:)=0
      wkb(:)=0d0

      do ir=nr-1,2,-1
      r=rvec(ir) 
      do ich=1,nch
      if (kch2(ich).lt.0d0) cycle !only for open channels
        barrier(ich)= (-l*(l+1)/r**2         ! centrifugal
     &          + kch2(ich)             ! energy
     &        - conv*dble(vcoup(ich,ich,ir))) ! potential
        if (barrier(ich).lt.0d0.and.rt_points(ich).eq.0) then !Compute WKB within forbidden region
             wkb(ich)=wkb(ich)+
     &        sqrt(abs(barrier(ich)))*h*log10(exp(1d0))
            if(.not. incvec(ich).and.wkb(ich).ge.tol_mag .and.
     &        rt_points(ich).eq.0) then
               rt_points(ich)=ir
            endif
            if(incvec(ich).and.wkb(ich).ge.tol_mag_inc .and.
     &        rt_points(ich).eq.0) then
               rt_points(ich)=ir
            endif
        endif
      enddo
      enddo

!For closed channels find the first point where the kinetic energy is positive
      do ich=1,nch
      if (kch2(ich).lt.0d0) then
      do ir=2,nr-1
      r=rvec(ir)
      if (ir.eq.nr-1) then
        rt_points(ich)=nr !Channel is always closed, do not bother with it 
        exit
      endif
      barrier(ich)= (-l*(l+1)/r**2         ! centrifugal
     &          + kch2(ich)             ! energy
     &        - conv*dble(vcoup(ich,ich,ir))) ! potential
        if (barrier(ich).ge.(-tol_closed*conv).and.rt_points(ich).eq.0) 
     &   then
            rt_points(ich)=max(2,ir-nint(abs(cutr)/h)) !set cut point CUTR before
            exit  
      endif
      enddo
      endif
      enddo

      irminaux=1000000
      do ich=1,nch
      if (rt_points(ich).eq.0) rt_points(ich)=2
      rturnv(ich)=max(rvec(rt_points(ich)),rturnv(ich))!Set turning point to WKB cut point if larger
      !if (incvec(ich)) irminaux=min(irminaux,rt_points(ich))! Fix irminaux to smallest turning point of incoming channels
      irminaux=min(irminaux,rt_points(ich))
      enddo

      write(444,*) 'irmin:',irmin

!Set all turning points to the minimum value, it seems different channels should have the same :(
      do ich=1,nch
            rturnv(ich)=minval(rturnv(:))
      enddo

      irmin=max(irmin,irminaux)

      if (2>1) then
      write(444,*) 'Turning points Jtot'
      do ich=1,nch
            irminaux=min(rt_points(ich),nr)
      write(444,*) ich,rt_points(ich),wkb(ich),rturnv(ich),irmin,
     & copen(ich),kch2(ich)
      enddo      
      endif

      write(*,*) 'BIG',big,' SMALL',small
      count_tolsma=0
c Start radial integration
      do ir=irmin,nr
      r=rvec(ir) ! (ir-1)*h    
      im=ir-1  ! r-h
      i0=ir    ! r
      ip=ir+1  ! r+h
      if (r.lt.1.0Q-6) r=1.0Q-6

!Check if renormalization is needed
      renorm=1Q0
! MODIFIED: abs(y0) etc. are now quad precision operations
      do ich=1,nch
        do is=1,nch
          renorm=max(renorm,abs(y0(ich,is)),abs(z0(ich,is)),
     &      abs(zm(ich,is)))
        enddo
      enddo
!      if (mod(ir,10) .eq. 0) write(*,*)'Maximum at',ir,':',renorm
      if (renorm.gt.big) then
      renorm=1Q0/renorm
      write(*,*)'Renormalizing y0 at ir,r=',ir,r,' renorm=',renorm
! MODIFIED: renormalizing complex*32 arrays z0, zm, y0
      do is=1,nch
      do ich=1,nch
      z0(ich,is)=z0(ich,is)*renorm
      zm(ich,is)=zm(ich,is)*renorm
      y0(ich,is)=y0(ich,is)*renorm
! MODIFIED: y remains complex*16, renorm is real*8 (double)
        do jr=irmin,ir-1
            y(ich,is,jr)=y(ich,is,jr)*renorm
        enddo
      enddo
      enddo
      endif

! Check maximum difference between points
      maxz0=1.0Q-50
      minz0=1.0Q50
      do is=1,nch
      do ich=1,nch
      maxz0=max(maxz0,abs(z0(ich,is)))
      if (abs(z0(ich,is)).gt.1.0Q0/big) then
       minz0=min(minz0,abs(z0(ich,is)))
      endif
      enddo
      enddo
      maxdiff=maxz0/minz0
      if (maxdiff.gt.1e35.and.maxz0.gt.1e18)then
      write(*,*)'Ratio between components at ir,r=',ir,r,
     & ' max=',maxz0,' min=',minz0,' ratio=',maxdiff

      write(*,*)'Renormalizing y0 at ir,r=',ir,r,' renorm=',1d0/maxz0
      renorm=1.0Q0/maxz0
      do is=1,nch
      do ich=1,nch
      z0(ich,is)=z0(ich,is)*renorm
      zm(ich,is)=zm(ich,is)*renorm
      y0(ich,is)=y0(ich,is)*renorm
        do jr=irmin,ir-1
            y(ich,is,jr)=y(ich,is,jr)*renorm
        enddo
      enddo
      enddo
      endif

      




c kinetic energy part
c Fresco: SMAT(ich) = -LL1(K)*RI2 + (-ECM(ich) + C) * H
! MODIFIED: s is now complex*32
      do ich=1,nch    
      l=ql(ich)
      s(ich)= h2*(-l*(l+1)/r**2         ! centrifugal
     &          + kch2(ich)             ! energy
     &        - conv*vcoup(ich,ich,ir)) ! potential
      enddo

!Correct for rturnv MGR
      do ich=1,nch
      if (r.lt.rturnv(ich)) s(ich)=zero
      enddo

! MODIFIED: y0, z0, s are complex*32, onec is complex*32
      do is=1,nch
      do ich=1,nch      
      y0(ich,is)= z0(ich,is)*
     X      (ONEC - S(ich) * (R12 - S(ich)*(ENA2 + S(ich)*ENA3)))
      enddo
      enddo

! MODIFIED: w0 is complex*32
      w0(:,:)=zero
     
       do is=1,nch
      w0(is,is)=(ONEC - S(is) * (R12 - S(is)*(ENA2 + S(is)*ENA3)))
      enddo 

!Correct for rturnv MGR
      do ich=1,nch
      if (r.lt.rturnv(ich))w0(ich,ich)=ONEC
      enddo

!      y0=matmul(w0,z0)

c non-diagonal part of V matrix 
      coupl(:,:)=-h2*conv*vcoup(:,:,ir)  
      do ich=1,nch
         coupl(ich,ich)=0d0
      enddo 

!  Correct for rturnv MGR
      do ich=1,nch
      do is=1,nch
         if (r.lt.rturnv(ich).or. r.lt.rturnv(is)) then
         coupl(ich,is)=zero
         if (ich.eq.is) then
! MODIFIED: Use quad precision constant
            y0(ich,ich)=1.0Q-20
         else
            y0(ich,is)=zero
         endif
         
         endif
      enddo
      enddo

! MODIFIED: w0 is complex*32
      w0=w0-r12*coupl
      
      if (2>1) then
      DO 24 ICH=1,NCH
	 DO 24 IS=1,NCH
! MODIFIED: C is complex*32
        C = COUPL(is,ich) * R12
! MODIFIED: y0, z0 are complex*32
	 if(C/=ZERO) y0(ich,1:nch) = y0(ich,1:nch) - C * z0(is,1:nch)
24    CONTINUE
        DO 34 ich=1,nch ! k
        V(ich,:)  = ZERO 
        DO 34 is=1,nch   ! j
	 C = COUPL(is,ich)
! MODIFIED: V, y0 are complex*32
 	 if(C/=ZERO) V(ich,1:nch) = V(ich,1:nch) + C * y0(is,1:nch)
34    CONTINUE
      else  !..........Idem using builtin matrix multiplication 
! MODIFIED: V is complex*32
        V(1:nch,1:nch)  = ZERO 
! MODIFIED: y0, z0 are complex*32
        y0=y0  - r12*matmul(coupl,z0)
! MODIFIED: v, y0 are complex*32
        v =   matmul(coupl,y0)
      endif

!! TEST 
      if (debug) then
! MODIFIED: yaux, w0, z0 are complex*32
      yaux=matmul(w0,z0)
      write(*,*)'y0 at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (y0(ich,is), is=1,min(15,nch))
      enddo
      endif

!!!!!!!!!!!!!!!!!!!!!!1


      if (debug) then
      write(*,*)'z0 at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (z0(ich,is), is=1,min(15,nch))
      enddo

      endif ! debug


c store solution 
! MODIFIED: This performs the requested precision demotion
! y (complex*16) <-- y0 (complex*32)
      y(:,:,ir)=y0(:,:)       

!     Compute nu factor from Tolsma CPC 40 233 (1986)
      if (ir.eq. irmin) write(447,*) '#R        Tolsma'
      count_tolsma=count_tolsma+1!Avoid using for the derivative points before and after orthogonalization
      tolsma=1d0
      if (count_tolsma.ge.5) then
      stolsma=s
      ytolsma(:,:,:) =y(:,:,ir-4:ir)
      call nutolsma(nch,stolsma,ytolsma,tolsma,H,ir)
      write(447,*) r,tolsma
      endif


c zm <- z0, wm <-w0
      do is=1,nch
      do ich=1,nch
! MODIFIED: ziv, z0 are complex*32
      ziv=z0(ich,is)
! MODIFIED: All variables here are complex*32. Use 2.0Q0 for quad constant.
      zpv=2.0Q0*ziv - zm(ich,is) - v(ich,is) - s(ich)*y0(ich,is)
! MODIFIED: zm, ziv, z0, zpv are complex*32
      zm(ich,is)=ziv
      z0(ich,is)=zpv
! MODIFIED: wm, w0 are complex*32
      wm(ich,is)=w0(ich,is)
      enddo  !ich
      enddo  !is


c Re-orthogonalize solution vectors ....................................
!      if (orto.and.(r.lt.rmorto).and.(ir.gt.irmin).and.(hort.gt.0)
!     & .and.(mod(ir-irmin,norto).eq.0).and.(ir.lt.nr-10)) then 
      if (orto .and. tolsma .lt.5e-4 .and. (ir.gt.irmin).and.
     & (ir.lt.nr-10)) then
      if (verb.ge.4) 
     &   write(*,'(5x,"->orthogonalizing at ir,r=",i4,1f7.2)')ir,r
  
         write(447,'(5x,"#->orthogonalizing at ir,r=",i4,1f7.2)')ir,r
       call cpu_time(ti)
! MODIFIED: Call to qrerwinz now passes complex*32 arrays
! NOTE: The subroutine qrerwinz MUST be modified to accept these.
!       y0qr=y0;z0qr=z0;zmqr=zm;w0qr=w0
       call qrerwinz_Q(y0,z0,zm,w0,nch,ir)      !
! QR factorization

!Permute columns and rows according to energy (lowest energy at last channel)

      sqr=s
      permute=0
      do ich=1,nch
      is=minloc(dble(sqr),1)
        permute(nch+1-ich)=is
        sqr(is)=onec*1e30
      enddo

!Now perform the permutation

      do ich=1,nch
      do is=1,nch
      qry0(ich,is)=y0(permute(ich),permute(is))
      qrz0(ich,is)=z0(permute(ich),permute(is))
      qrzm(ich,is)=zm(permute(ich),permute(is))
      qrw0(ich,is)=w0(permute(ich),permute(is))
      enddo
      enddo

       call qrerwinz_Q(qry0,qrz0,qrzm,qrw0,nch,ir)      ! QR factorization
       count_tolsma=0!Reset count for derivative

!Undo permutation

      do ich=1,nch
      do is=1,nch
      y0(permute(ich),permute(is))=qry0(ich,is)
      z0(permute(ich),permute(is))=qrz0(ich,is)
      zm(permute(ich),permute(is))=qrzm(ich,is)
      w0(permute(ich),permute(is))=qrw0(ich,is)
      enddo
      enddo


       y(:,:,ir)=y0(:,:) ! new (orthogonolized) solution at ir


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call cpu_time(tf)
       torto=torto+tf-ti

      endif !orto
!........................................................................


!!!!!!!!!!!!!!!!!!!!!!!!! TEST
      if (0>1) then
! MODIFIED: y remains complex*16, so this write statement is unchanged
       write(91,'(5x,1f7.3,3x,50g14.5)') r,(y(1,ich,ir),ich=1,nch)
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
      write(447,*)
      write(447,*)
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
      end subroutine schcc_erwin_cuts_Q

c *** ----------------------------------------
C ... Re-orthogonalize Z's by QR factorization quad precision
c *** ----------------------------------------
      subroutine qrerwinz_Q(yi,zi,zm,wi,n,ir)
      use nmrv, only: h,rvec
      implicit none
      logical    debug
c ... for zgeqrf
      integer info,lda,lwork,n,m,i,is,ich,ir
      complex*32 a(n,n),tau(n),work(2*n),rt(n,n)
      complex*32,dimension(n,n):: yi,zi,zm,wi
c ... for ztrsm
      complex*32 alpha
      integer ldb
      character diag,side,transa,uplo
c ...
      integer   :: nrhs, ipiv(n)
      character*1  trans
      EXTERNAL  ZGETRF_Q, ZGETRS_Q

c ---- TEST      
      complex*32 zaux(n,n),waux(n,n)
c .................................................
      lwork=2*n
      debug=.false.
     

c ... QR factorization ZI=Q.RT 
c     (Q=orthogonal matrix; R= triangular UPPER matrix)
      rt=zi
      call ZGEQRF_Q( N, N, RT, n, TAU, WORK, LWORK, INFO )
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
      
      call ZTRSM_Q(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,ZI,LDB)
      call ZTRSM_Q(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,ZM,LDB)


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
      call zgetrs_Q(trans,n,nrhs,waux,n,ipiv,zaux,ldb,info)

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
      call ZGETRF_Q(n, n, wi, n, ipiv, info)
      if (info /= 0) then
       write(*,*)'zgetrf returned error',info; stop
      endif

c ... Use IPIV & zgetrs to solve system  W*X= Y   (X=Z)
      trans='N';  ldb=n;  nrhs=n
      zi(:,:)= yi(:,:) 
      call zgetrs_Q(trans,n,nrhs,wi,n,ipiv,zi,ldb,info)

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


      SUBROUTINE ZGEQRF_Q( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*32         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGEQRF computes a QR factorization of a complex M-by-N matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the unitary matrix Q as a
*          product of min(m,n) elementary reflectors (see Further
*          Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is
*          the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEQR2_Q, ZLARFB_Q, ZLARFT_Q
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEQRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'ZGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZGEQRF', ' ', M, N, -1,
     $                 -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code initially
*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
*
*           Compute the QR factorization of the current block
*           A(i:m,i:i+ib-1)
*
            CALL ZGEQR2_Q( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL ZLARFT_Q( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H' to A(i:m,i+ib:n) from the left
*
               CALL ZLARFB_Q( 'Left', 'Conjugate transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
*
*     Use unblocked code to factor the last or only block.
*
      IF( I.LE.K )
     $   CALL ZGEQR2_Q( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,
     $                IINFO )
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of ZGEQRF
*
      END
      

      SUBROUTINE ZGEQR2_Q( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*32         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGEQR2 computes a QR factorization of a complex m by n matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(m,n) by n upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the unitary matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*32         ONE
      PARAMETER          ( ONE = ( 1.0Q+0, 0.0Q+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
      COMPLEX*32         ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF_Q, ZLARFG_Q
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEQR2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = 1, K
*
*        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
*
         CALL ZLARFG_Q( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,
     $                TAU( I ) )
         IF( I.LT.N ) THEN
*
*           Apply H(i)' to A(i:m,i+1:n) from the left
*
            ALPHA = A( I, I )
            A( I, I ) = ONE
            CALL ZLARF_Q( 'Left', M-I+1, N-I, A( I, I ), 1,
     $                  conjg( TAU( I ) ), A( I, I+1 ), LDA, WORK )
            A( I, I ) = ALPHA
         END IF
   10 CONTINUE
      RETURN
*
*     End of ZGEQR2
*
      END


      SUBROUTINE ZLARFB_Q( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*32         C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARFB applies a complex block reflector H or its transpose H' to a
*  complex M-by-N matrix C, from either the left or the right.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply H or H' from the Left
*          = 'R': apply H or H' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply H (No transpose)
*          = 'C': apply H' (Conjugate transpose)
*
*  DIRECT  (input) CHARACTER*1
*          Indicates how H is formed from a product of elementary
*          reflectors
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Indicates how the vectors which define the elementary
*          reflectors are stored:
*          = 'C': Columnwise
*          = 'R': Rowwise
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  K       (input) INTEGER
*          The order of the matrix T (= the number of elementary
*          reflectors whose product defines the block reflector).
*
*  V       (input) COMPLEX*16 array, dimension
*                                (LDV,K) if STOREV = 'C'
*                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
*                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*          if STOREV = 'R', LDV >= K.
*
*  T       (input) COMPLEX*16 array, dimension (LDT,K)
*          The triangular K-by-K matrix T in the representation of the
*          block reflector.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,K)
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          If SIDE = 'L', LDWORK >= max(1,N);
*          if SIDE = 'R', LDWORK >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*32         ONE
      PARAMETER          ( ONE = ( 1.0Q+0, 0.0Q+0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZCOPY_Q, ZGEMM_Q, ZLACGV_Q, ZTRMM_Q
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          conjg
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'C'
      ELSE
         TRANST = 'N'
      END IF
*
      IF( LSAME( STOREV, 'C' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1 )    (first K rows)
*                     ( V2 )
*           where  V1  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C1'
*
               DO 10 J = 1, K
                  CALL ZCOPY_Q( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV_Q( N, WORK( 1, J ), 1 )
   10          CONTINUE
*
*              W := W * V1
*
             CALL ZTRMM_Q( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2
*
               CALL ZGEMM_Q( 'Conjugate transpose', 'No transpose', N,
     $                        K, M-K, ONE, C( K+1, 1 ), LDC,
     $                        V( K+1, 1 ), LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL ZTRMM_Q( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2 * W'
*
                  CALL ZGEMM_Q( 'No transpose', 'Conjugate transpose',
     $                        M-K, N, K, -ONE, V( K+1, 1 ), LDV, WORK,
     $                        LDWORK, ONE, C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL ZTRMM_Q( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - conjg( WORK( I, J ) )
   20             CONTINUE
   30          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C1
*
               DO 40 J = 1, K
                  CALL ZCOPY_Q( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
*
*              W := W * V1
*
             CALL ZTRMM_Q( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2
*
               CALL ZGEMM_Q( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL ZTRMM_Q( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2'
*
              CALL ZGEMM_Q( 'No transpose', 'Conjugate transpose', M,
     $                        N-K, K, -ONE, WORK, LDWORK, V( K+1, 1 ),
     $                        LDV, ONE, C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL ZTRMM_Q( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
*
         ELSE
*
*           Let  V =  ( V1 )
*                     ( V2 )    (last K rows)
*           where  V2  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C2'
*
               DO 70 J = 1, K
                  CALL ZCOPY_Q( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV_Q( N, WORK( 1, J ), 1 )
   70          CONTINUE
*
*              W := W * V2
*
             CALL ZTRMM_Q( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1
*
               CALL ZGEMM_Q( 'Conjugate transpose', 'No transpose', N,
     $                        K, M-K, ONE, C, LDC, V, LDV, ONE, WORK,
     $                        LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL ZTRMM_Q( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1 * W'
*
                  CALL ZGEMM_Q( 'No transpose', 'Conjugate transpose',
     $                        M-K, N, K, -ONE, V, LDV, WORK, LDWORK,
     $                        ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL ZTRMM_Q( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V( M-K+1, 1 ), LDV, WORK,
     $                     LDWORK )
*
*              C2 := C2 - W'
*
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) -
     $                               conjg( WORK( I, J ) )
   80             CONTINUE
   90          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C2
*
               DO 100 J = 1, K
                  CALL ZCOPY_Q( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
*
*              W := W * V2
*
             CALL ZTRMM_Q( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1
*
              CALL ZGEMM_Q( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL ZTRMM_Q( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1'
*
               CALL ZGEMM_Q( 'No transpose', 'Conjugate transpose', M,
     $                        N-K, K, -ONE, WORK, LDWORK, V, LDV, ONE,
     $                        C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL ZTRMM_Q( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V( N-K+1, 1 ), LDV, WORK,
     $                     LDWORK )
*
*              C2 := C2 - W
*
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
*
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1  V2 )    (V1: first K columns)
*           where  V1  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C1'
*
               DO 130 J = 1, K
                  CALL ZCOPY_Q( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV_Q( N, WORK( 1, J ), 1 )
  130          CONTINUE
*
*              W := W * V1'
*
               CALL ZTRMM_Q( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2'
*
                  CALL ZGEMM_Q( 'Conjugate transpose',
     $                        'Conjugate transpose', N, K, M-K, ONE,
     $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
     $                        WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL ZTRMM_Q( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2' * W'
*
                  CALL ZGEMM_Q( 'Conjugate transpose',
     $                        'Conjugate transpose', M-K, N, K, -ONE,
     $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1
*
             CALL ZTRMM_Q( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - conjg( WORK( I, J ) )
  140             CONTINUE
  150          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C1
*
               DO 160 J = 1, K
                  CALL ZCOPY_Q( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
*
*              W := W * V1'
*
               CALL ZTRMM_Q( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2'
*
              CALL ZGEMM_Q( 'No transpose', 'Conjugate transpose', M,
     $                        K, N-K, ONE, C( 1, K+1 ), LDC,
     $                        V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL ZTRMM_Q( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2
*
            CALL ZGEMM_Q( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1
*
            CALL ZTRMM_Q( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
*
            END IF
*
         ELSE
*
*           Let  V =  ( V1  V2 )    (V2: last K columns)
*           where  V2  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C2'
*
               DO 190 J = 1, K
                  CALL ZCOPY_Q( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV_Q( N, WORK( 1, J ), 1 )
  190          CONTINUE
*
*              W := W * V2'
*
               CALL ZTRMM_Q( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V( 1, M-K+1 ), LDV, WORK,
     $                     LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1'
*
                  CALL ZGEMM_Q( 'Conjugate transpose',
     $                        'Conjugate transpose', N, K, M-K, ONE, C,
     $                        LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL ZTRMM_Q( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1' * W'
*
                  CALL ZGEMM_Q( 'Conjugate transpose',
     $                        'Conjugate transpose', M-K, N, K, -ONE, V,
     $                        LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
            CALL ZTRMM_Q( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) -
     $                               conjg( WORK( I, J ) )
  200             CONTINUE
  210          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C2
*
               DO 220 J = 1, K
                  CALL ZCOPY_Q( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
*
*              W := W * V2'
*
               CALL ZTRMM_Q( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V( 1, N-K+1 ), LDV, WORK,
     $                     LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1'
*
               CALL ZGEMM_Q( 'No transpose', 'Conjugate transpose', M,
     $                        K, N-K, ONE, C, LDC, V, LDV, ONE, WORK,
     $                        LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL ZTRMM_Q( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1
*
              CALL ZGEMM_Q( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
              CALL ZTRMM_Q( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of ZLARFB
*
      END

      SUBROUTINE ZLARFT_Q( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
      COMPLEX*32         T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARFT forms the triangular factor T of a complex block reflector H
*  of order n, which is defined as a product of k elementary reflectors.
*
*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*
*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*
*  If STOREV = 'C', the vector which defines the elementary reflector
*  H(i) is stored in the i-th column of the array V, and
*
*     H  =  I - V * T * V'
*
*  If STOREV = 'R', the vector which defines the elementary reflector
*  H(i) is stored in the i-th row of the array V, and
*
*     H  =  I - V' * T * V
*
*  Arguments
*  =========
*
*  DIRECT  (input) CHARACTER*1
*          Specifies the order in which the elementary reflectors are
*          multiplied to form the block reflector:
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Specifies how the vectors which define the elementary
*          reflectors are stored (see also Further Details):
*          = 'C': columnwise
*          = 'R': rowwise
*
*  N       (input) INTEGER
*          The order of the block reflector H. N >= 0.
*
*  K       (input) INTEGER
*          The order of the triangular factor T (= the number of
*          elementary reflectors). K >= 1.
*
*  V       (input/output) COMPLEX*16 array, dimension
*                               (LDV,K) if STOREV = 'C'
*                               (LDV,N) if STOREV = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i).
*
*  T       (output) COMPLEX*16 array, dimension (LDT,K)
*          The k by k triangular factor T of the block reflector.
*          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*          lower triangular. The rest of the array is not used.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  Further Details
*  ===============
*
*  The shape of the matrix V and the storage of the vectors which define
*  the H(i) is best illustrated by the following example with n = 5 and
*  k = 3. The elements equal to 1 are not stored; the corresponding
*  array elements are modified but restored on exit. The rest of the
*  array is not used.
*
*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*
*               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*                   ( v1  1    )                     (     1 v2 v2 v2 )
*                   ( v1 v2  1 )                     (        1 v3 v3 )
*                   ( v1 v2 v3 )
*                   ( v1 v2 v3 )
*
*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*
*               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*                   (     1 v3 )
*                   (        1 )
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*32         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0Q+0, 0.0Q+0 ),
     $                   ZERO = ( 0.0Q+0, 0.0Q+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      COMPLEX*32         VII
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMV_Q, ZLACGV_Q, ZTRMV_Q
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( LSAME( DIRECT, 'F' ) ) THEN
         DO 20 I = 1, K
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
*
*              general case
*
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
*
*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
*
                  CALL ZGEMV_Q( 'Conjugate transpose', N-I+1, I-1,
     $                        -TAU( I ), V( I, 1 ), LDV, V( I, I ), 1,
     $                        ZERO, T( 1, I ), 1 )
               ELSE
*
*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
*
                  IF( I.LT.N )
     $               CALL ZLACGV_Q( N-I, V( I, I+1 ), LDV )
                  CALL ZGEMV_Q( 'No transpose', I-1, N-I+1, -TAU( I ),
     $                        V( 1, I ), LDV, V( I, I ), LDV, ZERO,
     $                        T( 1, I ), 1 )
                  IF( I.LT.N )
     $               CALL ZLACGV_Q( N-I, V( I, I+1 ), LDV )
               END IF
               V( I, I ) = VII
*
*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
*
             CALL ZTRMV_Q( 'Upper', 'No transpose', 'Non-unit', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
*
*              general case
*
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
*
                     CALL ZGEMV_Q( 'Conjugate transpose', N-K+I, K-I,
     $                           -TAU( I ), V( 1, I+1 ), LDV, V( 1, I ),
     $                           1, ZERO, T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
*
                     CALL ZLACGV_Q( N-K+I-1, V( I, 1 ), LDV )
                  CALL ZGEMV_Q( 'No transpose', K-I, N-K+I, -TAU( I ),
     $                           V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO,
     $                           T( I+1, I ), 1 )
                     CALL ZLACGV_Q( N-K+I-1, V( I, 1 ), LDV )
                     V( I, N-K+I ) = VII
                  END IF
*
*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
*
            CALL ZTRMV_Q( 'Lower', 'No transpose', 'Non-unit', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
*
*     End of ZLARFT
*
      END


      SUBROUTINE ZLARF_Q( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*32         TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*32         C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARF applies a complex elementary reflector H to a complex M-by-N
*  matrix C, from either the left or the right. H is represented in the
*  form
*
*        H = I - tau * v * v'
*
*  where tau is a complex scalar and v is a complex vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*  To apply H' (the conjugate transpose of H), supply conjg(tau) instead
*  tau.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  V       (input) COMPLEX*16 array, dimension
*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*          The vector v in the representation of H. V is not used if
*          TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) COMPLEX*16
*          The value tau in the representation of H.
*
*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX*16 array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*32         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0Q+0, 0.0Q+0 ),
     $                   ZERO = ( 0.0Q+0, 0.0Q+0 ) )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMV_Q, ZGERC_Q
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  H * C
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C' * v
*
            CALL ZGEMV_Q( 'Conjugate transpose', M, N, ONE, C, LDC, V,
     $                  INCV, ZERO, WORK, 1 )
*
*           C := C - v * w'
*
            CALL ZGERC_Q( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C * v
*
            CALL ZGEMV_Q( 'No transpose', M, N, ONE, C, LDC, V, INCV,
     $                  ZERO, WORK, 1 )
*
*           C := C - w * v'
*
            CALL ZGERC_Q( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
*
*     End of ZLARF
*
      END


      SUBROUTINE ZLARFG_Q( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*32         ALPHA, TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*32         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARFG generates a complex elementary reflector H of order n, such
*  that
*
*        H' * ( alpha ) = ( beta ),   H' * H = I.
*             (   x   )   (   0  )
*
*  where alpha and beta are scalars, with beta real, and x is an
*  (n-1)-element complex vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a complex scalar and v is a complex (n-1)-element
*  vector. Note that H is not hermitian.
*
*  If the elements of x are all zero and alpha is real, then tau = 0
*  and H is taken to be the unit matrix.
*
*  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) COMPLEX*16
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) COMPLEX*16 array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) COMPLEX*16
*          The value tau.
*
*  =====================================================================
*
*     .. Parameters ..
      real*16   ONE, ZERO
      PARAMETER          ( ONE = 1.0Q+0, ZERO = 0.0Q+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      real*16   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      complex*32      ZLADIV_Q
      real*16 DLAMCH_Q, DLAPY3_Q, DZNRM2_Q
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, realpart, CMPLX, imagpart, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZDSCAL_Q, ZSCAL_Q
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = DZNRM2_Q( N-1, X, INCX )
      ALPHR = realpart( ALPHA )
      ALPHI = imagpart( ALPHA )
*
      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( DLAPY3_Q( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = DLAMCH_Q( 'S' ) / DLAMCH_Q( 'E' )
         RSAFMN = ONE / SAFMIN
*
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL ZDSCAL_Q( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = DZNRM2_Q( N-1, X, INCX )
            ALPHA = CMPLX( ALPHR, ALPHI,16 )
            BETA = -SIGN( DLAPY3_Q( ALPHR, ALPHI, XNORM ), ALPHR )
            TAU = CMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA,16 )
            ALPHA = ZLADIV_Q( CMPLX( ONE,zero,16 ), ALPHA-BETA )
            CALL ZSCAL_Q( N-1, ALPHA, X, INCX )
*
*           If ALPHA is subnormal, it may lose relative accuracy
*
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = CMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA,16 )
            ALPHA = ZLADIV_Q( CMPLX( ONE,zero,16 ), ALPHA-BETA )
            CALL ZSCAL_Q( N-1, ALPHA, X, INCX )
            ALPHA = BETA
         END IF
      END IF
*
      RETURN
*
*     End of ZLARFG
*
      END


      SUBROUTINE ZCOPY_Q(N,ZX,INCX,ZY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      complex*32 ZX(*),ZY(*)
*     ..
*
*  Purpose
*  =======
*
*     copies a vector, x, to a vector, y.
*     jack dongarra, linpack, 4/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX,IY
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          ZY(IY) = ZX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
   20 DO 30 I = 1,N
          ZY(I) = ZX(I)
   30 CONTINUE
      RETURN
      END

      SUBROUTINE ZGEMM_Q(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA
     & ,C,LDC)
*     .. Scalar Arguments ..
      complex*32 ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      complex*32 A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  Purpose
*  =======
*
*  ZGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC CONJG,MAX
*     ..
*     .. Local Scalars ..
      complex*32 TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL CONJA,CONJB,NOTA,NOTB
*     ..
*     .. Parameters ..
      complex*32 ONE
      PARAMETER (ONE= (1.0Q+0,0.0Q+0))
      complex*32 ZERO
      PARAMETER (ZERO= (0.0Q+0,0.0Q+0))
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
*     B  respectively are to be  transposed but  not conjugated  and set
*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
*     and the number of rows of  B  respectively.
*
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      CONJA = LSAME(TRANSA,'C')
      CONJB = LSAME(TRANSB,'C')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.CONJA) .AND.
     +    (.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.CONJB) .AND.
     +         (.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGEMM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (NOTB) THEN
          IF (NOTA) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
              DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      IF (B(L,J).NE.ZERO) THEN
                          TEMP = ALPHA*B(L,J)
                          DO 70 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
   70                     CONTINUE
                      END IF
   80             CONTINUE
   90         CONTINUE
          ELSE IF (CONJA) THEN
*
*           Form  C := alpha*conjg( A' )*B + beta*C.
*
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + conjg(A(L,I))*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
              DO 150 J = 1,N
                  DO 140 I = 1,M
                      TEMP = ZERO
                      DO 130 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  130                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  140             CONTINUE
  150         CONTINUE
          END IF
      ELSE IF (NOTA) THEN
          IF (CONJB) THEN
*
*           Form  C := alpha*A*conjg( B' ) + beta*C.
*
              DO 200 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 160 I = 1,M
                          C(I,J) = ZERO
  160                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 170 I = 1,M
                          C(I,J) = BETA*C(I,J)
  170                 CONTINUE
                  END IF
                  DO 190 L = 1,K
                      IF (B(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*conjg(B(J,L))
                          DO 180 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  180                     CONTINUE
                      END IF
  190             CONTINUE
  200         CONTINUE
          ELSE
*
*           Form  C := alpha*A*B'          + beta*C
*
              DO 250 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 210 I = 1,M
                          C(I,J) = ZERO
  210                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 220 I = 1,M
                          C(I,J) = BETA*C(I,J)
  220                 CONTINUE
                  END IF
                  DO 240 L = 1,K
                      IF (B(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*B(J,L)
                          DO 230 I = 1,M
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  230                     CONTINUE
                      END IF
  240             CONTINUE
  250         CONTINUE
          END IF
      ELSE IF (CONJA) THEN
          IF (CONJB) THEN
*
*           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
*
              DO 280 J = 1,N
                  DO 270 I = 1,M
                      TEMP = ZERO
                      DO 260 L = 1,K
                          TEMP = TEMP + conjg(A(L,I))*conjg(B(J,L))
  260                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  270             CONTINUE
  280         CONTINUE
          ELSE
*
*           Form  C := alpha*conjg( A' )*B' + beta*C
*
              DO 310 J = 1,N
                  DO 300 I = 1,M
                      TEMP = ZERO
                      DO 290 L = 1,K
                          TEMP = TEMP + conjg(A(L,I))*B(J,L)
  290                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  300             CONTINUE
  310         CONTINUE
          END IF
      ELSE
          IF (CONJB) THEN
*
*           Form  C := alpha*A'*conjg( B' ) + beta*C
*
              DO 340 J = 1,N
                  DO 330 I = 1,M
                      TEMP = ZERO
                      DO 320 L = 1,K
                          TEMP = TEMP + A(L,I)*conjg(B(J,L))
  320                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  330             CONTINUE
  340         CONTINUE
          ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
              DO 370 J = 1,N
                  DO 360 I = 1,M
                      TEMP = ZERO
                      DO 350 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  350                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  360             CONTINUE
  370         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of ZGEMM .
*
      END

      SUBROUTINE ZLACGV_Q( N, X, INCX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
*     ..
*     .. Array Arguments ..
      COMPLEX*32         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLACGV conjugates a complex vector of length N.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The length of the vector X.  N >= 0.
*
*  X       (input/output) COMPLEX*16 array, dimension
*                         (1+(N-1)*abs(INCX))
*          On entry, the vector of length N to be conjugated.
*          On exit, X is overwritten with conjg(X).
*
*  INCX    (input) INTEGER
*          The spacing between successive elements of X.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IOFF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG
*     ..
*     .. Executable Statements ..
*
      IF( INCX.EQ.1 ) THEN
         DO 10 I = 1, N
            X( I ) = CONJG( X( I ) )
   10    CONTINUE
      ELSE
         IOFF = 1
         IF( INCX.LT.0 )
     $      IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = CONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZLACGV
*
      END


      SUBROUTINE ZTRMM_Q(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*     .. Scalar Arguments ..
      complex*32 ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      complex*32 A(LDA,*),B(LDB,*)
*     ..
*
*  Purpose
*  =======
*
*  ZTRMM  performs one of the matrix-matrix operations
*
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
*
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
*
*  Arguments
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies B from
*           the left or right as follows:
*
*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*
*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain the matrix  B,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC CONJG,MAX
*     ..
*     .. Local Scalars ..
      complex*32 TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      complex*32 ONE
      PARAMETER (ONE= (1.0Q+0,0.0Q+0))
      complex*32 ZERO
      PARAMETER (ZERO= (0.0Q+0,0.0Q+0))
*     ..
*
*     Test the input parameters.
*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOCONJ = LSAME(TRANSA,'T')
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRMM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*A*B.
*
              IF (UPPER) THEN
                  DO 50 J = 1,N
                      DO 40 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              DO 30 I = 1,K - 1
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   30                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP*A(K,K)
                              B(K,J) = TEMP
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              B(K,J) = TEMP
                              IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                              DO 60 I = K + 1,M
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*A'*B   or   B := alpha*conjg( A' )*B.
*
              IF (UPPER) THEN
                  DO 120 J = 1,N
                      DO 110 I = M,1,-1
                          TEMP = B(I,J)
                          IF (NOCONJ) THEN
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              DO 90 K = 1,I - 1
                                  TEMP = TEMP + A(K,I)*B(K,J)
   90                         CONTINUE
                          ELSE
                              IF (NOUNIT) TEMP = TEMP*conjg(A(I,I))
                              DO 100 K = 1,I - 1
                                  TEMP = TEMP + conjg(A(K,I))*B(K,J)
  100                         CONTINUE
                          END IF
                          B(I,J) = ALPHA*TEMP
  110                 CONTINUE
  120             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = 1,M
                          TEMP = B(I,J)
                          IF (NOCONJ) THEN
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              DO 130 K = I + 1,M
                                  TEMP = TEMP + A(K,I)*B(K,J)
  130                         CONTINUE
                          ELSE
                              IF (NOUNIT) TEMP = TEMP*conjg(A(I,I))
                              DO 140 K = I + 1,M
                                  TEMP = TEMP + conjg(A(K,I))*B(K,J)
  140                         CONTINUE
                          END IF
                          B(I,J) = ALPHA*TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*B*A.
*
              IF (UPPER) THEN
                  DO 200 J = N,1,-1
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 170 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  170                 CONTINUE
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
  200             CONTINUE
              ELSE
                  DO 240 J = 1,N
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 210 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  210                 CONTINUE
                      DO 230 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 220 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  220                         CONTINUE
                          END IF
  230                 CONTINUE
  240             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
*
              IF (UPPER) THEN
                  DO 280 K = 1,N
                      DO 260 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = ALPHA*A(J,K)
                              ELSE
                                  TEMP = ALPHA*conjg(A(J,K))
                              END IF
                              DO 250 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = TEMP*A(K,K)
                          ELSE
                              TEMP = TEMP*conjg(A(K,K))
                          END IF
                      END IF
                      IF (TEMP.NE.ONE) THEN
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              ELSE
                  DO 320 K = N,1,-1
                      DO 300 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = ALPHA*A(J,K)
                              ELSE
                                  TEMP = ALPHA*conjg(A(J,K))
                              END IF
                              DO 290 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  290                         CONTINUE
                          END IF
  300                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = TEMP*A(K,K)
                          ELSE
                              TEMP = TEMP*conjg(A(K,K))
                          END IF
                      END IF
                      IF (TEMP.NE.ONE) THEN
                          DO 310 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  310                     CONTINUE
                      END IF
  320             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRMM .
*
      END

      SUBROUTINE ZGEMV_Q(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*     .. Scalar Arguments ..
      complex*32 ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
*     ..
*     .. Array Arguments ..
      complex*32 A(LDA,*),X(*),Y(*)
*     ..
*
*  Purpose
*  =======
*
*  ZGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
*
*     y := alpha*conjg( A' )*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      complex*32 ONE
      PARAMETER (ONE= (1.0Q+0,0.0Q+0))
      complex*32 ZERO
      PARAMETER (ZERO= (0.0Q+0,0.0Q+0))
*     ..
*     .. Local Scalars ..
      complex*32 TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
      LOGICAL NOCONJ
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC CONJG,MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +    .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGEMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
*
      NOCONJ = LSAME(TRANS,'T')
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      DO 50 I = 1,M
                          Y(I) = Y(I) + TEMP*A(I,J)
   50                 CONTINUE
                  END IF
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IY = KY
                      DO 70 I = 1,M
                          Y(IY) = Y(IY) + TEMP*A(I,J)
                          IY = IY + INCY
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
*
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 110 J = 1,N
                  TEMP = ZERO
                  IF (NOCONJ) THEN
                      DO 90 I = 1,M
                          TEMP = TEMP + A(I,J)*X(I)
   90                 CONTINUE
                  ELSE
                      DO 100 I = 1,M
                          TEMP = TEMP + conjg(A(I,J))*X(I)
  100                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  110         CONTINUE
          ELSE
              DO 140 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  IF (NOCONJ) THEN
                      DO 120 I = 1,M
                          TEMP = TEMP + A(I,J)*X(IX)
                          IX = IX + INCX
  120                 CONTINUE
                  ELSE
                      DO 130 I = 1,M
                          TEMP = TEMP + conjg(A(I,J))*X(IX)
                          IX = IX + INCX
  130                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  140         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of ZGEMV .
*
      END

      SUBROUTINE ZLACG_Q( N, X, INCX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
*     ..
*     .. Array Arguments ..
      COMPLEX*32         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLACGV conjugates a complex vector of length N.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The length of the vector X.  N >= 0.
*
*  X       (input/output) COMPLEX*16 array, dimension
*                         (1+(N-1)*abs(INCX))
*          On entry, the vector of length N to be conjugated.
*          On exit, X is overwritten with conjg(X).
*
*  INCX    (input) INTEGER
*          The spacing between successive elements of X.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IOFF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG
*     ..
*     .. Executable Statements ..
*
      IF( INCX.EQ.1 ) THEN
         DO 10 I = 1, N
            X( I ) = CONJG( X( I ) )
   10    CONTINUE
      ELSE
         IOFF = 1
         IF( INCX.LT.0 )
     $      IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = CONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZLACGV
*
      END

      SUBROUTINE ZTRMV_Q(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*     ..
*     .. Array Arguments ..
      complex*32 A(LDA,*),X(*)
*     ..
*
*  Purpose
*  =======
*
*  ZTRMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x,   or   x := conjg( A' )*x,
*
*  where x is an n element vector and  A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Arguments
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   x := A*x.
*
*              TRANS = 'T' or 't'   x := A'*x.
*
*              TRANS = 'C' or 'c'   x := conjg( A' )*x.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      complex*32 ZERO
      PARAMETER (ZERO= (0.0Q+0,0.0Q+0))
*     ..
*     .. Local Scalars ..
      complex*32 TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC CONJG,MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0) RETURN
*
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x := A*x.
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 10 I = 1,J - 1
                              X(I) = X(I) + TEMP*A(I,J)
   10                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX
                  DO 40 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 30 I = 1,J - 1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX + INCX
   30                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX + INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*A(I,J)
   50                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   60             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 70 I = N,J + 1,-1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX - INCX
   70                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX - INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
*
*        Form  x := A'*x  or  x := conjg( A' )*x.
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 110 J = N,1,-1
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 90 I = J - 1,1,-1
                              TEMP = TEMP + A(I,J)*X(I)
   90                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*conjg(A(J,J))
                          DO 100 I = J - 1,1,-1
                              TEMP = TEMP + conjg(A(I,J))*X(I)
  100                     CONTINUE
                      END IF
                      X(J) = TEMP
  110             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 140 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 120 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + A(I,J)*X(IX)
  120                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*conjg(A(J,J))
                          DO 130 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + conjg(A(I,J))*X(IX)
  130                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
  140             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 170 J = 1,N
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 150 I = J + 1,N
                              TEMP = TEMP + A(I,J)*X(I)
  150                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*conjg(A(J,J))
                          DO 160 I = J + 1,N
                              TEMP = TEMP + conjg(A(I,J))*X(I)
  160                     CONTINUE
                      END IF
                      X(J) = TEMP
  170             CONTINUE
              ELSE
                  JX = KX
                  DO 200 J = 1,N
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 180 I = J + 1,N
                              IX = IX + INCX
                              TEMP = TEMP + A(I,J)*X(IX)
  180                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*conjg(A(J,J))
                          DO 190 I = J + 1,N
                              IX = IX + INCX
                              TEMP = TEMP + conjg(A(I,J))*X(IX)
  190                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
  200             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRMV .
*
      END

      SUBROUTINE ZGERC_Q(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*     .. Scalar Arguments ..
      complex*32 ALPHA
      INTEGER INCX,INCY,LDA,M,N
*     ..
*     .. Array Arguments ..
      complex*32 A(LDA,*),X(*),Y(*)
*     ..
*
*  Purpose
*  =======
*
*  ZGERC  performs the rank 1 operation
*
*     A := alpha*x*conjg( y' ) + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Arguments
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      complex*32 ZERO
      PARAMETER (ZERO= (0.0Q+0,0.0Q+0))
*     ..
*     .. Local Scalars ..
      complex*32 TEMP
      INTEGER I,INFO,IX,J,JY,KX
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC CONJG,MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGERC ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*conjg(Y(JY))
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*conjg(Y(JY))
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
*
      RETURN
*
*     End of ZGERC .
*
      END


*  =====================================================================
      real*16 FUNCTION DLAMCH_Q( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*     .. Scalar Arguments ..
      real*16   A, B
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      real*16   ONE, ZERO
      PARAMETER          ( ONE = 1.0Q+0, ZERO = 0.0Q+0 )
*     ..
*     .. Local Scalars ..
      real*16   RND, EPS, SFMIN, SMALL, RMACH
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT,
     $                   MINEXPONENT, RADIX, TINY
*     ..
*     .. Executable Statements ..
*
*
*     Assume rounding, not chopping. Always.
*
      RND = ONE
*
      IF( ONE.EQ.RND ) THEN
         EPS = EPSILON(ZERO) * 0.5
      ELSE
         EPS = EPSILON(ZERO)
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         SFMIN = TINY(ZERO)
         SMALL = ONE / HUGE(ZERO)
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = EPS * RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = DIGITS(ZERO)
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = MINEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = tiny(zero)
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = MAXEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = HUGE(ZERO)
      ELSE
         RMACH = ZERO
      END IF
*
      DLAMCH_Q = RMACH
      RETURN
*
*     End of DLAMCH
*
      END

      real*16 FUNCTION DLAPY3_Q( X, Y, Z )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      real*16   X, Y, Z
*     ..
*
*  Purpose
*  =======
*
*  DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
*  unnecessary overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*  Z       (input) DOUBLE PRECISION
*          X, Y and Z specify the values x, y and z.
*
*  =====================================================================
*
*     .. Parameters ..
      real*16   ZERO
      PARAMETER          ( ZERO = 0.0Q+0 )
*     ..
*     .. Local Scalars ..
      real*16   W, XABS, YABS, ZABS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
*     W can be zero for max(0,nan,0)
*     adding all three entries together will make sure
*     NaN will not disappear.
         DLAPY3_Q =  XABS + YABS + ZABS
      ELSE
         DLAPY3_Q = W*SQRT( ( XABS / W )**2+( YABS / W )**2+
     $            ( ZABS / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY3
*
      END

      real*16 FUNCTION DZNRM2_Q(N,X,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      complex*32 X(*)
*     ..
*
*  Purpose
*  =======
*
*  DZNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DZNRM2 := sqrt( conjg( x' )*x )
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to ZLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
      real*16 ONE,ZERO
      PARAMETER (ONE=1.0Q+0,ZERO=0.0Q+0)
*     ..
*     .. Local Scalars ..
      real*16 NORM,SCALE,SSQ,TEMP
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,realpart,IMAGpart,SQRT
*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE
          SCALE = ZERO
          SSQ = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (realpart(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(realpart(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
              IF (imagpart(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(imagpart(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
*
      DZNRM2_Q = NORM
      RETURN
*
*     End of DZNRM2.
*
      END


      COMPLEX*32     FUNCTION ZLADIV_Q( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      COMPLEX*32         X, Y
*     ..
*
*  Purpose
*  =======
*
*  ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y
*  will not overflow on an intermediary step unless the results
*  overflows.
*
*  Arguments
*  =========
*
*  X       (input) COMPLEX*16
*  Y       (input) COMPLEX*16
*          The complex scalars X and Y.
*
*  =====================================================================
*
*     .. Local Scalars ..
      real*16   ZI, ZR
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          realpart, CMPLX, imagpart
*     ..
*     .. Executable Statements ..
*
      CALL DLADIV_Q( realpart(X), imagpart(X), realpart(Y), imagpart(Y), 
     $             ZR,ZI )
      ZLADIV_Q = CMPLX( ZR, ZI,16 )
*
      RETURN
*
*     End of ZLADIV
*
      END


      SUBROUTINE ZDSCAL_Q(N,DA,ZX,INCX)
*     .. Scalar Arguments ..
      real*16 DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      complex*32 ZX(*)
*     ..
*
*  Purpose
*  =======
*
*     scales a vector by a constant.
*     jack dongarra, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC CMPLX
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      IX = 1
      DO 10 I = 1,N
          ZX(IX) = CMPLX(DA,0.0q0,16)*ZX(IX)
          IX = IX + INCX
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
   20 DO 30 I = 1,N
          ZX(I) = CMPLX(DA,0.0q0,16)*ZX(I)
   30 CONTINUE
      RETURN
      END


      SUBROUTINE ZSCAL_Q(N,ZA,ZX,INCX)
*     .. Scalar Arguments ..
      complex*32 ZA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      complex*32 ZX(*)
*     ..
*
*  Purpose
*  =======
*
*     scales a vector by a constant.
*     jack dongarra, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      IX = 1
      DO 10 I = 1,N
          ZX(IX) = ZA*ZX(IX)
          IX = IX + INCX
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
   20 DO 30 I = 1,N
          ZX(I) = ZA*ZX(I)
   30 CONTINUE
      RETURN
      END


      SUBROUTINE DLADIV_Q( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      real*16   A, B, C, D, P, Q
*     ..
*
*  Purpose
*  =======
*
*  DLADIV performs complex division in  real arithmetic
*
*                        a + i*b
*             p + i*q = ---------
*                        c + i*d
*
*  The algorithm is due to Robert L. Smith and can be found
*  in D. Knuth, The art of Computer Programming, Vol.2, p.195
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*  C       (input) DOUBLE PRECISION
*  D       (input) DOUBLE PRECISION
*          The scalars a, b, c, and d in the above expression.
*
*  P       (output) DOUBLE PRECISION
*  Q       (output) DOUBLE PRECISION
*          The scalars p and q in the above expression.
*
*  =====================================================================
*
*     .. Local Scalars ..
      real*16   E, F
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( ABS( D ).LT.ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
*
      RETURN
*
*     End of DLADIV
*
      END


      SUBROUTINE ZTRSM_Q (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX*32         ALPHA
*     .. Array Arguments ..
      COMPLEX*32         A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  ZTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX*32         TEMP
*     .. Parameters ..
      COMPLEX*32         ONE
      PARAMETER        ( ONE  = ( 1.0Q+0, 0.0Q+0 ) )
      COMPLEX*32         ZERO
      PARAMETER        ( ZERO = ( 0.0Q+0, 0.0Q+0 ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZTRSM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*inv( A' )*B
*           or    B := alpha*inv( conjg( A' ) )*B.
*
            IF( UPPER )THEN
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 110, K = 1, I - 1
                           TEMP = TEMP - A( K, I )*B( K, J )
  110                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 120, K = 1, I - 1
                           TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  120                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/CONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  130             CONTINUE
  140          CONTINUE
            ELSE
               DO 180, J = 1, N
                  DO 170, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 150, K = I + 1, M
                           TEMP = TEMP - A( K, I )*B( K, J )
  150                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 160, K = I + 1, M
                           TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  160                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/CONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  170             CONTINUE
  180          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF( UPPER )THEN
               DO 230, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 190, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  190                CONTINUE
                  END IF
                  DO 210, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 220, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  220                CONTINUE
                  END IF
  230          CONTINUE
            ELSE
               DO 280, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 240, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  240                CONTINUE
                  END IF
                  DO 260, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 270, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  270                CONTINUE
                  END IF
  280          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*inv( A' )
*           or    B := alpha*B*inv( conjg( A' ) ).
*
            IF( UPPER )THEN
               DO 330, K = N, 1, -1
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/CONJG( A( K, K ) )
                     END IF
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
                  DO 310, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = CONJG( A( J, K ) )
                        END IF
                        DO 300, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  300                   CONTINUE
                     END IF
  310             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 320, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  320                CONTINUE
                  END IF
  330          CONTINUE
            ELSE
               DO 380, K = 1, N
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/CONJG( A( K, K ) )
                     END IF
                     DO 340, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  340                CONTINUE
                  END IF
                  DO 360, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = CONJG( A( J, K ) )
                        END IF
                        DO 350, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  350                   CONTINUE
                     END IF
  360             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 370, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  370                CONTINUE
                  END IF
  380          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of ZTRSM .
*
      END


      subroutine nutolsma(nch,s,yt,tolsma,h,ir)

      implicit none
      integer:: nch
      complex*16:: s(nch),yt(nch,nch,5),f(2*nch,nch)
      complex*16:: y0(nch,nch),yp0(nch,nch)
      real*8:: tolsma,norm,h
      integer:: ich,is,ir  
      real*8 :: svd(nch)
      complex*16 :: U(2*nch,2*nch),VT(nch,nch)
    
    ! LAPACK Variables
    ! WORK is Complex, RWORK is Real
      complex(8), allocatable :: WORK(:)
      real(8), allocatable :: RWORK(:)
      integer :: LWORK, INFO

!First build derivative and function
!      write(448,'(i5,2x,9999e14.6)') ir, s(1:nch)
!      write(449,'(i5,2x,9999e14.6)') ir, yt(1:nch,1:nch,1:5)
      do is=1,nch
      do ich=1,nch
            y0(ich,is)=yt(ich,is,5)
      yp0(ich,is)= (3d0*yt(ich,is,1)-16d0*yt(ich,is,2)+36d0*yt(ich,is,3)
     & -48d0*yt(ich,is,4)+25d0*yt(ich,is,5))/(12d0*h)
      enddo
      enddo

! Next build f matrix from Eq. (4.31)
      f=0d0
      do ich=1,nch
      do is=1,nch
            f(ich,is)=y0(ich,is)
            f(nch+ich,is)=-(0d0,1d0)*(s(ich))**(-0.5d0)*yp0(ich,is)
      enddo
      enddo

!      do ich=1,2*nch
!      write(447,'(9999g14.6)') ( f(ich,is), is=1,nch)
!      enddo

!Now normalize columns of f matrix
      do is=1,nch
      norm=0d0
      do ich=1,2*nch
            norm=norm+abs(f(ich,is))**2
      enddo
      if (norm.eq.0d0) cycle
      norm=sqrt(norm)
      do ich=1,2*nch
            f(ich,is)=f(ich,is)/norm
      enddo
      enddo
      


!Now compute single-value decomposition of f matrix using lapack 
      allocate(RWORK(5*nch))
 
    ! Query WORK size
      allocate(WORK(1))
    
    ! Call ZGESVD (Complex Double Precision SVD)
      call ZGESVD('A', 'A', 2*nch, nch, f, 2*nch, svd, U, 2*nch, VT, 
     & nch, WORK, -1, RWORK, INFO)
    
    ! WORK(1) is complex; take the real part to get the size
      LWORK = int(real(WORK(1)))
!      write(*,*)'NUTOLSMA: LWORK=',LWORK
      deallocate(WORK)
      allocate(WORK(LWORK))


      call ZGESVD('A', 'A', 2*nch, nch, f, 2*nch, svd, U, 2*nch, VT,
     &  nch, WORK, LWORK, RWORK, INFO)          

      if (INFO.ne.0) then
        write(*,*)'NUTOLSMA: Error in SVD calculation, INFO=',INFO
        tolsma=1e-10!Force orthogonalization
        return
      endif

      tolsma=svd(nch)/svd(1)

      end subroutine nutolsma



      subroutine schcc_ena_tolsma(nch,ecm,z12,incvec,ql,factor,dr,r0,
     & npt,wf,phase,smat,method,info,einc,icc)
      use xcdcc, only: smats
      use nmrv,only: nr,h,vcoup,conv,ech,debug,rmin,hort,rvec,rmort,cutr
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
      complex*16:: ytolsma(nch,nch,5),s(nch)
      real*8:: tolsma
      integer jr,count_tolsma
      complex*16:: qryp(nch,nch),qry(nch,nch,npt),sqr(nch)
      integer::  permute(nch)
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

      write(*,*) 'cutr',cutr      
      if (allocated(rvec)) deallocate(rvec)
      allocate(rvec(nr))
      do ir=1,nr
        rvec(ir)=r0+(ir-1)*dr
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius

! CHanged by AMoro to prevent undefined mod(ir-irmin,norto) below
!      norto=nint(hort/h)
      norto=max(nint(hort/h),1)

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
      rmin=max(h,rturn-abs(cutr)) ! CHECK h or r0 ?
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

      if (verb.ge.2) 
     & write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
     &  5x,"outer=",1f8.2," fm",
     &  5x,"R-min=",1f6.2," fm")') 
     &  minval(rturnv(:)), maxval(rturnv(:)), rmin      

c 
c Initial values: w(rmin) and w(rmin+h) 
c 
      if (irmin.lt.2) write(*,*)'wrong irmin=',irmin
      call zfun(nch,ql,kch2,irmin-1,ym,zm)
      wm=ym-(h2/12d0)*zm
      call zfun(nch,ql,kch2,irmin  ,y0,z0)
      w0=y0-(h2/12d0)*z0
 
      count_tolsma=0
      do ir=irmin,nr-1  
      r=rvec(ir) ! (ir-1)*h    
      im=ir-1  ! r-h
      i0=ir    ! r
      ip=ir+1  ! r+h
      call csh(nch,ql,kch2,i0,w0,z0,1) ! Cosh[Sqrt(sqrt(h2 F(x))]

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
!      if ((ir.ge.nr-8).or.(orto)) then      
        call w2y(wp,yp,ql,kch2,ir,nch)
        y(1:nch,1:nch,ip)=yp(1:nch,1:nch) 
!      endif

!Compute nu factor from Tolsma CPC 40 233 (1986)

      do ich=1,nch    
      l=ql(ich)
      s(ich)= h2*(-l*(l+1)/(r+h)**2         ! centrifugal
     &          + kch2(ich)             ! energy
     &        - conv*vcoup(ich,ich,ip)) ! potential
      enddo


      if (ir.eq. irmin) write(447,*) '#R        Tolsma'
      count_tolsma=count_tolsma+1!Avoid using for the derivative points before and after orthogonalization
      tolsma=1d0
      if (count_tolsma.gt.5) then
      ytolsma(:,:,:) =y(:,:,ip-4:ip)
      call nutolsma(nch,s,ytolsma,tolsma,h,ip)
      write(447,*) r,tolsma
      endif


c Re-orthogonalize solution vectors
      if (orto.and.(r.lt.rmorto).and.(ir.gt.irmin)
     & .and.tolsma.lt.5e-4) then 
        if (verb.ge.4) 
     &   write(*,*) 'orthogonalizing at ir,r=',ir,r

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
       
        sqr=s
      permute=0
      do ich=1,nch
      is=minloc(dble(sqr),1)
        permute(nch+1-ich)=is
        sqr(is)=1e30
      enddo

!Now perform the permutation

      do ich=1,nch
      do is=1,nch
      qryp(ich,is)=yp(permute(ich),permute(is))
      do jr=irmin,ip
      qry(ich,is,jr)=y(permute(ich),permute(is),jr)
      enddo
      enddo
      enddo

        call qrfac(qryp,qry,nch,ip,irmin,nr)      ! QR factorization
!       call qrfac(wp,w,nch,ip,irmin,nr)      ! QR factorization

      count_tolsma=0


!Undo the permutation

      do ich=1,nch
      do is=1,nch
      yp(permute(ich),permute(is))=qryp(ich,is)
      do jr=irmin,ip
      y(permute(ich),permute(is),jr)=qry(ich,is,jr)
      enddo
      enddo
      enddo

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
      wp=yp-(h2/12d0)*zp
      call zfun(nch,ql,kch2,i0,y0,z0)
      w0=y0-(h2/12d0)*z0


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


c
c wm <- w0 <- wp 
c
      wm(:,:)=w0(:,:) 
      w0(:,:)=wp(:,:) 

      enddo !ir
      if (ir.eq.nr) write(447,*)
      if (ir.eq.nr) write(447,*)
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

      call flush(6)
      end subroutine schcc_ena_tolsma

