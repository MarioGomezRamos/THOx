
      subroutine belam() !,jtot,partot
c ** Calculates B(E;lambda) distributions using THO pseudostates
      use globals !, only: written,mu12,pi,verb
      use constants, only: hc,pi,amu
      use hmatrix, only: hmatx
      use potentials, only: delta,maxlamb
      use sistema
      use wfs !, only: jc,parc,rvec,ql,qjc,qj,iord,cindex !,exc
      use channels
      use belambdamod
      use scattering
      implicit none
c      ------------------------------------------------------------------------
      logical :: energy,fail3!,ifbel
      integer :: i,np,ir,m,n,ncont,iil,is,itran,ik,jset,nex,ich
      integer:: pci,pcf
      integer, parameter:: igs=1
c     --------------------------------------------------------------------------
      real*8:: fival,res,res2,lambdar
      real*8:: Elam,Elamcore,BEl,BElcore,matel,besum,rms
      real*8:: Elamcorem,rotor
      real*8:: wid,qjcir,qjcr,cleb,sixj,coef
      real*8:: mec(maxcore,maxcore,maxlamb) 
      real*8:: econt,kcont,test,rme,geom,excore
      real*8:: jci,jcf,bec,alpha
      real*8:: mvc,mv, rm,conv,solap,rlam,dk
      real*8:: dkrdk,krel   !depends on wfcont and solap definitions
      real*8:: eps,aux
      real*8, allocatable:: edisc(:) ! discrete energies (PS's or Bins)
      real*8, allocatable:: uext(:,:)
      real*8, allocatable:: qlr(:)
      real*8, allocatable:: psh(:)
c     --------------------------------------------------------------------------
      complex*16:: Elamcont,zero,resc,mecc
      complex*16:: gaux(nr),faux(nr)!,test
      complex*16,allocatable:: smate(:,:),wfscat(:,:,:)
!      complex*16,allocatable :: ugs(:,:)
      integer coremodel,qc,ic,icp,nctrans
      CHARACTER*1 PSIGN(3)
      DATA PSIGN / '-','?','+' / 
c     --------------------------------------------------------------------------  
      namelist /belambda/ uwfgsfile,lambda,BElcore,ifbel,rms,
     &  coremodel,jset,emin,emax,nk
      namelist /coretrans/ic, icp, rme, qc  
c     --------------------------------------------------------------------------
      eps=1e-6
      zero=(0d0,0d0)
      BElcore=0d0
      Elamcorem=0d0
      rms=0d0
      coremodel=0
      nctrans=0
      uwfgsfile=""
      energy=.false. 
      jset=0
      read(kin,nml=belambda)
      if (.not.ifbel) return
      allocate(smate(nk,maxchan),psh(nk))
      smate(:,:)=0. ; psh(:)=0
c     ----------------- Core RME specified by user -----------------------------
      if (coremodel.eq.1) then
      mec=0d0
      write(*,*)'o Core reduced matrix elements:'
230   rme=0d0;
      read(kin,nml=coretrans)
      if (rme.eq.0) goto 235
      nctrans=nctrans+1
      jci=jc(ic)
      jcf=jc(icp)
      pci=parc(ic)
      pcf=parc(icp)
!      bec=rme**2/(2d0*jci+1)
!c CHECK THIS non-symmetric M(I->I')
      mec(ic,icp,qc)=rme !/sqrt(2*jcf+1)
      bec=rme**2*(2*jcf+1)/(2d0*jci+1)
      write(*,232)qc,jci,psign(pci+2),jcf,psign(pcf+2),bec
232   format(3x,"B(E",i1,";",1f3.1,a1," -> ",1f3.1,a1,")=",1f7.2)
      goto 230
      endif ! coremodel
c     --------------------------------------------------------------------------
        
235   write(*,*)    
      write(*,*)'  o Calculations of B(Elambda) in progress'

c     ------------------------------ GS WF  -----------------------------------
      if (uwfgsfile.eq."") then
      ncni=jpiset(igs)%nchan
   
      allocate(qjci(ncni),qli(ncni),qlir(ncni),qji(ncni),cindexi(ncni))
      allocate(ugs(nr,ncni))
      ugs(:,:)=0d0
      partoti        =jpiset(igs)%partot
      jtoti          =jpiset(igs)%jtot
      qli (1:ncni)   =jpiset(igs)%lsp(1:ncni)
      qlir(1:ncni)   =jpiset(igs)%lsp(1:ncni)
      qji (1:ncni)   =jpiset(igs)%jsp(1:ncni)
      qjci(1:ncni)   =jpiset(igs)%jc (1:ncni)
      cindexi(1:ncni)=jpiset(igs)%cindex(1:ncni)
      write(*,*)'The gs has j,parity=',jtoti,partoti,
     & 'and', ncni,' channel(s)'


        aux=0
        do ir=1,nr
        do ich=1,ncni
        ugs(ir,ich)=wfc(1,1,ich,ir)*rvec(ir)
        aux=aux+ugs(ir,ich)**2*dr
        enddo !ich
        enddo !ir
        write(*,*)'GS has norm=',aux
      else !....................................EXTERNAL GS WF
      open(20,file=uwfgsfile)
      read(20,246)ncni,eneri
      allocate(ugs(nr,ncni))
      ugs(:,:)=0d0
      write(*,*)'ncni,eneri=',ncni,eneri
!246	format('# ',i2)!,' Channels, Eigenvalue:',i2,' Enegy:',f8.4)
! AMoro v21b (we need to read energy too)
! 246	format(2x,i2)!,' Channels, Eigenvalue:',i2,' Enegy:',f8.4)
246   format(2x,i2,31x,1f8.4)
      read(20,247)jtoti,partoti
!247	format('# J:',f4.2,' parity:',i2)
247	format(4x,f4.2,8x,i2)
      read(20,*)
      write(*,251)ncni,jtoti,partoti

251   format('- Reading',i2,' chan. for gs: j=',f4.2,' partity=',i2)

      allocate(qjci(ncni),qli(ncni),qlir(ncni),qji(ncni)
     &,cindexi(ncni))
      do i=1,ncni
      read(20,245)incni,cindexi(incni),qjci(incni),qli(incni),
     &sni,qji(incni)
!245	format('# ',i2,' :',i2,f4.1,i2,f4.1,f4.1)
245	format(2x,i2,2x,i2,f4.1,i2,f4.1,f4.1)
!      write(*,*)'Core index=',cindexi(incni)
      enddo
      qlir(1:ncni)=qli(1:ncni)  ! initial (gs)
      if (sni.ne.sn) then
      write(*,*) 'FATAL ERROR: Electric 
     &transition do not allow changes in the spin'
      stop
      endif
      read(20,248)np
!248   format("# ",i5," points")
248   format(2x,i5)!," points")
      allocate(r(np),uext(np,ncni))
      do ir=1,np
      read(20,*)r(ir),(uext(ir,i),i=1,ncni)  !rlast would be r(np)
      enddo
      aux=sum(abs(uext(np,:))**2)
      if (aux.gt.1e-3) then
        write(*,*)' ** WARNING **: external function has norm of',
     &  aux,' at last point r=',r(np),' fm'
      endif
      write(99,*)'- gs function read'
      write(99,*)'- Please, make sure that gs goes to zero at rlast'
      write(99,*)'- First checks ok so:'
      write(*,'("- Calculating B(E",i1,"; gs->n)")')lambda
      alpha=2.
      do ir=1,nr
      if(rvec(ir).gt.r(np)) cycle
      do i=1,ncni
      ugs(ir,i)=fival(rvec(ir),r,uext(:,i),np,alpha)
      enddo
      enddo
      deallocate(uext)
      endif ! type of gs WF
c     -------------------------------------------------------------------------


c     Select FINAL j/pi set and check triangularity condition -----------------.  
      if ((jset.gt.jpsets).or.(jset.eq.0)) then
       write(*,*)'jset',jset,' not valid!'
       return
      endif
      partot=jpiset(jset)%partot
      jtot  =jpiset(jset)%jtot
      nex   =jpiset(jset)%nex
      nchan =jpiset(jset)%nchan
      write(*,*)'Final channel',jset,': j=',jtot, 'par=',partot,
     & ' chans:',nchan,' nex=',nex
      allocate(qlr(nchan))
      qjc(:)=0; 
      ql(:)=0; 
      qlr(:)=0; 
      qj(:)=0; 
      cindex(:)=0

!      allocate(qjc(nchan),ql(nchan),qlr(nchan),qj(nchan),cindex(nchan))
      ql (1:nchan)   =jpiset(jset)%lsp(1:nchan)
      qlr(1:nchan)   =jpiset(jset)%lsp(1:nchan)
      qj (1:nchan)   =jpiset(jset)%jsp(1:nchan)
      qjc(1:nchan)   =jpiset(jset)%jc (1:nchan)
      cindex(1:nchan)=jpiset(jset)%cindex(1:nchan)
      exc(1:nchan)   =jpiset(jset)%exc(1:nchan)
 
      qlr(1:nchan)=ql(1:nchan)  ! final 

      if (allocated(edisc)) deallocate(edisc)
      allocate(edisc(nex))
      edisc(1:nex)  =energ(jset,1:nex)
      allocate(mel(nex,nchan));  mel(:,:)=0d0
      
     
      if (partoti*partot.ne.(-1)**lambda) then
         write(*,'("!! B(E",i2,") does not satisfy parity")')lambda
         stop
      endif
      lambdar=lambda
      if (fail3(jtoti,lambdar,jtot)) then
        write(*,'("!! B(E",i2,") does not couple this Js")')lambda
        stop
      endif
c     --------------------------------------------------------------------------



!c  Here starts the calculation, we have to be sure that only allowed transitions are calculated
      Elam=0d0
      BEl=0d0
      besum=0d0
!     BElcore=0d0
      Elamcore=0d0 ! we have to change this

      select case(coremodel)
      case(0) ! rotor
        write(*,*)'rotor: jc(1),jc(2)=',jc(1),jc(2)
      if (BElcore.ne.0d0) then
        write(*,*)'  [Rotor model with BElcore=',BElcore,']'
!        Elamcorem=sqrt(BElcore*(2.*qjci(n)+1d0)/(2.*qjc(m)+1.))
        Elamcorem=sqrt(BElcore*(2.*jc(1)+1d0)/(2.*jc(2)+1.))
!         write(*,*)jc(1),jc(2)

!        Elamcorem=Elamcorem*delta/(abs(delta))                   !we have to add the other part of the sign yet (from the clebsch gordam)
        write(*,*)'- M(Elambda)_core=',Elamcorem,'  from experimental 
     &            B(Elambda)_core=',BElcore
        write(*,*)'WARNING: Not still full implemented and tested'
      else
        if (rms.ne.0d0) then
        write(*,*)'  [rotor model with delta=',delta,']'
        Elamcorem=3*zc*delta*(sqrt(5/3d0)*rms)**(lambda-1)/(4*pi)
        write(*,'(3x,"[including core excitation with Mn=",
     &   1f8.3," e.fm^lambda from rms and delta]")')  Elamcorem
        endif
        endif
       
       case(1) ! read reduced matrix elements for core transitions
       write(*,*)'  [User defined R.M.E. for core ]'

       case default
        write(*,*)'core model',coremodel,' not used!'
        stop
       end select 

c--------------------------------------------
c       Discrete distribution for B(Elambda)
c--------------------------------------------
        ncont=1
        do i=1,nex
        if (edisc(i).lt.0) then
          ncont=ncont+1 
        endif
        enddo
        write(*,*)'Final j/pi set has',ncont,' bound states'

        zeff=zv*(ac/(ac+av))**lambda+zc*(-av/(ac+av))**lambda
        write(*,*)'- Zeff,ac,av=',zeff,ac,av
        do i=1,nex
        Elam=0d0
        BEl=0d0
	 if (i.lt.ncont) then 
           wid=1. 
           else if (i.eq.ncont) then
           wid= 2.d0*edisc(ncont+1) -1.5d0*edisc(ncont)-
     &          0.5d0*edisc(2+ncont) ! change here!!!!!!!!!!!!1 
	    else if (i.eq.nex) then
           wid=-2.d0*edisc(i-1)+1.5d0*edisc(i)+0.5d0*edisc(i-2)
!dnex=-2.d0*ENEX(NSET-1)+1.5d0*ENEX(NSET)+0.5d0*ENEX(NSET-2)
	    else
           wid=(edisc(i+1)-edisc(i-1))/2d0
        endif

!------------------------------------------------------------------
        mvc=0d0
        itran=0
        do m=1,nchan
        do n=1,ncni
!c CORE contribution
        qjcir=qjci(n)
        qjcr =qjc(m)
        geom=(-1.)**(lambda+qj(m)+qjci(n)+jtot)*
     &       sqrt((2.*jtoti+1.)*(2.*qjc(m)+1.))*
     &       sixj(jtot,jtoti,lambdar,qjcir,qjcr,qj(m))

        if (.not.fail3(qjci(n),lambdar,qjc(m))) then
        if ((qj(m).eq.qji(n)).and.(qli(n).eq.ql(m))) then

!c Overlap between radial parts for n and m channels
!        faux(:)=rvec(:)*ugs(:,n)*wfeig(i,m,:) !the gs has already an r
        faux(1:nr)=rvec(1:nr)*ugs(1:nr,n)*wfc(jset,i,m,1:nr) !the gs has already an r
        call simc(faux,solap,1,nr,dr,nr)    !it will improve if we go only until rlast

!c Core reduced matrix element (model dependent) -------------------------------
        select case(coremodel)
        case(0) ! ----------->  rotor
        if (rms.eq.0.) then
        rme=Elamcorem*cleb(qjcir,0d0,lambdar,0d0,qjcr,0d0)*    ! this is only for K=0, but this is the same. (old version)
     &sqrt((2*qjcir+1)/(2*qjcr+1))*(-1)**(2*lambda)
        else
        rme=rotor(qjcir,qjcr,0d0,lambdar,rms,delta,zc)    ! this is only for K=0, but this is the same.
        endif
        
        case(1) !----------->  external M(E) 
        rme=mec(cindexi(n),cindex(m),lambda)
!         write(*,*)'rme:',cindexi(n),cindex(m),rme
        end select 

!        faux(:)=rvec(:)**(1)*ugs(:,n)*wfeig(iord(i),m,:) !the gs have already an r
!        call sim(faux,res,1,nr,dr,nr)    !it will improve if we go only until rlast
        Elamcore=geom*solap*rme
        Elam=Elam+Elamcore
        mvc=mvc+Elamcore
        itran=itran+1
        if (verb.ge.3) 
     &  write(99,'(4x,i3," Core(PS):",1f4.1,1i2,1f4.1,"->",
     &   1f4.1,1i2,1f4.1,5f9.3)') 
     & itran,qjcir,qli(n),qji(n),qjcr,ql(m),qj(m),
     & geom,solap,rme,Elamcore
        mel(i,m)=mel(i,m)+cmplx(Elamcore,0d0)
        Elamcore=0d0
        endif
        endif
c     --------------------------------------------------------------------------

c VALENCE contribution
      if ((.not.fail3(qji(n),lambdar,qj(m)))        !lamda is enough to connect the two states
     &.and.(cindexi(n).eq.cindex(m))) then          !we avoid different channels with same Jpi
         if (.not.fail3(qlir(n),lambdar,qlr(m))) then                
         
c < n | r^lambda | m > 
!        faux(:)=rvec(:)**(1+lambda)*ugs(:,n)*wfeig(i,m,:) !the gs have already an r
       faux(1:nr)=rvec(1:nr)**(1+lambda)
     &            *ugs(1:nr,n)*wfc(jset,i,m,1:nr) !the gs have already an r
       call simc(faux,rlam,1,nr,dr,nr)    !it will improve if we go only until rlast

       mv=zeff*rlam*
     & matel(lambdar,sn,qjc(m),jtoti,qji(n),qlir(n),jtot,qj(m),qlr(m))
       mvc=mvc+mv

! 2.0.4
!        mel(iord(i),m)=mel(iord(i),m)+cmplx(zeff*matel(
!     &lambdar,sn,qjc(m),jtoti,qji(n),qlir(n),jtot,qj(m),qlr(m))*
!     &res,0d0)+cmplx(Elamcore,0d0)

!        mel(iord(i),m)=mel(iord(i),m)+mvc   !it seems we have to add better mv or wait untill finish m cycle
!        write(*,*)'1) i,m=',i,m
        mel(i,m)=mel(i,m)+mv


         if (verb.ge.4) 
     &     write(99,'(4x,"Valence:",2f4.1,"->",2f4.1,1f9.4)') 
     &     qlir(n),qji(n),qlr(m),qj(m),mv
         endif
         endif
         enddo !n
        enddo  ! m

! AMoro 2.0.5
!        BEl=abs(Elam)**2*(2*jtot+1)/(2*jtoti+1)
        BEl=mvc**2*(2*jtot+1)/(2*jtoti+1)
        Elam=0d0
c        written(95)=.true.
        write(95,490) edisc(i),BEl/wid,BEl,wid
490     format(1x,f8.4,3f12.6)
        besum=besum+BEl
        write(*,300)edisc(i),lambda,i,BEl,2*lambda
300     format(4x,"Ex=",f8.4,"  B(E",i1,"; gs ->",i4,
     &")=",f8.4," e2.fm",i1)
!        test=0d0
        enddo !i=1,nex

	written(95)=.true.

        write(*,*)'- From Pseudo-States:'
        write(*,'(30x,"Total BE=",1f10.6," e2.fm",i1)')besum,2*lambda
        write(*,*)'- From Sum Rule: (only for SP excitation!)'
        call sumrule
        write(*,'(30x,"Total BE=",1f10.6," e2.fm",i1)')besr,2*lambda 
!        write(*,*)"(** sum rule only intended for sp excitations)"


c------------------------------------------------
c B(Elambda) from scattering states (version by AMM)
c------------------------------------------------ 
      written(96)=.true.
      rm=av*ac/(av+ac)
      conv=(2*amu/hc**2)*rm
      kmin=sqrt(2d0*mu12*abs(emin))/hc
      kmax=sqrt(2d0*mu12*abs(emax))/hc
      ili=1
      il =nchan

!!!!!!! NEED TO BE REVISED FOR OPEN CHANNELS WITH EXCITED CORE
!      write(*,*)'cont wf for jset,ili,il=',jset,ili,il
!      write(*,*)'emin,emax=',emin,emax
!      do iil=ili,il !open channels
!      call wfrange(jset,nchan,iil,emin,emax,nk,energy,
!     &             wfcont(:,iil,:,:),smate,psh)
!      wfcont(:,iil,:,:)=sqrt(2./pi)*wfcont(:,iil,:,:) ! set norm to delta(k-k')
!      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (nk.gt.1) then
         dk=(kmax-kmin)/(nk-1)
      else
         dk=(kmax-kmin)
      endif
      besum=0d0
      if (.not.allocated(dbde)) allocate(dbde(il,nk))
      dbde(:,:)=0.
      do iil=ili,il ! incoming channels
        if (allocated(wfscat)) deallocate(wfscat)
         allocate(wfscat(nk,nchan,nr))

         write(96,*)'# inc. channel=',iil
         excore =jpiset(jset)%exc(iil) ! core exc. energy for this inc. chan.
         call wfrange(jset,nchan,iil,emin,emax,nk,energy,
     &             wfscat(:,:,:),smate,psh)
         wfscat(:,:,:)=sqrt(2./pi)*wfscat(:,:,:) ! set norm to delta(k-k')

      do ik=1,nk
         Elamcont=zero
         Elamcore=0d0
         mecc=(0d0,0d0)
         BEl=0d0
         kcont=kmin+ (ik-1)*dk  !new kcont for schcc
         econt=(hc*kcont)**2/2/mu12

!         do iil=ili,il !open channels 
!           excore =jpiset(jset)%exc(iil) 
!            krel=sqrt(2d0*mu12*(econt-excore))/hc            
!	    if (krel.gt.eps) then
!		dkrdk=kcont/krel
!	    else
!		dkrdk=1d0
!            endif
!        Elamcont=zero
        itran=0
        do m=1,nchan 
        do n=1,ncni  ! init  channels
!c CORE contribution
        qjcir=qjci(n)
        qjcr =qjc(m)
        geom=(-1.)**(lambda+qj(m)+qjci(n)+jtot)
     &      *sqrt((2.*jtoti+1.)*(2.*qjc(m)+1.))
     &      *sixj(jtot,jtoti,lambdar,qjcir,qjcr,qj(m))

        if (.not.fail3(qjci(n),lambdar,qjc(m))) then
        if ((qj(m).eq.qji(n)).and.(qli(n).eq.ql(m))) then

!c Overlap between radial parts for n and m channels
!        gaux(:)=ugs(:,n)*wfcont(ik,iil,m,:)*sqrt(dkrdk) !only one incoming channel
        gaux(:)=ugs(:,n)*wfscat(ik,m,:) !only one incoming channel
        call simc(gaux(:),resc,1,nr,dr,nr)

!c Core reduced matrix element (model dependent)        
        select case(coremodel)
        case(0) ! ----------->  rotor
        if (rms.eq.0d0) then
        rme=Elamcorem*cleb(qjcir,0d0,lambdar,0d0,qjcr,0d0)*    ! this is only for K=0, but this is the same. (old version)
     &      sqrt((2*qjcir+1)/(2*qjcr+1))*(-1)**(2*lambda)
        else
            rme=rotor(qjcir,qjcr,0d0,lambdar,rms,delta,zc)    ! this is only for K=0, but this is the same.
        endif
        
        case(1) !----------->  external M(E) 
        rme=mec(cindexi(n),cindex(m),lambda)
        end select 

        Elamcore=geom*resc*rme*sqrt(pi/2)*
     &           (4*pi)/kcont
        mecc=geom*resc*rme*sqrt(pi/2)*
     &           (4*pi)/kcont
!commented by AMM (2.1)
!        Elamcont=Elamcont  + Elamcore
        Elamcont=Elamcont + mecc
        itran=itran+1
        if(verb.ge.3) then
        write(99,'(4x,i3,"Core(C):",1f4.1,1i2,1f4.1,"->",
     & 1f4.1,1i2,1f4.1,6f9.3)') 
     & itran,qjcir,qli(n),qji(n),qjcr,ql(m),qj(m),
     & geom,abs(resc),rme,Elamcore
        endif
        Elamcore=0d0

        endif
        endif

c VALENCE contribution
      if ((.not.fail3(qji(n),lambdar,qj(m)))        !lamda is enough to connect the two states
     &.and.(cindexi(n).eq.cindex(m))) then          !we avoid different channels with same Jpi
         if (.not.fail3(qlir(n),lambdar,qlr(m))) then                

       coef=zeff*matel(lambdar,sn,qjc(m),jtoti,qji(n),
     &                 qlir(n),jtot,qj(m),qlr(m))
         
c < n | r^lambda | m > 
       gaux(1:nr)=rvec(1:nr)**(lambda)
     &           *ugs(1:nr,n)*wfscat(ik,m,1:nr)
!     &           *ugs(1:nr,n)*wfcont(ik,iil,m,1:nr)*sqrt(dkrdk)
       call simc(gaux,resc,1,nr,dr,nr)   

       Elamcont = Elamcont  
     &          + coef*resc*sqrt(pi/2)*(4*pi)/kcont
         endif
         endif
        enddo !n
        enddo  !m

	 BEl=BEl+ (2*jtot+1)*abs(Elamcont)**2/(2*jtoti+1)
     &      *mu12*kcont/hc**2/(2*pi)**3
        dbde(iil,ik)=dbde(iil,ik)
     &      + (2*jtot+1)*abs(Elamcont)**2/(2*jtoti+1)
     &      *mu12*kcont/hc**2/(2*pi)**3
        Elam=0d0

        write(96,*) econt+excore,BEl,dbde(1,ik)
        besum=besum+BEl*hc**2*kcont/mu12*dk    !de=de/dk*dk
        enddo ! ik
        write(*,*)'inc,besum=',iil,besum
        write(96,*)'&'
        enddo ! iil
        write(*,*)'- From Continuum wfs:'
        write(*,'(30x,"Total BE=",1f10.6," e^2 fm^",i1)')besum,2*lambda

!c------------------------------------------------
!c        B(Elambda) from scattering states
!c------------------------------------------------        
!         if (.not.ifcont) return    !if we don't calculate continuum wf, we end here.
!        besum=0d0
!        do ik=1,nk
!         Elamcont=zero
!          kcont=kmin+(kmax-kmin)*(ik-1)/(nk-1)  !new kcont for schcc
!          econt=(hc*kcont)**2/2/mu12
!        do m=1,nchan
!          do n=1,ncni
!        qjcir=1d0*qjci(n)
!        qjcr=1d0*qjc(m)
!           do iil=ili,il
!	    krel=sqrt(2d0*mu12*(econt-exc(iil)))/hc
!		if (krel.gt.eps) then
!		  dkrdk=kcont/krel
!		  else
!		  dkrdk=1d0
!		endif
!c core contribution
!	   if (qjci(n).ne.qjc(m)) then
!             if (.not.fail3(qjci(n),lambdar,qjc(m))) then
!             if ((qj(m).eq.qji(n)).and.(qli(n).eq.ql(m))) then
!        if (rms.eq.0d0) then
!        Elamcore=Elamcorem*cleb(qjcir,0d0,lambdar,0d0,qjcr,0d0)*                    ! this is only for K=0, but this is the same. (old version)
!     &sqrt((2*qjcir+1)/(2*qjcr+1))*(-1)**(2*lambda)
!        else
!        Elamcore=rotor(qjcir,qjcr,0d0,lambdar,rms,delta,zc)                         ! this is only for K=0, but this is the same.
!        endif
!        Elamcore=Elamcore*(-1.)**(lambda+qj(m)+qjci(n)+jtot)*
!     &sqrt((2.*jtoti+1.)*(2.*qjc(m)+1.))*
!     &sixj(jtot,jtoti,lambdar,qjcir,qjcr,qj(m))
!              gaux(:)=ugs(:,n)*wfcont(ik,iil,m,:)*sqrt(dkrdk) !only one incoming channel
!              faux(:)=dreal(gaux(:))
!              coef=0d0
!              res=0d0
!              call sim(faux,res,1,nr,dr,nr)  
!               faux(:)=dimag(gaux(:))
!              res2=0d0
!              call sim(faux,res2,1,nr,dr,nr)    !it will improve if we go only until rlast
!               Elamcont=Elamcont+Elamcore*cmplx(res,res2)*sqrt(pi/2)*
!     &                  (4*pi)/kcont
!               Elamcore=0d0
!             endif
!             endif
!           else
!      if ((.not.fail3(qji(n),lambdar,qj(m)))                                !lamda is enough to connect the two states
!     &.and.(cindexi(n).eq.cindex(m))) then                                  !we avoid different channels with same Jpi
!              if (.not.fail3(qlir(n),lambdar,qlr(m))) then                  !I already checked that without these fail3s, the matel gives zero contribution in the not allowed cases
!             Elamcore=0d0
!             if ((qj(m).eq.qji(n)).and.(qli(n).eq.ql(m))) then
!        if (rms.eq.0d0) then
!        Elamcore=Elamcorem*cleb(qjcir,0d0,lambdar,0d0,qjcr,0d0)*                    ! this is only for K=0, but this is the same. (old version)
!     &sqrt((2*qjcir+1)/(2*qjcr+1))*(-1)**(2*lambda)
!        else
!        Elamcore=rotor(qjcir,qjcr,0d0,lambdar,rms,delta,zc)                         ! this is only for K=0, but this is the same.
!        endif
!        Elamcore=Elamcore*(-1.)**(lambda+qj(m)+qjci(n)+jtot)*
!     &sqrt((2.*jtoti+1.)*(2.*qjc(m)+1.))*
!     &sixj(jtot,jtoti,lambdar,qjcir,qjcr,qj(m))
!              gaux(:)=ugs(:,n)*wfcont(ik,iil,m,:)*sqrt(dkrdk) !the gs have already an r
!              faux(:)=dreal(gaux(:))
!              res=0d0
!              call sim(faux,res,1,nr,dr,nr)    !it will improve if we go only until rlast
!        Elamcont=Elamcont+cmplx(Elamcore*res,0d0)*sqrt(pi/2)*
!     &(4*pi)/kcont
!              faux(:)=dimag(gaux(:))
!              res=0d0
!              call sim(faux,res,1,nr,dr,nr)    !it will improve if we go only until rlast
!        Elamcont=Elamcont+cmplx(0d0,Elamcore*res)*sqrt(pi/2)*
!     &(4*pi)/kcont
!             else
!             Elamcore=0d0
!             endif
!c Valence
!              gaux(:)=rvec(:)**(lambda)*ugs(:,n)*wfcont(ik,iil,m,:)*
!     &sqrt(dkrdk) !only one incoming channel
!              faux(:)=dreal(gaux(:))
!              coef=0d0
!              coef=zeff*matel(lambdar,sn,qjc(m),jtoti,qji(n),
!     &qlir(n),jtot,qj(m),qlr(m))
!              res=0d0
!              call sim(faux,res,1,nr,dr,nr)  
!               faux(:)=dimag(gaux(:))
!              res2=0d0
!              call sim(faux,res2,1,nr,dr,nr)    !it will improve if we go only until rlast
!          Elamcont=Elamcont+coef*cmplx(res,res2)*sqrt(pi/2)*(4*pi)/kcont
!              endif
!           endif
!           endif
!           enddo
!          enddo
!          enddo
!	      BEl=(2*jtot+1)*abs(Elamcont)**2/(2*jtoti+1)
!          BEl=BEl*mu12*kcont/hc**2/(2*pi)**3
!          besum=besum+BEl*hc**2*kcont/mu12*(kmax-kmin)/(nk-1) !de/dk*dk
!          write(96,*)econt,BEl
!          enddo
!         written(96)=.true.
!         write(*,*)'- From Continuum wfs:'
!		 write(*,'(30x,"Total BE=",1f10.6," e^2 fm^",i1)')besum,2*lambda



c-----------------------------------------------------------------------------
c   B(Elambda) from discrete distribution convoluted with continuous wfs (AMM style)
c-----------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*)'NEED TO UPDATE SOLAPS HERE; RETURNING!!'; RETURN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      besum=0d0
      BEl=0d0

      do ik=1,nk
      BEl=0d0
      test=zero
      Elamcont=zero
      kcont=kmin+(kmax-kmin)*(ik-1)/(nk-1)  !new kcont for schcc          
      econt=(hc*kcont)**2/2/mu12
      do iil=ili,il	  
	krel=sqrt(2d0*mu12*(econt-exc(iil)))/hc
	if (krel.gt.eps) then
	  dkrdk=kcont/krel
	else
	  dkrdk=1d0
	endif
	  
       Elamcont=zero
       do i=1,nex
       do n=1,nchan
       Elamcont=Elamcont+solapmat(ik,iil,i,n)*
     &                 mel(i,n)*sqrt(dkrdk)
       enddo
       enddo
       Elamcont=Elamcont*sqrt(pi/2)*(4*pi)/kcont
      BEl=BEl+(2*jtot+1)*abs(Elamcont)**2/(2*jtoti+1)
       enddo
!          write(98,*)econt,BEl
      BEl=BEl*mu12*kcont/hc**2/(2*pi)**3
      besum=besum+BEl*hc**2*kcont/mu12*(kmax-kmin)/(nk-1) 
      write(97,*)econt,BEl
      enddo

      written(97)=.true.
      write(*,*)'- From Folding:'
      write(*,'(30x,"Total BE=",1f10.6," e^2 fm^",i1)')besum,2*lambda


!c-----------------------------------------------------------------------------
!c   B(Elambda) from discrete distribution convoluted with continuous wfs
!c-----------------------------------------------------------------------------
!         besum=0d0
!         BEl=0d0
!         do ik=1,nk
!         BEl=0d0
!         test=zero
!         Elamcont=zero
!          kcont=kmin+(kmax-kmin)*(ik-1)/(nk-1) 
!		  econt=(hc*kcont)**2/2/mu12
!         do i=1,nex
!            do n=1,nchan
!              do iil=ili,il
!		krel=sqrt(2d0*mu12*(econt-exc(iil)))/hc
!		if (krel.gt.eps) then
!		  dkrdk=kcont/krel
!		  else
!		  dkrdk=1d0
!		endif
!              Elamcont=Elamcont+solapmat(ik,iil,iord(i),n)*
!     &                 mel(iord(i),n)*dkrdk
!              enddo
!            enddo
!         enddo
!          Elamcont=Elamcont*sqrt(pi/2)*(4*pi)/kcont
!	  BEl=(2*jtot+1)*abs(Elamcont)**2/(2*jtoti+1)
!          BEl=BEl*mu12*kcont/hc**2/(2*pi)**3
!          besum=besum+BEl*hc**2*kcont/mu12*(kmax-kmin)/(nk-1) 
!          write(97,*)econt,BEl
!         enddo
!         written(97)=.true.
!         write(*,*)'- From Folding:'
!         write(*,'(30x,"Total BE=",1f10.6," e^2 fm^",i1)')besum,2*lambda
500      end 


        subroutine sumrule
        use constants, only: pi
        use sistema
        use wfs, only: dr,nr,rvec
        use channels
        use belambdamod
        implicit none
        integer:: n,m,c,q
        real*8:: cr,matel,threej,lambdar
        real*8:: cleb,qr,sum3j
        complex*16::faux(nr),res
        logical:: fail3



        lambdar=1d0*lambda
        besr=0d0
	do n=1,ncni
          do m=1,ncni
           if ((qjci(n).eq.qjci(m)).and.(cindexi(n).eq.cindexi(m))) then             !qjci is needed if cindexi takes into account only repeated Jpi
!           if (cindexi(n).eq.cindexi(m)) then                                    
           do c=0,2*lambda
           cr=1d0*c
           if (.not.fail3(qji(n),cr,qji(m))) then
              if (.not.fail3(qlir(n),cr,qlir(m))) then            
              faux(1:nr)=rvec(1:nr)**(2*lambda)*ugs(1:nr,n)*ugs(1:nr,m) 
              call simc(faux,res,1,nr,dr,nr)
                 
!              besr=besr+matel(cr,sn,qjci(m),jtoti,qji(n),
!     &qlir(n),jtoti,qji(m),qlir(m))*res*sqrt(2*cr+1)*
!     &threej(lambdar,lambdar,cr,0d0,0d0,0d0)**2
!     &*cleb(jtoti,jtoti,cr,0d0,jtoti,jtoti)

          sum3j=0d0
          do q=-c,c
            qr=1d0*q
            sum3j=sum3j+
     &threej(lambdar,lambdar,cr,-qr,qr,0d0)*(-1)**q
          enddo
              besr=besr+matel(cr,sn,qjci(m),jtoti,qji(n),          !calculating sum over q
     &qlir(n),jtoti,qji(m),qlir(m))*res*sqrt(2*cr+1)*
     &threej(lambdar,lambdar,cr,0d0,0d0,0d0)*(-1)**(2*c)*
     &sum3j*cleb(jtoti,jtoti,cr,0d0,jtoti,jtoti)

              endif
           endif
           enddo
           endif
          enddo
        enddo
        besr=besr*(2*lambda+1)/sqrt(4*pi)*(2*jtot+1)/(2*jtoti+1)
        besr=besr*zeff**2

        end


        function matel(lambda,sn,ic,jtoti,jni,li,jtotf,jnf,lf)
        use constants, only: pi
        implicit none
        real*8:: lambda,sn,ic,jtoti,jni,li,jtotf,jnf,lf   !they have to be all real at the end to fit the 3j,6j subroutines.
        real*8:: matel,sixj,threej
        matel=(-1)**(sn+ic+jtoti+jni+jnf)*sqrt((2*lambda+1)/(4*pi))

        matel=matel*sqrt((2*jtoti+1)*(2*jni+1)*(2*jnf+1)*
     &(2*lf+1)*(2*li+1))
        matel=matel*sixj(jtotf,jtoti,lambda,jni,jnf,ic)

        matel=matel*sixj(jnf,jni,lambda,li,lf,sn)

        matel=matel*threej(li,lf,lambda,0d0,0d0,0d0)

!        matel=(-1)**(jni-jnf)*matel !in case of coupling (I jn) J instead of (jn I) J; our coupling is ((l s) j I ) J
        end

