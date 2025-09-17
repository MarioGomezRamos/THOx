
      subroutine bmlam() !,jtot,partot
c ** Calculates B(M;lambda) distributions using THO pseudostates
      use globals !, only: written,mu12,pi,verb
      use constants, only: hc,pi,amu,e2
      use hmatrix, only: hmatx
      use potentials, only: maxlamb!,delta
      use sistema
      use wfs !, only: jc,parc,rvec,ql,qjc,qj,iord,cindex !,exc
      use channels
      use bmlambdamod
      use scattering
      use factorials
      implicit none
c      ------------------------------------------------------------------------
      logical :: energy,fail3,writewf!,ifbel
      integer :: i,np,ir,m,n,ncont,iil,is,itran,ik,jset,nex,ich
      integer:: pci,pcf,ibel!MGR 2022
      integer, parameter:: igs=1
      integer:: jseti,ni
c     --------------------------------------------------------------------------
      real*8:: fival,lambdar
      real*8:: mlam,mlamcore,Bml,Bmlcore,matml_l,bmsum,rms,delta
      real*8:: matml_sv,matml_sc
      real*8:: mlamcorem!,rotor
      real*8:: wid,qjcir,qjcr,cleb,sixj,coef
      real*8:: mmc(maxcore,maxcore,maxlamb) 
      real*8:: econt,kcont,test,rmm,geom,excore
      real*8:: jci,jcf,bmc,alpha
      real*8:: mvc,mv, rm,conv,solap,rlam,dk
      real*8:: dkrdk,krel   !depends on wfcont and solap definitions
      real*8:: eps,aux
      real*8, allocatable:: edisc(:) ! discrete energies (PS's or Bins)
      real*8, allocatable:: uext(:,:)
      real*8, allocatable:: qlr(:)
      real*8, allocatable:: psh(:)
      real*8:: capxs,photoxs,kphot,engs,sfactcap,sfactphoto,eta !MGR 2022
      real*8 massp
      real*8 gyrov,gyroc,mun2
c     --------------------------------------------------------------------------
      complex*16:: mlamcont,zero,resc,mmcc
      complex*16:: gaux(nr),faux(nr),solapc!,test
      complex*16,allocatable:: smate(:,:)
      complex*16,allocatable:: wfscat(:,:,:)
!      complex*16,allocatable :: ugs(:,:)
      integer coremodel,qc,ic,icp,nctrans
      CHARACTER*1 PSIGN(3)
      character*5 jpi

      
      DATA PSIGN / '-','?','+' / 
c     --------------------------------------------------------------------------  
      namelist /bmlambda/ uwfgsfile,lambda,Bmlcore,ifbml,rms,
     &  coremodel,jset,emin,emax,nk,engs,jseti,ni,delta,gyrov,gyroc
      namelist /coretrans/ic, icp, rmm, qc  
c     --------------------------------------------------------------------------
      pi=acos(-1.0)
      writewf=.false.
      ifcont=.false.
      call factorialgen(10)
      do ibel=1,100
      eps=1e-6
      zero=(0d0,0d0)
      Bmlcore=0d0
      mlamcorem=0d0
      rms=0d0
      delta=0d0
      engs=0d0
      coremodel=0
      nctrans=0
      uwfgsfile=""
      energy=.false. 
      jset=0
      lambda=-1 !MGR 2022
      jseti=1 ! default initial set
      ni=1    ! default initial state within initial set
      massp=1.0072d0
      mun2=0.01592234957868851303d0
      gyrov=0d0
      gyroc=0d0
      
      read(kin,nml=bmlambda)
      
      if (.not.ifbml) return
      
      if (lambda .eq. -1) return
!      write(0,*) 'nk',nk,'maxchan',maxchan
      allocate(smate(nk,maxchan),psh(nk))
      smate(:,:)=0. ; psh(:)=0
!      


!Extract gyromagnetic factor from Stone table (https://www.nndc.bnl.gov/nndc/stone_moments/nuclear-moments.pdf)
      !Valence
      if (abs(sn).gt.1e-3 .and. abs(gyrov).lt.1e-5) then
      write(*,*) 'Getting g-factor for valence from Stone table'
      write(*,*)'https://www.nndc.bnl.gov/nndc/stone_moments
     & /nuclear-moments.pdf'
      call stone_gyro(zv,av,gyrov)
      write(*,*) 'g-factor for valence set to:',gyrov
      endif
      
      !Core
      if (abs(qnc(1)%jc).gt.1e-3 .and. abs(gyroc).lt.1e-5) then
      write(*,*) 'Getting g-factor for core(g.s.) from Stone table'
      write(*,*)'https://www.nndc.bnl.gov/nndc/stone_moments
     & /nuclear-moments.pdf'
      call stone_gyro(zc,ac,gyroc)
      write(*,*) 'g-factor for core(g.s.) set to:',gyroc
      endif
      
      

c     ----------------- Core RME specified by user -----------------------------
      if (coremodel.eq.1) then
      mmc=0d0
      write(*,*)'o Core reduced matrix elements:'
230   rmm=0d0;
      read(kin,nml=coretrans)
      if (rmm.eq.0) goto 235
      nctrans=nctrans+1
      jci=jc(ic)
      jcf=jc(icp)
      pci=parc(ic)
      pcf=parc(icp)
!      bec=rme**2/(2d0*jci+1)
!c CHECK THIS non-symmetric M(I->I')
      mmc(ic,icp,qc)=rmm !/sqrt(2*jcf+1)
      bmc=rmm**2*(2*jcf+1)/(2d0*jci+1)
      write(*,232)qc,jci,psign(pci+2),jcf,psign(pcf+2),bmc
232   format(3x,"B(M",i1,";",1f3.1,a1," -> ",1f3.1,a1,")=",1f7.2)
      goto 230
      endif ! coremodel
c     --------------------------------------------------------------------------
        
235   write(*,*)    
      write(*,'(5x,//," ** CALCULATION OF B(MLAM) **",/ )')
      
   
      
      if ((jset.gt.jpsets).or.(jset.eq.0)) then
       write(*,'(5x,a,i3,a,i3,a)')'=> B(M) requested for set jset',jset,
     & ' but only', jpsets,' available!. Skipping'
      return
      endif
        
c     ------------------------------ GS WF  -----------------------------------
      write(*,'("- Calculating B(M",i1,"; gs->n)")')lambda
      write(95,'("#Calculating B(M",i1,"; gs->set ",i2,")")')lambda,jset
      write(96,'("#Calculating B(M",i1,"; gs->set ",i2,")")')lambda,jset

!      zeff=zv*(ac/(ac+av))**lambda+zc*(-av/(ac+av))**lambda
      zeff=massp*ac*av/(ac+av)*(zv/av**2*(ac/(av+ac))**(lambda-1)+
     & zc/ac**2*(-av/(av+ac))**(lambda-1))
      write(*,'(3x,"Effective charge=",1f9.4,"mu_N"/)')zeff

      if (uwfgsfile.eq."") then
      ncni=jpiset(jseti)%nchan
   
      if(.not. allocated(qjci))
     & allocate(qjci(ncni),qli(ncni),qlir(ncni),qji(ncni),cindexi(ncni))
      if(.not. allocated(ugs))
     & allocate(ugs(nr,ncni))
      ugs(:,:)=0d0
      partoti        =jpiset(jseti)%partot
      jtoti          =jpiset(jseti)%jtot
      qli (1:ncni)   =jpiset(jseti)%lsp(1:ncni)
      qlir(1:ncni)   =jpiset(jseti)%lsp(1:ncni)
      qji (1:ncni)   =jpiset(jseti)%jsp(1:ncni)
      qjci(1:ncni)   =jpiset(jseti)%jc (1:ncni)
      cindexi(1:ncni)=jpiset(jseti)%cindex(1:ncni)
      
      if (abs(engs).lt.1e-6)  engs=energ(jseti,ni)
      
      write(*,'(2x,"Initial state:",/, 5x,"Jset=",i3,3x, "J/pi=",a5,  
     & 3x,"Nb. chans:",i2,3x,"Energy:",1f8.3, "MeV")')
     &  jseti,jpi(jtoti,partoti),ncni,engs
!
!      write(*,*)'The initial state has j,parity,energy=',
!     & jtoti,partoti,engs,' and', ncni,' channel(s)'
        aux=0
        do ir=1,nr
        do ich=1,ncni
        ugs(ir,ich)=wfc(jseti,ni,ich,ir)*rvec(ir)
        aux=aux+ugs(ir,ich)**2*dr
        enddo !ich
        enddo !ir
        write(*,'(5x,"[ Norm=",1f10.5,"]")'),aux
      else !....................................EXTERNAL GS WF (ASSUMED REAL!!!!!)
      write(*,'(2x,"Initial WF from file: ",a)') uwfgsfile
      open(20,file=uwfgsfile)
      read(20,246)ncni,eneri
      allocate(ugs(nr,ncni))
      ugs(:,:)=0d0
!      write(*,*)'ncni,eneri=',ncni,eneri
!246	format('# ',i2)!,' Channels, Eigenvalue:',i2,' Enegy:',f8.4)
! AMoro v21b (we need to read energy too)
! 246	format(2x,i2)!,' Channels, Eigenvalue:',i2,' Enegy:',f8.4)
246   format(2x,i2,32x,1f8.4)
      read(20,247)jtoti,partoti
!247	format('# J:',f4.2,' parity:',i2)
247	format(4x,f4.2,8x,i2)
      read(20,*)
      write(*,251)ncni,jpi(jtoti,partoti)
251   format(3x,'-> reading',i2,' chans. for gs: j/pi=',a5)

      allocate(qjci(ncni),qli(ncni),qlir(ncni),qji(ncni)
     &,cindexi(ncni))
      do i=1,ncni
      read(20,245)incni,cindexi(incni),qjci(incni),qli(incni),
     &sni,qji(incni)
!245	format('# ',i2,' :',i2,f4.1,i2,f4.1,f4.1)
245	format(2x,i2,2x,i2,f4.1,i2,f4.1,f4.1)
!      write(*,*)'incni,Core index=',incni,cindexi(incni)
      enddo
      qlir(1:ncni)=qli(1:ncni)  ! initial (gs)
!      qlr(1:nchan)=jpiset(jset)%lsp(1:nchan)   ! final
!      qlr(1:nchan)=ql(1:nchan)  ! final 
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
      alpha=2.
      do ir=1,nr
      if(rvec(ir).gt.r(np)) cycle
      do i=1,ncni
      aux=fival(rvec(ir),r,uext(:,i),np,alpha)
      ugs(ir,i)=cmplx(aux,0.0)
      enddo
      enddo
      deallocate(uext)
      endif ! type of gs WF
c     -------------------------------------------------------------------------


c     Select FINAL j/pi set and check triangularity condition -----------------.  

      partot=jpiset(jset)%partot
      jtot  =jpiset(jset)%jtot
      nex   =jpiset(jset)%nex
      nchan =jpiset(jset)%nchan
      write(*,'(2x,"Final states:",/, 5x,"Jset=",i3,3x, "J/pi=",a5,  
     & 3x,"Nb. chans:",i2,3x,"Nb. states:",i3)')
     &  jset,jpi(jtot,partot),nchan,nex
!     & ' chans:',nchan,' nex=',nex
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
      if (allocated(mml)) deallocate(mml)
      allocate(mml(nex,nchan));  mml(:,:)=0d0
      
     
      if (-partoti*partot.ne.(-1)**lambda) then
         write(*,'("!! B(M",i2,") does not satisfy parity")')lambda
         stop
      endif
      lambdar=lambda
      if (fail3(jtoti,lambdar,jtot)) then
        write(*,'("!! B(M",i2,") does not couple this Js")')lambda
        stop
      endif
c     --------------------------------------------------------------------------



!c  Here starts the calculation, we have to be sure that only allowed transitions are calculated
      mlam=0d0
      Bml=0d0
      bmsum=0d0
!     BElcore=0d0
      mlamcore=0d0 ! we have to change this

      select case(coremodel)
      case(0) ! rotor
      if (Bmlcore.ne.0d0) then
        write(*,*)'  [Rotor model with BElcore=',Bmlcore,']'
!        Elamcorem=sqrt(BElcore*(2.*qjci(n)+1d0)/(2.*qjc(m)+1.))
         mlamcorem=sqrt(abs(Bmlcore)*(2.*jc(1)+1d0))!/(2.*jc(2)+1.))
!         write(*,*)jc(1),jc(2)
         mlamcorem=mlamcorem*Bmlcore/abs(Bmlcore) 
!        Elamcorem=Elamcorem*delta/(abs(delta))                   !we have to add the other part of the sign yet (from the clebsch gordam)
        write(*,'(3x,"[including core excitation with Mn=",
     &   1f8.3," e.fm^lambda from experimental]")')  mlamcorem
!      else
!        if (rms.ne.0d0) then
!        write(*,*)'  [rotor model with delta=',delta,']'
!        Elamcorem=3*zc*delta*(sqrt(5/3d0)*rms)**(lambda-1)/(4*pi)
!        write(*,'(3x,"[including core excitation with Mn=",
!     &   1f8.3," e.fm^lambda from rms and delta]")')  Elamcorem
!        endif
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
        write(*,'(5x,"[Final set has",i2," bound states ]",/)')ncont-1


        do i=1,nex
        mlam=0d0
        Bml=0d0
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

        call simc(faux,solapc,1,nr,dr,nr)    !it will improve if we go only until rlast
        solap=dble(solapc)
!c Core reduced matrix element (model dependent) -------------------------------
        select case(coremodel)
        case(0) ! ----------->  rotor
!        rme=Elamcorem*cleb(qjcir,0d0,lambdar,0d0,qjcr,0d0)*    ! this is only for K=0, but this is the same. (old version)
!     &sqrt((2*qjcir+1)/(2*qjcr+1))*(-1)**(2*lambda)
        
        rmm=0
        case(1) !----------->  external M(E) 
        rmm=mmc(cindexi(n),cindex(m),lambda)
!         write(*,*)'rme:',cindexi(n),cindex(m),rme
        end select 

!        faux(:)=rvec(:)**(1)*ugs(:,n)*wfeig(iord(i),m,:) !the gs have already an r
!        call sim(faux,res,1,nr,dr,nr)    !it will improve if we go only until rlast
        mlamcore=geom*solap*rmm
        mlam=mlam+mlamcore
        mvc=mvc+mlamcore
        itran=itran+1
        if (verb.ge.3) 
     &  write(99,'(4x,i3," Core(PS):",1f4.1,1i2,1f4.1,"->",
     &   1f4.1,1i2,1f4.1,5f9.3)') 
     & itran,qjcir,qli(n),qji(n),qjcr,ql(m),qj(m),
     & geom,solap,rmm,mlamcore
        mml(i,m)=mml(i,m)+cmplx(mlamcore,0d0)
        mlamcore=0d0
        endif
        endif
c     --------------------------------------------------------------------------

c VALENCE contribution
      if ((.not.fail3(qji(n),lambdar,qj(m)))        !lamda is enough to connect the two states
     &.and.(cindexi(n).eq.cindex(m))) then          !we avoid different channels with same Jpi
         if (.not.fail3(qlir(n),lambdar-1,qlr(m))) then                
         
c < n | r^lambda | m > 
!        faux(:)=rvec(:)**(1+lambda)*ugs(:,n)*wfeig(i,m,:) !the gs have already an r
       faux(1:nr)=rvec(1:nr)**(lambda)
     &            *ugs(1:nr,n)*wfc(jset,i,m,1:nr) !the gs have already an r
! AMoro March '22
!       call simc(faux,rlam,1,nr,dr,nr)    !it will improve if we go only until rlast
       call simc(faux,resc,1,nr,dr,nr)    !it will improve if we go only until rlast
       rlam=resc
       mv=zeff*rlam*
     & matml_l(lambdar,sn,qjc(m),jtoti,qji(n),qlir(n),jtot,qj(m),qlr(m))

       
       mv=mv+rlam*gyrov*(ac/(av+ac))**(lambda-1d0)*
     &matml_sv(lambdar,sn,qjc(m),jtoti,qji(n),qlir(n),jtot,qj(m),qlr(m))
     
       mv=mv+rlam*gyroc*(-av/(av+ac))**(lambda-1d0)*
     &matml_sc(lambdar,sn,qjc(m),jtoti,qji(n),qlir(n),jtot,qj(m),qlr(m))

       mvc=mvc+mv
! 2.0.4
!        mel(iord(i),m)=mel(iord(i),m)+cmplx(zeff*matel(
!     &lambdar,sn,qjc(m),jtoti,qji(n),qlir(n),jtot,qj(m),qlr(m))*
!     &res,0d0)+cmplx(Elamcore,0d0)

!        mel(iord(i),m)=mel(iord(i),m)+mvc   !it seems we have to add better mv or wait untill finish m cycle
!        write(*,*)'1) i,m=',i,m
        mml(i,m)=mml(i,m)+mv


         if (verb.ge.4) 
     &     write(99,'(4x,"Valence:",2f4.1,"->",2f4.1,1f9.4)') 
     &     qlir(n),qji(n),qlr(m),qj(m),mv
         endif
         endif
         enddo !n
        enddo  ! m

! AMoro 2.0.5
!        BEl=abs(Elam)**2*(2*jtot+1)/(2*jtoti+1)
        Bml=mvc**2*(2*jtot+1)/(2*jtoti+1)
        mlam=0d0
c        written(95)=.true.
        write(95,490) edisc(i),Bml/wid,Bml,wid
490     format(1x,f8.4,3f12.6)
        bmsum=bmsum+Bml
        write(99,*)'i,nex,edisc',i,nex,edisc(i),engs
!        apol = apol + e2*(8.*pi/9.)*Bel/(edisc(i)-engs)   
        write(*,300)edisc(i),lambda,i,Bml,2*lambda
300     format(4x,"Ex=",f8.4,"  B(M",i1,"; gs ->",i4,
     &")=",f8.4," e2.fm",i1)
!        test=0d0
        enddo !i=1,nex

	written(95)=.true.

        write(*,*)'- From Pseudo-States:'
        write(*,'(30x,"Total BM=",1f10.6," e2.fm",i1)')bmsum,2*lambda
!        write(*,*)'- From Sum Rule: (only for SP excitation!)'
!        write(*,'(1x,"- Polarizability=",1f10.6," fm^3")') apol
!        call sumrule
!        write(*,'(30x,"Total Bm=",1f10.6," e2.fm",i1)')besr,2*lambda 
!        write(*,*)"(** sum rule only intended for sp excitations)"


c------------------------------------------------
c B(Elambda) from scattering states (version by AMM)
c------------------------------------------------ 
      write(*,*)'Computation of B(Mlambda) from scattering states:'
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
      bmsum=0d0
      if (.not.allocated(dbde)) allocate(dbde(il,nk))
      dbde(:,:)=0.
      do iil=ili,il ! incoming channels
        
        if (allocated(wfscat)) deallocate(wfscat)
         allocate(wfscat(nk,nchan,nr))

         write(96,*)'# inc. channel=',iil
         excore =jpiset(jset)%exc(iil) ! core exc. energy for this inc. chan.
         call wfrange(jset,nchan,iil,emin,emax,nk,energy,
     &             wfscat(:,:,:),smate,psh,writewf)
         wfscat(:,:,:)=sqrt(2./pi)*wfscat(:,:,:) ! set norm to delta(k-k')

         write(96,*)'#E                               B(Ml)                   
     &   xs(photodis)                 xs(cap)                  S_fact_pd  
     &                  S_fact_cap             B(Ml)'
      do ik=1,nk
         mlamcont=zero
         mlamcore=0d0
         mmcc=(0d0,0d0)
         Bml=0d0
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
!        rme=Elamcorem*cleb(qjcir,0d0,lambdar,0d0,qjcr,0d0)*    ! this is only for K=0, but this is the same. (old version)
!     &      sqrt((2*qjcir+1)/(2*qjcr+1))*(-1)**(2*lambda)
        rmm=0d0
        case(1) !----------->  external M(E) 
        rmm=mmc(cindexi(n),cindex(m),lambda)
        end select 

!        Elamcore=geom*resc*rme*sqrt(pi/2)*
!     &           (4*pi)/kcont
        mmcc=geom*resc*rmm*sqrt(pi/2)*
     &           (4*pi)/kcont
!commented by AMM (2.1)
!        Elamcont=Elamcont  + Elamcore
        mlamcont=mlamcont + mmcc
        itran=itran+1
        if(verb.ge.3) then
        write(99,'(4x,i3,"Core(C):",1f4.1,1i2,1f4.1,"->",
     & 1f4.1,1i2,1f4.1,6f9.3)') 
     & itran,qjcir,qli(n),qji(n),qjcr,ql(m),qj(m),
     & geom,abs(resc),rmm,mmcc
        endif
        mlamcore=0d0

        endif
        endif

c VALENCE contribution
      if ((.not.fail3(qji(n),lambdar,qj(m)))        !lamda is enough to connect the two states
     &.and.(cindexi(n).eq.cindex(m))) then          !we avoid different channels with same Jpi
         if (.not.fail3(qlir(n),lambdar-1,qlr(m))) then                

!       coef=zeff*matel(lambdar,sn,qjc(m),jtoti,qji(n),
!     &                 qlir(n),jtot,qj(m),qlr(m))
        
        
         coef=zeff*
     & matml_l(lambdar,sn,qjc(m),jtoti,qji(n),qlir(n),jtot,qj(m),qlr(m))
       
       coef=coef+gyrov*(ac/(av+ac))**(lambda-1d0)*
     &matml_sv(lambdar,sn,qjc(m),jtoti,qji(n),qlir(n),jtot,qj(m),qlr(m))
     
       coef=coef+gyroc*(-av/(av+ac))**(lambda-1d0)*
     &matml_sc(lambdar,sn,qjc(m),jtoti,qji(n),qlir(n),jtot,qj(m),qlr(m))
         
                
!       write(777,*) coef,zeff,lambdar,sn,qjc(m),jtoti,qji(n),qlir(n),
!     &  jtot,qj(m),qlr(m),matml_l(lambdar,sn,qjc(m),jtoti,qji(n),qlir(n)
!     &,jtot,qj(m),qlr(m)),matml_sv(lambdar,sn,qjc(m),jtoti,qji(n),
!     &qlir(n),jtot,qj(m),qlr(m))*gyrov,matml_sc(lambdar,sn,qjc(m),jtoti,
!     &qji(n),qlir(n),jtot,qj(m),qlr(m))*gyroc
!      write(777,*)
         
         
c < n | r^lambda | m > 
       gaux(1:nr)=rvec(1:nr)**(lambda-1d0)
     &           *ugs(1:nr,n)*wfscat(ik,m,1:nr)
!     &           *ugs(1:nr,n)*wfcont(ik,iil,m,1:nr)*sqrt(dkrdk)
       call simc(gaux,resc,1,nr,dr,nr)   

       mlamcont = mlamcont  
     &          + coef*resc*sqrt(pi/2)*(4*pi)/kcont
         endif
         endif
        enddo !n
        enddo  !m

        Bml=Bml+ (2*jtot+1)*abs(mlamcont)**2/(2*jtoti+1)
     &      *mu12*kcont/hc**2/(2*pi)**3
        dbde(iil,ik)=dbde(iil,ik)
     &      + (2*jtot+1)*abs(mlamcont)**2/(2*jtoti+1)
     &      *mu12*kcont/hc**2/(2*pi)**3
        mlam=0d0
!        write(0,*)'engs=',engs, 'jtoti=',jtoti, 'excore=',excore
!        kphot=(econt+excore+engs)/hc!MGR 2022
        kphot=(econt+excore-engs)/hc !AMM2022

        eta=conv*zc*zv*e2/kcont/2
*
* Photo-absorption cross section (see e.g. Thompson Nunes (14.4.3))
*         
        photoxs=8*pi**3*(lambda+1d0)/(lambda*exp(2*dlfac2(2*lambda+1)))*
     &  kphot**(2d0*lambda-1d0)*Bml*mun2        ! fm^2
        photoxs=photoxs*10                    ! fm^2 -> mb
        sfactphoto=econt*exp(2*pi*eta)*photoxs
        
!        write(0,*)'ill,jc=',iil,jpiset(jset)%jc(iil)
        
*
* Radiative capture cross section
*        
        capxs=2*(2*jtoti+1d0)/(2*sn+1d0)/(2*jpiset(jset)%jc(iil)+1d0)*
     &  kphot**2/kcont**2*photoxs

*
* S-factor
*     
        sfactcap=econt*exp(2*pi*eta)*capxs!MGR 2022
        
        write(96,*) econt+excore,Bml,photoxs,capxs,sfactphoto,sfactcap,
     &  dbde(1,ik)
        bmsum=bmsum+Bml*hc**2*kcont/mu12*dk    !de=de/dk*dk
!        apol= apol +BEl*hc**2*kcont/mu12*dk/(econt-engs)*(8*pi/9.)*e2
        enddo ! ik
        write(*,'(3x, "For inc chan=",i3," B(Ml)=",1f8.3)') iil,bmsum
        write(96,*)'&'
        enddo ! iil
        write(*,*)'- From Continuum wfs:'
        write(*,'(30x,"Total BM=",1f10.6," e^2 fm^",i1)')bmsum,2*lambda
!        write(*,'(30x,"Polarizability=",1f10.6," fm^3")')apol
        write(95,*)'&'
        if (allocated(smate))deallocate(smate)
        if (allocated(psh))deallocate(psh)
        if (allocated(qlr))deallocate(qlr)
        if (allocated(dbde))deallocate(dbde)
        enddo !MGR 2022


c-----------------------------------------------------------------------------
c   B(Elambda) from discrete distribution convoluted with continuous wfs (AMM style)
c-----------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*)'NEED TO UPDATE SOLAPS HERE; RETURNING!!'; RETURN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      bmsum=0d0
      Bml=0d0

      do ik=1,nk
      Bml=0d0
      test=zero
      mlamcont=zero
      kcont=kmin+(kmax-kmin)*(ik-1)/(nk-1)  !new kcont for schcc          
      econt=(hc*kcont)**2/2/mu12
      do iil=ili,il	  
	krel=sqrt(2d0*mu12*(econt-exc(iil)))/hc
	if (krel.gt.eps) then
	  dkrdk=kcont/krel
	else
	  dkrdk=1d0
	endif
	  
       mlamcont=zero
       do i=1,nex
       do n=1,nchan
       mlamcont=mlamcont+solapmat(ik,iil,i,n)*
     &                 mml(i,n)*sqrt(dkrdk)
       enddo
       enddo
       mlamcont=mlamcont*sqrt(pi/2)*(4*pi)/kcont
      Bml=Bml+(2*jtot+1)*abs(mlamcont)**2/(2*jtoti+1)
       enddo
!          write(98,*)econt,BEl
      Bml=Bml*mu12*kcont/hc**2/(2*pi)**3
      bmsum=bmsum+Bml*hc**2*kcont/mu12*(kmax-kmin)/(nk-1) 
      write(97,*)econt,Bml
      enddo

      written(97)=.true.
      write(*,*)'- From Folding:'
      write(*,'(30x,"Total BE=",1f10.6," e^2 fm^",i1)')bmsum,2*lambda


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




        function matml_l(lambda,sn,ic,jtoti,jni,li,jtotf,jnf,lf)
        use constants, only: pi
        implicit none
        real*8:: lambda,sn,ic,jtoti,jni,li,jtotf,jnf,lf   !they have to be all real at the end to fit the 3j,6j subroutines.
        real*8:: matml_l,sixj,cleb
        matml_l=(-1)**(-sn+ic+jtoti-jni-jnf-li-lambda)*
     &   sqrt((2*lambda+1)*(2*lambda-1)*lf*(lf+1)/(4*pi))

        matml_l=matml_l*sqrt((2*jtoti+1)*(2*jni+1)*(2*jnf+1)*
     &  (2*lf+1)*(2*li+1))
        matml_l=matml_l*sixj(jtotf,jtoti,lambda,jni,jnf,ic)

        matml_l=matml_l*sixj(jnf,jni,lambda,li,lf,sn)
        
        matml_l=matml_l*sixj(lambda,1d0,lambda-1d0,lf,li,lf)

        matml_l=matml_l*cleb(li,0d0,lambda-1d0,0d0,lf,0d0)

!        matel=(-1)**(jni-jnf)*matel !in case of coupling (I jn) J instead of (jn I) J; our coupling is ((l s) j I ) J
        end

        function matml_sv(lambda,sn,ic,jtoti,jni,li,jtotf,jnf,lf)
        use constants, only: pi
        implicit none
        real*8:: lambda,sn,ic,jtoti,jni,li,jtotf,jnf,lf   !they have to be all real at the end to fit the 3j,6j subroutines.
        real*8:: matml_sv,sixj,cleb,u9
        matml_sv=(-1)**(lambda+jtoti+jnf-ic)*
     &   sqrt((2*lambda+1)*(2*lambda-1)*sn*(sn+1)/(4*pi))

        matml_sv=matml_sv*sqrt((2*jtoti+1)*(2*jni+1)*(2*jnf+1)*
     &  (2*sn+1)*(2*li+1))
        matml_sv=matml_sv*sixj(jtotf,jtoti,lambda,jni,jnf,ic)

        matml_sv=matml_sv*u9(li,lambda-1d0,lf,sn,1d0,sn,jni,lambda,jnf)

        matml_sv=matml_sv*cleb(li,0d0,lambda-1d0,0d0,lf,0d0)

!        matel=(-1)**(jni-jnf)*matel !in case of coupling (I jn) J instead of (jn I) J; our coupling is ((l s) j I ) J
        end
        

        function matml_sc(lambda,sn,ic,jtoti,jni,li,jtotf,jnf,lf)
        use constants, only: pi
        implicit none
        real*8:: lambda,sn,ic,jtoti,jni,li,jtotf,jnf,lf   !they have to be all real at the end to fit the 3j,6j subroutines.
        real*8:: matml_sc,sixj,cleb,u9
        matml_sc=(-1)**(lambda+1d0-jni+2*jnf-sn+lf)*
     &   sqrt((2*lambda+1)*(2*lambda-1)*sn*(sn+1)/(4*pi))

        matml_sc=matml_sc*sqrt((2*jtoti+1)*(2*jni+1)*(2*jnf+1)*
     &  (2*ic+1)*(2*li+1))
        matml_sc=matml_sc*sixj(lf,li,lambda,jni,jnf,sn)

        matml_sc=matml_sc*
     &   u9(jnf,ic,jtotf,jni,ic,jtoti,lambda-1d0,1d0,lambda)

        matml_sc=matml_sc*cleb(li,0d0,lambda-1d0,0d0,lf,0d0)

!        matel=(-1)**(jni-jnf)*matel !in case of coupling (I jn) J instead of (jn I) J; our coupling is ((l s) j I ) J
        end


       subroutine stone_gyro(z,a,gyro)
       
       implicit none
       
       logical :: file_exists
       character(len=100) :: filename
       integer i,ios
       real*8 z,a,gyro,kk,readgyro,readz,reada
       
       filename="nuclear_moments.csv"
       
       inquire(file=filename, exist=file_exists)

       if (file_exists) then
        open(10,file=filename,status="old")
        read(10,*)!header
        readgyro=0d0
        do i=1,4000
          read(10,*,iostat=ios) readz,reada,kk,kk,readgyro
          if (ios.ne.0) then
        write(*,*) 'Nucleus not found in stone table, g-factor set to 0'
          gyro=0d0
          return
          endif
          if (abs(readz-z).lt.0.01d0 .and. abs(reada-a).lt.0.01d0) then
          gyro=readgyro
          return
          endif
        enddo
       else
        print *, "File nuclear_moments.csv not foung g set to 0"
        gyro=0d0
       endif
       
       end subroutine
