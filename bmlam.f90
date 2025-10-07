
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
       
      real*8 z,a,gyro,gfact(0:200,0:200)
       
      gfact(   0,   1)=  -3.826085
      gfact(   1,   0)=   5.585695
      gfact(   1,   1)=   0.857438
      gfact(   1,   2)=   5.957925
      gfact(   2,   1)=  -4.254995
      gfact(   3,   3)=   0.822047
      gfact(   3,   4)=   2.170951
      gfact(   3,   5)=   0.826700
      gfact(   3,   6)=   2.292733
      gfact(   3,   8)=   2.445333
      gfact(   4,   5)=  -0.785200
      gfact(   5,   3)=   0.517750
      gfact(   5,   5)=   0.600215
      gfact(   5,   6)=   1.792433
      gfact(   5,   7)=   1.002720
      gfact(   5,   8)=   2.118533
      gfact(   5,   9)=   0.592500
      gfact(   5,  10)=   1.772667
      gfact(   5,  12)=   1.700000
      gfact(   6,   3)=   0.927600
      gfact(   6,   5)=  -0.642667
      gfact(   6,   6)=   0.030000
      gfact(   6,   7)=   1.404824
      gfact(   6,   8)=   0.273333
      gfact(   6,   9)=   2.640000
      gfact(   7,   5)=   0.457300
      gfact(   7,   6)=   0.644400
      gfact(   7,   7)=   0.403761
      gfact(   7,   8)=  -0.566378
      gfact(   7,   9)=   0.533333
      gfact(   7,  10)=   0.704000
      gfact(   8,   5)=   0.926067
      gfact(   8,   7)=   1.439020
      gfact(   8,   8)=   0.556000
      gfact(   8,   9)=  -0.757516
      gfact(   8,  10)=  -0.285000
      gfact(   8,  11)=   0.612780
      gfact(   8,  12)=   0.350000
      gfact(   9,   8)=   1.888520
      gfact(   9,   9)=   0.533333
      gfact(   9,  10)=   5.257736
      gfact(   9,  11)=   1.046675
      gfact(   9,  12)=   1.572000
      gfact(  10,   9)=  -3.770840
      gfact(  10,  10)=   0.540000
      gfact(  10,  11)=  -0.441198
      gfact(  10,  12)=   0.325000
      gfact(  10,  13)=  -0.432000
      gfact(  11,   9)=   0.184700
      gfact(  11,  10)=   1.890867
      gfact(  11,  11)=   0.582000
      gfact(  11,  12)=   1.478348
      gfact(  11,  13)=   0.422575
      gfact(  11,  14)=   1.473200
      gfact(  11,  15)=   0.950333
      gfact(  11,  16)=   1.558000
      gfact(  11,  17)=   2.426000
      gfact(  11,  18)=   1.632667
      gfact(  11,  19)=   1.041500
      gfact(  11,  20)=   1.536667
      gfact(  12,  11)=   0.357600
      gfact(  12,  12)=   0.510000
      gfact(  12,  12)=  -0.145000
      gfact(  12,  13)=  -0.342180
      gfact(  12,  14)=   0.500000
      gfact(  13,  12)=   1.458200
      gfact(  13,  13)=   0.560800
      gfact(  13,  14)=   1.456603
      gfact(  13,  15)=   1.080667
      gfact(  14,  13)=   0.342160
      gfact(  14,  14)=   0.550000
      gfact(  14,  15)=  -1.110580
      gfact(  14,  16)=   0.400000
      gfact(  14,  18)=  -0.080000
      gfact(  14,  19)=   0.806667
      gfact(  15,  14)=   2.469800
      gfact(  15,  16)=   2.263200
      gfact(  15,  16)=  -0.252400
      gfact(  16,  15)=   0.975860
      gfact(  16,  16)=   0.450000
      gfact(  16,  17)=   0.429214
      gfact(  16,  18)=   0.500000
      gfact(  16,  19)=   0.666667
      gfact(  17,  16)=   0.501333
      gfact(  17,  18)=   0.547916
      gfact(  17,  19)=   0.642735
      gfact(  17,  20)=   0.456082
      gfact(  17,  21)=   1.025000
      gfact(  18,  15)=  -1.446000
      gfact(  18,  17)=   0.422000
      gfact(  18,  18)=   0.055000
      gfact(  18,  19)=   0.533333
      gfact(  18,  21)=  -0.453714
      gfact(  18,  22)=  -0.100000
      gfact(  19,  17)=   0.274000
      gfact(  19,  18)=   0.135473
      gfact(  19,  19)=   0.457000
      gfact(  19,  20)=   0.260980
      gfact(  19,  21)=  -0.324525
      gfact(  19,  22)=   0.143247
      gfact(  19,  23)=  -0.571250
      gfact(  19,  24)=   0.108867
      gfact(  19,  25)=  -0.428000
      gfact(  19,  26)=   0.115600
      gfact(  19,  27)=  -0.525500
      gfact(  19,  28)=   3.866000
      gfact(  20,  19)=   0.681120
      gfact(  20,  20)=   0.533333
      gfact(  20,  21)=  -0.455652
      gfact(  20,  22)=  -0.095000
      gfact(  20,  23)=  -0.376371
      gfact(  20,  24)=  -0.300000
      gfact(  20,  25)=  -0.379257
      gfact(  20,  27)=  -0.394286
      gfact(  21,  20)=   1.551714
      gfact(  21,  22)=   1.320000
      gfact(  21,  23)=   1.280000
      gfact(  21,  24)=   1.358996
      gfact(  21,  25)=   0.757500
      gfact(  21,  26)=   1.525714
      gfact(  22,  21)=   0.242857
      gfact(  22,  23)=   0.027143
      gfact(  22,  24)=   0.500000
      gfact(  22,  25)=  -0.315392
      gfact(  22,  26)=   0.450000
      gfact(  22,  27)=  -0.315477
      gfact(  22,  28)=   1.350000
      gfact(  23,  23)=   0.546667
      gfact(  23,  25)=   0.503000
      gfact(  23,  26)=   1.277143
      gfact(  23,  27)=   0.557615
      gfact(  23,  28)=   1.471059
      gfact(  24,  25)=   0.190400
      gfact(  24,  26)=   0.600000
      gfact(  24,  27)=   0.266857
      gfact(  24,  28)=   1.500000
      gfact(  24,  29)=  -0.316360
      gfact(  24,  30)=   0.550000
      gfact(  25,  26)=   1.427320
      gfact(  25,  27)=   0.510367
      gfact(  25,  28)=   1.435429
      gfact(  25,  29)=   1.093967
      gfact(  25,  30)=   1.381280
      gfact(  25,  31)=   1.075533
      gfact(  26,  27)=  -0.257333
      gfact(  26,  28)=   1.700000
      gfact(  26,  29)=   1.080000
      gfact(  26,  30)=   0.610000
      gfact(  26,  31)=   0.180880
      gfact(  26,  32)=   0.450000
      gfact(  26,  33)=  -0.223867
      gfact(  27,  28)=   1.377714
      gfact(  27,  29)=   0.962500
      gfact(  27,  30)=   1.348571
      gfact(  27,  31)=   2.022000
      gfact(  27,  32)=   1.322000
      gfact(  27,  33)=   0.759800
      gfact(  28,  29)=  -0.531667
      gfact(  28,  30)=  -0.050000
      gfact(  28,  31)=   0.140000
      gfact(  28,  32)=   0.100000
      gfact(  28,  33)=  -0.500013
      gfact(  28,  34)=   0.340000
      gfact(  28,  35)=   0.300800
      gfact(  28,  36)=   0.450000
      gfact(  28,  37)=   0.276000
      gfact(  29,  31)=   0.609500
      gfact(  29,  32)=   1.426667
      gfact(  29,  33)=  -0.380000
      gfact(  29,  34)=   1.484804
      gfact(  29,  35)=  -0.217000
      gfact(  29,  36)=   1.587800
      gfact(  29,  37)=  -0.282000
      gfact(  30,  33)=  -0.187760
      gfact(  30,  34)=   0.400000
      gfact(  30,  35)=   0.307600
      gfact(  30,  36)=   0.250000
      gfact(  30,  37)=   0.350192
      gfact(  30,  38)=   0.500000
      gfact(  30,  39)=   0.257111
      gfact(  30,  40)=   0.300000
      gfact(  30,  41)=   0.233778
      gfact(  31,  35)=   0.505000
      gfact(  31,  36)=   1.233800
      gfact(  31,  37)=   0.011750
      gfact(  31,  38)=   1.344393
      gfact(  31,  39)=  -0.065000
      gfact(  31,  40)=   1.708180
      gfact(  31,  41)=  -0.044080
      gfact(  32,  35)=  -0.188667
      gfact(  32,  36)=   0.400000
      gfact(  32,  37)=   0.294000
      gfact(  32,  38)=   0.470000
      gfact(  32,  39)=   1.094000
      gfact(  32,  40)=   0.400000
      gfact(  32,  41)=  -0.195437
      gfact(  32,  42)=   0.435000
      gfact(  32,  43)=   1.020000
      gfact(  32,  44)=   0.420000
      gfact(  33,  35)=   0.230000
      gfact(  33,  36)=   0.632000
      gfact(  33,  37)=   0.526525
      gfact(  33,  38)=   0.669600
      gfact(  33,  39)=  -1.078300
      gfact(  33,  40)=   0.652000
      gfact(  33,  41)=  -0.798500
      gfact(  33,  42)=   0.959653
      gfact(  33,  43)=  -0.453000
      gfact(  33,  44)=   0.296000
      gfact(  34,  39)=   0.193333
      gfact(  34,  40)=  -0.180000
      gfact(  34,  41)=   0.268000
      gfact(  34,  42)=   0.400000
      gfact(  34,  43)=   1.070084
      gfact(  34,  44)=   0.400000
      gfact(  34,  45)=  -0.290857
      gfact(  34,  46)=   0.400000
      gfact(  34,  48)=   0.450000
      gfact(  35,  37)=   0.200000
      gfact(  35,  38)=   1.313333
      gfact(  35,  39)=   0.420000
      gfact(  35,  40)=   0.506667
      gfact(  35,  41)=   0.548210
      gfact(  35,  42)=   0.613333
      gfact(  35,  43)=   0.130000
      gfact(  35,  44)=   1.404267
      gfact(  35,  45)=   0.514000
      gfact(  35,  46)=   1.513708
      gfact(  35,  47)=   0.325400
      gfact(  35,  49)=   0.950000
      gfact(  36,  39)=  -0.212400
      gfact(  36,  41)=  -0.233200
      gfact(  36,  42)=   0.540000
      gfact(  36,  43)=   1.072000
      gfact(  36,  45)=  -0.259429
      gfact(  36,  47)=  -0.215704
      gfact(  36,  48)=  -0.246250
      gfact(  36,  49)=  -0.223333
      gfact(  36,  51)=  -0.407200
      gfact(  36,  53)=  -0.220000
      gfact(  36,  55)=  -0.233200
      gfact(  36,  57)=  -0.826000
      gfact(  36,  59)=  -0.820000
      gfact(  37,  39)=  -0.372623
      gfact(  37,  40)=   0.436312
      gfact(  37,  41)=   0.637250
      gfact(  37,  42)=   1.343160
      gfact(  37,  42)=   1.117778
      gfact(  37,  43)=  -0.083600
      gfact(  37,  44)=   1.373000
      gfact(  37,  45)=   0.554508
      gfact(  37,  46)=   0.569960
      gfact(  37,  47)=  -0.662058
      gfact(  37,  48)=   0.541192
      gfact(  37,  49)=  -0.846000
      gfact(  37,  50)=   1.834207
      gfact(  37,  51)=   0.254000
      gfact(  37,  52)=   1.589067
      gfact(  37,  53)=   0.538667
      gfact(  37,  54)=   1.454333
      gfact(  37,  56)=   0.564000
      gfact(  37,  57)=   0.499333
      gfact(  37,  58)=   0.533600
      gfact(  37,  59)=   0.733000
      gfact(  37,  60)=   1.227333
      gfact(  38,  39)=  -0.139200
      gfact(  38,  41)=  -0.316000
      gfact(  38,  43)=   1.086000
      gfact(  38,  44)=   0.400000
      gfact(  38,  45)=  -0.236857
      gfact(  38,  46)=   0.420000
      gfact(  38,  47)=  -0.222222
      gfact(  38,  48)=   0.275000
      gfact(  38,  49)=  -0.242844
      gfact(  38,  50)=   1.150000
      gfact(  38,  51)=  -0.458800
      gfact(  38,  53)=  -0.354000
      gfact(  38,  55)=  -0.317200
      gfact(  38,  57)=  -1.074000
      gfact(  38,  59)=  -0.996000
      gfact(  38,  60)=   0.380000
      gfact(  38,  61)=  -0.174000
      gfact(  39,  44)=   0.600000
      gfact(  39,  46)=   1.377778
      gfact(  39,  47)=   0.150000
      gfact(  39,  48)=   1.346667
      gfact(  39,  49)=   0.608750
      gfact(  39,  50)=  -0.274831
      gfact(  39,  51)=  -0.815000
      gfact(  39,  52)=   0.328200
      gfact(  40,  48)=  -0.226250
      gfact(  40,  49)=  -0.240000
      gfact(  40,  50)=   1.250000
      gfact(  40,  51)=  -0.521448
      gfact(  40,  52)=  -0.030000
      gfact(  40,  54)=  -0.260000
      gfact(  40,  55)=   0.452000
      gfact(  40,  57)=   0.391429
      gfact(  40,  59)=   0.280000
      gfact(  40,  60)=   0.260000
      gfact(  41,  48)=   1.381333
      gfact(  41,  49)=   0.620125
      gfact(  41,  50)=   1.406154
      gfact(  41,  51)=   3.068500
      gfact(  41,  52)=   1.371222
      gfact(  41,  54)=   1.364667
      gfact(  41,  55)=   0.829333
      gfact(  41,  56)=   1.367333
      gfact(  42,  47)=   0.790476
      gfact(  42,  48)=   1.100000
      gfact(  42,  49)=   0.839048
      gfact(  42,  50)=   1.412500
      gfact(  42,  51)=   0.945714
      gfact(  42,  52)=  -0.065000
      gfact(  42,  53)=  -0.365680
      gfact(  42,  54)=  -0.100000
      gfact(  42,  55)=  -0.373400
      gfact(  42,  56)=   0.350000
      gfact(  42,  57)=   0.750000
      gfact(  42,  58)=   0.350000
      gfact(  42,  60)=   0.420000
      gfact(  42,  62)=   0.200000
      gfact(  43,  49)=   0.809091
      gfact(  43,  50)=   1.404444
      gfact(  43,  51)=   0.731429
      gfact(  43,  52)=   1.320000
      gfact(  43,  53)=   0.727143
      gfact(  43,  56)=   1.263267
      gfact(  44,  49)=   0.854286
      gfact(  44,  50)=   1.353333
      gfact(  44,  51)=   0.344400
      gfact(  44,  52)=  -0.065000
      gfact(  44,  53)=  -0.314800
      gfact(  44,  54)=   0.400000
      gfact(  44,  55)=  -0.256400
      gfact(  44,  56)=   0.510000
      gfact(  44,  57)=  -0.287600
      gfact(  44,  58)=   0.370000
      gfact(  44,  59)=   0.137333
      gfact(  44,  60)=   0.410000
      gfact(  44,  61)=  -0.213333
      gfact(  45,  50)=   1.282353
      gfact(  45,  54)=   1.248889
      gfact(  45,  55)=   2.162000
      gfact(  45,  56)=   1.206667
      gfact(  45,  57)=   0.250000
      gfact(  45,  58)=  -1.768000
      gfact(  45,  60)=   1.260000
      gfact(  45,  61)=   2.575000
      gfact(  46,  50)=   1.371250
      gfact(  46,  55)=   0.264000
      gfact(  46,  56)=   0.410000
      gfact(  46,  57)=  -0.190909
      gfact(  46,  58)=   0.460000
      gfact(  46,  59)=  -0.256800
      gfact(  46,  62)=   0.360000
      gfact(  46,  64)=   0.310000
      gfact(  47,  54)=   1.266667
      gfact(  47,  55)=   0.920000
      gfact(  47,  56)=   1.277143
      gfact(  47,  57)=   0.783400
      gfact(  47,  58)=   0.202800
      gfact(  47,  59)=   2.900000
      gfact(  47,  60)=  -0.227140
      gfact(  47,  61)=   2.688400
      gfact(  47,  62)=   0.261120
      gfact(  47,  63)=   2.727100
      gfact(  47,  64)=  -0.292000
      gfact(  47,  65)=   0.027350
      gfact(  47,  66)=   0.318000
      gfact(  48,  52)=   1.237500
      gfact(  48,  54)=   1.287500
      gfact(  48,  55)=  -0.324000
      gfact(  48,  57)=  -0.295720
      gfact(  48,  58)=   0.400000
      gfact(  48,  59)=  -0.246022
      gfact(  48,  60)=   0.350000
      gfact(  48,  61)=  -0.331138
      gfact(  48,  62)=   0.285000
      gfact(  48,  63)=  -1.189772
      gfact(  48,  64)=   0.300000
      gfact(  48,  65)=  -1.244602
      gfact(  48,  66)=   0.290000
      gfact(  48,  67)=  -1.296852
      gfact(  48,  68)=   0.300000
      gfact(  49,  55)=   0.888000
      gfact(  49,  56)=   1.261111
      gfact(  49,  57)=   0.702286
      gfact(  49,  58)=   1.241111
      gfact(  49,  59)=   0.651571
      gfact(  49,  60)=   1.230667
      gfact(  49,  61)=   2.182500
      gfact(  49,  62)=   1.222889
      gfact(  49,  63)=   2.820000
      gfact(  49,  64)=   1.228644
      gfact(  49,  65)=   2.817000
      gfact(  49,  66)=   1.231289
      gfact(  49,  67)=   2.787600
      gfact(  49,  68)=   1.226444
      gfact(  49,  69)=   0.846200
      gfact(  49,  70)=   1.225556
      gfact(  49,  71)=   0.859000
      gfact(  49,  72)=   1.222667
      gfact(  49,  73)=   0.863600
      gfact(  49,  74)=   1.220222
      gfact(  49,  75)=   1.347667
      gfact(  49,  76)=   1.222667
      gfact(  49,  77)=   1.344667
      gfact(  49,  78)=   1.227111
      gfact(  50,  58)=  -0.040000
      gfact(  50,  59)=  -0.431600
      gfact(  50,  60)=   0.011667
      gfact(  50,  61)=   0.173714
      gfact(  50,  62)=   0.350000
      gfact(  50,  63)=  -1.758200
      gfact(  50,  65)=  -1.837660
      gfact(  50,  66)=  -0.150000
      gfact(  50,  67)=  -2.002080
      gfact(  50,  68)=   0.020000
      gfact(  50,  69)=  -2.094560
      gfact(  50,  70)=   0.011000
      gfact(  50,  71)=   0.465200
      gfact(  50,  72)=  -0.050000
      gfact(  50,  73)=  -0.249091
      gfact(  50,  74)=  -0.150000
      gfact(  50,  75)=  -0.245091
      gfact(  51,  61)=   0.274000
      gfact(  51,  63)=   0.573333
      gfact(  51,  64)=   1.384000
      gfact(  51,  65)=   0.905000
      gfact(  51,  66)=   1.372000
      gfact(  51,  67)=   2.470000
      gfact(  51,  68)=   1.380000
      gfact(  51,  69)=   2.300000
      gfact(  51,  70)=   1.345360
      gfact(  51,  71)=   0.950000
      gfact(  51,  72)=   0.728514
      gfact(  51,  73)=   0.400000
      gfact(  51,  74)=   0.751429
      gfact(  51,  75)=   0.160000
      gfact(  51,  76)=   0.770571
      gfact(  51,  77)=   0.162500
      gfact(  51,  78)=   0.797143
      gfact(  51,  80)=   0.825714
      gfact(  51,  82)=   0.857143
      gfact(  52,  63)=  -0.173455
      gfact(  52,  65)=  -0.314800
      gfact(  52,  67)=   0.500000
      gfact(  52,  68)=   0.390000
      gfact(  52,  69)=   0.162727
      gfact(  52,  70)=   0.330000
      gfact(  52,  71)=  -1.473896
      gfact(  52,  72)=   0.280000
      gfact(  52,  73)=  -1.777010
      gfact(  52,  74)=   0.310000
      gfact(  52,  75)=   0.423333
      gfact(  52,  76)=   0.250000
      gfact(  52,  77)=   0.468000
      gfact(  52,  78)=   0.290000
      gfact(  52,  79)=   0.464000
      gfact(  52,  80)=   0.783333
      gfact(  52,  82)=   0.846667
      gfact(  52,  83)=  -0.400000
      gfact(  53,  64)=   1.240000
      gfact(  53,  65)=   1.000000
      gfact(  53,  66)=   1.160000
      gfact(  53,  67)=   0.615000
      gfact(  53,  68)=   0.920000
      gfact(  53,  69)=   0.940000
      gfact(  53,  70)=   1.127200
      gfact(  53,  71)=   0.723000
      gfact(  53,  72)=   1.128400
      gfact(  53,  73)=   0.719000
      gfact(  53,  74)=   1.125308
      gfact(  53,  75)=  -0.180000
      gfact(  53,  76)=   0.748857
      gfact(  53,  77)=   0.669800
      gfact(  53,  78)=   0.783429
      gfact(  53,  79)=   0.772000
      gfact(  53,  80)=   0.816000
      gfact(  54,  63)=  -0.237520
      gfact(  54,  65)=  -0.261680
      gfact(  54,  67)=  -0.280400
      gfact(  54,  69)=  -0.300000
      gfact(  54,  70)=   0.230000
      gfact(  54,  71)=  -0.538000
      gfact(  54,  72)=   0.370000
      gfact(  54,  73)=  -1.006600
      gfact(  54,  74)=   0.410000
      gfact(  54,  75)=  -1.555952
      gfact(  54,  76)=   0.380000
      gfact(  54,  77)=   0.461000
      gfact(  54,  78)=   0.370000
      gfact(  54,  79)=   0.541933
      gfact(  54,  80)=   0.550000
      gfact(  54,  81)=   0.602133
      gfact(  54,  82)=   1.200000
      gfact(  54,  83)=  -0.276571
      gfact(  54,  85)=  -0.202667
      gfact(  54,  87)=   0.004000
      gfact(  54,  89)=  -0.183960
      gfact(  55,  63)=   1.938000
      gfact(  55,  64)=   1.213333
      gfact(  55,  65)=   1.935000
      gfact(  55,  66)=   0.513333
      gfact(  55,  67)=  -0.133300
      gfact(  55,  68)=   2.754000
      gfact(  55,  69)=   0.673000
      gfact(  55,  70)=   2.818000
      gfact(  55,  71)=   0.777000
      gfact(  55,  72)=   2.918000
      gfact(  55,  73)=   0.974000
      gfact(  55,  74)=   2.982000
      gfact(  55,  75)=   1.460000
      gfact(  55,  76)=   1.412000
      gfact(  55,  77)=   1.111000
      gfact(  55,  78)=   0.737721
      gfact(  55,  79)=   0.748425
      gfact(  55,  80)=   0.780686
      gfact(  55,  81)=   0.742200
      gfact(  55,  82)=   0.814657
      gfact(  55,  83)=   0.233333
      gfact(  55,  84)=   0.770286
      gfact(  55,  85)=   0.133895
      gfact(  55,  86)=   0.696571
      gfact(  55,  88)=   0.580000
      gfact(  55,  89)=  -0.546000
      gfact(  55,  90)=   0.522667
      gfact(  55,  91)=  -0.515000
      gfact(  56,  65)=   0.264000
      gfact(  56,  67)=  -0.272000
      gfact(  56,  69)=   0.354000
      gfact(  56,  71)=   0.166800
      gfact(  56,  73)=  -0.800000
      gfact(  56,  74)=   0.350000
      gfact(  56,  75)=   1.416226
      gfact(  56,  76)=   0.340000
      gfact(  56,  77)=   1.543340
      gfact(  56,  78)=   0.430000
      gfact(  56,  79)=   0.558627
      gfact(  56,  80)=   0.345000
      gfact(  56,  81)=   0.624913
      gfact(  56,  82)=   0.700000
      gfact(  56,  83)=  -0.278000
      gfact(  56,  85)=  -0.224667
      gfact(  56,  86)=   0.425000
      gfact(  56,  87)=   0.177200
      gfact(  56,  88)=   0.340000
      gfact(  56,  89)=  -0.114000
      gfact(  56,  90)=   0.280000
      gfact(  57,  76)=   1.363636
      gfact(  57,  78)=   0.000000
      gfact(  57,  80)=   0.770000
      gfact(  57,  81)=   0.742729
      gfact(  57,  82)=   0.795156
      gfact(  57,  83)=   0.243333
      gfact(  58,  68)=   1.000000
      gfact(  58,  76)=  -0.187000
      gfact(  58,  77)=  -0.069474
      gfact(  58,  78)=  -0.180000
      gfact(  58,  79)=   0.640000
      gfact(  58,  80)=  -0.170000
      gfact(  58,  81)=   0.706667
      gfact(  58,  82)=   0.950000
      gfact(  58,  83)=   0.311429
      gfact(  58,  84)=   0.210000
      gfact(  58,  85)=   0.286667
      gfact(  58,  88)=   0.240000
      gfact(  58,  90)=   0.370000
      gfact(  59,  77)=   0.575000
      gfact(  59,  80)=   1.200000
      gfact(  59,  82)=   1.710160
      gfact(  59,  83)=   0.117000
      gfact(  59,  84)=   0.771714
      gfact(  59,  85)=  -1.200000
      gfact(  60,  74)=   0.600000
      gfact(  60,  75)=  -0.173333
      gfact(  60,  76)=   1.100000
      gfact(  60,  77)=  -1.266000
      gfact(  60,  78)=  -0.174000
      gfact(  60,  79)=   0.604667
      gfact(  60,  80)=  -0.192000
      gfact(  60,  81)=   0.674667
      gfact(  60,  82)=   0.845000
      gfact(  60,  83)=  -0.304286
      gfact(  60,  84)=   0.160000
      gfact(  60,  85)=  -0.187429
      gfact(  60,  86)=   0.290000
      gfact(  60,  87)=   0.231200
      gfact(  60,  88)=   0.350000
      gfact(  60,  89)=   0.140400
      gfact(  60,  90)=   0.380000
      gfact(  61,  77)=   1.066667
      gfact(  61,  82)=   1.520000
      gfact(  61,  83)=   0.338000
      gfact(  61,  84)=   1.520000
      gfact(  61,  86)=   0.737143
      gfact(  61,  87)=   2.100000
      gfact(  61,  88)=   0.942857
      gfact(  61,  90)=   0.780000
      gfact(  62,  76)=   1.000000
      gfact(  62,  77)=  -1.060000
      gfact(  62,  78)=  -0.180000
      gfact(  62,  79)=  -1.480000
      gfact(  62,  80)=   0.060000
      gfact(  62,  81)=   0.673333
      gfact(  62,  82)=   0.750000
      gfact(  62,  83)=  -0.317143
      gfact(  62,  85)=  -0.232000
      gfact(  62,  86)=   0.255000
      gfact(  62,  87)=  -0.190771
      gfact(  62,  88)=   0.385000
      gfact(  62,  89)=  -0.144440
      gfact(  62,  90)=   0.400000
      gfact(  62,  91)=  -0.014000
      gfact(  62,  92)=   0.390000
      gfact(  62,  93)=   0.753333
      gfact(  63,  75)=   0.883333
      gfact(  63,  76)=   1.109091
      gfact(  63,  77)=   1.365000
      gfact(  63,  78)=   1.397600
      gfact(  63,  79)=   1.540000
      gfact(  63,  80)=   1.469200
      gfact(  63,  81)=   1.893000
      gfact(  63,  82)=   1.599600
      gfact(  63,  83)=   0.355250
      gfact(  63,  84)=   1.494400
      gfact(  63,  85)=   0.468000
      gfact(  63,  86)=   1.430400
      gfact(  63,  87)=   0.541600
      gfact(  63,  88)=   1.388680
      gfact(  63,  89)=  -0.646700
      gfact(  63,  90)=   0.612960
      gfact(  63,  91)=  -0.668333
      gfact(  63,  92)=   0.608000
      gfact(  63,  94)=   0.600000
      gfact(  63,  95)=   1.440000
      gfact(  63,  96)=   0.552000
      gfact(  64,  80)=   1.276000
      gfact(  64,  82)=   0.700000
      gfact(  64,  83)=   0.291429
      gfact(  64,  84)=  -0.017778
      gfact(  64,  85)=   0.251429
      gfact(  64,  87)=   0.220000
      gfact(  64,  88)=   0.480000
      gfact(  64,  89)=   0.253333
      gfact(  64,  90)=   0.480000
      gfact(  64,  91)=  -0.171467
      gfact(  64,  92)=   0.410000
      gfact(  64,  93)=  -0.226533
      gfact(  64,  94)=   0.390000
      gfact(  64,  95)=  -0.293333
      gfact(  64,  96)=   0.360000
      gfact(  65,  82)=   3.400000
      gfact(  65,  83)=  -0.875000
      gfact(  65,  84)=   2.700000
      gfact(  65,  85)=  -0.450000
      gfact(  65,  86)=   1.838000
      gfact(  65,  87)=  -0.290000
      gfact(  65,  88)=   1.376000
      gfact(  65,  89)=   0.533333
      gfact(  65,  90)=   1.340000
      gfact(  65, 101)=   0.566667
      gfact(  65,  92)=   1.340000
      gfact(  65,  93)=   0.586000
      gfact(  65,  94)=   1.342667
      gfact(  65,  95)=   0.596667
      gfact(  65,  96)=   1.466667
      gfact(  66,  81)=  -1.830000
      gfact(  66,  83)=  -0.034000
      gfact(  66,  85)=  -0.270000
      gfact(  66,  86)=   0.552381
      gfact(  66,  87)=  -0.223429
      gfact(  66,  89)=  -0.256667
      gfact(  66,  90)=   0.390000
      gfact(  66,  91)=  -0.200667
      gfact(  66,  92)=   0.360000
      gfact(  66,  93)=  -0.236000
      gfact(  66,  94)=   0.370000
      gfact(  66,  95)=  -0.192000
      gfact(  66,  96)=   0.345000
      gfact(  66,  97)=   0.269200
      gfact(  66,  98)=   0.340000
      gfact(  66,  99)=  -0.148571
      gfact(  67,  85)=  -0.510000
      gfact(  67,  86)=   1.238182
      gfact(  67,  87)=  -0.321500
      gfact(  67,  88)=   1.404000
      gfact(  67,  89)=   0.747500
      gfact(  67,  90)=   1.242857
      gfact(  67,  91)=   0.754000
      gfact(  67,  92)=   1.222857
      gfact(  67,  93)=   0.742000
      gfact(  67,  94)=   1.214286
      gfact(  67,  95)=   0.600000
      gfact(  67,  96)=   1.208571
      gfact(  67,  98)=   1.191429
      gfact(  67,  99)=   0.514286
      gfact(  68,  84)=  -0.075000
      gfact(  68,  85)=  -0.266857
      gfact(  68,  86)=   0.015364
      gfact(  68,  87)=  -0.191143
      gfact(  68,  88)=   0.400000
      gfact(  68,  89)=  -0.274667
      gfact(  68,  90)=   0.360000
      gfact(  68,  91)=  -0.202667
      gfact(  68,  92)=   0.320000
      gfact(  68,  93)=  -0.243333
      gfact(  68,  95)=   0.222800
      gfact(  68,  96)=   0.348500
      gfact(  68,  97)=   0.257200
      gfact(  68,  98)=   0.324500
      gfact(  68,  99)=  -0.161100
      gfact(  68, 100)=   0.310000
      gfact(  68, 101)=   1.040000
      gfact(  68, 102)=   0.316500
      gfact(  68, 103)=   0.263600
      gfact(  69,  87)=   0.200000
      gfact(  69,  88)=   0.952000
      gfact(  69,  89)=   0.020000
      gfact(  69,  90)=   1.368000
      gfact(  69,  91)=   0.160000
      gfact(  69,  92)=   0.685714
      gfact(  69,  93)=   0.068000
      gfact(  69,  94)=  -0.164000
      gfact(  69,  95)=   2.830000
      gfact(  69,  96)=  -0.278000
      gfact(  69,  97)=   0.046000
      gfact(  69,  98)=  -0.394000
      gfact(  69,  99)=   0.075667
      gfact(  69, 100)=  -0.462000
      gfact(  69, 101)=   0.246000
      gfact(  69, 102)=  -0.456000
      gfact(  70,  85)=  -0.240000
      gfact(  70,  87)=  -0.182571
      gfact(  70,  89)=  -0.147200
      gfact(  70,  91)=  -0.218000
      gfact(  70,  93)=  -0.249333
      gfact(  70,  95)=   0.191200
      gfact(  70,  97)=   0.249200
      gfact(  70,  99)=  -0.181429
      gfact(  70, 100)=   0.337000
      gfact(  70, 101)=   0.987340
      gfact(  70, 102)=  -0.575000
      gfact(  70, 103)=  -0.259200
      gfact(  70, 104)=   0.338000
      gfact(  70, 105)=   0.219429
      gfact(  70, 106)=   0.340000
      gfact(  71,  98)=   0.656286
      gfact(  71, 100)=   0.658571
      gfact(  71, 101)=   0.723250
      gfact(  71, 102)=   0.651429
      gfact(  71, 103)=   1.900000
      gfact(  71, 104)=   0.637800
      gfact(  71, 105)=   0.452714
      gfact(  71, 106)=   0.639714
      gfact(  72, 101)=   0.573913
      gfact(  72, 102)=   0.903333
      gfact(  72, 103)=  -0.248000
      gfact(  72, 104)=   0.315000
      gfact(  72, 105)=   0.226714
      gfact(  72, 106)=   0.240000
      gfact(  72, 107)=  -0.142422
      gfact(  72, 108)=   0.305000
      gfact(  73,  98)=   0.690000
      gfact(  73, 100)=   0.680000
      gfact(  73, 102)=   0.648571
      gfact(  73, 104)=   0.642857
      gfact(  73, 105)=   2.740000
      gfact(  73, 106)=   0.654000
      gfact(  73, 107)=   0.536111
      gfact(  73, 108)=   0.677286
      gfact(  73, 109)=   1.006667
      gfact(  73, 110)=   0.674000
      gfact(  74, 104)=   0.250000
      gfact(  74, 105)=   0.400000
      gfact(  74, 106)=   0.255000
      gfact(  74, 108)=   0.260000
      gfact(  74, 109)=   0.235570
      gfact(  74, 110)=   0.289000
      gfact(  74, 112)=   0.310000
      gfact(  74, 113)=   0.414000
      gfact(  75, 104)=   1.120000
      gfact(  75, 105)=   1.600000
      gfact(  75, 106)=   1.276000
      gfact(  75, 107)=   0.405714
      gfact(  75, 108)=   1.267200
      gfact(  75, 109)=   0.843000
      gfact(  75, 110)=   1.274840
      gfact(  75, 111)=   1.739000
      gfact(  75, 112)=   1.287880
      gfact(  75, 113)=   1.788000
      gfact(  76, 106)=   0.424000
      gfact(  76, 107)=   0.176000
      gfact(  76, 108)=  -1.200000
      gfact(  76, 110)=   0.280000
      gfact(  76, 111)=   0.129304
      gfact(  76, 112)=   0.290000
      gfact(  76, 113)=   0.439955
      gfact(  76, 114)=   0.350000
      gfact(  76, 115)=   0.213333
      gfact(  76, 116)=   0.395000
      gfact(  76, 117)=   0.486667
      gfact(  77, 106)=   0.944000
      gfact(  77, 107)=   0.139200
      gfact(  77, 108)=   1.042000
      gfact(  77, 109)=   0.776000
      gfact(  77, 110)=   0.627333
      gfact(  77, 111)=   0.302000
      gfact(  77, 112)=   0.086667
      gfact(  77, 113)=   0.010000
      gfact(  77, 114)=   0.100467
      gfact(  77, 115)=   0.481000
      gfact(  77, 116)=   0.109133
      gfact(  77, 117)=   0.390000
      gfact(  78, 105)=   1.020000
      gfact(  78, 106)=   0.280000
      gfact(  78, 107)=   0.172000
      gfact(  78, 108)=   0.270000
      gfact(  78, 109)=   0.272000
      gfact(  78, 110)=   0.290000
      gfact(  78, 111)=  -0.280667
      gfact(  78, 112)=   0.285000
      gfact(  78, 113)=  -0.334000
      gfact(  78, 114)=   0.285000
      gfact(  78, 115)=   1.206000
      gfact(  78, 116)=   0.300000
      gfact(  78, 117)=   1.219040
      gfact(  78, 118)=   0.295000
      gfact(  78, 119)=   1.020000
      gfact(  78, 120)=   0.315000
      gfact(  79, 104)=   0.788000
      gfact(  79, 105)=   0.414000
      gfact(  79, 106)=   0.868000
      gfact(  79, 107)=  -0.426667
      gfact(  79, 108)=   1.070000
      gfact(  79, 109)=  -0.070000
      gfact(  79, 110)=   0.988000
      gfact(  79, 111)=  -0.065000
      gfact(  79, 112)=   0.091267
      gfact(  79, 113)=  -0.010700
      gfact(  79, 114)=   0.093067
      gfact(  79, 115)=   0.076300
      gfact(  79, 116)=   0.099133
      gfact(  79, 117)=   0.290000
      gfact(  79, 118)=   0.097164
      gfact(  79, 119)=   0.320000
      gfact(  79, 120)=   0.174000
      gfact(  79, 121)=   0.491667
      gfact(  80, 101)=   1.014200
      gfact(  80, 103)=   1.048000
      gfact(  80, 105)=   1.018000
      gfact(  80, 107)=  -0.160615
      gfact(  80, 108)=  -0.168333
      gfact(  80, 109)=  -0.405733
      gfact(  80, 110)=  -0.208333
      gfact(  80, 111)=  -0.412000
      gfact(  80, 113)=  -0.418400
      gfact(  80, 115)=   1.082950
      gfact(  80, 116)=  -0.041429
      gfact(  80, 117)=   1.054749
      gfact(  80, 118)=   0.380000
      gfact(  80, 119)=   1.011771
      gfact(  80, 120)=   0.325000
      gfact(  80, 121)=  -0.373484
      gfact(  80, 122)=   0.390000
      gfact(  80, 123)=   0.339580
      gfact(  80, 124)=   0.450000
      gfact(  80, 125)=   1.201780
      gfact(  80, 126)=   1.090000
      gfact(  81, 106)=   3.100000
      gfact(  81, 107)=   0.069000
      gfact(  81, 108)=   0.861778
      gfact(  81, 109)=   0.127000
      gfact(  81, 110)=   3.176000
      gfact(  81, 111)=   0.100000
      gfact(  81, 112)=   3.182000
      gfact(  81, 113)=   0.070000
      gfact(  81, 114)=   3.160000
      gfact(  81, 115)=   0.036000
      gfact(  81, 116)=   3.160000
      gfact(  81, 117)=   0.000000
      gfact(  81, 118)=   3.200000
      gfact(  81, 119)=   0.020000
      gfact(  81, 120)=   3.210000
      gfact(  81, 121)=   0.030000
      gfact(  81, 122)=   3.244516
      gfact(  81, 123)=   0.045000
      gfact(  81, 124)=   3.276429
      gfact(  81, 125)=   0.854000
      gfact(  81, 126)=   3.752000
      gfact(  81, 127)=   0.058400
      gfact(  82, 109)=  -0.180308
      gfact(  82, 111)=  -0.176923
      gfact(  82, 112)=  -0.066667
      gfact(  82, 113)=  -0.173538
      gfact(  82, 114)=   0.098000
      gfact(  82, 115)=  -0.716667
      gfact(  82, 116)=   0.076000
      gfact(  82, 117)=  -0.716133
      gfact(  82, 118)=  -0.030000
      gfact(  82, 119)=   0.270120
      gfact(  82, 120)=   0.002000
      gfact(  82, 121)=   0.274560
      gfact(  82, 122)=   0.010000
      gfact(  82, 123)=   0.284680
      gfact(  82, 124)=   0.015000
      gfact(  82, 125)=   1.185166
      gfact(  82, 126)=   0.633333
      gfact(  82, 127)=  -0.327444
      gfact(  82, 128)=  -0.311667
      gfact(  82, 129)=  -0.311933
      gfact(  83, 116)=   1.022222
      gfact(  83, 118)=   1.066667
      gfact(  83, 119)=   0.980000
      gfact(  83, 120)=   0.892667
      gfact(  83, 121)=   0.720333
      gfact(  83, 122)=   1.023333
      gfact(  83, 123)=   0.726833
      gfact(  83, 124)=   0.906889
      gfact(  83, 125)=   0.926600
      gfact(  83, 126)=   0.913400
      gfact(  83, 127)=  -0.044510
      gfact(  83, 128)=   1.285714
      gfact(  83, 129)=   0.410000
      gfact(  83, 130)=   0.864444
      gfact(  84, 114)=   0.912500
      gfact(  84, 115)=   0.152308
      gfact(  84, 116)=   0.930000
      gfact(  84, 117)=   0.626667
      gfact(  84, 118)=   0.931250
      gfact(  84, 119)=   0.296000
      gfact(  84, 120)=   0.922500
      gfact(  84, 121)=   0.304000
      gfact(  84, 122)=   0.917500
      gfact(  84, 123)=   0.316000
      gfact(  84, 124)=   0.883333
      gfact(  84, 125)=   1.360000
      gfact(  84, 126)=   0.913333
      gfact(  84, 127)=  -0.050667
      gfact(  85, 122)=   0.300000
      gfact(  85, 123)=   0.269000
      gfact(  85, 124)=   0.952381
      gfact(  85, 125)=   0.890909
      gfact(  85, 126)=   0.910476
      gfact(  85, 127)=   0.540000
      gfact(  85, 132)=   0.844444
      gfact(  86, 117)=  -0.147692
      gfact(  86, 119)=   0.320800
      gfact(  86, 120)=   0.825000
      gfact(  86, 121)=   0.326400
      gfact(  86, 122)=   0.872500
      gfact(  86, 123)=   0.335520
      gfact(  86, 125)=   1.202000
      gfact(  86, 126)=   1.000000
      gfact(  86, 127)=   0.450476
      gfact(  86, 133)=  -0.176800
      gfact(  86, 135)=  -0.005714
      gfact(  86, 136)=   0.460000
      gfact(  86, 137)=   0.460000
      gfact(  86, 139)=  -0.198857
      gfact(  87, 120)=   0.864444
      gfact(  87, 121)=   0.678571
      gfact(  87, 122)=   0.877778
      gfact(  87, 123)=   0.733333
      gfact(  87, 124)=   0.888889
      gfact(  87, 125)=   0.924000
      gfact(  87, 126)=   0.893333
      gfact(  87, 127)=   0.510909
      gfact(  87, 133)=  -0.670000
      gfact(  87, 134)=   0.632000
      gfact(  87, 135)=   0.315000
      gfact(  87, 136)=   0.780000
      gfact(  87, 137)=   0.400000
      gfact(  87, 138)=   0.713333
      gfact(  87, 139)=   0.071200
      gfact(  87, 140)=   3.000000
      gfact(  87, 141)=  -0.380000
      gfact(  88, 121)=   0.346000
      gfact(  88, 123)=   0.351200
      gfact(  88, 124)=   0.887500
      gfact(  88, 125)=   1.226000
      gfact(  88, 126)=   0.885000
      gfact(  88, 127)=   0.733953
      gfact(  88, 133)=  -0.072000
      gfact(  88, 135)=   0.180667
      gfact(  88, 136)=   0.450000
      gfact(  88, 137)=  -1.467600
      gfact(  88, 139)=  -0.269333
      gfact(  88, 141)=   0.201200
      gfact(  89, 126)=   0.920000
      gfact(  89, 128)=   0.851111
      gfact(  89, 138)=   0.733333
      gfact(  90, 139)=   0.184000
      gfact(  91, 137)=   1.166667
      gfact(  91, 139)=   1.000000
      gfact(  91, 140)=   1.005000
      gfact(  91, 142)=   2.666667
      gfact(  92, 143)=  -0.108571
      gfact(  93, 144)=   1.256000
      gfact(  93, 146)=   0.800000
      gfact(  94, 145)=   0.406000
      gfact(  94, 147)=  -0.273200
      gfact(  95, 146)=   0.632000
      gfact(  95, 147)=   0.387900
      gfact(  96, 147)=   0.160000
      gfact(  96, 149)=   0.142857
      gfact(  96, 151)=   0.080000
      gfact(  97, 152)=   0.571429
      gfact(  99, 154)=   1.171429
      gfact(  99, 155)=   1.450000

      gyro=gfact(nint(z),nint(a-z))
       
      end subroutine
