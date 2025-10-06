c ----------------------------------------------------------------
c Double differential x-sections dsigma/dEx dW  (NEW VERSION) 
c ----------------------------------------------------------------
      subroutine d2sigma_alpha(nth,dth,thmin,thcut,emin,emax,ncont,
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
      real*8:: ti,tf,hbarc
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
      real*8  ::maxtab
      integer  :: imax,jmax
      real*8 :: polar(10)
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
      namelist /polarization/ polar
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
      hbarc=hc
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
      kcv=sqrt(2.d0*mucv*ecv)/hc    
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
      polar(:)=0d0
*     ---------------------------------------------------------------
      read(kin,nml=gridener)   !For alpha and beta this will be Ek (the v-C energy)
      read(kin,nml=gridthetac) !For alpha and beta this will be thcm (the vC-T theta) phicm=0
      read(kin,nml=gridthetav) !For alpha this will be thq (the v-C theta), for beta this will be the angle between kcm (vC-T) y q (v-C)
      read(kin,nml=gridphi)    !For alpha and beta this will be phiq (the v-C phi)
      read(kin,nml=polarization)
      if (sum(abs(polar)).lt.1e-2) polar(:)=1
*     ---------------------------------------------------------------
      iten=nint((Enup-Enlow)/dEn)+1
      dEn=(Enup-Enlow)/(iten-1)
*     ---------------------------------------------------------------
*     particle with specified energy defines phi=0 degrees
*     ---------------------------------------------------------------
      idet=1!MGR
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
      
      do ien=1,iten !MGR this energy will be the relative energy between core and valence
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
!$OMP& rnu,mult,xsig,nk,kinset,wkfac,haux,fxytab) 
!$OMP& PRIVATE(ii,maxtab,imax,jmax) 

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

       p1L=sqrt(2.d0*mucv*En)/hbarc

*     ---------------------------------------------------------------
*     the detected particle defines phi=0 degrees
*     ---------------------------------------------------------------
*     loop over azimuth detection angle          
      do 40 ip=1,itphim
      icount= icount+1
      sigphi(ip)=0.d0
      phid=phil+(ip-1)*dphi
      phi=phid*degrad  

      cospv=cos(phi)
      sinpv=sin(phi)

*     -----------------------------------------------------------------
*     construct remaining vectors
*     -----------------------------------------------------------------
!      do j=1,3
!MGR get wavevector modulus for p-t vector (reusing p2L) cm
      p2L=sqrt(kcmi**2+2d0*mupt/hbarc**2*(ebind-excore)-
     & mupt/mucv*p1L**2)
*     wavevector of cm of core and valence particles
      bkp(3)=p2L*costc
      bkp(1)=p2L*sintc
      bkp(2)=0d0
!*     wavevector of relative motion of core and valence particles
      !alpha
      kp(3)=p1L*costv
      kp(1)=p1L*sintv*cospv
      kp(2)=p1L*sintv*sinpv
      !beta
      !kp(3)=p2L*cos(tc+tv)
      !kp(1)=p2L*sin(tc+tv)*cospv
      !kp(2)=p2L*sintv*sinpv
      !enddo
!      eKb=hc**2*dot(bkp,bkp)/(2.d0*mupt)
!      eks=hc**2*dot(kp,kp)/(2.d0*mucv)
      eKb=kcmi**2/(2.d0*mupt)*hc**2+ebind-excore-En
      eks=En
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
!      xb=min(1d0,bkp(3)/ybo)
!      xb=max(-1d0,xb)
!      xb=acos(xb)
!      xbd=xb/degrad
      yb=abs(ybo-Kthr)
!      if(abs(bkp(1)).lt.1.d-6.and.abs(bkp(2)).lt.1.d-6) then
!      phKb=0.d0
!      else 
!      phKb=atan2(bkp(2),bkp(1))
!      endif
      xbd=tcd
      xb=tc
      phKb=0d0
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

      maxtab=-10d0
      imax=0
      jmax=0
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
      if (abs(sumps).gt. maxtab) then
        maxtab=abs(sumps)
        imax=ith+1
        jmax=ier
      endif
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
          write(99,'("fxyc=",20g12.5)') (fxyc(1:,ier),ier=1,6)
          write(99,*)' -> tmat=',tmat(iset,inc,iii)
      endif
      
      if (abs(tmat(iset,inc,iii)).gt.maxtab) then
       tmat(iset,inc,iii)=fxyc(imax,jmax)
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
      ampt2 =0. !!!! TESTING
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
      tmatsq=tmatsq+abs(ampt)**2*sqrt(polar(imo))
c2      tmatsq2=tmatsq2+abs(ampt2)**2*polar(imo)   !!! TESTING
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

      mult=mupt**2*mucv*m(kp)*m(bkp)/m(qp)/(2d0*pi*hc)**5/hc
      sigphi(ip)=sigphi(ip)+10.d0*mult*tmatsq
52    continue !sigphi(ip)=0.d0

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

