c *** ---------------------------------------------------
c     Build & solve CDCC equations
      subroutine cdcc(iexgs,ncc)
c *** ---------------------------------------------------
      use channels, only:jpiset,jpsets,jptset
      use xcdcc   , only: nex,jpch,exch,parch,iftrans,
     &            lamax,hcm,nmatch,elab,ecm,smats,
     &            jtmin,jtmax,rvcc,rmaxcc,nrcc,method,
     &            jump,jbord
      use nmrv,only: hort,rmort
      use sistema
      use globals, only: kin,debug,verb
      use constants, only: hc,amu,pi,e2
      use memory   , only: tcc
!#ifdef _OPENMP
!      use omp_lib
!#endif
      implicit none
c     ----------------------------------------------------------------
      logical:: dry,skip
      character*5 jpi 
      CHARACTER PARITY(3)
      character(len=80) numethod
      DATA PARITY / '-','?','+' / 
      integer:: i,iex,iexgs,ir,nfacmax,ncc
      integer:: l,lmin,lmax,icc,pargs,ijp,parp,partot
      integer:: ich,nch,icalc,nsets,ijt,njt,nchmax,inc,ijump,iblock
      integer:: nlag,ns ! R-matrix
c     ----------------------------------------------------------------
      real*8 :: jtot,jp,jtotmax,jblock
c *** Variables for CC (move to a module?) 
      real*8 :: jpgs,ex,kcmi,kcmf,mupt,ecmi,ecmf,etai,conv
      real*8 :: rturn
      real*8 :: xsr,xsinel,factor,xsrj,xsinelj
      real*8 start, end
c     OPENMP----------------------------------------------------------------
      integer:: num_threads,OMP_GET_NUM_THREADS,omp_get_thread_num
      logical:: omp_in_parallel
c     ---------------------------------------------------------------

c     ---------------------------------------------------------------
      complex*16 smat
c     ---------------------------------------------------------------
      namelist /reaction/ elab,jtmin,jtmax,jump,jbord,mp,mt,zt,skip
      namelist /numerov/ method,
     &         hcm,rmaxcc,rcc,hort,rmort,
     &         nlag,ns ! R-matrix
c     ---------------------------------------------------------------


c ... initialize variables --------------------------------------------------
      dry=.true.
      skip=.false.
      jtmin=-1; jtmax=-1
      jump=0; jbord=0; 
c ----------------------------------------------------------------------------
      write(*,'(//,5x,"******  REACTION CALCULATION ******")')
      read(kin,nml=reaction)
      if ((jtmin.lt.0).or.(skip)) then
        write(*,'(5x,"[No reaction calculation requested]",/)')
        return
      endif

      ecmi=elab*mt/(mp+mt)
      mupt=mp*mt/(mp+mt)*amu
      kcmi=sqrt(2*mupt*ecmi)/hc
      conv=(2*mupt/hc**2)
      etai=conv*zp*zt*e2/kcmi/2.
      zp  =zc+zv 
      if ((zp.lt.0).or.(zt.lt.0)) then 
          write(*,*)'Projectile or target charge not specified!'
          stop
      endif

      write(*,340) mp,zp,mt,zt !,mupt/amu
340   format(/,5x,"Projectile: Mass=",1f8.3," amu   Z=",1f6.1,3x,
     & /      ,5x,"Target:     Mass=",1f8.3," amu   Z=",1f6.1,/)
!     & /,"Reduced mass"," mu=",f8.3, ")

      write(*,'(5x,"J interval [min,max]=",2f6.1)') jtmin,jtmax
      if (maxval(jbord).gt.0) then
!        write(*,*)' J jumps defined:'
        do ijt=1,5
        if (jbord(ijt).eq.0) cycle
        write(*,345) jbord(ijt),jump(ijt)
345     format(8x,'o from J=',1f6.1,' in steps of', 1f6.1)
        enddo
      endif


c ... Numerov specifications
      method=0 
      hort  =0
      hcm   =0
      rmort =0
      nlag  =0; ns=0
      read(kin,nml=numerov)
      select case(method)
      case(0) 
       write(numethod,'(a60)') "Predictor-corrector Numerov"
      case(1)
       write(numethod,'(a60)')"ENA with 5 terms in Cosh[Sqrt[T]]"
      case(2)
       write(numethod,'(a60)')
     &      "ENA; 5 terms in Cosh[Sqrt[T]], only diag."
      case(3)
       write(numethod,'(a60)')"Raynal (ECIS)"
      case(4)
       write(numethod,'(a60)')"Modified ENA (fresco)"
      case(5)
       write(numethod,'(a60)')"R-matrix solution"
      end select 
      write(*,350) adjustl(numethod),hcm,rmaxcc
      if (abs(hort).gt.0) then
      write(*,'(7x,"[solutions orthogonalized every",f5.1,
     & " fm by QR algorithm]")') hort
350   format(5x, "-Method of solution: ",a60
     &     ,/,5x,"-Radial grid: Step",1f6.3,
     &     " fm and matching at",f8.2," fm")
      endif
c     ------------------------------------------------------------

      if (hcm.gt.1e-6) then 
        nrcc=nint(rmaxcc/hcm)+1
      else
       write(*,*) 'hcm too small!'; stop
      endif
      allocate(rvcc(nrcc))
      do ir=1,nrcc
        rvcc(ir)=dble(ir-1)*hcm
      enddo

      if (.not.iftrans) then !states & coupling pots read externally
          call read_states()
      endif

c ... BUILD CC SETS 
      nchmax=0
      jtotmax=0
      jpgs=jpch(iexgs)               ! angular momentum of proj gs
      pargs=parch(iexgs)             ! parity of proj gs
      lamax=nint(2*maxval(jpch))     ! Maximum multipole allowed Jp values
      write(*,'(/,3x,"=> there are",i5," states")') nex
      write(*,'(4x,"[Gs is state #",i2, 
     &  " and has spin/parity=",f3.1,a1," ]")')
     &  iexgs,jpgs, parity(pargs+2) !jpgs,pargs
      write(*,*)
!      lmax=nint(jtmax+maxval(jch))     
      if ((-1)**(2*jtmin)*(-1)**(jpgs).lt.0) then
         if (jtmin.eq.0.) then 
            jtmin=jpgs
         else
            jtmin=jtmin-0.5
         endif
      endif
      if (jtmin.lt.0) jtmin=jpgs
      njt  =jtmax-jtmin+1
      nsets=0
      if (debug) write(99,*)'Allocating memory for',2*njt,'CC sets'
      allocate(jptset(njt*2))
      icc=0
      jptset(:)%interp=.false.
      ijump=1
      do ijt=1,njt
      jtot=jtmin+dble(ijt-1)
      
      if (jtot.gt.jtotmax) jtotmax=jtot
      do i=1,5
        if ((jtot.eq.jbord(i)).and.(jump(i).gt.0)) then
        ijump=jump(i)
        jblock=jtot
        iblock=ijt
        write(*,*)'Start J-block at J=',Jtot,' J step=',ijump
        exit
        endif 
      enddo
      do partot=1,-1,-2
      nch=0     ! channel index within CC set icc
      icc=icc+1 ! CC set (total J/pi)
      if (verb.ge.3)
     &    write(*,400) icc,jtot,parity(partot+2)
400   format(/,4x,"CC set ",i3,' with J/pi=',1f5.1,a1)
      do iex=1,nex
      jp   =jpch(iex)
      parp =parch(iex)
      ex   =exch(iex)
      ecmf =ecmi + exch(iexgs)-ex
      if (ecmf> 0) then ! open channel    
        kcmf=sqrt(2*mupt*ecmf)/hc
      else 
        kcmf=sqrt(-2*mupt*ecmf)/hc
        write(*,*)'State',iex,' closed for Ecm=',Ecmi
      endif

      lmin=nint(abs(jtot-jp))
      lmax=nint(jtot+jp)
      do l=lmin,lmax
       if ((-1)**l*parp.ne.partot) cycle      
       nch=nch+1
!       if (nch.eq.1) icc=icc+1 
       if (verb.ge.3) write(*,500)  nch,l,jp,jtot,iex,ex,kcmf
500    format(3x," Ch.",i3, ": (",i3,f6.1,")",1f6.1,2x,
     &        "State: ",i3,5x,"Ex=",1f8.3," MeV;  kcm=",1f7.3)
       jptset(icc)%l(nch)  = l
       jptset(icc)%jp(nch) = jp
       jptset(icc)%idx(nch)= iex
       jptset(icc)%kcm(nch)= kcmf
       jptset(icc)%exc(nch)= ex - exch(iexgs)
      enddo !l
      enddo !iex
!      if (nch.gt.0) then
       jptset(icc)%nchan  =nch
       jptset(icc)%partot =partot
       jptset(icc)%jtot   =jtot
       if (nch.gt.nchmax) nchmax=nch
       if (mod(ijt-iblock,ijump).ne.0) jptset(icc)%interp=.true.
       if (mod(ijt-iblock,ijump).ne.0) write(*,*)'icc,Jp=',
     &     icc,Jtot,' interp'
       if (verb.gt.3) write(*,*)'[',nch,' channels]'
!      endif
      enddo !ipar
      enddo !icc

      dry=.false.
      ncc=icc
      write(*,'(/,8x,"=> Number of CC sets=",i4)') ncc  
      write(*,'(8x,  "=> Max number of chans =",i3)') nchmax  
      write(*,'(8x,  "=> Max L =",i4)') maxval(jptset(icc)%l(:))  
      write(*,'(8x,  "=> Max Jp =",f4.1)') maxval(jptset(icc)%jp(:))  
      write(*,'(8x,  "=> Max JTOT =",f5.1)') jtotmax 
      
!!!! TEMPORARY SOLUTION
      if (ncc.gt.2000) then
        write(*,*)'too many J/pi sets!.Increase dim of jptset'
        stop
      endif

      call flush(6)

c --------------------------------------------------------------------------
!!! CHECK FACTORIALS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!      nfacmax=2*nint(jtotmax+maxval(jptset(icc)%jp(:)))
      nfacmax=1000
      if (debug) then
      write(*,*)'Allocating factorials for n=',nfacmax
      endif

      call factorialgen(nfacmax)
c --------------------------------------------------------------------------


c ----------------------------------------------------------------------------
c Radial formfactors (coupling potentials) from this run or from external file
c (in both cases, they are interpolated to the specified radial grid for CC )
      if (iftrans) then
        call int_ff()  ! interpolate F(R)'s calculated in this run
      else
        call read_ff() ! read F(R)'s radial formfactors from ext file
      endif
      call flush(6)
c ---------------------------------------------------------------------------


c ***
c *** Solve CC eqns for ach J/pi set ***
c ***
      allocate(smats(ncc,nchmax,nchmax)) 
      call cpu_time(start)
!!#ifdef _OPENMP
!$      num_threads = 4
!$      call OMP_SET_NUM_THREADS(num_threads)
!!#endif

C$OMP PARALLEL
!$      write(*,*)'In parallel region (T)? ', OMP_IN_PARALLEL()
!$      write(*,*)'Threads allocated : ', OMP_GET_NUM_THREADS()
C$OMP END PARALLEL

!     do icc=1,ncc
      icc=0
C$OMP PARALLEL PRIVATE(nch,jtot,l,rturn)

C$OMP  DO
      do ijt=1,njt
      xsrj   =0      ! reaction  x-section
      xsinelj=0      ! inelastic x-section
      do partot=1,-1,-2
!$OMP CRITICAL
!        thread_id = omp_get_thread_num()
        jtot=jtmin+dble(ijt-1)
        icc=icc+1    
        nch=jptset(icc)%nchan   ! number of channels for this J/pi set
!$        print*,'thread:',omp_get_thread_num(),
!$     &  'icc,ncc,jtot,jptset()=',icc,ncc,jtot,jptset(icc)%jtot
!$OMP END CRITICAL      
        if (jtot.ne.jptset(icc)%jtot) then
          print*,'internal error: solvecc'
        endif
        if (partot.ne.jptset(icc)%partot) then
           print*,'internal error: solvecc'
        endif
!        jtot  =jptset(icc)%jtot 
!        partot=jptset(icc)%partot
        xsr   =0      ! reaction  x-section
        xsinel=0      ! inelastic x-section

        write(*,300) icc,jtot,parity(partot+2)
!300     format(/,2x," ***** CC SET",i4, " J/pi=",1f6.1,a1, "*****")
300     format(/,5x, 35("*")," CC SET",i4, " J/pi=",1f6.1,a1,2x,35("*")) 
        if (nch.eq.0) then 
          write(*,*) '  [no channels for this J/pi]'; cycle
        endif
        if (jptset(icc)%interp) then
          write(*,*) '  [ skipping => interpolated ]';cycle
!         goto 330
        endif 

c 1270 FORMAT( /,1X,118('#'),/,' #',116X,'#',/
c     X           ,' #',' Total SPIN and PARITY =',F8.1,1X,A1,',',I5,
c     X' channels,',I5,' in 1st block.  Rmin & Coul turning =',
c     X   2F7.1,' fm.  #'/ ' #',116X,'#',/1X,118('#')/)
        l=jptset(icc)%l(iexgs)
        rturn =(etai + SQRT(etai**2 + L*(L+1d0)))/kcmi
        write(*,310) nch,rturn
310     format(5x,"[",i3," channels",4x,"Elastic R-turn=",1f6.2," fm ]")
   
c ... construct < c| V | c'> matrix for this J/pi set from radial F(R)
        call makevcoup(icc,nch)  ! 


c ... Solve CC equations for each incoming channel 
        do inc=1,nch  ! incoming channels are those with iex=iexgs
          iex  = jptset(icc)%idx(inc)
          if (iex.ne.iexgs) cycle 
          call solvecc (icc,nch,inc,nrcc,nlag,ns)    
          
c ... Compute integrated cross sections
        do ich=1,nch
          kcmf  = jptset(icc)%kcm(ich)
          smat=smats(icc,inc,ich)*sqrt(kcmf/kcmi)  ! velocity factor
          factor=10*(pi/kcmi**2)*(2.*jtot+1)/(2*jpgs+1)
          if (jptset(icc)%idx(ich).eq.iexgs) then
            xsr=xsr+factor*(1-abs(smat)**2)
            xsrj=xsrj+factor*(1-abs(smat)**2)
          else 
            xsinel=xsinel  + factor*abs(smat)**2
            xsinelj=xsinelj+ factor*abs(smat)**2
          endif
        enddo ! ich
       enddo !inc
      
       write(*,320)xsr,xsinel,xsr-xsinel
320    format(5x,"o X-sections: Reac=", 1f8.3,5x,
     &   " Inel+BU=",1f8.3, 5x," Abs=",1f8.3," mb")

       enddo ! par
c       write(157,'(1x,1f6.1,a1,2x,3f12.6)') jtot, parity(partot+2),
c     &    xsrj-xsinelj,xsrj,xsinelj
       if (jptset(icc)%interp) cycle 
       write(157,'(1x,1f6.1,2x,3g14.5)') jtot, 
     &    xsrj-xsinelj,xsrj,xsinelj

       call flush(157)
       call flush(6) 
330    enddo !ijt
!      enddo ! ncc  -------------------------------------------------
C$OMP  END DO
C$OMP END PARALLEL


       call cpu_time(end)
       write(*,'(/,20x, "** All CC sets done in",1f12.2,"secs ** ")')
     & end-start
       tcc=end-start
!      call xsecs(kin,ncc,iexgs)
      end subroutine




c *** -------------------------------------------------------
c *** Solve CC eqns using Numerov method
c c *** -----------------------------------------------------
      subroutine solvecc(icc,nch,inc,nr,nlag,ns)
      use xcdcc   ,  only: hcm,elab,smats,rvcc,method
      use channels,  only: jptset
      use nmrv,      only: vcoup, ech,ql,hort,cutr
      use constants, only: hc,amu
      use sistema
      use globals  , only: verb,debug
!      use memory   , only: tcc
      implicit none
      integer    :: ir,n,nr,icc,nch,inc,nlag,ns
      complex*16 :: phase(nch),wf(nch,nr),smat(nch)
      real*8:: ecm
      real*8:: rmass,factor,rstart
!!!!!!!!!! TEST: delete me when I work 
      real*8 vcoul,rc
      logical:: test=.false., info=.true.
      real*8 ::v0,r0,a0,w0,ri,ai,a13,vaux,waux,ws,r
      real*8 ::ti,tf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      write(*,'(//,3x," Numerov integration:")') 
      debug=.false.
      rmass=mp*mt/(mp+mt)
      ecm=elab*mt/(mp+mt)
      factor=(2*amu/hc**2)*rmass


!      allocate(wf(nch,nmatch))
      if (allocated(ech))   deallocate(ech)
      if (allocated(ql))    deallocate(ql)
      allocate(ech(nch))
      allocate(ql(nch))

      ql(1:nch) =jptset(icc)%l(1:nch)  
      ech(1:nch)=jptset(icc)%exc(1:nch)


!!!!! TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (test) then
      write(*,*)'TESTING NUMEROV INTEGRATION'
c *** Channels and potentials
      ql(1:nch) =2
      ech(1:nch)=0.0
      inc=1
      zp=1; zt=6
      ecm=18.46;
      rc=1.25

      v0=-50; r0=1.2; a0=0.5
      w0=-10; ri=1.2; ai=0.5
c we do not consider the 1st point (r=rmin)
c becasuse in scatcc the 1st point is r=h
      rstart=hcm
      do ir=1,nr
        a13=12**.333333       
        r=ir*hcm 
        vaux=ws(r,v0,r0*a13,a0)
        waux=ws(r,w0,ri*a13,ai)
        vcoup(1,1,ir)=cmplx(vaux,waux) + vcoul(r,zp,zt,rc*a13) 
!        write(40,'(1f8.3,2x,50f12.8)') r,vcoup(1,1,ir)
      enddo      
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c Solve eqns
c-------------------------------------------------------
      if (debug) then
        write(*,*) 'Calling Numerov with:'
        write(*,*)'    Incoming channel=',inc
        write(*,*)'    (2*amu/hc**2)*rmass=',factor
        write(*,*)'    Ecm=',ecm
        write(*,*)'    l-values=',ql(1:nch)
        write(*,*)'    Channel energies:',ech(1:nch)
        write(*,*)'    Z1*Z2=',zp*zt
      if (.not.allocated(vcoup)) then
        write(*,*)'solvecc: vcoup not allocated!'
        stop
      endif
      endif
c -----------------------------------------------------
      call cpu_time(ti)
      rstart=rvcc(1)
      cutr=-10
      select case(method)
      case(0) ! predictor-corrector (Baylis & Peels)
        call schcc(nch,ecm,zp*zt,inc,ql,factor,hcm,
     &  rstart,nr,wf,phase,smat,info)

      case(1,2,3) ! Enhanced Numerov (Thorlacious & Cooper) / Raynal
        call schcc_ena(nch,ecm,zp*zt,inc,ql,factor,hcm,
     &  rstart,nr,wf,phase,smat,method,info)

      case(4)     ! Enhanced Numerov (Fresco version of T&C)
        call schcc_erwin(nch,ecm,zp*zt,inc,ql,factor,hcm,
     &  rstart,nr,wf,phase,smat,method,info)

      case(5)     ! R-matrix method (P. Desc. subroutine)
        call schcc_rmat(nch,ecm,zp*zt,inc,ql,factor,hcm,
     &  rstart,nr,wf,phase,smat,info,nlag,ns)

      case default
        write(*,*)'Method',method,' not valid!'
        stop
      end select

      call cpu_time(tf)
      write(*,'(5x,"[ CC solved in",1f8.3," sec ]")') tf-ti
!      tcc=tcc+finish-start
      smats(icc,inc,1:nch)=smat(1:nch)

      
!!!!! TEST!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (icc.eq.1) then 
      do ir=1,nr
      write(50,'(1f8.3,2x,50f12.8)') rvcc(ir),(wf(n,ir),n=1,min(2,nch))
      enddo
      write(50,*)'&'
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     Deallocate variables
      if (allocated(vcoup)) deallocate(vcoup)
      deallocate(ql)
      end subroutine



c *** --------------------------------------------------------------
c **  Scattering amplitudes & cross sections
c *** --------------------------------------------------------------
      subroutine xsecs(kin,ncc,iexgs)
      use xcdcc,    only: smats,elab,jpch,jtmin,jtmax,nex,exch,parch,
     &                    famps
      use channels, only: jptset,jpiset,jpsets
      use factorials
      use sistema
      use constants
      use parameters, only:  maxchan,maxsets
      use globals,  only: debug,written,verb
      use nmrv,     only: hort
      use wfs,      only: nr,dr,energ,rvec,wfr
      implicit none
      logical :: doublexs,jsets(maxsets)
c     ---------------------------------------------------------
      integer, parameter :: klog=99
      integer   , parameter:: kamp=136, kxs=300,kfam=137,
     &                        ksmat=170, !ksj=156,
     &                        nearfa=1
      integer :: istat,kin
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
      real*8  :: thmin,thmax,dth,cth,sth,s2,th,thrad
      real*8, allocatable :: pl(:,:),delci(:),delcf(:)
      real*8  :: rm,rmp,rmlp,lri,lrf
      real*8  :: factor, r1,r2,delif
      real*8  :: xsruth,xs,xstot,xsj,xr,xrj,sigex(nex),xsaux
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

      namelist/xsections/ thmin,thmax,dth,fileamp,doublexs,jsets,
     &                    ermin,ermax,ner,icore

c initialize -------------------------------------------------
      written(kamp)=.true.
      written(kxs)=.true.
      written(kfam)=.true.
!      written(ksj)=.false.
      written(kxs+1:kxs+min(nex,9))=.true.
      pi=acos(-1d0)
      thmin=0; thmax=0; dth=0; doublexs=.false.
      ermin=-1; ermax=-1
      jsets(:)=.true.        
      fileamp=""
      ner=0
c     -----------------------------------------------------------
      rewind(kin); 
      read(kin,nml=xsections)

      if (thmax.gt.thmin) then
        nth=nint((thmax-thmin)/dth)+1
      else
        write(*,*)'thmax<thmin!; Aborting'; stop
      endif

      if (ner.eq.0) then
         write(*,*)'NER=0!!';
      endif

!      write(*,*)'ncc=',ncc

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
      nmi    =2*jpgs+1
      nmfmax =2*maxval(jpiset(:)%jtot)+1
      allocate(famps(nex,nth,nmi*nmfmax))
      
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

!      write(kfam,301) iex,exc,jpf,parity(parch(iex)+2)
!      write(kfam,*) elab,nth
      write(kfam,'(4f6.1,i5,i2,1f10.3)') jpgs,jtarg,jpf,
     & jtarg,nth,nearfa,elab
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
      allocate(ampaux(0:lfmax,nmi,nmf),lfv(0:lfmax) )
      ampl(:,:,:)=0; 
      ampaux(:,:,:)=0; lfv(:)=0

      do icc=1,ncc
       partot=jptset(icc)%partot
       jtot  =jptset(icc)%jtot
       nch   =jptset(icc)%nchan
       if (jptset(icc)%interp) then
!          write(*,*)'skipping jtot=',jtot
          cycle 
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
       
       if (lf.gt.lfmax) stop 'internal error; lf >lfmax in famps!' 
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
c Satchler (4.58) for spin zero target:
        ampl(lf,im,imp)= ampl(lf,im,imp)  
     &                 + phci*phcf*r1*r2*c1*sqrt(2*lri+1)  
     &                 * (smat-kron(ni,nf))*YLMC2(lf,mlp) 
        ampaux(lf,im,imp)= ampaux(lf,im,imp)  
     &                   + phci*exp(-zz*(delcf(lf)-delcf(0)))
     &                   * r1*r2*c1*sqrt(2*lri+1)           
     &                   * (smat-kron(ni,nf))*YLMC2(lf,mlp) 
        lfv(lf)=lf

!        write(777,*)'lf,mlp,ylm=', lf,mlp,YLMC(lf,mlp),YLMC2(lf,mlp)

        if (debug) then
        write(klog,'(6x,"Contrib to AMP: ni,li,iex=",
     &  3i3," => nf,lf,iex=",3i3,6x,"S=",2f12.6," r1,r2=",2f12.5)')
     &  ni,li,1,nf,lf,iex,smat,r1,r2
       endif
       enddo ! imp
       enddo ! imp
       enddo ! nchf
       enddo ! nchi
       enddo ! icc (j/pi sets)
    

!!!!!!!!!!!!!!!!!!!!!!!!! INTERPOLATION !!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(gaux(lfmax+1),xv(lfmax+1))
      do icc=1,ncc
       partot=jptset(icc)%partot
       jtot  =jptset(icc)%jtot
       nch   =jptset(icc)%nchan
       if (nch.eq.0) cycle
       if (.not.jptset(icc)%interp) cycle

       do ni=1,nch
       if (jptset(icc)%idx(ni).ne.iexgs) cycle ! no inc waves 
       li   =jptset(icc)%l(ni)
       lri  =li
       jpi  =jptset(icc)%jp(ni)             ! should be just jpgs
       nmi  =2*jpi+1

       do nf=1,nch
       if (jptset(icc)%nchan.eq.0) cycle
       if (jptset(icc)%idx(nf).ne.iex) cycle 
       lf  =jptset(icc)%l(nf)
       lrf =lf
       jpf =jptset(icc)%jp(nf)
       nmf =2*jpf+1
       
       if (lf.gt.lfmax) stop 'internal error; lf >lfmax in famps!' 

       do im=1,nmi
        rm= -jpi + (im -1)
        do imp=1,nmf
        rmp=-jpf + (imp-1)
        mlp=nint(rm-rmp)
        rmlp=rm-rmp
        if (abs(mlp).gt.lf) cycle        
        phci =exp(zz*(delci(li)-delci(0)))
        phcf =exp(zz*(delcf(lf)-delcf(0)))
        gaux(1:lfmax+1)=ampaux(0:lfmax,im,imp)
        xv(1:lfmax+1)  =lfv(0:lfmax)
        write(*,'(10(i3,2f10.5))') (nint(xv(l)),gaux(l),l=1,lfmax)
        caux=cfival(lrf,xv,gaux,lfmax+1,alpha)
        write(*,*)'lf,caux=',lf,caux
        ampl(lf,im,imp)= ampl(lf,im,imp)+caux*(phcf*phcf)

       enddo ! imp
       enddo ! imp
       enddo ! nchf
       enddo ! nchi
       enddo ! icc (j/pi sets)
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c write angle-indep amplitudes ----------------------------------------
       write(kamp,300) iex,exc,jpf,parity(parch(iex)+2)
       do lf=0,lfmax
!225    WRITE(kamp,226) LF,((AMPL(lf,im,imp)*phcf,imp=1,nmf), im=1,nmi)
225    WRITE(kamp,226) lf,((ampl(lf,im,imp),imp=1,nmf), im=1,nmi)
226    FORMAT(' LP=',I5,1X,1P,10E12.4 / (10X,10E12.4),/)
       do im=1,nmi
        rm= -jpi + (im -1)
        do imp=1,nmf
        rmp=-jpf + (imp-1)
        enddo
       enddo
       enddo ! lf


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
       famps(iex,ith,iam)= fam(ith,im,imp)
       enddo !im
       enddo !imp
!       write(kfam,'(2x,"theta=",1f6.2)') th 
       write(kfam,'(2x,1f6.2)') th ! changed in v2.2e
       write(kfam,228) ((fam(ith,im,imp),imp=1,nmf),im=1,nmi)
228    format(1P,6E12.4)
!       write(kfam,*) ' '
       enddo !angles
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
        deallocate(ampl,fam,pl,delcf,ampaux,lfv,xv)
        if(allocated(gaux)) deallocate(gaux)
      enddo ! iex 


c *** Angle-integrated cross sections 
      sigex(:)=0
      xstot = 0.
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
       

       if (iex.eq.iexgs) then ! elastic channel 
         written(ksmat)=.true.
         smat=smats(icc,ni,nf)
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
790     format(8x,"#",i4,3x,1f5.1, a1,f10.3,"->",f12.6,' mb')
795     format(1x,1f10.3,2x,1f12.6,1x,i3,1f5.1,2x,i2)
      enddo

c Generate an approximate distribution dsigma/dE
      call budist(exch,nex,sigex)
      
c Double x-sections 
1000  if (doublexs) call d2sigma(nth,dth,thmin,ermin,ermax,
     &                           ner,icore,jsets,fileamp)
      if (allocated(famps)) deallocate(famps)
      end subroutine





c ----------------------------------------------------------------
c Double differential x-sections dsigma/dEx dW  (NEW VERSION) 
c ----------------------------------------------------------------
      subroutine d2sigma(nth,dth,thmin,emin,emax,ncont,
     &                   icore,jsets,fileamp)
      use xcdcc,    only: elab,jpch,nex,exch,parch,famps
      use channels, only: jpiset,jpsets,sn
      use sistema
      use constants
      use parameters, only:  maxchan,maxsets
      use globals,  only: debug,written,verb,kin
      use wfs,      only: nr,dr,energ,rvec,wfr
      implicit none
c     ----------------------------------------------------------------------- 
      character*40:: line,fileamp,sys*3
      logical :: jsets(maxsets),energy,zert
      integer , parameter:: kfam=137  ! file containing amplitudes
!      integer , parameter:: ncont=300 ! c-v relative energies
      real*8 , parameter:: eps=1e-3
      integer , parameter:: maxstat=300,ncleb=500
      integer :: ic,icore,ncont
      integer ::nsetmax,ir,iset,n,nchan,iam,ich,iecv,nth,iex,ith,mth
      integer ::nmi,nmf,nmfmax,nm,li,lf,im,imp,inc,nearfa
      integer,allocatable::ips(:,:)
c     -------------------------------------------------------------------------
      integer:: imu,isig,imo,inu,iang,iii,niag,nial,nchann
      integer:: ii,itc,iv,itv,itphim,itph,iflag,ien,iten,idet
      integer:: ip,iroots,iroot,iroot1n,iroot2n,j
      integer:: ilval,itphi,lmax,ind
c     -------------------------------------------------------------------------
      real*8:: kcv,dkrdk,Kthr,atarget
      real*8:: ji,jf,jci,jcf,raux,jac,jpi,jpf,jtarg
      real*8:: ecmi,ecmf,kcmi,kcmf,excore,th,thmin,dth
      real*8:: dec,ebind,ethr,emin,emax
      real*8:: facK,mv,mc,mucv,mupt,f2t,xstot,sumr
      real*8:: tmatsq,sc,mu,sig,mo,rnu,cleb,cgc(ncleb)
c     -------------------------------------------------------------------------
      real*8:: dot,m,plmrec,ylmc2
      real*8:: dphir,dtc,dtv,phiu,tcl,tvl,tcu,tvu,dtcr,dtvr
      real*8:: acc,degrad,mtot,totp
      real*8:: erel,erelmax,tcd,tc,costc,sintc,tvd,tv,costv,sintv
      real*8:: En,Enlow,Enup,dEn,p1L,p2L,phid,phil,dphi,phi,cospv,sinpv
      real*8:: cospc,sinpc,aq,bq,cq,dq,Ec,Ev,eKb,eks
      real*8:: qp(3),kcL(3),kvL(3),pcL(3),pvL(3),tem(3),cplm(0:10,21)
      real*8:: ptot(3),ktot(3),bkp(3),kp(3),sigphi(4000),xsig
      real*8:: mkp,co,phKb,phks,co2,si2,si,aleg,mult,xb,xbd,ybo,yb
c     -------------------------------------------------------------------------      
      real*8,allocatable:: angsr(:),xyt(:,:),ecv(:)
c     -------------------------------------------------------------------------
      complex*16,allocatable:: wfcont(:,:,:),fxyc(:,:),
     &                         ampaux(:,:,:),gaux(:),tmat(:,:,:)
      complex*16,allocatable,target:: gsolap(:,:,:,:)
      complex*16 :: resc,caux,faux,cfival
      complex*16,pointer:: overlap(:)
c     --------------------------------------------------------------------------
      complex*16:: ampt,xi,f2c,sumps
      complex*16:: phask1,phask,phas,ylm(0:10,21)
c ------------------------------------------------------------------------------
      namelist /framework/ sys,idet,atarget
      namelist /gridener/ Enlow,Enup,dEn
      namelist /gridthetac/ tcl,tcu,dtc
      namelist /gridthetav/ tvl,tvu,dtv
      namelist /gridphi/ phil,phiu,dphi
      nsetmax=maxval(jpiset(:)%nex) !max number of states within j/pi set
      energy=.true.
      mc=ac*amu
      mv=av*amu
      mt=mt*amu
      mp=mp*amu  !mc+mv
      mupt=mp*mt/(mp+mt)
      mucv=mc*mv/(mc+mv)
      nmi=nint(2*jpiset(1)%jtot+1)
      nmfmax =2*nint(maxval(jpiset(:)%jtot))+1

      write(*,'(//,3x, "** DOUBLE DIFFERENTIAL CROSS SECTIONS **",/ )')

!      write(*,*)'j/pi, nmi=',jpiset(1)%jtot,nmi
!      write(*,*)'ncont=',ncont
      write(*,*)'emin,max=',emin,emax

c     Read external f(theta) amplitudes
      if (.not.allocated(famps).and.(fileamp.ne."")) then
      write(*,'(3x," Reading f-amps from file ",a20)')fileamp
      open(kfam,file=fileamp,status='old')
!      rewind(kfam)
      read(kfam,*,err=900) jpi,jtarg,jpf,jtarg,mth,nearfa,elab
!      write(*,*) 'mth,nerfa,elab=',mth,nearfa,elab
      rewind(kfam)

      allocate(famps(maxstat,mth,nmi*nmfmax))
      iex=0
      do iset=1,jpsets  
      nex=jpiset(iset)%nex
      jpf=jpiset(iset)%jtot
      nmf=nint(2.*jpf+1.)
      nm=nmi*nmf
      do n=1,nex
      iex=iex+1
      if (iex.gt.maxstat) stop'increase dimension of famps' ! temporary solution
!      write(*,*)'expected jpf=',jpf
      read(kfam,*) jpi,jtarg,jpf,jtarg,mth,nearfa,elab
!      write(*,*) 'ie,jpf,mth=',iex,jpf,mth

      if (iex.eq.1) write(*,*) ' amplitudes are for elab=',elab
      if(nth.ne.mth) then
         write(*,*) 'mismatch in number of angles: specified:',nth,
     &              ' found:',mth
      endif 
      do ith=1,nth
      read(kfam,*) th  
!      write(*,*) 'th=',th   
      read(kfam,228) (famps(iex,ith,iam),iam=1,nm) 
228   format(1P,6E12.4)
      enddo !states
      enddo !ith
      enddo !iset
650   write(*,'(5x,"[f amplitudes read for ",i3," states]")')iex
      endif
c     -------------------------------------------------------------------------

      ecmi=elab*mt/(mp+mt)
      kcmi=sqrt(2*mupt*ecmi)/hc
      ebind=energ(1,1) ! exch(1)
      Ethr=Ecmi+ebind
      write(*,*)'Ecmi=',Ecmi, "Kcm=",kcmi, 'Ethr=',Ethr


c *** Assign global index to states of all j/pi sets iPS=iPS(iset,n)
      iex=0
      allocate(ips(jpsets,nsetmax))
      do iset=1,jpsets
      nex=jpiset(iset)%nex
      do n=1,nex
      iex=iex+1
      ips(iset,n)=iex
      if (verb.ge.3) write(*,*)'Set:',iset,' state=',n,"-> iex=",iex
      enddo
      enddo


c *** Sum discrete inel/breakup angular distributions      
      open(90,file='dsdw_ps.xs',status='unknown')
      do ith=1,nth
      raux=0
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nex=jpiset(iset)%nex
      jpf=jpiset(iset)%jtot
      nmf=nint(2*jpf+1)
      iam=0
!      if (ith.eq.1) write(*,*)'iset,nex,jpf,nmi,nmf=',
!     & iset,nex,jpf,nmi,nmf
      do im=1,nmi
      do imp=1,nmf
      iam=iam+1  
      do n=1,nex    
      iex=ips(iset,n)
      if (energ(iset,n).lt.0)    cycle  ! bound state 
      if (energ(iset,n).gt.Ethr) cycle  ! closed channel
      raux=raux+ 10.*abs(famps(iex,ith,iam))**2/nmi
      enddo ! states within j/pi set
      enddo !imp 
      enddo !im
      enddo !j/pi sets
      write(90,*) thmin + dth*(ith-1),raux
      enddo !ith
      close(90)
c ------------------------------------------------------------------------------
      
      
c *** Overlaps between PS's and Scattering states  -----------------------------
      if (allocated(gaux)) deallocate(gaux)
      allocate(gaux(nr))
      allocate(gsolap(jpsets,nsetmax,maxchan,ncont),ecv(ncont))
      gsolap(:,:,:,:)=0
      if (emax<0) emax=maxval(energ(:,:))
      if (emax>Ethr) emax=Ethr-0.001
      if (emin.le.0.) emin=0.01
      dec=(emax-emin)/dble(ncont-1)
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      nex =jpiset(iset)%nex
      do inc=1,nchan 
      if (inc.gt.maxchan) stop'inc too large in d2sigma!'
      jci    =jpiset(iset)%jc(inc)
      li     =jpiset(iset)%lsp(inc)
      excore =jpiset(iset)%exc(inc) ! core energy
      ic     =jpiset(iset)%cindex(inc)
      if (ic.ne.icore) then
         if (verb.ge.2)
     &   write(*,'(4x,"-> skipping core state",i3)')ic 
      cycle
      endif
!      if (emax.lt.excore) cycle !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST 2g
      if (allocated(wfcont)) deallocate(wfcont)
      allocate(wfcont(ncont,nchan,nr)); wfcont=0.
      write(*,*)'calling wfrange:emin,max=',emin,emax
      call wfrange_new(iset,nchan,inc,emin,emax,ncont,energy,wfcont)
      do n=1,nex
      if (energ(iset,n).lt.0)    cycle ! bound state
      if (energ(iset,n).gt.Ethr) cycle ! closed channel
      
      do iecv=1,ncont
      ecv(iecv)=emin+(iecv-1)*dec
      if (ecv(iecv).lt.eps) ecv(iecv)=eps
      kcv=sqrt(2.d0*mucv*ecv(iecv))/hc    
      sumr=0
      caux=0.
      do ich=1,nchan ! sum in final channels
      jcf =jpiset(iset)%jc(ich)
      lf  =jpiset(iset)%lsp(ich)
      jf  =jpiset(iset)%jc(ich)
      do ir=1,nr   
      gaux(ir)=wfr(iset,n,ich,ir)*wfcont(iecv,ich,ir)*rvec(ir)*
     &         (-1)**(li+jci+lf+jcf) !nfactl1+xfactI1+nfactl2+xfactI2)
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
     
c *** -------------- PRINT OVERLAPS FOR TESTING ----------------------------------
      if (verb.ge.3) then
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nex  =jpiset(iset)%nex
      nchan=jpiset(iset)%nchan
      
      do n=1,nex
      raux=0
      write(97,*)'# n -> Ex=',n,energ(iset,n)
      do iecv=1,ncont
       if (energ(iset,n).lt.0)    cycle ! bound state
       if (energ(iset,n).gt.Ethr) cycle ! closed channel
       kcv=sqrt(2.d0*mucv*ecv(iecv))/hc
       jac=mucv/(hc**2)/kcv
       write(97,111)ecv(iecv),
     &  (jac*abs(gsolap(iset,n,inc,iecv))**2,inc=1,nchan)
111    format(2x,1f12.6,2x,4g14.6)
       do inc=1,nchan
       raux=raux+jac*abs(gsolap(iset,n,inc,iecv))**2*dec
       enddo ! inc
      enddo !iecv
      write(97,*)'&'
      if (verb.ge.4) then
      write(*,*)'PS',n,' Norm solap=',raux*2/pi
      endif
      enddo ! n
      enddo !iset
      endif ! verb
*     -------------------------------------------------------------- 
*     compute and store CG coefficients
*     -------------------------------------------------------------- 
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
      do inu=1,2*li+1
      rnu=inu-li-1
      ind=ind+1
      cgc(ind)=cleb(li,rnu,sn,sig,ji,rnu+sig)*
     .cleb(ji,rnu+sig,jci,mu,jpf,rnu+sig+mu)
      enddo
      enddo
      enddo
      enddo
      enddo
      if(ind.gt.ncleb) then
      write(*,*),'increase ncleb up to',ind
      stop 
      endif
*     --------------------------------------------------------------- 
*     compute largest orbital angular momentum lmax
*     ---------------------------------------------------------------
      lmax=0
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      do inc=1,nchan 
      if(jpiset(iset)%lsp(inc).gt.lmax) then
      lmax=jpiset(iset)%lsp(inc)
      endif
      enddo
      enddo
*     ---------------------------------------------------------------
*     calculate and store coefficients for spherical harmonics
*     ---------------------------------------------------------------
      do ilval=0,lmax
      do inu=1,2*ilval+1
      cplm(ilval,inu)=ylmc2(ilval,inu-ilval-1)
      enddo
      enddo
*     --------------------------------------------------------------
      read(kin,nml=framework)
      write(77,'(a)') sys
      write(77,316) atarget,zt,elab/dble(nint(ac+av))
      write(77,316) dble(nint(ac)),zc
      write(77,316) dble(nint(av)),abs(ebind)
      write(77,316) elab
      write(77,*) idet
*     ---------------------------------------------------------------
      xi=(0.d0,1.d0)
      acc=1.d-6
      degrad=pi/180.d0
      mtot=mc+mv+mt
      Kthr=sqrt(2.d0*mupt*Ethr)/hc
      allocate(angsr(nth))
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
      dphi=(phiu-phil)/(itphi-1)
      dphir=dphi*degrad
*     ---------------------------------------------------------------
      write(77,307) itc,dtcr,tcl,tcu,dtc
      write(77,307) itv,dtvr,tvl,tvu,dtv
      write(77,307) iten,dEn,Enlow,Enup
      write(77,307) itphi,dphir,phil,phiu,dphi
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
*     ---------------------------------------------------------------
*     loop over core particle detection angle
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
*     -----------------------------------------------------------------
*     calculate second particle energy and momentum depending on idet
      do j=1,3
      if(idet.eq.1) then              
      tem(j)=ptot(j)-pcL(j)
      else
      tem(j)=ptot(j)-pvL(j)        
      endif
      enddo
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
      if(sys.eq.'lab') cq=cq-(elab+ebind-excore)
      if(sys.eq.'com') cq=cq-(ecmi+ebind-excore) 
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
      do j=1,3
      if(idet.eq.1) then
      kvL(j)=p2L*pvL(j)/hc                 
      kcL(j)=pcL(j)/hc                  
      else
      kcL(j)=p1L*pcL(j)/hc         
      kvL(j)=pvL(j)/hc              
      endif
      enddo
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
*     -----------------------------------------------------------------
*     increment maximum relative energy 
*     -----------------------------------------------------------------
      if(eks.gt.erelmax) erelmax=eks
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
*     -----------------------------------------------------------------
*     pick out nearest angles indices in array for interpolations
*     -----------------------------------------------------------------
      if (allocated(xyt)) deallocate(xyt)
      allocate(xyt(2,ncont))
      iang=nint((xbd-thmin)/dth)+1
      nial=max0(iang-2,1)
      niag=min0(nial+5,nth)
      nial=niag-5
      do iang=1,6
      xyt(1,iang)=angsr(nial+iang-1)
      enddo
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
      aleg=plmrec(ilval,inu-ilval-1,co,co2,si,si2)
      ylm(ilval,inu)=cplm(ilval,inu)*phas*aleg
      enddo
      enddo
      nchann=maxval(jpiset(:)%nchan)
      if (allocated(fxyc)) deallocate(fxyc)
      allocate(fxyc(10,ncont))
      if (allocated(tmat)) deallocate(tmat)
      allocate(tmat(jpsets,nchann,nmi*nmfmax))
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
      do inc=1,nchan
      do iii=1,iam
!      ic=jpiset(iset)%cindex(inc)
!      if (ic.ne.icore) cycle
!      excore =jpiset(iset)%exc(inc)
      do iecv=1,ncont 
      ecmf=ecmi-abs(ebind)-ecv(iecv)-excore
      if (ecmf.lt.0) then
      write(*,*) 'Ecm=',Ecmf,'Ecv=',Ecv(iecv)
      stop
      endif
      if (ecmf.lt.1e-4) ecmf=1e-4
      Kcmf=sqrt(2.d0*mupt*ecmf)/hc     
      f2t=-2.*pi*hc**2/mupt*sqrt(Kcmi/Kcmf) ! conversion factor f(theta) -> T-mat
      xyt(2,iecv)=Kcmf-Kthr
      do ith=0,5
      sumps=0.d0
      do n=1,nex
      if (energ(iset,n).lt.0)    cycle ! bound states
      if (energ(iset,n).gt.Ethr) cycle ! closed channels
      if(iii.eq.1.and.iecv.eq.1.and.ith.eq.0) then
      overlap=>gsolap(iset,n,inc,:)
      erel=hc**2*mkp**2/2.d0/mucv
      faux=cfival(erel,ecv,overlap,ncont,0.d0)
      endif
      iex=ips(iset,n)
      sumps=sumps+faux*f2t*famps(iex,nial+ith,iii)
      enddo ! PS's
      fxyc(ith+1,iecv)=sumps
      enddo ! theta
      enddo ! ecv
      tmat(iset,inc,iii)=f2c(xb,yb,xyt,fxyc,6,ncont,2,10,ncont)
      enddo ! channel
      enddo ! iii
      enddo ! jpset
*     ----------------------------------------------------------------- 
*     |t-matrix|**2 summed on spin projections / (2sp+1)
*     ----------------------------------------------------------------- 
      tmatsq=0.d0
      do imo=1,nmi
      mo=imo-jpiset(1)%jtot-1
      ind=0
      do imu=1,nint(2*sc+1)
      mu=imu-sc-1
      do isig=1,nint(2*sn+1)
      sig=isig-sn-1
      ampt=(0.d0,0.d0)
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      nex  =jpiset(iset)%nex
      jpf  =jpiset(iset)%jtot
      do inc=1,nchan  ! sum over open channels
!      ic  =jpiset(iset)%cindex(inc)
!      if (ic.ne.icore) cycle
      li     =jpiset(iset)%lsp(inc)
      ji     =jpiset(iset)%jsp(inc)
      jci    =jpiset(iset)%jc(inc)
      do inu=1,2*li+1
      rnu=inu-li-1
      ind=ind+1
      if(cgc(ind).lt.eps) cycle
      iam=nint((imo-1)*(2*jpf+1)+rnu+sig+mu+jpf+1)
      ampt=ampt+(-xi)**li*cgc(ind)*ylm(li,inu)*
     .tmat(iset,inc,iam)*exp(xi*nint(mo-rnu-sig-mu)*phKb)
      enddo
      enddo ! inc
      enddo ! j/pi
      tmatsq=tmatsq+abs(ampt)**2
      enddo
      enddo
      enddo
      tmatsq=(4.d0*pi/mkp)**2*tmatsq/dble(nmi)
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
      print 315,ien,iv,ii,xsig,En,iflag
      endif
*     -----------------------------------------------------------------
      write(77,315) ien,iv,ii,xsig,En
      call flush(77)
*     -----------------------------------------------------------------
*     close the angle (core and valence thetas) and energy loops
*     -----------------------------------------------------------------
30    continue      
20    continue
10    continue
      print 510
      print*,'  erelmax =',erelmax
      print 510
      return 
307   format(i5,f15.8,3f10.3)
315   format(3i5,d19.8,1x,f12.6,i5)
316   format(3f10.3)
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
*     l<=2 caluclated explicitly
*     l>2  calculated recursively
*-----------------------------------------------------------------------
      res=1.d0
      if(il.eq.0) return
      im=abs(imm)
*-----------------------------------------------------------------------
      if(il.eq.1) then
      if(im.eq.0) then
      res=c
      else if(im.eq.1)then
      res=s
      endif
      else if(il.eq.2) then
      if(im.eq.0) then
      res=(3.d0*c2-1.d0)/2.d0
      else if(im.eq.1) then
      res=3.d0*s*c
      else if(im.eq.2) then
      res=3.d0*s2
      endif
      elseif(im.gt.1) then
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
      if(enex(ia).lt.0) nb=nb+1
      enddo

      if (nset-nb.lt.3) return
      open(15,file='dsde_ps.xs') 
!      d1=2.d0*ENEX(2) - 1.5d0*ENEX(1) - 0.5d0*ENEX(3)
      d1=enex(nb+1)
      write(15,50) ENEX(nb+1)-be,sigset(nb+1)/d1,
     X	      enex(1)+d1/2.d0,d1 
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




c *** Calculate coupling matrix for CC set
      subroutine makevcoup(icc,nch)
      use xcdcc,     only: lamax,ffc,nrcc,nex,hcm,parch,rvcc
      use channels,  only: jptset
      use nmrv,      only: vcoup
      use sistema
      use globals,   only: verb, debug,nfacmax
      use memory
      implicit none
      logical   fail3
      integer   ir,nchans,nch,ni,nf,lam,icc,inc,ix
      real*8::  jtot,jpi,jpf,lri,lrf,lamr,ri,coef
      integer:: partot,pari,parf,idi,idf
      real*8::  r1,r2,r3,r13,x,zero,threej,sixj,zpt
      integer:: li,lf
      real*8 :: vcoul,frconv,ymem
      complex*16:: vaux,phc
      CHARACTER LNAME(0:14),PARITY(3)
      DATA PARITY / '-','?','+' / 
c     ...............................................................
      debug=.false.
       
!      write(*,*)'Coupling matrix for icc=',icc,' Nchan=',nch
!      write(*,'(5x,"zp,zt,rcc=",3f8.2)') zp,zt,rcc
      zero=0d0
      nchans =jptset(icc)%nchan   ! number of channels for this J/pi set
      partot =jptset(icc)%partot  ! parity of the set
      jtot   =jptset(icc)%jtot    ! total angular momentum of J/pi set
      if (nchans.ne.nch) then
       write(*,*)'Internal error in makevcoup: expected',nchans,
     &           'channels for CC set',icc,' but required',nch
       stop
      endif

      ymem=nch*nch*nrcc*lc16/1e6
      if (verb.ge.3) write(*,190) ymem
190   format(5x,"[ Coupling matrix require", 1f8.2," Mbytes ]")


      if (allocated(vcoup)) deallocate(vcoup)
      allocate(vcoup(nch,nch,1:nrcc))
      vcoup=0d0

      ix=0
      do ni=1,nch
      li =jptset(icc)%l(ni)
      lri=li
      jpi =jptset(icc)%jp(ni)
      idi =jptset(icc)%idx(ni)
      pari=parch(idi)
      do nf=ni,nch !Upper triangle available here only!!!
      lf=jptset(icc)%l(nf)
      lrf=lf
      jpf=jptset(icc)%jp(nf)
      idf =jptset(icc)%idx(nf)
      parf=parch(idf)
      r1=sqrt((2.*lri+1.)*(2.*lrf+1.)*(2.*jpi+1.)*(2.*jpf+1.))
     &  *(-1)**(jpi+jtot)

!      write(*,*)'makevcoup: jpi,jpf,jtot=',jpi,jpf,jtot,(-1)**(jpi+jtot)

      do lam=0,lamax
      lamr=lam
!      if (fail3(jpi,jpf,lamr)) cycle
      if (fail3(lri,lrf,lamr)) cycle
      if ((-1)**lam*pari.ne.parf) cycle

      r2=threej(lamr,lri,lrf,zero,zero,zero)
      r3=sixj(jpi,jpf,lamr,lrf,lri,jtot)
      r13=r1*r2*r3
      phc=(0d0,1d0)**(lrf-lri)
      ix=ix+1
      coef=(-1)**lamr*(2.*lamr+1.)*r13
      if (verb.ge.4) then     
      write(*,300) ix,ni,idi,li,jpi,parity(pari+2),
     &                nf,idf,lf,jpf,parity(parf+2),
     &             lam,r1,r2,r3,coef,ffc(idi,idf,lam,1)
300   format(5x,"IX=",i4," Ch. IEX L, Jpi=",3i3,f4.1,a1,
     &     " =>",3i3,f4.1,a1,3x,
     &     'with LAM=',i3,4x,'r1-3='3f8.3,
     &     2x,'COEF=',1f7.3,' F(1)=',2g10.4)
      endif
!      write(*,'("#r1,r2,r3,r =",5f8.3)') r1,r2,r3,r13
      do ir=1,nrcc      
!      ri=dble(ir-1)*hcm  ! CHECK!!!!!!!!!!!!
      ri=rvcc(ir)
      vaux= phc*coef*ffc(idi,idf,lam,ir)
      vcoup(ni,nf,ir)=vcoup(ni,nf,ir) + vaux


      if ((zp*zt.gt.1e-6).and.(ni.eq.nf).and.(li.eq.lf).
     &  and.(lam.eq.0)) then
        vcoup(ni,ni,ir)=vcoup(ni,ni,ir) + VCOUL(ri,zp,zt,Rcc)
      endif

      if (ni.ne.nf) vcoup(nf,ni,ir)=vcoup(ni,nf,ir)*conjg(phc)/phc
      
      if (debug) then 
      if (icc.eq.1) then
       if (ir.eq.1) write(97,'("#",2i3,2i4,i3,1g12.4)') 
     & ni,nf,idi,idf,lam,coef
        write(97,'(f8.3,10f12.5)')rvcc(ir),vcoup(ni,nf,ir),
     &   vcoup(nf,ni,ir)
        if (ir.eq.nrcc) write(97,*)'&'
      endif ! icc.eq.1 
      endif ! debug

      enddo !ir
      enddo !la
      enddo !m
      enddo !n


      if (.not.allocated(vcoup)) then
        write(*,*)'makevcoup: vcoup not allocated!'
        stop
      endif

      RETURN
      end subroutine 


      

c *** ---------------------------------------------------------------------
C     Use previously calcualted F(R)'s and interpolated in integration grid  
      subroutine int_ff()
c *** --------------------------------------------------------------------      
      use xcdcc,   only: hcm,ffc,lamax,nex,nrcc,jpch,Ff,
     &                   nrad3,rstep,parch,rmaxcc,rvcc
      use channels,only: jpiset,jpsets
      use wfs     ,only: idx
      use memory
      use globals ,only: verb
      implicit none
      integer :: ir
      integer :: m1,m2,nff,lam,par1,par2
      complex*16, pointer:: faux(:)
      real*8  :: rv(nrad3)
      real*8 :: jp1,jp2,r,ymem
      real*8, parameter:: alpha=0d0
      complex*16 cfival,caux
c     ----------------------------------------------------
     
      write(*,'(//,3x, "** INTERPOLATE FORMFACTORS **" )')

      lamax=2*maxval(jpiset(:)%jtot)
      write(*,*)'nrad3=',nrad3,' lamax=',lamax,' nrcc=',nrcc

      ymem=nex*nex*(lamax+1)*nrcc*lc16/1e6
      if (verb.ge.1) write(*,190) ymem
190   format(5x,"[ FF require", 1f8.2," Mbytes ]")

      if (allocated(ffc)) deallocate(ffc)
      allocate(ffc(nex,nex,0:lamax,1:nrcc))
      ffc=0
             
      do ir=1,nrad3
       rv(ir)= rstep + rstep*(ir-1) ! CHECK
      enddo     
      write(*,170) nrcc, hcm,rmaxcc,hcm
170   format(/,5x,"=>  Interpolating F(r) in grid with ",i4, ' points:',
     &           2x,"[ Rmin=",1f5.2,2x," Rmax=",
     &           1f6.1,1x," Step=",1f6.3, " fm ]",/)

      if (rmaxcc.gt.rv(nrad3)) then
        write(*,'(3x,"*** WARNING: Coupled eqns integrated up to",f6.1,
     &  1x,"fm, but formfactors calculated for R<",1f6.1," fm" )')
     &  rmaxcc,rv(nrad3) 
      endif 

      nff=0
      do m1=1,nex
      par1= parch(m1)
      jp1 = jpch(m1)
      do m2=m1,nex
      jp2 = jpch(m2)
      par2= parch(m2)
       do lam= nint(dabs(jp1-jp2)),nint(jp1+jp2)
       if (par1*par2*(-1)**lam<0) cycle
       nff=nff+1
       faux =>  Ff(m1,m2,lam,:)
c      interpolate & store in integration grid

       do ir=1,nrcc
!        r=hcm*ir ! grid used in scattcc
!        r=dble(ir-1)*hcm ! grid used in scattcc (Changed in v2.2)
        r=rvcc(ir)
        if (r.gt.rv(nrad3)) cycle ! DO NOT EXTRAPOLATE 
        caux=cfival(r,rv,faux,nrad3,alpha)
        ffc(m1,m2,lam,ir)=caux
       enddo ! ir
      enddo ! lam
      enddo !m2
      enddo !m1
  
      write(*,*) '=> there are', nff,' radial formfactors '
      nullify(faux)
      deallocate(ff)
      end subroutine



c *** ---------------------------------------------------------------
C     READ FORMFACTORS AND BUILD COUPLING MATRIX in integration grid  
      subroutine read_ff()
c *** ---------------------------------------------------------------      
      use xcdcc, only: hcm,nrcc,ffc,lamax,nex,nrcc,jpch,Ff,rmaxcc,rvcc
      use globals, only: debug,verb
      use memory
      implicit none
      complex*16, allocatable:: faux(:)
      real*8, allocatable:: rv(:)
      integer lam,m1,m2,npt,nff,ir
      real*8 factor,jt,ptr,ttr,rstep,rfirst,fscale
      real*8 r,x,y
      real*8 jpi,jpf,frconv,ymem
      real*8, parameter:: alpha=0d0
      integer,parameter :: kfr=4
      complex*16 cfival,caux
      character*40 comment
     
      write(*,'(//,3x, "** READ FORMFACTORS **" )')
      call flush(6)

      nff=0
      ymem=nex*nex*(lamax+1)*nrcc*lc16/1e6
      if (verb.ge.1) write(*,190) ymem
190   format(5x,"[ FF require", 1f8.2," Mbytes ]")

      if (allocated(ffc)) deallocate(ffc)
      allocate(ffc(nex,nex,0:lamax,1:nrcc))
      ffc=0

      write(*,170) nrcc, hcm,rmaxcc,hcm
170   format(/,5x,"=>  Interpolating F(r) in grid with ",i4, ' points:',
     &           2x,"[ Rmin=",1f5.2,2x," Rmax=",
     &           1f6.1,1x," Step=",1f6.3, " fm ]",/)


10    read(kfr,500,end=600) npt,rstep,rfirst,fscale,
     & lam,ptr,ttr,m2,m1,comment
 
!!!!!!!! CONVERSION FACTOR FROM FRESCO !!!!!!!!!!!!!
      jpi=jpch(m1)
      jpf=jpch(m2) 
      frconv=(2d0*ptr+1)*sqrt(2*jpi+1)*sqrt(2*jpf+1)
     &      *(-1)**ptr
      if (debug) then
      write(*,*)'o F(R) -> F(R) [fres]=',frconv
      endif
      if (verb.ge.4) write(*,'(i3," =>",i3, " with LAM=",i3,
     & 3x," o F(R) -> F(R) [fres]=",1f8.4)') 
     & m1,m2,lam,frconv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

500   format(i4,3f8.4,i4,2f4.0,2i4,a35)
      nff=nff+1
      if (nff.eq.1) allocate(faux(npt),rv(npt))
      do ir=1,npt
         read(kfr,*) x,y
         rv(ir)=rfirst+dble(ir-1)*rstep
         faux(ir)=fscale*cmplx(x,y)/frconv
!         write(98,'(3f16.6)') rv(ir),faux(ir)
      enddo !ir 
c     interpolate & store in integration grid
      do ir=1,nrcc
!        r=hcm*ir ! grid used in scattcc
!        r=dble(ir-1)*hcm    ! CHANGED in v2.2
        r=rvcc(ir)
        if (r.gt.rv(npt)) cycle ! DO NOT EXTRAPOLATE 
        caux=cfival(r,rv,faux,npt,alpha)
        ffc(m1,m2,lam,ir)=caux
!        write(98,'(3f16.6)') rv(ir),ffc(m1,m2,lam,ir)
      enddo ! ir
      goto 10

600   write(*,'(5x,"[",i6, " formfactors read from unit ",i2,"]")')
     &     nff,kfr
      deallocate(faux,rv)
      call flush(6)
      return
      end subroutine



c Read STATES specifications from states.fr file (fresco style)
      subroutine read_states()
      use xcdcc, only: exch,parch,jpch,nex
      implicit none
      namelist /states/jp,bandp,cpot,jt,copyt,bandt,ep,et,ptyp,ptyt
      integer, parameter:: kst=20
      integer ie,nst
      integer bandp,bandt,copyp,copyt,cpot,ptyp,ptyt
      real*8 jp,jt,ep,et
      CHARACTER PARITY(3)
      DATA PARITY / '-','?','+' / 

      write(*,*)
      write(*,*)' *** STATES FOR CC CALCULATIONS ***'
      
      open(kst,file='states.fr',status='old') 
      read(kst,*) nst
      write(*,'(3x,"Reading",i4," states from file states.fr")')
     &  nst
c let's consider all of them, BY NOW
      nex=nst
      if(allocated(exch)) then
        write(*,*)'read_states: trying to reallocate exch!'
        stop
      endif
      allocate(exch(nex))
      allocate(parch(nex))
      allocate(jpch(nex))

      do ie=1,nex
      ep=0; et=0; ptyp=0; ptyt=0;copyt=0; copyp=0;
      read(kst,nml=states)
!      write(*,nml=states)
      exch(ie)  =ep
      parch(ie) =ptyp
      jpch(ie)  =jp
      write(*,'(5x,i3,3x,1f6.3,3x,1f5.1,a1)') 
     &   ie,ep,jpch(ie),parity(ptyp+2)
      enddo !ie
      end subroutine


