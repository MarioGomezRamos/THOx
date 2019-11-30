c --------------------------------------------------------
c *** Modules
c -------------------------------------------------------
      module globals
        logical:: debug,written(399)
        integer:: kin,kout,verb,nfacmax
        real*8:: egs,mu12,kgs
      end module globals

      module memory
        integer :: lr8,lc16
        real*8  :: tzfun,torto,tcc,tcsh,tmatch,t3d
      end module memory

      module sistema
       integer:: nmult,ncoul
       integer:: partgs,ntex !MGR
       real*8:: zc,zv,zp,zt,rcc,mp,mt,jtgs
       character*25, namep, namet
       parameter(nmult=7)
!       real*8:: bel(0:nmult)
c 	real*8:: zv,zc
        real*8:: av,ac
c	real*8:: sn
      end module sistema


!      module lst
!        real*8 gamma, bosc
!      end module lst

!      module io
!        integer kin,kout
!        logical written(399)
!      end module io


      module trace
        logical:: wfprint(1:200)
        integer:: wfout(1:10)
        integer:: solapout(1:10)
        logical:: cdccwf
      end module trace

c    
c Some constants, as defined in FRESCO
c
      module constants
        real*8:: pi=acos(-1d0)
        real*8:: finec
        real*8:: hc
        parameter(finec=137.03599d0)
	parameter(hc=197.32705 ) ! MeV.fm
!       real*8:: e2=1.43997         ! MeV.fm
        real*8:: e2=hc/finec
        real*8:: amu=931.49432        ! MeV
!       real*8:: pi=2.*datan(-1.d0)
      end module constants
        

      module parameters
        integer, parameter:: maxcore=3    ! max number of core states
        integer, parameter:: maxchan=10    ! projectile channnel configs
        integer, parameter:: maxpauli=5   
        integer, parameter:: maxsets=100    ! max number of jp/parity sets
        integer           :: maxl=10
        integer, parameter:: maxeset=100
        integer, parameter:: maxch  =1000  ! max number of channels per JTOT/PI
      end module parameters


      module hmatrix
        use parameters
        logical xcdcc
        real*8,allocatable:: onm(:),hmat(:,:),ortmat(:,:)
        real*8,allocatable:: hmatx(:,:)
        real*8,allocatable:: gmat(:,:)
        real*8,allocatable:: eigv(:,:),idmat(:,:) 
        integer:: nrot=0,ieigen=0
      end module hmatrix


      module forbidden
        use parameters
        character*30:: paulifile(maxpauli)
        logical:: wfpread(maxpauli) 
        logical:: hindrance
        integer:: npauli,paul(maxpauli), paun(maxpauli)
        real*8:: pauj(maxpauli),eshift
        real*8,allocatable,target::wfpauli(:,:)
	real*8,allocatable,target:: ppauli(:,:),qpauli(:,:)
        real*8:: wblock(maxcore),pshift(maxpauli,maxcore)
        real*8:: pscale(maxpauli,maxcore)
        real*8:: uocup(maxpauli), vocup(maxpauli),eqsp(maxpauli)
      end module forbidden      


      module channels
        use parameters
        real*8 :: jtot,sn
        integer:: partot
        integer:: nchan,nsp,nce,bas2,nchmax
        integer:: ql(maxchan),qspl(maxchan)
        integer:: spindex(maxchan) ! identify the s.p. channel associated to this core/valence channel
        integer:: cindex(maxchan)
        integer:: qpar(maxchan)
        real*8 :: qj(maxchan), qspj(maxchan),qjc(maxchan),exc(maxchan)
        integer:: jpsets,inpsets
        integer:: indjset(maxsets)
 
      TYPE SPCHANNEL
          integer  l
          real*8:: j,jc,sn
      END TYPE
!      TYPE(SPCHANNEL):: spchan(maxsets,maxchan) 
      TYPE(SPCHANNEL),allocatable:: spchan(:,:) 


      TYPE PCHANNEL
          integer  bastype         ! basis type (1=THO, 2=CC Bins...)
          integer  partot          ! parity of composite
          real*8:: jtot            ! total spin of composite
          real*8:: sn              ! spin of valence (not really needed here)
          integer  nho             ! number of radial (PS) functions
          integer  bas2            ! digonalization procedure (see manual)
          integer  nsp             ! number of sp states retained after Hsp diagonalization
          integer  nchan           ! number of channels (i.e. c+v configurations) 
          integer  nchsp           ! number of valence configurations
          real*8:: exc(maxchan)    ! excitation energy of core for each channel
          real*8:: jsp(maxchan)    ! j of valence for each channel
          integer  lsp(maxchan)    ! l of valence for each channel
          integer  spindex(maxchan)! s.p. configuration index
          real*8:: jc (maxchan)    ! core spin for each channel
          integer  cindex(maxchan) ! core index for each channel
          integer  parc(maxchan)   ! core parity for each channel
          real*8   exmin,exmax     ! min, max eigenvalues to be retained
          real*8   wcut(maxchan)   ! minimum weight per chan to be retained in CDCC calcs
          integer  nex             ! number of energies retained 
          integer  inc             ! incoming channel (bins only)
          logical  complex         ! real/complex wfs?
!          type(spchannel) qnsp(maxchan)
      END TYPE
      TYPE(PCHANNEL),allocatable:: jpiset(:)

      TYPE PTCHANNEL
!          use parameters !, only: maxch
!          integer, parameter:: maxch=500
          integer partot           ! parity of system
!          integer:: mmm=60000
          real*8  jtot             ! total angular momentum (Jp+L+Jt)
          real*8  jt(6000)               ! spin of target MGR
          real*8  ext(6000)        ! energy of target MGR
          real*8  idt(6000)        !index of target MGR
          real*8  nchan            ! number of channels (p+t configurations)
          logical interp           ! true if this J/pi set will be obtained by interpolation
c         .................... can we allocate these vars dynamically?????????????
          integer l(6000)       ! p-t orbital angular momentum
          integer idx(6000)     ! index of projectile state
          real*8  jp(6000)      ! projectile spin
          real*8  jlp(6000)      ! J=L+Jp MGR
          real*8  exc(6000)     ! excitation energy of projectile
          real*8  kcm(6000)     ! wavenumber
          integer iex           ! index of excited state

      END TYPE
      TYPE(PTCHANNEL),allocatable:: jptset(:)
         
!MGR------------------------------------------------------------------------      
      TYPE TCHANNEL
          integer partarg        !Parity of target             
          real*8  jtarg          !Angular momentum of target
          real*8  extarg         !Excitation energy of target
          real*8 Krot           !Rotational projection of intrinsic angular momentum
          integer nphonon        !Vibrational number of phonons
          integer inc            !0 if not incoming channel, 1 if incoming channel
      END TYPE
      TYPE (TCHANNEL),allocatable:: tset(:)
      logical targdef
!-------------------------------------------------------------------------

c core state variables (deprecated; use qnc instead)
        real*8::  jc(maxcore)
        integer:: parc(maxcore) ! v2.0.5 
        integer:: nphon(maxcore)  ! v2.1b
        real*8 :: exc1(maxcore) 

c new derived type for core variables
        TYPE COREVAR  
         real*8 :: kband
         real*8 :: jc,exc
         integer:: parc,nphon
        END TYPE
        type(corevar) :: qnc(maxcore)       
        end module channels



      module wfs
        real*8 exmin,exmax
        integer:: ng,nr,hdim,nrint
        real*8::  rmin,rmax,dr,rlast,rint
        real*8,allocatable,target:: rvec(:),rweight(:) ! radial grid
        real*8,allocatable:: wchan(:,:),ebin(:)
        integer,allocatable:: iord(:)

        real*8,allocatable, target:: wftho(:,:,:) !n,l,r
        real*8,allocatable, target:: wftho2(:,:) !#(n,l),r !the order is first ncore and then each nho (l is then tied)
        real*8,pointer::wfaux(:,:)
        real*8,pointer::wfaux2(:,:)!,wfaux3(:)
        real*8,pointer::eaux(:)
        
        real*8,allocatable:: xref(:),uref(:)
        real*8,allocatable,target:: esp(:,:),wfsp(:,:,:) !sp energies and wfs
        real*8,allocatable,target:: wfeig(:,:,:) !final eigenfunctions (n,nchan,nr)
c       for all j/pi sets
        real*8,allocatable,target:: wfr(:,:,:,:)     !wfs for (jpiset,n,nchan,nr)
        complex*16,allocatable,target:: wfc(:,:,:,:) !wfs for (jpiset,n,nchan,nr)
        real*8,allocatable,target:: energ(:,:)       !energs for each jpiset
!        real*8,allocatable       :: wbin(:,:)        !energs for each jpiset
        real*8,allocatable,target:: idx(:,:)      
      end module wfs

      module xcdcc
      logical iftrans,realwf,rel
      integer nrad1,nquad,nex,numfilmax,nrad2,nrad3
      real*8 rin,dr, hin, rstep
      real*8, allocatable:: rvin(:),rvcc(:)
c frad changed to complex in v2.3
      complex*16,allocatable,target:: frad(:,:,:,:) !chann(:,:,:),
      real*8,allocatable,target:: xjp(:),xquad(:),wquad(:),
     &                            rquad(:),xrad1(:)
      complex*16,allocatable:: xintgn(:),xintgc(:)
      complex*16,allocatable:: potQKn(:,:,:,:),potQKc(:,:,:,:)
      complex*16,allocatable:: potQKn_targ(:,:,:,:,:), !MGR
     & potQKc_targ(:,:,:,:,:)                          !MGR
      complex*16:: potn,potc,fauxn,fauxc,caux,potnv,potcv,potncor,
     & potccor 
      complex*16,allocatable,target::Ff(:,:,:,:),Fc(:,:,:,:),Fn(:,:,:,:)
      complex*16,pointer:: ffr(:,:,:,:,:,:)
      parameter(nener=200,numfilmax=10) !,nchmax=10)
      integer, allocatable:: np(:)

c     Solving the CC (separate module??)
      integer:: method,lamax,nmatch,nrcc,lambmax,lambmin,lambpmax !MGR
      real*8,allocatable:: exch(:) ! ordered excitation energies for all JT/PI sets
      integer,allocatable:: parch(:) 
      real*8,allocatable:: jpch(:) 
      complex*16, allocatable:: ffc(:,:,:,:),smats(:,:,:),famps0(:,:,:)
      complex*16, allocatable:: famps(:,:,:,:) !MGR
      complex*16, allocatable:: ffcn(:,:,:,:),ffcc(:,:,:,:) !MGR
      real*8:: hcm,rmatch,elab,ecm,jtmin,jtmax,rmaxcc
      real*8:: jump(1:5),jbord(1:5) 

c ... to store bin information
      TYPE BIN_INFO
         logical :: ifhat
         integer :: nbins,nk(200)
         real*8 klow(200), kup(200), kmid(200),ebin(200),khat(200)
         real*8 wbin(200)
         real*8 emin, emax, kmin,kmax
         real*8 pshel(100,200)
         complex*16 wk(100,200)
      END TYPE
      TYPE(BIN_INFO), allocatable::binset(:) 

      !MGR lambdahipervector
      integer nlambhipr
      TYPE LAMBDAHIPER
        integer lambda,lambdap,q
      END TYPE LAMBDAHIPER
      
      TYPE(LAMBDAHIPER):: lambdahpr(10000)      

      end module xcdcc


c Scattering states of core+valence system
       module scattering
         integer:: nk,pcon,il,ilout,ili
         logical:: ifcont
         real*8:: emin,emax,kmin,kmax,eout
         complex*16,allocatable,target:: wfcont(:,:,:,:) !E,inchan,outchan,radius
         complex*16, allocatable,target:: solapmat(:,:,:,:)
         complex*16, allocatable:: solsum(:,:,:)
       end module scattering

       module coulwf
        real*8,allocatable,target:: wfcoul(:,:,:,:),wfcoulp(:,:,:,:)
        real*8,allocatable,target:: wfcoulpp(:,:,:,:)
       end module coulwf


       module potentials
         use parameters, only: maxl
         integer::ptype,lambda,pcmodel
         integer::cptype,lpot
         integer, parameter:: maxlamb=6
         logical:: laminc(0:maxlamb)
         real*8,allocatable,target::vcl(:,:)
         real*8,allocatable::vls(:,:),vss(:,:),vll(:,:),vlsc(:,:) !MGR
         real*8,allocatable::vtran(:,:,:,:,:)
         real*8,allocatable,target::vcou(:),vlcoup(:,:,:)
         complex*16,allocatable:: ccmat(:,:,:)
         real*8 :: beta,delta,kband ! fixed in 2.3b
       end module potentials



c for proj-target coupling potentials 
       module ptpots
         use sistema, only:nmult
         integer:: nr
         real*8, allocatable:: rfrag(:)
         real*8,pointer:: vfrag(:,:),vfragi(:,:),vcoup(:,:)
         real*8,allocatable,target::vcore(:,:),vcorei(:,:),vcorec(:,:)
         real*8,allocatable,target::vval(:,:),vvali(:,:),vvalc(:,:)
         REAL*8,target :: MELC(1:nmult),MELV(1:nmult)
         real*8:: MEL(1:nmult),MELT(1:nmult) !MGR
       end module ptpots

c--------------------------------------------------------------------------Added for laguerre cuadratures

      module factorials
        real*8:: dlfac(0:20500)
        real*8:: dl2fac(0:20500)
        real*8:: dlfac2(0:20500)
        real*8:: FACT(0:20500)
      end module factorials

c from DCE
      module factoriales
      real*8,allocatable:: fact(:)
      end module factoriales

      SUBROUTINE factorial(n)
      use factoriales
      implicit real*8 (a-h,o-z)
      allocate(fact(0:n))
      nfact=1
      fact(0)=dble(nfact)
      do j=1,n
      nfact=j*nfact
      fact(j)=dble(nfact)
      enddo
      return
      end


c--------------------------------------------------------------------------B(Elamda)

      module belambdamod
       character*40 :: uwfgsfile
        logical ifbel
        integer:: lambda,ncni,partoti,incni
        real*8,allocatable::qjci(:),qji(:),qlir(:),r(:)
!        real*8,allocatable:: ugs2(:,:)
        real*8,allocatable::dbde(:,:)
        complex*16,allocatable::ugs(:,:) 
        complex*16,allocatable:: mel(:,:)
        real*8:: zeff,besr,jtoti,sni,eneri
        integer,allocatable:: qli(:),cindexi(:)
      end module belambdamod


 
