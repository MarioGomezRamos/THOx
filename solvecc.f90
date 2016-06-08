c *** ---------------------------------------------------
c     Build & solve CDCC equations
      subroutine cdcc(iexgs,ncc)
c *** ---------------------------------------------------
      use channels, only:jpiset,jpsets,jptset
      use xcdcc   , only: nex,jpch,exch,parch,iftrans,
     &            lamax,hcm,nmatch,elab,ecm,smats,
     &            jtmin,jtmax,rvcc,rmaxcc,nrcc,method,
     &            jump,jbord
      use nmrv,only: hort,rmort,vcoup
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
      real*8 :: jtot,jp,jt,jtotmax,jblock
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
!      namelist /reaction/ elab,jtmin,jtmax,jump,jbord,mp,mt,zt,skip
      namelist /numerov/ method,
     &         hcm,rmaxcc,rcc,hort,rmort,
     &         jtmin,jtmax,jump,jbord,skip,
     &         nlag,ns ! R-matrix
c     ---------------------------------------------------------------


c ... initialize variables --------------------------------------------------
      dry=.true.
      skip=.false.
      jtmin=-1; jtmax=-1
      jump=0; jbord=0; 
c ----------------------------------------------------------------------------
!      write(*,'(//,5x,"******  REACTION CALCULATION ******")')
!      read(kin,nml=reaction)

      ecmi=elab*mt/(mp+mt)
      mupt=mp*mt/(mp+mt)*amu
      kcmi=sqrt(2*mupt*ecmi)/hc
      conv=(2*mupt/hc**2)
      etai=conv*zp*zt*e2/kcmi/2.
      zp  =zc+zv 

      write(*,'(//,5x,"******  CC SOLUTIONS ******")')


c ... Numerov specifications
      method=0 
      hort  =0
      hcm   =0
      rmort =0
      nlag  =0; ns=0
      read(kin,nml=numerov)

      if ((jtmin.lt.0).or.(skip)) then
        write(*,'(5x,"[No reaction calculation requested]",/)')
        return
      endif

      write(*,'(5x,"J interval [min,max]=",2f6.1)') jtmin,jtmax
      if (maxval(jbord).gt.0) then
!        write(*,*)' J jumps defined:'
        do ijt=1,5
        if (jbord(ijt).eq.0) cycle
        write(*,345) jbord(ijt),jump(ijt)
345     format(8x,'o from J=',1f6.1,' in steps of', 1f6.1)
        enddo
      endif


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
!      print*,'jtmin,jpgs,mod=',jtmin,jpgs,mod(2*(jtmin-jpgs),2.)  
      if (abs(mod(2*(jtmin-jpgs),2.)).gt.1e-3) then
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

c *** Channel energy
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
       jptset(icc)%jt(nch) = jt
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
        if (verb.ge.2) write(*,310) nch,rturn
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
!          if (jptset(icc)%idx(ich).eq.iexgs) then
         if (ich.eq.inc) then
            xsr=xsr+factor*(1-abs(smat)**2)
            xsrj=xsrj+factor*(1-abs(smat)**2)
          else 
            xsinel=xsinel  + factor*abs(smat)**2
            xsinelj=xsinelj+ factor*abs(smat)**2
          endif
        enddo ! ich
       enddo !inc
      if (allocated(vcoup)) deallocate(vcoup)

      
       write(*,320)xsr,xsinel,xsr-xsinel
320    format(5x,"o X-sections: Reac=", 1f8.3,5x,
     &   " Inel+BU=",1f8.3, 5x," Abs=",1f8.3," mb")

       enddo ! par
c       write(157,'(1x,1f6.1,a1,2x,3f12.6)') jtot, parity(partot+2),
c     &    xsrj-xsinelj,xsrj,xsinelj
       if (jptset(icc)%interp) cycle 
       write(156,'(1x,1f6.1,2x,3g14.5)') jtot, 
     &    xsrj-xsinelj,xsrj,xsinelj

       call flush(156)
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
!        write(0,*) 'Entering schcc_erwin Jtot',jptset(icc)%jtot,
!     &   jptset(icc)%partot
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
!      if (icc.eq.1) then 
!      do ir=1,nr
!      write(50,'(1f8.3,2x,50f12.8)') rvcc(ir),(wf(n,ir),n=1,min(2,nch))
!      enddo
!      write(50,*)'&'
!      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     Deallocate variables
!      if (allocated(vcoup)) deallocate(vcoup)
      deallocate(ql)
      end subroutine



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
c     .............................................................
      debug=.false.
c     .............................................................       
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
      li  =jptset(icc)%l(ni)
      lri =li
      jpi =jptset(icc)%jp(ni)
      idi =jptset(icc)%idx(ni)
      pari=parch(idi)
      do nf=ni,nch !Upper triangle available here only!!!
      lf=jptset(icc)%l(nf)
      lrf=lf
      jpf=jptset(icc)%jp(nf)
      idf =jptset(icc)%idx(nf)
      parf=parch(idf)
c1      r1=sqrt((2.*lri+1.)*(2.*lrf+1.)*(2.*jpi+1.)*(2.*jpf+1.))
c1     &  *(-1)**(jpi+jtot)

      r1=sqrt((2.*lri+1.)*(2.*lrf+1.)*(2.*jpi+1.)*(2.*jpf+1.))
     &  *(-1)**(jpf+jtot)


      do lam=0,lamax
      lamr=lam
!      if (fail3(jpi,jpf,lamr)) cycle
      if (fail3(lri,lrf,lamr)) cycle
      if ((-1)**lam*pari.ne.parf) cycle

c1      r2=threej(lamr,lri,lrf,zero,zero,zero)
c1      r3=sixj(jpi,jpf,lamr,lrf,lri,jtot)
      r2=threej(lamr,lrf,lri,zero,zero,zero)
      r3=sixj(jpf,jpi,lamr,lri,lrf,jtot)

      r13=r1*r2*r3
! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!
      phc=(0d0,1d0)**(lrf-lri)


      ix=ix+1
      coef=(-1)**lamr*(2.*lamr+1.)*r13
      frconv=(-1)**lamr*(2d0*lamr+1)*sqrt(2*jpi+1)*sqrt(2*jpf+1) ! FFC=FRCONV*FFR(fresco)
      if (verb.ge.4) then     
      write(*,300) ix,ni,idi,li,jpi,parity(pari+2),
     &                nf,idf,lf,jpf,parity(parf+2),
     &             lam,r1,r2,r3,coef,coef/frconv,phc !,ffc(idi,idf,lam,1)
300   format(5x,"IX=",i4," Ch. IEX L, Jpi=",3i3,f4.1,a1,
     &     " =>",3i3,f4.1,a1,3x,
     &     'with LAM=',i3,4x,'r1-3='3f8.3,
     &     2x,'COEF=',1f7.3,' COEF(FR)=',1g10.4,
     &  ' phc=',2f8.4)
!     &     2x,'COEF=',1f7.3,' F(1)=',2g10.4)
      endif
      do ir=1,nrcc      
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
       if (ir.eq.1) write(97,'("#",2i3,2i4,i3,2g12.4)') 
     & ni,nf,idi,idf,lam,coef,coef/frconv
        write(97,'(f8.3,10g16.4)')rvcc(ir),vcoup(ni,nf,ir),
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
C     Use previously calcualted F(R)'s and interpolate in radial grid
c     used to solve the CC equations  
      subroutine int_ff()
c *** --------------------------------------------------------------------      
      use xcdcc,   only: hcm,ffc,lamax,nex,nrcc,jpch,Ff,lambmax,
     &                   nrad3,rstep,parch,rmaxcc,rvcc
      use channels,only: jpiset,jpsets
      use wfs     ,only: idx
      use memory
      use globals ,only: verb
      implicit none
      integer :: ir,irr
      integer :: m1,m2,nff,lam,par1,par2
      complex*16, pointer:: faux(:)
      real*8  :: rv(nrad3),xpos,rfirst
      real*8 :: jp1,jp2,r,ymem
      real*8, parameter:: alpha=0d0
      complex*16 cfival,caux,fic,caux2,ffc4
c     ----------------------------------------------------
     
      write(*,'(//,3x, "** INTERPOLATE FORMFACTORS **" )')

      lamax=2*maxval(jpiset(:)%jtot)

c ... Memory requirements
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
170   format(/,5x,"=>  Interpolating F(R) in grid with ",i4, ' points:',
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
       do lam= nint(dabs(jp1-jp2)),min(nint(jp1+jp2),lambmax) !nint(jp1+jp2)
       if (par1*par2*(-1)**lam<0) cycle
       nff=nff+1
       faux =>  Ff(m1,m2,lam,:)
c      interpolate & store in integration grid
       do ir=1,nrcc
!        r=hcm*ir ! grid used in scattcc
!        r=dble(ir-1)*hcm ! grid used in scattcc (Changed in v2.2)
        r=rvcc(ir)
        if (r.gt.rv(nrad3)) cycle ! DO NOT EXTRAPOLATE 
!        caux=cfival(r,rv,faux,nrad3,alpha)
        xpos=(r-rstep)/rstep
        caux=FFC4(xpos,faux,nrad3)
        ffc(m1,m2,lam,ir)=caux
!        if (abs(caux).gt.1e20) then 
!          write(90,*)'# ** INT_FF: F(R)=',caux, 
!     &  'for ir,lam,m1,m2',ir,lam,m1,m2
!         write(90,'(1f8.3,2x,2g16.5)') (rv(irr),faux(irr),irr=1,nrad3)
!         write(90,*)'&'
!         stop
!        endif
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
      real*8 jpi,jpf,frconv,ymem,big,small
      real*8, parameter:: alpha=0d0
      integer,parameter :: kfr=4
      complex*16 cfival,caux
      character*40 comment
c     ...............................................................   
      big=huge(big)
      small=epsilon(small)
  
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
      frconv=(-1)**ptr*(2d0*ptr+1)*sqrt(2*jpi+1)*sqrt(2*jpf+1)
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
         if (abs(x).lt.small) x=small
         if (abs(y).lt.small) y=small
         rv(ir)=rfirst+dble(ir-1)*rstep
         faux(ir)=fscale*cmplx(x,y)/frconv
      enddo !ir 
c     interpolate & store in integration grid
      do ir=1,nrcc
        r=rvcc(ir)
        if (r.gt.rv(npt)) cycle ! DO NOT EXTRAPOLATE 
        caux=cfival(r,rv,faux,npt,alpha)
        ffc(m1,m2,lam,ir)=caux
      enddo ! ir
      goto 10

600   if (nff.gt.0) then
      write(*,'(5x,"[",i6, " formfactors read from unit ",i2,"]")')
     &     nff,kfr
      deallocate(faux,rv)
      else
      write(*,*)'No formfactors found in unit 4!';stop
      endif
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


