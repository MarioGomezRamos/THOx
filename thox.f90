c Calculates THO functions for a two-body system (core+valence) with core excitation
c
c by JA Lay and AM Moro (2011)
c

c v1.1 adding B(Elambda) 8/7/11
c v1.2 B(El) sum rule, named changed to thox on 15/6/2011 (original name:tho2frx)
c v1.4 September 2011. Added calculation of scattering states (bincc)   ! last stable version
c v1.5 B(El) for continuum states, B(El) folded added on 26/9/2011
c v1.5 Also core E2 contribution to B(El) and general makeup
c v1.6 Print PS wfs for XCDCC                                           ! minor mistakes but good for non multiple 0+'s
c v1.7 phaseshifts & widths 26/11/2011 (Not implemented yet)            ! totally crashed...
c v1.8 Pauli hindrance + Correction in B(El) for multiple 0+ cores 31/01/12

c v2.0 Stable, B(El), Solap and writting fort.79 errors corrected + Change in core index 3/02/12
c v2.0.3 core spin orbit (l.jc) implemented + read potentials from transition densities (fort.30)
c v2.0.5c change bincc by scatcc for scattering states calculation + Corrections in belambda
c         we get back to core index definitions previous to v2.0
c v2.1 general makeup 9/4/2013 !stable
c v2.1c Includes spin.spin central potential (PS and continuum wfs)
c v2.2 April 2014: allow for several j/pi sets on same input and calculates coupling potentials
c                  for (X)CDCC calculations, using DCE subroutine
c v2.2e Problems with reading fort.30 potentials solved, continue to read core+target potentials
c v2.3 AMM: CC bins included for the first time!
c           Change input ordering 

      program thox
      use globals
      use sistema
      use parameters
      use constants
      use hmatrix
      use trace
      use forbidden, only: npauli,ppauli,eshift,qpauli,hindrance
      use channels
      use belambdamod, only: ifbel
      use memory
      implicit none
      logical ehat
      integer parity,l,ic,iic,i,ir,ichsp,ncc,bastype,nk
      real*8 eps
     
      parameter (eps=1e-6)

      logical fail3,checkorth,ifphase,dummy,tres
      integer:: nfmax,mlst
      real*8 ::  lambda,r,norm
      integer :: al
      CHARACTER*1 BLANK,PSIGN(3)
      integer :: nset,nho, nchsp

!!! TEST 
      integer iset,inc,iexgs
      real*8 ecm


c Input namelists -------------------------
      namelist/system/ Zc,Ac,Zv,Av,egs,sn
      namelist/output/ wfout,checkorth,verb,ifphase,solapout,xcdcc

      DATA PSIGN / '-','?','+' /, BLANK / ' ' /
      
      interface 
      subroutine potdef(vc,qmax,rnc,vdef,betar)
        integer::qmax
        real*8:: vc(:),rnc,vdef(:),betar
      end subroutine
      end interface

      interface
      subroutine potdef2(rvec,nr,vc,qmax,rnc,vdef,betar)
      integer nr,qmax
      real*8 vc(:),rvec(:),vdef(:),betar,rnc
      end subroutine
      end interface

      interface 
        subroutine coulph(eta,cph,lmax)
        real*8 eta,cph(0:lmax) 
        integer lmax
        end subroutine
      end interface

      interface 
        subroutine cfival(lrf,lfv,gaux,nlf,alpha)
        real*8 lrf,lfv(:),gaux(:),alpha
        integer nlf
        end subroutine
      end interface

      interface 
        subroutine test_cont(nset,nchan,inc,ecm,wf)
        integer   :: nset,nchan,inc
        real*8    :: ecm
        complex*16:: wf(:,:)
        end subroutine
      end interface 
      
      interface
        subroutine wfrange(iset,nchan,inc,emin,emax,
     &             ncont,energy,wfcont)
        integer :: iset,nchan,inc,ncont
        logical :: energy
        complex*16 :: wfcont(:,:,:)
        real*8 emin,emax 
        end subroutine
      end interface
c     ------------------------------------------------------------

      debug=.false.

c *** Defined global constants
      call initialize() 
      write(*,'(50("*"))')
      write(*,'(" ***",8x,"THOx+DCE+CC code: version 2.3",8x, "***")')
      write(*,'(50("*"),/)')

c *** Print physical constants
      write(*,*)'***************** CONSTANTS **********************'
      write(*,'(" * hbarc=",1f9.5," MeV.fm",5x,
     & "e^2=",1f8.5," MeV.fm *")') hc,e2
      write(*,'(" * so, alpha= 1/",1f9.5,$)')hc/e2
      write(*,'(7x,"amu=",1f8.4, " MeV  *")') amu
      write(*,*)'**************************************************'


c *** Calculate memory sizes for memory allocation
      call memory_sizes()

c *** Calculate & store factorials
      call factorialgen(1000) 

c *** Read main input
      kin=5

c *** Pre-read input to set arrays dimensions
      call preread(kin)

c *** Actual read
      read(kin,nml=system)
      mu12=amu*ac*av/(ac+av)
      kgs=sqrt(2d0*mu12*abs(egs))/hc

      write(*,200)Ac,Zc,Av,Zv 
      write(*,'(/,1x,"- Reduced mass:",1f9.3," MeV/c2")') mu12
200   format(' - COMPOSITE:',/,5x, 'Core:    Ac=',1f8.4,2x,'Zc=',1f6.2,
     &      /,5x,'Valence: Av=',1f8.4,2x,'Zv=',f6.2)

       
c *** Read core states (last state followed an empty namelist)
      call readcorex(kin)
 
 
c *** Output trace
      wfprint(:)=.false.
      ifbel=.false.
      checkorth=.false.
      ifphase=.false.
      wfout(:)=0
      solapout(:)=0
      read(*,nml=output)
      do i=1,10
       if (wfout(i)>0) wfprint(wfout(i))=.true.
      enddo

c *** Radial grid (quadratures?)
      call readgrid(kin)

c *** Read potential parameters
      call readpot(kin)

c Pauli-forbidden states to be removed
      call pauli_forbidden(kin)

      do nset=1,jpsets
      dummy=.false.
      call read_jpiset(nset,bastype,nk,tres,ehat)

      jtot    =jpiset(nset)%jtot
      partot  =jpiset(nset)%partot
      nho     =jpiset(nset)%nho
      nchan   =jpiset(nset)%nchan
      nchsp   =jpiset(nset)%nchsp
      cindex(1:nchan)=jpiset(nset)%cindex(1:nchan)  !Lay: trying to recover fort.30

c these variables could be used to replace qspl, qspj, qjc(nsp)
!      qspl(:) =jpiset(nset)%lsp(:)  
!      qspj(:) =jpiset(nset)%jsp(:)
!      qjc (:) =jpiset(nset)%jc(:)
      bas2    =jpiset(nset)%bas2
      nsp     =jpiset(nset)%nsp


      select case(bastype)
      case(0,1) ! HO, THO -----------------------------------------
c *** Build and diagonalize s.p. Hamiltonian for each j/pi set
        if (bas2.eq.0) call hsp(nset,nchsp,nho)

c *** Build and diagonalize FULL VALENCE+CORE Hamiltonian 
        call fullham(nset,nho)

c *** Check orthonormality
        if (checkorth) call orthonorm

c *** Overlap between THO and scattering wfs
        call solap
c *** ------------------------------------------------------------
      case(2) ! Bins   
        call makebins(nset,nk,tres,ehat)
      end select 
      
      enddo ! loop on j/pi sets

c *** Continuum WFS
!      call continuum!(nchan)

c *** B(E;lambda)      
      call belam ! (kin,jtot,partot)

c *** Phase-shifts
      if (ifphase) call phase


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! TEST calculation of scattering wf for arbitrary j/pi and energy
!      iset=1
!      ecm= 1.
!      inc=1
!      call test_cont(iset,jpiset(iset)%nchan,inc,ecm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C *** Capture x-sections and S-factors (unfinished)
!!!      if (ifbel) call capture(jtot,partot)


c *** Define reaction (projectile, target, energy)
      call partition()

c *** Coupling potentials (includes DCE routines)
      if (.not. targdef) then
      write(*,*) 'using transition subroutine Target Spin=0'
      call transition 
      else
      write(*,*) 'using transition_targdef subroutine 
     & with target excitation'
      call transition_targdef
      endif

c *** Build & solve CDCC equations
      iexgs=1 !assume that the incident channel is the first state!!!!
      ncc  =0 
      call cdcc(iexgs,ncc)

c *** Scattering amplitudes and cross sections
      if (.not. targdef) then
      call xsecs(kin,ncc,iexgs)
      else
      call xsecs_tdef(kin,ncc,iexgs)
      endif
      

!      tripple=.true.
!      if (tripple) call three_body()

      call fkind(written,kout)

      call timings()

      end program


c *** ---------------------------------------
c *** Define PROJECTILE and TARGET properties 
c *** ----------------------------------------
      subroutine partition()
      use xcdcc
      use sistema
      use globals, only: kin
      use constants
      use channels, only: tset,targdef
      implicit none
      integer par,itex,part,nphon,inc,notgdef
      real*8 :: ecmi,mupt,conv,etai,kcmi,It,et,jt,K !MGR
      character*1::pchar !MGR
      namelist /reaction/ elab,namep,namet,mp,mt,zt,jt,par,ntex,nphon,
     & K,notgdef !MGR
      namelist /targstates/ et,It,part,nphon,K,inc !MGR
c           
      write(*,'(//,5x,"******  REACTION CALCULATION ******")')
      jt=0.0; zt=0.0; mt=0.0 
      ntex=0; nphon=0; K=0; inc=0;par=1 !MGR
      targdef=.false.                   !MGR
      notgdef=0                         !MGR
      
      read(kin,nml=reaction)
!MGR--------------------------------------------------------------------
      jtgs=jt
      partgs=par
      if (abs(jt).gt.0.2d0) targdef=.true.
      if(ntex.le.0) then
        allocate(tset(1))
        tset(:)%Krot=0d0
        tset(1)%partarg=par
        tset(1)%jtarg=jt
        tset(1)%extarg=0d0
        tset(1)%Krot=K
        tset(1)%nphonon=nphon
        tset(1)%inc=1
      else
        allocate(tset(ntex+1))
        tset(:)%Krot=0d0
        tset(1)%partarg=par
        tset(1)%jtarg=jt
        tset(1)%extarg=0d0
        tset(1)%Krot=K
        tset(1)%nphonon=nphon
        tset(1)%inc=1
        targdef=.true.
      endif
      
      if (ntex.gt.0) then
        do itex=1,ntex
          nphon=0; K=0; inc=0
          read(kin,nml=targstates)
          tset(itex+1)%partarg=part
          tset(itex+1)%jtarg=It
          tset(itex+1)%extarg=et
          tset(itex+1)%Krot=K
          tset(itex+1)%nphonon=nphon
          tset(itex+1)%inc=inc
          if (inc.ne.0) then
             tset(1)%inc=0 !Not fully implemented BEWARE
          endif
        enddo
      endif
      
      if (notgdef.gt.0) targdef=.false.
      
      write(*,*) ntex+1,' target states'
      do itex=1,ntex+1
      if (tset(itex)%partarg.eq.1) pchar='+'
      if (tset(itex)%partarg.eq.-1) pchar='-'
      write(*,*) 'State #',itex,': J/pi',tset(itex)%jtarg,pchar,
     & ' E:',tset(itex)%extarg
      enddo
!_----------------------------------------      

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
      end subroutine


      subroutine timings()
      use memory
      write(*,*)''
      write(*,*)'Timings:'
      write(*,*)' CC integration :',tcc,  ' secs'
      if (tzfun.gt.1e-4) write(*,*)' Calculation of zfun:',tzfun,' secs'
      if (torto.gt.1e-4) write(*,*)' Stabilization:',torto,' secs'
      if (tcsh.gt.1e-4) write(*,*)' Cosh :',tcsh,' secs'
      if (tmatch.gt.1e-4) write(*,*)' Matching :',tmatch,' secs'
      if (t3d.gt.1e-4) write(*,*)' Triple x-sections:',t3d,' secs'

      end subroutine


c *** Check memory required by real*8 and complex*16
      subroutine memory_sizes()
      use memory
      implicit none
      real*8 r8
      complex*16 c16
!      integer lr8, lc16
      inquire(iolength=lr8) r8
      write(*,*) '  A real*8 is ',lr8,' bytes'
      inquire(iolength=lc16) c16
      write(*,*) '  A complex*16 is ',lc16,' bytes'
      end






c
c Pre-read input to set dimensions
c
      subroutine preread(kin)
      use channels 
      use parameters, only: maxl
!      use potentials, only: vl0,vcp0
      implicit none
      logical :: tres,ehat
      integer::l,bastype,mlst,kin
      integer:: basold,parold,incold 
      real*8:: jn,exmin,exmax,jtold
      real*8:: bosc,gamma
      integer:: ichsp,ic,iset,nchsp,nfmax,nho,parity
      integer:: nk,nbins,inc
      
      
      namelist/jpset/ jtot,parity,lmax,
     &                bastype,nfmax,exmin,exmax,
     &                nho,bosc,     ! HO
     &                gamma,mlst,   !THO Amos 
     &                nsp,  ! sp eigenvalues to keep for full diag
     &                bas2,
     &                nk, nbins,inc,tres,ehat

      
      jpsets=0
      maxl=0
      
300   bastype=-1    
      tres=.false.
      inc=1
      read(kin,nml=jpset) 
      if (bastype.lt.0) goto 350
      jpsets=jpsets+1
      if (jpsets.eq.1) then
         jtold=jtot
         parold =parity
         basold=bastype
         incold=inc
      else
         if ((jtold.eq.jtot).and.(parity.eq.parold).and.
     x   (bastype.eq.basold).and.(inc.eq.incold)) then
         write(0,*)'Set',jpsets,' lecould be merged to previous one!'
         else
         jtold=jtot
         parold=parity
         basold=bastype
         incold=inc
         endif
      endif
      if (lmax.gt.maxl) maxl=lmax
      goto 300 

350   continue

      write(*,*)'' 
      if (jpsets.eq.0) then
        write(*,'(4x,"[No basis sets have been defined]")')
        rewind(kin)
        return
      endif 
      write(*,*)' Dimensions for memory allocation'
      write(*,*)'   o J/pi sets:', jpsets
      write(*,*)'   o Maximum internal l:', maxl

      write(*,*)''


      allocate(jpiset(jpsets)) ! sets, chans/set
      allocate(spchan(jpsets,maxchan))
      rewind(kin)
      end subroutine


c ***
c *** Read core states
c ***
      subroutine readcorex(kin)
       use channels
       implicit none
       integer:: kin
       integer:: nph,parity
       real*8 :: spin,ex
       character*1 BLANK,PSIGN(3)      
       DATA PSIGN / '-','?','+' /, BLANK / ' ' /

       namelist /corestates/ spin,parity,ex,nph


       nce=0 ! number of core states
       write(*,*) '- CORE states:'
50     ex=-100.
       parity=0
       nph=0
       read(kin,nml=corestates) 
       if (ex<-99) goto 100
       nce=nce+1   
       if (parity.eq.0) parity=+1
c backward compatibility (deprecated)
       jc(nce)=spin
       parc(nce)=parity
       exc1(nce)=ex
       nphon(nce)=nph
c new derived variable
       qnc(nce)%jc   = spin
       qnc(nce)%parc = parity
       qnc(nce)%exc  = ex
       qnc(nce)%nphon= nph

       write(*,230) nce,spin,psign(parity+2),ex
230    format(5x,"#",i2,": J/pi=",f4.1,a1,3x,'Ex=',f5.2,' MeV')
       goto 50

100    if (nce.eq.0) then
	write(*,*)'No core states found!'
       endif       
       end subroutine


c *** Read radial grid, nb of quadrature points,etc
       subroutine readgrid(kin)
       use globals, only:kgs
       use wfs,only:rvec,rmin,rmax,dr,rlast,rweight,ng,nr,rint
       implicit none
       integer ir,kin
       real*8:: eps,r,al,lambda
       parameter(eps=1e-6)
       namelist /grid/ ng, rmin,rmax,dr,rlast,rint

       ng=0
       rint=0.0
       write(*,*)
       write(*,*)' RADIAL GRID FOR INTERNAL WFS:'       
       read(kin,nml=grid)
       if ((rlast<0).or.(rlast>rmax)) rlast=rmax 
!       if (rlast<0) rlast=rmax 
       if (dr>0) then
           nr=ceiling((rmax-rmin)/dr)+1
         allocate(rvec(nr))
	   do ir=1, nr
              rvec(ir)=rmin+dble(ir-1)*dr
          enddo
        else
	   if (ng.eq.0) then
           write(0,*)'dr too small!. Aborting'; stop
	   endif
        endif
        if (ng>0) then
           write(*,*) '- Using', ng, 'quadrature points'
	   nr=ng
	   if (allocated(rvec)) deallocate(rvec)
	   allocate(rvec(nr))
	   allocate(rweight(nr))
	   rvec(:)=0d0
	   rweight(:)=0d0
c	   write(*,*)rmin,rmax,nr
	   call factorialgen(1000)
	   al=0
	   call qualag(rvec,rweight,nr,al)
c	   lambda=1
c	   lambda=rvec(nr)/rmax
	   lambda=0.6/kgs
	   do ir=1,nr
	      r=rvec(ir)
	      rweight(ir)=rweight(ir)/(r**al)*dexp(r)/lambda
	      rvec(ir)=rvec(ir)/lambda
	      write(20,*)rvec(ir),rweight(ir)
	   enddo
c	   dr=(rmax-rmin)/(nr-1)               !provisional
c	   do ir=1,nr
c              rvec(ir)=rmin+(ir-1)*dr          !provisional
c           enddo
        else
             write(*,170) nr, rmin,rmax,dr
170     format(1x,"- Using uniform grid with ",i4, ' points:',
     &           2x,"[ Rmin=",1f6.3,2x," Rmax=",
     &           1f6.1,1x," Step=",1f6.3, " fm ]",//)
         if (rint>0) write(0,*)'rint=',rint,' fm'
	endif
       end subroutine 


     



c *** ------------------------------------------------------
c *** Set value of global constants
c *** ------------------------------------------------------
        subroutine initialize()
          use globals
          use constants
          use xcdcc, only:  nex
          use channels, only: jpsets
          use wfs, only: rlast
          use memory
          implicit none
          real*8::cero,uno
          cero=0.
          uno=1.
          pi=acos(-1.)
          written(:)=.false.
          rlast=-1d0
          nex=0
          jpsets=0
          tzfun=0
          tcsh=0
          tmatch=0
          t3d=0
          return
        end subroutine initialize








c Write output files units and names
       subroutine fkind(written,KO)
       character*40 flkind(399)
       integer writf(399),nwrit
       logical written(399)
       flkind(:) = ' '
       flkind(20)='Val-core central potential for each l'
       flkind(22)='Val-core coupling potential for lambda=2'
       flkind(23)='Val-core coupling potential for lambda=3'
       flkind(40)='Valence-core potentials'
       flkind(44)='Fragment-target potentials'
       flkind(45)='Phase-shifts'
       flkind(50)='Full list of Cont. Wf.'
	   flkind(51)='Continuum Wavefunction'
       flkind(60:65)='wfs of Hsp eigenstates'
       flkind(70)='Full H matrix for each channel block'
       flkind(72)='Orthog. matrix for each channel block'
       flkind(74)='Full Hamiltonian'
       flkind(80)='Pauli projector for each chan block'
       flkind(101:120)='Channel WFs of H eigenvalues'
       flkind(76)='Single-particle Hamiltonians'
       flkind(78)='Channel decomposition of eigenfunctions'
       flkind(79)='List of full H matrix eigenvalues'
!       flkind(82)='Single-particle wfs for 1st channel'
       flkind(84)='WFS for XCDCC calculations'
       flkind(90)='Orthog. matrix for final wfs'
       flkind(95)='Discrete Distribution of B(Elambda)'
       flkind(97)='Discr. Dist. of B(Elambda) folded'
       flkind(97)='Discr. Dist. of B(Elambda) folded'
       flkind(96)='B(Elambda) disrib. from scatt. wfs'
       flkind(136)='Angle-independent scattering amplitudes'
       flkind(137)='Scattering amplitudes'
       flkind(156)='Inel/BU x-sections for each J'
       flkind(170)='Elastic S-matrix elements'

       flkind(201:249)='Momentum distribution of eigenstates'
       flkind(250)='Density of states vs energy'        
       flkind(251:299)='Momentum dist. vs energy'

       flkind(300)='All diff. x-sections'
       flkind(301:310)='Separate diff. x-sections'

	nwrit = 0
	do i=1,399
	 if(written(i)) then
		flkind(i) = trim(flkind(i))//'.'
	   nwrit = nwrit+1
	   writf(nwrit) = i
	 endif
	enddo
	write(*,990) (writf(i),flkind(writf(i)),i=1,nwrit)
990	format(/'  The following files have been created:',
     X	/(1x,2(i3,':',1x,a40)))
	return
	end




c ---------------------------------------------------------------
c Calculates overlap between scattering states and PS eigenstates
c ----------------------------------------------------------------
      subroutine solap
      use globals, only: written,verb,mu12
      use wfs
      use channels
      use constants, only: hc,pi
      use scattering
      use trace, only: solapout
      implicit none
      logical solprint(hdim)                
      integer n,m,ik,ie,ir,iil,iout
      real*8 :: dk,dkrdk(nchan),krel(nchan)
      complex*16 :: gaux(nr),caux
      complex*16:: zero,test2,resc
      real*8:: suma,total,suma2,test(il),kcont,econt
      complex*16, pointer:: sol

      if (.not.ifcont.or.nk.eq.0) return
        
      zero=(0d0,0d0)

      allocate(solapmat(nk,il,hdim,nchan))
      solapmat(:,:,:,:)=zero
      if (nk.gt.1) then 
          dk=(kmax-kmin)/(nk-1)
      else
          dk=(kmax-kmin)
      endif

      do ik=1,nk
        kcont=kmin+ (ik-1)*dk
        econt=(hc*kcont)**2/2/mu12
	 krel(1:nchan)=sqrt(2d0*mu12*(econt-exc(1:nchan)))/hc
	 dkrdk(1:nchan)=kcont/krel(1:nchan)

        test=0d0
        test2=0d0
        total=0d0
        do ie=1,hdim
        suma=0d0
	  do n=1,nchan
        do iil=ili,il
         sol=> solapmat(ik,iil,ie,n)
         gaux(:)=wfeig(ie,n,:)*conjg(wfcont(ik,iil,n,:))*rvec(:) 
         call simc(gaux,resc,1,nr,dr,nr)     


         sol=sol+resc
         suma=suma+abs(resc)**2

         test(iil)=test(iil)+abs(sol)**2  
        enddo ! iil
        enddo ! n
9           total=total+suma
        enddo ! ie (PS)

 
        test2=sum(test(ili:il))

        if (verb.ge.2) then
!           write(200,312) kcont,test2,(test(iil),iil=ili,il),total   !there are some differences between total & test2
!           write(250,312) econt,test2,(test(iil),iil=ili,il),total
           write(200,312) kcont,total,(test(iil),iil=ili,il) 
           write(250,312) econt,total,(test(iil),iil=ili,il)

        endif
312     format(1f10.4,2x,1f10.4,4x,10f10.4)  ! it only works when nchan=4
        enddo !ik

        if (verb.ge.2) then
        written(200)=.true.
        written(250)=.true.
        endif

        solprint(:)=.false.
        
        do ie=1,min(hdim,10)
        iout=solapout(ie)
        if ((iout>0).and.(iout.lt.hdim)) then
            solprint(iout)=.true.
            written(200+iout)=.true.
            written(250+iout)=.true.
        endif
!        write(*,*)'solapout=',solapout(ie)
        enddo


!-------------------- OVERLAPS AND NORMS ------------------
        total=0d0
	suma=0d0
	write(*,*)' - Printing Solaps with norms:'

        do ie=1,hdim
        total=0d0
		suma=0d0

        if (solprint(ie)) then
!        write(*,*)'ie=',ie
        do ik=1,nk
          kcont=kmin+(ik-1)*dk
          econt=(hc*kcont)**2/2/mu12
	   krel(1:nchan)=sqrt(2d0*mu12*(econt-exc(1:nchan)))/hc
	   dkrdk(1:nchan)=kcont/krel(1:nchan)
		  
         test(:)=0d0
        do n=1,nchan
           test(:)=test(:)+abs(solapmat(ik,:,ie,n))**2
        enddo !nchan

		do iil=ili,il

		  krel(iil)=sqrt(2d0*mu12*(econt-exc(iil)))/hc
		  if (krel(iil).gt.0.00001) then
		  dkrdk(iil)=kcont/krel(iil)
		  else
		  dkrdk(iil)=0d0
		  endif
!          dkrdk(:)=1d0 !test
		   test(iil)=test(iil)*dkrdk(iil)
		enddo !iil
		   suma=sum(test(ili:il))
		   total=total+suma*dk
		   
		   
        write(200+ie,313) kcont,(test(iil),iil=ili,il)
        write(250+ie,313) econt,(test(iil),iil=ili,il)
313     format(1f10.4,10f10.4)  ! it only works when nchan=4
        enddo

	write(*,*)'Norm(',ie,')=',total

        endif
        enddo !
! -------------------------- end of norms -----------------------

! Sum in outgoing channel (for 3-body observables) 
        allocate(solsum(nk,hdim,nchan))
        solsum(:,:,:)=zero
        do ik=1, nk
        do ie=1,hdim
        do iil=ili,il
        caux=0
        do n=1,nchan
        caux=caux+conjg(solapmat(nk,iil,ie,n))
        enddo
        solsum(ik,ie,iil)=caux
        enddo ! iil
        enddo ! ie
        enddo ! ik
! --------------------------------------------------------------
        end subroutine solap


  

        subroutine storecoulwf() !not working!!!!!!!!!!!!!!!!
        use globals
        use sistema
        use constants
        use wfs !remember to set an only:... to reduce memory
        use channels
        use coulwf
        use scattering
        implicit none
        real*8:: x,eta,k,k2
        integer:: in,ie,ir,ik
        integer:: ifail,m1
!        real*8:: f(lmax+1),g(lmax+1),fp(lmax+1),gp(lmax+1)
        real*8, dimension(0:lmax):: f,g,gp,fp
        real*8:: eps
        eps=1d-6
        ifail=0
!        allocate(wfcoul(nchan,hdim,nr,2))
        allocate(wfcoul(nchan,nk,nr,2))
        wfcoul(:,:,:,:)=0d0
!        allocate(wfcoulp(nchan,hdim,nr,2))
        allocate(wfcoulp(nchan,nk,nr,2))
        wfcoulp(:,:,:,:)=0d0

        do in=1,nchan
!        do ie=1,hdim
        do ik=1,nk
!        if (ebin(iord(ie)).lt.0d0) cycle
        
!        k2=2*mu12/hc**2*(ebin(iord(ie))-exc(in))
!        write(*,*)'Excitation energy:',exc(in)
!        k=sqrt(abs(k2))
         k=kmin+(ik-1)*(kmax-kmin)/(nk-1)
        k2=k-2*mu12/hc**2*exc(in)
        k=sqrt(abs(k2))
!        write(*,*)'k=',k
        eta=e2*zv*zc*mu12/hc**2/k

!        if ((in.eq.2).and.(ie.eq.45)) then
!        write(350,*)'#eta=',eta,'k=',k
!        write(350,*)'#l=',ql(in)
!        endif

        do ir=1,nr
        x=k*rvec(ir)
        if (x.eq.0d0) x=eps
        if (k2.lt.0d0)then
!        wfcoul(in,iord(ie),ir,1)=0d0        !change for whittaker functions
!        wfcoul(in,iord(ie),ir,2)=0d0
!        wfcoulp(in,iord(ie),ir,1)=0d0        !change for whittaker functions
!        wfcoulp(in,iord(ie),ir,2)=0d0
        wfcoul(in,ik,ir,1)=0d0        !change for whittaker functions
        wfcoul(in,ik,ir,2)=0d0
        wfcoulp(in,ik,ir,1)=0d0        !change for whittaker functions
        wfcoulp(in,ik,ir,2)=0d0

                                              !we also erase the bound part of the wavefunction
!        wfeig(iord(ie),in,ir)=0d0

        else
! AMoro
!        call coulfg(x,eta,0d0,1d0*lmax,f,g,fp,gp,1,0,ifail,m1)  !This was giving a problem
        CALL COUL90(x,eta,0d0,Lmax,f,g,FP,GP,0,IFAIL)
        IF (IFAIL.NE.0) THEN 
        WRITE(*,*) 'phase: IFAIL=',IFAIL
        ENDIF
		 
			 
			 
!        wfcoul(in,iord(ie),ir,1)=f(ql(in)+1)       
!        wfcoul(in,iord(ie),ir,2)=g(ql(in)+1)
!        wfcoulp(in,iord(ie),ir,1)=fp(ql(in)+1)    
!        wfcoulp(in,iord(ie),ir,2)=gp(ql(in)+1)
C commented by AMoro (2.0.5c)  
!        wfcoul(in,ik,ir,1)=f(ql(in)+1)       
!        wfcoul(in,ik,ir,2)=g(ql(in)+1)
!        wfcoulp(in,ik,ir,1)=fp(ql(in)+1)    
!        wfcoulp(in,ik,ir,2)=gp(ql(in)+1)   

        wfcoul(in,ik,ir,1)=f(ql(in))       
        wfcoul(in,ik,ir,2)=g(ql(in))
        wfcoulp(in,ik,ir,1)=fp(ql(in))    
        wfcoulp(in,ik,ir,2)=gp(ql(in))   

 
        endif
!        if ((in.eq.2).and.(ie.eq.45)) then
!        if (ir.eq.1) then
!        write(350,*)'#m1=',m1
!        else
!        write(350,'(6f12.4)')rvec(ir),wfcoul(in,iord(ie),ir,1),
!     &wfcoul(in,iord(ie),ir,2),wfcoulp(in,iord(ie),ir,1),
!     &wfcoulp(in,iord(ie),ir,2),x         
!        endif
!        endif
       
        enddo
        enddo
        enddo
        end subroutine storecoulwf
        
        subroutine phase     !not working
        use globals
        use sistema
        use constants
        use wfs !remember to set an only:... to reduce memory
        use channels
        use coulwf
!        use potentials, only:r0,at,ap
        use hmatrix, only:ortmat
        use scattering
        implicit none
        real*8:: faux,a,r,hnm,res,k
        real*8,allocatable :: gaux(:)
        integer:: ir,n,m,ie,p,q,ik
        complex*16,allocatable:: numh(:,:),nume(:,:)
        complex*16,allocatable:: denh(:,:),dene(:,:)

        complex*16,allocatable:: den(:,:,:),num(:,:,:),s(:,:,:)
        real*8,allocatable :: ia(:,:),ra(:,:),delta(:,:,:)
        real*8,allocatable,target:: wfcont2(:,:,:,:)

        real*8:: eps,test
        eps=1d-6

        allocate(delta(nchan,nchan,nk),gaux(nr))
        allocate(numh(nchan,nchan),nume(nchan,nchan))
        allocate(denh(nchan,nchan),dene(nchan,nchan))
        allocate(den(nchan,nchan,nk),num(nchan,nchan,nk))
        allocate(ia(nchan,nchan),ra(nchan,nchan))
        allocate(s(nchan,nchan,nk))
        allocate(wfcont2(nk,nchan,nchan,nr))
        wfcont2(:,:,:,:)=real(wfcont(:,:,:,:))

        write(*,*)' - Starting phaseshifts calculation'

!        a=1/r0/(at**(1./3)+ap**(1./3))
!        write(*,*)'a=',a,'r=',r0*(at**(1/3)+ap**(1/3))
!        write(*,*)r0,at,ap
        a=0.0001
        call storecoulwf()
        do ir=1,nr
        r=rvec(ir)
        if (rvec(ir).lt.eps) r=eps
        wfcoul(:,:,ir,:)=wfcoul(:,:,ir,:)*faux(a,r)
!        wfcoulp(:,:,ir,:)=wfcoulp(:,:,ir,:)*faux(a,r) !this is not true !in principle is not needed
!        write(351,'(8f12.4)')r,(wfcoul(n,4,ir,1),n=1,nchan),faux(a,r)
        enddo

        !Divided into two parts: <phi|H|f> and E<phi|f>

!        write(*,*)' - Starting phaseshifts calculation'
       
        den(:,:,:)=(0d0,0d0)
        num(:,:,:)=(0d0,0d0)
        s(:,:,:)=(0d0,0d0)




        do ik=1,nk
        k=kmin+(ik-1)*(kmax-kmin)/(nk-1)

        nume(:,:)=(0d0,0d0)
        dene(:,:)=(0d0,0d0)
        numh(:,:)=(0d0,0d0)
        denh(:,:)=(0d0,0d0)

        do ie=1,hdim
        if (ebin(ie).lt.0d0) cycle

        wfaux=>wfcoul(:,ik,:,1)  !first we calculate <wf_k|E-H|F_m>

        wfaux2=>wfeig(ie,:,:)
!        wfaux=>wfeig(iord(ie),:,:)
         do n=1,nchan

          do p=1,nchan
!           k=k-sqrt(2*mu12/(hc**2)*exc(p)) !not necesary: E= E_tot= hc*k**2/2m+ exc(n)
           do m=1,nchan
!          gaux(:)=wfeig(iord(ie),m,:)*wfcoul(n,iord(ie),:,1)*
!     &rvec(:)*rvec(:)
!          call sim(gaux,res,1,nr,dr,nr)
        numh(p,n)=numh(p,n)+solapmat(ik,n,ie,m)*hnm(m,m,p,p)
          if (p.eq.m) then
        nume(p,n)=nume(p,n)+
     &(hc*k)**2/2*mu12*solapmat(ik,n,ie,m)*ortmat(m,p)
          endif
!        write(*,'(8f12.4)')numh(n,m),nume(n,m),res*ebin(iord(ie)),res
          enddo
        num(p,n,ik)=nume(p,n)-numh(p,n)                   !we sum for <wf|E-H|F_m>
!          if (m.eq.1) then
!          write(*,*)'num=',num(1),nume(1,1)-numh(1,1)
!          endif
          enddo
         enddo

!         write(*,*)'numerator calculated'

!        write(*,*)

        wfaux=>wfcoul(:,ik,:,2)  !second we calculate <wf_k|E-H|G_n>
!        wfaux2=>wfeig(iord(ie),:,:)
         do n=1,nchan
          do p=1,nchan
!           k=k-sqrt(2*mu12/(hc**2)*exc(p))
           do m=1,nchan
        denh(p,n)=denh(p,n)+solapmat(ik,n,ie,m)*hnm(m,m,p,p)
!           write(*,*)'hnms=',hnm(m,m,p,p)
!           write(*,*)'solap=',solapmat(ik,n,iord(ie),m)
          if (p.eq.m) then
        dene(p,n)=dene(p,n)+
     &(hc*k)**2/2*mu12*solapmat(ik,n,ie,m)*ortmat(m,p)
          endif
          enddo
        den(p,n,ik)=dene(p,n)-denh(p,n)   !we sum for <wf|E-H|G_n>
!          if (n.eq.1) then
!          write(*,*)'den=',den(1),dene(1,1)-denh(1,1)
!          endif
          enddo
         enddo



!        delta(iord(ie),:,:)=atan(-(-numh(:,:)+nume(:,:))/
!     &(-denh(:,:)+dene(:,:)))*180/pi

!        write(352,'(8f12.4)')ebin(iord(ie)),
!     &(delta(iord(ie),1,m),m=1,nchan)


!        do n=1,nchan
!        do m=1,nchan
!        delta(iord(ie),n,m)=atan(-num(m)/den(n))*180/pi
!        enddo
!        enddo


!        write(350,'(8f12.4)')ebin(iord(ie)),
!     &(delta(iord(ie),1,m),m=1,nchan)

!        write(351,'(8f12.4)')ebin(iord(ie)),
!     &(delta(iord(ie),n,1),n=1,nchan)



!!!                         CHECK   <wf|E-H|wf>=0d0


!        test=0d0
!        den(:)=0d0
!        wfaux=>wfeig(iord(ie),:,:)
!         do n=1,nchan
!          do k=1,nchan
!          denh(k,n)=hnm(n,n,k,k)
!          if (k.eq.n) then
!          dene(k,n)=ebin(iord(ie))*ortmat(n,k)
!          endif
!          den(n)=den(n)+dene(k,n)-denh(k,n)   !we sum for <wf|E-H|G_n>
!!          if (n.eq.1) then
!!          write(*,*)'den=',den(1),dene(1,1)-denh(1,1)
!!          endif
!          enddo
!          test=test+den(n)
!         enddo
!
!          write(*,*)'test=',test

        enddo

!        write(*,*)'asignando real e imaginaria de A para A^-1'
        ra(:,:)=0d0
        ia(:,:)=0d0
        ra(:,:)=dreal(den(:,:,ik))
        ia(:,:)=aimag(den(:,:,ik))
!        write(*,*)ra,den(:,:,ik),den(1,1,40)
        if (ra(1,1).eq.0d0) then
        write(*,*)'ik=',ik,' skiped'
        delta(:,:,ik)=0d0
        else
!        write(*,*)'inverting A'
        call cmatin(ra,ia,nchan)
!        write(*,*)'done'
        den(:,:,ik)=cmplx(ra(:,:),ia(:,:))
!       we multiply A^(-1)*B
        do n=1,nchan
        do m=1,nchan
          do p=1,nchan
!        write(*,*)'here we go'
        s(n,m,ik)=s(n,m,ik)-den(n,p,ik)*num(p,m,ik)   ! is this k?
          enddo
        s(n,m,ik)=(cmplx(1d0,0d0)+cmplx(0d0,1d0)*s(n,m,ik))  ! to go from k to S
     &/(cmplx(1d0,0d0)-cmplx(0d0,1d0)*s(n,m,ik))
        enddo
        enddo
        
       delta(:,:,ik)=datan2(aimag(s(:,:,ik)),dreal(s(:,:,ik)))
       write(350,'(8f12.4)')(hc*k)**2/(2*mu12),(delta(1,m,ik),m=1,nchan)

        endif
c----------------------------------------------------------------------------------------------------------------------
!       directly from wave functions
        nume(:,:)=(0d0,0d0)
        dene(:,:)=(0d0,0d0)
        numh(:,:)=(0d0,0d0)
        denh(:,:)=(0d0,0d0)


        wfaux=>wfcoul(:,ik,:,1)  !first we calculate <wf_k|E-H|F_m>

         do n=1,nchan
        wfaux2=>wfcont2(ik,n,:,:)
          do p=1,nchan
        numh(p,n)=numh(p,n)+hnm(n,n,p,p)
          if (p.eq.n) then
        nume(p,n)=nume(p,n)+
     &(hc*k)**2/2*mu12*ortmat(n,p)
          endif

        num(p,n,ik)=nume(p,n)-numh(p,n)                   !we sum for <wf|E-H|F_m>
!          if (m.eq.1) then
!          write(*,*)'num=',num(1),nume(1,1)-numh(1,1)
!          endif
          enddo
         enddo

!         write(*,*)'numerator calculated'

!        write(*,*)

        wfaux=>wfcoul(:,ik,:,2)  !second we calculate <wf_k|E-H|G_n>
         do n=1,nchan
        wfaux2=>wfcont2(ik,n,:,:)
          do p=1,nchan

        denh(p,n)=denh(p,n)+hnm(n,n,p,p)
          if (p.eq.n) then
        dene(p,n)=dene(p,n)+
     &(hc*k)**2/2*mu12*ortmat(n,p)
          endif
        den(p,n,ik)=dene(p,n)-denh(p,n)   !we sum for <wf|E-H|G_n>
          enddo
         enddo

        ra(:,:)=0d0
        ia(:,:)=0d0
        ra(:,:)=dreal(den(:,:,ik))
        ia(:,:)=aimag(den(:,:,ik))

        if (ra(1,1).eq.0d0) then
        write(*,*)'ik=',ik,' skiped'
        delta(:,:,ik)=0d0
        else

        call cmatin(ra,ia,nchan)

        den(:,:,ik)=cmplx(ra(:,:),ia(:,:))

        do n=1,nchan
        do m=1,nchan
          do p=1,nchan

        s(n,m,ik)=s(n,m,ik)-den(n,p,ik)*num(p,m,ik)   ! is this k?
          enddo
!        s(n,m,ik)=(cmplx(1d0,0d0)+cmplx(0d0,1d0)*s(n,m,ik))  ! to go from k to S
!     &/(cmplx(1d0,0d0)-cmplx(0d0,1d0)*s(n,m,ik))
        enddo
        enddo
        
       delta(:,:,ik)=datan2(aimag(s(:,:,ik)),dreal(s(:,:,ik)))
       write(351,'(8f12.4)')(hc*k)**2/(2*mu12),(delta(1,m,ik),m=1,nchan)

        endif
        enddo        
        end subroutine phase

        function faux(a,r)  ! f(r)  !not used
         real*8:: faux
         real*8:: r ,a
         faux=1d0-dexp(-a*r**2)
        end function faux


!	Subroutine for calculating the cross section for capture reactions by electrical transitions.

      subroutine capture(jtot,partot)
      use constants, only:hc,e2,pi !finec=137.03599d0, hc=197.32705, e2=hc/finec, amu=931.49432
      use globals, only:mu12,egs
      use scattering, only: nk,kmin,kmax,emin,emax !il
      use sistema, only: zv,zc
      use channels, only: sn
!      use factorials
      use belambdamod, only: lambda,partoti,jtoti,eneri,dbde
      implicit none
!      real*8,allocatable,dimension(:):: crsect,kvec !dimension nk
      integer partot
      real*8  jtot
      real*8,dimension(nk):: crsect,kvec,kphot,ephot,enervec
      real*8,dimension(nk):: sommerf,se
      real*8 step
      real*8 doblefact
      integer i,il
!      WRITE(*,*)'AVER-1'
      
      if (partoti*(-1)**lambda.ne.partot) then
      print*
      write(*,*)'ERROR: parity violation.'
      write(*,*)'Parity in electrical transitions must change as'
      write(*,*)' (-1)**lambda.'
      print*
      crsect=0d0
      se=0d0
      goto 333
      endif
!			ANGULAR MOMENTUM RESTRICTION IS MISSING !!!

!      WRITE(*,*)'AVER0'
      write(*,*)'Capture x-section from Jt=',Jtoti,' to ',Jtot
      step=(kmax-kmin)/(nk-1)
      kvec(1)=kmin
      do i=2,nk
      kvec(i)=kvec(i-1)+step
      enddo 
      enervec=0.5*hc*hc*kvec*kvec/mu12

      ephot=(enervec-egs)
      kphot=sqrt(2d0*mu12*abs(ephot))/hc

      crsect=(lambda+1d0)/lambda*(2d0*pi)**3
      crsect=crsect/(doblefact(2*lambda+1))**2
      crsect=crsect*(ephot/hc)**(2d0*lambda-1d0)

      il=1 ! I'm not sure about it... is 'il' stored in a module so its value is il=1?
!      write(*,*)'AVER1'
      crsect(:)=crsect(:)*dbde(il,:)*e2 !e2 factor in order to obtain correct dimension unities
      crsect=crsect*2d0*(2d0*jtot+1d0)/(2d0*jtoti+1d0)/(2d0*sn+1d0)
      crsect=crsect*kphot**2/kvec**2 !capture cross section by electrical transition

      sommerf=zv*zc*e2*mu12/(hc)**2/kvec !Sommerfeld parameter
      se=enervec*exp(2d0*pi*sommerf)*crsect !astrophysical factor

!      write(*,*)'AVER2'
  333 write(600,*)'# core1+nucleon --> core2+photon cross section'
      write(600,*)'# energy c.m. (MeV) | cross section (fm^2)'
      write(601,*)'# core1+nucleon --> core2+photon cross section'
      write(601,*)'# energy c.m. (MeV) | astroph. factor (MeV fm^2)'
      do i=1,nk
      write(600,*)enervec(i),crsect(i),dbde(il,i)
      enddo
      do i=1,nk
      write(601,*)enervec(i),se(i)
      enddo
      print*
      write(*,*)'Capture cross section saved in "fort.600"'
      write(*,*)'astrophysical factor saved in "fort.601"'
      print*

      

      return
      end subroutine


      real*8 function doblefact(x) ! If 'factorials' module is used in subroutine cross, this is repetitive.

      implicit none
      integer x,a
      a=x
  310 if ((a-2).ge.2d0) then
      a=a-2
      x=x*a
      goto 310
      endif

      doblefact=x*1d0

      return
      end function



