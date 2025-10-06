c THOx: A CDCC code with core and target excitations 
c
c by A.M. Moro, J.A. Lay, M.Gomez-Ramos, R. de Diego (2011-2016)
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
c v2.4 AMM: 3-body observables with CC bins implemented
c v2.6 AMM: calculation of core+valence eigenphases
      program thox
      use globals
      use sistema
      use parameters
      use constants
      use hmatrix
      use trace
!      use forbidden, only: hindrance
      use channels
      use wfs, only: nst
      use belambdamod, only: ifbel
      use memory
      use potentials, only: reid93
! added by JLei 
      use cdccchannels
      use writecdccwf
      implicit none
      logical ehat
      integer :: lmax,i,ncc,bastype,nk
      real*8 eps
     
      parameter (eps=1e-6)

      logical checkort,ifphase,dummy,tres
!      real*8 ::  lambda,norm
      CHARACTER*1 BLANK,PSIGN(3)
      character*40 filename
      integer :: nset,nho,nodes
      logical :: changepot

!!! TEST 
      integer iset,inc,iexgs,nchsp
!      real*8 ecm


c Input namelists -------------------------
      namelist/system/ Zc,Ac,Zv,Av,egs,sn
      namelist/output/ wfout,checkort,verb,ifphase,solapout,
     &                 cdccwf,froverlaps,
     &                 xcdcc

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
        real*8 lrf,lfv(:),alpha
        complex*16:: gaux(:)
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
     &             ncont,energy,wfcont,smate,psh)
        integer :: iset,nchan,inc,ncont
        logical :: energy
        complex*16 :: wfcont(:,:,:),smate(:,:)
        real*8 emin,emax,psh(:)
        end subroutine
      end interface
c     ------------------------------------------------------------

      debug=.false.

c *** Defined global constants
      call initialize()
!      if (cdccwf) call alpha_cdcc_in()
      write(*,'(50("*"))')
      write(*,'(" ***",8x,"THOx+DCE+CC code: version 2.6",8x, "***")')
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
      wfprint(:) =.false.
      ifbel      =.false.
      checkort   =.false.
      ifphase    =.false.
      wfout(:)   =0
      solapout(:)=0
      cdccwf     =.false.
      froverlaps =0
      
      read(*,nml=output)
      do i=1,10
       if (wfout(i)>0) wfprint(wfout(i))=.true.
      enddo


      if (cdccwf) call alpha_cdcc_in()
      if (cdccwf) call alpha3b_in()
      if (cdccwf) call abar_cdcc()


c *** Radial grid (quadratures?)
      call readgrid(kin)

c *** Read potential parameters
      call readpot(kin)

c Pauli-forbidden states to be removed
      call pauli_forbidden(kin)

      nst=0 ! total number of selected eigenvalues
      do iset=1,inpsets
      dummy=.false.
      nset=iset
!      nset=indjset(iset)
!      call read_jpiset(iset,bastype,nk,tres,ehat,filename)
      call read_jpiset(nset,bastype,nk,tres,ehat,filename,nodes,
     & changepot)

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
      case(0,1,3,7) ! HO, THO cTHO -----------------------------------------
c *** Build and diagonalize s.p. Hamiltonian for each j/pi set
        if (bas2.eq.0) call hsp(nset,nchsp,nho)

c *** Build and diagonalize FULL VALENCE+CORE Hamiltonian 
        call fullham(nset,nho)
        
c *** Vertex functions
      if (froverlaps>0)  call vertex(nset,froverlaps)        

c *** Check orthonormality
        if (checkort) call orthonorm ! ham.f90

c *** Overlap between THO and scattering wfs
        call solap
c *** ------------------------------------------------------------
      case(2,4) ! Bins   
        call makebins(nset,nk,tres,ehat)

c *** Vertex functions
      if (froverlaps>0)  call vertex(nset,froverlaps)
	  

c *** ------------------------------------------------------------
      case(5) ! External   
       call readwfs(filename,nset,nchan)

      case(6) ! Eigcc   
       call pre_eigcc(nset,nchan,nodes,changepot)

      end select 
      
      enddo ! loop on j/pi sets
      
      call write_states()

c *** Continuum WFS (scatwf namelist)
!      call continuum!(nchan)
      call continuum_range


c *** B(E;lambda)      
      call belam ! (kin,jtot,partot)
      
c *** B(M;lambda)      
      call bmlam ! (kin,jtot,partot)      

c *** Phase-shifts
!      if (ifphase) call phase


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! TEST calculation of scattering wf for arbitrary j/pi and energy
!      iset=1
!      ecm= 1.
!      inc=1
!      call test_cont(iset,jpiset(iset)%nchan,inc,ecm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C *** Capture x-sections and S-factors (unfinished)
!!!      if (ifbel) call capture(jtot,partot)


!Make sure Reid93 potential is not used for target potentials
      

c *** Define reaction (projectile, target, energy)
      call partition()

c *** Coupling potentials (includes DCE routines)
      reid93=.false.
      if (.not. targdef) then
      write(*,'(5x,a)') 
     & '[using transition subroutine for spin-zero target]'
      call transition 
      else
      write(*,'(5x,a)') '[using transition_targdef subroutine 
     & with target excitation]'
      call transition_targdef
      endif

c *** Build & solve CDCC equations
      iexgs=1 !assume that the incident channel is the first state!!!!
      ncc  =0
      call recoupling()
      

      call cdcc(iexgs,ncc)
      if(cdccwf) call cdcc_wf_smoothie_out()


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
      character*5:: jpi
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

      write(*,'(5x,"Nb. of target states:",i3)') ntex+1
      do itex=1,ntex+1
      if (tset(itex)%partarg.eq.1) pchar='+'
      if (tset(itex)%partarg.eq.-1) pchar='-'
      write(*,'(5x,"o state #",i3,3x, "J/pi=",a5,3x, "Ex=",1f6.3)')
     &  itex, 
     & jpi(tset(itex)%jtarg,tset(itex)%partarg),tset(itex)%extarg
!      write(*,*) 'State #',itex,': J/pi',tset(itex)%jtarg,pchar,
!     & ' E:',tset(itex)%extarg
      enddo

      end subroutine


      subroutine timings()
      use memory
      write(*,*)''
      write(*,*)'Timings:'
      if (tcc.gt.1e-4) write(*,*)' CC integration :',tcc,  ' secs'
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
      use parameters, only: maxl,maxeset,maxchan
      use xcdcc, only: realwf
      use wfs, only: energ,wfc,nr
!      use potentials, only: vl0,vcp0
      implicit none
      logical :: tres,ehat,merge
      integer::l,lmin,lmax,bastype,mlst,kin,ng
      integer:: basold,parold,incold 
      real*8:: exmin,exmax,j,jtold,rmin,rmax,dr,rlast,rint
      real*8:: bosc,gamma,kband,wcut(1:maxchan),vscale,r1,rnmax
      integer:: ichsp,ic,iset,nchsp,nfmax,nho,parity
      integer:: nk,nbins,inc,nodes
      logical changepot
      character*40 filewf

      namelist /grid/ ng, rmin,rmax,dr,rlast,rint
            
      namelist/jpset/ jtot,parity,l,j,lmin,lmax,
     &                bastype,nfmax,exmin,exmax,
     &                nho,bosc,     ! HO
     &                gamma,mlst,   !THO Amos 
     &                nsp,  ! sp eigenvalues to keep for full diag
     &                bas2,
     &                nk, nbins,inc,tres,ehat,filewf,wcut,merge,
     &                vscale,
     &                r1,rnmax,nodes,changepot ! CG Basis

      
      jpsets=0; inpsets=0
      maxl=0
      realwf=.true.
      merge=.false.
      r1=0.0
      rnmax=0.0

      write(*,*)'Pre-reading input to set dimensions'

      read(kin,nml=grid)
      if (dr>0) then
           nr=ceiling((rmax-rmin)/dr)+1
      endif
      
300   bastype=-1    
      tres=.false.
      inc=1
      l=-1 ; j=-1; lmin=-1; lmax=-1;
      vscale=1;
      r1=0.0; rnmax=0.0
      read(kin,nml=jpset) 
      if (bastype.lt.0) goto 350
      if (l.lt.0) l=0;
      lmin=l;
      if (lmax.lt.lmin) lmax=lmin
!      if (bastype.eq.2) realwf=.false. ! complex bins
      inpsets=inpsets+1
      if (inpsets.eq.1) then
         jtold=jtot
         parold =parity
         basold=bastype
         incold=inc
         indjset(inpsets)=1
         jpsets=1
      else
         if ((jtold.eq.jtot).and.(parity.eq.parold).and.
     x     (bastype.eq.basold).and.(inc.eq.incold)) then
            write(0,*)'Set',jpsets,' could be merged with previous one!'
            jpsets=jpsets+1 ! temporary
         else
           jpsets=jpsets+1
           jtold=jtot
           parold=parity
           basold=bastype
           incold=inc
         endif
         indjset(inpsets)=jpsets
!         print*,'jpsets,inpsets=',jpsets,inpsets
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
      write(*,*)'   o J/pi sets:', inpsets,' merged into', jpsets
      write(*,*)'   o Maximum internal l:', maxl

      write(*,*)''


      allocate(jpiset(jpsets)) ! sets, chans/set
      allocate(spchan(jpsets,maxchan))

      if (.not.allocated(wfc)) then
        nchmax=0
        allocate(wfc(jpsets,maxeset,maxchan,nr))
        allocate(energ(jpsets,maxeset))
      endif

      write(*,*)' Pre-read finished'
!      close(kin)
      rewind(kin, err=360)
      return

360   write(*,*)'error rewindind std input file'
      return

      end subroutine


c ***
c *** Read core states
c ***
      subroutine readcorex(kin)
       use channels
       implicit none
       integer:: kin
       integer:: nph,parity
       real*8 :: spin,ex,kband
       character*1 BLANK,PSIGN(3)      
       DATA PSIGN / '-','?','+' /, BLANK / ' ' /

       namelist /corestates/ spin,parity,ex,nph,kband


       nce=0 ! number of core states
       write(*,*) '- CORE states:'
50     ex=-100.
       parity=0
       nph=0
       kband=0
       read(kin,nml=corestates) 
       if (ex<-99) goto 100
       nce=nce+1   
       if (parity.eq.0) parity=+1
       if (abs(spin-nint(spin)).gt.0.2d0 .and. kband.eq.0) then
       write(*,*) 'Changing kband=0 to kband=0.5'
       kband=0.5d0
       endif
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
       qnc(nce)%kband= kband

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
     &           1f6.1,1x," Step=",1f6.3, " fm ]",/)
         if ((rint>0).and.(rint.gt.rmax)) then
          write(*,172)rint, rmax
172       format(1x,"[ WFS calculated up to",1f6.1,
     & " fm and extrapolated up to ", 1f6.1," fm ]",//)
         endif
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
       character*40 flkind(600)
       integer writf(399),nwrit
       logical written(600)
       flkind(:) = ' '
       flkind(20)='Val-core central potential for each l'
       flkind(22)='Val-core coupling potential for lambda=2'
       flkind(23)='Val-core coupling potential for lambda=3'
       flkind(40)='Valence-core potentials'
       flkind(44)='Fragment-target potentials'
       flkind(45)='Phase-shifts'
       flkind(46)='Eigenphases'
       flkind(50)='Full list of Cont. Wf.'
       flkind(51)='Continuum Wavefunction'
       flkind(52)='Overlaps between PS and scattering states'
       flkind(53)='Overlaps between PS and scat.stat. (no head)'

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

       flkind(500)='Continuum WFS from SCATWF namelist'

	nwrit = 0
	do i=1,600
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



c     -----------------------------------------------------------
c     Write states & energies in Fresco format
c     ----------------------------------------------------------- 
       subroutine write_states()
       use wfs, only: energ,nst
       use channels, only: jpsets,jpiset
       use sistema
!      use hmatrix
       use globals
!      use trace
      implicit none
      integer, parameter:: kst=25
      integer i,ie,icopyt,ibandp,ibandt,icpot,n
      real*8:: et, jt,rjp,ep       
      open(kst,file="states2.fr",status='unknown')     
      i=0
      icopyt=0
      icpot=1
      ibandt=0
      et=0.d0
      jt=jtgs ! inherited from module 'sistema'
!      write(0,*)'energ(1,1)=',energ(1,1)
      write(kst,*) nst
      do n=1,jpsets
      rjp    =jpiset(n)%jtot    ! xjp(n)
      ibandp =jpiset(n)%partot  ! ipar(n)
      do ie=1,jpiset(n)%nex    ! np(n)
      i=i+1
      ep=energ(n,ie)-energ(1,1)
      if (i.eq.2) icopyt=1
      write(kst,11511) rjp,ibandp,ep,icpot
      if(icopyt==0) write(kst,1153) jt,ibandt,et
      if(icopyt.ne.0) write(kst,1154) icopyt
      enddo !ie
      enddo !n

11511 format(' &States jp=',f4.1,' ptyp=',i2,' ep=',f8.4,'  cpot=',i3)
 1153 format('         jt=',f4.1,' ptyt=',i2,' et=',f8.4,' /')
 1154 format('         copyt=',i2,' /')
 1155 format(' Skipping state #',i3,' J/pi',f4.1,i2,' and Ex=',f8.4)
c     ----------------------------------------------------------------
 
       close(kst)
       end subroutine 

