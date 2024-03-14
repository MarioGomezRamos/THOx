c  DCE code 
c  Calculates coupling potentials for 3-body CDCC calculations with core excitation
c  according to the formalism first proposed in PRC74, 014606 (2006) 
c
c  R. de Diego, A.M. Moro (2013-2014)
c  ----------------------------------------------------------------------------------

      subroutine transition
      use xcdcc
      use sistema
      use wfs, only:energ,idx,rvec
      use ptpots
      use globals
      use channels
      use memory
      use potentials, only: vcl
      implicit real*8 (a-h,o-z)
      logical skip,writeff
      integer qcmin,qcmax,kcmax,coups!,lambmax
      integer lgs
c     -------------------------------------------------------------------------------
      character*5:: jpi
      character*3:: jname
      CHARACTER PARITY(3)
      DATA PARITY / '-','?','+' / 
      character*40 filename,couplings,comment,fileampsd
c     -------------------------------------------------------------------------------
      real*8 factor,ptr,ttr,jt,iin,ifi,kband1,kband2
      real*8,allocatable:: xrad2(:),xrad3(:)
      real*8 :: qfactorc(0:50),qfactorn(0:50)
      
c changed to complex in v2.3
      complex*16,pointer:: ui(:),uf(:)
      complex*16 :: fprod
      complex*16 xi,pot,xsum, xsumn,xsumc
!      parameter(Kmax=6,xKrot=0.d0,nkmax=300)
      parameter(xKrot=0.d0,nkmax=300)
      
c  For adiabatic JT potential
      complex*16,allocatable:: xintad(:),vadnorm(:)
      complex*16 :: fauxad
      complex*16 :: xsumad
      real*8     ::xsumnad
      real*8,allocatable:: vclquad(:)
      
c     ------------------------------------------------------------------------------

! v0.6c
!      dimension energ(nener,numfilmax),!,xrad2(:),elapsed(2),
!     .np(numfilmax)
c      real*8, allocatable:: energ(:,:)
c      integer, allocatable:: np(:)
      TYPE GFAC  ! derived variable to store the geometrical factors
       INTEGER i,f,k,l,nq
       REAL*8:: rmatc,rmatn
      END TYPE
      TYPE(GFAC) PK(nkmax)
      interface
       subroutine potdef(vc,qmax,rnc,vdef,betar)
       integer::qmax
       real*8:: vc(:),rnc,vdef(:),betar
      end subroutine
      end interface

      interface
        subroutine extpot(filename,vr,vi,nr)
        integer::nr
        character*40 filename
        real*8 :: vr(:),vi(:)
       end subroutine
      end interface

! Modified by AMM in v2.2
!      namelist /trans/ mp,mt,zt,zp,rcc,namep,namet
      namelist /trans/ skip,rcc,writeff
      namelist /grid/ nquad,radmax,rmax,rextrap,rstep,rvecin,drvec,hin
      namelist /wf/ filename,eminf,emaxf
      namelist /coupling/ qcmin,qcmax,kcmax,lamin,lamax,coups,ncoul,
     &          qfactorn,qfactorc
      namelist /BEQ/ bel
!      interface
!      subroutine gauleg(x1,x2,x,w,n)
!      integer::n
!      real*8:: x1,x2,x(:),w(:)
!      end subroutine
!      end interface
      call cpu_time(start)
      debug=.false.
      xi=(0.d0,1.d0)
      kin=5
      kread=84
      kfr=4
      ken=8
      pi=acos(-1d0)
      writeff=.false.
      qfactorc(:)=1.0; qfactorn(:)=1.0; 
c     ----------------------------------------------------------------------------      
      write(*,'(/,5x,"*********** COUPLING POTENTIALS  *********** ")')
      
c ---------------------------------------------------------------------------------
c AMoro: partial preread to set dimensions
      rextrap=0
      skip=.false. 
      read(kin,nml=trans)

      if(jpsets.eq.0) then
      write(*,'(/,5x,"** No basis states have been defined",
     &    "=> skipping this section **" )') 
      return
      endif

      if (skip) then
      write(*,'(5x,"[No coupling potentials requested ]")')
      return
      endif

c *** Quantum number limits for transition potentials
      lamin=0
      lamax=-1
      qcmin=0
      coups=0
      qcmax=-1
      iftrans=.true.
      kcmax=-1
      read(kin,nml=coupling)
      if (qcmax.lt.0) then
       write(*,*)'[=> No coupling potentials requested]'
       iftrans=.false.
       return
      endif
!      write(*,'(/,2x,"o Couplings considered for formfactors:")')
      select case(coups)
      case (0)
       write(*,*) ' [ coups=0=> All transitions will be considered]'
      case (1)
       write(*,*) ' [ coups=1=> ONLY gs->cont couplings considered]'
      case (2)
       write(*,*) ' [ coups=2=> Diag nucl. + ALL Coulomb]'
       if (ncoul.ne.0) write(*,*) ' ( ncoul will be ignored )'
      case (3)
       write(*,*) ' [ coups=3=> V00 + gs->cont]'
   
      case (4)
       write(*,*) ' [ coups=4=> gs->cont (no diag) ]'
   
      end select
      
      select case(ncoul)
      case (0)
        write(*,*)' [ ncoul=0 => C+N couplings ]'
      case (1)
        write(*,*)' [ ncoul=1 => NUCLEAR couplings only ]'
      case (2)
        write(*,*)' [ ncoul=2 => COULOMB couplings only]'
      end select

!	  if (any(qfactorc(1:50)).ne.1.0)
!      write(*,*)"qfactorc=",qfactorc(0:5)
!      write(*,*)"qfactorc=",qfactorn(0:5)

c *** Radial grids ------------------------------------------------------------
      read(kin,nml=grid)
      coefc=av/mp
      coefv=-ac/mp ! = coefc-1
      coef=max(coefc,dabs(coefv))
      if (rstep.gt.1e-5) then 
      nrad2=nint((rmax-rstep)/rstep)+1
      if (rextrap.lt.rmax) rextrap=rmax
      nrad3=nint((rextrap-rstep)/rstep)+1
      allocate(xrad2(nrad2))
      allocate(xrad3(nrad3))
      else
      write(*,*)'radial step (rstep) too small!'; stop
      endif
      allocate(xquad(nquad),wquad(nquad),rquad(nquad))
      allocate(xintgn(nquad),xintgc(nquad)) ! AMoro
      allocate(xintad(nquad),vclquad(nquad),vadnorm(nquad)) ! adiabatic
      a=-1.d0
      b=1.d0
      call gauleg(a,b,xquad,wquad,nquad)

      rquad=radmax*0.5d0*(1.d0+xquad)
      rvecmax=coef*radmax+rmax
      nr=nint((rvecmax-rvecin)/drvec)

      allocate(rfrag(nr))
      allocate(rvin(nrin))
      do irvec=1,nr
      rfrag(irvec)=rvecin+drvec*irvec
      enddo
      do irvec=1,nrin
      rvin(irvec)=hin*irvec
      enddo

      if(.not.allocated(vcore))  allocate(vcore(1:nr,0:nmult))
      if(.not.allocated(vcorei)) allocate(vcorei(1:nr,0:nmult))
      if(.not.allocated(vcorec)) allocate(vcorec(1:nr,0:nmult))
      if(.not.allocated(vval))   allocate(vval(1:nr,0:nmult))
      if(.not.allocated(vvali))  allocate(vvali(1:nr,0:nmult))
      if(.not.allocated(vvalc))  allocate(vvalc(1:nr,0:nmult))
      kptype=0
      vfrag => vcore
      vfragi=> vcorei
      vcoup => vcorec
      write(*,'(/,2x,"** CORE-TARGET potential:")')
      call read_fragpot(kin,kptype,zc)
      melc(1:nmult)=mel(1:nmult)

      nullify(vcoup)
      kptype=1
      vfrag =>vval
      vfragi=>vvali
      vcoup =>vvalc
      write(*,'(/,2x,"** VALENCE-TARGET potential:")')
      call read_fragpot(kin,kptype,zv)
      melv(1:nmult)=mel(1:nmult)
      do irad2=1,nrad2
      xrad2(irad2)=rmax/dble(nrad2)*irad2
      enddo

c *** Interpolate projectile wfs at quadrature points
      call wf2quad()
      
c *** Interpolate projectile v+c potential at quadrature points
! Assuming single-channel and only central potential!!!!!!!!!!!!!!!!!!!
      alpha=0.0
      do iq=1,nquad
      raux=rquad(iq)
      lgs=jpiset(1)%lsp(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GENERALIZE FOR MULTICHANNEL!!!!!!!!!!!!
!      write(0,*) lgs,nr,vcl(lgs,1), size(vcl(lgs,:))
      faux=fival(raux,rvec,vcl(lgs,:),size(rvec),alpha) ! R(r)
!      faux=FFC4((raux-rmin)/dr,yvp,nr)
      vclquad(iq)=faux 
      write(95,'(1f8.3,2x,100g14.6)')raux,faux !,frad(jset,ie,iq,ich)
      enddo ! iq=quadrature point      


c *** Set limits for quantum numbers & multipoles ----------------------------------
      npe    =maxval(jpiset(:)%nex)
      xImax  =maxval(spchan(:,:)%jc)
      xJpmax =maxval(jpiset(:)%jtot)
      nIcmax =nint(2.d0*xImax)
!      nQmax=min(nmult,nIcmax) ! RDD

!!!! AMM: Iin, Iff, need to be checked
      nQmin=max(0,qcmin)
      nQmin=qcmin
      nQmax=min(qcmax,nIcmax)
      if (lamax.ge.0) then
          lambmax=min(lamax,nint(2.d0*xJpmax))
      else
          lambmax=nint(2.d0*xJpmax)
      endif
      
!! New in Mar/17
      lambmin=max(nint(dabs(xjp1-xjp2)),lamin)


!!! NEW IN v2.2g (CHECK!!)
      kmax=nqmax+lambmax
      if (kcmax.ge.0) kmax=min(kmax,kcmax)  

!      write(*,*)'qcmin,nqmin,Iin,Ifi=',qcmin,nqmin,Iin,Ifi
!      write(*,400) nqmin,nqmax,0,min(kmax,kcmax),0,lambmax
      write(*,400) nqmin,nqmax,0,kmax,lambmin,lambmax

400   format(/,2x,"o Quantum numbers limits:",/,5x,
     &  "- Core multipoles:",i2,"-",i1,/,5x,
     &  "- K:",i2,"-",i2,/,5x,
     &  "- Proj.-target multipoles:", i2,"-" i1)
c -------------------------------------------------------------------------------------


c *** START CALCULATION OF COUPLING POTENTIALS ---------------------------------------
      write(*,'(/,2x,"** FORMFACTORS **")') 
      fmem=nex*nex*(lambmax+1)*nrad3*lc16/1e6
      write(*,'(5x," [ formfactors need",1f7.1," Mbytes ]")') fmem

      allocate(potQKn(nquad,nrad2,0:nQmax,0:Kmax))
      allocate(potQKc(nquad,nrad2,0:nQmax,0:Kmax))
      allocate(Ff(nex,nex,0:lambmax,nrad3))
c commented by AMoro, to save memory
!      allocate(Fc(nex,nex,0:lambmax,nrad3))
!      allocate(Fn(nex,nex,0:lambmax,nrad3)) 
!      Fc=0.d0;  Fn=0.d0
      Ff=0.d0
      potQKn=0d0; potQKc=0d0
      rmat=0d0



       
      write(*,'(/,2x,"o V^{K,Q}(r,R) ",$)')
      call cpu_time(t1)
      call factorial(2*nQmax)

      do iquad=1,nquad
      do irad2=1,nrad2
      do nq=nqmin,nqmax
      do k=0,kmax
! AMoro: 
! ncoul=0: coulomb + nuclear
! ncoul=1: nuclear only
! ncoul=2: Coulomb only  
      potn=0d0; potc=0d0
      r1=rquad(iquad); ! internal coordinate (r)
      r2=xrad2(irad2)  ! p-t coordinate (R)
      if ((ncoul.eq.0).or.(ncoul.eq.1)) then  !nuclear part
        nc=1 ! nuclear part
        if(nq.eq.0) then                      !monopole
        potn= pot(r1,r2,nq,k,coefc,nc)
     &      + pot(r1,r2,nq,k,coefv,nc)
        else !      no core excitation for valence fragment! 
        potn=pot(r1,r2,nq,k,coefc,nc)
        endif
      endif
      if ((ncoul.eq.0).or.(ncoul.eq.2)) then ! Coulomb part
       nc=2 ! Coulomb part
       potc= pot(r1,r2,nq,k,coefc,nc)
     &     + pot(r1,r2,nq,k,coefv,nc)
      endif       

      potQKc(iquad,irad2,nq,k)=potc  
      potQKn(iquad,irad2,nq,k)=potn
      
      enddo !k 
      enddo !nq
      
      enddo ! nrad2 (R)
      enddo ! nquad (r)
      call cpu_time(t2)
      write(*,'(4x,"[",1f6.2," secs ]")') t2-t1
     
      rfirst=rstep 
      write(*,300) rfirst,rmax,rstep
300   format(/,5x,"[Radial grid: Rmin=",1f6.3," fm,", 
     &      " Rmax=",1f6.1," fm," 
     &      " Step=",1f6.3," fm]",/)

      if (rextrap.gt.rmax) then
       write(*,301) rextrap
301   format(5x,"[ Coulomb couplings extrapolated to",1f8.1," fm]",/)
      endif
      nff=0
      do m=1,jpsets
        nchann1=jpiset(m)%nchan ! nch(m)
        xjp1=jpiset(m)%jtot! xjp(m)
! TEst 11/OCT/16
!      do n=m,jpsets
      do n=1,jpsets
      if(realwf.and.n.lt.m) cycle
!!!!!!!!!!!!!!!!!!
        nchann2=jpiset(n)%nchan ! nch(n)
        xjp2=jpiset(n)%jtot     ! xjp(n)
      
!      do lc=nint(dabs(xjp1-xjp2)),min(lambmax,nint(xjp1+xjp2))
      do lc=lambmin,min(lambmax,nint(xjp1+xjp2))
        xlc=dble(lc)
        nk=0  !Nb. of P(IK) for each n,m,Lambda
        call cpu_time(t1)
        do i=1,nchann1
       xl1=jpiset(m)%lsp(i) 
       ic1=jpiset(m)%cindex(i)
       kband1=qnc(ic1)%kband
       Iin=jpiset(m)%jc(i) 

       do j=1,nchann2      
       xl2=jpiset(n)%lsp(j) 
       Ifi=jpiset(n)%jc(j) 
       ic2=jpiset(n)%cindex(j)
       kband2=qnc(ic2)%kband

! AMM: Assume rotor model and kband1=kband2
       if (kband1.ne.kband2) cycle
     
      do k=0,kmax
      do nq=max(nint(Iin-Ifi),nqmin),nqmax
      do l=0,nq ! small lambda
! CHECK AMM : i <-> f
!      call rmatel(m,n,xlc,i,j,k,nq,l,xjp1,xjp2,xKrot,rmatn,rmatc)
      call rmatel(n,m,xlc,j,i,k,nq,l,xjp2,xjp1,kband1,rmatn,rmatc)

      if (abs(rmatn).lt.1e-6.and.abs(rmatc).lt.1e-6) cycle
      nk=nk+1
      
      if (nk.gt.nkmax) then
        write(*,*)'NK > NKMAX=',NKMAX,' SO INCREASE NKMAX'
        stop
      endif
      pk(nk)%i     =i
      pk(nk)%f     =j
      pk(nk)%k     =k
      pk(nk)%nq    =nq  
      pk(nk)%l     =l    ! \lambda
      pk(nk)%rmatc =rmatc
      pk(nk)%rmatn =rmatn
      if ((m.eq.1).and.(n.eq.1).and.(verb.gt.0)) then
        if (nk.eq.1) then
        write(198,'("j=",1f4.2, " Jf=",1f4.2)')xjp1,xjp2
        write(198,'(7A7,2A10)')
     & "NK","Ch1","Ch2","LAM","K","QC","lambda","P(N)", "P(C)"
        endif
      endif
      if (verb.gt.0) write(198,1100) nk,i,j,lc,k,nq,l,
     & rmatn,rmatc,xl1,xl2
1100  format(7i7,2x,2f10.3,2x,2f4.1,2i5)
      enddo ! l (lambda)
      enddo ! nq (core multipole Q)
      enddo ! K
      enddo ! j (nchan)
      enddo ! i (nchan)


      if (nk.gt.0) then
      write(*,'(4x,a6,"-> ",a6," with LAM=",1i2,
     &   2x,"NK=",i3,2x, "non-zero P(K)  ",$)') 
     &   jpi(xjp1,jpiset(m)%partot),
     &   jpi(xjp2,jpiset(n)%partot),lc,nk
      else        
        if (debug) write(*,'(4x,"m=",i2,"-> n=",i2," with LAM=",1i2,
     &  2x,"(NO allowed transitions)")') m,n,lc
      endif


c ---------------------------------------------------------
c    Calculate formfactors F(R) for each |n,ie1> ->|m,ie2> 
c ---------------------------------------------------------
      if (nk.eq.0) cycle     
      if (debug) write(*,*)'radial integrals'
      do ie1=1,jpiset(m)%nex ! np(m)
      id1=idx(m,ie1)
      do ie2=1,jpiset(n)%nex ! np(n)
      id2=idx(n,ie2)
      nff=nff+1 !number of FF
      do irad2=1,nrad2
      do ik=1,NK
        i=pk(ik)%i 
		j=pk(ik)%f
        k=pk(ik)%k
		nq=pk(ik)%nq
		l=pk(ik)%l
        rmatc=pk(ik)%rmatc
		rmatn=pk(ik)%rmatn
!  radial integral: R(ik)=f*V(K,Q)*f'
        xsumn=(0.d0,0.d0)
        xsumc=(0.d0,0.d0)
        xsumad=0.
        xsumnad=0.
        if (abs(rmatc).lt.1e-5.and.abs(rmatn).lt.1e-5) goto 1200
      
!        ui=>frad(m,ie1,:,i)
!        uf=>frad(n,ie2,:,j)
!      call radquad(ui,uf,potQKn(:,irad2,nq,k))
      
      do iquad=1,nquad
      vmon=0d0
! AMM: frad(j) should be conjugate!! (fixed Sept 16)  
!      fprod=frad(m,ie1,iquad,i)*frad(n,ie2,iquad,j) ! i -> j = <j | V | i> 
      fprod=conjg(frad(n,ie2,iquad,j))*frad(m,ie1,iquad,i) ! i -> j = <j | V | i>
!      if (irad2.eq.1) 
!     & write(94,*) rquad(iquad),real(frad(m,ie1,iquad,i)), real(fprod)  
      if(l.eq.nq) then
      xintgn(iquad)=fprod*potQKn(iquad,irad2,nq,k)*xrad2(irad2)**l


      if ((m.eq.n).and.(ie1.eq.ie2).and.(n.eq.1).and.(k.eq.0)) then 
! AMM: For adiabatic JT potential
      xintad(iquad)=fprod*vclquad(iquad)*((
     &              potQKn(iquad,irad2,nq,k)+potQKc(iquad,irad2,nq,k))
     &              *xrad2(irad2)**l-VCOUL(xrad2(irad2),zp,zt,Rcc))
      vadnorm(iquad)=fprod*vclquad(iquad)
      endif

c subtract projectile-target monopole Coulomb
      if ((nq.eq.0).and.zp*zt.gt.1e-3.and.k.eq.0.and.ncoul.ne.1)
     & then
!     &  .and.(m.eq.n).and.(ie1.eq.ie2)) then 
       vmon=VCOUL(xrad2(irad2),zp,zt,Rcc)     
      endif

      xintgc(iquad)=fprod*
     & (potQKc(iquad,irad2,nq,k)*xrad2(irad2)**l-vmon)
      
      else 
!      write(0,*)'l,nq=',l,nq
      xintgn(iquad)=fprod*
     .potQKn(iquad,irad2,nq,k)*xrad2(irad2)**l*
     .(coefc*rquad(iquad))**(nq-l)

      xintgc(iquad)=fprod*
     .potQKc(iquad,irad2,nq,k)*xrad2(irad2)**l*
     .(coefc*rquad(iquad))**(nq-l)
      endif
      xsumn=xsumn + xintgn(iquad)*wquad(iquad)
      xsumc=xsumc + xintgc(iquad)*wquad(iquad)
      
! ADIABATIC      
      if ((m.eq.n).and.(ie1.eq.ie2).and.(n.eq.1).and.(k.eq.0)) then 
      xsumad=xsumad   +  xintad(iquad)*wquad(iquad)
      xsumnad=xsumnad + vadnorm(iquad)*wquad(iquad)
      !if (irad2.eq.1) write(0,*)iquad,vadnorm(iquad)
!      if (irad2.eq.1) write(0,*)iquad,xsumnad*0.5d0*radmax
      endif 
      
      enddo ! iquad
     

1200  fauxn=rmatn*xsumn*dsqrt(2.d0*dble(k)+1.d0)*0.5d0*radmax
      fauxc=rmatc*xsumc*dsqrt(2.d0*dble(k)+1.d0)*0.5d0*radmax
      
      
! ADIABATIC      
      if ((m.eq.n).and.(ie1.eq.ie2).and.(n.eq.1).and.(k.eq.0)) then 
!      write(0,*)'xsumna=',xsumnad
      xsumad=xsumad/xsumnad ! divide by norm <phi|Vbx| phi> 
      fauxad=rmatn*xsumad*dsqrt(2.d0*dble(k)+1.d0)!*0.5d0*radmax
      write(999,'(1f10.3,2g16.6)') xrad2(irad2),fauxad
      endif
      
! April 2019: only diagonal nuclear couplings for coups=2
      if (coups.eq.2.and.(ie1.ne.ie2.or.m.ne.n)) fauxn=0 

! March 2021      
      if (coups.eq.3.and.(ie1.ne.ie2.or.m.ne.n)) fauxn=0 
      
! Scale formfactors by qfactor()
!       write(*,*) lc,qfactorc(lc), qfactorn(lc),fauxc,fauxn
	   if (lc.lt.50) fauxc=fauxc*qfactorc(lc)
	   if (lc.lt.50) fauxn=fauxn*qfactorn(lc)
!       write(*,*) lc,qfactorc(lc), qfactorn(lc),fauxc,fauxn
 
!      Fc(id1,id2,lc,irad2)=Fc(id1,id2,lc,irad2)+fauxc    ! coulomb
!      Fn(id1,id2,lc,irad2)=Fn(id1,id2,lc,irad2)+fauxn    ! nuclear
      Ff(id1,id2,lc,irad2)=Ff(id1,id2,lc,irad2)+fauxn+fauxc ! total
      enddo ! ik
      enddo ! irad2 (R)
      enddo ! ie2
      enddo ! ie1
      call cpu_time(t2)
      write(*,'("(",1f6.2," secs )")')t2-t1
      call flush(6)
      enddo ! lc (\Lamda)
      enddo ! n (jpsets)
      enddo ! m (jpsets)

      call cpu_time(finish)
!      print*, 'Time = ',finish-start,'seconds'
      write(*,1220) nff,finish-start
1220  format(/,5x,"-> [",i7, " formfactors generated in ",
     & 1x,1f6.1," secs ]",/)



c    -----------------------------------------------------------
c    Write formfactors
c    ----------------------------------------------------------- 
      if (writeff) open(kfr,file="ff.fr",status='unknown')
      ttr=0.d0
      fscale=1.d0
      npa1=0
      write(*,*)' Output files:'
      select case(ncoul)
      case(0) ! nuc + coul
      write(*,*) ' ff.fr:  COUL+NUC couplings for Fresco'
!      ffr=>Ff(:,:,:,:,:,:)
      case(1) ! nuc
      write(*,*) ' ff.fr: NUCLEAR couplings for Fresco'
!      ffr=>Fn(:,:,:,:,:,:)
      case(2) ! coul
      write(*,*) ' ff.fr: COULOMB couplings for Fresco'
!      ffr=>Fc(:,:,:,:,:,:)
      end select 

      do if1=1,jpsets
      npa1=npa1+jpiset(if1)%nex ! np(if1)
      npa2=0
      xjp1 =jpiset(if1)%jtot
      ipar1=jpiset(if1)%partot
      do if2=1,jpsets
      npa2  =npa2+jpiset(if2)%nex ! np(if2)
      xjp2  =jpiset(if2)%jtot
      ipar2 =jpiset(if2)%partot

      do ie1=1,jpiset(if1)%nex  ! np(if1)
      do ie2=1,jpiset(if2)%nex  ! np(if2)
      e1=energ(if1,ie1)
      e2=energ(if2,ie2)

      m1=ie1+npa1-jpiset(if1)%nex !np(if1)
      m2=ie2+npa2-jpiset(if2)%nex !np(if2)

      select case(coups)
      case (0) ! all couplings (do nothing)
      case (1,3) ! gs->cont + diagonal
        if ((if1.ne.1).or.(ie1.ne.1)) then
          if ((if2.ne.if1).or.(ie2.ne.ie1)) cycle
        endif
        
      case(4) ! only gs-> cont NO DIAG
         if ((if1.ne.1).or.(ie1.ne.1)) cycle
         if ((if1.eq.if2).and.(ie2.eq.ie1)) cycle
     
        
!      case (3) ! V00 + gs->cont
!         if ((if1.ne.1).or.(ie1.ne.1)) cycle
!         if ((if2.ne.if1).or.(ie2.eq.ie1)) cycle
      end select

      if ((m1.gt.m2).and.realwf) cycle !!!! AMORO
c for fresco
      do lc= nint(dabs(xjp1-xjp2)),min(nint(xjp1+xjp2),lambmax)
      if (ipar1*ipar2*(-1)**lc<0) cycle
      ltr=lc
      ptr=lc
      
c \hat{Jp}*hat{Jp'}*(2*Lambda+1)*(-1)^Lambda
      factor=(2d0*ptr+1)*sqrt(2*xjp1+1)*sqrt(2*xjp2+1)*(-1)**ptr

      if (debug) then
      write(*,'(a,a,1f4.1,a,1f4.1,a,i3,2x,a,1g16.8)')
     & ' Conversion factor in F(r)[fres]=factor*F(r)',
     &  'ji=',xjp1,'-> jf=',xjp2,' LAM=',lc," =>  factor=",factor
      endif

      write(comment,'(2x,"<",i3,"|",i2,"|",i3,">")') m1,lc,m2
      if (writeff) 
     & write(kfr,500) nrad3,rstep,rfirst,fscale,ltr,ptr,ttr,m2,m1,
     &               comment   ! ORIG

!!! TEST APRIL
!      if (writeff) 
!     & write(888,500) nrad3,rstep,rfirst,fscale,ltr,ptr,ttr,m2,m1,
!     &               comment   ! ORIG


      if (verb.ge.4) write(120,'("# <",i3,"|",i2,"|",i3,">")') m1,lc,m2
      do irad=1,nrad2
      r2=xrad2(irad)
      fauxc=Ff(m1,m2,lc,irad)
      
      if ((coups.eq.3).and.(m1.eq.m2)) fauxc=Ff(1,1,lc,irad)
         
!      fauxn=0d0
!      if ((ncoul.eq.0).or.(ncoul.eq.1)) then
!       fauxn=Fn(m1,m2,lc,irad)
!      endif      
!      if ((ncoul.eq.0).or.(ncoul.eq.2)) then
!       fauxc=Fc(m1,m2,lc,irad)   
!      endif          
!      write(kfr,'(2x,1g16.10,2x,1g16.10)') factor*(fauxc+fauxn)
!      write(120,'(1x,1f8.3,2x,2g16.8)') r2,factor*(fauxc+fauxn)
       if (writeff) 
     &      write(kfr,'(2x,1g16.10,2x,1g16.10)') factor*fauxc
      if (verb.ge.4) 
     & write(120,'(1x,1f8.3,2x,2g16.8)') r2,factor*fauxc
      enddo !irad

c Check that imaginary part is small
      if (abs(aimag(fauxc)).gt.1e-5) then
      write(*,480) m1,xjp1,parity(ipar1+2),m2,xjp2,parity(ipar2+2),lc,
     & r2,fauxc
      endif
480   format('** Warning **: Large imaginary part for coupling: ',
     & 1i3,' (',1f4.1,a1,") -> ", 1i3," (",1f4.1,a1,")",
     & " for LAM=",1i3, ": R=",1f8.2," Fc(R)=",2g16.8)

c Extrapolate Coulomb formactors from R=Rmax to Rextrap
      if (rextrap.gt.rmax) then
      caux=real(fauxc)*rmax**(lc+1)
      
      do ir=nrad2+1,nrad3
      r2=rstep+rstep*(ir-1)
      fauxc=0
      if ((lc.gt.0).and.(ncoul.ne.1)) then
      fauxc=caux/r2**(lc+1)
      endif
     
      if (writeff) write(kfr,'(2x,1g16.10,2x,1g16.10)') factor*fauxc
      if (verb.ge.4) write(120,'(1x,1f8.3,2x,2g16.8)') r2,factor*fauxc
!      Fc(m1,m2,lc,ir)=fauxc 
      Ff(m1,m2,lc,ir)=fauxc ! total  
      enddo !nrad3
      endif ! rextrap> rmax
      
!!! TEST APRIL 18 
!      do ir=1,nrad3
!      r2=rstep+rstep*(ir-1)
!      fauxc=0
!      if ((lc.gt.0).and.(ncoul.ne.1)) then
!      fauxc=caux/r2**(lc+1)
!      endif
!      if (writeff) write(888,'(2x,1g16.10,2x,1g16.10)') factor*fauxc
!      enddo
!!!!!!!!!!! END TEST

      enddo !lc (next multipole)
     
      if (verb.ge.4) then
c ---------- Nuclear
!      write(10,800) '#','Jp1=',xjp1,'parity1=',ipar1,'Jp2=',xjp2,
!     & 'parity2=',ipar2,'ener1=',e1,'ener2=',e2
c ---------- Coulomb
!      write(11,800) '#','Jp1=',xjp1,'parity1=',ipar1,'Jp2=',xjp2,
!     & 'parity2=',ipar2,'ener1=',e1,'ener2=',e2
c ---------- nuclear + coulomb
      write(12,800) '#','Jp1=',xjp1,'parity1=',ipar1,'Jp2=',xjp2,
     & 'parity2=',ipar2,'ener1=',e1,'ener2=',e2
      endif


      i1=nint(dabs(xjp1-xjp2))
!      i2=nint(xjp1+xjp2)
      i2=min(nint(xjp1+xjp2),lambmax)
      if (verb.ge.4) then
      do irad=1,nrad3
!      r2=xrad2(irad)
      r2=rstep+rstep*dble(irad-1) 
      write(12,900) r2, (Ff(m1,m2,i,irad),i=i1,i2)
!      write(10,900) r2, (Fn(m1,m2,i,irad),i=i1,i2)
!      write(11,900) r2, (Fc(m1,m2,i,irad),i=i1,i2)
! TEST I4 -> I5
500   format(i4,3f8.4,i4,2f4.0,2i4,a35)
!500   format(i5,3f8.4,i4,2f4.0,2i4,a35)
800   format (a,2(a,(f8.4),2x,a,i3,2x),2(a,(f8.4),2x),/)
820   format (a,a5,a,a5,2(a,(f8.4),2x))
900   format (1f8.3,2x,100(e12.4,2x,e12.4,4x))
      enddo !irad
!      write(10,*) '&'
!      write(11,*) '&'
      write(12,*) '&'
      if (verb.ge.4)  write(120,*)'&'
      endif
      enddo 
      enddo
      enddo
      enddo
      write(*,*)
      if (verb.ge.2) then
!        write(*,*)'fort.10: nuclear formfactors'
!        write(*,*)'fort.11: Coulomb formfactors'
        write(*,*)'fort.12: nuclear + Coulomb formfactors'      
      endif ! verb

c     -----------------------------------------------------------
c     Write states & energies in Fresco format
c     ----------------------------------------------------------- 
      open(ken,file="states.fr",status='unknown')     
      i=0
      icopyt=0
      icpot=1
      ibandt=0
      et=0.d0
      jt=jtgs ! inherited from module 'sistema'
      write(ken,*) m1
      do n=1,jpsets
      rjp   =jpiset(n)%jtot    ! xjp(n)
      ibandp=jpiset(n)%partot  ! ipar(n)
      do ie=1,jpiset(n)%nex    ! np(n)
      i=i+1
      ep=energ(n,ie)-energ(1,1)
      if (i.eq.2) icopyt=1
      write(ken,11511) rjp,ibandp,ep,icpot
      if(icopyt==0) write(ken,1153) jt,ibandt,et
      if(icopyt.ne.0) write(ken,1154) icopyt
      enddo !ie
      enddo !n

11511 format(' &States jp=',f4.1,' ptyp=',i2,' ep=',f8.4,'  cpot=',i3)
 1153 format('         jt=',f4.1,' ptyt=',i2,' et=',f8.4,' /')
 1154 format('         copyt=',i2,' /')
 1155 format(' Skipping state #',i3,' J/pi',f4.1,i2,' and Ex=',f8.4)
c     ----------------------------------------------------------------


c      total = etime(elapsed)
c      print*,'total=',total,'user=',elapsed(1),
c     .'system=',elapsed(2)
      call flush(ken)
      call flush(kfr)
      close(ken); close(kfr)
!      deallocate(fn,fc)
      end



c *** Use projectile WFs from THOx *****************************
c     and interpolate them at quadrature points
c *** ---------------------------------------------------------
      subroutine wf2quad
      use sistema
      use channels, only: jpiset,jpsets,nchmax
      use wfs, only: wfr,wfc,energ,idx
      use xcdcc,only:nquad,rquad,frad,exch,parch,jpch,nex,elab
      use wfs, only:  wfeig,ebin,nr,rvec,dr,rmin
      use parameters, only: maxeset
      use trace, only: cdccwf
      implicit none
      integer:: jset,j,nchan,ne,nst,iq,ie,ich,iparity
c changed in v2.3
!      real*8,pointer:: yvp(:)
!      real*8 :: faux
      complex*16,pointer:: yvp(:)
      real*8      :: raux,xjtot
      complex*16  :: faux,caux,ffc4
      character*5 :: jpi
      integer,parameter:: alpha=0d0

!      write(*,*)'jpsets=',jpsets,' nquad=',nquad,' nchmax=',nchmax
!      allocate(frad(jpsets,nener,nquad,nchmax),chann(nchmax,4,jpsets))
      allocate(frad(jpsets,maxeset,nquad,nchmax))
      write(*,'(//,2x,"***  PROJECTILE WAVEFUNCTIONS *** ")')
      
      nst=0
      do jset=1,jpsets
      ne=jpiset(jset)%nex
      nst=nst+ne
      nchan=jpiset(jset)%nchan
      xjtot=jpiset(jset)%jtot
      iparity=jpiset(jset)%partot
      write(*,*)
      write(*,'(4x,"Set",i3,"  Jpi=",1x,a5,
     &    "=>",i3," state(s)" )')
     &    jset,jpi(xjtot,iparity),ne


      do ie=1,ne
      write(*,'(8x,"#",i3,2x,"Ex=",1f8.4)')
     & ie,energ(jset,ie)
      do ich=1,nchan
c changed in v2.3
!     yvp=>wfr(jset,ie,ich,:)
      yvp=>wfc(jset,ie,ich,:)

!      if(jset.eq.2) then
!      do iq=1,nr
!      write(95,'(1f8.3,2x,100g14.6)')rvec(iq),
!     &   rvec(iq)*wfc(jset,ie,ich,iq)
!      enddo
!      write(95,*)'&'
!      endif


      do iq=1,nquad
      raux=rquad(iq)
!      faux=cfival(raux,rvec,yvp,nr,alpha) ! R(r)
      faux=FFC4((raux-rmin)/dr,yvp,nr)
      frad(jset,ie,iq,ich)=faux*raux     ! u(r)=r*R(r)
!      write(95,'(1f8.3,2x,100g14.6)')raux,frad(jset,ie,iq,ich)
      enddo ! iq=quadrature point
      enddo ! ich=channel
      enddo ! ie = energy
!      write(*,'(5x,"Set:"i2," => ",i3," state(s) considered")')
!     &  jset,ne
      enddo ! jset=j/pi set
      write(*,'(/, 5x,"=>", i3,1x,"state(s) in ",i2," set(s)")')
     & nst,jpsets


c Write energies and internal WFS for full CDCC WF
      if (cdccwf) then
         write(85,'(5a,1f8.2,a)') '** CDCC WF for: ', 
     &   trim(namep),'+',trim(namet),'  at E=',elab,' MeV'
         write(85,'(a,i3)') 
     &   '# Nb. states: ', nst
         write(85,'(a,i5,2f8.3)') '# Radial grid:', 
     &   nr,dr,rmin
      endif

c     store excitation energies, parities for coupled-channels calculations
      nex=nst
      allocate(exch(nex))
      allocate(parch(nex))
      allocate(jpch(nex))
      allocate(idx(jpsets,maxval(jpiset(:)%nex)))
      
      nst=0
      do jset=1,jpsets
      ne=jpiset(jset)%nex
      do ie=1,ne
      nst=nst+1
      exch(nst)   = energ(jset,ie)
      parch(nst)  = jpiset(jset)%partot
      jpch(nst)   = jpiset(jset)%jtot
      idx(jset,ie)= nst

      if (cdccwf) then
       write(85,'(a,i3,1f10.5,1f5.1)') '# n, Ex, Jp:   ',
     &      nst,exch(nst),jpch(nst)
      do ich=1,jpiset(jset)%nchan
      write(85,'(6g14.6)')
     &   (rvec(iq)*wfc(jset,ie,ich,iq),iq=1,nr)
      enddo !nch
      endif !cdccwf?

      enddo ! ie
      enddo ! jset        

!      deallocate(wfr)
      return
      end subroutine

c nc=1: nuclear part
c nc=2: coulomb part
      complex*16 function pot(r1,r2,iq,k,coef,nc)
      use xcdcc,only:nquad,xquad,wquad
      implicit real*8 (a-h,o-z)
      complex*16 ss,func
      integer nc
      common/vari/nq,kk
      common/varr/xcoef,xr1,xr2
!      write(*,*)'pot: nc=',nc
      kk=k
      xcoef=coef
      xr1=r1
      xr2=r2
      nq=iq
      ind=mod(nquad,2)
      if(ind.eq.1) nquad=nquad-1
      ss=(0.d0,0.d0)
      do j=1,nquad/2
      dx=xquad(j)
      ss=ss+wquad(j)*(func(dx,nc)+func(-dx,nc))
      enddo
      pot=ss
      return
      end

c This is the integrand in eq.(22) of Summers et al
c Vq is the Q-multipole of the core-target or valence-target potential
c according to :
c
c Vct(\vec rc,\vec xi)=Sum_{Q,q} \hat{Q} V^Q(rc) C_{Qq}(\hat rc) C*_{Qq}(\hat xi)
c
c and likewise for the valence-target
c
c nc=1: nuclear part
c nc=2: coulomb part
      complex*16 function func(xx,nc)
      use ptpots, only:vcore,vcorei,vval,vvali,vcorec,vvalc,rfrag,nr
!      use wfs, only:rfrag,nr
      use constants, only:e2
      use sistema, only:zc,zv,zt,rcc,zp
      implicit real*8 (a-h,o-z)
      real*8 faux,zfrag,vcoul,raux,xpos
      complex*16 xi,vq
      common/vari/iq,k
      common/varr/coef,r1,r2 ! r,R ???
c     --------------------------------------------------------------
      xi=(0.d0,1.d0)
      alpha=0.d0
      func =0 

      xrvec=dsqrt((coef*r1)**2+r2**2-2.d0*r1*r2*coef*xx)
      if (coef.lt.0) then
         zfrag=zv
      else
         zfrag=zc
      endif
 
      if ((nc.eq.2).and.(zfrag.lt.1e-4)) return

      if(nc.eq.1) then                    ! nuclear part
      xpos=(xrvec-rfrag(1))/(rfrag(2)-rfrag(1)) 
      if(coef.lt.0.d0) then
!      vvalp=fival(xrvec,rfrag,vval(:,iq),nr,alpha)
        vvalp=FFR4(xpos,vval(:,iq),nr)
!       vvalip=fival(xrvec,rfrag,vvali(:,iq),nr,alpha)
       vvalip=ffr4(xpos,vvali(:,iq),nr)
       vq=vvalp + xi*vvalip
      else             
!        vcorep=fival(xrvec,rfrag,vcore(:,iq),nr,alpha)
!        vcoreip=fival(xrvec,rfrag,vcorei(:,iq),nr,alpha)
        vcorep=FFR4(xpos,vcore(:,iq),nr)
        vcoreip=FFR4(xpos,vcorei(:,iq),nr)
        vq=vcorep+xi*vcoreip
      endif
      func=Vq*fleg(xx,k)/(2.d0*xrvec**iq) !c.f. Eq. (22) of Summers 
      endif
 
      if((nc.eq.2).and.(zfrag.gt.0)) then      ! Coulomb part
      xpos=(xrvec-rfrag(1))/(rfrag(2)-rfrag(1)) 
      if (coef.lt.0) then ! valence-target
!        Vq=fival(xrvec,rfrag,vvalc(:,iq),nr,alpha)
        Vq=FFR4(xpos,vvalc(:,iq),nr)
       else                ! core-target
!        Vq=fival(xrvec,rfrag,vcorec(:,iq),nr,alpha)
        Vq=FFR4(xpos,vcorec(:,iq),nr)
!        if (iq.eq.2) write(199,*)xx, vq
       endif
       func=Vq*fleg(xx,k)/(2.d0*xrvec**iq) !c.f. Eq. (22) of Summers
      endif     
      RETURN
      end
