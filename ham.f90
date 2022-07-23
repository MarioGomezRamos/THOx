
c 
c *** Matrix elements for rotor model
c 
!      function rotor(i,f,k,lambda,rms,delta,zc)
!        use globals, only: pi
!        implicit none
!        real*8::rms,delta,rotor,cleb
!        real*8::i,f,lambda,k,zc
!        rotor=3*zc*delta*(sqrt(5./3.)*rms)**(lambda-1)/(4*pi)
!        rotor=rotor*sqrt((2.*i+1.)/(2.*f+1.))*(-1)**(2*lambda)
!        rotor=rotor*cleb(i,k,lambda,0d0,f,k)
!      end function


c 
c *** Build & diagonalize single-particle Hamiltonian matrix
c     for each channel {l,s,j} 
      subroutine hsp(nset,nch,nho)
      use globals
      use constants
      use wfs
      use channels
      use hmatrix
      implicit none
      real*8,allocatable:: gnm(:),faux(:),gaux(:)
      integer i,l,ich,nch,n,m,n1,m1,nprint,ifile,nset,nho
      parameter(nprint=11)
      real*8:: norm,r,rms,jn,vscale
      CHARACTER LNAME(0:14)
      DATA LNAME / 's','p','d','f','g','h','i','j','k','l','m',
     &                'n','o','p','q' /

!      write(*,'(/,2x,a40)')'SINGLE-PARTICLE Hamiltonian:'
      if (nho.eq.0) then
        write(*,*)'hsp: nho=0 ! Skipping this basis set'
        return
      endif
      
      jtot   =jpiset(nset)%jtot
      vscale= jpiset(nset)%vscale

      write(*,'(//,5x,"SINGLE-PARTICLE eigenstates:",/)')
      if (allocated(hmat))   deallocate(hmat)
      if (allocated(ortmat)) deallocate(ortmat)
      if (allocated(eigv))   deallocate(eigv)
      if (allocated(faux))   deallocate(faux,gaux)
      if (allocated(esp))    deallocate(esp)
      if (allocated(wfsp))   deallocate(wfsp)
      allocate(hmat(nho,nho)) ! s.p. Hamiltonian matrix
      allocate(ortmat(nho,nho)) !orthogonality Matrix
      allocate(eigv(nho,nho))
      allocate(esp(nch,nho)) !s.p. energies
      allocate(wfsp(nch,nho,nr)) !s.p. wfs
      allocate(faux(nr),gaux(nr))
      if (verb.ge.3) written(76)=.true.
      

      do ich=1,nch
        l=qspl(ich) 
        jn=qspj(ich)
!        write(*,*)'set,chan,l,jn',nset,ich,l,jn
!        write (*,*) spchan(nset,ich)%l,spchan(nset,ich)%j
c l, jn, etc could be also retrieved with the spchan() variable:
!        l =spchan(nset,ich)%l
!        jn=spchan(nset,ich)%j
	wfaux=>wftho(:,l,:)
        if (verb.ge.3) then
        write(76,'(/," # S.p. Hamiltonian for channel",i3,
     &  " and configuration", 1x,a1,1f3.1)')ich,lname(l),jn
        endif

        do n=1,nho
        do m=1,nho
        call hspnm(n,m,ich,vscale)  
        enddo !m
        if (written(76)) write(76,'(100f12.5)')(hmat(n,m),m=1,nho)
        enddo !n
! diagonalize single-particle Hamiltonian
        ieigen=0
        nrot=0
        eigv(:,:)=0.d0
        do n=1,nho
        eigv(n,n)=1.d0
        enddo
        call hdiag(hmat,nho,nho,ieigen,eigv,nrot)
      

! build sp basis for each ichan
      write(*,'(2x,"o S.p. configuration #",i2,1x,"=>",a1,f3.1)') 
     & ich,lname(l),jn
      do n=1,nho
         n1=nho-n+1
         faux=0d0
         esp(ich,n)=hmat(n1,n1)
         do m=1,nho
          m1=nho-m+1
          faux(:)=faux(:)+eigv(n1,m)*wfaux(m,:)
         enddo !m=1,nho
         wfsp(ich,n,1:nr)=faux(1:nr) 
             
c ** Calculate  & print  NORM and RMS
         norm=0.d0;rms=0d0
         do i=1,nr
          r=rmin+dr*dble(i-1)
          rms=rms+(faux(i)*r)**2*r**2*dr
          gaux(i)=(faux(i)*r)**2
          enddo
          call sim(gaux,norm,1,nr,dr,nr)
          rms=sqrt(rms)

         
          if ((n<nprint).and.(verb.gt.2))
     &        write(*,100) n,esp(ich,n),norm,rms
          if ((n.eq.nprint).and.(verb.gt.2))
     &        write(*,'(7x,"(further eigenstates not printed)",/)')
100   format(5x,"#",i2,2x,"ex=",1f10.4,2x,"Norm=",1f8.4," rms=",1f8.4)
        enddo !n=1,nho
      enddo !ich={lsj}

! PRINT WFS
       if (verb.ge.3) then
       do ich=1,min(nch,5)
         ifile=60+ich-1
	  written(ifile)=.true.
         write(ifile,*)'#S.p. wfs for channel:', ich
         do i=1,nr
         r=rmin+dr*dble(i-1)
         write(ifile,200) r,(r*wfsp(ich,n,i),n=1,min(3,nho))
         enddo
       enddo
       endif
200   format(1f8.3,2x,10f10.5)
      end    subroutine   
 


c *** ----------------------------------------------------
c     Calculates matrix elements for s.p. Hamiltonian:
c     hnm(n,m)=<m,nchanf|h|n,nchani>
c *** ---------------------------------------------------
      subroutine hspnm(n,m,ichan,vscale)
      use globals
      use sistema
      use potentials, only: vls,vss,vcl,vcou,vll,cptype,vtran
      use constants
      use wfs
      use channels
      use hmatrix
      implicit none
      real*8::h,fact,d2wf,xl,r,d1wfm,d1wfn
      real*8::deriv1,als,aspsp,all,vscale
      integer::l,i,n,m,ichan
      real*8::j,res
      real*8:: un,um,um1,um2,um3,um4,um5,aux,knm,vnm
c        real*8:: z1,z2,a1,a2,rc,rc2,aux6
c        real*8:: all,als,as0,j,vcou,vpot,deriv1,deriv2
      real*8,allocatable:: fmaux(:),fnaux(:),hnmaux(:)

      integer:: ic !jc
      real*8 :: jcore

      l=qspl(ichan)

      j    =qspj(ichan)
      jcore=qjc(ichan)
      ic   =cindex(ichan)
      

!      l =spchan(nset,ich)%l
!      jn=spchan(nset,ich)%j

      xl=l
      als=0.5*(j*(j+1)-xl*(xl+1.)-sn*(sn+1)) ! l.s
      all=xl*(xl+1.)                         ! l.l
      fact=hc*hc/2d0/mu12
      h=dr
       
      if (allocated(hnmaux)) deallocate(hnmaux)
      if (allocated(onm))    deallocate(onm)
      if (allocated(fnaux))  deallocate(fnaux)
      if (allocated(fmaux))  deallocate(fmaux)

      allocate(fmaux(nr),fnaux(nr))
      allocate(hnmaux(nr),onm(nr))
      do i=1,nr
        r=rvec(i)
        fnaux(i)=wfaux(n,i)*r
        fmaux(i)=wfaux(m,i)*r
      enddo
              
      rmax=rvec(nr)
  
	 
      do i=1,nr
           r=rvec(i)
           if (r.lt.1e-9) r=1e-9
           um=wfaux(m,i)*r
           un=wfaux(n,i)*r
           onm(i)=um*un
      
           d1wfm=deriv1(fmaux,h,nr,i)
           d1wfn=deriv1(fnaux,h,nr,i)
           knm=fact*d1wfn*d1wfm                !kinetic energy term
     &           + fact*xl*(xl+1d0)/r/r*un*um  !centrifugal term

	   vnm=  (vcou(i)+vscale*vcl(l,i)+als*vls(l,i)+all*vll(l,i))*un*um !coulomb+nucl(only s-o)  MGR l-dependence
!	   write(0,*)'ir,n,m,vnm=',i,r,vnm

           if (cptype.eq.5) then
!                 write(*,*)'hspnm: vscale',vscale
             if (partot.eq.1) then
 	     vnm= vnm+(vscale*vtran(ic,ic,0,i,1)/sqrt(2.*jcore+1))*un*um
             else
             vnm=vnm+(vscale*vtran(ic,ic,0,i,2)/sqrt(2.*jcore+1))*un*um
             endif
           endif

       hnmaux(i)=knm+vnm
c ***--------------------------------------------------------
        enddo
	 call sim(hnmaux,res,1,nr,h,nr)
       hmat(n,m)=res
	 call sim(onm,res,1,nr,h,nr)
       ortmat(n,m)=res

      deallocate(fnaux,fmaux,onm,hnmaux)
      end subroutine hspnm


c *** ----------------------------------------------------------
c     <nchani; n | H | nchanf; m > block in configuration space
c *** ----------------------------------------------------------
      subroutine hmat2(nchani,nchanf)
        use globals
        use constants
        use wfs
        use channels, only: ql,spindex,nsp,bas2,jtot
        use hmatrix
        implicit none
        integer:: n,m,nchani,nchanf,ndim
        integer:: li,lf
        real*8::  ji,jf,jci,jcf,pari,parf
        real*8::  hnm !,r,fact
        ndim=nsp
        written(70)=.false.
        written(72)=.false.

        if (verb.ge.2) then
        written(70)=.true.
        written(72)=.true.
        endif
 
	if (allocated(hmat)) deallocate(hmat)
        if (allocated(ortmat)) deallocate(ortmat)
        if (.not.allocated(hmat)) allocate(hmat(ndim,ndim))
        if (.not.allocated(ortmat)) allocate(ortmat(ndim,ndim))

c choose basis for block diagonalization
c THO basis
        if (bas2.eq.1) then  
	   wfaux=>wftho(:,ql(nchani),:)
	   wfaux2=>wftho(:,ql(nchanf),:)
        else if (bas2.eq.0) then
c single-particle eigenstates basis
           wfaux=>wfsp(spindex(nchani),:,:)
	   wfaux2=>wfsp(spindex(nchanf),:,:)
        else
           write(*,*) 'Basis type=',bas2,'not implemented!'
           stop
        endif


	if (written(70)) then
        write(70,700) nchani,nchanf
700     format("Full H matrix for nchani,nchanf=",i3,i3," ")
        write(72,720)nchani,nchanf
720     format("Orthogonality matrix for nchani,nchanf=",i3,i3," ")
        endif
       
         do n=1, ndim
         do m=1, ndim
             hmat(n,m)=hnm(n,nchani,m,nchanf)
         enddo
         if (written(70)) write(70,'(100f12.5)') (hmat(n,m),m=1,ndim)
	  if (written(72)) write(72,'(100f12.5)') (ortmat(n,m),m=1,ndim)
         enddo

      end subroutine hmat2


c Identity matrix
	subroutine identity(ndim)
	use hmatrix,only:idmat
	implicit none
	integer:: ndim,n,m
	if(allocated(idmat)) deallocate(idmat)
	allocate(idmat(ndim,ndim))
	idmat(:,:)=0d0
	do n=1,ndim
	do m=1,ndim
	if (n.eq.m) idmat(n,m)=1d0
	enddo
	enddo
	end subroutine



c *** ----------------------------------------------------
c     Matrix element hnm(n,m)=<m,nchanf|h|n,nchani>
c *** ---------------------------------------------------
      function hnm(n,nchani,m,nchanf)
      use globals
      use sistema
      use potentials, only: vls,vlsc,vss,vcl,vcou,vlcoup,
     &lpot,vll,vtran,cptype,maxlamb,  
     &lambda,pcmodel,laminc,vscale ! v2.1          
      use constants
      use channels
      use wfs
      use hmatrix
      use forbidden
      implicit none
      real*8:: hnm,res,wfpau
      real*8::h,fact,d2wf,xl,r,d1wfm,d1wfn
      real*8::deriv1,als,ass,alsc,corels,coefss,atens
      real*8 allp,aspsp,ccoef,crotor,gcoup,all,ctens
      integer::i,ip,n,nchani,m,nchanf,ichn,ichm
      integer li,lf
      real*8::ji,jf
      real*8:: un,um,um1,um2,um3,um4,um5,aux,knm,vnm
      real*8,allocatable:: fmaux(:),fnaux(:),hnmaux(:)
      real*8,allocatable:: fpaux(:),gpaux(:)
      real*8:: uno,cero
      parameter(uno=1d0, cero=0d0)
      integer:: ici,icf,il,ichsp,nphi,nphf
      real*8:: jci,jcf,lambdar,sixj,cleb,lfr,lir,cl
      real*8:: cmic,cmic2,raux
      character*10 chname
      real*8:: pfactor,p1,p2
      real*8:: ui,uf,vi,vf,ren
      real*8:: acoup(0:maxlamb),vcp(0:maxlamb)
      logical:: ifren=.false.
      real*8 :: ki,kf,kband
      real*8 :: allp2 !!!! TEST
         
      ichn=spindex(nchani)
      ichm=spindex(nchanf)

      li=ql(nchani); xl=li
      lf=ql(nchanf)
      ji=qj(nchani)
      jf=qj(nchanf)

      jci=qjc(nchani) 
      ici=cindex(nchani)
      nphi=nphon(ici)
      ki =qnc(ici)%kband

      jcf=qjc(nchanf)
      icf=cindex(nchanf)
      nphf=nphon(icf)
      kf =qnc(icf)%kband

      cl=1d0

      fact=hc*hc/2d0/mu12
      h=dr

      if (ki.eq.kf) kband=ki

C ***      Matrix elements of tensor operators           ****
c Spin-orbit for valence ------------------------------
      als=0.5*(ji*(ji+1.d0)-xl*(xl+1.d0)-sn*(sn+1d0))
!      if (abs(als).gt.0) write(0,*)'als=',als


c spin-orbit for core (added to v2.0.2 version by AMoro)
      alsc=corels(li,sn,ji,qjc(nchani),lf,sn,jf,qjc(nchanf),jtot)
!      if (abs(alsc).gt.0) write(0,*)'alsc=',alsc

c spin.spin term 
      ass= coefss(li,sn,ji,qjc(nchani),lf,sn,jf,qjc(nchanf),jtot)
      
c tensor
      atens=ctens(li,sn,ji,qjc(nchani),lf,sn,jf,qjc(nchanf),jtot)

c l.l term -----------------------------------------------
       all=xl*(xl+1.d0)

c particle-core coupling ---------------------------------
      acoup(:)=0d0
      do il=0,maxlamb
      if (.not.laminc(il)) cycle   
!      allp=ccoef(ql(nchanf),qj(nchanf),jcf,nphf,
!     &           ql(nchani),qj(nchani),jci,nphi,
!     &           jtot,kband,sn,il,pcmodel)          
!      allp=allp*sqrt(2d0*lambda+1d0)/sqrt(4d0*pi)!delta included in vcp(ir) 

      allp=gcoup(ql(nchanf),qj(nchanf),jcf,nphf,
     &           ql(nchani),qj(nchani),jci,nphi,
     &           jtot,kband,sn,il,pcmodel)          
      allp=allp*sqrt(2d0*lambda+1d0)/sqrt(4d0*pi)!delta included in vcp(ir) 

!      if(abs(allp).gt.1e-3) write(*,*)'hnm: pcmodel,il,allp=',
!     & ql(nchanf),qj(nchanf),jcf,
!     & ql(nchani),qj(nchani),jci,
!     & il,allp


      acoup(il)=allp
!      if ((nchani.eq.2).and.(nchanf.eq.2)) then
!      write(99,*) 'acoup:',il,lambda,allp,nchani,nchanf,jci,jcf,jtot,
!     & kband,sn
!      endif
      enddo

   
      if (allocated(hnmaux)) deallocate(hnmaux)
      if (allocated(onm))    deallocate(onm)
      if (allocated(fnaux))  deallocate(fnaux)
      if (allocated(fmaux))  deallocate(fmaux)
      if (allocated(fpaux))  deallocate(fpaux)
      if (allocated(gpaux))  deallocate(gpaux)
      allocate(hnmaux(nr),onm(nr))
      allocate(fmaux(nr),fnaux(nr))
      allocate(fpaux(nr),gpaux(nr))
      do i=1,nr
        r=rvec(i)
        fmaux(i)=wfaux2(m,i)*r
        fnaux(i)=wfaux(n,i)*r
      enddo
        
      fpaux(:)=0d0; gpaux(:)=0d0
	     

      do i=1,nr
        r=rvec(i)
        if (r.lt.1e-9) r=1e-9
        um=wfaux2(m,i)*r
        un=wfaux(n,i)*r
        onm(i)=um*un	    
        knm=0; vnm=0; vcp(:)=0d0
        do il=0,maxlamb
        if (.not.laminc(il)) cycle
	select case (lpot)
        case (0) ! mininum li,lf
            vcp(il)=vlcoup(il,min(li,lf),i)
        case (1) ! average
            vcp(il)=(vlcoup(il,li,i) +vlcoup(il,lf,i))/2d0
        case (2) ! maximum li,lf
            vcp(il)=vlcoup(il,max(li,lf),i)
        case (3) ! left
            vcp(il)=vlcoup(il,lf,i)
        case default
            write(*,*) 'lpot=',lpot,'not valid!'; stop
        end select
        enddo !il
!!$c ***-----------------   Method 1    ------------------------
!!$           d2wf=deriv2(fmaux,h,nptbin,i)            
!!$           if (r.lt.1e-9) r=1e-9
!!$           hnm(i)=-fact*d2wf           !kinetic energy term
!!$     &           + fact*xl*(xl+1d0)/r/r*um  !centrifugal term
!!$     &           + vcou(r)*um            !coulomb
!!$     &           + vpot(r)*um            !nuclear
!!$           hnm(i)=un*hnm(i) 
!!$!           if ((n.eq.1).and.(m.eq.1)) write(96,'(12f12.5)')r,um,d2wf,
!!$!     &      vcou(r),vpot(r), fact*xl*(xl+1d0)/r/r*um
!!$           
!!$c ***--------------------------------------------------------

       
c ***-----------------   Method 2    ------------------------   
           d1wfm=deriv1(fmaux,h,nr,i)
           d1wfn=deriv1(fnaux,h,nr,i)
c Diagonal term
           if (nchani.eq.nchanf) then       !diagonal terms
           if ((bas2.eq.0).and.(cptype.ne.5)) then
             ichsp=spindex(nchani)
!             if (i.eq.1) write(*,*) ichsp,n,esp(ichsp,n)
             vnm=0
             knm=un*um*esp(ichsp,n)
           do  ip=1,npauli
          if (li.eq.paul(ip).and.(ji.eq.pauj(ip))
     &  .and.(n.eq.paun(ip))) then
            if (eqsp(ip).ne.0) knm=un*um*eqsp(ip)
!            if (i.eq.1) write(*,*)n,li,ji,eqsp(ip)            
          endif
           enddo
           else ! bas2.ne.0
             knm=fact*d1wfn*d1wfm           !kinetic energy term
     &           + fact*xl*(xl+1d0)/r/r*un*um  !centrifugal term
             vnm= (
     &  vcou(i)        ! central Coulomb  
     &  + vcl(li,i)    ! central 
     &  + als*vls(li,i)   ! for ptype 5 it would sum vls anyway   MGR  l-dependence   
     &  + all*vll(li,i)   ! l.l term    MGR l-dependence
     &  )*un*um ! 
           endif  
           vnm=vnm + exc(nchani)*un*um  ! rotational energy of the core

           endif

c Non-diagonal
           vnm=vnm  
     &  + (alsc*vlsc(li,i)        ! core spin-orbit MGR l-dependence     
     &  + sum(acoup(:)*vcp(:))    ! particle-core coupling
     &    )*un*um

       if ((li.eq.lf).and.abs(ass*vss(li,i)).gt.0) then
!       if (i.eq.1) write(*,*)'spin.spin included! ass=',ass,vss(li,i)
             vnm=vnm + ass*vss(li,i)*un*um  ! spin-spin 
       endif


!old particle core
!     &  + allp*vcp      ! particle-core coupling

           if (cptype.eq.5) then
           lfr=1d0*lf
           lir=1d0*li
           do il=0,maxlamb
           lambdar=1d0*il
!           if (il.eq.0) then
!             cl=sqrt(4*pi)    !maybe we would need to divide
!           else
!             cl=1d0
!           endif
           raux=1.0
           if ((il.eq.0).and.(nchani.eq.nchanf)) raux=vscale
               if (partot.eq.1) then
           vnm=vnm+cmic2(lir,ji,jci,lfr,jf,jcf,jtot,sn,lambdar)
     &        *vtran(ici,icf,il,i,1)*un*um*raux
               else
            vnm=vnm+cmic2(lir,ji,jci,lfr,jf,jcf,jtot,sn,lambdar)
     &        *vtran(ici,icf,il,i,2)*un*um*raux
               endif

!           vnm=vnm+((-1d0)**(ji+jcfr+jtot)*sixj(jf,ji,lambdar
!     &,jcir,jcfr,jtot)*sqrt((2*jf+1.)*(2*ji+1.))*sixj(jf,ji,lambdar,
!     &lir,lfr,sn)*(-1d0)**(lambdar+lfr+ji+sn)*sqrt((2.*lf+1)*(2.*il+1.)
!     &/4/pi)*cleb(lambdar,0d0,lfr,0d0,lir,0d0)*
!     &vtran(jci,jcf,ici,icf,il,i))*un*um*cl     

           enddo ! lambda
           endif ! cptype.eq.5

           hnmaux(i)=knm+vnm

        enddo !ir




c ***--------------------------------------------------------

!!! Remove PAULI forbidden states ----------------   
        pfactor=1; p1=1; p2=1;
        do  ip=1,npauli
        if (li.eq.paul(ip).and.(ji.eq.pauj(ip))
     &  .and.(n.eq.paun(ip))) then

          p1=pscale(ip,ici)
 
         endif

         if (lf.eq.paul(ip).and.(jf.eq.pauj(ip)).
     &   and.(m.eq.paun(ip))) then
 
          p2=pscale(ip,icf)
 
          endif
        enddo !npauli
        hnmaux(:)=hnmaux(:)*p1*p2

c   v2.0.4
c   Renormalize non-diagonal couplings Quasi-Particle Core Model
      ui=1; uf=1; vi=0; vf=0; ifren=.false.
!      if (ichn.ne.ichm) then 
      if ((li.ne.lf).or.(ji.ne.jf).or.(n.ne.m).
     & or.(ici.ne.icf)) then
        do  ip=1,npauli
         if (li.eq.paul(ip).and.(ji.eq.pauj(ip))
     &  .and.(n.eq.paun(ip))) then
         ui=uocup(ip); vi=vocup(ip) ;
         ifren=.true.
         endif
         if (lf.eq.paul(ip).and.(jf.eq.pauj(ip))
     &  .and.(m.eq.paun(ip))) then
         uf=uocup(ip); vf=vocup(ip)
         ifren=.true.
         endif
        enddo !ip
         ren=(ui*uf-vi*vf)
!        if (abs(ren).lt.1d0) then
!        write(*,1500) n,li,ji,m,lf,jf,ui,vi,uf,vf,ren
1500    format("(",2i2,1f4.1,2x,")",
     &         "(",2i2,1f4.1,2x,")",4f5.1,":",1f6.3)
!        endif
        hnmaux(:)=hnmaux(:)*ren       
      endif 
              
!!!---------------------------------------------------------
        
        
	if(ng.eq.0) then !Simpson integration
          call sim(onm,res,1,nr,dr,nr)
          ortmat(n,m)=res
          call sim(hnmaux,hnm,1,nr,dr,nr)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!AÑADIR CUADRATURAS AQUÍ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else 
          write(*,*)' Integration by quadratures not yet implemented!'
          stop
        endif
      end function hnm


c *** Calculate du(r)/dr using five points derivative formula
      function deriv1(f,h,ndim,j)
        implicit none
        integer ndim,j,i
        real*8 f(ndim),h,deriv1
 

        if ((j.eq.1).or.(j.eq.2)) then
           deriv1=(-f(j+2)+4d0*f(j+1)-3d0*f(j))/2d0/h
        else if (j.eq.ndim-1) then
           deriv1=(3d0*f(j)-4d0*f(j-1)+f(j-2))/2d0/h
        else if (j.eq.ndim) then
           deriv1=0 !!!CHECK
        else ! five points formula
           deriv1=(f(j-2)-8*f(j-1)+8*f(j+1)-f(j+2))/h/12.
         end if
      end function deriv1


c *** Calculate d^2u(r)/dr^2 using five points derivative formula
c     f(ndim)=function to make derivative
c     h      =step
c     j      =point for derivative
      function deriv2(f,h,ndim,j)
        implicit none
        integer ndim,j,i
        real*8 f(ndim),h,deriv2
!        write(*,*) 'deriv'
!        do i=1,ndim
!           write(96,*) i,f(i)
!        enddo

        if (j.eq.1) then
           deriv2=(f(j+2)-2*f(j+1)+f(j))/h/h
        else if (j.eq.ndim) then
           deriv2=(f(j)-2*f(j-1)+f(j-2))/h/h           
        else if ((j.eq.2).or.(j.eq.ndim-1)) then
           deriv2=(f(j+1)-2*f(j)+f(j-1))/h/h
        else 
           deriv2=(f(j+1)-2*f(j)+f(j-1))/h/h ! O(h^2)
           deriv2=(-f(j+2)+16*f(j+1)-30*f(j)+16*f(j-1)-f(j-2))/(12*h*h)
        end if
      end function deriv2


c---------------------------------------------------------------------------------------------
c      coeficients for coupling with externally read potentials
	function cmic(lir,ji,jcir,lfr,jf,jcfr,jtot,sn,lambdar)
        use constants, only: pi
        implicit none
        real*8::lir,ji,jcir,lfr,jf,jcfr,jtot,sn,lambdar
        real*8::cmic,cl,cleb,sixj

           if (lambdar.eq.0d0) then
             cl=sqrt(4*pi)    !maybe we would need to divide
           else
             cl=1d0
           endif

        cmic=(-1d0)**(ji+jcfr+jtot)*sixj(jf,ji,lambdar
     &,jcir,jcfr,jtot)*sqrt((2*jf+1.)*(2*ji+1.))*sixj(jf,ji,lambdar,
     &lir,lfr,sn)*(-1d0)**(lambdar+lfr+ji+sn)*sqrt((2.*lfr+1)*
     &(2.*lambdar+1.)
     &/4/pi)*cleb(lambdar,0d0,lfr,0d0,lir,0d0)*cl

        end function cmic
c-----------------------------------------------------------------------------------------------------

	function cmic2(lir,ji,jcir,lfr,jf,jcfr,jtot,sn,lambdar)
        use constants, only: pi
        implicit none
        real*8::lir,ji,jcir,lfr,jf,jcfr,jtot,sn,lambdar
        real*8::cmic2,cl,cleb,sixj
        real*8::ph1,ph2,red

           if (lambdar.eq.0d0) then
             cl=sqrt(4*pi)    !maybe we would need to divide
           else
             cl=1d0
           endif


	ph2=(-1)**(lambdar+sn+jf)*(-1)**(lambdar-lir)
       red=ph2*sqrt((2d0*jf+1)*(2d0*lir+1.)*(2d0*lfr+1d0))
     &   *sqrt((2d0*lambdar+1)/4d0/pi)
     &   *sixj(ji,jf,lambdar,lfr,lir,sn)
     &   *cleb(lir,0d0,lambdar,0d0,lfr,0d0)/sqrt(2d0*lfr+1d0)


	ph1=(-1)**(jtot+jf+jcir)
	cmic2=ph1*sqrt((2d0*ji+1)!*(2d0*jcir+1d0)
     &       )*sixj(ji,jf,lambdar,jcfr,jcir,jtot)
     &       *red
!     &       *cleb(ri,kk,rlambda,zero,rf,kk)

        cmic2=cmic2*cl!*sqrt(4d0*pi/(2d0*lambdar+1)) ! 

        end function cmic2



c -----------------------------------------------------------
c <lf,jf,rf,Jtot,M_f|P_{lambda}(cos(theta))|li,ji,ri,Jtot,M_i>
c ------------------------------------------------------------
	function crotor(lf,jf,rf,li,ji,ri,jtot,kproj,sn,lambda)
        implicit none
	integer:: lambda,li,lf
	real*8:: ji,ri,jf,rf,jtot,sn,crotor
	real*8:: phase,zero,one,rlambda,kk,half,kproj
	real*8:: cleb,sixj,pi
	parameter(zero=0d0,one=1d0,half=0.5d0)

	rlambda=lambda
	kk=kproj
	pi=acos(-one)
 

	phase=(-1)**(jtot+jf+rf+rlambda)
 
c < | Y_lambda | > 
	crotor=sqrt((2d0*rlambda+1)*(2d0*ri+1)*(2d0*ji+1)/4d0/pi)*
     &         cleb(ji,half,rlambda,zero,jf,half)*
     &         sixj(ji,ri,jtot,rf,jf,rlambda)
     &       *cleb(ri,kk,rlambda,zero,rf,kk)


c < | P_lambda | > 
        crotor=crotor*sqrt(4d0*pi/(2d0*rlambda+1)) ! 
        end function


c -----------------------------------------------------------
c <li,ji,ri,Jtot,M_i|P_{lambda}(cos(theta))|lf,jf,rf,Jtot,M_f>
c 
	function gcoup(li,ji,ri,ni,lf,jf,rf,nf,
     &                  jtot,kproj,sn,lambda,model)
       implicit none
	integer::lambda,li,lf,ni,nf,model
	real*8:: ji,ri,jf,rf,jtot,sn,gcoup,vibrator
	real*8:: zero,one,rlambda,kk,half,red2, kproj
	real*8:: cleb,sixj,pi,ph1,ph2,rli,rlf,threej
	parameter(zero=0d0,one=1d0,half=0.5d0)
c      --------------------------------------------------------
	rlambda=lambda
	kk=kproj
	rli=li
	rlf=lf
	pi=acos(-one)

! <(lsj||Y_{lambda} || (l's)j'>  
	ph2=(-1)**(rlambda+sn+jf+2*rli)
       red2=ph2*sqrt((2d0*jf+1)*(2d0*rli+1.)*(2d0*rlf+1d0))
     &   *sqrt((2d0*rlambda+1)/4d0/pi)
     &   *sixj(ji,jf,rlambda,rlf,rli,sn)
     &   *threej(rli,rlambda,rlf,zero,zero,zero)
!     &   *cleb(rli,zero,rlambda,zero,rlf,zero)/sqrt(2d0*rlf+1d0)


 
	ph1=(-1)**(jtot+jf+ri)
c < j R; J || Y_lambda | j' R') J> 
	gcoup=ph1*sqrt((2d0*ji+1)*(2d0*ri+1d0))
     &        *sixj(ji,jf,rlambda,rf,ri,jtot)
     &        *red2


      select case (model)
      case(0) ! rotor
       gcoup=gcoup*cleb(ri,kk,rlambda,zero,rf,kk)
     &       *(-1)**(ri-rf)
      case(1) ! vibrator
!       if (abs(ccoef).gt.0.) write(*,*)'vib=',ccoef
       gcoup=gcoup*vibrator(ni,ri,nf,rf,lambda)
      case default
       write(*,*)'Core model',model,' not used!'
       stop
      end select

      gcoup=gcoup*sqrt(4d0*pi/(2d0*rlambda+1)) ! 
      end function



      function ccoef(lc,jc,rc,nc,l,j,r,n,jtot,kr,sn,lambda,model)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  subroutine for coupling coefficients:
c     <lc,jc,rc,I,M_i|Plambda(costheta)|l,j,r,I,M_i>  jc means j'
c 
c  l,sn,j=orbital, spin, and total angular momentum of valence
c  r, rc= initial (rhs), final (lhs) core spin
c  
c  kr=rotational band (projection of r on symmetry axis)
c  lambda=multipolarity
c  n,nc=nb of phonons

c  multiplied by (2j'+1)^(1/2) to be symmetric!!!
c  model=0 rotor
c       =1 vibrator
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	use globals
      implicit none
	real*8:: jc,j,jtot,sn,rc,r
	integer:: lc,l,lambda,model,n,nc
	real*8:: lcr,lr,rcr,rr,lambdar,krr,kr
	real*8:: fact,rme,ph,r6,r3 !,cg
	real*8:: ccoef,sixj,cleb
	real*8::zero,one,pi,vibrator
	parameter(zero=0d0,one=1d0)
	pi=acos(-1d0)
	ph=(-1)**(jtot+2*j+rc+lambda+lc+sn)
	fact=sqrt(2.*r+1)*sqrt(2*j+1)*sqrt(2.*l+1)*sqrt(2*jc+1)
	rr=r
	krr=kr
	lambdar=lambda
	rcr=1.*rc
        r6=sixj(j,rr,jtot,rcr,jc,lambdar)
	ccoef=ph*fact*r6
	lr=1.*l
	lcr=1.*lc
	rme=cleb(lambdar,zero,lr,zero,lcr,zero)
     &     *sixj(jc,j,lambdar,lr,lcr,sn)
	ccoef=ccoef*rme
      select case (model)
      case(0) ! rotor
       r3=cleb(rr,krr,lambdar,zero,rcr,krr)
       ccoef=ccoef*r3
!       write(*,*) j,rr,jtot,rcr,jc,lambdar,r6
      case(1) ! vibrator
!       if (rr.eq.rcr) ccoef=0d0
!       if (abs(ccoef).gt.0.) write(*,*)'vib=',ccoef
       ccoef=ccoef*vibrator(n,r,nc,rc,lambda)
      case default
       write(*,*)'Core model',model,' not used!'
       stop
      end select
	end

c Reduced matrix elements between phonon states
c
c < ni ji || Q(lam) || nf jf >
c
c in the vibrational model
c ni,nf= nb of phonons for initial/final states
c ji,jf= initial/final angular momentum
c lam  = multipolarity
c beta  is included in the radial part, so it is omitted here
c
c (uses Brink&Satchler convention for r.m.e.)
      function vibrator(ni,ji,nf,jf,lam)
      implicit none
      real*8 vibrator,ji,jf,rme,beta,lamr
      logical fail3
      integer ni,nf,lam
      vibrator=0.
      lamr=lam
     
      if (ni.eq.nf) return
      if (fail3(ji,lamr,jf)) return
      if ((ni.eq.1).and.(nf.eq.0)) rme=1./sqrt(2.*lamr+1)
      if ((nf.eq.1).and.(ni.eq.0)) rme=(-1)**lamr
!      write(*,*) 'vib:',ni,nf,ji,jf,lam,rme
      vibrator=rme
      end function
      



c spin.spin (borrowed from Fresco frxx4.f)
      function coefss(l,s1,j,s2,lp,s1p,jp,s2p,jt)
      IMPLICIT NONE
      integer:: l,lp,is,ns
      real*8:: s1,s2,s1p,s2p,j,jp,jt, coefss
      real*8:: z,rac,t,smin,smax,j2,t1,t2,s
      parameter(z=0.0)
      coefss=0d0
      
      IF(L.NE.LP) RETURN
      SMIN = MAX(ABS(L-JT),ABS(S1-S2))
      SMAX = MIN(    L+JT ,    S1+S2 )
      T = 0.0
!      DO 90 S=SMIN,SMAX
      NS=NINT(SMAX-SMIN)
      DO 90 IS=0,NS
      S=SMIN+IS
      T1 = SQRT((2*J+1)*(2*JP+1))*(2*S+1)
     X           * RAC(L+Z, S1,JT,S2,J,S)
     X           * RAC(LP+Z,S1,JT,S2,JP,S)
      T2 = (S*(S+1) - S1*(S1+1) - S2*(S2+1))*0.5
      T = T + T1 * T2 
 90   CONTINUE
      COEFSS = T
      RETURN
 95   RETURN
      END

c tensor (borrowed from Fresco frxx4.f)
c Fresco-> THOx
c racah -> RAC
c WIG9J -> U9
c cleb6 -> cleb ???
      function ctens(l,s1,j,s2,lp,s1p,jp,s2p,jt)
      IMPLICIT NONE
      integer:: l,lp,is,isp,ns,nsp
      real*8:: ctens
      real*8:: s1,s2,s1p,s2p,j,jp,jt, coefss
      real*8:: s,smin,smax,spmin,spmax,sp ! CHECK!!
      real*8:: z,rac,u9,cleb
      real*8:: t,j2,t1,t2,t3,t4,t5
      parameter(z=0.0)
      ctens=0d0
      SMIN = MAX(ABS(L-JT),ABS(S1-S2))
      SMAX = MIN(    L+JT ,    S1+S2 )
      T = 0.0
!      DO 80 S=SMIN,SMAX
      NS=NINT(SMAX-SMIN)
      DO 80 IS=0,NS
      S=SMIN+IS
      T1 = SQRT((2*J+1)*(2*S+1)) * RAC(L+Z, S1,JT,S2,J,S)
        SPMIN= MAX(ABS(LP-JT),ABS(S1-S2),ABS(S-2))
        SPMAX= MIN(    LP+JT ,    S1+S2 ,    S+2 )
!      DO 80 SP=SPMIN,SPMAX
      NSP=NINT(SPMAX-SPMIN)
      DO 80 ISP=0,NSP
      SP=SPMIN+ISP
      T2 = SQRT((2*JP+1)*(2*SP+1)) * RAC(LP+Z,S1,JT,S2, JP,SP)
      T3 = SQRT((2*L+1.)*(2*S+1.)) * (-1)**NINT(JT-L-SP)
     &      * RAC(L+Z,LP+Z,S,SP,2+Z,JT)
      T4 = SQRT(2*(2*LP+1)/(3.*(2*L+1.))) *CLEB(LP+Z,Z,2+Z,Z,L+Z,Z)
      T5 = SQRT((2*SP+1)*(2*2+1)*(2*S1+1)*(2*S2+1))
     &      * U9(S,   SP,   2+Z,
     &              S1,  S1,   1+Z,
     &              S2,  S2,   1+Z)
     &      * SQRT(S1*(S1+1) * S2*(S2+1))
      T = T + T1 * T2 * T3 * T4 * T5
 80    CONTINUE
!      write(*,*)'TENSOR: type,L,S1,J,S2,JT,LP,S1P,JP,S2P,LAM',
!     & L,S1,J,S2,JT,LP,S1P,JP,S2P,LAM
      CTENS = T
      RETURN
      END


c Core spin-orbit (OJO 2L.S) 
c  < (L,S1)J, S2; JT  /  2 l.s2  / (LP,S1P)JP, S2P; JT>
C           K1    K2                        K1p     K2p     projections
c
c rac(ria,rib,ric,rid,rie,rif)= w(abcd;ef)  (B&S)

      function corels(l,s1,j,s2,lp,s1p,jp,s2p,jt)
       implicit none
	integer:: l,lp,nj2,ij2
	real*8:: s1,s2,s1p,s2p,j,jp,jt, corels
	real*8:: z,rac,t,j2min,j2max,j2
       parameter(z=0.0)
       corels=0d0
       if (s1.ne.s1p) return
       if (l.ne.lp) return

      J2MIN = MAX(ABS(L-S2),ABS(S1-JT),ABS(LP-S2))
      J2MAX = MIN(    L+S2 ,    S1+JT ,    LP+S2 )
      T = 0.0
      NJ2=NINT(J2MAX-J2MIN)
      DO  IJ2=0,NJ2
      J2=J2MIN+IJ2
      T = T + (2*J2+1) * RAC(S1,L+Z ,JT,S2,J,J2)
     &                 * RAC(S1,LP+Z,JT,S2,JP,J2)
     &                 * (J2*(J2+1) - L*(L+1.) - S2*(S2+1))
!      write(300,'(5f8.4)') j2,RAC(S1,L+Z ,JT,S2,J,J2),
!     & RAC(S1,LP+Z,JT,S2,JP,J2),(J2*(J2+1) - L*(L+1.) - S2*(S2+1))
      ENDDO
      corels = T * SQRT((2*JP+1)*(2*J+1)) * (-1)**NINT(J-JP+LP-L)
      corels=corels/2.d0
      return
      end function 




c Build and diagonalize full Hamiltonian
      subroutine fullham(nset,nho)      
      use wfs
      use channels
      use sistema
      use forbidden
      use hmatrix
      use globals
      use trace
      use constants, only: hc
      implicit none
      integer:: i,ir,m,n,ip,n1,m1,inchani,inchanf,ninc
      integer:: incn,np,nho,nset
      real*8:: norm,r,norm2,rms,rl,rmstot,ex,excore
      real*8,parameter:: uno=1.0
      real*8 :: wcut(1:maxchan)
      real*8,allocatable:: faux(:),ebaux(:),vertex(:,:)
      real*8 :: raux,d2wf,deriv2
!      integer,allocatable:: iord(:)
      integer:: ifail
      real*8    , parameter:: zero=0.0
      CHARACTER(LEN=80) formato
c     --------------------------------------------------
      rmstot=0d0
      wcut(1:maxchan)=jpiset(nset)%wcut(1:maxchan)
!      write(*,*)'wcut=',wcut(1:nchan)

      if (nho.eq.0) then 
         write(*,*)'fullham: nho=0!'
         return
      endif
      if ((nsp.eq.0).or.(bas2.eq.1)) nsp=nho
       hdim=int(nchan*nsp)
      if ((nsp<nho).and.(bas2.eq.0)) then
      write(*,242) nho,nsp
242   format(7x,'WARNING: there are',i3, ' sp eigenstates but only',
     & i3,1x,'will be retained for full H')
      endif 
!      write(*,260) hdim,hdim
!260   format(" o Full Hamiltonian has dimension ",i3,"x",i3)

      if (allocated(hmatx)) deallocate(hmatx)
      allocate(hmatx(hdim,hdim))

      if (npauli>0) then
        write(*,*) 
        write(*,*) ' Pauli forbidden states removed with P projector'
        write(*,*) 
	call pauliproj(nsp)
      endif

      write(*,243) hdim,hdim
243   format(//,' o Building full Hamiltonian matrix (',i3,'x',i3,'):'
     &  ,$)
      do inchani=1,nchan
       do inchanf=inchani,nchan
        write(*,'(".",$)')			! Man! this is so cool!
	call hmat2(inchani,inchanf) ! nsp x nsp block

C Valence-Core coupling  already included, constructing full matrix
       do n=1,nsp
        do m=1,nsp
         hmatx(nsp*(inchani-1)+n,nsp*(inchanf-1)+m)=hmat(n,m) 
! this is not the best way, there should be a better way!!!!!
        enddo !m
       enddo !n
       enddo !inchani
      enddo !inchanf


      do n1=1,hdim
        do m1=1,n1
        hmatx(n1,m1)=hmatx(m1,n1)
        enddo
      enddo


      if (npauli>0) then
c     project out Pauli forbidden states using the projectors P=|forbidden><forbidden| and Q=1-P
c     as H -> QHQ + eshift*P  
c to recover hindrance, uncomment c1, c2, c3 (amoro)
      if (.not.hindrance) then
      hmatx=matmul(qpauli,hmatx)
      hmatx=matmul(hmatx,qpauli)        ! I recover qpauli although it would be the identity matrix if hindrance=false
      endif
!      hmatx=hmatx+eshift*ppauli

       hmatx=hmatx+ ppauli
 
      endif

      if (verb.ge.3) then
        call writemat(hmatx,hdim,74)
        written(74)=.true.
      endif
 

c *** Diagonalization 
      ieigen=0
      nrot=0
      if (allocated(eigv)) deallocate(eigv)
      allocate(eigv(hdim,hdim))
      eigv(:,:)=0.d0
      do n=1,hdim
      eigv(n,n)=1.d0
      enddo
      call hdiag(hmatx,hdim,hdim,ieigen,eigv,nrot)

      write(*,*)' '
      write(*,'(//,5x,"*** EIGENVALUES of full Hamiltonian *** ")')
c *** Construct eigenfunctions
      i=0;
      if (allocated(faux))  deallocate(faux)
      if (allocated(wfeig)) deallocate(wfeig)
      if (allocated(wchan)) deallocate(wchan)
      if (allocated(ebin))  deallocate(ebin)
      if (allocated(iord))  deallocate(iord)
      if (allocated(ebaux))  deallocate(ebaux)
      if (allocated(vertex))  deallocate(vertex)

      allocate(faux(nr))
      allocate(wfeig(hdim,nchan,nr))
      allocate(vertex(nchan,nr))

c *** we need to allocate this only for the first j/pi set
c     (changed in version 2.3 to complex)
!     if (.not.allocated(wfr)) then
      if (.not.allocated(wfc)) then
        nchmax=0
        allocate(wfc(jpsets,maxeset,maxchan,nr))
        allocate(energ(jpsets,maxeset))
      endif
!      allocate(idenerg(jpsets,maxeset))

      if (verb.ge.1) then
      written(78)=.true.
      endif
      if (verb.ge.0) then
      written(79)=.true.
      endif
      allocate(wchan(hdim,nchan))
      wchan(:,:)=0d0
      allocate(ebin(hdim),iord(hdim),ebaux(hdim))

!     Let's try ordering by external subroutine
      do i=1,hdim
        ebaux(i)=hmatx(hdim-i+1,hdim-i+1)
      enddo !i
      call m01daf(ebaux,1,hdim,'a',IORD,ifail)
      if (ifail.ne.0) then
         write(*,*)'m01daf: IFAIL=',ifail
      endif
      call m01zaf(IORD,1,hdim,ifail)    
      if (ifail.ne.0) then
         write(*,*)'m01zaf: IFAIL=',ifail
      endif
       
      ninc=0
      do i=1,hdim
	    ebin(i)=ebaux(iord(i))
        ex=ebin(i)
        rmstot=0d0
        n1=hdim-iord(i)+1
        if (verb.ge.1) then
           write(78, '("#",i2," Ex=",f10.5)') i,hmatx(n1,n1)
           if (wfprint(i)) write(194,'(i3,1f12.5)')nho, hmatx(n1,n1)      
           if (wfprint(i)) write(198,'(1f12.5)') hmatx(n1,n1)   
        endif
        if (verb.ge.0) then
           write(79,*)hdim,hmatx(n1,n1),1.2
        endif
        if ((hmatx(n1,n1)<0).or.wfprint(i)) then   
        write(100+i,246)nchan,i,hmatx(n1,n1)
!        write(200+i,246)nchan,i,hmatx(n1,n1)

246	   format('# ',i2,' Channels, Eigenvalue:',i2,' Energy: ',f8.4)
          write(100+i,247)jtot,partot
!          write(200+i,247)jtot,partot
247	   format('# J:',f4.2,' parity:',i2)
          write(100+i,*) '#  channel core Jc (l sn) jn'
!           write(200+i,*)'#  channel core Jc (l sn) jn'
        endif
 
        ip=0
        do inchanf=1,nchan
	 if (bas2.eq.0) then
	     wfaux=>wfsp(spindex(inchanf),:,:)
	 else if (bas2.eq.1) then
	     wfaux=>wftho(:,ql(inchanf),:)
	 endif
        faux=0d0
        norm=0d0
        do m=1,nsp
           ip=ip+1
           norm=norm+eigv(n1,ip)**2
           m1=hdim-ip+1
           faux(:)=faux(:)+eigv(n1,ip)*wfaux(m,:)
        enddo ! m
        wfeig(i,inchanf,:)=faux(:)

        call normfun(faux,faux,dr,rmin,nr,1,norm2,rms,rl)
        rmstot=rmstot+rms**2
        wchan(i,inchanf)=norm
        if ((hmatx(n1,n1)<0).or.wfprint(i)) then
        written(100+i)=.true.
        incn=inchanf
        write(100+i,245)incn,cindex(incn),qjc(incn),ql(incn),sn,qj(incn)
245	 format('# ',i2,' :',i2,f4.1,i2,f4.1,f4.1)
 
	 endif
        if (verb.ge.1) then
        write(78,'(5x,"Chan ",i2,":",f10.5,1x,"=",100f7.4)') inchanf,
     &      norm,(eigv(n1,(inchanf-1)*nsp+m)**2,m=1,nsp)
        endif
        enddo ! inchanf

        wfeig(i,:,:)=wfeig(i,:,:)*sign(uno,wfeig(i,1,1))

        ex=ebin(i)
c store wfs and energies for all j/pi sets
        if((ex.lt.exmin).or.(ex.gt.exmax)) cycle
        if (any(wchan(i,1:nchan).lt.wcut(1:nchan))) then   ! Added dec/16.
          print*,'i=',i, 'excluded due to small weight'; 
          cycle
        endif
          ninc=ninc+1
          energ(nset,ninc)=ex
          nst=nst +1 ! global counter
c changed in version 2.3 for complex bins
!         wfr(nset,ninc,1:nchan,:)=wfeig(i,1:nchan,:)
          wfc(nset,ninc,1:nchan,:)=wfeig(i,1:nchan,:)
!          idenerg(nset,ninc)=i
        
        write(formato,'( "(5x,",'' "#" '',",i3,1f10.4,2x,",
     & '' "[", '',I0,"f9.5",'' ,"]" '',",4x,",    
     & '' "rms=" '',",1f8.4)"
     & )' ) NCHAN
        rmstot=sqrt(rmstot)

       write(*,formato) i, hmatx(n1,n1), (wchan(i,m),m=1,nchan),rmstot
!       endif ! exmin < ex < exmax
       enddo !i=1,hdim (eigenvalues)

       jpiset(nset)%nex=ninc
c this is included in basis.f90, so it should be redundant here
       if (nchan.gt.nchmax) nchmax=nchan
       write(*,'(/,5x,"[",i3, " out of", i4,
     & " eigenvalues retained ]")')  ninc,hdim
     
     
!       if (any(froverlap(1:20)>0)>0) open(300,file='overlaps.fr')
       if (froverlaps>0) open(300,file='overlaps.fr')


! Write eigenfunctions 
        do i=1,hdim
        if (written(100+i).or.froverlaps>0) then
        if((ex.lt.exmin).or.(ex.gt.exmax)) cycle
        if (any(wchan(i,1:nchan).lt.wcut(1:nchan))) cycle
        if ((rlast.gt.0).and.(rlast.lt.rmax)) then 
            np=idint((rlast-rmin)/dr)
        else
            np=idint((rmax-rmin)/dr)
        endif
        if (written(100+i)) write(100+i,248)np
248     format("# ",i5," points")

! Compute & write vertex functions (march 2019) 
! Modified by Pedro (march 2022)
       do m=1,nchan
       excore =jpiset(nset)%exc(m) ! core energy
       
       do ir=1,np
        r=rmin+dr*dble(ir-1)
        !if (r.lt.1e-4) r=1e-4
        faux(ir)=r*wfeig(i,m,ir)
       enddo !ir
             
       do ir=1,np
       r=rmin+dr*dble(ir-1)
       if (r.lt.1e-5) then !Avoid indeterminate form
       if (ql(m).eq.0) then !l=0
       vertex(m,ir)=(ebin(i)-excore)*wfeig(i,m,ir)
     &   +(hc**2/2/mu12)*6*(wfeig(i,m,2)-wfeig(i,m,1))/(2*r+dr)/dr
!     &   +(hc**2/2/mu12)*2*(wfeig(i,m,ir+2)-wfeig(i,m,ir+1))/dr/dr
       else !l>0
       vertex(m,ir)=0
       endif 
       else !r>0
       d2wf=deriv2(faux,dr,np,ir)
!       write(99,*)rmin+dr*dble(ir-1),faux(ir),d2wf
       vertex(m,ir)=(ebin(i)-excore)*wfeig(i,m,ir)
     &   +(hc**2/2/mu12)*(d2wf/r-ql(m)*(ql(m)+1.)*wfeig(i,m,ir)/r/r)
       endif
       enddo !ir
       enddo !m
       if (written(100+i)) then      
        do ir=1,np
        r=rmin+dr*dble(ir-1)
        if (r> rlast) cycle
        write(100+i,'(1f8.3,2x,10g14.6)') r,(r*wfeig(i,m,ir),m=1,nchan)   
!        write(200+i,'(1f8.3,2x,10g14.6)') r,(vertex(m,ir),m=1,nchan)
        enddo !ir
        write(100+i,*)'& '
!        write(200+i,*)'& '



        if (verb.ge.4) then
        do ir=1,np
        r=rmin+dr*dble(ir-1)
        if (r> rlast) cycle
        write(100+i,'(1f8.3,2x,10g14.6)') r,(r*vertex(m,ir),m=1,nchan)
        enddo !ir
        write(100+i,*)'& '
        endif
!        write(200+i,*)'& '
	endif 
        endif !written(100+i)

! Write overlaps & vertex functions in fresco format
        if(froverlaps>0) then
        ex=ebin(i)
        write(0,*)'print overlap for i=',i, 'ex=',ebin(i)
        if((ex.lt.exmin).or.(ex.gt.exmax)) cycle
        if (any(wchan(i,1:nchan).lt.wcut(1:nchan))) then   
          print*,'i=',i, 'excluded due to small weight'; 
          cycle
        endif        
        
        do m=1,nchan
        if (froverlaps==1 .or. ebin(i).lt.0) then 
        write(300, '("#Single particle REAL wf & vertex for state",i2,
     &  2x,"Chan=",i3," E=",f10.5," MeV")') i, m,ebin(i)
        else ! complex
        write(300,'("#Single particle COMPLEX wf & vertex for state",i2,
     &  2x,"Chan=",i3," E=",f10.5," MeV")') i,m, ebin(i)
        endif
     
     
        write(300,*) np,dr,rmin
!        do m=1,nchan
        if (froverlaps==1 .or. ebin(i).lt.0) then ! real
          write(300,260) (wfeig(i,m,ir),ir=1,np)   
          write(300,260) (vertex(m,ir),ir=1,np)
        else                   ! complex   
          write(300,260) (cmplx(wfeig(i,m,ir),zero),ir=1,np)   
          write(300,260) (cmplx(vertex(m,ir),zero) ,ir=1,np)
        endif

        enddo ! m        
        endif !print overlaps

260    FORMAT(1P,6E12.4)
262    FORMAT(1P,12E12.4)
 

       enddo !i 
       call flush(6)


	   
! Write CDCC wfs	   
       if (xcdcc) call writewf(hdim,np)
       
c Free memory
      deallocate(hmatx,hmat)
      if (allocated(ortmat)) deallocate(ortmat)
      if (allocated(wftho)) deallocate(wftho)
	  end subroutine



       subroutine orthonorm
        use wfs, only: rmin,dr,nr,wfeig,hdim
        use channels, only: nchan
        use globals, only: written
        implicit none
        real*8,pointer:: faux(:),gaux(:)
        real*8:: norm,rms,rl
        real*8:: onm(hdim,hdim)
        integer:: i,j,k
        onm(:,:)=0d0
        write(*,*)
        write(*,*)' o Checking Orthogonality, see fort.90'
        write(90,*)'Othogonality matrix for final eigenfunctions'
        write(90,*)'Test for ensure the final orthogonality'
        do i=1,hdim
         do k=1,hdim
          do j=1,nchan
           faux=>wfeig(i,j,:)
           gaux=>wfeig(k,j,:)
           call normfun(faux,gaux,dr,rmin,nr,1,norm,rms,rl)
           onm(i,k)=onm(i,k)+norm  !we sum over the channels
          enddo
         enddo
        write(90,'(100f12.5)')(onm(i,k),k=1,hdim)
        enddo
        written(90)=.true.

        end 


c  --------------------------------
c  write wfs for XCDCC calculations
c  --------------------------------
      subroutine writewf(ndim,np)
      use hmatrix,only: hmatx
      use globals
      use sistema
      use wfs
      use channels
      implicit none
      integer kwf,ndim,n,ich,ir,m,ninc,np
      real*8 r,einc(ndim),ex
      parameter(kwf=84)
      written(kwf)=.true.

      ninc=0
      do n=1,ndim
      ex=hmatx(ndim-n+1,ndim-n+1)
      if((ex>exmin).and.(ex<exmax)) then
        ninc=ninc+1
        einc(ninc)=ex
      endif
      enddo

      write(kwf,230) jtot, partot,ninc, nchan, np, dr, rmin

230   format(1f4.1,2i3,i3,3x,i6,f6.3,f6.3)
240   format(5f10.6)
      write(kwf,*) '# Ic  (l  sn)  j '
      do ich=1,nchan
      write(kwf,245) qjc(ich),ql(ich),sn,qj(ich)
245   format(3x,f4.1,i2,f4.1,f4.1)
      enddo
c      u(r)
!      write(99,*) 'writewf: ndim,exmin,exmax=',ndim,exmin,exmax
      do n=1,ndim
!          ex=hmatx(ndim-n+1,ndim-n+1)
      	  ex=ebin(n)
          if((ex<exmin).or.(ex>exmax)) cycle
          write(kwf,*)  ex
	   do ir=1,np ! changed in 2.2 from nr to np
          r=rmin+dr*dble(ir-1)
          write(kwf,'(3x,10g14.6)') (r*wfeig(n,m,ir),m=1,nchan)       
	   enddo ! ir
      enddo ! n
      end subroutine
      
c  ----------------------------------------
c  write vertex functions for a given jpset
!  PS wfs from wfc(jset,i,ichan,ir)
!  < a | V | a'> is stored in  ccmat(1:nchan,1:nchan,1:nr)
c  -----------------------------------------
      subroutine vertex(jset)
      use wfs, only: energ,wfc,rmax,rmin,dr,rlast !,idenerg
      use potentials, only: ccmat
      use channels, only:  jpiset
      use globals, only: mu12
      use constants, only: hc
      implicit none
      integer:: nch,jset,ch,np,ni,i,j,ir
      character (len=9):: file_name

      ni=jpiset(jset)%nex
      nch=jpiset(jset)%nchan
      call coefmat(jset,nch)
      ccmat=ccmat*hc**2/mu12/2

      ch=1  !ONlY CHANEL 1

      if ((rlast.gt.0).and.(rlast.lt.rmax)) then 
            np=idint((rlast-rmin)/dr)
        else
            np=idint((rmax-rmin)/dr)
        endif
      
      write(file_name,'("vertex.",i0)') jset 
      open(525,file=trim(file_name))
      
      !wfc are converted to real
      do i=1, ni  !wfc states for this jset
      do j=1, nch !output chanel
      if (ch.GT.0 .and. ch.NE.j) cycle !if a chanel is selected
       write(525, '("#Single particle REAL wf & vertex for state",i2,2x,
     &  "Chan=",i3," E=",f10.5," MeV")') i,j,energ(jset,i)
        write(525,*)np, dr, rmin
        write(525,'(1P 6E12.4)') (dble(wfc(jset,i,j,ir)), ir=1,np)
        write(525,'(1P 6E12.4)') 
     & (dble(DOT_PRODUCT(ccmat(j,:,ir),wfc(jset,i,1:nch,ir))), ir=1,np)     
    
      enddo !j
      enddo !i
                 
      close(525)
      deallocate(ccmat)     

      end subroutine
      
      


