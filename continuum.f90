c -------------------------------------------------------------------------
c
c *** Angular matrix elements < (ls)j, I; J | V  | (ls)j, I; J>*CONV*phase
c 
c     V includes the radial part and any possible tensor part (eg L.S)
c     CONV=2*mu/h^2
c
c     phase=(-1)**NINT((JCI-JCF - ABS(JCI-JCF))/2.)
C           The above phase factors with JCORE etc
C           are there only because of definition of M(Ek) matrix element
c -------------------------------------------------------------------------
      subroutine coefmat(nchan)
      use wfs, only: nr,wfsp,rvec
      use channels,only:jtot,qj,spindex,ql,qj,qjc,cindex,nphon,sn,partot
      use sistema
      use globals
      use potentials
      use constants
      use forbidden ! added in v2.2f
      implicit none
      logical, save:: first=.true.
      integer:: ir      
      integer  ici,icf,il
      integer  nchani,nchanf,nchan,li,lf,nphi,nphf
c     -----------------------------------------------------------------------
      real*8:: fact,ji,jf,xl,jci,jcf,coefss,corels
      real*8:: all,ass,als,alsc,allp,ccoef,vcp,rm,conv,faux(nr),zero
      real*8:: lambdar,sixj,cleb,lfr,lir,cl
      real*8:: allp2, gcoup,big,small
c     -----------------------------------------------------------------------
c     Pauli blocking 
      integer ispi,ip
      integer:: lp,np
      real*8:: jp
      real*8,allocatable:: wfpau(:)       
c     -----------------------------------------------------------------------
!      big=huge(big)   
!      small=epsilon(small)
      big = huge(big)**0.8d0
      small=1./big
      zero=0d0

      if (npauli>0) then 
        write(*,*)npauli,' states to be removed!'
        if (allocated(wfpau)) deallocate(wfpau)
        allocate(wfpau(nr))
      endif

      allocate (ccmat(nchan,nchan,nr))
      ccmat(1:nchan,1:nchan,1:nr)=0d0
      fact=hc*hc/2d0/mu12
      rm=av*ac/(av+ac)
      conv=2*amu*rm/hc**2

      if (first) then
!        write(*,*)'coefmat: ql=',ql(1:nchan)
!        write(*,*)'coefmat: qj=',qj(1:nchan)
!        write(*,*)'coefmat: qjc=',qjc(1:nchan)
!        write(*,*)'coefmat: laminc=',laminc(:)
!         write(*,*)'coefmat: small=',small
!         write(*,*)'coefmat:conv=',conv
      endif 

      do nchani=1,nchan
      do nchanf=1,nchan
       als=0d0
       faux(1:nr)=0d0
       li=ql(nchani);  xl=li
       lf=ql(nchanf)

       ji=qj(nchani)
       jf=qj(nchanf)

       jci=qjc(nchani)
       jcf=qjc(nchanf)
        
       ici=cindex(nchani)
       icf=cindex(nchanf)

       nphi=nphon(ici)
       nphf=nphon(icf)

       if (nchani.eq.nchanf) then
          als=0.5*(ji*(ji+1)-li*(li+1.)-sn*(sn+1)) ! l.s
          all=xl*(xl+1.)                           ! l.l

          faux(1:nr)=
     &       vcou(1:nr)        ! coulomb
     &     + vcl(li,1:nr)      ! nuclear central
     &     + vls(1:nr)*als     ! spin-orbit
     &     + vll(1:nr)*all     ! l.l
 
!          do ir=1,nr
!            write(225,*)rvec(ir),faux(ir) ,vcl(li,ir)
!          enddo
        endif

c *** Pauli blocking ---------------------------------------------
c      if ((npauli.gt.0).and.(nchani.eq.nchanf).and.(li.eq.lf)
c     &   .and.(ji.eq.jf).and.(ici.eq.icf)) then
c        ispi=spindex(nchani)
c        write(*,*)'pauli: ispi,npauli=',ispi,npauli
c        do ip=1,npauli
c	lp=paul(ip) 
c        jp=pauj(ip)
c	np=paun(ip)
c        if (wfpread(ip)) then
c          wfpau(:)=wfpauli(ip,:)
c        else
c          wfpau(:)=wfsp(ispi,paun(ip),:) !Pauli forbidden 
c        endif
c        faux(1:nr)=faux(1:nr)+pshift(ip,ici)*abs(wfpau(1:nr))**2
c        enddo !ip
c      endif ! Pauli?
c -----------------------------------------------------------------


c spin.spin term 
      if (li.eq.lf) then
        ass= coefss(li,sn,ji,jci,lf,sn,jf,jcf,jtot)
        faux(1:nr)=faux(1:nr)+ ass*vss(li,1:nr)
      endif

c spin-orbit for core (added to v2.0.2 by AMoro)
      alsc=corels(li,sn,ji,qjc(nchani),lf,sn,jf,qjc(nchanf),jtot)
      faux(1:nr)=faux(1:nr)+ alsc*vlsc(1:nr)


c add particle-core coupling for each lambda
       do il=1,maxlamb
      if (.not.laminc(il)) cycle  
! matrix element of P(theta)
c	allp2=gcoup(ql(nchanf),qj(nchanf),qjc(nchanf),nphf,
c     &            ql(nchani),qj(nchani),qjc(nchani),nphi,
c     &            jtot,kband,sn,il,pcmodel)
c! matrix element of Y(theta)
c       allp2=allp2*sqrt(2d0*il+1)/sqrt(4*pi)
        
	allp=gcoup(li,ji,jci,nphi,lf,jf,jcf,nphf,
     &             jtot,kband,sn,il,pcmodel)
       allp=allp*sqrt(2d0*il+1)/sqrt(4*pi)

	allp2=ccoef(lf,jf,jcf,nphf,
     &              li,ji,jci,nphi,
     &            jtot,kband,sn,il,pcmodel)

c matrix element of Y(theta)
       allp2=allp2*sqrt(2d0*il+1)/sqrt(4*pi)

!       if (abs(allp).gt.1e-5) then 
!         write(99,'("chani=",i2," l,j,jci,nphi=",i2,2f4.1,i3,
!    &   " chanf=",i2," l,j,jc,nphf=",i2,2f4.1,i3
!     &   " lambda=",i2," sn=",1f3.1," il=",i3, 
!     &   "= > Geom=", 2f8.4)') nchani,li,ji,jci,nphi,nchanf,
!     &   lf,jf,jcf,nphf,
!     &   lambda,sn,il,allp,allp2
!       endif
!       all=li*(li+1)

c spin.spin term 
!      ass= coefss(li,sn,ji,jci,lf,sn,jf,jcf,jtot)

    
      do ir=1, nr
	 select case (lpot)
           case (0) ! mininum li,lf
             vcp=vlcoup(il,min(li,lf),ir)
           case (1) ! average
             vcp=(vlcoup(il,min(li,lf),ir) +vlcoup(il,lf,ir))/2d0
           case (2) ! maximum li,lf
              vcp=vlcoup(il,max(li,lf),ir)
           case (3) ! left
              vcp= vlcoup(il,lf,ir)
           case default
              write(*,*) 'lpot=',lpot,'not valid!'; stop
         end select

         if (abs(vcp*allp).lt.small) vcp=0.0 ! avoid underflow

         faux(ir)=faux(ir) + vcp*allp !!!+ ass*vss(ir)
!         write(225,*)ir,vcp*allp,faux(ir)
       enddo ! ir
       enddo !il

! For external coupling potentials
        if (cptype.eq.5) then
           do ir=1, nr
           lfr=1d0*lf
           lir=1d0*li
           do il=0,maxlamb
           lambdar=1d0*il
           if (il.eq.0) then
             cl=sqrt(4*pi)    !maybe we would need to divide
           else
             cl=1d0
           endif

            if (partot.eq.1) then
           faux(ir)=faux(ir)+((-1d0)**(ji+jcf+jtot)*sixj(jf,ji,lambdar
     &,jci,jcf,jtot)*sqrt((2*jf+1.)*(2*ji+1.))*sixj(jf,ji,lambdar,
     &lir,lfr,sn)*(-1d0)**(lambdar+lfr+ji+sn)*sqrt((2.*lf+1)*(2.*il+1.)
     &/4/pi)*cleb(lambdar,0d0,lfr,0d0,lir,0d0)*
     &vtran(ici,icf,il,ir,1))*cl
            else
           faux(ir)=faux(ir)+((-1d0)**(ji+jcf+jtot)*sixj(jf,ji,lambdar
     &,jci,jcf,jtot)*sqrt((2*jf+1.)*(2*ji+1.))*sixj(jf,ji,lambdar,
     &lir,lfr,sn)*(-1d0)**(lambdar+lfr+ji+sn)*sqrt((2.*lf+1)*(2.*il+1.)
     &/4/pi)*cleb(lambdar,0d0,lfr,0d0,lir,0d0)*
     &vtran(ici,icf,il,ir,2))*cl
            endif
!     &vtran(jci,jcf,ici,icf,il,ir))*cl
            enddo !il
        if (abs(faux(ir)).lt.small) faux(ir)=small ! avoid underflow
          enddo ! ir
        endif !cptype=5

!!!! For bincc !!!!!!!!!!!!!!!!!!!!!!!!!
!        ccmat(nchani,nchanf,1:nr)=cmplx(-faux(1:nr)*conv    
!     X      * (-1)**NINT((JCI-JCF - ABS(JCI-JCF))/2.),0d0)

! For "scatcc"    
! v2.1 I have changed the sign of ccmat so they have their physical sign
      do ir=1,nr
       if (abs(faux(ir)).lt.small) faux(ir)=0.0
       ccmat(nchani,nchanf,ir)=cmplx(faux(ir),zero)*conv
       enddo
       
      

!       ccmat(nchani,nchanf,1:nr)=cmplx(faux(1:nr)*conv 
!     X      * (-1)**NINT((JCI-JCF - ABS(JCI-JCF))/2.),0d0)
C           The above phase factors with JCORE etc
C           are there only because of definition of M(Ek) matrix element
!!!!!      
       enddo !inchani
      enddo !inchanf


      first=.false.

      if (verb.lt.4) return 

      do nchani=1,nchan
      write(143,*)'Coefmat for jtot,nchan=',jtot,nchan
      write(143,402) (real(ccmat(nchani,nchanf,1)),nchanf=1,nchan)
402   FORMAT( 1X,1P,14E10.2/(1X,24E9.2))
      enddo 

      end subroutine 



c ------------------------------------------------------------
c ** Calculate scattering states by differential integration
c ------------------------------------------------------------
      subroutine continuum_range
      use globals, only: mu12,egs,kin,written,verb
      use sistema
      use parameters, only: maxchan
      use constants
      use channels, only: jpiset,ql,jpsets
      use wfs     , only: nr,dr,rvec
      implicit none
c     ........................................................
      logical :: energy,ifcont,writewf
      character*5::jpi
      integer:: il, pcon ! backward compatibility
      integer inc,ilout,ili,jset,partot,nchan,nk,n,ir,ik
      real*8 :: emin,emax,eout,jtot,ecm
      real*8 :: kmin,dk,de,kmax,kcont
      real*8 :: r0,conv,rm
      real*8,allocatable :: psh_el(:)
!      real*8 :: deladd,deltai,deltap
      complex*16,allocatable:: smate(:,:),wfcont(:,:,:)
      namelist/scatwf/ emin, emax,ifcont,nk,inc,writewf,ilout,eout,jset,
     &                 energy,
     &                 il,pcon ! not used, backward compatibility only (AMM)

      inc=0; il=0
      jset=0; pcon=0
  
      writewf=.false.
      energy=.false.
      read(kin,nml=scatwf)
      if ((.not.ifcont).or.(nk.eq.0)) return
      if ((jset.gt.jpsets).or.(jset.eq.0)) then
       write(*,*)'jset',jset,' not valid!'
       return
      endif
      if (il.gt.0) then
       write(*,*)' IL deprecated; use INC instead!'
      endif


      partot=jpiset(jset)%partot
      jtot  =jpiset(jset)%jtot
      nchan =jpiset(jset)%nchan
      ql (1:nchan)   =jpiset(jset)%lsp(1:nchan)
      r0=rvec(1)
!      write(50,400) jpi(jtot,partot), nchan,inc

      rm=av*ac/(av+ac)
      conv=(2*amu/hc**2)*rm

      kmin=sqrt(conv*abs(emin))
      kmax=sqrt(conv*abs(emax))

      if (energy) then
        de=(emax-emin)/(nk-1)
      else
        dk=(kmax-kmin)/(nk-1)
      endif


      if(allocated(smate)) deallocate(smate)
      allocate(smate(nk,maxchan),psh_el(nk))
      if (allocated(wfcont)) deallocate(wfcont)
      allocate(wfcont(nk,nchan,nr))
!      deladd=0.
      call wfrange(jset,nchan,inc,emin,emax,nk,energy,wfcont,
     &     smate,psh_el)
     

c ... Write WFS
      if (writewf) then
      write(50,400) jpi(jtot,partot), nchan,inc,energy !, ql(1:nchan)
400   format("# Continuum WFS for Jp=",a5,"; Channels=",i2,
     & ";  Inc. chan.=",i2, ' Uniform Ecm=',l1)
      write(50,'(1x,i4,2x,2f8.4)') nr,dr,r0
      if (energy) then
       write(50,'(1x,i4,2x,2f8.4)') nk,de,emin
      else
       write(50,'(1x,i4,2x,2f8.4)') nk,dk,kmin
      endif

      do ik=1,nk
      if (energy) then
        ecm=emin+(ik-1)*de
      else
        kcont=kmin+ (ik-1)*dk  !new kcont for schcc
        ecm=(hc*kcont)**2/2/mu12
      endif
      write(50,*)'# Continuum wf with Ecm=',Ecm
      do ir=1,nr
      write(50,'(1f8.3,2x,50f12.8)')rvec(ir),
     &       (wfcont(ik,n,ir),n=1,nchan)
      enddo ! ir
      write(50,*)'& '
      enddo ! ik
      endif !writewf
      
      end subroutine 

c ------------------------------------------------------------
c ** Calculate scattering states by differential integration
c ------------------------------------------------------------
      subroutine continuum!(nchan)
      use globals, only: mu12,egs,kin,written,verb
      use constants
      use sistema
      use scattering
      use wfs, only:nr,dr,rvec
      use channels
      use parameters, only: maxchan
      use potentials, only: ccmat
      use factorials
!!!! scatcc
      use nmrv,only:nch,ech,vcoup
!!!!!
      implicit none
      logical info,writewf
      integer n,ir,ich,ichp,nskip,klog,method
      real*8:: r,r0,dk,kcont,econt,k,ecm, t1,t2
      real*8:: phase(nchan),aux,phadd(nchan)
      real*8,allocatable:: ph(:,:)
c  bincc variables
      logical tres,tdel,tkmat,tknrm
      integer:: nf   ! number of formfactors?
      integer:: isc,iil,ik
      integer:: maxn,nnk(nchan)
      real*8:: anc,k2(nchan),conv,rm,bphase(nk)
      complex*16:: y(nr,nchan),etap(nchan),wf(nchan,nr),smat(nchan)
      real*8:: fmscal
c  end bincc 
      namelist/scatwf/ emin, emax,ifcont,nk,il,ilout,writewf
      il=0
      ili=1
      klog=99
      read(kin,nml=scatwf)

      if ((.not.ifcont).or.(nk.eq.0)) return

      write(*,*)
      write(*,*)' SCATTERING WAVE FUNCTIONS:'
      write(*,*)
      if (il.eq.0) then
        write(*,*)'Max. incoming channel not specified: assuming nchan'
        il=nchan
      else 
        if (il.gt.nchan) then
        write(*,*)' # incoming channels bigger than #outcoming: 
     &assuming IL=nchan'
        il=nchan
      endif        

      if (il.lt.0) then
        write(*,*)' Max. incoming channel negative: computing only -il'
        il=-il
        ili=il
      endif
      endif 

      write(*,*)
      write(*,100)nk,emin,emax
100   format(3x,"Generating", i4, " continuum wf(s) for"
     &  ," [Emin=",1f6.3," Emax=",1f7.2," MeV]")

      call factorialgen(100)

      call coefmat(nchan)

      allocate(wfcont(nk,nchan,nchan,nr))
      wfcont(:,:,:,:)=0d0

      rm=av*ac/(av+ac)
      conv=(2*amu/hc**2)*rm
      kmin=sqrt(2d0*mu12*abs(emin))/hc
      kmax=sqrt(2d0*mu12*abs(emax))/hc


!! ! for bincc --------------------------------------------------
c1      nf=1;  isc=2      
c1      tres=.false.; tdel=.true.; tkmat=.false.; tknrm=.false.
c1      maxn=nr+1


c1      do n=1,nchan
c1       k2(n)=abs(egs)+exc(n) ! separation energy for each channel
c1       etap(n)= 2*e2*zv*zc*rm*amu/hc**2  !=2 k eta !! CHECK!!!!
c1       k2(n)=k2(n)*conv     
c1      enddo      
c1      do iil=ili,il
c1     call BINCC2(Y,ccmat,nf,NCHAN,k2,iil,conv,bphase,isc,kmin,kmax,nk,
c1     & etap,NR-1,dr,pcon,tres,tdel,tkmat,tknrm,lmax+1,maxn,anc)
c1      write(*,*)'il=',iil,' done'
c1      enddo
c1      RETURN !!!!!!! TEST BINCC 
c ----------------------------- end bincc call ---------------------
 

      allocate(vcoup(nchan,nchan,1:nr))
      allocate(ech(nchan))
      allocate(ph(nk,nchan))
      ech(1:nchan)=exc(1:nchan)
      vcoup(:,:,:)=0d0 
      vcoup(:,:,1:nr)=ccmat(:,:,1:nr)/conv ! changed in v2.2

      deallocate(ccmat)
      if (nk.gt.1) then
         DK = (KMAX - KMIN)/(NK-1)
      else
         DK=(KMAX - KMIN)
      endif
      do iil=ili,il
      nskip=0
      phadd=0
      write(*,'(5x,"o incoming channel:",i3,$)') iil
      call cpu_time(t1)
      do ik=1,nk
       K = KMIN + (IK-1)*DK
       ecm=k**2/conv-ech(iil)
!       write(*,'(".",$)')
       if (verb.ge.5) then 
        write(*,*)' **********************************'
        write(*,150) k**2/conv,ecm,il
150    format(4x,'Calling SCHCC for E=',1f7.3,' and Ecm=',1f7.3, 
     &  ' MeV with incident channel',i2)
       endif
!       call schcc(nchan,ecm,zv*zc,iil,ql,conv,dr,nr-1,wf,phase,smat) ! changed in v2.2
!       wfcont(ik,iil,1:nchan,2:nr)=wf(1:nchan,1:nr-1)
       r0=rvec(1)
       info=.false.
       method=4
c       call schcc(nchan,ecm,zv*zc,iil,ql,conv,dr,r0,nr,wf,phase,smat,
c     &            info)

       call schcc_erwin(nchan,ecm,zc*zv,iil,ql,conv,dr,r0,
     &                    nr,wf,phase,smat,method,info)

       wfcont(ik,iil,1:nchan,1:nr)=wf(1:nchan,1:nr)

!      write(47,'(1f10.3,3x,10g14.6)') ecm,(phase(n),n=1,nchan)
!      ph(ik,1:nchan)=phase(1:nchan)
!      write(46,'(1f10.3,3x,10g14.6)') ecm,(ph(ik,n),n=1,nchan) 


       if (ecm.gt.0) then
         written(45)=.true.
         if (ik.eq.1) then
          ph(1,1:nchan)=phase(1:nchan)
         else
          do n=1,nchan 
           aux=phase(n)+phadd(n)
           if (aux< ph(ik-1,n)-90) phadd(n)=phadd(n)+180
           if (aux> ph(ik-1,n)+90) phadd(n)=phadd(n)-180
           ph(ik,n)=phase(n)+phadd(n)
           enddo
         endif 
         write(45,'(1f10.3,3x,10g14.6)') ecm,
     &  (ph(ik,n),n=1,nchan) 


         if (iil.eq.ilout) then 
         written(50)=.true.
         write(50,*)'# Continuum wf with Ecm=',Ecm
         do ir=1,nr
         write(50,'(1f8.3,2x,50f12.8)')rvec(ir),
     &       (wfcont(ik,iil,n,ir),n=1,nchan)
         enddo 
         write(50,*)'& '
         endif
      else
        nskip=nskip+1
        if (nskip.le.5) then
        write(*,200) ecm
200     format(4x,'Skipping negative energy:',1f7.3, ' MeV')
        else if (nskip.eq.6) then
        write(*,*)"    (further not shown)"
        endif
      endif
      enddo ! ik
      call cpu_time(t2)
      write(*,'(3x,"(",1f6.3," secs )")') t2-t1
      enddo ! iil
      write(*,*)




!! AMoro: skip creepy part...
      RETURN 


        do iil=ili,il

         do ik=1,nk
          kcont=kmin+(kmax-kmin)*(ik-1)/(nk-1)
          econt=(hc*kcont)**2/2/mu12
! 	  if ((econt.gt.eout).and.(.not.written(51))) then
!          written(51)=.true.
!          write(51,*)'# Continuum wf with e=',econt,'ik=',ik
!          do ir=1,nr
!      write(51,'(1f8.3,2x,50f10.6)')rvec(ir),
!     &(wfcont(ik,iil,n,ir),n=1,nchan) 
!           enddo
!          endif
         enddo
!         endif
		 enddo



!!------------------------------------------------------------We are going to make the creepiest thing we can imagine to solve our problem with ik
!!------------------------------------------------------------I also use it to print out the function
!      ! first we calculate the number corresponding to the first k upon the barrier
!      dk=(kmax-kmin)/nk
!      do iil=ili,il
!
!         nnk(iil)=nint(sqrt(exc(iil)*conv)/dk)
!
!         do ik=nk,1,-1
!         if (ik.gt.nnk(iil)) then
!
!         wfcont(ik,iil,:,:)=wfcont(ik-nnk(iil),iil,:,:)
!
!         else
!         wfcont(ik,iil,:,:)=cmplx(0d0,0d0)
!         endif
!         enddo
!
!		 if (iil.eq.ilout) then
!         do ik=1,nk
!          kcont=kmin+(kmax-kmin)*ik/nk
!          econt=(hc*kcont)**2/2/mu12
!	  if ((econt.gt.eout).and.(.not.written(50))) then
!          written(50)=.true.
!          write(50,*)'# Continuum wf with e=',econt
!          do ir=1,nr
!      write(50,'(1f8.3,2x,50f10.6)')rvec(ir),
!     &(wfcont(ik,iil,n,ir),n=1,nchan) ! changed in v2.0.5
!!     &(dreal(wfcont(ik,iil,n,ir)),n=1,nchan)
!          enddo
!          endif
!         enddo
!         endif
!
!      enddo
      end subroutine         




