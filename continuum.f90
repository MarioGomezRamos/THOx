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
      subroutine coefmat(nset,nchan)
      use wfs, only: nr,wfsp,rvec
      use channels,only:jtot,qj,spindex,ql,qj,qjc,cindex,nphon,
     &                  sn,partot,qnc,jpiset
      use sistema
      use globals
      use potentials, only:ccmat,vtran,vlcoup,laminc,lpot,pcmodel,kband,
     &                   cptype,maxlamb,vcou,vss,vcl,vls,vlsc,vll,reid93 
      use constants
      use forbidden ! added in v2.2f
      implicit none
      logical, save:: first=.true.
      integer:: ir      
      integer  ici,icf,il
      integer  nchani,nchanf,nchan,li,lf,nphi,nphf,nset
c     -----------------------------------------------------------------------
      real*8:: fact,ji,jf,xl,jci,jcf,coefss,corels,vaux
      real*8:: all,ass,als,alsc,allp,ccoef,vcp,rm,conv,faux(nr),zero
      real*8:: lambdar,sixj,cleb,lfr,lir,cl
      real*8:: allp2, gcoup,big,small
      real*8 :: kbandi,kbandf
c     -----------------------------------------------------------------------
c     Pauli blocking 
      integer ispi,ip
      integer:: lp,np
      real*8:: jp,vscale
      real*8,allocatable:: wfpau(:)       

      character*3  reidname !MGR
      character*1 reidj
      real*8 :: reidpot(2,2),auxreid,vnm,r
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
      vscale= jpiset(nset)%vscale
             
      
!       if (first) then
!        write(*,*)'coefmat: ql=',ql(1:nchan)
!        write(*,*)'coefmat: nset,vscale=',nset,vscale
!        write(*,*)'coefmat: qj=',qj(1:nchan)
!       write(*,*)'coefmat: qjc=',qjc(1:nchan)
!       write(*,*)'coefmat: laminc=',laminc(:)
!         write(*,*)'coefmat: small=',small
!         write(*,*)'coefmat:conv=',conv
!        write(*,*)' partot=',partot
!        first=.false.
!      endif 

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

       kbandi=qnc(ici)%kband
       kbandf=qnc(icf)%kband

       if (kbandi.eq.kbandf) kband=kbandi

       if (nchani.eq.nchanf) then
          als=0.5*(ji*(ji+1)-li*(li+1.)-sn*(sn+1)) ! l.s
          all=xl*(xl+1.)                           ! l.l

          faux(1:nr)=
     &       vcou(1:nr)            ! coulomb
     &     + vcl(li,1:nr)*vscale   ! nuclear central
     &     + vls(li,1:nr)*als      ! spin-orbit MGR l-dependence
     &     + vll(li,1:nr)*all      ! l.l MGR l-dependence
 
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
      faux(1:nr)=faux(1:nr)+ alsc*vlsc(li,1:nr) !MGR l-dependence


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

      if (.not. reid93) then !MGR
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
       endif
       enddo !il
      

      if (reid93) then !MGR
      do ir=1, nr     
      r=rvec(ir)
        if (r.lt.1e-9) r=1e-9 
             vnm=0d0
             write(reidj,'(i1)') nint(jtot)
             reidname(3:3)=reidj
!             if(i.eq.1) write(*,*) 'Reid93 J', reidname(3:3)
!             write(0,*) 'Reid93 J',reidname(3:3),reidj,nint(jtot)
             if (abs(li-nint(jtot)).ne. 0) then !only S=1
             reidname(1:2)='3C'    
             if (nint(jtot).eq.0) reidname(1:2)='3P' 
!            if (ir.eq.1)  write(0,*) 'Reid93 J', reidname(1:3)   
             call rreid93(r,reidname,'PN',reidpot) 
              if (li.gt.lf) then
              vnm=reidpot(2,1)
!              if(i.eq.i) write(*,*) 'Reid93 ', reidname,'l',li,lf,vnm  
              else if (li.lt.lf) then
              vnm=reidpot(1,2)
!              if(i.eq.i) write(*,*) 'Reid93 ', reidname,'l',li,lf,vnm  
              else !li.eq.lf
               if (li.lt.nint(jtot) .or. nint(jtot).eq.0) then
               vnm=reidpot(1,1)
!               if(i.eq.i) write(*,*) 'Reid93 ', reidname,'l',li,lf,vnm  
               else
               vnm=reidpot(2,2)
!               if(i.eq.i) write(*,*) 'Reid93 ', reidname,'l',li,lf,vnm  
               endif
              endif
             else if (li.eq.nint(jtot) .and. lf.eq.nint(jtot)) then !S=0 and S=1
             select case(li)
             case(0)
             reidname(1:2)='1S'
             case(1)
             reidname(1:2)='1P'
             case(2)
             reidname(1:2)='1D'
             case(3)
             reidname(1:2)='1F'
             case(4)
             reidname(1:2)='1G'
             case(5)
             reidname(1:2)='1H'
             case(6)
             reidname(1:2)='1I'
             case(7)
             reidname(1:2)='1K'
             case(8)
             reidname(1:2)='1L'
             end select
!             if (ir.eq.1)  write(0,*) 'Reid93 J', reidname(1:3)  
             call rreid93(r,reidname,'PN',reidpot)
             !S=0
             vnm=sqrt((2d0*ji+1d0)*(2d0*jf+1d0))*                       &
     &       (-1d0)**(jf-ji)/2d0/(2d0*li+1d0)*reidpot(1,1) !Factors are Racah with S=0
             !S=1
             select case(li)
             case(0)
             reidname(1:2)='3S'
             cycle
             case(1)
             reidname(1:2)='3P'
             case(2)
             reidname(1:2)='3D'
             case(3)
             reidname(1:2)='3F'
             case(4)
             reidname(1:2)='3G'
             case(5)
             reidname(1:2)='3H'
             case(6)
             reidname(1:2)='3I'
             case(7)
             reidname(1:2)='3K'
             case(8)
             reidname(1:2)='3L'
             end select
!            if (ir.eq.1)  write(0,*) 'Reid93 J', reidname(1:3) 
             call rreid93(r,reidname,'PN',reidpot) 
              if (ji.gt.(jtot+0d0)) then !First racah
              auxreid=(-1d0)**(li+ji+0.5d0)/2d0*sqrt((li+1.5d0-ji)*      &
     &        (li+ji-0.5d0)/ji) !Racah from table 9.1 Varshalovich
              else
              auxreid=(-1d0)**(li+ji+1.5d0)/2d0*sqrt((-li+1.5d0+ji)*      &
     &        (li+ji+2.5d0)/(ji+1d0))
              endif
              if (jf.gt.(jtot+0d0)) then !Second racah
              auxreid=auxreid*(-1d0)**(lf+jf+0.5d0)                         &
     &        /2d0*sqrt((lf+1.5d0-jf)*(lf+jf-0.5d0)/jf) !Racah from table 9.1 Varshalovich
              else
              auxreid=auxreid*(-1d0)**(lf+jf+1.5d0)                            &
     & /2d0*sqrt((-lf+1.5d0+jf)*(lf+jf+2.5d0)/(jf+1d0))
              endif
              vnm=vnm+auxreid*reidpot(1,1)
             endif
             ccmat(nchani,nchanf,ir)=cmplx(vnm,zero)*conv
          enddo   
          endif

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
            vcp=vtran(ici,icf,il,ir,1)
           else
            vcp=vtran(ici,icf,il,ir,2)
           endif
           
! TEST Jul/22           
           if ((nchani.eq.nchanf).and.(ici.eq.icf).and.(il.eq.0)) then
              vcp=vcp*vscale
           endif    
           
           faux(ir)=faux(ir)+((-1d0)**(ji+jcf+jtot)*
     &          sixj(jf,ji,lambdar,jci,jcf,jtot)*
     &          sqrt((2*jf+1.)*(2*ji+1.))*
     &			sixj(jf,ji,lambdar,lir,lfr,sn)*
     &		    (-1d0)**(lambdar+lfr+ji+sn)*
     &			sqrt((2.*lf+1)*(2.*il+1.)/4/pi)*
     &			cleb(lambdar,zero,lfr,zero,lir,zero))*vcp*cl            

!            if (partot.eq.1) then            
!           faux(ir)=faux(ir)+((-1d0)**(ji+jcf+jtot)*sixj(jf,ji,lambdar
!     &,jci,jcf,jtot)*sqrt((2*jf+1.)*(2*ji+1.))*sixj(jf,ji,lambdar,
!     &lir,lfr,sn)*(-1d0)**(lambdar+lfr+ji+sn)*sqrt((2.*lf+1)*(2.*il+1.)
!     &/4/pi)*cleb(lambdar,0d0,lfr,0d0,lir,0d0)
!     &vtran(ici,icf,il,ir,1))*cl
!            else
!           faux(ir)=faux(ir)+((-1d0)**(ji+jcf+jtot)*sixj(jf,ji,lambdar
!     &,jci,jcf,jtot)*sqrt((2*jf+1.)*(2*ji+1.))*sixj(jf,ji,lambdar,
!     &lir,lfr,sn)*(-1d0)**(lambdar+lfr+ji+sn)*sqrt((2.*lf+1)*(2.*il+1.)
!     &/4/pi)*cleb(lambdar,0d0,lfr,0d0,lir,0d0)*
!     &vtran(ici,icf,il,ir,2))*cl
!            endif
            enddo !il
        if (abs(faux(ir)).lt.small) faux(ir)=small ! avoid underflow
          enddo ! ir
        endif !cptype=5

!!!! For bincc !!!!!!!!!!!!!!!!!!!!!!!!!
!        ccmat(nchani,nchanf,1:nr)=cmplx(-faux(1:nr)*conv    
!     X      * (-1)**NINT((JCI-JCF - ABS(JCI-JCF))/2.),0d0)

! For "scatcc"    
! v2.1 I have changed the sign of ccmat so they have their physical sign
      if (.not.reid93) then
      do ir=1,nr
       if (abs(faux(ir)).lt.small) faux(ir)=0.0
       ccmat(nchani,nchanf,ir)=cmplx(faux(ir),zero)*conv
       enddo
      endif
     

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
      use wfs     , only: nr,dr,rvec,energ,wfc
      use scattering, only: ifcont,method,nbas,ns
      implicit none
c     ........................................................
      logical :: energy,writewf,writesol
      character*5::jpi
      integer:: il, pcon ! deprecated, backward compatibility only
      integer inc,ilout,ili,jset,partot,nchan,nk,n,ir,ik, chan
      real*8 :: emin,emax,eout,jtot,ecv
      real*8 :: kmin,dk,de,kmax,kcv
      real*8 :: r0,conv,rm
      real*8,allocatable :: psh_el(:)
!      real*8 :: deladd,deltai,deltap

! ... Overlaps
	  integer :: ic,ich,li,lf,nex,nho
      real*8, parameter :: eps=1e-3
      real*8 :: sumr,jac,excore,jci,jcf,jf
      real*8, allocatable:: erel(:)
      complex*16,allocatable,target:: gsolap(:,:)
      complex*16 :: gaux(nr)
      complex*16 :: resc,caux,raux
      

      complex*16,allocatable:: smate(:,:),wfcont(:,:,:)
      namelist/scatwf/ emin, emax,ifcont,nk,inc,writewf,ilout,eout,jset,
     &                 energy, writesol,method,nbas,ns,
     &                 il,pcon ! not used, backward compatibility only (AMM)

      inc=0; il=0
      jset=0; pcon=0     
!      write(0,*)'method=',method
      method=4;
      nbas  =50
      ns    =1; 
  
      writewf  =.false.
      writesol =.false.
      energy   =.false.
      

      
      write(*,'(//,5x,"******  PROJECTILE SCATTERING STATES ******")')
      read(kin,nml=scatwf)
      
      if (method.eq.0) method=4
      
      if (method.eq.5) then
       write(*,*)' Using R-matrix with nbas=',nbas,' ns=',ns
      endif
      
      if ((.not.ifcont).or.(nk.eq.0)) return
      if ((jset.gt.jpsets).or.(jset.eq.0)) then
       write(*,*)'jset',jset,' not valid!'
       return
      endif
      if (il.gt.0) then
       write(*,*)' ** WARNING ** IL deprecated; use INC instead!'
      endif


      partot=jpiset(jset)%partot
      jtot  =jpiset(jset)%jtot
      nchan =jpiset(jset)%nchan
      ql (1:nchan)   =jpiset(jset)%lsp(1:nchan)
      r0=rvec(1)
      write(*,400) jpi(jtot,partot), nchan,inc,energy
      write(*,*) 

      rm=av*ac/(av+ac)
      conv=(2*amu/hc**2)*rm

      kmin=sqrt(conv*abs(emin))
      kmax=sqrt(conv*abs(emax))

      if (energy) then
        de=(emax-emin)/(nk-1)
      else
        dk=(kmax-kmin)/(nk-1)
      endif

! Store continuum energies in erel() grid
      allocate(erel(nk))
      do ik=1,nk
      if (energy) then
        ecv=emin+(ik-1)*de
      else
        kcv=kmin+ (ik-1)*dk  !new kcv for schcc
        ecv=(hc*kcv)**2/2/mu12
      endif
      erel(ik)= ecv
      enddo ! ik

      if(allocated(smate)) deallocate(smate)
      allocate(smate(nk,maxchan),psh_el(nk))
      if (allocated(wfcont)) deallocate(wfcont)
      allocate(wfcont(nk,nchan,nr))

      chan=inc
      do inc=1,nchan
      if (chan.ge.0.and.chan.ne.inc) cycle

!      deladd=0.
!      write(0,*)'continuum_range: calling wfrange';stop
      call wfrange(jset,nchan,inc,emin,emax,nk,energy,wfcont,
     &     smate,psh_el,writewf)
     

c ... Write WFS
      if (writewf) then
      write(50,400) jpi(jtot,partot), nchan,inc,energy !, ql(1:nchan)
400   format(5x,"# Continuum WFS for Jp=",a5,"; Channels=",i2,
     & " ;  Inc. chan.=",i2, ' Uniform Ecm=',l1)
      write(50,'(1x,i4,2x,2f8.4)') nr,dr,r0
      if (energy) then
       write(50,'(1x,i4,2x,2f8.4)') nk,de,emin
      else
       write(50,'(1x,i4,2x,2f8.4)') nk,dk,kmin
      endif

      do ik=1,nk   
      ecv=erel(ik)
      write(50,*)'# Continuum wf with Ecm=',ecv
      do ir=1,nr
      write(50,'(1f8.3,2x,50f12.8)')rvec(ir),
     &       (wfcont(ik,n,ir),n=1,nchan)
      enddo ! ir
      write(50,*)'& '
      enddo ! ik
      endif !writewf
      
      
! ... Compute ovelaps between scattering states and PS's (if defined)
      nho= jpiset(jset)%nho
      if (nho.eq.0) return ! no PS's in this jset
      nex =jpiset(jset)%nex             ! number of states in this jset
      if(.not. allocated(gsolap)) allocate(gsolap(nex,nk))
      gsolap(:,:)=0
      
!      if (emax<0)     emax=maxval(energ(:,:))
!      if (emin.le.0.) emin=0.01
!      if (emax>Ethr)  emax=Ethr-0.001
!      dec=(emax-emin)/dble(nk-1)
!      do jset=1,jpsets
!      if (.not.jsets(jset)) cycle
!      if (jpiset(jset)%nho.eq.0) cycle
!      nchan=jpiset(jset)%nchan
!      do inc=1,nchan 
      if (inc.gt.nchan) stop'inc too large in SCATWF!'
      jci    =jpiset(jset)%jc(inc)
      li     =jpiset(jset)%lsp(inc)
      excore =jpiset(jset)%exc(inc)     ! core energy
      ic     =jpiset(jset)%cindex(inc)  
!      if (ic.ne.icore) then
!         if (verb.ge.2)
!     &   write(*,'(4x,"-> skipping core state",i3)')ic 
!      cycle
!      endif

!      if (allocated(wfcont)) deallocate(wfcont)
!      if (ncont.gt.0) then
!        allocate(wfcont(ncont,nchan,nr)); wfcont=0.
!        call wfrange(jset,nchan,inc,emin,emax,ncont,energy,wfcont,
!     & smate,delta)
!      else
!       write(*,*)'*ERROR* ncont,nchan=',ncont,nchan
!      endif

      do n=1,nex
      if (energ(jset,n).lt.0)    cycle ! bound state
      do ik=1,nk
      ecv=erel(ik)
      if (ecv.lt.eps) ecv=eps
      kcv=sqrt(2.d0*mu12*ecv)/hc    
          
      sumr=0
      caux=0.
      do ich=1,nchan !  sum over outgoing channels
      jcf =jpiset(jset)%jc(ich)
      lf  =jpiset(jset)%lsp(ich)
      jf  =jpiset(jset)%jc(ich)
      do ir=1,nr   
      gaux(ir)=wfc(jset,n,ich,ir)*wfcont(ik,ich,ir)*rvec(ir)*
     &         (-1)**(li+jci+lf+jcf) ! AMM: I think this phase is always 1! 
      enddo !ir

      call simc(gaux,resc,1,nr,dr,nr)
      caux=caux+resc
      sumr = sumr+abs(caux)**2
      enddo ! ich 
      gsolap(n,ik)=caux
      enddo ! ik (c-v relative energy)
      enddo ! n  (PS within this j/pi set)
!      enddo ! inc  (incoming channel)
!      enddo ! jset (j/pi set)
     
c *** -------------- PRINT OVERLAPS  ----------------------------
      if (writesol) then
      written(52)=.true.
      written(53)=.true.
!      do jset=1,jpsets
!      if (.not.jsets(iset)) cycle
!      nex  =jpiset(iset)%nex
!      nchan=jpiset(iset)%nchan
      if (jpiset(jset)%nho.eq.0) return
      
      write(52,'(1x,"# Overlaps for set with J/pi=",a5)')
     &     jpi(jtot,partot)
      write(52,'(1x,"# THO basis with nho=",i3," functions")') nho
      write(52,'(1f8.3,2x,1f8.3,2x,i5,3x,L1,3x,i3)') 
     & emin,emax,nk,energy,nex
      
      do n=1,nex
      write(52,'(i3,1f10.4)') n , energ(jset,n)
      raux=0
      write(53,'(a,i3,a,i3,a,1f10.4)')'# Set:',jset,
     & ' n=',n,' Ex=',energ(jset,n)
      do ik=1,nk
       ecv=erel(ik)
       if (energ(jset,n).lt.0)    cycle ! bound state
       if (ecv.lt.eps) ecv=eps
       kcv=sqrt(2.d0*mu12*ecv)/hc
       jac=mu12/(hc**2)/kcv
       write(53,111)ecv, jac*abs(gsolap(n,ik))**2
       write(52,'(1f8.3,3x,2g14.6)')ecv,gsolap(n,ik)
!     &  (jac*abs(gsolap(n,ik))**2,inc=1,nchan)
111    format(2x,1f12.6,2x,10g14.6)
!       do inc=1,nchan
!       raux=raux+jac*abs(gsolap(n,inc,ik))**2*dec
!       enddo ! inc
      enddo !ik
!      if (verb.ge.3) then
!      write(97,'(3x,a,i3,a,i3,a,1f8.5)')'# -> Set:',jset, '  PS=',n,
!     & ' Norm solap=',raux*2/pi
!      endif
      write(53,*)'&'
      enddo ! n
!      enddo !iset
      endif ! verb
      
      enddo !inc
  
c-------------------------------------------------------------------------------
      
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
      integer n,ir,ich,ichp,nskip,klog !,method
      real*8:: r,r0,dk,kcont,econt,k,ecm, t1,t2
      real*8:: phase(nchan),aux,phadd(nchan)
      real*8,allocatable:: ph(:,:)
c  bincc variables
      logical tres,tdel,tkmat,tknrm
      integer:: nf   ! number of formfactors?
      integer:: isc,iil,ik,jset
      integer:: maxn,nnk(nchan)
      real*8:: anc,k2(nchan),conv,rm,bphase(nk)
      complex*16:: y(nr,nchan),etap(nchan),wf(nchan,nr),smat(nchan)
      real*8:: fmscal
c  end bincc 
      namelist/scatwf/ emin, emax,ifcont,nk,il,ilout,writewf,jset
      il=0
      ili=1
      klog=99
      method=4
      nbas= 50
      ns=   1
      read(kin,nml=scatwf)
      
      if (method.eq.0) method=4
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
      
      write(*,*)'R-matrix solution with nbas,ns=',nbas,ns

      write(*,*)
      write(*,100)nk,emin,emax
100   format(3x,"Generating", i4, " continuum wf(s) for"
     &  ," [Emin=",1f6.3," Emax=",1f7.2," MeV]")

      if (method.eq.5) write(*,*) "(using R-matrix routine)" 
      call factorialgen(100)

      call coefmat(jset,nchan)

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
!       method=4
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




