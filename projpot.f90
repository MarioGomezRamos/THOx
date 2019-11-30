      subroutine readpot(kin) 
c *** Read potential parameters and store the potential in the radial grid
c  The full CENTRAL potential will be assumed of the form
c
c    V(r)=Vcou(r)+ Vc(r)+ Vls*l.sv + Vss*sc.sv + Vll*l.l
c
c    Vcou(r) = coulomb central
c    Vcl(r)   = l-dependent nuclear central potential
c    Vls(r)  = spin-orbit potential for valence (s from valence)
c    Vlsc(r) = spin-orbit potential for core (sc from core)
c    Vss(r)  = spin-spin potential (sv.sc)
c    Vll(r)  = l.l potential          
c
c    The non-central potential is 
c     Sum_{lambda} Vlcoup(l,lambda,r) Y_{lm}(\hat r) 
c 
       use potentials
       use sistema
       use constants, only: pi
       use channels, only: qjc,cindex
       use globals, only: written,verb
       use wfs, only:rvec,nr,rmin,rmax,dr
       use parameters, only: maxl,maxcore!,maxchan
       implicit none
       logical:: lsderiv,uu
       integer kin,ir,ith,ncomp,mmultipole,l,p
       parameter(mmultipole=4) !!better in a module???
       real*8:: alpha,fival
       parameter(alpha=0.)
       real*8 a13,r,rc,rabs,rcp,rcpabs,raux,aaux
       real*8:: vcoul,ws,pt,gausspot,gausspotr2
       real*8:: vlsaux,wsso,gausder,dws,vllaux
       real*8:: vaux,vaux2,vssaux,sqwell
       real*8:: vcpaux,vdef(0:mmultipole)
       real*8:: YLMC,pl(4,4),test,theta,w
       integer np,nv,iv
       real*8, allocatable:: rvecp(:),vvec(:,:)
       integer jci,jcf,ici,icf
       real*8:: pfact(10),pgfact(2)
       real*8 :: vcp(0:maxl,nr), vc(0:maxl,nr)
       real*8 :: vcaux(nr),v0,v02
       character*15::potname
       real*8:: ap,at             ! masses for radius conversion
       real*8:: R0,A0        ! central WS
       real*8:: Vso(0:10),rso,aso       ! spin-orbit for valence MGR
       real*8:: Vsoc(0:10),rsoc,asoc       ! spin-orbit for core !MGR
       real*8:: Vss0(0:10),rss,ass       ! spin-spin
       real*8:: Vll0(0:10),rll,all       ! l.l !MGR
       real*8:: rc0               ! charge radius
       real*8:: rcp0
       real*8:: p0,acp !coupling
       real*8:: vl0(0:10),vcp0(0:10),normr(0:10)
       real*8:: th(40),wth(40)       
       real*8::fvol
       real*8:: vg1,rg1


       
       namelist /potential/ ptype,ap,at,
     & vl0,r0,a0,rc0,          !central
     & Vso,rso,aso,lsderiv,    !spin-orbit: l*sv
     & Vsoc,rsoc,asoc,         !spin-orbit: l*jc
     & Vss0,rss,ass,           !spin-spin : sv.sc
     & Vll0, rll,all,          !l.l
     & cptype,lpot,beta,delta,Vcp0,rcp0,acp,
     & np,nv,pfact,pgfact,
     & kband, lambda,pcmodel,normr

       pcmodel=0
       delta=0d0
       beta=0d0
       cptype=0
       ncomp=0 !number of potential components
       laminc(:)=.false.

		if (.not.allocated(vcl)) then
	  allocate(
     &    vcl(0:maxl,1:nr),     ! central, l-dependent
     &    vls(0:maxl,1:nr),            ! spin-orbit for valence MGR
     &    vlsc(0:maxl,1:nr),           ! spin-orbit for core MGR
     &    vss(0:maxl,1:nr),     ! spin-spin
     &    vll(0:maxl,1:nr),            ! l.l MGR
     &    vlcoup(0:maxlamb,0:maxl,1:nr),  ! coupling
     &    vcou(1:nr))           ! Coulomb monopole

          vcl(:,:)=0d0; vls(:,:)=0d0; vss(:,:)=0d0; vlsc(:,:)=0d0!MGR
          vcou(:)=0d0; vcaux(:)=0d0; vlcoup(:,:,:)=0d0;vll(:,:)=0d0!MGR
        endif

       vcp(:,:)=0d0

100    V0=0; vl0(:)=0d0; r0=0d0; a0=0d0; 
       Vso(:)=0d0; rso=0d0; aso=0d0; rc0=0d0!MGR
       Vsoc(:)=0; rsoc=0; rsoc=0!MGR
       Vcp0(:)=0d0; rcp0=0d0; acp=0d0!MGR
       Vss0(:)=0d0; rss=0d0; ass=0d0!MGR
       Vll0(:)=0d0; rll=0d0; all=0d0!MGR
       vcp(:,:)=0d0
       ap=0
       at=1
       if (cptype.ne.5) cptype=0
       normr(:)=1.
       delta=0d0
       beta=0d0

       lambda=2 
       kband=0
!       maxlamb=4
       written(40)=.false.
       pfact(:)=1d0
       pgfact(:)=1d0
       lsderiv=.true.

       read(kin,nml=potential)
       A13=Ap**.33333333 + At**.33333333
       rabs=r0*a13
       if (delta.eq.0d0) delta=beta*rabs

       if (ptype.eq.0) then
	if (ncomp.eq.0) then
	   write(*,*)'No potentials were found!'
           stop
        else
	  goto 200 !return !last component just read
        endif
       endif
       ncomp=ncomp+1
       write(*,'(3x,"------------- Potential component:",i2,
     &           "-------------" )') ncomp 
!       write(*,*)'number of potentials=',ncomp
       write(*,'(/," o Undeformed valence-core potential: ")')
       write(*,118) Rc0*a13
 118   format(3x,'COULOMB     : uniform charged sphere with Rc=',
     &        1f6.3,' fm')
!MGR--------------------------------------------------------------------
!In case the l dependence is ignored, take the given potential for all l
      if ((abs(Vso(0)).gt.0d0).and.(sum(abs(Vso(1:10))).lt.1e-12)) then
        do l=1,10
        Vso(l)=Vso(0)
        enddo
      endif
      
      if ((abs(Vsoc(0)).gt.0d0).and.(sum(abs(Vsoc(1:10))).lt.1e-12))then
        do l=1,10
        Vsoc(l)=Vsoc(0)
        enddo
      endif
      
      if((abs(Vll0(0)).gt.0d0).and.(sum(abs(Vll0(1:10))).lt.1e-12)) then
        do l=1,10
        Vll0(l)=Vll0(0)
        enddo
      endif
!MGR--------------------------------------------------------------------
       if (abs(vl0(0))>0d0) v0=vl0(0)
       write(*,'(3x,"NUCLEAR part:",$)')
       select case(ptype) 
c ------------------------------------------------------
        case(0)  ! Last potential already read
c ------------------------------------------------------
	   return
c ------------------------------------------------------
        case(1)  ! WS type
c ------------------------------------------------------
           write(*,*)'type=1 => Woods-Saxon'
           write(potname,*) "Woods-Saxon"
c ------------------------------------------------------          
           write(*,119) potname, Vl0(0:3),r0*a13,a0
           if (abs(maxval(Vso(0:3))).gt.0)write(*,122)Vso(0:3),
     &      rso*A13,aso !MGR

 119       format(7x,'o ',a15,"V0(l=0-3)=(",4f8.3,') MeV, R0=',
     &           f6.3,' fm, ar=',f6.3,' fm ')

 120       format(7x,'o Central   : V0(l=0-3)=(',4f8.3,') MeV, R0=',
     &           f6.3,' fm, ar=',f6.3,' fm ')
!     &           ' Rc=',f6.3,' fm ]')
 122       format(7x,'o Spin-orbit: Vso=(0:3)',4f8.3,' MeV, Rso=',f6.3, !MGR
     &           ' fm, aso=',f6.3,' fm ',/)
 124       format(7x,'o Spin-spin : Vss(l=0-3)=(',4f8.3,') MeV, Rss=',
     &          f6.3, ' fm, ass=',f6.3,' fm ',/)

c ------------------------------------------------------
        case(2) !  Posch-Teller
c ------------------------------------------------------
           write(*,*)'type=2 => Posch-Teller'
           write(*,130) Vl0(0:3),a0
 130       format(3x,'[V0(l=0-3)=(',4f8.3,') MeV, ar=',f6.3,
     &           ' fm ]',/)

c ------------------------------------------------------
        case(3) ! Gaussian 
c ------------------------------------------------------
             write(potname,*) "Gaussian"
             write(*,*)'type=3 => Gaussian'
             write(*,119) potname,Vl0(0:3),r0*a13,a0
         if (abs(maxval(Vso(0:3))).gt.0) write(*,122) Vso(0:3),
     &    rso*A13, aso!MGR
             if (maxval(abs(Vss0)).gt.0) 
     &          write(*,124) Vss0(0:3),rss*A13, ass

 140       format(3x,'[Vso(l=0-3)=(',4f8.3,') MeV  Rso=',f6.3,
     &           'aso=',f6.3,' fm  Rc0=',f6.3,' fm ]',/)

c ------------------------------------------------------
        case(30) ! Gaussian  x r2
c ------------------------------------------------------
             write(potname,*) "Gaussian"
             write(*,*)'type=30 => Gaussian x r2' 
             write(*,119) potname,Vl0(0:3),r0*a13,a0
         if (abs(maxval(Vso(0:3))).gt.0) write(*,122) Vso(0:3),
     &    rso*A13, aso!MGR
             if (maxval(abs(Vss0)).gt.0) 
     &          write(*,124) Vss0(0:3),rss*A13, ass

 142       format(3x,'[Vso(l=0-3)=(',4f8.3,') MeV  Rso=',f6.3,
     &           'aso=',f6.3,' fm  Rc0=',f6.3,' fm ]',/)



c ------------------------------------------------------
        case(4) !  Externally read
c ------------------------------------------------------
           write(*,*)'- type=4 => Read from fort.2x files'
           write(*,153)np, maxl
153        format(3x,i5,' points for maxl=',i3)

c ------------------------------------------------------
        case(5) !  Externally read
c ------------------------------------------------------
           write(*,*)'- type=5 => <I|V|I> read from fort.30 file'
           write(*,154)np,nv
154        format(3x,i5,' points & ',i3,' potentials')
           write(*,*)'- +WS s.o. with parametres:'
           write(*,122) Vso,rso*A13, aso


c ------------------------------------------------------
        case(10) !  Externally read
c -----------------------------------------------------

c ------------------------------------------------------
        case default
c ------------------------------------------------------
          write(*,*)'Potential type=',ptype,' not used'
          stop
        end select
 150    format(6f10.5)


c coupling component, if any
       select case(pcmodel)
       case(0)
        if (abs(beta).gt.1e-4.or.abs(delta).gt.1e-4) then
        write(*,'(" o Particle-core coupling:",$)')
        write(*,*)'  [rotor model with lambda=',lambda,']'
        laminc(lambda)=.true.
        endif
       case(1)
        write(*,'(" o Particle-core coupling:",$)')
        write(*,*)'  [vibrator model with lambda=',lambda,']'
        laminc(lambda)=.true.
        if (cptype.ne.1) then
         write(*,*)'     => Forcing cptype=1 (derivative)'
         cptype=1
        endif
       case default
        write(*,*)'Invalid particle-core model',pcmodel
        stop
       end select

       select case(cptype)
c ------------------------------------------------------
        case(1)  ! dV/dr type
c ------------------------------------------------------
	if (ptype.eq.1) then
          write(*,*)'- Coupling Pot. type=1 => deriv. W-S'
          if (Vcp0(0).eq.0) then
             write(*,141) delta
	   else
!	      if(rcp0.eq.0) rcp0=rc0
!	      if(acp.eq.0)  acp=a0
	      write(*,143)Vcp0(0),acp,rcp0
	      write(*,*)'- We are NOT using DELTA!!'
	   endif
 143      format(3x,'[Vcp=',1f8.3,' MeV acp=',f6.3,'  rcp0=',f6.3,
     &           ' fm ]',/)
 141      format(3x,'[delta=',1f8.3,' fm ]',/)
	else
	   write(*,*)'- Coupling Pot. type=1 => deriv. PT/Gaussian'
	   write(*,*)'- WARNING!!!! => not implemented yet'
        endif
c ------------------------------------------------------
        case(2,3)  ! Projection on multipoles
c ------------------------------------------------------
	  write(*,*)
     &       ' => Deform central potential & project on P_{lambda}'
         if (cptype.eq.3) write(*,*)
     &       '(monopole potential will be also recalculated!)'
           write(*,141) delta
c ------------------------------------------------------
        case(4)  ! V_2 type
c ------------------------------------------------------
	   write(*,*)'- Coupling Pot. type=2 => deformed fermi function'
	   write(*,*)'- WARNING!!!! => work in progress'                      !remember to erase!
	   write(*,144)beta
 144       format(3x,'[beta=',1f8.3,']',/)
	end select
	    

! ---------------------- END OF INFORMATIVE SECTION --------------------


c .. Monopole Coulomb ... 
       Rc=Rc0*A13
       if (ncomp.eq.1) then  !! coulomb potential should not be added to previous one
         do ir=1,nr
          r=rvec(ir)
	  vcou(ir)=VCOUL(r,zv,zc,Rc) 
         enddo
       endif

c .. External potential from fort.20 
        if (ptype.eq.4) then
        if (allocated(rvecp)) deallocate(rvecp) 
        if (allocated(vvec))  deallocate(vvec)
        allocate(rvecp(np),vvec(np,0:maxl))
        rvecp(:)=0d0
        vvec(:,:)=0d0
        open(20)
        do ir=1,np
!        read(20,'(12f10.4)')rvecp(ir),(vvec(ir,l),l=0,maxl)    !what happens if there is only a pot for all l's    
        read(20,*)rvecp(ir),(vvec(ir,l),l=0,maxl)    !what happens if there is only a pot for all l's    
        if (ir.lt.20) print*,rvecp(ir),vvec(ir,0)   
        enddo
        do l=0,maxl
        do ir=1,nr
        vaux=fival(rvec(ir),rvecp(:),vvec(:,l),np,alpha)  
!        write(*,*)rvec(ir),normr(l)*vaux
        vcl(l,ir)=vcl(l,ir) + normr(l)*vaux
        if (ir.eq.1)print*, 'norm,vaux=',normr(0),vaux
        enddo
        enddo
c .. External spin-orbit potential from fort.21
        uu = .false.
        inquire(file="fort.21",exist=uu)
        if (uu) then
        open(21)       
        rvecp(:)=0d0
        vvec(:,:)=0d0
        do ir=1,np
        read(21,*)rvecp(ir),vvec(ir,1) 
        enddo
        do ir=1,nr
        vls(0,ir)=fival(rvec(ir),rvecp(:),vvec(:,1),np,alpha) !MGR l-dependence not included here  
         do l=1,maxl!MGR
           vls(l,ir)=vls(0,ir)
         enddo!MGR
        enddo
        else
           write(*,*)'No spin orbit component found'
        endif
!---------------------------------------read from fort.22 coupling potential
        uu = .false.
        inquire(file="fort.22",exist=uu)
        if (uu) then
        open(22)
        rvecp(:)=0d0
        vvec(:,:)=0d0
        do ir=1,np
        read(22,'(12f10.4)')rvecp(ir),(vvec(ir,l),l=0,maxl)    !what happens if there is only a pot for all l's
        enddo
        do l=0,maxl
        do ir=1,nr
        vlcoup(lambda,l,ir)=fival(rvec(ir),rvecp(:),vvec(:,l),np,alpha)   
        enddo
        enddo
        else ! if file exists
           write(*,*)'No coupling potential found'
        endif
        endif
!-----------------------------------------



!---------------------------------------read from fort.30        Potentials from densities, <I|V|I'>
        if (ptype.eq.5) then
        cptype=5                     !as ptype is reset, cptype is the mark of being in ptype=5
!        write(*,*)qjc
!        write(*,*)maxval(nint(qjc)),maxval(nint(qjc)),
!     &maxval(cindex),maxval(cindex)
!        allocate(vtran(maxval(cindex),maxval(cindex),0:maxlamb,nr))
        allocate(vtran(maxcore,maxcore,0:maxlamb,nr,2))
        allocate(rvecp(np),vvec(np,1))
        rvecp(:)=0d0
        vvec(:,1)=0d0
        vtran(:,:,:,:,:)=0d0
        do iv=1,nv
        read(30,*)jci,jcf,ici,icf,l

        do ir=1,np
        read(30,*)rvecp(ir),vvec(ir,1)
        enddo
        do ir=1,nr
!        if (rvecp(np).gt.rvec(ir)) then
        DO P=1,2
        vtran(ici,icf,l,ir,p)=vtran(ici,icf,l,ir,p)+
     &fival(rvec(ir),rvecp(:),vvec(:,1),np,alpha)*pfact(iv)*pgfact(P)
        ENDDO

        write(31,*)rvec(ir),vtran(ici,icf,l,ir,1),vtran(ici,icf,l,ir,2)

        enddo
        vtran(icf,ici,l,:,:)=vtran(ici,icf,l,:,:)             !symmetric potentials

        write(31,*)'& '
        read(30,*)         !separation between different potentials, all in fort.30

        enddo
        
        write(*,*)'potentials read from fort.30: ok! cptype=',cptype
!        do ir=1,nr
!        write(400,*)rvec(ir),vtran(1,2,2,ir),vtran(2,1,2,ir)
!        enddo
 
        endif ! ptype=5

c ... Central and deformed formfactor ...
        rcp=rcp0*a13 
        do l=0,maxl
        v0=vl0(l)
        v02=vcp0(l)
        if((abs(v0).lt.1e-6).and.(abs(v02).lt.1e-6)) cycle
        do ir=1,nr
           r=rvec(ir) 
           vaux=0d0; vaux2=0d0
           select case(ptype)
c .. WS
           case(1)  ! WS type
	     if (abs(v0).gt.1e-6) vaux=ws(r,v0,rabs,a0)
             if (abs(v02).gt.1e-6) vaux2=ws(r,v02,rcp,acp)
c .. PT
           case(2) !  Posch-Teller
	     vaux=pt(r,v0,a0)
             vaux2=pt(r,vcp0(l),acp)
c .. Gaussian
           case(3) ! Gaussian 
 	      if (abs(v0).gt.1e-6) vaux=gausspot(r,v0,r0*a13,a0)
              if (abs(v02).gt.1e-6)vaux2=gausspot(r,v02,rcp*a13,acp)
c .. Gaussian x r2
           case(30) ! Gaussian x r2 
 	      if (abs(v0).gt.1e-6) vaux=gausspotr2(r,v0,r0*a13,a0)
              if (abs(v02).gt.1e-6)vaux2=gausspotr2(r,v02,rcp*a13,acp)
c .. Square well
           case(10) ! Square well 
 	      if (abs(v0).gt.1e-6) vaux=sqwell(r,v0,r0*a13)
              if (abs(v02).gt.1e-6)vaux2=sqwell(r,v02,rcp*a13)


	   end select
	    vc(l,ir)= vaux
!            if (ir.eq.1) 
!     &    print*,'adding central potential ptype,v=',ptype,vaux
           if (abs(v02).gt.1e-6) then
              vcp(l,ir)=vcp(l,ir)+ vaux2
           else
              vcp(l,ir)=vcp(l,ir)+ vaux
           endif
         enddo !ir
       enddo !l

c ... Spin-orbit for valence ....
       do l=0,maxl !MGR
       if ((abs(vso(l))).ge.1e-5) then!MGR       
       do ir=1,nr
         r=rvec(ir) 
         if (r.lt.1e-6) r=0.5*dr
         if (lsderiv) then  ! derivative form-factor
           select case(ptype)
           case(1,5)  ! WS type
                vlsaux=wsso(r,vso(l),rso*a13,aso)!MGR		
!            case(2) !  Posch-Teller
!		vlsaux=pt(r,v0,a0)
           case(3) ! Gaussian 
 	        vlsaux=gausder(r,vso(l),rso*a13,aso)!MGR
	   end select
          else ! non-derivative
           select case(ptype)
           case(1)  ! WS type
                vlsaux=ws(r,vso(l),rso*a13,aso)	!MGR
           case(2) !  Posch-Teller
		vlsaux=pt(r,vso(l),aso)
           case(3) ! Gaussian 
 	        vlsaux=gausspot(r,vso(l),rso*a13,aso)!MGR
	   end select            
          endif
	    vls(l,ir)=vls(l,ir)+vlsaux !MGR
        enddo !ir
        endif
        enddo !l MGR
        
c ... Spin-orbit for core ...

       do l=0,maxl!MGR
       if (abs(vsoc(l))>1e-5) then !MGR
       do ir=1,nr
           r=rvec(ir) 
           if (r.lt.1e-6) r=0.5*dr
           if (lsderiv) then ! derivative formfactor
            select case(ptype)
            case(1)  ! WS type
              vlsaux=wsso(r,vsoc(l),rsoc*a13,asoc)		!MGR
!            case(2) !  Posch-Teller
!		vaux=pt(r,v0,a0)
            case(3) ! Gaussian 
 	        vlsaux=gausder(r,vsoc(l),rsoc*a13,asoc)!MGR
	    end select
           else
            select case(ptype)
            case(1)  ! WS type
              vlsaux=ws(r,vsoc(l),rsoc*a13,asoc)		!MGR
            case(2) !  Posch-Teller
		vaux=pt(r,v0,asoc)
            case(3) ! Gaussian 
 	        vlsaux=gausspot(r,vsoc(l),rsoc*a13,asoc)!MGR
	    end select
           endif !derivative formfactor

	    vlsc(l,ir)=vlsc(l,ir)+vlsaux !MGR
        enddo !ir
        endif
        enddo !l MGR

! ...  l.l (only central!!!!)   ...
       do l=0,maxl!MGR
       if (abs(vll0(l))>1e-5) then !MGR
       do ir=1,nr
           r=rvec(ir) 
           if (r.lt.1e-6) r=0.5*dr
            select case(ptype)
            case(1)  ! WS type
              vllaux=ws(r,vll0(l),rll*a13,all)		!MGR
!            case(2) !  Posch-Teller
!		vaux=pt(r,v0,a0)
            case(3) ! Gaussian 
 	        vaux=gausspot(r,vll0(l),rll*a13,all)!MGR
	    end select
	    vll(l,ir)=vll(l,ir)+vllaux !MGR 
        enddo !ir
        endif
        enddo !l MGR

! ...  spin.spin (sv.jc)  
       do l=0,maxl
       v0=vss0(l)
       if (abs(v0)>1e-5) then
       do ir=1,nr
           r=rvec(ir) 
           if (r.lt.1e-6) r=0.5*dr
            select case(ptype)
            case(1)  ! WS type
              vssaux=ws(r,v0,rss*a13,ass)		
            case(2) !  Posch-Teller
              vssaux=pt(r,v0,ass)
            case(3) ! Gaussian 
 	      vssaux=gausspot(r,v0,rss*a13,ass)
	    end select
	    vss(l,ir)=vss(l,ir)+vssaux 
!	    if (abs(vss(l,ir))>0.1) write(0,*) vss(l,ir)  
        enddo !ir
        endif
        enddo ! l

      
    
c Coupling (deformed) potential
       do l=0,maxl
	if ((abs(vl0(l)).lt.1e-6).and.(abs(vcp0(l)).lt.1e-6)) cycle
        if (abs(vcp0(l)).lt.1e-5) then
            v0=vl0(l)
!            vc=>vcl(l,:)
            raux=rabs
            aaux=a0
        else
           v0=vcp0(l)
!           vc=>vcp(l)
           raux=rcp ! =rcp0*a13
           aaux=acp
        endif
        
        do ir=1,nr
        r=rvec(ir) 
        vcpaux=0.

        
        select case(cptype)
             case(0) ! no coupling
!                       write(*,*)'No coupling cptype=',cptype
             case(1) ! derivative formfactor
                vcpaux=-delta*dws(r,v0,raux,aaux)  

             case(2,3) ! Numerical projection on multipoles
!             if (ir.eq.1)  
!     &        write(99,'("l,ir,vcp(l,ir)=",2i3,5f12.6)') l,ir,vcp(l,1:5) 
!                 call potdef(vcl(l,:),lambda,r,vdef,delta)
                 call potdef(vcp(l,:),lambda,rvec(ir),vdef,delta)    ! when calling with r instead of rvec gives problem in mac gfortran
                 vcpaux=vdef(lambda)*sqrt(4d0*pi)
                 vcaux(ir)=vdef(0)
                 if ((cptype.eq.3).and.(v0.gt.1e-6)) 
     &           vc(l,ir)=vcaux(ir)
             case(-1) ! Idem 
                if (ir.eq.1) then
                if (beta.lt.1e-6) beta=delta/(r0*a13)
 
                th(:)=0.
                wth(:)=0.
!                call gauss(-1.0,1.,5,th,wth)   ! I don't know but it is not working
                call setmgl(40,th,wth)
                th(:)=(0.+pi)/2.0 + (pi-0.)/2.0*th(:)
                wth(:) = (pi-0.)/2.0*wth(:)
                call logfac3(100)
!                write(*,*)YLMC(2,0)    !with respect to a legendre function with values between -1 and 1
                fvol=0.
                test=0.
                write(*,*)'beta=',beta
                do ith=1,40
                theta=th(ith)
                call PLM(cos(th(ith)),2,0,2,pl)
                  w=wth(ith)
! Commented by AMM
!                 fvol=fvol+(1+beta*YLMC(2,0)*pl(3,1))**3*sin(theta)*w
                 fvol=fvol+(1+beta*YLMC(2,0)*pl(3,1))**3*sin(theta)*w
                 test=test+2*pi*(YLMC(2,0)*pl(3,1))**2*sin(theta)*w
                enddo
                fvol=fvol/2
                write(*,*)test,fvol,r0*a13,r0*a13/fvol**.333 !why test is not equal to 1?
                endif !ir=1
 

               case default
                  print*,'coupling potential',cptype,'not implemented!'
                  stop
                 !HEEEEEEEEEEEEERRRREEEEEEEEEEEEE!!!!!!
		end select
              vlcoup(lambda,l,ir)=vlcoup(lambda,l,ir)+vcpaux
!             if (ir.eq.5)  
!     &        write(99,'("l,ir,vcp(l,ir)=",2i3,5f12.6)') l,ir
!     &        ,vlcoup(lambda,l,1:5) 

              enddo !ir
!!!!  IMPORTANT!!!!
          if ((cptype.eq.3).and.(abs(vl0(l)).gt.1e-6)) then
c          write(*,*)'Monopole nuclear potential recalculated by quad!'
	      vc(l,:)=vcaux(:)
! JALay testing monopole recalculation
!          do ir=1,nr
!          write(700,*)rvec(ir),vcl(0,ir),vcaux(ir)
!          enddo
!---------------------------------------------------
          endif
c add monopole to previous potential
         if (abs(vl0(l)).gt.1e-6) then
           vcl(l,:)=vcl(l,:)+vc(l,:)
         endif
         enddo !maxl
       write(*,*)'--------------------------------------------------'
	goto 100 !read next potential component


200   continue

! TEMP
      do ir=1,nr
         write(40,*)rvec(ir),vcl(0,ir) 
      enddo


      if ((.not.written(40)).and.(verb.ge.2)) then
         written(40)=.true.
        close(20)
        close(21)
        close(22)
        close(23)
!        write(22)'#Coupling potential'
      do ir=1,nr
         write(40,*)rvec(ir),vcou(ir),vcl(0,ir) 
      if ((cptype.ne.5).and.(ptype.ne.4)) then
         written(20)=.true.

         write(20,'(12f10.4)')rvec(ir),(vcl(l,ir),l=0,maxl)
         write(21,'(12f10.4)')rvec(ir),(vls(l,ir),l=0,maxl)!MGR
         if (laminc(2)) then 
           written(22)=.true. 
           write(22,'(12f10.4)')rvec(ir),(vlcoup(2,l,ir),l=0,maxl)
         endif
         if (laminc(3)) then 
           written(23)=.true. 
           write(23,'(12f10.4)')rvec(ir),(vlcoup(3,l,ir),l=0,maxl)
         endif
      else
         write(20,'(12f10.4)')rvec(ir),vtran(1,1,0,ir,1)
     &,vtran(1,1,0,ir,2)
         write(21,'(12f10.4)')rvec(ir),(vls(l,ir),l=0,maxl)!MGR
         write(22,'(12f10.4)')rvec(ir),vtran(2,1,2,ir,1)
     &,vtran(1,2,2,ir,1),vtran(2,1,2,ir,2),vtran(1,2,2,ir,2)
         write(23,'(12f10.4)')rvec(ir),vtran(2,2,0,ir,1)
     &,vtran(2,2,2,ir,1),vtran(2,2,0,ir,2),vtran(2,2,2,ir,2)
      endif
      enddo
      endif
      end subroutine readpot


c *** Woods-Saxon (volume)
	function ws(r,v0,r0,a)
	implicit none
	real*8 r,v0,r0,a,ws,aux,small
        small=epsilon(small)
        ws=0d0
        if (abs(v0).lt.1e-6) return
	if (a>1e-6) then
          aux=exp(-(r-r0)/a)
          if (aux.lt.small) return
	  ws=v0/(1d0+1d0/aux)
	else
	  write(0,*)'WS: a too small!'
	endif
 	end function 

c *** Square well potential
	function sqwell(r,v0,r0)
	implicit none
	real*8 r,v0,r0,sqwell,aux,small
        small=epsilon(small)
        sqwell=0
        if (abs(v0).lt.1e-6) return
        if (r<r0) sqwell=v0
 	end function 



c *** Spin-orbit with WS derivative
c     This form factor will be then multiplied by l.s
c     We use fresco definiton for spin-orbit
        function wsso(r,v0,r0,a)
	implicit none
	real*8 r,v0,r0,a,wsso,aux,conls
	parameter(conls=2.0)
        if (r<1e-6) r=1e-6
	if (a>1e-6) then
          aux=exp(-(r-r0)/a)
	  wsso=-2d0*conls*v0*aux/(1d0+aux)**2/a/r	  
	else
	  write(0,*)'WS spin-orbit : a too small!'
         write(*,*)'v0,r0,a=',v0,r0,a
	endif

	end function 

c *** Spin-orbit with Gauss derivative
c     This form factor will be then multiplied by l.s
c     We use fresco definiton for spin-orbit
        function gausder(r,v0,r0,a)
	implicit none
	real*8 r,v0,r0,a,gausder,aux,conls,rh
	parameter(conls=2.0)
        if (r<1e-6) r=1e-6
	if (a>1e-6) then
          rh=(r-r0)/a
	  gausder=-exp(-rh**2)**rh*conls*v0/(r*a)	  
	else
	  write(0,*)'Gauss spin-orbit : a too small!'
	endif

	end function 

 
c *** WS derivative
       function dws(r,v0,r0,a)
	implicit none
	real*8 r,v0,r0,a,dws,aux
        if (r<1e-6) r=1e-6
	if (a>1e-6) then
          aux=exp(-(r-r0)/a)
	  dws=-v0*aux/(1d0+aux)**2/a	  
	else
	  write(0,*)'derivative WS: a too small!';stop
	endif

	end function 


c *** Gaussian 
	function gausspot(r,v0,r0,a)
	 implicit none
         real*8 small,big,aux,logsm,logbig
	 real*8 r,v0,r0,gausspot,a
         gausspot=0.0
         big=huge(big)
         small=epsilon(small)
         logsm=log(small)
         logbig=log(big)
         if (a.gt.1e-6) then
          aux=((r-r0)/a)**2
          if (aux.gt.logbig) return
	  gausspot=V0*exp(-(r-r0)**2/a**2)
         else
          write(*,*)'a too small in gausspot!'
          stop
         endif
         if (abs(gausspot).lt.small) gausspot=0
	end function
	
	
c *** Gaussian 
	function gausspotr2(r,v0,r0,a)
	 implicit none
         real*8 small,big,aux,logsm,logbig
	  real*8 r,v0,r0,gausspotr2,a
         gausspotr2=0.0
         big=huge(big)
         small=epsilon(small)
         logsm=log(small)
         logbig=log(big)
         if (a.gt.1e-6) then
          aux=((r-r0)/a)**2
          if (aux.gt.logbig) return
	     gausspotr2=V0*r**2*exp(-(r-r0)**2/a**2)
         else
          write(*,*)'a too small in gausspotr2!'
          stop
         endif
         if (abs(gausspotr2).lt.small) gausspotr2=0
	end function
         

c *** Yukawa
	function yukawa(r,v0,r0,a)
	implicit none
	real*8 r,v0,r0,yukawa,a
        if (r.lt.1e-4) r=1e-4
        if (a.gt.1e-6) then
	  yukawa=V0*exp(-(r-r0)/a)/r
        else
          write(*,*)'a too small in Yukawa!'
          stop
        endif
	end function



c Posch-Teller potential
      function pt(r,v0r,ar)
      implicit real*8 (a-h,o-z)
C     fOR THE DEUTERON     
C      aR=0.9407d0
C      v0R=-102.85d0
      PT=V0R/dcosh(AR*R)**2
      return
      end

c Coulomb potential
      FUNCTION VCOUL(R,z1,z2,Rc)
        use constants
        implicit none
        real*8 r,rc,rc2,aux,vcoul,z1,z2
        real*8 big,small
c       ..............................
        big=huge(big)
        small=epsilon(small)
        RC2=RC*2d0
        aux=e2*Z1*Z2
        vcoul=0
        if (z1*z2.lt.1e-4) return
        if (rc.lt.1e-6) rc=1e-6

        IF(R.LE.RC)THEN
        VCOUL=AUX*(3.-(R/RC)**2)/RC2
        ELSE
        VCOUL=AUX/R
        ENDIF
        if (vcoul.lt.small) vcoul=0 ! avoid underflow
        RETURN
      END


      subroutine potdef(vc,qmax,rnc,vdef,betar)
!      subroutine potdef(qmax,rnc,vdef,betar)
!        deform the form factor uc(l,r) by deformation def.lengths def(k),k=2,qmax
!     	 do numerical quadrature over angles
!
!	use factorials
!	use parameters
!      implicit real*8 (a-h,o-z)
!      use potentials
      use wfs, only: nr,rvec
      implicit none
      integer i,k,qmax,nu,mmultipole
!      integer ir,icall
      real*8 vc(nr)
      real*8 u,r,rsp,rnc,pi,sh,cns,alpha,fival,vcr,vdefr
      parameter(mmultipole=4,alpha=0d0)
      real*8 pl(0:mmultipole),c(0:mmultipole),sp(9),w(9),  
     &    vdef(0:mmultipole),def(0:mmultipole)
!!! TEST
      real*8 ws, p1,p2,p3,betar
!   
      data (sp(i),i=1,5) /-.96816024,-.8360311,-.61337143,-.32425342,0./  
     &,(w(i),i=1,5) /.08127439,.18064816,.26061070,.31234708,.33023936 /
       save rsp,sp,w,c
      pi=acos(-1d0)

!!!! TEST
      def(:)=0d0
      def(qmax)=betar
!!!!!

!!!! TEST
!!        write(99,*)'after sending:',rnc,vc(3)



!      if(icall.eq.1) then
!		it was never initialising c(k) so I changed this	!
!      if(l.le.0) then
! initialise
!      write(*,*)qmax,rnc,betar
      rsp = 1.0/sqrt(pi)
      do 2 i=1,4
      sp(5+i) =-sp(5-i)
2     w(5+i)  = w(5-i)
      do 15 k=0,qmax
15    c(k) = sqrt((2.*k+1.)/4.)
!	write(10,*) 'potdef initialised for qmax =',qmax
!      endif
      do 18 k=0,qmax
18    vdef(k) = 0.
         pl(0) = 1.
      do 50 nu=1,9
         u = sp(nu)
         pl(1) = u
         do 20 k=2,qmax
20       pl(k) = ((2*k-1)*u*pl(k-1) - (k-1)*pl(k-2))/dble(k)
      sh  = 0.0
      do 25 k=2,qmax
25       sh = sh + c(k)*rsp * pl(k) * def(k)
         r = rnc - sh
!!!!! TEST
!56.58  r0=2.483   a0=0.65  
         vcr=fival(r,rvec,vc,nr,alpha)
!         p1=ws(r,-54.239d0,2.483d0,0.65d0)
!         p1=ws(r,-56.00d0,2.483d0,0.65d0)
!         write(97,*) r,vcr,p1
!!!! 
	 do 40 k=0,qmax
!	    if(k.ge.2) then
!	    	if(def(k).eq.0d0) go to 40
!	    endif
            cns = w(nu) * pl(k) * c(k)
!          vdef(k)=vdef(k) + cns * vc(l,r,ib)
          vdef(k)=vdef(k) + cns * vcr
40      continue
50      continue
       vdefr=vdef(qmax)
	return
      end
