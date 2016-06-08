      subroutine readpot(kin,kpot,zfrag)
c *** Read fragment-target potential parameters and store the potential in the radial grid
c  Only CENTRAL potentials are allowed so far, i.e.:
c
c    V(r)=Vcou(r)+ Vc(r)
c
c    Vcou(r) = coulomb central
c    Vc(r)   = nuclear central potential
       use ptpots
       use sistema
       use wfs, only:rvec,nr
       use globals , only: written,verb
       implicit none
       integer kin,ir,ncomp,kpot,iq
       real*8 a13,r,rc,rabs,rabsi,zfrag
       real*8:: vaux,vcoul,ws,pt,gausspot,vlsaux,wsso,dws,vllaux
       real*8:: vdef(0:nmult),vcoulq,vcaux(nr),wcaux(nr)
       real*8,allocatable,target::vcl(:),vcli(:)
       character*40:: potfile
       character(len=10):: fragname
       integer::ptype,cptype,lambda    
	   integer np,nv,iv!,maxlamb    !Lay:added for read external potential   
	   real*8,allocatable:: veaux(:),weaux(:),raux(:)!Lay:added for read external potential   
	   real*8:: fival
       real*8:: ap,at             ! masses for radius conversion
       real*8:: V0,R0,A0
       real*8:: V0i,R0i,A0i
       real*8:: Vso,rso,aso       ! spin-orbit
       real*8:: Vss0,rss,ass       ! spin-spin
       real*8:: Vll0,rll,all       ! l.l
       real*8:: rc0               ! charge radius
       real*8:: beta,delta,Vcp0,rcp0,acp !coupling
       real*8:: betai,deltai,Vcp0i,rcp0i,acpi

       
       
c       namelist /potential/ ptype,ap,at,V0,vl0,r0,a0,rc0,
c     & Vso,rso,aso,
c     & Vss0,rss,ass,
c     & Vll0, rll,all,
c     & cptype,lpot,beta,delta,Vcp0,rcp0,acp
       namelist /corepotential/ ptype,ap,at,V0,r0,a0,rc0,
     & cptype,Vcp0,rcp0,acp,delta,beta,
     & V0i,r0i,a0i,
     & Vcp0i,rcp0i,acpi,deltai,betai,mel,
     & potfile,np,nv
       namelist /valencepotential/ ptype,ap,at,V0,r0,a0,rc0,
     & cptype,Vcp0,rcp0,acp,delta,beta,
     & V0i,r0i,a0i,
     & Vcp0i,rcp0i,acpi,deltai,betai,mel,
     & potfile


       ncomp=0 !number of potential components

       if (.not.allocated(vcl)) allocate(vcl(1:nr))
       if (.not.allocated(vcli)) allocate(vcli(1:nr))

       vcoup(:,:)=0d0;
       vfrag(:,:)=0.d0;    vfragi(:,:)=0.d0

c100    V0=0; vl0(:)=0d0; r0=0d0; a0=0d0;
c       Vso=0d0; rso=0d0; aso=0d0; rc0=0d0
c       Vcp0=0d0; rcp0=0d0; acp=0d0
c       Vss0=0d0; rss=0d0; ass=0d0
c       Vll0=0d0; rll=0d0; all=0d0
100    V0=0; r0=0d0; a0=0d0; rc0=0d0;
       Vcp0=0d0; rcp0=0d0; acp=0d0;
       V0i=0; r0i=0d0; a0i=0d0;
       Vcp0i=0d0; rcp0i=0d0; acpi=0d0
       mel=0d0;     
       delta=0d0
       beta=0d0
       deltai=0d0
       betai=0d0
       vcl(:)=0.d0;        
       vcli(:)=0.d0; 

     
c       lambda=2
       written(40)=.false.
       written(42)=.false.

       if(kpot.eq.0) then
        read(kin,nml=corepotential)
        fragname='core' 
       else
        read(kin,nml=valencepotential)
        fragname='valence'
       endif
       
       A13=Ap**.33333333 + At**.33333333
       rabs=r0*a13
       rabsi=r0i*a13
       if (delta.eq.0d0) delta=beta*rabs
       if (deltai.eq.0d0) deltai=betai*rabsi

       if (ptype.eq.0) then
	if (ncomp.eq.0) then
	   write(*,*)'No potentials were found!'
           stop
        endif
       endif
       ncomp=ncomp+1
!       write(*,*)'number of potentials=',ncomp
c       if (abs(vl0(0))>0d0) v0=vl0(0)
       if (ptype.ne.0) write(*, '(/,2x,"o Central potential",$)') 

       select case(ptype) ! Central potential
c ------------------------------------------------------
        case(0)  ! Last potential already read
c ------------------------------------------------------
        if (written(40)) then
         write(40,'("# Central potential for:",a10)') fragname
         do ir=1,nr
         if (rvec(ir).gt.10) exit
         write(40,2000)rvec(ir),vfrag(ir,0),vfragi(ir,0),vcoup(ir,0)
         enddo
         write(40,*)'& '
        endif

        if (written(42).and.(nmult.ge.2)) then
         write(42,'("# Deformed potential for:",a10)') fragname
         do ir=1,nr
         if (rvec(ir).gt.10) exit
         write(42,2000)rvec(ir),vfrag(ir,2),vfragi(ir,2),vcoup(ir,2)
         enddo
         write(42,*)'& '
        endif
2000    format(1x,1f7.3,2x,3f12.6)


	   return
c ------------------------------------------------------
        case(1,2)  ! WS / WS derivative types
c ------------------------------------------------------
          if (ptype.eq.1) write(*,*)'- type=1 => Woods-Saxon'
          if (ptype.eq.2) write(*,*)'- type=2 => Woods-Saxon derivative'
           write(*,120) V0,r0*a13,a0,V0i,R0i*a13,A0i,Rc0*a13
c           write(*,122) Vso,rso*A13, aso
 120       format(3x,'[V0=',1f7.2,' MeV, Rr=',f5.2,
     &           ' fm, ar=',f5.2,' fm, ',
     &           ' W0=',1f7.2,' MeV, Ri=',f5.2,
     &           ' fm, ai=',f5.2,' fm, '
     &           ' Rc=',f5.2,' fm ]')
c 122       format(3x,'[Vso=',1f8.3,' MeV, Rso=',f6.3,
c     &           ' fm, aso=',f6.3,' fm ]',/)
!c ------------------------------------------------------
!c        case(2) !  Posch-Teller
!c ------------------------------------------------------
!c           write(*,*)'- type=2 => Posch-Teller'
!c           write(*,130) V0,a0,rc0
!c 130       format(3x,'[V=',1f8.3,' MeV, ar=',f6.3,
!c     &           ' fm, rc0=',f6.3,' fm ]',/)
!c ------------------------------------------------------
!c        case(3) ! Gaussian 
!c ------------------------------------------------------
!c             write(*,*)'- type=3 => Gaussian'
!c             write(*,140) V0,r0,rc0
!c 140       format(3x,'[V=',1f8.3,' MeV  r0r=',f6.3,
!c     &           ' fm  rc0=',f6.3,' fm ]',/)
!c	    if (vg1>1e-3) then                                       
!c        	 write(*,142) vg1,rg1
!c 142       format(3x,'[V2=',1f8.3,' MeV  r0r=',f6.3,
!c     &           ' fm ]',/)
!c            endif

c ------------------------------------------------------
        case(3) 
c ------------------------------------------------------
             write(*,*)'- type=3 => Gaussian'
             write(*,140) V0,r0,rc0
 140       format(3x,'[V=',1f8.3,' MeV  r0r=',f6.3,
     &           ' fm  rc0=',f6.3,' fm ]',/)
c ------------------------------------------------------
        case(4)
c ------------------------------------------------------
           write(*,*)'- type=4 => Posch-Teller'
           write(*,160) V0,a0,rc0
 160       format(3x,'[V=',1f8.3,' MeV, ar=',f6.3,
     &           ' fm, rc0=',f6.3,' fm ]',/)

        case(5)
           write(*,*)'-type=5=> External potential'
		   write(*,*)'			Reading ',np,'points'
c ------------------------------------------------------
        case default
c ------------------------------------------------------
          write(*,*)'Potential type=',ptype,' not used'
          stop
        end select
 150    format(6f10.5)


c coupling component, if any
       if (abs(delta).lt.1e-4) goto 200
       write(*, '(/,2x,"o Deformed potential",$)')
       select case(cptype)
c ------------------------------------------------------
        case(1)  ! dV/dr type
c ------------------------------------------------------
	if (ptype.eq.1) then
           write(*,*)'- Coupling Pot. type=1 => deriv. W-S'
!	   write(*,*)'- WARNING!!!! => Possible 4pi factors not included'
           if (Vcp0.eq.0) then
           write(*,141) delta
	   else
	   if(rcp0.eq.0) rcp0=r0
           if(rcp0i.eq.0) rcp0i=r0i
	   if(acp.eq.0)  acp=a0
           if(acpi.eq.0)  acpi=a0i
	   write(*,143)Vcp0,acp,rcp0
	   write(*,*)'- We are NOT using DELTA!!'
	   endif
 143       format(3x,'[Vcp=',1f8.3,' MeV acp=',f6.3,'  rcp0=',f6.3,
     &           ' fm ]',/)
 141       format(6x,'[delta=',1f6.3,' fm ]',/)
	else
	   write(*,*)'- Coupling Pot. type=1 => deriv. PT/Gaussian'
	   write(*,*)'- WARNING!!!! => not implemented yet'
        endif
c ------------------------------------------------------
        case(2,3)  ! Projection on multipoles
c ------------------------------------------------------
	  write(*,*)
     &       ' => deform central potential & project on P_{lambda}'
         if (cptype.eq.3) write(*,*)
     &       '     (monopole potential will be also recalculated!)'
           write(*,141) delta
c ------------------------------------------------------
        case(4)  ! V_2 type
c ------------------------------------------------------
	   write(*,*)'- Coupling Pot. type=2 => deformed fermi function'
	   write(*,*)'- WARNING!!!! => work in progress'                      !remember to erase!
	   write(*,144)beta
 144     format(3x,'[beta=',1f8.3,']',/)
	end select
	    

c Coulomb potential
200   Rc=Rc0*A13
      if (zfrag.gt.1e-6) then
       if (maxval(abs(mel)).gt.1e-4) then
          write(*,'(2x,"o Coulomb deformation: Mn [e.fm^q]=",10f6.2)')
     &    mel(1:nmult)
          written(42)=.true.
       endif
	if (ncomp.eq.1) then  !! coulomb potential should not be added to previous one
        do ir=1,nr
        r=rvec(ir)
        vcoup(ir,0)=VCOUL(r,zfrag,zt,Rc) ! monopole
        do iq=1,nmult
        vcoup(ir,iq)=vcoulq(r,rc,zt,iq,mel(iq))  ! multipole
        enddo !iq
        enddo !ir
      endif
      endif

c Central NUCLEAR potential
        if (abs(v0).gt.1e-6) then
        written(40)=.true.
        do ir=1,nr
           r=rvec(ir)
            select case(ptype)
            case(1)  ! WS type
		vaux=ws(r,v0,rabs,a0)
            case(2)
		vaux=-4.d0*a0*dws(r,v0,rabs,a0)
            case(3) 
 	        vaux=gausspot(r,v0,r0)
            case(4)
		vaux=pt(r,v0,rabs,a0)
	    end select
              vcl(ir)=vcl(ir)+vaux
        enddo !ir
        endif
                

        if (abs(v0i).gt.1e-6) then
        written(40)=.true.
        do ir=1,nr
           r=rvec(ir)
            select case(ptype)
            case(1)
		vaux=ws(r,v0i,rabsi,a0i)
            case(2)
		vaux=-4.d0*a0i*dws(r,v0i,rabsi,a0i)
            case(3) 
 	        vaux=gausspot(r,v0i,r0i)
            case(4)
		vaux=pt(r,v0i,rabsi,a0i)
	    end select
              vcli(ir)=vcli(ir)+vaux
        enddo
        endif




c! Spin-orbit
c       if (abs(vso)>1e-5) then
c       do ir=1,nr
c           r=rvec(ir) 
c           if (r.lt.1e-6) r=0.5*dr
c            select case(ptype)
c            case(1)  ! WS type
c              vlsaux=wsso(r,vso,rso*a13,aso)		
c            case(2) !  Posch-Teller
c!		vaux=pt(r,v0,rabs,a0)
c            case(3) ! Gaussian 
c! 	        vaux=gausspot(r,v0,r0)
c	    end select
c	    vls(ir)=vls(ir)+vlsaux 
c!	    if (abs(vls(ir))>10000) write(0,*) vls(ir)  
c        enddo !ir
c        endif


c! l.l (only central!!!!) 
c       if (abs(vll0)>1e-5) then
c       do ir=1,nr
c           r=rvec(ir) 
c           if (r.lt.1e-6) r=0.5*dr
c            select case(ptype)
c            case(1)  ! WS type
c              vllaux=ws(r,vll0,rll*a13,all)		
c            case(2) !  Posch-Teller
c!		vaux=pt(r,v0,rabs,a0)
c            case(3) ! Gaussian 
c! 	        vaux=gausspot(r,v0,r0)
c	    end select
c	    vll(ir)=vll(ir)+vllaux 
c!	    if (abs(vll(ir))>1) write(0,*) vll(ir)  
c        enddo !ir
c        endif

            
    
c Non-central (deformed) part (REAL)
        vcaux(:)=vcl(:)
        if (abs(delta).gt.1e-4) then 
        do ir=1,nr
        r=rvec(ir)
        vdef(:)=0        
        select case(cptype)
             case(0) ! no coupling
             case(1) ! derivative formfactor
		 if (Vcp0.eq.0) then
                vdef(2)=-delta*dws(r,v0,  rabs,    a0)
		 else
                vdef(2)=-delta*dws(r,Vcp0,rcp0*a13,acp)
	        endif

              case(2,3) ! Numerical projection on multipoles
                 call potdef(vcl(:),nmult,r,vdef,delta)
            
              case default
                  print*,'coupling potential',cptype,'not implemented!'
                  stop
		end select

              if  (cptype.eq.3) vcaux(ir)=vdef(0)

              vfrag(ir,1:nmult)=vfrag(ir,1:nmult) + vdef(1:nmult)  
              enddo !ir
        endif ! delta>0

c Non-central (deformed) part (IMAGINARY)
        wcaux(:)=vcli(:)
        if (abs(deltai).gt.1e-4) then 
        written(42)=.true.
        do ir=1,nr
        r=rvec(ir)
        vdef=0
        select case(cptype)
             case(0) ! no coupling
             case(1) ! derivative formfactor
!!!!! CHECK THIS!!!!
		   if (Vcp0i.eq.0) then
		    vdef(2)=-deltai*dws(r,v0i,rabsi,a0i)
		   else
		    vdef(2)=-deltai*dws(r,Vcp0i,rcp0i*a13,acpi)
	          endif

             case(2,3) ! Numerical projection on multipoles
                 call potdef(vcli(:),nmult,r,vdef,deltai)
                				
             case default
                 print*,'coupling potential',cptype,'not implemented!'
                 stop
	      end select

             if  (cptype.eq.3) wcaux(ir)=vdef(0)
             vfragi(ir,1:nmult)=vfragi(ir,1:nmult) + vdef(1:nmult)    
             enddo !ir
        endif ! deltai > 0

       vfrag (1:nr,0)=vfrag (ir,0) + vcaux(1:nr)
       vfragi(1:nr,0)=vfragi(ir,0) + wcaux(1:nr)
	   
	   
	           if (ptype.eq.5) then !up to now only central part
!			 case(5) !read from external fort.32 !work in progress
                vcaux(:)=0d0
				wcaux(:)=0d0
		        allocate(veaux(np),weaux(np),raux(np))
			    veaux(:)=0.
				weaux(:)=0.
				read(32,*)
			    do ir=1,np	
				read(32,*)raux(ir),veaux(ir),weaux(ir)
				enddo
				do ir=1,nr
				vcaux(ir)=fival(rvec(ir),raux,veaux,np,0.)
				wcaux(ir)=fival(rvec(ir),raux,weaux,np,0.)
				enddo!          written(40)=.true.
				
			           vfrag (1:nr,0)=vfrag (ir,0) + vcaux(1:nr)
                       vfragi(1:nr,0)=vfragi(ir,0) + wcaux(1:nr)
					   
				vcaux(:)=0d0
				wcaux(:)=0d0
				veaux(:)=0.
				weaux(:)=0.
				read(32,*)             ! A second loop for deformed part
			    do ir=1,np	
				read(32,*)raux(ir),veaux(ir),weaux(ir)
				enddo
				do ir=1,nr
				vcaux(ir)=fival(rvec(ir),raux,veaux,np,0.)
				wcaux(ir)=fival(rvec(ir),raux,weaux,np,0.)
				enddo!
					   vfrag (1:nr,2)=vfrag (ir,2) + vcaux(1:nr)
                       vfragi(1:nr,2)=vfragi(ir,2) + wcaux(1:nr)
					   	!empirical squared root of 4pi to check   
					   
!          call extpot(potfile,vcl,vcli,nr) ! I am not sure if this is the best option
              endif
	   
	   
	   

	goto 100 !read next potential component
      return
      end subroutine readpot


c Read external potential
c Format: 
c 1st line: header
c 2nd line: npoints, rstep, rfirst
c Next lines: real, imag
      subroutine extpot(filename,vr,vi,nr)
        use wfs, only: rvec
        implicit none
        logical uu
        integer npt,n,nr
        real*8 vr(nr),vi(nr)
        integer,parameter:: kpot=50
        character*40 filename,header
        real*8:: fival,r,rstep,rfirst,x,y1,y2
        real*8,allocatable::xv(:),faux(:),gaux(:)
        real,parameter:: alpha=0

        uu = .false.
        inquire(file=filename,exist=uu)
        if (.not.uu) then
          write(0,*) 'Potential file:', filename,' not found!'
          stop
        endif
        write(*,'(8x, "Potential file:",a20)') filename 
        open(kpot,file=filename,status='old')
        read(kpot,*) header
        npt=0
20      read(kpot,*,end=200) x,y1,y2 
        npt=npt+1
        goto 20

200     if (npt.lt.1) goto 300        
        rewind(kpot)
        allocate(faux(npt),gaux(npt),xv(npt)) 
        faux(:)=0; gaux(:)=0; xv(:)=0      
        read(kpot,*) header        
!        read(kpot,*) npt,rstep,rfirst 
        do n=1,npt
           read(kpot,*,err=300)  x,y1,y2
           xv(n)=x
           faux(n)=y1
           gaux(n)=y2
        enddo

      write(*,'(3x,"=> read:",i3," points")') npt 
      write(*,250) xv(1),xv(npt),xv(2)-xv(1)
250   format(/,5x,"[Radial grid: Rmin=",1f6.3," fm,", 
     &      " Rmax=",1f6.1," fm," 
     &      " Step=",1f6.3," fm]",/)

        do n=1,nr
          r=rvec(n)
          vr(n)=fival(r,xv,faux,npt,alpha)
          vi(n)=fival(r,xv,gaux,npt,alpha)
        enddo

        deallocate(xv,faux,gaux)
        return
300     write(*,*)'error reading ext. potential !'
        stop
        return
      end subroutine 
         



c *** Woods-Saxon (volume)
	function ws(r,v0,r0,a)
	implicit none
	real*8 r,v0,r0,a,ws,aux
       ws=0
       if (abs(v0).lt.1e-6) return
	if (a<1e-6) then
         write(0,*)'WS: a too small!'
         stop
       endif
       aux=exp(-(r-r0)/a)
	ws=v0/(1d0+1d0/aux)
       return
 	end function 

c *** Spin-orbit with WS derivative
c     This form factor will be then multiplied by l.s
c     We use fresco definiton for spin-orbit
      function wsso(r,v0,r0,a)
	implicit none
	real*8 r,v0,r0,a,wsso,aux,conls
	parameter(conls=2.0)
       v0=0
       if (abs(v0)<1e-6) return
       if (r<1e-6) r=1e-6
	if (a>1e-6) then
         aux=exp(-(r-r0)/a)
	  wsso=-2d0*conls*v0*aux/(1d0+aux)**2/a/r	  
	else
	  write(0,*)'WS spin-orbit : a too small!'
	endif
       return
      end function 
 
c *** WS derivative
      function dws(r,v0,r0,a)
	implicit none
	real*8 r,v0,r0,a,dws,aux
       dws=0
       if (abs(v0).lt.1e-6) return
       if (r<1e-6) r=1e-6
	if (a>1e-6) then
        aux=exp(-(r-r0)/a)
	  dws=-v0*aux/(1d0+aux)**2/a	  
	else
	  write(0,*)'derivative WS: a too small!';stop
	endif
       return
	end function 


c *** Gaussian 
	function gausspot(r,v0,r0)
	 implicit none
	 real*8 r,v0,r0,gausspot
        gausspot=0
        if ((r0.lt.1e-6).or.(abs(v0).lt.1e-6))return
	 gausspot=V0*dexp(-r**2/r0**2)
        return
	end function
         

c Posch-Teller potential
C  For the deuteron, you may use:
C      aR=0.9407d0
C      v0=-102.85d0

      function pt(r,v0,ar)
      implicit real*8 (a-h,o-z)
      pt=0
      if ((abs(v0).lt.1e-6))return
      PT=V0/dcosh(AR*R)**2
      return
      end
c
c Monopole Coulomb potential between charges Z1 and Z2 
c 
        FUNCTION VCOUL(R,z1,z2,Rc)
          use constants,only: e2
          implicit none
          real*8 r,rc,rc2,aux,vcoul,z1,z2
          real*8 eps 
          parameter(eps=1e-3)
          RC2=RC*2d0
          aux=e2*Z1*Z2
          vcoul=0
          if (z1*z2.lt.1e-4) return
          if (r.lt.eps) r=eps

          IF(R.GT.RC)GO TO 1
          VCOUL=AUX*(3.-(R/RC)**2)/RC2
          RETURN
1         VCOUL=AUX/R
          RETURN
        END

c -------------------------------
c q>0 Coulomb potential   
c -------------------------------
      function vcoulq(r,rc,zt,q,mn)
      use constants, only: e2
      implicit none
      integer:: q
      real*8:: r,rc,vcoulq,rq,zt,mn
      real*8 eps,aux,pi
      parameter(eps=1e-3)
      vcoulq=0d0
      pi=acos(-1d0)
      rq=q
c  we use e^2 here since this includes thef 'e' 
c  coming from the Mn since this is 
c  typically given in units of e.fm^q
! old
!      aux=zt*e2/dsqrt(2d0*q+1)
! new 
      aux=mn*zt*e2*sqrt(4*pi)/(2d0*q+1)
      if (r.lt.eps) r=eps
      if (r.lt.rc) then
         vcoulq=aux*r**rq/rc**(2d0*q+1d0)
      else
         vcoulq=aux/r**(q+1d0)
      endif
      return
      end function


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
      real*8 u,r,rsp,rnc,pi,sh,cns,alpha,fival,vcr!,vdefr

!      parameter(mmultipole=4,alpha=0d0)
!      real*8 pl(0:mmultipole),c(0:mmultipole),sp(9),w(9),  
!     &    vdef(0:mmultipole),def(2:mmultipole)
      parameter(alpha=0d0)
      real*8, dimension(0:qmax):: pl,vdef,def,c
      real*8 sp(9),w(9) 
!!! TEST
      real*8 ws, p1,p2,p3,betar
!   
      data (sp(i),i=1,5) /-.96816024,-.8360311,-.61337143,-.32425342,0./  
     &,(w(i),i=1,5) /.08127439,.18064816,.26061070,.31234708,.33023936 /
!       save rsp,sp,w,c
       save rsp,sp,w
      pi=acos(-1d0)

!!!! TEST
      def(:)=0d0
      def(2)=betar
!!!!!

!      if(icall.eq.1) then
!		it was never initialising c(k) so I changed this	!
!      if(l.le.0) then
! initialise
      
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
         vcr=0d0
         if (r.le.rvec(nr)) then 
         vcr=fival(r,rvec,vc,nr,alpha)
         endif
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
c       vdefr=vdef(qmax)
	return
      end
