c Read sp configurations for Pauli-fobidden states
       subroutine pauli_forbidden(kin)
       use forbidden
       use parameters, only:maxcore
       use wfs,only:nr 
       implicit none
       integer kin,n,l
       real*8:: j,pfactor(maxcore),u,v,equasi
       real*8,pointer::wfaux(:)
       character*40:: wfname
       CHARACTER LNAME(0:14)
       DATA LNAME / 's','p','d','f','g','h','i','j','k','l','m',
     &                'n','o','p','q' /
       namelist /pauli/ n,l,j,wfname,wblock,hindrance,pfactor,
     &          u,v,equasi
c      ---------------------------------------------------
       eshift=1000
       pshift(:,:)=eshift
       pscale(:,:)=1
       write(*,*)'  '
       wfpread(:)=.false.
       hindrance=.false.
10     n=0; l=0;j=0.0
       wfname=""
       wblock(:)=eshift
       pfactor(:)=1
       u=1; v=0;
       equasi=0
c      -------------------------------------------------
       write(*,*)' PAULI forbidden sp states:'
       read(kin,nml=pauli)
       if (n.eq.0) goto 20
       npauli=npauli+1
       if (npauli.gt.maxpauli) then
         write(*,*)'Too many Pauli forbidden states!'
         write(*,*)'Increase maxpauli in modules.f90'
         stop
       endif
       paul(npauli)=l
       paun(npauli)=n
       pauj(npauli)=j
       pshift(npauli,:)=wblock(:)
! v2.0.4 ------------------------
       pscale(npauli,:)=pfactor(:)
       uocup(npauli)=u
       vocup(npauli)=v
       eqsp(npauli)=equasi
c ------------------------------
       write(*,30) npauli,n,lname(l),NINT(J*2.)
30     format(5x,"#" i2,":",1x,i2,a1,I1,'/2',$)
       if (wfname.eq."") then 
        write(*,*)' (calculated from Hsp diagonalization)'
       else
        if (.not.allocated(wfpauli)) allocate(wfpauli(maxpauli,nr))
        wfpread(npauli)=.true.
        paulifile(npauli)=wfname
        write(*,'(5x,"(read from file  ",a15,")")')wfname
	 wfaux=>wfpauli(npauli,:)
        call readwf(wfname,wfaux,nr)
       endif
	goto 10
20     if (npauli>0) then
       write(*,*) '    => So',npauli,' state(s) will be removed' 
       
!        write(*,*) '    => So',npauli,' state(s) will be removed'
       else
           write(*,*) '=> No states will be removed'
       endif
       write(*,*) ' '
       end subroutine


c +-------------------------------------------------------
c Read an external wavefuncion u(r) and interpolates in rvec() grid
c --------------------------------------------------------
      subroutine readwf(wfname,wfvec,ndim)
        use wfs, only: nr,dr,rvec
        use channels, only: bas2
        use parameters
        implicit none
        integer nptmax,ndim
        real*8:: eps,alpha
        parameter (nptmax=5000,eps=1e-3)    
        parameter(alpha=0d0)
        logical::uu
        character*40::wfname
        integer npt,ir
        real*8:: rstart,r,xmax,fival,res,dx
        real*8 ::yaux(nptmax),xaux(nptmax),faux(ndim)
        real*8:: wfvec(ndim),uaux
	 wfvec(:)=0d0
        faux(:)=0d0
        write(*,*) '- Reading WF from file: ',wfname
        uu = .false.
        inquire(file=wfname,exist=uu)
        if (.not.uu) then
           write(*,*)'Could not read file:',wfname
           write(*,*)'Aborting...'
           stop
        endif
        open(50,file=wfname,status='old')
        ir=1
	faux(:)=0d0
10      read(50,*,end=80)  xaux(ir), yaux(ir) 
!	faux(ir)=xaux(ir)**2*(yaux(ir))**2
!	gaux(ir)=xaux(ir)**4*(yaux(ir))**2
        ir=ir+1
        if (ir>nptmax) then
           write(*,*)'too many points in external function'  
           write(*,*)'increase variable nptmax'
           stop
        endif
        goto 10
80      npt=ir-1
        write(*,*)'read',npt,'points'
        xmax=xaux(npt)
	dx=xaux(2)-xaux(1)
 
c interpolate in rvec() grid (and avoid extrapolation)
!       allocate(uref(nr))
!       uref(1:nr)=0d0
       
       faux(:)=0d0
       do ir=1,nr
          r=rvec(ir)
          if (r.lt.eps) r=eps
          if (r.le.xmax) then
!            uref(ir)=fival(r,xaux,yaux,npt,alpha)
            uaux=fival(r,xaux,yaux,npt,alpha)
	     faux(ir)=uaux**2
            wfvec(ir)=uaux/r
!             write(90,*)r,uaux
          endif
       enddo!	call sim(faux,res,1,npt,dx,npt)
!	write(*,*)'    <wf_ref| r^2 | wf_ref>=',res
!	call sim(gaux,res,1,npt,dx,npt)
!	write(*,*)'    <wf_ref| r^4 | wf_ref>=',res

       call sim(faux,res,1,nr,dr,nr)
       write(*,*)'    Norm of WF=',res
       return
90      write(*,*)'End of File found while reading',wfname
      end subroutine readwf



c 
c *** Pauli projection operator for full matrix
c     P=Sum |forbidden> < forbidden| 
	subroutine pauliproj(ndim)
	use forbidden
	use globals, only: written
	use wfs !, only: wfsp,ql,qj,spindex,nr,rvec
        use channels
	implicit none
	integer:: ichi,ichf,ni,nf,ir,ncorei,ncoref
	integer:: lp,np
	real*8:: jp
	integer:: ip,n,m,ndim,ispi,ispf,li,lf,jci,jcf
	real*8:: ji,jf,un,um,r,up,res1,res2
	real*8,allocatable:: wfpau(:)
	real*8,allocatable:: fpaux(:),gpaux(:)
	real*8,allocatable:: fmaux(:),fnaux(:)
	real*8 iden(hdim,hdim)
	real*8::pblock(ndim,ndim),paux
	 
      
       if (allocated(fnaux)) deallocate(fnaux)
       if (allocated(fmaux)) deallocate(fmaux)
       if (allocated(fpaux)) deallocate(fpaux)
       if (allocated(gpaux)) deallocate(gpaux)
       if (allocated(wfpau)) deallocate(wfpau)
       if (allocated(ppauli)) deallocate(ppauli)
       if (allocated(qpauli)) deallocate(qpauli)
!	write(*,*)'pauli: nr,ndim,hdim',nr,ndim,hdim
	allocate(fpaux(nr),gpaux(nr))
	allocate(fnaux(nr),fmaux(nr))
	allocate(wfpau(nr))
	allocate(ppauli(hdim,hdim),qpauli(hdim,hdim))
	
	fnaux(:)=0d0; fmaux(:)=0d0
	ppauli(:,:)=0d0
	written(80)=.true.

	
      iden(:,:)=0d0 !Identity matrix
      write(*,*)'hdim=',hdim
	do n=1,hdim
	iden(n,n)=1d0
      enddo

      do ichi=1,nchan
       do ichf=1,nchan
	pblock(:,:)=0d0
       if (bas2.eq.1) then  
       	wfaux=>wftho(:,ql(ichi),:)
              wfaux2=>wftho(:,ql(ichf),:)
       else if (bas2.eq.0) then
       	wfaux=>wfsp(spindex(ichi),:,:)
       	wfaux2=> wfsp(spindex(ichf),:,:)
       else 
		write(*,*)'bas2=',bas2,' not implemented!'
              stop
       endif
       do n=1,ndim
	ni=(ichi-1)*ndim+n
	do m=1,ndim
	nf=(ichf-1)*ndim+m
!	write(*,*)ni,nf
	
       ppauli(ni,nf)=0d0; paux=0d0
       ispi=spindex(ichi)
       ispf=spindex(ichf)
       li=ql(ichi)	
       lf=ql(ichf)
       ji=qj(ichi)
       jf=qj(ichf)
       jci=qjc(ichi)
       jcf=qjc(ichf)
       ncorei=cindex(ichi)
       ncoref=cindex(ichf)
	
	
	do ip=1,npauli
	lp=paul(ip) 
        jp=pauj(ip)
	np=paun(ip)
        if (wfpread(ip)) then
          wfpau(:)=wfpauli(ip,:)
        else
          wfpau(:)=wfsp(ispi,paun(ip),:) !Pauli forbidden 
        endif

         fpaux(:)=0d0; gpaux(:)=0d0
        
	 if (li.ne.lp) cycle
         if (lf.ne.lp) cycle
	 if (ji.ne.jp) cycle
         if (jf.ne.jp) cycle
	 if (jci.ne.jcf) cycle
         if (ncorei.ne.ncoref) cycle

         do ir=1,nr
          r=rvec(ir)
          un=wfaux(n,ir)*r
          um=wfaux2(m,ir)*r
	   up=wfpau(ir)*r
          fpaux(ir)=fpaux(ir)+un*up ! <n|p>
	   gpaux(ir)=gpaux(ir)+um*up ! <p|m>
         enddo !ir
!         endif
	 call sim(fpaux,res1,1,nr,dr,nr)
        call sim(gpaux,res2,1,nr,dr,nr)
!!!!!! TEST PAULI HINDRANCE
!	 ppauli(ni,nf)=ppauli(ni,nf)+res1*res2
        paux=paux+res1*res2
	 ppauli(ni,nf)=ppauli(ni,nf)+pshift(ip,ncorei)*res1*res2
!         qpauli(ni,nf)=iden(ni,nf)-ppauli(ni,nf)
!         ppauli(ni,nf)=ppauli(ni,nf)*pshift(ip,ncorei)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 enddo !ip
        qpauli(ni,nf)=iden(ni,nf)-paux
	 pblock(n,m)=ppauli(ni,nf)
        
!!!!!! TEST HINDRANCE
!!         ppauli(ni,nf)=ppauli(ni,nf)*pshift(ip,ncorei)
!!!!!!!!!
        enddo !m
        enddo !n
	 write(80,*)'# Pauli projector for block',ichi,ichf
        call writemat(pblock,ndim,80)
	 enddo !ichf
        enddo !ichi
        
c write block
!	 write(80,*)'# Pauli projector'
!         call writemat(ppauli,hdim,80)
         !ip=1,npauli
	
c Q= 1 - P = 1 - Sum |forbidden> <forbidden|
!	 qpauli=iden-ppauli
	end subroutine 


c 
c *** Pauli projection operator for {nchani,nchanf} block
c     P=Sum |forbidden> < forbidden| 
	subroutine pauliproj2(nchani,nchanf,ndim)
	use forbidden
	use globals, only: written
	use wfs !, only: wfsp,ql,qj,spindex,nr,rvec
        use channels
	implicit none
	integer:: nchani,nchanf
	integer:: ip,n,m,ndim,ichn,ichm,li,lf,ir
	real*8:: ji,jf,un,um,r,up,res1,res2
	real*8,allocatable:: wfpau(:)
	real*8,allocatable:: fpaux(:),gpaux(:)
	real*8,allocatable:: fmaux(:),fnaux(:)
	real*8 iden(ndim,ndim)
	 
       ichn=spindex(nchani)
       ichm=spindex(nchanf)
       li=ql(nchani)	
       lf=ql(nchanf)
       ji=qj(nchani)
       jf=qj(nchanf)
       if (allocated(fnaux)) deallocate(fnaux)
       if (allocated(fmaux)) deallocate(fmaux)
       if (allocated(fpaux)) deallocate(fpaux)
       if (allocated(gpaux)) deallocate(gpaux)
       if (allocated(wfpau)) deallocate(wfpau)
       if (allocated(ppauli)) deallocate(ppauli)
       if (allocated(qpauli)) deallocate(qpauli)
!	write(*,*)'pauli: nr,ndim',nr,ndim
	allocate(fpaux(nr),gpaux(nr))
	allocate(fnaux(nr),fmaux(nr))
	allocate(wfpau(nr))
	allocate(ppauli(ndim,ndim),qpauli(ndim,ndim))
	
	
	fnaux(:)=0d0; fmaux(:)=0d0
	ppauli(:,:)=0d0
	iden(1:ndim,1:ndim)=0d0

	do ip=1,npauli
        wfpau(:)=wfsp(ichn,paun(ip),:) !Pauli forbidden WF

	 do n=1,ndim
	  do m=1,ndim
	  if (n.eq.m) iden(n,n)=1d0
         fpaux(:)=0d0; gpaux(:)=0d0
         ppauli(n,m)=0d0
        if (li.eq.paul(ip).and.(ji.eq.pauj(ip))
     &  .and.(lf.eq.paul(ip)).and.(jf.eq.pauj(ip)) 
     &  .and.(nchani.eq.nchanf)) then
!     &  .and.(n.eq.paun(ip))) then
	  do ir=1,nr
          r=rvec(ir)
          un=wfaux(n,ir)*r
          um=wfaux2(m,ir)*r
	  up=wfpau(ir)*r
          fpaux(ir)=fpaux(ir)+un*up ! <n|p>
	   gpaux(ir)=gpaux(ir)+um*up ! <p|m>
          enddo !ir
         endif

	 call sim(fpaux,res1,1,nr,dr,nr)
        call sim(gpaux,res2,1,nr,dr,nr)
	 ppauli(n,m)=ppauli(n,m)+res1*res2
        enddo !n
        enddo !m

c write block
	 written(80)=.true.
	 write(80,*)'# Pauli projector for ichan block',nchani,nchanf
         call writemat(ppauli,ndim,80)
        enddo !ip=1,npauli
	
c Q= 1 - P = 1 - Sum |forbidden> <forbidden|
	 qpauli=iden-ppauli
	

	end subroutine 


