c 
c *** Calculate valence/core configurations compatible with Jtot/pi
c     nchan={Ic,lsj} configurations  
      subroutine read_jpiset(iset,bastype,nk,tres,ehat,filewf)
      use channels
      use globals,only: kin
      use wfs    ,only: exmin,exmax
      use potentials, only: vscale
      implicit none
      logical fail3,tres,ehat,merge
      integer l,lmin,lmax,bastype,mlst
      real*8:: xl,j,jn,jcore,ex
      real*8:: bosc,gamma,acg,r1
      real*8:: eta
      integer:: ichsp,ic,iset,nset,nchsp,nfmax,nho,parity
      integer:: nk,nbins,inc
      integer,save:: prevset
      CHARACTER*1 BLANK,PSIGN(3)      
      character*40 filewf
      real*8 :: wcut(maxchan) !,vscale
      DATA PSIGN / '-','?','+' /, BLANK / ' ' /



      namelist/jpset/ jtot,parity,l,j,lmin,lmax,
     &                bastype,nfmax,exmin,exmax,
     &                nho,bosc,     ! HO
     &                gamma,mlst,   !THO Amos 
     &                eta,          !cTHO Lay
     &                nsp,  ! sp eigenvalues to keep for full diag
     &                bas2,
     &                nk,    ! nb of continuum wfs to construct each bin
     &                nbins, ! nb of bins
     &                inc,   ! incoming channel (bins)
     &                tres,  ! weight by T-matrix 
     &                ehat,  ! (logical, default T) to use mean bin energies (otherwise midpoint), 
     &                filewf, ! external file for wfs
     &                wcut,   ! mininum weight per channel to be retained (default 1) 
     &                vscale, ! scaling factor for v-core potential
     &                acg,r1  !CG	 
!     &                realcc ! if TRUE, calculate real multichannel states instead of scat. states


c Initialize variables and assign default values
      nfmax=0; exmin=-1e-20; exmax=1e20; inc=1; tres=.false.; nho=0;
      ehat=.false. ; filewf=""
      wcut(1:maxchan) = 0.0
      lmin=-1; lmax=-1; l=-1; j=-1;
      gamma=0d0
      eta=0d0
      vscale=1.
  	nho=0

      read(kin,nml=jpset) 
      nset=indjset(iset)
      if (l.lt.0) l=0;
      lmin=l;
      if (lmax.lt.lmin) lmax=lmin

c To be done!
c      if (iset.eq.1) then
c         prevset=1
c     else if(nset.eq.prevset) then
c          write(*,*)'(merged with previous set)'
c      else
c         prevset=nset
c      endif

      partot=parity !composite parity
      jpiset(nset)%nho     = nho
      jpiset(nset)%nsp     = nsp
      jpiset(nset)%partot  = parity
      jpiset(nset)%jtot    = jtot
      jpiset(nset)%exmin   = exmin
      jpiset(nset)%exmax   = exmax
      jpiset(nset)%inc     = inc
      jpiset(nset)%bastype = bastype
      jpiset(nset)%vscale  = vscale
 
!      write(0,*)'set',nset,'scaled by',vscale

      write(*,'(/)') 
      write(*,*)  '******************************************'
      write(*,220) nset, jtot,psign(parity+2)
220   format(5x," BASIS SET:",i3," WITH TOTAL J/pi=",f4.1,a1)
!      write(*,*)' lmax restricted to=',lmax
      write(*,*)  '******************************************'


c determine allowed single-particle configurations 
      call spconf(nset,nchsp,lmin,lmax,j)
      jpiset(nset)%nchsp =nchsp
      jpiset(nset)%lsp(maxchan)=0; 
      
      write(*,*)'  '
      write(*,*)'   CORE+VALENCE channels:' 
      nchan=0
!      write(99,*)'read_jpi: nce,nchsp=',nce,nchsp
      do ic=1,nce
!        ex    = exc1(ic)
!        jcore = jc(ic)
         ex    = qnc(ic)%exc    
         jcore= qnc(ic)%jc   
      do ichsp=1,nchsp
        l     = qspl(ichsp); xl=l
        jn    = qspj(ichsp)
!        write(*,*) 'chan=',ichsp,l,jn,jcore
        if ((-1)**l*parc(ic).ne.partot) cycle      
        if (fail3(jn,jcore,jtot))       cycle
        if (fail3(xl,jn,sn))            cycle
        nchan=nchan+1
        if (nchan.gt.maxchan) then 
          write(*,*)"Too many channels!"
          write(*,*)" Increase maxchan in modules.f90"
          stop
        endif

        jpiset(nset)%spindex(nchan)=ichsp
        jpiset(nset)%cindex (nchan)=ic
        jpiset(nset)%lsp(nchan)    =l
        jpiset(nset)%jsp(nchan)    =jn
        jpiset(nset)%jc(nchan)     =jcore
        jpiset(nset)%exc(nchan)    =ex
        jpiset(nset)%parc(nchan)   =parc(ic)
        jpiset(nset)%bas2          =bas2
        jpiset(nset)%lsp(nchan)   =l  
!         write(0,*)'jpiset',nset,l
        jpiset(nset)%jsp(nchan)   =jn
        jpiset(nset)%jc(nchan)    =jcore
        jpiset(nset)%wcut(1:nchan)=wcut(1:nchan)
!        if (bastype.ne.1) jpiset(nset)%nex =
!     &     jpiset(nset)%nex+nbins
         if (bastype.ne.1) jpiset(nset)%nex =nbins


c deprecated variables
        spindex(nchan)=ichsp  
        cindex(nchan)=ic
        ql(nchan)=l
        qj(nchan)=jn
        qjc(nchan)=jcore
        exc(nchan)=ex               !I redefined it to have the excitation energy for each nchan
        qpar(nchan)=parc(ic)

      write(*,240) nchan,jcore,l,sn,jn,ex,ichsp,ic
240	format(5x,'Channel #',i2,' : Ic=',f4.1,' (',i2,f4.1,')',f4.1
     &,' Ex=',f6.3,' MeV with sp config #',i1,' and core #',i1)
       enddo !ichsp
      enddo !ic

      if (nchan.eq.0) then 
          write(*,*)' No core+valece channels found!'
          stop
      endif

      jpiset(nset)%nchan=nchan
      if (nchan.gt.nchmax) nchmax=nchan
      write(*,*)'  '

      call genbasis(bastype,nfmax,nho,bosc,acg,r1,
     & gamma,mlst,nsp,bas2,lmax)

      end subroutine


c
c determine & store single-particle channels compatible with jtot & lmax
c
      subroutine spconf(iset,nchsp,lmin,lmax,j)
!      use wfs, only:  ql,qspl,qspj,spindex,jc,parc,spchan,nce
      use channels
      implicit none
      logical:: fail3
      real*8 :: xl,jnmin,jnmax,jcore,j,jn
      integer:: l,lmin,lmax,njn,ic,ijn,iset
      integer,intent(out):: nchsp
      write(*,*)'         '
      write(*,*)'SINGLE-PARTICLE configurations:'
      write(*,'(3x,"[ l val-core truncated at lmin,lmax=",i2,1x,i2,"]")'
     &  ) lmin,lmax
     
      nchsp=0
	  	   
      do l=lmin,lmax
       xl=l
       jnmin=abs(xl-sn)
       jnmax=xl+sn
       njn=nint(jnmax-jnmin)
       do ijn=0,njn
       jn=jnmin+ijn
       if ((j.ge.0).and.(jn.ne.j)) cycle
       do ic=1,nce
          jcore=jc(ic)
          if ((-1)**l*parc(ic).ne.partot) cycle      
          if(fail3(xl,jn,sn))cycle
          if(fail3(jn,jcore,jtot))cycle
          nchsp=nchsp+1
	      qspl(nchsp)=l
          qspj(nchsp)=jn

          spchan(iset,nchsp)%l  =l
          spchan(iset,nchsp)%sn =sn
          spchan(iset,nchsp)%j  =jn
          spchan(iset,nchsp)%jc =jcore

          write(*,235) nchsp,l,sn,jn
!          goto 236  ! I change this go to by a cycle
          exit
235	format(7x,'S.p. configuration #',i2,' (',i2,f4.1,')',f4.1)
       enddo !ic
236    enddo !jn
       enddo !l 
       write(*,*)'   '
       end subroutine




c -----------------------------------------------
c ** Choose basis and generate basis functions
c -----------------------------------------------
      subroutine genbasis(bastype,nfmax,nho,bosc,acg,r1,
     & gamma,mlst,nsp,bas2,lmax)
       use wfs,only : wftho,nr,rvec
       use globals
       use sistema
       use constants   
       implicit none
       integer l,mlst,nfmax,bastype,nho,bas2,lmax,nsp,n,i,unit_out
       real*8 ::bosc,gamma,eta,acg,r1 

 
c *** THO basis
       if (allocated(wftho)) deallocate(wftho)
       select case(bastype)
c      ----------------------------------------------------
       case(0) !HO 
         write(*,*)' ** RADIAL BASIS **'
		 write(*,*)'TESTIIIIIIIIIIIIIIIIIIIIIIING HO'
		 
		 !Vamos a crear un archivo externo para que se impriman las funciones de la Base !!!Añadido 15/05/24
		 
         write(unit_out,251) bosc,nho
251       format(" - HO basis with b=",1f5.2,2x,
     &           "and N=",i3," states")		
	 
		 open(unit_out, file='FuncionesBaseHO.dat', status='replace') !!!Añadido 15/05/24
		 

c        if (dr>eps) then
c           nr=(rmax-rmin)/dr+1                             !!!! we may have problems with definition of nr !!!!!
c        else
c           write(0,*)'dr too small!. Aborting'; stop
c        endif

       allocate(wftho(nho,0:lmax,nr))
       do l=0,lmax
       call hobasis(l,bosc,nho)
	     ! Imprime las funciones de la base en el archivo externo
		 !write(unit_out, '(A, 3X, A)') "R", "U(R)"
	     !write(unit_out, '(A)') "----------"
         do i = 1, nr
            do n = 1, nho
                ! Imprime r y u(r) en dos columnas separadas
                write(unit_out, '(2E12.6, A)') rvec(i), wftho(n, l, i)
			enddo
         enddo
	    enddo
		
		! Cierra el archivo
		close(unit_out)
	   
c -------------------------------------------------------------	   
	   ! Base Creada por Daniel Arjona:
  
       case(3) !Complex Gaussian Basis (CG)

         write(*,*)' ** COMPLEX GAUSSIAN BASIS **'
		 !write(*,*)'TESTIIIIIIIIIIIIIIIIIIIIIIING CGO - BASIS.F90'
		 
		 !Vamos a crear un archivo externo para que se impriman las funciones de la Base !!!Añadido 15/05/24
		 
		write(0,252) acg,r1,nho 
252       format(" - CG basis with acg,r1=",
     & 2f5.2,2x,"and N=",i3," states")

		 open(unit_out, file='FuncionesBaseCG.dat', status='replace') !!!Añadido 15/05/24
		 
         
       

       allocate(wftho(nho,0:lmax,nr))
       do l=0,lmax
       call cgbasis(l,acg,r1,nho)
	     ! Imprime las funciones de la base en el archivo externo
		 !write(unit_out, '(A, 3X, A)') "R", "U(R)"
	     !write(unit_out, '(A)') "----------"
		 do n = 1, nho
			do i = 1, nr
            !do n = 1, nho
                ! Imprime r y u(r) en dos columnas separadas
                write(unit_out,*) rvec(i), wftho(n, l, i)
			enddo
			write(unit_out,*)'&'
         enddo
	    enddo
		
		! Cierra el archivo
		close(unit_out)	   
   


c -------------------------------------------------------------
       case(1) ! THO 
        write(*,*)' ** RADIAL BASIS **'

        if (gamma<1e-6) then
!        gamma=bosc*(8*mu12*abs(egs)/hc/hc)**(1d0/4d0)
          write(*,*)'( gamma=0 ! => aborting )';stop
        endif
        write(*,250) mlst,bosc,gamma,nho
250     format(" - THO basis with m=",i2,2x,
     &           "b=",1f5.2,2x,"fm, gamma=",1f6.3,1x,
     &           "and N=",i3," states")

       if (nr.eq.0) then 
        write(*,*)'nr=0 in genbasis()!!'; stop
       endif

       allocate(wftho(nho,0:lmax,nr))
       do l=0,lmax
       call thobasis(l,mlst,bosc,gamma,nho)
       enddo

       if (bas2.eq.0) then
	  write(*,*)' [bas2=0 => Use Hsp eigenstates to diagonalize
     &full H]'
       else if  (bas2.eq.1) then
	  write(*,*)'  [bas2=1 => Use THO basis to diagonalizd full H]'
        else
       write(*,*)' bas2=',bas2,' not implemented!'
       endif

c      -----------------------------------------------------
       case(7)  ! cTHO
c      -----------------------------------------------------


        if (gamma<1e-6) then
!          write(*,*)'Im here'
         gamma=bosc*(8*mu12*abs(egs)/hc/hc)**(1d0/4d0)       !8???
         write(*,*)'( gamma=0 => will use Amos formula )'
         write(*,*)'gamma=',gamma
         endif
        if (abs(egs)<1e-6) then
!        gamma=bosc*(8*mu12*abs(egs)/hc/hc)**(1d0/4d0)
          write(*,*)'( egs=0 ! => aborting )';stop
        endif
        if (abs(eta)<1e-6) then
         eta=zv*zc*e2/dsqrt(2.*mu12*abs(egs)/hc/hc)*mu12/hc/hc
        endif
         write(*,*)'eta=',eta
       allocate(wftho(nho,0:lmax,nr))
       do l=0,lmax
         call cthobasis(l,mlst,bosc,gamma,eta,nho)	
       enddo

          write(0,52) mlst,bosc,gamma,eta,nho
52        format(" - Charged THO basis with m=",i2,2x,
     &           "b=",1f5.2,2x,"fm, gamma=",1f6.3,2x,
     &           " somm=",1f6.3," and N=",i3," states")

       if (bas2.eq.0) then
	  write(*,*)' bas2=0 => Use Hsp eigenstates to diagonalize full H'
       else if  (bas2.eq.1) then
	  write(*,*)' bas2=1 => Use THO basis to diagonalizd full H'
        else
       write(*,*)' bas2=',bas2,' not implemented!'
       endif


c -------------------------------------------------------------
c -------------------------------------------------------------
       case(2) ! CC bins
       write(*,*)' ** BIN wfs ** '

c -------------------------------------------------------------
       case(4) ! Real CC Bins
       write(*,*)' ** Real CC BIN wfs ** '

c -------------------------------------------------------------
       case(5) ! External 
       write(*,*)' ** External wfs ** '

c -------------------------------------------------------------
       case default ! undefined bastype
       write(*,*)' Bastype=',bastype,' not valid. Stopping!'
       stop
       end select 

      end subroutine



 	
