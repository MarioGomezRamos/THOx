c 
c *** Calculate valence/core configurations compatible with Jtot/pi
c     nchan={Ic,lsj} configurations  
      subroutine read_jpiset(iset,bastype,nk,tres,ehat,filewf)
      use channels
      use globals,only: kin
      use wfs    ,only: exmin,exmax
      implicit none
      logical fail3,tres,ehat
      integer l,integer,bastype,mlst
      real*8:: xl,jn,jcore,ex
      real*8:: bosc,gamma
      integer:: ichsp,ic,iset,nchsp,nfmax,nho,parity
      integer:: nk,nbins,inc
      CHARACTER*1 BLANK,PSIGN(3)      
      character*40 filewf
      DATA PSIGN / '-','?','+' /, BLANK / ' ' /



      namelist/jpset/ jtot,parity,lmax,
     &                bastype,nfmax,exmin,exmax,
     &                nho,bosc,     ! HO
     &                gamma,mlst,   !THO Amos 
     &                nsp,  ! sp eigenvalues to keep for full diag
     &                bas2,
     &                nk,    ! nb of continuum wfs to construct each bin
     &                nbins, ! nb of bins
     &                inc,   ! incoming channel (bins)
     &                tres,  ! weight by T-matrix 
     &                ehat,  ! (logical, default T) to use mean bin energies (otherwise midpoint), 
     &                filewf  ! external file for wfs
!     &                realcc ! if TRUE, calculate real multichannel states instead of scat. states


      nfmax=0; exmin=-1e-20; exmax=1e20; inc=1; tres=.false.; nho=0;
      ehat=.false. ; filewf=""
      read(kin,nml=jpset) 
      partot=parity !composite parity
      jpiset(iset)%nho     = nho
      jpiset(iset)%nsp     = nsp
      jpiset(iset)%partot  = parity
      jpiset(iset)%jtot    = jtot
      jpiset(iset)%exmin   = exmin
      jpiset(iset)%exmax   = exmax
      jpiset(iset)%inc     = inc
      jpiset(iset)%bastype = bastype

      write(*,'(/)') 
      write(*,*)  '******************************************'
      write(*,220) iset, jtot,psign(parity+2)
220   format(5x," BASIS SET:",i3," WITH TOTAL J/pi=",f4.1,a1)
!      write(*,*)' lmax restricted to=',lmax
      write(*,*)  '******************************************'


c determine allowed single-particle configurations 
      call spconf(iset,nchsp)
      jpiset(iset)%nchsp =nchsp
      
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

        jpiset(iset)%spindex(nchan)=ichsp
        jpiset(iset)%cindex (nchan)=ic
        jpiset(iset)%lsp(nchan)    =l
        jpiset(iset)%jsp(nchan)    =jn
        jpiset(iset)%jc(nchan)     =jcore
        jpiset(iset)%exc(nchan)    =ex
        jpiset(iset)%parc(nchan)   =parc(ic)
        jpiset(iset)%bas2          =bas2

        jpiset(iset)%lsp(nchan)   =l  
        jpiset(iset)%jsp(nchan)   =jn
        jpiset(iset)%jc(nchan)    =jcore
        if (bastype.ne.1) jpiset(iset)%nex =nbins

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

      jpiset(iset)%nchan=nchan
      if (nchan.gt.nchmax) nchmax=nchan
      write(*,*)'  '

      call genbasis(bastype,nfmax,nho,bosc,
     & gamma,mlst,nsp,bas2,lmax)

      end subroutine


c
c determine & store single-particle channels compatible with jtot & lmax
c
      subroutine spconf(iset,nchsp)
!      use wfs, only:  ql,qspl,qspj,spindex,jc,parc,spchan,nce
      use channels
      implicit none
      logical:: fail3
      real*8 :: xl,jnmin,jnmax,jcore,jn
      integer:: l,njn,ic,ijn,iset
      integer,intent(out):: nchsp
      write(*,*)'         '
      write(*,*)'SINGLE-PARTICLE configurations:'
      write(*,'(3x,"[ l val-core truncated at lmax=",i2,"]")') lmax
      nchsp=0
	  	   
      do l=0,lmax
       xl=l
       jnmin=abs(xl-sn)
       jnmax=xl+sn
       njn=nint(jnmax-jnmin)
       do ijn=0,njn
       jn=jnmin+ijn
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
      subroutine genbasis(bastype,nfmax,nho,bosc,
     & gamma,mlst,nsp,bas2,lmax)
       use wfs,only : wftho,nr
       use globals
       implicit none
       integer l,mlst,nfmax,bastype,nho,bas2,lmax,nsp
       real*8 ::bosc,gamma

 
c *** THO basis
!       read(kin,nml=basis) 
       if (allocated(wftho)) deallocate(wftho)
       select case(bastype)
c      ----------------------------------------------------
       case(0) !HO 
         write(*,*)' ** RADIAL BASIS **'
         write(0,251) bosc,nho
251       format(" - HO basis with b=",1f5.2,2x,
     &           "and N=",i3," states")
c        if (dr>eps) then
c           nr=(rmax-rmin)/dr+1                             !!!! we may have problems with definition of nr !!!!!
c        else
c           write(0,*)'dr too small!. Aborting'; stop
c        endif

       allocate(wftho(nho,0:lmax,nr))
       do l=0,lmax
       call hobasis(l,bosc,nho)
       enddo
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
	  write(*,*)' bas2=0 => Use Hsp eigenstates to diagonalize full H'
       else if  (bas2.eq.1) then
	  write(*,*)' bas2=1 => Use THO basis to diagonalizd full H'
        else
       write(*,*)' bas2=',bas2,' not implemented!'
       endif


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



 	
