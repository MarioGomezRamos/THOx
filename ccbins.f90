c     Construct continuum single or multi-channel bins
      subroutine makebins(iset,nk,tres,ifhat)
      use globals, only: mu12,egs,kin,written,verb
      use constants
      use sistema
      use wfs, only:nr,dr,rvec,wfc,energ,wbin,rint,rmin
      use channels, only: jpiset,jpsets,nchmax,ql
      use parameters, only: maxchan,maxeset
      use potentials, only: ccmat
      use xcdcc, only: binset
      implicit none
c     ----------------------------------------------
      logical :: energy,tres,debug,ifhat
      integer :: iset,inc,nbins,nchan,ib,ich,ik,ir,nk
      real*8  :: emin,emax,ei,ef,ebin,excore,ehat,khatsq
      real*8  :: r,bnorm,raux,r2aux,bcoef,chnorm(maxchan),rms
      real*8  :: ddk,ki,kf,kmin,kmax,kstep,kmid,k
      real*8  :: psh_el(nk),deltap, deladd,deltai
      real*8  :: fconv,yint

      complex*16, parameter::iu=(0.,1.)
      complex*16,allocatable,target:: wfcont(:,:,:)
      complex*16,pointer       :: wfbin(:)
      complex*16 :: smate(nk,maxchan),tfac,yfac
      CHARACTER(LEN=80) formato
c  bincc variables
      logical tdel,tkmat,tknrm,bincc
      integer:: nf   ! number of formfactors?
      integer:: isc,il,pcon,lmax
      integer:: maxn,nnk(maxchan)
      real*8:: anc,k2(maxchan),conv,bphase(nk)
      complex*16:: y(nr-1,maxchan),etap(maxchan)
      real*8:: fmscal
      complex*16:: ccoup(maxchan,maxchan,nr)
c  end bincc 
c *** ----------------------------------------------
      fconv=hc**2/2/mu12 ! E=fconv*k^2

c *** TEMP
      energy=.false.
      debug=.false.
      bincc=.false.
c *** 

      emin  =jpiset(iset)%exmin
      emax  =jpiset(iset)%exmax
      inc   =jpiset(iset)%inc
      nbins =jpiset(iset)%nex
      nchan =jpiset(iset)%nchan
      excore=jpiset(iset)%exc(inc)
      
      if (inc.gt.nchan) then
        write(*,*)'inc=',inc,'but only',nchan,' in this j/pi set'
      endif

      kmax=sqrt(2*mu12*emax)/hc
      if (emin.lt.0) emin=0
      kmin=sqrt(2*mu12*emin)/hc
      if (nbins.gt.0) then
       kstep=(kmax-kmin)/nbins
      else
       write(*,*)'makebins: nbins=0!'; stop
      endif

c *** Construct and store bin wfs
c *** We need to allocate this only for the first j/pi set
      if (.not.allocated(wfc)) then
        nchmax=0
        allocate(wfc(jpsets,maxeset,maxchan,nr))
        wfc(:,:,:,:)=(0.,0.); 
      endif
      if (.not.allocated(energ)) then
           allocate(energ(jpsets,maxeset))
           energ(:,:)=0.
      endif
      if (.not.allocated(wbin)) then 
           allocate(wbin(jpsets,maxeset))
           allocate(binset(jpsets)) 
           wbin(:,:)=0.
      endif
      
      binset(iset)%emin=emin
      binset(iset)%emax=emax
      binset(iset)%kmin=kmin
      binset(iset)%kmax=kmax
      binset(iset)%nbins=nbins
      binset(iset)%ifhat=ifhat
      deladd=0.

      do ib=1,nbins
      tfac=1.
      chnorm=0.
      bnorm=0.
      ki=kmin+(ib-1.)*kstep
      kf=ki+kstep
!      bcoef=sqrt(2./pi/(kf-ki))
      bcoef=sqrt(2./pi)
      ddk=(kf-ki)/(nk-1)   !!!! CHECK!!!!!!!!!!!!
      kmid=(kf+ki)/2.
      ei=fconv*ki**2
      ef=fconv*kf**2
!!      ebin=fconv*kmid**2  
      wbin(iset,ib) =ef-ei     
c < phi |H | phi> 
      khatsq=(kf-ki)**2/12.d0
      khatsq=khatsq+(kf+ki)**2/4.d0
      ehat=fconv*khatsq
!!! !!!!!CHECK BEST DEFINITION
      if (ifhat) then
         ebin=ehat
      else
         ebin=(ei+ef)/2.
      endif

      energ(iset,ib)=ebin + excore  !!! CHECK!!!!       

      binset(iset)%khat(ib)=sqrt(khatsq)
      binset(iset)%klow(ib)=ki
      binset(iset)%kup(ib) =kf
      binset(iset)%kmid(ib)=kmid
      binset(iset)%ebin(ib)=ebin
      binset(iset)%nk(ib)  =nk
!      write(*,*)'ib=',ib,'ei-ef=',ei,ef,' tres=',tres
         
      if (ei.lt.1e-3) ei=1e-3
      if (allocated(wfcont)) deallocate(wfcont)
      allocate(wfcont(nk,nchan,nr)); wfcont=0.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (.NOT.BINCC) THEN 
      call wfrange(iset,nchan,inc,ei,ef,nk,energy,wfcont,smate,psh_el)

c make phase-shifts continuous (not really needed for different bins!)
      deltai=psh_el(1)
      if (ib.eq.1) then
          deltap=deltai
      else
         if(deltai<deltap-90) deladd=deladd+180
         if(deltai>deltap+90) deladd=deladd-180
      endif
      deltap=psh_el(nk)
      psh_el(:)=psh_el(:) + deladd
     
!!!!!!!!!!!!!!!!!!!!!!!! !TEST
!      if (ib.eq.1) then
!        do ik=25,25 ! nk
!        do ir=1,nr
!        r=rvec(ir)
!        write(97,'(1f8.3,3x,20g14.6)')r,
!     &  (exp(-iu*psh_el(ik)*pi/180)*wfcont(ik,ich,ir),ich=1,nchan)
!        enddo ! ir 
!        write(97,*)'&'
!        enddo ! ik
!      endif
!!!!!!!!!!!!!!!!!!!!!!!! TEST

c k-average
      raux=0 ; r2aux=0; chnorm(:)=0; bnorm=0
      do ich=1,nchan
      nullify(wfbin)
      wfbin=> wfc(iset,ib,ich,:)
      
      wfbin=0
      yint=0.
      do ik=1,nk
      k=ki+(ik-1)*ddk
      yfac=exp(-iu*psh_el(ik)*pi/180.)
      tfac=(smate(ik,inc)-1.)/(0.,2.)
      if (tres) yfac=conjg(tfac)
      binset(iset)%pshel(ib,ik)=psh_el(ik)*pi/180.
!!! TEST
!      if (ich.eq.1)then 
!      write(460,'(1f6.3,2x,10f12.5)') fconv*k**2,tfac    
      write(450,'(1f6.3,2x,10f12.5)') fconv*k**2,psh_el(ik)*pi/180
      if (ik.eq.nk) write(450,*)'&'
!      endif 
!!! TEST
      wfbin(:)= wfbin(:) + bcoef*yfac*wfcont(ik,ich,:)*ddk
      YINT = YINT + ABS(YFAC)**2 * ddK
      enddo ! ik
      wfbin(:)=wfbin(:)/sqrt(yint)

c compute channel normalization
      do ir=1,nr
      r=rvec(ir)
      raux=raux   + abs(wfbin(ir))**2*dr
      r2aux=r2aux + (abs(wfbin(ir))*r)**2*dr
      chnorm(ich)=chnorm(ich)+abs(wfbin(ir))**2*dr
      enddo !ir
      bnorm=bnorm + chnorm(ich)
!      write(*,*)'Bin:',ib,' Chan=',ich,' Norm=',chnorm(ich)
c     ----------
      enddo ! ich (channels)
      rms=sqrt(r2aux/bnorm)
!      write(formato,'( "(5x,",'' "#" '',", 3i3,2x, ",'' "Ecm=[", '',
!     & "2f8.4,",'' "Norm=" '',
!     & ", 1f8.4,",'' "[", '',I0,"f8.4",'' ,"]" '',",4x,",    
!     & '' "rms=" '',",1f8.4)"
!     & )' ) NCHAN


!       write(formato,'( "(5x,",'' "#" '',",i3,1f10.4,2x,",
!     & '' "[", '',I0,"f8.4",'' ,"]" '',",4x,",    
!     & '' "rms=" '',",1f8.4)"
!     & )' ) NCHAN
!      write(*,*)formato
!      write(*,formato) ib, inc,nchan,ei,ef,bnorm,chnorm(1:nchan),rms

      ELSE ! ...... TEST WITH BINCC FROM FRESCO ......................
      nf=1;  isc=2; il=inc ; pcon=2
      tres=.false.; tdel=.true.; tkmat=.false.; tknrm=.false.
      y(:,:)=0.0
      if ((rint.gt.0).and.(rint.lt.rvec(nr))) then
          maxn=ceiling((rint-0.0)/dr)+1 ! bincc assumes rmin=0!
      else
          maxn=nr
      endif

      call coefmat(nchan)
      ccmat(:,:,:)=-ccmat(:,:,:)  ! needed for bincc2 convention
      do ich=1,nchan
       ql(ich)    =jpiset(iset)%lsp(ich)  
       k2(ich)=abs(egs)+jpiset(iset)%exc(ich) ! separation energy for each channel
       conv=1/fconv     
       k2(ich)=k2(ich)*conv     
       etap(ich)= 2*e2*zv*zc*mu12/hc**2  !=2 k eta 
      enddo      
      lmax=maxval(ql)
!      write(*,*)'bincc: nchan,inc,lmax,conv=',nchan,inc,lmax,conv
      call BINCC2(Y,ccmat,nf,NCHAN,k2,il,conv,bphase,isc,ki,kf,nk,
     & etap,maxn-2,dr,pcon,tres,tdel,tkmat,tknrm,lmax+1,maxn,anc)
 
      do ik=1,nk
      k=ki+(ik-1)*ddk
      write(450,'(1f6.3,2x,10f12.5)')fconv*k**2,bphase(ik)*180./pi
      enddo

      do ich=1,nchan
      wfc(iset,ib,ich,1:maxn-1)=y(1:maxn-1,ich)
      enddo
      deallocate(ccmat)
      ENDIF 
!.........................................................................


      write(*,100) ib,ebin,ei,ef,nchan,inc,bnorm,rms
100   format(5x,"Bin #",i3,":", 3x,
     & "Emid=",1f8.3,2x,"Erel=[",1f8.4,1x,"-",1f8.4,"] MeV;",
     & i2," chan(s) (Inc=",i3,")",3x,
     & "Norm=",1f8.4,3x,"Rms=",1f8.4)   


c  write bin wfs
      if (debug)
     &write(98,101) ib,ebin,ei,ef,nchan,inc,bnorm,rms
101   format(5x,"#Bin ",i3,":", 3x,
     & "Emid=",1f8.3,2x,"Erel=[",1f8.4,1x,"-",1f8.4,"] MeV;",
     & i2," chan(s) (Inc=",i3,")",3x,
     & "Norm=",1f8.4,3x,"Rms=",1f8.4) 
      do ir=1,nr
        r=rvec(ir)
        if (debug)write(98,'(1f8.3,3x,50g14.6)')r,
     &  (wfc(iset,ib,ich,ir),ich=1,nchan)
      enddo ! ir
      if (debug)write(98,*)'&'

      if (bnorm.gt.1.2) then 
        write(*,*)'** ERROR ** Bin',ib,' gave Norm=',bnorm
        stop
      endif


c  u(r) -> R(r) for compatibility with the rest of the program
      do ir=1,nr
      r=rvec(ir)
      if (r.lt.1e-5) r=1e-5
      wfc(iset,ib,:,ir)=wfc(iset,ib,:,ir)/r
      enddo !ir

      enddo ! ib
      deallocate(wfcont)
      end subroutine
