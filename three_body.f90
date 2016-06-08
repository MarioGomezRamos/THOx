*     ---------------------------------------------------------------      
      subroutine three_body
      use globals, only: egs
      use scattering
      use wfs
      use channels
      use parameters, only: maxchan
      implicit real*8 (a-h,k-m,o-z)
*     --------------------------------------------------------------- 
*     JAT October 1999 (FN+IJT CDCC for three-body observables)
*     Modified February 2004 for use of closed channels
*     Minor modifications and checks by JAT, ANU 2008.
*     Mods/corrections for kinematic roots criteria - April 2010 (IJT)
*     --------------------------------------------------------------- 
*     max number of: bins, k-steps/bin, cm angles, total 
*     excited states, different stacks, number of bins per stack.
*     ---------------------------------------------------------------
      parameter (ibg=150,nkm=201,nanm=1009,nexg=150,nst=55,nbst=32)
      parameter (ncleb=2000,immmp=80)
      parameter (ncemax=5,nchmaxx=5,nchmaxx2=25,nerel=100)
*     ---------------------------------------------------------------
      real*8 binls(ibg),binjs(ibg),emids(ibg),kminJs(ibg),kmaxJs(ibg)
      real*8 phases(ibg,nkm),pars(ibg),bindk(ibg)
      real*8 khat(ibg),ehat(ibg),exen(ibg)
      real*8 exjs(nexg),par(nexg),eex(nexg),angsd(nanm),angsr(nanm)
      real*8 chjs(nst),chls(nst),chpar(nst),chjis(nst),chpari(nst)
      real*8 kminJ(nst),kmaxJ(nst),xyt(2,nbst+1),Ko(nst,nbst+1)
      real*8 cgc(ncleb),cplm(0:10,21),sigmad(0:10),sigphi(4000)
      real*8 xIce(ncemax),xece(ncemax),chann(3,nchmaxx)
      real*8, allocatable:: dsdew(:,:)
*     ---------------------------------------------------------------
*     real*8 sigcs(nst,nbst,nanm),qbig(3)
*     ---------------------------------------------------------------
*     THREE COMPONENT ARRAYS TO STORE VECTORS
*     z-axis in incident beam direction
*     ---------------------------------------------------------------    
      real*8 qp(3),kcL(3),kvL(3),pcL(3),pvL(3),tem(3)
      real*8 ptot(3),ktot(3),bkp(3),kp(3),qv(3),qc(3)
*     --------------------------------------------------------------- 
*     last index in fam (immmp) is max(2*sp+1)*(2*J*+1) dimension
*     --------------------------------------------------------------- 
      complex*16 ylm(0:10,21),phas,phask,phask1,amp,fin(nst,immmp)  
      complex*16 phiph(-21:21),pph,f2c,phbin(nst),phasil(0:10)
      complex*16 fam(nexg,nanm,immmp),i,fxyc(10,nbst+1)
      complex*16 xscat(nchmaxx2),fauxi(nr),caux,zero,resc,suma
      complex*16,allocatable:: wf(:,:),sumi(:)
      complex*16, pointer:: over
*     --------------------------------------------------------------- 
      integer iiex(nst,nbst),iibin(nst,nbst),iiist(nst),ibused(nst)
      integer nks(ibg),ibin(nexg),iich(nst),ifinc(nst),isc(ibg)
      integer ipar(ncemax)
      integer, allocatable:: ips(:,:)
*     --------------------------------------------------------------- 
      character ijtfile*30,ansa*12,outfile*20,sys*3,head1*40,head2*40
      logical zert,tes1,tes2
*     ---------------------------------------------------------------
      real*8 faclog(500)
      common/clebma/faclog
      call factor
*     --------------------------------------------------------------- 
*     Factorials print diagnostic
*     --------------------------------------------------------------- 
*     do ii=0,5
*      print*,ii,exp(faclog(ii+1))
*     enddo
*     ---------------------------------------------------------------
*     constants (IJT/FRESCO values)
*     ---------------------------------------------------------------
      acc=1.d-6
      i=(0.d0,1.d0)
      pi=4.d0*datan(1.d0)
      degrad=pi/180.d0
      zero=(0d0,0d0)
*     ---------------------------------------------------------------
*     FRESCO use may involve unitmass parameter - needs to be read
*     ---------------------------------------------------------------
      read*,idcore
      read*,unitmass
      hbarc=197.32705d0
      fsc=1.d0/137.03599d0
      amu=931.49432d0
      amu=unitmass*amu
      hbarc2=hbarc*hbarc
      h2m=hbarc2/(2.d0*amu)
      coul=hbarc*fsc
      print 400
      print*,' unitmass used is ',unitmass
      print*,' hbar**2/2*amu value ',real(h2m)
      print*,' Coulomb parameter is',real(coul)
      print 400
      read*,nord
      print*,' order of 2D-Aitken interpolation is ',nord
*     ---------------------------------------------------------------
      print 400
      print*,' output file for cross sections '
      read'(a)',outfile
      print'(a,a)','  output will be written to file: ',outfile
      print 400
      open(77,file=outfile,status='unknown')
*     --------------------------------------------------------------- 
*     frame for specification of detected energies and observables 
      print*,' lab or com system '
      read'(a)',sys
      write(77,'(a)') sys
      print*,' coordinate system is set to: ',sys
      print 400
*     --------------------------------------------------------------- 
      print*,' target mass and charge and incident energy/nucleon '
      read*,mt,zt,eoa
      print 316,mt,zt,eoa
      write(77,316) mt,zt,eoa
*     projectile constituents and binding energy
      print*,' core particle mass and charge'
      read*,mc,zc
      print 316,mc,zc
      write(77,316) mc,zc
      print*,' valence particle mass and sep energy (>0) '
      read*,mv,ebind
      print 316,mv,ebind
      write(77,316) mv,ebind
*     ---------------------------------------------------------------
*     projectile incident lab energy (MeV)
      ep=eoa*(mc+mv)
      write(77,316) ep
      print 303,'  projectile of mass ',mc+mv
      print 303,'  projectile lab energy = ',ep
      print 400  
*     ---------------------------------------------------------------
      print*,' detected energy is for: core(1) or valence(2) '
      read*,idet
      write(77,*) idet
      if(idet.eq.1) then
       print*,' idet=1: core particle energy is specified '
       else
       print*,' idet=2: valence particle energy specified ' 
      endif
      print 400  
*     ---------------------------------------------------------------
      print*,' detected energy range: lower:upper:step '
      read*,Enlow,Enup,dEn
      iten=nint((Enup-Enlow)/dEn)+1
*     adjust dEn in case of partial number of steps
      dEn=(Enup-Enlow)/(iten-1)
      print 316,Enlow,Enup,dEn
*     ---------------------------------------------------------------
*     particle with specified energy defines phi=0 degrees
*     ---------------------------------------------------------------
      if(idet.eq.1) then
       cospc=1.d0
       sinpc=0.d0
      else
       cospv=1.d0
       sinpv=0.d0
      endif
*     ---------------------------------------------------------------           
      print*,' core theta range: lower upper and step (degrees)'
      read*,tcl,tcu,dtc
      itc=nint((tcu-tcl)/dtc)+1
*     adjust dtc in case of partial number of steps
      dtc=(tcu-tcl)/(itc-1)
      print 316,tcl,tcu,dtc
      dtcr=dtc*degrad
      write(77,307) itc,dtcr,tcl,tcu,dtc
*     ---------------------------------------------------------------           
      print*,' valence theta range: lower upper and step (degrees)'
      read*,tvl,tvu,dtv
      itv=nint((tvu-tvl)/dtv)+1
*     adjust dtv in case of partial number of steps
      dtv=(tvu-tvl)/(itv-1)
      print 316,tvl,tvu,dtv
      dtvr=dtv*degrad
      write(77,307) itv,dtvr,tvl,tvu,dtv
*     ---------------------------------------------------------------           
      write(77,307) iten,dEn,Enlow,Enup 
*     ---------------------------------------------------------------           
 112  print*,' phi range: lower upper and step (degrees)         '
      print*,' ================================================= '
      print*,' if phi(upper) > phi(lower) then these are usually '
      print*,' expected to be 360 and 0 with a specified step if '
      print*,' to be used correctly later to integrate over one  '
      print*,' or more particle solid angles                     '
      print*,' ================================================= '
      read*,phil,phiu,dphi
      itphi=nint((phiu-phil)/dphi)+1
*     adjust dphi in case of partial number of steps
      dphi=(phiu-phil)/(itphi-1)
      print 316,phil,phiu,dphi
*     ---------------------------------------------------------------
      if((-1)**nint((phiu-phil)/(2.d0*dphi)).lt.0) then
       print*,' not even number of intervals (0-180) so problem if'
       print*,' plan to integrate (uses Simpson): reenter (y/n?) '
       read'(a)',ansa
       if(ansa.eq.'y') goto 112
      endif
      if(itphi.gt.1) print*,' phi values will be integrated over '
      dphir=dphi*degrad
      write(77,307) itphi,dphir,phil,phiu,dphi
      print 400
*     ---------------------------------------------------------------
*     here start the modifications to read CDCC breakup amplitudes
*     from FRESCO. Detailed conventions specified elsewhere.
*     --------------------------------------------------------------- 
      print*,'  input file for CDCC amplitudes '
      read'(a)',ijtfile
      print'(a,a)','   CDCC input will be read from file: ',ijtfile
      open(88,file=ijtfile,status='unknown')
*     ---------------------------------------------------------------
*     strictly should compare inputs with those read already but
*     user is assumed to be alive to this risk!!
*     --------------------------------------------------------------- 
      read(88,'(a,a)') head1,head2
      read(88,*) ep,ebind
      read(88,'(4f8.4)') mp,mt,mc,mv
      read(88,'(4f8.4)') zp,zt,zc,zv
      read(88,'(a)') ansa
*     --------------------------------------------------------------- 
*     spins sx: st is assumed zero; compute statistical factors issx
*     --------------------------------------------------------------- 
      read(88,'(4f8.4)') sp,st,sc,sv
      read(88,*) xxp,xxt,xxc,xxv
      ipar(1)=xxc
      xIce(1)=sc
      xece(1)=0.d0
      issp=nint(2*sp+1)
      issc=nint(2*sc+1)
      issv=nint(2*sv+1)
*     --------------------------------------------------------------- 
      print 400
      print '(3x,a,/,3x,a)',head1,head2
      print 400
      print*,'   mass:      charge:   spin:     parity ---- target'
      print 304,mt,zt,st,xxt
      print*,'   mass:      charge:   spin:     parity ---- core'
      print 304,mc,zc,sc,xxc
      print*,'   mass:      charge:   spin:     parity ---- valence'
      print 304,mv,zv,sv,xxv
      print*,'   mass:      charge:   spin:     parity ---- proj'
      print 304,mp,zp,sp,xxp
*     --------------------------------------------------------------- 
      print*,'  lab energy and separation energy (>0) '
      print 304,ep, ebind
      print 400
*     ebind is now defined negative
      ebind=-ebind
*     --------------------------------------------------------------- 
*     read back bin and angles information: put angles in arrays 
*     --------------------------------------------------------------- 
      read(88,*) nbins,nkmaxJ,nexb
      read(88,*) nangl,thmin,thinc
      read(88,*) nce
      do ii=1,nce
      read(88,*) iparce,xjce,ece
      ipar(ii+1)=iparce
      xIce(ii+1)=xjce
      xece(ii+1)=ece
      enddo
*     --------------------------------------------------------------- 
      print*,'  There are ',nangl,' angles in steps of',thinc
      print 400
      if(nangl.gt.nanm) then
       print*,'  nanm is ',nanm,' so too small '
       print 400
       stop 
      endif
*     --------------------------------------------------------------- 
      do iang=1,nangl
       angsd(iang)=thmin+(iang-1)*thinc
       angsr(iang)=angsd(iang)*degrad
      enddo
*     --------------------------------------------------------------- 
      print*,'  There are ',nexb,' excited states '
      print 400
      if(nexb.gt.nexg) then
       print*,'  nexg is ',nexg,' so too small '
       print 400
       stop 
      endif
*     --------------------------------------------------------------- 
      print*,'  There are ',nbins,' bins '
      if(nbins.gt.ibg) then
       print*,'  ibg is ',ibg,' so too small '
       print 400
       stop 
      endif
*     --------------------------------------------------------------- 
      if(nkmaxJ.gt.nkm) then
       print*,'  nkm is ',nkm,' so too small '
       print 400
       stop 
      endif
*     --------------------------------------------------------------- 
*     loop over all bins reading limits, energies and nuclear phases 
*     Also compute bindk(i)=1/sqrt(deltak) factors for later scaling
*     --------------------------------------------------------------- 
ccc      print*,'   i  l    j     emid        kminJ        kmaxJ      nks'
ccc     +      ,'  isc'
ccc      do ibi=1,nbins
ccc       read(88,*) binls(ibi),binjs(ibi),emids(ibi),kminJs(ibi),
ccc     + kmaxJs(ibi),nks(ibi),nko,isc(ibi)
ccc       print 309, ibi,binls(ibi),binjs(ibi),emids(ibi),kminJs(ibi),
ccc     + kmaxJs(ibi),nks(ibi),isc(ibi)
ccc       read(88,305) (phases(ibi,ik), ik=1,nks(ibi))
ccc       bindk(ibi)=1.d0/sqrt(kmaxJs(ibi)-kminJs(ibi))
ccc       pars(ibi)=(-1.d0)**nint(binls(ibi))
ccc      enddo
      print*,'   i  l    j     emid        kminJ        kmaxJ      nks'
     +      ,'  isc'
      do ibi=1,nbins
       read(88,*) xjex,npex,nchann,nchu,emids(ibi),kminJs(ibi),
     + kmaxJs(ibi),nks(ibi),nko,isc(ibi)
       read(88,*) ((chann(ind1,ind2),ind1=1,3),ind2=1,nchann)
       binls(ibi)=chann(1,nchu)
       binjs(ibi)=chann(2,nchu)
       ice=chann(3,nchu)
       nch2=nchann*nchann
       print 309, ibi,binls(ibi),binjs(ibi),emids(ibi),kminJs(ibi),
     + kmaxJs(ibi),nks(ibi),isc(ibi)
       nksl=nks(ibi)
       do ik=1,nksl
       read(88,'(2f10.6)') phases(ibi,ik),xk
       read(88,'(10f25.6)') (xscat(i1),i1=1,nch2)
       enddo
       bindk(ibi)=1.d0/sqrt(kmaxJs(ibi)-kminJs(ibi))
       pars(ibi)=(-1.d0)**nint(binls(ibi))
      enddo
      print 400
*     --------------------------------------------------------------- 
      mucv=mc*mv/(mc+mv)*amu         
      print*,'   i  l    j     emid        ehat        exen '
      do ibi=1,nbins
       khatsq=(kmaxJs(ibi)-kminJs(ibi))**2/12.d0
       khatsq=khatsq+(kmaxJs(ibi)+kminJs(ibi))**2/4.d0
       khat(ibi)=sqrt(khatsq)
       if(isc(ibi).eq.2.or.isc(ibi).eq.0) then
       ehat(ibi)=hbarc2*khatsq/(2.d0*mucv)
       else if(isc(ibi).eq.12) then
        kvalsq=(kmaxJs(ibi)**5-kminJs(ibi)**5)/5.d0
        kvalsq=kvalsq/(kmaxJs(ibi)-kminJs(ibi))/khatsq
        ehat(ibi)=hbarc2*kvalsq/(2.d0*mucv)
       endif
       exen(ibi)=ehat(ibi)-ebind
       print 309, ibi,binls(ibi),binjs(ibi),emids(ibi),ehat(ibi),
     + exen(ibi)
      enddo
      print 400
*     --------------------------------------------------------------- 
*     loop over excited states, reading properties and bin labels
*     --------------------------------------------------------------- 
      print*,'  There are originally ',nexb,' excited states '
      print*,'    i  #m      parity       eex   ibin     Deltak'
      iex=1
      nbnd=0
      nclo=0
      do iix=1,nexb
       read(88,*,end=555) exjs(iex),par(iex),eex(iex),ibin(iex)
       iii=issp*nint(2*exjs(iex)+1)
*     --------------------------------------------------------------- 
       if(iex.gt.1.and.ibin(iex).gt.0) then
        ijump=ibin(iex)-(ibin(iex-1)+1)
        if(ijump.gt.0) then
         print*,'  there were',ijump,' closed channels'
         nclo=nclo+ijump
        endif
       endif
*     --------------------------------------------------------------
       if(ibin(iex).gt.0) then
        Deltak=kmaxJs(ibin(iex))-kminJs(ibin(iex))
        print 317,iex,iii,par(iex),eex(iex),ibin(iex),Deltak
       else
        print*,'  bound state with iex= ',iex,' discarded '
        nbnd=nbnd+1
       endif
*     --------------------------------------------------------------- 
       if(iii.gt.immmp) then
        print*,'  immmp is ',immmp,' so too small '
        print 400
        stop 
       endif
*     --------------------------------------------------------------- 
*      loop over all centre of mass angles reading amplitudes
*     --------------------------------------------------------------- 
       do iang=1,nangl
        read(88,306) (fam(iex,iang,iam),iam=1,iii)
*     --------------------------------------------------------------- 
*       diagnostic print if needed
*     --------------------------------------------------------------- 
*       if(iang.eq.1) print 308, (fam(iex,iang,iam),iam=1,iii)
       enddo
       if(ibin(iex).gt.0) then
        iex=iex+1
       endif
      enddo
  555 continue
      numread=iex-1
      ijump=nexb-numread-nbnd-nclo
      if(ijump.gt.0) then
       print*,'  there were',ijump,' closed channels'
       nclo=nclo+ijump
      endif
      print 400
      if(nexb.ne.numread) then
       print*,'  there were',nbnd,' bound states  '
       print*,'  there were',nclo,' closed states '
       print*,'  so now ',numread,' excited states '
       nexb = numread
      endif
      close(88)
      print 400
      nsetmax=maxval(jpiset(:)%nex)
      allocate(overlapmat(hdimtot,hdimtot,nchmaxx))
      allocate(oversum(hdimtot,hdimtot))
      allocate(gsolap(jpsets,nsetmax,maxchan,nq))
      overlapmat(:,:,:)=zero
      oversum=zero
*     --------------------------------------------------------------- 
*     print*,'  breakup amplitude m,mp elements counter check  '
*     spp=2
*     iam=0
*     do im=1,issp
*      mm=im-sp-1
*      do imp=1,nint(2*spp+1)
*       iam=iam+1
*       mmp=imp-spp-1
*       print'(4i6)',iam,nint((mm+sp)*(2*spp+1)+mmp+spp+1),nint(mm),
*    +  nint(mmp)
*      enddo
*     enddo
*     --------------------------------------------------------------- 
*     if charged core and valence particles, add Coulomb phases
*     this modified to be correct in the case of closed channels
*     loops on retained nexb excited states rather than bin index
*     --------------------------------------------------------------- 
ccc      if(ie(zv*zc,0.d0).gt.0) then
ccc       print*,'  Calculating Coulomb phases  '
ccc       etacon=coul*zv*zc*mucv/hbarc2
ccc       do iex=1,nexb
ccc        ibi=ibin(iex)
ccc        ilval=nint(binls(ibi))
ccc        dkb=(kmaxJs(ibi)-kminJs(ibi))/(nks(ibi)-1)
ccc        do ik=1,nks(ibi)
ccc         ak=kminJs(ibi)+(ik-1)*dkb
ccc         eta=etacon/ak
ccc         call sigma(ilval,eta,sigmad)
ccc         phases(ibi,ik)=phases(ibi,ik)+sigmad(ilval)
ccc        enddo
ccc       enddo
ccc       print 400
ccc      endif
*     --------------------------------------------------------------- 
*     loop again over excited states sorting and counting
*     initialise counters for stacks and channels in each stack
*     --------------------------------------------------------------- 
      ist=1
      ich=0
      ajval=exjs(1)
      alval=binls(ibin(1))
      ajival=binjs(ibin(1))
      apar=par(1)
*     --------------------------------------------------------------- 
*     loop over all excited states (ie(x,y) checks real equality)
*     --------------------------------------------------------------- 
      do iex=1,nexb
       ibi=ibin(iex)
       tes1=(ie(exjs(iex),ajval).eq.0.and.ie(par(iex),apar).eq.0)
       tes2=(ie(binjs(ibi),ajival).eq.0.and.ie(binls(ibi),alval).eq.0)
       if(tes1.and.tes2) then
        ich=ich+1
       else
        ist=ist+1
        ich=1
        ajval=exjs(iex)
        alval=binls(ibi)
        ajival=binjs(ibi)
        apar=par(iex)
       endif
       iiex(ist,ich)=iex
       iibin(ist,ich)=ibi
       iich(ist)=ich
*     --------------------------------------------------------------- 
       if(ich.gt.nbst) then
	print*,'  ich has reached ',ich
        print*,'  nbst is ',nbst,' so too small '
        print 400
        stop 
       endif
*     --------------------------------------------------------------- 
       chjs(ist)=exjs(iex)
       chls(ist)=binls(ibi)
       chjis(ist)=binjs(ibi)
       chpar(ist)=par(iex)
       chpari(ist)=pars(ibi)
*     --------------------------------------------------------------- 
*      for lowest bin take care of the bin boundaries
*     --------------------------------------------------------------- 
       if(ich.eq.1) then
        if(kminJs(ibi).lt.0.03d0) kminJs(ibi)=0.d0
        kmaxJs(ibi)=kminJs(ibin(iex+1))
        kminJ(ist)=kminJs(ibi)
       endif
       kmaxJ(ist)=kmaxJs(ibi)
      enddo
      istmax=ist
*     --------------------------------------------------------------- 
      print*,' so',nexb,' excited states in',istmax,' stacks'
*     --------------------------------------------------------------- 
      if(istmax.gt.nst) then
       print*,'  nst is ',nst,' so too small '
       print 400
       stop 
      endif
*     --------------------------------------------------------------- 
*     loop over stacks with diagnostic printouts of included bins
*     --------------------------------------------------------------- 
      do ist=1,istmax
       ifinc(ist)=1
       print*,' -------- stack',ist,' -------------------------'
       print*,'     J       parity      l         j       '
       print 304,chjs(ist),chpari(ist),chls(ist),chjis(ist)
       print 311,'   k-range',kminJ(ist),kmaxJ(ist),'  in',iich(ist),
     + ' bins'
*     --------------------------------------------------------------- 
*     loop over bins in stacks
*     --------------------------------------------------------------- 
       do ich=1,iich(ist)
        iii=iibin(ist,ich)
        print 310,'   ich=',ich,' ibin=',iii,kminJs(iii),kmaxJs(iii),
     +  nint(binls(iii)),iiex(ist,ich)
       enddo
      enddo
      print 400
*     --------------------------------------------------------------- 
*     eikonal Ylm constraint
*     --------------------------------------------------------------- 
      ifeik=0
*     read*,ifeik
*     if(ifeik.eq.1) then
*       print*,'  eikonal Ylm combinations '
*       print 400
*     endif
*     --------------------------------------------------------------- 
*     include all stacks in calculation or not
*     --------------------------------------------------------------- 
      read*,isubset
      if(isubset.gt.0) then
       read*,(ifinc(ist),ist=1,istmax)
       print*,'  stacks included in calculations '
       do ist=1,istmax
        print*,'  ist=',ist,'  ifinc(ist)=',ifinc(ist)
       enddo
       print 400
      else
       read*,igugg
      endif
*     --------------------------------------------------------------- 
*     relative energy cutuff
*     --------------------------------------------------------------- 
      ifuecut=0
*     print *,'     relative energy cutuff'
*     read*,ifuecut,uecut
      if(ifuecut.eq.0) uecut=1.d6
*     if(ifuecut.eq.1) then
*       print*,'  energy cutoff is ',uecut
*       print 400
*     endif
*     --------------------------------------------------------------- 
*     calculate differential cross sections for each bin state (mb/sr)
*     this with IJT original normalisations
*     --------------------------------------------------------------- 
*     do ist=1,istmax
*      spp=chjs(ist)
*      do ich=1,iich(ist)
*       iix=iiex(ist,ich)
*       do iang=1,nangl
*        sum=0.d0
*        do im=1,issp
*         do imp=1,nint(2*spp+1)
*          iii=nint((im-1)*(2*spp+1)+imp) 
*          sum=sum+abs(fam(iix,iang,iii))**2
*         enddo
*        enddo
*        sigcs(ist,ich,iang)=sum/issp*10.d0
*       enddo
*      enddo
*     enddo
*     --------------------------------------------------------------- 
*     print differential cross sections for each bin state
*     these agree with fresco output (channels 16 and 20*) 7/08
*     print the first 10 bins only
*     --------------------------------------------------------------- 
*     ijtcount=0
*     do ist=1,istmax
*      do ich=1,iich(ist)
*       ijtcount=ijtcount+1
*       if(ijtcount.gt.10) go to 937
*       iout=10+(ist-1)*iich(ist)+ich
*       do iang=1,nangl
*        write(iout,*) angsd(iang),sigcs(ist,ich,iang)
*       enddo
*      close(iout)
* 937  enddo
*     enddo
*     --------------------------------------------------------------- 
*     print cross sections summed over all bin states on fort.99
*     --------------------------------------------------------------- 
*     do iang=1,nangl
*      sum=0.d0
*      do ist=1,istmax
*       do ich=1,iich(ist)
*       sum=sum+sigcs(ist,ich,iang)
*       enddo
*      enddo
*      write(99,*) angsd(iang),sum
*     enddo
*     --------------------------------------------------------------- 
*     renormalise all two-body kinematics amplitudes to T-matrix
*     compute appropriate wave numbers from bin excitation energies
*     must be done according to read in eex(iex)
*     --------------------------------------------------------------- 
      print*,'  Renormalising (T-matrix) for all excited states '
      mupt=mp*mt/(mp+mt)*amu
      Ecmin=mt/(mt+mp)*ep                    
      Kinc=sqrt(2.d0*mupt*Ecmin)/hbarc
      con=mupt/(2.d0*pi*hbarc2)
      con=con*con
      do iex=1,nexb
       iii=issp*nint(2*exjs(iex)+1)
       Ecmout=Ecmin-eex(iex)
       Kout=sqrt(2.d0*mupt*Ecmout)/hbarc  
       renor=sqrt(con*Kout/Kinc)
       do iang=1,nangl
        do iam=1,iii
         fam(iex,iang,iam)=fam(iex,iang,iam)/renor
        enddo
       enddo
*     print 312,'   iex=',iex,eex(iex),Ecmin,Ecmout,Kout
      enddo
      print 400
*     --------------------------------------------------------------- 
*     fam now has T-matrix normalisations 
*     --------------------------------------------------------------- 
*     c.m energy and wavenumber at breakup threshold
*     store cm wave numbers for each stack of bins (interpolation)
*     --------------------------------------------------------------- 
      Ethr=Ecmin+ebind
      Kthr=sqrt(2.d0*mupt*Ethr)/hbarc
*     print*,'                            Kout      |Kout-Kthr| '
      do ist=1,istmax
       Ko(ist,1)=Kthr
*      print 313,'   ist=',ist,' ich=',0,' iix=',0,
*    +    Ko(ist,1),abs(Ko(ist,1)-Kthr)
       do ich=1,iich(ist)
        iix=iibin(ist,ich)
        Ecmout=Ethr-eex(iix)
        Ko(ist,ich+1)=sqrt(2.d0*mupt*Ecmout)/hbarc    
*       print 313,'   ist=',ist,' ich=',ich,' iix=',iix,
*    +    Ko(ist,ich+1),abs(Ko(ist,ich+1)-Kthr)
       enddo
*     print 400
      enddo
*     --------------------------------------------------------------- 
*     calculate differential cross sections for each bin with two-
*     body phase space included
*     --------------------------------------------------------------- 
*     print*,'  Calculating summed stack cross sections '
*     do ist=1,istmax
*     print*,'  ---------------'
*      spp=chjs(ist)
*      do ich=1,iich(ist)
*       ibi=iibin(ist,ich)
*       iix=iiex(ist,ich)
*       Kout=Ko(ist,ich+1) 
*       Ecmout=Ethr-ehat(ibi)
*       Kout=sqrt(2.d0*mupt*Ecmout)/hbarc  
*       do iang=1,nangl
*        sum=0.d0
*        do im=1,issp
*         do imp=1,nint(2*spp+1)
*          iii=nint((im-1)*(2*spp+1)+imp) 
*          sum=sum+abs(fam(iix,iang,iii))**2
*         enddo
*        enddo
*        sigcs(ist,ich,iang)=con*Kout/Kinc*sum/issp*10.d0
*       enddo
*      print 313,'   ist=',ist,'  ich=',ich,'  iix=',iix,
*    + ehat(ibi)-ebind,Kout
*      enddo
*     enddo
*     print 400
*     --------------------------------------------------------------- 
*     print differential cross sections for each bin state
*     --------------------------------------------------------------- 
*     ijtcount=0
*     do ist=1,istmax 
*      do ich=1,iich(ist)
*       ijtcount=ijtcount+1
*       if(ijtcount.gt.10) goto 761
*       iout=30+(ist-1)*iich(ist)+ich
*       do iang=1,nangl
*        write(iout,*) angsd(iang),sigcs(ist,ich,iang)
*       enddo
*      close(iout)
*      enddo
*761  enddo
*     --------------------------------------------------------------- 
*     dump differential cross section summed over all bins on fort.98
*     --------------------------------------------------------------- 
*     do iang=1,nangl 
*      sum=0.d0
*      do ist=1,istmax
*       do ich=1,iich(ist)
*       sum=sum+sigcs(ist,ich,iang)
*       enddo
*      enddo
*      write(98,*) angsd(iang),sum
*     enddo
*     --------------------------------------------------------------- 
*     store cm wavenumber deviations from threshold in each stack
*     --------------------------------------------------------------- 
      do ist=1,istmax
       Ko(ist,1)=0.d0
       do ich=2,iich(ist)+1
        Ko(ist,ich)=abs(Ko(ist,ich)-Kthr)    
       enddo
      enddo


ccc double differential cross section (E_rel,theta_K)
      
*     ------------------------------------
*     overlap between scattering states and PS eigenstates
*     ------------------------------------
      mc=mc*amu
      mv=mv*amu
      mt=mt*amu
      mp=mc+mv
      mtot=mc+mv+mt
      gamma=mt/(mc+mt)
      alfa =mv/(mc+mv)
      mupt=mp*mt/mtot
      muct=mc*mt/(mc+mt)
      mucv=mc*mv/mp
      muvct=mv*(mc+mt)/mtot
*     ---------------------------------------------------------------    
*     projectile wavenumber in cm frame
      ecm=mt/mtot*ep                    
      Kcm=sqrt(2.d0*mupt*ecm)/hbarc        
      Ethr=Ecmin+ebind
ccc      qthr=sqrt(2.d0*mupt*Ethr)/hbarc
      derel=Ethr/dble(nerel)
      do iecv=1,nerel
      erel=(iecv-1)*derel
      indexx=0
      do isett=1,jpsets
      nchnn=jpiset(isett)%nchan
      if (allocated(wf)) deallocate(wf)
      allocate(wf(nchnn,nr))
      do iil=1,nchnn
      xfactI1=jpiset(isett)%jc(iil)
      if(xfactI1.ne.xIce(idcore)) cycle
      nfactl1=jpiset(isett)%lsp(iil)
      indexx=indexx+1
      call test_cont(isett,nchnn,iil,erel,wf)
      nex=jpiset(isett)%nex
      do iie=1,nex
      caux=zero
      do iilp=1,nchnn
      nfactl2=jpiset(isett)%lsp(n)
      xfactI2=jpiset(isett)%jc(n)
      do j=1,nr
      fauxi(j)=wfeign(iie,n,j)*wf(n,j)*rvec(j)*
     .(-1)**(nfactl1+xfactI1+nfactl2+xfactI2)
      enddo
      call simc(fauxi,resc,1,nr,dr,nr)
      caux=caux+resc
      enddo ! iilp
      gsolap(isett,iil,iie,iecv)=caux
      enddo ! iie (PS)
      enddo ! iil
      enddo ! isett
      enddo ! iq
     
ccc   related index inside the stacks iPS=iPS(iset,n)
      iex=0
      allocate(ips(jpsets,nsetmax))
      do isett=1,jpsets
      nex=jpiset(isett)%nex
      do n=1,nex
      iex=iex+1
      ips(isett,n)=iex
      enddo
      enddo
      
      
ccc   x section dsdew
      allocate(sumi(nangl),dsdew(nerel,nangl))
      faux=0.5/(pi**3)/hbarc**6*mupt**2*mucv/issp
      do iecv=1,nerel
      erel=(iecv-1)*derel
      kcv=sqrt(2.d0*mucv*erel)/hbarc
      ecmf=ecmin-dabs(ebind)-erel-xece(idcore)
      Kpt=sqrt(2.d0*mupt*ecmf)/hbarc     
      facK=kpt/kinc/kcv
      do isett=1,jpsets
      nex=jpiset(isett)%nex
      jpf=jpiset(isett)%jtot
      isspf=nint(2.*jpf+1.)
      iam=0
      do imi=1,issp
      do imf=1,isspf
      iam=iam+1
      nchnn=jpiset(isett)%nchan
      do iil=1,nchnn
      xfactI1=jpiset(isett)%jc(iil)
      if(xfactI1.ne.xIce(idcore)) cycle
      lcv=jpiset(isett)%lsp(iil)
      jcv=jpiset(isett)%jsp(iil)
      sumi(:)=0.d0
      do n=1,nex
      if (energ(isett,n).lt.0) cycle    ! bound states
      if (energ(isett,n).gt.Ethr) cycle ! closed channels
      iex=ips(isett,n)
      sumi(:)=sumi(:)+gsolap(isett,i,iil,iecv)*fam(iex,:,iam)
      enddo
      do iang=1,nangl
      dsdew(iecv,iang)=faux*facK*abs(sumi(iang))**2
      enddo
      enddo ! iil
      enddo ! imf
      enddo ! imi
      enddo ! isett
      enddo ! iecv (valence-core relative energy)  
      
      

*     --------------------------------------------------------------- 
*     Calculate angular momentum coefficients 
*     --------------------------------------------------------------- 
      inz=0
*     print 400
      print*,'  Calculating angular momentum coefficients '
      print 400
      print*,'               J       parity      l         j       '
*     --------------------------------------------------------------- 
*     loop on asymptotic z-component spin labels first, mu, sigma, m
*     --------------------------------------------------------------- 
      imaxmmp=0
      ifl=0
      iam=0
      do imu=1,issc
       mu=imu-sc-1.d0
       do isig=1,issv
        sig=isig-sv-1.d0
        do im=1,issp
         mm=im-sp-1.d0
         ifl=ifl+1
         do ist=1,istmax
          ajval=chjs(ist) 
          alval=chls(ist) 
          issl=nint(2*alval+1)
          aji =chjis(ist) 
*     --------------------------------------------------------------- 
*      print stack quantum numbers
*     --------------------------------------------------------------- 
          if(ifl.eq.1) then
           print 314,'   ist=',ist,ajval,chpari(ist),alval,aji
          endif
*     --------------------------------------------------------------- 
*     remaining free z-projection, nu
*     --------------------------------------------------------------- 
          do inu=1,issl
           iam=iam+1
           rnu=inu-alval-1.d0
           ajim=rnu+sig
           mmp=mu+ajim
           if(abs(ajim).gt.aji+0.1d0.or.abs(mmp).gt.ajval+0.1d0) then
            cgc(iam)=0.d0
           else
            inz=inz+1
            c1=cleb(alval,rnu,sv,sig,aji,ajim)
            c2=cleb(aji,ajim,sc,mu,ajval,mmp)
            cgc(iam)=c1*c2
            if(abs(cgc(iam)).lt.1.d-10) inz=inz-1
*     --------------------------------------------------------------- 
*     store the maximum value of |M-M'| for later phi phases
*     --------------------------------------------------------------- 
            iii=nint(abs(mm-mmp))
            if(iii.gt.imaxmmp) imaxmmp=iii
           endif
          enddo
         enddo
        enddo
       enddo
      enddo
*     --------------------------------------------------------------- 
*     print some simple diagnostics 
*     --------------------------------------------------------------- 
      print 400
      print*,' Total number of amplitudes, non zero ones, max|M-M"|'
      print'(3(a,i4))','   iam =',iam,'    inz =',inz,
     +     '    maxmmp =',imaxmmp
      print 400
*     --------------------------------------------------------------- 
      if(iam.gt.ncleb) then
       print*,'  ncleb is ',ncleb,' so too small '
       print 400
       stop 
      endif
*     --------------------------------------------------------------- 
*     Calculate and store coefficients for spherical harmonics
*     Compute largest bin orbital angular momentum needed
*     --------------------------------------------------------------- 
      almax=0.d0
      do ist=1,istmax
       if(chls(ist).gt.almax) almax=chls(ist)
      enddo
      print*,'  Calculating spherical harmonic coefficients '
      print*,'  Maximum internal l-value is',nint(almax)
      print 400
*     --------------------------------------------------------------- 
*     test angles for P_{lm}(0)
*     --------------------------------------------------------------- 
      co=0.d0
      co2=co*co
      si=1.d0
      si2=si*si
*     --------------------------------------------------------------- 
      do ilval=0,nint(almax)
       ali=ilval
*     --------------------------------------------------------------- 
*      (after some pain) the correct relative partial wave phase is: 
*     --------------------------------------------------------------- 
       phasil(ilval)=(-i)**ilval 
*     --------------------------------------------------------------- 
       c1=sqrt((2.d0*ali+1.d0)/(4.d0*pi))
       do inu=1,2*ilval+1
        rnu=abs(inu-ilval-1)
        c2=faclog(nint(ali-rnu+1))-faclog(nint(ali+rnu+1))
        cplm(ilval,inu)=c1*exp(c2/2.d0)
*     --------------------------------------------------------------- 
*     eikonal condition on the m-substates or orbital momentum
*     --------------------------------------------------------------- 
*       if(ifeik.eq.1) then
*        if((-1)**(nint(ilval+rnu)).lt.0.d0) cplm(ilval,inu)=0.d0
*       endif
*     --------------------------------------------------------------- 
       enddo
       do inu=2*ilval+1,1,-1
        iii=inu-ilval-1
        if(iii.gt.0) cplm(ilval,inu)=(-1)**iii*cplm(ilval,inu)  
*     --------------------------------------------------------------- 
*       test values for Plm: corrected Arfken formulae for P_{lm}(0)
*     --------------------------------------------------------------- 
        iij=abs(iii)
        if((-1)**(ilval+iij).lt.0) then
         pan=0.d0
        else
         pan=(-1)**((ilval-iij)/2)/(2.d0**ilval)
         pan=pan*exp(faclog(ilval+iij+1))
         pan=pan/exp(faclog((ilval+iij)/2+1)+faclog((ilval-iij)/2+1))
        endif
*     --------------------------------------------------------------- 
*       diagnostic print of Y_{lm} coefficients and P_{lm}(0)
*     --------------------------------------------------------------- 
        print'(2i4,3f15.7)',ilval,iii,cplm(ilval,inu),
     +   plmrec(ilval,iii,co,co2,si,si2),pan
       enddo
       print 400
      enddo
*     --------------------------------------------------------------- 
*     order of the Aitken 2D-interpolation
*     --------------------------------------------------------------- 
*     nord=2
*     --------------------------------------------------------------- 
*     masses, mass ratios, and reduced masses
      mc=mc*amu
      mv=mv*amu
      mt=mt*amu
      mp=mc+mv
      mtot=mc+mv+mt
      gamma=mt/(mc+mt)
      alfa =mv/(mc+mv)
      mupt=mp*mt/mtot
      muct=mc*mt/(mc+mt)
      mucv=mc*mv/mp         
      muvct=mv*(mc+mt)/mtot
*     ---------------------------------------------------------------    
*     projectile wavenumber in cm frame
      ecm=mt/mtot*ep                    
      Kcm=sqrt(2.d0*mupt*ecm)/hbarc
*     print*,' Kcm = ',Kcm, ' ecm = ',ecm
*     ---------------------------------------------------------------
*     specify total incident momentum (momentum of c.m.) in MeV/c
      if(sys.eq.'lab') totp=sqrt(2.d0*mp*ep)
*     if detector positions refer to c.m. frame 
      if(sys.eq.'com') totp=0.d0
*     print*,' ktot = ',totp/hbarc
*     ---------------------------------------------------------------
*     set up projectile cm K vector (qp) (in z-direction)
*     set up ptot vector and ktot (also in z-direction)             
      do j=1,2      
       qp(j)=0.d0
       ptot(j)=0.d0       
       ktot(j)=0.d0
      enddo
      qp(3)=Kcm
      ptot(3)=totp
      ktot(3)=totp/hbarc
      erelmax=0.d0
*     ---------------------------------------------------------------
*     loop over core particle detection angle
      do 10 ic=1,itc
       tcd=tcl+(ic-1)*dtc
       tc=tcd*degrad
       costc=cos(tc)
       sintc=sin(tc)
*     ---------------------------------------------------------------
*     loop over valence particle detection angle
      do 20 iv=1,itv
       tvd=tvl+(iv-1)*dtv
       tv=tvd*degrad
       costv=cos(tv)
       sintv=sin(tv)
*     ---------------------------------------------------------------
*     check if either theta is zero since no need to do multiple phi
*     calculations (if wanted, itphi>1) in this geometry
*     ---------------------------------------------------------------
      itphim=itphi/2+1
      iflag=0
      zert=(abs(tcd).le.acc.or.abs(tvd).le.acc)
      if(itphi.gt.1.and.zert) then
       itphim=1
       iflag=1
      endif
*     ---------------------------------------------------------------
*     loop over detected energy of chosen particle 
      iiiflag=0
      do 30 ien=1,iten
      En=Enlow+(ien-1)*dEn
      if(idet.eq.1) then
       p1L=sqrt(2.d0*mc*En)
      else
       p2L=sqrt(2.d0*mv*En)       
      endif
*     ---------------------------------------------------------------
*     the detected particle defines phi=0 degrees
*     ---------------------------------------------------------------
*     loop over azimuth detection angle          
      do 40 ip=1,itphim
       sigphi(ip)=0.d0
       phid=phil+(ip-1)*dphi
       phi=phid*degrad  
       if(idet.eq.1) then
        cospv=cos(phi)
        sinpv=sin(phi)
       else
        cospc=cos(phi)
        sinpc=sin(phi)
       endif
*     -----------------------------------------------------------------    
*     core particle momentum UNIT vector in chosen frame             
      pcL(1)=sintc*cospc
      pcL(2)=sintc*sinpc
      pcL(3)=costc
*     valence particle momentum UNIT vector in chosen frame
      pvL(1)=sintv*cospv
      pvL(2)=sintv*sinpv
      pvL(3)=costv
      do j=1,3
       if(idet.eq.1) then
        pcL(j)=p1L*pcL(j)
       else
        pvL(j)=p2L*pvL(j)
       endif
      enddo
*     -----------------------------------------------------------------
*     calculate second particle energy and momentum depending on idet
      do j=1,3
       if(idet.eq.1) then              
        tem(j)=ptot(j)-pcL(j)
       else
        tem(j)=ptot(j)-pvL(j)        
       endif
      enddo
*     coefficients in quadratic for second particle momentum
      if(idet.eq.1) then
       aq=(1.d0/mv+1.d0/mt)/2.d0
       bq=-dot(tem,pvL)/mt
       cq=(dot(tem,tem)/mt+dot(pcL,pcL)/mc)/2.d0
      else
       aq=(1.d0/mc+1.d0/mt)/2.d0
       bq=-dot(tem,pcL)/mt
       cq=(dot(tem,tem)/mt+dot(pvL,pvL)/mv)/2.d0
      endif
      if(sys.eq.'lab') cq=cq-(ep+ebind)
      if(sys.eq.'com') cq=cq-(ecm+ebind) 
*     -----------------------------------------------------------------
*     discriminant of the quadratic equation
*     -----------------------------------------------------------------
      dq=bq*bq-4.d0*aq*cq
*     -----------------------------------------------------------------
*     possible there are no roots for light targets and some angles
*     and energies in the lab frame
*     -----------------------------------------------------------------
      iroots=0
      if(dq.lt.0.d0) goto 52
*     -----------------------------------------------------------------
*     otherwise proceeed to the two roots and check for positive values
*     -----------------------------------------------------------------
      iroots=2
      iroot=1
*     -----------------------------------------------------------------
*     finalise momentum vector of second particle
*     -----------------------------------------------------------------
      iroot1n=0
      if(idet.eq.1) then
       p2L=(-bq+sqrt(dq))/(2.d0*aq)
       if(p2L.le.0.d0) iroot1n=1
      else
       p1L=(-bq+sqrt(dq))/(2.d0*aq)          
       if(p1L.le.0.d0) iroot1n=1
      endif
      if(iroot1n.gt.0) iroots=1
*     -----------------------------------------------------------------
   53 if(iroot.eq.2.or.iroot1n.gt.0) then
       iroot2n=0
       if(idet.eq.1) then
        p2L=(-bq-sqrt(dq))/(2.d0*aq)
        if(p2L.le.0.d0) iroot2n=1
       else
        p1L=(-bq-sqrt(dq))/(2.d0*aq)          
        if(p1L.le.0.d0) iroot2n=1
       endif
       if(iroot2n.gt.0) then
	iroots=iroots-1
	goto 52
       endif
      endif
*     -----------------------------------------------------------------
*     core and valence particle energies and momenta in chosen frame
*     -----------------------------------------------------------------
      Ec=(p1L*p1L)/(2.d0*mc)
      Ev=(p2L*p2L)/(2.d0*mv)
      do j=1,3
       if(idet.eq.1) then
        kvL(j)=p2L*pvL(j)/hbarc                 
        kcL(j)=pcL(j)/hbarc                  
       else
        kcL(j)=p1L*pcL(j)/hbarc         
        kvL(j)=pvL(j)/hbarc              
       endif
      enddo
*     -----------------------------------------------------------------
*     construct remaining vectors
*     -----------------------------------------------------------------
      do j=1,3
*     wavevector of cm of core and valence particles
       bkp(j)=kcL(j)+kvL(j)-(mc+mv)/mtot*ktot(j)
*     wavevector of relative motion of core and valence particles
       kp(j)=mc/(mc+mv)*kvL(j)-mv/(mc+mv)*kcL(j)
*     -----------------------------------------------------------------
*     wavevector for core-target relative motion
*      qc(j)=(1.d0-alfa*gamma)*bkp(j)-gamma*kp(j)
*     wavevector for valence-(core+target) relative motion
*      qv(j)=alfa*bkp(j)+kp(j)
*     -----------------------------------------------------------------
      enddo
*     -----------------------------------------------------------------
*     qbig is mv/(mc+mv)*(K-K') in the cm frame
*     -----------------------------------------------------------------
*     do j=1,3
*      qbig(j)=mv/(mc+mv)*(qp(j)-bkp(j))
*     enddo
*     -----------------------------------------------------------------
*     important parameters for RCJ model test set up here, first Q
*     -----------------------------------------------------------------
*     qq=m(qbig)
*     -----------------------------------------------------------------
*     cosine of angle between big Q and k' for RCJ test
*     core+valence particles relative energy in their cm frame
*     -----------------------------------------------------------------
*     xx=dot(qbig,kp)/qq/m(kp)
*     -----------------------------------------------------------------
*     diagnostic checks for energies in cm frame
*     -----------------------------------------------------------------
*     qc2=dot(qc,qc)
*     qv2=dot(qv,qv)
*     ener=hbarc2*qv2/(2.d0*muvct)   
*     ener=ener+hbarc2*qc2/(2.d0*muct)
*     print*,ener,ecm+ebind,iroot
*     print*,En,eks,acos(xx)*180.d0/pi
*     -----------------------------------------------------------------
      eKb=hbarc2*dot(bkp,bkp)/(2.d0*mupt)
      eks=hbarc2*dot(kp,kp)/(2.d0*mucv)
*     -----------------------------------------------------------------
*     increment maximum relative energy 
*     -----------------------------------------------------------------
      if(eks.gt.erelmax) erelmax=eks
*     -----------------------------------------------------------------
*     CALCULATION OF T-MATRIX <bkp,kp|T|qp> STARTS HERE
*     -----------------------------------------------------------------
*     ybo = Kout_cm, xb =theta_cm (rad), xbd (degrees)
*     -----------------------------------------------------------------
*     compute polar angles (th,ph) of cm wave vector
*     -----------------------------------------------------------------
      ybo=m(bkp)
      xb=min(1d0,bkp(3)/ybo)
      xb=max(-1d0,xb)
      xb=acos(xb)
      if(abs(bkp(1)).lt.1.d-6.and.abs(bkp(2)).lt.1.d-6) then
       phKb=0.d0
      else 
       phKb=atan2(bkp(2),bkp(1))
      endif
      xbd=xb/degrad
*     -----------------------------------------------------------------
*     pick out nearest angles indices in array for interpolations
*     -----------------------------------------------------------------
      iang=int((xbd-thmin)/thinc)+1
      nial=max0(iang-2,1)
      niag=min0(nial+5,nangl)
      nial=niag-5
      do iang=1,6
       xyt(1,iang)=angsr(nial+iang-1)
      enddo
*     -----------------------------------------------------------------
*     compute polar angles (theta,phi) of relative wave vector
*     -----------------------------------------------------------------
      mkp=m(kp)
      co=kp(3)/mkp
      if(abs(kp(1)).lt.1.d-6.and.abs(kp(2)).lt.1.d-6) then
       phks=0.d0
      else 
       phks=atan2(kp(2),kp(1))
      endif
      co2=co*co
      si2=1.d0-co2
c	if(si2 .lt.0.) print *,kp,mkp,co,co2,si2
      si=sqrt(abs(si2))
*     ------------------------------------
*     overlap between scattering states and PS eigenstates
*     ------------------------------------
      energy=0.5d0*hbarc2*mkp**2/mucv-xece(idcore)!-dabs(egs)
      indexx=0
      do isett=1,jpsets
      nchnn=jpiset(isett)%nchan
      if (allocated(wf)) deallocate(wf)
      allocate(wf(nchnn,nr))
      do iil=1,nchnn
      xfactI1=jpiset(isett)%jc(iil)
      if(xfactI1.ne.xIce(idcore)) cycle
      nfactl1=jpiset(isett)%lsp(iil)
      indexx=indexx+1
      call test_cont(isett,nchnn,iil,energy,wf)
      do iie=1,hdimtot
      do n=1,nchnn
      nfactl2=jpiset(isett)%lsp(n)
      xfactI2=jpiset(isett)%jc(n)
      over=> overlapmat(indexx,iie,n)
      do j=1,nr
      fauxi(j)=wfeign(iie,n,j)*wf(n,j)*rvec(j)*
     .(-1)**(nfactl1+xfactI1+nfactl2+xfactI2)
      enddo
      call simc(fauxi,resc,1,nr,dr,nr)
      over=resc
      enddo ! n
      enddo ! iie (PS)
      enddo ! iil
      enddo ! isett
      do iie=1,hdimtot
      do ii=1,indexx
      caux=zero
      do n=1,nchnn
      caux=caux+overlapmat(ii,iie,n)
      enddo
      oversum(iie,ii)=caux
      enddo ! ii
      enddo ! iie
ccc      if(ic.eq.1.and.iv.eq.1.and.ip.eq.1) then
ccc      do inn=1,hdimtot
ccc      write(1234,911) energy,oversum(inn,1:indexx)
ccc      enddo
ccc      write(1234,*) '&'
ccc911   format (1f8.3,2x,100(e12.4,2x,e12.4,4x))
ccc      endif
*     -----------------------------------------------------------------
*     check if relative momentum less than any of the stack maxima
*     --------------------------------------------------------------- 
      iisum=0
      do ist=1,istmax
       iiist(ist)=1
       if(mkp.gt.kmaxJ(ist)+0.05d0.or.mkp.lt.kminJ(ist).or. ! ****
     *    eks.gt.uecut) then
ccc       if(eks.gt.uecut) then
ccc        phbin(ist)=(1.d0,0.d0)
       iiist(ist)=0
       iiiflag=1
       else ! ****
        iisum=iisum+iiist(ist)
*     -----------------------------------------------------------------
*     Phase shift for the bin state at this relative energy
*     -----------------------------------------------------------------
ccc        do ich=1,iich(ist)
ccc         ibi=iibin(ist,ich)
ccc         if(mkp.ge.kminJs(ibi).and.mkp.le.kmaxJs(ibi)) then
*     -----------------------------------------------------------------
*        record of bin, that mkp value is in, for this stack
*     -----------------------------------------------------------------
*         ibused(ist)=ibi
ccc          dkb=(kmaxJs(ibi)-kminJs(ibi))/(nks(ibi)-1)
ccc          ik=nint((mkp-kminJs(ibi))/dkb+1)
*     -----------------------------------------------------------------
ccc          phbin(ist)=exp(i*phases(ibi,ik))
*     -----------------------------------------------------------------
ccc         endif
ccc        enddo
       endif ! ****
      enddo
*     --------------------------------------------------------------- 
*     if iisum=0 then no bins defined at this relative energy
*     --------------------------------------------------------------- 
      if(iisum.eq.0) then
       tmatsq=0.d0
      else
*     --------------------------------------------------------------- 
*     calculate the spherical harmonics for this relative wave vector
*     --------------------------------------------------------------- 
      phask1=(1.d0,0.d0)
      phask=exp(i*phks)
      do ilval=0,nint(almax)
       phask1=phask1/phask
       phas=phask1
       do inu=1,2*ilval+1
        phas=phas*phask
        aleg=plmrec(ilval,inu-ilval-1,co,co2,si,si2)
        ylm(ilval,inu)=cplm(ilval,inu)*phas*aleg
        ylm(ilval,inu)=ylm(ilval,inu)*phasil(ilval)
       enddo
      enddo
*     --------------------------------------------------------------- 
*     calculate the Kcm phi phases (m-m' is always an integer)
*     --------------------------------------------------------------- 
      pph=exp(i*phKb)
      phiph(0)=(1.d0,0.d0)
      do iii=1,imaxmmp
       phiph(iii)=phiph(iii-1)*pph
       phiph(-iii)=1.d0/phiph(iii)
      enddo
*     --------------------------------------------------------------- 
*     interpolate the scattering amplitudes for each ist, m and m'
*     --------------------------------------------------------------- 
      yb=abs(ybo-Kthr)
      do ist=1,istmax
*     --------------------------------------------------------------- 
*      see if this is an s-state stack
*     --------------------------------------------------------------- 
       issstat=0
       if(nint(chls(ist)).eq.0) issstat=1
*     --------------------------------------------------------------- 
       iam=issp*nint(2*chjs(ist)+1)
       if(iiist(ist).eq.0.or.ifinc(ist).eq.0) then
        do iii=1,iam
         fin(ist,iii)=(0.d0,0.d0)
        enddo
       else
        nny=iich(ist)+1
        xyt(2,1)=0.d0
        do iang=1,6
         fxyc(iang,1)=(0.d0,0.d0)
        enddo
*      -------------------------------------------------------------- 
*      set up angle (x) v Kcm (y) array fxyc for 2D-interpolation
*      --------------------------------------------------------------
        do iii=1,iam
         do ich=2,nny
ccc          ibi=iibin(ist,ich-1)
ccc          iex=iiex(ist,ich-1)
*      -------------------------------------------------------------- 
          if(issstat.eq.1) then
           nnu=nny-1
ccc           c1=bindk(ibi)/khat(ibi)
           xyt(2,ich-1)=Ko(ist,ich)
           do iang=0,5
            suma=zero
            do iex=1,hdimtot
            suma=suma+
     .      fam(iex,nial+iang,iii)*oversum(iex,ist)
            enddo
            fxyc(iang+1,ich-1)=suma/mkp
           enddo
          else
           nnu=nny
ccc           if(isc(ibi).eq.2.or.isc(ibi).eq.0) then
ccc            c1=bindk(ibi)/mkp
ccc           else if(isc(ibi).eq.12) then
ccc            c1=bindk(ibi)/khat(ibi)
ccc           endif
           xyt(2,ich)=Ko(ist,ich)
           do iang=0,5
            suma=zero
            do iex=1,hdimtot
            suma=suma+
     .      fam(iex,nial+iang,iii)*oversum(iex,ist)
            enddo
            fxyc(iang+1,ich+0)=suma/mkp
           enddo
          endif
*      -------------------------------------------------------------- 
         enddo
c         xmin=minval(xyt(:,:)) ! ****
c         xmax=maxval(xyt(:,:)) ! ****
c         if (xmin.lt.xb.or.xmax.gt.xb.or.xmin.lt.yb.or.xmax.gt.yb) then ! ****
c         print*, 'Fresco:',xb,yb ! ****
c         print*, 'Jeff:',xmin,xmax ! ****
c         stop ! ****
c         endif ! ****
         fin(ist,iii)=f2c(xb,yb,xyt,fxyc,6,nnu,nord,10,nbst+1)
*      -------------------------------------------------------------- 
         if(issstat.eq.1) then
          fin(ist,iii)=fin(ist,iii)
         else
          fin(ist,iii)=fin(ist,iii)
         endif
*      -------------------------------------------------------------- 
        enddo
       endif
      enddo
*     --------------------------------------------------------------- 
*     calculate |t-matrix|**2 summed on spin projections / (2sp+1)
*     --------------------------------------------------------------- 
*     loop on asymptotic z-component spin labels, mu, sigma, m
*     --------------------------------------------------------------- 
      sum=0.d0
      iam=0
      do imu=1,issc
       mu=imu-sc-1
       do isig=1,issv
        sig=isig-sv-1
        do im=1,issp
         mm=im-sp-1
         amp=(0.d0,0.d0)
         do ist=1,istmax
          spp=chjs(ist) 
          ilval=nint(chls(ist))
          do inu=1,2*ilval+1
           iam=iam+1
           if(iiist(ist).eq.1) then
            cg=cgc(iam)*ifinc(ist)
            if(abs(cg).gt.1.d-10) then
             rnu=inu-ilval-1
             mmp=mu+rnu+sig 
             iii=nint((im-1)*(2*spp+1)+mmp+spp+1)
             imp=nint(mm-mmp)
             amp=amp+cg*fin(ist,iii)*phiph(imp)
     +          *ylm(ilval,inu)
            endif
           endif
          enddo
         enddo
         sum=sum+abs(amp)**2
        enddo
       enddo
      enddo
      tmatsq=sum*((4.d0*pi)**2)/issp
      endif
*     -----------------------------------------------------------------
*     CALCULATION OF |T-MATRIX|**2 COMPLETED (is in tmatsq)
*     -----------------------------------------------------------------
*     CALCULATE PHASE SPACE
*     common bits independent of particle detected
*     -----------------------------------------------------------------
      mult=(mupt*mv/(2.d0*pi*hbarc2)**2)*m(kvL)/m(qp)
      mult=mult*mc*m(kcL)*hbarc/(2.d0*pi*hbarc)**3
*     -----------------------------------------------------------------
*     vector (k[detected]-ktot) for phase space factor
*     -----------------------------------------------------------------
      do j=1,3
       if(idet.eq.1) then
        tem(j)=kcL(j)-ktot(j)
       else
        tem(j)=kvL(j)-ktot(j)  
       endif
      enddo
*     -----------------------------------------------------------------
*     detected particle dependent parts
*     -----------------------------------------------------------------
      if(idet.eq.1) then
       mult=mult*abs((mt/((mv+mt)+mv*dot(tem,kvL)/dot(kvL,kvL))))
      else
       mult=mult*abs((mt/((mc+mt)+mc*dot(tem,kcL)/dot(kcL,kcL))))       
      endif
*     -----------------------------------------------------------------
*     COMPLETES PHASE SPACE CALCULATION (is in mult)
*     calculate the triple differential cross section (10 for mb)
*     -----------------------------------------------------------------
      sigphi(ip)=sigphi(ip)+10.d0*mult*tmatsq
   52 if(iroots.eq.0) sigphi(ip)=0.d0
*     -----------------------------------------------------------------
*     prepare to go around again at this phi if two roots required
*     -----------------------------------------------------------------
      if(iroots.eq.2.and.iroot.eq.1) then
       iroot=2
       goto 53
      endif
*     -----------------------------------------------------------------
*     close the phi loop
*     -----------------------------------------------------------------
   40 continue
*     -----------------------------------------------------------------
*     integrate over phi angles NOW if more than one phi angle
*     -----------------------------------------------------------------
      if(itphi.eq.1) then
       sig=sigphi(1)
      else if(iflag.eq.0) then
       call simpson(sigphi,sig,1,itphim,dphir)
       sig=2.d0*sig
      else
       sig=sigphi(1)*(phiu-phil)*degrad
      endif
      if(ien.eq.iten) then
       print 315,ien,iv,ic,sig,En,iflag
      endif
*     -----------------------------------------------------------------
      write(77,315) ien,iv,ic,sig,En
      call flush(77)
*     if(ien.eq.1) write(99,*) '# text',tcd,tvd
*     write(99,*) En,sig,real(mkp),real(eks)
*     -----------------------------------------------------------------
*     close the angle (core and valence thetas) and energy loops
*     -----------------------------------------------------------------
   30 continue      
   20 continue
   10 continue
      print 400
      print*,'  erelmax =',erelmax
      print 400
*     --------------------------------------------------------------- 
*     collected formatting statements 
*     --------------------------------------------------------------- 
*     789012345678901234567890123456789012345678901234567890123456789012
*     --------------------------------------------------------------- 
  303 format(a,f9.4,a,f9.4)
  304 format(5f10.3)
  305 format(10f8.4)
  306 format(6e12.4)
  307 format(i5,f15.8,3f10.3)
  308 format(1p,6e12.4)
  309 format(i5,2f5.2,3e12.5,2i4)
  310 format(a,i3,a,i3,2f9.4,2i5)
  311 format(a,2f10.5,a,i4,a)
  312 format(a,i4,4f10.6)
  313 format(a,i3,a,i3,a,i3,3f11.6)
  314 format(a,i3,5f10.3)
  315 format(3i5,d19.8,1x,f12.6,i5)
  316 format(5f10.3)
  317 format(i6,i4,2f12.4,i4,f12.4)
  318 format(i5,2f5.2,4e12.5)
  400 format('  ------------------------------------------------------',
     +'-----')
*     --------------------------------------------------------------- 
      return
      end
*----------------------------------------------------------------------
      real*8 function dot(a,b)
*     scalar product of two vectors
      real*8 a(3),b(3)
      dot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end
*----------------------------------------------------------------------
      real*8 function m(a)
*     modulus of vector
      real*8 a(3)
      m=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      return
      end
*----------------------------------------------------------------------
      subroutine simpson(fa,res,m,n,h)
      implicit real*8(a-h,o-z)
      dimension fa(91),dq(91)
      do i=m,n
       dq(i)=fa(i)
      enddo
      rq1=dq(m+1)
      rq2=dq(m+2)
      i=m+3
   98 continue
      if(i.ge.n) go to 99
      rq1=rq1+dq(i)
      rq2=rq2+dq(i+1)
      i=i+2
      go to 98
   99 continue
      res=0.33333333333d0*h*(dq(m)+4.d0*rq1+2.d0*rq2-dq(n))
      return
      end
*----------------------------------------------------------------------
      real*8 function pleg(ill,imm,c,c2,s,s2)
      implicit real*8(a-h,o-z)
*-----------------------------------------------------------------------
* Associated Legendre functions for l<=6
*-----------------------------------------------------------------------
      pleg=1.d0
      if(ill.eq.0) return
      im=abs(imm)
*-----------------------------------------------------------------------
      if(ill.eq.1) then
       if(im.eq.0) then
        pleg=c
       else if(im.eq.1)then
        pleg=s
       endif
*-----------------------------------------------------------------------
      else if(ill.eq.2) then
       if(im.eq.0) then
       pleg=(3.d0*c2-1.d0)/2.d0
       else if(im.eq.1) then
        pleg=3.d0*s*c
       else if(im.eq.2) then
        pleg=3.d0*s2
       endif
*-----------------------------------------------------------------------
      else if(ill.eq.3) then
       if(im.eq.0) then
        pleg=(5.d0*c2-3.d0)*c/2.d0
       else if(im.eq.1) then
        pleg=s*(15.d0*c2-3.d0)/2.d0
       else if(im.eq.2) then
        pleg=15.d0*s2*c
       else if(im.eq.3) then
        pleg=15.d0*s2*s
       endif
*-----------------------------------------------------------------------
      else if(ill.eq.4) then
       if(im.eq.0) then
        pleg=(35.d0*c2*c2-30.d0*c2+3.d0)/8.d0
       else if(im.eq.1) then
        pleg=s*(35.d0*c2-15.d0)*c/2.d0
       else if(im.eq.2) then
        pleg=s2*(105.d0*c2-15.d0)/2.d0
       else if(im.eq.3) then
        pleg=105.d0*s*s2*c
       else if(im.eq.4) then
        pleg=105.d0*s2*s2
       endif
*-----------------------------------------------------------------------
      else if(ill.eq.5) then
       if(im.eq.0) then
        pleg=(63.d0*c2*c2-70.d0*c2+15.d0)*c/8.d0
       else if(im.eq.1) then
        pleg=s*(315.d0*c2*c2-210.d0*c2+15.d0)/8.d0
       else if(im.eq.2) then
        pleg=s2*(315.d0*c2-105.d0)*c/2.d0
       else if(im.eq.3) then
        pleg=s*s2*(945.d0*c2-105.d0)/2.d0
       else if(im.eq.4) then
        pleg=945.d0*s2*s2*c
       else if(im.eq.5) then
        pleg=945.d0*s2*s2*s
       endif
*-----------------------------------------------------------------------
      else if(ill.eq.6) then
       if(im.eq.0) then
        pleg=(231.d0*c2*c2*c2-315.d0*c2*c2+105.d0*c2-5.d0)/16.d0
       else if(im.eq.1) then
        pleg=21.d0*s*c*(33.d0*c2*c2-30.d0*c2+5.d0)/8.d0
       else if(im.eq.2) then
        pleg=(33.d0*c2*c2*c2-51.d0*c2*c2+19.d0*c2-1.d0)
        pleg=-105.d0*pleg/8.d0
       else if(im.eq.3) then
        pleg=315.d0*s*s2*c*(11.d0*c2-3.d0)/2.d0
       else if(im.eq.4) then
        pleg=(11.d0*c2*c2*c2-23.d0*c2*c2+13.d0*c2-1.d0)
        pleg=945.d0*pleg/2.d0
       else if(im.eq.5) then
        pleg=10395.d0*s2*s2*s*c
       else if(im.eq.6) then
        pleg=10395.d0*s2*s2*s2
       endif
*-----------------------------------------------------------------------
      endif
      return
      end
*-----------------------------------------------------------------------
      recursive real*8 function plmrec(ill,imm,c,c2,s,s2)
      implicit real*8(a-h,o-z)
*-----------------------------------------------------------------------
*     Associated Legendre functions
*     ECS 24th Nov 2008
*     l<=2 caluclated explicitly
*     l>2  calculated recursively
*-----------------------------------------------------------------------
      plmrec=1.d0
      if(ill.eq.0) return
      im=abs(imm)
*-----------------------------------------------------------------------
      if(ill.eq.1) then
       if(im.eq.0) then
        plmrec=c
       else if(im.eq.1)then
        plmrec=s
       endif
      else if(ill.eq.2) then
       if(im.eq.0) then
        plmrec=(3.d0*c2-1.d0)/2.d0
       else if(im.eq.1) then
        plmrec=3.d0*s*c
       else if(im.eq.2) then
        plmrec=3.d0*s2
       endif
      elseif(im.gt.1) then
       plmrec=(1.d0/sqrt(1.d0-c2))*
     &            (2*(im-1)*c*plmrec(ill,im-1,c,c2,s,s2)) - 
     &        ((ill+im-1)*(ill-im+2))*plmrec(ill,im-2,c,c2,s,s2)
      elseif(ill.gt.2) then
       plmrec=(1.d0/real(ill-im))*
     . ((2*ill-1)*c*plmrec(ill-1,im,c,c2,s,s2) 
     &        - (ill+im-1)*plmrec(ill-2,im,c,c2,s,s2))
*-----------------------------------------------------------------------
      endif
      return
      end
*----------------------------------------------------------------------
      integer function ie(a,b) 
      real*8 a,b
      ie=1
      if(abs(a-b).lt.1.d-5) ie=0
      return
      end
*----------------------------------------------------------------------
      subroutine factor
      implicit real*8(a-h,o-z)
      common/clebma/faclog(500)
      faclog(1)=0.0d0
      faclog(2)=0.0d0
      fn=1.0d0
      do i=3,500
       fn=fn+1.0d0
       faclog(i)=faclog(i-1)+log(fn)
      enddo
      return
      end
*----------------------------------------------------------------------
ccc      real*8 function cleb(ria,rid,rib,rie,ric,rif)
ccc      implicit real*8(a-h,o-z)
ccc      common/clebma/faclog(500)
ccc      ia=2.d0*(ria+.0001d0)
ccc      ib=2.d0*(rib+.0001d0)
ccc      ic=2.d0*(ric+.0001d0)
ccc      id=int(sign(1.d0,rid)*2.d0*(abs(rid)+.0001d0))
ccc      ie=int(sign(1.d0,rie)*2.d0*(abs(rie)+.0001d0))
ccc      if=int(sign(1.d0,rif)*2.d0*(abs(rif)+.0001d0))
ccc      wwww=-1.0d0
ccc      cleb=0.0d0
ccc      if(id+ie-if) 7000,105,7000
ccc  105 k1=ia+ib+ic
ccc      if((-1)**k1) 7000,107,107
ccc  107 if(.not.((id.eq.0).and.(ie.eq.0))) go to 110
ccc      k1=k1/2
ccc      if((-1)**k1) 7000,110,110
ccc  110 k1=ia+ib-ic
ccc      k2=ic-iabs(ia-ib)
ccc      k3=min0(k1,k2)
ccc      if(k3) 7000,130,130
ccc  130 if((-1)**(ib+ie)) 7000,7000,140
ccc  140 if((-1)**(ic+if)) 7000,7000,150
ccc  150 if(ia-iabs (id)) 7000,152,152
ccc  152 if(ib-iabs (ie)) 7000,154,154
ccc  154 if(ic-iabs (if)) 7000,160,160
ccc  160 if(ia) 7000,175,165
ccc  165 if(ib) 7000,175,170
ccc  170 if(ic) 7000,180,250
ccc  175 cleb=1.0d0
ccc      go to 7000
ccc  180 fb=float(ib+1)
ccc      cleb=((wwww)**((ia-id)/2))/sqrt(fb)
ccc      go to 7000
ccc  250 fc2=ic+1
ccc      iabcp=(ia+ib+ic)/2+1
ccc      iabc=iabcp-ic
ccc      icab=iabcp-ib
ccc      ibca=iabcp-ia
ccc      iapd=(ia+id)/2+1
ccc      iamd=iapd-id
ccc      ibpe=(ib+ie)/2+1
ccc      ibme=ibpe-ie
ccc      icpf=(ic+if)/2+1
ccc      icmf=icpf-if
ccc      vvv=0.5d0
ccc      sqfclg=vvv*(log(fc2)-faclog(iabcp+1)
ccc     1      +faclog(iabc)+faclog(icab)+faclog(ibca)
ccc     2      +faclog(iapd)+faclog(iamd)+faclog(ibpe)
ccc     3      +faclog(ibme)+faclog(icpf)+faclog(icmf))
ccc      nzmic2=(ib-ic-id)/2
ccc      nzmic3=(ia-ic+ie)/2
ccc      nzmi= max0(0,nzmic2,nzmic3)+1
ccc      nzmx= min0(iabc,iamd,ibpe)
ccc      if(nzmx.lt.nzmi) go to 7000
ccc      s1=(wwww)**(nzmi-1)
ccc      do 400 nz=nzmi,nzmx
ccc      nzm1=nz-1
ccc      nzt1=iabc-nzm1
ccc      nzt2=iamd-nzm1
ccc      nzt3=ibpe-nzm1
ccc      nzt4=nz-nzmic2
ccc      nzt5=nz-nzmic3
ccc      termlg=sqfclg-faclog(nz)-faclog(nzt1)-faclog(nzt2)
ccc     1           -faclog(nzt3)-faclog(nzt4)-faclog(nzt5)
ccc      ssterm=s1*exp (termlg)
ccc      cleb=cleb+ssterm
ccc  400 s1=-s1
ccc 7000 return
ccc      end
*----------------------------------------------------------------------
      complex*16 function f2c(xb,yb,xyt,fxyt,nnx,nny,nord,mmx,mmy)
      implicit real*8 (a-h,o-z)
      complex*16 fxyt(mmx,mmy),fxy(33,33),cm1,cj3
      real*8 x(33),y(33),xyb(2),xyt(2,mmy)
      integer nxy(2),nbg(2)
      nnn=nord+1
      nxy(1)=nnx
      nxy(2)=nny
      xyb(1)=xb
      xyb(2)=yb
*     ---------------------------------------------------------------
*     loop on x and y grid finding sub-matrix needed
*     ---------------------------------------------------------------
      do i=1,2
       if(xyb(i).lt.xyt(i,1)) then
        num=1
       else if(xyb(i).gt.xyt(i,nxy(i))) then
        num=nxy(i)-nnn+1
       else
        min=1
        max=nxy(i)
        num=nxy(i)/2
 55     if(max-min.lt.4) goto 70
        if(xyt(i,num)-xyb(i)) 60,82,61
 60     min=num
        num=num+(max-min+1)/2
        goto 55
 61     max=num
        num=num-(max-min+1)/2
        goto 55
 70     num=max
 71     if(xyt(i,num)-xyb(i)) 82,82,81
 81     num=num-1
        goto 71
*     ---------------------------------------------------------------
*     change/correction July 2008
*     ---------------------------------------------------------------
 82     num=max0(1,num-nord/2)
        num=min0(num,nxy(i)-nord)
       endif
       nbg(i)=num-1
      enddo
*     ---------------------------------------------------------------
      nx=nbg(1)
      ny=nbg(2)
      do i=1,nnn
       x(i)=xyt(1,nx+i)-xb
       y(i)=xyt(2,ny+i)-yb
       do j=1,nnn
        fxy(i,j)=fxyt(nx+i,ny+j)
       enddo
      enddo
*     ---------------------------------------------------------------
      do ii=2,nnn
       iii=ii-1
       do jj=ii,nnn
        cj2=x(jj)-x(iii)
        cj3=fxy(iii,iii)*x(jj)-fxy(jj,iii)*x(iii)
        do mm=ii,nnn
        cm1=fxy(iii, mm)*x(jj)-fxy(jj, mm)*x(iii)
        fxy(jj,mm)=(cj3*y(mm)-cm1*y(iii))/(y(mm)-y(iii))/cj2 ! ****
        enddo
       enddo
      enddo
      f2c=fxy(nnn,nnn)
      return
      end
*----------------------------------------------------------------------
      subroutine sigma(lmax,eta,sigmad)                          
      implicit real*8(a-h,o-z)                                          
      dimension sigmad(0:10)                                            
      if(eta.ge.10.) go to 20                                           
      eta2=eta*eta                                                      
      t1=eta2+16.0                                                      
      t2=t1*t1                                                          
      t4=t2*t2                                                          
      t6=eta/6.0                                                        
      sigma0=-(eta/(12.0*t1))*(1.0+(eta2-48.0)/(30.0*t2)+(eta2**2-160.0*
     1eta2+1280.0)/(t4*105.0))-eta+3.0*t6*log(t1)+3.5*atan (1.5*t6)-at  
     2an (eta)-atan (3.0*t6)-atan (2.0*t6)                              
      go to 25                                                          
   20 t1=1.0/eta                                                        
      t2=t1*t1                                                          
      t3=t1*t2                                                          
      t5=t3*t2                                                          
      t7=t5*t2                                                          
      t9=t7*t2                                                          
      sigma0=0.7853981634d0+eta*log(eta)-eta                            
     1      -(0.08333333333d0*t1+0.00277777777d0*t3                     
     2       +0.00079365079d0*t5+0.00059523810d0*t7                     
     3       +0.00084175084d0*t9)                                       
   25 mn=sigma0/6.2831853072d0                                          
      sigma0=sigma0-6.2831853072d0*mn                                   
      sigmad(1)=sigma0                                                  
      do l=2,lmax+1                                                    
       fl1=l-1                                                         
       sigmad(l)=sigmad(l-1)+atan (eta/fl1)                            
      enddo 
      do l=0,lmax
       sigmad(l)=sigmad(l+1)
      enddo                                                         
      return                                                            
      end                                                               
