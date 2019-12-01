*     ---------------------------------------------------------------      
      program coreval2
      implicit real*8 (a-h,k-z)
*     JAT May 1998 
*     ---------------------------------------------------------------
      character outfile*30,infile*30,sys*3,ansc,ansv,anse,ansa
*     array for cross sections/indices/energy:val theta:core theta
      integer iemax,ivmax,icmax,ntmp,ipmax   
      parameter(iemax=821,ivmax=400,icmax=iemax,ntmp=iemax)
      parameter (ipmax=600)
      real*8 temp(ntmp),sig(iemax,ivmax,icmax,ipmax)
      real*8 raux, xstot
      logical apt
*     constants
      hbarc=197.34d0
      amu=938.0d0               
      pi=4.d0*datan(1.d0)
*     --------------------------------------------------------------- 
 117  print*,' input file for cross sections:  '
      read'(a)',infile
 119  open(77,file=infile,status='unknown')
      rewind 77
*     --------------------------------------------------------------- 
*     read back essential run information from file and print 
*     --------------------------------------------------------------- 
      read(77,'(a)') sys
      print*,' coordinate system is set to: ',sys
      read(77,304) mt,zt,eoa
  304 format(5f10.3)
      print*,' target mass: charge: incident energy/nucleon '
      print 304,mt,zt,eoa
      read(77,304) mc
      print*,' core particle mass  '
      print 304,mc
      read(77,304) mv,ebind
      print*,' valence particle mass and sep energy '
      print 304,mv,ebind
      print 303,' projectile of mass ',mc+mv
  303 format(1x,a,f9.4,a,f9.4)
      read(77,304) ep
      print 303,' projectile lab energy = ',ep
      read(77,*) idet
      if(idet.eq.1) then
      print*,' idet=1: core particle energy is specified '
      else
      print*,' idet=2: valence particle energy specified ' 
      endif
*     --------------------------------------------------------------- 
*     read back stored thetas, energy and phi information in file 
*     --------------------------------------------------------------- 
      read(77,307) itc,dtcr,tcl,tcu,dtc
      read(77,307) itv,dtvr,tvl,tvu,dtv
      read(77,307) iten,den,enlow,enup 
      read(77,307) itp,dphr,phl,phu,dph
  307 format(i5,f15.8,3f10.3)
*     --------------------------------------------------------------- 
      if(itc.gt.icmax) then
       print*,'  itc is ',itc
       print*,'  icmax is ',icmax,' so too small '
       stop 
      endif
      if(itv.gt.ivmax) then
       print*,'  itv is ',itv
       print*,'  ivmax is ',ivmax,' so too small '
       stop 
      endif
      if(iten.gt.iemax) then
       print*,'  iten is ',iten
       print*,'  iemax is ',iemax,' so too small '
       stop 
      endif
*     --------------------------------------------------------------- 
*     read back ALL cross section information and close file
*     --------------------------------------------------------------- 
      xstot=0
      do 10 ic=1,itc
      do 10 iv=1,itv
      do 10 ien=1,iten
      do 10 ip=1,itp
      read(77,*) itcl,itcl,itcl,itcl,raux
      sig(ien,iv,ic,ip)=raux
      xstot=xstot + raux*dtcr*sin(dtcr)*dtvr*sin(dtvr)*den*dphr
   10 continue
      close(77)
      
      write(*,*)'Total x-section=',xstot
      
*     --------------------------------------------------------------- 
*     check for any spike in energy distributions
*     --------------------------------------------------------------- 
c      do ic=1,itc
c      do iv=1,itv
c       do ien=3,iten-2
c        grad1=abs(sig(ien+1,iv,ic)-sig(ien-1,iv,ic))/2.d0
c        grad=(sig(ien,iv,ic)-sig(ien-1,iv,ic))
c        grad=abs(grad-(sig(ien+1,iv,ic)-sig(ien,iv,ic)))
c        if(grad1.gt.0.d0) then
c         grad=grad/grad1
c        else
c         grad=0.d0
c        endif
c*       write(25,*) ien,grad
c        if(grad.gt.5.d0) then
c         sigav=2.d0*(sig(ien+1,iv,ic)+sig(ien-1,iv,ic))/3.d0
c         sigav=sigav-(sig(ien+2,iv,ic)+sig(ien-2,iv,ic))/6.d0
c         sig(ien,iv,ic)=sigav
c        endif
c       enddo
c      enddo
c      enddo


! AMM; check total integrated cross section

      








*     ---------------------------------------------------------------    
*     what to do with these now depends on what is stored ?
*     first notify regarding the phi status
*     ---------------------------------------------------------------   
      print*
      print*,' ================================================='
      if(itp.eq.1) then
      print 330,' file contains sigmas for a single phi =',phl
      print*,' SO NO INTEGRATIONS OVER SOLID ANGLES POSSIBLE'
  330 format(1x,a,f10.4)
      ansc='n'
      ansv='n'
      else      
      print*,' file contains sigmas already integrated over phi '
      print 333,phl,phu,dph,' degrees '
      print*,' MUST INTEGRATE OVER AT LEAST ONE SOLID ANGLE '
  333 format(3f10.4,a)
      ansc='y'
      ansv='y'
      endif
      print*,' ================================================='
      print*
*     ---------------------------------------------------------------   
*     sort out core (charged) particle angles 
*     ---------------------------------------------------------------   
      itcl=1
      itcu=itc
      if(itc.gt.1) then
       print 332,' file contains sigmas for',itc,' core thetas'
       print 333,tcl,tcu,dtc,' degrees '
  332  format(1x,a,i4,a)
       print*,' use full angular range above? (y/n?) '
       read'(a)',ansa
       if(ansa.eq.'n') then
  111   print*,' input lower and upper theta_(core) limits: '
        read*,tcln,tcun
        if(tcln.lt.tcl.or.tcun.gt.tcu) then
         print*,' these are not in stated range '
         go to 111
        endif
        itcl=nint((tcln-tcl)/dtc)+1
        itcu=itcl+nint((tcun-tcln)/dtc)
        if((-1)**(itcu-itcl).lt.0) then
         print*,' not even number of intervals so problem if '
         print*,' plan to integrate (uses Simpson): reenter (y/n?) '
         read'(a)',ansa
         if(ansa.eq.'y') goto 111
        endif
        tcl=tcln
        tcu=tcun
       endif
      else
       print*,' file contains sigmas for single theta_(core)'
       ansc='n'
      endif
*     ---------------------------------------------------------------   
      if(ansc.eq.'y') then
      print*,' integrate on dOmega_(core) (y/n?)'
      read'(a)',ansc
       if(ansc.eq.'y'.and.(-1)**(itcu-itcl).lt.0) then
       print*,' need even number of intervals for Simpsons rule '
       print*,' sorry, you were warned !!'
       print*
       goto 117
       endif
      endif
*     ---------------------------------------------------------------   
      print*,' so sigma for theta_(core) in range '
      print 333,tcl,tcu,dtc,' degrees '
      if(ansc.ne.'y') print*,' without solid angle integration '
      if(ansc.eq.'y') print*,' with solid angle integration '
      print*,' ============================================='
*     ---------------------------------------------------------------   
*     sort out valence particle angles 
*     ---------------------------------------------------------------   
  888 itvl=1
      itvu=itv
      if(itv.gt.1) then
       print 332,' file contains sigmas for',itv,' val thetas'
       print 333,tvl,tvu,dtv,' degrees '
       print*,' use full angular range above? (y/n?) '
       read'(a)',ansa
       if(ansa.eq.'n') then
  112   print*,' input lower and upper theta_(val) limits: '
        read*,tvln,tvun
        if (tvln.lt.tvl.or.tvun.gt.tvu) then
         print*,' these are not in stated range '
         go to 112
        endif
        itvl=nint((tvln-tvl)/dtv)+1
        itvu=itvl+nint((tvun-tvln)/dtv)
        if((-1)**(itvu-itvl).lt.0) then
         print*,' not even number of intervals so problem if '
         print*,' plan to integrate (uses Simpson): reenter (y/n?) '
         read'(a)',ansa
         if(ansa.eq.'y') goto 112
        endif
        tvl=tvln
        tvu=tvun
       endif
      else
       print*,' file contains sigmas for single theta_(val)'
       ansv='n'
      endif
*     ---------------------------------------------------------------   
      if(ansv.eq.'y') then
      print*,' integrate on dOmega_(val) (y/n?)'
      read'(a)',ansv
       if(ansv.eq.'y'.and.(-1)**(itvu-itvl).lt.0) then
       print*,' need even number of intervals for Simpsons rule '
       print*,' sorry, you were warned !!'
       print*
       goto 117
       endif
      endif
*     ---------------------------------------------------------------   
      print*,' so sigma for theta_(val) in range '
      print 333,tvl,tvu,dtv,' degrees '
      if(ansv.ne.'y') print*,' without solid angle integration '
      if(ansv.eq.'y') print*,' with solid angle integration '
      print*,' ============================================='
*     ---------------------------------------------------------------   
      if(itp.gt.1.and.ansc.eq.'n'.and.ansv.eq.'n') then
      print*,' ============================================='
      print*,' further calculation does not make sense !! '
      print*,' sigmas input are phi integrated and so should '
      print*,' be integrated over at least one domega '
      print*,' ============================================='
      stop 99
      endif
*     ---------------------------------------------------------------   
*     sort out energy range and integration
*     --------------------------------------------------------------- 
      anse='y'  
      itel=1
      iteu=iten
      if(iten.gt.1) then
       print 332,' file contains sigmas for',iten,' energies'
       print 333,enlow,enup,den,' MeV '
       print*,' use full energy range above? (y/n?) '
       read'(a)',ansa
       if(ansa.ne.'y') then
  113   print*,' lower and upper energy limits : '
        read*,enlown,enupn
        if (enlown.lt.enlow.or.enupn.gt.enup) then
         print*,' these are not in stated range '
         go to 113
        endif
        itel=nint((enlown-enlow)/den)+1
        iteu=itel+nint((enupn-enlown)/den)
        if((-1)**(iteu-itel).lt.0) then
         print*,' not even number of intervals so problem if '
         print*,' plan to integrate (uses Simpson): reenter (y/n?) '
         read'(a)',ansa
         if(ansa.eq.'y') goto 113
        endif
        enlow=enlown
        enup=enupn
       endif 
      else
       print*,' file contains sigmas for single energy'
       anse='n'
      endif
*     --------------------------------------------------------------- 
      if(anse.eq.'y') then
      print*,' integrate on energy (y/n?)'
      read'(a)',anse
       if(anse.eq.'y'.and.(-1)**(iteu-itel).lt.0) then
       print*,' need even number of intervals for Simpsons rule '
       print*,' sorry, you were warned !!'
       print*
       goto 117
       endif
      endif
*     --------------------------------------------------------------- 
      print*,' so sigma for energies in range '
      print 333,enlow,enup,den,' MeV '
      if(anse.ne.'y') print*,' without energy integration '
      if(anse.eq.'y') print*,' with energy integration '
      print*,' ============================================='
*     ---------------------------------------------------------------
*     if anse=y integrate over the selected energies
*     ---------------------------------------------------------------
      if(anse.eq.'y') then
       do 20 ic=itcl,itcu
       do 20 iv=itvl,itvu
       do 30 ien=itel,iteu
        temp(ien)=sig(ien,iv,ic)
   30  continue
        call sim(temp,res,itel,iteu,den)
        sig(1,iv,ic)=res
   20  continue
       itel=1
       iteu=1
      endif
*     ---------------------------------------------------------------
*     if ansc=y integrate over the selected theta_(core)
*     ---------------------------------------------------------------
      if(ansc.eq.'y') then
!       write(0,*)'ok3: iteu,itvu,itcu=',iteu,itvu,itcu
       do 40 ien=itel,iteu
       do 40 iv=itvl,itvu
        do 50 ic=itcl,itcu
        tcd=tcl+(ic-itcl)*dtc
        temp(ic)=sig(ien,iv,ic)*sin(tcd*pi/180.d0)
   50   continue
        call sim(temp,res,itcl,itcu,dtcr)
        sig(ien,iv,1)=res
   40  continue
       itcl=1
       itcu=1
      endif
      
*     ---------------------------------------------------------------
*     if ansv=y integrate over the selected theta_(val)
*     ---------------------------------------------------------------
      if(ansv.eq.'y') then
       do 60 ien=itel,iteu
       do 60 ic=itcl,itcu
        do 70 iv=itvl,itvu
        tvd=tvl+(iv-itvl)*dtv
        temp(iv)=sig(ien,iv,ic)*sin(tvd*pi/180.d0)
  70    continue
        call sim(temp,res,itvl,itvu,dtvr)
        sig(ien,1,ic)=res
  60   continue
       itvl=1
       itvu=1
      endif
*     ---------------------------------------------------------------
*     if ansv=y and ansc=y then extra factor of 2*pi for integration
*     over the second phi angle (cylindrical symmetry)
*     ---------------------------------------------------------------
      if(ansv.eq.'y'.and.ansc.eq.'y') then
       do 80 ien=itel,iteu
       sig(ien,1,1)=sig(ien,1,1)*2.d0*pi
   80  continue
      endif
  331 format(f12.4,d18.7)
*     ---------------------------------------------------------------
*     arrange printout now (to external file if required)
*     ---------------------------------------------------------------
      apt=(1.eq.0)
      print*,' write output to external file (y/n?)  '
      read'(a)',ansa
      if(ansa.eq.'y') then
      print*,' output filename:  '
      read'(a)',outfile
      apt=(1.eq.1)
      open(88,file=outfile,status='unknown')
      rewind 88
      endif
*     ---------------------------------------------------------------
*     if both core and valence particle solid angles integrated
*     ---------------------------------------------------------------
      if(ansc.eq.'y'.and.ansv.eq.'y') then
       if(anse.eq.'y') then
        print*, 'sigma = ',sig(1,1,1),' (mb)'
        if(apt) write(88,*) 'sigma = ',sig(1,1,1),' (mb)'
       else
        print*,'       E        dsigma/dE (mb/MeV) '
        print*
        do 90 ien=itel,iteu
! AMoro
!        print 331,enlow+(ien-itel)*den,sig(ien,1,1)
!        print 331,enlow+(ien-itel)*den,sig(ien,1,1)
        if(apt) write(88,331) enlow+(ien-itel)*den,sig(ien,1,1)
   90   continue
       endif
      endif
*     ---------------------------------------------------------------
*     if only core particle solid angle integrated
*     ---------------------------------------------------------------
      if(ansc.eq.'y'.and.ansv.ne.'y') then    
       if(anse.eq.'y') then
        print*,'    theta_v     dsigma/dOmega_v  '
        print*
        do 100 iv=itvl,itvu
        print 331,tvl+(iv-itvl)*dtv,sig(1,iv,1)
        if(apt) write(88,331) tvl+(iv-itvl)*dtv,sig(1,iv,1)
  100   continue
       else
        do 110 iv=itvl,itvu
        print*,'theta_v = ',tvl+(iv-itvl)*dtv
        if(apt) write(88,*) 'theta_v = ',tvl+(iv-itvl)*dtv
        print*,'       E        dsigma/dEdOmega_v (mb/MeV/sr) '
        print*
        do 120 ien=itel,iteu
        print 331,enlow+(ien-itel)*den,sig(ien,iv,1)
        if(apt) write(88,331) enlow+(ien-itel)*den,sig(ien,iv,1)
  120   continue
        print*
  110   continue
       endif
      endif
*     ---------------------------------------------------------------
*     if only valence particle solid angle integrated
*     ---------------------------------------------------------------
      if(ansc.ne.'y'.and.ansv.eq.'y') then       
       if(anse.eq.'y') then
        print*,'  #  theta_c     dsigma/dOmega_c  '
        print*
        do 220 ic=itcl,itcu
        print 331,tcl+(ic-itcl)*dtc,sig(1,1,ic)
        if(apt) write(88,331) tcl+(ic-itcl)*dtc,sig(1,1,ic)
  220   continue
       else
        write(88,'(a,i3,a,i4)')'nel',iteu-itel+1,' nang',itcu-itcl+1
        do 230 ic=itcl,itcu
        print*,'# theta_c = ',tcl+(ic-itcl)*dtc
        if(apt) write(88,*) '# theta_c = ',tcl+(ic-itcl)*dtc
        print*,'       E        dsigma/dEdOmega_c (mb/MeV/sr) '
        print*
        do 240 ien=itel,iteu
*     ---------------------------------------------------------------
*     if momentum distribution (test lines)
       enn=enlow+(ien-itel)*den
       mval=sqrt(2.d0*mc*amu*enn)
       mcon=sqrt(2.d0*enn/mc/amu)
       print 331,mval,mcon*sig(ien,1,ic)
!       if(apt) write(88,331) mval,mcon*sig(ien,1,ic)
       write(89,331) mval,mcon*sig(ien,1,ic)
*     ---------------------------------------------------------------
        print 331,enlow+(ien-itel)*den,sig(ien,1,ic)
        if(apt) write(88,331) enlow+(ien-itel)*den,sig(ien,1,ic)
        write(99,'(1x,f8.3,2x,f12.6)') tcl+(ic-itcl)*dtc,sig(ien,1,ic)
*     ---------------------------------------------------------------
  240   continue
        print*
  230   continue
       endif
      endif
*     ---------------------------------------------------------------
*     if no solid integrations carried out and so one phi angle
*     ---------------------------------------------------------------
      if(ansc.ne.'y'.and.ansv.ne.'y'.and.itp.eq.1) then       
*     ---------------------------------------------------------------
*     if energy variable is integrated
*     ---------------------------------------------------------------
       if(anse.eq.'y') then
*     ---------------------------------------------------------------
*      if several valence particle thetas
*     ---------------------------------------------------------------
        if (itv.gt.1) then
         do 430 ic=itcl,itcu
          print*,'# theta_c = ',tcl+(ic-itcl)*dtc
         if(apt) write(88,*) '# theta_c = ',tcl+(ic-itcl)*dtc
         print*,'   theta_v    dsigma/dOmega_vdOmega_c (mb/sr/sr) '
         print*
         do 440 iv=itvl,itvu
         print 331,tvl+(iv-itvl)*dtv,sig(1,iv,ic)
         if(apt) write(88,331) tvl+(iv-itvl)*dtv,sig(1,iv,ic)
  440    continue
         print*
  430    continue
        else
*     ---------------------------------------------------------------
*      if one valence particle theta
*     ---------------------------------------------------------------
         do 489 ic=itcl,itcu
          print 331,tcl+(ic-itcl)*dtc,sig(1,1,ic)
          if(apt) write(88,331) tcl+(ic-itcl)*dtc,sig(1,1,ic)
  489    continue
        endif
*     ---------------------------------------------------------------
*     if energy variable not integrated and one valence theta
*     ---------------------------------------------------------------
       else if(itv.eq.1) then
        do 438 ic=itcl,itcu
        print*,'# theta_c = ',tcl+(ic-itcl)*dtc
        if(apt) write(88,*) '# theta_c = ',tcl+(ic-itcl)*dtc
        print*,'   E          dsigma/dE dOmega_vdOmega_c  '
        print*
        do 247 ien=itel,iteu
        print 331,enlow+(ien-itel)*den,sig(ien,1,ic)
        if(apt) write(88,331) enlow+(ien-itel)*den,sig(ien,1,ic)
  247   continue
        print*
  438   continue
       endif
      endif
*     ---------------------------------------------------------------
*     calculation complete or calculate other observables?
*     ---------------------------------------------------------------
      if(apt) close(88)
      print*
      print*,' finished (y/n?)  '
      read'(a)',ansa
! Changed by AMoro
!      if(ansa.ne.'y') goto 117
      if(ansa.ne.'y') goto 119
      stop
      end
*     ---------------------------------------------------------------
      subroutine sim(fa,res,m,n,h)
      implicit real*8(a-h,o-z)
      dimension fa(*),dq(801)
      do 90 i=m,n
      dq(i)=fa(i)
   90 continue
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
