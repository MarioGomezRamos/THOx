*     ---------------------------------------------------------------      
      program coreval2
      implicit real*8 (a-h,k-z)
*     JAT May 1998 
*     ---------------------------------------------------------------
      character outfile*30,infile*30,sys*3,ansc,ansv,anse,ansa
*     array for cross sections/indices/energy:val theta:core theta
      integer iemax,ivmax,icmax,ntmp,ipmax   
      parameter(iemax=621,ivmax=200,icmax=iemax,ntmp=iemax)
      parameter (ipmax=300)
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

      
      end program
