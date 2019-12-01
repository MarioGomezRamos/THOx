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
      real*8 temp(ntmp) ! sig(iemax,ivmax,icmax,ipmax)
      real*8, allocatable:: sig(:,:,:,:)
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
  
      write(*,*)'Dims:',itc,itv,iten,itp
      allocate(sig(itc,itv,iten,itp))
  
*     --------------------------------------------------------------- 
*     read back ALL cross section information and close file
*     --------------------------------------------------------------- 
      xstot=0
      icount=0
      do 10 ic=1,itc
      tc=tcl+(ic-1)*dtc
      tcr=tc*pi/180.
      do 10 iv=1,itv
      tv= tvl + (iv-1)*dtv
      tvr=tv*pi/180.
      do 10 ien=1,iten
      ed= enlow + (ien-1)*den
      
      do 10 ip=1,itp
      phi=phl + (ip-1)*dph
      phir=phi*pi/180.
      icount=icount + 1 
      read(77,*) i1,i2,i3,i4,raux
!      write(*,'(5i7,1g14.6)') icount,ien,iv,ic,ip,raux
      sig(ic,iv,ien,ip)=raux
      xstot=xstot + raux*dtcr*sin(tcr)*dtvr*sin(tvr)*den*dphr
   10 continue
      xstot=xstot*2.0
      close(77)
      
      write(*,*)'Total x-section=',xstot
      
*
*     dsigma/dphi
*
      xstot=0 
      do ip=1,itp
      phi=phl + (ip-1)*dph
      xstot=0
      do ic=1,itc
      tc=tcl+(ic-1)*dtc
      tcr=tc*pi/180.
      do iv=1,itv
      tv= tvl + (iv-1)*dtv
      tvr=tv*pi/180.
      do ien=1,iten
      ed= enlow + (ien-1)*den
      raux=sig(ic,iv,ien,ip)

      xstot =xstot + 2.*raux*dtcr*sin(tcr)*dtvr*sin(tvr)*den
!            write(*,*)ip,ic,iv,ien,raux,xstot
!      xstot=xstot + 2.*raux*dtcr*sin(tcr)*dtvr*sin(tvr)*den*dphr

      enddo ! ien
      enddo ! iv
      enddo ! ic
      write(98,*) phi, xstot
      enddo  ! ip (phi)
      
!      write(*,*)'Total x-section=',xstot
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
