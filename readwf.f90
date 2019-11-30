      subroutine readwfs(filename, iset,nch) 
      use wfs, only: rvec,nr,wfc,energ
      use parameters, only: maxchan
      use channels, only: jpiset
      implicit none
      logical:: realwf
      integer:: ich,ir,i
      integer:: iset,nch,nfun,npt
      character*40 filename,line
      character*10 type
      real*8:: r
      real*8,parameter:: alpha=0.0
      real*8:: energy, x(maxchan), y(maxchan)
      real*8,allocatable:: rv(:)
      complex*16,allocatable:: ucwf(:,:)
      complex*16:: cfival ,caux

      open(20,file=filename)
      realwf=.true.
      read(20,'(a)') line
      read(line,*) nfun, npt, type
      if(type(1:1)=='c'.or.type(1:1)=='C') realwf=.false.

      write(*,'(5x,a,i3,a,a10,a,i3,a)') 'Reading', nfun,
     & ' wf(s) from file: ', filename ,
     & 'with ',nch, ' chan(s)'
      jpiset(iset)%nex=nfun
 
      allocate(rv(npt))
      allocate(ucwf(npt,nch))
  
      do i=1,nfun
      ucwf(:,:)=0.
      read(20,*) energy
      jpiset(iset)%exc(i)=energy
      energ(iset,i)=energy
 
      write(*,'(5x,a,1f8.4,a)') '-> Ext. function with E=',energy, 
     & ' MeV'
      do ir=1,npt
      if (realwf) then    ! real
         read(20,*) rv(ir),(x(ich),ich=1,nch)
!      if (nch.eq.1) then  
!         read(20,*) rv(ir), x(1) 
         ucwf(ir,1:nch)=cmplx(x(1:nch),0.0)
      else                ! complex
         read(20,*) rv(ir),(x(ich),y(ich),ich=1,nch)
         ucwf(ir,1:nch)=cmplx(x(1:nch),y(1:nch)) 
      endif
      enddo !ir

c interpolate in THOx grid    
      do ir=1,nr
      r=rvec(ir)
      if (r.lt.1e-4) r=1e-4
      if(r.gt.rv(npt)) cycle
      do ich=1,nch
      wfc(iset,i,ich,ir)=cfival(r,rv(:),ucwf(:,ich),npt,alpha)/r
      enddo !ich
!      write(99,*)r,real(r*wfc(iset,i,1,ir))
      enddo !ir
      deallocate(ucwf)
      enddo !i

      

      end subroutine

