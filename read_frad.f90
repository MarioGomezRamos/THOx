!      subroutine read_frad(l,kin,nfun,emin,enerm,energ,np)
      subroutine read_frad(l,kin,nfun,emin,enerm)
      use xcdcc,only:nch,nchmax,nrad1,nquad,rin,dr,rquad,frad,
     .xrad1,nener,numfilmax,
     .energ,np
      implicit real*8 (a-h,o-z)
      real*8,pointer:: yvp(:)
      real*8,allocatable,target:: yv(:,:)
!      dimension energ(nener,numfilmax),np(numfilmax)
      allocate(yv(nrad1,nchmax),xrad1(nrad1))
      do irad1=1,nrad1
      xrad1(irad1)=rin+dr*(irad1-1)
      enddo
!      print*,'read_frad: starting'
      i=0
      alpha=0.d0
      do n=1,nfun
      read(kin,*,end=9) ener
      if((ener.gt.enerm)) goto 9
      if (ener.lt.emin) then
         write(*,'(8x,"skipping energy Ex=",1f8.4)') ener  
         do j=1,nrad1
         read(kin,*) (yv(j,k),k=1,nch(l))
         enddo
         cycle
      endif
      i=i+1
      energ(i,l)=ener
      write(*,'(8x,"State",i3,2x,"Ex=",1f8.4)') i,energ(i,l)
      do j=1,nrad1
      read(kin,*) (yv(j,k),k=1,nch(l))
      enddo
      do k=1,nch(l)
      yvp=>yv(:,k)
      do j=1,nquad
      frad(l,i,j,k)=fival(rquad(j),xrad1,yvp,nrad1,alpha)
      enddo
      enddo
      enddo
9     continue !new function
      deallocate(yv,xrad1)
!      print*,'read_frad: finished'
      np(l)=i
      if (i.gt.0) write(*,'(/)')
      return
      end
