
! ------------------------------------------------------------ 
      orto =.false.
      if (orto.and.(r.lt.rturn).and.(mod(ir,5).eq.0)) then 
      a(:,:)=0
      do i=1,nch
      a(i,i)=1
      do j=i+1,nch
        c1=0; c2=0;
        do k=1,nch
        c1=c1+conjg(yp(k,i))*yp(k,j)
        enddo !k
        do k=1,nch
        c2=c2+abs(yp(k,i))**2
        enddo !k
        a(i,j)=-c1/c2
      enddo !j
      enddo !i

      if (debug) then
      write(95,*)'Ortog. matrix at ir=',ir
      do ich=1,nch
        write(95,'(5x,50g14.5)') 
     &  (a(ich,is),is=1,nch)
      enddo
      endif 


      do i=1,ir+1
      y(:,:,i)=matmul(y(:,:,i),a(:,:))
      enddo ! i
      yp=y(:,:,ip)
      y0=y(:,:,i0)
      
      endif ! orto
! ------------------------------------------------------------


c
c orthogonalize solution
c 
      orto =.false.

      if (orto.and.(r.lt.rturn)) then    

        if (debug.and.(mod(ir,10).eq.0).or.(ir.eq.irmin)) then
        write(95,*)'y(is).y(is-1) (before G-S) for ir',ir
        do is=2,nch
        write(95,'(5x,"is=",i2,3x,50g14.5)') is,
     &     dot_product(yp(:,is),yp(:,is-1))
        enddo
        endif !debug


!       y(:,1,ir+1)= yp(:,1) ! leave first column unchanged
       do is=2,nch
       yaux=0
       do isp=1,is-1
       yaux(:)=yaux(:) + yp(:,isp)
     &        *dot_product(yp(:,isp),yp(:,is))/
     &         dot_product(yp(:,isp),yp(:,isp))
       enddo !isp


      yp(:,is)=yp(:,is) - yaux(:)
      enddo !is

      call zfun(nch,ql,kch2,ip,yp,zp)
      endif

