c *** ----------------------------------------------
c     Build THO basis with analytical LST transformation 
c *** -----------------------------------------------
      subroutine thobasis(l,m,bosc,gam,nho)
        use wfs, only: wftho,rmin,rmax,dr,rvec,nr,rweight
        implicit none
        logical:: debug=.false.
        integer:: l,nho,n,ir,m,kbout,nmin,ibas
        real*8:: r,wf,wfho,nu,bosc,gam,eps
!        real*8:: rstart,rmax,dr
        real*8:: sr,dsr,lst,dlst
        real*8:: norm
	real*8:: test
        parameter(eps=1e-6)
        nu=1d0/bosc**2d0
	kbout=10


!        open(50,file='lst.out',status='unknown')                !uncoment to check local scale transformation

!        write(*,*)'Building THO wfs and writing in fort.10'
!        write(*,*)'m=',m
c        if (dr>eps) then
c           nr=(rmax-rmin)/dr+1
c        else
c           write(0,*)'hobasis:dr too small!. Aborting'; stop
c        endif
c	write(*,*)'nr=',nr
	test=0.

! TEST !!!!!!!!!!!!!!!!!!!!!!
        nmin=0
        ibas=0
        do n=nmin,nho+nmin-1
           ibas=ibas+1
!        do n=0,nho-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           norm=0
           do ir=1,nr
!              r=rmin+(ir-1)*dr  !old definition, we now use rvec
	      r=rvec(ir)
!	      write(*,*)ir,rvec(ir)
              if (r<eps) r=eps
              sr=lst(m,gam,r)
!              if (debug)write(50,*) r,sr                     !uncoment to check local scale transformation
              dsr=dlst(m,gam,r)  
              call ho3d(nu,n, l, sr,wfho)
            
              wf=wfho *sr*sqrt(dsr)/r
              wftho(ibas,l,ir)=wf
              if(debug) write(kbout,*) r,r*wf !,sr,dsr
	      if (.not.allocated(rweight)) then
              norm=norm+dr*(wf*r)**2
	      else
	      norm=norm+rweight(ir)*(wf*r)**2
	      endif
           enddo
           if (debug)write(kbout,*)'#Norm=',norm
           if (debug)write(kbout,*)'&'
	      test=test+(norm-1)**2
!        write(*,*)'l,n,norm=',l,n,norm
        enddo
        close(50)
	   
      end subroutine thobasis





c *** -------------------------------------------
c  LST proposed by Karatagladis, Amos and Giraud,
c  PRC71, 064601 (2005)
c *** -------------------------------------------
      real*8 function lst(m,g,r)
        implicit none
        integer m
        real*8:: g,r,one,eps,xm
        parameter(eps=1e-6)
        one=1d0
	xm=m
        if (r<eps) r=eps
        lst=((1d0/r)**xm + (1d0/g/r**0.5d0)**xm)**(-1d0/xm)
      end function lst



c *** -----------------------------------------------
c  Derivative of the LST of Karatagladis et al
c *** -----------------------------------------------
       real*8 function dlst(m,g,r)
         implicit none
         integer:: m
         real*8:: g,r,c1,c2,eps,xm
         parameter(eps=1e-6)
	 xm=m

         if (r<eps) r=eps

         c1=1d0/((1d0/r)**xm + (1d0/g/dsqrt(r))**xm)**(1d0/xm+1d0)
         c2=(1d0/r)**(xm+1d0)+(1d0/2d0/g**xm)*(1/dsqrt(r))**(xm+2d0)
         dlst=c1*c2
       end function dlst

!---------------------------------------------------------------
!       CTHO basis by Lay
!---------------------------------------------------------------

       subroutine cthobasis(l,m,bosc,gam,eta,nho)
        use wfs, only: wftho,rmin,rmax,dr,rvec,nr,rweight
       implicit none
       integer:: l,nho,n,ir,m,kbout,ibas,nmin
       real*8:: r,wf,wfho,nu,bosc,gam,eps,eta
       real*8:: rstart,sr,dsr,clst,dclst
       real*8:: norm,test
       parameter(eps=1e-6)
       nu=1d0/bosc**2d0
       kbout=10
       rstart=rmin
       open(50,file='clst.out',status='unknown')
       open(51,file='dclst.out',status='unknown')
       !        write(*,*)'Building THO wfs and writing in fort.10'
       !        write(*,*)'m=',m
       if (dr>eps) then
       nr=(rmax-rstart)/dr+1
       else
       write(0,*)'hobasis:dr too small!. Aborting'; stop
       endif



        nmin=0
        ibas=0
        do n=nmin,nho+nmin-1
           ibas=ibas+1
!        do n=0,nho-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           norm=0
           do ir=1,nr
!              r=rmin+(ir-1)*dr  !old definition, we now use rvec
	      r=rvec(ir)
!	      write(*,*)ir,rvec(ir)
              if (r<eps) r=eps
              sr=clst(m,bosc,gam,eta,r)
!              if (n.eq.0)write(50,*) r,sr                     !uncoment to check local scale transformation
              dsr=dclst(m,bosc,gam,eta,r)
              call ho3d(nu,n, l, sr,wfho)
            
              wf=wfho *sr*sqrt(dsr)/r
              wftho(ibas,l,ir)=wf
!              write(kbout,*) r,r*wf !,sr,dsr
	      if (.not.allocated(rweight)) then
              norm=norm+dr*(wf*r)**2
	      else
	      norm=norm+rweight(ir)*(wf*r)**2
	      endif
           enddo
!           write(kbout,*)'#Norm=',norm
!           write(kbout,*)'&'
	      test=test+(norm-1)**2
!        write(*,*)'l,n,norm=',l,n,norm
        enddo
        close(50)

       end subroutine cthobasis




c *** -------------------------------------------
c  LST proposed by Karatagladis, Amos and Giraud,
c  PRC71, 064601 (2005)  +   Charged addition
c *** -------------------------------------------
       real*8 function clst(m,b,g,e,r)
       implicit none
       integer m
       real*8:: b,g,e,r,one,eps,xm,k,c1
       parameter(eps=1e-6)
       one=1.d0
       xm=m
       if (r<eps) r=eps

                   k=g*g/b/b/2

       c1=g*dsqrt(r+e/k*dlog(1.d0+2.d0*k*r))

         clst=((one/r)**xm+(one/c1)**xm)**(-one/xm)

!       clst=((1d0/r)**xm +(1d0/g/
!     &dsqrt(r+2*e*(b/g)**2d0*dlog((g/b)**2d0*r+1d0))**xm))**(-1d0/xm)
       end function clst

c *** -----------------------------------------------
c  Derivative of the LST of Karatagladis et al +   Charged addition
c *** -----------------------------------------------
       real*8 function dclst(m,b,g,e,r)
       implicit none
       integer:: m
       real*8:: b,g,e,r,c1,c2,c3,eps,xm,k
       parameter(eps=1e-6)
       xm=m

       if (r<eps) r=eps

         k=g*g/b/b/2

       c1=g*dsqrt(r+e/k*dlog(1.d0+2.d0*k*r))
       c2=g*g/2.*(1.+2.*e/(1+2*k*r))
       c3=(r**(-xm)+c1**(-xm))**(-(xm+1)/xm)
       dclst=(r**(-1.-xm)+c2*c1**(-2.-xm))*c3
       end function dclst



