c *** ----------------------------------------------
c     Build HO basis
c *** -----------------------------------------------

      subroutine hobasis(l,bosc,nho)
        use wfs, only: wftho,rmin,rmax,dr,rvec,nr
        implicit none
        integer:: l,nho,n,ir,kbout,ibas,nmin
        real*8:: r,wf,nu,bosc,eps
!        real*8:: rmin,rmax,dr
        real*8:: norm
        parameter(eps=1e-5)
        nu=1d0/bosc**2d0
	kbout=10
!        write(*,'("Building HO wfs and writing in fort.",i2)'),kbout
c        if (dr>eps) then
c           nr=(rmax-rmin)/dr+1
c        else
c           write(0,*)'hobasis:dr too small!. Aborting'; stop
c        endif

! TEST        do n=0,nho-1
        nmin=0
        ibas=0
        do n=nmin,nho+nmin-1
           ibas=ibas+1
           norm=0
           do ir=1,nr
c              r=rmin+(ir-1)*dr !old definition of r, now we use rvec
	      r=rvec(ir)        !new definition
              call ho3d(nu,n, l, r,wf)
              wftho(ibas,l,ir)=wf
!              write(kbout,*) r,r*wf
              norm=norm+dr*(wf*r)**2
           enddo
           write(kbout,*)'#Norm=',norm
           write(kbout,*)'&'
        enddo
      end subroutine hobasis



c *** -------------------------------------------------
c     Calculates HO wavefunctions (radial part) in 3D
c     nu=(mu*w/hbar)=1/bosc^2, where bosc is the oscillator 
c     parameter (typically in fm)
c *** ------------------------------------------------
      subroutine ho3d(nu,n, l, r,wf)
        implicit none
        integer l,n,p
        real*8:: x,r,tp,fact,norma,normho,wf,nu,lag
        real*8:: laguerre

        x=nu*r**2
        lag=laguerre(n,l,x)
        norma=normho(nu,n,l)
!        if (norma<1e-6) then
!           write(*,*)'ho3d: Norm=0!!!for  nu,n,l',nu,n,l
!        endif
        wf=norma*x**(l/2d0)*dexp(-0.5*x)*lag
      end subroutine ho3d
 

c *** ------------------------------------------------
c     Generalized Laguerre function L(n,l+1/2,x)
c *** -----------------------------------------------
      function laguerre(n,l,x)
        implicit none
        integer n,l,p
        real*8:: x,eps,aux,laguerre,fact,aux1,aux2,alpha
        parameter(eps=1e-6)

        alpha=l+0.5d0
        aux1=1d0
        aux2=-x+alpha+1d0


        if (n.eq.0) then
           aux=1d0
        else if (n.eq.1) then
           aux=aux2
        else
           do p=2,n
              aux=((2*p-1+alpha-x)*aux2-(p-1+alpha)*aux1)/p
              aux1=aux2
              aux2=aux
           enddo
        endif
         laguerre=aux


      end function laguerre


c *** ------------------------------------------------
c     Generalized Laguerre function L(n,l+1/2,x)
c *** -----------------------------------------------
      function laguerre2(n,l,x)
        implicit none
        integer n,l,p
        real*8:: x,eps,tp,laguerre2,fact
        parameter(eps=1e-6)

        if(x<eps) x=eps
        tp=(-x)**n/fact(n)

        do p=1,n
!           write(*,*)'p=',p
           tp=tp-(n+l+1.5-p)*(n+1-p)/p/x*tp
        enddo
        laguerre2=tp
      end function laguerre2


c *** --------------------------------------------------------
c     Calculates normalization coefficient of HO wavefunctions: N(n,l) 
c *** --------------------------------------------------------      
      function normho(nu,n,l)
        use factorials, only: dlfac
        implicit none
        real*8:: nu,pi,normho,fact,lnorm,logfac
        integer:: n,l,m,p
        pi=acos(-1d0)

        lnorm=(n+l+1)*log(2d0)+
     &        0.5d0*log(nu**(3d0/2d0)/sqrt(pi)) +
     &        0.5d0*(dlfac(n)+dlfac(n+l)-dlfac(2*n+2*l+1))
!     &        0.5d0*(logfac(n)+logfac(n+l)-logfac(2*n+2*l+1))

!        write(*,*) 'normho: n,l ->',n,l,logfac(n+l),dlfac(n+l)
        normho=dexp(lnorm)

!        if(normho<1e-6) then00
!           write(*,'(a,1f8.4,2i4,a,1g12.6)')
!     &    'Internal ERROR: nu,n,l=',nu,n,l,
!     &    ' Norm=',normho
!           stop
!        endif

      end function normho


c *** --------------------------------------------------------
c     Calculates normalization coefficient of HO wavefunctions: N(n,l) 
c *** --------------------------------------------------------      
      function normho2(nu,n,l)
        implicit none
        real*8:: nu,pi,normho2,fact
        integer:: n,l,m,p
        pi=acos(-1d0)

!!$c recursive formula
!!$        normho=2*(nu**3/pi)**0.25d0
!!$        do p=1,l ! recurrence: N(0,0) -> N(0,l)
!!$           write(*,*)'l=0 -> l'
!!$           normho=(2/(2*p+1))**0.5d0*normho
!!$        enddo
!!$        do m=1,n ! recurrence: N(0,l) -> N(n,l)
!!$           normho=(m/(m+l+0.5d0))**0.5d0*normho
!!$        enddo

        normho2=2**(n+l+1)*sqrt(nu**(3d0/2d0)/sqrt(pi))*
     c       sqrt(fact(n)*fact(n+l)/fact(2*n+2*l+1))
        if(normho2<1e-6) then
           write(*,*)'nu,n,l,normho2',nu,n,l,normho2
           write(*,*)'fact2(n)',fact(n)
           write(*,*)'fact2(n+l)=',fact(n+l)
           write(*,*)'fact2(2n+2l+1)=',fact(2*n+2*l+1)
           stop
        endif

      end function normho2


c ----------------------------------------------
c  factorial
c -----------------------------------------------
      function fact2(n)
        implicit none
        integer p,n
        real*8 fact2
        fact2=1
        do p=1,n
           fact2=fact2*p
        enddo
        if (fact2<1e-5) write(*,*)'n=',n,'fact2=0' 
      end function fact2





c -----------------------------------------------
c  factorial calcualted with Gamma function
c -----------------------------------------------
      function fact(n)
        implicit none
        integer n
        real*8 fact, dgamma,x
        x=dfloat(n+1)
        fact=dgamma(x)
      end function fact

c *** ---------------------------------------------
c Factorial in recursive form 
c *** ---------------------------------------------
      RECURSIVE FUNCTION factorialr(n) RESULT(nfact)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER :: nfact
        IF(n > 0) THEN
           nfact = n * factorialr(n-1)
        ELSE
           nfact = 1
        END IF
      END FUNCTION factorialr




