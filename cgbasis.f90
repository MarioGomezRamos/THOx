! *** ----------------------------------------------
!     Build CG basis
! *** -----------------------------------------------
      	
	!SUBRUTINA PRINCIPAL -  Obtendra los valores de las distintas funciones de onda en funcion de varios parametros
	
      subroutine cgbasis(l,r1,rnmax,nho)
	  
		!DECLARACION DE VARIABLES:
		
        use wfs, only: wftho,dr,rvec,nr !Revisar modulo wfs de Fortran
        implicit none
        integer:: nmin,ibas,jbas,nho,l,n,ir,kbout
        real*8:: r,wf,acg,r1,rnmax,wf1,wf2,nhomitad,norm1,norm2
!		real*8,allocatable:: wfaux(:,:),faux(:),gaux(:)
!		real*8,allocatable:: norm(:,:)
!		real*8,allocatable:: eigv(:,:)
	real*8:: wfaux(nho,nr),faux(nr),gaux(nr)
	real*8:: norm(nho,nho)
	real*8:: eigv(nho,nho)
	real*8 :: kk1
	real*8 :: kk2
	integer:: inousable
	real*8:: aux
	logical:: debug=.false.
	kk1 = 1
	kk2 = 1
	nhomitad = nho/2
	kbout=10
	wf=1 
        nmin=1
        ibas=0
        if (r1.lt.1e-6) then
          write(*,*) 'rmin too small for Gaussian basis: rmin=',r1
        endif

        if (rnmax.lt.1e-6) then
          write(*,*) 'rnmax too small for Gaussian basis: rnmax=',
     &               rnmax
        endif


        acg=(rnmax/r1)**(1./nhomitad)
        write(0,*)'cgbasis: acg=',acg
        do n=nmin,nhomitad 
           norm1=0
	   norm2=0
		    		
    	ibas=ibas+1 
            do ir=1,nr
	        r=rvec(ir)        
			  
		call cg3d(n,l,r,acg,r1,wf1,wf2)
		wfaux(ibas,ir)=wf1
		norm1=norm1+dr*(wf1*r)**2 
	
	    enddo
	    write(99,*)'&'		
	    ibas=ibas+1
	    do ir=1,nr
	        r=rvec(ir)        
		call cg3d(n,l,r,acg,r1,wf1,wf2)
		wfaux(ibas,ir)=wf2
		norm2=norm2+dr*(wf2*r)**2
				
           enddo     
	   write(99,*)'&'
           write(kbout,*)'#Norm=',norm1,norm2 
           write(kbout,*)'&'
        enddo


	
	do ibas=1,nho
	do jbas=ibas,nho
		faux(:) = wfaux(ibas,:)
		gaux(:) = wfaux(jbas,:)
		call normfun(faux,gaux,dr,rvec(1),nr,l,aux,kk1,kk2)
		norm(ibas,jbas) = aux
		if (ibas .ne. jbas) norm(jbas,ibas)=norm(ibas,jbas)
	enddo
	enddo
	

	
!SUBROUTINE HDIAG(H,N,NSMAX,IEGEN,U,NR). Me devuelve eigv (U) que sera una matriz nho,nho que contiene los coeficientes de los autovectores
	call  HDIAG(norm,nho,nho,0,eigv,inousable) 
	

        if (debug) write(0,*)'Norm eigenvalues:'
	do ibas=1,nho
	if (debug) write(0,*)'ibas,eigenvalue=',ibas,norm(ibas,ibas)

	enddo
	      
	do ibas=1,nho
	faux=0
        do jbas=1,nho
	    do ir=1,nr
	    aux=1./sqrt(norm(ibas,ibas))
            faux(ir)=faux(ir)+aux*eigv(ibas,jbas)*wfaux(jbas,ir)
	    enddo
        enddo
        wftho(ibas,l,:)=faux(:)
	enddo	

	do ibas=1,nho
		do jbas=ibas,nho
		faux(:) = wftho(ibas,l,:)
		gaux(:) = wftho(jbas,l,:)
		call normfun(faux,gaux,dr,rvec(1),nr,l,aux,kk1,kk2)
		norm(ibas,jbas) = aux
!		write(0,*)'ibas,jbas,sol=',ibas,jbas,aux
		if (ibas .ne. jbas) norm(jbas,ibas)=norm(ibas,jbas)
		enddo
	enddo	
	
	do ibas=1,nho
	write(30,*)(norm(ibas,jbas),jbas=1,nho)
	enddo	
	
	
		
      end subroutine cgbasis



! *** -------------------------------------------------
!     Calculates CG wavefunctions (radial part) in 3D
! *** ------------------------------------------------
      subroutine cg3d(n,l,r,acg,r1,wf1,wf2)
        implicit none
        integer l,n
        real*8:: pi,r,norma,normCG,wf,acg,alfa,eta1,eta2,rn,v,wf1,wf2
		real*8:: r1,rmax
		!complex*16 :: z1
		
!		r1=1
		pi=acos(-1d0)
		!z1 = (0,1)
		!z1 = dcmplx(0.d0, 1.d0)
		rn = r1*acg**(n-1) 
		if (rn.lt.1e-6) rn=1e-6
		v = 1./(rn**2)
		alfa = pi/2. 
		!alfa = 0
		!eta1 = (1+z1*alfa)*v
		!eta2 = (1-z1*alfa)*v
		
              norma=normcg(n,l,v)
		!write(*,*)n,r,acg,norma
		

		!Funciones Expresión Exponenciales
        !wf1=norma*r**(l)*((dexp(-eta1*r**(2)) + dexp(-eta2*r**(2)))/2)
		!wf2=norma*r**(l)*((dexp(-eta1*r**(2)) - dexp(-eta2*r**(2)))/z1*2)
		
		!Funciones Expresión Trigonometricas
		wf1=norma*r**(l)*(dexp(-v*r**(2)))*(cos(alfa*v*r**(2)))
		wf2=norma*r**(l)*(dexp(-v*r**(2)))*(sin(alfa*v*r**(2)))
		
		write(99,*)r,wf1,wf2
				
      end subroutine cg3d
 


! *** --------------------------------------------------------
!     Calculates normalization coefficient of CG wavefunctions: 
! *** --------------------------------------------------------      
      function normcg(n,l,v)
        use factorials, only: FACT,dlfac2
        implicit none
        real*8:: pi,acg,v,normcg, param_1, param_2 
        real*8:: param_numerador, param_denominador,param_fact 
        integer:: n,l
!	REAL :: param_numerador
!	REAL :: param_denominador
	real*8:: r1 !Luego lo meteremos como argumento
		
		!r1=1
        pi=acos(-1d0)
		!rn = r1*acg**(n-1)
		!if (rn.lt.1e-6) rn=1e-6
		!v = 1/(rn**2)
		param_1 = 2**(l+2)
		param_2 = 2*v
		param_fact = dexp(dlfac2(2*l+1))
		param_numerador = param_1 * param_2**(l+(3d0/2d0))
		param_denominador = (sqrt(pi) * param_fact)
		
			
		normcg = sqrt(param_numerador/param_denominador)
		
		write(105,*)n,l,v,normcg
		
	
		!lnorm = sqrt((param_1 * param_2**(l+(3d0/2d0))) / (sqrt(pi) * param_dlfac))
	 
	    !lnorm = sqrt(((2**(l+2))*((2*v)**(l+(3/2))))/(sqrt(pi)*(dlfac(2*l+1))))
		
		!lnorm = sqrt((2**(l+2)*(2*v)**(l+(3/2))/(sqrt(pi)*dlfac(2*l+1))
		


      end function normcg


! *** --------------------------------------------------------





