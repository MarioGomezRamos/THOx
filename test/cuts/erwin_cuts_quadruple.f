subroutine schcc_erwin_cuts_Q(nch,ecm,z12,incvec,ql,factor,dr,r0,
     & npt,wf,phase,smat,method,info,einc,icc)
      use xcdcc, only: smats
      use nmrv,only: nr,h,vcoup,conv,ech,debug,rmin,hort,rvec,cutr
      use globals, only: verb
      use factorials
      use constants , only: e2
      use memory
      use trace , only: cdccwf
      implicit none
      
      integer:: method
      logical:: copen(nch),orto,raynal,info
      integer:: nch,klog,npt,ql(nch)
      integer:: ir,irmin,is,isp,ich,inc,l,lmax
      integer:: i,j,k,im,ip,i0,norto,icc
c     .........................................................
      real*8     :: ymem,ti,tf
      real*8 , parameter:: pi=acos(-1d0)
      real*8 :: ecm,tkch,eta,z12,r12,ena2,ena3
      real*8 :: dr,r,r0,rm,h2,factor,einc,kinc
      real*8 :: rturn,rmorto,rmatch
      real*8 :: kch(nch),kch2(nch),rturnv(nch),etav(nch),barrier(nch)
      real*8 :: start,finish,aux,big,small,tmin,tmax
      real*8 :: tol_mag,tol_mag_inc
      real*8, dimension(nch):: wkb
      integer:: rt_points(nch)
      integer:: irminaux
c     .........................................................
      complex*16, intent(out):: wf(nch,nr)
      
! MODIFIED: y remains complex*16 as requested
      complex*16 :: y(nch,nch,npt)

! MODIFIED: y0, z0, zm promoted to complex*32
      complex*32, dimension(nch,nch):: y0
      complex*32, dimension(nch,nch):: z0,zm
      
! MODIFIED: v moved from complex*16 to complex*32
      complex*16, dimension(nch,nch):: coupl
      complex*32, dimension(nch,nch):: v
      
      complex*16 :: phase(nch),smat(nch)
      complex*16 :: a(nch,nch)

! MODIFIED: s, ziv, zpv promoted to complex*32
      complex*32 :: s(nch),ziv,zpv
! MODIFIED: c1, c2, c promoted to complex*32
      complex*32 :: c1,c2,c
      
! MODIFIED: onec and zero are now complex*32 parameters
      complex*32,parameter:: onec=(1.0Q0, 0.0Q0),zero=(0.0Q0, 0.0Q0)
c ... for QR factorization
! MODIFIED: w0, wm promoted to complex*32
      complex*16 :: HDIAG(nch)
      complex*32 :: w0(nch,nch),wm(nch,nch)
c ... TEST (delete after debugging)
! MODIFIED: yaux promoted to complex*32
      complex*32:: yaux(nch,nch)
      
      logical incvec(nch)
      real*16:: renorm
      integer jr
c     ---------------------------------------------------------
      R12 = 1D0/12D0
      ENA2 = 2D0/5D0 * R12**2
      ENA3 = - 4D0/35D0 * R12**3
!
! MODIFIED: This line is redundant with the parameter, but kept for logic
      ONEC = (1.0Q0, 0.0Q0) 
      TMAX = 20.
      TMIN = -125.
      big=huge(big)**(0.8d0)
      small=epsilon(small)**0.8d0
     
      h=dr
      h2=h*h
      conv=factor
      copen(:)=.true.
      nr=npt
      klog=99
      rmin=0
      lmax=maxval(ql)
      orto=.false.
c ... Initialize some variables 
      
! MODIFIED: y remains complex*16 and is initialized with 0d0 (double)
      y=0d0      
      coupl=0
c ... Set TRUE for debugging
      debug=.false.
c ..............................................................................
      if (nch.lt.1) then 
        write(*,*)'Scatcc_erwin: nch=',nch; stop
      endif

      if (ecm.lt.1e-6) then
        write(*,*)'Scatcc_erwin: Ecm too small! =',ecm
        stop
      endif 
 
      if (allocated(rvec)) deallocate(rvec)
      allocate(rvec(nr))

      do ir=1,nr
        rvec(ir)=r0+(ir-1)*dr
      enddo
      rm      =rvec(nr)    !Max radius for integration
      rmatch  =rvec(nr-2)  !Matching radius if derivative is used
     
! CHanged by AMoro to prevent undefined mod(ir-irmin,norto) below
!      norto=nint(hort/h)
      norto=max(nint(hort/h),1)   

      if (abs(hort).gt.0.) orto=.true.
      if (info.and.(verb.ge.3)) then
        ymem=nch*nch*(npt+2)*lc16/1e6
        write(*,190) ymem
190     format(5x,"[ WFS require", 1f8.2," Mbytes ]")
        ymem=nch*nch*lc16/1e6
!        write(*,*)'Auxiliary arrays require 8x',
!     &  ymem,'=',ymem*8, 'Mbyes'
!        write(*,*)' Biggest real=',big
      endif 

      if (ecm.lt.0) then 
      write(*,200) ecm
200   format(4x,'Skipping negative energy:',1f7.3, ' MeV')
      wf(:,:)=0d0
      return
      endif

      rmorto=0
      kinc=sqrt(conv*ecm)
      do ich=1,nch
      aux=ecm+einc-ech(ich)
      kch2(ich)=conv*aux

      if (aux.gt.0) then
   
           copen(ich)=.true.
        if (debug) write(99,300) ich,aux,"open"
        kch(ich)=sqrt(conv*aux)
        etav(ich)=conv*z12*e2/kch(ich)/2.
        rturnv(ich) =(etav(ich)+SQRT(etav(ich)**2 
     &              + ql(ich)*(ql(ich)+1d0)))/kch(ich) 
        if(rturnv(ich).gt.rmorto) rmorto=rturnv(ich)
      else
        copen(ich)=.false.
        if (debug) write(99,300) ich,aux,"closed"
        kch(ich)  =sqrt(-conv*aux)
        rturnv(ich)=rm
        rmorto     =rm
      endif
300   format(3x,"Channel",i2," Final Ecm=",1f8.3," MeV ->",a6)
      enddo


c classical turning point for incoming channel
      eta=conv*z12*e2/kinc/2.
      l=1000
      do ich=1,nch
      if (incvec(ich)) then
      l=min(ql(ich),l)
      endif
      enddo
      RTURN =(eta+SQRT(eta**2 + L*(L+1d0)))/kinc  
 
      call factorialgen(lmax+1)
      call cpu_time(start)
c 
c Starting values and initial integration point
c 
!      rmin=h ;       irmin=1
!      rmin=max(h,rturn-10.) ! CHECK h or r0 ?
      rmin  =max(h,minval(rturnv(:)-abs(cutr)))
      irmin=nint(rmin/h)+1

      if (debug) write(95,*)'Starting value ir=,',irmin,'r=',rvec(irmin)
      
! MODIFIED: Initializing complex*32 variables with quad 'zero'
      z0(:,:)=zero;
      zm(:,:)=zero; 
      w0(:,:)=zero; wm(:,:)=zero
      y0(:,:)=zero; v(:,:)=zero 
      do ich=1,nch
        l=ql(ich)
!!!!!
      CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ZI(L,L) = H**(QNF(9,L)+1) / EXP(0.5 * FACT(QNF(9,L)+1))
! MODIFIED: Use quad precision for constant
        z0(ich,ich)= 1.0Q-10  !TEEEEEEEEEEEEEEEEEEEEEESSSSST!
!        z0(ich,ich)= H**(l+1.) / EXP(0.5 * FACT(l+1.))  
!         write(*,*)'ich,z0=:',  ich,z0(ich,ich)   
!!!!!! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
      
      if (info.and.(verb.ge.2)) then
      write(*,'(5x,"o Classical turning point: inner=",1f6.2," fm",
     &  5x,"outer=",1f8.2," fm",
     &  5x,"R-min=",1f6.2," fm")') 
     &  minval(rturnv(:)), maxval(rturnv(:)), rmin      
      endif

!Displace rturnv by cutr MGR
      do ich=1,nch
            rturnv(ich)=rturnv(ich)-abs(cutr)
      enddo


      if (irmin.lt.2) write(*,*)'** ERWIN: wrong irmin=',irmin
 
! Use WKB approximation to find cut point.
! Set a tolerance, which is the number of orders of magnitude we allow our wf to
!
! decrease within the classically forbidden region and add all decreases until we reach that tolerance
      tol_mag=8d0
      tol_mag_inc=8d0
      rt_points(:)=0
      wkb(:)=0d0
      do ir=nr-1,2,-1
      r=rvec(ir) 
      do ich=1,nch
      if (kch2(ich).lt.0d0)cycle
        barrier(ich)= (-l*(l+1)/r**2         ! centrifugal
     &          + kch2(ich)             ! energy
     &        - conv*dble(vcoup(ich,ich,ir))) ! potential
        if (barrier(ich).lt.0d0.and.rt_points(ich).eq.0) then !Compute WKB within forbidden region
             wkb(ich)=wkb(ich)+
     &        sqrt(abs(barrier(ich)))*h*log10(exp(1d0))
            if(.not. incvec(ich).and.wkb(ich).ge.tol_mag .and.
     &        rt_points(ich).eq.0) then
               rt_points(ich)=ir
            endif
            if(incvec(ich).and.wkb(ich).ge.tol_mag_inc .and.
     &        rt_points(ich).eq.0) then
               rt_points(ich)=ir
            endif
        endif
      enddo
      enddo

      irminaux=1000000
      do ich=1,nch
      if (rt_points(ich).eq.0) rt_points(ich)=2
      rturnv(ich)=max(rvec(rt_points(ich)),rturnv(ich))!Set turning point to WKB cut point if larger
      if (incvec(ich)) irminaux=min(irminaux,rt_points(ich))! Fix irminaux to smallest turning point of incoming channels
      enddo

      irmin=max(irmin,irminaux)

      if (2>1) then
      write(444,*) 'Turning points Jtot'
      do ich=1,nch
            irminaux=rt_points(ich)
      write(444,*) ich,rt_points(ich),wkb(ich),rvec(irminaux),irmin
      enddo      
      endif

      write(*,*) 'BIG',big,' SMALL',small
c Start radial integration
      do ir=irmin,nr
      r=rvec(ir) ! (ir-1)*h    
      im=ir-1  ! r-h
      i0=ir    ! r
      ip=ir+1  ! r+h
      if (r.lt.1e-10) r=1e-10

!Check if renormalization is needed
      renorm=1Q0
! MODIFIED: abs(y0) etc. are now quad precision operations
      do ich=1,nch
        do is=1,nch
          renorm=max(renorm,abs(y0(ich,is)),abs(z0(ich,is)),
     &      abs(zm(ich,is)))
        enddo
      enddo
!      if (mod(ir,10) .eq. 0) write(*,*)'Maximum at',ir,':',renorm
      if (renorm.gt.big) then
      renorm=1d0/renorm
      write(*,*)'Renormalizing y0 at ir,r=',ir,r,' renorm=',renorm
! MODIFIED: renormalizing complex*32 arrays z0, zm, y0
      do is=1,nch
      do ich=1,nch
      z0(ich,is)=z0(ich,is)*renorm
      zm(ich,is)=zm(ich,is)*renorm
      y0(ich,is)=y0(ich,is)*renorm
! MODIFIED: y remains complex*16, renorm is real*8 (double)
        do jr=irmin,ir-1
            y(ich,is,jr)=y(ich,is,jr)*renorm
        enddo
      enddo
      enddo
      endif

c kinetic energy part
c Fresco: SMAT(ich) = -LL1(K)*RI2 + (-ECM(ich) + C) * H
! MODIFIED: s is now complex*32
      do ich=1,nch    
      l=ql(ich)
      s(ich)= h2*(-l*(l+1)/r**2         ! centrifugal
     &          + kch2(ich)             ! energy
     &        - conv*vcoup(ich,ich,ir)) ! potential
      enddo

!Correct for rturnv MGR
      do ich=1,nch
      if (r.lt.rturnv(ich)) s(ich)=0d0
      enddo

! MODIFIED: y0, z0, s are complex*32, onec is complex*32
      do is=1,nch
      do ich=1,nch      
      y0(ich,is)= z0(ich,is)*
     X      (ONEC - S(ich) * (R12 - S(ich)*(ENA2 + S(ich)*ENA3)))
      enddo
      enddo

! MODIFIED: w0 is complex*32
      w0(:,:)=zero
     
       do is=1,nch
      w0(is,is)=(ONEC - S(is) * (R12 - S(is)*(ENA2 + S(is)*ENA3)))
      enddo 

!Correct for rturnv MGR
      do ich=1,nch
      if (r.lt.rturnv(ich))w0(ich,ich)=ONEC
      enddo

!      y0=matmul(w0,z0)

c non-diagonal part of V matrix 
      coupl(:,:)=-h2*conv*vcoup(:,:,ir)  
      do ich=1,nch
         coupl(ich,ich)=0d0
      enddo 

!  Correct for rturnv MGR
      do ich=1,nch
      do is=1,nch
         if (r.lt.rturnv(ich).or. r.lt.rturnv(is)) then
         coupl(ich,is)=zero
         if (ich.eq.is) then
! MODIFIED: Use quad precision constant
            y0(ich,ich)=1.0Q-10
         else
            y0(ich,is)=zero
         endif
         
         endif
      enddo
      enddo

! MODIFIED: w0 is complex*32
      w0=w0-r12*coupl
      
      if (2>1) then
      DO 24 ICH=1,NCH
	 DO 24 IS=1,NCH
! MODIFIED: C is complex*32
        C = COUPL(is,ich) * R12
! MODIFIED: y0, z0 are complex*32
	 if(C/=ZERO) y0(ich,1:nch) = y0(ich,1:nch) - C * z0(is,1:nch)
24    CONTINUE
        DO 34 ich=1,nch ! k
        V(ich,:)  = ZERO 
        DO 34 is=1,nch   ! j
	 C = COUPL(is,ich)
! MODIFIED: V, y0 are complex*32
 	 if(C/=ZERO) V(ich,1:nch) = V(ich,1:nch) + C * y0(is,1:nch)
34    CONTINUE
      else  !..........Idem using builtin matrix multiplication 
! MODIFIED: V is complex*32
        V(1:nch,1:nch)  = ZERO 
! MODIFIED: y0, z0 are complex*32
        y0=y0  - r12*matmul(coupl,z0)
! MODIFIED: v, y0 are complex*32
        v =   matmul(coupl,y0)
      endif

!! TEST 
      if (debug) then
! MODIFIED: yaux, w0, z0 are complex*32
      yaux=matmul(w0,z0)
      write(*,*)'y0 at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (y0(ich,is), is=1,min(15,nch))
      enddo
      endif

!!!!!!!!!!!!!!!!!!!!!!1


      if (debug) then
      write(*,*)'z0 at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (z0(ich,is), is=1,min(15,nch))
      enddo

      endif ! debug


c store solution 
! MODIFIED: This performs the requested precision demotion
! y (complex*16) <-- y0 (complex*32)
      y(:,:,ir)=y0(:,:)       


c zm <- z0, wm <-w0
      do is=1,nch
      do ich=1,nch
! MODIFIED: ziv, z0 are complex*32
      ziv=z0(ich,is)
! MODIFIED: All variables here are complex*32. Use 2.0Q0 for quad constant.
      zpv=2.0Q0*ziv - zm(ich,is) - v(ich,is) - s(ich)*y0(ich,is)
! MODIFIED: zm, ziv, z0, zpv are complex*32
      zm(ich,is)=ziv
      z0(ich,is)=zpv
! MODIFIED: wm, w0 are complex*32
      wm(ich,is)=w0(ich,is)
      enddo  !ich
      enddo  !is


c Re-orthogonalize solution vectors ....................................
      if (orto.and.(r.lt.rmorto).and.(ir.gt.irmin).and.(hort.gt.0)
     & .and.(mod(ir-irmin,norto).eq.0).and.(ir.lt.nr-10)) then 
      if (verb.ge.4) 
     &   write(*,'(5x,"->orthogonalizing at ir,r=",i4,1f7.2)')ir,r
  
       call cpu_time(ti)
! MODIFIED: Call to qrerwinz now passes complex*32 arrays
! NOTE: The subroutine qrerwinz MUST be modified to accept these.
       call qrerwinz(y0,z0,zm,w0,nch,ir)      !
! QR factorization

! MODIFIED: Store the *new* orthogonalized y0
! This again performs the complex*32 -> complex*16 demotion
       y(:,:,ir)=y0(:,:) ! new (orthogonolized) solution at ir


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (debug) then
      write(*,*)'y0 (after) at ir,r=',ir,r
      do ich=1,min(15,nch)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (y0(ich,is), is=1,min(15,nch))
      enddo

        write(*,*)'|yp x yp| (after QR)'
        do is=1,nch
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(y0(:,is),y0(:,ich)))
     &      / abs(dot_product(y0(:,is),y0(:,is)))
     &      / abs(dot_product(y0(:,ich),y0(:,ich))), ich=1,nch) 
        enddo
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call cpu_time(tf)
       torto=torto+tf-ti

      endif !orto
!........................................................................


!!!!!!!!!!!!!!!!!!!!!!!!! TEST
      if (0>1) then
! MODIFIED: y remains complex*16, so this write statement is unchanged
       write(91,'(5x,1f7.3,3x,50g14.5)') r,(y(1,ich,ir),ich=1,nch)
      endif       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         DO 46 IT=1,NEQS
!         DO 43 K=1,NEQS
!            ZIV = ZI(IT,K)
!            ZPV = ZIV + ZIV - ZM(IT,K) - V(IT,K) - SMAT(K) * FI(IT,K)
!            ZM(IT,K) = ZIV
!            ZI(IT,K) = ZPV
!43            CONTINUE
!46    CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      enddo !ir
      call cpu_time(finish)


c 
c Match y with asymptotic solution -> S-matrix
c            
      call cpu_time(start)
!      call matching(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! gauss5 + derivative
!      call matching2(ecm,z12,nch,ql,lmax,inc,rmatch,y,wf,phase,smat,info) ! LAPACK 
      do ich=1,nch
      if (incvec(ich)) then
      inc=ich
      write(445,*) 'Incoming channel=',inc
      call matching3(ecm,z12,nch,ql,lmax,inc,nr,y,wf,phase,smat,info) ! gauss5 + 2 points
      endif
            smats(icc,inc,1:nch)=smat(1:nch) 
      enddo
       call cpu_time(finish)
       tmatch=tmatch+finish-start

      if (debug) write(klog,'(2f10.3,3x,10g14.6)') ecm,z12,ql(1),
     &  (phase(ich),ich=1,nch) 


      call flush(6)
      end subroutine schcc_erwin_cuts
