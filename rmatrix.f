      module rmat_mod
      implicit real*8(a,b,d-h,o-z)
      real*8 rmax0
      allocatable wle(:),xle(:),tc(:,:,:),blo0(:,:),blo1(:,:),blo2(:,:),
     1 q1(:,:),q2(:,:)
      end

      subroutine rmat_ini(nr,ns,rmax,zrma)
c This subroutine initializes the R-matrix calculation. 
c It computes:
c  the roots and weights of the Legendre quadrature
c  the matrix elements of the kinetic energy and of the Bloch operator
      use rmat_mod
      implicit real*8(a,b,d-h,o-z)
      dimension zrma(ns*nr)
      if(allocated(wle))deallocate(wle)
      if(allocated(xle))deallocate(xle)
      if(allocated(tc))deallocate(tc)
      if(allocated(blo0))deallocate(blo0)
      if(allocated(blo1))deallocate(blo1)
      if(allocated(blo2))deallocate(blo2)
      if(allocated(q1))deallocate(q1)
      if(allocated(q2))deallocate(q2)
      allocate (tc(nr,nr,ns),blo0(nr,nr),blo1(nr,nr),blo2(nr,nr),
     1 wle(nr),xle(nr),q1(nr,ns),q2(nr,ns))
      rmax0=rmax/ns
c  roots and weights of the Legendre quadrature
      call legzo(nr,xle,wle)
c  matrix elements of the kinetic energy and of the Bloch operator
      do 31 is=1,ns
      do 30 i1=1,nr
      xi=xle(i1)
      xi2=xi*(1-xi)
      if(is.eq.1)then
      xx=4*nr*(nr+1)+3+(1-6*xi)/xi2
      tc(i1,i1,1)=xx/(3*xi2)
      blo0(i1,i1)=1/xi2
      else
      xlb=xi/(1-xi)*(nr*(nr+1)-1/(1-xi))
      xla=(1-xi)/xi*(-nr*(nr+1)+1/xi)
      tc(i1,i1,is)=(nr*nr+nr+6-2/xi2)/(3*xi2)+xlb-xla
      blo1(i1,i1)=(1-xi)/xi
      blo2(i1,i1)=xi/(1-xi)
      end if
      do 30 i2=1,i1-1
      xj=xle(i2)
      xj2=xj*(1-xj)
      if(is.eq.1)then
      xx=nr*(nr+1)+1+(xi+xj-2*xi*xj)/(xi-xj)**2-1/(1-xi)-1/(1-xj)
      tc(i1,i2,1)=xx/sqrt(xi2*xj2)
      blo0(i1,i2)=1/sqrt(xi2*xj2)
      else
      yy=sqrt(xj2/xi2**3)*(2*xi*xj+3*xi-xj-4*xi**2)/(xj-xi)**2
      xlb=sqrt(xi*xj/(1-xi)/(1-xj))*(nr*(nr+1)-1/(1-xj))
      xla=sqrt((1-xi)*(1-xj)/xi/xj)*(-nr*(nr+1)+1/xj)
      tc(i1,i2,is)=yy+xlb-xla
      blo1(i1,i2)=sqrt((1-xi)*(1-xj)/xi/xj)
      blo2(i1,i2)=sqrt(xi*xj/(1-xi)/(1-xj))
      end if
      if(mod(i1+i2,2).eq.1)then
      tc(i1,i2,is)=-tc(i1,i2,is)
      blo0(i1,i2)=-blo0(i1,i2)
      blo1(i1,i2)=-blo1(i1,i2)
      blo2(i1,i2)=-blo2(i1,i2)
      end if
      tc(i2,i1,is)=tc(i1,i2,is)
      blo0(i2,i1)=blo0(i1,i2)
      blo1(i2,i1)=blo1(i1,i2)
      blo2(i2,i1)=blo2(i1,i2)
   30 continue
      if(is.eq.1)then
      q2(:,1)=1/sqrt(xle(:)*(1-xle(:)))
      q1(:,1)=0
      else
      q2(:,is)=sqrt(xle(:)/(1-xle(:)))
      q1(:,is)=-1/q2(:,is)
      end if
      if(mod(nr,2).eq.1)q2(1:nr,is)=-q2(1:nr,is)
      q1(1:nr:2,is)=-q1(1:nr:2,is)
      q2(1:nr:2,is)=-q2(1:nr:2,is)
   31 continue
      tc=tc/rmax0**2
      blo0=blo0/rmax0
      blo1=blo1/rmax0
      blo2=blo2/rmax0
      q1=q1/sqrt(rmax0)
      q2=q2/sqrt(rmax0)
      nn=0
      do is=1,ns
      zrma(nn+1:nn+nr)=((is-1)+xle(1:nr))*rmax0
      nn=nn+nr
      end do
      return
      end


      subroutine rmatrix(nch,lval,qk,eta,rmax,nr,ns,cpot,cu,
     1 ncp1,ndim,nopen,twf,cf,nwf1,nwf2,nc,nvc,ncp2,cpnl)
c Computes the R-matrix and the associated collision matrix
c Inputs:
c  nch=number of channels
c  lval=array containing the ell values of all channels
c  qk=array containing the wave numbers (<0 if the channel is closed)
c  eta=array containing the Sommerfeld parameters
c  rmax=channel radius a
c  nr=number of Lagrange function per interval
c  ns=number of intervals
c  cpot(ncp1,ndim,ndim)= complex array with the potentials
c  ncp1=first dimension of cpot
c  ndim=2nd and 3rd dimensions of cpot
c  twf=logical variable (.true. to determine the wave function)
c  nwf1=first dimension of cf
c  nwf2=second dimension of cf
c  nc=number of entrance channels
c  nvc=array containing the references to the entrance channels
c  ncp2=first dimension of cpnl (must be ncp2=0 if no n-l potential is present)
c  cpnl(ncp2,ndim,ndim)=complex array with the non-local potentials
c
c Output:
c  cu(ndim,ndim)=complex array with the collision matrix
c  npoen=number of open channels
c  cf(nwf1,nwf2,nc)=complex array with the wave function
      use rmat_mod
      implicit real*8(a,b,d-h,o-z)
      implicit complex*16(c)
! AMoro: lval changed from default to integer
      integer lval(nch)
      integer nr,ns,ncp1,ndim,nopen,nwf1,nwf2,nc,nvc,ncp2
!      dimension lval(nch),qk(nch),eta(nch),cpot(ncp1,ndim,ndim),
      dimension qk(nch),eta(nch),cpot(ncp1,ndim,ndim),

     1 cu(ndim,ndim),cf(nwf1,nwf2,nc),nvc(nc),cx(3)
      dimension cpnl(ncp2,ndim,ndim)
      logical twf,tnl
      allocatable ch(:,:,:),crma(:,:,:,:),co(:),cop(:),
     1 cz(:,:),crma0(:,:,:),crma2(:,:,:),fc(:),dfc(:),gc(:),dgc(:),
     2 cfp(:,:),npo(:,:),xc(:,:)
      tnl=ncp2.ne.0
      if(tnl.and.ns.ne.1)then
      print*,'ns must be ns=1 for non-local potentials'
      stop
      end if
      lmax=maxval(lval)+1
      ntot=nch
      ntot2=ntot*nr
      allocate (fc(lmax),dfc(lmax),gc(lmax),dgc(lmax))
      allocate (ch(ntot2,ntot2,ns),crma(ntot,ntot,3,ns),co(ntot),
     1 cop(ntot),cz(ntot,ntot),crma2(ntot,ntot,ns),
     2 crma0(ntot,ntot,0:ns),cfp(ntot,0:ns),xc(ntot,2),npo(ntot,2))
      nopen=0
      nclo=0
      cu=0
      cfp=0

c Stores Coulomb and Whittaker functions
      do 2 i1=1,nch
      l=lval(i1)
      ll=l+1
      if(qk(i1).gt.0)then
      nopen=nopen+1
      npo(nopen,1)=i1
      xl=l
      call coulfg(qk(i1)*rmax,eta(i1),xl,xl,fc,gc,dfc,dgc,1,0,ifail)
      co(nopen)=dcmplx(gc(ll),fc(ll))
      cop(nopen)=dcmplx(dgc(ll),dfc(ll))*qk(i1)*rmax
      else
      nclo=nclo+1
      npo(nclo,2)=i1
      qk2=-qk(i1)
      ie=0
      call whit(eta(i1),rmax,qk2,-qk2**2,l,fc,dfc,ie)
      xc(i1,1)=rmax*dfc(l+1)/fc(l+1)
      xc(i1,2)=fc(l+1)
      end if
    2 continue

c Calculation of matrix C
      rma0=0
      ms=-nr
      ch=0
      do 21 is=1,ns
      rmax1=(is-1)*rmax0
      rmax2=is*rmax0
      ms=ms+nr
      m1=0
      do 5 i1=1,ntot
      ch(m1+1:m1+nr,m1+1:m1+nr,is)=tc(:,:,is)
      fac=lval(i1)*(lval(i1)+1)
      do 8 ir=1,nr
      xx=rmax1+xle(ir)*rmax0
      ch(m1+ir,m1+ir,is)=ch(m1+ir,m1+ir,is)+fac/xx**2
      if(qk(i1).gt.0)ch(m1+ir,m1+ir,is)=ch(m1+ir,m1+ir,is)-qk(i1)**2
      if(qk(i1).lt.0)ch(m1+ir,m1+ir,is)=ch(m1+ir,m1+ir,is)+qk(i1)**2
    8 continue
      if(qk(i1).lt.0)then
      if(is.eq.1)ch(m1+1:m1+nr,m1+1:m1+nr,is)=
     1 ch(m1+1:m1+nr,m1+1:m1+nr,is)-xc(i1,1)*blo0/rmax2
      if(is.ne.1)ch(m1+1:m1+nr,m1+1:m1+nr,is)=
     1 ch(m1+1:m1+nr,m1+1:m1+nr,is)-xc(i1,1)*(blo2/rmax2-blo1/rmax1)
      end if
      m2=0
      do 6 i2=1,ntot
      ii=(is-1)*nr*nr
      do 7 ir=1,nr
      ch(m1+ir,m2+ir,is)=ch(m1+ir,m2+ir,is)+cpot(ms+ir,i1,i2)
      if(.not.tnl)go to 7
      do irp=1,nr
      ii=ii+1
      ch(m1+ir,m2+irp,is)=ch(m1+ir,m2+irp,is)+sqrt(wle(ir)*wle(irp))*
     1 rmax0*cpnl(ii,i1,i2)
      end do
    7 continue
      m2=m2+nr
    6 continue
      m1=m1+nr
    5 continue

c Inverse of matrix C
      call cminv_sym(ch(1,1,is),ntot2,ntot2)
c Calculation of the R matrix
      m1=0
      do 15 i1=1,ntot
      m2=0
      do 16 i2=1,ntot
      cx=0
      do 17 ir1=1,nr
      do 17 ir2=1,nr
      cx(1)=cx(1)+ch(m1+ir1,m2+ir2,is)*q2(ir1,is)*q2(ir2,is)
      cx(2)=cx(2)+ch(m1+ir1,m2+ir2,is)*q1(ir1,is)*q2(ir2,is)
      cx(3)=cx(3)+ch(m1+ir1,m2+ir2,is)*q1(ir1,is)*q1(ir2,is)
   17 continue
      crma(i1,i2,:,is)=cx
      m2=m2+nr
   16 continue
      m1=m1+nr
   15 continue
      cz=crma(:,:,3,is)+crma0(:,:,is-1)*rmax1
      call cminv_sym(cz,ntot,ntot)
      crma0(:,:,is)=matmul(transpose(crma(:,:,2,is)),cz)
      crma0(:,:,is)=matmul(crma0(:,:,is),crma(:,:,2,is))
      crma0(:,:,is)=crma(:,:,1,is)-crma0(:,:,is)
      crma0(:,:,is)=crma0(:,:,is)/rmax2
   21 continue

c Calculation of matrices Z_I and Z_O
      do 20 m1=1,nopen
      i1=npo(m1,1)
      do 20 m2=1,nopen
      i2=npo(m2,1)
      crma2(m1,m2,1:ns)=crma0(i1,i2,1:ns)*sqrt(qk(i1)/qk(i2))
   20 continue
      do 22 m1=1,nopen
      cz(1:nopen,m1)=-crma2(1:nopen,m1,ns)*cop(m1)
      cu(1:nopen,m1)=-crma2(1:nopen,m1,ns)*conjg(cop(m1))
      cz(m1,m1)=cz(m1,m1)+co(m1)
      cu(m1,m1)=cu(m1,m1)+conjg(co(m1))
   22 continue
      call cminv_nsym(cz,nopen,ntot)
c Calculation of the collision matrix
      cu(1:nopen,1:nopen)=matmul(cz(1:nopen,1:nopen),
     1 cu(1:nopen,1:nopen))

c Calculation of the wave function (if twf=.true.)
      if(.not.twf)return
c Loop over the entrance channels
      do 100 ic=1,nc
      nvex=nvc(ic)
      i2=npo(nvex,1)
      do 101 m1=1,nopen
      i1=npo(m1,1)
      cfp(i1,ns)=-cu(m1,nvex)*cop(i1)*sqrt(qk(i2)/qk(i1))
      if(i1.eq.i2)cfp(i1,ns)=cfp(i1,ns)+conjg(cop(i1))
  101 continue
      cfp(:,ns)=cfp(:,ns)/rmax
      do 102 is=ns-1,1,-1
      rin=is*rmax0
      cz=crma(:,:,3,is+1)+rin*crma0(:,:,is)
      call cminv_nsym(cz,nch,nch)
      cz=matmul(cz,crma(:,:,2,is+1))
      cfp(:,is)=matmul(cz,cfp(:,is+1))
  102 continue
      do 110 is=1,ns
      m1=0
      do 111 i1=1,nch
      do 111 ir1=1,nr
      m1=m1+1
      cy=0
      m2=0
      do 112 i2=1,nch
      do 112 ir2=1,nr
      m2=m2+1
      cy=cy+ch(m1,m2,is)*(q2(ir2,is)*cfp(i2,is)-
     1 q1(ir2,is)*cfp(i2,is-1)) 
  112 continue
      mm1=(is-1)*nr+ir1
      cf(mm1,i1,ic)=cy
  111 continue
  110 continue
c Calculation of bound-state amplitudes
      do 120 m1=1,nclo
      i1=npo(m1,2)
      cy=sum(cf((ns-1)*nr+1:ns*nr,i1,ic)*q2(:,ns))
      cu(m1+nopen,ic)=cy/xc(i1,2)
  120 continue
      do 122 i1=1,nch
      mm=0
      do 123 is=1,ns
      cf(mm+1:mm+nr,i1,ic)=cf(mm+1:mm+nr,i1,ic)/sqrt(rmax0*wle(1:nr))
      mm=mm+nr
  123 continue
  122 continue
  100 continue
      return
      end

	      subroutine cminv_sym(c,n,mdim)
c This subroutine computes the inverse of a symmetric complex matrix
c  zsytrf and  zsytri are included in the LAPACK library
c C=matrix to be inverted 
c mdim: dimension of matrix C as declared in the program from which
c       cminv_nsym is called
c n: order of matrix C
c
c 2 options are available:
c   1) use of the subroutine CMATINV included in the package
c   2) use of LAPACK subroutines
c LAPACK subroutines are faster and should be used for optimal performances
      implicit complex*16(c)
      allocatable cwork(:),inde(:)
      dimension c(mdim,mdim)
      allocate(inde(mdim))
c Option 1: use of subroutine CMATINV (included in the package)
      call cmatinv(c,mdim,n,inde,nerror,cdet)
c Option 2: use of the LAPACK subroutines
c     lwork=64*mdim
c     allocate(cwork(lwork))
c     call zsytrf('U',n,c,mdim,inde,cwork,lwork,info)
c     call zsytri('U',n,c,mdim,inde,cwork,info2)
c     do i=1,n
c     c(i,1:i)=c(1:i,i)
c     enddo
      return
      end

      subroutine cminv_nsym(c,n,mdim)
c This subroutine computes the inverse of a non-symmetric complex matrix
c  zgetrf and  zgetri are included in the LAPACK library
c C=matrix to be inverted 
c mdim: dimension of matrix C as declared in the program from which
c       cminv_nsym is called
c n: order of matrix C
c
c 2 options are available:
c   1) use of the subroutine CMATINV included in the package
c   2) use of LAPACK subroutines
c LAPACK subroutines are faster and should be used for optimal performances
      implicit complex*16(c)
      allocatable cwork(:),inde(:)
      dimension c(mdim,mdim)
      allocate(inde(mdim))
c Option 1: use of subroutine CMATINV (included in the package)
      call cmatinv(c,mdim,n,inde,nerror,cdet)
c Option 2: use of the LAPACK subroutines
c     lwork=64*mdim
c     allocate(cwork(lwork))
c     call zgetrf(n,n,c,mdim,inde,info)
c     call zgetri(n,c,mdim,inde,cwork,lwork,info2)
      return
      end

        SUBROUTINE LEGZO(N,X,W)
C
C       =========================================================
C       Purpose : Compute the zeros of Legendre polynomial Pn(x)
C                 in the interval [0,1], and the corresponding
C                 weighting coefficients for Gauss-Legendre
C                 integration
C       Input :   n    --- Order of the Legendre polynomial
C       Output:   X(n) --- Zeros of the Legendre polynomial
C                 W(n) --- Corresponding weighting coefficients
C       =========================================================
C
C Author: J. M. Jin
C Downloaded from http://jin.ece.illinois.edu/routines/routines.html
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION X(N),W(N)
        data pi/3.1415926535898d0/,one/1/
        N0=(N+1)/2
        DO 45 NR=1,N0
c          Z=COS(pi*(NR-0.25D0)/N)
           Z=COS(pi*(NR-0.25D0)/(N+0.5d0))
10         Z0=Z
           P=1
           DO 15 I=1,NR-1
15            P=P*(Z-X(I))
           F0=1
           IF (NR.EQ.N0.AND.N.NE.2*INT(N/2)) Z=0
           F1=Z
           DO 20 K=2,N
              PF=(2-one/K)*Z*F1-(1-one/K)*F0
              PD=K*(F1-Z*PF)/(1-Z*Z)
              F0=F1
20            F1=PF
           IF (Z.EQ.0) GO TO 40
           FD=PF/P
           Q=0
           DO 35 I=1,NR-1
              WP=1
              DO 30 J=1,NR-1
                 IF (J.NE.I) WP=WP*(Z-X(J))
30            CONTINUE
35            Q=Q+WP
           GD=(PD-Q*FD)/P
           Z=Z-FD/GD
           IF (ABS(Z-Z0).GT.ABS(Z)*1.0D-15) GO TO 10
40         X(NR)=Z
           X(N+1-NR)=-Z
           W(NR)=2/((1-Z*Z)*PD*PD)
45         W(N+1-NR)=W(NR)
        x(1:n)=(1+x(n:1:-1))/2
        w(1:n)=w(n:1:-1)/2
        RETURN
        END

C                                                                       ABNK0330
      SUBROUTINE COULFG(XX,ETA1,XLMIN,XLMAX, FC,GC,FCP,GCP,
     *                  MODE1,KFN,IFAIL)
C  REVISED #5 IJT WITH L-T ALGORITHMN FOR CONTINUED FRACTIONS,
C  AND IFAIL > 0 FOR AVOIDED EXPONENT CHECKS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
C                                                                      C
C  A. R. BARNETT           MANCHESTER  MARCH   1981                    C
C                                                                      C
C  ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395           C
C                 + 'RCWFF'      IN    CPC 11 (1976) 141-142           C
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314           C
C  THIS VERSION WRITTEN UP       IN    CPC 27 (1982) 147-166           C
C                                                                      C
C  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), C
C   AND REAL LAMBDA(XLMIN) .GT. -1 FOR INTEGER-SPACED LAMBDA VALUES    C
C   THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER   C
C   EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF     C
C   THE DIRAC EQUATION ,ALSO SPHERICAL & CYLINDRICAL BESSEL EQUATIONS  C
C                                                                      C
C  FOR A RANGE OF LAMBDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,    C
C  STARTING ARRAY ELEMENT IS M1 = MAX0(IDINT(XLMIN+ACCUR),0) + 1       C
C      SEE TEXT FOR MODIFICATIONS FOR INTEGER L-VALUES                 C
C                                                                      C
C  IF 'MODE' = 1  GET F,G,F',G'   FOR INTEGER-SPACED LAMBDA VALUES     C
C            = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN    C
C            = 3      F               CALL TO AT LEAST LENGTH (1)      C
C  IF 'KFN'  = 0 REAL        COULOMB FUNCTIONS ARE RETURNED            C
C            = 1 SPHERICAL   BESSEL      "      "     "                C
C            = 2 CYLINDRICAL BESSEL      "      "     "                C
C  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT                          C
C                                                                      C
C  PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'    C
C   IN OSCILLATING REGION X .GE. ETA1 + SQRT(ETA1**2 + XLM(XLM+1))     C
C   COULFG IS CODED FOR REAL*8 ON IBM OR EQUIVALENT  ACCUR = 10**-16   C
C   USE AUTODBL + EXTENDED PRECISION ON HX COMPILER  ACCUR = 10**-33   C
C   FOR MANTISSAS OF 56 & 112 BITS. FOR SINGLE PRECISION CDC (48 BITS) C
C   REASSIGN DSQRT=SQRT ETC.  SEE TEXT FOR COMPLEX ARITHMETIC VERSION  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION    FC(*),GC(*),FCP(*),GCP(*)
      LOGICAL      ETANE0,XLTURN
      COMMON       /STEE / PACCQ,NFP,NPQ,IEXP,M1
C***  COMMON BLOCK IS FOR INFORMATION ONLY.  NOT REQUIRED IN CODE
C***  COULFG HAS CALLS TO: DSQRT,DABS,DMOD,IDINT,DSIGN,DFLOAT,DMIN1
      DATA ZERO,ONE,TWO,TEN2,ABORT /0.0D0, 1.0D0, 2.0D0, 1.0D2, 2.0D4/
      DATA HALF,TM30,BIG / 0.5D0, 1.0D-30, 1.0D+100/
      DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 D0/
C *** THIS CONSTANT IS  DSQRT(TWO/PI):  USE Q0 FOR IBM REAL*16: D0 FOR
C ***  REAL*8 & CDC DOUBLE P:  E0 FOR CDC SINGLE P; AND TRUNCATE VALUE.
C
                        ACCUR = 1.0D-16
C ***            CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
      MODE  = 1
      IF(MODE1 .EQ. 2 .OR. MODE1 .EQ. 3 ) MODE = MODE1
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      ETA   = ETA1
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
                 ETANE0  = ETA .NE. ZERO
      ACC   = ACCUR * 10D0
      ACC4  = ACC*TEN2*TEN2
      ACCH  = DSQRT(ACC)
C ***    TEST RANGE OF XX, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE
C
      IF(XX .LE. ACCH)                          GO TO 100
      X     = XX
      XLM   = XLMIN
      IF(KFN .EQ. 2)  XLM = XLM - HALF
      IF(XLM .LE. -ONE .OR. XLMAX .LT. XLMIN)   GO TO 105
      E2MM1 = ETA*ETA + XLM*XLM + XLM
      XLTURN= X*(X - TWO*ETA) .LT. XLM*XLM + XLM
      DELL  = XLMAX - XLMIN + ACC
      IF(DABS(DMOD(DELL,ONE)) .GT. 2*ACC) WRITE(6,2040)XLMAX,XLMIN,DELL
      LXTRA = IDINT(DELL)
      XLL   = XLM + DFLOAT(LXTRA)
C ***       LXTRA IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
C ***       XLL  IS MAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS
C ***         DETERMINE STARTING ARRAY ELEMENT (M1) FROM XLMIN
      M1  = MAX0(IDINT(XLMIN + ACC),0) + 1
      L1  = M1 + LXTRA
C
C ***    EVALUATE CF1  =  F   =  FPRIME(XL,ETA,X)/F(XL,ETA,X)
C
      XI  = ONE/X
      FCL = ONE
      PK  = XLL + ONE
      PX  = PK  + ABORT
      F   =  ETA/PK + PK*XI
         IF(DABS(F).LT.TM30) F = TM30
         D = ZERO
         C = F
C
C ***   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1
C
    4 PK1   = PK + ONE
        EK  = ETA / PK
        RK2 = ONE + EK*EK
        TK  = (PK + PK1)*(XI + EK/PK1)
        D   =  TK - RK2 * D
        C   =  TK - RK2 / C
         IF(DABS(C).LT.TM30) C = TM30
         IF(DABS(D).LT.TM30) D = TM30
         D = ONE/D
         DF = D * C
         F  = F * DF
            IF(D .LT. ZERO) FCL = - FCL
         PK = PK1
                          IF( PK .GT. PX ) GO TO 110
      IF(DABS(DF-ONE) .GE. ACC)             GO TO 4
                  NFP = PK - XLL - 1
      IF(LXTRA .EQ. 0)                          GO TO 7
C
C *** DOWNWARD RECURRENCE TO LAMBDA = XLM. ARRAY GC,IF PRESENT,STORES RL
C
      FCL = FCL/BIG
      FPL = FCL*F
      IF(MODE .EQ. 1) FCP(L1) = FPL
                      FC (L1) = FCL
      XL  = XLL
      RL  = ONE
      EL  = ZERO
      DO 6  LP = 1,LXTRA
         IF(ETANE0) EL = ETA/XL
         IF(ETANE0) RL = DSQRT(ONE + EL*EL)
         SL    =  EL  + XL*XI
         L     =  L1  - LP
         FCL1  = (FCL *SL + FPL)/RL
         FPL   =  FCL1*SL - FCL *RL
         FCL   =  FCL1
         FC(L) =  FCL
         IF(MODE .EQ. 1) FCP(L)  = FPL
         IF(MODE .NE. 3 .AND. ETANE0) GC(L+1) = RL
	 if(abs(FCL).gt.BIG) then
		do 55 LP1=L,M1+LXTRA
      		IF(MODE .EQ. 1) FCP(LP1) = FCP(LP1)*1d-20
 55                    		FC (LP1) = FC(LP1)*1d-20
		FCL=FC(L)
		FPL=FPL*1d-20
		endif
    6 XL = XL - ONE
      IF(FCL .EQ. ZERO) FCL = ACC
      F  = FPL/FCL
C ***    NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
C ***    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
C ***    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM
C
    7 IF( XLTURN ) CALL JWKB(X,ETA,DMAX1(XLM,ZERO),FJWKB,GJWKB,IEXP)
      IF( IEXP .GT. 1 .OR. GJWKB .GT. ONE/(ACCH*TEN2))  GO TO 9
          XLTURN = .FALSE.
      TA =  TWO*ABORT
      PK =  ZERO
      WI =  ETA + ETA
      P  =  ZERO
      Q  =  ONE - ETA*XI
      AR = -E2MM1
      AI =  ETA
      BR =  TWO*(X - ETA)
      BI =  TWO
      DR =  BR/(BR*BR + BI*BI)
      DI = -BI/(BR*BR + BI*BI)
      DP = -XI*(AR*DI + AI*DR)
      DQ =  XI*(AR*DR - AI*DI)
    8 P     = P  + DP
         Q  = Q  + DQ
         PK = PK + TWO
         AR = AR + PK
         AI = AI + WI
         BI = BI + TWO
         D  = AR*DR - AI*DI + BR
         DI = AI*DR + AR*DI + BI
         C  = ONE/(D*D + DI*DI)
         DR =  C*D
         DI = -C*DI
         A  = BR*DR - BI*DI - ONE
         B  = BI*DR + BR*DI
         C  = DP*A  - DQ*B
         DQ = DP*B  + DQ*A
         DP = C
         IF(PK .GT. TA)                         GO TO 120
      IF(DABS(DP)+DABS(DQ).GE.(DABS(P)+DABS(Q))*ACC)   GO TO 8
                      NPQ   = PK/TWO
                      PACCQ = HALF*ACC/DMIN1(DABS(Q),ONE)
                      IF(DABS(P) .GT. DABS(Q)) PACCQ = PACCQ*DABS(P)
C
C *** SOLVE FOR FCM = F AT LAMBDA = XLM,THEN FIND NORM FACTOR W=W/FCM
C
      GAM = (F - P)/Q
            IF(Q .LE. ACC4*DABS(P))             GO TO 130
      W   = ONE/DSQRT((F - P)*GAM + Q)
            GO TO 10
C *** ARRIVE HERE IF G(XLM) .GT. 10**6 OR IEXP .GT. 70 & XLTURN = .TRUE.
    9 W   = FJWKB
      GAM = GJWKB*W
      P   = F
      Q   = ONE
C
C *** NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
C
   10                     ALPHA = ZERO
          IF(KFN  .EQ. 1) ALPHA = XI
          IF(KFN  .EQ. 2) ALPHA = XI*HALF
                          BETA  = ONE
          IF(KFN  .EQ. 1) BETA  = XI
          IF(KFN  .EQ. 2) BETA  = DSQRT(XI)*RT2DPI
      FCM  = DSIGN(W,FCL)*BETA
           FC(M1)  = FCM
                      IF(MODE .EQ. 3)           GO TO 11
           IF(.NOT. XLTURN)   GCL =  FCM*GAM
           IF(      XLTURN)   GCL =  GJWKB*BETA
           IF( KFN .NE. 0 )   GCL = -GCL
           GC(M1)  = GCL
           GPL =  GCL*(P - Q/GAM) - ALPHA*GCL
                      IF(MODE .EQ. 2)           GO TO 11
           GCP(M1) = GPL
           FCP(M1) = FCM*(F - ALPHA)
   11 IF(LXTRA .EQ. 0 ) RETURN
C *** UPWARD RECURRENCE FROM GC(M1),GCP(M1)  STORED VALUE IS RL
C *** RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
C ***    XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
         W    = BETA*W/DABS(FCL)
         MAXL = L1 - 1
      DO 12 L = M1,MAXL
                      IF(MODE .EQ. 3)           GO TO 12
                      XL = XL + ONE
         IF(ETANE0)   EL = ETA/XL
         IF(ETANE0)   RL = GC(L+1)
                      SL = EL + XL*XI
         GCL1     = ((SL - ALPHA)*GCL - GPL)/RL
		IF(ABS(GCL1).gt.BIG) GO TO 140
         GPL      =   RL*GCL -  (SL + ALPHA)*GCL1
         GCL      = GCL1
         GC(L+1)  = GCL1
                      IF(MODE .EQ. 2)           GO TO 12
         GCP(L+1) = GPL
         FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
   12 FC(L+1)     = W* FC(L+1)
      RETURN
C
C ***    ERROR MESSAGES
C
  100 IFAIL = -1
      WRITE(6,2000) XX,ACCH
 2000 FORMAT(' FOR XX = ',1P,D12.3,' TRY SMALL-X  SOLUTIONS',
     *' OR X NEGATIVE'/ ,' SQUARE ROOT ACCURACY PARAMETER =  ',D12.3/)
      RETURN
  105 IFAIL = -2
      WRITE (6,2005) XLMAX,XLMIN,XLM
 2005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES:XLMAX,XLMIN,XLM = ',
     *1P,3D15.6/)
      RETURN
  110 IFAIL = -3
      WRITE (6,2010) ABORT,F ,DF,PK,PX,ACC
 2010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/
     *' F,DF,PK,PX,ACCUR =  ',1P,5D12.3//)
      RETURN
  120 IFAIL = -4
      WRITE (6,2020) ABORT,P,Q,DP,DQ,ACC
 2020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/
     *' P,Q,DP,DQ,ACCUR =  ',1P,4D17.7,D12.3//)
      RETURN
  130 IFAIL = -5
      WRITE (6,2030) P,Q,ACC,DELL,LXTRA,M1
 2030 FORMAT(' FINAL Q.LE.DABS(P)*ACC*10**4 , P,Q,ACC = ',1P,3D12.3,4X,
     *' DELL,LXTRA,M1 = ',D12.3,2I5 /)
      RETURN
 2040 FORMAT(' XLMAX - XLMIN = DELL NOT AN INTEGER ',1P,3D20.10/)
  140 IFAIL=L1-L
      RETURN
      END
C
      SUBROUTINE JWKB_bis(XX,ETA1,XL,FJWKB,GJWKB,IEXP)                      ABNK0331
      implicit real*8(a-h,o-z)
c     REAL*8          XX,ETA1,XL,FJWKB,GJWKB,DZERO                      ABNK0332
C *** COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS    FOR XL.GE. 0 ABNK0333
C *** AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554       ABNK0334
C *** CALLS DMAX1,SQRT,ALOG,EXP,ATAN2,FLOAT,INT        BARNETT FEB 1981 ABNK0335
      DATA   ZERO,HALF,ONE,SIX,TEN/ 0, 0.5d0, 1, 6, 10 /                ABNK0336
c     DATA  DZERO, RL35, ALOGE  /0,35, 0.43429 45 E0 /                  ABNK0337
      DATA  DZERO, RL35 /0,35 /                                         ABNK0337
      aloge=log(ten)
      X     = XX                                                        ABNK0339
      ETA   = ETA1                                                      ABNK0339
      GH2   = X*(ETA + ETA - X)                                         ABNK0340
      XLL1  = MAX(XL*XL + XL,DZERO)                                     ABNK0341
      IF(GH2 + XLL1 .LE. ZERO) RETURN                                   ABNK0342
       HLL  = XLL1 + SIX/RL35                                           ABNK0343
       HL   = SQRT(HLL)                                                 ABNK0344
       SL   = ETA/HL + HL/X                                             ABNK0345
       RL2  = ONE + ETA*ETA/HLL                                         ABNK0346
       GH   = SQRT(GH2 + HLL)/X                                         ABNK0347
       PHI  = X*GH - HALF*( HL*LOG((GH + SL)**2/RL2) - LOG(GH) )        ABNK0348
          IF(ETA .NE. ZERO) PHI = PHI - ETA*ATAN2(X*GH,X - ETA)         ABNK0349
      PHI10 = -PHI*ALOGE                                                ABNK0350
      IEXP  =  INT(PHI10)                                               ABNK0351
      IF(IEXP .GT. 70) GJWKB = TEN**(PHI10 - IEXP)                      ABNK0352
      IF(IEXP .LE. 70) GJWKB = EXP(-PHI)                                ABNK0353
      IF(IEXP .LE. 70) IEXP  = 0                                        ABNK0354
      FJWKB = HALF/(GH*GJWKB)                                           ABNK0355
      RETURN                                                            ABNK0356
      END                                                               ABNK0357

  
      SUBROUTINE WHIT_bis(HETA,R,XK,E,LL,F,FD,IE)
C
C     CALCULATES  WHITTAKER  FUNCTION  WL(K,R)  WITH
C     ASYMPTOTIC  FORM  EXP(-(KR + ETA(LOG(2KR)))
C     E  IS  NEGATIVE
C     If IE = 0, allowed to return result e**IE larger than Whittaker,
C                for the IE value returned.
C     If IE > 0, must scale results by that amount.
C
C     Author: I.J. Thompson, downloaded from http://www.fresco.org.uk
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(LL+1),FD(LL+1) ,T(12),S(7)
      data fpmax/1.0d290/
      L = LL+1
C              NOW L = NO. OF VALUES TO FIND
      EE=-1.0
      AK=XK
      ETA=HETA
      LP1=L+1
      RHO=AK*R
	S(:) = 0
      IF(L-50)1,1,2
    1 LM=60
      GO TO 3
    2 LM=L+10
    3 LMP1=LM+1
      IS=7
      PJE=30.0*RHO+1.0
      H=max(INT(PJE),4)
      H=RHO/H
!	write(147,111) R,RHO,H
!111	format(3f10.6)
      RHOA=10.0*(ETA+1.0)
      IF(RHOA-RHO)13,13,14
   13 IFEQL=1
      RHOA=RHO
      GO TO 15
   14 IFEQL=0
   15 PJE=RHOA/H+0.5
      RHOA=H*INT(PJE)
      IF(IFEQL)16,16,18
   16 IF(RHOA-RHO-1.5*H)17,18,18
   17 RHOA=RHO+2.0*H
   18 IF(EE)55,55,19
   19 STOP 'WHIT'
   27 A=2.0-10.0/12.0*H*H*EE
      B=1.0/6.0*H*ETA
      C=1.0+1.0/12.0*H*H*EE
      M1=INT(RHOA/H-0.5)
      M2=INT(RHO/H-1.5)
      T(2)=B/FLOAT(M1+1)
      T(3)=B/FLOAT(M1)
      JS=M1
      DO 29 IS=M2,M1
      DO 28 I=1,6
      S(I)=S(I+1)
   28 CONTINUE
      T(1)=T(2)
      T(2)=T(3)
      T(3)=B/FLOAT(JS-1)
      S(7)=((A+10.0*T(2))*S(6)-(C-T(1))*S(5))/(C-T(3))
      JS=JS-1
      IF(ABS(S(7)).LE.FPMAX) GO TO 29
       DO 285 I=2,7
  285   S(I) = S(I) / FPMAX
   29 CONTINUE
      T(1)=S(4)
      T(2)=(1.0/60.0*(S(1)-S(7))+0.15*(S(6)-S(2))+0.75*(S(3)-S(5)))/H
      GO TO 60
   55 C=1.0/RHOA
      A=1.0
      B=1.0-C*ETA
      F(1)=A
      FD(1)=B
      DO 56 M=1,26
      D=0.5*(ETA+FLOAT(M-1))*(ETA+FLOAT(M))*C/FLOAT(M)
      A=-A*D
      B=-B*D-A*C
      F(1)=F(1)+A
      FD(1)=FD(1)+B
   56 CONTINUE
      A=-ETA*LOG(2.0*RHOA)-RHOA
      FPMINL = -LOG(FPMAX)
      if(IE.eq.0.and.A.LT.FPMINL) IE = INT(FPMINL-A)
      A=EXP(A+IE)
      F(1)=A*F(1)
c      FD(1)=A*FD(1)
      FD(1)=A*FD(1) * (-1d0 - 2*ETA/(RHOA))
      IF(IFEQL)57,57,61
   57 S(IS)=F(1)
      IF(IS-7)27,58,27
   58 IS=6
      RHOA=RHOA+H
      GO TO 55
   60 F(1)=T(1)
      FD(1)=T(2)
   61 C=1.0/RHO
      DO 63 M=1,L-1
      A=ETA/FLOAT(M)
      B=A+C*FLOAT(M)
      F(M+1)=(B*F(M)-FD(M))/(A+1.0)
      FD(M+1)=(A-1.0)*F(M)-B*F(M+1)
   63 CONTINUE
      DO 65 M=1,L
      FD(M)=AK*FD(M)
   65 CONTINUE
      RETURN
      END

      SUBROUTINE CMATINV (A,IDIM1,N1,INDE,NERROR,DETERM)
      IMPLICIT complex*16(A-H,O-Z)
      real*8 pivot2
C
C         Inversion of matrix A (overwritten on exi)
C
      DIMENSION A(IDIM1*IDIM1),INDE(IDIM1)
      DETER=1
      N=N1
      IEMAT=N
      KDIM=IDIM1
      NMIN1=N-1
C        THE ROUTINE DOES ITS OWN EVALUATION FOR DOUBLE SUBSCRIPTING OF
C        ARRAY A.
      IPIVC=1-KDIM
C        MAIN LOOP TO INVERT THE MATRIX
      DO 11 MAIN=1,N
      PIVOT=0
      IPIVC=IPIVC+KDIM
C        SEARCH FOR NEXT PIVOT IN COLUMN MAIN.
      IPIVC1=IPIVC+MAIN-1
      IPIVC2=IPIVC +NMIN1
      DO 2 I1=IPIVC1,IPIVC2
      IF(ABS(A(I1))-ABS(PIVOT)) 2,2,1
    1 PIVOT=A(I1)
      LPIV=I1
    2 CONTINUE
C        IS PIVOT DIFFERENT FROM ZERO
      pivot2=pivot
      IF(PIVOT2) 3,15,3
C        GET THE PIVOT-LINE INDICATOR AND SWAP LINES IF NECESSARY
    3 ICOL=LPIV-IPIVC+1
      INDE(MAIN)=ICOL
      IF(ICOL-MAIN) 6,6,4
C        COMPLEMENT THE DETERMINANT
    4 DETER=-DETER
C        POINTER TO LINE PIVOT FOUND
      ICOL=ICOL-KDIM
C        POINTER TO EXACT PIVOT LINE
      I3=MAIN-KDIM
      DO 5 I=1,IEMAT
      ICOL=ICOL+KDIM
      I3=I3+KDIM
      SWAP=A(I3)
      A(I3)=A(ICOL)
    5 A(ICOL)=SWAP
C        COMPUTE DETERMINANT
    6 DETER=DETER*PIVOT
      PIVOT=1/PIVOT
C        TRANSFORM PIVOT COLUMN
      I3=IPIVC+NMIN1
      DO 7 I=IPIVC,I3
    7 A(I)=-A(I)*PIVOT
      A(IPIVC1)=PIVOT
C        PIVOT ELEMENT TRANSFORMED
C
C        NOW CONVERT REST OF THE MATRIX
      I1=MAIN-KDIM
C        POINTER TO PIVOT LINE ELEMENTS
      ICOL=1-KDIM
C        GENERAL COLUMN POINTER
      DO 10 I=1,IEMAT
      ICOL=ICOL+KDIM
      I1=I1+KDIM
C        POINTERS MOVED
      IF(I-MAIN) 8,10,8
C        PIVOT COLUMN EXCLUDED
    8 JCOL=ICOL+NMIN1
      SWAP=A(I1)
      I3=IPIVC-1
      DO 9 I2=ICOL,JCOL
      I3=I3+1
    9 A(I2)=A(I2)+SWAP*A(I3)
      A(I1)=SWAP*PIVOT
   10 CONTINUE
   11 CONTINUE
C        NOW REARRANGE THE MATRIX TO GET RIGHT INVERS
      DO 14 I1=1,N
      MAIN=N+1-I1
      LPIV=INDE(MAIN)
      IF(LPIV-MAIN) 12,14,12
   12 ICOL=(LPIV-1)*KDIM+1
      JCOL=ICOL+NMIN1
      IPIVC=(MAIN-1)*KDIM+1-ICOL
      DO 13 I2=ICOL,JCOL
      I3=I2+IPIVC
      SWAP=A(I2)
      A(I2)=A(I3)
   13 A(I3)=SWAP
   14 CONTINUE
      DETERM=DETER
      NERROR=0
      RETURN
   15 NERROR=MAIN
      DETERM=DETER
      RETURN
      END
