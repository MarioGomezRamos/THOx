****PRE_EIGCC************************************************************
      subroutine pre_eigcc(nset,nchan,nodes)
      use globals, only: mu12,egs,kin,written,verb
      use constants
      use sistema
      use wfs, only:nr,dr,rvec,wfc,energ,rint,rmin,rlast,rmax
      use channels, only:jtot,sn,partot,nchmax,ql,cindex,qj,qjc,jpiset,
     &  jpsets
      use parameters, only: maxchan,maxeset
      use potentials, only: ccmat
      use xcdcc, only: binset,realwf
      use nmrv, only : nopen
      implicit none
c     ----------------------------------------------
      logical :: energy,tres,debug,ifhat,writewf
      integer, parameter:: kfr=20
      integer :: nset,inc,nbins,nchan,ib,ibp,ich,ik,ir,nk,li,nodes
      integer :: bastype,wftype,inchanf,incn,m,np
c     .....................................................
      real*8  :: emin,emax,ei,ef,ebin,excore,ehat,khatsq,wbin
      real*8  :: r,bnorm,raux,r2aux,bcoef,chnorm(maxchan),rms
      real*8  :: ddk,ki,kf,kmin,kmax,kstep,kmid,k
      real*8  :: fconv,yint,rm,ener
      real*8 ,   allocatable  :: cph(:)

      complex*16, parameter::iu=(0.,1.)
      complex*16 :: phc
      
      CHARACTER(LEN=80) formato
c  bincc variables (Fresco)
      logical tdel,tkmat,tknrm,bincc
      integer:: nf   ! number of formfactors?
      integer:: isc,il,pcon,lmax
      integer:: maxn,nnk(maxchan)
      real*8:: anc,k2(maxchan),conv
      complex*16:: y(nr-1,maxchan)
      real*8 etap(maxchan)
      real*8:: fmscal
      complex*16:: ccoup(maxchan,maxchan,nr)
      
      integer maxiter,ifail,nods
      real*8 eps
c  end bincc 
c *** ----------------------------------------------
       fconv=hc**2/2/mu12 ! E=fconv*k^2
      
      pcon=2;maxiter=300;eps=1e-5
      y(:,:)=0.0
      if ((rint.gt.0).and.(rint.lt.rvec(nr))) then
          maxn=ceiling((rint-0.0)/dr)+1 ! bincc assumes rmin=0!
      else
          maxn=nr
      endif

      inc     =jpiset(nset)%inc
      nchan   =jpiset(nset)%nchan
      
       call coefmat(nset,nchan)
       ccmat(:,:,:)=-ccmat(:,:,:)
      
       do ich=1,nchan
       ql(ich)    =jpiset(nset)%lsp(ich)  
       k2(ich)=abs(egs)+jpiset(nset)%exc(ich) ! separation energy for each channel
       conv=1/fconv     
       k2(ich)=k2(ich)*conv     
       etap(ich)= 2*e2*zv*zc*mu12/hc**2  !=2 k eta 
      enddo      
      lmax=maxval(ql)
      
      if (.not.allocated(wfc)) then
        nchmax=0
        allocate(wfc(jpsets,maxeset,maxchan,nr))
        allocate(energ(jpsets,maxeset))
      endif
      
      do ich=1,nchan
      k2(ich)=k2(ich)-conv*egs
      enddo
      
      ener=egs
      
      call EIGCC(y,CCMAT,ETAP,1d0,nr-1,K2,conv,ener,nr,inc,NODES,nr,dr   
     &  ,nchan,PCON,2*nchan+1,maxiter,EPS,IFAIL,lmax+1)
      
      
      energ(nset,1)=-ener
      write(*,*) 'Energy',ener,egs
      
      do ich=1,nchan
      wfc(nset,1,ich,1:maxn-1)=y(1:maxn-1,ich)
      enddo
      deallocate(ccmat)
      
      write(101,246)nchan,1,-ener
!        write(200+i,246)nchan,i,hmatx(n1,n1)

      

246	 format('# ',i2,' Channels, Eigenvalue:',i2,' Energy: ',f8.4)
       write(101,247)jtot,partot
!          write(200+i,247)jtot,partot
247	   format('# J:',f4.2,' parity:',i2)
        write(101,*) '#  channel core Jc (l sn) jn'
!           write(200+i,*)'#  channel core Jc (l sn) jn'
        
        do inchanf=1,nchan
        incn=inchanf
         
        write(101,245)incn,cindex(incn),qjc(incn),ql(incn),sn,qj(incn)
245	 format('# ',i2,' :',i2,f4.1,i2,f4.1,f4.1)
 
        enddo ! inchanf

       if (nchan.gt.nchmax) nchmax=nchan
       
     
     
         if ((rlast.gt.0).and.(rlast.lt.rmax)) then 
            np=idint((rlast-rmin)/dr)
        else
            np=idint((rmax-rmin)/dr)
        endif


! Write eigenfunctions 
         write(101,248)np
248     format("# ",i5," points")


       
! Compute & write vertex functions (march 2019) 
! Modified by Pedro (march 2022)
       
       do ir=1,np
        r=rmin+dr*dble(ir-1)
       enddo !ir
             
        do ir=1,np
        r=rmin+dr*dble(ir-1)
        if (r> rlast) cycle
!        write(100+i,'(1f8.3,2x,10g14.6)') r,(r*wfeig(i,m,ir),m=1,nchan)   
        write(101,'(1f8.3,2x,10f14.6)') r,
     &   (dble(wfc(nset,1,m,ir)),m=1,nchan)   
        
!        write(200+i,'(1f8.3,2x,10g14.6)') r,(vertex(m,ir),m=1,nchan)
        enddo !ir
        write(101,*)'& '
!        write(200+i,*)'& '



!        if (verb.ge.4) then
!        do ir=1,np
!        r=rmin+dr*dble(ir-1)
!        if (r> rlast) cycle
!        write(100+i,'(1f8.3,2x,10g14.6)') r,(r*vertex(m,ir),m=1,nchan)
!        enddo !ir
!        write(100+i,*)'& '
!        endif


      
      
      
      end subroutine





*****EIGCC**************************************************************
      SUBROUTINE EIGCC(PSI,CCMAT,ETAP,AB,MR,KAP2,THETA,P,
     &  MAXN,MC,NODES,NP,H,M,PCON,MM2,MAXC,EPS,IFAIL,LMX1)
		use factorials
c! AMoro ---------------------
!     	use io
c!	use drier
       use globals
       use channels, only: ql
	! NP here is MINT
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 PSI(NP+1,M),F(NP+1,M,M),                                  &
     &       MAT(MM2,MM2)
      real*8 RMAT(MM2,MM2)
      real*8 KAP2(M),ETAP(M)
      COMPLEX*16 ccmat(m,m,maxn)
      INTEGER PCON,COUNT,BC,CCP,CC
      complex*16 CENT(M),ZI(M,M),ZM(M,M),ZP(M,M)
     &             ,COUT(M),COUTP(M),COUPL(M,M)
      REAL*8 ETA,K,WL(LMX1),WLD(LMX1)
      LOGICAL SING
C
      ko=6
      fpmax = huge(acc8)**0.8d0
      COUNT = 1
      CCP = 1
      BC = 0
      PP=P
      N = NP-1
      RN = (NP-1)*H
      HP = H*H
      R12= 1./12.
      HP12=HP/12.
      NR = 2*M+1
      IFAIL = 0
      SMALL = 1D0 / FPMAX
      SMALLR = SQRT(SMALL)
!      write(*,*) 'smallr',smallr,'fpmax',fpmax
      SMALLQ = SQRT(SMALLR)
C     IF(NP.GT.MAXN.OR.NR.GT.MM2) STOP 101
      IF(NP.GT.MAXN.OR.NR.GT.MM2) CALL ABEND(8)
      IF(PCON.GE.3) WRITE(KO,207) NODES,MC
  207 FORMAT(/'    Parameters      Mismatch   Nodes ->',I2,' in ch',i3/)
102   DO 21 J=1,M
      CENT(J) = Ql(J) * (Ql(J)+1)
      K = SQRT(ABS(KAP2(J) + THETA*P)) * AB
!      write(*,*) 'K',k
      ETA  = 0.5*ETAP(J)/K
      L = Ql(J)
      IE = 0
      CALL WHIT(ETA,RN+H,K,E,L,WL,WLD,IE)
      COUTP(J) = WL(L+1)
      CALL WHIT(ETA,RN,K,E,L,WL,WLD,IE)
21    COUT(J) = WL(L+1)
      I   =(NP + 1) *3/4
         IMAX = NP/2
         DEL = -SQRT(FPMAX)
22    I=I-1
!      write(778,*) I,dble(CCMAT(MC,MC,I))
      UDIAG = -KAP2(MC) - THETA*P
23    UDIAG = UDIAG + CCMAT(MC,MC,I)
      DEN = -CENT(MC)/((I-1)*H)**2 + UDIAG
         IF(DEL.LT.DEN) THEN
            DEL = DEN
            IMAX = I
            ENDIF
      IF(DEN.LT.0 .AND. I.GT.10) GO TO 22
!      write(778,*)
!      write(778,*)
      MAM = I
      IF(I.EQ.10) MAM = IMAX
      MAP = MAM+1
      IF(COUNT.EQ.1 .AND. PCON.GE.7)
     &WRITE(*,*)MC,Ql(MC),K,ETA,RN+H,COUTP(MC),RN,COUT(MC),MAM,MAP
     &             ,CENT(MC)
103   DO 30 J=1,M
      DO 25 IT=1,M
      ZM(IT,J) = 0.0
25    ZI(IT,J) = 0.0
      ZM(J,J) = COUTP(J)
30    ZI(J,J) = COUT(J)
C
C      outer integration,  zi from np to map
      NT = -1
      NF = NP
      NO = NT * (MAP-NP) + 1
40    DO 60 III=1,NO
         I = NF + (III-1)*NT
         II= (I + MIN0(NT,0)) + 1
         RRI= 1.0/(II-1.)**2
      DO 42 IT=1,M
      DO 42 J=1,M
42    F(I,J,IT) = ZI(IT,J) * (1. + CENT(J)*RRI*R12 )
       DO 45 J=1,M
       DO 45 L=1,M
         C = 0.0
            T = CCMAT(L,J,II)
         C = C + T 
425      CONTINUE
         IF(L.EQ.J) C = C - KAP2(J) - THETA * P
       COUPL(L,J) = C
       IF(C.EQ.0) GO TO 44
       C = C * HP12
          DO 43 IT=1,M
43         F(I,L,IT) = F(I,L,IT) - C * ZI(IT,J)
44    CONTINUE
45    CONTINUE
      DO 49 IT=1,M
         DO 49 L=1,M
49       MAT(IT,L) = 0.0
      DO 54 L=1,M
      DO 54 J=1,M
      C = COUPL(L,J)
      IF(C.EQ.0.0) GO TO 54
      C = C * HP
      DO 53 IT=1,M
53      MAT(IT,L) = MAT(IT,L) + C * F(I,J,IT)
!        write(*,*) 'IT',IT, 'count',count,'C',c,'F',F(I,:,:)
!53      write(*,*) 'IT',IT, 'count',count,'MAT',mat
54    CONTINUE
      DO 55 J=1,M
      DO 55 IT=1,M
      ZP(IT,J) = 2*ZI(IT,J) - ZM(IT,J) - MAT(IT,J)
     &                      + F(I,J,IT) * CENT(J) * RRI
      ZM(IT,J) = ZI(IT,J)
55    ZI(IT,J) = ZP(IT,J)
C                              now check for incipient overflows:
      DO 59 IT=1,M
      C = 0.
      DO 57 J=1,M
57    C = MAX(C,ABS(F(I,J,IT)))
      IF(C .LT. FPMAX) GO TO 59
      C = SMALLR
      DO 58 J=1,M
         ZI(IT,J) = ZI(IT,J) * C
         ZM(IT,J) = ZM(IT,J) * C
         DO 58 L=1,III
            I = NF + (L-1)*NT
58       F(I,J,IT) = F(I,J,IT) * C
59    CONTINUE
!      write(777,*) I,dble(ZI)
60    CONTINUE
!      write(777,*)
!      write(777,*)
C
      NT = -NT
      IF(NT.LE.0) GO TO 3
      DO 65 J=1,M
      DO 64 IT=1,M
      ZM(IT,J) = 0.0
64    ZI(IT,J) = 0.0
65    ZI(J,J) = 1E-10 * H**(Ql(J)+1) / EXP(0.5 * FACT(Ql(J)+1))
C      inner integration,  zi from 1 to mam
      NF = 1
      NO = NT * (MAM-1) + 1
      GO TO 40
3     COUNT = COUNT + 1
C   now calc. derivatives at matching pt. inner(zm) & outer(zp)
      DO 80 J=1,M
         MAT(J,NR) = 0.0
         MAT(J+M,NR) = 0.0
      DO 75 IT=1,M
      DEL =147.0*F(MAM,J,IT)-360.0*F(MAM-1,J,IT)+450.0*F(MAM-2,J,IT)
     1 -400.0*F(MAM-3,J,IT)+225.0*F(MAM-4,J,IT)
     2 -72*F(MAM-5,J,IT)+10.0*F(MAM-6,J,IT)
      ZM(IT,J) = 1  * DEL / (60.0 * H)
      DEL =147.0*F(MAP,J,IT)-360.0*F(MAP+1,J,IT)+450.0*F(MAP+2,J,IT)
     1 -400.0*F(MAP+3,J,IT)+225.0*F(MAP+4,J,IT)
     2 -72*F(MAP+5,J,IT)+10.0*F(MAP+6,J,IT)
      ZP(IT,J) =(-1)* DEL / (60.0 * H)
!      WRITE(KO,*) j,it,zm(it,j)/f(mam,j,it)-zp(it,j)/f(map,j,it)
      MAT(J,IT) = F(MAM,J,IT)
      MAT(J,IT+M) = -F(MAP,J,IT)
      MAT(J+M,IT) = ZM(IT,J)
75    MAT(J+M,IT+M) = - ZP(IT,J)
80    CONTINUE
      DO 78 IT=1,MM2
78    MAT(M+MC,IT) = 0.0
      MAT(M+MC,M+MC) = 1.0
      MAT(M+MC,NR)  = 1.0
c     do 81 j=1,nr
c     WRITE(KO,*) j,(mat(j,it),it=1,nr)
c81        continue
       RMAT=dble(MAT)
       CALL GAUSSr(2*M,MM2,RMAT,SING,DET,SMALL,.false.)
       IF(SING) GO TO 200
       DO 88 L=1,M
          DO I=1,NP
          PSI(I,L) = 0.0
          enddo
       DO 88 IT=1,M
       DO 87 I=2,MAM
87     PSI(I,L) = PSI(I,L) + F(I-1,L,IT) * RMAT(IT,NR)
       DO 88 I=MAP,NP
88     PSI(I,L) = PSI(I,L) + F(I,L,IT) * RMAT(IT+M,NR)
       DERIN = 0.0
       DEROUT= 0.0
       DO 90 IT=1,M
       DERIN = DERIN + dble(ZM(IT,MC)) * RMAT(IT,NR)
       DEROUT= DEROUT+ dble(ZP(IT,MC)) * RMAT(IT+M,NR)
90     CONTINUE
       NCO = 1
       DO 91 I=3,NP
C 91     IF(PSI(I,MC)*PSI(I-1,MC).LT.0.) NCO = NCO + 1
       IF(dble(PSI(I,MC)).GT.0 .AND. dble(PSI(I-1,MC)).LT.0.) NCO=NCO+1
91     IF(dble(PSI(I,MC)).LT.0 .AND. dble(PSI(I-1,MC)).GT.0.) NCO=NCO+1
       DEN = 0.0
       DO  96 J=1,M
       DO  96 L=1,M
       DO 92 IMAX=NP,1,-1
92          IF(ABS(CCMAT(L,J,IMAX)).GT.SMALLR .AND.
     X         ABS(PSI(IMAX,L)).GT.SMALLQ) GO TO 925
925         DO  93 I=2,IMAX
!            write(778,*) I,dble(CCMAT(J,L,I))
93          DEN = DEN + dble(PSI(I,J))*dble(CCMAT(J,L,I))*dble(PSI(I,L))
!         write(778,*)
!         write(778,*)
         IF(THETA.EQ.0.0 .OR. J.NE.L) GO TO 96
         DO 95 I=2,NP
95       DEN = DEN - dble(PSI(I,J))*THETA*dble(PSI(I,L))
96    CONTINUE
       DEL =(dble(DEROUT) - dble(DERIN)) / abs(PSI(MAP,MC))
       IF(PCON.GE.3) WRITE(KO,217) P,DEL,NCO,MAM
C 217    format(e14.6,e15.3,i6)
  217    FORMAT(1X,F13.6,F15.8,2I6)
       POLD = P
       IF(COUNT.GT.MAXC) GO TO 210
       IF(abs(P).lt.1d-4) GO TO 210
!       write(*,*) 'NCO',NCO,'NODES',NODES
       IF(NCO.NE.NODES) GO TO 6
       IF(ABS(DEL).LT.EPS) GO TO 7
         IF(DEN.EQ.0.0) GO TO 220
       Q =DEL * abs(PSI(MAP,MC))**2 / (DEN * H)
 !      write(*,*) 'Q',q,p
       IF(P.LT.Q) Q=P
       IF(P.LT.-Q) Q = -0.5*P
       P = P+Q
       BC = 1
       Q = 2*Q
99     IF(THETA.NE.0.0) GO TO 102
       if(P>POLD*1.8) GO TO 102  ! recalculate matching point MAM for big potl changes!
       if(P<POLD*0.6) GO TO 102  ! recalculate matching point MAM for big potl changes!
       GO TO 103
6     IF(BC.NE.0) GO TO 61
      CC = 1
      IF(NODES.LT.NCO) CC = -1
      IF(DEN.LT.0.0) CC = -CC
      IF(CC.LT.0) CCP = -IABS(CCP)
      IF(CCP.GT.0) PP = P
      IF(CCP.LT.0) PP = 0.5*PP
      P = PP*CC + P
      GO TO 99
61    Q = 0.5*Q
      P = P - 0.5*Q
      GO TO 99
7     Q = 0.0
      DO 700 J=1,M
      DO 700 I=1,NP
700   Q = Q + PSI(I,J)**2
         PIC = 1.0
      Q = SIGN(PIC,dble(PSI(3,MC))) / SQRT(Q * H)
      DO 710 J=1,M
      DO 710 I=1,NP
710   PSI(I,J) = PSI(I,J) * Q
      IF(PCON.EQ.0) RETURN
      IF(MOD(PCON,2).EQ.0)      RETURN
      DO 996 J=1,M
!      WRITE(KO,992) J,Ql(J),0.0,0
!  992 FORMAT(/' For Channel #',I3,' with l =',I2,',  j =',F5.1,
!     &  ', around core state #',I3,             //,'   R    Y(R)')
  996 WRITE(KO,998) ((I-1)*H,PSI(I,J),I=1,NP,MR)
  998 FORMAT(6(1X,0P,F6.2,2X,1P,E11.4,', '))
      RETURN
200   IFAIL = 1
      WRITE(KO,202) DET
202   FORMAT(' ***** FAIL IN EIGCC : DETERMINANT =',1P,2E12.4)
      STOP
210   IFAIL = 2
      WRITE(KO,212) DEL,MAXC
212   FORMAT(' ***** FAIL IN EIGCC : DISCREPANCY =',1P,E12.4,
     &       ' EVEN AFTER',I3,' ITERATIONS'/)
      STOP
220   IFAIL = 3
      WRITE(KO,222)
222   FORMAT(' ***** FAIL IN EIGCC : VARIABLE POTENTIAL PART ZERO'/)
      STOP
      END
      SUBROUTINE GAUSSR(N,NR,A,SING,DET,EPS,SHOW)
C
C    solve by Gaussian elimination sum(j): A(i,j).P(j) = A(i,N+1)
C             where A(j) is left in A(j,N+1)
C
c	use io
       IMPLICIT REAL*8(A-H,O-Z)
       PARAMETER(M = 1)
       REAL*8 A(NR,N+M),DET,RA
       LOGICAL SING,SHOW
       SING  = .FALSE.
       ko=6
      NPM = N + M
C      DO 201 I=1,N
C201   IF(SHOW) WRITE(KO,402) (            A(I,J)  ,J=1,NPM)
      DET = 1
      DO 9 K = 1, N
      DET = DET * A(K,K)
      IF (ABS (A(K,K)) .GT. EPS ) GO TO 5
         SING = .TRUE.
         WRITE(KO,3) K,DET
    3    FORMAT(//' THE MATRIX IS SINGULAR AT',I3,'  determinant is ',
     &  E16.8/)
      DO 401 I=1,N
401   IF(SHOW) WRITE(KO,402) (            A(I,J)  ,J=1,NPM)
402   FORMAT( 1X,20F6.3/(1X,20F6.3))
         RETURN
    5 KP1 = K + 1
         RA = 1.0/A(K,K)
      DO 6 J = KP1, NPM
    6 A(K,J) = A(K,J) * RA
      A(K,K) = 1
      DO 9 I = 1, N
      IF (I .EQ. K  .OR. ABS (A(I,K)) .EQ. 0) GO TO 9
Cdir$ ivdep
         DO 8 J = KP1, NPM
    8    A(I,J) = A(I,J) - A(I,K)*A(K,J)
         A(I,K) = 0
    9 CONTINUE
      IF(SHOW) WRITE(KO,15 ) DET
15    FORMAT(/' The determinant is ',E16.8)
      RETURN
      END
      
C  SUBROUTINE NEEDED ON NON-IBM MACHINES:
      SUBROUTINE ABEND(IRC)
      STOP
      END
