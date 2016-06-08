

      SUBROUTINE HDIAG(H,N,NSMAX,IEGEN,U,NR)
C     MIHDI3,FORTRAN 22 DIAGONALIZATION OF A REAL SYMMETRIC MATRIX BY
C        THE JACOBI METHOD.
C     MAY 19, 1959
C     CALLING SEQUENCE FOR DIAGONALIZATION
C            CALL  HDIAG(H, N, IEGEN, U, NR)
C     IEGEN MUST BE SET UNEQUAL TO ZERO IF ONLY EIGENVALUES ARE
C             TO BE COMPUTED.
C
C     N IS THE ORDER OF THE MATRIX.  H.
C            WHERE H IS THE ARRAY TO BE DIAGONALIZED.
C     IEGEN MUST BE SET EQUAL TO ZERO IF EIGENVALUES AND EIGENVECTORS
C           ARE TO BE COMPUTED.
C
C     U IS THE UNITARY MATRIX USED FOR FORMATION OF THE EIGENVECTORS.
C
C     NR IS THE NUMBER OF ROTATIONS.
C
C     A DIMENSION STATEMENT MUST BE INSERTED IN THE SUBROUTINE.
C     DIMENSION H(N,N), U(N,N), X(N), IQ(N)
C
C
C     THE SUBROUTINE OPERATES ONLY ON THE ELEMENTS OF H THAT ARE TO THE
C             RIGHT OF THE MAIN DIAGONAL. THUS, ONLY A TRIANGULAR
C             SECTION NEED BE STORED IN THE ARRAY H.
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION H(NSMAX,NSMAX), U(NSMAX,NSMAX), X(NSMAX), IQ(NSMAX)
      ONE = 1d0
      TWO = 2d0
      IF(IEGEN) 15,10,15
   10 DO 14 I = 1,N
      DO 14 J = 1,N
      IF(I - J) 12,11,12
   11 U(I,J) = 1.0
      GO TO 14
   12 U(I,J) = 0.
   14 CONTINUE
C
   15 NR=0
      IF(N - 1) 1000,1000,17
C
C     SCAN FOR LARGEST OFF DIAGONAL ELEMENT IN EACH ROW
C     X(I) CONTAINS LARGEST ELEMENT IN ITH ROW
C     IQ(I) HOLDS SECOND SUBSCRIPT DEFINING POSITION OF ELEMENT
C
   17 NMI1 = N - 1
      DO 30 I = 1,NMI1
      X(I) = 0.
      IPL1 = I + 1
      DO 30 J = IPL1,N
      IF(X(I) - ABS (H(J,I))) 20,20,30
   20 X(I) = ABS (H(J,I))
      IQ(I) = J
   30 CONTINUE
C
C     SET INDICATOR FOR SHUT-OFF.RAP=2**-27,NR = NO. OF ROTATIONS
      RAP=7.450580596d-9
      HDTEST=1.0d100
C
C     FIND MAXIMUM OF X(I) S FOR PIVOT ELEMENT AND
C     TEST FOR END OF PROBLEM
C
   40 DO 70 I = 1,NMI1
      IF(I - 1) 60,60,45
   45 IF(XMAX-X(I)) 60,70,70
   60 XMAX = X(I)
      JPIV = IQ(I)
      IPIV = I
   70 CONTINUE
C
C     IS MAX. X(I) EQUAL TO ZERO, IF LESS THAN HDTEST, REVISE HDTEST
      IF(XMAX) 1000,1000,80
   80 IF (HDTEST) 90,90,85
   85 IF(XMAX - HDTEST) 90,90,148
   90 HDIMIN = ABS (H(1,1))
      DO 110 I = 2,N
      IF(HDIMIN - ABS (H(I,I))) 110,110,100
  100 HDIMIN = ABS (H(I,I))
  110 CONTINUE
C
      HDTEST = HDIMIN*RAP
C
C     RETURN IF MAX.H(J,I)LESS THAN(2**-27)ABSF(H(K,K)-MIN)
      IF (HDTEST - XMAX) 148,1000,1000
  148 NR = NR + 1
C
C     COMPUTE TANGENT, SINE AND COSINE, H(I,I),H(J,J)
  150 TANG = SIGN (TWO,(H(IPIV,IPIV)-H(JPIV,JPIV)))*H(JPIV,IPIV)/(ABS (H
     1(IPIV,IPIV)-H(JPIV,JPIV))+SQRT ((H(IPIV,IPIV)-H(JPIV,JPIV))**2+4.0
     2*H(JPIV,IPIV)**2))
      COSINE = 1.0/SQRT (1.0+TANG**2)
      SINE = TANG*COSINE
      HII = H(IPIV,IPIV)
      H(IPIV,IPIV) = COSINE**2*(HII+TANG*(2.*H(JPIV,IPIV)+TANG*H(JPIV,JP
     1IV)))
      H(JPIV,JPIV)=COSINE**2*(H(JPIV,JPIV)-TANG*(2.*H(JPIV,IPIV)-TANG*H
     1II))
      H(JPIV,IPIV)= 0.
C
C     PSEUDO RANK THE EIGENVALUES
C     ADJUST SINE AND COS FOR COMPUTATION OF H(IK) AND U(IK)
      IF(H(IPIV,IPIV) -H(JPIV,JPIV)) 152,153,153
  152 HTEMP = H(IPIV,IPIV)
      H(IPIV,IPIV) = H(JPIV,JPIV)
      H(JPIV,JPIV) = HTEMP
      HTEMP = SIGN (ONE, -SINE)*COSINE
      COSINE = ABS (SINE)
      SINE = HTEMP
  153 CONTINUE
C
C     INSPECT THE IQS BETWEEN I + 1 AND N-1 TO DETERMINE
C     WHETHER A NEW MAXIMUM VALUE SHOULD BE COMPUTED SINCE
C     THE PRESENT MAXIMUM IS IN THE I OR J ROW.
C
      DO 350 I = 1,NMI1
      IF(I-IPIV) 210,350,200
  200 IF(I-JPIV) 210,350,210
  210 IF(IQ(I)-IPIV) 230,240,230
  230 IF(IQ(I)-JPIV)350,240,350
  240 K = IQ(I)
  250 HTEMP=H(K,I)
      H(K,I) = 0.
      IPL1 = I + 1
      X(I) = 0.
C
C     SEARCH IN DEPLETED ROW FOR NEW MAXIMUM
C
      DO 320 J = IPL1,N
      IF(X(I)-ABS (H(J,I))) 300,300,320
  300 X(I) = ABS (H(J,I))
      IQ(I) = J
  320 CONTINUE
      H(K,I) = HTEMP
  350 CONTINUE
C
      X(IPIV) = 0.
      X(JPIV) = 0.
C
C     CHANGE THE OTHER ELEMENTS OF H
C
      DO 530 I = 1,N
C
      IF(I-IPIV) 370,530,420
  370 HTEMP = H(IPIV,I)
      H(IPIV,I) = COSINE*HTEMP + SINE * H(JPIV,I)
      IF(X(I) - ABS (H(IPIV,I)))380,390,390
  380 X(I) = ABS (H(IPIV,I))
      IQ(I) = IPIV
  390 H(JPIV,I) = -SINE*HTEMP + COSINE*H(JPIV,I)
      IF(X(I) - ABS (H(JPIV,I))) 400,530,530
  400 X(I) = ABS (H(JPIV,I))
      IQ(I) = JPIV
      GO TO 530
C
  420 IF(I-JPIV) 430,530,480
  430 HTEMP = H(I,IPIV)
      H(I,IPIV) = COSINE*HTEMP + SINE*H(JPIV,I)
      IF(X(IPIV) - ABS (H(I,IPIV)) ) 440,450,450
  440 X(IPIV) = ABS (H(I,IPIV))
      IQ(IPIV) = I
  450 H(JPIV,I) = -SINE*HTEMP + COSINE*H(JPIV,I)
      IF(X(I) - ABS ( H(JPIV,I)) ) 400,530,530
C
  480 HTEMP = H(I,IPIV)
      H(I,IPIV) = COSINE*HTEMP + SINE*H(I,JPIV)
      IF(X(IPIV) - ABS ( H(I,IPIV)) ) 490,500,500
  490 X(IPIV) = ABS (H(I,IPIV))
      IQ(IPIV) = I
  500 H(I,JPIV) = -SINE*HTEMP + COSINE*H(I,JPIV)
      IF(X(JPIV) - ABS ( H(I,JPIV))) 510,530,530
  510 X(JPIV) = ABS (H(I,JPIV))
      IQ(JPIV) = I
  530 CONTINUE
C
C     TEST FOR COMPUTATION OF EIGENVECTORS
C
      IF(IEGEN) 40,540,40
  540 DO 550 I = 1,N
      HTEMP = U(IPIV,I)
      U(IPIV,I) = COSINE*HTEMP + SINE*U(JPIV,I)
  550 U(JPIV,I) = -SINE*HTEMP + COSINE*U(JPIV,I)
      GO TO 40
 1000 RETURN
      END




      subroutine cmatin(a, b, n)
************************************************************************
c
c
c *** finds the inverse of complex matrix h=(rea,aima)=(a,b)
c
c *** code provided by norman bardsley to kwon+tabakin s prog bopit
c *** n=no of elements in a assigned(and in the sub=dimension)
c *** matrix a,b(real+imag parts) destroyed by call
c *** the other arguments are dummy just to save space

      implicit real*8 (a-h, o-z)
      integer, allocatable :: indexa(:),indexb(:),ipivot(:)
      real*8 :: a(n,n),b(n,n)

      equivalence (irow,jrow), (icolum,jcolum)

      allocate(indexa(n),indexb(n),ipivot(n))

c ***
c *** detr=1.0
c *** rhl add to define tempi for real potential

      tempi = 0.

c *** deti=0.0

      do 10 j = 1,n
  10     ipivot(j) = 0

      do 70 i = 1,n
         xmax = 0.0
         do 35 j = 1,n
            if (ipivot(j) - 1) 15, 35, 15

  15        do 30 k=1,n
               if (ipivot(k) - 1) 20, 30, 90
  20           if (abs(a(j,k)) .lt. 1.0e-20) go to 30
               xjk = abs(a(j,k)) + abs(b(j,k))
               if (xmax - xjk)  25,30,30
  25           irow = j
               icolum = k
               xmax = xjk
  30        continue

  35  continue

      ipivot(icolum) = ipivot(icolum) + 1
      if (irow .eq. icolum) go to 50

c *** detr=-deti
c *** deti=-deti

      do 45 l = 1,n
         swap = a(irow,l)
         swapi = b(irow,l)
         a(irow,l) = a(icolum,l)
         b(irow,l) = b(icolum,l)
         b(icolum,l) = swapi
  45     a(icolum,l) = swap

  50  indexa(i) = irow
      indexb(i) = icolum
      pivotr = a(icolum,icolum)
      pivoti = b(icolum,icolum)

c *** temp=detr*pivotr-deti*pivoti
c *** deti=detr*pivoti+deti*pivotr
c *** detr=temp

      a(icolum,icolum) = 1.0
      b(icolum,icolum) = 0.0
      if (pivoti .eq. 0.0) tempr = 1.0/pivotr
      if (pivoti .ne. 0.0) then
         tempr =  pivotr/(pivotr * pivotr + pivoti * pivoti)
         tempi = -pivoti/(pivotr * pivotr + pivoti * pivoti)
      endif

      do 55 l = 1,n
         temp = a(icolum,l) * tempr - b(icolum,l) * tempi
         b(icolum,l) = a(icolum,l) * tempi + b(icolum,l) * tempr
  55     a(icolum,l) = temp

      do 70 l1 = 1,n
         if (l1 - icolum) 60, 70, 60
  60     tempa = a(l1,icolum)
         tempb = b(l1,icolum)
         a(l1,icolum) = 0.0
         b(l1,icolum) = 0.0

         do 65 l = 1,n
            b(l1,l) = b(l1,l) - a(icolum,l) * tempb
     $                - b(icolum,l) * tempa
            a(l1,l) = a(l1,l) - a(icolum,l) * tempa
     $                + b(icolum,l) * tempb
  65     continue
  70  continue

      do 85 i = 1,n
         l = n + 1 - i
         if (indexa(l) - indexb(l)) 75, 85, 75
  75     jrow = indexa(l)
         jcolum = indexb(l)

         do 80 k = 1,n
            swap = a(k,jrow)
            swapi = b(k,jrow)
            a(k,jrow) = a(k,jcolum)
            b(k,jrow) = b(k,jcolum)
            a(k,jcolum) = swap
            b(k,jcolum) = swapi
  80     continue
  85  continue

            
  90  deallocate(indexa,indexb,ipivot)
      return
      end
