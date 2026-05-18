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
      IF (IEGEN .EQ. 0) THEN
         DO I = 1,N
            DO J = 1,N
               IF (I .EQ. J) THEN
                  U(I,J) = 1.0
               ELSE
                  U(I,J) = 0.
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C
      NR=0
      IF (N - 1 .LE. 0) THEN
         RETURN
      ENDIF
C
C     SCAN FOR LARGEST OFF DIAGONAL ELEMENT IN EACH ROW
C     X(I) CONTAINS LARGEST ELEMENT IN ITH ROW
C     IQ(I) HOLDS SECOND SUBSCRIPT DEFINING POSITION OF ELEMENT
C
      NMI1 = N - 1
      DO I = 1,NMI1
         X(I) = 0.
         IPL1 = I + 1
         DO J = IPL1,N
            IF (X(I) .LT. ABS (H(J,I))) THEN
               X(I) = ABS (H(J,I))
               IQ(I) = J
            ENDIF
         ENDDO
      ENDDO
C
C     SET INDICATOR FOR SHUT-OFF.RAP=2**-27,NR = NO. OF ROTATIONS
      RAP=7.450580596d-9
      HDTEST=1.0d100
C
C     FIND MAXIMUM OF X(I) S FOR PIVOT ELEMENT AND
C     TEST FOR END OF PROBLEM
C
   40 DO I = 1,NMI1
         IF (I .EQ. 1 .OR. XMAX .LT. X(I)) THEN
            XMAX = X(I)
            JPIV = IQ(I)
            IPIV = I
         ENDIF
      ENDDO
C
C     IS MAX. X(I) EQUAL TO ZERO, IF LESS THAN HDTEST, REVISE HDTEST
      IF (XMAX .LE. 0.0) THEN
         RETURN
      ENDIF
      IF (HDTEST .LE. 0.0 .OR. XMAX .LE. HDTEST) THEN
   90    HDIMIN = ABS (H(1,1))
         DO I = 2,N
            IF (HDIMIN .LT. ABS (H(I,I))) THEN
               HDIMIN = ABS (H(I,I))
            ENDIF
         ENDDO
         HDTEST = HDIMIN*RAP
         IF (HDTEST .GE. XMAX) THEN
            RETURN
         ENDIF
      ENDIF
      NR = NR + 1
C
C     COMPUTE TANGENT, SINE AND COSINE, H(I,I),H(J,J)
      TANG = SIGN (TWO,(H(IPIV,IPIV)-H(JPIV,JPIV)))*H(JPIV,IPIV)/(ABS (H
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
      IF (H(IPIV,IPIV) .LT. H(JPIV,JPIV)) THEN
         HTEMP = H(IPIV,IPIV)
         H(IPIV,IPIV) = H(JPIV,JPIV)
         H(JPIV,JPIV) = HTEMP
         HTEMP = SIGN (ONE, -SINE)*COSINE
         COSINE = ABS (SINE)
         SINE = HTEMP
      ENDIF
C
C     INSPECT THE IQS BETWEEN I + 1 AND N-1 TO DETERMINE
C     WHETHER A NEW MAXIMUM VALUE SHOULD BE COMPUTED SINCE
C     THE PRESENT MAXIMUM IS IN THE I OR J ROW.
C
      DO I = 1,NMI1
         IF (I .EQ. IPIV .OR. I .EQ. JPIV) THEN
            CYCLE
         ENDIF
         IF (IQ(I) .EQ. IPIV .OR. IQ(I) .EQ. JPIV) THEN
            K = IQ(I)
            HTEMP=H(K,I)
            H(K,I) = 0.
            IPL1 = I + 1
            X(I) = 0.
C
C     SEARCH IN DEPLETED ROW FOR NEW MAXIMUM
C
            DO J = IPL1,N
               IF (X(I) .LT. ABS (H(J,I))) THEN
                  X(I) = ABS (H(J,I))
                  IQ(I) = J
               ENDIF
            ENDDO
            H(K,I) = HTEMP
         ENDIF
      ENDDO
C
      X(IPIV) = 0.
      X(JPIV) = 0.
C
C     CHANGE THE OTHER ELEMENTS OF H
C
      DO I = 1,N
         IF (I .EQ. IPIV .OR. I .EQ. JPIV) THEN
            CYCLE
         ENDIF
         IF (I .LT. IPIV) THEN
            HTEMP = H(IPIV,I)
            H(IPIV,I) = COSINE*HTEMP + SINE * H(JPIV,I)
            IF (X(I) .LT. ABS (H(IPIV,I))) THEN
               X(I) = ABS (H(IPIV,I))
               IQ(I) = IPIV
            ENDIF
            H(JPIV,I) = -SINE*HTEMP + COSINE*H(JPIV,I)
            IF (X(I) .LT. ABS (H(JPIV,I))) THEN
               X(I) = ABS (H(JPIV,I))
               IQ(I) = JPIV
            ENDIF
         ELSE IF (I .LT. JPIV) THEN
            HTEMP = H(I,IPIV)
            H(I,IPIV) = COSINE*HTEMP + SINE*H(JPIV,I)
            IF (X(IPIV) .LT. ABS (H(I,IPIV))) THEN
               X(IPIV) = ABS (H(I,IPIV))
               IQ(IPIV) = I
            ENDIF
            H(JPIV,I) = -SINE*HTEMP + COSINE*H(JPIV,I)
            IF (X(I) .LT. ABS ( H(JPIV,I))) THEN
               X(I) = ABS (H(JPIV,I))
               IQ(I) = JPIV
            ENDIF
         ELSE
            HTEMP = H(I,IPIV)
            H(I,IPIV) = COSINE*HTEMP + SINE*H(I,JPIV)
            IF (X(IPIV) .LT. ABS ( H(I,IPIV))) THEN
               X(IPIV) = ABS (H(I,IPIV))
               IQ(IPIV) = I
            ENDIF
            H(I,JPIV) = -SINE*HTEMP + COSINE*H(I,JPIV)
            IF (X(JPIV) .LT. ABS ( H(I,JPIV))) THEN
               X(JPIV) = ABS (H(I,JPIV))
               IQ(JPIV) = I
            ENDIF
         ENDIF
      ENDDO
C
C     TEST FOR COMPUTATION OF EIGENVECTORS
C
      IF (IEGEN .EQ. 0) THEN
         DO I = 1,N
            HTEMP = U(IPIV,I)
            U(IPIV,I) = COSINE*HTEMP + SINE*U(JPIV,I)
            U(JPIV,I) = -SINE*HTEMP + COSINE*U(JPIV,I)
         ENDDO
      ENDIF
      GO TO 40
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

      do j = 1,n
         ipivot(j) = 0
      enddo

      do i = 1,n
         xmax = 0.0
         do j = 1,n
            if (ipivot(j) .ne. 1) then
               do k=1,n
                  if (ipivot(k) .eq. 1) then
                     ! do nothing
                  else if (ipivot(k) .gt. 1) then
                     go to 90
                  else
                     if (abs(a(j,k)) .ge. 1.0e-20) then
                        xjk = abs(a(j,k)) + abs(b(j,k))
                        if (xmax .lt. xjk) then
                           irow = j
                           icolum = k
                           xmax = xjk
                        endif
                     endif
                  endif
               enddo
            endif
         enddo

         ipivot(icolum) = ipivot(icolum) + 1
         if (irow .ne. icolum) then
            do l = 1,n
               swap = a(irow,l)
               swapi = b(irow,l)
               a(irow,l) = a(icolum,l)
               b(irow,l) = b(icolum,l)
               b(icolum,l) = swapi
               a(icolum,l) = swap
            enddo
         endif

         indexa(i) = irow
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

         do l = 1,n
            temp = a(icolum,l) * tempr - b(icolum,l) * tempi
            b(icolum,l) = a(icolum,l) * tempi + b(icolum,l) * tempr
            a(icolum,l) = temp
         enddo

         do l1 = 1,n
            if (l1 .ne. icolum) then
               tempa = a(l1,icolum)
               tempb = b(l1,icolum)
               a(l1,icolum) = 0.0
               b(l1,icolum) = 0.0

               do l = 1,n
                  b(l1,l) = b(l1,l) - a(icolum,l) * tempb
     1                      - b(icolum,l) * tempa
                  a(l1,l) = a(l1,l) - a(icolum,l) * tempa
     1                      + b(icolum,l) * tempb
               enddo
            endif
         enddo
      enddo

      do i = 1,n
         l = n + 1 - i
         if (indexa(l) .ne. indexb(l)) then
            jrow = indexa(l)
            jcolum = indexb(l)

            do k = 1,n
               swap = a(k,jrow)
               swapi = b(k,jrow)
               a(k,jrow) = a(k,jcolum)
               b(k,jrow) = b(k,jcolum)
               a(k,jcolum) = swap
               b(k,jcolum) = swapi
            enddo
         endif
      enddo

            
   90 deallocate(indexa,indexb,ipivot)
      return
      end
