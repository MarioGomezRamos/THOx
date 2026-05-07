! This file contains deprecated/unused subroutines extracted from the codebase.
! They are not currently compiled with the project.

! From whittaker.f90
      SUBROUTINE PHASES(ETA,LLMAX,SIGMA)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ETA,SIGMA(LLMAX)
      INTEGER MAX
      COMPLEX*16 S,CLGAMM,Z
      Z=CMPLX( 1D0,ETA)
      S=CLGAMM(Z)
      SIGMA(1)=S*(0D0,-1D0)
      IF (LLMAX.EQ.1) RETURN
      DO 1 MAX=2,LLMAX
    1 SIGMA(MAX)= ATAN2(ETA,(MAX-1D0))+SIGMA(MAX-1)
      RETURN
      END


! From nag2.f90
      subroutine f02agf(a,ia,n,rr,ri,vr,ivr,vi,ivi,intger,ifail)
c     mark 13 re-issue. nag copyright 1988.
c     mark 14a revised. ier-685 (dec 1989).
c
c     eigenvalues and eigenvectors of real unsymmetric matrix
c     1st august 1971
c
c     .. parameters ..
      character*6       srname
      parameter         (srname='f02agf')
c     .. scalar arguments ..
      integer           ia, ifail, ivi, ivr, n
c     .. array arguments ..
      double precision  a(ia,n), ri(n), rr(n), vi(ivi,n), vr(ivr,n)
      integer           intger(n)
c     .. local scalars ..
      double precision  c, d, machep, max, sum, term
      integer           i, ib, isave, j, k, l
c     .. local arrays ..
      character*1       p01rec(1)
c     .. external functions ..
      double precision  x02ajf
      integer           p01abf, x02bhf
      external          x02ajf, p01abf, x02bhf
c     .. external subroutines ..
      external          f01akf, f01apf, f01atf, f01auf, f02aqf
c     .. intrinsic functions ..
      intrinsic         abs, sqrt
c     .. executable statements ..
      isave = ifail
      ifail = 1
      machep = x02ajf()
      ib = x02bhf()
      call f01atf(n,ib,a,ia,k,l,rr)
      call f01akf(n,k,l,a,ia,intger)
      call f01apf(n,k,l,intger,a,ia,vr,ivr)
      call f01auf(n,k,l,n,rr,vr,ivr)
      call f02aqf(n,1,n,machep,a,ia,vr,ivr,rr,ri,intger,ifail)
      if (ifail.eq.0) go to 20
      ifail = p01abf(isave,ifail,srname,0,p01rec)
      return
   20 do 140 i = 1, n
         if (ri(i).eq.0.0d0) go to 60
         if (ri(i).gt.0.0d0) go to 100
         do 40 j = 1, n
            vr(j,i) = vr(j,i-1)
            vi(j,i) = -vi(j,i-1)
   40    continue
         go to 140
   60    do 80 j = 1, n
            vi(j,i) = 0.0d0
   80    continue
         go to 140
  100    do 120 j = 1, n
            vi(j,i) = vr(j,i+1)
  120    continue
  140 continue
      do 280 i = 1, n
         sum = 0.0d0
         max = 0.0d0
         do 180 j = 1, n
            if (abs(vr(j,i)).le.max) go to 160
            max = abs(vr(j,i))
  160       if (abs(vi(j,i)).le.max) go to 180
            max = abs(vi(j,i))
  180    continue
         do 200 j = 1, n
            vr(j,i) = vr(j,i)/max
            vi(j,i) = vi(j,i)/max
  200    continue
         max = 0.0d0
         do 240 j = 1, n
            term = vr(j,i)**2 + vi(j,i)**2
            sum = sum + term
            if (term.le.max) go to 220
            max = term
            c = vr(j,i)
            d = -vi(j,i)
  220       continue
  240    continue
         sum = sum*(c**2+d**2)
         sum = sqrt(sum)
         do 260 j = 1, n
            term = vr(j,i)
            vr(j,i) = (vr(j,i)*c-vi(j,i)*d)/sum
            vi(j,i) = (d*term+c*vi(j,i))/sum
  260    continue
  280 continue
      return
      end


! From nag2.f90
      subroutine f06paf( trans, m, n, alpha, a, lda, x, incx,
     $                   beta, y, incy )
c     mark 12 release. nag copyright 1986.
c     axp4 version for vector machines
c     .. entry points ..
      entry      dgemv ( trans, m, n, alpha, a, lda, x, incx,
     $                   beta, y, incy )
c     .. scalar arguments ..
      double precision   alpha, beta
      integer            incx, incy, lda, m, n
      character*1        trans
c     .. array arguments ..
      double precision   a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  dgemv  performs one of the matrix-vector operations
c
c     y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,
c
c  where alpha and beta are scalars, x and y are vectors and a is an
c  m by n matrix.
c
c  parameters
c  ==========
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   y := alpha*a*x + beta*y.
c
c              trans = 't' or 't'   y := alpha*a'*x + beta*y.
c
c              trans = 'c' or 'c'   y := alpha*a'*x + beta*y.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - double precision.
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - double precision array of dimension ( lda, n ).
c           before entry, the leading m by n part of the array a must
c           contain the matrix of coefficients.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, m ).
c           unchanged on exit.
c
c  x      - double precision array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
c           before entry, the incremented array x must contain the
c           vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - double precision.
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - double precision array of dimension at least
c           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
c           before entry with beta non-zero, the incremented array y
c           must contain the vector y. on exit, y is overwritten by the
c           updated vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c  -- do-loops unrolled on 20-november-1986.
c     peter mayes, nag central office.
c
c
c     .. parameters ..
      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     .. local scalars ..
      double precision   temp, temp1, temp2, temp3, temp4
      integer            i, info, iy, j, jx, kx, ky, lenx, leny, m4, n4
c     .. external subroutines ..
      external           f06aaz
c     .. intrinsic functions ..
      intrinsic          max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.(trans.eq.'n' .or. trans.eq.'n').and.
     $         .not.(trans.eq.'t' .or. trans.eq.'t').and.
     $         .not.(trans.eq.'c' .or. trans.eq.'c')      )then
         info = 1
      else if( m.lt.0 )then
         info = 2
      else if( n.lt.0 )then
         info = 3
      else if( lda.lt.max( 1, m ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      else if( incy.eq.0 )then
         info = 11
      end if
      if( info.ne.0 )then
         call f06aaz( 'f06paf/dgemv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     set  lenx  and  leny, the lengths of the vectors x and y, and set
c     up the start points in  x  and  y.
c
      if( (trans.eq.'n' .or. trans.eq.'n') )then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( lenx - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( leny - 1 )*incy
      end if
c
c     start the operations. in this version the inner loops are all
c     equivalent to axpy operations.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, leny
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, leny
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, leny
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, leny
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      jx = kx
      if( (trans.eq.'n' .or. trans.eq.'n') )then
c
c        form  y := alpha*a*x + y.
c
         if( incy.eq.1 )then
c**** u n r o l l   t o   d e p t h   4 ********************************
            n4 = 4*( n/4 )
            do 60, j = 1, n4, 4
               temp1 = alpha*x( jx )
               temp2 = alpha*x( jx + incx )
               temp3 = alpha*x( jx + 2*incx )
               temp4 = alpha*x( jx + 3*incx )
               if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.zero.or.
     $             temp4.ne.zero )then
                  do 50, i = 1, m
                     y( i ) = ( ( ( ( y( i ) + temp1*a( i, j ) )
     $                        + temp2*a( i, j + 1 ) )
     $                        + temp3*a( i, j + 2 ) )
     $                        + temp4*a( i, j + 3 ) )
   50             continue
               end if
               jx = jx + 4*incx
   60       continue
c**** clean-up loop ****************************************************
            do 80, j = n4 + 1, n, 1
               temp = alpha*x( jx )
               if( temp.ne.zero )then
                  do 70, i = 1, m
                     y( i ) = y( i ) + temp*a( i, j )
   70             continue
               end if
               jx = jx + incx
   80       continue
         else
c**** u n r o l l   t o   d e p t h   4 ********************************
            n4 = 4*( n/4 )
            do 100, j = 1, n4, 4
               temp1 = alpha*x( jx )
               temp2 = alpha*x( jx + incx )
               temp3 = alpha*x( jx + 2*incx )
               temp4 = alpha*x( jx + 3*incx )
               if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.zero.or.
     $             temp4.ne.zero )then
                  iy = ky
                  do 90, i = 1, m
                     y( iy ) = ( ( ( ( y( iy ) + temp1*a( i, j ) )
     $                         + temp2*a( i, j + 1 ) )
     $                         + temp3*a( i, j + 2 ) )
     $                         + temp4*a( i, j + 3 ) )
                     iy = iy + incy
   90             continue
               end if
               jx = jx + 4*incx
  100       continue
c**** clean-up loop ****************************************************
            do 120, j = n4 + 1, n, 1
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy = ky
                  do 110, i = 1, m
                     y( iy ) = y( iy ) + temp*a( i, j )
                     iy = iy + incy
  110             continue
               end if
               jx = jx + incx
  120       continue
         end if
      else
c
c        form  y := alpha*a'*x + y.
c
         if( incy.eq.1 )then
c**** u n r o l l   t o   d e p t h   4 ********************************
            m4 = 4*( m/4 )
            do 140, j = 1, m4, 4
               temp1 = alpha*x( jx )
               temp2 = alpha*x( jx + incx )
               temp3 = alpha*x( jx + 2*incx )
               temp4 = alpha*x( jx + 3*incx )
               if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.zero.or.
     $             temp4.ne.zero )then
                  do 130, i = 1, n
                     y( i ) = ( ( ( ( y( i ) + temp1*a( j, i ) )
     $                        + temp2*a( j + 1, i ) )
     $                        + temp3*a( j + 2, i ) )
     $                        + temp4*a( j + 3, i ) )
  130             continue
               end if
               jx = jx + 4*incx
  140       continue
c**** clean-up loop ****************************************************
            do 160, j = m4 + 1, m, 1
               temp = alpha*x( jx )
               if( temp.ne.zero )then
                  do 150, i = 1, n
                     y( i ) = y( i ) + temp*a( j, i )
  150             continue
               end if
               jx = jx + incx
  160       continue
         else
c**** u n r o l l   t o   d e p t h   4 ********************************
            m4 = 4*( m/4 )
            do 180, j = 1, m4, 4
               temp1 = alpha*x( jx )
               temp2 = alpha*x( jx + incx )
               temp3 = alpha*x( jx + 2*incx )
               temp4 = alpha*x( jx + 3*incx )
               if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.zero.or.
     $             temp4.ne.zero )then
                  iy = ky
                  do 170, i = 1, n
                     y( iy ) = ( ( ( ( y( iy ) + temp1*a( j, i ) )
     $                         + temp2*a( j + 1, i ) )
     $                         + temp3*a( j + 2, i ) )
     $                         + temp4*a( j + 3, i ) )
                     iy = iy + incy
  170             continue
               end if
               jx = jx + 4*incx
  180       continue
c**** clean-up loop ****************************************************
            do 200, j = m4 + 1, m, 1
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy = ky
                  do 190, i = 1, n
                     y( iy ) = y( iy ) + temp*a( j, i )
                     iy = iy + incy
  190             continue
               end if
               jx = jx + incx
  200       continue
         end if
      end if
c
      return
c
c     end of f06paf (dgemv ).
c
      end


! From nag2.f90
      subroutine f06pjf( uplo, trans, diag, n, a, lda, x, incx )
c     mark 12 release. nag copyright 1986.
c     axp4 version for vector machines
c     .. entry points ..
      entry      dtrsv ( uplo, trans, diag, n, a, lda, x, incx )
c     .. scalar arguments ..
      integer            incx, lda, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      double precision   a( lda, * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  dtrsv  solves one of the systems of equations
c
c     a*x = b,   or   a'*x = b,
c
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular matrix.
c
c  no test for singularity or near-singularity is included in this
c  routine. such tests must be performed before calling this routine.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the equations to be solved as
c           follows:
c
c              trans = 'n' or 'n'   a*x = b.
c
c              trans = 't' or 't'   a'*x = b.
c
c              trans = 'c' or 'c'   a'*x = b.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  a      - double precision array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular matrix and the strictly lower triangular part of
c           a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular matrix and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced either, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c  x      - double precision array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element right-hand side vector b. on exit, x is overwritten
c           with the solution vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c  -- do-loops unrolled on 20-november-1986.
c     peter mayes, nag central office.
c
c
c     .. parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      double precision   temp1, temp2, temp3, temp4
      integer            i, info, ix, j, jx, kx, n4
      logical            nounit
c     .. external subroutines ..
      external           f06aaz
c     .. intrinsic functions ..
      intrinsic          max, mod
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.(uplo .eq.'u' .or. uplo .eq.'u').and.
     $         .not.(uplo .eq.'l' .or. uplo .eq.'l')      )then
         info = 1
      else if( .not.(trans.eq.'n' .or. trans.eq.'n').and.
     $         .not.(trans.eq.'t' .or. trans.eq.'t').and.
     $         .not.(trans.eq.'c' .or. trans.eq.'c')      )then
         info = 2
      else if( .not.(diag .eq.'u' .or. diag .eq.'u').and.
     $         .not.(diag .eq.'n' .or. diag .eq.'n')      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( lda.lt.max( 1, n ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      end if
      if( info.ne.0 )then
         call f06aaz( 'f06pjf/dtrsv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      nounit = (diag.eq.'n' .or. diag.eq.'n')
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the inner loops are all
c     equivalent to axpy operations.
c
      if( (trans.eq.'n' .or. trans.eq.'n') )then
c
c        form  x := inv( a )*x.
c
         if( (uplo.eq.'u' .or. uplo.eq.'u') )then
            if( incx.eq.1 )then
c**** u n r o l l   t o   d e p t h  4 ********** l o o p  a ***********
               n4 = mod( n, 4 ) + 1
               do 20, j = n, n4, -4
c**** 4 by 4 triangle at bottom of block *******************************
                  if( nounit )
     $               x( j ) = x( j )/a( j, j )
                  x( j - 1 ) = x( j - 1 ) - x( j )*a( j - 1, j )
                  if( nounit )
     $               x( j - 1 ) = x( j - 1 )/a( j - 1, j - 1 )
                  x( j - 2 ) = x( j - 2 ) - x( j )*a( j - 2, j ) -
     $                         x( j - 1 )*a( j - 2, j - 1 )
                  if( nounit )
     $               x( j - 2 ) = x( j - 2 )/a( j - 2, j - 2 )
                  x( j - 3 ) = x( j - 3 ) - x( j )*a( j - 3, j ) -
     $                         x( j - 1 )*a( j - 3, j - 1 ) - x( j - 2 )
     $                         *a( j - 3, j - 2 )
                  if( nounit )
     $               x( j - 3 ) = x( j - 3 )/a( j - 3, j - 3 )
c**** unrolled main loop ***********************************************
                  temp1 = x( j )
                  temp2 = x( j - 1 )
                  temp3 = x( j - 2 )
                  temp4 = x( j - 3 )
                  if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.
     $                zero.or.temp4.ne.zero )then
                     do 10, i = j - 4, 1, -1
                        x( i ) = ( ( ( ( x( i ) - temp1*a( i, j ) )
     $                           - temp2*a( i, j - 1 ) )
     $                           - temp3*a( i, j - 2 ) )
     $                           - temp4*a( i, j - 3 ) )
   10                continue
                  end if
   20          continue
c**** left-overs on top left corner ************************************
               if( n4.ge.2 )then
                  if( nounit )
     $               x( n4 - 1 ) = x( n4 - 1 )/a( n4 - 1, n4 - 1 )
               end if
               if( n4.ge.3 )then
                  x( n4 - 2 ) = x( n4 - 2 ) - x( n4 - 1 )
     $                          *a( n4 - 2, n4 - 1 )
                  if( nounit )
     $               x( n4 - 2 ) = x( n4 - 2 )/a( n4 - 2, n4 - 2 )
               end if
               if( n4.ge.4 )then
                  x( n4 - 3 ) = x( n4 - 3 ) - x( n4 - 1 )
     $                          *a( n4 - 3, n4 - 1 ) - x( n4 - 2 )
     $                          *a( n4 - 3, n4 - 2 )
                  if( nounit )
     $               x( n4 - 3 ) = x( n4 - 3 )/a( n4 - 3, n4 - 3 )
               end if
            else
               jx = kx + ( n - 1 )*incx
c**** u n r o l l   t o   d e p t h  4 ********** l o o p  b ***********
               n4 = mod( n, 4 ) + 1
               do 40, j = n, n4, -4
c**** 4 by 4 triangle at bottom of block *******************************
                  if( nounit )
     $               x( jx ) = x( jx )/a( j, j )
                  x( jx - incx ) = x( jx - incx ) - x( jx )
     $                             *a( j - 1, j )
                  if( nounit )
     $               x( jx - incx ) = x( jx - incx )/a( j - 1, j - 1 )
                  x( jx - 2*incx ) = x( jx - 2*incx ) - x( jx )
     $                               *a( j - 2, j ) - x( jx - incx )
     $                               *a( j - 2, j - 1 )
                  if( nounit )
     $               x( jx - 2*incx ) = x( jx - 2*incx )
     $                                  /a( j - 2, j - 2 )
                  x( jx - 3*incx ) = x( jx - 3*incx ) - x( jx )
     $                               *a( j - 3, j ) - x( jx - incx )
     $                               *a( j - 3, j - 1 ) -
     $                               x( jx - 2*incx )*a( j - 3, j - 2 )
                  if( nounit )
     $               x( jx - 3*incx ) = x( jx - 3*incx )
     $                                  /a( j - 3, j - 3 )
c**** unrolled main loop ***********************************************
                  temp1 = x( jx )
                  temp2 = x( jx - incx )
                  temp3 = x( jx - 2*incx )
                  temp4 = x( jx - 3*incx )
                  if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.
     $                zero.or.temp4.ne.zero )then
                     ix = jx - 3*incx
                     do 30, i = j - 4, 1, -1
                        ix = ix - incx
                        x( ix ) = ( ( ( ( x( ix ) - temp1*a( i, j ) )
     $                            - temp2*a( i, j - 1 ) )
     $                            - temp3*a( i, j - 2 ) )
     $                            - temp4*a( i, j - 3 ) )
   30                continue
                  end if
                  jx = jx - 4*incx
   40          continue
c**** left-overs on top left corner ************************************
               if( n4.ge.2 )then
                  if( nounit )
     $               x( jx ) = x( jx )/a( n4 - 1, n4 - 1 )
               end if
               if( n4.ge.3 )then
                  x( jx - incx ) = x( jx - incx ) - x( jx )
     $                             *a( n4 - 2, n4 - 1 )
                  if( nounit )
     $               x( jx - incx ) = x( jx - incx )/a( n4 - 2, n4 - 2 )
               end if
               if( n4.ge.4 )then
                  x( jx - 2*incx ) = x( jx - 2*incx ) - x( jx )
     $                               *a( n4 - 3, n4 - 1 ) -
     $                               x( jx - incx )*a( n4 - 3, n4 - 2 )
                  if( nounit )
     $               x( jx - 2*incx ) = x( jx - 2*incx )
     $                                  /a( n4 - 3, n4 - 3 )
               end if
            end if
         else
            if( incx.eq.1 )then
c**** u n r o l l   t o   d e p t h  4 ********** l o o p  c ***********
               n4 = 4*( n/4 )
               do 60, j = 1, n4, 4
c**** 4 by 4 triangle at bottom of block *******************************
                  if( nounit )
     $               x( j ) = x( j )/a( j, j )
                  x( j + 1 ) = x( j + 1 ) - x( j )*a( j + 1, j )
                  if( nounit )
     $               x( j + 1 ) = x( j + 1 )/a( j + 1, j + 1 )
                  x( j + 2 ) = x( j + 2 ) - x( j )*a( j + 2, j ) -
     $                         x( j + 1 )*a( j + 2, j + 1 )
                  if( nounit )
     $               x( j + 2 ) = x( j + 2 )/a( j + 2, j + 2 )
                  x( j + 3 ) = x( j + 3 ) - x( j )*a( j + 3, j ) -
     $                         x( j + 1 )*a( j + 3, j + 1 ) - x( j + 2 )
     $                         *a( j + 3, j + 2 )
                  if( nounit )
     $               x( j + 3 ) = x( j + 3 )/a( j + 3, j + 3 )
c**** unrolled main loop ***********************************************
                  temp1 = x( j )
                  temp2 = x( j + 1 )
                  temp3 = x( j + 2 )
                  temp4 = x( j + 3 )
                  if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.
     $                zero.or.temp4.ne.zero )then
                     do 50, i = j + 4, n
                        x( i ) = ( ( ( ( x( i ) - temp1*a( i, j ) )
     $                           - temp2*a( i, j + 1 ) )
     $                           - temp3*a( i, j + 2 ) )
     $                           - temp4*a( i, j + 3 ) )
   50                continue
                  end if
   60          continue
c**** left-overs on top left corner ************************************
               if( n4 + 1.le.n )then
                  if( nounit )
     $               x( n4 + 1 ) = x( n4 + 1 )/a( n4 + 1, n4 + 1 )
               end if
               if( n4 + 2.le.n )then
                  x( n4 + 2 ) = x( n4 + 2 ) - x( n4 + 1 )
     $                          *a( n4 + 2, n4 + 1 )
                  if( nounit )
     $               x( n4 + 2 ) = x( n4 + 2 )/a( n4 + 2, n4 + 2 )
               end if
               if( n4 + 3.le.n )then
                  x( n4 + 3 ) = x( n4 + 3 ) - x( n4 + 1 )
     $                          *a( n4 + 3, n4 + 1 ) - x( n4 + 2 )
     $                          *a( n4 + 3, n4 + 2 )
                  if( nounit )
     $               x( n4 + 3 ) = x( n4 + 3 )/a( n4 + 3, n4 + 3 )
               end if
            else
               jx = kx
c**** u n r o l l   t o   d e p t h  4 ********** l o o p  d ***********
               n4 = 4*( n/4 )
               do 80, j = 1, n4, 4
c**** 4 by 4 triangle at bottom of block *******************************
                  if( nounit )
     $               x( jx ) = x( jx )/a( j, j )
                  x( jx + incx ) = x( jx + incx ) - x( jx )
     $                             *a( j + 1, j )
                  if( nounit )
     $               x( jx + incx ) = x( jx + incx )/a( j + 1, j + 1 )
                  x( jx + 2*incx ) = x( jx + 2*incx ) - x( jx )
     $                               *a( j + 2, j ) - x( jx + incx )
     $                               *a( j + 2, j + 1 )
                  if( nounit )
     $               x( jx + 2*incx ) = x( jx + 2*incx )
     $                                  /a( j + 2, j + 2 )
                  x( jx + 3*incx ) = x( jx + 3*incx ) - x( jx )
     $                               *a( j + 3, j ) - x( jx + incx )
     $                               *a( j + 3, j + 1 ) -
     $                               x( jx + 2*incx )*a( j + 3, j + 2 )
                  if( nounit )
     $               x( jx + 3*incx ) = x( jx + 3*incx )
     $                                  /a( j + 3, j + 3 )
c**** unrolled main loop ***********************************************
                  temp1 = x( jx )
                  temp2 = x( jx + incx )
                  temp3 = x( jx + 2*incx )
                  temp4 = x( jx + 3*incx )
                  if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.
     $                zero.or.temp4.ne.zero )then
                     ix = jx + 3*incx
                     do 70, i = j + 4, n
                        ix = ix + incx
                        x( ix ) = ( ( ( ( x( ix ) - temp1*a( i, j ) )
     $                            - temp2*a( i, j + 1 ) )
     $                            - temp3*a( i, j + 2 ) )
     $                            - temp4*a( i, j + 3 ) )
   70                continue
                  end if
                  jx = jx + 4*incx
   80          continue
c**** left-overs on top left corner ************************************
               if( n4 + 1.le.n )then
                  if( nounit )
     $               x( jx ) = x( jx )/a( n4 + 1, n4 + 1 )
               end if
               if( n4 + 2.le.n )then
                  x( jx + incx ) = x( jx + incx ) - x( jx )
     $                             *a( n4 + 2, n4 + 1 )
                  if( nounit )
     $               x( jx + incx ) = x( jx + incx )/a( n4 + 2, n4 + 2 )
               end if
               if( n4 + 3.le.n )then
                  x( jx + 2*incx ) = x( jx + 2*incx ) - x( jx )
     $                               *a( n4 + 3, n4 + 1 ) -
     $                               x( jx + incx )*a( n4 + 3, n4 + 2 )
                  if( nounit )
     $               x( jx + 2*incx ) = x( jx + 2*incx )
     $                                  /a( n4 + 3, n4 + 3 )
               end if
            end if
         end if
      else
c
c        form  x := inv( a' )*x.
c
         if( (uplo.eq.'u' .or. uplo.eq.'u') )then
            if( incx.eq.1 )then
c**** u n r o l l   t o   d e p t h  4 ********** l o o p  e ***********
               n4 = 4*( n/4 )
               do 100, j = 1, n4, 4
c**** 4 by 4 triangle at bottom of block *******************************
                  if( nounit )
     $               x( j ) = x( j )/a( j, j )
                  x( j + 1 ) = x( j + 1 ) - x( j )*a( j, j + 1 )
                  if( nounit )
     $               x( j + 1 ) = x( j + 1 )/a( j + 1, j + 1 )
                  x( j + 2 ) = x( j + 2 ) - x( j )*a( j, j + 2 ) -
     $                         x( j + 1 )*a( j + 1, j + 2 )
                  if( nounit )
     $               x( j + 2 ) = x( j + 2 )/a( j + 2, j + 2 )
                  x( j + 3 ) = x( j + 3 ) - x( j )*a( j, j + 3 ) -
     $                         x( j + 1 )*a( j + 1, j + 3 ) - x( j + 2 )
     $                         *a( j + 2, j + 3 )
                  if( nounit )
     $               x( j + 3 ) = x( j + 3 )/a( j + 3, j + 3 )
c**** unrolled main loop ***********************************************
                  temp1 = x( j )
                  temp2 = x( j + 1 )
                  temp3 = x( j + 2 )
                  temp4 = x( j + 3 )
                  if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.
     $                zero.or.temp4.ne.zero )then
                     do 90, i = j + 4, n
                        x( i ) = ( ( ( ( x( i ) - temp1*a( j, i ) )
     $                           - temp2*a( j + 1, i ) )
     $                           - temp3*a( j + 2, i ) )
     $                           - temp4*a( j + 3, i ) )
   90                continue
                  end if
  100          continue
c**** left-overs on top left corner ************************************
               if( n4 + 1.le.n )then
                  if( nounit )
     $               x( n4 + 1 ) = x( n4 + 1 )/a( n4 + 1, n4 + 1 )
               end if
               if( n4 + 2.le.n )then
                  x( n4 + 2 ) = x( n4 + 2 ) - x( n4 + 1 )
     $                          *a( n4 + 1, n4 + 2 )
                  if( nounit )
     $               x( n4 + 2 ) = x( n4 + 2 )/a( n4 + 2, n4 + 2 )
               end if
               if( n4 + 3.le.n )then
                  x( n4 + 3 ) = x( n4 + 3 ) - x( n4 + 1 )
     $                          *a( n4 + 1, n4 + 3 ) - x( n4 + 2 )
     $                          *a( n4 + 2, n4 + 3 )
                  if( nounit )
     $               x( n4 + 3 ) = x( n4 + 3 )/a( n4 + 3, n4 + 3 )
               end if
            else
               jx = kx
c**** u n r o l l   t o   d e p t h  4 ********** l o o p  d ***********
               n4 = 4*( n/4 )
               do 120, j = 1, n4, 4
c**** 4 by 4 triangle at bottom of block *******************************
                  if( nounit )
     $               x( jx ) = x( jx )/a( j, j )
                  x( jx + incx ) = x( jx + incx ) - x( jx )
     $                             *a( j, j + 1 )
                  if( nounit )
     $               x( jx + incx ) = x( jx + incx )/a( j + 1, j + 1 )
                  x( jx + 2*incx ) = x( jx + 2*incx ) - x( jx )
     $                               *a( j, j + 2 ) - x( jx + incx )
     $                               *a( j + 1, j + 2 )
                  if( nounit )
     $               x( jx + 2*incx ) = x( jx + 2*incx )
     $                                  /a( j + 2, j + 2 )
                  x( jx + 3*incx ) = x( jx + 3*incx ) - x( jx )
     $                               *a( j, j + 3 ) - x( jx + incx )
     $                               *a( j + 1, j + 3 ) -
     $                               x( jx + 2*incx )*a( j + 2, j + 3 )
                  if( nounit )
     $               x( jx + 3*incx ) = x( jx + 3*incx )
     $                                  /a( j + 3, j + 3 )
c**** unrolled main loop ***********************************************
                  temp1 = x( jx )
                  temp2 = x( jx + incx )
                  temp3 = x( jx + 2*incx )
                  temp4 = x( jx + 3*incx )
                  if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.
     $                zero.or.temp4.ne.zero )then
                     ix = jx + 3*incx
                     do 110, i = j + 4, n
                        ix = ix + incx
                        x( ix ) = ( ( ( ( x( ix ) - temp1*a( j, i ) )
     $                            - temp2*a( j + 1, i ) )
     $                            - temp3*a( j + 2, i ) )
     $                            - temp4*a( j + 3, i ) )
  110                continue
                  end if
                  jx = jx + 4*incx
  120          continue
c**** left-overs on top left corner ************************************
               if( n4 + 1.le.n )then
                  if( nounit )
     $               x( jx ) = x( jx )/a( n4 + 1, n4 + 1 )
               end if
               if( n4 + 2.le.n )then
                  x( jx + incx ) = x( jx + incx ) - x( jx )
     $                             *a( n4 + 1, n4 + 2 )
                  if( nounit )
     $               x( jx + incx ) = x( jx + incx )/a( n4 + 2, n4 + 2 )
               end if
               if( n4 + 3.le.n )then
                  x( jx + 2*incx ) = x( jx + 2*incx ) - x( jx )
     $                               *a( n4 + 1, n4 + 3 ) -
     $                               x( jx + incx )*a( n4 + 2, n4 + 3 )
                  if( nounit )
     $               x( jx + 2*incx ) = x( jx + 2*incx )
     $                                  /a( n4 + 3, n4 + 3 )
               end if
            end if
         else
            if( incx.eq.1 )then
c**** u n r o l l   t o   d e p t h  4 ********** l o o p  g ***********
               n4 = mod( n, 4 ) + 1
               do 140, j = n, n4, -4
c**** 4 by 4 triangle at bottom of block *******************************
                  if( nounit )
     $               x( j ) = x( j )/a( j, j )
                  x( j - 1 ) = x( j - 1 ) - x( j )*a( j, j - 1 )
                  if( nounit )
     $               x( j - 1 ) = x( j - 1 )/a( j - 1, j - 1 )
                  x( j - 2 ) = x( j - 2 ) - x( j )*a( j, j - 2 ) -
     $                         x( j - 1 )*a( j - 1, j - 2 )
                  if( nounit )
     $               x( j - 2 ) = x( j - 2 )/a( j - 2, j - 2 )
                  x( j - 3 ) = x( j - 3 ) - x( j )*a( j, j - 3 ) -
     $                         x( j - 1 )*a( j - 1, j - 3 ) - x( j - 2 )
     $                         *a( j - 2, j - 3 )
                  if( nounit )
     $               x( j - 3 ) = x( j - 3 )/a( j - 3, j - 3 )
c**** unrolled main loop ***********************************************
                  temp1 = x( j )
                  temp2 = x( j - 1 )
                  temp3 = x( j - 2 )
                  temp4 = x( j - 3 )
                  if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.
     $                zero.or.temp4.ne.zero )then
                     do 130, i = j - 4, 1, -1
                        x( i ) = ( ( ( ( x( i ) - temp1*a( j, i ) )
     $                           - temp2*a( j - 1, i ) )
     $                           - temp3*a( j - 2, i ) )
     $                           - temp4*a( j - 3, i ) )
  130                continue
                  end if
  140          continue
c**** left-overs on top left corner ************************************
               if( n4.ge.2 )then
                  if( nounit )
     $               x( n4 - 1 ) = x( n4 - 1 )/a( n4 - 1, n4 - 1 )
               end if
               if( n4.ge.3 )then
                  x( n4 - 2 ) = x( n4 - 2 ) - x( n4 - 1 )
     $                          *a( n4 - 1, n4 - 2 )
                  if( nounit )
     $               x( n4 - 2 ) = x( n4 - 2 )/a( n4 - 2, n4 - 2 )
               end if
               if( n4.ge.4 )then
                  x( n4 - 3 ) = x( n4 - 3 ) - x( n4 - 1 )
     $                          *a( n4 - 1, n4 - 3 ) - x( n4 - 2 )
     $                          *a( n4 - 2, n4 - 3 )
                  if( nounit )
     $               x( n4 - 3 ) = x( n4 - 3 )/a( n4 - 3, n4 - 3 )
               end if
            else
               jx = kx + ( n - 1 )*incx
c**** u n r o l l   t o   d e p t h  4 ********** l o o p  b ***********
               n4 = mod( n, 4 ) + 1
               do 160, j = n, n4, -4
c**** 4 by 4 triangle at bottom of block *******************************
                  if( nounit )
     $               x( jx ) = x( jx )/a( j, j )
                  x( jx - incx ) = x( jx - incx ) - x( jx )
     $                             *a( j, j - 1 )
                  if( nounit )
     $               x( jx - incx ) = x( jx - incx )/a( j - 1, j - 1 )
                  x( jx - 2*incx ) = x( jx - 2*incx ) - x( jx )
     $                               *a( j, j - 2 ) - x( jx - incx )
     $                               *a( j - 1, j - 2 )
                  if( nounit )
     $               x( jx - 2*incx ) = x( jx - 2*incx )
     $                                  /a( j - 2, j - 2 )
                  x( jx - 3*incx ) = x( jx - 3*incx ) - x( jx )
     $                               *a( j, j - 3 ) - x( jx - incx )
     $                               *a( j - 1, j - 3 ) -
     $                               x( jx - 2*incx )*a( j - 2, j - 3 )
                  if( nounit )
     $               x( jx - 3*incx ) = x( jx - 3*incx )
     $                                  /a( j - 3, j - 3 )
c**** unrolled main loop ***********************************************
                  temp1 = x( jx )
                  temp2 = x( jx - incx )
                  temp3 = x( jx - 2*incx )
                  temp4 = x( jx - 3*incx )
                  if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.
     $                zero.or.temp4.ne.zero )then
                     ix = jx - 3*incx
                     do 150, i = j - 4, 1, -1
                        ix = ix - incx
                        x( ix ) = ( ( ( ( x( ix ) - temp1*a( j, i ) )
     $                            - temp2*a( j - 1, i ) )
     $                            - temp3*a( j - 2, i ) )
     $                            - temp4*a( j - 3, i ) )
  150                continue
                  end if
                  jx = jx - 4*incx
  160          continue
c**** left-overs on top left corner ************************************
               if( n4.ge.2 )then
                  if( nounit )
     $               x( jx ) = x( jx )/a( n4 - 1, n4 - 1 )
               end if
               if( n4.ge.3 )then
                  x( jx - incx ) = x( jx - incx ) - x( jx )
     $                             *a( n4 - 1, n4 - 2 )
                  if( nounit )
     $               x( jx - incx ) = x( jx - incx )/a( n4 - 2, n4 - 2 )
               end if
               if( n4.ge.4 )then
                  x( jx - 2*incx ) = x( jx - 2*incx ) - x( jx )
     $                               *a( n4 - 1, n4 - 3 ) -
     $                               x( jx - incx )*a( n4 - 2, n4 - 3 )
                  if( nounit )
     $               x( jx - 2*incx ) = x( jx - 2*incx )
     $                                  /a( n4 - 3, n4 - 3 )
               end if
            end if
         end if
      end if
c
      return
c
c     end of f06pjf (dtrsv ).
c
      end


! From transition.f90
      subroutine transition
      use xcdcc
      use sistema
      use wfs, only:energ,idx,rvec
      use ptpots
      use globals
      use channels
      use memory
      use potentials, only: vcl
      implicit real*8 (a-h,o-z)
      logical skip,writeff
      integer qcmin,qcmax,kcmax,coups!,lambmax
      integer lgs
c     -------------------------------------------------------------------------------
      character*5:: jpi
      character*3:: jname
      CHARACTER PARITY(3)
      DATA PARITY / '-','?','+' / 
      character*40 filename,couplings,comment,fileampsd
c     -------------------------------------------------------------------------------
      real*8 factor,ptr,ttr,jt,iin,ifi,kband1,kband2
      real*8,allocatable:: xrad2(:),xrad3(:)
      real*8 :: qfactorc(0:50),qfactorn(0:50)
      
c changed to complex in v2.3
      complex*16,pointer:: ui(:),uf(:)
      complex*16 :: fprod
      complex*16 xi,pot,xsum, xsumn,xsumc
!      parameter(Kmax=6,xKrot=0.d0,nkmax=300)
      parameter(xKrot=0.d0,nkmax=300)
      
c  For adiabatic JT potential
      complex*16,allocatable:: xintad(:),vadnorm(:)
      complex*16 :: fauxad
      complex*16 :: xsumad
      real*8     ::xsumnad
      real*8,allocatable:: vclquad(:)
      
c     ------------------------------------------------------------------------------

! v0.6c
!      dimension energ(nener,numfilmax),!,xrad2(:),elapsed(2),
!     .np(numfilmax)
c      real*8, allocatable:: energ(:,:)
c      integer, allocatable:: np(:)
      TYPE GFAC  ! derived variable to store the geometrical factors
       INTEGER i,f,k,l,nq
       REAL*8:: rmatc,rmatn
      END TYPE
      TYPE(GFAC) PK(nkmax)
      interface
       subroutine potdef(vc,qmax,rnc,vdef,betar)
       integer::qmax
       real*8:: vc(:),rnc,vdef(:),betar
      end subroutine
      end interface

      interface
        subroutine extpot(filename,vr,vi,nr)
        integer::nr
        character*40 filename
        real*8 :: vr(:),vi(:)
       end subroutine
      end interface

! Modified by AMM in v2.2
!      namelist /trans/ mp,mt,zt,zp,rcc,namep,namet
      namelist /trans/ skip,rcc,writeff
      namelist /grid/ nquad,radmax,rmax,rextrap,rstep,rvecin,drvec,hin
      namelist /wf/ filename,eminf,emaxf
      namelist /coupling/ qcmin,qcmax,kcmax,lamin,lamax,coups,ncoul,
     &          qfactorn,qfactorc
      namelist /BEQ/ bel
!      interface
!      subroutine gauleg(x1,x2,x,w,n)
!      integer::n
!      real*8:: x1,x2,x(:),w(:)
!      end subroutine
!      end interface
      call cpu_time(start)
      debug=.false.
      xi=(0.d0,1.d0)
      kin=5
      kread=84
      kfr=4
      ken=8
      pi=acos(-1d0)
      writeff=.false.
      qfactorc(:)=1.0; qfactorn(:)=1.0; 
c     ----------------------------------------------------------------------------      
      write(*,'(/,5x,"*********** COUPLING POTENTIALS  *********** ")')
      
c ---------------------------------------------------------------------------------
c AMoro: partial preread to set dimensions
      rextrap=0
      skip=.false. 
      read(kin,nml=trans)

      if(jpsets.eq.0) then
      write(*,'(/,5x,"** No basis states have been defined",
     &    "=> skipping this section **" )') 
      return
      endif

      if (skip) then
      write(*,'(5x,"[No coupling potentials requested ]")')
      return
      endif

c *** Quantum number limits for transition potentials
      lamin=0
      lamax=-1
      qcmin=0
      coups=0
      qcmax=-1
      iftrans=.true.
      kcmax=-1
      read(kin,nml=coupling)
      if (qcmax.lt.0) then
       write(*,*)'[=> No coupling potentials requested]'
       iftrans=.false.
       return
      endif
!      write(*,'(/,2x,"o Couplings considered for formfactors:")')
      select case(coups)
      case (0)
       write(*,*) ' [ coups=0=> All transitions will be considered]'
      case (1)
       write(*,*) ' [ coups=1=> ONLY gs->cont couplings considered]'
      case (2)
       write(*,*) ' [ coups=2=> Diag nucl. + ALL Coulomb]'
       if (ncoul.ne.0) write(*,*) ' ( ncoul will be ignored )'
      case (3)
       write(*,*) ' [ coups=3=> V00 + gs->cont]'
   
      case (4)
       write(*,*) ' [ coups=4=> gs->cont (no diag) ]'
   
      end select
      
      select case(ncoul)
      case (0)
        write(*,*)' [ ncoul=0 => C+N couplings ]'
      case (1)
        write(*,*)' [ ncoul=1 => NUCLEAR couplings only ]'
      case (2)
        write(*,*)' [ ncoul=2 => COULOMB couplings only]'
      end select

!	  if (any(qfactorc(1:50)).ne.1.0)
!      write(*,*)"qfactorc=",qfactorc(0:5)
!      write(*,*)"qfactorc=",qfactorn(0:5)

c *** Radial grids ------------------------------------------------------------
      read(kin,nml=grid)
      coefc=av/mp
      coefv=-ac/mp ! = coefc-1
      coef=max(coefc,dabs(coefv))
      if (rstep.gt.1e-5) then 
      nrad2=nint((rmax-rstep)/rstep)+1
      if (rextrap.lt.rmax) rextrap=rmax
      nrad3=nint((rextrap-rstep)/rstep)+1
      allocate(xrad2(nrad2))
      allocate(xrad3(nrad3))
      else
      write(*,*)'radial step (rstep) too small!'; stop
      endif
      allocate(xquad(nquad),wquad(nquad),rquad(nquad))
      allocate(xintgn(nquad),xintgc(nquad)) ! AMoro
      allocate(xintad(nquad),vclquad(nquad),vadnorm(nquad)) ! adiabatic
      a=-1.d0
      b=1.d0
      call gauleg(a,b,xquad,wquad,nquad)

      rquad=radmax*0.5d0*(1.d0+xquad)
      rvecmax=coef*radmax+rmax
      nr=nint((rvecmax-rvecin)/drvec)
      nrin=nint(rvecmax/hin)
      allocate(rfrag(nr))
      allocate(rvin(nrin))
      do irvec=1,nr
      rfrag(irvec)=rvecin+drvec*irvec
      enddo
      do irvec=1,nrin
      rvin(irvec)=hin*irvec
      enddo

      if(.not.allocated(vcore))  allocate(vcore(1:nr,0:nmult))
      if(.not.allocated(vcorei)) allocate(vcorei(1:nr,0:nmult))
      if(.not.allocated(vcorec)) allocate(vcorec(1:nr,0:nmult))
      if(.not.allocated(vval))   allocate(vval(1:nr,0:nmult))
      if(.not.allocated(vvali))  allocate(vvali(1:nr,0:nmult))
      if(.not.allocated(vvalc))  allocate(vvalc(1:nr,0:nmult))
      kptype=0
      vfrag => vcore
      vfragi=> vcorei
      vcoup => vcorec
      write(*,'(/,2x,"** CORE-TARGET potential:")')
      call read_fragpot(kin,kptype,zc)
      melc(1:nmult)=mel(1:nmult)

      nullify(vcoup)
      kptype=1
      vfrag =>vval
      vfragi=>vvali
      vcoup =>vvalc
      write(*,'(/,2x,"** VALENCE-TARGET potential:")')
      call read_fragpot(kin,kptype,zv)
      melv(1:nmult)=mel(1:nmult)
      do irad2=1,nrad2
      xrad2(irad2)=rmax/dble(nrad2)*irad2
      enddo

c *** Interpolate projectile wfs at quadrature points
      call wf2quad()
      
c *** Interpolate projectile v+c potential at quadrature points
! Assuming single-channel and only central potential!!!!!!!!!!!!!!!!!!!
      alpha=0.0
      do iq=1,nquad
      raux=rquad(iq)
      lgs=jpiset(1)%lsp(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GENERALIZE FOR MULTICHANNEL!!!!!!!!!!!!
!      write(0,*) lgs,nr,vcl(lgs,1), size(vcl(lgs,:))
      faux=fival(raux,rvec,vcl(lgs,:),size(rvec),alpha) ! R(r)
!      faux=FFC4((raux-rmin)/dr,yvp,nr)
      vclquad(iq)=faux 
      write(95,'(1f8.3,2x,100g14.6)')raux,faux !,frad(jset,ie,iq,ich)
      enddo ! iq=quadrature point      


c *** Set limits for quantum numbers & multipoles ----------------------------------
      npe    =maxval(jpiset(:)%nex)
      xImax  =maxval(spchan(:,:)%jc)
      xJpmax =maxval(jpiset(:)%jtot)
      nIcmax =nint(2.d0*xImax)
!      nQmax=min(nmult,nIcmax) ! RDD

!!!! AMM: Iin, Iff, need to be checked
      nQmin=max(0,qcmin)
      nQmin=qcmin
      nQmax=min(qcmax,nIcmax)
      if (lamax.ge.0) then
          lambmax=min(lamax,nint(2.d0*xJpmax))
      else
          lambmax=nint(2.d0*xJpmax)
      endif
      
!! New in Mar/17
      lambmin=max(nint(dabs(xjp1-xjp2)),lamin)


!!! NEW IN v2.2g (CHECK!!)
      kmax=nqmax+lambmax
      if (kcmax.ge.0) kmax=min(kmax,kcmax)  

!      write(*,*)'qcmin,nqmin,Iin,Ifi=',qcmin,nqmin,Iin,Ifi
!      write(*,400) nqmin,nqmax,0,min(kmax,kcmax),0,lambmax
      write(*,400) nqmin,nqmax,0,kmax,lambmin,lambmax

400   format(/,2x,"o Quantum numbers limits:",/,5x,
     &  "- Core multipoles:",i2,"-",i1,/,5x,
     &  "- K:",i2,"-",i2,/,5x,
     &  "- Proj.-target multipoles:", i2,"-" i1)
c -------------------------------------------------------------------------------------


c *** START CALCULATION OF COUPLING POTENTIALS ---------------------------------------
      write(*,'(/,2x,"** FORMFACTORS **")') 
      fmem=nex*nex*(lambmax+1)*nrad3*lc16/1e6
      write(*,'(5x," [ formfactors need",1f7.1," Mbytes ]")') fmem

      allocate(potQKn(nquad,nrad2,0:nQmax,0:Kmax))
      allocate(potQKc(nquad,nrad2,0:nQmax,0:Kmax))
      allocate(Ff(nex,nex,0:lambmax,nrad3))
c commented by AMoro, to save memory
!      allocate(Fc(nex,nex,0:lambmax,nrad3))
!      allocate(Fn(nex,nex,0:lambmax,nrad3)) 
!      Fc=0.d0;  Fn=0.d0
      Ff=0.d0
      potQKn=0d0; potQKc=0d0
      rmat=0d0



       
      write(*,'(/,2x,"o V^{K,Q}(r,R) ",$)')
      call cpu_time(t1)
      call factorial(2*nQmax)

      do iquad=1,nquad
      do irad2=1,nrad2
      do nq=nqmin,nqmax
      do k=0,kmax
! AMoro: 
! ncoul=0: coulomb + nuclear
! ncoul=1: nuclear only
! ncoul=2: Coulomb only  
      potn=0d0; potc=0d0
      r1=rquad(iquad); ! internal coordinate (r)
      r2=xrad2(irad2)  ! p-t coordinate (R)
      if ((ncoul.eq.0).or.(ncoul.eq.1)) then  !nuclear part
        nc=1 ! nuclear part
        if(nq.eq.0) then                      !monopole
        potn= pot(r1,r2,nq,k,coefc,nc)
     &      + pot(r1,r2,nq,k,coefv,nc)
        else !      no core excitation for valence fragment! 
        potn=pot(r1,r2,nq,k,coefc,nc)
        endif
      endif
      if ((ncoul.eq.0).or.(ncoul.eq.2)) then ! Coulomb part
       nc=2 ! Coulomb part
       potc= pot(r1,r2,nq,k,coefc,nc)
     &     + pot(r1,r2,nq,k,coefv,nc)
      endif       

      potQKc(iquad,irad2,nq,k)=potc  
      potQKn(iquad,irad2,nq,k)=potn
      
      enddo !k 
      enddo !nq
      
      enddo ! nrad2 (R)
      enddo ! nquad (r)
      call cpu_time(t2)
      write(*,'(4x,"[",1f6.2," secs ]")') t2-t1
     
      rfirst=rstep 
      write(*,300) rfirst,rmax,rstep
300   format(/,5x,"[Radial grid: Rmin=",1f6.3," fm,", 
     &      " Rmax=",1f6.1," fm," 
     &      " Step=",1f6.3," fm]",/)

      if (rextrap.gt.rmax) then
       write(*,301) rextrap
301   format(5x,"[ Coulomb couplings extrapolated to",1f8.1," fm]",/)
      endif
      nff=0
      do m=1,jpsets
        nchann1=jpiset(m)%nchan ! nch(m)
        xjp1=jpiset(m)%jtot! xjp(m)
! TEst 11/OCT/16
!      do n=m,jpsets
      do n=1,jpsets
      if(realwf.and.n.lt.m) cycle
!!!!!!!!!!!!!!!!!!
        nchann2=jpiset(n)%nchan ! nch(n)
        xjp2=jpiset(n)%jtot     ! xjp(n)
      
!      do lc=nint(dabs(xjp1-xjp2)),min(lambmax,nint(xjp1+xjp2))
      do lc=lambmin,min(lambmax,nint(xjp1+xjp2))
        xlc=dble(lc)
        nk=0  !Nb. of P(IK) for each n,m,Lambda
        call cpu_time(t1)
        do i=1,nchann1
       xl1=jpiset(m)%lsp(i) 
       ic1=jpiset(m)%cindex(i)
       kband1=qnc(ic1)%kband
       Iin=jpiset(m)%jc(i) 

       do j=1,nchann2      
       xl2=jpiset(n)%lsp(j) 
       Ifi=jpiset(n)%jc(j) 
       ic2=jpiset(n)%cindex(j)
       kband2=qnc(ic2)%kband

! AMM: Assume rotor model and kband1=kband2
       if (kband1.ne.kband2) cycle
     
      do k=0,kmax
      do nq=max(nint(Iin-Ifi),nqmin),nqmax
      do l=0,nq ! small lambda
! CHECK AMM : i <-> f
!      call rmatel(m,n,xlc,i,j,k,nq,l,xjp1,xjp2,xKrot,rmatn,rmatc)
      call rmatel(n,m,xlc,j,i,k,nq,l,xjp2,xjp1,kband1,rmatn,rmatc)

      if (abs(rmatn).lt.1e-6.and.abs(rmatc).lt.1e-6) cycle
      nk=nk+1
      
      if (nk.gt.nkmax) then
        write(*,*)'NK > NKMAX=',NKMAX,' SO INCREASE NKMAX'
        stop
      endif
      pk(nk)%i     =i
      pk(nk)%f     =j
      pk(nk)%k     =k
      pk(nk)%nq    =nq  
      pk(nk)%l     =l    ! \lambda
      pk(nk)%rmatc =rmatc
      pk(nk)%rmatn =rmatn
      if ((m.eq.1).and.(n.eq.1).and.(verb.gt.0)) then
        if (nk.eq.1) then
        write(198,'("j=",1f4.2, " Jf=",1f4.2)')xjp1,xjp2
        write(198,'(7A7,2A10)')
     & "NK","Ch1","Ch2","LAM","K","QC","lambda","P(N)", "P(C)"
        endif
      endif
      if (verb.gt.0) write(198,1100) nk,i,j,lc,k,nq,l,
     & rmatn,rmatc,xl1,xl2
1100  format(7i7,2x,2f10.3,2x,2f4.1,2i5)
      enddo ! l (lambda)
      enddo ! nq (core multipole Q)
      enddo ! K
      enddo ! j (nchan)
      enddo ! i (nchan)


      if (nk.gt.0) then
      write(*,'(4x,a6,"-> ",a6," with LAM=",1i2,
     &   2x,"NK=",i3,2x, "non-zero P(K)  ",$)') 
     &   jpi(xjp1,jpiset(m)%partot),
     &   jpi(xjp2,jpiset(n)%partot),lc,nk
      else        
        if (debug) write(*,'(4x,"m=",i2,"-> n=",i2," with LAM=",1i2,
     &  2x,"(NO allowed transitions)")') m,n,lc
      endif


c ---------------------------------------------------------
c    Calculate formfactors F(R) for each |n,ie1> ->|m,ie2> 
c ---------------------------------------------------------
      if (nk.eq.0) cycle     
      if (debug) write(*,*)'radial integrals'
      do ie1=1,jpiset(m)%nex ! np(m)
      id1=idx(m,ie1)
      do ie2=1,jpiset(n)%nex ! np(n)
      id2=idx(n,ie2)
      nff=nff+1 !number of FF
      do irad2=1,nrad2
      do ik=1,NK
        i=pk(ik)%i 
		j=pk(ik)%f
        k=pk(ik)%k
		nq=pk(ik)%nq
		l=pk(ik)%l
        rmatc=pk(ik)%rmatc
		rmatn=pk(ik)%rmatn
!  radial integral: R(ik)=f*V(K,Q)*f'
        xsumn=(0.d0,0.d0)
        xsumc=(0.d0,0.d0)
        xsumad=0.
        xsumnad=0.
        if (abs(rmatc).lt.1e-5.and.abs(rmatn).lt.1e-5) goto 1200
      
!        ui=>frad(m,ie1,:,i)
!        uf=>frad(n,ie2,:,j)
!      call radquad(ui,uf,potQKn(:,irad2,nq,k))
      
      do iquad=1,nquad
      vmon=0d0
! AMM: frad(j) should be conjugate!! (fixed Sept 16)  
!      fprod=frad(m,ie1,iquad,i)*frad(n,ie2,iquad,j) ! i -> j = <j | V | i> 
      fprod=conjg(frad(n,ie2,iquad,j))*frad(m,ie1,iquad,i) ! i -> j = <j | V | i>
!      if (irad2.eq.1) 
!     & write(94,*) rquad(iquad),real(frad(m,ie1,iquad,i)), real(fprod)  
      if(l.eq.nq) then
      xintgn(iquad)=fprod*potQKn(iquad,irad2,nq,k)*xrad2(irad2)**l


      if ((m.eq.n).and.(ie1.eq.ie2).and.(n.eq.1).and.(k.eq.0)) then 
! AMM: For adiabatic JT potential
      xintad(iquad)=fprod*vclquad(iquad)*((
     &              potQKn(iquad,irad2,nq,k)+potQKc(iquad,irad2,nq,k))
     &              *xrad2(irad2)**l-VCOUL(xrad2(irad2),zp,zt,Rcc))
      vadnorm(iquad)=fprod*vclquad(iquad)
      endif

c subtract projectile-target monopole Coulomb
      if ((nq.eq.0).and.zp*zt.gt.1e-3.and.k.eq.0.and.ncoul.ne.1)
     & then
!     &  .and.(m.eq.n).and.(ie1.eq.ie2)) then 
       vmon=VCOUL(xrad2(irad2),zp,zt,Rcc)     
      endif

      xintgc(iquad)=fprod*
     & (potQKc(iquad,irad2,nq,k)*xrad2(irad2)**l-vmon)
      
      else 
!      write(0,*)'l,nq=',l,nq
      xintgn(iquad)=fprod*
     .potQKn(iquad,irad2,nq,k)*xrad2(irad2)**l*
     .(coefc*rquad(iquad))**(nq-l)

      xintgc(iquad)=fprod*
     .potQKc(iquad,irad2,nq,k)*xrad2(irad2)**l*
     .(coefc*rquad(iquad))**(nq-l)
      endif
      xsumn=xsumn + xintgn(iquad)*wquad(iquad)
      xsumc=xsumc + xintgc(iquad)*wquad(iquad)
      
! ADIABATIC      
      if ((m.eq.n).and.(ie1.eq.ie2).and.(n.eq.1).and.(k.eq.0)) then 
      xsumad=xsumad   +  xintad(iquad)*wquad(iquad)
      xsumnad=xsumnad + vadnorm(iquad)*wquad(iquad)
      !if (irad2.eq.1) write(0,*)iquad,vadnorm(iquad)
!      if (irad2.eq.1) write(0,*)iquad,xsumnad*0.5d0*radmax
      endif 
      
      enddo ! iquad
     

1200  fauxn=rmatn*xsumn*dsqrt(2.d0*dble(k)+1.d0)*0.5d0*radmax
      fauxc=rmatc*xsumc*dsqrt(2.d0*dble(k)+1.d0)*0.5d0*radmax
      
      
! ADIABATIC      
      if ((m.eq.n).and.(ie1.eq.ie2).and.(n.eq.1).and.(k.eq.0)) then 
!      write(0,*)'xsumna=',xsumnad
      xsumad=xsumad/xsumnad ! divide by norm <phi|Vbx| phi> 
      fauxad=rmatn*xsumad*dsqrt(2.d0*dble(k)+1.d0)!*0.5d0*radmax
      write(999,'(1f10.3,2g16.6)') xrad2(irad2),fauxad
      endif
      
! April 2019: only diagonal nuclear couplings for coups=2
      if (coups.eq.2.and.(ie1.ne.ie2.or.m.ne.n)) fauxn=0 

! March 2021      
      if (coups.eq.3.and.(ie1.ne.ie2.or.m.ne.n)) fauxn=0 
      
! Scale formfactors by qfactor()
!       write(*,*) lc,qfactorc(lc), qfactorn(lc),fauxc,fauxn
	   if (lc.lt.50) fauxc=fauxc*qfactorc(lc)
	   if (lc.lt.50) fauxn=fauxn*qfactorn(lc)
!       write(*,*) lc,qfactorc(lc), qfactorn(lc),fauxc,fauxn
 
!      Fc(id1,id2,lc,irad2)=Fc(id1,id2,lc,irad2)+fauxc    ! coulomb
!      Fn(id1,id2,lc,irad2)=Fn(id1,id2,lc,irad2)+fauxn    ! nuclear
      Ff(id1,id2,lc,irad2)=Ff(id1,id2,lc,irad2)+fauxn+fauxc ! total
      enddo ! ik
      enddo ! irad2 (R)
      enddo ! ie2
      enddo ! ie1
      call cpu_time(t2)
      write(*,'("(",1f6.2," secs )")')t2-t1
      call flush(6)
      enddo ! lc (\Lamda)
      enddo ! n (jpsets)
      enddo ! m (jpsets)

      call cpu_time(finish)
!      print*, 'Time = ',finish-start,'seconds'
      write(*,1220) nff,finish-start
1220  format(/,5x,"-> [",i7, " formfactors generated in ",
     & 1x,1f6.1," secs ]",/)



c    -----------------------------------------------------------
c    Write formfactors
c    ----------------------------------------------------------- 
      if (writeff) open(kfr,file="ff.fr",status='unknown')
      ttr=0.d0
      fscale=1.d0
      npa1=0
      write(*,*)' Output files:'
      select case(ncoul)
      case(0) ! nuc + coul
      write(*,*) ' ff.fr:  COUL+NUC couplings for Fresco'
!      ffr=>Ff(:,:,:,:,:,:)
      case(1) ! nuc
      write(*,*) ' ff.fr: NUCLEAR couplings for Fresco'
!      ffr=>Fn(:,:,:,:,:,:)
      case(2) ! coul
      write(*,*) ' ff.fr: COULOMB couplings for Fresco'
!      ffr=>Fc(:,:,:,:,:,:)
      end select 

      do if1=1,jpsets
      npa1=npa1+jpiset(if1)%nex ! np(if1)
      npa2=0
      xjp1 =jpiset(if1)%jtot
      ipar1=jpiset(if1)%partot
      do if2=1,jpsets
      npa2  =npa2+jpiset(if2)%nex ! np(if2)
      xjp2  =jpiset(if2)%jtot
      ipar2 =jpiset(if2)%partot

      do ie1=1,jpiset(if1)%nex  ! np(if1)
      do ie2=1,jpiset(if2)%nex  ! np(if2)
      e1=energ(if1,ie1)
      e2=energ(if2,ie2)

      m1=ie1+npa1-jpiset(if1)%nex !np(if1)
      m2=ie2+npa2-jpiset(if2)%nex !np(if2)

      select case(coups)
      case (0) ! all couplings (do nothing)
      case (1,3) ! gs->cont + diagonal
        if ((if1.ne.1).or.(ie1.ne.1)) then
          if ((if2.ne.if1).or.(ie2.ne.ie1)) cycle
        endif
        
      case(4) ! only gs-> cont NO DIAG
         if ((if1.ne.1).or.(ie1.ne.1)) cycle
         if ((if1.eq.if2).and.(ie2.eq.ie1)) cycle
     
        
!      case (3) ! V00 + gs->cont
!         if ((if1.ne.1).or.(ie1.ne.1)) cycle
!         if ((if2.ne.if1).or.(ie2.eq.ie1)) cycle
      end select

      if ((m1.gt.m2).and.realwf) cycle !!!! AMORO
c for fresco
      do lc= nint(dabs(xjp1-xjp2)),min(nint(xjp1+xjp2),lambmax)
      if (ipar1*ipar2*(-1)**lc<0) cycle
      ltr=lc
      ptr=lc
      
c \hat{Jp}*hat{Jp'}*(2*Lambda+1)*(-1)^Lambda
      factor=(2d0*ptr+1)*sqrt(2*xjp1+1)*sqrt(2*xjp2+1)*(-1)**ptr

      if (debug) then
      write(*,'(a,a,1f4.1,a,1f4.1,a,i3,2x,a,1g16.8)')
     & ' Conversion factor in F(r)[fres]=factor*F(r)',
     &  'ji=',xjp1,'-> jf=',xjp2,' LAM=',lc," =>  factor=",factor
      endif

      write(comment,'(2x,"<",i3,"|",i2,"|",i3,">")') m1,lc,m2
      if (writeff) 
     & write(kfr,500) nrad3,rstep,rfirst,fscale,ltr,ptr,ttr,m2,m1,
     &               comment   ! ORIG

!!! TEST APRIL
!      if (writeff) 
!     & write(888,500) nrad3,rstep,rfirst,fscale,ltr,ptr,ttr,m2,m1,
!     &               comment   ! ORIG


      if (verb.ge.4) write(120,'("# <",i3,"|",i2,"|",i3,">")') m1,lc,m2
      do irad=1,nrad2
      r2=xrad2(irad)
      fauxc=Ff(m1,m2,lc,irad)
      
      if ((coups.eq.3).and.(m1.eq.m2)) fauxc=Ff(1,1,lc,irad)
         
!      fauxn=0d0
!      if ((ncoul.eq.0).or.(ncoul.eq.1)) then
!       fauxn=Fn(m1,m2,lc,irad)
!      endif      
!      if ((ncoul.eq.0).or.(ncoul.eq.2)) then
!       fauxc=Fc(m1,m2,lc,irad)   
!      endif          
!      write(kfr,'(2x,1g16.10,2x,1g16.10)') factor*(fauxc+fauxn)
!      write(120,'(1x,1f8.3,2x,2g16.8)') r2,factor*(fauxc+fauxn)
       if (writeff) 
     &      write(kfr,'(2x,1g16.10,2x,1g16.10)') factor*fauxc
      if (verb.ge.4) 
     & write(120,'(1x,1f8.3,2x,2g16.8)') r2,factor*fauxc
      enddo !irad

c Check that imaginary part is small
      if (abs(aimag(fauxc)).gt.1e-5) then
      write(*,480) m1,xjp1,parity(ipar1+2),m2,xjp2,parity(ipar2+2),lc,
     & r2,fauxc
      endif
480   format('** Warning **: Large imaginary part for coupling: ',
     & 1i3,' (',1f4.1,a1,") -> ", 1i3," (",1f4.1,a1,")",
     & " for LAM=",1i3, ": R=",1f8.2," Fc(R)=",2g16.8)

c Extrapolate Coulomb formactors from R=Rmax to Rextrap
      if (rextrap.gt.rmax) then
      caux=real(fauxc)*rmax**(lc+1)
      
      do ir=nrad2+1,nrad3
      r2=rstep+rstep*(ir-1)
      fauxc=0
      if ((lc.gt.0).and.(ncoul.ne.1)) then
      fauxc=caux/r2**(lc+1)
      endif
     
      if (writeff) write(kfr,'(2x,1g16.10,2x,1g16.10)') factor*fauxc
      if (verb.ge.4) write(120,'(1x,1f8.3,2x,2g16.8)') r2,factor*fauxc
!      Fc(m1,m2,lc,ir)=fauxc 
      Ff(m1,m2,lc,ir)=fauxc ! total  
      enddo !nrad3
      endif ! rextrap> rmax
      
!!! TEST APRIL 18 
!      do ir=1,nrad3
!      r2=rstep+rstep*(ir-1)
!      fauxc=0
!      if ((lc.gt.0).and.(ncoul.ne.1)) then
!      fauxc=caux/r2**(lc+1)
!      endif
!      if (writeff) write(888,'(2x,1g16.10,2x,1g16.10)') factor*fauxc
!      enddo
!!!!!!!!!!! END TEST

      enddo !lc (next multipole)
     
      if (verb.ge.4) then
c ---------- Nuclear
!      write(10,800) '#','Jp1=',xjp1,'parity1=',ipar1,'Jp2=',xjp2,
!     & 'parity2=',ipar2,'ener1=',e1,'ener2=',e2
c ---------- Coulomb
!      write(11,800) '#','Jp1=',xjp1,'parity1=',ipar1,'Jp2=',xjp2,
!     & 'parity2=',ipar2,'ener1=',e1,'ener2=',e2
c ---------- nuclear + coulomb
      write(12,800) '#','Jp1=',xjp1,'parity1=',ipar1,'Jp2=',xjp2,
     & 'parity2=',ipar2,'ener1=',e1,'ener2=',e2
      endif


      i1=nint(dabs(xjp1-xjp2))
!      i2=nint(xjp1+xjp2)
      i2=min(nint(xjp1+xjp2),lambmax)
      if (verb.ge.4) then
      do irad=1,nrad3
!      r2=xrad2(irad)
      r2=rstep+rstep*dble(irad-1) 
      write(12,900) r2, (Ff(m1,m2,i,irad),i=i1,i2)
!      write(10,900) r2, (Fn(m1,m2,i,irad),i=i1,i2)
!      write(11,900) r2, (Fc(m1,m2,i,irad),i=i1,i2)
! TEST I4 -> I5
500   format(i4,3f8.4,i4,2f4.0,2i4,a35)
!500   format(i5,3f8.4,i4,2f4.0,2i4,a35)
800   format (a,2(a,(f8.4),2x,a,i3,2x),2(a,(f8.4),2x),/)
820   format (a,a5,a,a5,2(a,(f8.4),2x))
900   format (1f8.3,2x,100(e12.4,2x,e12.4,4x))
      enddo !irad
!      write(10,*) '&'
!      write(11,*) '&'
      write(12,*) '&'
      if (verb.ge.4)  write(120,*)'&'
      endif
      enddo 
      enddo
      enddo
      enddo
      write(*,*)
      if (verb.ge.2) then
!        write(*,*)'fort.10: nuclear formfactors'
!        write(*,*)'fort.11: Coulomb formfactors'
        write(*,*)'fort.12: nuclear + Coulomb formfactors'      
      endif ! verb

c     -----------------------------------------------------------
c     Write states & energies in Fresco format
c     ----------------------------------------------------------- 
      open(ken,file="states.fr",status='unknown')     
      i=0
      icopyt=0
      icpot=1
      ibandt=0
      et=0.d0
      jt=jtgs ! inherited from module 'sistema'
      write(ken,*) m1
      do n=1,jpsets
      rjp   =jpiset(n)%jtot    ! xjp(n)
      ibandp=jpiset(n)%partot  ! ipar(n)
      do ie=1,jpiset(n)%nex    ! np(n)
      i=i+1
      ep=energ(n,ie)-energ(1,1)
      if (i.eq.2) icopyt=1
      write(ken,11511) rjp,ibandp,ep,icpot
      if(icopyt==0) write(ken,1153) jt,ibandt,et
      if(icopyt.ne.0) write(ken,1154) icopyt
      enddo !ie
      enddo !n

11511 format(' &States jp=',f4.1,' ptyp=',i2,' ep=',f8.4,'  cpot=',i3)
 1153 format('         jt=',f4.1,' ptyt=',i2,' et=',f8.4,' /')
 1154 format('         copyt=',i2,' /')
 1155 format(' Skipping state #',i3,' J/pi',f4.1,i2,' and Ex=',f8.4)
c     ----------------------------------------------------------------


c      total = etime(elapsed)
c      print*,'total=',total,'user=',elapsed(1),
c     .'system=',elapsed(2)
      call flush(ken)
      call flush(kfr)
      close(ken); close(kfr)
!      deallocate(fn,fc)
      end


! From scatcc.f90
      subroutine test_cont(nset,nchan,inc,ecm,wf)
      use nmrv,only:nch,ech,vcoup,hort
      use constants
      use channels, only: jpiset, ql,qspl, qj, qspj,qjc,exc,cindex
      use sistema
      use potentials, only: ccmat
      use wfs, only: nr,dr,rvec
      implicit none
c     ------------------------------------------------------------
      logical :: info
      integer ir,n,nchan,method
      integer nset,inc,partot
      real*8  ecm,z12,jtot,rm,factor,r0,vscale
      complex*16 :: phase(nchan),smat(nchan)
      complex*16 :: wf(nchan,nr)
c     ------------------------------------------------------------
      real*8 :: eph(nchan)
c     ------------------------------------------------------------
      rm=av*ac/(av+ac)
      factor=(2*amu/hc**2)*rm

      jtot    =jpiset(nset)%jtot
      partot  =jpiset(nset)%partot
      vscale  =jpiset(nset)%vscale
      
      
      if (nchan.ne.jpiset(nset)%nchan) then
         write(*,*)'Wrong number of channels passed to test_cont!'
      endif
      nch=nchan

      if (allocated(vcoup)) deallocate(vcoup)
      if (allocated(ech))   deallocate(ech)
      allocate(vcoup(nch,nch,1:nr)) ! changed in v2.2

      vcoup(:,:,:)=0d0 
      allocate(ech(nch))
      ql(1:nchan)    =jpiset(nset)%lsp(1:nchan)  
      ech(1:nch)     =jpiset(nset)%exc(1:nch)
      qj(1:nchan)    =jpiset(nset)%jsp(1:nch)
      qjc(1:nchan)   =jpiset(nset)%jc(1:nch)
      cindex(1:nchan)=jpiset(nset)%cindex(1:nch)

      call coefmat(nset,nch)
      vcoup(:,:,1:nr)=ccmat(:,:,1:nr)/factor

      r0=rvec(1)
      method=4       ! enhanced Numerov as used in Fresco
      hort  =0       ! no stabilization
      info  =.false. ! silent output
!     call schcc(nch,ecm,zv*zc,inc,ql,factor,dr,r0,nr-1,wf,phase,smat)
!      call schcc(nch,ecm,zv*zc,inc,ql,factor,dr,r0,
!     & nr,wf,phase,smat,info)
      call schcc_erwin(nch,ecm,zc*zv,inc,ql,factor,dr,r0,
     & nr,wf,phase,smat,method,info,eph)
c to get continuum wfs 'normalized' as <k|k'> = delta(k-k')
!      wf(:,:)=wf(:,:)*sqrt(2./pi)

!      if ((ecm.gt.0.5).and.(ecm.lt.1)) then 
!      write(500,*)'# WF wf for ecm,ql(:)=',ecm,ql(1:nchan)
!      do ir=1,nr
!         write(500,'(1f8.3,2x,50f12.8)') rvec(ir),
!     &       (real(wf(n,ir)),n=1,nch)
!      enddo
!      write(500,*)'&'
!      endif
      deallocate(ccmat,vcoup)
      end subroutine


! From scatcc.f90
      subroutine zfun_clean(nch,ql,kch2,ir,ym,zm)
      use nmrv,only   : conv,vcoup,h,debug,rvec
      use memory, only: tzfun
      implicit none
      integer:: l,ir,ich,is,nch,ql(nch)
      real*8 :: r,kch2(nch)
      real*8 :: ti,tf
      real,parameter:: eps=1e-6
      complex*16,intent(in)  :: ym(nch,nch)
      complex*16,intent(out) :: zm(nch,nch)
      complex*16             :: gmat(nch,nch)
      debug=.false.
      gmat(:,:)=0d0
      call cpu_time(ti)
!      r=h*ir
      r=rvec(ir)

      if (r.lt.eps) r=eps
      do ich=1,nch
      l=ql(ich)
      gmat(ich,ich)= kch2(ich)-l*(l+1)/r**2  
      enddo
      
      gmat(1:nch,1:nch)=gmat(1:nch,1:nch)
     &                 -conv*vcoup(1:nch,1:nch,ir)
      zm=-matmul(gmat,ym)

      call cpu_time(tf)
      tzfun=tzfun+(tf-ti)
      return
      end


! From scatcc.f90
      subroutine qrerwin(yi,zi,zm,wi,wm,y,n,ir,irmin,nr)
      use nmrv, only: h,rvec
      implicit none
      logical    debug
      integer    nr
      complex*16 y(n,n,nr)
c ... for zgeqrf
      integer info,lda,lwork,n,m,i,is,ich,ir,irmin
      complex*16 a(n,n),tau(n),work(2*n),rt(n,n)
      complex*16,dimension(n,n):: yi,ym,zi,zm,wi,wm
c ... for ztrsm
      complex*16 alpha
      integer ldb
      character diag,side,transa,uplo
c ...
      integer   :: nrhs, ipiv(n)
      character*1  trans
      EXTERNAL  ZGETRF, ZGETRS

c ---- TEST      
      integer ndim,lworkx
      complex*16 aux(3,3),taux(3),raux(3,3),workx(6)
      complex*16 zaux(n,n),waux(n,n)
c .................................................
      lwork=2*n
      debug=.false.
     

c ... QR factorization YI=Q.R 
c     (Q=orthogonal matrix; R= triangular UPPER matrix)
      rt=yi
      call ZGEQRF( N, N, RT, n, TAU, WORK, LWORK, INFO )
      if (info.ne.0)  then
         write(*,*)'zgeqrf failed; info=',info
      endif


c ... Apply transformation backwards
c       Y(old)= Y(new) R 
c     Solve:
c        X    *  RT  = Y(old)
c     and make: Y(new)= X 
      
c ... X*op( A ) = alpha*B
      side  ='R' ! A acts on the right
      uplo  ='U' ! R is upper triangle matrix
      transa='N' ! no transpose
      diag  ='N' !??????????
      m=n
      alpha =1d0
      lda   =n   ! leading dim of A=RT
      ldb   =n   ! leading dim of B=Y(old)

c Apply transformation to YI,YM
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,YI,LDB)
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,YM,LDB)
!      y(:,:,ir  )=yi(:,:) ! new (orthogonolized) solution at ir
!      y(:,:,ir-1)=ym(:,:) ! idem at ir-1
      
 
c Apply the same transformation to ZI,ZM
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,ZI,LDB)
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,ZM,LDB)

      if (debug) then
      write(*,*)'yi(new) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(yi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

        write(*,*)'|yp x yp| (after QR)'
        do is=1,n
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(yi(:,is),yi(:,ich)))
     &      / abs(dot_product(yi(:,is),yi(:,is)))
     &      / abs(dot_product(yi(:,ich),yi(:,ich))), ich=1,n) 

        enddo
       endif !debug



      RETURN  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      do i=irmin,ir
      yi(:,:)=y(:,:,i)
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,RT,LDA,YI,LDB)
      y(:,:,i)=yi(:,:) ! override 
      enddo ! ir
      yi(:,:)=y(:,:,ir  ) ! new (orthogonolized) solution at ir
      ym(:,:)=y(:,:,ir-1) ! idem at ir-1






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c ... Recalculate ZI functions for new Y solutions applying the inverse transformation
c       YI(new) = W ZI(new)   -> ZI(new)
c       YM(new) = W ZM(new)   -> ZM(new) 
c     Use Lapack to solve:
c       A * X  =  alpha*B
c       W * Z  =  1    *Y

      write(*,*)'Zi (old) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

c ... First, compute LU factorization of matrix W to get IPIV
      call ZGETRF(n, n, wi, n, ipiv, info)
      if (info /= 0) then 
       print*, 'zgetrf returned error',info
      endif

c ... Second, use IPIV & zgetrs to solve system  W*Z= Y  -> Z 
      trans='N';  ldb=n;  nrhs=n
      zi(:,:)= yi(:,:) 
      call zgetrs(trans,n,nrhs,wi,n,ipiv,zi,ldb,info)

      if (info /= 0) then
        print*,'zgetrs failed at qrerwin!'
      endif

      write(*,*)'Zi (new) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

c ... Repeat for  YM(new) = W ZM(new) 
      call ZGETRF(n, n, wm, n, ipiv, info)
      if (info /= 0) stop 'zgetrf returned error'
      zm(:,:)= ym(:,:) 

      call zgetrs(trans,n,nrhs,wm,n,ipiv,zm,ldb,info)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      Return

c ......................... TEST ........................

      if (1>2) then
      write(*,*)'R (triangular) matrix:'
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(RT(ich,is)), is=1,min(15,n))
      enddo

      write(*,*)'yi(old) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(y(ich,is,ir)), is=1,min(15,n))
      enddo
      write(*,*)' '

      endif


! TEST THAT Y=W*Z
      if (1>2) then 
      write(*,*)'y(old) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(y(ich,is,ir)), is=1,min(15,n))
      enddo
      write(*,*)' '

      write(*,*)'yi(old) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(yi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

!! TEST that Z=W^-1*Y
      write(*,*)'Zi (arg) at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zi(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

      waux=wi
      call ZGETRF(n, n, waux, n, ipiv, info)
      if (info /= 0) then
       write(*,*)'zgetrf returned error',info; stop
      endif

c ... Use IPIV & zgetrs to solve system  W*X= Y   (X=Z)
      trans='N';  ldb=n;  nrhs=n
      zaux(:,:)= yi(:,:) 
      call zgetrs(trans,n,nrhs,waux,n,ipiv,zaux,ldb,info)

      write(*,*)'Zi = W^-1*Y at ir,r=',ir
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(zaux(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

      if (info /= 0) then
          stop 'Matrix inversion failed!'
      else
       write(*,*)'zgetrs info=',info
       write(*,*)'ipiv=',ipiv
      endif


!!! TEST  Y= W*Z
      YM=MATMUL(WM,ZM)
      write(*,*)'W*ZM at ir-1=',ir-1
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(ym(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '

      ym=y(:,:,ir-1)
      write(*,*)'y(old) at ir-1',ir-1
      do ich=1,min(15,n)
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(ym(ich,is)), is=1,min(15,n))
      enddo
      write(*,*)' '
      endif

      if (1>2) then
        write(*,*)'|yi x yi| (before QR)'
        do is=1,n
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(y(:,is,ir),y(:,ich,ir)))
     &      / abs(dot_product(y(:,is,ir),y(:,is,ir)))
     &      / abs(dot_product(y(:,ich,ir),y(:,ich,ir))), 
     &     ich=1,n) 

        enddo
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (0>1) then
      aux(1,1)=12; aux(1,2)=-51; aux(1,3)=4
      aux(2,1)=6;  aux(2,2)=167; aux(2,3)=-68
      aux(3,1)=-4; aux(3,2)=24; aux(3,3)=-41
      lda=3
      lworkx=6
      raux=aux
      call ZGEQRF( lda, lda, raux, lda, TAUX, WORKX, LWORKX, INFO )
      if (info.ne.0)  then
         write(*,*)'zgeqrf failed; info=',info
      endif

      write(*,*)'R matrix:'
      do ich=1,3
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(raux(ich,is)), is=1,lda)
      enddo

c ... X*op( A ) = alpha*B
      side  ='R'
      uplo  ='U'
      transa='N'
      diag  ='N' !??????????
      m=n
      alpha =1d0
      lda   =3
      ldb   =3
      m     =lda
      call ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,m,ALPHA,RAUX,LDA,AUX,LDB)

      write(*,*)'Q matrix:'
      do ich=1,3
        write(*,'(5x,i3,50g14.5)') ich,
     &  (real(aux(ich,is)), is=1,lda)
      enddo


        write(*,*)'|yp x yp| (before QR)'
        do is=1,3
        write(*,'(5x,50g14.5)') 
     &    ( abs(dot_product(aux(:,is),aux(:,ich)))
     &      / abs(dot_product(aux(:,is),aux(:,is)))
     &      / abs(dot_product(aux(:,ich),aux(:,ich))), ich=1,3) 

        enddo

      stop
      endif

c .......................................................
      end subroutine


! From scatcc.f90
	SUBROUTINE POTWF2(COUPL,ZI,FI,V,NEQS,MAXB,MAXCH,R12)
	IMPLICIT NONE
	REAL*8 R12
	INTEGER NEQS,MAXB,K,J,IT,MAXCH
	COMPLEX*16 COUPL(MAXB,MAXB),ZI(MAXB,MAXCH),FI(MAXB,MAXCH),
     X		   V(MAXB,MAXCH),C,ZERO
	PARAMETER(ZERO = (0d0,0d0))
C					VECTORISED MATRIX MULTIPLIES
C					ALLOWS SKIPS IF ZERO COUPLING

         DO 24 K=1,NEQS
	     DO 24 J=1,NEQS
             C = COUPL(K,J) * R12
	      if(C/=ZERO) FI(1:NEQS,K) = FI(1:NEQS,K) - C * ZI(1:NEQS,J)
 24       CONTINUE

         DO 34 K=1,NEQS
	    V(:,K)  = ZERO
            DO 34 J=1,NEQS
	      C = COUPL(K,J)
	      if(C/=ZERO) V(1:NEQS,K) = V(1:NEQS,K) + C * FI(1:NEQS,J)
34       CONTINUE
	return
	end


! From scatcc.f90
      subroutine matching(ecm,z12,nch,ql,lmax,inc,
     & rmatch,y,wf,phase,smat,show)
      use nmrv, only: nr,mu,ech,h,conv
      use constants, only : e2,pi
      use globals, only: verb,debug
      implicit none
      logical :: sing,show
      integer ich,is,inc,ir,klog
      integer nch,ql(nch),l,linc,lmax,ifail,ie,nmatch

      real*8, dimension(0:lmax):: f,g,gp,fp
      real*8, dimension(1:lmax+1):: f1,fp1
      real*8 :: kron,rmatch
      real*8 :: ecm,tkch,eta,krm,z12,kch(nch)
      real*8 :: phr,cph(0:lmax) !!!! changed to lmax check !!!!
      real*8 :: fpmax,small,acc8,flux,magn(nr)

      complex*16:: wf(nch,nr)
      complex*16:: yd(nch,nch),yaux(5)
      complex*16:: y(nch,nch,nr)
      complex*16:: ch,chd,mat(2*nch,2*nch+1),det
      complex*16:: ci,test,yfac,tmat
      complex*16:: phase(nch),phc
      complex*16:: smat(nch),delta,scoul,anc,svel

      acc8 = epsilon(acc8);
      fpmax = huge(acc8)**0.8d0 
      SMALL=1D0/FPMAX
      ci=(0d0,1d0)
      pi=acos(-1d0) 
      debug=.false.
!      show=.false.
      klog=99
      magn(:)=0      

      if (nch.eq.0) then
       write(*,*)'Matching: nch=0! Abort'; stop     
      endif
      mat(:,:)=0d0
!      mat(1:nch,1:nch)=y(1:nch,1:nch,nr)
      mat(1:nch,1:nch)=y(1:nch,1:nch,nr-2) ! matrix solution at r=rmatch
      nmatch=nr-2                        ! matching done at ir=nr-2 

c
c derivatives at matching point
c 
      if (debug) write(99,*)'Derivatives at  rmatch=',rmatch
      do ich=1,nch
      do is=1,nch
!      yaux(1:5)= y(ich,is,nr-2:nr+2)
      yaux(1:5)= y(ich,is,nr-4:nr) ! use last 5 points for derivative
!      yd(ich,is)=deriv1(yaux,h,nr+2,nr)
      yd(ich,is)=(yaux(1)
     &   -8*yaux(2)+8*yaux(4)-yaux(5))/(h*12)
      mat(nch+ich,is)= yd(ich,is)
c3      mat(nch+ich,is)= y(ich,is,nr-1)
      enddo
      if (debug) write(99,'(5x,50f12.6)') (yd(ich,is),is=1,nch) 
      enddo
      
      linc=ql(inc)
      do ich=1,nch
      l=ql(ich) 
      tkch=ecm+ech(inc)-ech(ich) 

      kch(ich)=sqrt(conv*abs(tkch))
      krm=kch(ich)*rmatch
      eta=conv*z12*e2/kch(ich)/2.

      if (tkch.gt.0d0) then   ! open channel
       call coulph(eta,cph,linc)
       phc=exp(ci*cph(linc))

!       write(190,'(6f10.4)') tkch,kch(ich),eta,krm,cph(linc)
       call coul90(krm,eta,0d0,l,f,g,fp,gp,0,IFAIL)


       if (ifail.ne.0) then 
       write(*,*) 'coul90: ifail=',ifail; stop
       endif
       ch  = dcmplx(g(l),f(l))  * (0.,.5)   ! outgoing (i/2) H(+)
       chd = dcmplx(gp(l),fp(l)) * (0.,.5) ! (i/2) H'(+)
c3      call coul90(kch(ich)*(rm-h),eta,0d0,l,f,g,fp,gp,0,IFAIL)
c3      chd = dcmplx(g(l),f(l)) * (0.,.5) 
c lhs
      mat(ich,nch+ich)    = ch
      mat(nch+ich,nch+ich)= chd*kch(ich)
c rhs
      mat(ich,2*nch+1)    =-kron(ich,inc)*conjg(ch) 
      mat(nch+ich,2*nch+1)=-kron(ich,inc)*
     &                    conjg(chd)*kch(ich) 
      else                      ! closed channel
       IE = 0  
       call whit(eta,rmatch,kch(ich),tkch,l,f1,fp1,ie)
       ch = -f1(l+1)  ! Whittaker
       chd= -fp1(l+1) ! derivative 
!       write(0,*) 'l,k,f',l,kch(ich),f1(l+1)
c lhs
      mat(ich,nch+ich)    = -ch
      mat(nch+ich,nch+ich)= -chd*kch(ich)
c rhs
      mat(ich,2*nch+1)    =0. 
      mat(nch+ich,2*nch+1)=0. 
      endif
c2       mat(ich,nch+ich)    = -(g(l)+ci*f(l))
c2       mat(nch+ich,nch+ich)= -(gp(l)+ci*fp(l))*kch(ich)
c3        mat(ich,nch+ich)    = ch
c3        mat(nch+ich,nch+ich)    = chd
c2      mat(ich,2*nch+1)    =kron(ich,inc)*f(l)
c2      mat(nch+ich,2*nch+1)=kron(ich,inc)*fp(l)*kch(ich)
c3      mat(ich,2*nch+1)    =-kron(ich,inc)*conjg(ch)
c3      mat(nch+ich,2*nch+1)=-kron(ich,inc)*conjg(chd)

      enddo
      call GAUSS5(2*nch,2*nch,mat,1,sing,det,small,debug)
      if (debug) write(klog,'(5x,"Coefficients:",100f8.3)') 
     & (mat(ich,2*nch+1),ich=1,nch)

!      IF(SING) then
!      DO 605 Ir=1,nr
!605   IF(MAGN(Ir).NE.0.0) MAGN(Ir) = LOG10(MAGN(Ir))
!      WRITE(*,610) (MAGN(Ir),Ir=1,nr)
!610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
!      endif

c
c S-matrix and phase-shifts
c 
      flux=0d0
      phase(1:nch)=0.
      if (show.and.(verb.ge.1)) 
!     &  write(*,'(8x,"S-matrix (wo/with velocity factors):")') 
     &   write(*,'(a40,6x,a45)')
     &   "S-matrix (WITHOUT velocity factors)", 
     &   "S-matrix (WITH velocity factors)"
      do ich=1,nch
      tkch=ecm+ech(inc)-ech(ich) 
c2      tmat=mat(nch+ich,2*nch+1)
c2      smat(ich)=(tmat-1)*(0.,-0.5)
      if (tkch.gt.0) then                                   ! open channel
!commented v2.1b
!      smat(ich)=mat(nch+ich,2*nch+1)*ci**(QL(INC)-QL(ICH))
       smat(ich)=mat(nch+ich,2*nch+1)
       svel=smat(ich)*sqrt(kch(ich)/kch(inc)) ! S-matrix with velocity factors
       phase(ich)=(0.,-0.5) * LOG(mat(nch+ich,2*nch+1))
       test=(0.,-0.5) * LOG(mat(nch+ich,2*nch+1)) 
       if (((show).and.(ich.eq.inc)).or.(verb.ge.1)) then
!       write(*,500) ich,smat(ich),svel,abs(svel),test*180/pi
       write(*,500) ich,smat(ich),svel,abs(svel) !,test*180/pi
       endif
       flux=flux + abs(svel)**2
500    format(10x,"Chan. #",i3,5x,
     &    "S=(",1f10.6,",",1f10.6,")", 5x, 
     &    "S=(",1f10.6,",",1f10.6,") ->  |S|=",1f9.6)
!     &    ,"delta=(",1f8.3,",",1f8.3,")")

      else                                                 ! closed channel
        smat(ich)=mat(nch+ich,2*nch+1)
        write(*,500) ich,smat(ich),svel,abs(svel),test*180/pi
        if (verb.gt.3) then
        anc=mat(nch+ich,2*nch+1)
        if (debug)write(*,510) ich,anc,abs(anc)
510     format(5x,"Channel",i3,3x,
     &    "ANC=(",1g12.6,",",1g12.6,") -> |ANC|=",1g12.6)
      endif
      endif
      enddo !ich

      phase(1:nch)=phase(1:nch)*180/pi
!      write(45,'(1f10.3,3x,10g14.6)') ecm,
!    &  (phase(ich),ich=1,nch) 
       
      write(*,'(10x,"=> Unitarity=",1f12.8)')flux


c
c scattering wavefunctions
c
      wf(1:nch,1:nr)=(0d0,0d0)
      if(abs(smat(inc)).gt.1e-20) then
        phr = (0.,-0.5) * LOG(smat(inc))
        yfac=exp(-ci*phr)
      endif

      do ich=1,nch
      do is=1,nch
      wf(ich,1:nr)=wf(ich,1:nr)+mat(is,2*nch+1)*y(ich,is,1:nr) 
!      wf(ich,1:nr)=wf(ich,1:nr)+xmat(is,1)*y(ich,is,1:nr)
      enddo
      enddo

!      wf(:,:)=dsqrt(2d0/pi)*phc*yfac*wf(:,:)
ccc      wf(:,:)=dsqrt(2d0/pi)*yfac*wf(:,:)
! TEST
!      wf(:,:)=dsqrt(2d0/pi)*wf(:,:)
! TO MATCH WITH FRESCO FORT.17 , NO FACTORS ARE NEEDED HERE!


c      do ir=1,nr
c      write(51,'(2x,1f6.3,2x,10f14.8)')ir*h,(wf(ich,ir),ich=1,nch)
c      enddo

      end subroutine


! From scatcc.f90
      subroutine matching2(ecm,z12,nch,ql,lmax,inc,
     & rmatch,y,wf,phase,smat)
      use nmrv, only: nr,mu,ech,h,conv
      use constants, only : e2,pi
      use globals, only: verb
      implicit none
!!!!!!!!!!!!!!!!!!!  LAPACK    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer  nch
      integer  iter,info
      integer,  parameter:: nrhs=1
!      integer          ::  n=2*nch
      integer    ipiv(2*nch)
      double precision :: rwork(2*nch) 
      complex    swork(2*nch*(2*nch+nrhs))
      complex*16 mat(2*nch,2*nch),bmat(2*nch,nrhs),xmat(2*nch,nrhs),
     &           work(2*nch)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      logical :: debug !sing,show
      integer ich,is,inc,ir,klog
      integer ql(nch),l,linc,lmax,ifail,ie

      real*8, dimension(0:lmax):: f,g,gp,fp
      real*8, dimension(1:lmax+1):: f1,fp1
      real*8 :: kron,rmatch
      real*8 :: ecm,tkch,eta,krm,z12,kch(nch)
      real*8 :: phr,cph(0:lmax) !!!! changed to lmax check !!!!
      real*8 :: fpmax,small,acc8,flux,magn(nr)

      complex*16:: wf(nch,nr)
      complex*16:: yd(nch,nch),yaux(5)
      complex*16:: y(nch,nch,nr+2)
      complex*16:: ch,chd,det
      complex*16:: ci,test,yfac,tmat
      complex*16:: phase(nch),phc
      complex*16:: smat(nch),delta,scoul,anc,svel

      acc8 = epsilon(acc8);
      fpmax = huge(acc8)**0.8d0 
      SMALL=1D0/FPMAX
      ci=(0d0,1d0)
      pi=acos(-1d0) 
      debug=.false.
      klog=99
      magn(:)=0      

      if (nch.eq.0) then
       write(*,*)'Matching: nch=0! Abort'; stop     
      endif
      mat(:,:)=0d0
      mat(1:nch,1:nch)=y(1:nch,1:nch,nr-2)

c
c derivatives at end point
c 
      if (debug) write(99,*)'Derivatives at end point rm=',rmatch
      do ich=1,nch
      do is=1,nch
!      yaux(1:5)= y(ich,is,nr-2:nr+2)
      yaux(1:5)= y(ich,is,nr-4:nr) ! use last 5 points for derivative
!      yd(ich,is)=deriv1(yaux,h,nr+2,nr)
      yd(ich,is)=(yaux(1)-
     &   8*yaux(2)+8*yaux(4)-yaux(5))/(h*12)
      mat(nch+ich,is)= yd(ich,is)
c3      mat(nch+ich,is)= y(ich,is,nr-1)
      enddo
      if (debug) write(99,'(5x,50f12.6)') (yd(ich,is),is=1,nch) 
      enddo
      
      linc=ql(inc)
      do ich=1,nch
      l=ql(ich) 
      tkch=ecm+ech(inc)-ech(ich) 

      kch(ich)=sqrt(conv*abs(tkch))
      krm=kch(ich)*rmatch
      eta=conv*z12*e2/kch(ich)/2.

      if (tkch.gt.0d0) then   ! open channel
       call coulph(eta,cph,linc)
       phc=exp(ci*cph(linc))
       call coul90(krm,eta,0d0,l,f,g,fp,gp,0,IFAIL)


       if (ifail.ne.0) then 
       write(*,*) 'coul90: ifail=',ifail; stop
       endif
       ch  = dcmplx(g(l),f(l))  * (0.,.5)   ! outgoing (i/2) H(+)
       chd = dcmplx(gp(l),fp(l)) * (0.,.5) ! (i/2) H'(+)
c3      call coul90(kch(ich)*(rm-h),eta,0d0,l,f,g,fp,gp,0,IFAIL)
c3      chd = dcmplx(g(l),f(l)) * (0.,.5) 
c lhs
      mat(ich,nch+ich)    = ch
      mat(nch+ich,nch+ich)= chd*kch(ich)
c rhs
      bmat(ich,1)      =-kron(ich,inc)*conjg(ch) 
      bmat(nch+ich, 1) =-kron(ich,inc)*conjg(chd)*kch(ich) 
      else                      ! closed channel
       IE = 0  
       call whit(eta,rmatch,kch(ich),tkch,l,f1,fp1,ie)
       ch = -f1(l+1) !* (0.,.5) ! Whittaker*i/2
       chd= -fp1(l+1)!* (0.,.5) ! derivative 
c lhs
      mat(ich,nch+ich)    = -ch
      mat(nch+ich,nch+ich)= -chd*kch(ich)
c rhs
      bmat(ich,    1) =0. 
      bmat(nch+ich,1) =0. 
      endif
      enddo


c LAPACK SOLVER
      call zcgesv(2*nch,nrhs,mat,2*nch,ipiv,bmat,2*nch,xmat,2*nch,
     &           work,swork,rwork,iter,info)

c
c S-matrix and phase-shifts
c 
      flux=0d0
      phase(1:nch)=0.
      if (verb.gt.0)
!     &  write(*,'(8x,"S-matrix (wo/with velocity factors):")') 
     &   write(*,'(a40,6x,a45)')
     &   "S-matrix (WITHOUT velocity factors)", 
     &   "S-matrix (WITH velocity factors)"
      do ich=1,nch
      tkch=ecm+ech(inc)-ech(ich) 
c2      tmat=mat(nch+ich,2*nch+1)
c2      smat(ich)=(tmat-1)*(0.,-0.5)
      if (tkch.gt.0) then                                   ! open channel
      smat(ich)=xmat(nch+ich,1)
      svel=smat(ich)*sqrt(kch(ich)/kch(inc)) ! S-matrix with velocity factors
      phase(ich)=(0.,-0.5) * LOG(smat(ich))
      test=(0.,-0.5) * LOG(smat(ich)) 
      if (((verb.eq.0).and.(ich.eq.inc)).or.(verb.ge.1)) then
       write(*,500) ich,smat(ich),svel,abs(svel),test*180/pi
      endif
      flux=flux + abs(svel)**2
500   format(10x,"Chan. #",i3,5x,
     &    "S=(",1f10.6,",",1f10.6,")", 5x, 
     &    "S=(",1f10.6,",",1f10.6,") ->  |S|=",1f9.6,3x, 
     &    "delta=(",1f8.4,",",1f8.4,")")

      else                                                  ! closed channel
        if (verb.gt.3) then
        anc=xmat(nch+ich,1)
        smat(ich)=anc
        if (debug)write(*,510) ich,anc,abs(anc)
510     format(5x,"Channel",i3,3x,
     &    "ANC=(",1g12.6,",",1g12.6,") -> |ANC|=",1g12.6)
      endif
      endif
      enddo !ich
      phase(1:nch)=phase(1:nch)*180/pi
       
      write(*,'(10x,"=> Unitarity=",1f12.8)')flux
      
c
c scattering wavefunctions
c
      wf(1:nch,1:nr)=(0d0,0d0)
      if(abs(smat(inc)).gt.1e-20) then
        phr = (0.,-0.5) * LOG(smat(inc))
        yfac=exp(-ci*phr)
      endif

      do ich=1,nch
      do is=1,nch
      wf(ich,1:nr)=wf(ich,1:nr)+xmat(is,1)*y(ich,is,1:nr)
      enddo
      enddo

!      wf(:,:)=dsqrt(2d0/pi)*phc*yfac*wf(:,:)
ccc      wf(:,:)=dsqrt(2d0/pi)*yfac*wf(:,:)
! TEST
!      wf(:,:)=dsqrt(2d0/pi)*wf(:,:)
! TO MATCH WITH FRESCO FORT.17 , NO FACTORS ARE NEEDED HERE!

!      do ir=1,nr
!      write(51,'(2x,1f6.3,2x,10f14.8)')ir*h,(wf(ich,ir),ich=1,nch)
!      enddo

      end subroutine


! From scatcc.f90
      subroutine match_real(ecm,z12,nch,nop,ql,lmax,inc,n1,y,wf,show)
      use nmrv, only: nr,mu,ech,h,conv,rvec
      use constants, only : e2,pi
      use globals, only: verb,debug
      implicit none
      logical :: sing,show
      integer ich,icl,is,inc,ir,klog,nd,n1,n2
      integer nch,nop,ql(nch),l,linc,lmax,ifail,ie
      integer ncl,iop,sgn,i,j
c     --------------------------------------------------------------
      real*8, dimension(0:lmax):: f,g,gp,fp
      real*8, dimension(1:lmax+1):: f1,fp1
      real*8 :: kron,r1,r2,r,whitp,raux,d1,d2,recdet,det
      real*8 :: ecm,tkch,eta,krm1,krm2,z12,kch(nch)
      real*8 :: phr,cph(0:lmax) !!!! changed to lmax check !!!!
      real*8 :: fpmax,small,acc8,flux,magn(nr)
      real*8:: yd(nch,nch),yaux(5),y1,y2
      real*8:: y(nch,nch,nr),yop(nch,nop,nr),yor(nch,nop,nr)
      real*8:: ccoef,dcoef,gaux(nch,nr),amod,rcoef,eph(nch)
      real*8:: cmat(nch-nop,nch),dmat(nch-nop,nch)
      real*8:: caux(nch-nop,nch-nop),cinv(nch-nop,nch-nop)
      real*8:: kmat(nch-nop,nop)
      real*8:: subnmat(nop,nop),nmat(nop,nop)
c     ----------------------------------------------------------------
      complex*16, parameter::iu=cmplx(0.,1.)
      complex*16:: wf(nch,nch,nr) 
      complex*16:: ch1,ch2,chp1,chp2,zaux
      complex*16:: ci,test,cj
      complex*16:: afas,cnorm,acoef
      complex*16:: anc
      complex*16:: amat(nop,nch),bmat(nop,nop)
c     ---------------------------------------------------------------
!      debug=.true.
      show=.false.

      acc8 = epsilon(acc8);
      fpmax = huge(acc8)**0.8d0 
      SMALL=1D0/FPMAX
      pi=acos(-1d0) 
      klog=99
      magn(:)=0      
      nd=6
      n2=n1-nd
      ncl=nch-nop
 
!      print*,'Match_real for Ecm=',ecm,' nop=',nop

      if (nch.eq.0) then
       write(*,*)'Matching: nch=0! Abort'; stop     
      endif

      amat(:,:)=0d0; cmat(:,:)=0; dmat(:,:)=0
      linc=ql(inc)

      do ich=1,nch
        tkch=ecm+ech(inc)-ech(ich) 
        kch(ich)=sqrt(conv*abs(tkch))
      enddo

      do is=1,nch
      do ich=1,nch
      l=ql(ich) 
      tkch=ecm+ech(inc)-ech(ich) 
!      kch(ich)=sqrt(conv*abs(tkch))
      eta=conv*z12*e2/kch(ich)/2.

      r1=rvec(n1)
      r2=rvec(n2)

      krm1=kch(ich)*r1
      krm2=kch(ich)*r2

      y1=y(ich,is,n1)
      y2=y(ich,is,n2)
   
      IF (tkch.gt.0d0) THEN   ! open channel
       call coulph(eta,cph,linc)
!       phc=exp(ci*cph(linc))
!       write(190,'(6f10.4)') tkch,kch(ich),eta,krm1,cph(linc)
       call coul90(krm1,eta,0d0,l,f,g,fp,gp,0,IFAIL)
       if (ifail.ne.0) then 
       write(*,*) 'coul90: ifail=',ifail; stop
       endif
       ch1  = dcmplx(g(l),f(l))     !  H(+)(r2)
       call coul90(krm2,eta,0d0,l,f,g,fp,gp,0,IFAIL)
       ch2 = dcmplx(g(l),f(l))      !  H(+)(r2)

!       afas=(y1*ch2-y2*ch1)/(y1*conjg(ch2)-y2*conjg(ch1))
!       amod=(y1*ch2-y2*ch1)/(y1*conjg(ch2)-y2*conjg(ch1))

!       acoef=amod*sqrt(afas)

! Alt method
       acoef=(y1*ch2-y2*ch1)/(conjg(ch2)*ch1-conjg(ch1)*ch2)
       acoef=acoef*2*iu*sqrt(kch(ich)*kch(inc))
       amat(ich,is)=acoef

!       print*,'Sol=',is,' Chan=',ich,'(open) a=',amod*sqrt(afas),acoef
!       acoef=(y2*conjg(ch1)-y1*conjg(ch2))/
!     &  (ch2*conjg(ch1)-ch1*conjg(ch2))
!     &  *2*iu*sqrt(kch(ich)*kch(inc))
!       print*,'Sol=',is,' Chan=',ich,'(open) a=',amod*sqrt(afas),acoef

      ELSE !......................................... closed channel
       IE = 0  
       call whit(eta,r1,kch(ich),tkch,l,f1,fp1,ie)
       ch1 = f1(l+1) ! Whittaker at r1
       call whit(eta,r2,kch(ich),tkch,l,f1,fp1,ie)
       ch2=  f1(l+1) ! Whittaker at r2  
 
       chp1=whitp(eta,r1,kch(ich)) ! exp[+k r1]
       chp2=whitp(eta,r2,kch(ich)) ! exp[+k r2]

       ccoef=(y1*ch2-y2*ch1)/(chp1*ch2-chp2*ch1)*kch(ich)
       dcoef=(y1*chp2-y2*chp1)/(ch1*chp2-ch2*chp1)*kch(ich)

       icl=ich-nop
       cmat(icl,is)=ccoef
       dmat(icl,is)=dcoef
       
!      print*,'Sol=',is,' Chan=',ich,' (closed) c,d=',ccoef,dcoef
       ENDIF ! OPEN/CLOSED CHANNELS      
      enddo ! ich
      enddo ! is

c ... Open channels solutions
      yop(:,1:nop,1:nr)=y(:,1:nop,1:nr)
!      do iop=1,nop
!       write(96,*)'#Open regular solution ',iop,nop,' Ecm=',ecm
!       do ir=1,nr
!       r=rvec(ir)
!       write(96,'(1f8.3,3x,50g14.6)')r,(yop(ich,iop,ir),ich=1,nch)
!       enddo !ir 
!       write(96,*)'&'
!      enddo !iop


c ... Regularize in case of closed channels
      if (ncl.gt.0) then
        cinv(1:ncl,1:ncl)=cmat(1:ncl,nop+1:nch)
        caux(:,:)=0.0
        call cmatin(cinv,caux,ncl)
! TEST      
!        caux=matmul(cinv,cmat(1:ncl,nop+1:nch))
!        write(*,*)'cinv.c=',caux 

        kmat=-matmul(cinv,cmat(1:ncl,1:nop))

c     Regular solutions
       do iop=1,nop
       gaux(:,1:nr)=0.0
       do icl=1,ncl
        gaux(:,1:nr)= gaux(:,1:nr)
     &              + kmat(icl,iop)*y(:,nop+icl,1:nr) ! CHECH!!!
       enddo ! closed channels
       yop(:,iop,1:nr)=yop(:,iop,1:nr)+gaux(:,1:nr)
      enddo ! open channels 
      endif ! open/closed channels?

c B-matrix
!      print*,'B-matrix: ncl=',ncl
      do iop=1,nop
      do ich=1,nop
       zaux=0
       do icl=1,ncl
       zaux=zaux+kmat(icl,iop)*amat(ich,nop+icl)
       enddo !icl
       bmat(ich,iop)=amat(ich,iop)+zaux
!       write(*,*)'iop,ich,b=',iop,ich,bmat(ich,iop)
      enddo !ich
      enddo !iop

c Overlap matrix
      nmat(:,:)=0
      do iop=1,nop
      do ich=1,nop
    
      nmat(iop,ich)=(sum(conjg(bmat(iop,1:nop))*bmat(ich,1:nop))
     +             +sum(bmat(iop,1:nop)*conjg(bmat(ich,1:nop))))/2.0
!      do   j=1,nop  
!      nmat(iop,ich)=nmat(iop,ich)
!     &             +0.5*(conjg(bmat(iop,j))*bmat(ich,j)
!     &             +bmat(iop,j)*conjg(bmat(ich,j)))
!      enddo ! j
!      write(*,*)'iop,ich,N=',iop,ich,nmat(iop,ich)
      enddo !ich
      enddo !iop

c orthogonalize and normalize solutions using Gram-Schmidt 
      if (1>2) then
! GS in determinant form
      rcoef=1./sqrt(nmat(1,1))
      yor(:,:,:)=0.
      yor(1:nch,1,:)=yop(1:nch,1,:)*rcoef
      sgn=+1
      do iop=2,nop
        d1=recdet(nmat(1:iop,1:iop),iop)
        d2=recdet(nmat(1:iop-1,1:iop-1),iop-1)
        do is=1,iop
        subnmat(1:iop-1,1:is-1)   =nmat(1:iop-1,1:is-1)
        subnmat(1:iop-1,is:iop-1) =nmat(1:iop-1,is+1:iop) 
        det=recdet(subnmat(1:iop-1,1:iop-1),iop-1)
        yor(1:nch,iop,:)=yor(1:nch,iop,:)+sgn*det*yop(1:nch,is,:)
        sgn=-sgn
        enddo !is
        yor(1:nch,iop,:)=yor(1:nch,iop,:)/sqrt(d1*d2)
       
      if (debug) then
      write(97,*)'# Regular solution ',iop,'of',nop,' Ecm=',ecm
      do ir=1,nr
      r=rvec(ir)
      write(97,'(1f8.3,3x,50g14.6)')r,(yor(ich,iop,ir),ich=1,nch)
      enddo !ir 
      write(97,*)'&'
      endif
      enddo !iop
 
      ELSE 
! Conventional GS form 
      rcoef=1./sqrt(nmat(1,1))
!      yor(1:nch,1:nop,1:nr)=yop(1:nch,1:nop,1:nr)
      wf(1:nch,1,:)=yop(1:nch,1,:)*rcoef
      if (debug) then
      write(97,'(a,i2,a,i2,a,1f8.3)')'# Regular solution ',1,
     & 'of',nop,' Ecm=',ecm
      do ir=1,nr
      r=rvec(ir)
      write(97,'(1f8.3,3x,50g14.6)')r,(wf(ich,1,ir),ich=1,nch)
      enddo !ir 
      write(97,*)'&'
      endif
      do iop=2,nop
        wf(1:nch,iop,:)=yop(1:nch,iop,:)
        do is=1,iop-1
        wf(:,iop,:)=wf(:,iop,:)-nmat(iop,is)*yop(:,is,:)
        enddo ! is
c norm
        raux=0
        do i=1,iop
        do j=1,iop
        ci=-nmat(iop,i)/nmat(i,i)
        cj=-nmat(iop,j)/nmat(j,j)
        if (i.eq.iop) ci=1.
        if (j.eq.iop) cj=1.
        raux=raux + ci*cj*nmat(i,j)
        enddo !j
        enddo !i
        wf(:,iop,:)=wf(:,iop,:)/sqrt(raux)
!!!!!!!!
      if (debug) then
      write(97,'(a,i2,a,i2,a,1f8.3)')'# Regular solution ',
     &  iop,'of',nop,' Ecm=',ecm
      do ir=1,nr
      r=rvec(ir)
      write(97,'(1f8.3,3x,50g14.6)')r,(wf(ich,iop,ir),ich=1,nch)
      enddo !ir 
      write(97,*)'&'
      endif
!!!!!!!!
      enddo ! iop
      endif
      end subroutine


! From scatcc.f90
      subroutine gsorto(yp,a,nch)
      implicit none
      integer nch,ic,icp,n,i,j
      complex*16 yp(nch,nch),a(nch,nch),yaux(nch)

      a(:,:)=0;
      a(1,1)=1
      do j=2,nch      ! columns
         a(j,j)=1
         yaux(:)=0
         do i=1,j-1   ! rows
         a(i,j)=-dot_product(yp(:,i),yp(:,j))
     &          /dot_product(yp(:,i),yp(:,i))
         yaux(:)=yaux(:) +  yp(:,i)*a(i,j)
         enddo ! i
         yp(:,j)=yp(:,j)+ yaux(:)
      enddo !j
      return
      end subroutine


! From scatcc.f90
      subroutine mgsorto(yp,a,nch)
      implicit none
      integer nch,ic,icp,n,i,j
      complex*16 yp(nch,nch),s(nch,nch),a(nch,nch)
!!! TEST
      character*1 jobvs,sort
      logical select
      integer lda,ldvs,sdim,lwork,info
      complex*16 w(nch),vs(nch,nch),work(2*nch)
      double precision rwork(nch)
      logical bwork(nch)
!!!!!
      a(:,:)=0; a(1,1)=1
      do j=2,nch      ! columns
         a(j,j)=1
!         yaux(:)=0
         do i=1,j-1   ! rows
         a(i,j)=-dot_product(yp(:,i),yp(:,j))
     &          /dot_product(yp(:,i),yp(:,i))
!         yaux(:)=yaux(:) +  yp(:,i)*a(i,j)
         yp(:,j)=yp(:,j)-yp(:,i)*a(i,j)
         enddo ! i
!         yp(:,j)=yp(:,j)+ yaux(:)
      enddo !j

      RETURN
c overlap matrix 
      do ic=1,nch
      do icp=1,nch
         s(ic,icp)=dot_product(yp(:,ic),yp(:,icp))
      enddo ! icp
      enddo ! ic


      lwork=2*nch
      lda=nch
      ldvs=nch
      jobvs='N'
      sort='N'
      call ZGEES(JOBVS, SORT, SELECT, NCH, A, LDA, SDIM, W, VS,
     $                  LDVS, WORK, LWORK, RWORK, BWORK, INFO )
*
      if (info.eq.0) then
      write(*,'(5(i3,2x,1g12.6,3x))') (n,w(n),n=1,nch)
      else
         write(*,*)'error: info=',info
      endif
            
      end subroutine


! From scatcc.f90
      SUBROUTINE ESMOOV(NC,QUAD,QUAD1,PHZ,EN,NUME,NFTA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C***********************************************************************
C     Downloaded from:  http://www.tampa.phys.ucl.ac.uk/rmat/old/eigenp.f
C
C     ESMOOV was formerly part of ETABLE and carries out smoothing
C     of eigenphases using second order finite differences
C
C***********************************************************************
C
      DIMENSION QUAD(NC),QUAD1(NC),PHZ(NC,NUME),EN(NUME)
      DATA PI/3.1415926535897932D+00/,ZERO/0.D0/,ONE/1.0D0/
C
C     INITIALIZE ARRAY QUAD1 WITH THE MULTIPLES OF PI WHICH ARE
C     TO BE ADDED TO THE EIGENPHASES AT THE FIRST ENERGY
C
      scale = one
      ieav = 0
      DO 3 I=1,NC
      QUAD1(I)=QUAD(I)
    3 CONTINUE
C
C     CHECK CONTINUITY OF EIGENPHASE PARAMETERS AS A FUNCTION OF
C     SCATTERING ENERGY BY THE USE OF SECOND-ORDER DIVIDED DIFFERENCES
C
      DO 190 NO=1,NC
C
C     INITIATE SMOOTHING ; CONTINUITY OF FIRST THREE ENERGY POINTS
C     IS ASSUMED
C
      Q1=QUAD1(NO)*PI
C
      II=3
      E1=EN(1)
      E2=EN(2)
      E3=EN(3)
      P1=PHZ(NO,1)+Q1
      P2=PHZ(NO,2)+Q1
      P3=PHZ(NO,3)+Q1
      PHZ(NO,1)=P1
      PHZ(NO,2)=P2
      PHZ(NO,3)=P3
C
C     FORM DIVIDED DIFFERENCES
C
      FD1=(P2-P1)/(E2-E1)
      FD2=(P3-P2)/(E3-E2)
      SD2=(FD2-FD1)/(E3-E1)
C
C     BEGIN ITERATION OVER ENERGY POINTS
C
    5 II=II+1
      IF(II .GT. NUME) GO TO 190
      Q1=QUAD1(NO)*PI
      P2=P3
      P3=PHZ(NO,II)+Q1
C
      E1=E2
      E2=E3
      E3=EN(II)
      FD1=FD2
      FD2=(P3-P2)/(E3-E2)
      SD1=SD2
      SD2=(FD2-FD1)/(E3-E1)
      EPS1=SCALE*DABS(SD1)
C
C     CHECK CONTINUITY
C
      IF(DABS(SD2-SD1) .GT. EPS1) GO TO 7
C
C     CONTINUITY CHECKS CORRECTLY
C
    6 PHZ(NO,II)=P3
      GO TO 5
C
C                     **************************
C
C     CORRECTION SEQUENCE :
C
C     (1) ADD QUADRANT CORRECTION OF PI
C
    7 P3=P3+PI
      FD2=(P3-P2)/(E3-E2)
      SD2=(FD2-FD1)/(E3-E1)
      IF(DABS(SD2-SD1) .GT. EPS1) GO TO 9
C
      QUAD1(NO)=QUAD1(NO)+ONE
      GO TO 6
C
C     (2) SUBTRACT QUADRANT CORRECTION OF PI
C
    9 P3=P3-PI-PI
      FD2=(P3-P2)/(E3-E2)
      SD2=(FD2-FD1)/(E3-E1)
      IF(DABS(SD2-SD1) .GT. EPS1) GO TO 11
C
      QUAD1(NO)=QUAD1(NO)-ONE
      GO TO 6
C
C     (3) CHECK WHETHER CHANNELS HAVE BEEN SWAPPED
C
   11 NO1=NO+1
      IF(NO1 .GT. NC) GO TO 20
C
      DO 18 NOX=NO1,NC
      IPT=1
      NOP=NOX
C
      P3=PHZ(NOP,II)+Q1
   12 FD2=(P3-P2)/(E3-E2)
      SD2=(FD2-FD1)/(E3-E1)
C
C     JUMP OUT IF SEARCH SUCCEEDS
C
      IF(DABS(SD2-SD1) .LE. EPS1) GO TO 25
C
      if(ipt.eq.1) then
C
C     TRY ADDING PI
C
        P3=P3+PI
        IPT=2
        GO TO 12
C
C     TRY SUBTRACTING PI
C
      else if(ipt.eq.2) then
        P3=P3-PI-PI
        IPT=3
        GO TO 12
      endif
C
   18 CONTINUE
C
C     SMOOTHING FAILURE
C
   20 WRITE(NFTA,22)NO,II,E3
   22 FORMAT(' CONTINUITY ERROR ENCOUNTERED IN CHANNEL ',I3,
     1       ' EIGENPHASE AT ENERGY E(',I3,') =',D16.8)
      P3=PHZ(NO,II)+Q1
      GO TO 5
C
C     SWAP EIGENPHASE COLUMNS NO AND NOP BEGINNING AT ENERGY POINT II
C
   25 DO 27 IEP=II,NUME
      SAVE=PHZ(NO,IEP)
      PHZ(NO,IEP)=PHZ(NOP,IEP)
      PHZ(NOP,IEP)=SAVE
   27 CONTINUE
C
C     REVERSE QUAD1 VALUES
C
      SAVE=QUAD1(NOP)
      QUAD1(NOP)=QUAD1(NO)
      QUAD1(NO)=SAVE
      if(ipt.eq.2) then
        QUAD1(NO)=QUAD1(NO)+ONE
      else if(ipt.eq.3) then
        QUAD1(NO)=QUAD1(NO)-ONE
      endif
C
      PHZ(NO,II)=P3
C
      WRITE(NFTA,32)NO,NOP,II
   32 FORMAT(' CHANNELS NO =',I3,' AND NOP =',I3,' SWAPPED BEGINNING',
     1       ' AT ENERGY POINT IEN =',I4)
      GO TO 5
C
  190 CONTINUE
      RETURN
      END


! From scatcc.f90
      SUBROUTINE ZLACG_Q( N, X, INCX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
*     ..
*     .. Array Arguments ..
      COMPLEX*32         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLACGV conjugates a complex vector of length N.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The length of the vector X.  N >= 0.
*
*  X       (input/output) COMPLEX*16 array, dimension
*                         (1+(N-1)*abs(INCX))
*          On entry, the vector of length N to be conjugated.
*          On exit, X is overwritten with conjg(X).
*
*  INCX    (input) INTEGER
*          The spacing between successive elements of X.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IOFF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG
*     ..
*     .. Executable Statements ..
*
      IF( INCX.EQ.1 ) THEN
         DO 10 I = 1, N
            X( I ) = CONJG( X( I ) )
   10    CONTINUE
      ELSE
         IOFF = 1
         IF( INCX.LT.0 )
     $      IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = CONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZLACGV
*
      END


! From ham.f90
	subroutine identity(ndim)
	use hmatrix,only:idmat
	implicit none
	integer:: ndim,n,m
	if(allocated(idmat)) deallocate(idmat)
	allocate(idmat(ndim,ndim))
	idmat(:,:)=0d0
	do n=1,ndim
	do m=1,ndim
	if (n.eq.m) idmat(n,m)=1d0
	enddo
	enddo
	end subroutine


! From fragpot.f90
      subroutine potdef3(vc,qmax,rnc,vdef,betar)
      use wfs, only: nr,rvec
      implicit none
      integer i,k,qmax,nu,mmultipole
!      integer ir,icall
      real*8 vc(nr)
      real*8 u,r,rsp,rnc,pi,sh,cns,alpha,fival,vcr!,vdefr
      parameter(alpha=0d0)
      real*8, dimension(0:qmax):: pl,vdef,def,c
      real*8 sp(9),w(9) 
!!! TEST
      real*8 ws, p1,p2,p3,betar
!   
      data (sp(i),i=1,5) /-.96816024,-.8360311,-.61337143,-.32425342,0./  &
     &,(w(i),i=1,5) /.08127439,.18064816,.26061070,.31234708,.33023936 /
!       save rsp,sp,w,c
       save rsp,sp,w
      pi=acos(-1d0)

!!!! TEST
      def(:)=0d0
      def(2)=betar
!!!!!      
      rsp = 1.0/sqrt(pi)
      do 2 i=1,4
      sp(5+i) =-sp(5-i)
2     w(5+i)  = w(5-i)
      do 15 k=0,qmax
15    c(k) = sqrt((2.*k+1.)/4.)
!	write(10,*) 'potdef initialised for qmax =',qmax
!      endif
      do 18 k=0,qmax
18    vdef(k) = 0.
         pl(0) = 1.
      do 50 nu=1,9
         u = sp(nu)
         pl(1) = u
         do 20 k=2,qmax
20       pl(k) = ((2*k-1)*u*pl(k-1) - (k-1)*pl(k-2))/dble(k)
      sh  = 0.0
      do 25 k=2,qmax
25       sh = sh + c(k)*rsp * pl(k) * def(k)
         r = rnc - sh
         vcr=0d0
         if (r.le.rvec(nr)) then 
         vcr=fival(r,rvec,vc,nr,alpha)
         endif
	 do 40 k=0,qmax
            cns = w(nu) * pl(k) * c(k)
          vdef(k)=vdef(k) + cns * vcr
40      continue
50      continue
	return
      end


! From solvecc.f90
      subroutine write_cdcc_wf(icc,wf,nch,nr)
      use channels, only: jptset
      implicit none
      integer, parameter:: kwf=85
      integer :: ich,nch,nr,ir,icc,partot
      integer :: lsp,jsp,jc,parc,lpt,idx
      real*8 :: jtot,exc,jp,ex,jlp,kcm
      complex*16 wf(nch,nr)

!     nch   =jptset(icc)%nchan  
      partot=jptset(icc)%partot
      jtot  =jptset(icc)%jtot  

      write(kwf,350) icc,nch,jtot,partot
350   format(2x,"CC set:", i4, " Chans:",i4,
     &  " JTOT=",1f6.1, " Parity:",i2)

      do ich=1,nch
       idx   =jptset(icc)%idx(ich)
       lpt   =jptset(icc)%l  (ich)
       jp    =jptset(icc)%jp (ich)
       jlp   =jptset(icc)%jlp(ich)
       ex    =jptset(icc)%exc(ich)
       kcm   =jptset(icc)%kcm(ich)
  
      write(kwf,354) ich,idx,lpt,jlp,jp,ex,kcm
354   format(2x," Channel:",i4, " IEX: ",i3,
     &       " QN:",i3,1f5.1,1f5.1, 
     &       " Ex=",1f6.3," Kcm=",1f6.3)
!      write(kwf,'(6g16.6)') (wf(nch,ir),ir=1,nr)
      enddo !ich
      end subroutine


! From solvecc.f90
      subroutine solvecc(icc,nch,inc,nr,nlag,ns)
      use xcdcc   ,  only: hcm,elab,smats,rvcc,method
      use channels,  only: jptset
      use nmrv,      only: vcoup, ech,ql,hort,cutr
      use constants, only: hc,amu
      use sistema
      use globals  , only: verb,debug
!      use memory   , only: tcc
      implicit none
      integer    :: ir,n,nr,icc,nch,inc,nlag,ns
      complex*16 :: phase(nch),wf(nch,nr),smat(nch)
      real*8:: ecm
      real*8:: rmass,factor,rstart
!!!!!!!!!! TEST: delete me when I work 
      real*8 vcoul,rc
      logical:: test=.false., info=.true.
      real*8 ::v0,r0,a0,w0,ri,ai,a13,vaux,waux,ws,r
      real*8 ::ti,tf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      write(*,'(//,3x," Numerov integration:")') 
      debug=.false.
      rmass=mp*mt/(mp+mt)
      ecm=elab*mt/(mp+mt)
      factor=(2*amu/hc**2)*rmass


!      allocate(wf(nch,nmatch))
      if (allocated(ech))   deallocate(ech)
      if (allocated(ql))    deallocate(ql)
      allocate(ech(nch))
      allocate(ql(nch))

      ql(1:nch) =jptset(icc)%l(1:nch)  
      ech(1:nch)=jptset(icc)%exc(1:nch)+jptset(icc)%ext(1:nch)


!!!!! TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (test) then
      write(*,*)'TESTING NUMEROV INTEGRATION'
c *** Channels and potentials
      ql(1:nch) =2
      ech(1:nch)=0.0
      inc=1
      zp=1; zt=6
      ecm=18.46;
      rc=1.25

      v0=-50; r0=1.2; a0=0.5
      w0=-10; ri=1.2; ai=0.5
c we do not consider the 1st point (r=rmin)
c becasuse in scatcc the 1st point is r=h
      rstart=hcm
      do ir=1,nr
        a13=12**.333333       
        r=ir*hcm 
        vaux=ws(r,v0,r0*a13,a0)
        waux=ws(r,w0,ri*a13,ai)
        vcoup(1,1,ir)=cmplx(vaux,waux) + vcoul(r,zp,zt,rc*a13) 
!        write(40,'(1f8.3,2x,50f12.8)') r,vcoup(1,1,ir)
      enddo      
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c Solve eqns
c-------------------------------------------------------
      if (debug) then
        write(*,*) 'Calling Numerov with:'
        write(*,*)'    Incoming channel=',inc
        write(*,*)'    (2*amu/hc**2)*rmass=',factor
        write(*,*)'    Ecm=',ecm
        write(*,*)'    l-values=',ql(1:nch)
        write(*,*)'    Channel energies:',ech(1:nch)
        write(*,*)'    Z1*Z2=',zp*zt
      if (.not.allocated(vcoup)) then
        write(*,*)'solvecc: vcoup not allocated!'
        stop
      endif
      endif
c -----------------------------------------------------
      call cpu_time(ti)
      rstart=rvcc(1)
      cutr=-10
      select case(method)
      case(0) ! predictor-corrector (Baylis & Peels)
        call schcc(nch,ecm,zp*zt,inc,ql,factor,hcm,
     &  rstart,nr,wf,phase,smat,info)

      case(1,2,3) ! Enhanced Numerov (Thorlacious & Cooper) / Raynal
        call schcc_ena(nch,ecm,zp*zt,inc,ql,factor,hcm,
     &  rstart,nr,wf,phase,smat,method,info)

      case(4)     ! Enhanced Numerov (Fresco version of T&C)
        call schcc_erwin(nch,ecm,zp*zt,inc,ql,factor,hcm,
     &  rstart,nr,wf,phase,smat,method,info)

      case(5)     ! R-matrix method (P. Desc. subroutine)
        write(0,*)'calling scchcc_rmat'
        write(*,*) '=== R-matrix parameters: nlag=',nlag,' ns=',ns
        call schcc_rmat(nch,ecm,zp*zt,inc,ql,factor,hcm,
     &  rstart,nr,wf,phase,smat,info,nlag,ns)

      case default
        write(*,*)'Method',method,' not valid!'
        stop
      end select

      call cpu_time(tf)
      if (verb.ge.3) 
     &write(*,'(5x,"[ CC solved in",1f8.3," sec ]")') tf-ti
!      tcc=tcc+finish-start
      smats(icc,inc,1:nch)=smat(1:nch)

c     Deallocate variables
!      if (allocated(vcoup)) deallocate(vcoup)
      deallocate(ql)
      end subroutine


! From solvecc.f90
      subroutine int_ff()
c *** --------------------------------------------------------------------      
      use xcdcc,   only: hcm,ffc,lamax,nex,nrcc,jpch,Ff,lambmax,
     &                   nrad3,rstep,parch,rmaxcc,rvcc,realwf
      use channels,only: jpiset,jpsets
      use wfs     ,only: idx
      use memory
      use globals ,only: verb
      implicit none
      integer :: ir,irr
      integer :: m1,m2,nff,lam,par1,par2
      complex*16, pointer:: faux(:)
      real*8  :: rv(nrad3),xpos,rfirst
      real*8 :: jp1,jp2,r,ymem
      real*8, parameter:: alpha=0d0
      complex*16 cfival,caux,fic,caux2,ffc4
c     ----------------------------------------------------
     
      write(*,'(//,3x, "** INTERPOLATE FORMFACTORS **" )')

      lamax=2*maxval(jpiset(:)%jtot)

c ... Memory requirements
      ymem=nex*nex*(lamax+1)*nrcc*lc16/1e6
      if (verb.ge.1) write(*,190) ymem
190   format(5x,"[ FF require", 1f8.2," Mbytes ]")

      if (allocated(ffc)) deallocate(ffc)
      allocate(ffc(nex,nex,0:lamax,1:nrcc))
      ffc=0
      do ir=1,nrad3
       rv(ir)= rstep + rstep*(ir-1) ! CHECK
      enddo     
 
      write(*,170) nrcc, hcm,rmaxcc,hcm
170   format(/,5x,"=>  Interpolating F(R) in grid with ",i10, ' points:',
     &           2x,"[ Rmin=",1f5.2,2x," Rmax=",
     &           1f6.1,1x," Step=",1f6.3, " fm ]",/)

      if (rmaxcc.gt.rv(nrad3)) then
        write(*,'(3x,"*** WARNING: Coupled eqns integrated up to",f6.1,
     &  1x,"fm, but formfactors calculated for R<",1f6.1," fm" )')
     &  rmaxcc,rv(nrad3) 
      endif 

      nff=0
      do m1=1,nex
      par1= parch(m1)
      jp1 = jpch(m1)
! Changed in v2.4 !!!!!!!!!!!!!!!!!!!!
!      do m2=m1,nex
      do m2=1,nex
      if (realwf.and.(m2.lt.m1)) cycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      jp2 = jpch(m2)
      par2= parch(m2)
       do lam= nint(dabs(jp1-jp2)),min(nint(jp1+jp2),lambmax) !nint(jp1+jp2)
       if (par1*par2*(-1)**lam<0) cycle
       nff=nff+1
       faux =>  Ff(m1,m2,lam,:)
c      interpolate & store in integration grid
       do ir=1,nrcc
!        r=hcm*ir ! grid used in scattcc
!        r=dble(ir-1)*hcm ! grid used in scattcc (Changed in v2.2)
        r=rvcc(ir)
        if (r.gt.rv(nrad3)) cycle ! DO NOT EXTRAPOLATE 
!        caux=cfival(r,rv,faux,nrad3,alpha)
        xpos=(r-rstep)/rstep
        caux=FFC4(xpos,faux,nrad3)
        ffc(m1,m2,lam,ir)=caux
!        if (abs(caux).gt.1e20) then 
!          write(90,*)'# ** INT_FF: F(R)=',caux, 
!     &  'for ir,lam,m1,m2',ir,lam,m1,m2
!         write(90,'(1f8.3,2x,2g16.5)') (rv(irr),faux(irr),irr=1,nrad3)
!         write(90,*)'&'
!         stop
!        endif
       enddo ! ir
      enddo ! lam
      enddo !m2
      enddo !m1
  
      write(*,*) '=> there are', nff,' radial formfactors '
      nullify(faux)
      deallocate(ff)
      end subroutine


! From solvecc.f90
      subroutine int_ff_tdef()
c *** --------------------------------------------------------------------      
      use xcdcc,   only: hcm,ffcn,ffcc,lamax,nex,nrcc,jpch,Fc,nlambhipr,
     &                   nrad3,rstep,parch,rmaxcc,rvcc,Fn
      use channels,only: jpiset,jpsets
      use wfs     ,only: idx
      use memory
      use globals ,only: verb
      implicit none
      integer :: ir,irr
      integer :: m1,m2,nff,lam,par1,par2
      complex*16, pointer:: faux1(:),faux2(:)
      real*8  :: rv(nrad3)
      real*8 :: jp1,jp2,r,ymem,yffc
      real*8, parameter:: alpha=0d0
      complex*16 cfival,caux,ffc4
c     ----------------------------------------------------
     
      write(*,'(//,3x, "** INTERPOLATE FORMFACTORS TARG DEF **" )')


      ymem=nex*nex*(nlambhipr)*nrcc*lc16/1e6
      if (verb.ge.1) write(*,190) ymem
190   format(5x,"[ FF require", 1f8.2," Mbytes ]")

      if (allocated(ffcn)) deallocate(ffcn)
      if (allocated(ffcc)) deallocate(ffcc)
      allocate(ffcn(nex,nex,nlambhipr,1:nrcc),
     & ffcc(nex,nex,nlambhipr,1:nrcc))
      ffcn=0
      ffcc=0
             
      do ir=1,nrad3
       rv(ir)= rstep + rstep*(ir-1) ! CHECK
      enddo     
      write(*,170) nrcc, hcm,rmaxcc,hcm
170   format(/,5x,"=>  Interpolating F(R) in grid with ",i4, ' points:',
     &           2x,"[ Rmin=",1f5.2,2x," Rmax=",
     &           1f6.1,1x," Step=",1f6.3, " fm ]",/)

      if (rmaxcc.gt.rv(nrad3)) then
        write(*,'(3x,"*** WARNING: Coupled eqns integrated up to",f6.1,
     &  1x,"fm, but formfactors calculated for R<",1f6.1," fm" )')
     &  rmaxcc,rv(nrad3) 
      endif 

      nff=0
      do m1=1,nex
      par1= parch(m1)
      jp1 = jpch(m1)
      do m2=m1,nex
      jp2 = jpch(m2)
      par2= parch(m2)
       do lam= 1,nlambhipr !nint(jp1+jp2)
       nff=nff+1
       faux1 =>  Fc(m1,m2,lam,:)
       faux2 =>  Fn(m1,m2,lam,:)
c      interpolate & store in integration grid

       do ir=1,nrcc
!        r=hcm*ir ! grid used in scattcc
!        r=dble(ir-1)*hcm ! grid used in scattcc (Changed in v2.2)
        r=rvcc(ir)
        if (r.gt.rv(nrad3)) cycle ! DO NOT EXTRAPOLATE 
!        caux=cfival(r,rv,faux1,nrad3,alpha)
        yffc=(r-rv(1))/(rv(2)-rv(1))
        caux=ffc4(yffc,faux1,nrad3)
        ffcc(m1,m2,lam,ir)=caux
!        caux=cfival(r,rv,faux2,nrad3,alpha)
        yffc=(r-rv(1))/(rv(2)-rv(1))
        caux=ffc4(yffc,faux2,nrad3)
        ffcn(m1,m2,lam,ir)=caux
!        if (abs(caux).gt.1e20) then 
!          write(90,*)'# ** INT_FF: F(R)=',caux, 
!     &  'for ir,lam,m1,m2',ir,lam,m1,m2
!         write(90,'(1f8.3,2x,2g16.5)') (rv(irr),faux(irr),irr=1,nrad3)
!         write(90,*)'&'
!         stop
!        endif
       enddo ! ir
      enddo ! lam
      enddo !m2
      enddo !m1
  
      write(*,*) '=> there are', nff,' radial formfactors '
      nullify(faux1,faux2)
      deallocate(fc,fn)
      end subroutine


! From continuum.f90
      subroutine continuum!(nchan)
      use globals, only: mu12,egs,kin,written,verb
      use constants
      use sistema
      use scattering
      use wfs, only:nr,dr,rvec
      use channels
      use parameters, only: maxchan
      use potentials, only: ccmat
      use factorials
!!!! scatcc
      use nmrv,only:nch,ech,vcoup
!!!!!
      implicit none
      logical info,writewf
      integer n,ir,ich,ichp,nskip,klog !,method
      real*8:: r,r0,dk,kcont,econt,k,ecm, t1,t2
      real*8:: phase(nchan),aux,phadd(nchan)
      real*8,allocatable:: ph(:,:)
c  bincc variables
      logical tres,tdel,tkmat,tknrm
      integer:: nf   ! number of formfactors?
      integer:: isc,iil,ik,jset
      integer:: maxn,nnk(nchan)
      real*8:: anc,k2(nchan),conv,rm,bphase(nk)
      complex*16:: y(nr,nchan),etap(nchan),wf(nchan,nr),smat(nchan)
      real*8:: fmscal
c  end bincc 
      namelist/scatwf/ emin, emax,ifcont,nk,il,ilout,writewf,jset
      il=0
      ili=1
      klog=99
      method=4
      nbas= 50
      ns=   1
      read(kin,nml=scatwf)
      
      if (method.eq.0) method=4
      if ((.not.ifcont).or.(nk.eq.0)) return

      write(*,*)
      write(*,*)' SCATTERING WAVE FUNCTIONS:'
      write(*,*)
      if (il.eq.0) then
        write(*,*)'Max. incoming channel not specified: assuming nchan'
        il=nchan
      else 
        if (il.gt.nchan) then
        write(*,*)' # incoming channels bigger than #outcoming: 
     &assuming IL=nchan'
        il=nchan
      endif        

      if (il.lt.0) then
        write(*,*)' Max. incoming channel negative: computing only -il'
        il=-il
        ili=il
      endif
      endif 
      
      write(*,*)'R-matrix solution with nbas,ns=',nbas,ns

      write(*,*)
      write(*,100)nk,emin,emax
100   format(3x,"Generating", i4, " continuum wf(s) for"
     &  ," [Emin=",1f6.3," Emax=",1f7.2," MeV]")

      if (method.eq.5) write(*,*) "(using R-matrix routine)" 
      call factorialgen(100)

      call coefmat(jset,nchan)

      allocate(wfcont(nk,nchan,nchan,nr))
      wfcont(:,:,:,:)=0d0

      rm=av*ac/(av+ac)
      conv=(2*amu/hc**2)*rm
      kmin=sqrt(2d0*mu12*abs(emin))/hc
      kmax=sqrt(2d0*mu12*abs(emax))/hc


!! ! for bincc --------------------------------------------------
c1      nf=1;  isc=2      
c1      tres=.false.; tdel=.true.; tkmat=.false.; tknrm=.false.
c1      maxn=nr+1


c1      do n=1,nchan
c1       k2(n)=abs(egs)+exc(n) ! separation energy for each channel
c1       etap(n)= 2*e2*zv*zc*rm*amu/hc**2  !=2 k eta !! CHECK!!!!
c1       k2(n)=k2(n)*conv     
c1      enddo      
c1      do iil=ili,il
c1     call BINCC2(Y,ccmat,nf,NCHAN,k2,iil,conv,bphase,isc,kmin,kmax,nk,
c1     & etap,NR-1,dr,pcon,tres,tdel,tkmat,tknrm,lmax+1,maxn,anc)
c1      write(*,*)'il=',iil,' done'
c1      enddo
c1      RETURN !!!!!!! TEST BINCC 
c ----------------------------- end bincc call ---------------------
 

      allocate(vcoup(nchan,nchan,1:nr))
      allocate(ech(nchan))
      allocate(ph(nk,nchan))
      ech(1:nchan)=exc(1:nchan)
      vcoup(:,:,:)=0d0 
      vcoup(:,:,1:nr)=ccmat(:,:,1:nr)/conv ! changed in v2.2

      deallocate(ccmat)
      if (nk.gt.1) then
         DK = (KMAX - KMIN)/(NK-1)
      else
         DK=(KMAX - KMIN)
      endif
      do iil=ili,il
      nskip=0
      phadd=0
      write(*,'(5x,"o incoming channel:",i3,$)') iil
      call cpu_time(t1)
      do ik=1,nk
       K = KMIN + (IK-1)*DK
       ecm=k**2/conv-ech(iil)
!       write(*,'(".",$)')
       if (verb.ge.5) then 
        write(*,*)' **********************************'
        write(*,150) k**2/conv,ecm,il
150    format(4x,'Calling SCHCC for E=',1f7.3,' and Ecm=',1f7.3, 
     &  ' MeV with incident channel',i2)
       endif
!       call schcc(nchan,ecm,zv*zc,iil,ql,conv,dr,nr-1,wf,phase,smat) ! changed in v2.2
!       wfcont(ik,iil,1:nchan,2:nr)=wf(1:nchan,1:nr-1)
       r0=rvec(1)
       info=.false.
!       method=4
c       call schcc(nchan,ecm,zv*zc,iil,ql,conv,dr,r0,nr,wf,phase,smat,
c     &            info)

       call schcc_erwin(nchan,ecm,zc*zv,iil,ql,conv,dr,r0,
     &                    nr,wf,phase,smat,method,info)

       wfcont(ik,iil,1:nchan,1:nr)=wf(1:nchan,1:nr)

!      write(47,'(1f10.3,3x,10g14.6)') ecm,(phase(n),n=1,nchan)
!      ph(ik,1:nchan)=phase(1:nchan)
!      write(46,'(1f10.3,3x,10g14.6)') ecm,(ph(ik,n),n=1,nchan) 


       if (ecm.gt.0) then
         written(45)=.true.
         if (ik.eq.1) then
          ph(1,1:nchan)=phase(1:nchan)
         else
          do n=1,nchan 
           aux=phase(n)+phadd(n)
           if (aux< ph(ik-1,n)-90) phadd(n)=phadd(n)+180
           if (aux> ph(ik-1,n)+90) phadd(n)=phadd(n)-180
           ph(ik,n)=phase(n)+phadd(n)
           enddo
         endif 
         write(45,'(1f10.3,3x,10g14.6)') ecm,
     &  (ph(ik,n),n=1,nchan) 


         if (iil.eq.ilout) then 
         written(50)=.true.
         write(50,*)'# Continuum wf with Ecm=',Ecm
         do ir=1,nr
         write(50,'(1f8.3,2x,50f12.8)')rvec(ir),
     &       (wfcont(ik,iil,n,ir),n=1,nchan)
         enddo 
         write(50,*)'& '
         endif
      else
        nskip=nskip+1
        if (nskip.le.5) then
        write(*,200) ecm
200     format(4x,'Skipping negative energy:',1f7.3, ' MeV')
        else if (nskip.eq.6) then
        write(*,*)"    (further not shown)"
        endif
      endif
      enddo ! ik
      call cpu_time(t2)
      write(*,'(3x,"(",1f6.3," secs )")') t2-t1
      enddo ! iil
      write(*,*)




!! AMoro: skip creepy part...
      RETURN 


        do iil=ili,il

         do ik=1,nk
          kcont=kmin+(kmax-kmin)*(ik-1)/(nk-1)
          econt=(hc*kcont)**2/2/mu12
! 	  if ((econt.gt.eout).and.(.not.written(51))) then
!          written(51)=.true.
!          write(51,*)'# Continuum wf with e=',econt,'ik=',ik
!          do ir=1,nr
!      write(51,'(1f8.3,2x,50f10.6)')rvec(ir),
!     &(wfcont(ik,iil,n,ir),n=1,nchan) 
!           enddo
!          endif
         enddo
!         endif
		 enddo



!!------------------------------------------------------------We are going to make the creepiest thing we can imagine to solve our problem with ik
!!------------------------------------------------------------I also use it to print out the function
!      ! first we calculate the number corresponding to the first k upon the barrier
!      dk=(kmax-kmin)/nk
!      do iil=ili,il
!
!         nnk(iil)=nint(sqrt(exc(iil)*conv)/dk)
!
!         do ik=nk,1,-1
!         if (ik.gt.nnk(iil)) then
!
!         wfcont(ik,iil,:,:)=wfcont(ik-nnk(iil),iil,:,:)
!
!         else
!         wfcont(ik,iil,:,:)=cmplx(0d0,0d0)
!         endif
!         enddo
!
!		 if (iil.eq.ilout) then
!         do ik=1,nk
!          kcont=kmin+(kmax-kmin)*ik/nk
!          econt=(hc*kcont)**2/2/mu12
!	  if ((econt.gt.eout).and.(.not.written(50))) then
!          written(50)=.true.
!          write(50,*)'# Continuum wf with e=',econt
!          do ir=1,nr
!      write(50,'(1f8.3,2x,50f10.6)')rvec(ir),
!     &(wfcont(ik,iil,n,ir),n=1,nchan) ! changed in v2.0.5
!!     &(dreal(wfcont(ik,iil,n,ir)),n=1,nchan)
!          enddo
!          endif
!         enddo
!         endif
!
!      enddo
      end subroutine         


! From xsections_alpha.f90
      subroutine d2sigma_alpha(nth,dth,thmin,thcut,emin,emax,ncont,
     &                   icore,jsets,fileamp,doublexs,triplexs,phixs)
      use xcdcc,    only: elab,jpch,nex,exch,parch,famps0,binset,rel
      use channels, only: jpiset,jpsets,nchmax,sn,qnc
      use sistema
      use constants
      use parameters, only:  maxchan,maxsets,maxeset
      use globals,  only: debug,written,verb,kin
      use wfs,      only: nr,dr,energ,rvec,wfc !,wbin
      use memory,   only: t3d,lr8
      use scattering, only: method,nbas,ns
      implicit none

      include "omp_lib.h"

c     ----------------------------------------------------------------------- 
      character*40:: line,fileamp
      logical :: writewf=.false.
      logical, parameter:: energy=.true.   ! NEEDED FOR FFC4 INTERPOLATION!! 
      logical :: jsets(maxsets),doublexs,triplexs,phixs
      integer , parameter:: kfam=137  ! file containing amplitudes
      integer , parameter:: eps=1e-3
      integer , parameter:: maxstat=400
      integer :: ic,icore,ncont,idet,ier
      integer ::nsetmax,ir,iset,n,nchan,iam,ich,iecv,nth,iex,ith,mth
      integer ::nmi,nmf,nmfmax,nm,li,lf,im,imp,inc,nearfa
      integer,allocatable::ips(:,:)
c     -------------------------------------------------------------------------
      real*8:: kpt,kcv,krel,dkrdk,wbin
      real*8:: ji,jf,jci,jcf,raux,xsaux,jac,jpi,jpf,jtarg
      real*8:: ecmi,ecmf,kcmi,kcmf,excore,th,thmin,dth,thcut
      real*8:: dec,ebind,ecv,ethr,emin,emax,exc
      real*8:: facK,mv,mc,mucv,mupt,f2t,xstot,sumr
      real*8,allocatable::dsdew(:,:),dsdew_b(:,:,:)
      real*8:: delta(ncont)
      real*8:: ti,tf,hbarc
c     -------------------------------------------------------------------------
      complex*16,allocatable:: wfcont(:,:,:),ampaux(:,:,:) !,gaux(:)
      complex*16,allocatable,target:: gsolap(:,:,:,:)
      complex*16 :: smate(ncont,maxchan),gaux(nr),haux(ncont)
c     -------------------------------------------------------------------------
c     Triple x-sections
c     -------------------------------------------------------------------------
      logical zert,fail3,wrt,kinset(maxsets)
      character*3 :: sys ! 'com' or 'lab'
c    ......................................................................
      integer imu,imo,ival,j,iii,ii,iflag,iroot,iroot1n,iroots,iroot2n
      integer iten,isig,ip,inu,iang,ilval,itphi,iv,itv,itphim,itc,lmax
      integer ind,ien,nchann,niag,nial,ibin,nieg,niel,ik,maxne,nk
      integer, parameter:: ncleb=2000, nord=2
c     .....................................................................
      real*8,allocatable:: angsr(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      real*8,allocatable:: xyt(:,:)
!      real*8 xyt(2,ncont)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8:: Kthr,dk,dkb,kmax,a,b
      real*8:: dot,m,plmrec,ylmc2,plm_nr
      real*8:: dphir,dtc,dtv,phiu,tcl,tvl,tcu,tvu,dtcr,dtvr
      real*8:: acc,degrad,mtot,totp
      real*8:: erel(ncont),ki,kf,ebin
      real*8:: erelmax,tcd,tc,costc,sintc,tvd,tv,costv,sintv
      real*8:: En,Enlow,Enup,dEn,p1L,p2L,phid,phil,dphi,phi,cospv,sinpv
      real*8:: cospc,sinpc,aq,bq,cq,dq,Ec,Ev,eKb,eks
      real*8:: qp(3),kcL(3),kvL(3),pcL(3),pvL(3),tem(3)
      real*8:: ptot(3),ktot(3),bkp(3),kp(3),sigphi(4000),xsig
      real*8:: mkp,co,phKb,phks,co2,si2,si,aleg,mult,xb,xbd,ybo,yb
      real*8:: tmatsq,sc,mu,sig,mo,rnu,cleb,rli
      real*8:: cplm(0:10,21),cgc(ncleb)
c     ..........................................................................
      complex*16 :: resc,caux,tmat_aux,faux,cfival,ffc4
      complex*16:: ampt,xi,f2c,sumps,phbin(maxsets),wkfac(maxsets)
      complex*16:: phask1,phask,phas,ylm(0:10,21)
      complex*16,allocatable:: tmat(:,:,:)
      complex*16,pointer:: overlap(:)
      complex*16,allocatable:: fxyc(:,:)
!      complex*16:: fxyc(10,ncont)
      complex*16,allocatable:: fv(:)
c  *** Relativistic kinematics
! included now in xcdcc module
!      logical  :: rel=.true.
      real*8   :: sinv,etcm,pcmi,pcmf,mpx
      real*8  ::maxtab
      integer  :: imax,jmax
      real*8 :: polar(10)
c     ..........................................................................

!!! TESTING
      real*8 tmatsq2,ampt2
      integer icount,if2c
      complex*16,allocatable:: fin(:,:,:)
      complex*16:: bindk(maxsets,maxeset)
      real*8:: xtab(6),ytab(6),kcmf_cont
      complex*16:: fxytab(6,6),fint2d,fint2dd,fint2db_jin

!MGR---------------------------------------------------------------------
      real*8,allocatable :: xs3body(:,:,:), energ3(:)
! AMM July 19
      real*8,allocatable :: xs3body_phi(:,:,:,:), energ4(:,:,:,:)
!------------------------------------------------------------------------
c ------------------------------------------------------------------------------
      namelist /framework/ sys,idet ! ,atarget not needed
      namelist /gridener/ Enlow,Enup,dEn
      namelist /gridthetac/ tcl,tcu,dtc
      namelist /gridthetav/ tvl,tvu,dtv
      namelist /gridphi/ phil,phiu,dphi
      namelist /polarization/ polar
c ------------------------------------------------------------------------------
      nsetmax=maxval(jpiset(:)%nex) !max number of states within j/pi set
      allocate(fv(nsetmax))
      nmfmax =2*nint(maxval(jpiset(:)%jtot))+1
      mc=ac*amu
      mv=av*amu
      mt=mt*amu
      mp=mp*amu  
      mupt=mp*mt/(mp+mt)
      mucv=mc*mv/(mc+mv)
      hbarc=hc
      nmi=nint(2*jpiset(1)%jtot+1)

      write(*,'(//,3x, "** DOUBLE DIFFERENTIAL CROSS SECTIONS **",/ )')

c     If no amplitudes have been calculated before, we try to read them
c     from an external file, specified by the user 
      if (.not.allocated(famps0).and.(fileamp.ne."")) then
      if (fileamp.eq."") then
          write(*,*)'No amplitudes calculated and no file specified';
          stop
      endif
      write(*,'(3x," Reading f-amps from file ",a20)')fileamp
      open(kfam,file=fileamp,status='old')
      read(kfam,'(a)',err=900) line 
!      write(*,'(a)') line
      read(line,*) jpi,jtarg,jpf,jtarg,mth,nearfa,elab     
!      read(kfam,*,err=900) jpi,jtarg,jpf,jtarg,mth,nearfa,elab
!     & ,a,b
!        write(*,*) 'ie,jpf,mth,elab,a=',iex,jpf,mth,elab,a
      rewind(kfam)

      allocate(famps0(maxstat,mth,nmi*nmfmax))
      iex=0
      famps0(:,:,:)=0
      do iset=1,jpsets  
      nex=jpiset(iset)%nex
      jpf=jpiset(iset)%jtot
      nmf=nint(2.*jpf+1.)
      nm=nmi*nmf
      write(*,'(a,i3,a,i3,a,1f4.1)') ' Expected: iset=',iset, 
     & '  nex=', nex, ' jpf=',jpf
      do n=1,nex
      iex=iex+1
      if (iex.gt.maxstat) stop'increase dimension of famps0' ! temporary solution
!      write(*,*)'expected jpf=',jpf
      read(kfam,'(a)',err=900) line
      read(line,*) jpi,jtarg,jpf,jtarg,mth,nearfa,elab    
!      read(kfam,*) jpi,jtarg,jpf,jtarg,mth,nearfa,elab
!     & ,a,b
      write(*,'(a,i3,a,i3,a,1f4.1)') ' Found: ie=',iex, 
     & ' n=',n,' jpf=',jpf

!      if (iex.eq.1) write(*,*) ' amplitudes are for elab=',elab
      if(nth.ne.mth) then
         write(*,*) 'mismatch in number of angles: specified:',nth,
     &              ' found:',mth
      endif 
      do ith=1,nth
      read(kfam,*) th  
      read(kfam,228) (famps0(iex,ith,iam),iam=1,nm) 
228   format(1P,6E12.4)
      enddo !states
      enddo !ith
      write(*,*)' ' 
      enddo !iset
650   write(*,'(5x,"[f amplitudes read for ",i3," states and",i5,
     & " angles ]")')iex,nth
      endif
c     -------------------------------------------------------------------------
 

      if (.not.rel) then       ! non-relativistic
        ecmi=elab*mt/(mp+mt)
        mupt=mp*mt/(mp+mt)!*amu
        kcmi=sqrt(2*mupt*ecmi)/hc
        write(*,*)' ** USING NON-RELATIVISTIC KINEMATICS **'
      else                     ! relativistic
        write(*,*)' ** USING RELATIVISTIC KINEMATICS **'
        sinv=((mp+mt))**2+2*mt*elab
        etcm=sqrt(sinv)          ! total relativistic energy
        ecmi=etcm-mp-mt       ! kinetic energy in initial channel
        Pcmi=sqrt((sinv-mp**2-mt**2)**2-4.*(mp*mt)**2)/2/sqrt(sinv)
        Kcmi=pcmi/hc
        write(*,'(3x,a,1f10.3,a)')
     & 'Relativistic invariant mass: sqrt(s)=',
     & sqrt(sinv), ' MeV'
        write(*,'(3x,"Ecmi=",1f8.3, "MeV   Kcmi=",1f8.3," fm-1" )') 
     &  Ecmi, kcmi
      endif

      ebind=energ(1,1) ! exch(1)          
!      Kcm=sqrt(2.d0*mupt*ecm)/hc        
      Ethr=Ecmi+ebind
      Kthr=sqrt(2.d0*mupt*Ethr)/hc

!      write(99,*) 'Ethr=',Ethr, ' Kthr=',Kthr
!      write(99,'(a,1f8.3,a,3f10.5)')'Elab=',elab,
!     & ' mp,mc,mv=',mp/amu,mc/amu,mv/amu
!      write(99,*)'Ecmi=',Ecmi, "Kcm=",kcmi, 'Ethr=',Ethr


c *** Assign global index to states of all j/pi sets iPS=iPS(iset,n)
      iex=0
      allocate(ips(jpsets,nsetmax))
      do iset=1,jpsets
      nex=jpiset(iset)%nex
      do n=1,nex
      iex=iex+1
      ips(iset,n)=iex
      if ((verb.ge.3).and.(iset.gt.1).and.(n.eq.1)) write(*,*)
      if (verb.ge.3) write(*,*)'Set:',iset,' state=',n,"-> iex=",iex
      enddo
      enddo


c *** Sum discrete breakup angular distributions    
      if (any(jpiset(1:jpsets)%nho>0)) then  
         allocate(dsdew_b(jpsets,nsetmax,nth))
         dsdew_b(:,:,:)=0
      endif
      
      open(90,file='dsdw_bu.xs',status='unknown')
      do ith=1,nth
      xstot=0     
      do iset=1,jpsets
!      if (jsets(iset)) write(*,*)'include j/pi set', iset
      if (.not.jsets(iset)) cycle
!      if(jpiset(iset)%nho>0) cycle ! not bins (PS)
!      if (.not.allocated(dsdew_b)) allocate(dsdew_b(jpsets,nsetmax,nth))
      nex=jpiset(iset)%nex
      jpf=jpiset(iset)%jtot
      nmf=nint(2*jpf+1)
      do n=1,nex  
      xsaux=0
      iex=ips(iset,n)
      
      if (energ(iset,n).lt.0)    cycle  ! bound state 
      if (energ(iset,n).gt.Ethr) cycle  ! closed channel
      iam=0
      raux=0
      do im=1,nmi
      do imp=1,nmf
      iam=iam+1  
      raux=10.*abs(famps0(iex,ith,iam))**2/nmi
      xsaux = xsaux  + raux 
      xstot = xstot  + raux
      enddo !imp 
      enddo !im
      
      if(jpiset(iset)%nho.eq.0)then
          wbin=binset(iset)%wbin(n)
          dsdew_b(iset,n,ith)=  dsdew_b(iset,n,ith) + xsaux/wbin
      endif
      enddo ! states within j/pi set
      enddo !j/pi sets
      write(90,*) thmin + dth*(ith-1),xstot
      enddo !ith
      call flush(90)
      close(90)
c ----------------------------------------------------------------------


c Print ds/dE for bins !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (any(jpiset(1:jpsets)%nho>0)) then
      open(92,file='dsde_bin.xs')
      write(*,*)' ds/dE for bins sets: thcut=',thcut
      do iset=1,jpsets
      if(jpiset(iset)%nho>0) cycle ! not bins (PS)
      nex=jpiset(iset)%nex
      do n=1,nex    
      iex=ips(iset,n)
      exc  = energ(iset,n)
      xsaux=0
	  do ith=1,nth
      raux=0.
      th = thmin + dth*(ith-1)
      if (th.gt.thcut) cycle      
      xsaux=xsaux + dsdew_b(iset,n,ith)*2*pi*sin(th*pi/180.)*dth*pi/180.   
      enddo !ith
      write(92,*) exc,xsaux,iex
      enddo !n within set
      write(92,*)'&'
	  enddo ! iset
      endif ! are there bin sets?



      
c *** Overlaps between PS's and Scattering states  -----------------------------
      allocate(gsolap(jpsets,nsetmax,maxchan,ncont))
      gsolap(:,:,:,:)=0
      if (emax<0)     emax=maxval(energ(:,:))
      if (emax>Ethr)  emax=Ethr-0.001
      if (emin.le.0.) emin=0.01
      dec=(emax-emin)/dble(ncont-1)
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      if (jpiset(iset)%nho.eq.0) cycle
      nchan=jpiset(iset)%nchan
      nex =jpiset(iset)%nex
      do inc=1,nchan 
      if (inc.gt.maxchan) stop'inc too large in d2sigma!'
      jci    =jpiset(iset)%jc(inc)
      li     =jpiset(iset)%lsp(inc)
      excore =jpiset(iset)%exc(inc) ! core energy
      ic     =jpiset(iset)%cindex(inc)
      if (ic.ne.icore) then
         if (ic.eq.abs(icore)) then
         call solapread(gsolap(iset,:,inc,:),jpiset(iset)%partot,
     & inc,nsetmax,ncont,jpiset(iset)%jtot,dec,emin+excore)
         else if (verb.ge.2) then
         write(*,'(4x,"-> skipping core state",i3)')ic 
         endif
      cycle
      endif
      
!          write(0,*)'method=',method
      method=4;
      nbas  =50
      ns    =1; 

      if (allocated(wfcont)) deallocate(wfcont)
      if (ncont.gt.0) then
        allocate(wfcont(ncont,nchan,nr)); wfcont=0.
        call wfrange(iset,nchan,inc,emin,emax,ncont,energy,wfcont,
     & smate,delta,writewf)
      else
       write(*,*)'*ERROR* ncont,nchan=',ncont,nchan
      endif

!!!! CHECK Continuum WFS
!      write(505,*)'#WF for Ecm;J=',emin,jpiset(iset)%jtot
!      do ir=1,nr
!         write(505,111)rvec(ir),wfcont(1,1,ir)
!      enddo
!!!!!!!!!!!!!!!!!!!!!!!!

      do n=1,nex
      if (energ(iset,n).lt.0)    cycle ! bound state
      if (energ(iset,n).gt.Ethr) cycle ! closed channel
         
      erel(:)=0
      do iecv=1,ncont
      ecv=emin+(iecv-1)*dec
      erel(iecv)=ecv
      if (ecv.lt.eps) ecv=eps
      kcv=sqrt(2.d0*mucv*ecv)/hc    
      sumr=0
      caux=0.
      do ich=1,nchan ! sum in final channels
      jcf =jpiset(iset)%jc(ich)
      lf  =jpiset(iset)%lsp(ich)
      jf  =jpiset(iset)%jc(ich)
      do ir=1,nr   
      gaux(ir)=wfc(iset,n,ich,ir)*wfcont(iecv,ich,ir)*rvec(ir)*
     &         (-1)**(li+jci+lf+jcf) ! AMM: I think this phase is always 1! 
      enddo !ir

      call simc(gaux,resc,1,nr,dr,nr)
      caux=caux+resc
      sumr = sumr+abs(caux)**2
      enddo ! ich 
      gsolap(iset,n,inc,iecv)=caux
      enddo ! iecv (c-v relative energy)
      enddo ! n  (PS within this j/pi set)
      enddo ! inc  (incoming channel)
      enddo ! iset (j/pi set)
      icore=abs(icore)
     
c *** -------------- PRINT OVERLAPS FOR TESTING ----------------------------------
      if (verb.ge.3) then
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nex  =jpiset(iset)%nex
      nchan=jpiset(iset)%nchan
      if (jpiset(iset)%nho.eq.0) cycle
      do n=1,nex
      raux=0
      write(97,'(a,i3,a,i3,a,1f10.4)')'# Set:',iset,
     & ' n=',n,' Ex=',energ(iset,n)
      do iecv=1,ncont
       ecv=emin+(iecv-1)*dec
       if (energ(iset,n).lt.0)    cycle ! bound state
       if (energ(iset,n).gt.Ethr) cycle ! closed channel
       if (ecv.lt.eps) ecv=eps
       kcv=sqrt(2.d0*mucv*ecv)/hc
       jac=mucv/(hc**2)/kcv
       write(97,111)ecv,
     &  (jac*abs(gsolap(iset,n,inc,iecv))**2,inc=1,nchan)
111    format(2x,1f12.6,2x,10g14.6)
       do inc=1,nchan
       raux=raux+jac*abs(gsolap(iset,n,inc,iecv))**2*dec
       enddo ! inc
      enddo !iecv
      if (verb.ge.3) then
      write(97,'(3x,a,i3,a,i3,a,1f8.5)')'# -> Set:',iset, '  PS=',n,
     & ' Norm solap=',raux*2/pi
      endif
      write(97,*)'&'
      enddo ! n
      enddo !iset
      endif ! verb
c-------------------------------------------------------------------------------


c *** Compute DOUBLE x-sections dsigma/(dEdW) FOR SELECTED CORE STATE ICORE
      if (ncont.gt.0) then
        allocate(dsdew(ncont,nth))
      else  
         write(*,*)'*ERROR* ncont,nchan=',ncont,nchan
      endif
        
      dsdew=0
      raux=0.5/(pi**3)/hc**6*mupt**2*mucv/nmi
      raux=raux*10 ! fm^2 -> mb

      excore=qnc(icore)%exc
      write(*,'(3x,a,i3,a,1f7.3,a)')
     & '=> Double diff. x-sections for final core state',
     &  icore,' with energy=',excore, ' MeV'

      do iecv=1,ncont
      ecv=emin+(iecv-1)*dec
      if (ecv.lt.1e-3) ecv=1e-3
      kcv=sqrt(2.d0*mucv*ecv)/hc  ! above threshold!!!
      do iset=1,jpsets
      if (jpiset(iset)%nho.eq.0) cycle ! bins
!      if (jsets(iset)) write(*,*) 'include j/pi set',iset
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      nex  =jpiset(iset)%nex
      jpf  =jpiset(iset)%jtot
      nmf  =nint(2.*jpf+1.)
      iam=0
      do im =1,nmi
      do imp=1,nmf
      iam=iam+1
      nchan=jpiset(iset)%nchan
      do inc=1,nchan  ! sum over open channels
      ic  =jpiset(iset)%cindex(inc)
      if (ic.ne.icore) cycle
      li     =jpiset(iset)%lsp(inc)
      ji     =jpiset(iset)%jsp(inc)
      jci    =jpiset(iset)%jc(inc)
      excore =jpiset(iset)%exc(inc) ! core energy 
      ic     =jpiset(iset)%cindex(inc)

      ecmf   =ecmi-dabs(ebind)-ecv-excore ! p-t c.m. energy in final channel !! CHEEEEEEEEEEEEEEEEECK
      if (ecmf.lt.0) then
          write(*,*)'Ecm=',Ecmf, 'Ecv=',Ecv
          stop
      endif
      if (ecmf.lt.1e-4) ecmf=1e-4

      if (.not.rel) then ! non-relativistic
        Kcmf=sqrt(2.d0*mupt*ecmf)/hc   
        if ((ith.eq.1).and.(iecv.eq.1))
     &     write(*,*)'Kcmf(nr)=',kcmf    
      else               ! relativistic
        mpx=mp+abs(ebind)+ecv+excore
        Pcmf=sqrt((sinv-mpx**2-mt**2)**2-4.*(mpx*mt)**2)/2/sqrt(sinv)
        Kcmf=pcmf/hc
!        if ((ith.eq.1).and.(iecv.eq.1))
!         write(*,*)'Kcmf(rel)=',kcmf,sqrt(2.d0*mupt*ecmf)/hc
      endif

     

      Kcmf=sqrt(2.d0*mupt*ecmf)/hc 
      facK=kcmf/kcmi/kcv
! commented v26
!      f2t=-2.*pi*hc**2/mupt*sqrt(Kcmi/Kcmf) ! conversion factor f(theta) -> T-mat
!
      do ith=1,nth  
      sumps=0d0
      do n=1,nex ! PS states within this j/pi set
      if (energ(iset,n).lt.0)    cycle ! bound states
      if (energ(iset,n).gt.Ethr) cycle ! closed channels
      iex=ips(iset,n)
!      if((iecv.eq.1).and.(imp.eq.1)) 
!     &   write(*,*)'inc,iset,n,iex=',inc,iset,n,iex
c v26 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ecmf=ecmi-dabs(ebind)-energ(iset,n)-excore
      Kcmf=sqrt(2.d0*mupt*ecmf)/hc   
      f2t=-2.*pi*hc**2/mupt*sqrt(Kcmi/Kcmf) ! conversion factor f(theta) -> T-mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      tmat_aux=f2t*famps0(iex,ith,iam)
      sumps=sumps+gsolap(iset,n,inc,iecv)*tmat_aux
!      if (ith.eq.1) write(*,*) iset,n,inc,iecv,
!     & abs(gsolap(iset,n,inc,iecv))
      enddo ! PS's
      dsdew(iecv,ith)= dsdew(iecv,ith)
     &               + raux*facK*abs(sumps)**2
      enddo ! ith
      enddo ! inc
      enddo ! imf
      enddo ! imi
      enddo ! j/pi 
      enddo ! iecv (valence-core relative energy)  


c Integrate in ENERGY, to get dsigma/dOmega for core state ICORE
      open(90,file='dsdw_conv.xs',status='unknown')
      open(91,file='dsdwe_conv.xs',status='unknown')
      open(92,file='dsdwe_conv.gnu',status='unknown')
      write(91,'("nel ", i5,2x, "ang ", i5)') ncont,nth 
      do ith=1,nth
      raux=0.
      th = thmin + dth*(ith-1)
      write(91,'(a,1f8.4)') '#thcm=',th
      do iecv=1,ncont
      ecv=emin+(iecv-1)*dec
      raux=raux + dsdew(iecv,ith)*dec
      write(91,'(1f10.4,1g14.6)') ecv, dsdew(iecv,ith)      
      write(92,'(2f10.4,1g14.6)') th,ecv, dsdew(iecv,ith)      
      enddo !iecv
      write(92,*)' ' 
      write(90,*) th, raux
      enddo !ith
      call flush(90)
      close(90); close(91); close(92)



c Integrate in ANGLE, to get dsigma/dE  for core state ICORE
      open(93,file='dsde_conv.xs',status='unknown')
      xstot=0
      do iecv=1,ncont
      ecv=emin+(iecv-1)*dec
      raux=0.
      do ith=1,nth
      th = thmin + dth*(ith-1)
      if (th.gt.thcut) cycle
      raux=raux + 2*pi*dsdew(iecv,ith)*
     &     sin(th*pi/180.)*dth*pi/180.   
      enddo !iecv
      write(93,*) ecv +excore, raux
      xstot=xstot + raux*dec
      enddo !ith
      write(*,*)' Integrated dsigma/dEx=',xstot
      call flush(93)
      close(93)
      
      

c---------------------------------------------------------------------------------
c     Start calculation of triple differential cross sections 
c     (inspired in cdcn code by JAT)
c---------------------------------------------------------------------------------
      if (.not.triplexs) return
      write(*,'(//,3x, "** TRIPLE DIFFERENTIAL CROSS SECTIONS **",/ )')

!      open(99,file='thox.log')
      
      read(kin,nml=framework)
      write(77,'(a)') sys
      write(77,316) mt/amu,zt,elab/(ac+av)
      write(77,316) ac,zc
      write(77,316) av,abs(ebind)
      write(77,316) elab
      write(77,*) idet
      
      if (phixs) then
      write(777,'(a)') sys
      write(777,316) mt/amu,zt,elab/(ac+av)
      write(777,316) ac,zc
      write(777,316) av,abs(ebind)
      write(777,316) elab
      write(777,*) idet
      end if      
      
*     ---------------------------------------------------------------
      call cpu_time(ti)
      xi=(0.d0,1.d0)
      acc=1.d-6
      degrad=pi/180.d0
      mtot=mc+mv+mt
!!!! REDEFINE Mp as IN TESTN 
      mp=mc+mv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      excore=qnc(icore)%exc
      write(*,'(5x,a,i3,a,1f7.3,a)')
     & '=> Triple x-sections for final core state',
     &  icore,' with energy=',excore, ' MeV'


      allocate(angsr(nth))
      nchann=maxval(jpiset(:)%nchan)
      if (allocated(tmat)) deallocate(tmat)
      allocate(tmat(jpsets,nchann,nmi*nmfmax))
      allocate(fin(jpsets,nchann,nmi*nmfmax)) ! TESTING

*     --------------------------------------------------------------- 
      do iang=1,nth
      angsr(iang)=thmin+(iang-1)*dth
      angsr(iang)=angsr(iang)*degrad
      enddo
      polar(:)=0d0
*     ---------------------------------------------------------------
      read(kin,nml=gridener)   !For alpha and beta this will be Ek (the v-C energy)
      read(kin,nml=gridthetac) !For alpha and beta this will be thcm (the vC-T theta) phicm=0
      read(kin,nml=gridthetav) !For alpha this will be thq (the v-C theta), for beta this will be the angle between kcm (vC-T) y q (v-C)
      read(kin,nml=gridphi)    !For alpha and beta this will be phiq (the v-C phi)
      read(kin,nml=polarization)
      if (sum(abs(polar)).lt.1e-2) polar(:)=1
*     ---------------------------------------------------------------
      iten=nint((Enup-Enlow)/dEn)+1
      dEn=(Enup-Enlow)/(iten-1)
*     ---------------------------------------------------------------
*     particle with specified energy defines phi=0 degrees
*     ---------------------------------------------------------------
      idet=1!MGR
      if(idet.eq.1) then
      cospc=1.d0
      sinpc=0.d0
      else
      cospv=1.d0
      sinpv=0.d0
      endif
*     ---------------------------------------------------------------           
      itc=nint((tcu-tcl)/dtc)+1
      dtc=(tcu-tcl)/(itc-1)
      dtcr=dtc*degrad
*     ---------------------------------------------------------------           
      itv=nint((tvu-tvl)/dtv)+1
      dtv=(tvu-tvl)/(itv-1)
      dtvr=dtv*degrad
*     ---------------------------------------------------------------           
      itphi=nint((phiu-phil)/dphi)+1
      itphim=itphi/2+1 ! AMM
      dphi=(phiu-phil)/(itphi-1)
      dphir=dphi*degrad
*     ---------------------------------------------------------------
      write(77,307) itc,dtcr,tcl,tcu,dtc
      write(77,307) itv,dtvr,tvl,tvu,dtv
      write(77,307) iten,dEn,Enlow,Enup
      write(77,307) itphi,dphir,phil,phiu,dphi
      
      if (phixs) then
      write(777,307) itc,dtcr,tcl,tcu,dtc
      write(777,307) itv,dtvr,tvl,tvu,dtv
      write(777,307) iten,dEn,Enlow,Enup
      write(777,307) itphim,dphir,phil,phil+(itphim-1)*dphi,dphi
      end if
      
*     ---------------------------------------------------------------
*     specify total incident momentum (momentum of c.m.) in MeV/c
      if(sys.eq.'lab') totp=sqrt(2.d0*mp*elab)
*     if detector positions refer to c.m. frame 
      if(sys.eq.'com') totp=0.d0
*     ---------------------------------------------------------------
*     set up projectile cm K vector (qp) (in z-direction)
*     set up ptot vector and ktot (also in z-direction)             
      do j=1,2      
      qp(j)=0.d0
      ptot(j)=0.d0       
      ktot(j)=0.d0
      enddo
      qp(3)=Kcmi
      ptot(3)=totp
      ktot(3)=totp/hc
      erelmax=0.d0
!      write(99,*)'mp,elab,ptot(3)=',mp/amu,elab,totp
*     -------------------------------------------------------------- 
*     compute and store CG coefficients
*     -------------------------------------------------------------- 
      sc=qnc(icore)%jc
      ind=0
      do imu=1,nint(2*sc+1)
      mu=imu-sc-1
      do isig=1,nint(2*sn+1)
      sig=isig-sn-1
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      jpf  =jpiset(iset)%jtot
      do inc=1,nchan
      li     =jpiset(iset)%lsp(inc)
      ji     =jpiset(iset)%jsp(inc)
      jci    =jpiset(iset)%jc(inc)
      ic  =jpiset(iset)%cindex(inc)
      if (ic.ne.icore) cycle
      rli=li
      if (fail3(rli,sn,ji)) cycle
      if (fail3(ji,jci,jpf)) cycle
      do inu=1,2*li+1
      rnu=inu-li-1
      if (abs(rnu+sig).gt.ji.or.abs(rnu+sig+mu).gt.jpf) cycle
      ind=ind+1
      cgc(ind)=cleb(rli,rnu,sn,sig,ji,rnu+sig)*
     .cleb(ji,rnu+sig,jci,mu,jpf,rnu+sig+mu)
      enddo
      enddo
      enddo
      enddo
      enddo
      if(ind.gt.ncleb) then
      write(*,*)'increase ncleb up to',ind
      stop 
      endif
*     --------------------------------------------------------------- 
*     compute largest orbital angular momentum (lmax), maximum relative 
*     momentum (kmax), and maximum nb of states per set (maxne)
*     ---------------------------------------------------------------
!      write(99,'(/,4(a,1f9.5,4x),/)')'Ecmi=',Ecmi," Kcmi=",Kcmi, 
!     &"Kthr=",Kthr, 'Ethr=',Ethr

      lmax=0
      kmax=0
      maxne=0
!      write(*,*)'energ(1,1)=',energ(1,1)
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      nex=jpiset(iset)%nex
      jpf=jpiset(iset)%jtot
      nmf=nint(2.*jpf+1.)
      iam=nmi*nmf 
      if (nex.gt.maxne) maxne=nex
      do n=1,nex     
      ebin=energ(iset,n) 
      if (ebin.lt.0) cycle ! bound state
      iex=ips(iset,n)
!      write(99,*)'iset,n=',iset,n
!      write(99,308) (famps0(iex,1,iii),iii=1,iam)
308   format(1p,6e12.4)
      if (jpiset(iset)%nho.eq.0) then ! bin set
!        ebin=hc**2*(binset(iset)%khat(n))**2/(2.d0*mucv)
!        energ(iset,n)=ebin
        if (binset(iset)%kmax.gt.kmax) 
     &  kmax=binset(iset)%kmax
      else                            ! PS set
        kmax=sqrt(2*mucv*emax)/hc
      endif
      ecmf=ecmi-abs(ebind)-ebin
      Kcmf=sqrt(2.d0*mupt*ecmf)/hc    
      f2t=-2.*pi*hc**2/mupt*sqrt(Kcmi/Kcmf) ! conversion factor f(theta) -> T-mat
      famps0(iex,:,:)=f2t*famps0(iex,:,:)      ! convert f(theta) -> T(theta)
      if (ecmf.lt.1e-4) ecmf=1e-4
 3130 format(a,i3,a,i3,a,1f11.6,a,2f11.5,a,1f11.5,a,1f11.5)
      enddo ! n
      do inc=1,nchan 
      if(jpiset(iset)%lsp(inc).gt.lmax) then
      lmax=jpiset(iset)%lsp(inc)
      endif
      enddo !inc
      enddo !jpsets

      write(*,'(3x,a,i3,a)') 'There are',jpsets,' sets with:'
      write(*,'(5x,a,1f7.3, a,1f8.3,a)') 'o Maximum relative momentum=',
     & kmax, ' fm-1 and relative energy:',hc**2*kmax**2/2/mucv,' MeV'
      write(*,'(5x,a,i3)') 'o Maximum orbital ang. momentum=',lmax
      write(*,'(5x,a,i3)') 'o Max. number of states=',maxne

      if (allocated(xyt)) deallocate(xyt)
      allocate(xyt(2,maxne))
      if (allocated(fxyc)) deallocate(fxyc)
      allocate(fxyc(10,maxne))

*     ---------------------------------------------------------------
*     calculate and store coefficients for spherical harmonics
*     ---------------------------------------------------------------
      do ilval=0,lmax
      do inu=1,2*ilval+1
      cplm(ilval,inu)=ylmc2(ilval,inu-ilval-1)
      enddo
      enddo
*     ---------------------------------------------------------------
*     loop over core particle detection angle
      wrt=.false. ; icount=0; if2c=0 !MGR
      
      write(0,*) '- Triple diff xsections need ',
     &           itc*itv*iten*lr8/1e6,'MB' !MGR  
      allocate(xs3body(iten,itv,itc),energ3(iten))
!      allocate(energ4(itphim*itc*itv*iten))

      xs3body=0d0
!AMM 
      if (phixs) then
       itphim=itphi/2+1
       write(0,*) '- Triple diff xsections WITHOUT phi int. need ',
     &           itc*itv*iten*itphim*lr8/1e6,'MB' !AMM
      allocate(xs3body_phi(iten,itv,itc,itphim))
      allocate(energ4(iten,itv,itc,itphim))
      xs3body_phi=0.0
      end if
      
      do ien=1,iten !MGR this energy will be the relative energy between core and valence
        energ3(ien)=Enlow+(ien-1)*dEn
      enddo
!      write(0,*) 'Ener', energ3(:)
    
!$OMP PARALLEL DO FIRSTPRIVATE(tcd,tc,costc,sintc,iv,tvd,tv,costv,
!$OMP& sintv,itphim,iflag,zert,ien,En,p1L,p2L,ip,sigphi,phid,phi,cospv,
!$OMP& sinpv,cospc,sinpc,pcL,pvL,tem,aq,bq,cq,dq,iroots,iroot,iroot1n,
!$OMP& iroot2n,Ec,Ev,kvL,kcL,bkp,kp,ekb,eks,
!$OMP& tmatsq,tmatsq2,ybo,xb,xbd,yb,phkb,xyt,iang,nial,mkp,co,phks,co2,
!$OMP& si2,si,phask1,phask,phas,ilval,j,inu,aleg,ylm,fxyc,tmat,
!$OMP& fin,iset,nchan,jpf,nex,nmf,iam,inc,fv,n,iii,ic,excore,ier,nieg,
!$OMP& niel,iecv,ecv,ecmf,kcmf,th,sumps,iex,ibin,ki,kf,dk,dkb,ik,phbin,
!$OMP& if2c,icount,imo,mo,ind,imu,mu,isig,sig,ampt,ampt2,li,ji,jci,rli,
!$OMP& rnu,mult,xsig,nk,kinset,wkfac,haux,fxytab) 
!$OMP& PRIVATE(ii,maxtab,imax,jmax) 

!      write(*,*)'itc,itv,iten=',itc,itv,iten
      
      do 10 ii=1,itc
      tcd=tcl+(ii-1)*dtc
      tc=tcd*degrad
      costc=cos(tc)
      sintc=sin(tc)
*     ---------------------------------------------------------------
*     loop over valence particle detection angle
      do 20 iv=1,itv
      tvd=tvl+(iv-1)*dtv
      tv=tvd*degrad
      costv=cos(tv)
      sintv=sin(tv)
*     ---------------------------------------------------------------
*     check if either theta is zero since no need to do multiple phi
*     calculations (if wanted, itphi>1) in this geometry
*     ---------------------------------------------------------------
      itphim=itphi/2+1
      iflag=0
      zert=(abs(tcd).le.acc.or.abs(tvd).le.acc)
      if(itphi.gt.1.and.zert) then
      itphim=1
      iflag=1
      endif
*     ---------------------------------------------------------------
*     loop over detected energy of chosen particle
      do 30 ien=1,iten
      En=Enlow+(ien-1)*dEn

       p1L=sqrt(2.d0*mucv*En)/hbarc

*     ---------------------------------------------------------------
*     the detected particle defines phi=0 degrees
*     ---------------------------------------------------------------
*     loop over azimuth detection angle          
      do 40 ip=1,itphim
      icount= icount+1
      sigphi(ip)=0.d0
      phid=phil+(ip-1)*dphi
      phi=phid*degrad  

      cospv=cos(phi)
      sinpv=sin(phi)

*     -----------------------------------------------------------------
*     construct remaining vectors
*     -----------------------------------------------------------------
!      do j=1,3
!MGR get wavevector modulus for p-t vector (reusing p2L) cm
      p2L=sqrt(kcmi**2+2d0*mupt/hbarc**2*(ebind-excore)-
     & mupt/mucv*p1L**2)
*     wavevector of cm of core and valence particles
      bkp(3)=p2L*costc
      bkp(1)=p2L*sintc
      bkp(2)=0d0
!*     wavevector of relative motion of core and valence particles
      !alpha
      kp(3)=p1L*costv
      kp(1)=p1L*sintv*cospv
      kp(2)=p1L*sintv*sinpv
      !beta
      !kp(3)=p2L*cos(tc+tv)
      !kp(1)=p2L*sin(tc+tv)*cospv
      !kp(2)=p2L*sintv*sinpv
      !enddo
!      eKb=hc**2*dot(bkp,bkp)/(2.d0*mupt)
!      eks=hc**2*dot(kp,kp)/(2.d0*mucv)
      eKb=kcmi**2/(2.d0*mupt)*hc**2+ebind-excore-En
      eks=En
!     write(99,'(5i3,20g14.5)') iv,ii,ien,ip,iroot,ekb,eks,bkp(1:3)
      phbin(:)=1
*     -----------------------------------------------------------------
*     increment maximum relative energy 
*     -----------------------------------------------------------------
      if(eks.gt.erelmax) erelmax=eks

!      if (wrt)
!     & write(99,'(5i3,20g12.5)') iv,ic,ien,ip,iroot,ekb,eks,bkp(1:3)
      tmatsq=0; tmatsq2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ((eks.gt.emax).or.(eks.lt.emin)) goto 500  ! NEEDS CHECKING (AMM)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*     -----------------------------------------------------------------
*     CALCULATION OF T-MATRIX <bkp,kp|T|qp> STARTS HERE
*     -----------------------------------------------------------------
*     ybo = Kout_cm, xb =theta_cm (rad), xbd (degrees)
*     -----------------------------------------------------------------
*     compute polar angles (th,ph) of cm wave vector
*     -----------------------------------------------------------------
      ybo=m(bkp)
!      xb=min(1d0,bkp(3)/ybo)
!      xb=max(-1d0,xb)
!      xb=acos(xb)
!      xbd=xb/degrad
      yb=abs(ybo-Kthr)
!      if(abs(bkp(1)).lt.1.d-6.and.abs(bkp(2)).lt.1.d-6) then
!      phKb=0.d0
!      else 
!      phKb=atan2(bkp(2),bkp(1))
!      endif
      xbd=tcd
      xb=tc
      phKb=0d0
      if(wrt) 
     & write(99,'(4i3," phikb,xbd=",20g12.5)') iv,ii,ien,ip,phKb,xbd  
*     -----------------------------------------------------------------
*     pick out nearest angles indices in array for interpolations
*     -----------------------------------------------------------------
      xyt(:,:)=0
!      iang=nint((xbd-thmin)/dth)+1
!      if (wrt) write(99,*)'thmin,thinc=',thmin,dth
      iang=int((xbd-thmin)/dth)+1 ! CHECK
      nial=max0(iang-2,1);  
!      if (wrt)write(99,*) 'iang,nial=',iang,nial
      niag=min0(nial+5,nth)
      nial=niag-5
      do iang=1,6
      xyt(1,iang)=angsr(nial+iang-1)
c v26
      xtab(iang)=angsr(nial+iang-1)
      enddo
      if (wrt)write(99,*) 'iang,nial,niag,angsr=',iang,nial,niag
      if (wrt)write(99,*)'angrs=',angsr(1),angsr(nth),angsr(nial+iang-1)
      
*     -----------------------------------------------------------------
*     compute polar angles (theta,phi) of relative wave vector
*     -----------------------------------------------------------------
      mkp=m(kp)
      co=kp(3)/mkp
      if(abs(kp(1)).lt.1.d-6.and.abs(kp(2)).lt.1.d-6) then
      phks=0.d0
      else 
      phks=atan2(kp(2),kp(1))
      endif
      co2=co*co
      si2=1.d0-co2
      si=sqrt(abs(si2))
      if (wrt)write(99,'(4i3," phks=",20g12.5)') iv,ii,ien,ip,phKs
*     --------------------------------------------------------------- 
*     calculate the spherical harmonics for this relative wave vector
*     ---------------------------------------------------------------
      phask1=(1.d0,0.d0)
      phask=exp(xi*phks)
      do ilval=0,lmax
      phask1=phask1/phask
      phas=phask1
      do inu=1,2*ilval+1
      phas=phas*phask
c AMM: Use explicit forms for PLM (l<9); otherwise, original recursive formula
       if (ilval.le.10) then
       aleg=plm_nr(ilval,inu-ilval-1,co,co2,si,si2)
       else
       aleg=plmrec(ilval,inu-ilval-1,co,co2,si,si2)
       endif
!      ylm(ilval,inu)=ylmc2(ilval,inu-ilval-1)*phas*aleg
      ylm(ilval,inu)=cplm(ilval,inu)*phas*aleg
!     write(99,'(6i3,20g12.5)') iv,ii,ien,ip,ilval,inu,
!     & ylmc2(ilval,inu-ilval-1),cplm(ilval,inu)
      enddo
      enddo
      fxyc(:,:)=0; fxytab=0
      tmat(:,:,:)=0; fin=0;
      kinset(:)=.false.
      wkfac(:)=0
*     -------------------------------------------- 
*     interpolate the sum \sum_i gsolap(i)*tmat(i)
*     --------------------------------------------
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      jpf=jpiset(iset)%jtot
      nex=jpiset(iset)%nex
      nmf=nint(2.*jpf+1.)
      iam=nmi*nmf 
      if (iam.lt.0) stop 'iam<0! '
      do inc=1,nchan
      ic=jpiset(iset)%cindex(inc)
      if (ic.ne.icore) cycle

      IF ((jpiset(iset)%nho).gt.0) then ! PS's -> need overlaps
        if (eks.gt.emax) then
        fv(1:nex)=0
        else
        do n=1,nex
        haux(1:ncont)=gsolap(iset,n,inc,1:ncont)
        fv(n)=ffc4((eks-emin)/dec,haux,ncont)
        enddo
        endif
      ELSE                              ! BINs
        if (inc.ne.jpiset(iset)%inc) cycle  ! CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ik=0
! determine and store bin w(k) factor for this j/pi set
      do n=1,nex
      ki=binset(iset)%klow(n)
      kf=binset(iset)%kup(n)
      dk=kf-ki
      nk=binset(iset)%nk(n)
      dkb=(kf-ki)/(nk - 1)

      if ((mkp.ge.ki).and.(mkp.le.kf)) then
        ik=nint((mkp-ki)/dkb+1)     
        if (ik.gt.nk) stop 'ik too large'
        if (ik.le.0) stop 'ik <=0!!'
        kinset(iset)=.true.
!        phbin(iset)=exp(xi*binset(iset)%pshel(n,ik))
        wkfac(iset)=binset(iset)%wk(n,ik)
       if (wrt)write(99,'(a,1f6.2,i3,4f12.6)')'eks,iset,phbin=',
     &   eks,iset,phbin(iset)
       if (wrt)write(99,*)'iset,ik,wkfac=',mkp,iset,ik,wkfac(iset)
      endif !mkp within this bin
      enddo !n 

      endif  ! PS or Bins?

c Interpolate m-dependent amplitudes
      do iii=1,iam

c *** PS's .......................................................
      if ((jpiset(iset)%nho).gt.0) then 
      ier=nint((eks-emin)/dec)+1
      niel=max0(ier-2,1)
      nieg=min0(niel+5,ncont)
      niel=nieg-5
      if (wrt) write(99,*)'eks=',eks, 'ier=',ier,'niel=',niel

      maxtab=-10d0
      imax=0
      jmax=0
      do ier=1,6
      iecv=niel+ier-1 
      ecv=emin+(iecv-1)*dec 
      xyt(2,ier)=emin+(iecv-1)*dec 
c v26
      ytab(ier)=emin+(iecv-1)*dec
!      if (wrt) write(99,*)'ier,ecv=',ier,ecv
       ecmf=ecmi-abs(ebind)-ecv-excore
!      if (excore.gt.0) stop' excore!'
!      if (ecmf.lt.0) then
!          write(*,*) 'Ecm=',Ecmf,'Ecv=',Ecv; stop
!      endif
!      if (ecmf.lt.1e-4) ecmf=1e-4
       Kcmf_cont=sqrt(2.d0*mupt*ecmf)/hc     
c1      f2t=-2.*pi*hc**2/mupt*sqrt(Kcmi/Kcmf) ! conversion factor f(theta) -> T-mat
c1 ------------------------------------
      do ith=0,5
      sumps=0.d0
      do n=1,nex
      if (energ(iset,n).lt.0)    cycle ! bound states
      if (energ(iset,n).gt.Ethr) cycle ! closed channels
      iex=ips(iset,n)
! famps0() already transformed to T-matrix
!      sumps=sumps+fv(n)*f2t*famps0(iex,nial+ith,iii)

!       ecmf=ecmi-abs(ebind)-energ(iset,n)
!       Kcmf=sqrt(2.d0*mupt*ecmf)/hc  
!       raux=sqrt(kcmf/kcmf_cont)
!       sumps=sumps+fv(n)*raux*famps0(iex,nial+ith,iii)
      sumps=sumps+fv(n)*famps0(iex,nial+ith,iii)


!!!!

      enddo !  n (PS's)
      if (abs(sumps).gt. maxtab) then
        maxtab=abs(sumps)
        imax=ith+1
        jmax=ier
      endif
      fxyc(ith+1,ier)=sumps
c v26
      fxytab(ith+1,ier)=sumps
!      if (wrt) write(99,*)'ith,ier,fxy=',ith,ier,sumps
      enddo !  ith (theta)
      enddo ! ier
c v26
!      tmat(iset,inc,iii)=f2c(xb,eks,xyt,fxytab,6,6,nord,6,6)
       tmat(iset,inc,iii)=f2c(xb,eks,xyt,fxyc,6,6,nord,10,maxne)
!      tmat(iset,inc,iii)=fint2dd(xtab,ytab,fxytab,xb,eks,6,6,6)
!      tmat(iset,inc,iii)=fint2d(xtab,ytab,fxytab,xb,eks,6,6,nord,6)

!             ijcmsec=fint2db2(ecmb,thcm,xsecmb,iecm,jthcm,necm+1,
!     +       nthcm+1,1.0d0)
!       tmat(iset,inc,iii)=fint2db_jin(xtab,ytab,fxytab,xb,eks,6,6,1.0d0)
     

      if (wrt) then
          write(99,*)'iii=',iii,'xb,yb=',xb,eks
          write(99,'("xyt(1,:)=",20g12.5)') (xyt(1,ith),ith=1,6)
          write(99,'("xyt(2,:)=",20g12.5)') (xyt(2,ier),ier=1,6)
          write(99,'("fxyc=",20g12.5)') (fxyc(1:,ier),ier=1,6)
          write(99,*)' -> tmat=',tmat(iset,inc,iii)
      endif
      
      if (abs(tmat(iset,inc,iii)).gt.maxtab) then
       tmat(iset,inc,iii)=fxyc(imax,jmax)
      endif

      ELSE 
c *** BINs  .......................................................................
      ibin=0
      if (kinset(iset)) then 
!      if (abs(wkfac(iset)).lt.1e-6) 
!     &    write(99,*)'iset,wkfac=',iset,wkfac(iset)
!       if (eks < binset(iset)%emin-5) cycle
!      if (eks > binset(iset)%emax+5) cycle
!      if (mkp < binset(iset)%klow(1))   cycle
!      if (mkp > binset(iset)%kup(nex)) cycle

      do n=1,nex
!      ki=binset(iset)%klow(n)
!      kf=binset(iset)%kup(n)
!      dk=kf-ki
!      nk=binset(iset)%nk(n)
      ecmf=ecmi-abs(ebind)-energ(iset,n)
      if(wrt) write(99,'(a,i3,6f12.6)')'n,k,bindk=',
     & n,binset(iset)%khat(n),wkfac(iset) !phbin(iset)/sqrt(dk)
      if (ecmf.lt.0) then
          write(*,*) 'Ecm=',Ecmf,'Ecv=',Ecv; stop
      endif
      if (ecmf.lt.1e-4) ecmf=1e-4
      Kcmf=sqrt(2.d0*mupt*ecmf)/hc     
c1      f2t=-2.*pi*hc**2/mupt*sqrt(Kcmi/Kcmf) ! conversion factor f(theta) -> T-mat
      xyt(2,n)=abs(Kcmf-Kthr)
!      call flush(99)
      do ith=0,5
        iex=ips(iset,n)
        fxyc(ith+1,n)=famps0(iex,nial+ith,iii)
     &               *sqrt(pi/2.)*wkfac(iset) ! includes phase e^(-i delta_el)
!     &             *sqrt(pi/2.)/sqrt(dk)  * phbin(iset)
      enddo !ith

      if(wrt)write(99,'(a,4i4,14g12.6)') 'iset,n,nial,ik,',
     &    iset,n,nial,ik,fxyc(1:6,n)

      if(wrt)write(99,'(a,i3,20e12.4)') '3 famp,iii=',iii,
     &  (famps0(iex,nial+ith,iii),ith=0,5)


      enddo !n 
      if2c=if2c+1
      tmat(iset,inc,iii)=f2c(xb,yb,xyt,fxyc,6,nex,nord,10,maxne) ! CHECK DIMS HERE!!!!!!!!!!!!!!!!
c1      fin(iset,inc,iii)=f2c(xb,yb,xyt,fxyc2,6,nex,nord,10,50) ! TESTING
      if (wrt) then
          write(99,*)'iii=',iii,'xb,yb=',xb,yb
          write(99,'("xyt(1,:)=",20g12.5)') (xyt(1,ith),ith=1,6)
          write(99,'("xyt(2,:)=",20g12.5)') (xyt(2,ier),ier=1,nex)
          write(99,'("fxyc(1,:)=",12g12.5)') fxyc(1,1:nex)
          write(99,'("fxyc(2,:)=",12g12.5)') fxyc(2,1:nex)
          write(99,*)
          endif
!          if(abs(tmat(iset,inc,iii)).gt.1e5) then
!           write(*,*)'Big tmat for iset,inc,iii=',iset,inc,iii,
!     &     tmat(iset,inc,iii)
!          endif
      else 
         tmat(iset,inc,iii)=0
      endif ! k in this set?
         
      endif ! Bins or PS's     
c    ...............................................................................
      enddo ! iii
      enddo ! inc channel
      enddo ! jpset
*     ----------------------------------------------------------------- 
*     |t-matrix|**2 summed on spin projections / (2sp+1)
*     ----------------------------------------------------------------- 
      tmatsq=0.d0
      tmatsq2=0. !!TESTING
      do imo=1,nmi
      mo=imo-jpiset(1)%jtot-1
      ind=0
      do imu=1,nint(2*sc+1)
      mu=imu-sc-1
      do isig=1,nint(2*sn+1)
      sig=isig-sn-1
      ampt=(0.d0,0.d0)
      ampt2 =0. !!!! TESTING
      do iset=1,jpsets
      if (.not.jsets(iset)) cycle
      nchan=jpiset(iset)%nchan
      nex  =jpiset(iset)%nex
      jpf  =jpiset(iset)%jtot
      do inc=1,nchan  ! sum over open channels
      ic  =jpiset(iset)%cindex(inc)
      if (ic.ne.icore) cycle
      li     =jpiset(iset)%lsp(inc)
      ji     =jpiset(iset)%jsp(inc)
      jci    =jpiset(iset)%jc(inc)
      rli=li
      do inu=1,2*li+1
      rnu=inu-li-1
      if (fail3(rli,sn,ji)) cycle
      if (fail3(ji,jci,jpf)) cycle
      if (abs(rnu+sig).gt.ji.or.abs(rnu+sig+mu).gt.jpf) cycle
      ind=ind+1
!      if (abs(cgc(ind).lt.eps)) cycle
      iam=nint((imo-1)*(2*jpf+1)+rnu+sig+mu+jpf+1)
      if (iam.lt.0) then
          write(*,*)'solvecc error: iam<0!';
          write(*,*)'iam,imo,jpf,rnu,sig,mu=',iam,imo,jpf,rnu,sig,mu
          write(*,*)'shape(tmat)=',shape(tmat)
          stop
      endif

!!! CHECK
!      raux=cleb(rli,rnu,sn,sig,ji,rnu+sig)*
!     .cleb(ji,rnu+sig,jci,mu,jpf,rnu+sig+mu)
!      if (abs(raux-cgc(ind)).gt.1e-3) then
!      write(99,*)'cgc mismatch!: rli,rnu,sig,mu,jpf=',rli,rnu,sig,mu,jpf     
!      endif
!AMMoro Sept 16
!      ampt=ampt+(-xi)**li*cgc(ind)*ylm(li,inu)*  ! NOTE phbin() inserted here!!!!!
      ampt=ampt+(-xi)**li*cgc(ind)*ylm(li,inu)* ! phbin(iset)*  ! NOTE phbin() inserted here!!!!!
     .tmat(iset,inc,iam)*exp(xi*nint(mo-rnu-sig-mu)*phKb)

!TESTING
c1      ampt2=ampt2+(-xi)**li*cgc(ind)*ylm(li,inu)*phbin(iset)*
c1     & fin(iset,inc,iam)*exp(xi*nint(mo-rnu-sig-mu)*phKb)
!!
      enddo ! inu
      enddo ! inc
      enddo ! j/pi
      tmatsq=tmatsq+abs(ampt)**2*sqrt(polar(imo))
c2      tmatsq2=tmatsq2+abs(ampt2)**2*polar(imo)   !!! TESTING
      enddo ! isigma
      enddo ! imu
      enddo ! imo
c2      tmatsq2=tmatsq2*((2.d0*pi)**3)/dble(nmi)
      tmatsq=(4.d0*pi/mkp)**2*tmatsq/dble(nmi)
c1      if (icount.lt.100)
c1     & write(*,*)'tmatsq=',tmatsq,tmatsq2,tmatsq2/tmatsq
      

500   continue
*     -----------------------------------------------------------------
*     Phase Space Factor
*     -----------------------------------------------------------------

      mult=mupt**2*mucv*m(kp)*m(bkp)/m(qp)/(2d0*pi*hc)**5/hc
      sigphi(ip)=sigphi(ip)+10.d0*mult*tmatsq
52    continue !sigphi(ip)=0.d0

*     -----------------------------------------------------------------
*     close the phi loop
*     -----------------------------------------------------------------
40    continue
*     -----------------------------------------------------------------
*     integrate over phi angles NOW if more than one phi angle
*     -----------------------------------------------------------------
      if(itphi.eq.1) then
      xsig=sigphi(1)
      else if(iflag.eq.0) then
      call sim(sigphi,xsig,1,itphim,dphir,itphim)
      xsig=2.d0*xsig
      else
      xsig=sigphi(1)*(phiu-phil)*degrad
      endif
      if(ien.eq.iten) then
!      print 315,ien,iv,ii,xsig,En,iflag
      endif
*     -----------------------------------------------------------------
       xs3body(ien,iv,ii)=xsig

! AMM July 2019
       if (phixs) then
       do ip=1,itphim
       xs3body_phi(ien,iv,ii,ip)=sigphi(ip)
       enddo !ip
       endif
       
!      write(77,315) ien,iv,ii,xsig,En
!      STOP
!      call flush(77)
*     -----------------------------------------------------------------
*     close the angle (core and valence thetas) and energy loops
*     -----------------------------------------------------------------
!      energ4(ien)=ev
!      write(0,*)ec,energ3(ien),ev
30    continue      
20    continue
10    continue
!$OMP END PARALLEL DO
      print 510
      print*,'  erelmax =',erelmax
      write(*,*)'icount=',icount,' if2c=',if2c
      write(*,*)'Dims=',itc,itv,iten,itphim
      print 510
!MGR--------------------------------------------------------------------
      icount=0
      do ii=1,itc
        do iv=1,itv
          do ien=1,iten
            write(77,315) ien,iv,ii,xs3body(ien,iv,ii),energ3(ien)
            do ip=1,itphim
            icount=icount+1
            if (phixs) write(777,318) ien,iv,ii,ip,
     &       xs3body_phi(ien,iv,ii,ip),energ3(ien),energ4(ien,iv,ii,ip)
            enddo !ip
          enddo
        enddo
      enddo    
      deallocate(xs3body,energ3)
      
      if (phixs) deallocate(xs3body_phi,energ4)
      write(*,*)'xs3body_phi has',icount, 'elements'
!-----------------------------------------------------------------------      
      call cpu_time(tf)
      t3d=tf-ti
      return 
307   format(i5,f15.8,3f10.3)
315   format(3i5,d19.8,1x,f12.6,i5)
316   format(3f10.3)
318   format(4i5,d19.8,1x,2f12.6,i5)
510   format('  ------------------------------------------------------',
     +'-----')
900   write(*,*)'Error reading f-amps file!'; stop
      end subroutine


! From xsections.f90
      subroutine budist(enex,nset,sigset)
      implicit real*8(a-h,o-z)
      parameter (nmax=500)
      integer nb
      real*8  enex(1:nset),sigset(1:nset)
      real*8 d1,di,dnex,be
      be=0.0
!c     -----------------------------------------------
      nb=0
      do ia=1,nset
!      write(*,*)'ia=',ia
      if(enex(ia).lt.0) nb=nb+1 ! bound states
      enddo

!      print*,'budist: nb=',nb

      if (nset-nb.lt.3) return
      open(15,file='dsde_ps.xs') 
!      d1=2.d0*ENEX(2) - 1.5d0*ENEX(1) - 0.5d0*ENEX(3)
      d1=enex(nb+1)
!      print*,'d1=',d1
      write(15,50) ENEX(nb+1)-be,sigset(nb+1)/d1,enex(1)+d1/2.d0 
      write(18,50) ENEX(nb+1)-be,sigset(nb+1),                          &
     &     enex(nb+1)+d1/2.d0 

!c     Interpolation for 1<=IA<NEX-1
      do IA=nb+2,nset-1
         if (enex(ia).lt.enex(ia-1)) write(15,*)'&'
         di=(ENEX(IA+1)-ENEX(IA-1))/2.d0
         write(15,50) ENEX(IA)-be,sigset(IA)/di,                        &
     &           enex(ia)+di/2.d0,di
         write(18,50) ENEX(IA)-be,sigset(IA),                           &
     &           enex(ia)+di/2.d0
      enddo
!c     Interpolation for IA=NSET
      dnex=-2.d0*ENEX(NSET-1)+1.5d0*ENEX(NSET)+0.5d0*ENEX(NSET-2)
      write(15,50) ENEX(NSET)-be,sigset(NSET)/dnex,                     &
     &        enex(nset)+dnex/2.d0,dnex
      write(18,50) ENEX(NSET)-be,sigset(NSET),                          &
     &        enex(nset)+dnex/2.d0

 50   format(6f12.4)
      sig=0
      do i=1,nset
         sig=sig+sigset(i)
      enddo
      write(*,*) ' => Total inel+bu x-section=',sig
      write(*,*)  '------------------------------------'
      write(15,*) '&'
      write(18,*) '&'
      close(15)
      return
      end


! From rmat_solvers.F90
subroutine solve_rmatrix(cmat, B_vector, nch, nlag, normfac, Rmat, solver_type, K_pure)
  implicit none
  integer, intent(in) :: nch, nlag
  real(dp), intent(in) :: normfac
  complex(dp), intent(in) :: cmat(nch*nlag, nch*nlag), B_vector(nlag)
  complex(dp), intent(out) :: Rmat(nch, nch)
  integer, intent(in), optional :: solver_type
  real(dp), intent(in), optional :: K_pure(nlag, nlag)

  integer :: stype

  stype = 1  ! Default to dense LAPACK
  if (present(solver_type)) stype = solver_type

  select case (stype)
  case (1)
    call solve_rmatrix_dense(cmat, B_vector, nch, nlag, normfac, Rmat)
  case (2)
    call solve_rmatrix_mixed(cmat, B_vector, nch, nlag, normfac, Rmat)
  case (3)
    if (present(K_pure)) then
      call solve_rmatrix_woodbury(cmat, B_vector, nch, nlag, normfac, Rmat, K_pure)
    else
      call solve_rmatrix_woodbury(cmat, B_vector, nch, nlag, normfac, Rmat)
    end if
  case (4)
    call solve_rmatrix_gpu(cmat, B_vector, nch, nlag, normfac, Rmat)
  case (5)
    call solve_rmatrix_tf32(cmat, B_vector, nch, nlag, normfac, Rmat)
  case default
    call solve_rmatrix_dense(cmat, B_vector, nch, nlag, normfac, Rmat)
  end select

end subroutine solve_rmatrix


! From rmat_solvers.F90
subroutine rmatrix_propagation(cmat, B0, B1, nch, nlag, normfac, ins, rmax, Rmat0, Rmat)
  implicit none
  integer, intent(in) :: nch, nlag, ins
  real(dp), intent(in) :: normfac, rmax
  complex(dp), intent(in) :: cmat(nch*nlag, nch*nlag), B0(nlag), B1(nlag)
  complex(dp), intent(in) :: Rmat0(nch, nch)
  complex(dp), intent(out) :: Rmat(nch, nch)

  complex(dp) :: R00(nch, nch), R01(nch, nch), R10(nch, nch), R11(nch, nch)
  complex(dp) :: inv_matrix(nch, nch)
  complex(dp), allocatable :: A_work(:,:), X_vector(:,:)
  integer, allocatable :: IPIV(:)
  integer :: ich, ichp, ir, ntotal, INFO

  ntotal = nch * nlag

  allocate(A_work(ntotal, ntotal))
  allocate(X_vector(ntotal, 2*nch))
  allocate(IPIV(ntotal))

  A_work = cmat

  ! Setup RHS for both B0 and B1
  X_vector = (0.0_dp, 0.0_dp)
  do ich = 1, nch
    do ir = 1, nlag
      X_vector((ich-1)*nlag + ir, ich) = B0(ir)
      X_vector((ich-1)*nlag + ir, ich + nch) = B1(ir)
    end do
  end do

  call ZGESV(ntotal, 2*nch, A_work, ntotal, IPIV, X_vector, ntotal, INFO)
  if (INFO /= 0) then
    write(*,*) 'ERROR in propagation ZGESV, INFO =', INFO
    Rmat = (0.0_dp, 0.0_dp)
    deallocate(A_work, X_vector, IPIV)
    return
  end if

  ! Compute the four R-matrices
  R00 = (0.0_dp, 0.0_dp)
  R01 = (0.0_dp, 0.0_dp)
  R10 = (0.0_dp, 0.0_dp)
  R11 = (0.0_dp, 0.0_dp)

  do ichp = 1, nch
    do ich = 1, nch
      do ir = 1, nlag
        R00(ich, ichp) = R00(ich, ichp) + B0(ir) * X_vector(ir + (ich-1)*nlag, ichp) * normfac
        R01(ich, ichp) = R01(ich, ichp) + B0(ir) * X_vector(ir + (ich-1)*nlag, ichp + nch) * normfac
        R10(ich, ichp) = R10(ich, ichp) + B1(ir) * X_vector(ir + (ich-1)*nlag, ichp) * normfac
        R11(ich, ichp) = R11(ich, ichp) + B1(ir) * X_vector(ir + (ich-1)*nlag, ichp + nch) * normfac
      end do
    end do
  end do

  ! Propagate: Rmat = (R11 - R10 * inv(R00 + (ins-1)*rmax*Rmat0) * R01) / (ins*rmax)
  inv_matrix = R00 + (ins - 1.0_dp) * rmax * Rmat0

  deallocate(IPIV)
  allocate(IPIV(nch))
  call ZGESV(nch, nch, inv_matrix, nch, IPIV, R01, nch, INFO)
  if (INFO /= 0) then
    write(*,*) 'ERROR in propagation inversion, INFO =', INFO
    Rmat = (0.0_dp, 0.0_dp)
    deallocate(A_work, X_vector, IPIV)
    return
  end if

  Rmat = (R11 - matmul(R10, R01)) / (ins * rmax)

  deallocate(A_work, X_vector, IPIV)

end subroutine rmatrix_propagation


! From rmat_solvers.F90
subroutine compute_smatrix(nch, Zmat_I, Zmat_O, Smat)
  implicit none
  integer, intent(in) :: nch
  complex(dp), intent(in) :: Zmat_I(nch, nch), Zmat_O(nch, nch)
  complex(dp), intent(out) :: Smat(nch, nch)

  complex(dp) :: Z_work(nch, nch)
  integer :: IPIV(nch), INFO

  Z_work = Zmat_O
  Smat = Zmat_I

  call ZGESV(nch, nch, Z_work, nch, IPIV, Smat, nch, INFO)
  if (INFO /= 0) then
    write(*,*) 'ERROR in S-matrix computation, INFO =', INFO
    Smat = (0.0_dp, 0.0_dp)
  end if

end subroutine compute_smatrix


! From gpu_solver_interface.F90
  subroutine gpu_solver_finalize(ierr)
    implicit none
    integer, intent(out) :: ierr

    integer(c_int) :: c_ierr

    if (.not. gpu_initialized) then
       ierr = 0
       return
    end if

    c_ierr = gpu_solver_finalize_c()
    ierr = c_ierr
    gpu_initialized = .false.
  end subroutine gpu_solver_finalize


! From gpu_solver_interface.F90
  subroutine gpu_solve_mixed(A, B, n, nrhs, max_refine, tol, ierr)
    implicit none
    complex(dp), intent(in) :: A(n, n)
    complex(dp), intent(inout) :: B(n, nrhs)
    integer, intent(in) :: n, nrhs, max_refine
    real(dp), intent(in) :: tol
    integer, intent(out) :: ierr

    integer(c_int) :: c_n, c_nrhs, c_max_refine, c_info, c_ret
    real(c_double) :: c_tol
    integer :: init_err

    ! Auto-initialize GPU on first call
    if (.not. gpu_initialized) then
       call gpu_solver_init(0, init_err)
    end if

    ! If GPU not available, fall back to CPU
    if (.not. gpu_available) then
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
       return
    end if

    c_n = n
    c_nrhs = nrhs
    c_max_refine = max_refine
    c_tol = tol

    c_ret = gpu_solve_mixed_c(A, B, c_n, c_nrhs, c_max_refine, c_tol, c_info)
    ierr = c_info

    if (c_ret /= 0) then
       write(*,*) "GPU solve failed, falling back to CPU"
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
    end if
  end subroutine gpu_solve_mixed


! From gpu_solver_interface.F90
  subroutine gpu_solve_double(A, B, n, nrhs, ierr)
    implicit none
    complex(dp), intent(in) :: A(n, n)
    complex(dp), intent(inout) :: B(n, nrhs)
    integer, intent(in) :: n, nrhs
    integer, intent(out) :: ierr

    integer(c_int) :: c_n, c_nrhs, c_info, c_ret

    ! If GPU not available, fall back to CPU
    if (.not. gpu_available) then
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
       return
    end if

    c_n = n
    c_nrhs = nrhs

    c_ret = gpu_solve_double_c(A, B, c_n, c_nrhs, c_info)
    ierr = c_info

    if (c_ret /= 0) then
       write(*,*) "GPU solve failed, falling back to CPU"
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
    end if
  end subroutine gpu_solve_double


! From gpu_solver_interface.F90
  subroutine gpu_get_info(info_str)
    implicit none
    character(len=*), intent(out) :: info_str

    character(kind=c_char) :: c_str(256)
    integer(c_int) :: c_len
    integer :: i

    c_len = 256
    call gpu_get_info_c(c_str, c_len)

    ! Convert C string to Fortran
    info_str = ''
    do i = 1, 255
       if (c_str(i) == c_null_char) exit
       info_str(i:i) = c_str(i)
    end do
  end subroutine gpu_get_info


! From gpu_solver_interface.F90
  subroutine gpu_multi_init(ngpu, ierr)
    implicit none
    integer, intent(inout) :: ngpu
    integer, intent(out) :: ierr

    integer(c_int) :: c_ngpu, c_ret

    c_ngpu = ngpu
    c_ret = gpu_multi_init_c(c_ngpu)
    ngpu = c_ngpu
    ierr = c_ret
  end subroutine gpu_multi_init


! From gpu_solver_interface.F90
  subroutine gpu_multi_finalize(ierr)
    implicit none
    integer, intent(out) :: ierr

    integer(c_int) :: c_ret

    c_ret = gpu_multi_finalize_c()
    ierr = c_ret
  end subroutine gpu_multi_finalize


! From gpu_solver_interface.F90
  subroutine gpu_solve_multi(A, B, n, nrhs, ierr)
    implicit none
    complex(dp), intent(in) :: A(n, n)
    complex(dp), intent(inout) :: B(n, nrhs)
    integer, intent(in) :: n, nrhs
    integer, intent(out) :: ierr

    integer(c_int) :: c_n, c_nrhs, c_info, c_ret

    c_n = n
    c_nrhs = nrhs

    c_ret = gpu_solve_multi_c(A, B, c_n, c_nrhs, c_info)
    ierr = c_info

    if (c_ret /= 0) then
       write(*,*) "Multi-GPU solve failed, falling back to CPU"
       call cpu_fallback_solve(A, B, n, nrhs, ierr)
    end if
  end subroutine gpu_solve_multi


! From utils.f90
      subroutine sbesjh(x,lmax,xj,xjp,xh1,xh1p,ifail)
c ***                                                       i.j.thompson
c ***                                                       31 may 1985.
c ***  complex spherical bessel functions from l=0 to l=lmax
c ***    for x in the upper half plane ( im(x) > -3)
c ***
c ***    xj(l)   = j/l(x)          regular solution: xj(0)=sin(x)/x
c ***    xjp(l)  = d/dx j/l(x)
c ***    xh1(l)  = h(1)/l(x)       irregular hankel function:
c ***    xh1p(l) = d/dx h(1)/l(x)            xh1(0) = j0(x) + i. y0(x)
c ***                                               =(sin(x)-i.cos(x))/x
c ***                                               = -i.exp(i.x)/x
c ***  using complex cf1, and trigonometric forms for l=0 solutions.
c ***
      implicit complex*16 (a-h,o-z)
      parameter (limit=20000)
      dimension xj(0:lmax),xjp(0:lmax),xh1(0:lmax),xh1p(0:lmax)
      real*8 zero,one,accur,tm30,absc
      data zero,one/ 0.0d0,1.0d0 /, accur /1.0d-12/, tm30 / 1d-30 /,
     #     ci / (0d0,1d0) /
      absc(w) = abs(real(w)) + abs(dimag(w))
      ifail= -1
      if(absc(x).lt.accur .or. dimag(x).lt.-3.0) go to 5
      xi = one/x
      w  = xi + xi
      pl = lmax*xi
      f = pl + xi
      b  = f + f + xi
      d  = zero
      c  = f
      do 1 l=1,limit
      d  = b - d
      c  = b - one/c
         if(absc(d).lt.tm30) d = tm30
         if(absc(c).lt.tm30) c = tm30
      d = one / d
      del= d * c
      f = f * del
      b = b + w
    1 if(absc(del-one).lt.accur) go to 2
        ifail = -2
        go to 5
c
    2 xj(lmax)   = tm30
      xjp(lmax)  = f * xj(lmax)
c
c *** downward recursion to l=0 (n.b.  coulomb functions)
c
      do 3 l = lmax-1,0,-1
      xj(l) = pl*xj(l+1) + xjp(l+1)
      xjp(l)= pl*xj(l)   - xj(l+1)
    3 pl = pl - xi
c *** calculate the l=0 bessel functions
      xj0  = xi * sin(x)
      xh1(0) = exp(ci*x) * xi * (-ci)
      xh1p(0)= xh1(0) * (ci - xi)
c
c *** rescale xj, xjp,  converting to spherical bessels.
c *** recur   xh1,xh1p             as spherical bessels.
c
        w = one/xj(0)
         pl = xi
      do 4  l = 0,lmax
      xj(l)  =  xj0*(w*xj(l))
      xjp(l) =  xj0*(w*xjp(l)) - xi*xj(l)
         if(l.eq.0) go to 4
      xh1(l) = (pl-xi) * xh1(l-1) - xh1p(l-1)
         pl = pl + xi
      xh1p(l)=- pl     * xh1(l)   + xh1(l-1)
   4  continue
      ifail = 0
      return
    5      write(6,10) ifail
   10      format( 'sbesjh : ifail = ',i4)
      return
      end


! From utils.f90
      subroutine gauss(a,b,npoint,xri,wri)
c*****************************************************************
      implicit real*8(a-h,o-z)
      dimension xg(npoint),wg(npoint),xri(npoint),wri(npoint)
      if(npoint.le.0) return
      call setmgl(npoint,xg,wg)
      do 20 j=1,npoint
      xri(j) = (a+b)/2.0 + (b-a)/2.0*xg(j)
      wri(j) = (b-a)/2.0*wg(j)
   20 continue
      return
      end


! From utils.f90
        SUBROUTINE LPMN(MM,M,N,X,PM,PD)
C
C       =====================================================
C       Purpose: Compute the associated Legendre functions 
C                Pmn(x) and their derivatives Pmn'(x)
C       Input :  x  --- Argument of Pmn(x)
C                m  --- Order of Pmn(x),  m = 0,1,2,...,n
C                n  --- Degree of Pmn(x), n = 0,1,2,...,N
C                mm --- Physical dimension of PM and PD
C       Output:  PM(m,n) --- Pmn(x)
C                PD(m,n) --- Pmn'(x)
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PM(0:MM,0:N),PD(0:MM,0:N)
        DO 10 I=0,N
        DO 10 J=0,M
           PM(J,I)=0.0D0
10         PD(J,I)=0.0D0
        PM(0,0)=1.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 I=1,N
              PM(0,I)=X**I
15            PD(0,I)=0.5D0*I*(I+1.0D0)*X**(I+1)
           DO 20 J=1,N
           DO 20 I=1,M
              IF (I.EQ.1) THEN
                 PD(I,J)=1.0D+300
              ELSE IF (I.EQ.2) THEN
                 PD(I,J)=-0.25D0*(J+2)*(J+1)*J*(J-1)*X**(J+1)
              ENDIF
20         CONTINUE
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XQ=DSQRT(LS*(1.0D0-X*X))
        XS=LS*(1.0D0-X*X)
        DO 30 I=1,M
30         PM(I,I)=-LS*(2.0D0*I-1.0D0)*XQ*PM(I-1,I-1)
        DO 35 I=0,M
35         PM(I,I+1)=(2.0D0*I+1.0D0)*X*PM(I,I)
        DO 40 I=0,M
        DO 40 J=I+2,N
           PM(I,J)=((2.0D0*J-1.0D0)*X*PM(I,J-1)-
     &             (I+J-1.0D0)*PM(I,J-2))/(J-I)
40      CONTINUE
        PD(0,0)=0.0D0
        DO 45 J=1,N
45         PD(0,J)=LS*J*(PM(0,J-1)-X*PM(0,J))/XS
        DO 50 I=1,M
        DO 50 J=I,N
           PD(I,J)=LS*I*X*PM(I,J)/XS+(J+I)
     &             *(J-I+1.0D0)/XQ*PM(I-1,J)
50      CONTINUE
        RETURN
        END


! From utils.f90
      SUBROUTINE GAUSS5_Q(N,NR,A,M,SING,DET,EPS,SHOW)
C
C    SOLVE BY GAUSSIAN ELIMINATION SUM(J): A(I,J).P(J) = A(I,N+1)
C             WHERE P(J) IS LEFT IN MAT(J,N+1)
C
      IMPLICIT REAL*16(A-H,O-Z)
      parameter(KO=99)
      COMPLEX*32 A(NR,N+M),DET,RA
      LOGICAL SING,SHOW
      SING  = .FALSE.
      ZERO = 0.0Q0
      ONE = 1.0Q0
      NPM = N + M
      DO 201 I=1,N
201   IF(SHOW) WRITE(KO,402) (A(I,J)  ,J=1,NPM)
      DET = ZERO
      DO 9 K = 1, N
      IF (abs(A(K,K)) .NE. ZERO ) GO TO 5
         SING = .TRUE.
         WRITE(KO,3) K,DET/LOG(10.0Q0)
    3  FORMAT(//' THE MATRIX IS SINGULAR AT',I3,', Log10 determinant is
     &  ',2E16.8/)
      DO 401 I=1,N
401   IF(SHOW) WRITE(KO,402) (            A(I,J)  ,J=1,NPM)
C402   FORMAT( 1X,20F6.1/(1X,20F6.3))
402   FORMAT( 1X,1P,14E9.1/(1X,24E9.1))
         RETURN
    5 KP1 = K + 1
      DET = DET + LOG(A(K,K))
         RA = ONE/A(K,K)
      DO 6 J = KP1, NPM
    6 A(K,J) = A(K,J) * RA
      A(K,K) = ONE
      DO 9 I = 1, N
      IF (I .EQ. K  .OR. ABS (A(I,K)) .EQ. 0) GO TO 9
CDIR$ IVDEP
         DO 8 J = KP1, NPM
    8    A(I,J) = A(I,J) - A(I,K)*A(K,J)
         A(I,K) = ZERO
    9 CONTINUE
      IF(SHOW) WRITE(KO,15 ) DET/LOG(10.0Q0)
15    FORMAT(/' Log10 determinant is ',2F10.5)
      IF(SHOW) WRITE(KO,402) (A(I,N+1),I=1,N)
      RETURN
      END


! From utils.f90
      subroutine spline2(n,x,y,b,c,d)
      implicit real*8(a-h,o-z)
      dimension x(n),y(n),b(n),c(n),d(n)
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
c
      nm1 = n-1
      if (n.lt.2) return
      if (n.lt.3) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1)=x(2)-x(1)
      c(2)=(y(2)-y(1))/d(1)
      do 10 i =2,nm1
      d(i)=x(i+1)-x(i)
      b(i)=2.d0*(d(i-1)+d(i))
      c(i+1)=(y(i+1)-y(i))/d(i)
      c(i)=c(i+1)-c(i)
10    continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1)=-d(1)
      b(n)=-d(n-1)
      c(1)=0.d0
      c(n)=0.d0
      if(n.eq.3) go to 15
      c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
      c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
      c(1)=c(1)*d(1)**2/(x(4)-x(1))
      c(n)=-c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
15    do 20 i=2,n
      t=d(i-1)/b(i-1)
      b(i)=b(i)-t*d(i-1)
      c(i)=c(i)-t*c(i-1)
20    continue
c
c  back substitution
c
      c(n)=c(n)/b(n)
      do 30 ib=1,nm1
      i=n-ib
      c(i)=(c(i)-d(i)*c(i+1))/b(i)
30    continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n)=(y(n)-y(nm1))/d(nm1)+d(nm1)*(c(nm1)+2.d0*c(n))
      do 40 i=1,nm1
      b(i)=(y(i+1)-y(i))/d(i)-d(i)*(c(i+1)+2.d0*c(i))
      d(i)=(c(i+1)-c(i))/d(i)
      c(i)=3.d0*c(i)
40    continue
      c(n)=3.d0*c(n)
      d(n)=d(n-1)
      return
c
50    b(1)=(y(2)-y(1))/(x(2)-x(1))
      c(1)=0.d0
      d(1)=0.d0
      b(2)=b(1)
      c(2)=0.d0
      d(2)=0.d0
      return
      end 


! From pauli.f90
	subroutine pauliproj2(nchani,nchanf,ndim)
	use forbidden
	use globals, only: written
	use wfs !, only: wfsp,ql,qj,spindex,nr,rvec
        use channels
	implicit none
	integer:: nchani,nchanf
	integer:: ip,n,m,ndim,ichn,ichm,li,lf,ir
	real*8:: ji,jf,un,um,r,up,res1,res2
	real*8,allocatable:: wfpau(:)
	real*8,allocatable:: fpaux(:),gpaux(:)
	real*8,allocatable:: fmaux(:),fnaux(:)
	real*8 iden(ndim,ndim)
	 
       ichn=spindex(nchani)
       ichm=spindex(nchanf)
       li=ql(nchani)	
       lf=ql(nchanf)
       ji=qj(nchani)
       jf=qj(nchanf)
       if (allocated(fnaux)) deallocate(fnaux)
       if (allocated(fmaux)) deallocate(fmaux)
       if (allocated(fpaux)) deallocate(fpaux)
       if (allocated(gpaux)) deallocate(gpaux)
       if (allocated(wfpau)) deallocate(wfpau)
       if (allocated(ppauli)) deallocate(ppauli)
       if (allocated(qpauli)) deallocate(qpauli)
!	write(*,*)'pauli: nr,ndim',nr,ndim
	allocate(fpaux(nr),gpaux(nr))
	allocate(fnaux(nr),fmaux(nr))
	allocate(wfpau(nr))
	allocate(ppauli(ndim,ndim),qpauli(ndim,ndim))
	
	
	fnaux(:)=0d0; fmaux(:)=0d0
	ppauli(:,:)=0d0
	iden(1:ndim,1:ndim)=0d0

	do ip=1,npauli
        wfpau(:)=wfsp(ichn,paun(ip),:) !Pauli forbidden WF

	 do n=1,ndim
	  do m=1,ndim
	  if (n.eq.m) iden(n,n)=1d0
         fpaux(:)=0d0; gpaux(:)=0d0
         ppauli(n,m)=0d0
        if (li.eq.paul(ip).and.(ji.eq.pauj(ip))
     &  .and.(lf.eq.paul(ip)).and.(jf.eq.pauj(ip)) 
     &  .and.(nchani.eq.nchanf)) then
!     &  .and.(n.eq.paun(ip))) then
	  do ir=1,nr
          r=rvec(ir)
          un=wfaux(n,ir)*r
          um=wfaux2(m,ir)*r
	  up=wfpau(ir)*r
          fpaux(ir)=fpaux(ir)+un*up ! <n|p>
	   gpaux(ir)=gpaux(ir)+um*up ! <p|m>
          enddo !ir
         endif

	 call sim(fpaux,res1,1,nr,dr,nr)
        call sim(gpaux,res2,1,nr,dr,nr)
	 ppauli(n,m)=ppauli(n,m)+res1*res2
        enddo !n
        enddo !m

c write block
	 written(80)=.true.
	 write(80,*)'# Pauli projector for ichan block',nchani,nchanf
         call writemat(ppauli,ndim,80)
        enddo !ip=1,npauli
	
c Q= 1 - P = 1 - Sum |forbidden> <forbidden|
	 qpauli=iden-ppauli
	

	end subroutine 


! From read_frad.f90
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


! From rmatrix_hp.F90
subroutine wf_print(nch, lval, qk, eta, rmax, nr, ns, cu, &
                    ndim, nopen, cf, nwf1, nwf2, zrma, iv, nom, npoin, h, cwftab)
  implicit none
  integer, intent(in) :: nch, nr, ns, ndim, nopen, nwf1, nwf2, nom, iv, npoin
  integer, intent(in) :: lval(nch)
  real*8, intent(in) :: qk(nch), eta(nch), rmax, zrma(ns*nr), h
  complex*16, intent(in) :: cu(ndim, ndim), cf(nwf1, nwf2, nom)
  complex*16, intent(out) :: cwftab(npoin)

  integer :: nsr, nop, nclo, iw, i, ll, ifail, jw
  real*8 :: r, wfr, wfi, xl, yp1, yp2
  real*8 :: xfr(ns*nr), xfi(ns*nr), y2r(ns*nr), y2i(ns*nr)
  real*8 :: fc(500), dfc(500), gc(500), dgc(500)
  complex*16 :: co, cfx

  yp1 = 1.0d30
  yp2 = 1.0d30

  nsr = ns * nr
  nop = 0
  nclo = nopen

  do iw = 1, iv
    if (qk(iw) > 0) nop = nop + 1
    if (qk(iw) < 0) nclo = nclo + 1
  end do

  xfr(1:nsr) = real(cf(1:nsr, iv, nom))
  xfi(1:nsr) = aimag(cf(1:nsr, iv, nom))

  call spline(zrma, xfr, nsr, yp1, yp2, y2r)
  call spline(zrma, xfi, nsr, yp1, yp2, y2i)

  do i = 1, npoin
    r = i * h
    if (r <= rmax) then
      call splint(zrma, xfr, y2r, nsr, r, wfr)
      call splint(zrma, xfi, y2i, nsr, r, wfi)
      cwftab(i) = dcmplx(wfr, wfi)
    else
      ll = lval(iv)
      if (qk(iv) > 0) then
        xl = lval(iv)
        call coulfg(qk(iv)*r, eta(iv), xl, xl, fc, gc, dfc, dgc, 1, 0, ifail)
        co = dcmplx(gc(ll+1), fc(ll+1)) * sqrt(qk(nom)/qk(iv))
        cfx = -cu(nop, nom) * co
        if (iv == nom) cfx = cfx + conjg(co)
      else
        jw = 0
        call whit(eta(iv), r, -qk(iv), -qk(iv)**2, ll, fc, dfc, jw)
        cfx = cu(nclo, nom) * fc(ll+1)
      end if
      cwftab(i) = cfx
    end if
  end do

end subroutine wf_print


! From ceigen.F
	subroutine CEigensystem(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	Complex A(ldA,*), U(ldU,*), d(*)

	integer p, q, j
	Real red, off, thresh
	Complex delta, t, s, invc, sx, sy, tx, ty
	Complex x, y
	Complex ev(2,MAXDIM)

	integer sweep
	common /nsweeps/ sweep

	Real Sq
	Complex c
	Sq(c) = Re(c*Conjugate(c))

	if( n .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = A(p,p)
	  d(p) = ev(2,p)
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = 0
	  enddo
	  U(p,p) = 1
	enddo

	red = .01D0/n**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, n
	    do p = 1, q - 1
	      off = off + Sq(A(p,q)) + Sq(A(q,p))
	    enddo
	  enddo
	  if( .not. off .gt. EPS ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
	      off = Sq(A(p,q)) + Sq(A(q,p))
	      if( sweep .gt. 4 .and. off .lt.
     &              EPS*(Sq(ev(2,p)) + Sq(ev(2,q))) ) then
	        A(p,q) = 0
	        A(q,p) = 0
	      else if( off .gt. thresh ) then
	        delta = A(p,q)*A(q,p)
	        x = .5D0*(ev(2,p) - ev(2,q))
	        y = sqrt(x**2 + delta)
	        t = x - y
	        s = x + y
	        if( Sq(t) .lt. Sq(s) ) t = s

	        t = 1/t
	        delta = delta*t
	        ev(1,p) = ev(1,p) + delta
	        ev(2,p) = d(p) + ev(1,p)
	        ev(1,q) = ev(1,q) - delta
	        ev(2,q) = d(q) + ev(1,q)

	        invc = sqrt(delta*t + 1)
	        s = t/invc
	        t = t/(invc + 1)
	        sx = s*A(p,q)
	        ty = t*A(p,q)
	        sy = s*A(q,p)
	        tx = t*A(q,p)

	        do j = 1, n
	          x = A(j,p)
	          y = A(j,q)
	          A(j,p) = x + sy*(y - ty*x)
	          A(j,q) = y - sx*(x + tx*y)
	          x = A(p,j)
	          y = A(q,j)
	          A(p,j) = x + sx*(y - tx*x)
	          A(q,j) = y - sy*(x + ty*y)
	        enddo

	        A(p,q) = 0
	        A(q,p) = 0

	        do j = 1, n
	          x = U(p,j)
	          y = U(q,j)
	          U(p,j) = x + sx*(y - tx*x)
	          U(q,j) = y - sy*(x + ty*y)
	        enddo
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in CEigensystem"

1	if( sort .eq. 0 ) return

* sort the eigenvalues by their real part

	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
	    if( sort*(Re(t) - Re(d(q))) .gt. 0 ) then
	      j = q
	      t = d(q)
	    endif
	  enddo

	  if( j .ne. p ) then
	    d(j) = d(p)
	    d(p) = t
	    do q = 1, n
	      x = U(p,q)
	      U(p,q) = U(j,q)
	      U(j,q) = x
	    enddo
	  endif
	enddo
	end


! From ceigen.F
	subroutine HEigensystem(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	Complex A(ldA,*), U(ldU,*)
	Real d(*)

	integer p, q, j
	Real red, off, thresh
	Real t, delta, invc, s
	Complex x, y, Apq
	Real ev(2,MAXDIM)

	integer sweep
	common /nsweeps/ sweep

	Real Sq
	Complex c
	Sq(c) = Re(c*Conjugate(c))

	if( n .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = Re(A(p,p))
	  d(p) = ev(2,p)
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = 0
	  enddo
	  U(p,p) = 1
	enddo

	red = .04D0/n**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, n
	    do p = 1, q - 1
	      off = off + Sq(A(p,q))
	    enddo
	  enddo
	  if( .not. off .gt. SYM_EPS ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
	      Apq = A(p,q)
	      off = Sq(Apq)
	      if( sweep .gt. 4 .and. off .lt.
     &              SYM_EPS*(ev(2,p)**2 + ev(2,q)**2) ) then
	        A(p,q) = 0
	      else if( off .gt. thresh ) then
	        t = .5D0*(ev(2,p) - ev(2,q))
	        t = 1/(t + sign(sqrt(t**2 + off), t))

	        delta = t*off
	        ev(1,p) = ev(1,p) + delta
	        ev(2,p) = d(p) + ev(1,p)
	        ev(1,q) = ev(1,q) - delta
	        ev(2,q) = d(q) + ev(1,q)

	        invc = sqrt(delta*t + 1)
	        s = t/invc
	        t = delta/(invc + 1)

	        do j = 1, p - 1
	          x = A(j,p)
	          y = A(j,q)
	          A(j,p) = x + s*(Conjugate(Apq)*y - t*x)
	          A(j,q) = y - s*(Apq*x + t*y)
	        enddo

	        do j = p + 1, q - 1
	          x = A(p,j)
	          y = A(j,q)
	          A(p,j) = x + s*(Apq*Conjugate(y) - t*x)
	          A(j,q) = y - s*(Apq*Conjugate(x) + t*y)
	        enddo

	        do j = q + 1, n
	          x = A(p,j)
	          y = A(q,j)
	          A(p,j) = x + s*(Apq*y - t*x)
	          A(q,j) = y - s*(Conjugate(Apq)*x + t*y)
	        enddo

	        A(p,q) = 0

	        do j = 1, n
	          x = U(p,j)
	          y = U(q,j)
	          U(p,j) = x + s*(Apq*y - t*x)
	          U(q,j) = y - s*(Conjugate(Apq)*x + t*y)
	        enddo
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in HEigensystem"

1	if( sort .eq. 0 ) return

* sort the eigenvalues

	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
	    if( sort*(t - d(q)) .gt. 0 ) then
	      j = q
	      t = d(q)
	    endif
	  enddo

	  if( j .ne. p ) then
	    d(j) = d(p)
	    d(p) = t
	    do q = 1, n
	      x = U(p,q)
	      U(p,q) = U(j,q)
	      U(j,q) = x
	    enddo
	  endif
	enddo
	end


! From ceigen.F
	subroutine SVD(m, n, Ao,ldA, d, Vo,ldV, Wo,ldW, sort)
	implicit none
	integer m, n, ldA, ldV, ldW, sort
	Complex Ao(ldA,*), Vo(ldV,*), Wo(ldW,*)
	Real d(*)

	integer nx, nm, p, q, px, qx, j, rev, pi(MAXDIM)
	Real red, off, thresh
	Real t, dv, dw, xv, xw, invc
	Complex App, Apq, Aqp, Aqq
	Complex x, y, sv, sw, tv, tw, f
	Complex VW(MAXDIM,MAXDIM,0:2)

* note: for better cache efficiency, the Vx, Wx arrays
* contain the *transpose* of the transformation matrices
	Complex V(MAXDIM,MAXDIM)
	Complex W(MAXDIM,MAXDIM)
	Complex A(MAXDIM,MAXDIM)
	equivalence (VW(1,1,0), V)
	equivalence (VW(1,1,1), W)
	equivalence (VW(1,1,2), A)

	integer sweep
	common /nsweeps/ sweep

	Real Sq
	Complex c
	Sq(c) = Re(c*Conjugate(c))

	nx = max(m, n)

	if( nx .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, nx
	  do q = 1, nx
	    V(q,p) = 0
	    W(q,p) = 0
	    A(q,p) = 0
	  enddo
	  V(p,p) = 1
	  W(p,p) = 1
	enddo

	rev = ibits(m - n, 15, 1)
	if( rev .eq. 1 ) then
	  do p = 1, n
	    do q = 1, m
	      A(p,q) = Ao(q,p)
	    enddo
	  enddo
	else
	  do p = 1, n
	    do q = 1, m
	      A(q,p) = Ao(q,p)
	    enddo
	  enddo
	endif

	red = .01D0/nx**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, nx
	    do p = 1, q - 1
	      off = off + Sq(A(p,q)) + Sq(A(q,p))
	    enddo
	  enddo
	  if( .not. off .gt. EPS ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, nx
	    do p = 1, q - 1
	      px = p
	      qx = q
	      if( Sq(A(p,p)) + Sq(A(q,q)) .lt.
     &            Sq(A(p,q)) + Sq(A(q,p)) ) then
	        px = q
	        qx = p
	      endif

	      App = A(px,p)
	      Aqq = A(qx,q)
	      Apq = A(px,q)
	      Aqp = A(qx,p)
	      off = Sq(Apq) + Sq(Aqp)
	      if( sweep .gt. 4 .and.
     &            off .lt. EPS*(Sq(App) + Sq(Aqq)) ) then
	        A(px,q) = 0
	        A(qx,p) = 0
	      else if( off .gt. thresh ) then
	        xv = Re((App - Aqq)*Conjugate(App + Aqq))
	        xw = Re((Apq - Aqp)*Conjugate(Apq + Aqp))
	        dv = .5D0*(xv + xw)
	        dw = .5D0*(xv - xw)

	        tv = Conjugate(App)*Aqp + Aqq*Conjugate(Apq)
	        tw = Conjugate(App)*Apq + Aqq*Conjugate(Aqp)
c	        t = sqrt(dv**2 + Sq(tv))
	        t = sqrt(dw**2 + Sq(tw))

	        xv = min(abs(dv + t), abs(dw + t))
	        xw = min(abs(dv - t), abs(dw - t))
	        if( xv + xw .gt. DBL_EPS ) then
	          t = sign(t, xv - xw)
	          tv = tv/(dv + t)
	          tw = tw/(dw + t)
	        else
	          tv = 0
	          tw = Apq/App
	        endif

	        invc = sqrt(1 + Sq(tv))
	        sv = tv/invc
	        tv = tv/(invc + 1)

	        invc = sqrt(1 + Sq(tw))
	        sw = tw/invc
	        tw = tw/(invc + 1)

	        do j = 1, nx
	          x = A(j,p)
	          y = A(j,q)
	          A(j,p) = x + Conjugate(sw)*(y - tw*x)
	          A(j,q) = y - sw*(x + Conjugate(tw)*y)
	          x = A(px,j)
	          y = A(qx,j)
	          A(p,j) = x + Conjugate(sv)*(y - tv*x)
	          A(q,j) = y - sv*(x + Conjugate(tv)*y)
	        enddo

	        A(p,p) = invc*(App + Conjugate(sv)*(Aqp - tv*App))
	        A(q,p) = 0
	        A(p,q) = 0
	        A(q,q) = invc*(Aqq - sv*(Apq + Conjugate(tv)*Aqq))

	        do j = 1, nx
	          x = V(j,px)
	          y = V(j,qx)
	          V(j,p) = x + sv*(y - Conjugate(tv)*x)
	          V(j,q) = y - Conjugate(sv)*(x + tv*y)
	        enddo

	        do j = 1, nx
	          x = W(j,p)
	          y = W(j,q)
	          W(j,p) = x + sw*(y - Conjugate(tw)*x)
	          W(j,q) = y - Conjugate(sw)*(x + tw*y)
	        enddo
	        goto 2
	      endif

	      if( p .ne. px ) then
	        do j = 1, nx
	          x = A(p,j)
	          A(p,j) = A(q,j)
	          A(q,j) = x
	        enddo

	        do j = 1, nx
	          x = V(j,p)
	          V(j,p) = V(j,q)
	          V(j,q) = x
	        enddo
	      endif

2	      continue
	    enddo
	  enddo
	enddo

	print *, "Bad convergence in SVD"

1	continue

	nm = min(m, n)

* make the diagonal elements nonnegative

	do p = 1, nm
	  d(p) = abs(A(p,p))
	  if( d(p) .gt. DBL_EPS .and. d(p) .ne. Re(A(p,p)) ) then
	    f = A(p,p)/d(p)
	    do q = 1, nm
	      W(q,p) = W(q,p)*f
	    enddo
	  endif
	enddo

* sort the singular values

	do p = 1, nm
	  pi(p) = p
	enddo

	do p = 1, nm
	  j = p
	  t = d(p)
	  if( sort .ne. 0 ) then
	    do q = p + 1, nm
	      if( sort*(t - d(q)) .gt. 0 ) then
	        j = q
	        t = d(q)
	      endif
	    enddo
	  endif

	  d(j) = d(p)
	  d(p) = t

	  q = pi(j)
	  pi(j) = pi(p)

	  do j = 1, m
	    Vo(p,j) = VW(j,q,rev)
	  enddo
	  do j = 1, n
	    Wo(p,j) = VW(j,q,1-rev)
	  enddo
	enddo
	end


! From ceigen.F
	subroutine TakagiFactor(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	Complex A(ldA,*), U(ldU,*)
	Real d(*)

	integer p, q, j
	Real red, off, thresh
	Real sqp, sqq, t, invc
	Complex f, x, y
	Complex ev(2,MAXDIM)

	integer sweep
	common /nsweeps/ sweep

	Real Sq
	Complex c
	Sq(c) = Re(c*Conjugate(c))

	if( n .gt. MAXDIM ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = A(p,p)
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = 0
	  enddo
	  U(p,p) = 1
	enddo

	red = .04D0/n**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, n
	    do p = 1, q - 1
	      off = off + Sq(A(p,q))
	    enddo
	  enddo
	  if( .not. off .gt. SYM_EPS ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
	      off = Sq(A(p,q))
	      sqp = Sq(ev(2,p))
	      sqq = Sq(ev(2,q))
	      if( sweep .gt. 4 .and.
     &            off .lt. SYM_EPS*(sqp + sqq) ) then
	        A(p,q) = 0
	      else if( off .gt. thresh ) then
	        t = .5D0*abs(sqp - sqq)
	        if( t .gt. 0 ) then
	          f = sign(1D0, sqp - sqq)*
     &              (ev(2,q)*Conjugate(A(p,q)) +
     &               Conjugate(ev(2,p))*A(p,q))
	        else
	          f = 1
	          if( sqp .ne. 0 ) f = sqrt(ev(2,q)/ev(2,p))
	        endif
	        t = t + sqrt(t**2 + Sq(f))
	        f = f/t

	        ev(1,p) = ev(1,p) + A(p,q)*Conjugate(f)
	        ev(2,p) = A(p,p) + ev(1,p)
	        ev(1,q) = ev(1,q) - A(p,q)*f
	        ev(2,q) = A(q,q) + ev(1,q)

	        t = Sq(f)
	        invc = sqrt(t + 1)
	        f = f/invc
	        t = t/(invc*(invc + 1))

	        do j = 1, p - 1
	          x = A(j,p)
	          y = A(j,q)
	          A(j,p) = x + (Conjugate(f)*y - t*x)
	          A(j,q) = y - (f*x + t*y)
	        enddo

	        do j = p + 1, q - 1
	          x = A(p,j)
	          y = A(j,q)
	          A(p,j) = x + (Conjugate(f)*y - t*x)
	          A(j,q) = y - (f*x + t*y)
	        enddo

	        do j = q + 1, n
	          x = A(p,j)
	          y = A(q,j)
	          A(p,j) = x + (Conjugate(f)*y - t*x)
	          A(q,j) = y - (f*x + t*y)
	        enddo

	        A(p,q) = 0

	        do j = 1, n
	          x = U(p,j)
	          y = U(q,j)
	          U(p,j) = x + (f*y - t*x)
	          U(q,j) = y - (Conjugate(f)*x + t*y)
	        enddo
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    A(p,p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in TakagiFactor"

1	continue

* make the diagonal elements nonnegative

	do p = 1, n
	  d(p) = abs(A(p,p))
	  if( d(p) .gt. DBL_EPS .and. d(p) .ne. Re(A(p,p)) ) then
	    f = sqrt(A(p,p)/d(p))
	    do q = 1, n
	      U(p,q) = U(p,q)*f
	    enddo
	  endif
	enddo

	if( sort .eq. 0 ) return

* sort the eigenvalues

	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
	    if( sort*(t - d(q)) .gt. 0 ) then
	      j = q
	      t = d(q)
	    endif
	  enddo

	  if( j .ne. p ) then
	    d(j) = d(p)
	    d(p) = t
	    do q = 1, n
	      x = U(p,q)
	      U(p,q) = U(j,q)
	      U(j,q) = x
	    enddo
	  endif
	enddo
	end


! From rmatrix.f
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


! From thox.f90
        subroutine cfival(lrf,lfv,gaux,nlf,alpha)
        real*8 lrf,lfv(:),alpha
        complex*16:: gaux(:)
        integer nlf
        end subroutine


! From thox.f90
        subroutine test_cont(nset,nchan,inc,ecm,wf)
        integer   :: nset,nchan,inc
        real*8    :: ecm
        complex*16:: wf(:,:)
        end subroutine


! From utils/pwl_interp_2d.f90
subroutine pwl_interp_2d ( nxd, nyd, xd, yd, zd, ni, xi, yi, zi )

!*****************************************************************************80
!
!! PWL_INTERP_2D: piecewise linear interpolant to data defined on a 2D grid.
!
!  Discussion:
!
!    Thanks to Adam Hirst for pointing out an error in the formula that
!    chooses the interpolation triangle, 04 February 2018.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NXD, NYD, the number of X and Y data values.
!
!    Input, real ( kind = 8 ) XD(NXD), YD(NYD), the sorted X and Y data.
!
!    Input, real ( kind = 8 ) ZD(NXD,NYD), the Z data.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), YI(NI), the coordinates of the
!    interpolation points.
!
!    Output, real ( kind = 8 ) ZI(NI), the value of the interpolant.
!
  implicit none

  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nxd
  integer ( kind = 4 ) nyd

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) det
  real ( kind = 8 ) dxa
  real ( kind = 8 ) dxb
  real ( kind = 8 ) dxi
  real ( kind = 8 ) dya
  real ( kind = 8 ) dyb
  real ( kind = 8 ) dyi
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_huge
  integer ( kind = 4 ) r8vec_bracket5
  real ( kind = 8 ) xd(nxd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yd(nyd)
  real ( kind = 8 ) yi(ni)
  real ( kind = 8 ) zd(nxd,nyd)
  real ( kind = 8 ) zi(ni)

  do k = 1, ni
!
!  For interpolation point (xi(k),yi(k)), find data intervals I and J so that:
!
!    xd(i) <= xi(k) <= xd(i+1),
!    yd(j) <= yi(k) <= yd(j+1).
!
!  But if the interpolation point is not within a data interval, 
!  assign the dummy interpolant value zi(k) = infinity.
!
    i = r8vec_bracket5 ( nxd, xd, xi(k) )
    if ( i == -1 ) then
      zi(k) = r8_huge ( )
      cycle
    end if

    j = r8vec_bracket5 ( nyd, yd, yi(k) )
    if ( j == -1 ) then
      zi(k) = r8_huge ( )
      cycle
    end if
!
!  The rectangular cell is arbitrarily split into two triangles.
!  The linear interpolation formula depends on which triangle 
!  contains the data point.
!
!    (I,J+1)--(I+1,J+1)
!      |\       |
!      | \      |
!      |  \     |
!      |   \    |
!      |    \   |
!      |     \  |
!    (I,J)---(I+1,J)
!
    if ( yi(k) < yd(j+1) &
      + ( yd(j) - yd(j+1) ) * ( xi(k) - xd(i) ) / ( xd(i+1) - xd(i) ) ) then

      dxa = xd(i+1) - xd(i)
      dya = yd(j)   - yd(j)

      dxb = xd(i)   - xd(i)
      dyb = yd(j+1) - yd(j)

      dxi = xi(k)   - xd(i)
      dyi = yi(k)   - yd(j)

      det = dxa * dyb - dya * dxb

      alpha = ( dxi * dyb - dyi * dxb ) / det
      beta =  ( dxa * dyi - dya * dxi ) / det
      gamma = 1.0D+00 - alpha - beta

      zi(k) = alpha * zd(i+1,j) + beta * zd(i,j+1) + gamma * zd(i,j)

    else

      dxa = xd(i)   - xd(i+1)
      dya = yd(j+1) - yd(j+1)

      dxb = xd(i+1) - xd(i+1)
      dyb = yd(j)   - yd(j+1)

      dxi = xi(k)   - xd(i+1)
      dyi = yi(k)   - yd(j+1)

      det = dxa * dyb - dya * dxb

      alpha = ( dxi * dyb - dyi * dxb ) / det
      beta =  ( dxa * dyi - dya * dxi ) / det
      gamma = 1.0D+00 - alpha - beta

      zi(k) = alpha * zd(i,j+1) + beta * zd(i+1,j) + gamma * zd(i+1,j+1)

    end if

  end do

  return
end


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
