!http://www.feynarts.de/diag/

* CEigensystem.F
* diagonalization of a complex n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 9 Aug 11 th

#include "diag-f.h"


************************************************************************
** CEigensystem diagonalizes a general complex n-by-n matrix.
** Input: n, A = n-by-n matrix
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^-1.

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

* HEigensystem.F
* diagonalization of a Hermitian n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 9 Aug 11 th

#include "diag-f.h"


************************************************************************
** HEigensystem diagonalizes a Hermitian n-by-n matrix.
** Input: n, A = n-by-n matrix, Hermitian
** (only the upper triangle of A needs to be filled).
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^+ = U A U^-1 with U unitary.

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

* SEigensystem.F
* diagonalization of a complex symmetric n-by-n matrix using
* the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 9 Aug 11 th

#include "diag-f.h"


************************************************************************
** SEigensystem diagonalizes a complex symmetric n-by-n matrix.
** Input: n, A = n-by-n matrix, complex symmetric
** (only the upper triangle of A needs to be filled).
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^T = U A U^-1 with U U^T = 1.

	subroutine SEigensystem(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	Complex A(ldA,*), U(ldU,*), d(*)

	integer p, q, j
	Real red, off, thresh
	Complex delta, t, invc, s
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
	      delta = A(p,q)
	      off = Sq(delta)
	      if( sweep .gt. 4 .and. off .lt.
     &              SYM_EPS*(Sq(ev(2,p)) + Sq(ev(2,q))) ) then
	        A(p,q) = 0
	      else if( off .gt. thresh ) then
	        x = .5D0*(ev(2,p) - ev(2,q))
	        y = sqrt(x**2 + delta**2)
	        t = x - y
	        s = x + y
	        if( Sq(t) .lt. Sq(s) ) t = s

	        t = delta/t
	        delta = delta*t
	        ev(1,p) = ev(1,p) + delta
	        ev(2,p) = d(p) + ev(1,p)
	        ev(1,q) = ev(1,q) - delta
	        ev(2,q) = d(q) + ev(1,q)

	        invc = sqrt(t**2 + 1)
	        s = t/invc
	        t = t/(invc + 1)

	        do j = 1, p - 1
	          x = A(j,p)
	          y = A(j,q)
	          A(j,p) = x + s*(y - t*x)
	          A(j,q) = y - s*(x + t*y)
	        enddo

	        do j = p + 1, q - 1
	          x = A(p,j)
	          y = A(j,q)
	          A(p,j) = x + s*(y - t*x)
	          A(j,q) = y - s*(x + t*y)
	        enddo

	        do j = q + 1, n
	          x = A(p,j)
	          y = A(q,j)
	          A(p,j) = x + s*(y - t*x)
	          A(q,j) = y - s*(x + t*y)
	        enddo

	        A(p,q) = 0

	        do j = 1, n
	          x = U(p,j)
	          y = U(q,j)
	          U(p,j) = x + s*(y - t*x)
	          U(q,j) = y - s*(x + t*y)
	        enddo
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in SEigensystem"

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

* SVD.F
* singular value decomposition of an m-by-n matrix
* this file is part of the Diag library
* last modified 9 Aug 11 th

#include "diag-f.h"


************************************************************************
** SVD performs a singular value decomposition.
** Input: m, n, A = m-by-n matrix.
** Output: d = nm-vector of singular values,
** V = nm-by-m left transformation matrix,
** W = nm-by-n right transformation matrix, nm = min(m, n),
** these fulfill diag(d) = V^* A W^+.

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

* TakagiFactor.F
* computes the Takagi factorization of a complex symmetric matrix
* code adapted from the "Handbook" routines
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 9 Aug 11 th

#include "diag-f.h"


************************************************************************
** TakagiFactor factorizes a complex symmetric n-by-n matrix
** Input: n, A = n-by-n matrix, complex symmetric
** (only the upper triangle of A needs to be filled).
** Output: d = vector of diagonal values, U = transformation matrix
** these fulfill diag(d) = U^* A U^+ with U unitary.

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

