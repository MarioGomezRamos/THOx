      subroutine f01akf(n,k,l,a,ia,intger)
c     mark 2 release. nag copyright 1972
c     mark 4 revised.
c     mark 4.5 revised
c     mark 8 revised. ier-248 (jun 1980).
c     mark 11 revised. vectorisation (jan 1984).
c     mark 11.5(f77) revised. (sept 1985.)
c     mark 12 revised. extended blas (june 1986)
c
c     dirhes
c     august 1st, 1971 .
c     given the unsymmetric matrix, a, stored in the array a(n,n),
c     this subroutine reduces the sub-matrix of order l - k + 1,
c     which starts at the element a(k,k) and finishes at the
c     element a(l,l), to hessenberg form, h, by the direct
c     method(an = nh). the matrix h is overwritten on a with
c     details of the transformations (n) stored in the remaining
c     triangle under h and in elements k to l of the array
c     intger(n).
c     1st august 1971
c
c     .. scalar arguments ..
      integer           ia, k, l, n
c     .. array arguments ..
      double precision  a(ia,n)
      integer           intger(n)
c     .. local scalars ..
      double precision  x, y
      integer           i, j, k1, m
c     .. external subroutines ..
      external          dgemv, dtrsv
c     .. intrinsic functions ..
      intrinsic         abs
c     .. executable statements ..
      k1 = k + 1
      if (k1.gt.l) return
      do 140 j = k1, n
         m = j
         x = 0.0d0
         if (j.gt.l) go to 120
         do 20 i = j, l
            if (abs(a(i,j-1)).le.abs(x)) go to 20
            x = a(i,j-1)
            m = i
   20    continue
         intger(j) = m
         if (m.eq.j) go to 80
c        interchange rows and columns of a.
         do 40 i = k, n
            y = a(m,i)
            a(m,i) = a(j,i)
            a(j,i) = y
   40    continue
         do 60 i = 1, l
            y = a(i,m)
            a(i,m) = a(i,j)
            a(i,j) = y
   60    continue
   80    if (x.ne.0.0d0 .and. j.lt.l) then
            do 100 i = j + 1, l
               a(i,j-1) = a(i,j-1)/x
  100       continue
            call dgemv('n',l,l-j,1.0d0,a(1,j+1),ia,a(j+1,j-1),1,1.0d0,
     *                 a(1,j),1)
         end if
  120    call dtrsv('l','n','u',j-k,a(k+1,k),ia,a(k+1,j),1)
         if (j.lt.l) call dgemv('n',l-j,j-k,-1.0d0,a(j+1,k),ia,a(k+1,j),
     *                          1,1.0d0,a(j+1,j),1)
  140 continue
      return
      end
      subroutine f01apf(n,low,iupp,intger,h,ih,v,iv)
c     mark 2 release. nag copyright 1972
c     mark 4 revised.
c     mark 4.5 revised
c     mark 5c revised
c     mark 11.5(f77) revised. (sept 1985.)
c
c     dirtrans
c     form the matrix of accumulated transformations in the array
c     v(n,n) from the information left by subroutine f01akf
c     below the upper hessenberg matrix, h, in the array h(n,n)
c     and in the integer array intger(n).
c     1st august 1971
c
c     .. scalar arguments ..
      integer           ih, iupp, iv, low, n
c     .. array arguments ..
      double precision  h(ih,n), v(iv,n)
      integer           intger(n)
c     .. local scalars ..
      double precision  x
      integer           i, i1, ii, j, low1, m
c     .. executable statements ..
      do 40 i = 1, n
         do 20 j = 1, n
            v(i,j) = 0.0d0
   20    continue
         v(i,i) = 1.0d0
   40 continue
      low1 = low + 1
      if (low1.gt.iupp) return
      do 120 ii = low1, iupp
         i = low1 + iupp - ii
         i1 = i - 1
         if (low1.gt.i1) go to 80
         do 60 j = low1, i1
            v(i,j) = h(i,j-1)
   60    continue
   80    m = intger(i)
         if (m.eq.i) go to 120
         do 100 j = low1, iupp
            x = v(m,j)
            v(m,j) = v(i,j)
            v(i,j) = x
  100    continue
  120 continue
      return
      end
      subroutine f01atf(n,ib,a,ia,low,lhi,d)
c     mark 2 release. nag copyright 1972
c     mark 4 revised.
c     mark 4.5 revised
c     mark 11.5(f77) revised. (sept 1985.)
c
c     balance
c     reduce the norm of a(n,n) by exact diagonal similarity
c     transformations stored in d(n).
c     december 1st.,1971
c
c     .. scalar arguments ..
      integer           ia, ib, lhi, low, n
c     .. array arguments ..
      double precision  a(ia,n), d(n)
c     .. local scalars ..
      double precision  b2, c, f, g, r, s
      integer           i, j, jj, k, l
      logical           noconv
c     .. external subroutines ..
      external          f01atz
c     .. intrinsic functions ..
      intrinsic         abs, dble
c     .. executable statements ..
      b2 = ib*ib
      l = 1
      k = n
   20 if (k.lt.1) go to 100
c     search for rows isolating an eigenvalue and push them down
      j = k + 1
      do 60 jj = 1, k
         j = j - 1
         r = 0.0d0
         do 40 i = 1, k
            if (i.eq.j) go to 40
            r = r + abs(a(j,i))
   40    continue
         if (r.eq.0.0d0) go to 80
   60 continue
      go to 100
   80 call f01atz(k,a,ia,d,k,l,n,j)
      k = k - 1
      go to 20
c     search for columns isolating an eigenvalue and push them
c     left.
  100 if (l.gt.k) go to 180
      do 140 j = l, k
         c = 0.0d0
         do 120 i = l, k
            if (i.eq.j) go to 120
            c = c + abs(a(i,j))
  120    continue
         if (c.eq.0.0d0) go to 160
  140 continue
      go to 180
  160 call f01atz(l,a,ia,d,k,l,n,j)
      l = l + 1
      go to 100
c     now balance the submatrix in rows l through k.
  180 low = l
      lhi = k
      if (l.gt.k) go to 220
      do 200 i = l, k
         d(i) = 1.0d0
  200 continue
  220 noconv = .false.
      if (l.gt.k) go to 420
      do 400 i = l, k
         c = 0.0d0
         r = 0.0d0
         do 240 j = l, k
            if (j.eq.i) go to 240
            c = c + abs(a(j,i))
            r = r + abs(a(i,j))
  240    continue
         g = r/dble(ib)
         f = 1.0d0
         s = c + r
  260    if (c.ge.g) go to 280
         f = f*dble(ib)
         c = c*b2
         go to 260
  280    g = r*dble(ib)
  300    if (c.lt.g) go to 320
         f = f/dble(ib)
         c = c/b2
         go to 300
  320    if (((c+r)/f).ge.(0.95d0*s)) go to 400
         g = 1.0d0/f
         d(i) = d(i)*f
         noconv = .true.
         if (l.gt.n) go to 360
         do 340 j = l, n
            a(i,j) = a(i,j)*g
  340    continue
  360    do 380 j = 1, k
            a(j,i) = a(j,i)*f
  380    continue
  400 continue
  420 if (noconv) go to 220
      return
      end
      subroutine f01auf(n,low,lhi,m,d,z,iz)
c     mark 2 release. nag copyright 1972
c     mark 4 revised.
c     mark 4.5 revised
c     mark 11.5(f77) revised. (sept 1985.)
c
c     balbak
c     backward transformation of a set of right-hand eigenvectors
c     of a balanced matrix into the eigenvectors of the original
c     matrix from which the balanced matrix was derived by a call
c     of subroutine f01atf.
c     december 1st.,1971
c
c     .. scalar arguments ..
      integer           iz, lhi, low, m, n
c     .. array arguments ..
      double precision  d(n), z(iz,m)
c     .. local scalars ..
      double precision  s
      integer           i, ii, j, k, lhi1, low1
c     .. executable statements ..
      if (low.gt.lhi) go to 60
      do 40 i = low, lhi
         s = d(i)
c        left-hand eigenvectors are back transformed if the
c        foregoing statement is replaced by s=1/d(i)
         do 20 j = 1, m
            z(i,j) = z(i,j)*s
   20    continue
   40 continue
   60 i = low
      low1 = low - 1
      if (low1.lt.1) go to 120
      do 100 ii = 1, low1
         i = i - 1
         k = d(i)
         if (k.eq.i) go to 100
         do 80 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
   80    continue
  100 continue
  120 lhi1 = lhi + 1
      if (lhi1.gt.n) return
      do 160 i = lhi1, n
         k = d(i)
         if (k.eq.i) go to 160
         do 140 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
  140    continue
  160 continue
      return
      end
      subroutine f02aqf(n,low,upp,machep,h,ih,vecs,ivecs,wr,wi,cnt,
     *                  ifail)
c     mark 13 re-issue. nag copyright 1988.
c
c     hqr2
c     finds the eigenvalues and eigenvectors of a real matrix
c     which has been reduced to upper hessenberg form in the array
c     h(n,n) with the accumulated transformations stored in
c     the array vecs(n,n). the real and imaginary parts of the
c     eigenvalues are formed in the arrays wr, wi(n) and the
c     eigenvectors are formed in the array vecs(n,n) where
c     only one complex vector, corresponding to the root with
c     positive imaginary part, is formed for a complex pair. low
c     and upp are two integers produced in balancing where
c     eigenvalues are isolated in positions 1 to low-1 and upp+1
c     to n. if balancing is not used low=1, upp=n. machep is the
c     relative machine precision. the subroutine will fail if
c     all eigenvalues take more than 30*n iterations.
c     1st december 1971
c
c     .. parameters ..
      character*6       srname
      parameter         (srname='f02aqf')
c     .. scalar arguments ..
      double precision  machep
      integer           ifail, ih, ivecs, low, n, upp
c     .. array arguments ..
      double precision  h(ih,n), vecs(ivecs,n), wi(n), wr(n)
      integer           cnt(n)
c     .. local scalars ..
      double precision  norm, p, q, r, ra, s, sa, t, u, vi, vr, w, x, y,
     *                  z
      integer           en, en2, i, i1, ii, isave, itn, its, j, jj, k,
     *                  kk, l, ll, low1, m, m2, m3, mm, na, na1, nhs,
     *                  upp1
      logical           notlas
c     .. local arrays ..
      character*1       p01rec(1)
c     .. external functions ..
      double precision  a02abf, x02ajf, x02akf
      integer           p01abf
      external          a02abf, x02ajf, x02akf, p01abf
c     .. external subroutines ..
      external          a02acf, dgemv
c     .. intrinsic functions ..
      intrinsic         abs, min, dble, sqrt
c     .. executable statements ..
      isave = ifail
c     compute matrix norm
      norm = 0.0d0
      k = 1
      do 40 i = 1, n
         do 20 j = k, n
            norm = norm + abs(h(i,j))
   20    continue
         k = i
   40 continue
      nhs = n*(n+1)/2 + n - 1
c     isolated roots
      if (low.le.1) go to 80
      j = low - 1
      do 60 i = 1, j
         wr(i) = h(i,i)
         wi(i) = 0.0d0
         cnt(i) = 0
   60 continue
   80 if (upp.ge.n) go to 120
      j = upp + 1
      do 100 i = j, n
         wr(i) = h(i,i)
         wi(i) = 0.0d0
         cnt(i) = 0
  100 continue
  120 en = upp
      t = 0.0d0
      itn = 30*n
  140 if (en.lt.low) go to 880
      its = 0
      na = en - 1
c     look for single small sub-diagonal element
  160 if (low+1.gt.en) go to 200
      low1 = low + 1
      do 180 ll = low1, en
         l = en + low1 - ll
         s = abs(h(l-1,l-1)) + abs(h(l,l))
         if (s.lt.x02akf()/x02ajf()) s = norm/dble(nhs)
         if (abs(h(l,l-1)).le.machep*s) go to 220
  180 continue
  200 l = low
  220 x = h(en,en)
      if (l.eq.en) go to 740
      y = h(na,na)
      w = h(en,na)*h(na,en)
      if (l.eq.na) go to 760
      if (itn.le.0) go to 1500
c     form shift
      if ((its.ne.10) .and. (its.ne.20)) go to 280
      t = t + x
      if (low.gt.en) go to 260
      do 240 i = low, en
         h(i,i) = h(i,i) - x
  240 continue
  260 s = abs(h(en,na)) + abs(h(na,en-2))
      x = 0.75d0*s
      y = x
      w = -0.4375d0*s**2
  280 its = its + 1
      itn = itn - 1
c     look for two consecutive small sub-diagonal elements
      if (l.gt.en-2) go to 320
      en2 = en - 2
      do 300 mm = l, en2
         m = l + en2 - mm
         z = h(m,m)
         r = x - z
         s = y - z
         p = (r*s-w)/h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - z - r - s
         r = h(m+2,m+1)
         s = abs(p) + abs(q) + abs(r)
         p = p/s
         q = q/s
         r = r/s
         if (m.eq.l) go to 320
         if ((abs(h(m,m-1))*(abs(q)+abs(r))).le.(machep*abs(p)
     *       *(abs(h(m-1,m-1))+abs(z)+abs(h(m+1,m+1))))) go to 320
  300 continue
  320 m2 = m + 2
      if (m2.gt.en) go to 360
      do 340 i = m2, en
         h(i,i-2) = 0.0d0
  340 continue
  360 m3 = m + 3
      if (m3.gt.en) go to 400
      do 380 i = m3, en
         h(i,i-3) = 0.0d0
  380 continue
  400 if (m.gt.na) go to 720
      do 700 k = m, na
         notlas = k .ne. na
         if (k.eq.m) go to 420
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0
         if (notlas) r = h(k+2,k-1)
         x = abs(p) + abs(q) + abs(r)
         if (x.eq.0.0d0) go to 700
         p = p/x
         q = q/x
         r = r/x
  420    s = sqrt(p**2+q**2+r**2)
         if (p.lt.0.0d0) s = -s
         if (k.ne.m) go to 440
         if (l.ne.m) h(k,k-1) = -h(k,k-1)
         go to 460
  440    h(k,k-1) = -s*x
  460    p = p + s
         x = p/s
         y = q/s
         z = r/s
         q = q/p
         r = r/p
c        row modification
         if (notlas) go to 500
         do 480 j = k, n
            p = h(k,j) + q*h(k+1,j)
            h(k+1,j) = h(k+1,j) - p*y
            h(k,j) = h(k,j) - p*x
  480    continue
         go to 540
  500    do 520 j = k, n
            p = h(k,j) + q*h(k+1,j) + r*h(k+2,j)
            h(k+2,j) = h(k+2,j) - p*z
            h(k+1,j) = h(k+1,j) - p*y
            h(k,j) = h(k,j) - p*x
  520    continue
  540    j = en
         if (k+3.lt.en) j = k + 3
c        column modification
         if (notlas) go to 580
         do 560 i = 1, j
            p = x*h(i,k) + y*h(i,k+1)
            h(i,k+1) = h(i,k+1) - p*q
            h(i,k) = h(i,k) - p
  560    continue
         go to 620
  580    do 600 i = 1, j
            p = x*h(i,k) + y*h(i,k+1) + z*h(i,k+2)
            h(i,k+2) = h(i,k+2) - p*r
            h(i,k+1) = h(i,k+1) - p*q
            h(i,k) = h(i,k) - p
  600    continue
c        accumulate transformations
  620    if (low.gt.upp) go to 700
         if (notlas) go to 660
         do 640 i = low, upp
            p = x*vecs(i,k) + y*vecs(i,k+1)
            vecs(i,k+1) = vecs(i,k+1) - p*q
            vecs(i,k) = vecs(i,k) - p
  640    continue
         go to 700
  660    do 680 i = low, upp
            p = x*vecs(i,k) + y*vecs(i,k+1) + z*vecs(i,k+2)
            vecs(i,k+2) = vecs(i,k+2) - p*r
            vecs(i,k+1) = vecs(i,k+1) - p*q
            vecs(i,k) = vecs(i,k) - p
  680    continue
  700 continue
  720 go to 160
c     one root found
  740 wr(en) = x + t
      h(en,en) = wr(en)
      wi(en) = 0.0d0
      cnt(en) = its
      en = na
      go to 140
c     two roots found
  760 p = (y-x)/2.0d0
      q = p**2 + w
      z = sqrt(abs(q))
      x = x + t
      h(en,en) = x
      h(na,na) = y + t
      cnt(en) = -its
      cnt(na) = its
      if (q.lt.0.0d0) go to 840
c     real pair
      if (p.lt.0.0d0) z = p - z
      if (p.gt.0.0d0) z = p + z
      wr(na) = x + z
      wr(en) = wr(na)
      if (z.ne.0.0d0) wr(en) = x - w/z
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      x = h(en,na)
      r = a02abf(x,z)
      p = x/r
      q = z/r
c     row modification
      do 780 j = na, n
         z = h(na,j)
         h(na,j) = q*z + p*h(en,j)
         h(en,j) = q*h(en,j) - p*z
  780 continue
c     column modification
      do 800 i = 1, en
         z = h(i,na)
         h(i,na) = q*z + p*h(i,en)
         h(i,en) = q*h(i,en) - p*z
  800 continue
c     accumulate transformations
      do 820 i = low, upp
         z = vecs(i,na)
         vecs(i,na) = q*z + p*vecs(i,en)
         vecs(i,en) = q*vecs(i,en) - p*z
  820 continue
      go to 860
c     complex pair
  840 wr(na) = x + p
      wr(en) = x + p
      wi(na) = z
      wi(en) = -z
  860 en = en - 2
      go to 140
c     all roots found now backsubstitute
  880 if (norm.eq.0.0d0) go to 1480
      norm = norm*machep
c     backsubstitution
      do 1340 kk = 1, n
         en = n + 1 - kk
         p = wr(en)
         q = wi(en)
         na = en - 1
         if (q.ne.0.0d0) go to 1120
c        real vector
         h(en,en) = 1.0d0
         if (na.lt.1) go to 1340
         do 1100 ii = 1, na
            i = na + 1 - ii
            i1 = i - 1
            w = h(i,i) - p
            r = h(i,en)
            if (wi(i).ge.0.0d0) go to 900
            z = w
            s = r
            go to 1100
  900       if (wi(i).gt.0.0d0) go to 1020
c           modification to stop overflow
            if (w.ne.0.0d0) go to 940
            if (abs(r).lt.10.0d0*norm) go to 960
            r = -r
            do 920 j = 1, en
               h(j,en) = h(j,en)*norm
  920       continue
            go to 980
  940       r = -r/w
            go to 980
  960       r = -r/norm
  980       h(i,en) = r
            if (i1.eq.0) go to 1100
            do 1000 j = 1, i1
               h(j,en) = h(j,en) + h(j,i)*r
 1000       continue
            go to 1100
c           solve real equations
 1020       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i)-p)**2 + wi(i)**2
            t = (x*s-z*r)/q
            h(i,en) = t
            if (abs(x).gt.abs(z)) go to 1040
            r = (-s-y*t)/z
            go to 1060
 1040       r = (-r-w*t)/x
 1060       h(i+1,en) = r
            if (i1.eq.0) go to 1100
            do 1080 j = 1, i1
               h(j,en) = (h(j,en)+h(j,i+1)*r) + h(j,i)*t
 1080       continue
 1100    continue
c        end real vector
         go to 1340
 1120    if (q.gt.0.0d0) go to 1340
c        complex vector associated with lambda=p-i*q
         if (abs(h(en,na)).le.abs(h(na,en))) go to 1140
         r = q/h(en,na)
         s = -(h(en,en)-p)/h(en,na)
         go to 1160
 1140    call a02acf(0.0d0,-h(na,en),h(na,na)-p,q,r,s)
 1160    h(en,na) = 0.0d0
         h(en,en) = 1.0d0
         h(na,na) = r
         h(na,en) = s
         if (na.lt.2) go to 1340
         na1 = na - 1
         do 1180 j = 1, na1
            h(j,en) = h(j,en) + h(j,na)*s
            h(j,na) = h(j,na)*r
 1180    continue
         do 1320 ii = 1, na1
            i = 1 + na1 - ii
            i1 = i - 1
            w = h(i,i) - p
            ra = h(i,na)
            sa = h(i,en)
            if (wi(i).ge.0.0d0) go to 1200
            z = w
            r = ra
            s = sa
            go to 1320
 1200       if (wi(i).eq.0.0d0) go to 1280
c           solve complex equations
            x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i)-p)**2 + wi(i)**2 - q**2
            vi = (wr(i)-p)*2.0d0*q
            if ((vr.eq.0.0d0) .and. (vi.eq.0.0d0))
     *          vr = machep*norm*(abs(w)+abs(q)+abs(x)+abs(y)+abs(z))
            call a02acf(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi,t,u)
            if (abs(x).le.abs(z)+abs(q)) go to 1220
            r = (-ra-w*t+q*u)/x
            s = (-sa-w*u-q*t)/x
            go to 1240
 1220       call a02acf(-r-y*t,-s-y*u,z,q,r,s)
 1240       h(i,na) = t
            h(i,en) = u
            h(i+1,na) = r
            h(i+1,en) = s
            if (i1.eq.0) go to 1320
            do 1260 j = 1, i1
               h(j,na) = (h(j,na)+h(j,i+1)*r) + h(j,i)*t
               h(j,en) = (h(j,en)+h(j,i+1)*s) + h(j,i)*u
 1260       continue
            go to 1320
 1280       call a02acf(-ra,-sa,w,q,r,s)
            h(i,na) = r
            h(i,en) = s
            if (i1.eq.0) go to 1320
            do 1300 j = 1, i1
               h(j,na) = h(j,na) + h(j,i)*r
               h(j,en) = h(j,en) + h(j,i)*s
 1300       continue
 1320    continue
c        end complex vector
 1340 continue
c     end backsubstitution
c     vectors of isolated roots
      low1 = low - 1
      upp1 = upp + 1
      do 1420 j = 1, n
         m = min(j,low1)
         if (m.lt.1) go to 1380
         do 1360 i = 1, m
            vecs(i,j) = h(i,j)
 1360    continue
 1380    if (upp1.gt.j) go to 1420
         do 1400 i = upp1, j
            vecs(i,j) = h(i,j)
 1400    continue
 1420 continue
c     multiply by transformation matrix to give
c     vectors of original full matrix
      do 1460 jj = low, n
         j = low + n - jj
         m = min(j,upp)
         do 1440 i = low, upp
            vecs(i,j) = vecs(i,m)*h(m,j)
 1440    continue
         m = m - 1
         if (m+1.ge.low) call dgemv('n',upp-low+1,m-low+1,1.0d0,
     *                              vecs(low,low),ivecs,h(low,j),1,
     *                              1.0d0,vecs(low,j),1)
 1460 continue
 1480 ifail = 0
      return
 1500 ifail = p01abf(isave,1,srname,0,p01rec)
      return
      end
      subroutine f01atz(m,a,ia,d,k,l,n,j)
c     mark 2 release. nag copyright 1972
c     mark 4 revised.
c     mark 11.5(f77) revised. (sept 1985.)
c     nag copyright 1975
c     mark 4.5 revised
c
c     auxiliary routine called by f01atf.
c     interchanges elements 1 to k of columns j and m,
c     and elements l to n of rows j and m.
c     .. scalar arguments ..
      integer           ia, j, k, l, m, n
c     .. array arguments ..
      double precision  a(ia,n), d(n)
c     .. local scalars ..
      double precision  f
      integer           i
c     .. executable statements ..
      d(m) = j
      if (j.eq.m) go to 60
      do 20 i = 1, k
         f = a(i,j)
         a(i,j) = a(i,m)
         a(i,m) = f
   20 continue
      if (l.gt.n) go to 60
      do 40 i = l, n
         f = a(j,i)
         a(j,i) = a(m,i)
         a(m,i) = f
   40 continue
   60 return
      end
      subroutine f06aaz ( srname, info )
c     mark 12 release. nag copyright 1986.
c     .. scalar arguments ..
      integer            info
      character*13       srname
c     ..
c
c  purpose
c  =======
c
c  f06aaz  is an error handler for the level 2 blas routines.
c
c  it is called by the level 2 blas routines if an input parameter is
c  invalid.
c
c  parameters
c  ==========
c
c  srname - character*13.
c           on entry, srname specifies the name of the routine which
c           called f06aaz.
c
c  info   - integer.
c           on entry, info specifies the position of the invalid
c           parameter in the parameter-list of the calling routine.
c
c
c  auxiliary routine for level 2 blas.
c
c  written on 20-july-1986.
c
c     .. local scalars ..
      integer            ifail
      character*80       rec (1)
c     .. external functions ..
      integer            p01abf
      external           p01abf
c     ..
c     .. executable statements ..
      write (rec (1),99999) srname, info
      ifail = 0
      ifail = p01abf (ifail, -1, srname(1:6), 1, rec)
c
      return
c
99999 format ( ' ** on entry to ', a13, ' parameter number ', i2,
     $         ' had an illegal value' )
c
c     end of f06aaz.
c
      end
      integer function p01abf(ifail,ierror,srname,nrec,rec)
c     mark 11.5(f77) release. nag copyright 1986.
c     mark 13 revised. ier-621 (apr 1988).
c     mark 13b revised. ier-668 (aug 1988).
c
c     p01abf is the error-handling routine for the nag library.
c
c     p01abf either returns the value of ierror through the routine
c     name (soft failure), or terminates execution of the program
c     (hard failure). diagnostic messages may be output.
c
c     if ierror = 0 (successful exit from the calling routine),
c     the value 0 is returned through the routine name, and no
c     message is output
c
c     if ierror is non-zero (abnormal exit from the calling routine),
c     the action taken depends on the value of ifail.
c
c     ifail =  1: soft failure, silent exit (i.e. no messages are 
c                 output)
c     ifail = -1: soft failure, noisy exit (i.e. messages are output)
c     ifail =-13: soft failure, noisy exit but standard messages from
c                 p01abf are suppressed
c     ifail =  0: hard failure, noisy exit
c
c     for compatibility with certain routines included before mark 12
c     p01abf also allows an alternative specification of ifail in which
c     it is regarded as a decimal integer with least significant digits
c     cba. then
c
c     a = 0: hard failure  a = 1: soft failure
c     b = 0: silent exit   b = 1: noisy exit
c
c     except that hard failure now always implies a noisy exit.
c
c     s.hammarling, m.p.hooper and j.j.du croz, nag central office.
c
c     .. scalar arguments ..
      integer                 ierror, ifail, nrec
      character*(*)           srname
c     .. array arguments ..
      character*(*)           rec(*)
c     .. local scalars ..
      integer                 i, nerr
      character*72            mess
c     .. external subroutines ..
      external                p01abz, x04aaf, x04baf
c     .. intrinsic functions ..
      intrinsic               abs, mod
c     .. executable statements ..
      if (ierror.ne.0) then
c        abnormal exit from calling routine
         if (ifail.eq.-1 .or. ifail.eq.0 .or. ifail.eq.-13 .or.
     *       (ifail.gt.0 .and. mod(ifail/10,10).ne.0)) then
c           noisy exit
            call x04aaf(0,nerr)
            do 20 i = 1, nrec
               call x04baf(nerr,rec(i))
   20       continue
            if (ifail.ne.-13) then
               write (mess,fmt=99999) srname, ierror
               call x04baf(nerr,mess)
               if (abs(mod(ifail,10)).ne.1) then
c                 hard failure
                  call x04baf(nerr,
     *                     ' ** nag hard failure - execution terminated'
     *                        )
                  call p01abz
               else
c                 soft failure
                  call x04baf(nerr,
     *                        ' ** nag soft failure - control returned')
               end if
            end if
         end if
      end if
      p01abf = ierror
      return
c
99999 format (' ** abnormal exit from nag library routine ',a,': ifail',
     *  ' =',i6)
      end
      subroutine p01abz
c     mark 11.5(f77) release. nag copyright 1986.
c
c     terminates execution when a hard failure occurs.
c
c     ******************** implementation note ********************
c     the following stop statement may be replaced by a call to an
c     implementation-dependent routine to display a message and/or
c     to abort the program.
c     *************************************************************
c     .. executable statements ..
      stop
      end
      double precision function a02abf(xxr,xxi)
c     nag copyright 1975
c     mark 4.5 revised
c     mark 5c revised
c     mark 11.5(f77) revised. (sept 1985.)
c
c     returns the absolute value of a complex number via routine
c     name
c
c     .. scalar arguments ..
      double precision                 xxi, xxr
c     .. local scalars ..
      double precision                 h, one, xi, xr, zero
c     .. intrinsic functions ..
      intrinsic                        abs, sqrt
c     .. data statements ..
      data                             zero/0.0d0/, one/1.0d0/
c     .. executable statements ..
c
      xr = abs(xxr)
      xi = abs(xxi)
      if (xi.le.xr) go to 20
      h = xr
      xr = xi
      xi = h
   20 if (xi.ne.zero) go to 40
      a02abf = xr
      return
   40 h = xr*sqrt(one+(xi/xr)**2)
      a02abf = h
      return
      end
      subroutine a02acf(xxr,xxi,yyr,yyi,zr,zi)
c     mark 2a release.  nag copyright 1973
c     mark 4.5 revised
c     mark 5c revised
c     mark 11.5(f77) revised. (sept 1985.)
c
c     divides one complex number by a second
c
c     .. scalar arguments ..
      double precision  xxi, xxr, yyi, yyr, zi, zr
c     .. local scalars ..
      double precision  a, h, one
c     .. intrinsic functions ..
      intrinsic         abs
c     .. data statements ..
      data              one/1.0d0/
c     .. executable statements ..
c
      if (abs(yyr).le.abs(yyi)) go to 20
      h = yyi/yyr
      a = one/(h*yyi+yyr)
      zr = (xxr+h*xxi)*a
      zi = (xxi-h*xxr)*a
      return
   20 h = yyr/yyi
      a = one/(h*yyr+yyi)
      zr = (h*xxr+xxi)*a
      zi = (h*xxi-xxr)*a
      return
      end
      double precision function x02ajf()
c     mark 12 release. nag copyright 1986.
c
c     returns  (1/2)*b**(1-p)  if rounds is .true.
c     returns  b**(1-p)  otherwise
c
      double precision x02con
!      data x02con /z'3ca0000000000001' /
c     .. executable statements ..
!      x02ajf = x02con
      x02ajf = epsilon(x02con)
      return
      end
      integer function x02bhf()
c     mark 12 release. nag copyright 1986.
c
c     returns the model parameter, b.
c
c     .. executable statements ..
      x02bhf =     2
      x02bhf =     radix(x02bhf)
      return
      end
      double precision function x02akf()
c     mark 12 release. nag copyright 1986.
c
c     returns  b**(emin-1)  (the smallest positive model number)
c
      double precision x02con
!      data x02con /z'0010000000000000' /
c     .. executable statements ..
!      x02akf = x02con
      x02akf = tiny(x02con)
      return
      end
      subroutine x04aaf(i,nerr)
c     mark 7 release. nag copyright 1978
c     mark 7c revised ier-190 (may 1979)
c     mark 11.5(f77) revised. (sept 1985.)
c     if i = 0, sets nerr to current error message unit number
c     (stored in nerr1).
c     if i = 1, changes current error message unit number to
c     value specified by nerr.
c
c     .. scalar arguments ..
      integer           i, nerr
c     .. local scalars ..
      integer           nerr1
c     .. save statement ..
      save              nerr1
c     .. data statements ..
      data              nerr1/0/
c     .. executable statements ..
      if (i.eq.0) nerr = nerr1
      if (i.eq.1) nerr1 = nerr
      return
      end
      subroutine x04baf(nout,rec)
c     mark 11.5(f77) release. nag copyright 1986.
c
c     x04baf writes the contents of rec to the unit defined by nout.
c
c     trailing blanks are not output, except that if rec is entirely
c     blank, a single blank character is output.
c     if nout.lt.0, i.e. if nout is not a valid fortran unit identifier,
c     then no output occurs.
c
c     .. scalar arguments ..
      integer           nout
      character*(*)     rec
c     .. local scalars ..
      integer           i
c     .. intrinsic functions ..
      intrinsic         len
c     .. executable statements ..
      if (nout.ge.0) then
c        remove trailing blanks
         do 20 i = len(rec), 2, -1
            if (rec(i:i).ne.' ') go to 40
   20    continue
c        write record to external file
   40    write (nout,fmt=99999) rec(1:i)
      end if
      return
c
99999 format (a)
      end
