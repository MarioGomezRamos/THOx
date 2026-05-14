import re

with open("belam.f90", "r") as f:
    content = f.read()

# Discrete block
old_discrete = """c SIMULTANEOUS Core and Valence contribution --------------------------
      if (lambda.gt.1) then
        do ikterm = 1, lambda - 1
           rkt = dble(ikterm)
           rlv = dble(lambda - ikterm)
           if (fail3(qjci(n), rkt, qjc(m))) cycle
           if (fail3(qji(n), rlv, qj(m))) cycle
           if (fail3(qlir(n), rlv, qlr(m))) cycle
           select case(coremodel)
           case(0) 
             if (rms.ne.0d0) then
               rkcorem = 3d0*zc*delta*(sqrt(5d0/3d0)*rms)**(ikterm-1)/
      &                   (4d0*pi)
               rkcore  = rkcorem*cleb(qjci(n), 0d0, rkt, 0d0, qjc(m), 
      &                   0d0) * sqrt((2*qjci(n)+1)/(2*qjc(m)+1)) * 
      &                   (-1)**(2*ikterm)
             else
               rkcore = 0d0
             endif
           case(1)
             rkcore = mec(cindexi(n), cindex(m), ikterm)
           end select
           if (abs(rkcore).lt.1d-10) cycle
           faux(1:nr) = rvec(1:nr)**(1 + lambda - ikterm) 
      &                * ugs(1:nr,n) * wfc(jset,i,m,1:nr)
           call simc(faux, resc, 1, nr, dr, nr)
           rlam = resc
           term_pref = sqrt(4d0 * pi) * dfactratio(lambda, ikterm) 
      &                 * ((-av)/(ac+av))**(lambda - ikterm)
           term_mix = matmix(rkt, rlv, lambdar, sn, 
      &                       qlir(n), qji(n), jtoti, qjci(n),
      &                       qlr(m), qj(m), jtot, qjc(m))
           Elamcore = term_pref * term_mix * rlam * rkcore
           Elam = Elam + Elamcore
           mvc  = mvc + Elamcore
           mel(i,m) = mel(i,m) + cmplx(Elamcore, 0d0)
        enddo
      endif
c ----------------------------------------------------------------------------"""

new_discrete = """c SIMULTANEOUS Core and Valence contribution --------------------------
        do ikterm = 1, maxlamb
           do ilv = 1, maxlamb
           rkt = dble(ikterm)
           rlv = dble(ilv)
           if (fail3(rkt, rlv, lambdar)) cycle
           if (fail3(qjci(n), rkt, qjc(m))) cycle
           if (fail3(qji(n), rlv, qj(m))) cycle
           if (fail3(qlir(n), rlv, qlr(m))) cycle
           select case(coremodel)
           case(0) 
             if (rms.ne.0d0) then
               rkcorem = 3d0*zc*delta*(sqrt(5d0/3d0)*rms)**(ikterm-1)/
      &                   (4d0*pi)
               rkcore  = rkcorem*cleb(qjci(n), 0d0, rkt, 0d0, qjc(m), 
      &                   0d0) * sqrt((2*qjci(n)+1)/(2*qjc(m)+1)) * 
      &                   (-1)**(2*ikterm)
             else
               rkcore = 0d0
             endif
           case(1)
             rkcore = mec(cindexi(n), cindex(m), ikterm)
           end select
           if (abs(rkcore).lt.1d-10) cycle
           faux(1:nr) = rvec(1:nr)**(1 + ilv) 
      &                * ugs(1:nr,n) * wfc(jset,i,m,1:nr)
           call simc(faux, resc, 1, nr, dr, nr)
           rlam = resc
           if (ikterm + ilv == lambda) then
             term_pref = sqrt(4d0 * pi) * dfactratio(lambda, ikterm) 
      &                 * ((-av)/(ac+av))**(ilv)
           else
             term_pref = 1.0d0
           endif
           term_mix = matmix(rkt, rlv, lambdar, sn, 
      &                       qlir(n), qji(n), jtoti, qjci(n),
      &                       qlr(m), qj(m), jtot, qjc(m))
           Elamcore = term_pref * term_mix * rlam * rkcore
           Elam = Elam + Elamcore
           mvc  = mvc + Elamcore
           mel(i,m) = mel(i,m) + cmplx(Elamcore, 0d0)
           enddo
        enddo
c ----------------------------------------------------------------------------"""

old_cont = """c SIMULTANEOUS Core and Valence contribution --------------------------
      if (lambda.gt.1) then
        do ikterm = 1, lambda - 1
           rkt = dble(ikterm)
           rlv = dble(lambda - ikterm)
           if (fail3(qjci(n), rkt, qjc(m))) cycle
           if (fail3(qji(n), rlv, qj(m))) cycle
           if (fail3(qlir(n), rlv, qlr(m))) cycle
           select case(coremodel)
           case(0) 
             if (rms.ne.0d0) then
               rkcorem = 3d0*zc*delta*(sqrt(5d0/3d0)*rms)**(ikterm-1)/
      &                   (4d0*pi)
               rkcore  = rkcorem*cleb(qjci(n), 0d0, rkt, 0d0, qjc(m), 
      &                   0d0) * sqrt((2*qjci(n)+1)/(2*qjc(m)+1)) * 
      &                   (-1)**(2*ikterm)
             else
               rkcore = 0d0
             endif
           case(1)
             rkcore = mec(cindexi(n), cindex(m), ikterm)
           end select
           if (abs(rkcore).lt.1d-10) cycle
           gaux(1:nr) = rvec(1:nr)**(lambda - ikterm) 
      &                * ugs(1:nr,n) * wfscat(ik,m,1:nr)
           call simc(gaux, resc, 1, nr, dr, nr)
           term_pref = sqrt(4d0 * pi) * dfactratio(lambda, ikterm) 
      &                 * ((-av)/(ac+av))**(lambda - ikterm)
           term_mix = matmix(rkt, rlv, lambdar, sn, 
      &                       qlir(n), qji(n), jtoti, qjci(n),
      &                       qlr(m), qj(m), jtot, qjc(m))
           write(0,*)'term_mix=',term_mix
           Elamcont = Elamcont + term_pref * term_mix * rkcore * resc 
      &                * sqrt(pi/2) * (4*pi) / kcont
        enddo
      endif
c ----------------------------------------------------------------------------"""

new_cont = """c SIMULTANEOUS Core and Valence contribution --------------------------
        do ikterm = 1, maxlamb
           do ilv = 1, maxlamb
           rkt = dble(ikterm)
           rlv = dble(ilv)
           if (fail3(rkt, rlv, lambdar)) cycle
           if (fail3(qjci(n), rkt, qjc(m))) cycle
           if (fail3(qji(n), rlv, qj(m))) cycle
           if (fail3(qlir(n), rlv, qlr(m))) cycle
           select case(coremodel)
           case(0) 
             if (rms.ne.0d0) then
               rkcorem = 3d0*zc*delta*(sqrt(5d0/3d0)*rms)**(ikterm-1)/
      &                   (4d0*pi)
               rkcore  = rkcorem*cleb(qjci(n), 0d0, rkt, 0d0, qjc(m), 
      &                   0d0) * sqrt((2*qjci(n)+1)/(2*qjc(m)+1)) * 
      &                   (-1)**(2*ikterm)
             else
               rkcore = 0d0
             endif
           case(1)
             rkcore = mec(cindexi(n), cindex(m), ikterm)
           end select
           if (abs(rkcore).lt.1d-10) cycle
           gaux(1:nr) = rvec(1:nr)**(ilv) 
      &                * ugs(1:nr,n) * wfscat(ik,m,1:nr)
           call simc(gaux, resc, 1, nr, dr, nr)
           if (ikterm + ilv == lambda) then
             term_pref = sqrt(4d0 * pi) * dfactratio(lambda, ikterm) 
      &                 * ((-av)/(ac+av))**(ilv)
           else
             term_pref = 1.0d0
           endif
           term_mix = matmix(rkt, rlv, lambdar, sn, 
      &                       qlir(n), qji(n), jtoti, qjci(n),
      &                       qlr(m), qj(m), jtot, qjc(m))
           write(0,*)'term_mix=',term_mix
           Elamcont = Elamcont + term_pref * term_mix * rkcore * resc 
      &                * sqrt(pi/2) * (4*pi) / kcont
           enddo
        enddo
c ----------------------------------------------------------------------------"""

if old_discrete in content:
    content = content.replace(old_discrete, new_discrete)
    print("Replaced discrete block")
else:
    print("Failed to find discrete block")

if old_cont in content:
    content = content.replace(old_cont, new_cont)
    print("Replaced continuum block")
else:
    print("Failed to find continuum block")

with open("belam.f90", "w") as f:
    f.write(content)
