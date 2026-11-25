      program trip2lab
c     -----------------------------------------------------------------
c     trip2lab.f
c
c     Transforms a triple-differential breakup cross section
c     d(5)sigma/(dEc dOmega_c dOmega_v)  [tabulated, in the CM frame of
c     a+A, as a function of (Ec,theta_c,theta_v,phi_v), with c="core"
c     particle b, v="valence" particle x, a=b+x -- exactly the format
c     written/read by husk-phi.f] for the reaction a+A -> b+x+A, into
c     the fully exclusive LAB-frame differential cross section
c
c         d(5)sigma / (dEb_lab dOmega_b_lab dEx_lab dOmega_x_lab)
c
c     for a chosen incoming LAB beam: reverse=.true. (default) means
c     A (the target of the original normal-kinematics calculation) is
c     the incoming beam here (inverse kinematics, matching an actual
c     heavy-beam experiment); reverse=.false. means a (=b+x, the
c     original projectile) is the incoming beam instead (i.e. the
c     lab frame coincides with the frame the original calculation
c     itself used). elabbeam (&control) is that beam's LAB kinetic
c     energy, expressed as a function of
c         (theta_b_lab, Eb_lab, theta_x_lab, Ex_lab, Dphi_lab)
c     where Dphi_lab is the azimuthal angle of x relative to b in the
c     lab frame.
c
c     ---------------------------------------------------------------
c     Method: stratified (deterministic Monte-Carlo-like) sampling,
c     the same technique used for the two-body cm2lab_am.f code, now
c     applied independently to each of the two detected fragments.
c
c     Key simplification used here: the CM->LAB boost is a pure
c     Galilean boost along the beam (z) axis (vcm is a fixed vector,
c     the same for every fragment). Such a boost changes only the
c     LONGITUDINAL momentum component of a particle; the TRANSVERSE
c     components -- and therefore the particle's own azimuthal angle
c     -- are exactly unchanged. Since phi_c=0 by the convention used
c     to tabulate the CM file (phi_v is measured relative to the
c     a-b/beam-core plane), this means:
c
c         Dphi_lab = phi_v   (identically, no transformation needed)
c
c     and the two fragments' (theta,E) lab kinematics transform
c     completely independently of one another (b's boost depends only
c     on b's own CM (Ec,theta_c); x's only on x's own CM
c     (Ev,theta_v)) -- even though their weight must still be
c     deposited jointly, since Ec and Ev are correlated per input
c     grid cell (both stored directly per cell, exactly as clarified:
c     the file already gives the fully energy-momentum-consistent
c     (Ec,Ev) pair for each (ic,iv,ien,ip), so no conservation
c     equation needs to be solved and no branch ambiguity arises).
c
c     Approximation made: within one input grid cell, Ec and Ev are
c     held fixed at that cell's tabulated values while theta_c,
c     theta_v and phi_v are sub-sampled continuously across the cell.
c     This is the same pragmatic choice used before (Ecm held at bin
c     center while only angles are sub-sampled) -- adequate as long
c     as the input grid is reasonably fine, and easy to refine later
c     (interpolating Ec/Ev across neighbouring energy-grid points) if
c     found to be insufficient.
c
c     Normalization convention (matching husk-phi.f's own xstot check):
c     the tabulated raux is d(5)sigma/(dEc dtheta_c dtheta_v dphi_v)
c     (angles in RADIANS, not yet weighted by sin(theta)), and the
c     total cross section is obtained via
c       sigma_tot = 2*integral{ raux * 2*pi*sin(tc)*dtc
c                                *sin(tv)*dtv*dEc*dphi_v }
c     where the single 2*pi accounts for the trivial (unobserved,
c     because irrelevant by axial symmetry) overall rotation of the
c     whole event about the beam axis (i.e. the azimuth of particle c
c     itself), and the overall factor 2 in husk-phi.f's own xstot
c     line reflects that only phi_v in [0,pi] is tabulated (by
c     up-down/mirror symmetry about the reaction plane) -- see the
c     ISYMPHI note below for how that is handled here.
c
c     Output: a flat, sparse table (only non-negligible bins) of
c         theta_b_lab  Eb_lab  theta_x_lab  Ex_lab  Dphi_lab  sigma
c     with sigma = d(5)sigma/(dEb_lab dOmega_b_lab dEx_lab dOmega_x_lab)
c     in the same convention as the input (i.e. dOmega_b_lab means
c     2*pi*sin(theta_b_lab)*dtheta_b_lab -- trivial phi_b_lab
c     integrated out -- while dOmega_x_lab means simply
c     sin(theta_x_lab)*dtheta_x_lab*dDphi_lab, matching the input's
c     phi_v convention).
c     -----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 mt,zt,mc,mv,ep,ebind,mthird,mbeam
      character*80 filecm,sys*3,line
      logical isymphi,reverse
      integer pair

      real*8, allocatable:: sig(:,:,:,:),ecarr(:,:,:,:),evarr(:,:,:,:)
      real*8, allocatable:: hist(:,:,:,:,:)
      real*8, allocatable:: hist3th(:),hist3e(:)

c     ---------------------------------------------------------------
c     &control namelist (read from trip2lab.in) -- input variables:
c
c       filecm    - path to the CM triple-differential file, in
c                   husk-phi.f's format (masses, idet, grids, then
c                   the (raux,ec,ev) data, per (theta_c,theta_v,
c                   E,phi_v) grid cell)
c       elabbeam  - LAB kinetic energy (MeV) of whichever particle
c                   is chosen as the incoming beam via "reverse"
c       reverse   - .true. (default): A (the target of the original
c                   normal-kinematics calculation) is the incoming
c                   LAB beam (inverse kinematics)
c                   .false.: a (=b+x, the original projectile) is
c                   the incoming LAB beam instead (the same particle
c                   roles as the original calculation)
c       pair      - which two particles the output angle/energy
c                   variables (theta_x_lab,Ex_lab in the code, i.e.
c                   the arrays/limits named thxlab*/exlab*) refer to:
c                   1 (default): particles b and x
c                   2: particles b and A (A's kinematics are then
c                   obtained from momentum conservation, p_A =
c                   -(p_b+p_x), so they depend on theta_c, theta_v
c                   AND phi_v jointly -- unlike x's own kinematics,
c                   which depend only on theta_v)
c
c       thblabmin,thblabmax,thblabinc - output grid (deg) for b's own
c                   LAB polar angle, theta_b_lab
c       eblabmin,eblabmax,eblabinc   - output grid (MeV) for b's own
c                   LAB kinetic energy, Eb_lab
c       thxlabmin,thxlabmax,thxlabinc - output grid (deg) for the
c                   second particle's (x or A, per "pair") LAB polar
c                   angle
c       exlabmin,exlabmax,exlabinc  - output grid (MeV) for the
c                   second particle's LAB kinetic energy -- note the
c                   very different physical scale needed for x
c                   (tens of MeV) vs A (close to the full beam
c                   energy, since the heavy beam particle keeps most
c                   of its kinetic energy through the breakup)
c       dphlabmin,dphlabmax,dphlabinc - output grid (deg) for
c                   Dphi_lab, the second particle's azimuth relative
c                   to b's (b's own azimuth is 0 by convention)
c
c       th3labmin,th3labmax,th3labinc - output grid (deg) for the
c                   "third" particle's (A if pair=1, x if pair=2)
c                   LAB polar angle. This particle is NOT part of
c                   the main 5-D table (it is fully redundant given
c                   the other two particles' kinematics plus
c                   momentum conservation) but is written as two
c                   extra marginal files (angle-only and energy-
c                   only), derived at each sub-sample via lab-frame
c                   momentum conservation, for convenience -- see
c                   trip2lab_dsdth3.out / trip2lab_dsdE3.out below
c       e3labmin,e3labmax,e3labinc - output grid (MeV) for the third
c                   particle's LAB kinetic energy (same physical
c                   scale caveat as exlabmin/exlabmax above)
c
c       nsubc,nsubv,nsubphi - number of stratified sub-samples per
c                   input CM cell, in theta_c, theta_v and phi_v
c                   respectively (see comments below); more sub-
c                   samples cost more runtime but resolve the LAB
c                   mapping more accurately -- check convergence by
c                   increasing them and confirming the total cross
c                   section (printed at the end) is stable
c       isymphi   - .true. (default) if the CM file only tabulates
c                   phi_v in [0,180] (up-down mirror symmetry about
c                   the reaction plane, as in husk-phi.f); the
c                   [180,360) half is then reconstructed by mirroring
c       wtol      - input CM cells with a total weight (mb) at or
c                   below this are skipped entirely (0.d0 default:
c                   no cells skipped other than exact zeros)
c     ---------------------------------------------------------------
      namelist /control/ filecm,elabbeam,reverse,pair,
     x   thblabmin,thblabmax,thblabinc,
     x   eblabmin,eblabmax,eblabinc,
     x   thxlabmin,thxlabmax,thxlabinc,
     x   exlabmin,exlabmax,exlabinc,
     x   dphlabmin,dphlabmax,dphlabinc,
     x   th3labmin,th3labmax,th3labinc,
     x   e3labmin,e3labmax,e3labinc,
     x   nsubc,nsubv,nsubphi,isymphi,wtol

      pi=4.d0*datan(1.d0)
      degrad=pi/180.d0
      amu=938.0d0

c     ---------------------------------------------------------------
c     defaults, then read the run's own control namelist
c     ---------------------------------------------------------------
      nsubc=4
      nsubv=4
      nsubphi=3
      isymphi=.true.
      wtol=0.d0
      reverse=.true.
      pair=1
      open(unit=10,file='trip2lab.in',status='old')
      read(10,nml=control)
      close(10)
      open(unit=99,file='trip2lab.log',status='unknown')

c     ---------------------------------------------------------------
c     read the CM triple-differential file (husk-phi.f format)
c     ---------------------------------------------------------------
      open(77,file=filecm,status='old')
      read(77,'(a)') sys
      write(*,*) ' coordinate system (as stored): ',sys
      write(99,*) ' coordinate system (as stored): ',sys
      read(77,304) mt,zt,eoa
  304 format(5f10.3)
      write(*,*) ' target(A) mass, charge, projectile(a) E/A: ',
     x  mt,zt,eoa
      write(99,*) ' target(A) mass, charge, projectile(a) E/A: ',
     x  mt,zt,eoa
      read(77,304) mc
      write(*,*) ' core (b) mass: ',mc
      write(99,*) ' core (b) mass: ',mc
      read(77,304) mv,ebind
      write(*,*) ' valence (x) mass, separation energy: ',mv,ebind
      write(99,*) ' valence (x) mass, separation energy: ',mv,ebind
      read(77,304) ep
      write(*,*) ' projectile (a) LAB energy (normal kinematics): ',ep
      write(99,*) ' projectile (a) LAB energy (normal kinematics): ',ep
      read(77,*) idet
      if(idet.eq.1) then
       write(*,*) ' idet=1: tabulated energy is Ec (core/b)'
       write(99,*) ' idet=1: tabulated energy is Ec (core/b)'
      else
       write(*,*) ' idet=2: tabulated energy is Ev (valence/x)'
       write(99,*) ' idet=2: tabulated energy is Ev (valence/x)'
      endif

      read(77,307) itc,dtcr,tcl,tcu,dtc
      read(77,307) itv,dtvr,tvl,tvu,dtv
      read(77,307) iten,den,enlow,enup
      read(77,307) itp,dphr,phl,phu,dph
  307 format(i5,f15.8,3f10.3)
      write(*,*) ' grid dims (theta_c,theta_v,E,phi_v): ',
     x  itc,itv,iten,itp
      write(99,*) ' grid dims (theta_c,theta_v,E,phi_v): ',
     x  itc,itv,iten,itp

      allocate(sig(itc,itv,iten,itp))
      allocate(ecarr(itc,itv,iten,itp))
      allocate(evarr(itc,itv,iten,itp))

      xstot=0.d0
      do ic=1,itc
       tc=tcl+(ic-1)*dtc
       tcr=tc*degrad
       do iv=1,itv
        tv=tvl+(iv-1)*dtv
        tvr=tv*degrad
        do ien=1,iten
         do ip=1,itp
          phi=phl+(ip-1)*dph
          read(77,*) i1,i2,i3,i4,raux,ecv,evv
          sig(ic,iv,ien,ip)=raux
          ecarr(ic,iv,ien,ip)=ecv
          evarr(ic,iv,ien,ip)=evv
          term=raux*2.d0*pi*sin(tcr)*dtc*degrad
     x                *sin(tvr)*dtv*degrad*den*dph*degrad
c         proper mirror-symmetric phi weight: phi=0/180 are their own
c         mirror image under the reaction-plane reflection and must
c         be counted once, not twice, when isymphi folds [0,180] out
c         to the full [0,360) range (husk-phi.f's own xstot check
c         instead doubles everything uniformly, including these
c         boundary points, which slightly overestimates the total
c         whenever the boundary phi values carry non-negligible
c         cross section -- fixed here for a meaningful comparison)
          phiwt=2.d0
          if(isymphi.and.(phi.eq.0.d0.or.phi.eq.180.d0)) phiwt=1.d0
          if(.not.isymphi) phiwt=1.d0
          xstot=xstot + term*phiwt
         enddo
        enddo
       enddo
      enddo
      close(77)
      write(*,*) ' Total CM cross section (input, check) = ',xstot,
     x  ' mb'
      write(99,*) ' Total CM cross section (input, check) = ',xstot,
     x  ' mb'

      ma=mc+mv
      if(reverse) then
c      A (target of the original normal-kinematics calculation) is
c      the incoming LAB beam here (inverse kinematics). Since
c      theta_c,theta_v are ALWAYS tabulated relative to a's direction
c      (regardless of reverse), and momentum conservation means A
c      moves in the OPPOSITE direction to a in the CM frame, the
c      frame in which A is a forward-moving beam is reached by first
c      re-orienting into "A-forward" coordinates (z -> -z, and
c      correspondingly y -> -y to keep a right-handed system, x
c      unchanged) and THEN adding the boost -- not by simply adding
c      the boost in the file's original z direction. zsign=-1
c      implements that re-orientation for the CM-frame components
c      below, before vcm (always added with its own natural +sign)
c      is applied.
       vcm=sqrt(2.d0*mt*amu*elabbeam)/((mt+ma)*amu)
       zsign=-1.d0
       write(*,*) ' reverse=T: A is the incoming LAB beam, v_cm = ',
     x  vcm,'  (A lab energy = ',elabbeam,' MeV)'
       write(99,*) ' reverse=T: A is the incoming LAB beam, v_cm = ',
     x  vcm,'  (A lab energy = ',elabbeam,' MeV)'
      else
c      a (=b+x, the original projectile) is the incoming LAB beam
c      here -- the same particle roles and z-axis as the original
c      calculation, so no re-orientation is needed (zsign=+1)
       vcm=sqrt(2.d0*ma*amu*elabbeam)/((mt+ma)*amu)
       zsign=1.d0
       write(*,*) ' reverse=F: a is the incoming LAB beam, v_cm = ',
     x  vcm,'  (a lab energy = ',elabbeam,' MeV)'
       write(99,*) ' reverse=F: a is the incoming LAB beam, v_cm = ',
     x  vcm,'  (a lab energy = ',elabbeam,' MeV)'
      endif

c     ---------------------------------------------------------------
c     beam momentum (along the +z_lab axis, by construction of the
c     boost above -- always +z regardless of reverse) and the mass
c     of the "third" particle (not part of the chosen pair): A if
c     pair=1 (output pair b,x), x if pair=2 (output pair b,A). Used
c     below to derive the third particle's LAB kinematics via
c     momentum conservation: p_third = p_beam - p_b - p_(other
c     already-computed particle).
c     ---------------------------------------------------------------
      if(reverse) then
       mbeam=mt
      else
       mbeam=ma
      endif
      pbeamz=sqrt(2.d0*mbeam*amu*elabbeam)
      if(pair.eq.1) then
       mthird=mt
      else
       mthird=mv
      endif

c     ---------------------------------------------------------------
c     output grid setup
c     ---------------------------------------------------------------
      nthb=nint((thblabmax-thblabmin)/thblabinc)+1
      neb =nint((eblabmax -eblabmin )/eblabinc )+1
      nthx=nint((thxlabmax-thxlabmin)/thxlabinc)+1
      nex =nint((exlabmax -exlabmin )/exlabinc )+1
      nph =nint((dphlabmax-dphlabmin)/dphlabinc)+1

      rbytes=8.d0*dble(nthb)*dble(neb)*dble(nthx)*dble(nex)*dble(nph)
      write(*,*) ' output grid points: ',nthb,neb,nthx,nex,nph
      write(99,*) ' output grid points: ',nthb,neb,nthx,nex,nph
      write(*,*) ' dense 5-D array would need ',rbytes/1.d9,
     x  ' GB -- allocating it now'
      write(99,*) ' dense 5-D array would need ',rbytes/1.d9,
     x  ' GB -- allocating it now'
      allocate(hist(nthb,neb,nthx,nex,nph))
      hist=0.d0

      nth3=nint((th3labmax-th3labmin)/th3labinc)+1
      ne3 =nint((e3labmax -e3labmin )/e3labinc )+1
      allocate(hist3th(nth3))
      allocate(hist3e(ne3))
      hist3th=0.d0
      hist3e=0.d0

      write(*,*) ' '
      write(99,*) ' '
      write(*,*) '** stratified CM->LAB sampling: nsubc=',nsubc,
     x  ' nsubv=',nsubv,' nsubphi=',nsubphi,' sub-samples per cell'
      write(99,*) '** stratified CM->LAB sampling: nsubc=',nsubc,
     x  ' nsubv=',nsubv,' nsubphi=',nsubphi,' sub-samples per cell'

c     ---------------------------------------------------------------
c     main sampling loop
c     ---------------------------------------------------------------
      do ic=1,itc
       tc=tcl+(ic-1)*dtc
       do iv=1,itv
        tv=tvl+(iv-1)*dtv
        do ien=1,iten
         do ip=1,itp
          phi=phl+(ip-1)*dph

          raux=sig(ic,iv,ien,ip)
          if(raux.eq.0.d0) cycle
          ecv=ecarr(ic,iv,ien,ip)
          evv=evarr(ic,iv,ien,ip)
          if(ecv.le.0.d0.or.evv.le.0.d0) cycle

          tcr=tc*degrad
          cellwt=raux*2.d0*pi*sin(tcr)*dtc*degrad
     x               *sin(tv*degrad)*dtv*degrad*den*dph*degrad
          if(cellwt.le.wtol) cycle

          vprimec=sqrt(2.d0*ecv/(mc*amu))
          vprimev=sqrt(2.d0*evv/(mv*amu))

          wtsub=cellwt/dble(nsubc*nsubv*nsubphi)

c         clamp each cell's sub-sampling range to the physical
c         [0,180] boundary -- otherwise a boundary cell (theta=0,
c         theta=180, phi=0 or phi=180) would spill sub-samples into
c         unphysical negative/>180 territory, which for theta wrongly
c         flips the implicit phi_c=0 convention, and for phi (via
c         isymphi) wrongly duplicates part of the mirror deposit
          tc_lo=max(tc-0.5d0*dtc,0.d0)
          tc_hi=min(tc+0.5d0*dtc,180.d0)
          tv_lo=max(tv-0.5d0*dtv,0.d0)
          tv_hi=min(tv+0.5d0*dtv,180.d0)
          ph_lo=max(phi-0.5d0*dph,phl)
          ph_hi=min(phi+0.5d0*dph,phu)

c         --- sub-sample theta_c: get (theta_b_lab,Eb_lab) list ---
          do isc=1,nsubc
           tc_s=tc_lo+(isc-0.5d0)*(tc_hi-tc_lo)/dble(nsubc)
           tcr_s=tc_s*degrad
           vlx=vcm+zsign*vprimec*cos(tcr_s)
           vly=      vprimec*sin(tcr_s)
           thb_s=atan2(vly,vlx)*180.d0/pi
           if(thb_s.lt.0.d0) thb_s=thb_s+360.d0
           eb_s=0.5d0*mc*amu*(vlx*vlx+vly*vly)
           ithb=nint((thb_s-thblabmin)/thblabinc)+1
           ieb =nint((eb_s -eblabmin )/eblabinc )+1
           if(ithb.lt.1.or.ithb.gt.nthb) cycle
           if(ieb .lt.1.or.ieb .gt.neb ) cycle

c          --- sub-sample theta_v ---
           do isv=1,nsubv
            tv_s=tv_lo+(isv-0.5d0)*(tv_hi-tv_lo)/dble(nsubv)
            tvr_s=tv_s*degrad

c           --- sub-sample phi_v: the "second particle" (x if
c           --- pair=1, A if pair=2) is computed fresh at each
c           --- (theta_v,phi_v) sub-sample, since A's kinematics
c           --- need the full vector sum of b's and x's momenta and
c           --- therefore depend on phi_v too (unlike x's own
c           --- kinematics, which depend only on theta_v). Dphi_lab
c           --- is defined as the second particle's azimuth minus
c           --- b's (=0 by convention); for reverse=F this is simply
c           --- phi_v for pair=1 (and the natural atan2 result for
c           --- pair=2), while for reverse=T the A-forward
c           --- re-orientation's y-flip (see zsign above) negates it
c           --- (handled explicitly below for pair=1, and already
c           --- built into vAx/vAy for pair=2). Regardless of pair
c           --- or reverse, the isymphi mirror of whatever Dphi_lab
c           --- was computed is always (360-Dphi_lab): reflecting
c           --- phi_v -> -phi_v flips the sign of every transverse
c           --- y-component (b's, x's, and hence A's, since A's is a
c           --- linear combination of the other two), which for any
c           --- angle computed as atan2(y,x) with x unchanged and y
c           --- sign-flipped is exactly equivalent to negating that
c           --- angle -- so the *trigger* for mirroring is still
c           --- based on the input phi_s (below), but the mirrored
c           --- deposit is always at (360-Dphi_lab).
            wtphi=wtsub
            if(isymphi.and.(phi.eq.phl.or.phi.eq.phu)) wtphi=0.5d0*wtsub
            do isp=1,nsubphi
             phi_s=ph_lo+(isp-0.5d0)*(ph_hi-ph_lo)/dble(nsubphi)
             phir_s=phi_s*degrad

             if(pair.eq.1) then
c             second particle = x: depends only on theta_v
              vlx2=vcm+zsign*vprimev*cos(tvr_s)
              vly2=      vprimev*sin(tvr_s)
              thx_s=atan2(vly2,vlx2)*180.d0/pi
              if(thx_s.lt.0.d0) thx_s=thx_s+360.d0
              ex_s=0.5d0*mv*amu*(vlx2*vlx2+vly2*vly2)
c             Dphi_lab = x's azimuth minus b's (=0). The y-axis flip
c             (zsign=-1, i.e. reverse=T) negates every azimuthal
c             angle measured from it, phi_v included.
              dphi_s=phi_s
              if(reverse) dphi_s=mod(360.d0-phi_s,360.d0)
             else
c             second particle = A: momentum conservation,
c             p_A = -(p_b+p_x). The transverse-x component is
c             unaffected by the A-forward re-orientation (x is kept
c             fixed by construction); the transverse-y and
c             longitudinal-z components both pick up the zsign flip
c             before vcm (always added with its own natural +sign)
c             is applied to the (now re-oriented) z component.
              vAx=-(mc*vprimec*sin(tcr_s)
     x              +mv*vprimev*sin(tvr_s)*cos(phir_s))/mt
              vAy=zsign*(-(mv*vprimev*sin(tvr_s)*sin(phir_s))/mt)
              vAz=vcm+zsign*(-(mc*vprimec*cos(tcr_s)
     x                 +mv*vprimev*cos(tvr_s))/mt)
              thx_s=atan2(sqrt(vAx*vAx+vAy*vAy),vAz)*180.d0/pi
              if(thx_s.lt.0.d0) thx_s=thx_s+360.d0
              ex_s=0.5d0*mt*amu*(vAx*vAx+vAy*vAy+vAz*vAz)
              dphi_s=atan2(vAy,vAx)*180.d0/pi
              if(dphi_s.lt.0.d0) dphi_s=dphi_s+360.d0
             endif

c            --- third particle (A if pair=1, x if pair=2), derived
c            --- by lab-frame momentum conservation from the beam's
c            --- known total momentum and the two already-computed
c            --- particles -- no ambiguity, since both of the other
c            --- two full 3-vectors are already in hand (unlike
c            --- solving for one particle from just its direction)
             pbx=mc*amu*vly
             pby=0.d0
             pbz=mc*amu*vlx
             if(pair.eq.1) then
              p2x=mv*amu*vprimev*sin(tvr_s)*cos(phir_s)
              p2y=mv*amu*zsign*vprimev*sin(tvr_s)*sin(phir_s)
              p2z=mv*amu*vlx2
             else
              p2x=mt*amu*vAx
              p2y=mt*amu*vAy
              p2z=mt*amu*vAz
             endif
             p3x=-pbx-p2x
             p3y=-pby-p2y
             p3z=pbeamz-pbz-p2z
             p3mag=sqrt(p3x*p3x+p3y*p3y+p3z*p3z)
             th3_s=atan2(sqrt(p3x*p3x+p3y*p3y),p3z)*180.d0/pi
             if(th3_s.lt.0.d0) th3_s=th3_s+360.d0
             e3_s=p3mag*p3mag/(2.d0*mthird*amu)

             ith3=nint((th3_s-th3labmin)/th3labinc)+1
             ie3 =nint((e3_s -e3labmin )/e3labinc )+1
             if(ith3.ge.1.and.ith3.le.nth3)
     x        hist3th(ith3)=hist3th(ith3)+wtphi
             if(ie3.ge.1.and.ie3.le.ne3)
     x        hist3e(ie3)=hist3e(ie3)+wtphi
c            the mirror copy (isymphi, interior phi_s) has the same
c            transverse MAGNITUDE and z-component for the third
c            particle (only the sign of every y-component flips,
c            which a magnitude/atan2(sqrt(..),z) does not see), so
c            th3_s/e3_s are unchanged -- deposit the same bins again
             if(isymphi.and.phi_s.gt.0.d0.and.phi_s.lt.180.d0) then
              if(ith3.ge.1.and.ith3.le.nth3)
     x         hist3th(ith3)=hist3th(ith3)+wtphi
              if(ie3.ge.1.and.ie3.le.ne3)
     x         hist3e(ie3)=hist3e(ie3)+wtphi
             endif

             ithx=nint((thx_s-thxlabmin)/thxlabinc)+1
             iex =nint((ex_s -exlabmin )/exlabinc )+1
             if(ithx.lt.1.or.ithx.gt.nthx) cycle
             if(iex .lt.1.or.iex .gt.nex ) cycle

             iph=nint((dphi_s-dphlabmin)/dphlabinc)+1
             if(iph.ge.1.and.iph.le.nph) then
              hist(ithb,ieb,ithx,iex,iph)=
     x          hist(ithb,ieb,ithx,iex,iph)+wtphi
             endif

             if(isymphi.and.phi_s.gt.0.d0.and.phi_s.lt.180.d0) then
              phim=360.d0-dphi_s
              iphm=nint((phim-dphlabmin)/dphlabinc)+1
              if(iphm.ge.1.and.iphm.le.nph) then
               hist(ithb,ieb,ithx,iex,iphm)=
     x           hist(ithb,ieb,ithx,iex,iphm)+wtphi
              endif
             endif
            enddo  ! isp
           enddo  ! isv
          enddo  ! isc

         enddo  ! ip
        enddo  ! ien
       enddo  ! iv
      enddo  ! ic

c     ---------------------------------------------------------------
c     single differentials w.r.t. b's own lab angle/energy, obtained
c     by summing the raw (mb) histogram over all the OTHER four
c     dimensions -- a useful cross-check against, e.g., an inclusive
c     single-particle calculation of b alone
c     ---------------------------------------------------------------
      open(951,file='trip2lab_dsdthb.out',status='unknown')
      write(951,'(a)') '#   thb_lab      dsigma/dOmega_b'
      xchk1=0.d0
      do ithb=1,nthb
       thb=thblabmin+(ithb-1)*thblabinc
       thbr_lo=max((thb-0.5d0*thblabinc)*degrad,0.d0)
       thbr_hi=min((thb+0.5d0*thblabinc)*degrad,pi)
       domb=2.d0*pi*(cos(thbr_lo)-cos(thbr_hi))
       if(domb.le.0.d0) cycle
       wsum=0.d0
       do ieb=1,neb
        do ithx=1,nthx
         do iex=1,nex
          do iph=1,nph
           wsum=wsum+hist(ithb,ieb,ithx,iex,iph)
          enddo
         enddo
        enddo
       enddo
       write(951,'(f10.4,1x,d16.8)') thb,wsum/domb
       xchk1=xchk1+wsum
      enddo
      close(951)

      open(952,file='trip2lab_dsdEb.out',status='unknown')
      write(952,'(a)') '#   Eb_lab       dsigma/dEb'
      xchk2=0.d0
      do ieb=1,neb
       eb=eblabmin+(ieb-1)*eblabinc
       wsum=0.d0
       do ithb=1,nthb
        do ithx=1,nthx
         do iex=1,nex
          do iph=1,nph
           wsum=wsum+hist(ithb,ieb,ithx,iex,iph)
          enddo
         enddo
        enddo
       enddo
       write(952,'(f10.4,1x,d16.8)') eb,wsum/eblabinc
       xchk2=xchk2+wsum
      enddo
      close(952)

c     ---------------------------------------------------------------
c     single differentials w.r.t. the SECOND particle's (x if
c     pair=1, A if pair=2) own lab angle/energy -- same pattern as
c     the b marginals above, marginalizing over (thb,Eb,Dphi) instead
c     ---------------------------------------------------------------
      open(956,file='trip2lab_dsdthx.out',status='unknown')
      write(956,'(a)') '#   thx_lab      dsigma/dOmega_x'
      xchk6=0.d0
      do ithx=1,nthx
       thx=thxlabmin+(ithx-1)*thxlabinc
       thxr_lo=max((thx-0.5d0*thxlabinc)*degrad,0.d0)
       thxr_hi=min((thx+0.5d0*thxlabinc)*degrad,pi)
       domxfull=2.d0*pi*(cos(thxr_lo)-cos(thxr_hi))
       if(domxfull.le.0.d0) cycle
       wsum=0.d0
       do iex=1,nex
        do ithb=1,nthb
         do ieb=1,neb
          do iph=1,nph
           wsum=wsum+hist(ithb,ieb,ithx,iex,iph)
          enddo
         enddo
        enddo
       enddo
       write(956,'(f10.4,1x,d16.8)') thx,wsum/domxfull
       xchk6=xchk6+wsum
      enddo
      close(956)

      open(957,file='trip2lab_dsdEx.out',status='unknown')
      write(957,'(a)') '#   Ex_lab       dsigma/dEx'
      xchk7=0.d0
      do iex=1,nex
       ex=exlabmin+(iex-1)*exlabinc
       wsum=0.d0
       do ithx=1,nthx
        do ithb=1,nthb
         do ieb=1,neb
          do iph=1,nph
           wsum=wsum+hist(ithb,ieb,ithx,iex,iph)
          enddo
         enddo
        enddo
       enddo
       write(957,'(f10.4,1x,d16.8)') ex,wsum/exlabinc
       xchk7=xchk7+wsum
      enddo
      close(957)

c     ---------------------------------------------------------------
c     2-D marginal: d2sigma/dEb dEx (b's energy vs the second
c     particle's energy -- e.g. proton vs 238U energy for pair=2),
c     integrated over both particles' angles and Dphi_lab. This is
c     the classic energy-correlation ("kinematic locus") plot.
c     ---------------------------------------------------------------
      open(953,file='trip2lab_dsdEbdEx.out',status='unknown')
      write(953,'(a)') '#   Eb_lab       Ex_lab       d2sigma/dEbdEx'
      xchk3=0.d0
      do ieb=1,neb
       eb=eblabmin+(ieb-1)*eblabinc
       do iex=1,nex
        ex=exlabmin+(iex-1)*exlabinc
        wsum=0.d0
        do ithb=1,nthb
         do ithx=1,nthx
          do iph=1,nph
           wsum=wsum+hist(ithb,ieb,ithx,iex,iph)
          enddo
         enddo
        enddo
        write(953,'(2f11.4,1x,d16.8)') eb,ex,wsum/(eblabinc*exlabinc)
        xchk3=xchk3+wsum
       enddo
      enddo
      close(953)

c     ---------------------------------------------------------------
c     single differentials w.r.t. the THIRD particle's (A if pair=1,
c     x if pair=2) own lab angle/energy. Unlike the b/second-particle
c     marginals above, hist3th/hist3e are not summed from the main
c     5-D histogram -- they were accumulated directly during the
c     sampling loop (the third particle is not one of the 5-D
c     table's axes at all, only ever derived via momentum
c     conservation), so here they just need normalizing and writing.
c     ---------------------------------------------------------------
      open(954,file='trip2lab_dsdth3.out',status='unknown')
      write(954,'(a)') '#   th3_lab      dsigma/dOmega_3'
      xchk4=0.d0
      do ith3=1,nth3
       th3=th3labmin+(ith3-1)*th3labinc
       th3r_lo=max((th3-0.5d0*th3labinc)*degrad,0.d0)
       th3r_hi=min((th3+0.5d0*th3labinc)*degrad,pi)
       dom3=2.d0*pi*(cos(th3r_lo)-cos(th3r_hi))
       if(dom3.le.0.d0) cycle
       write(954,'(f10.4,1x,d16.8)') th3,hist3th(ith3)/dom3
       xchk4=xchk4+hist3th(ith3)
      enddo
      close(954)

      open(955,file='trip2lab_dsdE3.out',status='unknown')
      write(955,'(a)') '#   E3_lab       dsigma/dE3'
      xchk5=0.d0
      do ie3=1,ne3
       e3=e3labmin+(ie3-1)*e3labinc
       write(955,'(f10.4,1x,d16.8)') e3,hist3e(ie3)/e3labinc
       xchk5=xchk5+hist3e(ie3)
      enddo
      close(955)

      write(*,*) ' '
      write(99,*) ' '
      write(*,*) ' trip2lab_dsdthb.out : dsigma/dOmega_b vs'//
     x  ' theta_b_lab [mb/sr]  (integrated over Eb,thx,Ex,Dphi)'
      write(99,*) ' trip2lab_dsdthb.out : dsigma/dOmega_b vs'//
     x  ' theta_b_lab [mb/sr]  (integrated over Eb,thx,Ex,Dphi)'
      write(*,*) '   total from this table = ',xchk1,
     x  ' mb  (should match the overall total below)'
      write(99,*) '   total from this table = ',xchk1,
     x  ' mb  (should match the overall total below)'
      write(*,*) ' trip2lab_dsdEb.out  : dsigma/dEb vs Eb_lab'//
     x  ' [mb/MeV]  (solid-angle integrated over thb, and over'//
     x  ' thx,Ex,Dphi)'
      write(99,*) ' trip2lab_dsdEb.out  : dsigma/dEb vs Eb_lab'//
     x  ' [mb/MeV]  (solid-angle integrated over thb, and over'//
     x  ' thx,Ex,Dphi)'
      write(*,*) '   total from this table = ',xchk2,
     x  ' mb  (should match the overall total below)'
      write(99,*) '   total from this table = ',xchk2,
     x  ' mb  (should match the overall total below)'
      write(*,*) ' trip2lab_dsdthx.out : dsigma/dOmega_x vs'//
     x  ' theta_x_lab [mb/sr]  (second particle, integrated'//
     x  ' over Eb,thb,Ex,Dphi)'
      write(99,*) ' trip2lab_dsdthx.out : dsigma/dOmega_x vs'//
     x  ' theta_x_lab [mb/sr]  (second particle, integrated'//
     x  ' over Eb,thb,Ex,Dphi)'
      write(*,*) '   total from this table = ',xchk6,
     x  ' mb  (should match the overall total below)'
      write(99,*) '   total from this table = ',xchk6,
     x  ' mb  (should match the overall total below)'
      write(*,*) ' trip2lab_dsdEx.out  : dsigma/dEx vs Ex_lab'//
     x  ' [mb/MeV]  (second particle, solid-angle integrated'//
     x  ' over thx, and over thb,Eb,Dphi)'
      write(99,*) ' trip2lab_dsdEx.out  : dsigma/dEx vs Ex_lab'//
     x  ' [mb/MeV]  (second particle, solid-angle integrated'//
     x  ' over thx, and over thb,Eb,Dphi)'
      write(*,*) '   total from this table = ',xchk7,
     x  ' mb  (should match the overall total below)'
      write(99,*) '   total from this table = ',xchk7,
     x  ' mb  (should match the overall total below)'
      write(*,*) ' trip2lab_dsdEbdEx.out : d2sigma/dEbdEx vs'//
     x  ' (Eb_lab,Ex_lab) [mb/MeV/MeV]  (integrated over'//
     x  ' thb,thx,Dphi)'
      write(99,*) ' trip2lab_dsdEbdEx.out : d2sigma/dEbdEx vs'//
     x  ' (Eb_lab,Ex_lab) [mb/MeV/MeV]  (integrated over'//
     x  ' thb,thx,Dphi)'
      write(*,*) '   total from this table = ',xchk3,
     x  ' mb  (should match the overall total below)'
      write(99,*) '   total from this table = ',xchk3,
     x  ' mb  (should match the overall total below)'
      write(*,*) ' trip2lab_dsdth3.out : dsigma/dOmega_3 vs'//
     x  ' theta_3_lab [mb/sr]  (third particle, derived by'//
     x  ' momentum conservation)'
      write(99,*) ' trip2lab_dsdth3.out : dsigma/dOmega_3 vs'//
     x  ' theta_3_lab [mb/sr]  (third particle, derived by'//
     x  ' momentum conservation)'
      write(*,*) '   total from this table = ',xchk4,
     x  ' mb  (should match the overall total below)'
      write(99,*) '   total from this table = ',xchk4,
     x  ' mb  (should match the overall total below)'
      write(*,*) ' trip2lab_dsdE3.out  : dsigma/dE3 vs E3_lab'//
     x  ' [mb/MeV]  (third particle, derived by momentum'//
     x  ' conservation)'
      write(99,*) ' trip2lab_dsdE3.out  : dsigma/dE3 vs E3_lab'//
     x  ' [mb/MeV]  (third particle, derived by momentum'//
     x  ' conservation)'
      write(*,*) '   total from this table = ',xchk5,
     x  ' mb  (should match the overall total below)'
      write(99,*) '   total from this table = ',xchk5,
     x  ' mb  (should match the overall total below)'
      open(950,file='trip2lab.out',status='unknown')
      write(950,'(a)') '#   thb_lab      Eb_lab      thx_lab'//
     x  '      Ex_lab      Dphi_lab      d5sigma'
      xstotlab=0.d0
      do ithb=1,nthb
       thb=thblabmin+(ithb-1)*thblabinc
       thbr_lo=max((thb-0.5d0*thblabinc)*degrad,0.d0)
       thbr_hi=min((thb+0.5d0*thblabinc)*degrad,pi)
       domb=2.d0*pi*(cos(thbr_lo)-cos(thbr_hi))
       if(domb.le.0.d0) cycle
       do ieb=1,neb
        eb=eblabmin+(ieb-1)*eblabinc
        do ithx=1,nthx
         thx=thxlabmin+(ithx-1)*thxlabinc
         thxr_lo=max((thx-0.5d0*thxlabinc)*degrad,0.d0)
         thxr_hi=min((thx+0.5d0*thxlabinc)*degrad,pi)
         domxth=cos(thxr_lo)-cos(thxr_hi)
         if(domxth.le.0.d0) cycle
         do iex=1,nex
          ex=exlabmin+(iex-1)*exlabinc
          do iph=1,nph
           dphv=dphlabmin+(iph-1)*dphlabinc
           w=hist(ithb,ieb,ithx,iex,iph)
           if(w.le.0.d0) cycle
           domx=domxth*dphlabinc*degrad
           binvol=domb*eblabinc*domx*exlabinc
           sigma=w/binvol
           write(950,'(5f11.4,1x,d16.8)')
     x       thb,eb,thx,ex,dphv,sigma
           xstotlab=xstotlab+w
          enddo
         enddo
        enddo
       enddo
      enddo
      close(950)

      write(*,*) ' '
      write(99,*) ' '
      write(*,*) ' Total LAB cross section (from histogram sum) = ',
     x  xstotlab,' mb  (compare to input check above)'
      write(99,*) ' Total LAB cross section (from histogram sum) = ',
     x  xstotlab,' mb  (compare to input check above)'

      close(99)
      stop
      end
