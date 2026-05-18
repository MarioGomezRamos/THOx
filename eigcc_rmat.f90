!***********************************************************************
! EIGCC_RMAT  –  Coupled-channel bound-state solver on a Lagrange mesh
!
! Method: R-matrix method on a Lagrange-Legendre mesh
!   M. Hesse, J.-M. Sparenberg, F. Van Raemdonck, D. Baye
!   Nuclear Physics A 640 (1998) 37-51
!***********************************************************************
subroutine eigcc_rmat(PSI, CCMAT, ETAP, KAP2, THETA, P, E_fixed,  &
                      MAXN, MC, NODES, NP, H, M, NLAG, NS, RMAX, &
                      MAXC, EPS, IFAIL, POT, CONV)
  use channels, only: ql
  use rmat_mod,  only: tc, blo0, xle, wle, rmax0, q2
  implicit none

  !----- arguments -------------------------------------------------------
  integer,    intent(in)    :: M, MAXN, MC, NODES, NP, NLAG, NS, MAXC
  real*8,     intent(in)    :: H, RMAX, EPS, THETA, E_fixed, CONV
  real*8,     intent(in)    :: ETAP(M), POT(MAXN)
  real*8,     intent(inout) :: KAP2(M), P
  complex*16, intent(in)    :: CCMAT(M, M, MAXN)
  real*8,     intent(out)   :: PSI(NP, M)
  integer,    intent(out)   :: IFAIL

  integer              :: i1, i2, ir, irp, ich, ich2, iter, nco, lmax1, irg, info_lap, nodes_ich, nb, info2
  integer              :: nscan, isc, br_found, j
  integer              :: bc_iter, MAX_BC_IT
  real*8               :: fac, kap_c, eta_c, xc1_c, rn, xx, bval, P0, P1, fold, fnew, norm, aux, det_s
  real*8               :: xi, xi2, xj, xj2, dP
  real*8, allocatable  :: svals_m(:), VT_m(:,:), wk_m(:), q2_sol(:,:), integrand(:), lag_vals(:)
  real*8, allocatable  :: secmat_copy(:,:)
  real*8, allocatable  :: wl(:), wld(:), cmat(:,:), Rmat(:,:), secmat(:,:), kap_vec(:), xc1_vec(:)
  integer, allocatable :: ipiv2(:), ipiv_loc(:)
  real*8               :: cmat_min, cmat_max, tmpv
  real*8               :: a, b, fa, fb, fc, c, d, e, s, fs, tol1, xm, pbr, qbr, rbr
  real*8               :: ER_MIN, ER_MAX  ! Search range (adaptive based on problem type)
  real*8, dimension(12) :: samples
  real*8 :: sval
  real*8, allocatable :: xc1_old(:)
  real*8               :: tol_bc, max_bc_diff, E_old
  real*8               :: dummy_u(1,1)
  integer              :: ie_w
  real*8, parameter :: MAX_ESTEP = 100d0  ! Limit energy step size
  logical              :: is_potential_scaling_search

  !----- check NS=1 -------------------------------------------------------
  if (NS /= 1) then
    write(*,*) '[eigcc_rmat] WARNING: NS>1 not supported; using NS=1'
  end if

  ! Determine search range based on problem type
  is_potential_scaling_search = abs(THETA) < 1d-5
  if (is_potential_scaling_search) then
    ! Searching for potential scaling factor (around p ≈ 1.0)
    ER_MIN = 0.1d0
    ER_MAX = 5.0d0
    write(*,'(a)') '  [eigcc_rmat] Potential scaling search: range [0.1, 5.0]'
  else
    ! Searching for energy
    ER_MIN = -200d0
    ER_MAX = 200d0
    write(*,'(a)') '  [eigcc_rmat] Energy search: range [-200, 200] MeV'
  end if

  !----- check NS=1 -------------------------------------------------------
  if (NS /= 1) then
    write(*,*) '[eigcc_rmat] WARNING: NS>1 not supported; using NS=1'
  end if

  nb    = M * NLAG            ! total basis size
  lmax1 = maxval(ql(1:M)) + 2
  allocate(wl(lmax1), wld(lmax1), cmat(nb, nb), Rmat(M, M), &
           secmat(M, M), kap_vec(M), xc1_vec(M), ipiv2(M))

  IFAIL = 0
  write(*,'(a,i3,a,i3,a,i5)') &
    '  [eigcc_rmat] nch=',M,' nlag=',NLAG,' nb=',nb

  ! Initial seeds
  P0 = P - 0.1d0
  P1 = P + 0.1d0
  if (abs(P0) < 1d-12) P0 = -1d-3
  if (abs(P1) < 1d-12) P1 =  1d-3

  call secular_fn(P0, fold); write(*,'(a,1p,2e14.6)') '    seed0: E,f=',P0,fold
  call secular_fn(P1, fnew); write(*,'(a,1p,2e14.6)') '    seed1: E,f=',P1,fnew

  ! ------------------------------------------------------------------
  ! Bracketing scan and Brent root-finding (robust)
  ! ------------------------------------------------------------------
  ! Outer iteration over boundary parameters (Bc) as in Hesse et al.
  MAX_BC_IT = 8
  tol_bc = 1d-6

  ! IMPORTANT: Initialize boundary condition parameters from initial energy guess
  ! This ensures the secular function det(Rmat - 1/xc1) has the correct behavior
  call update_xc1_vec(P)
  write(*,'(a)') '  [eigcc_rmat] Initial boundary condition parameters computed'
  do isc = 1, M
    write(*,'(a,i3,a,1p,e14.6)') '    xc1_vec(',isc,')=',xc1_vec(isc)
  end do

  E_old = P
  allocate(xc1_old(M))
  do bc_iter = 1, MAX_BC_IT
    write(*,'(a,i3)') '  [eigcc_rmat] Bc-iter=', bc_iter
    write(*,'(a)') '  [eigcc_rmat] Scanning energy range for sign change...'
    nscan = 400
    br_found = 0
    xc1_old = xc1_vec
  a = ER_MIN
  call secular_fn(a, fa)
  do isc = 1, nscan
    b = ER_MIN + isc * (ER_MAX - ER_MIN) / nscan
    call secular_fn(b, fb)
    if (fa * fb <= 0d0) then
      br_found = 1
      write(*,'(a,1p,e14.6,a,e14.6)') '    Bracket: [', a, ', ', b, ']  f=[', fa, ', ', fb, ']'
      exit
    end if
    a = b
    fa = fb
  end do

  if (br_found == 1) then
    ! Brent's method initialization
    c = a
    call secular_fn(a, fa)
    call secular_fn(b, fb)
    fc = fa
    d = b - a
    e = d
    do iter = 1, MAXC
      if (fb * fc > 0d0) then
        c = a
        fc = fa
        d = b - a
        e = d
      end if
      if (abs(fc) < abs(fb)) then
        a = b
        b = c
        c = a
        fa = fb
        fb = fc
        fc = fa
      end if
      tol1 = 2d0 * 1d-15 * abs(b) + 0.5d0 * EPS
      xm = 0.5d0 * (c - b)
      if (abs(xm) <= tol1 .or. fb == 0d0) exit
      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        ! Attempt inverse quadratic interpolation
        s = fb / fa
        if (a == c) then
          pbr = 2d0 * xm * s
          qbr = 1d0 - s
        else
          qbr = fa / fc
          rbr = fb / fc
          pbr = s * (2d0 * xm * qbr * (qbr - rbr) - (b - a) * (rbr - 1d0))
          qbr = (qbr - 1d0) * (rbr - 1d0) * (s - 1d0)
        end if
        if (pbr > 0d0) qbr = -qbr
        pbr = abs(pbr)
        ! Accept interpolation only if it's within bounds
        if (2d0 * pbr < 3d0 * xm * qbr - abs(tol1 * qbr) .and. pbr < abs(0.5d0 * e * qbr)) then
          e = d
          d = pbr / qbr
        else
          d = xm
          e = d
        end if
      else
        d = xm
        e = d
      end if
      ! Move by at least tol1
      if (abs(d) > tol1) then
        s = b + d
      else
        s = b + sign(tol1, xm)
      end if
      call secular_fn(s, fs)
      a = b
      fa = fb
      if (fa * fs < 0d0) then
        b = s
        call secular_fn(b, fb)
      else
        c = s
        fc = fs
      end if
      write(*,'(a,i3,a,1p,e14.6,a,e12.4)') '    brent iter=',iter,'  E=',b,'  f=',fb
      if (abs(fb) < EPS) exit
    end do
    P = b
    if (iter > MAXC) IFAIL = 2
    ! Update boundary parameters from the found energy and check convergence
    call update_xc1_vec(P)
    max_bc_diff = 0d0
    do isc = 1, M
      max_bc_diff = max(max_bc_diff, abs(xc1_vec(isc) - xc1_old(isc)))
    end do
    write(*,'(a,1p,e14.6,a,1p,e14.6)') '  [eigcc_rmat] after Bc-update: E=', P, '  max|dBc|=', max_bc_diff
    if (max_bc_diff < tol_bc .and. abs(P - E_old) < EPS) then
      write(*,'(a)') '  [eigcc_rmat] Bc and E converged; breaking outer loop'
      exit
    end if
    E_old = P
    ! continue outer loop (re-run root-finder with updated xc1_vec)
  else
    write(*,'(a)') '  [eigcc_rmat] No bracket found; printing sample secular_fn values and falling back to controlled secant'
    ! Print sample secular function values for diagnosis
    samples = (/ -50d0, -20d0, -10d0, -5d0, -2.5d0, -2.25d0, -2.2d0, -2.1d0, -1.5d0, -1.0d0, 0d0, 10d0 /)
    write(*,'(a)') '    Sample E      f(E)'
    do isc = 1, size(samples)
      call secular_fn(samples(isc), sval)
      write(*,'(a,1p,e14.6,a,1p,e14.6)') '    ', samples(isc), '  ', sval
    end do
    ! Fallback: limited secant from initial seeds
    do iter = 1, MAXC
      if (abs(fnew - fold) < 1d-20 * max(abs(fnew),1d0)) then
        write(*,*) '[eigcc_rmat] stagnation'
        IFAIL = 2; exit
      end if
      dP = fnew * (P1 - P0) / (fnew - fold)
      if (abs(dP) > MAX_ESTEP) then
        dP = sign(MAX_ESTEP, dP)
      end if
      P = P1 - dP
      if (abs(P) < 1d-12) P = 1d-4
      fold = fnew; P0 = P1; P1 = P
      call secular_fn(P, fnew)
      write(*,'(a,i3,a,1p,e14.6,a,e12.4)') '    iter=',iter,'  E=',P,'  f=',fnew
      if (abs(fnew) < EPS) exit
    end do
    if (iter > MAXC) IFAIL = 2
  end if

  end do  ! end bc_iter loop
  deallocate(xc1_old)

  call secular_fn(P, fnew)
  allocate(svals_m(M), VT_m(M, M), wk_m(5*M))
  allocate(secmat_copy(M,M))
  secmat_copy = secmat
  call DGESVD('N', 'A', M, M, secmat_copy, M, svals_m, &
               dummy_u, 1, VT_m, M, wk_m, 5*M, info_lap)
  write(*,'(a,i3)') '  [eigcc_rmat] DGESVD info=', info_lap
  write(*,'(a)') '  [eigcc_rmat] singular values:'
  do i1 = 1, M
    write(*,'(a,i3,a,1p,e14.6)') '    svals(', i1, ')=', svals_m(i1)
  end do
  write(*,'(a)') '  [eigcc_rmat] VT_m (rows):'
  do i1 = 1, min(M,5)
    write(*,'(a,i3,a)') '    VT row ', i1, ':'
    do j = 1, M
      write(*,'(1p,e14.6,1x)', advance='no') VT_m(i1,j)
    end do
    write(*,*)
  end do
  deallocate(secmat_copy)
  
  allocate(lag_vals(NLAG))
  call build_cmat(P)
  allocate(q2_sol(nb, M), ipiv_loc(nb))
  q2_sol = 0d0
  do ich = 1, M
    do ir = 1, NLAG
      q2_sol((ich-1)*NLAG + ir, ich) = ((-1d0)**(NLAG-ir) / sqrt(RMAX * xle(ir)*(1d0-xle(ir))))
    end do
  end do
  call DGESV(nb, M, cmat, nb, ipiv_loc, q2_sol, nb, info_lap)
  write(*,'(a,i3)') '  [eigcc_rmat] DGESV(nb,M) info=', info_lap
  write(*,'(a)') '  [eigcc_rmat] sample q2_sol first rows:'
  do i1 = 1, min(5, nb)
    write(*,'(a,i4,a)') '    row', i1, ':'
    do j = 1, M
      write(*,'(1p,e14.6,1x)', advance='no') q2_sol(i1,j)
    end do
    write(*,*)
  end do

  write(*,'(a,1p,e14.6)') '  [eigcc_rmat] Rmat(1,1)=', Rmat(1,1)
  write(*,'(a,1p,e14.6)') '  [eigcc_rmat] secmat(1,1)=', secmat(1,1)
  write(*,'(a,1p,e14.6)') '  [eigcc_rmat] xc1_vec(1)=', xc1_vec(1)
  if (abs(xc1_vec(1)) > 1d-30) write(*,'(a,1p,e14.6)') '  [eigcc_rmat] 1/xc1_vec(1)=', 1d0/xc1_vec(1)

  ! Print cmat stats
  cmat_min = 1d300
  cmat_max = -1d300
  do i1 = 1, nb
    do j = 1, nb
      tmpv = cmat(i1,j)
      if (tmpv < cmat_min) cmat_min = tmpv
      if (tmpv > cmat_max) cmat_max = tmpv
    end do
  end do
  write(*,'(a,1p,e14.6,a,1p,e14.6)') '  [eigcc_rmat] cmat min=', cmat_min, ' max=', cmat_max

  PSI = 0d0
  do ir = 1, NP
    rn = (ir-1) * H
    xx = rn / RMAX
    if (xx < 1d-10) xx = 1d-10
    if (xx > 1d0-1d-10) xx = 1d0-1d-10
    do irp = 1, NLAG
      bval = 1d0
      do i2 = 1, NLAG
        if (i2 /= irp) bval = bval * (xx - xle(i2)) / (xle(irp) - xle(i2))
      end do
      bval = bval * sqrt(xx*(1d0-xx) / (xle(irp)*(1d0-xle(irp))))
      do ich2 = 1, M
        do ich = 1, M
          PSI(ir, ich) = PSI(ir, ich) + VT_m(ich2, ich) * q2_sol((ich-1)*NLAG + irp, ich2) * bval
        end do
      end do
    end do
  end do

  allocate(integrand(NP))
  norm = 0d0
  do ich = 1, M
    integrand(1:NP) = PSI(1:NP, ich)**2
    call sim(integrand, aux, 1, NP, H, NP)
    norm = norm + aux
  end do
  deallocate(integrand)
  write(*,'(a,1p,e14.6)') '    computed norm (integral)=', norm
  if (norm > 1d-30) PSI = PSI / sqrt(norm)

  ! Diagnostic prints: sample PSI values and stats
  write(*,'(a)') '  [eigcc_rmat] PSI diagnostics:'
  write(*,'(a,1p,e14.6)') '    PSI norm after normalization=', sum(PSI(:,1)**2)
  write(*,'(a)') '    sample PSI at radii:'
  if (NP >= 1) then
    write(*,'(a,i6,2x,1p,e14.6)') '      r idx=', 1, PSI(1,1)
  end if
  if (NP >= 2) then
    write(*,'(a,i6,2x,1p,e14.6)') '      r idx=', NP/2, PSI((NP/2),1)
  end if
  write(*,'(a,i6,2x,1p,e14.6)') '      r idx=', NP, PSI(NP,1)


  do ich = 1, M
    call count_nodes(PSI(:,ich), NP, H, nodes_ich)
    write(*,'(a,i3,a,i3)') '  [eigcc_rmat] Channel ', ich, ' nodes=', nodes_ich
  end do

  deallocate(q2_sol, ipiv_loc, svals_m, VT_m, wk_m, lag_vals)
  deallocate(wl, wld, cmat, Rmat, secmat, kap_vec, xc1_vec, ipiv2)

contains

  subroutine build_cmat(E_trial)
    real*8, intent(in) :: E_trial
    integer :: m1, m2, ir2
    real*8  :: kap_loc, xc_loc, r_ir

    cmat = 0d0
    m1 = 0
    do i1 = 1, M
      do ir = 1, NLAG
        xi = xle(ir)
        xi2 = xi * (1d0 - xi)
        r_ir = xi * RMAX
        xx = ((4d0*NLAG*(NLAG+1d0) + 3d0)*xi2 + 1d0 - 6d0*xi) / (3d0 * xi2**2)
        cmat(m1+ir, m1+ir) = xx / RMAX**2
        do ir2 = 1, ir - 1
          xj = xle(ir2)
          xj2 = xj * (1d0 - xj)
          xx = ((-1d0)**(ir-ir2) / sqrt(xi2*xj2)) * (2d0 / (xi - xj)**2)
          cmat(m1+ir, m1+ir2) = xx / RMAX**2
          cmat(m1+ir2, m1+ir) = cmat(m1+ir, m1+ir2)
        end do
        fac = dble(ql(i1)) * (dble(ql(i1)) + 1d0)
        kap_loc = -E_trial * CONV
        if (r_ir < 1d-10) r_ir = 1d-10
        cmat(m1+ir, m1+ir) = cmat(m1+ir, m1+ir) + fac/r_ir**2 + kap_loc
        
        if (abs(THETA) < 1d-5 .and. i1 == MC) then
          irg = nint(r_ir / H) + 1
          if (irg < 1)    irg = 1
          if (irg > MAXN) irg = MAXN
          cmat(m1+ir, m1+ir) = cmat(m1+ir, m1+ir) + (E_trial - 1.0d0) * POT(irg)
        end if
      end do

      ! Note: boundary parameters `xc1_vec` are computed in an outer iteration
      ! via `update_xc1_vec` (they must be kept fixed while searching for the
      ! root). Here we only leave placeholders; `xc1_vec` will be used later
      ! when forming `secmat` in `secular_fn`.
      kap_vec(i1) = 0d0

      m2 = 0
      do i2 = 1, M
        do ir = 1, NLAG
          irg = nint(xle(ir) * RMAX / H) + 1
          if (irg < 1)    irg = 1
          if (irg > MAXN) irg = MAXN
          cmat(m1+ir, m2+ir) = cmat(m1+ir, m2+ir) + dble(CCMAT(i1, i2, irg))
        end do
        m2 = m2 + NLAG
      end do
      m1 = m1 + NLAG
    end do
  end subroutine build_cmat

  subroutine update_xc1_vec(E_bound)
    real*8, intent(in) :: E_bound
    integer :: ii, ie_w_loc, irg_loc
    real*8 :: kap_loc_loc, eta_loc, xc1_c_loc, r_ir_loc
    do ii = 1, M
      kap_loc_loc = KAP2(ii)
      if (abs(THETA) > 1d-5) kap_loc_loc = kap_loc_loc + THETA * E_bound
      kap_loc_loc = sqrt(abs(kap_loc_loc))
      kap_loc_loc = max(kap_loc_loc, 1d-10)
      eta_loc = 0.5d0 * ETAP(ii) / kap_loc_loc
      ie_w_loc = 0
      call WHIT(eta_loc, RMAX, kap_loc_loc, -kap_loc_loc**2, ql(ii), wl, wld, ie_w_loc)
      if (abs(wl(ql(ii)+1)) > 1d-30) then
        xc1_c_loc = RMAX * wld(ql(ii)+1) / wl(ql(ii)+1)
      else
        xc1_c_loc = -RMAX * kap_loc_loc
      end if
      kap_vec(ii) = kap_loc_loc
      xc1_vec(ii) = xc1_c_loc
    end do
  end subroutine update_xc1_vec

  subroutine secular_fn(E_trial, fval)
    real*8, intent(in)  :: E_trial
    real*8, intent(out) :: fval
    real*8 :: work_mat(nb, nb), q2_sol(nb, M), qi
    integer :: ipiv_loc(nb), info_loc, ich1, ich2
    call build_cmat(E_trial)
    work_mat = cmat
    q2_sol = 0d0
    do ich = 1, M
      do ir = 1, NLAG
        qi = ((-1d0)**(NLAG-ir) / sqrt(RMAX * xle(ir)*(1d0-xle(ir))))
        q2_sol((ich-1)*NLAG + ir, ich) = qi
      end do
    end do
    call DGESV(nb, M, work_mat, nb, ipiv_loc, q2_sol, nb, info_loc)
    if (info_loc /= 0) then
      fval = 1d10; return
    end if
    Rmat = 0d0
    do ich1 = 1, M
      do ich2 = 1, M
        do ir = 1, NLAG
          qi = ((-1d0)**(NLAG-ir) / sqrt(RMAX * xle(ir)*(1d0-xle(ir))))
          Rmat(ich1, ich2) = Rmat(ich1, ich2) + qi * q2_sol((ich1-1)*NLAG + ir, ich2)
        end do
      end do
    end do
    Rmat = Rmat / RMAX
    secmat = Rmat
    do ich = 1, M
      if (abs(xc1_vec(ich)) < 1d-30) then
         secmat(ich,ich) = Rmat(ich,ich) - 1d20
      else
         secmat(ich,ich) = Rmat(ich,ich) - 1d0/xc1_vec(ich)
      endif
    end do
    if (M == 1) then
      fval = secmat(1, 1)
    else
      call DGETRF(M, M, secmat, M, ipiv2, info2)
      fval = 1d0
      do ich = 1, M
        fval = fval * secmat(ich, ich)
      end do
    end if
  end subroutine secular_fn

  subroutine count_nodes(wf, n, dx, nodes)
    real*8, intent(in) :: wf(n), dx
    integer, intent(in) :: n
    integer, intent(out) :: nodes
    integer :: i
    nodes = 0
    do i = 1, n-1
      if (wf(i) * wf(i+1) < 0d0) nodes = nodes + 1
    end do
  end subroutine count_nodes

end subroutine eigcc_rmat

subroutine pre_eigcc_rmat(nset, nchan, nodes, changepot)
  use globals,    only: mu12, egs, kin, written, verb
  use constants
  use sistema
  use wfs,        only: nr, dr, rvec, wfc, energ, rint, rmin, rlast, rmax
  use channels,   only: jtot, sn, partot, nchmax, ql, cindex, qj, qjc, &
                        jpiset, jpsets, qspstot, lscoup
  use parameters, only: maxchan, maxeset
  use potentials, only: ccmat, vcl
  use xcdcc,      only: binset, realwf
  use nmrv,       only: nopen
  use rmat_mod
  implicit none
  integer, intent(in)  :: nset, nodes
  integer, intent(out) :: nchan
  logical, intent(in)  :: changepot
  integer :: inc, ich, ir, np_out, nlag, ns_loc, ifail
  real*8  :: fconv, ener, conv, rau, norm, aux, rms
  real*8  :: etap(maxchan), k2(maxchan), theta_loc, p_loc, ener_loc
  real*8, allocatable :: zrma(:), psi_rmat(:,:), integrand(:)

  nlag   = jpiset(nset)%nlag
  if (nlag <= 0) nlag = 40
  ns_loc = 1
  fconv = hc**2 / 2d0 / mu12
  inc   = jpiset(nset)%inc
  nchan = jpiset(nset)%nchan
  conv  = 1d0 / fconv
  if (changepot) then
    theta_loc = 0d0
    p_loc     = 1.0d0
    ener_loc  = egs
  else
    theta_loc = -conv
    p_loc     = egs
    ener_loc  = 0d0
  end if
  call coefmat(nset, nchan)
  do ich = 1, nchan
    ql(ich)   = jpiset(nset)%lsp(ich)
    k2(ich) = (jpiset(nset)%exc(ich) - ener_loc) * conv
    etap(ich) = 2d0 * e2 * zv * zc * mu12 / hc**2
  end do
  if (.not. allocated(wfc)) then
    nchmax = 0
    allocate(wfc(jpsets, maxeset, maxchan, nr))
    allocate(energ(jpsets, maxeset))
  end if
  allocate(zrma(ns_loc * nlag))
  call rmat_ini(nlag, ns_loc, rmax, zrma)
  deallocate(zrma)
  write(*,'(/,a,i3,a,i3,a,i3,a,f8.3,a)') &
    '  [pre_eigcc_rmat] nch=',nchan,' nlag=',nlag,' ns=',ns_loc, &
    '  rmax=',rmax,' fm'
  allocate(psi_rmat(nr, nchan))
  call eigcc_rmat(psi_rmat, ccmat, etap, k2, theta_loc, p_loc, &
                  ener_loc, nr, inc, nodes, nr, dr, nchan,     &
                  nlag, ns_loc, rmax, 300, 1d-6, ifail,        &
                  vcl(ql(inc), 1:nr) * conv, conv)
  if (ifail == 1) then
    write(*,*) '[pre_eigcc_rmat] singular matrix – stopping'; stop
  end if
  if (ifail == 2) write(*,*) '[pre_eigcc_rmat] WARNING: not fully converged'
  jpiset(nset)%nex = 1
  energ(nset, 1) = -p_loc
  allocate(integrand(nr))
  norm = 0d0; rms = 0d0
  do ich = 1, nchan
    do ir = 1, nr-1
      integrand(ir) = psi_rmat(ir,ich)**2
    end do
    call sim(integrand, aux, 1, nr-1, dr, nr-1)
    write(*,'(a,i3,a,5f6.2,a,f10.6)') '  Channel',ich,              &
      ' (l,s,j,I,J)',dble(ql(ich)),sn,qj(ich),qjc(ich),jtot,'  Norm:',aux
    norm = norm + aux
    do ir = 1, nr-1
      integrand(ir) = rvec(ir)**2 * psi_rmat(ir,ich)**2
    end do
    call sim(integrand, rau, 1, nr-1, dr, nr-1)
    rms = rms + rau
  end do
  write(*,'(a,f10.6)') '  Total norm   :', norm
  write(*,'(a,f10.6)') '  sqrt(rms) fm :', sqrt(rms / max(norm,1d-30))
  deallocate(integrand)
  if (changepot) then
    write(*,'(a,f10.5)') '  Pot-scaling p=', p_loc
  else
    ener = p_loc
    energ(nset,1) = ener
    write(*,'(a,f10.5,a)') '  Converged energy: ',ener,' MeV'
  end if
  do ich = 1, nchan
    do ir = 1, nr
      wfc(nset,1,ich,ir) = psi_rmat(ir,ich) / (rvec(ir) + 1d-10)
      if (ir == 1 .and. ql(ich) /= 0) wfc(nset,1,ich,ir) = 0d0
    end do
  end do
  deallocate(ccmat)
  if (nchan > nchmax) nchmax = nchan
  np_out = int((min(rlast,rmax) - rmin) / dr)
  if (np_out < 2) np_out = int((rmax-rmin)/dr)
  write(101,'(a,i2,a,f8.4)') &
    '# ',nchan,' Channels, Eigenvalue: 1  Energy: ',-energ(nset,1)
  write(101,'(a,f4.2,a,i2)') '# J:',jtot,'  parity:',partot
  if (lscoup) then
    write(101,*) '#  channel l (sn Jc) S'
    do ich=1,nchan
      write(101,'(a,i2,a,i2,3f5.1)') '# ',ich,' :',ql(ich),sn,qjc(ich),qspstot(ich)
    end do
  else
    write(101,*) '#  channel core Jc (l sn) jn'
    do ich=1,nchan
      write(101,'(a,i2,a,i2,4f5.1)') '# ',ich,' :',cindex(ich),qjc(ich),dble(ql(ich)),sn,qj(ich)
    end do
  end if
  write(101,'(a,i5,a)') '# ',np_out,' points'
  do ir = 2, np_out
    rau = rmin + dr*dble(ir-1)
    if (rau > rlast) cycle
    write(101,'(1f8.3,2x,10g14.6)') rau,(wfc(nset,1,ich,ir),ich=1,nchan)
  end do
  write(101,*) '& '
  deallocate(psi_rmat)
end subroutine pre_eigcc_rmat
