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
  
  ! Adaptive bracketing variables
  integer :: nbracket, maxbracket, ibrac
  real*8  :: Emin, Emax, Estep, ftest, ftmp
  integer :: iter_bis, maxiter_bis
  real*8  :: Emid, fmid
  real*8, parameter :: EPS_BIS = 1d-10  ! Bisection tolerance

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

  ! IMPORTANT: Initialize boundary condition parameters from initial energy guess
  ! This ensures the secular function det(Rmat - 1/xc1) has the correct behavior
  call update_xc1_vec(P)
  write(*,'(a)') '  [eigcc_rmat] Initial boundary condition parameters computed'
  do isc = 1, M
    write(*,'(a,i3,a,1p,e14.6)') '    xc1_vec(',isc,')=',xc1_vec(isc)
  end do

  ! Initial seeds
  P0 = P - 0.1d0
  P1 = P + 0.1d0
  if (abs(P0) < 1d-12) P0 = -1d-3
  if (abs(P1) < 1d-12) P1 =  1d-3

  call secular_fn(P0, fold)
  call secular_fn(P1, fnew)

  ! Root Search Phase (Scan + Bisection)
  nbracket = 0
  Emin = ER_MIN
  Emax = ER_MAX
  maxbracket = 100
  

  call secular_fn(Emin, fold)
  do ibrac = 1, maxbracket
    Estep = (Emax - Emin) / maxbracket
    P0 = Emin + (ibrac-1) * Estep
    P1 = Emin + ibrac * Estep
    
    call secular_fn(P0, fold)
    call secular_fn(P1, fnew)
    
    
    if (fold * fnew <= 0d0) then
      nbracket = 1
      ! Bracket found
      exit
    end if
  end do
  
  if (nbracket == 0) then
    write(*,'(a)') '[eigcc_rmat] ERROR: No root bracket found in energy range'
    write(*,'(a,1p,e14.6,a,1p,e14.6,a)') '  Scan range: [', ER_MIN, ',', ER_MAX, ']'
    IFAIL = 1
    return
  end if
  

  maxiter_bis = 50
  do iter_bis = 1, maxiter_bis
    Emid = 0.5d0 * (P0 + P1)
    call secular_fn(Emid, fmid)
    
    if (abs(P1 - P0) < EPS_BIS) then

      P = Emid
      exit
    end if
    
    if (fold * fmid <= 0d0) then
      P1 = Emid
      fnew = fmid
    else
      P0 = Emid
      fold = fmid
    end if
  end do
  
  ! Wavefunction reconstruction
  block
    real*8 :: cmat_wb(nb, nb), w_ev(nb), work_ev(3*nb), coefficients(nb)
    integer :: info_ev, ich, ir, ir2, irp, i2, m1, ie_w_loc
    real*8 :: bval, xx, rn, r_trans
    real*8 :: u_boundary, kap_loc_loc, eta_loc
    real*8, allocatable :: wl_rmax(:), wld_rmax(:), wl_rn(:), wld_rn(:)
    
    if (is_potential_scaling_search) then
      call update_xc1_vec(E_fixed)
    else
      call update_xc1_vec(P)
    end if
    
    call build_cmat(P)
    cmat_wb = cmat
    
    do ich = 1, M
      m1 = (ich-1)*NLAG
      do ir = 1, NLAG
        do ir2 = 1, NLAG
          cmat_wb(m1+ir, m1+ir2) = cmat_wb(m1+ir, m1+ir2) - xc1_vec(ich) * blo0(ir, ir2) / RMAX
        end do
      end do
    end do
    
    call DSYEV('V', 'U', nb, cmat_wb, nb, w_ev, work_ev, 3*nb, info_ev)

    if (info_ev == 0) then
      coefficients = cmat_wb(:, min(NODES + 1, nb))  ! Column corresponding to state with NODES nodes
    else
      coefficients = 0d0
      do ich = 1, M
        coefficients((ich-1)*NLAG + 1) = 1.0d0
      end do
    end if
    
    PSI = 0d0
    r_trans = RMAX * xle(NLAG)
    
    do ich = 1, M
      ! Compute boundary value u(RMAX)
      u_boundary = 0d0
      do irp = 1, NLAG
        u_boundary = u_boundary + coefficients((ich-1)*NLAG + irp) * q2(irp, 1)
      end do
      
      ! Setup Whittaker function parameters
      if (is_potential_scaling_search) then
        kap_loc_loc = KAP2(ich) + THETA * E_fixed
      else
        kap_loc_loc = KAP2(ich) + THETA * P
      end if
      kap_loc_loc = sqrt(abs(kap_loc_loc))
      if (kap_loc_loc < 1d-10) kap_loc_loc = 1d-10
      eta_loc = 0.5d0 * ETAP(ich) / kap_loc_loc
      
      allocate(wl_rmax(ql(ich)+1), wld_rmax(ql(ich)+1))
      ie_w_loc = 0
      call WHIT(eta_loc, RMAX, kap_loc_loc, -kap_loc_loc**2, ql(ich), wl_rmax, wld_rmax, ie_w_loc)
      
      allocate(wl_rn(ql(ich)+1), wld_rn(ql(ich)+1))
      
      do ir = 1, NP
        rn = (ir-1) * H
        if (ir == 1) then
          PSI(ir, ich) = 0d0
          cycle
        end if
        if (rn < r_trans) then
          ! Interior Lagrange-mesh expansion
          xx = rn / RMAX
          if (xx < 1d-10) xx = 1d-10
          if (xx > 1d0-1d-10) xx = 1d0-1d-10
          do irp = 1, NLAG
            bval = 1d0
            do i2 = 1, NLAG
              if (i2 /= irp) bval = bval * (xx - xle(i2)) / (xle(irp) - xle(i2))
            end do
            bval = bval * sqrt(xx*(1d0-xx) / (xle(irp)*(1d0-xle(irp))))
            PSI(ir, ich) = PSI(ir, ich) + coefficients((ich-1)*NLAG + irp) * bval / sqrt(RMAX * wle(irp))
          end do
        else
          ! Exterior matched Whittaker tail
          ie_w_loc = 0
          call WHIT(eta_loc, rn, kap_loc_loc, -kap_loc_loc**2, ql(ich), wl_rn, wld_rn, ie_w_loc)
          if (abs(wl_rmax(ql(ich)+1)) > 1d-30) then
            PSI(ir, ich) = u_boundary * wl_rn(ql(ich)+1) / wl_rmax(ql(ich)+1)
          else
            PSI(ir, ich) = 0d0
          end if
        end if
      end do
      
      deallocate(wl_rmax, wld_rmax, wl_rn, wld_rn)
    end do
  end block

  allocate(integrand(NP))
  norm = 0d0
  do ich = 1, M
    integrand(1:NP) = PSI(1:NP, ich)**2
    call sim(integrand, aux, 1, NP, H, NP)
    norm = norm + aux
  end do
  deallocate(integrand)

  if (norm > 1d-30) then
    PSI = PSI / sqrt(norm)
    ! Enforce a single global sign (same as eigcc.f) based on the first channel
    if (PSI(min(3, NP), 1) < 0d0) then
      PSI = -PSI
    end if
  end if




  do ich = 1, M
    call count_nodes(PSI(:,ich), NP, H, nodes_ich)

  end do

  if (allocated(q2_sol)) deallocate(q2_sol)
  if (allocated(ipiv_loc)) deallocate(ipiv_loc)
  if (allocated(svals_m)) deallocate(svals_m)
  if (allocated(VT_m)) deallocate(VT_m)
  if (allocated(wk_m)) deallocate(wk_m)
  if (allocated(lag_vals)) deallocate(lag_vals)
  if (allocated(wl)) deallocate(wl)
  if (allocated(wld)) deallocate(wld)
  if (allocated(cmat)) deallocate(cmat)
  if (allocated(Rmat)) deallocate(Rmat)
  if (allocated(secmat)) deallocate(secmat)
  if (allocated(kap_vec)) deallocate(kap_vec)
  if (allocated(xc1_vec)) deallocate(xc1_vec)
  if (allocated(ipiv2)) deallocate(ipiv2)

contains

  subroutine build_cmat(E_trial)
    real*8, intent(in) :: E_trial
    integer :: m1, m2, ir2, ir, i1, i2, idx_low, idx_high
    real*8  :: kap_loc, xc_loc, r_ir, fac, t, v_val

    cmat = 0d0
    m1 = 0
    do i1 = 1, M
      do i2 = 1, M
        m2 = (i2-1) * NLAG
        if (i1 == i2) then
          do ir = 1, NLAG
            do ir2 = 1, NLAG
              cmat(m1+ir, m1+ir2) = tc(ir, ir2, 1)
            end do
          end do
          fac = dble(ql(i1)) * (dble(ql(i1)) + 1d0)
          if (is_potential_scaling_search) then
            kap_loc = KAP2(i1)
            if (abs(THETA) > 1d-5) kap_loc = kap_loc + THETA * E_fixed
          else
            kap_loc = KAP2(i1)
            if (abs(THETA) > 1d-5) kap_loc = kap_loc + THETA * E_trial
          end if
          do ir = 1, NLAG
            r_ir = xle(ir) * RMAX
            if (r_ir < 1d-10) r_ir = 1d-10
            cmat(m1+ir, m1+ir) = cmat(m1+ir, m1+ir) + fac/r_ir**2 + kap_loc
          end do
        end if
        do ir = 1, NLAG
          r_ir = xle(ir) * RMAX
          idx_low = int(r_ir / H) + 1
          idx_high = idx_low + 1
          if (idx_low < 1)    idx_low = 1
          if (idx_low > MAXN) idx_low = MAXN
          if (idx_high < 1)    idx_high = 1
          if (idx_high > MAXN) idx_high = MAXN
          if (idx_low == idx_high) then
            v_val = CCMAT(i1, i2, idx_low)
          else
            t = (r_ir - (idx_low - 1) * H) / H
            v_val = (1d0 - t) * CCMAT(i1, i2, idx_low) + t * CCMAT(i1, i2, idx_high)
          end if
          if (is_potential_scaling_search .and. i1 == i2 .and. i1 == MC) then
            cmat(m1+ir, m2+ir) = cmat(m1+ir, m2+ir) + E_trial * v_val
          else
            cmat(m1+ir, m2+ir) = cmat(m1+ir, m2+ir) + v_val
          end if
        end do
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
    real*8 :: cmat_wb(nb, nb), w_ev(nb), work_ev(3*nb)
    integer :: info_ev, ich, ir, ir2, m1
    
    if (is_potential_scaling_search) then
      call update_xc1_vec(E_fixed)
    else
      call update_xc1_vec(E_trial)
    end if
    
    call build_cmat(E_trial)
    cmat_wb = cmat
    
    do ich = 1, M
      m1 = (ich-1)*NLAG
      do ir = 1, NLAG
        do ir2 = 1, NLAG
          cmat_wb(m1+ir, m1+ir2) = cmat_wb(m1+ir, m1+ir2) - xc1_vec(ich) * blo0(ir, ir2) / RMAX
        end do
      end do
    end do
    
    call DSYEV('N', 'U', nb, cmat_wb, nb, w_ev, work_ev, 3*nb, info_ev)
    if (info_ev /= 0) then
      fval = 1d10
    else
      fval = w_ev(min(NODES + 1, nb))
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
  real*8  :: fconv, ener, conv, rau, norm, aux, rms, rmat_boundary
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
  if (rint > 0d0 .and. rint < rmax) then
    rmat_boundary = rint
    write(*,'(/,a,f8.3,a,f8.3,a)') &
      '  [pre_eigcc_rmat] Extrapolating bound states: R-matrix boundary rint =', rint, &
      ' fm, extrapolating to rmax =', rmax, ' fm'
  else
    rmat_boundary = rmax
  end if
  allocate(zrma(ns_loc * nlag))
  call rmat_ini(nlag, ns_loc, rmat_boundary, zrma)
  deallocate(zrma)
  write(*,'(/,a,i3,a,i3,a,i3,a,f8.3,a)') &
    '  [pre_eigcc_rmat] nch=',nchan,' nlag=',nlag,' ns=',ns_loc, &
    '  rmat_boundary=',rmat_boundary,' fm'
  allocate(psi_rmat(nr, nchan))
  call eigcc_rmat(psi_rmat, ccmat, etap, k2, theta_loc, p_loc, &
                  ener_loc, nr, inc, nodes, nr, dr, nchan,     &
                  nlag, ns_loc, rmat_boundary, 300, 1d-6, ifail, &
                  vcl(ql(inc), 1:nr) * conv, conv)
  if (ifail == 1) then
    write(*,*) '[pre_eigcc_rmat] singular matrix – stopping'; stop
  end if
  if (ifail == 2) write(*,*) '[pre_eigcc_rmat] WARNING: not fully converged'
  jpiset(nset)%nex = 1
  energ(nset, 1) = -p_loc
  allocate(integrand(nr))
  norm = 0d0; rms = 0d0
  block
    real*8, allocatable :: comp_norms(:), comp_rmss(:), comp_ancs(:)
    allocate(comp_norms(nchan), comp_rmss(nchan), comp_ancs(nchan))
    do ich = 1, nchan
      do ir = 1, nr-1
        integrand(ir) = psi_rmat(ir,ich)**2
      end do
      call sim(integrand, aux, 1, nr-1, dr, nr-1)
      write(*,'(a,i3,a,5f6.2,a,f10.6)') '  Channel',ich,              &
        ' (l,s,j,I,J)',dble(ql(ich)),sn,qj(ich),qjc(ich),jtot,'  Norm:',aux
      norm = norm + aux
      comp_norms(ich) = aux
      
      do ir = 1, nr-1
        integrand(ir) = rvec(ir)**2 * psi_rmat(ir,ich)**2
      end do
      call sim(integrand, rau, 1, nr-1, dr, nr-1)
      rms = rms + rau
      comp_rmss(ich) = sqrt(rau / max(aux, 1d-30))

      block
        real*8 :: kap_loc_loc, eta_loc, u_boundary, anc_val
        real*8, allocatable :: wl_rmax(:), wld_rmax(:)
        integer :: ie_w_loc, ir_boundary
        
        ! Compute ANC
        kap_loc_loc = k2(ich) + theta_loc * p_loc
        kap_loc_loc = sqrt(abs(kap_loc_loc))
        if (kap_loc_loc < 1d-10) kap_loc_loc = 1d-10
        eta_loc = 0.5d0 * etap(ich) / kap_loc_loc
        
        allocate(wl_rmax(ql(ich)+1), wld_rmax(ql(ich)+1))
        ie_w_loc = 0
        call WHIT(eta_loc, rmat_boundary, kap_loc_loc, -kap_loc_loc**2, ql(ich), wl_rmax, wld_rmax, ie_w_loc)
        
        ir_boundary = nint(rmat_boundary / dr) + 1
        if (ir_boundary < 1) ir_boundary = 1
        if (ir_boundary > nr) ir_boundary = nr
        u_boundary = psi_rmat(ir_boundary, ich)
        
        if (abs(wl_rmax(ql(ich)+1)) > 1d-30) then
          anc_val = u_boundary / wl_rmax(ql(ich)+1)
        else
          anc_val = 0d0
        end if
        comp_ancs(ich) = anc_val
        deallocate(wl_rmax, wld_rmax)
      end block
    end do

    write(*,'(/,a)') '  ===================================================================================='
    write(*,'(a)')   '                       Component Properties (bastype=8 R-matrix)'
    write(*,'(a)')   '  ===================================================================================='
    write(*,'(a)')   ' Chan      l      s      j      I      J       Norm       RMS (fm)     ANC (fm^-1/2)'
    write(*,'(a)')   '  ------------------------------------------------------------------------------------'
    do ich = 1, nchan
      write(*,'(i5,5f7.2,f11.6,f15.4,es18.5)') &
        ich, dble(ql(ich)), sn, qj(ich), qjc(ich), jtot, comp_norms(ich), comp_rmss(ich), comp_ancs(ich)
    end do
    write(*,'(a)')   '  ===================================================================================='
    deallocate(comp_norms, comp_rmss, comp_ancs)
  end block
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
      if (ir == 1) then
        wfc(nset,1,ich,ir) = 0d0
      else
        wfc(nset,1,ich,ir) = psi_rmat(ir,ich) / rvec(ir)
      end if
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
    write(101,'(1f8.3)', advance='no') rau
    do ich = 1, nchan
      write(101,'(1x,1p,e14.6)', advance='no') real(wfc(nset,1,ich,ir))
    end do
    write(101,*) ''
  end do
  write(101,*) '& '
  deallocate(psi_rmat)
end subroutine pre_eigcc_rmat
