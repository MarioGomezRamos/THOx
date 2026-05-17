! PROPOSED FIX FOR eigcc_rmat root-finding algorithm
! 
! Replace the root-finding section (lines 47-67) with this improved version
! that properly brackets the secular function before using secant method

! ============================================================================
! IMPROVED ROOT-FINDING: Bisection with Secant Method Fallback
! ============================================================================

  ! Adaptive bracketing phase
  integer :: nbracket, maxbracket, ibrac
  real*8  :: Emin, Emax, Estep, ftest, ftmp
  real*8, parameter :: ERANGE_MIN = -100d0, ERANGE_MAX = 1d0  ! Bound state range
  
  ! Initialize and search for bracket
  nbracket = 0
  Emin = ERANGE_MIN
  Emax = ERANGE_MAX
  maxbracket = 100
  
  ! Phase 1: Scan for sign change in secular function
  write(*,'(a)') '  [eigcc_rmat] Scanning for root bracket...'
  
  call secular_fn(Emin, fold)
  do ibrac = 1, maxbracket
    Estep = (Emax - Emin) / maxbracket
    P0 = Emin + (ibrac-1) * Estep
    P1 = Emin + ibrac * Estep
    
    call secular_fn(P0, fold)
    call secular_fn(P1, fnew)
    
    ! Check for sign change
    if (fold * fnew < 0d0) then
      nbracket = 1
      write(*,'(a,1p,2e14.6,a,2e14.6)') &
        '  [eigcc_rmat] Bracket found: E in [',P0,',',P1,'] with f=[',fold,',',fnew,']'
      exit
    end if
  end do
  
  if (nbracket == 0) then
    write(*,'(a)') '[eigcc_rmat] ERROR: No root bracket found in energy range'
    write(*,'(a,1p,2e14.6)') '  Scan range: [', ERANGE_MIN, ',', ERANGE_MAX, ']'
    write(*,'(a)') '  Suggestions:'
    write(*,'(a)') '    1. Check coupling matrix CCMAT is properly initialized'
    write(*,'(a)') '    2. Verify potential parameters are physical'
    write(*,'(a)') '    3. Check that bound state energy is within expected range'
    IFAIL = 1
    return
  end if
  
  ! Phase 2: Bisection root refinement (guaranteed convergence)
  write(*,'(a)') '  [eigcc_rmat] Bisection phase...'
  
  integer :: iter_bis, maxiter_bis, iter
  real*8  :: Emid, fmid
  real*8, parameter :: EPS_BIS = 1d-10  ! Bisection tolerance
  
  maxiter_bis = 50
  
  do iter_bis = 1, maxiter_bis
    Emid = 0.5d0 * (P0 + P1)
    call secular_fn(Emid, fmid)
    
    ! Check convergence
    if (abs(P1 - P0) < EPS_BIS) then
      write(*,'(a,i3,a,1p,e14.6)') &
        '  [eigcc_rmat] Bisection converged at iter=',iter_bis,'  E=',Emid
      P = Emid
      exit
    end if
    
    ! Update bracket
    if (fold * fmid < 0d0) then
      P1 = Emid
      fnew = fmid
    else
      P0 = Emid
      fold = fmid
    end if
  end do
  
  if (iter_bis > maxiter_bis) then
    write(*,'(a)') '[eigcc_rmat] WARNING: Bisection exceeded max iterations'
    write(*,'(a,1p,e14.6,a,e14.6)') '  Final bracket: [',P0,',',P1,']'
  end if
  
  ! Phase 3: Refine with secant method (faster convergence near root)
  write(*,'(a)') '  [eigcc_rmat] Secant refinement phase...'
  
  P1 = P0
  P0 = P0 - 0.1d0
  call secular_fn(P0, fold)
  call secular_fn(P1, fnew)
  
  do iter = 1, MAXC
    ! Stagnation check
    if (abs(fnew - fold) < 1d-20 * max(abs(fnew), 1d0)) then
      write(*,'(a,i3)') '[eigcc_rmat] Stagnation after ', iter, ' iterations'
      IFAIL = 2
      exit
    end if
    
    ! Secant step
    P = P1 - fnew * (P1 - P0) / (fnew - fold)
    
    ! Safeguard against wild extrapolation
    if (P < ERANGE_MIN .or. P > ERANGE_MAX) then
      write(*,'(a,1p,e14.6)') &
        '[eigcc_rmat] WARNING: Secant step out of range, E=',P
      P = 0.5d0 * (P0 + P1)  ! Bisect instead
    end if
    
    if (abs(P) < 1d-12) P = 1d-4
    
    fold = fnew
    P0 = P1
    P1 = P
    call secular_fn(P, fnew)
    
    write(*,'(a,i3,a,1p,e14.6,a,e12.4)') &
      '    iter=',iter,'  E=',P,'  f=',fnew
    
    if (abs(fnew) < EPS) exit
  end do
  
  if (iter > MAXC) then
    write(*,'(a,i3)') '[eigcc_rmat] WARNING: Exceeded max iterations', MAXC
    IFAIL = 2
  end if
  
  ! Energy range validation
  if (P > ERANGE_MAX .or. P < ERANGE_MIN) then
    write(*,'(a,1p,e14.6,a,1p,2e14.6)') &
      '[eigcc_rmat] WARNING: Converged energy out of physical range: E=',P, &
      '  acceptable: [',ERANGE_MIN,',',ERANGE_MAX,']'
    IFAIL = 2
  end if
  
  ! Final convergence check
  call secular_fn(P, fnew)
  if (abs(fnew) > EPS) then
    write(*,'(a,1p,2e14.6)') &
      '[eigcc_rmat] WARNING: Final |f|=',abs(fnew),' exceeds tolerance',EPS
    IFAIL = 2
  end if

! ============================================================================
! END OF PROPOSED FIX
! ============================================================================
