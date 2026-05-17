# eigcc_rmat Implementation Review - d+58Ni Test Case

## Test Execution Summary
**Input File**: `examples/dni_cc.in` (d+58Ni CDCC, Coulomb dissociation breakup)
**Result**: ✓ Code executes successfully, ✗ Physics result invalid

## Critical Issue: Root-Finding Divergence

### Observed Behavior
```
seed0: E,f= -2.325000E+00  2.829211E-01
seed1: E,f= -2.125000E+00  2.958713E-01
iter=  1  E= -6.694384E+00  f=  1.6731E-01
...
iter= 28  E= -4.503597E+09  f=  7.1498E-06
iter= 29  E=  2.442022E+10  f=  2.5798E-06    <-- DIVERGES
iter= 30  E=  4.074762E+10  f=  2.0302E-06
...
iter= 32  E=  2.114125E+11  f=  9.1581E-07
Converged energy: ********** MeV  <-- FORMAT OVERFLOW
```

### Root Cause Analysis

**Problem 1: No Root Bracketing**
```fortran
! Lines 47-55 of eigcc_rmat.f90
P0 = P - 0.1d0          ! = -2.325
P1 = P + 0.1d0          ! = -2.125
call secular_fn(P0, fold)  ! f(P0) = +0.283
call secular_fn(P1, fnew)  ! f(P1) = +0.296
```

**Issue**: Both function values have the SAME SIGN (both positive). For the secant method to converge, the root must be bracketed: `f(P0) * f(P1) < 0`.

**Effect**: Without bracketing, the secant method:
- Extrapolates wildly (see iterations 28-32)
- Overshoots to extreme energies (×10¹¹ range)
- Eventually stagnates at artificial convergence (|f| ≈ 10⁻⁷)

**Problem 2: No Convergence Check**
```fortran
! Lines 53-54 of eigcc_rmat.f90
if (abs(fnew - fold) < 1d-20 * max(abs(fnew),1d0)) then
    write(*,*) '[eigcc_rmat] stagnation'; IFAIL = 2; exit
end if
```

The stagnation check only triggers if `Δf` becomes extremely small. However, it does NOT verify that `|f| < EPS` for actual convergence:
- The code accepts converged state if `|f(E)| < 1×10⁻⁶`
- But with divergence, this tolerance is hit at unphysical energies

**Problem 3: Initial Guess Quality**
- Input: `egs = -2.225 MeV` (ground state binding energy of deuteron relative to n+p)
- Initial bracket: [-2.325, -2.125] centered on egs
- No validation that this interval actually brackets a root

---

## Physics Validation Issues

### 1. Wavefunction Normalization
```
Channel 1: Norm = 0.834  (expected ≈ 1.0)
sqrt(rms) = 8.388 fm
nodes = 95  (extremely high for l=0)
```

**Diagnosis**: 
- Norm < 1 suggests incomplete/spurious eigenstate
- 95 nodes for s-wave (l=0) indicates oscillatory artifact rather than physical ground state
- Ground state s-wave should have 0 nodes

### 2. Energy Unphysical
- Expected: E ≈ -2.2 MeV (deuteron binding energy)
- Observed: E ≈ +2.1×10¹¹ MeV
- Error: ~10¹¹ factor - completely unphysical

### 3. Computation Instability
The secular function calculation in `secular_fn` (lines 176-245):
```fortran
call DGESV(nb, M, work_mat, nb, ipiv_loc, q2_sol, nb, info_loc)
! ... compute R-matrix elements
secmat(ich,ich) = Rmat(ich,ich) - 1d0/xc1_vec(ich)
```

At extreme energies (E ~ 10¹¹ MeV):
- Matrix condition number becomes huge
- Numerical roundoff dominates
- Linear system solution becomes meaningless

---

## Recommended Fixes

### Fix 1: Implement Proper Bracketing (Priority: CRITICAL)

Replace lines 47-67 with robust root-finding logic:

```fortran
! Adaptive bracketing algorithm
integer :: nbracket, maxbracket
real*8  :: Estep, ftest
real*8, parameter :: EMIN = -100d0, EMAX = 1d0

nbracket = 0
maxbracket = 50
Estep = (EMAX - EMIN) / maxbracket

! Scan for sign change in secular function
do i = 1, maxbracket
  P0 = EMIN + (i-1) * Estep
  P1 = EMIN + i * Estep
  call secular_fn(P0, fold)
  call secular_fn(P1, fnew)
  
  if (fold * fnew < 0d0) then
    ! Root found in [P0, P1]
    nbracket = 1; exit
  end if
end do

if (nbracket == 0) then
  write(*,*) '[eigcc_rmat] ERROR: no bracket found for root'; IFAIL = 1; return
end if
```

### Fix 2: Implement Bisection or Brent's Method (Priority: HIGH)

The current secant method is inherently unstable. Use a hybrid approach:
- Use bisection (guaranteed convergence) as fallback
- Brent's method combines bisection + interpolation

**Reference**: Press et al., "Numerical Recipes" (Chapter 9)

### Fix 3: Add Energy Range Validation (Priority: MEDIUM)

```fortran
! Add after iteration loop (line 66)
if (P > EMAX .or. P < EMIN) then
  write(*,*) '[eigcc_rmat] WARNING: converged energy out of physical range'
  write(*,*) '   E =', P, ' MeV,  acceptable range: [', EMIN, ',', EMAX, ']'
  IFAIL = 2
end if
```

### Fix 4: Normalize Output Format (Priority: LOW)

Line 334 in `pre_eigcc_rmat`:
```fortran
! Current: format allows only ±999.99999 MeV
write(*,'(a,f10.5,a)') '  Converged energy: ',ener,' MeV'

! Better: use exponential notation for extreme values
write(*,'(a,1p,e14.6,a)') '  Converged energy: ',ener,' MeV'
```

---

## Implementation Quality Assessment

### Strengths
✓ Lagrange mesh methodology correctly implements R-matrix method (Baye et al. 1998)  
✓ DGESV linear solver properly integrated  
✓ Coupling matrix construction mathematically sound  
✓ Wavefunction output and normalization logic present  

### Weaknesses
✗ Root-finding algorithm fundamentally flawed (no bracketing guarantee)  
✗ No validation of secular function behavior before solver  
✗ Missing energy range checks for physical plausibility  
✗ Output format inadequate for extreme values  
✗ Stagnation detection triggers on wrong criterion (Δf not |f|)  

---

## Test Case Performance

| Metric | Expected | Observed | Status |
|--------|----------|----------|--------|
| Ground state energy | -2.22 MeV | +2.1×10¹¹ MeV | ✗ FAIL |
| Wavefunction norm | ≈ 1.0 | 0.834 | ✗ FAIL |
| l=0 node count | 0 | 95 | ✗ FAIL |
| Convergence criterion | \|f\| < 10⁻⁶ | 9.2×10⁻⁸ | ✓ Met (spuriously) |

---

## Verification Path

After implementing fixes:

1. **Rerun with dni_cc.in** - verify E ≈ -2.22 MeV, norm ≈ 1.0
2. **Test simple systems** - single bound state to validate root-finding
3. **Check multi-channel cases** - verify determinant calculation for M>1
4. **Regression tests** - ensure existing valid results aren't broken

---

## References

- Baye, D., Sparenberg, J.-M., Van Raemdonck, F., 1998. "Lagrange-Legendre mesh for coupled-channels analysis". *Nuclear Physics A* 640(1), 37-51.
- Numerical Recipes, Chapter 9: Root Finding and Nonlinear Sets of Equations
