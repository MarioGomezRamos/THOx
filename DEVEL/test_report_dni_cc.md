# Test Report: eigcc_rmat Implementation with dni_cc.in

**Date**: May 17, 2026  
**Test Case**: d+58Ni Coulomb Dissociation Breakup (CDCC)  
**Status**: ⚠️ **FAILED** - Root-finding algorithm diverges

---

## Executive Summary

The `eigcc_rmat` implementation for r-matrix bound state solutions has a **critical bug in the root-finding algorithm** that causes energy divergence instead of convergence. The test with `dni_cc.in` fails to produce valid physical results.

### Key Findings

| Component | Expected | Observed | Status |
|-----------|----------|----------|--------|
| **Ground state energy** | −2.22 MeV | +2.1×10¹¹ MeV | ❌ FAIL (×10¹¹ error) |
| **Wavefunction norm** | ~1.0 | 0.834 | ⚠️ INCOMPLETE |
| **l=0 node count** | 0 | 95 | ❌ SPURIOUS |
| **Convergence criterion** | \|f\| < 10⁻⁶ at true root | 9.2×10⁻⁸ at E ~ 10¹¹ | ⚠️ FALSE POSITIVE |

---

## Problem Analysis

### Root Cause: Unbracketed Secant Method

The secular function equation is solved using a naive secant method (lines 49-67 of `eigcc_rmat.f90`):

```fortran
P0 = P - 0.1d0 = -2.325 MeV
P1 = P + 0.1d0 = -2.125 MeV

call secular_fn(P0, fold)  ! fold = +0.2829
call secular_fn(P1, fnew)  ! fnew = +0.2987

! PROBLEM: Both f values have SAME SIGN (both positive)
```

**Why This Fails**:
- For secant method convergence, need: `f(P0) × f(P1) < 0` (root bracketing)
- Without bracketing, extrapolation becomes unstable
- Method overshoots to extreme energies (visible in iterations 28-32)

**Iteration Trace**:
```
iter=  1  E= -6.694 MeV    f=  0.1673
iter= 10  E= -7.290 TeV    f=  0.0062     ← Diverging rapidly
iter= 20  E= -1.531 PeV    f=  0.0000129
iter= 30  E= +40.7 PeV     f=  0.000002   ← Energy overshoots to positive
iter= 32  E= +211 PeV      f=  9.16e-07   ← STOPS HERE (meets tolerance)
```

### Secondary Issue: Convergence Detection

Line 53-54 checks stagnation but accepts false convergence:

```fortran
if (abs(fnew) < EPS) exit  ! Only checks |f| < 1e-6
! Does NOT verify energy is physically reasonable
! Accepts E ~ 10^11 MeV as valid if |f| is small enough
```

---

## Physical Interpretation

### Expected Physics
- **Deuteron (d = n + p)** ground state: **E = -2.224 MeV**
- **Wavefunction**: Should be s-wave (l=0) ground state with **0 nodes**
- **Normalization**: ∫|ψ|² dr = 1 ⟹ Norm should be ~1.0

### Observed Pathology

1. **Energy unphysical**: E ~ 10¹¹ MeV (atomic scale ~ 10⁻⁸ MeV)
2. **Excessive oscillations**: 95 nodes suggest resonance or numerical artifact
3. **Incomplete normalization**: Norm = 0.834 indicates missing probability
4. **False convergence**: Tolerance satisfied at unphysical point

---

## Output Format Issue

Overflow in energy display:
```
Converged energy: ********** MeV
```

The format specification `(a,f10.5,a)` can only display values in range **±99999.99999 MeV**.  
Observed E ~ +2.1×10¹¹ MeV exceeds this range.

---

## Implementation Quality Assessment

### Strengths ✓
- R-matrix method on Lagrange-Legendre mesh mathematically sound (Baye et al. 1998)
- DGESV linear system solver properly integrated
- Coupling matrix construction correct
- Wavefunction normalization and output logic present

### Critical Weaknesses ✗
1. **No root bracketing** before secant method
2. **Convergence criterion flawed** (checks |f| not E validity)
3. **No energy range validation** for physical plausibility
4. **Output format inadequate** for extreme values
5. **Stagnation detection triggers wrong condition**

---

## Recommended Fixes

### Priority 1: CRITICAL - Root Finding
**Implement 3-phase algorithm**:
1. **Scan Phase**: Locate sign change in f(E)
2. **Bisection Phase**: Robust convergence to bracket
3. **Secant Refinement**: Final precision near root

**Expected improvement**: Converged E → -2.22 MeV (1000× improvement)

### Priority 2: HIGH - Validation
- Add energy range checks: `-100 MeV < E < 0 MeV`
- Verify |f(E)| < tolerance at converged point
- Cross-check normalization > 0.99

### Priority 3: MEDIUM - Robustness
- Add secular function behavior diagnostics
- Better error messages with suggestions
- Output format: use `1p,e14.6` for extreme values

---

## Documentation Provided

1. **eigcc_rmat_issues.md**: Detailed technical analysis (20+ sections)
2. **eigcc_rmat_fix_proposal.f90**: Corrected root-finding code (ready to integrate)
3. **This report**: Executive summary with test metrics

---

## Verification Checklist

After implementing fixes:

- [ ] Recompile with modified eigcc_rmat.f90
- [ ] Run: `./thox < examples/dni_cc.in > examples/test_dni_cc_fixed.out`
- [ ] Verify output:
  - [ ] Energy: -2.2 ± 0.1 MeV
  - [ ] Norm: 0.99 ≤ norm ≤ 1.01
  - [ ] Nodes: 0 (for l=0)
  - [ ] Format: Displays numeric energy value
- [ ] Compare convergence history: Should monotonically approach root
- [ ] Test multi-channel case (M > 1): Check determinant calculation
- [ ] Regression test: Ensure existing valid results unchanged

---

## References

- Baye, D., et al. (1998). "Coupled-channel analysis of low-energy scattering using Lagrange mesh." *Nuclear Physics A* 640(1), 37-51.
- Press, W.H., et al. (2007). "Numerical Recipes". 3rd ed., Chapter 9: Root Finding.
- THOx Code Structure: See [README.md](README.md) for project overview

---

## Next Steps for User

1. **Review** the detailed analysis in `DEVEL/eigcc_rmat_issues.md`
2. **Integrate** the fix from `DEVEL/eigcc_rmat_fix_proposal.f90` into `eigcc_rmat.f90`
3. **Retest** with dni_cc.in and verify converged energy ≈ -2.22 MeV
4. **Validate** on multi-channel test cases
5. **Consider** additional safeguards for matrix conditioning at extreme energies

