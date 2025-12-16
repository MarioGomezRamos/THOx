# R-matrix vs Numerov Method Debugging Report

## Summary

This report documents the investigation into why the R-matrix method (method=5,6) gives different results compared to the Numerov method (method=4) in THOx coupled-channel calculations.

## Key Finding

**The R-matrix method completely ignores the absorption (imaginary part of optical potential)**, while the Numerov method correctly handles it.

### Evidence from Be11+Au test:

| Method | |S| (elastic) | Absorption |
|--------|-------------|------------|
| Numerov (method=4) | 0.055 | 1.667 mb |
| R-matrix (method=5) | 0.9999 | 0.000 mb |

The same input file, same potentials - but completely different results.

## Tests Performed

### 1. Be11+Au System (test/Be11/)
- Input: `be11au_simple.in`
- Target: 197Au (Z=79), heavy nucleus with strong Coulomb
- Core potential: V0=-113 MeV, V0i=-50 MeV
- Valence potential: V0=-46.56 MeV, V0i=-5 MeV

**Results:**
- Numerov: Strong absorption, |S|≈0.055, unitarity≈0.016
- R-matrix: No absorption, |S|≈1.0, Abs=0

### 2. d+58Ni System (test/)
- Input: `dni_numerov.in`, `dni_rmatrix.in`
- Target: 58Ni (Z=28)
- Surface derivative absorption: V0i=-50 MeV

**Results:**
- Both methods give |S|=1.0, Abs=0
- This suggests a potential issue with how THOx handles absorption in general for this system

## Issues Discovered

### Issue 1: R-matrix ignores imaginary potential
The R-matrix solver receives the coupling potential `cpot` which should be complex (contains imaginary absorption part). However, the S-matrix always shows |S|≈1.0, indicating no absorption effect.

Verified that:
- `cpot` is declared as `complex*16`
- The imaginary part exists in fort.44 output file
- The `ch` matrix in R-matrix should be complex

### Issue 2: nlag convergence problem
When testing R-matrix with different `nlag` values (Lagrange points per interval):
- ns=1, nlag=100: |S|≈0.27
- ns=1, nlag=200: |S|≈1.13
- ns=1, nlag=400: |S|≈2.07
- ns=1, nlag=600: |S|≈2.65

**Results diverge instead of converging!** This is abnormal behavior.

### Issue 3: d+58Ni shows no absorption even with Numerov
When testing with d+58Ni system:
- Even with V0i=-50 MeV (strong absorption)
- Even with Numerov method
- |S|=1.0 everywhere, Abs=0

This suggests either:
1. The fragment-target potential imaginary part is not being used in CC calculations
2. Or there's a configuration issue with the input file

## Files Modified/Created

1. `test/Be11/be11au_simple.in` - Simplified Be11+Au test
2. `test/dni_numerov.in` - d+58Ni test with Numerov
3. `test/dni_rmatrix.in` - d+58Ni test with R-matrix

## Code Locations Investigated

- `/Users/jinlei/Desktop/code/THOx/scatcc.f90` - Main CC solver
  - Line 15: `vcoup` declared as `complex*16`
  - Line 685: `schcc_rmat` subroutine (Pierre's R-matrix)
  - Line 6097: `schcc_rmat_MGR` subroutine
  - Line 6286: `schcc_rmat_hp_MGR` subroutine (HPRMAT)

- `/Users/jinlei/Desktop/code/THOx/solvecc.f90` - CC equation driver
  - Line 673: calls `schcc_rmat`
  - Line 773: calls `schcc_rmat_MGR`

## Possible Root Causes

1. **R-matrix Hamiltonian construction**: The imaginary part of the potential may not be correctly incorporated into the R-matrix Hamiltonian matrix.

2. **Sign convention**: Numerov uses `-conv*vcoup` while R-matrix may use `+conv*vcoup`. A sign test was attempted but made results worse.

3. **THOx formfactor generation**: The imaginary part may be lost during formfactor interpolation before reaching CC solver.

## Next Steps

1. Add debug prints in `schcc_rmat` to verify `cpot` has non-zero imaginary part
2. Check if the `ch` matrix (R-matrix Hamiltonian) contains imaginary components
3. Investigate why d+58Ni shows no absorption even with Numerov - this may reveal a separate bug in THOx
4. Compare with standalone HPRMAT code which is known to work correctly

## Conclusion

The R-matrix implementation in THOx appears to have a fundamental issue with handling absorptive (complex) potentials. The absorption effect is completely missing, resulting in |S|=1.0 instead of |S|<1 as expected for optical model calculations. Additionally, the nlag parameter shows divergent behavior instead of convergence, suggesting numerical instability in the R-matrix solver.
