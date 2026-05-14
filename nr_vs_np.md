# Walkthrough - Debugging TELP Runtime Errors

I have resolved the runtime memory out-of-bounds error that occurred during the scattering calculation with the TELP potential enabled.

## Changes Made

### Debugging and Root Cause Identification
- Identified that the `wf` (wavefunction) and `y` (solution matrix) arrays in several solver subroutines in `scatcc.f90` were being dimensioned using the module variable `nr` from the `nmrv` module.
- Discovered that `nr` often contains a stale value from previous grid initializations when a subroutine is entered.
- When the subroutine updates `nr` to a new (larger) value `npt` passed as an argument, the subsequent loops exceed the bounds of the arrays that were fixed at entry based on the old `nr`.
- This led to the reported error: `Index '1202' of dimension 2 of array 'wf' above upper bound of 1201`.

### Memory Safety Fixes
- **Updated `scatcc.f90`**: Modified the declarations of `wf`, `y`, and other related arrays in the following subroutines to use the explicitly passed size argument (`npt` or `n1`) instead of the module variable `nr`:
    - `schcc`
    - `schcc_ena`
    - `schcc_erwinrc`
    - `schcc_MGR`
    - `schcc_erwin_MGR`
    - `schcc_ena_MGR`
    - `schcc_erwin_cuts`
    - `schcc_erwin_cuts_Q`
    - `schcc_ena_tolsma`
    - `matching3`
    - `matching3_eph`
    - `matching4`
    - `match_real`
- This ensures that the array bounds correctly match the actual number of points used in the calculation, regardless of the previous state of the module variables.

### Verification Results
- Compiled the code using `gfortran` with bounds checking enabled.
- Ran the example `examples/dni_cc.in`.
- **Outcome**: The calculation completed successfully without any runtime errors.
- **TELP Potential**: Verified that the `telp_pot.dat` file is correctly generated and contains valid radial data for the entire grid.

## Final Status
The single-pass, efficient TELP calculation is now stable and memory-safe across all integrated solvers.
