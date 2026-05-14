# Simultaneous Core and Valence Excitation for $B(E\lambda)$

This plan provides the mathematical framework and implementation strategy to integrate the "Simultaneous Excitation" cross-terms into the Electric Multipole transition strength calculation.

Currently, the $B(E\lambda)$ operator treats core excitation and valence excitation as strictly separate one-body operators. Consequently, matrix elements where BOTH core and valence components change states evaluate to zero due to angular orthogonality. 

Theoretical analysis of the project documentation (Eq. 830) confirms that shifting the total multipole operator into relative coordinates yields coherent mixed product terms arising from the solid harmonic expansion:
$$ \hat{\mathcal{M}}_{mixed} = \sum_{k=1}^{\lambda-1} \sqrt{\frac{4\pi (2\lambda+1)!}{(2k+1)!(2(\lambda-k)+1)!}} \left(-\frac{A_v}{A_{tot}}\right)^{\lambda-k} \cdot r^{\lambda-k} [ Y_{\lambda-k}(\hat{r}) \otimes \hat{\mathcal{M}}_{core}(Ek) ]_{\lambda} $$
These terms explicitly connect states where both the core and valence spaces change configuration simultaneously (crucial for octupole $\lambda=3$ and higher multipoles).

## User Review Required

> [!IMPORTANT]
> **Core Model Expansion**: The simultaneous excitation routine will look up core matrix elements $\langle I_f \| \hat{\mathcal{M}}_{core}(Ek) \| I_i \rangle$ for intermediate multipoles $k \in [1, \lambda-1]$. 
> - **CoreModel=1 (RME)**: Requires the user to have supplied matrix elements for those intermediate multipolarities in the user-defined matrix elements array. If missing, they revert to zero.
> - **CoreModel=0 (Rotor)**: Will automatically calculate the internal core moments using the existing scaling relation $M(Ek) \propto \delta (rms)^{k-1}$ applied to the requested intermediate $k$.

## Proposed Changes

### 1. Extend Mathematical Utilities

#### [MODIFY] [belam.f90](file:///home/amoro/cloud/programs/thox/github/belam.f90)
Create a new helper function `matmix` at the bottom of the file to evaluate the coupled angular matrix element of the tensor product of operators acting on Valence and Core spaces.
The 9-j formulation ensures consistent angular conservation during simultaneous redistribution:
```fortran
function matmix(k, lambdav, lambda, sn, li, jni, jtoti, ici, lf, jnf, jtotf, icf)
  ! Incorporates:
  ! 1. 9-j coupling between Valence and Core space dynamics.
  ! 2. 6-j coupling reduction for the valence operator into the inert-spin frame.
  ! 3. Standard 3-j spherical harmonic elements for Valence parity change.
```

### 2. Implement Matrix Element Accumulation

#### [MODIFY] [belam.f90](file:///home/amoro/cloud/programs/thox/github/belam.f90)
Introduce a tertiary accumulation loop within the discrete and continuum channel iterations for coherent addition of the mixed terms.

1.  **Discrete States Loop** (inside `do m=1,nchan`, `do n=1,ncni`):
    -   Detect transition candidates that may vary both core and valence quantum numbers.
    -   Apply selection guards for coupling: `|ji-jf| <= lambdav <= ji+jf` AND `|Ii-If| <= k <= Ii+If`.
    -   Loop over $k=1 \dots \lambda-1$:
        -   Obtain core moment $M(Ek)$ from `mec` (Model 1) or scale `Elamcorem` (Model 0).
        -   Perform radial integral $\langle u_f | r^{\lambda-k} | u_i \rangle$.
        -   Evaluate Angular 9j element via `matmix`.
        -   Accumulate to cumulative matrix element `Elam`.
2.  **Continuum Scat-States Loop**:
    -   Replicate the identical inner logic block into the iterative continuous iteration for parity.

## Verification Plan

### Automated Tests
-   Compile functionality check via `make clean && make`.
-   Execute standard dipole test-case to guarantee that simultaneous contributions naturally default to zero without overhead or regression for $\lambda=1$.

### Manual Verification
-   Enable Octupole scenario ($\lambda=3$) with Quadrupole core configuration to empirically observe population of the simultaneous transition channel strength, which is mathematically prevented by the current legacy guards.
