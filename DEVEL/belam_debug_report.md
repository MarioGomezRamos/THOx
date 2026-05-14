# Debug Report: Simultaneous Valence + Core Excitation in `belam.f90`

The simultaneous (mixed) term appears in **two places**:
- **Discrete** B(Eλ): lines 415–459
- **Continuum** B(Eλ): lines 631–673

Both blocks are structurally identical; every bug appears twice.

---

## Bug 1 — Critical: Spurious contributions when `ikterm + ilv ≠ lambda`

### Location
- Discrete: lines 444–449
- Continuum: lines 659–664

### Code (identical in both sections)
```fortran
if (ikterm + ilv == lambda) then
   term_pref = sqrt(4d0 * pi) * dfactratio(lambda, ikterm) &
             * ((-av)/(ac+av))**(ilv)
else
   term_pref = 1.0d0    ! <-- BUG
endif
```

### Why it's wrong
The mixed term arises from the expansion of `r_rel^λ Y_λ(r_rel)` into products of
core (`r_c^k Y_k`) and valence (`r_v^lv Y_lv`) operators. **The only valid expansion
terms satisfy `k + lv = λ` exactly** — this is a property of the multipole expansion
in CM coordinates, not a triangle condition.

The `fail3` guard above admits all `(k, lv)` satisfying the *triangle* condition
`|k−lv| ≤ λ ≤ k+lv`. So for `λ=2, k=1, lv=2` the triangle passes but `k+lv=3≠2`.
The code then falls into the `else` branch and sets `term_pref = 1.0d0`, adding a
completely unphysical contribution.

### Fix
Replace the `else` clause with `cycle`:

```fortran
if (ikterm + ilv == lambda) then
   term_pref = sqrt(4d0 * pi) * dfactratio(lambda, ikterm) &
             * ((-av)/(ac+av))**(ilv)
else
   cycle   ! k+lv ≠ λ: these pairs are absent from the multipole expansion
endif
```

---

## Bug 2 — Stray debug `write` left in the continuum loop

### Location
- Continuum only: line 668

### Code
```fortran
write(0,*)'term_mix=',term_mix
```

### Why it's wrong
This is inside the quadruple loop `do iil / do ik / do m / do n` and the double
`do ikterm / do ilv` loop.  For any non-trivial run it prints **millions of lines to
stderr**, making the calculation unusable. It was clearly left from a debugging session.

### Fix
Remove the line entirely.

---

## Bug 3 — Guard condition `lambda.gt.0` is too permissive

### Location
- Discrete: line 416
- Continuum: line 632

### Code
```fortran
if (lambda.gt.0) then
   do ikterm = 1, maxlamb
     do ilv = 1, maxlamb
```

### Why it matters
For `λ=1` (E1), both `ikterm` and `ilv` start at 1, so `ikterm+ilv ≥ 2 > 1 = λ`.
After fixing Bug 1, **every iteration would cycle immediately** — the block does zero
useful work for E1.  With Bug 1 still present it is actively harmful: the pair
`(k=1, lv=1)` satisfies the triangle condition for `λ=1` and triggers `term_pref=1.0`.

### Fix
Change to `lambda.gt.1` (mixed term only exists for `λ ≥ 2`):

```fortran
if (lambda.gt.1) then
```

---

## Bug 4 — Missing `cindexi(n) == cindex(m)` guard in the mixed term

### Location
- Discrete: ~line 415 (outer context)
- Continuum: ~line 631

### Why it matters
The **pure valence** contribution correctly guards:
```fortran
if (.not.fail3(qji(n),lambdar,qj(m)) .and. (cindexi(n).eq.cindex(m))) then
```
The `cindexi(n).eq.cindex(m)` condition avoids coupling channels that share the same
J/π quantum numbers but belong to different core states.

The **mixed term** has no such guard. While the core selection rules (`fail3` on the
core angular momenta) partially cover this, it is possible for `rkcore` to be
non-zero even when the core states are unphysical partners (e.g. in `coremodel=1`
where `mec` can be set for any pair). Without the guard, matrix elements are computed
between unphysical channel combinations.

### Fix
Add the guard before the `ikterm/ilv` loops:

```fortran
! Simultaneous Core and Valence contribution
if (lambda.gt.1 .and. cindexi(n).eq.cindex(m)) then
  do ikterm = 1, maxlamb
    ...
```

> [!NOTE]  
> This guard may be intentionally absent if the physics requires coupling channels with
> different core indices (e.g. fragmented core states). Confirm with your physics intent
> before applying. The pure valence block immediately above applies it, so consistency
> suggests it should be here too.

---

## Summary Table

| # | Severity | Location | Description |
|---|----------|----------|-------------|
| 1 | **Critical** | Lines 444–449 & 659–664 | `term_pref = 1.0d0` for `k+lv≠λ` → spurious contributions |
| 2 | **Critical** | Line 668 | Stray `write(0,*)` inside inner loop → unusable output |
| 3 | Medium | Lines 416 & 632 | Guard `lambda.gt.0` should be `lambda.gt.1` |
| 4 | Medium | Lines 415 & 631 | Missing `cindexi(n).eq.cindex(m)` guard |

---

## Recommended Patch (both sections — discrete shown, continuum identical)

```diff
- if (lambda.gt.0) then
+ if (lambda.gt.1 .and. cindexi(n).eq.cindex(m)) then
    do ikterm = 1, maxlamb
     do ilv = 1, maxlamb
       ...
       if (abs(rkcore).lt.1d-10) cycle
       faux(1:nr) = rvec(1:nr)**(1 + ilv) &
                  * ugs(1:nr,n) * wfc(jset,i,m,1:nr)
       call simc(faux, resc, 1, nr, dr, nr)
       rlam = resc
       if (ikterm + ilv == lambda) then
          term_pref = sqrt(4d0 * pi) * dfactratio(lambda, ikterm) &
                    * ((-av)/(ac+av))**(ilv)
       else
-         term_pref = 1.0d0
+         cycle
       endif
       term_mix = matmix(...)
       Elamcore = term_pref * term_mix * rlam * rkcore
       Elam = Elam + Elamcore
       ...
     enddo
    enddo
  endif
```

For the **continuum** block, additionally remove line 668:

```diff
- write(0,*)'term_mix=',term_mix
```
