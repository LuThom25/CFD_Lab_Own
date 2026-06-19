# PCG_SSOR Pressure Solver — Strong Scaling Analysis

**Branch:** `project_precondioned_cg_analysis`
**Date:** 2026-06-19
**Solvers compared:** SOR_STANDARD · CG_STANDARD · PCG_SSOR

---

## 1. Overview

This analysis extends the CG strong scaling study by adding a Preconditioned Conjugate Gradient solver
with a Red-Black SSOR preconditioner (PCG_SSOR). All four test cases from the previous study are retained,
and every case is run with three MPI decompositions each (serial / 1×1 / 2×2 / 4×1 or 1×4).

### Test Cases

| Case | Grid | t_end | Decompositions |
|------|------|------:|----------------|
| Rayleigh-Bénard convection | 40×18 | 10 000 | 1×1, 2×2, 4×1 |
| Fluid Trap (small)         | 100×50 | 2 000 | 1×1, 2×2, 4×1 |
| Fluid Trap (large)         | 200×100 | 2 000 | 1×1, 2×2, 4×1 |
| Lid-Driven Cavity          | 150×150 | 50 | serial, 1×1, 2×2, 1×4 |

---

## 2. Algorithm: PCG_SSOR

The Preconditioned Conjugate Gradient method replaces the standard search direction `d = r` with
`d = M⁻¹r`, where `M` is a symmetric positive definite preconditioner approximating `A`.

### Why Red-Black SSOR?

PCG requires `M = Mᵀ` (symmetric preconditioner) for guaranteed convergence. Standard SSOR with MPI
domain decomposition produces domain-dependent lower triangular factors `L`, making `M` non-symmetric
across decompositions. **Red-Black coloring** uses global coordinates `(i+j) % 2` — consistent across
all MPI ranks — which guarantees `M = Mᵀ` regardless of the decomposition.

### Preconditioner Application (`apply_ssor`)

Each call to `M⁻¹r` performs four passes (two forward, two backward):

```
Initialize z = 0
Forward red:   z_red   = (r_red   + offdiag(z_black)) / d_ii  →  communicate(z)
Forward black: z_black = (r_black + offdiag(z_red))   / d_ii  →  communicate(z)
Backward black: z_black = (r_black + offdiag(z_red))  / d_ii  →  communicate(z)
Backward red:   z_red   = (r_red   + offdiag(z_black)) / d_ii →  communicate(z)
```

This yields a symmetric Gauss-Seidel preconditioner. The backward black pass does not change `z_black`
(since `z_red` is unchanged at that point), but is retained for symmetry of the communication pattern.

### PCG Iteration (per timestep)

```
r ← Δₕp − RS         (warm start from current p)
z ← M⁻¹r             (apply_ssor)
d ← z
ρ = (r, z)            (M⁻¹-weighted inner product)

for k = 0, 1, ...:
    communicate(d)
    q ← −Δₕd
    α = ρ / (d, q)
    p ← p + α·d
    r ← r − α·q
    z ← M⁻¹r          (apply_ssor)
    ρ_new = (r, z)
    β = ρ_new / ρ
    d ← z + β·d        ← key difference from plain CG: uses z, not r
    ρ ← ρ_new
```

**Theoretical advantage:** Condition number κ(M⁻¹A) = O(√N) instead of O(N) for plain CG,
reducing iterations from O(N^½) to O(N^¼).

---

## 3. Wall-Time Results

### 3.1 Rayleigh-Bénard (40×18)

| Config | SOR | CG | PCG_SSOR |
|--------|----:|---:|--------:|
| 1×1    | 21.7s | 19.8s | 24.5s |
| 2×2    |  8.9s |  9.5s | 18.4s |
| 4×1    |  9.8s |  9.7s | 17.3s |

→ See `pictures/plot_rb_walltime.pdf`, `plot_rb_speedup.pdf`

### 3.2 Fluid Trap 100×50

| Config | SOR | CG | PCG_SSOR |
|--------|----:|---:|--------:|
| 1×1    | 5.2s | 7.3s | **4.6s** |
| 2×2    | 2.3s | 3.7s | **2.2s** |
| 4×1    | 2.7s | 3.3s | **2.4s** |

→ See `pictures/plot_fts100_walltime.pdf`, `plot_fts100_speedup.pdf`

**PCG_SSOR is the fastest solver on this case.**

### 3.3 Fluid Trap 200×100

| Config | SOR | CG | PCG_SSOR |
|--------|----:|---:|--------:|
| 1×1    | 27.2s | 68.1s | 55.0s |
| 2×2    |  8.7s | 37.1s | 27.6s |
| 4×1    | 10.0s | 38.4s | 30.2s |

→ See `pictures/plot_ft_walltime.pdf`, `plot_ft_speedup.pdf`

PCG improves over CG but cannot beat SOR for this geometry/resolution.

### 3.4 Lid-Driven Cavity 150×150

| Config | SOR | CG | PCG_SSOR |
|--------|----:|---:|--------:|
| Serial |  185s | 833s | — |
| 1×1    |  185s | 873s | 966s |
| 2×2    |   65s | 313s | 457s |
| 1×4    |   61s | 290s | 409s |

→ See `pictures/plot_ldc_walltime.pdf`, `pictures/plot_ldc_speedup.pdf`

---

## 4. Iteration Counts

| Case | CG avg | PCG avg | Reduction |
|------|-------:|--------:|----------:|
| Rayleigh-Bénard 40×18   | 0.6 | 0.5 | 1.2× |
| Fluid Trap 100×50       | 1.4 | 0.7 | 2.0× |
| Fluid Trap 200×100      | 1.5 | 0.9 | 1.7× |
| LDC 150×150             | 2.0 | 0.6 | **3.3×** |

→ See `pictures/plot_ldc_all_solvers_1_1.pdf` and corresponding plots for all cases.

The preconditioner demonstrably reduces iteration counts. For LDC, PCG converges in 3.3× fewer
iterations than CG. However, each PCG iteration is significantly more expensive.

---

## 5. Why PCG_SSOR Underperforms Despite Fewer Iterations

### 5.1 Communication overhead per iteration

| Solver | `communicate_field` calls/iter | `reduce_sum` calls/iter |
|--------|-------------------------------:|------------------------:|
| SOR    | 1                              | 1 (residual)            |
| CG     | 1                              | 2                       |
| PCG_SSOR | **4** (in apply_ssor, called twice per iter = 8 total) | **4** |

Each `communicate_field` is a blocking MPI halo exchange. PCG needs `apply_ssor` once before the loop
and once per iteration, totalling **8 communicate_field + 4 reduce_sum** per outer iteration vs.
CG's **1 + 2**. The communication overhead outweighs the iteration savings for all tested grid sizes.

### 5.2 Break-even analysis

For PCG to beat CG in wall time, the iteration reduction must exceed the communication overhead factor:

```
Required iteration reduction > (PCG comm cost) / (CG comm cost) ≈ 8 / 1 = 8×
Observed reduction (best case, LDC): 3.3×
```

The theoretical maximum reduction for SSOR preconditioning on 2D Poisson is O(N^¼ / N^½) = O(N^{-¼}).
For a 150×150 grid (N ≈ 22 500), this yields at most ~7× — barely sufficient, and only achievable
with optimal ω and a perfect preconditioner.

### 5.3 Scaling behavior

| Case | SOR speedup 1×1→4R | CG speedup | PCG speedup |
|------|-------------------:|-----------:|------------:|
| RB   | 2.2–2.4×           | 2.0–2.1×   | 1.3–1.4× |
| FTs  | 1.9–2.3×           | 2.0–2.2×   | 1.9–2.1× |
| FT   | 2.7–3.1×           | 1.8×       | 1.8–2.0× |
| LDC  | 2.8–3.0×           | 2.7–2.9×   | **2.1–2.4×** |

PCG scales worst due to its higher communication cost per iteration:
more MPI calls per iteration → more time spent waiting at barriers → lower parallel efficiency.

---

## 6. Known Implementation Note

The backward-black sweep in `apply_ssor` (3rd of 4 passes) does not change `z_black` because
`z_red` has not been updated yet at that point — it reads the same values as the forward-black pass.
This sweep, together with its `communicate_field`, is therefore a **no-op in terms of numerical
result**, but it preserves the symmetric communication structure. Removing it would reduce communicates
from 8 to 6 per outer iteration and save ~25% of the MPI overhead, without changing the computed
preconditioned direction `z`.

---

## 7. Conclusions

1. **PCG_SSOR converges in fewer iterations** in all cases. The preconditioner is mathematically
   correct (M = Mᵀ guaranteed by global red-black coloring) and PCG converges.

2. **PCG_SSOR wins only for small grids** (Fluid Trap 100×50), where the iteration reduction
   outweighs the communication overhead.

3. **For large grids with MPI**, the 4× higher communication cost per `apply_ssor` call dominates.
   PCG does not beat SOR or CG in wall time for Rayleigh-Bénard, 200×100 Fluid Trap, or 150×150 LDC.

4. **SOR_STANDARD remains the fastest solver** for all tested cases in MPI parallel runs. Its low
   iteration cost per step and single halo exchange per iteration outperform the theoretically
   superior CG-based solvers on these grid sizes.

5. **Potential improvements** (not implemented):
   - Remove redundant 3rd sweep in `apply_ssor` (6 comm/iter instead of 8)
   - Use ω > 1 (over-relaxation) for stronger preconditioning with same communication pattern
   - Replace convergence check from `||r||` to `(r,z)` to eliminate one `reduce_sum` per iteration
   - Polynomial (Chebyshev) preconditioner: same iteration reduction, zero extra MPI calls

---

## 8. Pictures Index

### Iterations per timestep
| File | Content |
|------|---------|
| `plot_rb_sor_iters.pdf` | SOR iterations — Rayleigh-Bénard all decomps |
| `plot_rb_cg_iters.pdf`  | CG iterations — Rayleigh-Bénard all decomps |
| `plot_rb_pcg_iters.pdf` | PCG iterations — Rayleigh-Bénard all decomps |
| `plot_rb_all_solvers_1_1.pdf` | All 3 solvers, 1×1 — Rayleigh-Bénard |
| `plot_rb_all_solvers_2_2.pdf` | All 3 solvers, 2×2 — Rayleigh-Bénard |
| `plot_rb_all_solvers_4_1.pdf` | All 3 solvers, 4×1 — Rayleigh-Bénard |
| `plot_fts100_*_iters.pdf` | Same for Fluid Trap 100×50 |
| `plot_fts100_all_solvers_*.pdf` | Per-decomp comparison, Fluid Trap 100×50 |
| `plot_ft_*_iters.pdf` | Same for Fluid Trap 200×100 |
| `plot_ft_all_solvers_*.pdf` | Per-decomp comparison, Fluid Trap 200×100 |
| `plot_ldc_*_iters.pdf` | Same for LDC 150×150 |
| `plot_ldc_all_solvers_*.pdf` | Per-decomp comparison, LDC 150×150 |

### Wall time & Speedup
| File | Content |
|------|---------|
| `plot_rb_walltime.pdf`     | 3-solver bar chart — Rayleigh-Bénard |
| `plot_rb_speedup.pdf`      | Speedup lines — Rayleigh-Bénard |
| `plot_fts100_walltime.pdf` | 3-solver bar chart — Fluid Trap 100×50 |
| `plot_fts100_speedup.pdf`  | Speedup lines — Fluid Trap 100×50 |
| `plot_ft_walltime.pdf`     | 3-solver bar chart — Fluid Trap 200×100 |
| `plot_ft_speedup.pdf`      | Speedup lines — Fluid Trap 200×100 |
| `plot_ldc_walltime.pdf`    | 3-solver bar chart — LDC 150×150 |
| `plot_ldc_speedup.pdf`     | Speedup lines — LDC 150×150 |
