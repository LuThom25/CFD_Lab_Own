# CFD Lab Worksheet 1 — Master Guide: Theory, Numerics & Code

> A comprehensive reference bridging the worksheet theory, the governing equations, and the C++ implementation. Read this alongside the code to understand both the *why* and the *how*.

---

## 1. What Problem Are We Solving?

We simulate **2D incompressible viscous flow** in a **lid-driven cavity** (LDC): a unit square domain filled with fluid, where the top wall moves at velocity $U_\text{lid} = 1\,\text{m/s}$ and the three remaining walls are fixed. The moving lid drags the fluid, creating a large recirculating vortex.

The LDC is the canonical CFD benchmark because:
- Its geometry is trivially simple (uniform Cartesian grid)
- Its physics is rich (vortex structure, corner eddies, Re-dependent transitions)
- Exact Ghia et al. (1982) reference data exists for validation

---

## 2. Governing Equations

### 2.1 Navier-Stokes Equations (Incompressible, 2D)

Two **momentum equations** (one per direction):

$$\frac{\partial u}{\partial t} + \frac{\partial(u^2)}{\partial x} + \frac{\partial(uv)}{\partial y} = -\frac{\partial p}{\partial x} + \nu\left(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}\right) + g_x \tag{1}$$

$$\frac{\partial v}{\partial t} + \frac{\partial(uv)}{\partial x} + \frac{\partial(v^2)}{\partial y} = -\frac{\partial p}{\partial y} + \nu\left(\frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2}\right) + g_y \tag{2}$$

One **continuity equation** (enforces incompressibility):

$$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0 \tag{3}$$

| Term | Physical meaning |
|------|-----------------|
| $\partial(\cdot)/\partial t$ | Local acceleration |
| $\partial(u^2)/\partial x + \partial(uv)/\partial y$ | Convective (nonlinear) transport |
| $-\partial p/\partial x$ | Pressure gradient force |
| $\nu\nabla^2 u$ | Viscous diffusion |
| $g_x, g_y$ | External body forces (gravity) |
| $\nu = \mu/\rho$ | Kinematic viscosity |

### 2.2 Knowns vs. Unknowns

| Category | Quantities |
|----------|-----------|
| **Knowns** (inputs) | Domain size $(L_x, L_y)$, viscosity $\nu$, lid velocity $U_\text{lid}$, initial fields $u_0=v_0=p_0=0$, boundary conditions |
| **Unknowns** (to solve) | Velocity field $u(x,y,t)$, $v(x,y,t)$, pressure field $p(x,y,t)$ at every grid cell and time step |

The Reynolds number $\text{Re} = U_\text{lid} \cdot L / \nu = 1/\nu$ (with $U=L=1$) is the key dimensionless parameter controlling flow behaviour.

---

## 3. Spatial Discretisation: The Staggered Grid

### 3.1 Why Staggered?

On a **collocated** grid (all variables at cell centres), the pressure-velocity coupling leads to *checkerboard instabilities* — spurious pressure oscillations that do not appear in the velocity divergence. A **staggered grid** avoids this by placing each variable at a different location:

```
  j+1  ·    v(i,j+1)   ·
       |        ↑       |
  j    u(i-1,j)→  p(i,j)  →u(i,j)
       |        ↑       |
  j-1  ·    v(i,j)     ·
              i        i+1
```

| Variable | Location | Indices |
|----------|----------|---------|
| $u(i,j)$ | Right face of cell $(i,j)$ | $i=0\ldots i_\text{max}$, $j=1\ldots j_\text{max}$ |
| $v(i,j)$ | Top face of cell $(i,j)$ | $i=1\ldots i_\text{max}$, $j=0\ldots j_\text{max}$ |
| $p(i,j)$ | Cell centre | $i=1\ldots i_\text{max}$, $j=1\ldots j_\text{max}$ |

### 3.2 Ghost Cells

The matrices have size $(i_\text{max}+2)\times(j_\text{max}+2)$. The outermost ring of cells ($i=0$, $i=i_\text{max}+1$, $j=0$, $j=j_\text{max}+1$) are **ghost cells** — they lie outside the physical domain and are used exclusively to enforce boundary conditions without special-casing the interior stencils.

**Example — no-slip bottom wall** ($j=0$ ghost layer):

The physical boundary is at $j=0$. No $u$-values lie on this horizontal boundary, so the no-slip condition $u=0$ is imposed by averaging:

$$u_{i,0} = -u_{i,1} \quad \Rightarrow \quad \frac{u_{i,0}+u_{i,1}}{2} = 0$$

The ghost cell value mirrors the interior value with a sign flip.

### 3.3 Finite Difference Stencils

**Laplacian** (central differences, 2nd order):

$$\nabla^2 A\big|_{i,j} = \frac{A_{i+1,j} - 2A_{i,j} + A_{i-1,j}}{\delta x^2} + \frac{A_{i,j+1} - 2A_{i,j} + A_{i,j-1}}{\delta y^2}$$

**Code:** `Discretization::laplacian(A, i, j)` in `src/Discretization.cpp`.

---

## 4. The Donor-Cell Convection Scheme

Standard central differences for convection terms ($\gamma=0$) are 2nd-order accurate but can produce non-physical oscillations at high Re. The **donor-cell (upwind) scheme** ($\gamma=1$) is stable but only 1st-order. The worksheet uses a blended version:

$$\frac{\partial(u^2)}{\partial x}\bigg|_{i,j} \approx \frac{u_e^2 - u_w^2}{\delta x} + \frac{\gamma}{\delta x}\left(\frac{|u_e|(u_i - u_{i+1})}{2} - \frac{|u_w|(u_{i-1} - u_i)}{2}\right)$$

where $u_e = (u_{i,j}+u_{i+1,j})/2$, $u_w = (u_{i-1,j}+u_{i,j})/2$.

- $\gamma = 0$: pure central differences (2nd order, less stable)
- $\gamma = 1$: full donor-cell upwind (1st order, maximum stability)
- $\gamma = 0.5$: practical compromise used in the base case

**Code:** `Discretization::convection_u()` and `convection_v()` — Equations (4) and (5) from the worksheet.

---

## 5. The Chorin Projection Algorithm

The algorithm splits each time step into three phases:

```
WHILE t < t_end:
  ┌─────────────────────────────────────────────────────┐
  │ Step 1: Apply velocity BCs (ghost cells, moving lid) │
  │         → fields.applyVelocity() for all boundaries  │
  ├─────────────────────────────────────────────────────┤
  │ Step 2: Compute intermediate fluxes F, G (Eq. 9,10) │
  │         F = u + dt*(ν∇²u - conv_u + gx)             │
  │         G = v + dt*(ν∇²v - conv_v + gy)             │
  │         → Fields::calculate_fluxes()                 │
  ├─────────────────────────────────────────────────────┤
  │ Step 3: Apply flux BCs (F=u on walls)                │
  ├─────────────────────────────────────────────────────┤
  │ Step 4: Compute RHS of pressure Poisson (Eq. 11)     │
  │         RS = (1/dt)*[(F_i-F_{i-1})/dx+(G_j-G_{j-1})/dy] │
  │         → Fields::calculate_rs()                     │
  ├─────────────────────────────────────────────────────┤
  │ Step 5: SOR pressure solve (Eq. 18)                  │
  │         WHILE iter < itermax AND res > eps:          │
  │           p = SOR_iteration(p, RS)                   │
  │           apply pressure BCs (Neumann)               │
  │         → SOR::solve()                               │
  ├─────────────────────────────────────────────────────┤
  │ Step 6: Velocity correction (Eq. 7,8)                │
  │         u^{n+1} = F - dt*(p_{i+1}-p_i)/dx           │
  │         v^{n+1} = G - dt*(p_{j+1}-p_j)/dy           │
  │         → Fields::calculate_velocities()             │
  ├─────────────────────────────────────────────────────┤
  │ Step 7: Compute adaptive dt for next step (Eq. 12,13)│
  │         dt = τ·min(dt_visc, dx/|u_max|, dy/|v_max|) │
  │         → Fields::calculate_dt()                     │
  ├─────────────────────────────────────────────────────┤
  │ Step 8: Advance time, write VTK output if due        │
  └─────────────────────────────────────────────────────┘
```

**Why this works:** Steps 2–4 produce an intermediate velocity field $(F,G)$ that satisfies the momentum equation but not necessarily the continuity equation. Step 5 finds a pressure correction that, when applied in Step 6, makes the updated velocities divergence-free (Eq. 3). This is the essence of the projection method.

---

## 6. Intermediate Fluxes F and G (Explicit Euler)

$$F_{i,j} = u_{i,j} + \delta t \left[\nu \nabla^2 u_{i,j} - \frac{\partial(u^2)}{\partial x}\bigg|_{i,j} - \frac{\partial(uv)}{\partial y}\bigg|_{i,j} + g_x\right] \tag{9}$$

$$G_{i,j} = v_{i,j} + \delta t \left[\nu \nabla^2 v_{i,j} - \frac{\partial(uv)}{\partial x}\bigg|_{i,j} - \frac{\partial(v^2)}{\partial y}\bigg|_{i,j} + g_y\right] \tag{10}$$

This is **explicit** in time: the right-hand side uses only values from the current time step $n$. No system of equations needs to be solved here.

**Code:** `Fields::calculate_fluxes()` in `src/Fields.cpp`.

---

## 7. Pressure Poisson Equation (PPE)

Substituting Equations (7,8) into the discrete continuity equation yields the **Pressure Poisson Equation**:

$$\frac{p_{i+1,j} - 2p_{i,j} + p_{i-1,j}}{\delta x^2} + \frac{p_{i,j+1} - 2p_{i,j} + p_{i,j-1}}{\delta y^2} = \underbrace{\frac{1}{\delta t}\left(\frac{F_{i,j} - F_{i-1,j}}{\delta x} + \frac{G_{i,j} - G_{i,j-1}}{\delta y}\right)}_{\text{RS}_{i,j}} \tag{11}$$

This is a **linear system** $\mathbf{A}\mathbf{p} = \mathbf{b}$, where $\mathbf{A}$ is the discrete Laplacian operator (sparse, symmetric). It is solved **implicitly** (the full pressure field at the new time step must be determined simultaneously).

### 7.1 Why SOR?

The system $\mathbf{A}\mathbf{p} = \mathbf{b}$ has $i_\text{max} \times j_\text{max}$ unknowns. Rather than a direct solver (expensive), we use **Successive Over-Relaxation (SOR)**, a stationary iterative method:

$$p_{i,j}^{(it+1)} = (1-\omega)\,p_{i,j}^{(it)} + \frac{\omega}{2\left(\frac{1}{\delta x^2}+\frac{1}{\delta y^2}\right)} \left(\frac{p_{i+1,j}^{(it)} + p_{i-1,j}^{(it+1)}}{\delta x^2} + \frac{p_{i,j+1}^{(it)} + p_{i,j-1}^{(it+1)}}{\delta y^2} - \text{RS}_{i,j}\right) \tag{18}$$

The **relaxation factor** $\omega \in (0,2)$ accelerates convergence:
- $\omega = 1$: Gauss-Seidel (no acceleration)
- $1 < \omega < 2$: over-relaxation (faster convergence for well-conditioned systems)
- $\omega \to 2$: instability / divergence

### 7.2 The Null-Space Problem (Pure Neumann BCs)

All boundary conditions for pressure are **Neumann** (zero gradient): $\partial p / \partial n = 0$. This means the system $\mathbf{A}\mathbf{p} = \mathbf{b}$ is **singular** — pressure is only determined up to an additive constant. Two independent mechanisms prevent convergence in the naive implementation:

1. **Fredholm incompatibility**: The discrete divergence $\sum \text{RS} \neq 0$ due to floating-point flux imbalances. The Fredholm alternative requires $\sum b = 0$ for a singular system $\mathbf{Ap} = \mathbf{b}$ to have *any* solution.
2. **Null-space drift**: Even if Fix 1 holds, each SOR sweep accumulates a constant offset in $p$ along the null space. This drift prevents the residual from decaying.

Without these fixes:
- The SOR residual $\|\nabla^2 p - \text{RS}\|$ **never reaches** `eps = 0.001`
- The solver always exhausts `itermax` iterations with avg\_res ≈ 1.2

### 7.3 Solution: Fredholm Compatibility + Zero-Mean Projection

**Fix 1 — RHS compatibility** (`Fields::calculate_rs()`, branch `ws1_further_extensions_improved_SOR`):  
After computing all RS values, subtract the mean to enforce $\sum \text{RS} = 0$:
```cpp
double sum = 0.0;
for (auto cell : cells) sum += _RS(cell->i(), cell->j());
const double mean = sum / static_cast<double>(N);
for (auto cell : cells) _RS(cell->i(), cell->j()) -= mean;
```

**Fix 2 — Zero-mean pressure projection** (`SOR::solve()`, after each sweep):  
Subtract the mean pressure from all fluid cells after every SOR iteration. This projects the iterate onto the subspace orthogonal to the null space, where the Poisson operator is invertible and the residual decays geometrically:
```cpp
double p_sum = 0.0;
for (auto cell : cells) p_sum += field.p(cell->i(), cell->j());
const double p_mean = p_sum / static_cast<double>(N);
for (auto cell : cells) field.p(cell->i(), cell->j()) -= p_mean;
```

**Result** (50×50, Re=100, $\omega=1.7$, itermax=100, $\varepsilon=10^{-3}$):

| Metric | Before fix | After fix |
|--------|-----------|-----------|
| avg SOR iterations | 100 (always) | **4.8** |
| avg residual | 1.230 | **3.4×10⁻³** |
| Converges to $\varepsilon$? | Never | ✅ (steps 2 onward) |

The fix is physically harmless: velocity fields depend only on pressure *gradients*, not absolute levels. The zero-mean normalization only removes the undetermined constant.

**Code:** `SOR::solve()` in `src/PressureSolver.cpp`, `calculate_rs()` in `src/Fields.cpp`.

---

## 8. Velocity Correction

Once the new pressure $p^{n+1}$ is known, the velocities are corrected:

$$u_{i,j}^{n+1} = F_{i,j} - \frac{\delta t}{\delta x}\left(p_{i+1,j}^{n+1} - p_{i,j}^{n+1}\right) \tag{7}$$

$$v_{i,j}^{n+1} = G_{i,j} - \frac{\delta t}{\delta y}\left(p_{i,j+1}^{n+1} - p_{i,j}^{n+1}\right) \tag{8}$$

**Code:** `Fields::calculate_velocities()`.

---

## 9. Stability & Adaptive Time Stepping

Three conditions must hold simultaneously for stability:

$$2\nu\,\delta t < \frac{(\delta x)^2(\delta y)^2}{(\delta x)^2 + (\delta y)^2}, \quad |u_\text{max}|\,\delta t < \delta x, \quad |v_\text{max}|\,\delta t < \delta y \tag{12}$$

The first is the **viscous (diffusion) stability condition**, the latter two are **CFL conditions**.

Adaptive time stepping selects:

$$\delta t = \tau \cdot \min\!\left(\frac{1}{2\nu}\left(\frac{1}{\delta x^2}+\frac{1}{\delta y^2}\right)^{-1},\ \frac{\delta x}{|u_\text{max}|},\ \frac{\delta y}{|v_\text{max}|}\right) \tag{13}$$

For the base case (50×50, $\nu=0.01$, $\tau=0.5$):

$$\delta t_\text{visc} = \frac{0.02^2}{4 \times 0.01} = 0.010 \text{ s}, \quad \delta t = 0.5 \times 0.010 = \mathbf{0.005\text{ s}}$$

The viscous limit (0.010 s) is **more restrictive** than the CFL limit ($\delta x / U_\text{lid} = 0.020$ s), so it dominates.

**Code:** `Fields::calculate_dt()` — returns `_dt` unchanged when `_tau ≤ 0` (fixed dt mode).

---

## 10. Boundary Conditions

### 10.1 No-Slip Walls (Fixed Walls)

Velocity zero on the boundary is enforced through ghost cells. For the **bottom wall** ($j=0$):

| Field | Condition | Ghost cell value |
|-------|-----------|-----------------|
| $u$ (horizontal) | $u = 0$ on boundary → average of ghost + interior = 0 | $u_{i,0} = -u_{i,1}$ |
| $v$ (vertical) | $v_{i,0} = 0$ directly (on boundary face) | — |
| $p$ | Neumann: $\partial p/\partial n = 0$ | $p_{i,0} = p_{i,1}$ |

**Code:** `FixedWallBoundary::applyVelocity()` and `applyPressure()` in `src/Boundary.cpp`.

### 10.2 Moving Lid (Top Wall)

The top wall moves at $U_\text{lid} = 1\,\text{m/s}$. The $u$-velocity at the boundary must equal $U_\text{lid}$:

$$\frac{u_{i,j_\text{max}} + u_{i,j_\text{max}+1}}{2} = U_\text{lid} \quad \Rightarrow \quad u_{i,j_\text{max}+1} = 2U_\text{lid} - u_{i,j_\text{max}}$$

**Code:** `MovingWallBoundary::applyVelocity()`.

---

## 11. Code Variable Mapping

| Physical quantity | C++ member | Class | Notes |
|------------------|-----------|-------|-------|
| $u_{i,j}$ | `_U(i,j)` | `Fields` | staggered right face |
| $v_{i,j}$ | `_V(i,j)` | `Fields` | staggered top face |
| $p_{i,j}$ | `_P(i,j)` | `Fields` | cell centre |
| $F_{i,j}$ | `_F(i,j)` | `Fields` | x-momentum flux |
| $G_{i,j}$ | `_G(i,j)` | `Fields` | y-momentum flux |
| $\text{RS}_{i,j}$ | `_RS(i,j)` | `Fields` | PPE right-hand side |
| $\nu$ | `_nu` | `Fields` | kinematic viscosity |
| $\delta t$ | `_dt` | `Fields` | current time step |
| $\tau$ | `_tau` | `Fields` | safety factor (≤0 → fixed dt) |
| $\gamma$ | `_gamma` | `Discretization` | donor-cell coefficient |
| $\delta x$, $\delta y$ | `_dx`, `_dy` | `Discretization` | (static members) |
| $\omega$ | `_omega` | `SOR` | relaxation factor |
| `itermax`, `eps` | `_max_iter`, `_tolerance` | `Case` | SOR stopping criteria |

All matrices use column-major storage: `_container[num_cols * j + i]`.

---

## 12. Key Insights & Common Pitfalls

| Insight | Explanation |
|---------|-------------|
| **Pressure is relative** | Pure Neumann BCs → only pressure *differences* are physical. Absolute value drifts but gradients stay correct. Per-frame rescaling is necessary for meaningful visualisation. |
| **SOR fix: Fredholm + zero-mean** | Two-step fix (branch `ws1_further_extensions_improved_SOR`): (1) subtract mean RS before iteration to satisfy Fredholm compatibility; (2) subtract mean pressure after each sweep to project out null-space drift. Result: avg\_sor drops from 100 → 4.8, avg\_res from 1.23 → 3.4×10⁻³. |
| **Viscous stability dominates** | At Re=100, $\delta t_\text{visc}=0.010 < \delta t_\text{CFL}=0.020$. Adaptive stepping automatically selects a safe margin of $\tau \cdot \delta t_\text{visc}=0.005$. |
| **Ghost cells before F/G** | `applyVelocity()` must be called *before* `calculate_fluxes()` so that boundary-adjacent stencils access correct ghost values. |
| **Fixed dt + fine grid = divergence** | For $\delta t=0.05$ and $\delta x < 0.05$, CFL>1. Instability grows exponentially. Use adaptive dt for all fine-grid runs. |
| **High Re → more SOR iterations** | The pressure equation becomes more ill-conditioned as Re increases. At Re≥1000, itermax=100 is insufficient for full convergence. |
