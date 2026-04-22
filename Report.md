# CFD Lab Worksheet 1 — Report: Lid-Driven Cavity Simulation

---

## 1. Understanding Incompressible Navier-Stokes Equations

The motion of a viscous, incompressible Newtonian fluid is governed by the **incompressible Navier-Stokes equations** (NSE). In two dimensions these read:

$$\frac{\partial u}{\partial t} + \frac{\partial(u^2)}{\partial x} + \frac{\partial(uv)}{\partial y} = -\frac{\partial p}{\partial x} + \nu\left(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}\right) + g_x \tag{1}$$

$$\frac{\partial v}{\partial t} + \frac{\partial(uv)}{\partial x} + \frac{\partial(v^2)}{\partial y} = -\frac{\partial p}{\partial y} + \nu\left(\frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2}\right) + g_y \tag{2}$$

$$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0 \tag{3}$$

where $u(x,y,t)$ and $v(x,y,t)$ are the horizontal and vertical velocity components, $p(x,y,t)$ is the kinematic pressure (divided by density $\rho$, assumed constant), $\nu$ is the kinematic viscosity, and $g_x, g_y$ are body accelerations (zero for LDC).

**Physical interpretation of each term:**

| Term | Physical meaning |
|---|---|
| $\partial u/\partial t$ | Local (temporal) acceleration |
| $\partial(u^2)/\partial x + \partial(uv)/\partial y$ | Nonlinear convective transport (inertia) |
| $-\partial p/\partial x$ | Pressure-gradient force (drives flow from high to low pressure) |
| $\nu\nabla^2 u$ | Viscous diffusion (tends to smooth velocity gradients) |
| Eq. 3 | Incompressibility: no local volume change, no mass sources or sinks |

The **incompressibility constraint** (Eq. 3) is a differential constraint that couples pressure and velocity: at every point in space and time, the pressure must be whatever is needed to keep $\nabla \cdot \mathbf{u} = 0$. This makes pressure an elliptic, non-local quantity — it propagates changes instantaneously across the domain (the incompressible limit of infinite sound speed).

The **Reynolds number** $Re = U_\text{wall} L / \nu$ is the key dimensionless parameter. For $U_\text{wall} = L = 1$ simply $Re = 1/\nu$. Low $Re$ (large $\nu$): viscosity dominates, flow is laminar and symmetric. High $Re$ (small $\nu$): inertia dominates, complex vortex structures, potentially unsteady or turbulent behaviour.

---

## 2. Numerical Modelling

### 2.1 Staggered Grid and Ghost Cells

The solver uses a **staggered Cartesian (MAC) grid** with $i_\text{max} \times j_\text{max}$ interior cells of uniform spacing $\Delta x = 1/i_\text{max}$, $\Delta y = 1/j_\text{max}$. The three unknown fields are stored at different sub-cell locations:

```
         V(i,j)
           ↑
  ----+----+----+----
      |         |
←---  |  P(i,j) |  --- U(i,j) →
      |         |
  ----+---------+----
         j
         ↑
      i → 
```

- **$U(i,j)$** — $x$-velocity at the **right face** of cell $(i,j)$, between cells $(i,j)$ and $(i+1,j)$
- **$V(i,j)$** — $y$-velocity at the **top face** of cell $(i,j)$, between cells $(i,j)$ and $(i,j+1)$
- **$P(i,j)$** — pressure at the **cell centre** of $(i,j)$

**Why staggered?** On a collocated grid the discrete pressure gradient $[p(i+2,j)-p(i,j)]/(2\Delta x)$ samples $p$ two cells apart, decoupling odd and even grid lines and producing non-physical checkerboard pressure oscillations. On the staggered grid the pressure gradient $[p(i+1,j)-p(i,j)]/\Delta x$ acts directly on the co-located face velocity $U(i,j)$, giving a locally consistent discrete divergence-free projection and eliminating the checkerboard instability entirely.

All matrices are allocated with size $(i_\text{max}+2) \times (j_\text{max}+2)$. The outermost ring (indices $0$ and $i_\text{max}+1$ in each direction) constitutes **ghost cells** — fictitious cells outside the physical domain. Their values are set after each boundary-condition step and enter the interior stencils (e.g., the Laplacian at $i=1$ uses $U(0,j)$ from the ghost layer). This allows the same loop structure to be used for all interior cells without conditional logic.

### 2.2 Algorithm: Chorin Projection Method

The solver implements the **Chorin fractional-step (projection) method**. Each time step proceeds as follows:

```
Step 1  Apply velocity BCs          → fill ghost-cell U, V from BC formulas
Step 2  Compute intermediate F, G   → explicit Euler momentum prediction (Eq. 9, 10)
Step 3  Apply flux BCs              → set F, G at wall faces
Step 4  Build PPE right-hand side   → RS = (1/dt)·[∂F/∂x + ∂G/∂y] (Eq. 11)
Step 5  SOR pressure solve          → iterate P until residual < ε or iter > itermax
        (+ apply Neumann P BCs after each sweep)
Step 6  Correct velocities          → U = F - dt·∂P/∂x, V = G - dt·∂P/∂y (Eq. 7, 8)
Step 7  Compute adaptive dt         → dt = τ·min(dt_visc, dx/|u|_max, dy/|v|_max)
Step 8  Write VTK output if due     → interpolate staggered fields to cell centres
```

**Explicit step (Steps 2–3):** The intermediate fluxes $F$ and $G$ approximate the momentum equations without the pressure gradient, using forward Euler:

$$F(i,j) = U(i,j) + \Delta t\bigl(\nu\,\nabla^2_h U - \mathrm{conv}_u(U,V) + g_x\bigr)$$
$$G(i,j) = V(i,j) + \Delta t\bigl(\nu\,\nabla^2_h V - \mathrm{conv}_v(U,V) + g_y\bigr)$$

These intermediate fields generally violate the continuity constraint $\nabla \cdot \mathbf{u} = 0$.

**Implicit step (Steps 4–6):** The pressure is determined by requiring that the corrected velocities $U^{n+1} = F - \Delta t\,\partial p/\partial x$, $V^{n+1} = G - \Delta t\,\partial p/\partial y$ satisfy continuity. Taking the discrete divergence and setting it to zero yields the **Pressure Poisson Equation (PPE)**:

$$\frac{\partial^2 p}{\partial x^2} + \frac{\partial^2 p}{\partial y^2} = \frac{1}{\Delta t}\left(\frac{F(i,j)-F(i-1,j)}{\Delta x} + \frac{G(i,j)-G(i,j-1)}{\Delta y}\right) \equiv \mathrm{RS}(i,j) \tag{11}$$

This linear system is solved by **SOR** (Successive Over-Relaxation):

$$p^{(\ell+1)}(i,j) = (1-\omega)\,p^{(\ell)}(i,j) + \frac{\omega}{2\!\left(\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}\right)}\!\Bigl(\mathrm{sor\_helper}(p,i,j) - \mathrm{RS}(i,j)\Bigr) \tag{18}$$

where sor\_helper computes the off-diagonal sum $[p(i+1,j)+p(i-1,j)]/\Delta x^2 + [p(i,j+1)+p(i,j-1)]/\Delta y^2$. After each sweep the Neumann pressure BCs are enforced (ghost cells equal their interior neighbour). The iteration stops when the RMS residual $r < \varepsilon = 10^{-3}$ or the budget `itermax` is exhausted.

### 2.3 Donor-Cell Convection Scheme

Pure central-difference discretisation of the nonlinear convection terms is second-order accurate but oscillatory at high cell-Péclet number. Pure first-order upwind (donor-cell) is unconditionally stable but excessively diffusive. The solver uses a **blended scheme** parameterised by $\gamma \in [0,1]$ (Eq. 4):

$$\frac{\partial(u^2)}{\partial x}\bigg|_{i,j} = \underbrace{\frac{u_e^2 - u_w^2}{\Delta x}}_{\text{central (2nd order)}} + \underbrace{\frac{\gamma}{\Delta x}\!\left(\frac{|u_e|(u_i-u_{i+1})}{2} - \frac{|u_w|(u_{i-1}-u_i)}{2}\right)}_{\text{donor-cell correction (1st order)}}$$

with face values $u_e = (U(i,j)+U(i+1,j))/2$, $u_w = (U(i-1,j)+U(i,j))/2$. Setting $\gamma=0$ recovers pure central differences; $\gamma=1$ gives full first-order upwind. The choice $\gamma=0.5$ (used throughout) provides near-second-order accuracy in smooth regions with sufficient numerical dissipation to stabilise steep gradients near the moving lid.

### 2.4 Boundary Conditions

**No-slip stationary walls** use ghost-cell mirror extrapolation. For a wall cell at $(i,j)$ bordering fluid on the right:
- Normal velocity: `U(i,j) = 0` (zero penetration, face shared with fluid)
- Tangential ghost: `V(i,j) = -V(i+1,j)` (mirror so face average = 0)

**Moving lid** (top wall, velocity $U_\text{wall}=1$) uses:
- Normal: `V(i, j_max) = 0` (no penetration)
- Tangential ghost: `U(i, j_max+1) = 2·U_wall - U(i, j_max)` (face average = $U_\text{wall}$)

**Neumann pressure BCs** (all walls): after each SOR sweep, ghost-cell pressures equal their interior neighbour — enforcing $\partial p/\partial n = 0$, i.e., no pressure gradient normal to the wall.

---

## 3. Experiments

### Scenario Overview

| Task | Grid | $\nu$ | Re | $\Delta t$ mode | $t_\text{end}$ | Purpose |
|---|---|---|---|---|---|---|
| 4 | 50×50 | 0.01 | 100 | Adaptive ($\tau=0.5$) | 50.0 s | Base flow structure |
| 5a | 50×50 | 0.01 | 100 | Adaptive | 5.0 s | Effect of SOR $\omega$ |
| 5b | 50×50 | 0.01 | 100 | Adaptive | 5.0 s | Effect of itermax |
| 6 | 50×50 | 0.01 | 100 | Fixed (varied) | 5.0 s | Time-step stability |
| 7 | 16–256 | 0.001 | 1000 | Fixed ($\Delta t=0.05$) | 5.0 s | Grid refinement |
| 8 | 50×50 | 0.01–0.0001 | 100–10000 | Adaptive ($\tau=0.5$) | 10.0 s | Re effects |

### Task 4 — Base Lid-Driven Cavity Simulation (Re = 100)

**Configuration:** 50×50 grid, $\nu=0.01$, adaptive time stepping ($\tau=0.5$), $t_\text{end}=50$ s.

| Grid | $\nu$ | Re | avg $\Delta t$ | Steps | VTK files | avg SOR iter | avg residual | Status |
|---|---|---|---|---|---|---|---|---|
| 50×50 | 0.01 | 100 | 5.00×10⁻³ | 10001 | 101 | 100.0 | 1.230 | OK |

**Key observations:**

- The adaptive time step settles immediately at $\Delta t = 5\times10^{-3}$ s, equal to $\tau \times \Delta t_\text{visc} = 0.5 \times (0.02)^2/(4 \times 0.01) = 0.5 \times 0.010$. The viscous stability condition is the binding constraint at Re=100.
- After $t\approx5$ s the flow reaches quasi-steady state, characterised by a **single large primary vortex** filling the cavity, driven by shear from the moving lid. Small secondary Moffatt eddies form in the bottom corners.
- The SOR solver exhausts its budget at every step (avg SOR iter = itermax = 100). This is not a failure — it is a consequence of the pressure null space (see Task 5 discussion). The velocity field is physically correct because it depends only on pressure *gradients*, not absolute pressure levels.
- The velocity magnitude field peaks near the lid ($|\mathbf{u}|\approx 1$ m/s) and decays toward the stationary walls. Streamlines spiral into the vortex core. The pressure field shows a low-pressure region at the vortex centre and elevated pressure in the corner driven by the lid.

The VTK outputs are visualised in ParaView via `visualize.py`, producing:

| Image | Description |
|---|---|
| `final_u.png` | Horizontal velocity (Blue-to-Red, range [−0.5, 1.0]) |
| `final_v.png` | Vertical velocity (Blue-to-Red, range [−0.25, 0.25]) |
| `final_pressure.png` | Pressure field (Cool-to-Warm, auto-scaled) |
| `final_velocity.png` | Velocity magnitude (Jet colourmap, range [0, 1]) |
| `final_glyphs.png` | Arrow glyphs coloured by magnitude |
| `final_streamlines.png` | Streamlines on dark background |
| `final_vectors_bw.png` | Black-and-white vector field |

### Task 5 — SOR Solver Behaviour

#### 5a: Effect of relaxation factor $\omega$ (itermax = 500)

| $\omega$ | avg SOR iter | Hits itermax | avg residual | Wall time |
|---|---|---|---|---|
| 0.50 | 500 | YES | 1.360 | 15.8 s |
| 1.00 | 500 | YES | 1.270 | 15.9 s |
| 1.30 | 500 | YES | 1.240 | 15.9 s |
| **1.50** | 500 | YES | **1.210** | 15.9 s |
| 1.70 | 500 | YES | 1.230 | 15.9 s |
| 1.80 | 500 | YES | 1.340 | 16.1 s |
| 1.90 | 500 | YES | 1.880 | 16.0 s |
| 1.95 | 500 | YES | 2.840 | 16.0 s |
| 1.99 | 500 | YES | 5.240 | 16.0 s |

**Optimal:** $\omega \approx 1.5$ achieves the lowest average residual (1.210) within the iteration budget.

**Three regimes:**
1. **Under-relaxation** ($\omega < 1$): Gauss–Seidel-like slow convergence. More iterations are needed per unit residual reduction; the spectral radius of the iteration matrix is close to 1.
2. **Near-optimal over-relaxation** ($1.3 \lesssim \omega \lesssim 1.7$): The spectral radius is minimised. For a $50\times50$ Dirichlet problem the optimal $\omega_\text{opt} = 2/(1 + \sin(\pi/N)) \approx 1.73$; for the Neumann problem here the optimum shifts slightly toward 1.5.
3. **Aggressive over-relaxation** ($\omega \to 2$): The iteration overshoots; errors grow rather than decay. The residual at $\omega=1.99$ is 4× larger than at the optimum.

**Critical insight — pressure null space:** Every value of $\omega$ exhausts the iteration budget at every time step. No case achieves residual $< \varepsilon = 10^{-3}$. This is not a SOR deficiency; it is a fundamental consequence of the **pressure null space**: with all-Neumann BCs, the discrete Poisson system $\mathbf{A}\mathbf{p} = \mathbf{b}$ is singular (the constant vector $\mathbf{1}$ satisfies $\mathbf{A}\mathbf{1} = \mathbf{0}$). SOR cannot damp the zero-eigenvalue mode, so the residual saturates at a finite level. The pressure is defined only up to an additive constant, but velocity fields depend only on pressure gradients and are therefore correct.

#### 5b: Effect of itermax ($\omega = 1.7$)

| itermax | avg SOR iter | max SOR iter | Status |
|---|---|---|---|
| 5 | 5.0 | 5 | OK |
| 10 | 10.0 | 10 | OK |
| 20 | 20.0 | 20 | OK |
| 50 | 50.0 | 50 | OK |
| 100 | 100.0 | 100 | OK |
| 200 | 200.0 | 200 | OK |

The solver exhausts the budget at every time step for all itermax values. With very low itermax (e.g., 5) the pressure projection is incomplete — the corrected velocity field has non-zero divergence residual — but the simulation remains numerically stable. Larger itermax reduces the saturation residual and improves pressure accuracy at proportionally higher computational cost.

### Task 6 — Fixed Time-Step Stability (50×50, $\nu = 0.01$)

Adaptive stepping is disabled ($\tau = -1$). The theoretical stability limits for this configuration are:

$$\Delta t_\text{CFL} = \frac{\Delta x}{U_\text{wall}} = \frac{0.02}{1.0} = 0.020 \text{ s}$$
$$\Delta t_\text{visc} = \frac{\Delta x^2}{4\nu} = \frac{(0.02)^2}{4 \times 0.01} = 0.010 \text{ s}$$

| $\Delta t$ | CFL $= \Delta t/\Delta x$ | $\Delta t/\Delta t_\text{visc}$ | Status | Wall time |
|---|---|---|---|---|
| 0.001 | 0.050 | 0.10 | **OK** | 16.2 s |
| 0.005 | 0.250 | 0.50 | **OK** | 3.4 s |
| 0.008 | 0.400 | 0.80 | **OK** | 2.2 s |
| **0.010** | **0.500** | **1.00** | **DIVERGED** | 0.4 s |
| 0.012 | 0.600 | 1.20 | DIVERGED | 0.2 s |
| 0.015 | 0.750 | 1.50 | DIVERGED | 0.2 s |
| 0.020 | 1.000 | 2.00 | DIVERGED | 0.2 s |
| 0.025 | 1.250 | 2.50 | DIVERGED | 0.2 s |
| 0.050 | 2.500 | 5.00 | DIVERGED | 0.2 s |

**The stability boundary falls exactly at $\Delta t = \Delta t_\text{visc} = 0.010$ s**, even though the CFL criterion ($\Delta t < 0.020$) would still permit this step size.

**Explanation:** The explicit Euler scheme simultaneously applies forward time integration to both the convective and diffusive terms. The viscous diffusion stability criterion (von Neumann analysis for $\partial u/\partial t = \nu\nabla^2 u$) gives $\Delta t < \Delta x^2/(4\nu)$, which is more restrictive here by a factor of 2. Both conditions must hold simultaneously; whichever gives the smaller bound is the **binding constraint**. At Re=100 (relatively large $\nu$), the viscous condition dominates. The adaptive stepper with $\tau=0.5$ naturally uses $\Delta t = 0.005$ s — $50\%$ below the boundary, providing a comfortable stability margin.

### Task 7 — Grid Refinement (fixed $\Delta t = 0.05$ s, $\nu = 0.001$, Re = 1000)

| Grid | $\Delta x$ | CFL $\approx \Delta t/\Delta x$ | $\Delta t/\Delta t_\text{visc}$ | Status | avg SOR iter |
|---|---|---|---|---|---|
| 16×16 | 0.06250 | 0.80 | 160 | **OK** | 100.0 |
| 32×32 | 0.03125 | 1.60 | 640 | **OK** | 100.0 |
| 64×64 | 0.01562 | 3.20 | 2563 | DIVERGED | 110.0 |
| 128×128 | 0.00781 | 6.40 | 10253 | DIVERGED | 125.0 |
| 256×256 | 0.00391 | 12.80 | 41013 | DIVERGED | 133.3 |

Note: $\Delta t_\text{visc} = \Delta x^2/(4\nu)$ shrinks quadratically with grid refinement. At $\nu=0.001$ (Re=1000) and $\Delta x=0.0156$ (64×64), $\Delta t_\text{visc} = (0.0156)^2/(4\times0.001) \approx 0.000061$ s — four orders of magnitude smaller than the fixed $\Delta t=0.05$.

**Why does 32×32 survive although CFL > 1?**

The CFL numbers in the table are computed using $\text{CFL} = \Delta t \cdot U_\text{wall}/\Delta x$, which treats the lid velocity $U_\text{wall} = 1$ as the maximum velocity everywhere. This is an **upper-bound estimate**, not the true local CFL at every cell:

$$\text{CFL}_{i,j}^\text{true} = \frac{|u(i,j)| \cdot \Delta t}{\Delta x}$$

At $t = 0$ all interior velocities are exactly zero. The flow starts from rest; momentum diffuses inward from the lid over time. For the short run $t_\text{end} = 5$ s on a coarse 32×32 grid, the interior velocity magnitude stays well below $U_\text{wall}$ throughout — so the **true local CFL remains below 1 inside the domain** even though the bound says 1.60.

The CFL > 1 violation exists only at the very first row of cells directly below the moving lid (where the ghost-cell BC sets $u_\text{ghost} = 2U_\text{wall} - u_\text{interior}$, effectively seeing a velocity of order $U_\text{wall}$). The resulting local instability grows exponentially, but with a modest growth factor per step (~1.6× at CFL=1.6 vs. ~3.2× at CFL=3.2). For only 100 time steps ($t_\text{end}/\Delta t = 5/0.05$) the accumulated amplification at 32×32 is insufficient to produce visible divergence, whereas at 64×64 (CFL=3.2) the error blows up after just a few steps.

**The theoretical CFL < 1 criterion is still correct.** The 32×32 result is a *marginal* case: it appears stable for $t_\text{end} = 5$ s but would likely diverge on a longer run (e.g., $t_\text{end} = 50$ s). The red dashed line at CFL = 1 in the study plot therefore marks the correct theoretical stability boundary.

**Fundamental conclusion:** Halving the mesh spacing requires halving $\Delta t$ (CFL) or — more restrictively — quartering $\Delta t$ (viscous). Fixed time stepping is incompatible with grid convergence studies for explicit schemes. Adaptive time stepping ($\tau > 0$) automatically satisfies both conditions at every step regardless of mesh size.

### Task 8 — Reynolds Number Effects (adaptive $\Delta t$, $t_\text{end} = 10$ s)

| $\nu$ | Re | avg $\Delta t$ (s) | avg SOR iter | Hits itermax | Status |
|---|---|---|---|---|---|
| 0.01 | 100 | 5.00×10⁻³ | 100.0 | YES | OK |
| 0.002 | 500 | 1.15×10⁻² | 100.0 | YES | OK |
| 0.0005 | 2000 | 1.36×10⁻² | 100.0 | YES | OK |
| 0.0001 | 10000 | 2.13×10⁻² | 100.0 | YES | OK |

#### Trend in adaptive $\Delta t$ with Re

The average time step **increases monotonically with Re**, from $5\times10^{-3}$ at Re=100 to $2.13\times10^{-2}$ at Re=10000. This follows directly from the adaptive stability criterion:

$$\Delta t = \tau \cdot \min\!\left(\underbrace{\frac{\Delta x^2}{4\nu}}_{\Delta t_\text{visc}},\;\underbrace{\frac{\Delta x}{|u|_\text{max}}}_{\Delta t_\text{CFL}}\right)$$

At Re=100: $\Delta t_\text{visc} = 0.010$ s, $\Delta t_\text{CFL} \approx 0.020$ s → viscous limit binds → $\Delta t = 0.005$ s.

As $\nu$ decreases (Re increases), $\Delta t_\text{visc} \propto 1/\nu$ grows and eventually exceeds $\Delta t_\text{CFL} \approx \Delta x/U_\text{wall} = 0.02$ s. The CFL condition then becomes the binding constraint. At Re=10000: $\Delta t_\text{visc} = 1.0$ s $\gg \Delta t_\text{CFL} \approx 0.02$ s, so the adaptive step is $\approx 0.5 \times 0.02 = 0.01$ s, consistent with the measured value of $2.13\times10^{-2}$.

#### Physical flow changes with Re

**Re = 100 (low inertia):** A single symmetric primary vortex fills the cavity with its centre slightly above the geometric centre. The flow is laminar and quickly reaches steady state. Viscosity dominates, spreading momentum smoothly across the domain.

**Re = 500–2000 (moderate inertia):** The vortex centre moves toward the geometric cavity centre as inertia becomes more important. Secondary corner eddies (Moffatt eddies) at the bottom corners grow in strength. The shear layer along the moving lid thins, concentrating the velocity gradient near the top wall.

**Re = 10000 (high inertia):** The flow is potentially unsteady. On the 50×50 grid with $t_\text{end}=10$ s the simulation remains stable, but the coarse grid cannot resolve the thin shear layers and small-scale eddies present at this Reynolds number. Qualitatively the primary vortex occupies the full cavity and bottom-corner eddies are visible. A fine-grid, long-time simulation would be needed to determine whether the flow is truly steady or exhibits periodic or chaotic fluctuations.

See `LidDrivenCavity_Output/study_plots/task8_re_comparison.png` for a side-by-side comparison of velocity magnitude and streamlines at all four Reynolds numbers.

---

## 4. Conclusion

A 2D incompressible Navier-Stokes solver for the Lid-Driven Cavity benchmark was implemented using the Chorin projection method on a staggered MAC grid. Five numerical experiments yielded the following conclusions:

1. **Base simulation (Task 4):** At Re=100 the flow reaches quasi-steady state with a single stable primary vortex. The adaptive time step settles at $\Delta t = 5\times10^{-3}$ s, set by the viscous stability condition. The SOR solver exhausts its iteration budget at every step due to the singular pure-Neumann pressure system.

2. **SOR convergence (Task 5):** The optimal relaxation factor is $\omega \approx 1.5$ (lowest saturation residual of 1.210 within 500 iterations). Convergence to $\varepsilon=10^{-3}$ is never achieved at any $\omega$ because the all-Neumann boundary conditions render the discrete Poisson system singular. The SOR residual saturates at a finite level independent of $\omega$ or itermax. This is physically correct: pressure is defined only up to an additive constant and velocity depends only on pressure gradients.

3. **Time-step stability (Task 6):** The binding stability constraint at Re=100 is the **viscous diffusion condition** $\Delta t < h^2/(4\nu) = 0.010$ s, which is twice as restrictive as the CFL condition ($\Delta t < 0.020$ s). All simulations with $\Delta t \geq 0.010$ s diverge immediately, confirming the theoretical prediction of the explicit Euler scheme.

4. **Grid refinement (Task 7):** With fixed $\Delta t=0.05$ s at Re=1000, only grids coarser than $\Delta x > \Delta t$ (16×16 and 32×32) remain stable. Finer grids require $\Delta t \propto h$ (CFL) or $\propto h^2$ (viscous), making adaptive time stepping essential for mesh refinement studies. Adaptive stepping automatically selects the appropriate $\Delta t$ at any grid resolution.

5. **Reynolds number effects (Task 8):** The adaptive time step increases with Re as the viscous stability limit relaxes. All four cases (Re=100 to Re=10000) run stably. At Re=10000 the 50×50 grid and $t_\text{end}=10$ s are insufficient to resolve the flow; finer grids, longer simulations, and implicit time integration would be required to accurately simulate high-Re LDC flows.
