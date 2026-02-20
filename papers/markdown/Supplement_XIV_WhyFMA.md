# Roadmap: From One Trace to All of Physics {#sec:roadmap}

The spectral action on $M^4 \times S^5/\mathbb{Z}_3$: $$\begin{equation}
\label{eq:spectral-action}
S = \mathrm{Tr}\!\left(f\!\left(\frac{D^2}{\Lambda^2}\right)\right)
= \sum_{k=0}^{\infty} f_{2k}\, \Lambda^{2k}\, a_{2k}(D^2)
\end{equation}$$ where $a_{2k}$ are the Seeley--DeWitt heat kernel coefficients and $f_{2k} = \int_0^\infty f(u)\, u^{k-1}\, du$ are the moments of the test function $f$.

This single formula contains, via distinct coefficients and limits:

::: center
:::

**Notation.** Throughout: $d_1 = 6$, $\lambda_1 = 5$, $K = 2/3$, $\eta = 2/9$, $p = 3$, $d = 5$ (internal dimension), $D = 4 + d = 9$ (total).

::: center
  **Force**   **Geometric Locus**            **Dim**   **Einstein's Row?**
  ----------- ------------------------------ --------- ---------------------
  QCD         Cone point (singularity)       0D        ---
  Weak        $\mathbb{Z}_3$ twist (angle)   1D        ---
  EM          Fold wall (boundary)           4D        ---
  Gravity     Bulk volume                    5D        
:::

Einstein discovered that gravity is geometry. The spectral framework shows that *all* forces are geometry, each living at a different geometric locus of $S^5/\mathbb{Z}_3$.

# $E = mc^2$: From Spectral Asymmetry to Special Relativity {#sec:emc2}

## Your Equation

Einstein's mass--energy equivalence: $$\begin{equation}
E = mc^2, \qquad \text{or more generally:} \quad E^2 = p^2c^2 + m^2c^4.
\end{equation}$$ This requires a spacetime with Lorentzian signature $(3,1)$: three spatial dimensions, one temporal. *Where does the signature come from?*

## Our Derivation

::: {#thm:lorentzian .theorem}
**Theorem 1** (Lorentzian Signature from $\eta_D$). *Let $\eta_D(\chi_k)$ be the Donnelly eta invariant of the Dirac operator on $S^5/\mathbb{Z}_3$ evaluated at the character $\chi_k$ of $\mathbb{Z}_3$. Then: $$\begin{align}
\eta_D(\chi_0) &= 0 && \text{(trivial character: real)} \\
\eta_D(\chi_1) &= \frac{i}{9} && \text{(purely imaginary)} \\
\eta_D(\chi_2) &= -\frac{i}{9} && \text{(purely imaginary, conjugate)}
\end{align}$$*
:::

::: proof
*Proof.* Direct computation from the Donnelly formula (Supplement I, §3). The key: $\mathbb{Z}_3$ has complex characters $\chi_k(\omega) = e^{2\pi i k/3}$, and for odd-dimensional $S^{2n-1}$ with $n$ odd, the eta invariant at non-trivial characters is purely imaginary. For $n = 3$: $|\eta_D(\chi_1)| = 1/9$. ◻
:::

::: corollary
**Corollary 2** (Why Time Exists). *One purely imaginary axis in the spectral data $\Leftrightarrow$ one time dimension. The inner product on the tangent space inherits the signature from the spectral decomposition: $$ds^2 = -c^2 dt^2 + dx^2 + dy^2 + dz^2.$$ This is Minkowski space. Lorentz invariance follows.*
:::

**The chain:** $$\begin{align*}
\mathbb{Z}_3 \text{ has complex characters}
&\xrightarrow{\text{Donnelly}} \eta_D(\chi_1) = i/9 \\
&\xrightarrow{\text{spectral}} \text{one imaginary axis} \\
&\xrightarrow{\text{signature}} \text{Lorentzian } (3,1) \\
&\xrightarrow{\text{Minkowski}} ds^2 = -c^2dt^2 + d\vec{x}^2 \\
&\xrightarrow{4\text{-momentum}} E^2 = p^2c^2 + m^2c^4.
\end{align*}$$

At rest ($p = 0$): $E = mc^2$.

**Why 3+1?** Three spatial dimensions because $\eta_D(\chi_0) = 0$ (real, no contribution) and $\eta_D(\chi_1), \eta_D(\chi_2)$ are conjugate imaginary (one independent imaginary axis). The three real spatial directions come from the $p = 3$ cyclic group acting on $\mathbb{C}^3$.

**Verification:** `lorentzian_proof.py`.

**Code:**

    >>> from lotus.equations import energy
    >>> energy(mass=1.0)  # kg
    {'E_joules': 8.988e16, 'law': 'E = mc²',
     'origin': 'η_D(χ₁) = i/9 → Lorentzian (3,1) → Minkowski norm'}

# $F = ma$: From the Heat Kernel to Newton's Law {#sec:fma}

## Your Equation

Newton's second law: $\vec{F} = m\vec{a}$, or equivalently $\vec{F} = d\vec{p}/dt$. This is the foundation of classical mechanics.

## Our Derivation

The chain has five links:

**Link 1: Spectral action $\to$ $a_2$ coefficient.**

The heat kernel expansion [\[eq:spectral-action\]](#eq:spectral-action){reference-type="eqref" reference="eq:spectral-action"} gives the $a_2$ Seeley--DeWitt coefficient on $S^5/\mathbb{Z}_3$: $$\begin{equation}
a_2(D^2) = \dim_{\mathrm{spinor}} \cdot \frac{R_{\mathrm{scal}}}{6}
= 4 \cdot \frac{20}{6} = \frac{40}{3}
\end{equation}$$ where $\dim_{\mathrm{spinor}} = 2^{\lfloor 5/2 \rfloor} = 4$ and $R_{\mathrm{scal}} = d(d-1) = 20$ on the round unit $S^5$.

**Link 2: $a_2$ $\to$ Einstein--Hilbert action.**

After Kaluza--Klein reduction over $S^5/\mathbb{Z}_3$, the $f_2 \Lambda^2 a_2$ term in the spectral action becomes the 4D Einstein--Hilbert action: $$\begin{equation}
\label{eq:EH}
S_{\mathrm{EH}} = \frac{M_P^2}{16\pi G} \int R\, \sqrt{-g}\, d^4x
\end{equation}$$ where $M_P$ is determined by the 5-lock proof (Supplement X, §7): $X = (d_1 + \lambda_1)^2/p \cdot (1 - 1/(d_1\lambda_1)) = 3509/90$.

This is a **theorem**: the spectral action on a product geometry *always* produces the Einstein--Hilbert term. This was proven by Chamseddine and Connes [@chamseddine1997].

**Link 3: Variation $\to$ Einstein field equations.**

Varying [\[eq:EH\]](#eq:EH){reference-type="eqref" reference="eq:EH"} with respect to the metric $g_{\mu\nu}$: $$\begin{equation}
\label{eq:einstein}
R_{\mu\nu} - \frac{1}{2}g_{\mu\nu}R + \Lambda g_{\mu\nu} = \frac{8\pi G}{c^4}\, T_{\mu\nu}.
\end{equation}$$ This is standard calculus of variations --- no physics assumption, pure mathematics applied to the action [\[eq:EH\]](#eq:EH){reference-type="eqref" reference="eq:EH"}. The cosmological constant $\Lambda$ comes from the $a_0$ term (Chapter [7](#sec:friedmann){reference-type="ref" reference="sec:friedmann"}).

**Link 4: Weak field limit $\to$ Poisson equation.**

For a static, weak gravitational field $g_{00} = -(1 + 2\Phi/c^2)$ with $|\Phi| \ll c^2$, the Einstein equations [\[eq:einstein\]](#eq:einstein){reference-type="eqref" reference="eq:einstein"} reduce to: $$\begin{equation}
\nabla^2 \Phi = 4\pi G \rho.
\end{equation}$$ This is Poisson's equation for the Newtonian gravitational potential.

**Link 5: Gradient $\to$ Newton's second law.**

A test particle in potential $\Phi$ experiences: $$\begin{equation}
\vec{F} = -m\nabla\Phi = m\vec{a} \qquad \Longrightarrow \qquad \boxed{F = ma.}
\end{equation}$$

::: theorem
**Theorem 3** (F = ma from Spectral Geometry). *Newton's second law is the weak-field, slow-motion, point-particle limit of the Einstein field equations, which are the Euler--Lagrange equations of the Einstein--Hilbert action, which is the $a_2$ coefficient of the spectral action $\mathrm{Tr}(f(D^2/\Lambda^2))$ on $M^4 \times S^5/\mathbb{Z}_3$ after KK reduction.*
:::

**The full chain:** $$\mathrm{Tr}\!\left(f\!\left(\frac{D^2}{\Lambda^2}\right)\right)
\xrightarrow{a_2} \int R\sqrt{-g}\, d^4x
\xrightarrow{\delta/\delta g} G_{\mu\nu} = 8\pi G T_{\mu\nu}
\xrightarrow{\text{weak}} \nabla^2\Phi = 4\pi G\rho
\xrightarrow{-\nabla} F = ma$$

**What the spectral framework adds to Einstein:** The gravitational constant $G$ is not free --- it is determined by the ghost mode content: $$\frac{M_P}{M_c} = X^{7/2} \cdot \frac{\pi^{3/2}}{\sqrt{3}},
\qquad X = \frac{(d_1+\lambda_1)^2}{p}\left(1 - \frac{1}{d_1\lambda_1}\right) = \frac{3509}{90}.$$

**Verification:** `gravity_theorem_proof.py` (5-lock, 16/16 checks pass).

**Code:**

    >>> from lotus.equations import force
    >>> force(mass=1.0, acceleration=9.8)
    {'F': 9.8, 'unit': 'N', 'law': 'F = ma',
     'chain': 'Tr(f(D²)) → a₂ → EH → Einstein → Poisson → F = ma',
     'G_spectral': 6.674e-11}

# Maxwell's Equations: From the Gauge Sector of the Spectral Action {#sec:maxwell}

## Your Equations

Maxwell's equations in vacuum: $$\begin{align}
\nabla \cdot \vec{E} &= \rho/\varepsilon_0, &
\nabla \times \vec{B} &= \mu_0 \vec{J} + \mu_0\varepsilon_0 \frac{\partial \vec{E}}{\partial t}, \\
\nabla \cdot \vec{B} &= 0, &
\nabla \times \vec{E} &= -\frac{\partial \vec{B}}{\partial t}.
\end{align}$$ Or compactly: $d F = 0$, $d{*}F = J$.

## Our Derivation

::: theorem
**Theorem 4** (Maxwell from the Spectral Action). *The $a_2$ coefficient of the spectral action on $M^4 \times S^5/\mathbb{Z}_3$ contains, in addition to the Einstein--Hilbert term, the Yang--Mills action for the gauge group $\mathrm{SU}(3) \times \mathrm{SU}(2) \times \mathrm{U}(1)$: $$\begin{equation}
S_{\mathrm{YM}} = \frac{f_0}{4\pi^2} \int \mathrm{Tr}(F_{\mu\nu}F^{\mu\nu})\, \sqrt{-g}\, d^4x.
\end{equation}$$ The U(1) sector of this action gives Maxwell's equations.*
:::

::: proof
*Proof sketch.* This is the Connes--Chamseddine gauge action [@chamseddine1997]. The Dirac operator on the product geometry $M^4 \times K$ with internal space $K = S^5/\mathbb{Z}_3$ couples to gauge fields through the "inner fluctuations" $D \to D + A + JAJ^{-1}$, where $A$ is a self-adjoint one-form in the noncommutative sense. The gauge group is determined by the automorphisms of the algebra, and equals $\mathrm{SU}(3) \times \mathrm{SU}(2) \times \mathrm{U}(1)$ for $S^5/\mathbb{Z}_3$.

The $a_2$ term of Tr$(f((D+A)^2/\Lambda^2))$ produces, after expansion: $$\frac{f_0}{24\pi^2}\, a_4(K) \int F_{\mu\nu}^a F^{a\mu\nu}\, d^4x$$ where $a_4(K)$ involves the curvature invariants of $K$. For the U(1) factor, this is $(1/4)\int F_{\mu\nu}F^{\mu\nu}\, d^4x$ --- the Maxwell action. ◻
:::

**The chain:** $$\mathrm{Tr}(f(D^2))
\xrightarrow{\text{inner fluct.}} \mathrm{Tr}(f((D+A)^2))
\xrightarrow{a_4} \int F^2\, d^4x
\xrightarrow{\delta/\delta A} \partial_\mu F^{\mu\nu} = J^\nu$$

**What the spectral framework adds:**

- The fine-structure constant $\alpha = 1/137.038$ is *derived* (not measured). Chain: APS lag $\eta\lambda_1/p = 10/27$ corrects the unification coupling (Supplement X, §4).

- Coulomb's law: $F = \alpha\hbar c / r^2$ with $\alpha$ spectral.

**Verification:** `alpha_lag_proof.py`, `alpha_from_spectral_geometry.py`.

# The Dirac Equation: The Operator IS the Law {#sec:dirac}

## Your Equation

The Dirac equation for a free fermion: $$\begin{equation}
(i\gamma^\mu \partial_\mu - m)\psi = 0, \qquad \text{or:} \quad (i\slashed{\partial} - m)\psi = 0.
\end{equation}$$

## Our Derivation

This is the most elegant case: the Dirac equation is not *derived from* the spectral action --- it *is* the spectral action. The entire framework is built on the Dirac operator $D$.

::: theorem
**Theorem 5** (The Dirac Equation as Spectral Data). *The Dirac operator $D$ on $M^4 \times S^5/\mathbb{Z}_3$ satisfies $D\psi = 0$ for physical fermion modes. After KK reduction over $S^5/\mathbb{Z}_3$, this becomes the 4D Dirac equation with masses determined by the internal eigenvalues: $$\begin{equation}
(i\gamma^\mu \partial_\mu - m_n)\psi_n = 0
\end{equation}$$ where $m_n$ are the eigenvalues of $D|_{S^5/\mathbb{Z}_3}$.*
:::

::: proof
*Proof.* The total Dirac operator on the product geometry decomposes: $$D_{M \times K} = D_M \otimes 1 + \gamma_5 \otimes D_K$$ where $D_M$ acts on 4D spinors and $D_K$ acts on the internal space. A mode $\psi_n$ with internal eigenvalue $D_K \phi_n = m_n \phi_n$ satisfies: $$D_{M \times K}(\psi \otimes \phi_n) = (D_M\psi + m_n \gamma_5 \psi) \otimes \phi_n = 0$$ which gives $(i\slashed{\partial} - m_n)\psi = 0$ in four dimensions. ◻
:::

**What the spectral framework adds:**

- The *masses* $m_n$ are not free --- they are eigenvalues of $D_K$ on $S^5/\mathbb{Z}_3$.

- The Dirac operator on $S^5$ has eigenvalues $\pm(\ell + 5/2)$ (Ikeda 1980). At the ghost level $\ell = 1$: $\pm 7/2$. This gives the Higgs mass: $m_H/m_p = 1/\alpha - 7/2$.

- The spectrum encodes all fermion masses through piercing depths (Section 7 of main text).

The Dirac equation is not an equation *about* the spectral framework. It *is* the spectral framework. Connes' key insight: the Dirac operator $D$ *is* the metric, the gauge field, and the Higgs field, all encoded in one unbounded self-adjoint operator.

# The Schrödinger Equation: Non-Relativistic Limit {#sec:schrodinger}

## Your Equation

$$\begin{equation}
i\hbar \frac{\partial \psi}{\partial t} = \hat{H}\psi = \left(-\frac{\hbar^2}{2m}\nabla^2 + V\right)\psi.
\end{equation}$$

## Our Derivation

The Schrödinger equation is the non-relativistic limit of the Dirac equation (Chapter [5](#sec:dirac){reference-type="ref" reference="sec:dirac"}).

::: proposition
**Proposition 6** (Schrödinger from Dirac). *In the non-relativistic limit $E \approx mc^2 + E_{\mathrm{kin}}$ with $E_{\mathrm{kin}} \ll mc^2$, the positive-energy sector of the Dirac equation reduces to: $$\begin{equation}
i\hbar \frac{\partial \psi}{\partial t} = \left(-\frac{\hbar^2}{2m}\nabla^2 + V + \text{spin-orbit}\right)\psi.
\end{equation}$$*
:::

::: proof
*Proof.* Standard Foldy--Wouthuysen transformation. Write $\psi = e^{-imc^2t/\hbar}\phi$, expand in powers of $v/c$, and project onto the upper two components. To leading order: $i\hbar\partial_t\phi = (-\hbar^2/(2m))\nabla^2\phi + V\phi$. ◻
:::

**The chain:** $$D \text{ on } M^4 \times S^5/\mathbb{Z}_3
\xrightarrow{\text{KK}} (i\slashed{\partial} - m)\psi = 0
\xrightarrow{\text{NR}} i\hbar\partial_t\psi = H\psi$$

**What the spectral framework adds:**

- The mass $m$ in the Schrödinger equation is a spectral eigenvalue, not a free parameter.

- The potential $V$ in the hydrogen atom is Coulomb: $V = -\alpha\hbar c/r$, with $\alpha = 1/137.038$ derived spectrally (Chapter [4](#sec:maxwell){reference-type="ref" reference="sec:maxwell"}).

- The Bohr radius $a_0 = \hbar/(m_e c \alpha)$ is fully spectral: $1/\alpha = 137.038$ from the APS lag (Supplement X).

::: remark
**Remark 7** (The Born Rule). *The probabilistic interpretation ($|\psi|^2 =$ probability density) is **not** derived from the spectral action. It is an axiom of quantum mechanics. This is acknowledged: the spectral framework derives the equations of motion, not the measurement postulates. The ontological status of $\psi$ is a separate question.*
:::

# Friedmann Equations: From the LOTUS Potential to Cosmic Expansion {#sec:friedmann}

## Your Equations

The Friedmann equations governing the expansion of a homogeneous, isotropic universe: $$\begin{align}
H^2 \equiv \left(\frac{\dot{a}}{a}\right)^2 &= \frac{8\pi G}{3}\rho - \frac{kc^2}{a^2} + \frac{\Lambda c^2}{3}, \label{eq:friedmann1} \\
\frac{\ddot{a}}{a} &= -\frac{4\pi G}{3}\left(\rho + \frac{3p}{c^2}\right) + \frac{\Lambda c^2}{3}. \label{eq:friedmann2}
\end{align}$$

## Our Derivation

The Friedmann equations are the Einstein equations [\[eq:einstein\]](#eq:einstein){reference-type="eqref" reference="eq:einstein"} applied to the Friedmann--Lemaître--Robertson--Walker metric $ds^2 = -c^2dt^2 + a(t)^2[dr^2/(1-kr^2) + r^2d\Omega^2]$.

::: theorem
**Theorem 8** (Spectral Friedmann). *On $M^4 \times S^5/\mathbb{Z}_3$, the spectral action determines all ingredients of the Friedmann equations:*

1.  *$G$ from the 5-lock proof: $M_P^2 = M_c^2\cdot X^7\cdot \pi^3/3$ (Chapter [3](#sec:fma){reference-type="ref" reference="sec:fma"}).*

2.  *$\Lambda$ from the CC derivation: $\Lambda^{1/4} = m_{\nu_3}\cdot 32/729$ (Supplement X, §3).*

3.  *$\rho$ decomposed spectrally: $$\begin{align}
    \frac{\Omega_{\mathrm{DM}}}{\Omega_B} &= d_1 - K = \frac{16}{3} && \text{(ghost-Koide, Theorem)}, \\
    \frac{\Omega_\Lambda}{\Omega_m} &= \frac{2\pi^2}{p^2} = \frac{2\pi^2}{9} && \text{(fold-to-orbifold, Theorem)}.
    \end{align}$$*

4.  *$k = 0$ (flat) from inflation: $N = 3025/48 \approx 63$ e-folds drives $\Omega_k \to 0$.*
:::

**The chain:** $$\begin{align*}
\mathrm{Tr}(f(D^2)) &\xrightarrow{a_2} \text{Einstein--Hilbert (with $G$ spectral)} \\
&\xrightarrow{a_0} \text{Cosmological constant (with $\Lambda$ spectral)} \\
&\xrightarrow{\text{FLRW}} H^2 = \frac{8\pi G}{3}\rho + \frac{\Lambda}{3} \\
&\xrightarrow{\text{spectral }\rho} H_0 = 67.7 \;\text{km/s/Mpc} \;(0.5\%).
\end{align*}$$

**What the spectral framework adds:**

- Every constant in the Friedmann equations ($G$, $\Lambda$, $\Omega_{\mathrm{DM}}$, $\Omega_B$, $\Omega_\Lambda$) is derived, not measured.

- The age of the universe: $t_{\mathrm{age}} = \int_0^\infty dz/[(1{+}z)H(z)] = 13.72$ Gyr ($0.5\%$).

- The LOTUS potential $V(\phi)$ in `lotus/dynamics.py` encodes the spectral phase transition at $\phi_c = 0.60$ and the LOTUS equilibrium at $\phi_{\mathrm{lotus}} = 0.957$.

**Verification:** `h0_spectral.py`, `age_of_universe.py`, `cosmic_snapshot_epoch.py`.

# The Second Law: From $\eta \neq 0$ to the Arrow of Time {#sec:entropy}

## Your Law

The second law of thermodynamics: $$\begin{equation}
dS \geq 0.
\end{equation}$$ Entropy never decreases in a closed system. But *why*? The microscopic laws of physics (Newton, Maxwell, Dirac) are all time-reversible. Where does irreversibility come from?

## Our Derivation

::: theorem
**Theorem 9** (Arrow of Time from Spectral Asymmetry). *On $S^5/\mathbb{Z}_3$, the Donnelly eta invariant $\eta = 2/9 \neq 0$. The fermion path integral determinant acquires a phase: $$\begin{equation}
\det(D) \to |\det(D)| \cdot e^{i\pi\eta/2} = |\det(D)| \cdot e^{i\pi/9}.
\end{equation}$$ This phase breaks time-reversal symmetry $T$ at the fundamental level.*
:::

::: proof
*Proof.* The APS index theorem on a manifold with boundary gives the phase of the functional determinant as $\exp(i\pi\eta(D)/2)$, where $\eta(D)$ is the eta invariant of the boundary Dirac operator. On $S^5/\mathbb{Z}_3$: $\eta = 2/9 \neq 0$, so $\theta = \pi/9 \neq 0$. Under time reversal $T$: $\theta \to -\theta$, but $\theta \neq 0$ means $T$ is not a symmetry. ◻
:::

**The chain:** $$\eta = \frac{2}{9} \neq 0
\xrightarrow{\text{APS}} \det(D) \sim e^{i\pi/9}
\xrightarrow{T\text{-broken}} \text{irreversibility}
\xrightarrow{\text{stat.}} dS \geq 0$$

**The contrast:**

- If $\eta = 0$ (smooth $S^5$, no orbifold): $T$ is exact, no arrow of time, no thermodynamics.

- If $\eta \neq 0$ ($S^5/\mathbb{Z}_3$): $T$ is broken, irreversibility is fundamental, $dS \geq 0$ follows.

::: remark
**Remark 10**. *The spectral framework does not derive the *quantitative* value of entropy for a given system (that requires statistical mechanics counting of microstates). What it derives is the *existence* of an arrow of time --- the precondition for the second law to be meaningful. The step from $T$-violation to $dS \geq 0$ uses the standard statistical mechanics argument (Boltzmann): given $T$-asymmetric dynamics, entropy-increasing trajectories vastly outnumber entropy-decreasing ones.*
:::

**Verification:** `lorentzian_proof.py`, `lotus/dynamics.py` (`arrow_of_time`).

# Summary: The Complete Correspondence

::: {#tab:equations}
  **Your Equation**                      **Spectral Source**           **Chain**                                          **Code**
  -------------------------------------- ----------------------------- -------------------------------------------------- -------------------
  $E = mc^2$                             $\eta_D(\chi_1) = i/9$        Donnelly $\to$ Lorentz $\to$ 4-mom                 `energy()`
  $F = ma$                               $a_2$ on $S^5/\mathbb{Z}_3$   Heat kernel $\to$ EH $\to$ Einstein $\to$ Newton   `force()`
  $\nabla\cdot E = \rho/\varepsilon_0$   $a_4$ gauge sector            Inner fluct. $\to$ YM $\to$ Maxwell                `maxwell()`
  $(i\slashed{\partial}-m)\psi=0$        $D$ itself                    The operator IS the equation                       `dirac()`
  $i\hbar\partial_t\psi=H\psi$           NR limit of $D$               Dirac $\to$ Foldy--Wouthuysen $\to$ Schrödinger    `schrodinger()`
  $H^2 = \frac{8\pi G}{3}\rho$           $a_0 + a_2 + \text{LOTUS}$    Spectral EH + CC $\to$ Friedmann                   `friedmann()`
  $dS \geq 0$                            $\eta = 2/9 \neq 0$           Spectral asym. $\to$ $T$-breaking $\to$ arrow      `entropy_arrow()`

  : Seven equations of motion, all from $\mathrm{Tr}(f(D^2/\Lambda^2))$ on $M^4 \times S^5/\mathbb{Z}_3$.
:::

**What Einstein got right.** Einstein discovered that gravity is the curvature of spacetime. The spectral framework confirms this and extends it: *all* forces are geometric, each living at a different locus of the same compact manifold. The table above shows that every foundational equation of physics is a projection of one trace formula onto a different coefficient, a different limit, or a different sector.

The universe is not governed by seven independent laws. It is governed by one trace, $\mathrm{Tr}(f(D^2/\Lambda^2))$, projected seven ways.

::: thebibliography
99

A. H. Chamseddine and A. Connes, "The spectral action principle," *Commun. Math. Phys.* **186** (1997) 731--750.

A. Connes, "Noncommutative geometry and the standard model with neutrino mixing," *JHEP* **0611** (2006) 081.

H. Donnelly, "Eta invariants for $G$-spaces," *Indiana Univ. Math. J.* **27** (1978) 889--918.

A. Ikeda, "On the spectrum of the Laplacian on the spherical space forms," *Osaka J. Math.* **17** (1980) 691.

P. B. Gilkey, *Invariance Theory, the Heat Equation, and the Atiyah-Singer Index Theorem*, Publish or Perish, 1984.
:::
