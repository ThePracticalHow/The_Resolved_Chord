*The fold field $\phi$ parameterizes the deformation from smooth $S^5$ to $S^5/\mathbb{Z}_3$. The 77 predictions of the main paper are the equilibrium values of $\phi$-dependent functions. This supplement derives the **dynamics**: the LOTUS Lagrangian, the cosmological timeline, the arrow of time, and the inflaton--Higgs unification.*

# Introduction

The main text fixes $\phi = \phi_{\mathrm{lotus}}$ and obtains 77 predictions from spectral geometry. Those predictions are *equilibrium* values---the outcome when the fold field has settled. Here we derive the *dynamics*: how $\phi$ evolves, what drives inflation, when the spectral phase transition occurs, and why the arrow of time emerges from spectral asymmetry.

**Key result.** One field $\phi$, one potential $V(\phi)$, one history. At high energy: $\phi$ drives Starobinsky inflation. At low energy: $\phi$ is the Higgs. The canonical field $H = v_{\max} \phi$ *is* the Higgs field.

# The LOTUS Lagrangian

::: definition
**Definition 1** (Fold field). *The fold field $\phi \in [0,1]$ parameterizes the deformation from smooth $S^5$ ($\phi = 0$) to the orbifold $S^5/\mathbb{Z}_3$ ($\phi = 1$). The equilibrium value is $$\begin{equation}
\phi_{\mathrm{lotus}} = 1 - \frac{\alpha(d_1 + \lambda_1 + K)}{2} = 0.9574,
\end{equation}$$ where $\alpha \approx 1/137$, $d_1 = 6$, $\lambda_1 = 5$, $K = 2/3$.*
:::

::: definition
**Definition 2** (LOTUS Lagrangian). *$$\begin{equation}
\mathcal{L}(\phi) = \frac{v_{\max}^2}{2}\left(\frac{d\phi}{dt}\right)^2 - V(\phi),
\end{equation}$$ with potential $$\begin{equation}
V(\phi) = \frac{\lambda_H}{4}\, v_{\max}^4\, (\phi^2 - \phi_{\mathrm{lotus}}^2)^2.
\end{equation}$$ The canonical field is $H = v_{\max} \phi$; this *is* the Higgs field.*
:::

::: remark
**Remark 1**. *At $\phi = \phi_{\mathrm{lotus}}$, $V = 0$ (tree-level minimum). The quartic coupling $\lambda_H$ is fixed by the Higgs mass prediction $m_H/m_p = 1/\alpha - 7/2$.*
:::

# $\phi$-Dependent Universe

All spectral-derived quantities become $\phi$-dependent. At $\phi_{\mathrm{lotus}}$, they recover the 77 predictions exactly.

## Coupling and asymmetry

$$\begin{equation}
\alpha(\phi) = \frac{2(1-\phi)}{d_1 + \lambda_1 + K}, \qquad
\eta(\phi) = \frac{d_1}{p^n}\left(\frac{\phi}{\phi_{\mathrm{lotus}}}\right)^3,
\end{equation}$$ with $\eta$ normalized so that $\eta(\phi_{\mathrm{lotus}}) = 2/9$.

## Effective parameters

$$\begin{equation}
K(\phi) = 1 - (1-K)\phi^2, \qquad
d_{1,\mathrm{eff}}(\phi) = d_1 \phi^2, \qquad
G(\phi) = \lambda_1 \eta(\phi).
\end{equation}$$

::: theorem
**Theorem 1** (Recovery at equilibrium). *At $\phi = \phi_{\mathrm{lotus}}$, all 77 predictions of the main paper are recovered exactly: $\alpha = 1/137.038$, $\eta = 2/9$, and the full dictionary of masses, mixings, CKM, PMNS, gravity, and CC.*
:::

# The Cosmological Timeline: Three Epochs

## Epoch 1: Inflation

The spectral action $R^2$ term drives Starobinsky inflation [@starobinsky1980]. $$\begin{equation}
N = \frac{3025}{48} \approx 63 \quad \text{(e-folds)},
\end{equation}$$ matching CMB observations ($n_s \approx 0.968$). The inflaton is $\phi$ at high energy.

## Epoch 2: Spectral Phase Transition

At $\phi_c = 0.60$, the spectral phase transition occurs:

- Ghost decoupling: modes at $\ell = 1$ decouple from the low-energy spectrum.

- Baryogenesis: $\eta_B \propto \alpha^4 \eta$ from spectral asymmetry.

- Dark matter freeze-out: $\Omega_{\mathrm{DM}}/\Omega_B = 16/3$ from $\mathbb{Z}_3$ counting.

## Epoch 3: Electroweak Settlement

$\phi$ settles at $\phi_{\mathrm{lotus}}$. All 77 predictions lock in. The Higgs VEV $v = v_{\max} \phi_{\mathrm{lotus}}$ sets the electroweak scale.

# The Arrow of Time

::: proposition
**Proposition 1** (Spectral arrow). *$\eta(0) = 0$ (smooth $S^5$, no asymmetry) and $\eta(\phi_{\mathrm{lotus}}) = 2/9$. The growth of spectral asymmetry $\eta(\phi)$ from 0 to $2/9$ *is* the arrow of time.*
:::

::: remark
**Remark 2** (Three consequences).

1.  ***CP violation:** $\eta \neq 0$ implies T-violation in the spectral action.*

2.  ***Baryogenesis:** The asymmetry $\eta$ seeds the baryon asymmetry $\eta_B$.*

3.  ***T-symmetry breaking:** $\eta$ cannot be reversed without reversing the fold deformation.*
:::

# The Cosmological Constant as One-Loop Correction

::: theorem
**Theorem 2** (Tree-level cancellation). *$V(\phi_{\mathrm{lotus}}) = 0$ at tree level. The $\mathbb{Z}_3$ equidistribution cancels heavy KK modes in the one-loop effective potential.*
:::

::: proposition
**Proposition 2** (One-loop CC). *The surviving contribution comes from the lightest mode $m_{\nu_3}$: $$\begin{equation}
V_{\mathrm{1-loop}} = \Lambda_{\mathrm{CC}} = \left(m_{\nu_3} \cdot \frac{32}{729}\right)^4.
\end{equation}$$ Numerically: $\Lambda^{1/4} \approx 2.22$ meV ($1.4\%$ of observed CC).*
:::

# The Lorentzian Signature Theorem

::: theorem
**Theorem 3** (Lorentzian Signature from Spectral Asymmetry). *$\eta_D(\chi_1) = i/9$ is purely imaginary (Donnelly, $n=3$ odd, $\mathbb{Z}_3$ complex characters). The imaginary direction maps to the time direction (Osterwalder--Schrader reconstruction). $\mathbb{C}$ has one imaginary axis $\Rightarrow d_{\mathrm{time}}=1 \Rightarrow$ signature $(3,1)$. The time dimension count: $d_{\mathrm{time}} = \dim(\mathrm{center}\, U(3)/\mathbb{Z}_3) = 1$.*
:::

::: remark
**Remark 3** (Status). *Theorem. The proof chain: (1) $\eta_D = i/9$ purely imaginary \[Donnelly\]; (2) $\mathbb{Z}_3$ characters complex because $\omega = e^{2\pi i/3} \neq \bar\omega$ \[algebra\]; (3) Wick rotation maps imaginary to time \[Osterwalder--Schrader\]; (4) $\dim_{\mathrm{Im}}(\mathbb{C}) = 1 \Rightarrow d_{\mathrm{time}} = 1$ \[algebra\]; (5) uniqueness selects $p=3$ (complex) over $p=2$ (real, no time) \[Theorem\]. Verification: `lorentzian_proof.py`. The Connes--Chamseddine spectral action [@connes1996] all align.*
:::

# The Inflaton--Higgs Unification

::: theorem
**Theorem 4** (One field, one potential, one history).

- ***High energy:** $\phi$ at large values drives Starobinsky $R^2$ inflation.*

- ***Low energy:** $\phi$ at $\phi_{\mathrm{lotus}}$ is the Higgs field $H = v_{\max}\phi$.*

- ***Single potential:** $V(\phi) = (\lambda_H/4) v_{\max}^4 (\phi^2 - \phi_{\mathrm{lotus}}^2)^2$.*

*The inflaton and the Higgs are the same degree of freedom at different epochs.*
:::

# Verification Scripts

The following Python scripts verify the dynamics:

- `lotus_dynamics.py` --- $\phi$-dependent functions, equilibrium check

- `lotus_eom.py` --- Equations of motion, phase transition

- `lotus_arrow.py` --- $\eta(\phi)$ evolution, arrow of time

- `lotus_cc_oneloop.py` --- One-loop CC from $m_{\nu_3}$

- `lotus_signature.py` --- Lorentzian signature verification

- `lorentzian_proof.py` --- Full proof: $\eta_D = i/9 \Rightarrow$ signature $(3,1)$

- `sheet_music_spectral.py` --- Two-stave score: spatial eigenvalues (masses) + temporal eigenvalues (decay rates). Tests: $\tau_n = 899$ s ($2.3\%$), $\tau_{\pi^\pm} = 2.70 \times 10^{-8}$ s ($3.5\%$), $\tau_\mu = 2.19 \times 10^{-6}$ s ($0.5\%$). CKM matrix identified as the temporal channel.

# Provenance Table

  **Result**                                 **Source**                        **Verification**                   **Status**
  ------------------------------------------ --------------------------------- ---------------------------------- -------------
  $\phi_{\mathrm{lotus}} = 0.9574$           $\alpha(d_1{+}\lambda_1{+}K)/2$   Exact                              Theorem
  LOTUS Lagrangian                           Higgs potential + fold            EOM                                Definition
  $\eta(\phi)$ arrow                         Donnelly $\eta$, $\phi^3$         Script                             Proposition
  $N = 63$ e-folds                           $R^2$ Starobinsky                 CMB $n_s$                          Theorem
  $\phi_c = 0.60$ transition                 Ghost decoupling                  Script                             Definition
  $\Lambda^{1/4} = m_{\nu_3} \cdot 32/729$   One-loop, $Z_3$ cancel            $1.4\%$                            Theorem
  Inflaton = Higgs                           Same $\phi$ field                 Unification                        Theorem
  Sheet Music (temporal)                     $\mathrm{Im}(\eta_D) = 1/9$       $\tau_n$, $\tau_\pi$, $\tau_\mu$   Framework

  : Provenance map for Supplement XIII.

::: thebibliography
99

H. Donnelly, "Eta invariants for $G$-spaces," *Indiana Univ. Math. J.* **27** (1978) 889--918.

M. Atiyah, V. K. Patodi, and I. M. Singer, "Spectral asymmetry and Riemannian geometry," *Math. Proc. Cambridge Phil. Soc.* **77** (1975) 43--69.

A. A. Starobinsky, "A new type of isotropic cosmological model without singularity," *Phys. Lett. B* **91** (1980) 99--102.

A. Connes and A. H. Chamseddine, "The spectral action principle," *Commun. Math. Phys.* **186** (1997) 731--750.
:::
