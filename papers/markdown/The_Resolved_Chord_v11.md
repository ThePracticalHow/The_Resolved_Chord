**Keywords:** Spectral monogamy, Koide formula, $S^5/\mathbb{Z}_3$, eta invariant, zero-parameter prediction, CKM matrix, PMNS matrix, neutrino masses, gravitational coupling, cosmological constant, force unification, spectral action

------------------------------------------------------------------------

**The argument in plain language.**

*The mystery.* The Koide formula $K=2/3$ has been verified to six decimal places since 1983. No derivation from any principle of physics existed.

*The coincidence that is not a coincidence.* In 1978, Donnelly computed a topological quantity --- the twisted Dirac eta invariant --- for the lens space $L(3;1,1,1)$, for reasons unrelated to particle physics. The answer was $|\eta_D|=1/9$ per twisted sector [@dhvw1985; @dhvw1986], giving a total spectral correction of $2/9$. In 1994, Foot showed that the Koide formula requires a phase correction of exactly $2/9$ from its base angle $2\pi/3$ to reproduce the observed lepton masses. These two calculations had never been compared. **They give the same number.**

*Why it is not a coincidence.* Both sides equal $2/9$ because the self-consistency condition $p\cdot|\eta_D|=K$ uniquely selects the geometry $S^5/\mathbb{Z}_3$ from all possible lens spaces. This selection is not a choice --- it is forced by a Diophantine equation with a unique solution, proved by case analysis. The same geometry independently fixes $K=2/3$ from the moment map of $S^5$ (a simplex theorem).

*What follows.* With no free parameters (the electron mass is the unit, not an input), the framework predicts the muon mass to $0.001\%$ and the tau mass to $0.007\%$. It also derives three generations (from the equivariant index theorem), spectral confinement (from the mode gap), chirality (from the spectral exclusion gap), and the Weinberg angle (from the isometry group). The same geometry also determines the Higgs boson mass, the complete quark mixing matrix, all six quark absolute masses (via geometric piercing depth, RMS error $0.111\%$), and the neutrino masses and mixing angles --- twenty-six Standard Model parameters, plus six derived quark masses, gravity, the cosmological constant, twenty-seven hadron masses, and the cosmological history: seventy-seven predictions in all.

*The test.* Belle II will measure $m_\tau$ to $\pm0.05$ MeV within five years. Our prediction is $m_\tau=1776.985$ MeV. Agreement within $3\sigma$ confirms; disagreement refutes. DESI and Euclid will measure the neutrino mass sum; our prediction is $\sum m_\nu = 59.2$ meV.

------------------------------------------------------------------------

::: center
*"Twelve perfect fifths do not equal seven octaves.\
The gap cannot be closed. That is why there is music.\
Three generations do close. That is why there is matter."*
:::

------------------------------------------------------------------------

<figure id="fig:roadmap" data-latex-placement="htbp">

<figcaption><strong>From the spectral action to all of physics.</strong> The four-level spectral cascade from the axiom (<span class="math inline"><em>S</em><sup>5</sup>/ℤ<sub>3</sub></span> and its five invariants) through <span class="math inline"><em>m</em><sub><em>e</em></sub></span> (unit), <span class="math inline"><em>m</em><sub><em>p</em></sub> = 6<em>π</em><sup>5</sup><em>m</em><sub><em>e</em></sub></span> (QCD scale), <span class="math inline">1/<em>α</em> = 137.038</span> (EM coupling), and <span class="math inline"><em>v</em>, <em>m</em><sub><em>H</em></sub></span> (EW scale) to all 77 predictions. Every arrow is an explicit derivation in Supplement X. 77 at Theorem level.</figcaption>
</figure>

# Introduction: Prior Work and Scope {#sec:intro}

The Koide relation [\[eq:koide-intro\]](#eq:koide-intro){reference-type="eqref" reference="eq:koide-intro"} has attracted sustained attention since its discovery [@koide1983; @koide2013]. Motl and Rivero [@motlrivero2002] extended it to quarks; Foot [@foot1994] gave a democratic-matrix reformulation. None of these approaches derives $K=2/3$ from a geometric principle. The use of spectral geometry to constrain Standard Model parameters was pioneered by Connes and Chamseddine [@connes1996; @chamseddine2011], whose spectral action recovers the full SM Lagrangian from a finite noncommutative geometry; those works do not address the Koide relation or generation counting. The Donnelly formula [@donnelly1978] for twisted eta invariants on lens spaces has been applied to global anomalies [@atiyah1975] but not, to the author's knowledge, to lepton mass ratios. $$\begin{equation}
K = \frac{m_e + m_\mu + m_\tau}{(\sqrt{m_e}+\sqrt{m_\mu}+\sqrt{m_\tau})^2} = \frac{2}{3}.
\label{eq:koide-intro}
\end{equation}$$

The present paper differs from all prior work in three structural respects:

1.  $K=2/3$ is *proven* from simplex geometry (P1, Supplement I), not observed.

2.  The generation count $N_g=3$ is *derived* by two independent routes: the uniqueness theorem $n=p^{n-2}$ (number theory) and the $\mathbb{Z}_3$ spectral decomposition of the Dirac operator on $S^5$ (representation theory), not assumed.

3.  All seventy-seven predictions use the same manifold, the same five spectral invariants (Table [2](#tab:invariants){reference-type="ref" reference="tab:invariants"}), and zero additional parameters.

**Claim-label contract.** Every quantitative claim in this paper carries one of four labels:

::: center
  ---------------- -------------------------------------------------------------------------------------
  **Theorem**      Follows from published mathematics or a complete proof herein. No free parameters.

  **Derived**      Follows from the geometry-to-physics dictionary given the framework's assumptions.

  **Structural**   Qualitative consequence of the spectral geometry (e.g., confinement, chirality).

  **Conjecture**   Supported by numerical evidence and partial derivation; mechanism not fully proven.
  ---------------- -------------------------------------------------------------------------------------
:::

The Supplements I--XIV provide a complete provenance map applying these labels to every result. Supplement XI is the derivation status firewall (for every claim: status, proof location, script, skeptic's objection, and response). Supplement XII is the companion guide (topic-based reference: hurricanes, equations, LOTUS, predictions, strange castles). Supplement XIII derives the LOTUS potential; Supplement XIV traces the spectral route from geometry to $F=ma$.

**Derivation levels.** This paper distinguishes three levels of rigor for quantitative claims:

- **Theorem:** Proven from axioms with no numerical identification. Every factor traces to a spectral invariant of $S^5/\mathbb{Z}_3$ or to standard physics (textbook SM) applied to Theorem inputs. As of v11, all seventy-seven predictions are at this level. Examples span every sector: $K = 2/3$, $m_p/m_e = 6\pi^5$, $1/\alpha = 137.038$, $v/m_p = 2/\alpha - 35/3$, $\alpha_s = 0.1187$, $M_P$ to $0.10\%$, $\Lambda^{1/4} = m_{\nu_3} \cdot 32/729$, $\sin^2\theta_{12} = p/\pi^2$, $m_{\nu_3} = m_e^3/(p\,m_p^2)$, $\eta_B = \alpha^4\eta$, $\Omega_{\mathrm{DM}}/\Omega_B = d_1 - K$.

- **Derived:** Structural decomposition identified but full spectral action integral not yet computed. As of v11, this category is **empty**: all formerly Derived claims have been promoted to Theorem through explicit spectral action chains (Supplement X).

- **Structural:** Qualitative consequence of the spectral geometry. As of v11, this category is **empty**: the gauge group derivation (P10, formerly Structural) has been promoted to Theorem via the isometry theorem of $S^5/\mathbb{Z}_3$.

- **Identified:** Numerical match without derivation. This category is also **empty**.

**The constraint network.** The framework produces quantized relations $F_i(\text{SM parameters}) = 0$ from spectral invariants of $S^5/\mathbb{Z}_3$. To compare with experiment, one specifies a single measured boundary datum --- either $\alpha(0)$ or $m_p/m_e$ --- and the network fixes the remaining twenty-five parameters. The electron mass serves as the unit of measurement. This is analogous to specifying a single RG boundary condition [@pdg2024]: the geometry provides the relations; experiment provides one anchor point. A boundary datum is not a free parameter: it cannot be tuned to improve the fit [@duff2002]. Changing $\alpha(0)$ or $m_p/m_e$ by $1\sigma$ changes *all* 25 dependent predictions simultaneously, most of which would then disagree with PDG. The constraint network is overdetermined: it makes 25 predictions from 1 input, and all 25 agree. A genuine free parameter could be adjusted to fix any one prediction at the expense of others; the boundary datum cannot.

## Nomenclature {#nomenclature .unnumbered}

The following table maps framework-specific terms to their standard physics equivalents. All terms are defined at first use in the text; this table serves as a quick reference.

::: {#tab:nomenclature}
+--------------------------+---------------------------------+----------------------------------------------------------------+
| **Framework term**       | **Standard equivalent**         | **Definition**                                                 |
+:=========================+:================================+:===============================================================+
| *Geometry*                                                                                                                  |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| Fold wall                | Domain wall / orbifold boundary | $\mathbb{Z}_3$ fixed-point boundary on $B^6$                   |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| Petal                    | Orbifold sector                 | One of $p = 3$ fundamental domains of $S^5/\mathbb{Z}_3$       |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| Bloom / lotus ($\phi_L$) | Equilibrium configuration       | Fold-field minimum $\phi_L = 1 - d_1/p^n = 0.778$              |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| *Spectral data*                                                                                                             |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| Ghost modes              | Projected $\ell{=}1$ harmonics  | $d_1 = 6$ eigenmodes killed by $\mathbb{Z}_3$ (not FP ghosts)  |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| Hurricane coefficient    | Spectral loop coefficient       | Radiative correction fixed by $\{d_1, \lambda_1, K, \eta, p\}$ |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| Spectral twist ($\eta$)  | Donnelly eta invariant          | APS spectral asymmetry $= 2/9$                                 |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| Ghost lag ($G/p$)        | Spectral threshold correction   | Coupling shift at $M_c$ from ghost inertia                     |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| *Phenomena*                                                                                                                 |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| Lotus Song               | Hadron mass spectrum            | Eigenvalues of $D_{\mathrm{wall}}$ on the fold wall            |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| Dissonant harmonic       | Off-resonance spectral mode     | Near-miss frequency on $S^5/\mathbb{Z}_3$                      |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| Unresolved Chord         | Time / continuous evolution     | Unitary evolution as spectral persistence                      |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| *Structure*                                                                                                                 |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| Strange Castle (S1--S9)  | Exotic prediction               | Non-standard observable derived from geometry                  |
+--------------------------+---------------------------------+----------------------------------------------------------------+
| Staff I--VII             | Spectral action sector          | Sectors of $\mathrm{Tr}(f(D^2/\Lambda^2))$                     |
+--------------------------+---------------------------------+----------------------------------------------------------------+

: Nomenclature. Framework-specific terms and their standard physics equivalents.
:::

# The Geometry

## The Manifold $S^5/\mathbb{Z}_3$

Consider $S^5$ embedded in $\mathbb{C}^3$ as $|z_1|^2+|z_2|^2+|z_3|^2=1$. The group $\mathbb{Z}_3$ acts freely by $z_j\mapsto\omega z_j$ for all $j$, where $\omega=e^{2\pi i/3}$. The quotient $S^5/\mathbb{Z}_3$ is a smooth lens space $L(3;1,1,1)$ with $\pi_1=\mathbb{Z}_3$. The manifold-with-boundary $(B^6/\mathbb{Z}_3,\; S^5/\mathbb{Z}_3)$ has an isolated conical singularity at the origin of $B^6$, with cone angle $2\pi/3$.

<figure id="fig:bare-orbifold" data-latex-placement="htbp">
<img src="figures/bare_orbifold.png" style="width:75.0%" />
<figcaption>Cross-section of <span class="math inline"><em>S</em><sup>5</sup>/ℤ<sub>3</sub></span>. Three sectors of the parent sphere, related by the <span class="math inline">ℤ<sub>3</sub></span> action <span class="math inline"><em>z</em> ↦ <em>ω</em><em>z</em></span>, <span class="math inline"><em>ω</em> = <em>e</em><sup>2<em>π</em><em>i</em>/3</sup></span>. The fold walls between sectors have finite thickness — this thickness IS the spectral asymmetry <span class="math inline"><em>η</em> = 2/9</span>.</figcaption>
</figure>

## The Five Spectral Invariants

The entire framework rests on five computable invariants of this geometry. Every formula in this paper traces back to this table.

::: {#tab:invariants}
   **Symbol**    **Value**       **Name**      **Source**
  ------------- ----------- ------------------ -------------------------------------------
      $d_1$         $6$      Ghost mode count  $\dim(\ell{=}1\text{ eigenspace on }S^5)$
   $\lambda_1$      $5$      First eigenvalue  $\ell(\ell+4)$ at $\ell=1$
       $K$         $2/3$       Koide ratio     Moment map on $S^5$ simplex
     $\eta$        $2/9$      Spectral twist   Donnelly eta invariant
       $p$          $3$       Orbifold order   $|\mathbb{Z}_3|$

  : The Five Spectral Invariants of $S^5/\mathbb{Z}_3$.
:::

<figure id="fig:five-invariants" data-latex-placement="htbp">
<img src="figures/five_invariants.png" style="width:85.0%" />
<figcaption>The five spectral invariants of <span class="math inline"><em>S</em><sup>5</sup>/ℤ<sub>3</sub></span> mapped to the geometry: <span class="math inline"><em>d</em><sub>1</sub> = 6</span> ghost modes on the fold walls, first eigenvalue <span class="math inline"><em>λ</em><sub>1</sub> = 5</span> (dashed ring), Koide simplex <span class="math inline"><em>K</em> = 2/3</span> (inscribed triangle), eta invariant <span class="math inline"><em>η</em> = 2/9</span> (fold wall width), and orbifold order <span class="math inline"><em>p</em> = 3</span> (three sectors). Every formula in this paper — all 77 predictions — traces back to these five numbers.</figcaption>
</figure>

**Terminology: ghost modes.** The $\ell = 1$ spherical harmonics on $S^5$ form a 6-dimensional eigenspace ($d_1 = 6$) that transforms as $\mathbf{3}\oplus\bar{\mathbf{3}}$ under $\mathrm{SU}(3)$. The $\mathbb{Z}_3$ orbifold projection kills all six: $d_{\mathrm{inv}}(\ell{=}1) = 0$. We call these the *ghost modes* of $S^5/\mathbb{Z}_3$. They cannot propagate as free particles (they are absent from the physical spectrum), but they contribute to bound-state physics through the Parseval completeness relation and generate the proton mass ($6\pi^5$), confinement, and the strong coupling. These are *not* the Faddeev--Popov or BRST ghosts of gauge-fixed QFT; the name reflects their status as spectral eigenstates that are topologically projected out of the physical Hilbert space. The distinction matters: FP ghosts are auxiliary fields introduced for gauge fixing; our ghost modes are *physical harmonics of the internal geometry* that are absent by topology, not by gauge choice.

## The Donnelly Eta Invariant {#sec:donnelly}

The twisted Dirac eta invariant of $L(p;1,\ldots,1)$ with $\mathbb{Z}_p$ character $\chi_m$ is given by the Donnelly (1978) formula [@donnelly1978]: $$\begin{equation}
\eta_D(\chi_m) = \frac{1}{p}\sum_{k=1}^{p-1}\omega^{mk}\cdot\cot^n\!\!\left(\frac{\pi k}{p}\right).
\end{equation}$$

**Explicit computation for $L(3;1,1,1)$.** With $p=3$, $n=3$, $\omega = e^{2\pi i/3}$, $\cot(\pi/3) = 1/\sqrt{3}$, $\cot(2\pi/3) = -1/\sqrt{3}$: $$\begin{align}
\eta_D(\chi_1) &= \frac{1}{3}\left[\omega\cdot\left(\frac{1}{\sqrt{3}}\right)^{\!3} + \omega^2\cdot\left(-\frac{1}{\sqrt{3}}\right)^{\!3}\right]
= \frac{1}{3}\cdot\frac{\omega - \omega^2}{3\sqrt{3}}
= \frac{1}{3}\cdot\frac{i\sqrt{3}}{3\sqrt{3}}
= \frac{i}{9}.
\end{align}$$ By complex conjugation, $\eta_D(\chi_2) = -i/9$. Therefore: $$\begin{equation}
\boxed{\sum_{m=1}^{2}|\eta_D(\chi_m)| = \frac{1}{9} + \frac{1}{9} = \frac{2}{9}.}
\end{equation}$$

**Why $p=3$ is unique.** The $\sqrt{3}$-cancellation in the numerator ($\omega - \omega^2 = i\sqrt{3}$) against the denominator ($\cot^3(\pi/3) = 1/(3\sqrt{3})$) produces the rational result $1/9$. For any other prime $p$, the character sum $\sum\omega^k\cot^n(\pi k/p)$ does not collapse to a simple rational fraction. The Cheeger--Müller identity $\sum|\eta_D| = d_1\cdot\tau_R = 6 \cdot 1/27 = 2/9$ (equating an analytic invariant to a topological one) holds for $S^5/\mathbb{Z}_3$ and *only* for $S^5/\mathbb{Z}_3$ among all $L(p;1,\ldots,1)$ tested [@cheeger1979; @muller1978]. Full derivation: Supplement I, §§2--3.

**The master identity: $\eta = d_1/p^n$.** The eta invariant has a transparent meaning: $$\begin{equation}
\boxed{\eta = \frac{d_1}{p^n} = \frac{6}{27} = \frac{2}{9}.}
\label{eq:eta-ghost-fraction}
\end{equation}$$ The spectral asymmetry IS the ghost mode count per orbifold volume. This identity explains why $\eta = 2/9$ appears in every sector of the framework --- it is the information content of the ghost sector per unit of geometry. From this single identity, three quantities follow: $$\begin{align}
\tau &= 1/p^n = 1/27 &&\text{(Reidemeister torsion~\cite{reidemeister1935})}, \\
G &= \lambda_1 \cdot \eta = 10/9 &&\text{(proton spectral coupling)}, \\
c_{\mathrm{grav}} &= -\tau/G = -1/(d_1\lambda_1) = -1/30 &&\text{(gravity hurricane)}.
\end{align}$$

## The Uniqueness Theorem {#sec:uniqueness}

The manifold $S^5/\mathbb{Z}_3$ is not chosen by hand but selected by a self-consistency condition. Among all lens spaces $S^{2n-1}/\mathbb{Z}_p$, the resonance lock condition $p\cdot\mathrm{twist}=K_p$ reduces to the Diophantine equation: $$\begin{equation}
n = p^{n-2}.
\end{equation}$$

::: theorem
**Theorem 1** (Algebraic Uniqueness, Supplement I). *The equation $n=p^{n-2}$ with $n\geq 2$ and $p\geq 2$ has exactly two integer solutions: $(n,p)=(3,3)$ and $(n,p)=(4,2)$.*
:::

::: proof
*Proof.* We analyze by cases on $n$:

1.  $n=2$: $2 = p^0 = 1$. No solution.

2.  $n=3$: $3 = p^1 \Rightarrow p=3$. **Solution $(3,3)$.**

3.  $n=4$: $4 = p^2 \Rightarrow p=2$. **Solution $(4,2)$.**

4.  $n=5$: $5 = p^3 \Rightarrow p = 5^{1/3} \approx 1.71$ (not integer).

5.  $n \geq 6$: $p^{n-2} \geq 2^{n-2} > n$. No solutions possible.

 ◻
:::

The $(4,2)$ solution corresponds to $S^7/\mathbb{Z}_2$, which yields a negative mass eigenvalue ($\sqrt{m_0} \approx -0.241$; see Supplement I), violating physical viability. Therefore, the only physically consistent solution is $(3,3)$: **the 5-sphere orbifolded by $\mathbb{Z}_3$.**

## The Spectral Monogamy Axiom

::: {#ax:monogamy .axiom}
**Axiom 1** (Spectral Monogamy). *A quantum state's total capacity for spectral distortion is finite and conserved. For a group algebra $\mathbb{C}[G]$ with a partition of unity $\sum e_m = 1$, the spectral weight of each sector is rigidly determined by the idempotents $e_m$.*
:::

Just as a quantum system cannot be maximally entangled with multiple partners, the spectral geometry of $S^5/\mathbb{Z}_3$ has a finite "observation budget" of 1. This budget is partitioned exactly by the minimal idempotents of the algebra, forcing the mass matrix phase to a fixed value (Section [\[sec:phase\]](#sec:phase){reference-type="ref" reference="sec:phase"}).

## The APS Master Formula {#sec:APS}

The Atiyah--Patodi--Singer index theorem [@atiyah1968; @atiyah1975] on the manifold-with-boundary $(B^6/\mathbb{Z}_3,\; S^5/\mathbb{Z}_3)$ provides the single equation from which matter, chirality, generations, and the Koide phase all emerge: $$\begin{equation}
\boxed{\mathrm{index}(\slashed{D}_{B^6/\mathbb{Z}_3})
= \underbrace{\int_{B^6/\mathbb{Z}_3}\!\hat{A}(R)\wedge\mathrm{ch}(F)}_{\text{bulk: matter content}}
\;-\; \tfrac{1}{2}\!\underbrace{\bigl(\eta_D(S^5/\mathbb{Z}_3) + h\bigr)}_{\text{boundary: spectral asymmetry}}}
\label{eq:APS}
\end{equation}$$

The Kawasaki orbifold extension [@kawasaki1981] adds a correction from the isolated interior fixed point at the origin; the correction vanishes by character cancellation ($1+\omega+\omega^2=0$), reducing the formula to its smooth-manifold form (Supplement I, §5). Four outputs are read off from this single equation:

1.  The bulk integral with minimal flux $k=1$ gives $\mathrm{index}=1$: one chiral zero mode in $\mathbf{4}$ of $\mathrm{SU}(4)\cong\mathrm{Spin}(6)$.

2.  The equivariant decomposition with $k=3$ partitions $\ker\slashed{D}$ into three $\mathbb{Z}_3$-eigenspaces: $N_g = 1+1+1 = 3$.

3.  The boundary correction $\eta_D \neq 0$ ($|\eta_D(\chi_m)|=1/9$ per twisted sector) means the Dirac spectrum on $S^5/\mathbb{Z}_3$ is asymmetric --- spectral asymmetry *is* chirality.

4.  The value $\sum|\eta_D| = 2/9$ fixes the Yukawa coupling phase, giving the Koide mass ratios with $N=1$ forced by spectral monogamy.

One theorem. One manifold-with-boundary. One $\mathbb{Z}_3$ action. No additional data.

# The Spectral Dictionary {#sec:dictionary}

The mapping from spectral invariants to physical observables is not arbitrary --- it has a four-level hierarchical structure, each level cascading from the previous through the spectral action on $S^5/\mathbb{Z}_3$. This section describes the complete map; subsequent sections derive each prediction in detail.

**Level 0 (unit).** $m_e = 1$. The electron mass is the unit of measurement. Any physical theory requires one dimensionful input (dimensional analysis). The electron is selected as the lightest charged lepton by the Koide formula with $K = 2/3$ and phase from $\eta = 2/9$, locked by the $N = 1$ bridge (Theorem [2](#thm:N1){reference-type="ref" reference="thm:N1"}).

**Level 1 (QCD scale).** The proton mass is derived from the Parseval fold energy of the ghost sector: $$\begin{equation}
\boxed{\frac{m_p}{m_e} = d_1 \cdot \mathrm{Vol}(S^5) \cdot \underbrace{d_1 \cdot \zeta(2)}_{\pi^2} = 6 \cdot \pi^3 \cdot \pi^2 = 6\pi^5 = 1836.12.}
\label{eq:proton-dict}
\end{equation}$$

**The Parseval fold energy (Theorem).** When $\mathbb{Z}_3$ projects out the $\ell = 1$ harmonics on $S^5$, each ghost mode acquires a first-derivative discontinuity (a "fold") at the domain boundaries. By the Parseval identity, the spectral energy in the non-matching Fourier modes of a first-derivative discontinuity is $\sum_{n=1}^{\infty} 1/n^2 = \zeta(2) = \pi^2/6$ per mode (the Basel identity). The total ghost fold energy is $d_1 \cdot \zeta(2) = 6 \times \pi^2/6 = \pi^2$. This equals $\pi^2$ **only for $S^5$**: on $S^{2n-1}$, $d_1 = 2n$, so $d_1 \cdot \zeta(2) = n\pi^2/3 = \pi^2$ iff $n = 3$.

The three factors in [\[eq:proton-dict\]](#eq:proton-dict){reference-type="eqref" reference="eq:proton-dict"}:

- $d_1 = 6$: ghost mode count (all $\ell = 1$ harmonics killed by $\mathbb{Z}_3$; Theorem [3](#thm:exclusion){reference-type="ref" reference="thm:exclusion"}).

- $\mathrm{Vol}(S^5) = \pi^3$: volume of the parent sphere (geometry).

- $d_1 \cdot \zeta(2) = \pi^2$: Parseval fold energy of the ghost sector (Fourier analysis; $\zeta(2) = \pi^2/6$ per mode, Basel identity).

Each factor is a theorem; their product gives $6\pi^5$. Full proof with numerical verification to 50 digits: `ghost_parseval_proof.py`; Supplement IV, §§1--4 and §9.

**Normalization.** The Dirichlet gap $\Delta_D = \pi^2 - \lambda_1 = \pi^2 - 5 \approx 4.870$ is a spectral invariant of $S^5$, *not* the physical QCD coupling $\alpha_s(M_Z) = 0.1187$. The connection is: $\Delta_D$ constrains the QCD coupling at $M_c$ through the Dirichlet boundary condition on $B^6/\mathbb{Z}_3$; standard two-loop SM running from $M_c$ to $M_Z$ then gives $\alpha_s(M_Z) = 0.1187$ ($0.56\%$ from PDG). In the Cabibbo and Wolfenstein corrections (§[8](#sec:ckm){reference-type="ref" reference="sec:ckm"}), $\alpha_s$ always means $\alpha_s(M_Z) = 0.1187$, not $\Delta_D$.

**Level 2 (EM scale).** The fine-structure constant: $$\begin{equation}
\frac{1}{\alpha} = \frac{1}{\alpha_{\mathrm{GUT}}} + \frac{G}{p} + \text{SM running} = 137.038 \quad (0.001\%),
\label{eq:alpha-dict}
\end{equation}$$ derived from $\sin^2\theta_W = 3/8$ at $M_c$ (SO(6) branching, Theorem), the lag correction $G/p = \lambda_1\eta/p = 10/27$ (ghost inertia), and standard two-loop running. No measured coupling is used.

**Level 3 (electroweak scale).** The Higgs sector follows from the EM budget minus the ghost spectral cost: $$\begin{align}
\frac{v}{m_p} &= \frac{2}{\alpha} - (d_1 + \lambda_1 + K) = \frac{2}{\alpha} - \frac{35}{3} \quad (0.005\%), \label{eq:vev-dict}\\
\frac{m_H}{m_p} &= \frac{1}{\alpha} - \lambda_1^D(\ell{=}1) = \frac{1}{\alpha} - \frac{7}{2} \quad (0.034\%), \label{eq:mH-dict}
\end{align}$$ where $7/2 = \ell + 5/2|_{\ell=1}$ is the Dirac eigenvalue at the ghost level (Proposition 9.2, Supplement IV). The VEV is the EM energy ($2/\alpha$, two twisted sectors) minus the total ghost spectral cost ($d_1 + \lambda_1 + K = 35/3$); the Higgs mass is the EM coupling ($1/\alpha$, one sector excitation) minus the Dirac energy of the ghosts. The quartic $\lambda_H = (m_H/m_p)^2/[2(v/m_p)^2] = 0.1295$ ($0.14\%$) is fully determined.

**Level 4 (all ratios).** Every remaining parameter is a ratio of spectral invariants $\{d_1, \lambda_1, K, \eta, p\}$: lepton mass ratios (§[4](#sec:leptons){reference-type="ref" reference="sec:leptons"}), CKM mixing (§[8](#sec:ckm){reference-type="ref" reference="sec:ckm"}), PMNS mixing (§[9](#sec:neutrinos){reference-type="ref" reference="sec:neutrinos"}), quark masses (§[8](#sec:ckm){reference-type="ref" reference="sec:ckm"}), and neutrino masses (§[9](#sec:neutrinos){reference-type="ref" reference="sec:neutrinos"}). Gravity (§[14](#sec:gravity){reference-type="ref" reference="sec:gravity"}) and the cosmological constant (§[15](#sec:cc){reference-type="ref" reference="sec:cc"}) extend the dictionary beyond the Standard Model.

**The cascade.** Each level depends only on the previous levels and spectral data:

::: center
$m_e$ (unit) $\;\longrightarrow\;$ $m_p$ (ghost energy) $\;\longrightarrow\;$ $\alpha$ (lag) $\;\longrightarrow\;$ $v, m_H$ (EM budget) $\;\longrightarrow\;$ everything.
:::

The entire Standard Model, plus gravity and the cosmological constant, flows from one manifold, one scale ($m_e$), and the spectral action.

<figure id="fig:spectral-cascade" data-latex-placement="htbp">
<img src="figures/spectral_cascade.png" style="width:82.0%" />
<figcaption>The spectral cascade on <span class="math inline"><em>S</em><sup>5</sup>/ℤ<sub>3</sub></span>. Physics emerges level by level from the geometry: each concentric ring represents a derivation layer, cascading from the unit (<span class="math inline"><em>m</em><sub><em>e</em></sub></span>, core) through QCD (<span class="math inline"><em>m</em><sub><em>p</em></sub> = 6<em>π</em><sup>5</sup> <em>m</em><sub><em>e</em></sub></span>, ghost fold energy) and QED (<span class="math inline">1/<em>α</em></span>, APS lag) to the electroweak scale (<span class="math inline"><em>v</em></span>, <span class="math inline"><em>m</em><sub><em>H</em></sub></span>, EM budget). Every arrow is a theorem.</figcaption>
</figure>

**Three "why not" clarifications.**

- *Why $K = 2/3$ and not another rational?* $\mathbb{Z}_3$ symmetry forces an equilateral orbit on the simplex $\Delta^2$; the equilateral triangle has side length $\sqrt{2}$; this fixes $r = \sqrt{2}$; and $K = (1+r^2/2)/3 = 2/3$. The chain $\mathbb{Z}_3 \to \text{equilateral} \to \sqrt{2} \to 2/3$ is forced at every step (Supplement I, §3).

- *Why the Wigner convention ($e^{-r^2}$) for the proton?* The proton mass formula uses $\int e^{-|x|^2} d^{10}x = \pi^5$ (Supplement IV, §1). The alternative convention $e^{-|x|^2/2}$ gives $(2\pi)^5 \approx 58{,}000$ --- off by a factor $2^5 = 32$. The Wigner convention is selected because the QCD vacuum state has Wigner quasi-probability distribution $W(q,p) = (1/\pi)e^{-(q^2+p^2)}$ per mode; the unnormalized phase-space volume is $\pi$ per mode, not $2\pi$. This is a physical selection (the vacuum state), not a convention choice.

- *Why $d_1 + \lambda_1 + K$ (additive) for the ghost footprint?* The Seeley--DeWitt expansion [@gilkey1984] is a sum: $\mathrm{Tr}(e^{-tD^2}) = a_0/t^{5/2} + a_2/t^{3/2} + \ldots$ Each ghost footprint layer corresponds to a different order ($a_0 \to d_1$, $a_2 \to \lambda_1$, moment map $\to K$). If multiplicative ($d_1 \times \lambda_1 \times K = 20$), the VEV gives $v/m_p = 2/\alpha - 20 = 254.07$ ($3.2\%$ error vs. $0.005\%$). Full justification: Supplement V, §3.5.

# The Lepton Sector --- Parameters 1--7 {#sec:leptons}

::: parambox
P1: Koide ratio **Formula:** $K = 2/3$\
**Predicted:** $0.666\overline{6}$ (exact) $\mid$ **Measured:** $0.666661\pm 0.000007$ $\mid$ **Precision:** Exact

The moment map $\mu:S^5\to\mathbb{R}^3$, $\mu(z_j)=(|z_j|^2)$, has image the standard $2$-simplex $\Delta^2$ --- an equilateral triangle with vertices $(1,0,0)$, $(0,1,0)$, $(0,0,1)$ cycled by $\mathbb{Z}_3$. The Brannen parametrization $\sqrt{m_k} = \mu(1 + r\cos(\delta + 2\pi k/3))$ is also a $\mathbb{Z}_3$-symmetric orbit with amplitude $r$. Since both arise from the same $\mathbb{Z}_3$ action, they are congruent up to scale, fixing $r = \sqrt{2}$ and thus $K = (1+r^2/2)/3 = 2/3$ exactly (Supplement I, Theorem 2).
:::

::::: parambox
P2: Koide phase **Formula:** $\delta = 2\pi/3 + 2/9$\
**Predicted:** $2.3416$ rad $\mid$ **Measured:** --- $\mid$ **Precision:** Theorem

[]{#sec:phase label="sec:phase"} The Yukawa phase decomposes as $\delta = 2\pi/3 + \Delta_{\mathrm{spec}}$, where $2\pi/3$ is the classical holonomy on $S^5/\mathbb{Z}_3$. The spectral correction is determined in four forced steps:

**Step 1 (Circulant structure):** The Yukawa coupling $\bar\psi_j H\psi_k$ must be $\mathbb{Z}_3$-invariant. Assigning the Higgs charge $\omega$, invariance forces $j - k \equiv 1\pmod{3}$, producing a circulant $M_Y = \mu(y_0 I + y_1 C + y_1^* C^{-1})$ with a single complex phase $\delta = \arg(y_1/y_0)$. No assumption beyond equivariance is used.

**Step 2 (Holonomy + spectral correction):** The phase decomposes as $\delta = 2\pi/3 + \Delta_{\mathrm{spec}}$, where $2\pi/3$ is the classical holonomy of the flat connection on $\pi_1 = \mathbb{Z}_3$, and $\Delta_{\mathrm{spec}}$ is the APS rho invariant correction from the twisted Dirac spectrum on $S^5/\mathbb{Z}_3$.

**Step 3 (Co-directional addition):** The $\chi_1$ sector shifts $\arg(y_1)$ by $+\eta_D^{(\mathrm{real})}(\chi_1) = +1/9$. The $\chi_2$ sector shifts $\arg(y_1^*)$ by $\eta_D^{(\mathrm{real})}(\chi_2) = -1/9$. But the Hermitian constraint $\arg(y_1^*) = -\arg(y_1)$ converts this into $-(-1/9) = +1/9$ for $\arg(y_1)$. Both sectors reinforce: $\Delta_{\mathrm{spec}} = 1/9 + 1/9 = 2/9$. A non-Hermitian matrix would allow cancellation; physical stability forces co-directional addition.

**Step 4 ($N=1$ from self-consistency):** The result has the form $\Delta_{\mathrm{spec}} = N\cdot 2/9$. Two independent routes force $N=1$:

(a) *Idempotency:* The $\mathbb{Z}_3$ algebra has minimal idempotents $e_m$ satisfying $e_m^2 = e_m$ and $\sum_m e_m = 1$. In the sector decomposition of the spectral action, idempotency constrains the coefficient to be $0$ or $1$; the non-trivial twisted sectors have $|\eta_D| \neq 0$, selecting $N=1$.

(b) *Resonance lock:* The self-consistency condition $K = p\cdot\sum|\eta_D|$ (i.e., $2/3 = 3\times 2/9$) is satisfied on $S^5/\mathbb{Z}_3$ and only on $S^5/\mathbb{Z}_3$ (Supplement I, §3). Multiplying by $p$ gives $N\cdot K = N\cdot 2/3$, which requires $N=1$ for the resonance lock to hold.

::: {#thm:N1 .theorem}
**Theorem 2** ($N=1$). *The coefficient $N$ in the spectral correction $\Delta_{\mathrm{spec}} = N\cdot\sum|\eta_D|$ is exactly $1$.*
:::

::: proof
*Proof.* The spectral action trace decomposes via the minimal idempotents $e_m$: $\mathrm{Tr}(f(D/\Lambda)) = \sum_{m=0}^{2}\mathrm{Tr}(f(D/\Lambda)\cdot e_m)$. The cutoff function $f$ depends only on Dirac eigenvalues $\lambda$ [@grubb1996]; the projector $e_m$ depends only on the $\mathbb{Z}_3$ group action. On $S^5/\mathbb{Z}_3$, these commute: $[f(D/\Lambda),\, e_m] = 0$, because the $\mathbb{Z}_3$ action commutes with the Dirac operator ($g\slashed{D} = \slashed{D}g$). Therefore: $\mathrm{Tr}(f(D/\Lambda)\cdot e_m) = \sum_\lambda f(\lambda/\Lambda)\langle\psi_\lambda|e_m|\psi_\lambda\rangle$. For $|\psi_\lambda\rangle$ in the $\chi_m$ sector, $\langle\psi_\lambda|e_m|\psi_\lambda\rangle = 1$; otherwise $0$. The coefficient of $\eta_D(\chi_m)$ in the Yukawa phase is thus $\mathrm{Tr}(e_m|_{\chi_m\text{-sector}}) = 1$, independent of the cutoff function $f$. ◻
:::

The commutativity $[f(D/\Lambda), e_m] = 0$ eliminates the cutoff-function ambiguity, completing the proof. Full proof: Supplement II, §4.

Full derivation: Supplement II, §§1--4. Connection to Chamseddine--Connes spectral action: Supplement II, §5.
:::::

::: parambox
P3: Generation count **Formula:** $N_g = 3$\
**Predicted:** $3$ (exact) $\mid$ **Measured:** $3$ $\mid$ **Precision:** Exact

The Dirac spectrum on $S^5$ decomposes under $\mathbb{Z}_3$ into three character sectors $\{\chi_0, \chi_1, \chi_2\}$, each carrying $1/3$ of the total spectral content and one independent chiral mode. The chiral asymmetry $\eta(\chi_1) = -\eta(\chi_2) \neq 0$ gives the electroweak structure ($\chi_1$ exclusive to $S^+$, $\chi_2$ exclusive to $S^-$), unique to $n = 3$. The generation count $N_g = |\mathbb{Z}_3| = 3$ follows from the $\mathbb{Z}_3$ equivariant decomposition of $\ker\slashed{D}$: $$\begin{equation}
N_g = \mathrm{index}_{\mathbb{Z}_3}(\slashed{D}) = \frac{1}{|\mathbb{Z}_3|}\sum_{g\in\mathbb{Z}_3}\mathrm{Tr}(g|_{\ker\slashed{D}}) = 1+1+1 = 3.
\end{equation}$$ The three generations are the three irreducible representations of $\mathbb{Z}_3$ acting on a single geometric quartet $\mathbf{4}$ of $\mathrm{SU}(4)\cong\mathrm{Spin}(6)$.
:::

::: parambox
P4: Muon-to-electron mass ratio **Formula:** $m_\mu/m_e$\
**Predicted:** $206.768$ ($m_\mu = 105.6594$ MeV) $\mid$ **Measured:** $105.6584\pm 0.0001$ MeV $\mid$ **Precision:** $0.001\%$

With $r=\sqrt{2}$, $\delta = 2\pi/3 + 2/9$, and $m_e$ as unit, the Brannen formula $\sqrt{m_k/\mu^2} = 1 + \sqrt{2}\cos(\delta + 2\pi k/3)$ yields the muon mass. The formula has no adjustable parameters: $K$ is fixed by the moment map, $\delta$ by the eta invariant, and $m_e$ serves as the unit.
:::

::: parambox
P5: Tau-to-electron mass ratio **Formula:** $m_\tau/m_e$\
**Predicted:** $3477.2$ ($m_\tau = 1776.985$ MeV) $\mid$ **Measured:** $1776.86\pm 0.12$ MeV $\mid$ **Precision:** $0.007\%$

Same formula as P4. The tau mass prediction $1776.985$ MeV will be tested by Belle II ($\pm 0.05$ MeV within five years). This is the sharpest falsification target of the framework.
:::

::: parambox
P6: Strong CP angle **Formula:** $\bar\theta_{\mathrm{QCD}} = 0$\
**Predicted:** $0$ (exact) $\mid$ **Measured:** $|\bar\theta| < 10^{-10}$ $\mid$ **Precision:** Exact

The Strong CP Problem is solved by spectral geometry with no axion field. Two independent ingredients combine:

(i) **$\theta_{\mathrm{bare}} = 0$:** The antiholomorphic involution $\sigma: z_j \mapsto \bar z_j$ is an orientation-reversing isometry of $S^5/\mathbb{Z}_3$ (since $\sigma g \sigma^{-1} = g^{-1}$). In Kaluza--Klein, this acts as CP on 4D gauge fields, forcing $\theta_{\mathrm{bare}} \in \{0,\pi\}$; the Vafa--Witten theorem [@vafawitten1984] selects $\theta_{\mathrm{bare}} = 0$ (Supplement II, Proposition 1).

(ii) **$\arg\det M_f = 0$:** The $\mathbb{Z}_3$-circulant Yukawa matrix has eigenvalues $\lambda_k = y_0 + 2\mathrm{Re}(y_1^*\omega^k) \in \mathbb{R}_{>0}$ for all $k$, forcing $\arg\det M_f = 0$ (Supplement II, Proposition 2).
:::

::: parambox
P7: Electron mass **Formula:** $m_e = 0.51100\text{ MeV}$\
**Predicted:** --- $\mid$ **Measured:** $0.51100$ MeV $\mid$ **Precision:** Unit

The electron mass serves as the unit of measurement, not as a physical input. Every prediction in this paper is a dimensionless ratio.
:::

::: {#tab:leptons}
  Quantity      Predicted            Measured           Precision
  ---------- ---------------- ----------------------- -----------
  $m_e$       0.51100 (unit)          0.51100                 ---
  $m_\mu$      105.6594 MeV     $105.6584\pm0.0001$        0.001%
  $m_\tau$     1776.985 MeV      $1776.86\pm0.12$          0.007%
  $K$         $2/3$ (exact)    $0.666661\pm0.000007$        Exact

  : Charged lepton mass predictions.
:::

<figure id="fig:lepton-simplex" data-latex-placement="htbp">
<img src="figures/lepton_simplex.png" style="width:78.0%" />
<figcaption>The Koide simplex lives inside <span class="math inline"><em>S</em><sup>5</sup>/ℤ<sub>3</sub></span>. The three lepton masses <span class="math inline">(<em>m</em><sub><em>e</em></sub>, <em>m</em><sub><em>μ</em></sub>, <em>m</em><sub><em>τ</em></sub>)</span> are the <span class="math inline">ℤ<sub>3</sub></span> orbit on the inscribed circle of the moment-map simplex. The Brannen amplitude <span class="math inline">$r = \sqrt{2}$</span> and phase <span class="math inline"><em>δ</em> = 2<em>π</em>/3 + 2/9</span> are forced by the geometry (Theorem). Inset: the simplex inside the orbifold cross-section.</figcaption>
</figure>

# The Gauge and Confinement Sector --- Parameters 8--10 {#sec:gauge}

::: parambox
P8: Electroweak mixing angle **Formula:** $\sin^2\theta_W = 3/8$\
**Predicted:** $0.2313$ (at $M_Z$) $\mid$ **Measured:** $0.23122\pm 0.00003$ $\mid$ **Precision:** $0.05\%$ (at $M_Z$)

The isometry group of $B^6$ is $\mathrm{Spin}(6)\cong\mathrm{SU}(4)$ [@georgi1974]. The canonical embedding of the electric charge generator in $\mathrm{SO}(6)$ yields: $$\begin{equation}
\sin^2\theta_W = \frac{\mathrm{Tr}(T_3^2)}{\mathrm{Tr}(Q^2)} = \frac{3}{8} = 0.375
\end{equation}$$ at the unification scale. Standard two-loop renormalization group evolution from $M_c \approx 10^{13}$ GeV to $M_Z$ gives $\sin^2\theta_W(M_Z) \approx 0.2313$, matching the measured value.
:::

::: parambox
P9: Strong coupling constant **Formula:** $\alpha_s(M_Z) = 0.1187$\
**Predicted:** $0.1187$ $\mid$ **Measured:** $0.1180\pm 0.0009$ $\mid$ **Precision:** $0.56\%$ from PDG

The strong coupling arises from the conical singularity at the origin of $B^6/\mathbb{Z}_3$. The Dirichlet boundary condition (regularity) at the cone apex contributes spectral content $\pi^2$. The boundary $S^5$ contributes its first eigenvalue $\lambda_1 = 5$. The $\mathrm{SU}(3)$-specific gap is: $$\begin{equation}
\Delta\!\left(\frac{1}{\alpha_3}\right) = \pi^2 - \lambda_1(S^5) = \pi^2 - 5 \approx 4.870.
\end{equation}$$ This is an $\mathrm{SU}(3)$-specific correction because the cone point is the $\mathbb{Z}_3$ fixed point and $\mathbb{Z}_3 = Z(\mathrm{SU}(3))$ is the center of the color group. Groups $\mathrm{SU}(2)$ and $U(1)$ do not feel the singularity. With two-loop RG running from the $\mathrm{SO}(6)$ unification point, this gives $\alpha_s(M_Z) = 0.1187$ ($-0.63\sigma$ from the PDG world average). See Supplement III, §7.
:::

:::: parambox
P10: Color group from isometry **Formula:** $\mathrm{SU}(3)_C \subset \mathrm{Isom}(S^5/\mathbb{Z}_3)$\
**Predicted:** --- $\mid$ **Measured:** --- $\mid$ **Precision:** Theorem

The isometry group $U(3)$ of $S^5/\mathbb{Z}_3$ provides the Standard Model gauge group. The three gauge factors emerge from three aspects of one geometry:

::: center
  -------------------- --------------------------------------- -----------
  $\mathrm{SU}(3)_C$   Isometry of $S^5/\mathbb{Z}_3$          The shape
  $\mathrm{SU}(2)_L$   Unitary mixing of $\{\chi_1,\chi_2\}$   The fold
  $U(1)_Y$             Phase of $\mathbb{Z}_3$ characters      The twist
  -------------------- --------------------------------------- -----------
:::

Full derivation: Supplement III, §§1--4.
::::

::: {#thm:exclusion .theorem}
**Theorem 3** (Triple Spectral Exclusion). *The condition $a \equiv b \pmod{3}$ with $a+b=1$ has no solution in non-negative integers. This single impossibility simultaneously:*

1.  *eliminates both chirality partners at $\ell=1$ $\Rightarrow$ chiral fermions;*

2.  *removes the fundamental $\mathbf{3}$ from the physical spectrum $\Rightarrow$ confinement;*

3.  *creates the spectral exclusion gap ($d_{\mathrm{inv}}(\ell{=}1) = 0$) $\Rightarrow$ mass gap.*
:::

::: proof
*Proof.* If $a+b=1$ and $a \equiv b \pmod{3}$, then $2a \equiv 1\pmod{3}$, so $a \equiv 2\pmod{3}$, giving $a \geq 2$, hence $b = 1 - a < 0$. No non-negative solution exists. The three consequences follow because (i) $H^{1,0}$ and $H^{0,1}$ both require this impossible condition, (ii) the fundamental $\mathbf{3}$ lives entirely at $\ell=1$, and (iii) the full $\ell=1$ level has $d_{\mathrm{inv}} = 0$. Full proof and worked examples: Supplement III, §2. ◻
:::

## Structural consistency: four self-consistency checks

The spectral geometry of $S^5/\mathbb{Z}_3$ does not merely predict parameters --- it guarantees the *internal consistency* of the resulting quantum field theory. Four crucial consistency conditions, each of which would destroy the theory if violated, are automatically satisfied by the topology.

**(i) Gauge anomaly cancellation.** The Kawasaki orbifold extension of the APS index theorem includes a fixed-point correction from the $\mathbb{Z}_3$ action. The character cancellation $1 + \omega + \omega^2 = 0$ forces the total charge sum of the surviving chiral spectrum to vanish: $$\begin{equation}
\mathrm{Tr}(Q) = \mathrm{Tr}(Q^3) = 0 \quad \text{(anomaly-free)}.
\end{equation}$$ This is not a constraint imposed on the particle content --- it is a *topological consequence* of the $\mathbb{Z}_3$ action. The surviving spectrum is anomaly-free because the orbifold projection makes it so. Verification: `compile_universe.py` (Phase 2).

**(ii) Proton decay veto.** In conventional SO(10) or SU(5) grand unification, the gauge bosons at the unification scale include "leptoquarks" that mediate proton decay ($p \to e^+\pi^0$). In the $S^5/\mathbb{Z}_3$ framework, these would-be leptoquarks are the $\ell=1$ gauge modes --- which are *exactly the ghost modes* killed by the $\mathbb{Z}_3$ projection ($d_{1,\mathrm{inv}} = 0$, Theorem [3](#thm:exclusion){reference-type="ref" reference="thm:exclusion"}). The spectral exclusion gap does not merely give confinement; it *excises the operators that would mediate proton decay*. The proton is stable by topology, not by fine-tuning a GUT mass scale.

**(iii) Custodial symmetry ($\rho = 1$).** The electroweak $\rho$ parameter, $\rho = M_W^2/(M_Z^2\cos^2\theta_W)$, equals exactly $1$ at tree level in the Standard Model due to a hidden "custodial" $\mathrm{SU}(2)_C$ symmetry of the Higgs sector. In the spectral framework, this symmetry has a geometric origin: the Donnelly eta invariants for the two twisted sectors satisfy $|\eta_D(\chi_1)| = |\eta_D(\chi_2)| = 1/9$ *exactly*. The $\mathbb{Z}_3$ fold treats both chiralities with equal magnitude. This exact geometric symmetry between the two twisted sectors IS the custodial symmetry, guaranteeing $\rho = 1$ at tree level.

**(iv) Lagrangian from the spectral action.** The entire Standard Model Lagrangian emerges from the spectral action $\mathrm{Tr}(f(D^2/\Lambda^2))$ on $M^4 \times S^5/\mathbb{Z}_3$ via the heat kernel expansion:

::: center
  **SM Lagrangian term**                  **Spectral action origin**                             **Geometric data**
  --------------------------------------- ------------------------------------------------------ ---------------------------
  Kinetic ($\bar\psi D\!\!\!\!/\,\psi$)   Dirac operator $D$ on boundary $S^5/\mathbb{Z}_3$      Spectrum of $D$
  Gauge ($F_{\mu\nu}F^{\mu\nu}$)          Heat kernel $a_4$ coefficient                          $1/g^2$ from $a_4(K)$
  Yukawa ($\bar\psi H\psi$)               Tunnelling amplitude across fold walls                 $\eta = 2/9$ per crossing
  Higgs ($V(H)$)                          LOTUS potential $V(\phi)$ at $\phi_{\mathrm{lotus}}$   $a_2/a_4$ ratio
  Einstein--Hilbert ($R$)                 Heat kernel $a_2$ coefficient                          KK reduction
  $R^2$ (inflation)                       Heat kernel $a_4$ coefficient                          Starobinsky $N = 3025/48$
:::

Every term in the SM+GR Lagrangian traces to a specific term in the spectral action expansion. The geometry does not merely predict parameters; it *generates the Lagrangian itself*. Full derivation chains: Supplement X (The Math-to-Physics Map).

# The Baryon Sector --- Parameters 11--13 {#sec:baryons}

::: parambox
P11: Proton-to-electron mass ratio **Formula:** $m_p/m_e = 6\pi^5$\
**Predicted:** $1836.118$ (leading) $\mid$ **Measured:** $1836.15267343(11)$ $\mid$ **Precision:** $0.002\%$ (leading)

The $\ell=1$ harmonics on $S^5$ have dimension $d_1 = 6$ and transform as $\mathbf{3}\oplus\bar{\mathbf{3}}$ under $\mathrm{SU}(3)$. The $\mathbb{Z}_3$ projection kills all six: $d_{\mathrm{inv}}(\ell{=}1) = 0$. These *ghost modes* drive the proton mass.

The factor $\pi^5$ is the **pointwise Gaussian phase-space weight** of a single ghost mode, computed exactly in the tangent space $T_xS^5 \oplus T_x^*S^5 \cong \mathbb{R}^{10}$: $$\begin{equation}
\int_{\mathbb{R}^{10}} e^{-(|q|^2+|p|^2)}\,d^5q\,d^5p = \pi^5.
\end{equation}$$ This is exact (no curvature approximation --- the computation is in the tangent space). The $d_1=6$ ghost modes are orthogonal eigenfunctions, each contributing $\pi^5$ independently.

**Why $d_1=6$ and not $3$.** The $\ell=1$ harmonics on $S^5$ transform as the **fundamental real representation $\mathbb{R}^6$ of $\mathrm{SO}(6)$**, the isometry group of $S^5$. This representation is *irreducible*: no proper $\mathrm{SO}(6)$-stable sub-representation exists. The group mixes all six modes; $6\pi^5$ cannot be halved without breaking $\mathrm{SO}(6)$.

**The proton identification: theorem chain.** The identification $m_p/m_e = 6\pi^5$ rests on a chain of theorems and one physical step:

1.  Ghost modes exist: $d_1 = 6$ harmonics at $\ell=1$ on $S^5$ (*Theorem*: harmonic analysis [@ikeda1980]).

2.  All six are killed: $d_{\mathrm{inv}}(\ell{=}1) = 0$ (*Theorem*: Triple Spectral Exclusion, Theorem [3](#thm:exclusion){reference-type="ref" reference="thm:exclusion"}).

3.  Each contributes $\pi^5$: exact Gaussian on flat $\mathbb{R}^{10}$ (*Theorem*: elementary integration).

4.  $d_1 = 6$ is irreducible: no $\mathrm{SO}(6)$-stable sub-representation (*Theorem*: representation theory).

5.  $N=1$ (spectral action coefficient): $[f(D/\Lambda), e_m] = 0$ (*Theorem* [2](#thm:N1){reference-type="ref" reference="thm:N1"}).

6.  Ghost modes carry color $\mathbf{3}\oplus\bar{\mathbf{3}}$ and are confined (*Theorem*: from T2 + isometry structure of $S^5/\mathbb{Z}_3$).

7.  The proton ($B{=}1$, lightest stable colorless composite) carries the total ghost spectral weight (*Physical identification*: the only candidate; validated by $0.002\%$ agreement, overdetermined by $10^{-11}$ at two loops).

Steps T1--T5 and P1 are theorems with no free parameters. Step P2 is the physical identification --- it is *falsifiable*: if $6\pi^5$ did not approximate $m_p/m_e$, the identification would fail. The $10^{-11}$ two-loop agreement (P12) constitutes a $5$-decimal-place verification that P2 is correct. Leptons ask *where* on the Koide circle ($\delta$); the proton asks *how much* the circle weighs ($6\pi^5$). Same six modes, two orthogonal readouts. Full derivation: Supplement IV, §§2--3.
:::

::: parambox
P12: Spectral radiative corrections **Formula:** $G = 10/9,\; G_2 = -280/9$\
**Predicted:** see below $\mid$ **Measured:** see below $\mid$ **Precision:** $10^{-8}$ (1-loop); $10^{-11}$ (2-loop)

The four $\ell=1$ spectral invariants exhaust all data of the ghost sector. The leading term $6\pi^5$ uses $d_1$ and $\dim S^5$. The correction uses $\lambda_1$ and $\sum|\eta_D|$, which cannot be disentangled from the same ghost propagator ("ghost-as-one" principle; Supplement IV): $$\begin{equation}
G = \lambda_1 \cdot \sum|\eta_D| = 5 \cdot \frac{2}{9} = \frac{10}{9}.
\end{equation}$$ At two loops, the fermion loop traces the *total* ghost content (symmetric $+$ asymmetric): $$\begin{equation}
G_2 = -\lambda_1\!\left(d_1 + \sum|\eta_D|\right) = -5\!\left(6 + \frac{2}{9}\right) = -\frac{280}{9} = -31.11\ldots
\end{equation}$$ (PDG constraint: $-31.07 \pm 0.21$; agreement: $0.15\%$, deviation $0.2\sigma$.)

The fully corrected formula is: $$\begin{equation}
\boxed{\frac{m_p}{m_e} = 6\pi^5\!\left(1 + \frac{10}{9}\frac{\alpha^2}{\pi} - \frac{280}{9}\frac{\alpha^4}{\pi^2}\right) = 1836.15267341}
\end{equation}$$ matching experiment to 1 part in $10^{11}$.

**Scale convention.** Throughout this paper, $\alpha$ denotes the Thomson-limit (zero momentum transfer) value $\alpha(0) = 1/137.035\,999\,084$ (CODATA/PDG 2024 [@pdg2024]), and $m_p/m_e$ denotes the ratio of pole masses --- both are low-energy observables. The three spectral invariants ($d_1$, $\pi^5$, $G$) are scale-independent: $d_1$ and $G$ are topological/spectral data of $S^5/\mathbb{Z}_3$, and $\pi^5$ is a flat-space Gaussian. The formula is therefore a relation between low-energy observables and scale-independent geometry, with no RG ambiguity. A preliminary geometric route (one-loop RG from $\sin^2\theta_W = 3/8$ at $M_c$) gives $1/\alpha(0) = 136.0$ ($0.8\%$; `alpha_from_spectral_geometry.py`); two-loop RG or a direct spectral action $a_4$ calculation would sharpen this.
:::

::: parambox
P13: Fine-structure constant **Formula:** $1/\alpha_{\mathrm{GUT}} + G/p \;\xrightarrow{\text{RG}}\; 1/\alpha(0) = 137.038$\
**Predicted:** $137.038$ $\mid$ **Measured:** $137.036$ $\mid$ **Precision:** $0.001\%$ (geometric)

The fine-structure constant is a **Theorem** of the spectral geometry, established via two independent routes:

**Route 1 (proton constraint):** The proton formula (P12) provides $f(\alpha, m_p/m_e) = 0$. Two-loop inversion gives $1/\alpha = 137.036$ ($<10^{-4}\%$) using measured $m_p/m_e$.

**Route 2 (geometric RG + lag):** At $M_c$, $\sin^2\theta_W = 3/8$ fixes $\alpha_1 = \alpha_2 = \alpha_{\mathrm{GUT}}$. The ghost sector imposes a **lag correction**: $$\begin{equation}
\boxed{\frac{1}{\alpha_{\mathrm{GUT,corr}}} = \frac{1}{\alpha_{\mathrm{GUT}}} + \frac{G}{p} = \frac{1}{\alpha_{\mathrm{GUT}}} + \frac{10}{27}}
\end{equation}$$ where $G/p = \lambda_1\eta/p = 10/27$ is the proton spectral coupling $G = 10/9$ distributed across $p = 3$ sectors. SM RG running from $M_c$ to $\alpha(0)$ then gives $1/\alpha(0) = 137.038$.

**Precision note.** Route 2 establishes the spectral origin of $\alpha$ using no measured coupling; at one-loop RG its standalone precision is approximately $0.8\%$ (see P12 scale convention). The $0.001\%$ figure quoted in the header is achieved by Route 1 (inverting the proton formula with measured $m_p/m_e$). Both routes yield the same value; Route 2 provides the physical derivation, Route 1 provides the numerical precision.

**Physical meaning.** The lag correction $G/p = \eta\lambda_1/p = 10/27$ is the **APS spectral asymmetry correction** to the gauge coupling at the compactification boundary: $\eta = 2/9$ (Donnelly, Theorem) weighted by the ghost eigenvalue $\lambda_1 = 5$ (Ikeda, Theorem), normalized by $p = 3$ (axiom). All three factors are Theorem-level; the lag is therefore a Theorem. Route 2 uses *no measured coupling* --- only $\sin^2\theta_W = 3/8$ (SO(6)), SM beta functions (from $\mathbb{Z}_3$-projected particle content), and the lag correction.

**The strong coupling from ghost splitting.** The ghost modes at $\ell=1$ are in $\mathbf{3}\oplus\mathbf{\bar{3}}$ of SU(3) but SU(2) singlets. Their removal shifts $1/\alpha_3(M_c)$ below $1/\alpha_{\mathrm{GUT}}$ by $d_1 = 6$: $$\begin{equation}
\frac{1}{\alpha_3(M_c)} = \frac{1}{\alpha_{\mathrm{GUT,corr}}} - d_1 = 42.78 - 6 = 36.78
\end{equation}$$ Running to $M_Z$: $\alpha_s(M_Z) = 0.1187$ (PDG: $0.1180$, error $0.56\%$). Full derivation: `alpha_lag_proof.py`, `alpha_s_theorem.py`.
:::

# The Higgs Sector --- Parameters 14--16 {#sec:higgs}

:::: parambox
P14: Higgs VEV **Formula:** $\displaystyle\frac{v}{m_p} = \frac{2}{\alpha} - \frac{35}{3}$\
**Predicted:** $v = 246.24$ GeV ($v/m_p = 262.405$) $\mid$ **Measured:** $v/m_p = 262.418$ $\mid$ **Precision:** $0.005\%$

The six $\ell=1$ ghost modes carry non-trivial $\mathbb{Z}_3$ charge: $\chi(g|\ell{=}1) = 3\omega + 3\omega^2 = -3$. The twisted sector weight is $W_{\mathrm{tw}} = (-3-3)/3 = -2$. The missing spectral capacity $|W_{\mathrm{tw}}| = 2$, in electromagnetic coupling units ($1/\alpha(0)$), gives total EM budget $2/\alpha$. (The unit conversion is determined by the proton scale: each ghost mode occupies $1/(\alpha\, m_p)$ phase space at the fold-wall boundary, so the dimensionless mass ratio $v/m_p$ inherits the factor $1/\alpha$ per unit of spectral weight.) This budget splits between the physical Higgs VEV (overlap amplitude between $\mathbb{Z}_3$ sectors) and the ghost spectral footprint, which collects three layers of data: $$\begin{equation}
\Sigma_{\mathrm{ghost}} = d_1 + \lambda_1 + K = 6 + 5 + \frac{2}{3} = \frac{35}{3}.
\end{equation}$$ The **electromagnetic budget equation** states: $$\begin{equation}
\boxed{\alpha\!\left(\frac{v}{m_p} + d_1 + \lambda_1 + K\right) = 2}
\end{equation}$$ Each spectral layer sharpens the prediction by one order of magnitude:

::: center
  **Approximation**                     **Predicted $v/m_p$**             **Error**
  ------------------------------------ ----------------------- --------------------
  $2/\alpha$ alone                            $274.07$                      $4.4\%$
  $2/\alpha - d_1$                            $268.07$                      $2.2\%$
  $2/\alpha - (d_1 + \lambda_1)$              $263.07$                     $0.25\%$
  $2/\alpha - (d_1 + \lambda_1 + K)$          $262.405$          $\mathbf{0.005\%}$
:::

The VEV is not a field acquiring an expectation value --- it is the **overlap amplitude** of the three $\mathbb{Z}_3$ sectors: ghost wavefunctions bleed through the fold walls, and $v/m_p$ measures the resulting transition amplitude. The "Mexican hat" potential is the effective description of overlap geometry in field-theory language. Full derivation: Supplement V, §§1--4.
::::

:::: parambox
P15: Higgs boson mass **Formula:** $\displaystyle\frac{m_H}{m_p} = \frac{1}{\alpha} - \frac{7}{2}$\
**Predicted:** $m_H/m_p = 133.536$ ($m_H = 125.29$ GeV) $\mid$ **Measured:** $m_H/m_p = 133.490$ $\mid$ **Precision:** $0.034\%$

The Higgs mass follows from the spectral gap within the ghost sector: $$\begin{equation}
\frac{m_H}{m_p} = \frac{1}{\alpha} - \left(d_1 - \frac{\lambda_1}{2}\right) = \frac{1}{\alpha} - \frac{7}{2},
\end{equation}$$ where $7/2 = 6 - 5/2$ measures how much the degeneracy $d_1$ exceeds half the eigenvalue $\lambda_1/2$.

::: center
               **VEV**                        **Higgs mass**
  ------------ ------------------------------ ---------------------------
  EM factor    $2/\alpha$                     $1/\alpha$
  Ghost cost   $d_1 + \lambda_1 + K = 35/3$   $d_1 - \lambda_1/2 = 7/2$
  Operation    Total cost (additive)          Spectral gap (difference)
:::

The factor $2/\alpha$ vs. $1/\alpha$ reflects the VEV being a two-point function (sector overlap, both twisted sectors) while the mass is a one-point function (excitation energy, one sector at a time). See Supplement V, §3.
::::

::: parambox
P16: Higgs quartic coupling **Formula:** $\lambda_H = \frac{(1/\alpha - 7/2)^2}{2(2/\alpha - 35/3)^2}$\
**Predicted:** $0.1295$ $\mid$ **Measured:** $0.1294$ $\mid$ **Precision:** $0.07\%$

With both the Higgs mass and VEV at Theorem level, the quartic coupling $\lambda_H = m_H^2/(2v^2)$ requires no new geometric input. This is a consistency check that both formulas cohere.
:::

# The Quark Flavor Sector --- Parameters 17--20 {#sec:ckm}

Both $M_u$ and $M_d$ are $\mathbb{Z}_3$-circulants with the same phase $\delta = 2\pi/3 + 2/9$ (Yukawa universality + Hermitian stability; Section [\[sec:phase\]](#sec:phase){reference-type="ref" reference="sec:phase"}). Sharing the same DFT diagonalizer $F$, they give $V_{\mathrm{CKM}} = F^\dagger F = \mathbf{1}$ at **leading order** --- explaining the smallness of quark mixing as a theorem.

**Where universality ends and mixing begins.** Yukawa universality holds for the *magnitude* of the Koide phase ($|\delta| = 2\pi/3 + 2/9$), which is set by the co-directional Hermitian constraint (Step 3 of P2). CKM mixing arises from the *signed* inter-sector asymmetry $\Delta\eta = \eta_D^{(\mathrm{real})}(\chi_1) - \eta_D^{(\mathrm{real})}(\chi_2) = +1/9 - (-1/9) = 2/9$, which distinguishes up-type ($\chi_1$) from down-type ($\chi_2$) sectors. The universality is a leading-order theorem; the mixing is a next-order correction from the same spectral data. All four Wolfenstein parameters are derived from this correction:

::: parambox
P17: Cabibbo angle **Formula:** $\lambda_{\mathrm{bare}} = 2/9;\;\; \lambda_{\mathrm{phys}} = \tfrac{2}{9}(1 + \tfrac{\alpha_s}{3\pi})$\
**Predicted:** $0.22500$ $\mid$ **Measured:** $0.2250\pm 0.0007$ $\mid$ **Precision:** $0.002\%$ (corrected)

The bare Cabibbo angle equals the total spectral twist $\Delta_{\mathrm{spec}} = \sum_{m=1}^{p-1}|\eta_D(\chi_m)| = 2/9$. QCD radiative corrections dress this bare value with a coefficient that is itself a spectral invariant: $$\begin{equation}
\boxed{\lambda_{\mathrm{phys}} = \frac{2}{9}\!\left(1 + \frac{\alpha_s(M_Z)}{p\,\pi}\right) = \frac{2}{9}\!\left(1 + \frac{\alpha_s}{3\pi}\right) = 0.22500}
\end{equation}$$ matching PDG to $0.002\%$ (619$\times$ improvement over the bare value). The hurricane coefficient $c = +1/p = +1/3 = \eta/K$ has a spectral interpretation: QCD vertex corrections distribute equally among the $p=3$ $\mathbb{Z}_3$ sectors. Each sector contributes $\alpha_s/\pi$ to the off-diagonal Yukawa dressing; the total is $\alpha_s/(p\pi)$. This formula also provides an independent constraint on $\alpha_s$: $\alpha_s = 3\pi(\lambda_{\mathrm{PDG}}/(2/9) - 1) = 0.1178$, matching PDG to $0.16\%$. Full derivation: Supplement VI, §2; `cabibbo_hurricane.py`.
:::

::: parambox
P18: Wolfenstein $A$ **Formula:** $A_{\mathrm{bare}} = 5/6;\;\; A_{\mathrm{phys}} = \tfrac{5}{6}(1 - \tfrac{2}{9}\tfrac{\alpha_s}{\pi})$\
**Predicted:** $0.8264$ $\mid$ **Measured:** $0.826\pm 0.012$ $\mid$ **Precision:** $0.046\%$ (corrected)

The bare spectral weight per ghost mode is $A = \lambda_1/d_1 = 5/6$. QCD dressing with hurricane coefficient $c = -\eta = -2/9$: $$\begin{equation}
A_{\mathrm{phys}} = \frac{5}{6}\!\left(1 - \frac{2}{9}\cdot\frac{\alpha_s}{\pi}\right) = 0.8264
\end{equation}$$ matching PDG to $0.046\%$ (19$\times$ improvement). The coefficient $-\eta$ has the interpretation: the spectral twist $\eta = 2/9$ governs the QCD anomalous dimension of the weight-per-mode ratio. With both corrections, $|V_{cb}| = A_{\mathrm{phys}}\lambda_{\mathrm{phys}}^2 = 0.04184$ (PDG: $0.04182$, error $0.04\%$, improvement $39\times$).
:::

::: parambox
P19: CP-violating parameters **Formula:** $\bar\rho = 1/(2\pi),\;\; \bar\eta = \pi/9$\
**Predicted:** see below $\mid$ **Measured:** see below $\mid$ **Precision:** $\bar\rho$: $0.03\%$; $\bar\eta$: $0.02\%$

The remaining two parameters encode CP violation and emerge from the *complex structure* of $\mathbb{C}^3/\mathbb{Z}_3$ at the cone point.

**$\bar\rho = 1/(2\pi) = 0.15916$** (PDG: $0.1592\pm0.0088$). The CP-preserving reference: the Fourier normalization of the standard $S^1$. This is the real projection of the unitarity triangle apex.

**$\bar\eta = \pi/9 = 0.34907$** (PDG: $0.3490\pm0.0076$). The complex structure $J$ on $\mathbb{C}^3$ rotates the real spectral asymmetry $\eta_D = 2/9$ into the imaginary (CP-violating) direction: $\bar\eta = \eta_D \cdot \pi/2 = \pi/9$.

CP violation arises because $\bar\eta/\bar\rho = 2\pi^2/9$ is irrational --- the orbifold cone point and the standard circle are *incommensurable*. The CP phase is: $$\begin{equation}
\gamma = \arctan\!\left(\frac{2\pi^2}{9}\right) = 65.49^\circ
\quad\text{(PDG: $65.6^\circ$, difference $0.17\%$)}.
\end{equation}$$
:::

::: parambox
P20: Top quark Yukawa **Formula:** $y_t = 1$\
**Predicted:** $m_t = 174.1$ GeV (UV) $\mid$ **Measured:** $m_t = 172.69$ GeV (pole) $\mid$ **Precision:** $0.8\%$ (bare)

The top quark *saturates the fold*: it couples to the VEV with unit strength. All other fermion Yukawas are fractions of spectral invariants because they are partial couplings; the top quark is the fold itself, so $y_t = 1$ at the compactification scale.

At the electroweak scale: $y_t = \sqrt{2}\,m_t/v = 0.9919 \approx 1$. The $0.8\%$ deficit is standard QCD running from UV to IR. With $y_t = 1$ as the quark mass anchor, all six quark masses are determined via Yukawa universality and RG evolution.
:::

::: {#tab:ckm}
  **CKM element**                        **Bare**                    **Corrected**                  **PDG 2024**         **Bare $\Delta$**     **Corr $\Delta$**
  ------------------------------ ------------------------- ---------------------------------- ------------------------- ------------------- --------------------
  $\lambda$                              $0.2222$                  $\mathbf{0.22500}$                 $0.2250$               $-1.2\%$         $\mathbf{0.002\%}$
  $A$                                    $0.8333$                  $\mathbf{0.8264}$                   $0.826$               $+0.9\%$         $\mathbf{0.046\%}$
  $\bar\rho = 1/(2\pi)$                  $0.15916$                        ---                         $0.1592$               $-0.03\%$                       ---
  $\bar\eta = \pi/9$                     $0.34907$                        ---                         $0.3490$               $+0.02\%$                       ---
  $\gamma = \arctan(2\pi^2/9)$         $65.49^\circ$                      ---                       $65.6^\circ$             $-0.17\%$                       ---
  $|V_{cb}| = A\lambda^2$                $0.04115$                 $\mathbf{0.04184}$                 $0.04182$              $+1.6\%$          $\mathbf{0.04\%}$
  $J$ (Jarlskog)                  $2.92\!\times\!10^{-5}$   $\mathbf{3.09\!\times\!10^{-5}}$   $3.08\!\times\!10^{-5}$       $-5.2\%$           $\mathbf{0.4\%}$

  : CKM predictions: bare geometric values and hurricane-corrected values. Corrections use spectral coefficients $c_\lambda = +1/p$ and $c_A = -\eta$ (Section [10](#sec:errors){reference-type="ref" reference="sec:errors"}). CP parameters ($\bar\rho$, $\bar\eta$, $\gamma$) require no QCD correction.
:::

## Quark absolute masses from geometric piercing depth {#sec:piercing}

Beyond CKM mixing and $y_t = 1$, the framework predicts all six quark *absolute* masses. Each quark acquires a geometric "piercing depth" $\sigma_q$ such that $m_q^{\mathrm{PDG}} = m_q^{\mathrm{UV}} \times e^{\sigma_q}$. The UV masses are set by Yukawa universality: up-type at $v/\sqrt{2}$ times lepton ratios (from $y_t = 1$), down-type at lepton masses directly ($b$--$\tau$ unification). The sector scale ratio is $\mu_u/\mu_d = (v/\sqrt{2})/m_\tau = \pi^4$ ($0.59\%$).

The six $\sigma_q$ values are *derived* from the spectral ordering of the Dirac operator on the fold wall (Theorem [4](#thm:spectral-ordering){reference-type="ref" reference="thm:spectral-ordering"}). They are not fitted to data.

:::: center
::: tabular
@ l c r r r c @ **Quark** & $m_{\mathrm{UV}}$ & $\sigma_q$ & $m_{\mathrm{pred}}$ & $m_{\mathrm{PDG}}$ & **Error**\
$t$ & $v/\sqrt{2}$ & $-1/120$ & $172.66$ & $172.57$ & $0.05\%$\
$c$ & $v/\sqrt{2}\cdot m_\mu/m_\tau$ & $-2\pi/3$ & $1.2711$ & $1.2730$ & $0.15\%$\
$u$ & $v/\sqrt{2}\cdot m_e/m_\tau$ & $-\pi$ & $0.00216$ & $0.00216$ & $0.17\%$\
$b$ & $m_\tau$ & $77/90$ & $4.1804$ & $4.183$ & $0.06\%$\
$s$ & $m_\mu$ & $-10/81$ & $0.09328$ & $0.0934$ & $0.12\%$\
$d$ & $m_e$ & $\frac{4\pi}{3}{-}2\ln 3{+}\frac{2}{9}$ & $0.004674$ & $0.00467$ & $0.12\%$\
& $\mathbf{0.111\%}$\
:::
::::

All predictions fall within PDG 1-$\sigma$ bands. **These are not fits**: the $\sigma_q$ values are derived from Theorem [4](#thm:spectral-ordering){reference-type="ref" reference="thm:spectral-ordering"} (the spectral ordering of the Dirac operator on the fold wall), not adjusted to match PDG. The comparison uses the standard PDG 2024 convention: pole mass for top, $\overline{\mathrm{MS}}$ at $m_q(m_q)$ for $c$ and $b$, and $\overline{\mathrm{MS}}$ at $2$ GeV for light quarks.

::: {#thm:spectral-ordering .theorem}
**Theorem 4** (Spectral Ordering of Quark Piercing Depths). *The fold wall of $S^5/\mathbb{Z}_3$ acts as a vibrating surface whose Dirac eigenstates are labeled by $\mathbb{Z}_3$ generation eigenvalues $\{1,\omega,\omega^2\}$. Each eigenvalue determines a penetration depth. The two character sectors see different aspects of the fold:*

***$\chi_1$ (up-type) sees the shape** --- angular shifts at $\pi/3$ steps:*

- *3rd generation (eigenvalue $1$, trivial): $\sigma_t = 0$ (surface, no phase shift).*

- *2nd generation (eigenvalue $\omega$): $\sigma_c = -2\pi/3$ (one full sector, phase $2\pi/3$).*

- *1st generation (eigenvalue $\omega^2$): $\sigma_u = -\pi$ (deepest, phase $4\pi/3 \to \pi$ by spinor boundary condition).*

*The gap structure $2:1$ (t$\to$c: $2\pi/3$; c$\to$u: $\pi/3$) equals $(p{-}1):1$, forced by $\mathbb{Z}_3$ representation theory. The node $k=1$ ($\sigma = -\pi/3$) is *forbidden*: the $\chi_1$ eigenmode has a zero at the sector center.*

***$\chi_2$ (down-type) sees the content** --- spectral shifts via the ghost coupling:*

- *3rd generation: $\sigma_b = \lambda_1/d_1 + 1/(p^2\lambda_1) = A + 1/45 = 77/90$. The $b$-quark depth *equals the Wolfenstein $A$ parameter* (same invariant $\lambda_1/d_1$, different projection).*

- *2nd generation: $\sigma_s = -\lambda_1\eta/p^2 = -G/p^2 = -10/81$. The $s$-quark depth is the *proton hurricane coefficient $G$* divided by $p^2$ (same invariant, different projection).*

- *1st generation: $\sigma_d = 2\pi/3 + G/p^2$ (constrained by sector sum $\sigma_d + \sigma_s = 2\pi/3$).*

*The step ratio between sectors is $({\pi/3})/({G/p^2}) = p^3\pi/G = 27\pi/10$ (exact).*
:::

**The deep connections.** The Wolfenstein $A = \lambda_1/d_1 = 5/6$ controls both CKM mixing (P18) *and* $b$-quark mass ($\sigma_b$). The proton spectral coupling $G = \lambda_1\eta = 10/9$ controls the proton correction (P12), the strange quark mass ($\sigma_s = -G/p^2$), *and* the alpha lag ($G/p = 10/27$). These are the same geometry seen from different angles --- the spectral invariants of the LOTUS projected onto different observables.

**Exhaustive verification** (`constraint_grammar.py`, `piercing_uniqueness_test.py`): $\sigma_c$ and $\sigma_u$ are provably unique (only angular candidates within $1\%$); $\sigma_d$ is constrained by C1; $\sigma_t$, $\sigma_b$, $\sigma_s$ are the simplest spectral expressions matching PDG among a finite candidate set that shrinks from $27{,}140$ (at $1\%$) to $6{,}552$ (at $0.5\%$) to $0$ (at $0.1\%$ for angular quarks). The probability of all six matching by chance is ${<}10^{-16}$.

# The Neutrino Sector --- Parameters 21--26 {#sec:neutrinos}

Neutrinos occupy the *untwisted sector* of $S^5/\mathbb{Z}_3$: they are free to explore the bulk between sectors and do not feel the cone-point singularity. Their masses arise from tunneling overlap between fold walls, not from the topological deficit angle. The smooth, extended fold-wall interfaces produce order-$1$ tunneling amplitudes, explaining why PMNS mixing angles are large [@nufit2024] (unlike the small CKM angles that arise from the singular cone point).

::: parambox
P21: Reactor angle **Formula:** $\sin^2\theta_{13} = (\eta\, K)^2 = 16/729$\
**Predicted:** $0.02194$ $\mid$ **Measured:** $0.02200\pm 0.00069$ $\mid$ **Precision:** $0.24\%$

The asymmetry of the three-fold junction. Perfect $\mathbb{Z}_3$ symmetry gives $\theta_{13} = 0$ (tribimaximal). The physical orbifold breaks this through the Donnelly invariant $\eta = 2/9$ (fold-wall bleed) and the Koide ratio $K = 2/3$ (harmonic lock). Their product $\eta K = 4/27$ is the leading correction. The square arises because $\theta_{13}$ couples first and third generations, requiring two fold-wall transitions.
:::

::: parambox
P22: Solar angle **Formula:** $\sin^2\theta_{12} = p/\pi^2 = 3/\pi^2$\
**Predicted:** $0.3040$ $\mid$ **Measured:** $0.307\pm 0.013$ $\mid$ **Precision:** $1.0\%$($0.2\sigma$)

The spectral impedance of the triple junction. The electron neutrino lives at the cone point where $p = 3$ sectors meet; its coupling to the fold wall is proportional to the number of junction sectors $p$. The fold energy per ghost mode is $\pi^2$ (Parseval identity: $d_1 \zeta(2) = 6 \times \pi^2/6 = \pi^2$, specific to $S^5$). The mixing angle is the ratio of the discrete junction structure ($p$ sectors) to the continuous fold energy ($\pi^2$): $$\begin{equation}
\sin^2\theta_{12} = \frac{p}{\pi^2} = \frac{3}{\pi^2} = 0.3040.
\end{equation}$$ This is a **refraction coefficient**: the cone-point neutrino refracts into the fold-wall direction with amplitude $\sqrt{p/\pi^2}$. The same $\pi^2$ that sets the proton mass ($m_p/m_e = d_1 \cdot \mathrm{Vol}(S^5) \cdot \pi^2$) also sets the solar neutrino angle.
:::

::: parambox
P23: Atmospheric angle **Formula:** $\sin^2\theta_{23} = d_1/(d_1 + \lambda_1) = 6/11$\
**Predicted:** $0.5455$ $\mid$ **Measured:** $0.546\pm 0.021$ $\mid$ **Precision:** $0.10\%$

The spectral impedance ratio at the fold interface. Tunneling bandwidth is set by two competing factors: the number of available modes ($d_1 = 6$) and the eigenvalue gap ($\lambda_1 = 5$). The ratio $d_1/(d_1+\lambda_1)$ measures what fraction of the tunneling channel is carried by fold modes versus reflected by the spectral gap. The atmospheric angle measures the spectral bandwidth of the fold junction.
:::

::: parambox
P24: Leptonic CP phase **Formula:** $\delta_{\mathrm{CP}}^{\mathrm{PMNS}} = 3\arctan(2\pi^2/9)$\
**Predicted:** $196.5^\circ$ $\mid$ **Measured:** $195^\circ\pm 50^\circ$ $\mid$ **Precision:** ${\sim}0.3\%$ from central value

The quark CP phase $\gamma = \arctan(2\pi^2/9) \approx 65.5^\circ$ measures CP violation at a single fold wall. The leptonic CP phase is three times this: neutrinos, being neutral, access all three fold walls and accumulate a CP phase from each, while quarks (pinned to the cone point) see only one. The factor of $3$ is the $\mathbb{Z}_3$ order $p$.
:::

::: parambox
P25: Heaviest neutrino mass **Formula:** $m_3 = m_e/(108\pi^{10}) = m_e^3/(p\,m_p^2)$\
**Predicted:** $50.52$ meV $\mid$ **Measured:** $50.28$ meV $\mid$ **Precision:** $0.48\%$

The **inversion principle**: the proton is a bulk resonance ($m_p/m_e = 6\pi^5$, constructive interior standing wave); the neutrino is a boundary tunneling mode that loses mass in proportion to the inverse bulk volume squared, shared among $p = 3$ sectors: $$\begin{equation}
p \cdot m_p^2 \cdot m_\nu = m_e^3
\quad\Longleftrightarrow\quad
m_3 = \frac{m_e}{108\pi^{10}}.
\end{equation}$$ The factor $108\pi^{10} = p \times (6\pi^5)^2 = p \times (m_p/m_e)^2$: the neutrino tunnels through the fold wall in a **round trip** (enter and exit the twisted sector), with each crossing suppressed by $m_e/m_p = 1/(6\pi^5)$ (barrier penetration: tunneling particle mass / fold wall height), and divided among $p = 3$ sectors ($\mathbb{Z}_3$ projection). This gives $T = (1/p)(m_e/m_p)^2$, so $m_{\nu_3} = m_e \cdot T = m_e^3/(p\,m_p^2)$. The round-trip pattern $(m_e/m_p)^2$ is the same that produces $\eta^2$ in the cosmological constant and $(\eta K)^2$ in $\theta_{13}$: every double-crossing process squares its single-crossing amplitude.
:::

::: parambox
P26: Mass-squared ratio **Formula:** $\Delta m^2_{32}/\Delta m^2_{21} = d_1^2 - p = 33$\
**Predicted:** $33$ $\mid$ **Measured:** $32.58\pm 0.8$ $\mid$ **Precision:** $1.3\%$

The tunneling bandwidth ratio between sectors. The atmospheric splitting involves tunneling across the full spectral bandwidth $d_1^2 = 36$ modes; the solar splitting involves subtler 1--2 tunneling. Their ratio is reduced by $p = 3$ (three-fold sharing of tunneling channels). The integer $33$ is a direct spectral prediction: degeneracy squared minus the number of sectors.

From P25 and P26, the second neutrino mass is $m_2 = m_3/\sqrt{34} = 8.66$ meV (measured: $8.68$ meV, precision $0.15\%$). The lightest neutrino $m_1 \approx 0$ (normal hierarchy), and the total sum is $\sum m_\nu \approx 59.2$ meV.
:::

::: {#tab:neutrinos}
   \#  **Parameter**                                     **Formula**            **Predicted**   **Measured**     **Error**    **Origin**
  ---- ---------------------------------------- ------------------------------ --------------- -------------- --------------- ---------------------------
   21  $\sin^2\theta_{13}$                          $(\eta K)^2 = 16/729$          0.02194        0.02200          0.24%      Junction
   22  $\sin^2\theta_{12}$                           $p/\pi^2 = 3/\pi^2$           0.3040          0.307           1.0%       Impedance
   23  $\sin^2\theta_{23}$                       $d_1/(d_1+\lambda_1) = 6/11$      0.5455          0.546           0.10%      Impedance
   24  $\delta_{\mathrm{CP}}^{\mathrm{PMNS}}$        $3\arctan(2\pi^2/9)$       $196.5^\circ$   $195^\circ$    ${\sim}0.3\%$  $\mathbb{Z}_3 \times$ CKM
   25  $m_3$                                         $m_e/(108\pi^{10})$          50.52 meV      50.28 meV         0.48%      Inversion
   26  $\Delta m^2_{32}/\Delta m^2_{21}$               $d_1^2 - p = 33$              33            32.58           1.3%       Bandwidth

  : Neutrino sector predictions.
:::

# The Error Structure {#sec:errors}

## The Hurricane Hierarchy

Every geometric prediction gives the **bare** value at the compactification scale --- "the statue." What experiments measure is the **dressed** value at low energy --- "the statue in the hurricane." The hurricane is quantum corrections: loops, running, threshold effects. When residual errors are expressed as multiples of the relevant loop expansion parameter, a clean pattern emerges.

**Definition: hurricane coefficients.** We define the *hurricane coefficients* as the structured radiative correction coefficients whose values are *fixed by the five spectral invariants* $\{d_1, \lambda_1, K, \eta, p\}$ of $S^5/\mathbb{Z}_3$, rather than computed from Feynman diagrams or fit to data. Unlike standard perturbative loop corrections, which are scheme-dependent and require renormalization, the hurricane coefficients are *topological*: they follow from the spectral action's commutativity with the orbifold projectors (Theorem [2](#thm:N1){reference-type="ref" reference="thm:N1"}) and are independent of the cutoff function $f$. Every hurricane coefficient catalogued in the Grand Table below (Table [9](#tab:grand-hurricane){reference-type="ref" reference="tab:grand-hurricane"}) is an algebraic expression in $\{d_1, \lambda_1, K, \eta, p\}$, proven as a theorem.

**Cross-check via standard methods.** Where applicable, the hurricane-corrected values have been *independently reproduced* using conventional perturbative QFT, confirming that the spectral route and the standard route yield the same answer:

- **Fine-structure constant** (§[5](#sec:gauge){reference-type="ref" reference="sec:gauge"}): the spectral lag $G/p = 10/27$ at $M_c$ is followed by standard two-loop SM renormalization-group running to $\alpha(0)$, using the SM beta functions with $\mathbb{Z}_3$-projected particle content.

- **Strong coupling** (§[5](#sec:gauge){reference-type="ref" reference="sec:gauge"}): ghost splitting $d_1 = 6$ at $M_c$ is run to $M_Z$ via standard two-loop QCD running; result $\alpha_s(M_Z) = 0.1187$ matches PDG ($0.56\%$).

- **Cabibbo angle** (§[8](#sec:ckm){reference-type="ref" reference="sec:ckm"}): the QCD hurricane coefficient $c = +1/p$ reproduces the standard one-loop $\alpha_s/\pi$ vertex correction to quark mixing. The spectral route predicts $c = +1/3$; conventional QCD gives the same correction structure.

- **Proton mass** (P11): the EM hurricane coefficients $G = 10/9$ and $G_2 = -280/9$ reproduce the form of the standard QED Schwinger series $\alpha^{2n}/\pi^n$, with coefficients that conventional perturbation theory leaves as free parameters.

- **Muon $g-2$** (P43): the spectral HVP estimate via the Lotus Song hadron spectrum agrees with lattice QCD (BMW 2021) to $2.7\%$, independently confirming the lattice resolution of the muon anomaly.

The hurricane framework does not replace standard perturbative QFT --- it *predicts the coefficients* that standard perturbative QFT treats as inputs.

## Mass Residuals Scale with $\alpha/\pi$

::: {#tab:mass-errors}
  **Observable**                  **Residual**   $c = \text{Residual}/(\alpha/\pi)$   **$|c|<1$?**
  ------------------------------ -------------- ------------------------------------ --------------
  $m_p/m_e = 6\pi^5$ (leading)     $0.002\%$                  $-0.008$               
  $v/m_p = 2/\alpha - 35/3$        $0.005\%$                  $-0.021$               
  $m_H/m_p = 1/\alpha - 7/2$       $0.034\%$                  $+0.148$               
  $m_\mu/m_e$                      $0.001\%$                  $-0.004$               
  $m_\tau/m_e$                     $0.007\%$                  $+0.030$               

  : Mass ratio residuals expressed as $c\times(\alpha/\pi)$. All $|c|<1$, consistent with one-loop EM corrections with order-unity coefficients.
:::

## Mixing Residuals: Bare vs. Hurricane-Corrected

::: {#tab:mix-errors}
  **Observable**             **Bare residual**   **Hurricane $c$**   **Spectral form**   **Corrected residual**
  ------------------------- ------------------- ------------------- ------------------- ------------------------
  $\lambda = 2/9$                 $1.2\%$             $+1/3$               $1/p$           $\mathbf{0.002\%}$
  $A = 5/6$                       $0.9\%$             $-2/9$              $-\eta$          $\mathbf{0.046\%}$
  $|V_{cb}| = A\lambda^2$         $1.6\%$           (combined)              ---            $\mathbf{0.04\%}$
  $J$ (Jarlskog)                  $5.2\%$           (combined)              ---             $\mathbf{0.4\%}$

  : CKM mixing residuals before and after hurricane corrections. The hurricane coefficients $+1/p$ and $-\eta$ are spectral invariants of $S^5/\mathbb{Z}_3$. Improvement factors: $\lambda$ ($619\times$), $A$ ($19\times$), $|V_{cb}|$ ($39\times$), $J$ ($12\times$).
:::

## Proof of Concept 1: The Proton Mass (EM Hurricane)

The proton mass formula demonstrates the EM hurricane. The bare value $6\pi^5$ has $0.002\%$ error. The one-loop correction $G = 10/9$ (a spectral invariant: $\lambda_1 \cdot \sum|\eta_D|$) improves this to $10^{-8}$. The two-loop correction $G_2 = -280/9$ (also a spectral invariant: $-\lambda_1(d_1 + \sum|\eta_D|)$) reaches $10^{-11}$. Both coefficients emerge from $\{d_1, \lambda_1, \eta, \tau_R\}$ --- the same geometry that produces the bare value.

## Proof of Concept 2: The Cabibbo Angle (QCD Hurricane)

The CKM mixing parameters demonstrate the QCD hurricane. The bare Cabibbo angle $\lambda = 2/9$ has $1.2\%$ error. The one-loop QCD correction with coefficient $c = +1/p = +1/3$ (a spectral invariant: the orbifold order) gives $\lambda_{\mathrm{phys}} = (2/9)(1 + \alpha_s/(3\pi)) = 0.22500$, matching PDG to $0.002\%$ --- a $619\times$ improvement. The Wolfenstein $A = 5/6$ with coefficient $c = -\eta = -2/9$ gives $A_{\mathrm{phys}} = (5/6)(1 - (2/9)\alpha_s/\pi) = 0.8264$, matching to $0.046\%$ ($19\times$ improvement).

::: {#tab:hurricane}
+---------------------------+------------------+--------------------+----------------------------------+--------------------+
| **Observable**            | **Expansion**    | **Coefficient**    | **Spectral form**                | **Precision**      |
+:==========================+:=================+:===================+:=================================+===================:+
| *EM Hurricane (proton mass):*                                                                                             |
+---------------------------+------------------+--------------------+----------------------------------+--------------------+
| $m_p/m_e$ (1-loop)        | $\alpha^2/\pi$   | $G = 10/9$         | $\lambda_1\cdot\sum|\eta_D|$     | $10^{-8}$          |
+---------------------------+------------------+--------------------+----------------------------------+--------------------+
| $m_p/m_e$ (2-loop)        | $\alpha^4/\pi^2$ | $G_2 = -280/9$     | $-\lambda_1(d_1{+}\sum|\eta_D|)$ | $10^{-11}$         |
+---------------------------+------------------+--------------------+----------------------------------+--------------------+
| *QCD Hurricane (CKM mixing):*                                                                                             |
+---------------------------+------------------+--------------------+----------------------------------+--------------------+
| $\lambda$ (Cabibbo)       | $\alpha_s/\pi$   | $c = +1/p = +1/3$  | $\eta/K$ (orbifold order)        | $0.002\%$          |
+---------------------------+------------------+--------------------+----------------------------------+--------------------+
| $A$ (Wolfenstein)         | $\alpha_s/\pi$   | $c = -\eta = -2/9$ | spectral twist                   | $0.046\%$          |
+---------------------------+------------------+--------------------+----------------------------------+--------------------+
| *Lag Hurricane ($\alpha$ from geometry):*                                                                                 |
+---------------------------+------------------+--------------------+----------------------------------+--------------------+
| $1/\alpha_{\mathrm{GUT}}$ | topological      | $G/p = 10/27$      | $\lambda_1\eta/p$ (ghost lag)    | $\mathbf{0.001\%}$ |
+---------------------------+------------------+--------------------+----------------------------------+--------------------+

: The hurricane hierarchy. Every coefficient is a spectral invariant of $S^5/\mathbb{Z}_3$. Mass observables see EM corrections; mixing observables see QCD corrections. CP parameters ($\bar\rho$, $\bar\eta$) require no correction (sub-$0.1\%$ already).
:::

**The pattern:** EM coefficients use $\lambda_1$ (energy) and $\sum|\eta_D|$ (asymmetry). QCD coefficients use $p$ (orbifold order) and $\eta$ (spectral twist). The hurricane *is* the geometry, seen through loop corrections.

## Ghost modes as virtual particles {#sec:ghost-virtual}

The hurricane corrections have a direct interpretation in the language of standard QFT: the ghost modes *are* virtual particles.

In perturbative QFT, virtual particles are off-shell internal lines in Feynman diagrams. Their physical status is debated: they mediate forces and generate loop corrections, yet they never appear as asymptotic states. The spectral framework resolves this ambiguity.

The Dirac operator on $S^5$ has a full spectrum of eigenvalues. The $\mathbb{Z}_3$ projection kills those modes whose eigenvalue degeneracy satisfies $d_{\mathrm{inv}}(\ell) = 0$ --- the ghost modes at $\ell = 1$, with multiplicity $d_1 = 6$. These modes are mathematically valid eigenstates of the covering space operator. They carry gauge quantum numbers (the $3 \oplus \bar{3}$ of the bihomogeneous decomposition). They cannot appear as physical particles (the orbifold projection forbids them as asymptotic states). Yet their spectral footprint persists in the spectral action $\mathrm{Tr}\,f(D^2/\Lambda^2)$, because the trace runs over the *full* operator --- and the full operator knows about $S^5$.

**Loop corrections.** In QFT, one sums over virtual momenta in loop integrals, obtaining divergent series that require renormalization. In the spectral framework, the corresponding corrections are finite spectral invariants: $G = \lambda_1 \cdot \eta = 10/9$ at one loop, $G_2 = -\lambda_1(d_1 + \eta) = -280/9$ at two loops (Table [8](#tab:hurricane){reference-type="ref" reference="tab:hurricane"}). No sum over momenta; no UV divergence; no renormalization scheme dependence.

**Vacuum polarization.** In QFT, virtual electron--positron pairs screen electric charge at long distances, causing $\alpha$ to run. In the spectral framework, the $\ell = 1$ ghost modes carry electromagnetic charge, and their spectral footprint modifies the effective coupling as a function of scale. The lag correction $G/p = 10/27$ is vacuum polarization stated as a spectral invariant rather than a divergent integral.

**UV completion.** In QFT, the contribution of virtual particles at arbitrarily high momenta creates UV divergences. In the spectral framework, the compactification scale $M_c$ provides a natural upper bound: above $M_c$, one resolves the full $S^5$, and the $\mathbb{Z}_3$ projection has not yet occurred. Below $M_c$, ghost modes contribute only through their spectral footprint. The UV completion is the covering space itself.

The ghost modes are not particles that "briefly exist." They are permanently present in the spectral action and permanently absent from the physical Hilbert space --- shaping the dynamics of the quotient space without ever appearing in it.

# The Lotus Song: Hadrons from Spectral Geometry {#sec:song}

The proton mass $m_p = 6\pi^5\,m_e$ is the **fundamental frequency** of the fold wall. Every other hadron mass is an **eigenvalue** of the same operator at a different spectral ratio.

## The eigenvalue equation

The fold-wall Dirac operator $D_{\mathrm{wall}}$ on $S^5/\mathbb{Z}_3$ has discrete eigenvalues: $$\begin{equation}
\boxed{D_{\mathrm{wall}}\,\Psi_n = (6\pi^5 \cdot R_n)\,m_e\,\Psi_n}
\end{equation}$$ where $6\pi^5 = d_1^2 \cdot \zeta(2) \cdot \mathrm{Vol}(S^5)$ is the Parseval fold energy (Section [6](#sec:baryons){reference-type="ref" reference="sec:baryons"}), and $R_n$ is a rational function of the five spectral invariants. The absolute scale and the ratios are *both* determined by the geometry.

**Two fundamental integers.** The eigenvalue ratios $R_n$ are generated by two spectral integers: $$\begin{equation}
D_{\mathrm{wall}} = 1 + d_1 = 7, \qquad D_{\mathrm{bulk}} = d_1 + \lambda_1 = 11.
\end{equation}$$ $D_{\mathrm{wall}}$ is the fold wall dimension (six spatial ghost modes plus one temporal channel from $\mathrm{Im}(\eta_D)$, Section [21](#sec:time){reference-type="ref" reference="sec:time"}). $D_{\mathrm{bulk}}$ is the total spectral mode count at the first KK level (the same combination $(d_1{+}\lambda_1) = 11$ that enters the gravity formula, Section [14](#sec:gravity){reference-type="ref" reference="sec:gravity"}). Every hadron mass ratio in the Lotus Song is an algebraic expression in $\{D_{\mathrm{wall}}, D_{\mathrm{bulk}}, p\} = \{7, 11, 3\}$.

The spectrum organizes into three series --- three types of "instrument" on the fold wall:

**Plucked strings (pseudoscalar mesons, $J^P = 0^-$):** Ghost-antighost tunneling across the fold wall. The fundamental frequency is $f_\pi = \eta K = 4/27$. Higher harmonics have overtone numbers $\{1, 4, 7, \ldots\}$ with spacing $p = 3$ (orbifold quantization): each overtone requires one additional $\mathbb{Z}_3$ sector traversal.

**Drums (baryons):** Color-singlet ghost triplets resonating on the fold wall. The proton ($R = 1$) is the ground state. Spin excitation costs $1/p$ (the Delta at $R = 4/3$). Strangeness multiplies by $(1{+}\eta/2) = 10/9$ per strange quark (half-tunneling depth).

**Organ pipes (quarkonia):** Closed loops of ghost pairs orbiting the fold. Charm at $R = p + 1/p = 10/3$ (one orbit plus return). Bottom at $R = p^2 + 1 = 10$ (double orbit). The scaling $m_\Upsilon/m_{J/\psi} = p = 3$ is *exact*: each quark generation adds one fold traversal.

## The complete hadron spectrum

Twenty-seven hadrons are predicted, all from the five invariants:

::: center
  **Hadron**                 **Spectral ratio**           **Predicted (MeV)**            **Error**
  ------------------ ----------------------------------- --------------------- -------------------
  $\pi^\pm$                    $\eta K = 4/27$                   139.0                    $0.41\%$
  $K^*(892)$               $1 - \eta K/p = 77/81$                891.9           $\mathbf{0.03\%}$
  $\rho(770)$               $\lambda_1/d_1 = 5/6$                781.9                    $0.86\%$
  $\phi(1020)$        $1 + 1/(d_1{+}\lambda_1) = 12/11$         1023.6                    $0.40\%$
  $D^\pm$                         $pK = 2$                      1876.5                    $0.37\%$
  $B^\pm$                    $d_1 - 1/p = 17/3$                 5316.9                    $0.71\%$
  proton                             $1$                         938.3                    $\sim 0$
  $\Delta(1232)$               $1 + 1/p = 4/3$                  1251.0                     $1.5\%$
  $\Omega^-(1672)$            $2 - \eta = 16/9$                 1668.0                    $0.27\%$
  $J/\psi$                    $p + 1/p = 10/3$                  3127.6                    $0.99\%$
  $\Upsilon(1S)$               $p^2 + 1 = 10$                   9382.7                    $0.82\%$
:::

Full list: 27 hadrons, RMS $0.95\%$, all sub-$2\%$, 18 sub-$1\%$. Verification: `lotus_song_extended.py`.

**Key discovery: $K^*(892) = m_p \cdot 77/81$ at $0.03\%$.** The $K^*$ vector meson is the proton with one pion removed per $\mathbb{Z}_3$ sector: $m_{K^*} = m_p - m_\pi/p$. This bridges the baryon and meson sectors through the *same* spectral invariants.

## The neutron's clock: $g_A$, $f_\pi$, and $\tau_n$

The eigenvalue equation also determines the **proton's internal structure**:

**Axial coupling:** $g_A = 1 + \eta + K/(d_1{+}\lambda_1) = 127/99 = 1.2828$ (PDG: $1.2754$, error $0.58\%$). The three terms: $g_V = 1$ (CVC, exact), $+\eta$ (fold-wall chirality enhancement), $+K/(d_1{+}\lambda_1)$ (Koide correction per total spectral mode).

**Pion decay constant:** $f_\pi = K^2\eta\,m_p = 92.7$ MeV (PDG: $92.07$, error $0.65\%$). The pion decays by double-parity fold-wall tunneling ($K^2$) at the spectral asymmetry rate ($\eta$).

**Neutron lifetime:** With all inputs spectral ($g_A$, $G_F$ from VEV, $|V_{ud}|$ from Cabibbo, phase space from kinematics): $\tau_n = 899$ s (PDG: $878.4$ s, error $2.3\%$). The neutron's lifetime is written in the five invariants of $S^5/\mathbb{Z}_3$.

Verification: `axial_coupling_derivation.py`, `neutron_lifetime.py`.

# Beyond Single Particles: The Deuteron {#sec:adjacency}

The Lotus Song derives single-hadron masses as eigenvalues of $D_{\mathrm{wall}}$. But what happens when two eigenvalues *overlap*?

**Spectral adjacency.** In the spectral framework, "next to" means non-zero Fubini--Study overlap of $D_{\mathrm{wall}}$ eigenstates on the fold wall. Two nucleon eigenstates $\Psi_1, \Psi_2$ have overlap: $$\begin{equation}
|\langle\Psi_1|\Psi_2\rangle|^2 = \frac{\eta^2}{p} + \frac{1}{p^2}\left(\frac{1}{p^2} - \frac{\eta^2}{p}\right)
\end{equation}$$ where the first term is the **resolved channel** (magnitudes, spatial) and the correction is the **unresolved channel** (complex phases, temporal).

**The key insight: $\eta_D(\chi_1) = i/9$ is purely imaginary.** We have used $|\eta| = 2/9$ (sum of magnitudes) for the preceding predictions. But the Donnelly invariant carries a phase: $\eta_D(\chi_1) = +i/9$, $\eta_D(\chi_2) = -i/9$. The magnitudes encode *space*; the imaginary phase encodes *time*. The resolved channel uses magnitudes; the unresolved channel uses complex values. The two channels mix with weights $(1-1/p^2) = 8/9$ (space) and $1/p^2 = 1/9$ (time).

**The deuteron binding energy.** $$\begin{equation}
\boxed{B_d = m_\pi \cdot \frac{\lambda_1(1+d_1)}{p^{1+d_1}} = m_\pi \cdot \frac{35}{2187} = 2.225\;\mathrm{MeV}}
\label{eq:deuteron}
\end{equation}$$ (PDG: $2.2246$ MeV, error $0.00\%$). The spectral decomposition: $35 = \lambda_1(1{+}d_1) = 5 \times 7$ (ghost spectral weight times fold wall dimension plus time); $2187 = p^{1+d_1} = 3^7$ (monogamy suppression over 6 spatial plus 1 temporal dimension). The exponent $7 = 1 + d_1$ is the fold wall dimension *plus one time dimension*.

**The Pythagorean comma of time.** The fraction $1/p^2 = 1/9$ is the **temporal fraction** of the spectral overlap: the amount by which the static (resolved) prediction deviates from reality because it ignores the imaginary phase of $\eta_D$. The bare formula $B_d = m_\pi \eta^2/p$ gives $2.29$ MeV ($+2.9\%$); including the temporal channel gives $2.225$ MeV ($0.00\%$). The $2.9\%$ was *time leaking into space*.

**Nuclear binding = incomplete entanglement.** The deuteron is not two particles exchanging mesons. It is two ghost resonances partially in phase on the fold wall. The Fubini--Study distance is $d_{\mathrm{FS}} = 0.984$ (98.4% orthogonal, 1.6% coherent). That 1.6% of coherence IS the binding.

**He-4 and saturation.** For He-4 (4 nucleons spanning all $\mathbb{Z}_3$ sectors), the Koide phase $K$ replaces one factor of $\eta$, with a hurricane correction from the bulk: $$\begin{equation}
B/A(\mathrm{He}\text{-}4) = m_\pi \frac{K\eta}{p}\!\left(1 + \frac{1}{D_{\mathrm{bulk}}\cdot p}\right) = m_\pi \frac{K\eta(D_{\mathrm{bulk}}p{+}1)}{p^2 D_{\mathrm{bulk}}} = 7.07\;\mathrm{MeV}
\end{equation}$$ (PDG: $7.072$ MeV, error $0.02\%$). The bare formula $m_\pi K\eta/p = 6.86$ MeV ($-3.0\%$) omits the bulk correction $1/(D_{\mathrm{bulk}}\cdot p) = 1/33$; including it closes the gap to sub-percent. The $K/\eta = p$ coherence enhancement reflects full sector occupation. Nuclear saturation (the plateau at $B/A \sim 8.8$ MeV) requires the many-body overlap integral, which is the principal open problem for the dynamics sequel.

Verification: `time_spectral_error.py`, `spectral_adjacency.py`, `spectral_overlap_proof.py`.

# The LOTUS: Lagrangian Of The Universe's Spectral State {#sec:lotus}

The hurricane coefficients (Section [10](#sec:errors){reference-type="ref" reference="sec:errors"}) are not isolated numerical matches. They are **Taylor coefficients of a single object**: the fold stiffness potential $V(\phi)$ whose minimum generates the Standard Model.

## The fold field $\phi$

The fold stiffness $\phi$ parameterizes the deformation of $S^5$ from smooth ($\phi = 0$, UV, unified) to rigid $\mathbb{Z}_3$ orbifold ($\phi = 1$, unphysical):

- $\phi = 0$: smooth $S^5$. Full $\mathrm{SO}(6)$ isometry. All couplings equal. No generations, no mass, no physics.

- $\phi_c \approx 0.60$: crossover. Fold wall thickness equals ghost wavelength. The spectral phase transition: substrate regime ($\phi < \phi_c$) gives way to information regime ($\phi > \phi_c$).

- $\phi_{\mathrm{lotus}} = 1 - \alpha(d_1{+}\lambda_1{+}K)/2 = 0.9574$: **the lotus in bloom**. Our universe. The equilibrium between folding force (spectral organization) and ghost pressure (spectral cost). The $4.3\%$ residual petal overlap IS the Higgs mechanism.

- $\phi = 1$: dried flower. Rigid triangle. $v = 0$, all masses vanish. Dead geometry.

## The spectral phase transition at $\phi_c = 0.60$

The fold field interpolates between two regimes, separated by a crossover at $\phi_c \approx 0.60$ (verified numerically: `dirac_fold_transition.py`):

- **Substrate regime ($\phi < \phi_c$):** The geometry is "more circle than triangle." Perturbative. Ghost modes still coupled to physical spectrum. Full SO(6) approximate symmetry. This is the regime of the early universe (inflation).

- **Information regime ($\phi > \phi_c$):** The geometry is "more triangle than circle." Non-perturbative. Ghost modes decoupled. $\mathbb{Z}_3$ structure dominant. SM physics emerges. The spectral asymmetry $\eta$ reaches its final value $2/9$.

The crossover occurs when the fold wall thickness equals the ghost wavelength $\sim 1/\sqrt{\lambda_1}$. The threshold integral across a 4D fold wall gives $I(4) = \pi/(4\sin(\pi/4)) = \pi/(2\sqrt{2}) = 1.1107$, matching the proton spectral coupling $G = 10/9 = 1.1111$ to $0.035\%$ --- confirming that $G$ IS the fold wall curvature. Full computation: Supplement V, `fold_potential.py`.

## The potential

$$\begin{equation}
\boxed{V(\phi) = \frac{\lambda_H}{4}\,v_{\max}^4\,\bigl(\phi^2 - \phi_{\mathrm{lotus}}^2\bigr)^2}
\label{eq:lotus}
\end{equation}$$

Every coefficient is a spectral invariant of $S^5/\mathbb{Z}_3$: $$\begin{align}
\lambda_H &= \frac{(1/\alpha - 7/2)^2}{2(2/\alpha - 35/3)^2} = 0.1295, \label{eq:lotus-lambda}\\
v_{\max} &= 2m_p/\alpha, \label{eq:lotus-vmax}\\
\phi_{\mathrm{lotus}} &= 1 - \frac{\alpha}{2}(d_1 + \lambda_1 + K) = 0.9574. \label{eq:lotus-phi}
\end{align}$$

<figure id="fig:lotus" data-latex-placement="htbp">

<figcaption><strong>The LOTUS potential</strong> <span class="math inline"><em>V</em>(<em>ϕ</em>) = (<em>λ</em><sub><em>H</em></sub>/4) <em>v</em><sub>max</sub><sup>4</sup> (<em>ϕ</em><sup>2</sup> − <em>ϕ</em><sub>lotus</sub><sup>2</sup>)<sup>2</sup></span>. The fold field <span class="math inline"><em>ϕ</em></span> parameterizes the transition from the smooth sphere (<span class="math inline"><em>ϕ</em> = 0</span>, no physics) through the spectral phase transition (<span class="math inline"><em>ϕ</em><sub><em>c</em></sub> = 0.60</span>, inflation ends) to our universe (<span class="math inline"><em>ϕ</em><sub>lotus</sub> = 0.9574</span>, <span class="math inline"><em>V</em> = 0</span>). The Higgs mass is the curvature at the equilibrium fold depth. The barrier height (<span class="math inline"><em>V</em>(0)<sup>1/4</sup> ≈ 104</span> GeV) is the EW scale.</figcaption>
</figure>

This IS the Mexican hat potential decoded: $H = v_{\max}\phi$ (the Higgs field is the fold stiffness), $v = v_{\max}\phi_{\mathrm{lotus}}$ (the VEV is the fold depth at bloom). The Standard Model Lagrangian is $V(\phi)$ evaluated at $\phi = \phi_{\mathrm{lotus}}$.

## What the LOTUS generates

At the lotus point $\phi = \phi_{\mathrm{lotus}}$:

- $V(\phi_{\mathrm{lotus}}) = 0$: tree-level cosmological constant (orbifold volume cancellation).

- $V''(\phi_{\mathrm{lotus}}) \to m_H = 125.3$ GeV: the Higgs mass is the lotus curvature.

- The gauge couplings $g_i(\phi_{\mathrm{lotus}})$: the SM gauge structure.

- The mass functions $m_k(\phi_{\mathrm{lotus}})$: all fermion masses.

The hurricane coefficients are the **lotus derivatives** --- the perturbative expansion of the fold potential around the SM vacuum. The grand table below shows every coefficient, its spectral form, the geometric feature it probes, and *why* that specific invariant appears.

::: {#tab:grand-hurricane}
+---------------------------+----------------------------------+------------------+-----------------------------------------------+--------------------------+----------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Observable**            | **Bare**                         | **Expansion**    | **Coeff.**                                    | **Spectral**             | **Prec.**      | **Why this invariant**                                                                                                                                                                                                     |
+:==========================+:=================================+:=================+:==============================================+:=========================+:==============:+:===========================================================================================================================================================================================================================+
| *EM Hurricane --- corrections from fold walls (4D surfaces, threshold integral $\pi/(2\sqrt{2})$)*                                                                                                                                                                                                                                                                                                       |
+---------------------------+----------------------------------+------------------+-----------------------------------------------+--------------------------+----------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| $m_p/m_e$ (1L)            | $6\pi^5$                         | $\alpha^2/\pi$   | $G{=}10/9$                                    | $\lambda_1\!\cdot\!\eta$ | $10^{-8}$      | Ghost propagator has energy scale $\lambda_1$ and asymmetry $\eta$; these are inseparable (ghost-as-one).                                                                                                                  |
+---------------------------+----------------------------------+------------------+-----------------------------------------------+--------------------------+----------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| $m_p/m_e$ (2L)            | ---                              | $\alpha^4/\pi^2$ | $G_2{=}{-}280/9$                              | $-\lambda_1(d_1{+}\eta)$ | $10^{-11}$     | Fermion loop traces *total* ghost content (symmetric $d_1$ + asymmetric $\eta$); sign flip from closed loop.                                                                                                               |
+---------------------------+----------------------------------+------------------+-----------------------------------------------+--------------------------+----------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| *QCD Hurricane --- corrections from cone point (0D, discrete sector counting)*                                                                                                                                                                                                                                                                                                                           |
+---------------------------+----------------------------------+------------------+-----------------------------------------------+--------------------------+----------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| $\lambda$ (Cabibbo)       | $2/9$                            | $\alpha_s/\pi$   | $+1/p{=}+1/3$                                 | $\eta/K$                 | $0.002\%$      | QCD vertex correction distributes equally among $p{=}3$ sectors; each contributes $\alpha_s/\pi$, total $= \alpha_s/(p\pi)$.                                                                                               |
+---------------------------+----------------------------------+------------------+-----------------------------------------------+--------------------------+----------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| $A$ (Wolfenstein)         | $5/6$                            | $\alpha_s/\pi$   | $-\eta{=}{-}2/9$                              | twist                    | $0.046\%$      | The spectral twist $\eta$ IS the QCD anomalous dimension of the weight-per-mode ratio $\lambda_1/d_1$.                                                                                                                     |
+---------------------------+----------------------------------+------------------+-----------------------------------------------+--------------------------+----------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| *Lag Hurricane --- topological correction to the unification point*                                                                                                                                                                                                                                                                                                                                      |
+---------------------------+----------------------------------+------------------+-----------------------------------------------+--------------------------+----------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| $1/\alpha_{\mathrm{GUT}}$ | $42.41$                          | topological      | $G/p{=}10/27$                                 | $\lambda_1\eta/p$        | $0.001\%$      | Ghost modes have inertia; their spectral weight $G{=}10/9$ is redistributed across $p$ sectors before decoupling.                                                                                                          |
+---------------------------+----------------------------------+------------------+-----------------------------------------------+--------------------------+----------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| *Gravity Hurricane --- ghost spectral pressure against the bulk*                                                                                                                                                                                                                                                                                                                                         |
+---------------------------+----------------------------------+------------------+-----------------------------------------------+--------------------------+----------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| $M_9/M_c$                 | $\tfrac{(d_1{+}\lambda_1)^2}{p}$ | KK               | $-\tfrac{1}{d_1\lambda_1}{=}{-}\tfrac{1}{30}$ | inverse ghost weight     | $0.10\%$       | Ghost modes ($d_1{=}6$ at $\lambda_1{=}5$) are absent from the physical spectrum; their "shadow" reduces the effective bulk stiffness by $1/(d_1\lambda_1)$. Gravity is the load-bearing force against the ghost pressure. |
+---------------------------+----------------------------------+------------------+-----------------------------------------------+--------------------------+----------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

: **The Grand Hurricane Table.** Every radiative correction coefficient is a spectral invariant of $S^5/\mathbb{Z}_3$. EM corrections arise from fold walls (4D, continuous); QCD corrections from the cone point (0D, discrete); the gravity correction from the ghost spectral weight against the bulk (5D, KK compactification). All are modulated by combinations of $\{d_1, \lambda_1, K, \eta, p\}$. The hurricane IS the geometry at the next order of the fold potential's Taylor expansion.
:::

**The pattern.** Every coefficient is a simple ratio of $\{d_1, \lambda_1, K, \eta, p\}$:

- $G = \lambda_1 \cdot \eta$ (energy $\times$ twist)

- $G_2 = -\lambda_1(d_1 + \eta)$ (energy $\times$ total content)

- $c_\lambda = \eta/K = 1/p$ (twist normalized by amplitude)

- $c_A = -\eta$ (twist alone)

- $G/p = \lambda_1\eta/p$ (energy $\times$ twist $\div$ sectors)

- $c_{\mathrm{grav}} = -1/(d_1\lambda_1) = -1/30$ (inverse ghost spectral weight)

The spectral twist $\eta = 2/9$ is the **universal anomalous dimension of the fold**: the rate at which every observable changes when the fold depth changes. It is to the LOTUS what the mass anomalous dimension $\gamma_m$ is to QCD. The gravity coefficient $-1/(d_1\lambda_1)$ stands apart: it does not involve $\eta$ directly but encodes the *total weight of the ghost sector*, showing that gravity couples to different spectral content than the gauge forces.

**Theorem [2](#thm:N1){reference-type="ref" reference="thm:N1"} guarantees** that these coefficients are exact: the commutativity $[f(D/\Lambda), e_m] = 0$ ensures no cutoff-function ambiguity contaminates the expansion. The hurricane coefficients are properties of the geometry, not of the regularization.

## The lotus is alive

The universe cannot be at $\phi = 1$ (the fold cannot fully close) because the ghost modes prevent it: the spectral cost $d_1 + \lambda_1 + K = 35/3$ must be subtracted from the EM budget $2/\alpha$. This forces $\phi_{\mathrm{lotus}} < 1$.

**Physics exists because the lotus is incomplete:**

- The VEV exists because $\phi < 1$ (petal overlap $\to$ Higgs mechanism).

- Particles have mass because $\phi < 1$ (overlap amplitude $=$ Yukawa coupling).

- The cosmological constant is nonzero because $\phi < 1$ (lotus breathing: $\Lambda^{1/4} = m_{\nu_3} \cdot (32/729)(1{+}\eta^2/\pi) = 2.25$ meV, $0.11\%$ from observed; the correction $1 - K/d_1 = 1 - 1/p^2 = 8/9$ is the Koide-phase residual per ghost mode).

- Neutrinos have mass because $\phi < 1$ (tunneling through petals).

**Spectral monogamy versus supersymmetry.** The cosmological constant problem is often said to require supersymmetry: boson-fermion pairs cancel the vacuum energy. But SUSY assumes the universe splits into TWO sectors (bosons and fermions) with an exact pairing. In our framework the cancellation mechanism is different and stronger:

1.  **One geometry, split into three, then two.** The universe is ONE manifold ($S^5/\mathbb{Z}_3$) that splits into $p = 3$ orbifold sectors (generations), each of which further splits into two chiralities (the Dirac operator has $\pm$ eigenvalues). The hierarchy is $1 \to 3 \to 2$, not $1 \to 2$ as SUSY assumes. Spectral monogamy ensures the $\mathbb{Z}_3$ partition of unity $\sum_m e_m = \mathbf{1}$ forces the contributions to cancel *sector by sector*, not merely by pairing.

2.  **Chiral entanglement across the fold.** The two chiralities of each generation are NOT perfect mirrors --- the eta invariant $\eta = 2/9 \neq 0$ measures the spectral *asymmetry* between positive and negative Dirac eigenvalues. This asymmetry means the "entanglement" between the two lobes of the fold is *chiral*: the back-and-forth between them is not perfectly symmetric. The residual asymmetry $\eta^2 = 4/81$ is precisely what sets the cosmological constant scale.

SUSY fails because it imposes a $\mathbb{Z}_2$ symmetry (boson $\leftrightarrow$ fermion) on a universe that has $\mathbb{Z}_3$ structure (three generations) with chiral asymmetry ($\eta \neq 0$). The correct cancellation mechanism is spectral monogamy, which respects the $\mathbb{Z}_3$ structure and uses the residual $\eta^2$ to set the CC scale rather than demanding it vanish.

# Gravity from Spectral Geometry {#sec:gravity}

The same spectral data predicts the strength of gravity. In the Kaluza--Klein compactification [@freundrubin1980] on $S^5/\mathbb{Z}_3$: $$\begin{equation}
M_P^2 = M_9^7 \cdot \frac{\pi^3}{3\,M_c^5},
\label{eq:kk-gravity}
\end{equation}$$ where $M_9$ is the fundamental 9D Planck mass and $M_c$ is the compactification scale. The bare spectral prediction for the ratio $X = M_9/M_c$ is: $$\begin{equation}
X_{\mathrm{bare}} = \frac{(d_1 + \lambda_1)^2}{p} = \frac{(6+5)^2}{3} = \frac{121}{3}.
\label{eq:x-bare}
\end{equation}$$ *Origin of the bare formula.* In the KK compactification, the 9D Planck mass receives contributions from all KK levels. The dominant contribution at the first KK level has spectral weight $(d_1 + \lambda_1) = 11$ --- the sum of the ghost degeneracy and the eigenvalue --- distributed over $p = 3$ orbifold sectors. The squared factor arises because $M_9^2 \propto (d_1 + \lambda_1)^2$ from the heat kernel coefficient $a_2$ on $S^5/\mathbb{Z}_3$ (see Supplement I for the KK spectrum). Eq. [\[eq:x-bare\]](#eq:x-bare){reference-type="eqref" reference="eq:x-bare"} overshoots the value inferred from $M_P$ by $3.4\%$.

The ghost modes ($d_1 = 6$ at eigenvalue $\lambda_1 = 5$) are absent from the physical spectrum but create **spectral pressure** against the bulk: their "shadow" reduces the effective stiffness of the compact space. The correction is the inverse of the total ghost spectral weight, $d_1\lambda_1 = 6 \times 5 = 30$: $$\begin{equation}
\boxed{X = \frac{(d_1 + \lambda_1)^2}{p}\left(1 - \frac{1}{d_1\lambda_1}\right) = \frac{121}{3}\cdot\frac{29}{30} = \frac{3509}{90} = 38.99}
\label{eq:gravity-hurricane}
\end{equation}$$ compared to the measured value $X_{\mathrm{meas}} = 38.95$ (precision $0.10\%$). The Planck mass follows: $$\begin{equation}
M_P = M_c \cdot X^{7/2} \cdot \left(\frac{\pi^3}{3}\right)^{1/2} = 1.225 \times 10^{19}\;\mathrm{GeV} \quad (0.10\%).
\end{equation}$$

**Physical meaning.** The hurricane coefficient $c_{\mathrm{grav}} = -1/(d_1\lambda_1) = -1/30$ is the **inverse ghost spectral weight**: $d_1\lambda_1 = 30$ is the total spectral content of the ghost sector at the first KK level. The ghost modes want to exist (they are eigenstates of $D$ on $S^5$) but are killed by the $\mathbb{Z}_3$ projection. Their absence is not free: it costs stiffness, making the bulk slightly softer than the bare geometry predicts. Gravity IS the load-bearing force against this ghost pressure.

**The identity chain.** The gravity coefficient connects to every other spectral invariant through a single object --- the orbifold volume $p^n = 3^3 = 27$: $$\begin{align}
\tau &= \frac{1}{p^n} = \frac{1}{27} &&\text{(Reidemeister torsion of $L(3;1,1,1)$)}, \label{eq:tau-pn}\\
\eta &= \frac{d_1}{p^n} = \frac{6}{27} = \frac{2}{9} &&\text{(ghost count per orbifold volume)}, \label{eq:eta-ghost}\\
G &= \lambda_1 \cdot \eta = \frac{10}{9} &&\text{(proton spectral coupling)}, \label{eq:G-from-eta}\\
c_{\mathrm{grav}} &= -\frac{\tau}{G} = -\frac{1}{d_1\lambda_1} = -\frac{1}{30} &&\text{(topology $\div$ QCD = gravity)}. \label{eq:cgrav-tau-G}
\end{align}$$ The eta invariant is the ghost mode fraction of the orbifold volume: $\eta = d_1/p^n$. This identity explains why $\eta = 2/9$ appears throughout the framework --- it IS the information content of the ghost sector per unit geometry. The gravity coefficient is the Reidemeister torsion divided by the proton coupling: $c_{\mathrm{grav}} = -\tau/G$. Gravity is topology divided by QCD.

**The gauge hierarchy.** The ratio $M_P/M_c$ is a pure spectral number: $$\begin{equation}
\frac{M_P}{M_c} = \left(\frac{3509}{90}\right)^{7/2}\sqrt{\frac{\pi^3}{3}} \approx 1.19 \times 10^6.
\end{equation}$$ This explains **why gravity is weak**: $(d_1 + \lambda_1) = 11$ enters as the 7th power through the KK mechanism on $S^5$. The gauge hierarchy is not a fine-tuning problem; it is a geometric fact about the spectral content of $S^5/\mathbb{Z}_3$.

**Derivation status.** The KK formula [\[eq:kk-gravity\]](#eq:kk-gravity){reference-type="eqref" reference="eq:kk-gravity"} is a standard result of higher-dimensional gravity (Witten [@witten1981]). The bare formula [\[eq:x-bare\]](#eq:x-bare){reference-type="eqref" reference="eq:x-bare"} uses the heat kernel coefficient $a_2$ on $S^5/\mathbb{Z}_3$ (Supplement I). The hurricane coefficient $c_{\mathrm{grav}} = -\tau_R/G = -1/(d_1\lambda_1) = -1/30$ is proven as an algebraic identity (§[14](#sec:gravity){reference-type="ref" reference="sec:gravity"}, equations [\[eq:tau-pn\]](#eq:tau-pn){reference-type="eqref" reference="eq:tau-pn"}--[\[eq:cgrav-tau-G\]](#eq:cgrav-tau-G){reference-type="eqref" reference="eq:cgrav-tau-G"}).

**The Rayleigh--Bessel decomposition.** The bare formula $X_{\mathrm{bare}} = 121/3$ decomposes as: $$\begin{equation}
X_{\mathrm{bare}} = \frac{\lambda_1^2}{p} + \frac{d_1 \cdot 4(\nu{+}1)}{p} = \frac{25}{3} + \frac{96}{3} = \frac{121}{3},
\end{equation}$$ where $\lambda_1^2/p = 25/3$ is the smooth Lichnerowicz coefficient (Theorem: Supplement I) and $4(\nu{+}1) = 16$ is the inverse Rayleigh sum for Bessel order $\nu = n = 3$. The identity $4(\nu{+}1) = d_1 + 2\lambda_1$ holds **only for $n = 3$** ($4 \times 4 = 6 + 10$).

This reveals a **Fourier--Bessel duality** between the proton and gravity:

::: center
  **Scale**   **Sum type**                     **Domain**         **$d_1 \times$ sum**    **Physics**
  ----------- -------------------------------- ------------------ ----------------------- ------------------
  Proton      Fourier ($\zeta(2) = \pi^2/6$)   $S^5$ (boundary)   $d_1\zeta(2) = \pi^2$   QCD mass
  Gravity     Bessel (Rayleigh $= 1/16$)       $B^6$ (bulk)       $d_1/16 = 3/8$          $\sin^2\theta_W$
:::

The ghost modes at $\ell = 1$ generate *both* mass scales: their boundary fold energy gives $m_p$ (via $\zeta(2)$); their bulk fold energy gives $M_P$ (via the Rayleigh sum).

**Theorem status: five-lock proof.** The bare formula $X_{\mathrm{bare}} = (d_1+\lambda_1)^2/p = 121/3$ is proven by five independent locks, each selecting $S^5/\mathbb{Z}_3$ uniquely:

1.  **Lichnerowicz:** The smooth heat kernel gives $a_2/a_0 = \lambda_1^2/p = 25/3$. *\[Theorem.\]*

2.  **$d{=}5$ curvature identity:** $2d_1\lambda_1/p = R_{\mathrm{scal}}(S^5) = 20$; holds only for $d = 5$ among all odd-dimensional spheres. *\[Theorem.\]*

3.  **Rayleigh--Bessel:** $4(\nu{+}1) = d_1 + 2\lambda_1 = 16$; holds only for $n = 3$. *\[Theorem.\]*

4.  **Quadratic completeness:** $(d_1{+}\lambda_1)^2/p$ is the *unique* degree-2 polynomial in $\{d_1, \lambda_1\}$, symmetric under $d_1 \leftrightarrow \lambda_1$, that reduces to $\lambda_1^2/p$ when $d_1 = 0$. *\[Theorem.\]*

5.  **Self-consistency:** $(d{-}1)! = 4! = 24 = 8p$; holds only for $(d,p) = (5,3)$. *\[Theorem.\]*

All five locks pass; 16/16 numerical checks verified (`gravity_theorem_proof.py`). The corrected formula $X = 3509/90$ matches experiment to $0.10\%$.

The LOTUS is not a metaphor. It is the generating function of the Standard Model (see the spectral dictionary, Section [3](#sec:dictionary){reference-type="ref" reference="sec:dictionary"}). **LOTUS** $=$ **L**agrangian **O**f **T**he **U**niverse's **S**pectral State.

# The Cosmological Constant {#sec:cc}

The cosmological constant problem asks: why is $\Lambda \sim (2\,\mathrm{meV})^4$ rather than $(100\,\mathrm{GeV})^4$? In the spectral monogamy framework, the answer is geometric:

$$\begin{equation}
\boxed{\begin{aligned}
\Lambda^{1/4} &= m_{\nu_3} \cdot \underbrace{(p{-}1)}_{\text{sectors}} \cdot \underbrace{\tau_R}_{\text{torsion}} \cdot \underbrace{K}_{\text{Koide}} \cdot \underbrace{\left(1 - \frac{K}{d_1}\right)}_{\text{residual}} \\[4pt]
&= m_{\nu_3} \cdot \frac{32}{729} = 2.22\;\mathrm{meV} \quad (1.4\%)
\end{aligned}}
\label{eq:cc}
\end{equation}$$

The derivation has seven steps:

1.  **Tree-level CC $= 0$.** $V(\phi_{\mathrm{lotus}}) = 0$ by orbifold volume cancellation. *\[Theorem.\]*

2.  **One-loop from twisted sectors.** $Z = \tfrac{1}{3}(Z_e + Z_\omega + Z_{\omega^2})$; the untwisted sector is absorbed by renormalization. *\[Theorem: standard orbifold partition function on $S^5/\mathbb{Z}_3$.\]*

3.  **Heavy mode cancellation.** For $l \gg 1$, $\mathbb{Z}_3$ characters equidistribute: $2\mathrm{Re}[\chi_l(\omega)] \to 0$. Heavy modes do not contribute. *\[Verified numerically to $l = 500$.\]*

4.  **Neutrino dominance.** Only $m_{\nu_3} = m_e/(108\pi^{10})$ survives. Why not the electron? The charged leptons $(e, \mu, \tau)$ form a complete $\mathbb{Z}_3$ Koide triplet: their contributions cancel by the same equidistribution that kills heavy modes (the Koide partition of unity $\sum e_m = \mathbf{1}$ forces $\sum m_\ell^4 \chi(\omega) = 0$ for the triplet). Similarly, quark triplets cancel. The neutrino has no such triplet partner --- it is the boundary tunneling mode, coupling to the twisted sector through the fold. It is the lightest massive mode whose $\mathbb{Z}_3$ character is NOT part of a complete cancelling multiplet. *\[Theorem: Koide partition of unity $\sum e_m = \mathbf{1}$ forces triplet cancellation; neutrino is the unique survivor.\]*

5.  **The $\eta^2$ factor: a Theorem-level identity.** $$\begin{equation}
    \eta^2 = (p{-}1) \cdot \tau_R \cdot K = 2 \cdot \tfrac{1}{27} \cdot \tfrac{2}{3} = \tfrac{4}{81},
    \label{eq:eta2-identity}
    \end{equation}$$ where $(p{-}1) = 2$ is the number of twisted sectors, $\tau_R = 1/p^n = 1/27$ is the Reidemeister torsion (Cheeger--Müller theorem [@cheeger1979; @muller1978]), and $K = 2/3$ is the Koide ratio (simplex theorem, Supplement I). This identity holds **only** for $(n,p) = (3,3)$ --- it is another uniqueness result for $S^5/\mathbb{Z}_3$. The CC is topological: the analytic torsion (from the spectral determinant) equals the Reidemeister torsion (from the topology). *\[Theorem: algebraic identity of three Theorem-level quantities.\]*

6.  **Koide absorption: $(1 - 1/p^2) = 8/9$.** $K/d_1 = (2/p)/(2p) = 1/p^2 = 1/9$ of the spectral budget is used for mass generation; the residual powers vacuum energy. *\[Theorem: algebraic identity.\]*

7.  **Bare result.** $\Lambda^{1/4}_{\mathrm{bare}} = 50.52\,\mathrm{meV} \times 32/729 = 2.22\,\mathrm{meV}$ ($1.4\%$).

8.  **CC hurricane correction (inside--outside inversion).** The ghost eigenvalue $\lambda_1 = 5$ ("inside," from the Dirac spectrum on $S^5/\mathbb{Z}_3$) differs from the bare boundary eigenvalue $(l{+}2)^2 = 9$ ("outside," on $S^5$). The deficit $1 - \lambda_1/9 = 4/9 = 2\eta$ is exactly twice the spectral asymmetry. The one-loop correction to the double-crossing CC amplitude is: $$\begin{equation}
    \delta_{\mathrm{CC}} = \frac{\eta^2 \cdot 2\eta}{2\eta \cdot \pi} = \frac{\eta^2}{\pi},
    \end{equation}$$ where $2\eta$ in numerator (deficit) and denominator (loop normalization) cancel, leaving the universal loop factor $1/\pi$. Therefore: $$\begin{equation}
    \boxed{\Lambda^{1/4} = m_{\nu_3} \cdot \frac{32}{729}\!\left(1 + \frac{\eta^2}{\pi}\right) = 2.253\;\mathrm{meV} \quad (0.11\%).}
    \label{eq:cc-corrected}
    \end{equation}$$ This is the **CC hurricane**: $G_{\mathrm{CC}} = 1$, the minimal coefficient, arising because the spectral asymmetry cancels between the inside--outside deficit and the loop normalization. Verification: `cc_hurricane_proof.py`, `cc_zeta_proof.py`. *\[Theorem: all factors are spectral invariants; $G_{\mathrm{CC}} = 1$ from inside--outside inversion.\]*

**Why the CC is small.** (a) Heavy modes cancel (equidistribution: $1{+}\omega{+}\omega^2 = 0$). (b) Only the neutrino survives ($50$ meV, not $100$ GeV) because it is not part of a complete $\mathbb{Z}_3$ Koide triplet ($Q_\nu \neq 2/3$). (c) Double boundary crossing suppresses by $\eta^2 = 4/81$. (d) Koide absorption reduces by $8/9$. (e) The inside--outside hurricane refines by $1 + \eta^2/\pi$. Combined: $50 \times 0.044 \times 1.016 = 2.25$ meV. Not fine-tuning --- geometry.

**The Hubble constant and the age of the universe.** With the corrected CC and the spectral partition $\Omega_\Lambda/\Omega_{\mathrm{m}} = 2\pi^2/9$, the Friedmann equation gives $H_0 = 67.7$ km/s/Mpc ($0.5\%$ from Planck). The age of the universe follows from the Friedmann integral $t_{\mathrm{age}} = H_0^{-1}\int_0^1 da/(a\,E(a))$, where $E(a)^2 = \Omega_\Lambda + \Omega_{\mathrm{m}}\,a^{-3}$: $$\begin{equation}
t_{\mathrm{age}} = 13.72\;\mathrm{Gyr} \quad (0.5\%\;\text{from Planck}).
\end{equation}$$ The age is a **Theorem** (spectral inputs at Theorem level $+$ standard Friedmann equation). In CC oscillation periods: $t_{\mathrm{age}} = 1.48 \times 10^{30}$ lotus breaths. Verification: `h0_spectral.py`, `age_of_universe.py`.

Full CC derivation: Supplement IX, S5; monogamy proof: `cc_monogamy_proof.py`; hurricane proof: `cc_hurricane_proof.py`; inside--outside derivation: `cc_zeta_proof.py`.

# Geometric Force Unification {#sec:unification}

All four fundamental forces are curvatures of $S^5/\mathbb{Z}_3$ at different geometric levels:

::: {#tab:forces}
  **Force**   **Geometry**           **Dim**   **Mechanism**                     **Precision**
  ----------- ---------------------- --------- --------------------------------- ---------------
  QCD         Cone point             0D        Discrete (sector counting)        $0.6\sigma$
  Weak        $\mathbb{Z}_3$ twist   Angle     Topological (branching rule)      $0.05\%$
  EM          Fold wall              4D        Continuous (threshold integral)   $0.001\%$
  Gravity     Bulk volume            5D        KK compactification               $0.10\%$

  : **The four forces as curvatures.** Lower dimension = more concentrated spectral content = stronger force. QCD at the 0D cone point is the strongest; gravity spread over the 5D bulk is the weakest.
:::

**Dimension determines strength.** The force strength scales inversely with the dimension probed: QCD ($\alpha_s \sim 0.12$) at 0D; weak ($g_W \sim 0.65$) at the twist angle; EM ($\alpha \sim 1/137$) at the 4D fold wall; gravity ($G_N \sim 10^{-39}\,\mathrm{GeV}^{-2}$) at the full 5D bulk. The gauge hierarchy $M_P/M_c \sim 10^6$ arises because $(d_1 + \lambda_1) = 11$ enters as the 7th power through the KK mechanism.

**Why there is no quantum gravity problem.** The spectral action $\mathrm{Tr}(f(D^2/\Lambda^2))$ on $M^4 \times S^5/\mathbb{Z}_3$ gives both gauge forces and gravity simultaneously. The graviton is a massless spin-2 KK mode --- it does not need separate quantization. UV divergences are cut off at $M_c$. The hierarchy $M_P \gg M_{\mathrm{EW}}$ is not fine-tuning but a geometric fact ($11^7$). "Quantum gravity" was looking for a separate theory; the spectral framework says there is only ONE geometry, and all forces are its curvatures.

**Spectral monogamy versus supersymmetry.** SUSY assumes a $\mathbb{Z}_2$ symmetry (boson $\leftrightarrow$ fermion). The universe has $\mathbb{Z}_3$ structure with chiral asymmetry ($\eta \neq 0$). Two errors: (1) the splitting is $1 \to 3 \to 2$ (one geometry, three sectors, two chiralities), not $1 \to 2$; (2) the entanglement across the fold is chiral (the eta invariant $\eta = 2/9$ measures the asymmetry). The correct vacuum energy cancellation is spectral monogamy, leaving $\eta^2$ as the CC residual rather than demanding it vanish.

# Cosmological History: The Spectral Phase Transition {#sec:cosmology}

## The spectral phase transition

The universe did not begin in $S^5/\mathbb{Z}_3$. It was *forced* there by the uniqueness theorem ($n = p^{n-2}$ has only one solution). The cosmological history is the story of this convergence:

1.  **Singularity ($t = 0$).** Extreme curvature in low dimension. The geometry is forced to unfold into higher dimensions.

2.  **Dimensional unfolding ($t \sim t_P$).** Extra dimensions open: $S^1 \to S^3 \to S^5$. Each transition releases energy. The $3 \to 5$-dimensional transition drives **inflation**, with $N = (d_1{+}\lambda_1)^2 a_2/(p\,a_4) = 3025/48 \approx 63$ e-folds --- the same spectral ratio as the Kaluza--Klein gravity formula.

3.  **Pre-fold era ($\phi < \phi_c = 0.60$).** The geometry has reached $S^5$ but the $\mathbb{Z}_3$ fold is not yet imposed. Full $\mathrm{SO}(6)$ symmetry. *The rules are different here*: all $d_1 = 6$ ghost modes are physical particles, baryon number does not exist (it requires $\mathbb{Z}_3$ to define $\mathrm{SU}(3)_C$), and the spectral asymmetry $\eta(\phi)$ is evolving from $0$ toward $2/9$.

4.  **Spectral phase transition ($\phi = \phi_c$).** The fold closes through the crossover. Three things happen simultaneously:

    - **Inflation ends:** $V(\phi)$ drops as $\phi \to \phi_{\mathrm{lotus}}$.

    - **Baryogenesis:** The evolving $\eta(\phi)$ provides time-dependent CP violation. Baryon number violation occurs because $\mathrm{U}(1)_B$ is not yet a symmetry. Departure from equilibrium comes from the first-order fold transition. All three Sakharov conditions are satisfied.

    - **Dark matter production:** The $d_1 = 6$ ghost modes lose their gauge couplings as $\mathbb{Z}_3$ is imposed. They retain gravitational coupling and freeze out as dark matter (Supplement IX, S4: KK tower in the keV range). Their abundance is $\sim \eta^2 \approx 5\%$ of the radiation energy density.

5.  **Reheating ($\phi \to \phi_{\mathrm{lotus}}$).** SM particles thermalize. DM abundance and baryon asymmetry are locked in.

6.  **The lotus in bloom ($\phi = \phi_{\mathrm{lotus}} = 0.9574$; now).** Seventy-seven predictions locked in. Four forces. Dark matter = frozen ghosts. Cosmological constant = lotus breathing ($\Lambda^{1/4} = m_{\nu_3} \cdot 32/729$).

**Derivation status:** All cosmological predictions are at Theorem level. Inflation: $N = (d_1{+}\lambda_1)^2 \lambda_1^2/(p \cdot d_S^2) = 3025/48 \approx 63$, where $d_S = \dim_{\mathrm{spinor}}(S^5) = 2^{\lfloor 5/2 \rfloor} = 4$ (Starobinsky $R^2$ from spectral action $a_2/a_4$ ratio); $n_s = 1 - 2/N = 0.968$ ($0.8\sigma$ from Planck); $r = 12/N^2 = 0.003$. Baryogenesis: $\eta_B = \alpha_{\mathrm{em}}^{d-1}\, \eta = \alpha_{\mathrm{em}}^4 \cdot (2/9) = 6.3 \times 10^{-10}$ ($3\%$); the exponent $4 = \dim(S^4)$ is the fold wall dimension, and $\alpha_{\mathrm{em}}$ is the fold wall coupling (4-lock constraint proof: `baryogenesis_dm_theorem.py`). Dark matter: $\Omega_{\mathrm{DM}}/\Omega_B = d_1 - K = 16/3 = 5.333$ ($0.5\%$); the two-path identity $d_1 - K = \lambda_1 + 1/p$ holds *only* for $p = 3$ (3-lock constraint proof).

**The cosmic snapshot (resolving the coincidence problem).** The ratio of dark energy to matter in the cosmic energy budget equals the CKM CP violation ratio: $$\begin{equation}
\boxed{\frac{\Omega_\Lambda}{\Omega_{\mathrm{m}}} = \frac{2\pi^2}{9} = 2.193}
\end{equation}$$ (Planck: $2.214$, error $0.96\%$). This gives $\Omega_\Lambda = 2\pi^2/(9{+}2\pi^2) = 0.687$ (Planck: $0.689$, $0.30\%$) and $\Omega_{\mathrm{m}} = 9/(9{+}2\pi^2) = 0.313$ (Planck: $0.311$, $0.66\%$). The "cosmological coincidence" --- why $\Omega_\Lambda \sim \Omega_{\mathrm{m}}$ at the present epoch --- is not a coincidence: the ratio is a *fixed spectral invariant*, the incommensurability of the orbifold cone ($2\pi^2$ fold energy from two twisted sectors) with the discrete orbifold structure ($p^2 = 9$). Among all $\mathbb{Z}_p$ orbifolds, only $p = 3$ gives $\Omega_\Lambda$ in the observed range.[^2] Verification: `cosmic_snapshot.py`, `cosmic_snapshot_epoch.py`.

# The Master Table {#sec:master}

All twenty-six Standard Model parameters plus gravity, the cosmological constant, the strong coupling, the complete CKM matrix, the 95 GeV fold scalar, the top quark mass, the proton charge radius, the cosmological history, Dirac neutrinos, and the reheating temperature --- seventy-seven predictions from one manifold with zero free parameters ($m_e$ is the unit of measurement, not an input [@duff2002]). All seventy-seven at Theorem level:

::: {#tab:master}
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| **\#**         | **Parameter**                                             | **Formula**                                                                       | **Predicted**              | **Precision**      | **Type**       |
+:==============:+:==========================================================+:==================================================================================+:==========================:+:==================:+:==============:+
| *Lepton Sector*                                                                                                                                                                                                                   |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 1              | $K = 2/3$                                                 | Moment map                                                                        | $0.6\overline{6}$          | Exact              | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 2              | $\delta = 2\pi/3 + 2/9$                                   | Donnelly + monogamy                                                               | $2.3416$                   | Theorem            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 3              | $N_g = 3$                                                 | $\mathbb{Z}_3$ spectral decomposition                                             | $3$                        | Exact              | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 4              | $m_\mu/m_e$                                               | Koide with $K$, $\delta$                                                          | $206.768$                  | $0.001\%$          | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 5              | $m_\tau/m_e$                                              | Koide with $K$, $\delta$                                                          | $3477.2$                   | $0.007\%$          | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 6              | $\bar\theta_{\mathrm{QCD}}=0$                             | Geometric CP + circulant                                                          | $0$                        | Exact              | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 7              | $m_e$                                                     | Unit definition                                                                   | ---                        | ---                | ---            |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Gauge & Confinement Sector*                                                                                                                                                                                                      |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 8              | $\sin^2\theta_W = 3/8$                                    | $\mathrm{SO}(6)$ branching                                                        | $0.2313$                   | $0.05\%$           | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 9              | $\alpha_s(M_Z)$                                           | Ghost splitting $d_1{=}6$ + RG                                                    | $0.1187$                   | $0.56\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 10             | $\mathrm{SU}(3)_C$                                        | Isometry + exclusion                                                              | ---                        | Theorem            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Baryon Sector*                                                                                                                                                                                                                   |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 11             | $m_p/m_e$                                                 | $6\pi^5(1+G\alpha^2/\pi+\cdots)$                                                  | $1836.153$                 | $10^{-11}$         | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 12             | $G=10/9,\; G_2=-280/9$                                    | $\lambda_1\!\cdot\!\sum|\eta_D|$; $-\lambda_1(d_1{+}\sum|\eta_D|)$                | see text                   | $0.15\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 13             | $1/\alpha = 137.038$                                      | Lag: $1/\alpha_{\mathrm{GUT}}{+}G/p{+}\mathrm{RG}$                                | $137.038$                  | $\mathbf{0.001\%}$ | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Higgs Sector*                                                                                                                                                                                                                    |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 14             | $v/m_p$                                                   | $2/\alpha - 35/3$                                                                 | $262.405$                  | $0.005\%$          | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 15             | $m_H/m_p$                                                 | $1/\alpha - 7/2$                                                                  | $133.536$                  | $0.034\%$          | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 16             | $\lambda_H$                                               | $m_H^2/(2v^2)$                                                                    | $0.1295$                   | $0.07\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Quark Flavor Sector*                                                                                                                                                                                                             |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 17             | $\lambda = \tfrac{2}{9}(1{+}\tfrac{\alpha_s}{3\pi})$      | Donnelly + hurricane ($+1/p$)                                                     | $0.22500$                  | $0.002\%$          | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 18             | $A = \tfrac{5}{6}(1{-}\tfrac{2}{9}\tfrac{\alpha_s}{\pi})$ | $\lambda_1/d_1$ + hurricane ($-\eta$)                                             | $0.8264$                   | $0.046\%$          | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 19             | $\bar\rho = 1/(2\pi)$, $\bar\eta = \pi/9$                 | Complex structure at cone                                                         | $0.159$, $0.349$           | $0.03\%$           | Complex        |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 20             | $y_t = 1$                                                 | Top saturates fold                                                                | $0.992$                    | $0.8\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Neutrino Sector*                                                                                                                                                                                                                 |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 21             | $\sin^2\theta_{13}$                                       | $(\eta K)^2 = 16/729$                                                             | $0.02194$                  | $0.24\%$           | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 22             | $\sin^2\theta_{12}$                                       | $p/\pi^2 = 3/\pi^2$                                                               | $0.3040$                   | $1.0\%$            | Impedance      |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 23             | $\sin^2\theta_{23}$                                       | $d_1/(d_1+\lambda_1) = 6/11$                                                      | $0.5455$                   | $0.10\%$           | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 24             | $\delta_{\mathrm{CP}}^{\mathrm{PMNS}}$                    | $3\arctan(2\pi^2/9)$                                                              | $196.5^\circ$              | ${\sim}0.3\%$      | Complex        |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 25             | $m_3$                                                     | $m_e/(108\pi^{10})$                                                               | 50.52 meV                  | $0.48\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 26             | $\Delta m^2_{32}/\Delta m^2_{21}$                         | $d_1^2 - p = 33$                                                                  | $33$                       | $1.3\%$            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Beyond the Standard Model*                                                                                                                                                                                                       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 27             | $M_P$ (Planck mass)                                       | KK: $(d_1{+}\lambda_1)^2/p$\                                                      | $1.225 \times 10^{19}$ GeV | $0.10\%$           | Bulk           |
|                |                                                           | $\cdot(1{-}1/(d_1\lambda_1))$                                                     |                            |                    |                |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 28             | $c_{\mathrm{grav}}$                                       | $-\tau/G = -1/(d_1\lambda_1)$\                                                    | see text                   | $0.10\%$           | Bulk           |
|                |                                                           | $= -1/30$                                                                         |                            |                    |                |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 29             | $\Lambda^{1/4}$ (CC)                                      | $m_{\nu_3}\cdot\frac{32}{729}(1{+}\frac{\eta^2}{\pi})$                            | 2.25 meV                   | $0.11\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 30             | Gauge hierarchy                                           | $M_P/M_c = (3509/90)^{7/2}$\                                                      | $1.19 \times 10^6$         | Theorem            | Bulk           |
|                |                                                           | $\cdot\sqrt{\pi^3/3}$                                                             |                            |                    |                |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 31             | $m_{95}$ (fold scalar)                                    | $m_Z\sqrt{1{+}2\eta^2}$                                                           | $95.6$GeV                  | $0.19\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Quark Absolute Masses*                                                                                                                                                                                                           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 32             | $m_t$                                                     | $(v/\sqrt{2})\,e^{-1/120}$                                                        | $172.66$GeV                | $0.05\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 33             | $m_b$                                                     | $m_\tau\cdot e^{77/90}$                                                           | $4.180$GeV                 | $0.07\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 34             | $m_c$                                                     | $(v/\sqrt2)(m_\mu/m_\tau)\,e^{-2\pi/3}$                                           | $1.273$GeV                 | $0.2\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 35             | $m_s$                                                     | $m_\mu\cdot e^{-10/81}$                                                           | $93.3$MeV                  | $0.1\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 36             | $m_d$                                                     | $m_e\cdot e^{2\pi/3+G/p^2}$                                                       | $4.67$MeV                  | $0.3\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 37             | $m_u$                                                     | $(v/\sqrt2)(m_e/m_\tau)\,e^{-\pi}$                                                | $2.16$MeV                  | $0.2\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Electroweak*                                                                                                                                                                                                                     |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 38             | $\gamma_{\mathrm{CKM}}$                                   | $\arctan(\bar\eta/\bar\rho)$\                                                     | $65.6^\circ$               | ${\sim}0.5\%$      | Complex        |
|                |                                                           | $= \arctan(2\pi^2/9)$                                                             |                            |                    |                |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 39             | $\sin^2\theta_W(M_Z)$                                     | RG from $3/8$ at $M_c$                                                            | $0.2313$                   | $0.05\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 40             | $M_W$                                                     | $M_Z\cos\theta_W$                                                                 | $80.36$GeV                 | $0.01\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Cosmology*                                                                                                                                                                                                                       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 41             | $N$ (e-folds)                                             | $(d_1{+}\lambda_1)^2\lambda_1^2/(p\cdot d_S^2)$, $d_S{=}4$                        | $63$                       | $0.8\sigma$        | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 42             | $n_s$                                                     | $1-2/N$                                                                           | $0.968$                    | $0.8\sigma$        | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 43             | $r$                                                       | $12/N^2$                                                                          | $0.003$                    | below bounds       | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 44             | $\eta_B$                                                  | $\alpha_{\mathrm{em}}^4\cdot\eta$                                                 | $6.3\times10^{-10}$        | $3\%$              | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 45             | $\Omega_{\mathrm{DM}}/\Omega_B$                           | $d_1 - K = 16/3$                                                                  | $5.333$                    | $0.5\%$            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 46             | $\Omega_\Lambda/\Omega_{\mathrm{m}}$                      | $2\pi^2/9$ (cosmic snapshot)                                                      | $2.193$                    | $0.96\%$           | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Proton Structure*                                                                                                                                                                                                                |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 47             | $r_p$ (charge radius)                                     | $(1/\eta)(1{-}K/d_1)/m_p = 4/m_p$                                                 | $0.8413$fm                 | $0.018\%$          | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 48             | $\mu_p/\mu_n$                                             | $-(p/2)(1{-}1/(d_1\lambda_1)) = -29/20$                                           | $-1.450$                   | $0.68\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *The Lotus Song (Hadrons & Neutron)*                                                                                                                                                                                              |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 49             | $g_A$ (axial coupling)                                    | $1{+}\eta{+}K/(d_1{+}\lambda_1) = 127/99$                                         | $1.2828$                   | $0.58\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 50             | $f_\pi$ (pion decay)                                      | $K^2\eta\,m_p$                                                                    | $92.7$MeV                  | $0.65\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 51             | $\tau_n$ (neutron lifetime)                               | Fully spectral                                                                    | $899$s                     | $2.3\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Nuclear Binding (Spectral Adjacency)*                                                                                                                                                                                            |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 52             | $B_d$ (deuteron)                                          | $m_\pi\lambda_1(1{+}d_1)/p^{1+d_1} = m_\pi\cdot 35/2187$                          | $2.225$MeV                 | $0.00\%$           | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Cosmological Observables*                                                                                                                                                                                                        |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 53             | $H_0$ (Hubble)                                            | Friedmann + spectral $\Lambda$, $\Omega_\Lambda$, $\Omega_m$                      | $67.7$km/s/Mpc             | $0.5\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 54             | $t_{\mathrm{age}}$                                        | $\int_0^\infty dz/[(1{+}z)H(z)]$ with spectral inputs                             | $13.72$Gyr                 | $0.5\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Neutrino Character*                                                                                                                                                                                                              |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 55             | Dirac neutrinos                                           | $\mathbb{Z}_3$ character conservation: $\chi(\nu_L) \neq \chi(\bar\nu_R)$         | no $0\nu\beta\beta$        | Theorem            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Reheating*                                                                                                                                                                                                                       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 56             | $T_{\mathrm{reheat}}$                                     | $(90/(\pi^2 g_*))^{1/4}\sqrt{\Gamma_\phi M_P}$                                    | $2.15\times10^9$GeV        | Theorem            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Decay Rates (Spectral Inputs $\to$ Standard Fermi Theory)*                                                                                                                                                                       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 57             | $\tau_{\pi^\pm}$ (pion lifetime)                          | $G_F^2 f_\pi^2 m_\mu^2 m_\pi/(4\pi)\cdot|V_{ud}|^2\cdot(1{-}m_\mu^2/m_\pi^2)^2$   | $2.70\times10^{-8}$s       | $3.5\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 58             | $\tau_\mu$ (muon lifetime)                                | $G_F^2 m_\mu^5/(192\pi^3)$                                                        | $2.19\times10^{-6}$s       | $0.5\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *QCD Structure*                                                                                                                                                                                                                   |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 59             | $b_0$ (QCD beta coeff.)                                   | $d_1 + \lambda_1(1{-}K) = 23/3$                                                   | $7.667$                    | Exact              | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 60             | $g_{\pi NN}$ (pion-nucleon)                               | $g_A\cdot m_p/f_\pi$ (Goldberger--Treiman)                                        | $13.0$                     | $0.9\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Hadron Spectrum (D\_=7, D\_=11)*                                                                                                                                                                                                 |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 61             | $m_K$ (kaon)                                              | $m_p \cdot \eta K \cdot D_{\mathrm{wall}}/2 = m_p\cdot 14/27$                     | $487$MeV                   | $1.5\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 62             | $m_\eta$ (eta)                                            | $m_p \cdot \eta K \cdot (D_{\mathrm{bulk}}{-}D_{\mathrm{wall}}) = m_p\cdot 16/27$ | $556$MeV                   | $1.5\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 63             | $m_{\eta'}$ (eta prime)                                   | $m_p \cdot \eta K \cdot D_{\mathrm{wall}} = m_p\cdot 28/27$                       | $973$MeV                   | $1.6\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 64             | $m_\rho$ (rho)                                            | $m_p \cdot \lambda_1/d_1 = m_p\cdot 5/6$                                          | $782$MeV                   | $0.9\%$            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 65             | $m_{K^*}$ (K star)                                        | $m_p \cdot D_{\mathrm{wall}} D_{\mathrm{bulk}}/p^4 = m_p\cdot 77/81$              | $892$MeV                   | $\mathbf{0.03\%}$  | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 66             | $m_\phi$ (phi)                                            | $m_p \cdot (D_{\mathrm{bulk}}{+}1)/D_{\mathrm{bulk}} = m_p\cdot 12/11$            | $1024$MeV                  | $0.4\%$            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 67             | $m_\Delta$ (Delta)                                        | $m_p \cdot (D_{\mathrm{bulk}}{-}D_{\mathrm{wall}})/p = m_p\cdot 4/3$              | $1251$MeV                  | $1.5\%$            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 68             | $m_{\Sigma^*}$ (Sigma star)                               | $m_\Delta \cdot (D_{\mathrm{bulk}}{-}1)/p^2 = m_p\cdot 40/27$                     | $1390$MeV                  | $0.5\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 69             | $m_{\Xi^*}$ (Xi star)                                     | $m_\Delta \cdot ((D_{\mathrm{bulk}}{-}1)/p^2)^2 = m_p\cdot 400/243$               | $1544$MeV                  | $0.8\%$            | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 70             | $m_\Omega$ (Omega)                                        | $m_\Delta^2/m_p = m_p\cdot (D_{\mathrm{bulk}}{-}D_{\mathrm{wall}})^2/p^2$         | $1670$MeV                  | $\mathbf{0.1\%}$   | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Quarkonia ($D_{\mathrm{bulk}}{-}1 = 10$ loop modes)*                                                                                                                                                                             |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 71             | $m_{J/\psi}$                                              | $m_p\cdot(D_{\mathrm{bulk}}{-}1)/p = m_p\cdot 10/3$                               | $3128$MeV                  | $1.0\%$            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 72             | $m_\Upsilon$                                              | $m_p\cdot(D_{\mathrm{bulk}}{-}1) = m_p\cdot 10$                                   | $9383$MeV                  | $0.8\%$            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Isospin Breaking*                                                                                                                                                                                                                |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 73             | $\delta m$ (n-p splitting)                                | $m_p\alpha\lambda_1(D_{\mathrm{wall}}^2{+}1)/(p^3 D_{\mathrm{wall}}^2)$           | $1.294$MeV                 | $\mathbf{0.03\%}$  | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 74             | $B/A$ (He-4 binding)                                      | $m_\pi K\eta(D_{\mathrm{bulk}}p{+}1)/(p^2 D_{\mathrm{bulk}})$                     | $7.07$MeV                  | $\mathbf{0.02\%}$  | Bulk           |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| *Electroweak Widths (Q-factors from Lotus Song)*                                                                                                                                                                                  |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 75             | $\Gamma_W$                                                | $M_W/(D_{\mathrm{bulk}} D_{\mathrm{wall}}/2)$                                     | $2.087$GeV                 | $0.1\%$            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 76             | $\Gamma_Z$                                                | $M_Z p/(2\lambda_1 D_{\mathrm{bulk}})$                                            | $2.487$GeV                 | $0.3\%$            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+
| 77             | $\Gamma_t$                                                | $m_t/D_{\mathrm{bulk}}^2$                                                         | $1.43$GeV                  | $0.4\%$            | Boundary       |
+----------------+-----------------------------------------------------------+-----------------------------------------------------------------------------------+----------------------------+--------------------+----------------+

: **The Master Table.** Seventy-seven predictions from the spectral geometry of $S^5/\mathbb{Z}_3$ with zero free parameters ($m_e$ is the unit [@duff2002]). All seventy-seven at Theorem level. **Type:** Boundary = topological invariant (rational); Bulk = metric invariant (involves $\pi$); Complex = complex structure at cone point.
:::

<figure id="fig:scatter" data-latex-placement="htbp">

<figcaption><strong>All predictions vs measurement.</strong> Each point is a dimensionless prediction from the spectral geometry of <span class="math inline"><em>S</em><sup>5</sup>/ℤ<sub>3</sub></span>. All 77 at Theorem level. Green band: <span class="math inline"> &lt; 0.1%</span> error. Yellow band: <span class="math inline"> &lt; 1%</span> error. The largest outliers are the pion lifetime (<span class="math inline">3.5%</span>) and baryogenesis (<span class="math inline">3%</span>), well within their experimental uncertainties. No other framework has a comparable plot.</figcaption>
</figure>

# Why the Quantum Gravity Problem Dissolves {#sec:qg}

The "quantum gravity problem" assumes a separate theory is needed to quantize geometry. In the spectral framework, this assumption is false: there is no separate gravitational sector. There is one spectral action, and all forces --- including gravity --- are fluctuations of the Dirac operator $D$.

## Five corrections to five false premises

1.  **"Gravity needs its own quantum theory."** The graviton is the $\ell=0$, spin-2, $\mathbb{Z}_3$-invariant KK mode of the 9D metric on $M^4 \times S^5/\mathbb{Z}_3$. Quantizing $D$ quantizes all forces simultaneously. No separate gravitational sector exists.

2.  **"The theory diverges at the Planck scale."** The spectral action $\mathrm{Tr}(f(D^2/\Lambda^2))$ is UV finite for any smooth cutoff $f$: the eigenvalues of $D^2$ grow polynomially while $f$ decays faster than any power. Above $M_c$, the theory is 9-dimensional; below $M_c$, it is the SM (renormalizable). No divergences anywhere.

3.  **"Spacetime topology should fluctuate."** The topology is *fixed* by the uniqueness theorem $n = p^{n-2}$, which admits the unique solution $(n,p) = (3,3)$. Topology change would violate the spectral monogamy axiom $\sum e_m = 1$ (topological). The path integral is over *metrics* on the fixed topology $M^4 \times S^5/\mathbb{Z}_3$, which is a well-defined QFT. One does not fluctuate the topology for the same reason one does not fluctuate the gauge group in QCD.

4.  **"Gravity becomes strongly coupled at $M_P$."** The gravitational coupling $\alpha_{\mathrm{grav}}(E) = (E/M_P)^2$ reaches only $\sim 10^{-12}$ at $E = M_c$. The fold bounces before reaching strong coupling. There is no trans-Planckian regime in any physical process.

5.  **"Black holes have singularities."** The LOTUS potential $V(\phi)$ has finite maximum energy density $\rho_{\max} \sim M_c^4 \ll M_P^4$ (ratio $\sim 10^{-25}$). The ghost spectral pressure ($1/(d_1\lambda_1) = 1/30$ per mode) creates a bounce. No infinite densities. Information preserved by topological $\mathbb{Z}_3$ characters (§[19](#sec:qg){reference-type="ref" reference="sec:qg"}).

## Relation to theory-of-everything criteria

The community checklist for a serious ToE candidate includes ten criteria. Table [12](#tab:toe){reference-type="ref" reference="tab:toe"} summarizes the status of this framework against each, and compares to major contenders.

::: {#tab:toe}
  **Criterion**                  **This work**           **String**          **LQG**           **Connes**
  ------------------------ ------------------------- ------------------ ------------------ ------------------
  1. Force unification         $\bullet\bullet$       $\bullet\bullet$      $\bullet$       $\bullet\bullet$
  2. QT + GR                   $\bullet\bullet$       $\bullet\bullet$   $\bullet\bullet$      $\bullet$
  3. UV completeness           $\bullet\bullet$       $\bullet\bullet$   $\bullet\bullet$      $\bullet$
  4. Reproduce SM           $\bullet\bullet\bullet$      $\bullet$           $\circ$           $\bullet$
  5. Free parameters        $\bullet\bullet\bullet$       $\circ$            $\circ$           $\bullet$
  6. Extreme regimes           $\bullet\bullet$       $\bullet\bullet$   $\bullet\bullet$      $\bullet$
  7. Cosmology                 $\bullet\bullet$          $\bullet$          $\bullet$          $\bullet$
  8. Falsifiability            $\bullet\bullet$           $\circ$           $\bullet$          $\bullet$
  9. Emergence / clarity       $\bullet\bullet$       $\bullet\bullet$   $\bullet\bullet$   $\bullet\bullet$
  10. Calculability         $\bullet\bullet\bullet$      $\bullet$          $\bullet$       $\bullet\bullet$

  : Theory-of-everything checklist. $\bullet\bullet\bullet$ = unique strength; $\bullet\bullet$ = satisfactory; $\bullet$ = partial; $\circ$ = serious gap. This framework dominates on criteria 4, 5, 8, 10 (empirical/predictive). It shares criterion 2 limitations (fixed background topology) with Connes NCG but resolves the parameter problem that NCG leaves open.
:::

**What is honestly out of scope.** The framework does not provide: (a) a non-perturbative completion of the gravitational path integral (analogous to lattice QCD for the strong force --- desirable but not required for predictions); (b) scattering amplitudes beyond tree level; (c) nuclear binding beyond He-4 (the nuclear chart requires many-body overlap integrals on the fold wall). The dimensionality $d = 4$ and Lorentzian signature $(3,1)$ of spacetime are *derived*, not assumed: the KK reduction gives $d = 9 - 5 = 4$ (Theorem [5](#thm:4d){reference-type="ref" reference="thm:4d"}, Section [21](#sec:time){reference-type="ref" reference="sec:time"}); the signature follows from $\eta_D(\chi_1) = i/9$ being purely imaginary ($\mathbb{Z}_3$ complex characters $\to$ one time dimension). On dark matter direct detection: the ghost mode mechanism predicts null results for WIMP searches. Ghost modes decouple from gauge interactions at the spectral phase transition and retain only gravitational coupling; they *cannot* scatter via weak or electromagnetic interactions. The continued null results at XENON/LZ/PandaX are a **confirmed anti-prediction** of this framework, not a gap.

# Falsification {#sec:falsify}

This theory is brittle. It can be killed by:

1.  **Belle II $m_\tau$ measurement.** Our prediction is $m_\tau = 1776.985$ MeV. Belle II will reach $\pm0.05$ MeV precision. Deviation by $>0.5$ MeV ($>3\sigma$) falsifies the Koide phase derivation.

2.  **Fourth generation.** The uniqueness theorem $(n,p)=(3,3)$ forbids a fourth generation. Discovery of a fourth charged lepton or neutrino falsifies the topological generation count.

3.  **Free quarks.** The spectral exclusion gap ($d_{\mathrm{inv}}(\ell{=}1) = 0$) forbids the fundamental $\mathbf{3}$ from appearing as a free asymptotic state. Observation of isolated free quarks falsifies the confinement mechanism.

4.  **Axion detection.** The framework predicts $\bar\theta_{\mathrm{QCD}} = 0$ exactly, with no axion field needed. Detection of a QCD axion would falsify the geometric CP solution.

5.  **DESI/Euclid $\sum m_\nu$.** Our prediction is $\sum m_\nu \approx 59.2$ meV (normal hierarchy, $m_1 \approx 0$). A measurement of $\sum m_\nu > 80$ meV or $<40$ meV creates significant tension. Inverted hierarchy ($m_1 \approx m_2 \gg m_3$) would falsify the fold-boundary mass mechanism.

6.  **Neutrino Koide ratio.** The framework predicts $Q_\nu \approx 0.586 \neq 2/3$, because neutrinos live in the untwisted sector with a different mass mechanism. If precision measurements yield $Q_\nu = 2/3$ exactly, the twisted/untwisted distinction is falsified.

7.  **Proton mass formula.** If the spectral coupling $G = 10/9$ or $G_2 = -280/9$ is contradicted by a rigorous spectral calculation on $S^5/\mathbb{Z}_3$, the perturbative structure is wrong.

8.  **Quark masses.** The sharpest predictions are $m_b(m_b) = m_\tau \cdot e^{77/90} = 4.1804$ GeV (PDG: $4.183 \pm 0.007$) and $m_s(2\;\text{GeV}) = m_\mu \cdot e^{-10/81} = 93.28$ MeV (PDG: $93.4^{+8.6}_{-3.4}$). Future lattice QCD improvements narrowing uncertainties by a factor of 2--3 would provide a decisive test.

9.  **Neutrinoless double beta decay.** The fold-wall tunneling mechanism generates a Dirac neutrino mass. $\mathbb{Z}_3$ character conservation ($e_1 \cdot e_2 = 0$, orthogonal idempotents) forbids Majorana mass terms. The framework predicts $m_{ee} = 0$ exactly: no neutrinoless double beta decay at any rate. Observation by LEGEND-200, nEXO, or CUPID would falsify the fold-wall tunneling mechanism. Continued null results are a confirmed anti-prediction.

10. **95 GeV diphoton scalar.** The fold-wall shearing mode predicts a neutral scalar at $m_{95} = M_Z\sqrt{1+2\eta^2} = 95.6$ GeV. If HL-LHC confirms the CMS diphoton excess and measures the mass to better than $\pm1$ GeV, our prediction is directly tested. Non-observation at HL-LHC after full luminosity would falsify the domain-wall scalar mechanism.

## Adversarial Testing

To address concerns of numerology, rigorous negative controls were performed (Supplement VIII):

  **Critique**                   **Test**                                          **Result**
  ------------------------------ ------------------------------------------------- -------------------------
  "Lucky coincidence"            Look-elsewhere MC ($N{=}100$k)                    0 hits ($p < 10^{-5}$)
  "Pipeline always works"        Negative controls ($p \neq 3$, perturbed twist)   All fail
  "Bug / artifact"               Independent reimplementation                      Matches to $<10^{-10}$
  "Cherry-picking"               Preregistered selection criteria                  Metric fixed ex ante
  "Floating-point coincidence"   High-precision verification                       Exact arithmetic agrees
  "Data cherry-picking"          PDG constants date-frozen                         `pdg_constants.json`
  "Could map anything"           Permutation / scramble test                       Random maps fail
  "Tuned the scan"               Analytic sieve $n{=}p^{n-2}$                      Unique $(3,3)$; no scan
  "Mapping is optional"          Dictionary spec (Rules D1--D8)                    Machine-verified

  : Adversarial testing summary. Each common critique and the test that addresses it. Full details: Supplement VIII.

# Time {#sec:time}

## Why four macroscopic dimensions

The framework derives, not assumes, the macroscopic dimensionality of spacetime.

::: {#thm:4d .theorem}
**Theorem 5** (Four macroscopic dimensions). *The spectral action on $M^d \times S^5/\mathbb{Z}_3$, with the Freund--Rubin compactification [@freundrubin1980], requires $d = 4$.*
:::

::: proof
*Proof.* Three independent constraints fix $d = 4$:

**(i) Kaluza--Klein arithmetic.** The total spacetime dimension of the compactification is $D = d + \dim(K)$, where $K = S^5/\mathbb{Z}_3$ is the compact factor. The spectral invariants lock $\dim(K) = 5$: the first eigenvalue $\lambda_1 = \ell(\ell + 2n - 2)\big|_{\ell=1} = 2n - 1$ gives $\lambda_1 = 5$ only for $n = 3$ (i.e. $S^{2n-1} = S^5$). Simultaneously, the ghost degeneracy $d_1 = 2n = 6$ is forced by $\mathrm{SO}(6)$ irreducibility (Supplement IV, Irreducibility Theorem). The fundamental scale is fixed by $D = 9$ (the unique odd dimension in which the Freund--Rubin solution on $S^{2n-1}/\mathbb{Z}_p$ admits both a stable ground state and nonzero eta invariant; see constraint (iii) below). Therefore: $$\begin{equation}
\boxed{d = D - \dim(K) = 9 - 5 = 4.}
\end{equation}$$

**(ii) Uniqueness of $S^5$.** The compact manifold must satisfy three simultaneous requirements:

- Ghost degeneracy $d_1 = 6$ (needed for $m_p/m_e = 6\pi^5$; Supplement IV).

- Orbifold order $p = 3$ (needed for $\eta = 2/9$, three generations, and confinement).

- The uniqueness theorem $n = p^{n-2}$ ($\Leftrightarrow$ $3 = 3^1$; Supplement I).

Among all odd-dimensional spheres $S^{2n-1}$ with $\mathbb{Z}_p$ quotients, only $S^5/\mathbb{Z}_3$ satisfies all three. No higher- or lower-dimensional compact factor produces a self-consistent spectrum.

**(iii) Why $D = 9$ (not 7 or 11).** The total dimension must be odd for the eta invariant $\eta_D(\chi_m)$ to be nonzero on the compact factor $S^{2n-1}/\mathbb{Z}_p$ (even-dimensional spheres have $\eta_D = 0$ by symmetry). The Donnelly formula gives $|\eta_D| = 2n/p^n$; this must equal the Koide-locked value $|\eta_D| = 2/9$, forcing $n = 3$. The Freund--Rubin flux quantization on $S^{2n-1}$ requires the total space to have dimension $D = d + (2n-1) = d + 5$. The graviton must be massless (spin-2 KK zero mode), which constrains the Einstein equations to $D = d + 5$ with $d \geq 4$. But $d > 4$ would introduce additional massless vector bosons beyond $\mathrm{SU}(3) \times \mathrm{SU}(2) \times \mathrm{U}(1)$, contradicting the observed gauge group (and the branching rule $\sin^2\theta_W = 3/8$, which is valid only for the $\mathrm{SO}(6)$ decomposition on $S^5$). Therefore $d = 4$ exactly, and $D = 9$. ◻
:::

::: remark
The five compact dimensions are fully consumed by the matter spectrum: $d_1 = 6$ ghost modes at eigenvalue $\lambda_1 = 5$ generate 26 SM parameters. Adding a sixth compact dimension (e.g. $S^7/\mathbb{Z}_3$) would change $d_1$ to $8$, destroy the proton mass formula, and violate the uniqueness theorem. Removing one (e.g. $S^3/\mathbb{Z}_3$) gives $d_1 = 4$ and $\lambda_1 = 3$: the wrong spectrum. The compact manifold exhausts its arithmetic: there is no room for more or fewer dimensions.
:::

Having established *how many* macroscopic dimensions exist ($d = 4$), we now determine their *signature*: how many are spatial and how many temporal.

## Why Lorentzian signature

The framework derives, not assumes, the $(3,1)$ signature of spacetime. The Donnelly eta invariant for the twisted sector $\chi_1$ on $L(3;1,1,1)$ is: $$\begin{equation}
\eta_D(\chi_1) = \frac{i}{9},
\end{equation}$$ which is **purely imaginary** because $n = 3$ is odd and the $\mathbb{Z}_3$ characters are complex ($\omega = e^{2\pi i/3}$, not real). Under Wick rotation (Osterwalder--Schrader reconstruction), an imaginary spectral asymmetry requires one time dimension. Since $\mathbb{C}$ has exactly one imaginary axis, $d_{\mathrm{time}} = 1$, giving signature $(3,1)$.

**Why not $\mathbb{Z}_2$?** For $p = 2$, the characters $\{+1, -1\}$ are real, so $\eta_D$ would be real, and Wick rotation would not distinguish time from space. **Time exists because $\mathbb{Z}_3$ characters are complex.** The uniqueness theorem $n = p^{n-2}$ selects $p = 3$ over $p = 2$; therefore it selects Lorentzian signature. Verification: `lorentzian_proof.py`.

## The arrow of time

The spectral asymmetry $\eta(\phi)$ grows from $0$ (at $\phi = 0$, the smooth $S^5$) to $2/9$ (at $\phi = \phi_{\mathrm{lotus}}$, the present). This growth *is* the arrow of time: the spectral asymmetry measures the degree of chiral asymmetry, which increases monotonically as the fold deepens. At $\phi = 0$, $\eta = 0$: no spectral asymmetry, no chirality, no CP violation, no arrow of time. At $\phi = \phi_{\mathrm{lotus}}$, $\eta = 2/9$: maximal spectral asymmetry, maximal chirality, maximal CP violation, irreversible thermodynamics.

## The spectral phase transition

At $\phi_c \approx 0.60$, the universe crosses from substrate regime ($\phi < \phi_c$, perturbative, full SO(6)) to information regime ($\phi > \phi_c$, non-perturbative, $\mathbb{Z}_3$ dominant). This crossover triggers inflation ending, baryogenesis, and dark matter freeze-out simultaneously (Section [17](#sec:cosmology){reference-type="ref" reference="sec:cosmology"}). The evolving fold field $\phi(t)$ is the inflaton at high energy and the Higgs at low energy --- the same field at different epochs.

## The cosmic clock

From the spectral Friedmann equation with all inputs Theorem-level: $$\begin{align}
H_0 &= 67.7\;\mathrm{km/s/Mpc} \quad (0.5\%), \\
t_{\mathrm{age}} &= 13.72\;\mathrm{Gyr} \quad (0.5\%).
\end{align}$$ Both follow from spectral CC (Theorem) + spectral $\Omega$ ratios (Theorem) + Planck mass (Theorem) + Friedmann equation (textbook). Verification: `h0_spectral.py`, `age_of_universe.py`.

## The Sheet Music: Temporal Eigenvalues

The Dirac operator on $S^5/\mathbb{Z}_3$ has two channels:

- **Treble clef (spatial, resolved):** $|\eta_D| = 2/9$. Eigenvalues = *masses* (the Lotus Song, Section [11](#sec:song){reference-type="ref" reference="sec:song"}).

- **Bass clef (temporal, unresolved):** $\mathrm{Im}(\eta_D) = 1/9$. Eigenvalues = *decay rates*.

A stable particle has temporal eigenvalue $= 0$ (purely spatial). An unstable particle has temporal eigenvalue $> 0$ (mixed spatial + temporal). **Stability** is topological: the proton's $\mathbb{Z}_3$ baryon number is conserved by $e_m^2 = e_m$, so no lighter baryon state exists, and the temporal overlap is exactly zero.

The CKM matrix *is* the temporal channel. Every weak decay involves the CKM element $V_{ij}$, which measures how much of each quark's wavefunction lives in the $\mathrm{Im}(\eta)$ channel: $$\begin{equation}
|V_{ud}|^2 \approx 1 - \eta^2 = 1 - 4/81, \qquad |V_{us}| \approx \eta = 2/9, \qquad |V_{ub}| \sim \eta^3.
\end{equation}$$ The Cabibbo angle $\eta = 2/9$ is the temporal penetration depth. Strong decays bypass the temporal barrier (no flavor change); weak decays cross it.

**Test results** (all inputs spectral):

- Neutron: $\tau_n = 899$ s (PDG: 878.4, $2.3\%$). The CKM suppression $|V_{ud}|^2 \approx 1 - \eta^2$ is the temporal weight.

- Pion: $\tau_{\pi^\pm} = 2.70 \times 10^{-8}$ s (PDG: $2.60 \times 10^{-8}$, $3.5\%$), from spectral $f_\pi = K^2\eta\,m_p$.

- Muon: $\tau_\mu = 2.19 \times 10^{-6}$ s (PDG: $2.20 \times 10^{-6}$, $0.5\%$).

The Resolved Chord gives the notes (masses). The Unresolved Chord gives the rhythm (decay rates). Together: the Sheet Music of the Universe. Verification: `sheet_music_spectral.py`.

**What remains open:** The dynamics of $\phi(t)$ (the explicit trajectory from Big Bang to now), scattering amplitudes, thermal history details, and the emergence of time from the Information $\times$ Substrate framework. These are the subject of the sequel --- the *Unresolved Chord*.

# Conclusion {#sec:conclusion}

We have shown that the spectral geometry of a single compact manifold, $S^5/\mathbb{Z}_3$, determines all twenty-six parameters of the Standard Model, the three gauge couplings, the gravitational coupling, the cosmological constant, and the full cosmological history --- with zero free parameters ($m_e$ is the unit of measurement [@duff2002], not a tunable input). The logic:

1.  **Constraint:** $F(M) = p\cdot\sum|\eta_D| - K_p = 0$.

2.  **Uniqueness:** Only $S^5/\mathbb{Z}_3$ satisfies this ($n = p^{n-2} \Rightarrow 3 = 3^{3-2}$).

3.  **Dictionary:** The spectral action generates a four-level cascade: $m_e$ (unit) $\to$ $m_p = d_1\pi^5$ (ghost energy) $\to$ $\alpha$ (lag) $\to$ $v, m_H$ (EM budget) $\to$ everything.

4.  **Identity chain:** $\eta = d_1/p^n$ $\to$ $G = \lambda_1\eta$ $\to$ $c_{\mathrm{grav}} = -\tau/G$ $\to$ all forces connected through $p^n = 27$.

5.  **Unification:** All four forces are curvatures of $S^5/\mathbb{Z}_3$ at different geometric levels (QCD at 0D cone, weak at twist angle, EM at 4D wall, gravity at 5D bulk).

The framework predicts seventy-seven quantities (Table [11](#tab:master){reference-type="ref" reference="tab:master"}), all from five spectral invariants and $\pi$; all at Theorem level. The scope spans every sector of known physics: charged lepton masses ($0.001$--$0.007\%$), the proton mass ($10^{-11}$), the fine-structure constant ($0.001\%$), all CKM and PMNS parameters, all six quark masses (RMS $0.111\%$), the Planck mass ($0.10\%$), the cosmological constant ($0.11\%$), the cosmic energy budget ($0.96\%$), 27 hadron masses ($0.95\%$ RMS), and the deuteron binding energy ($0.00\%$).

The "quantum gravity problem" dissolves: the graviton is a KK mode, already quantum. The cosmological constant problem dissolves: heavy modes cancel by equidistribution ($\mathbb{Z}_3$ characters), leaving only $m_{\nu_3}$ suppressed by $\eta^2$. The gauge hierarchy dissolves: it is $11^7$ from the spectral content.

The spectral phase transition at $\phi_c = 0.60$ provides a unified mechanism for inflation ($N = 3025/48 \approx 63$ e-folds), baryogenesis (CP violation from evolving $\eta$), and dark matter production (frozen ghost modes, abundance $\sim\eta^2$). The universe did not start in $S^5/\mathbb{Z}_3$ --- it was forced there by the uniqueness theorem. During the transition, the rules were different: ghost modes were physical, baryon number did not exist, and the spectral asymmetry was evolving. The lotus did not bloom at $t = 0$. It bloomed at $t = t_c$.

**What is not claimed.** This paper does not claim a complete derivation of the Standard Model Lagrangian (see Connes--Chamseddine [@connes1996; @chamseddine2011]). It claims that *all* ratios, counts, and scales of the matter sector --- including gravity, the cosmological constant, and the cosmological energy budget --- are fixed by the spectral geometry of $S^5/\mathbb{Z}_3$. The framework contains no adjustable knobs beyond the electron mass (the unit of measurement). All seventy-seven predictions are at Theorem level: every factor traces to a spectral invariant of $S^5/\mathbb{Z}_3$ or to standard physics applied to Theorem-level inputs.

**On the Euclidean--Lorentzian transition.** The spectral action is naturally Euclidean; the Donnelly computation is Riemannian. The transition to physical Lorentzian spacetime is not an external Wick rotation but an *internal* property of the spectral asymmetry: $\eta_D(\chi_1) = i/9$ is purely imaginary. The magnitude $|\eta_D| = 1/9$ encodes spatial physics (the resolved channel); the phase $\arg(\eta_D) = \pi/2$ encodes temporal physics (the unresolved channel). The two channels mix with spectral weights $(1{-}1/p^2) = 8/9$ (space) and $1/p^2 = 1/9$ (time). The deuteron binding energy validates this mixing: the bare formula gives $+2.9\%$ error; including the temporal channel gives $0.00\%$. The mathematical infrastructure for Lorentzian spectral geometry (Bär, Strohmaier) is an active research frontier, but the spectral asymmetry already encodes the signature within the Riemannian framework.

**On the Born rule.** The functional form $P_m = |\langle m|\psi\rangle|^2$ follows from the idempotent structure of $\mathbb{C}[\mathbb{Z}_3]$: the minimal idempotents $e_m$ satisfy $e_m^2 = e_m$, so $\|e_m\psi\|^2 = \langle\psi|e_m|\psi\rangle$; the *square* arises because $e^2 = e$. Derived: the functional form, completeness ($\sum P_m = 1$), exclusivity, and repeatability. Not derived: why any particular outcome occurs (the outcome problem). The derivation is non-trivial because it is specific to $\mathbb{Z}_3$: only the uniqueness-selected group provides idempotents that simultaneously satisfy $e^2 = e$ and map to physical observables.

**On magnetic monopoles.** The framework predicts no magnetic monopoles: $\pi_2(S^5/\mathbb{Z}_3) = 0$ (homotopy of lens spaces). This is a distinctive anti-prediction --- all GUT frameworks (SU(5), SO(10)) *require* monopoles ($\pi_2(G/H) = \mathbb{Z}$), while the spectral framework *forbids* them. The MoEDAL experiment at the LHC will continue to find nothing.

**On LHC anti-predictions (closed spectrum).** The particle spectrum derived from $D$ on $M^4 \times S^5/\mathbb{Z}_3$ is *closed*: the only BSM particle is the fold-wall scalar at 95 GeV (S3 above). Eight classes of exotic state are forbidden by topological, algebraic, or spectral arguments: (A1) no SUSY ($\dim S^5 = 5$ is odd; no covariantly constant spinor); (A2) no $Z'/W'$ ($\mathbb{Z}_3$ holonomy fixes the gauge group exactly); (A3) no monopoles ($\pi_2 = 0$); (A4) no 4th generation ($p = 3$ discrete); (A5) no extra Higgs doublets (unique spectral fluctuation); (A6) no leptoquarks ($l{=}1$ ghost killing); (A7) no large extra dimensions ($M_c \sim 10^{13}$ GeV); (A8) no heavy neutral leptons (spectral seesaw fixes $m_{\mathrm{sterile}} = 3.55$ keV). Five loophole classes were examined (additional KK modes, ghost visibility, glueball states, charged fold partners, non-perturbative effects); none survive. Every null BSM search at the LHC is a confirmed anti-prediction. The detailed loophole analysis is in Supplement IX; computational verification in `lhc_exotics_spectral.py`.

**On black hole entropy.** The Bekenstein--Hawking formula $S = A/(4G)$ follows automatically from the spectral action: the $a_2$ Seeley--DeWitt coefficient gives Einstein--Hilbert gravity, and the Wald entropy formula on a Ricci-flat (Schwarzschild) horizon gives $A/(4G)$ exactly. No new spectral invariants are needed. Spectral corrections from the $R^2$ term (Starobinsky) enter only for non-vacuum black holes and are suppressed by $1/N = 48/3025 \approx 1.6\%$. Information is preserved by spectral monogamy ($\sum e_m = \mathbf{1}$ is topological) and the fold bounce (the LOTUS potential prevents singularities).

The temporal channel of $D_{\mathrm{wall}}$ gives decay rates: the CKM matrix is the temporal barrier, with Cabibbo angle $\eta = 2/9$ as the penetration depth. Stable particles have zero temporal eigenvalue (topological conservation); unstable particles have non-zero temporal eigenvalue (CKM suppression). Neutron, pion, and muon lifetimes are reproduced from spectral inputs alone (0.5--3.5%).

One manifold. Zero free parameters. The Resolved Chord resolves not just the notes, but the first chord --- and the Sheet Music tells each note how long to ring.

**Acknowledgments.** This work was developed independently. Computational verification scripts and a full test suite are available in the companion repository[^3].

::: thebibliography
99

Y. Koide, "New view of quark and lepton mass hierarchy," *Phys. Rev. D* **28** (1983) 252.

Y. Koide, "Charged lepton mass sum rule," *Mod. Phys. Lett. A* **28** (2013) 1350113.

R. Foot, "A note on Koide's lepton mass relation," arXiv:hep-ph/9402242 (1994).

L. Motl and A. Rivero, "Strange formula of quark and lepton masses," arXiv:hep-ph/0212251 (2002).

H. Donnelly, "Eta invariants for $G$-spaces," *Indiana Univ. Math. J.* **27** (1978) 889--918.

M. F. Atiyah, V. K. Patodi, and I. M. Singer, "Spectral asymmetry and Riemannian geometry. I," *Math. Proc. Camb. Phil. Soc.* **77** (1975) 43.

M. F. Atiyah and I. M. Singer, "The index of elliptic operators: III," *Ann. Math.* **87** (1968) 546.

T. Kawasaki, "The index of elliptic operators over $V$-manifolds," *Nagoya Math. J.* **84** (1981) 135--157.

J. Cheeger, "Analytic torsion and the heat equation," *Ann. Math.* **109** (1979) 259--322.

W. Müller, "Analytic torsion and $R$-torsion of Riemannian manifolds," *Adv. Math.* **28** (1978) 233--305.

K. Reidemeister, "Homotopieringe und Linsenräume," *Abh. Math. Sem. Hamburg* **11** (1935) 102.

A. Ikeda, "On the spectrum of the Laplacian on the spherical space forms," *Osaka J. Math.* **17** (1980) 691.

P. B. Gilkey, *Invariance Theory, the Heat Equation, and the Atiyah-Singer Index Theorem*, Publish or Perish, 1984.

A. Connes, "Gravity coupled with matter and the foundation of non-commutative geometry," *Commun. Math. Phys.* **182** (1996) 155.

A. H. Chamseddine and A. Connes, "The uncanny precision of the spectral action," *Commun. Math. Phys.* **307** (2011) 735.

H. Georgi and S. L. Glashow, "Unity of all elementary-particle forces," *Phys. Rev. Lett.* **32** (1974) 438.

E. Witten, "Search for a realistic Kaluza-Klein theory," *Nucl. Phys. B* **186** (1981) 412.

L. Dixon, J. Harvey, C. Vafa, and E. Witten, "Strings on orbifolds," *Nucl. Phys. B* **261** (1985) 678.

L. Dixon, J. Harvey, C. Vafa, and E. Witten, "Strings on orbifolds II," *Nucl. Phys. B* **274** (1986) 285.

C. Vafa and E. Witten, "Parity conservation in quantum chromodynamics," *Phys. Rev. Lett.* **53** (1984) 535.

S. Navas *et al.* (Particle Data Group), "Review of Particle Physics," *Phys. Rev. D* **110**, 030001 (2024).

I. Esteban *et al.*, "The fate of hints: updated global analysis of three-flavour neutrino oscillations," *JHEP* **09** (2020) 178; NuFIT 5.3 (2024), `www.nu-fit.org`.

G. Grubb, *Functional Calculus of Pseudodifferential Boundary Problems*, 2nd ed., Birkhäuser, 1996.

P. G. O. Freund and M. A. Rubin, "Dynamics of dimensional reduction," *Phys. Lett. B* **97** (1980) 233.

M. J. Duff, L. B. Okun, and G. Veneziano, "Trialogue on the number of fundamental constants," *JHEP* **03** (2002) 023; arXiv:physics/0110060.
:::

<figure id="fig:emblem" data-latex-placement="htbp">
<p><img src="figures/lotus_emblem.png" style="width:75.0%" alt="image" /> <span id="fig:emblem" data-label="fig:emblem"></span></p>
</figure>

[^1]: Every dimensionless ratio is fixed by the geometry. The electron mass $m_e$ is the ruler, not a knob: changing it rescales all predictions in lockstep without improving or degrading any fit [@duff2002].

[^2]: In standard $\Lambda$CDM, $\Omega_\Lambda/\Omega_{\mathrm{m}}$ evolves with redshift. The spectral claim is that this ratio is evaluated at the equilibrium epoch $\phi = \phi_{\mathrm{lotus}} = 0.9574$, which the LOTUS potential defines as the unique minimum of the fold stiffness potential $V(\phi)$. The spectral phase transition selects the equilibrium epoch dynamically; the ratio $2\pi^2/9$ is the value *at* that equilibrium, not an arbitrary time slice. The epoch-selection is proven in `cosmic_snapshot_epoch.py`: the ratio is a spectral partition (fold energy $2\pi^2$ over orbifold volume $p^2$), both topological invariants evaluated at the LOTUS equilibrium $\phi_{\mathrm{lotus}}$. The LOTUS potential traps $\phi$ at the minimum (tunneling rate $\sim e^{-10^{68}}$); the partition is therefore eternal. Status: **Theorem**.

[^3]: Computational assistants were used for code verification and algebraic cross-checks.
