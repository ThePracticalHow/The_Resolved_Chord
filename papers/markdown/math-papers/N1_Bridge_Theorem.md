# Introduction {#sec:intro}

A recurring situation in spectral geometry and mathematical physics: an operator $D$ acts on a Hilbert space $\mathcal{H}$, a finite group $G$ acts by symmetries, and the spectrum decomposes into sectors labeled by the irreducible representations of $G$. A spectral quantity --- a trace, a zeta value, a heat coefficient, an eta invariant --- is computed for each sector. The question arises:

::: center
*Does the per-sector quantity depend on the regularization scheme?*
:::

If it does, the decomposition is an artifact; if it does not, it is a geometric invariant. The purpose of this paper is to provide a single, general theorem that settles this question for all finite isometric group actions, and to demonstrate its reach through a catalog of applications.

The theorem itself is elementary --- three steps of functional calculus --- and experts in operator algebra may regard it as "well known." However, we argue that its value lies not in the proof but in the *generality of application*. Like a wrench, the tool is simple; the art is knowing which bolts it fits. We exhibit five such bolts from different areas of mathematics and physics, and one explicit counterexample showing that the result has content: it fails when the isometry condition is violated.

# Setup and Notation {#sec:setup}

::: {#def:setup .definition}
**Definition 1** (Spectral triple with finite group action). *Let $(\mathcal{H}, D)$ be a spectral datum consisting of:*

- *A separable Hilbert space $\mathcal{H}$.*

- *A self-adjoint operator $D$ with compact resolvent (hence discrete spectrum $\{\lambda_n\}_{n \in \mathbb{N}}$ with $|\lambda_n| \to \infty$).*

*Let $G$ be a finite group acting on $\mathcal{H}$ by unitary operators $\rho: G \to \mathcal{U}(\mathcal{H})$ satisfying the **isometry condition**: $$\begin{equation}
\rho(g)\, D = D\, \rho(g) \qquad \text{for all } g \in G.
\label{eq:commutation}
\end{equation}$$*
:::

::: {#def:idempotents .definition}
**Definition 2** (Minimal central idempotents). *The group algebra $\mathbb{C}[G]$ decomposes as a direct sum of matrix algebras, one for each irreducible representation $\pi_m$ of $G$ ($m = 0, 1, \ldots, r-1$ where $r$ is the number of irreducible representations). The **minimal central idempotents** are: $$\begin{equation}
e_m = \frac{\dim \pi_m}{|G|} \sum_{g \in G} \overline{\chi_m(g)}\, \rho(g),
\label{eq:idempotent}
\end{equation}$$ where $\chi_m$ is the character of $\pi_m$. These satisfy: $$\begin{align}
e_m^2 &= e_m, \label{eq:idem-square}\\
e_m e_{m'} &= 0 \quad\text{for } m \neq m', \label{eq:idem-orthogonal}\\
\sum_{m=0}^{r-1} e_m &= \mathbf{1}. \label{eq:partition-of-unity}
\end{align}$$*
:::

::: definition
**Definition 3** (Sector Hilbert space and restricted operator). *The $m$-th sector Hilbert space is $\mathcal{H}_m = e_m \mathcal{H}$. The restricted operator is $D_m = D|_{\mathcal{H}_m}$.*
:::

# The Main Theorem {#sec:main}

::: {#thm:N1 .theorem}
**Theorem 4** (Coefficient One). *Let $(\mathcal{H}, D, G, \rho)$ be as in Definition [1](#def:setup){reference-type="ref" reference="def:setup"}, satisfying the isometry condition [\[eq:commutation\]](#eq:commutation){reference-type="eqref" reference="eq:commutation"}. Let $f: \mathbb{R} \to \mathbb{C}$ be a Borel function such that $f(D)$ is trace-class[^1]. Then for each minimal central idempotent $e_m$: $$\begin{equation}
\boxed{\mathrm{Tr}(f(D)\, e_m) = \mathrm{Tr}(f(D_m)).}
\label{eq:main}
\end{equation}$$ In particular, the coefficient of every per-sector spectral quantity (eta invariant, heat coefficient, zeta value) is exactly $1$, independent of $f$.*
:::

::: proof
*Proof.* The proof proceeds in three steps.

**Step 1: $f(D)$ commutes with $e_m$.**

Since $G$ acts by isometries, $\rho(g) D = D \rho(g)$ for all $g \in G$. By the functional calculus for self-adjoint operators, for any bounded Borel function $f$: $$\begin{equation}
\rho(g)\, f(D) = f(D)\, \rho(g) \qquad \text{for all } g \in G.
\label{eq:f-commutes-g}
\end{equation}$$ *Proof of [\[eq:f-commutes-g\]](#eq:f-commutes-g){reference-type="eqref" reference="eq:f-commutes-g"}:* Let $D = \int \lambda\, dE(\lambda)$ be the spectral decomposition, so $f(D) = \int f(\lambda)\, dE(\lambda)$. Since $\rho(g)$ commutes with $D$, it commutes with every spectral projection $E(B)$ for Borel sets $B \subset \mathbb{R}$ (standard result: see [@reed-simon], Theorem VIII.5). Therefore $\rho(g)$ commutes with $\int f(\lambda)\, dE(\lambda) = f(D)$.

Since $e_m$ is a $\mathbb{C}$-linear combination of the $\rho(g)$ (equation [\[eq:idempotent\]](#eq:idempotent){reference-type="eqref" reference="eq:idempotent"}), it follows that: $$\begin{equation}
[f(D),\, e_m] = 0.
\label{eq:f-commutes-em}
\end{equation}$$

**Step 2: The trace decomposes.**

Since $e_m$ is a projection ($e_m^2 = e_m$) and commutes with $f(D)$, the operator $f(D)\, e_m = e_m\, f(D)\, e_m$ acts on $\mathcal{H}_m = e_m\mathcal{H}$. Let $\{|\psi_n^{(m)}\rangle\}$ be an orthonormal eigenbasis of $D_m$ in $\mathcal{H}_m$ with eigenvalues $\{\lambda_n^{(m)}\}$. Then: $$\begin{align}
\mathrm{Tr}(f(D)\, e_m)
&= \sum_n \langle \psi_n^{(m)} | f(D)\, e_m | \psi_n^{(m)} \rangle
\label{eq:trace-expand}\\
&= \sum_n \langle \psi_n^{(m)} | f(D) | \psi_n^{(m)} \rangle
\qquad \text{(since $e_m |\psi_n^{(m)}\rangle = |\psi_n^{(m)}\rangle$)}
\label{eq:em-acts}\\
&= \sum_n f(\lambda_n^{(m)})
\qquad \text{(since $|\psi_n^{(m)}\rangle$ is an eigenstate of $D$)}
\label{eq:eigenvalue}\\
&= \mathrm{Tr}(f(D_m)).
\label{eq:final}
\end{align}$$

**Step 3: The coefficient is $1$.**

Equation [\[eq:final\]](#eq:final){reference-type="eqref" reference="eq:final"} shows $\mathrm{Tr}(f(D)\, e_m) = \mathrm{Tr}(f(D_m))$ with coefficient exactly $1$, for any trace-class $f(D)$. No normalization, no $f$-dependent prefactor. 0◻ ◻
:::

# Why the Result Is Not Trivial: A Counterexample {#sec:counterexample}

One might suspect that Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} is vacuous --- that "the coefficient is always $1$" regardless of assumptions. We show this is false: when the isometry condition [\[eq:commutation\]](#eq:commutation){reference-type="eqref" reference="eq:commutation"} is violated, the coefficient becomes $f$-dependent.

::: {#ex:counter .example}
**Example 5** (Non-isometric $\mathbb{Z}_2$ action on $S^1$). *Let $\mathcal{H} = L^2(S^1)$ with the standard Dirac operator $D_0 = -i\,d/d\theta$, spectrum $\{n : n \in \mathbb{Z}\}$. Let $\mathbb{Z}_2$ act by the reflection $\theta \mapsto -\theta$, which is an isometry of the round $S^1$. In this case, $[\rho(\sigma), D_0] = 0$ and Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} applies: all spectral quantities decompose with coefficient $1$.*

*Now deform the metric: let $D_\epsilon = -i\,h(\theta)^{-1}\, d/d\theta$ where $h(\theta) = 1 + \epsilon \cos\theta$ ($|\epsilon| < 1$). The reflection $\theta \mapsto -\theta$ preserves $h$ (since $h(-\theta) = h(\theta)$), so it is still an isometry and $[\rho(\sigma), D_\epsilon] = 0$. Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} still applies.*

*But consider instead the *non-isometric* action $\theta \mapsto \theta + \pi$ (translation by half-period) with the *same* deformed metric $h(\theta) = 1 + \epsilon\cos\theta$. Since $h(\theta + \pi) = 1 - \epsilon\cos\theta \neq h(\theta)$ for $\epsilon \neq 0$, this action does *not* preserve the metric. Therefore $[\rho(\sigma), D_\epsilon] \neq 0$, and: $$\begin{equation}
\mathrm{Tr}(f(D_\epsilon)\, e_{\mathrm{even}}) \neq \mathrm{Tr}(f((D_\epsilon)_{\mathrm{even}}))
\quad\text{for generic } f.
\end{equation}$$ Explicitly: the eigenvalues of $D_\epsilon$ are $\mu_n = n + \epsilon\, c_n + O(\epsilon^2)$ where the perturbation coefficients $c_n$ depend on the parity of $n$ asymmetrically. The even-sector trace $\sum_{n\,\mathrm{even}} f(\mu_n)$ differs from $\mathrm{Tr}(f(D_\epsilon)\, e_{\mathrm{even}})$ by $O(\epsilon)$ corrections whose sign and magnitude depend on $f$. The "coefficient" is no longer $1$ but a function of the cutoff.*
:::

::: remark
**Remark 6** (The isometry condition has teeth). *The counterexample shows that Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} is genuinely a theorem about isometric group actions. The word "isometry" cannot be weakened to "diffeomorphism," "conformal map," or "homeomorphism." The wrench fits only isometric bolts.*
:::

# Application 1: Spectral Action on Orbifolds {#sec:spectral-action}

::: {#cor:spectral-action .corollary}
**Corollary 7** (Cutoff independence of sector corrections). *In the Connes--Chamseddine spectral action [@connes1996; @chamseddine2011], the bosonic action on a Riemannian orbifold $M/G$ is $S = \mathrm{Tr}(f(D/\Lambda))$ where $f$ is a cutoff function and $\Lambda$ a scale. The sector decomposition $$\begin{equation}
S = \sum_{m=0}^{r-1} S_m, \qquad S_m = \mathrm{Tr}(f(D_m/\Lambda)),
\end{equation}$$ holds with coefficient $1$ per sector, independent of the choice of $f$.*
:::

This is the original motivation for the theorem. Applied to $S^5/\mathbb{Z}_3$: the three $\mathbb{Z}_3$ sectors ($\chi_0$, $\chi_1$, $\chi_2$) contribute their respective eta invariants $\eta_D(\chi_m)$ to the Yukawa coupling phase with coefficient exactly $1$, regardless of whether $f$ is a sharp cutoff, a smooth Schwartz function, or a heat kernel. This eliminates a potential source of scheme dependence in the derivation of the Koide phase $\delta = 2\pi/3 + 2/9$ from the spectral geometry of $S^5/\mathbb{Z}_3$ [@leng2026].

::: remark
**Remark 8**. *The Seeley--DeWitt expansion $\mathrm{Tr}(f(D/\Lambda)) \sim \sum_k f_k \Lambda^{d-k} a_k(D)$ involves "moments" $f_k = \int_0^\infty f(u)\, u^{k/2-1}\, du$ that *do* depend on $f$. The theorem does not say these moments are universal --- it says the *sector decomposition* $a_k = \sum_m a_k^{(m)}$ is universal. The moments multiply the total; the decomposition multiplies with coefficient $1$.*
:::

# Application 2: Heat Kernel Asymptotics {#sec:heat}

::: {#cor:heat .corollary}
**Corollary 9** (Heat trace decomposition). *Let $D$ be the Dirac operator on a compact Riemannian orbifold $M/G$. The heat trace $$\begin{equation}
K(t) = \mathrm{Tr}(e^{-tD^2}) = \sum_{n} e^{-t\lambda_n^2}
\end{equation}$$ decomposes by $G$-representation sector as: $$\begin{equation}
K(t) = \sum_{m=0}^{r-1} K_m(t), \qquad K_m(t) = \mathrm{Tr}(e^{-tD_m^2}),
\end{equation}$$ with coefficient $1$ per sector, for all $t > 0$.*
:::

::: proof
*Proof.* Apply Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} with $f(x) = e^{-tx^2}$, which is a bounded Borel function making $f(D)$ trace-class (since $D$ has compact resolvent and $e^{-tx^2}$ decays rapidly). 0◻ ◻
:::

**Consequence for Seeley--DeWitt coefficients.** The small-$t$ expansion $K_m(t) \sim \sum_{k \geq 0} a_k^{(m)}\, t^{(k-d)/2}$ gives Seeley--DeWitt coefficients $a_k^{(m)}$ for each sector. These are the building blocks of one-loop effective actions on orbifold backgrounds. Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} guarantees: $$\begin{equation}
a_k = \sum_{m=0}^{r-1} a_k^{(m)},
\end{equation}$$ with no renormalization of the per-sector contributions. In particular, the equivariant Euler characteristic, the equivariant signature, and the equivariant $\hat{A}$-genus all decompose with coefficient $1$.

::: remark
**Remark 10** (One-loop determinants in string theory). *On orbifold string backgrounds $\mathcal{M}/G$, the one-loop partition function $Z = (\det D^2)^{-1/2}$ factorizes by twisted sector. The coefficient $1$ theorem ensures that no sector receives an anomalous weight --- the twisted-sector contributions to the vacuum energy are exactly $\log\det(D_m^2)$ with no cutoff-dependent normalization. This is implicit in standard orbifold CFT calculations [@dixon1985] but is not usually stated as a general theorem.*
:::

# Application 3: Casimir Energy on Orbifolds {#sec:casimir}

::: {#cor:casimir .corollary}
**Corollary 11** (Casimir energy decomposition). *The regularized vacuum energy of a quantum field on $M/G$ is $$\begin{equation}
E_{\mathrm{Cas}} = \frac{1}{2}\zeta_D'(0), \qquad
\zeta_D(s) = \sum_{\lambda_n > 0} \lambda_n^{-2s},
\end{equation}$$ where the sum is over positive eigenvalues of $D$. This decomposes as: $$\begin{equation}
E_{\mathrm{Cas}} = \sum_{m=0}^{r-1} E_{\mathrm{Cas}}^{(m)}, \qquad
E_{\mathrm{Cas}}^{(m)} = \frac{1}{2}(\zeta_{D_m})'(0),
\end{equation}$$ with coefficient $1$ per sector, independent of the regularization prescription.*
:::

::: proof
*Proof.* The zeta function $\zeta_D(s) = \mathrm{Tr}(D^{-2s}\, \Pi_+)$ (where $\Pi_+$ projects onto positive eigenvalues) decomposes by sector via Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} applied to $f(x) = |x|^{-2s}$ in the region of absolute convergence $\mathrm{Re}(s) > d/2$, then extended to $s = 0$ by analytic continuation. The continuation preserves the coefficient because it acts independently in each sector (the sectors are spectrally disjoint). 0◻ ◻
:::

::: remark
**Remark 12** (Randall--Sundrum and orbifold GUTs). *In Randall--Sundrum models [@randall1999] and orbifold GUT compactifications [@hall2002], the Casimir energy on $S^1/\mathbb{Z}_2$ stabilizes the extra dimension. The $\mathbb{Z}_2$-even and $\mathbb{Z}_2$-odd sectors contribute independently to the Casimir force, and the coefficient $1$ theorem guarantees that no regularization artifact contaminates the even/odd decomposition. This is physically important: the hierarchy between the Planck and TeV scales depends on the *ratio* of even to odd Casimir contributions, and a scheme-dependent coefficient would destroy the prediction.*
:::

# Application 4: Equivariant Index Theory {#sec:index}

::: {#cor:index .corollary}
**Corollary 13** (Equivariant APS index). *Let $(M, \partial M)$ be a compact Riemannian manifold with boundary, $G$ a finite group of isometries, and $D$ the Dirac operator with APS boundary conditions. The equivariant index in the $m$-th sector is: $$\begin{equation}
\mathrm{ind}_m(D) = \mathrm{Tr}(\gamma_5\, e_m) = \mathrm{Tr}(\gamma_5|_{\mathcal{H}_m}),
\end{equation}$$ with coefficient $1$, independent of the regularization used to define the index.*
:::

::: proof
*Proof.* The index can be expressed as $\mathrm{ind}_m = \mathrm{Tr}(\gamma_5\, e^{-tD^2}\, e_m)$ for any $t > 0$ (McKean--Singer formula). Since $G$ acts by isometries, $\gamma_5$ commutes with the $G$-action (the chirality grading is preserved by orientation-preserving isometries). By Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"}, the trace decomposes with coefficient $1$. Taking $t \to 0^+$ gives the index. 0◻ ◻
:::

**Connection to generation counting.** In the companion paper [@leng2026aps], the equivariant APS index on $B^6/\mathbb{Z}_3$ is computed to give $N_g = 1 + 1 + 1 = 3$ --- one chiral zero mode per $\mathbb{Z}_3$ sector. The coefficient $1$ theorem is the reason this count does not depend on the choice of heat-kernel regulator $t$ or any other regularization parameter. The generation count $N_g = 3$ is a topological invariant precisely *because* the equivariant decomposition has coefficient $1$.

Without the isometry condition, one could imagine an anomalous weighting where sector 1 receives weight $1 + \epsilon$ and sector 2 receives weight $1 - \epsilon$ for some $\epsilon$ depending on the regularization. The theorem rules this out.

# Application 5: Crystallographic Point Groups {#sec:crystal}

::: {#cor:crystal .corollary}
**Corollary 14** (Band structure decomposition). *Let $\mathcal{H} = L^2(\mathbb{R}^3/\Lambda)$ be the Hilbert space of a crystalline solid with lattice $\Lambda$, and let $G$ be the point group of the crystal (a finite subgroup of $\mathrm{O}(3)$). The Hamiltonian $H = -\nabla^2 + V(x)$, where $V$ is $G$-invariant, satisfies $[\rho(g), H] = 0$. The density of states per symmetry channel decomposes as: $$\begin{equation}
\rho(\epsilon) = \sum_{m} \rho_m(\epsilon), \qquad
\rho_m(\epsilon) = \mathrm{Tr}(\delta(\epsilon - H_m)),
\end{equation}$$ with coefficient $1$, independent of any smearing or broadening prescription.*
:::

::: proof
*Proof.* Apply Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} with $f(x) = \delta(\epsilon - x^2)$ (distributional limit of smooth approximations), or equivalently with $f(x) = \chi_{[\epsilon, \epsilon+d\epsilon]}(x^2)$. In practice, one uses a Lorentzian or Gaussian broadening $f_\sigma(x) = (\pi\sigma)^{-1}(\sigma^2 + (x^2-\epsilon)^2)^{-1}$; the theorem says the per-irrep decomposition is independent of $\sigma$. 0◻ ◻
:::

::: remark
**Remark 15** (Tight-binding models). *In tight-binding models on molecules or clusters with point-group symmetry $G$, the molecular orbitals classify by irreps of $G$. The total energy $E_{\mathrm{tot}} = \sum_m E_m$ decomposes by irrep with coefficient $1$. This is exploited routinely in computational chemistry (e.g., symmetry-adapted basis sets), but the underlying mathematical guarantee is precisely Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"}. When symmetry-breaking perturbations are introduced (e.g., Jahn--Teller distortions), the isometry condition [\[eq:commutation\]](#eq:commutation){reference-type="eqref" reference="eq:commutation"} is violated and the clean decomposition is lost --- consistent with the counterexample in §[4](#sec:counterexample){reference-type="ref" reference="sec:counterexample"}.*
:::

# Non-Examples and Obstructions {#sec:obstructions}

The hypotheses of Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} --- finite group, compact resolvent, isometry condition --- are sharp. We catalog the failure modes.

::: {#nonex:infinite .nonexample}
**Non-Example 16** (Non-finite groups). *Let $G = \mathrm{U}(1)$ act on $L^2(S^1)$ by rotation. The group algebra $\mathbb{C}[\mathrm{U}(1)]$ is infinite-dimensional and has no minimal central idempotents in the algebraic sense. The decomposition into Fourier modes (irreps of $\mathrm{U}(1)$) still works spectrally, but the idempotent construction [\[eq:idempotent\]](#eq:idempotent){reference-type="eqref" reference="eq:idempotent"} involves an integral over $G$ rather than a finite sum. Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} extends to compact groups via the Peter--Weyl theorem, but the proof requires the additional hypothesis that the multiplicity spaces are finite-dimensional (automatic for compact $G$ on compact $M$, but subtle for non-compact $M$).*
:::

::: {#nonex:continuous .nonexample}
**Non-Example 17** (Continuous spectrum). *Let $D = -i\,d/dx$ on $L^2(\mathbb{R})$ (continuous spectrum, no compact resolvent). Let $\mathbb{Z}_2$ act by $x \mapsto -x$. The even/odd decomposition of $L^2(\mathbb{R})$ is well-defined, but $\mathrm{Tr}(f(D)\, e_{\mathrm{even}})$ is not defined for most $f$ because $f(D)$ is not trace-class. The theorem requires compact resolvent precisely to ensure trace-class regularity.*
:::

::: {#nonex:anti .nonexample}
**Non-Example 18** (Anti-unitary actions). *Time reversal $T$ is anti-unitary: $T(c|\psi\rangle) = \bar{c}\, T|\psi\rangle$. The idempotent construction [\[eq:idempotent\]](#eq:idempotent){reference-type="eqref" reference="eq:idempotent"} uses $\mathbb{C}$-linear combinations of group elements, which fails for anti-unitary operators. For systems with time-reversal symmetry, the relevant decomposition is by *real* or *quaternionic* representations (Dyson's threefold way), and the coefficient $1$ result must be replaced by the appropriate Kramers multiplicity.*
:::

::: {#nonex:diffeo .nonexample}
**Non-Example 19** (Non-isometric diffeomorphisms). *A diffeomorphism $\phi: M \to M$ that is not an isometry does not commute with the Laplacian or Dirac operator. The pullback $\phi^*$ acts on functions but $[\phi^*, \Delta] \neq 0$ in general (the Laplacian depends on the metric, which $\phi$ does not preserve). The counterexample in §[4](#sec:counterexample){reference-type="ref" reference="sec:counterexample"} is a concrete instance. In such cases, the "coefficient" in the sector decomposition becomes a function of the regularization, and per-sector spectral quantities are not intrinsic.*
:::

# Connection to Equivariant K-Theory {#sec:ktheory}

The decomposition $\mathcal{H} = \bigoplus_{m} e_m\mathcal{H}$ is a manifestation of the equivariant K-theory of the algebra of observables.

::: {#prop:ktheory .proposition}
**Proposition 20** (K-theoretic interpretation). *Let $\mathcal{A} = C(M)$ (continuous functions on $M$) with the $G$-action by pullback. The equivariant K-group $K_G^0(M)$ decomposes as: $$\begin{equation}
K_G^0(M) \cong \bigoplus_{m=0}^{r-1} K^0(M/G) \otimes R(G)_m,
\end{equation}$$ where $R(G)_m$ is the $m$-th component of the representation ring. The Chern character $\mathrm{ch}: K_G^0(M) \to H_G^*(M; \mathbb{Q})$ intertwines the idempotent decomposition: $\mathrm{ch}(e_m \cdot [D]) = e_m \cdot \mathrm{ch}([D])$. Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} is the *operator-trace shadow* of this K-theoretic decomposition.*
:::

::: remark
**Remark 21** (Baum--Connes for finite groups). *For finite groups, the Baum--Connes assembly map $\mu: K_*^G(\underline{E}G) \to K_*(C^*_r G)$ is an isomorphism [@baum-connes1994]. The coefficient $1$ in Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} reflects the fact that the assembly map for finite groups is an *exact* isomorphism --- no correction terms, no anomalous dimensions, no renormalization. For infinite discrete groups, the assembly map may fail to be surjective (the Baum--Connes conjecture), and the clean coefficient $1$ decomposition would not hold in general. The finiteness of $G$ is therefore not merely a technical convenience but a reflection of the exactness of the assembly map.*
:::

# Generality of the Result {#sec:general}

::: remark
**Remark 22** (Summary of scope). *Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"} holds for:*

1.  *Any finite group $G$ (not just cyclic groups).*

2.  *Any self-adjoint operator with compact resolvent (not just the Dirac operator).*

3.  *Any bounded Borel function $f$ (not just smooth cutoffs).*

4.  *Any Riemannian manifold $M$ on which $G$ acts by isometries.*

*The only requirement is the isometry condition [\[eq:commutation\]](#eq:commutation){reference-type="eqref" reference="eq:commutation"}: the group action must commute with the operator. For the Dirac operator on a Riemannian manifold, this is equivalent to requiring that $G$ act by isometries --- a geometric condition, not an analytic one.*
:::

# Summary {#sec:summary}

The coefficient one theorem (Theorem [4](#thm:N1){reference-type="ref" reference="thm:N1"}) is a three-line proof with five applications and one counterexample:

::: center
  **Application**                                                                                   **Eliminates**                                    **Section**
  ------------------------------------------------------------------------------------------------- ------------------------------------------------- ----------------------------------------------------------------------------------
  Spectral action (orbifolds)                                                                       Cutoff dependence in sector corrections           §[5](#sec:spectral-action){reference-type="ref" reference="sec:spectral-action"}
  Heat kernel asymptotics                                                                           Ambiguity in per-sector Seeley--DeWitt coeff.     §[6](#sec:heat){reference-type="ref" reference="sec:heat"}
  Casimir energy                                                                                    Regularization artifact in sector forces          §[7](#sec:casimir){reference-type="ref" reference="sec:casimir"}
  Equivariant APS index                                                                             Regulator dependence in generation counting       §[8](#sec:index){reference-type="ref" reference="sec:index"}
  Band structure (crystals)                                                                         Broadening dependence in per-irrep DOS            §[9](#sec:crystal){reference-type="ref" reference="sec:crystal"}
  Counterexample (§[4](#sec:counterexample){reference-type="ref" reference="sec:counterexample"})   Shows isometry condition is necessary             §[4](#sec:counterexample){reference-type="ref" reference="sec:counterexample"}
  Non-examples (§[10](#sec:obstructions){reference-type="ref" reference="sec:obstructions"})        Infinite $G$, continuous spectrum, anti-unitary   §[10](#sec:obstructions){reference-type="ref" reference="sec:obstructions"}
  K-theory (§[11](#sec:ktheory){reference-type="ref" reference="sec:ktheory"})                      Connects to Baum--Connes assembly map             §[11](#sec:ktheory){reference-type="ref" reference="sec:ktheory"}
:::

The theorem is the wrench. The applications are the bolts. The tool is simple; the art is knowing where it fits.

::: thebibliography
99

A. Connes, "Gravity coupled with matter and the foundation of non-commutative geometry," *Commun. Math. Phys.* **182** (1996) 155.

A. H. Chamseddine and A. Connes, "The uncanny precision of the spectral action," *Commun. Math. Phys.* **307** (2011) 735.

H. Donnelly, "Eta invariants for $G$-spaces," *Indiana Univ. Math. J.* **27** (1978) 889--918.

M. Reed and B. Simon, *Methods of Modern Mathematical Physics*, Vol. I: Functional Analysis, Academic Press, 1980.

L. Dixon, J. Harvey, C. Vafa, and E. Witten, "Strings on orbifolds," *Nucl. Phys. B* **261** (1985) 678--686.

L. Randall and R. Sundrum, "A large mass hierarchy from a small extra dimension," *Phys. Rev. Lett.* **83** (1999) 3370.

L. J. Hall and Y. Nomura, "Gauge unification in higher dimensions," *Phys. Rev. D* **64** (2001) 055003.

P. Baum, A. Connes, and N. Higson, "Classifying space for proper actions and $K$-theory of group $C^*$-algebras," *Contemporary Mathematics* **167** (1994) 241--291.

J. Leng, "The Resolved Chord: The Theorem of Everything," v10 (2026).

J. Leng, "Equivariant APS index on $B^6/\mathbb{Z}_3$ and the emergence of three generations," (2026).

J. Leng, "Eta invariants, Reidemeister torsion, and a ghost-mode identity on the lens space $L(3;1,1,1)$," (2026).
:::

[^1]: *This holds for heat kernels $f(x) = e^{-tx^2}$, resolvents $f(x) = (x^2 + m^2)^{-s}$ with $\mathrm{Re}(s)$ sufficiently large, smooth cutoffs with sufficient decay, and more generally any $f$ for which $\sum_n |f(\lambda_n)| < \infty$. The eta function case $f(x) = \mathrm{sign}(x)|x|^{-s}$ is handled by analytic continuation.*
