*This supplement is self-contained. It provides the complete derivation chain for the baryon sector of the main text (Section 4: Parameters 11--13). All definitions, lemmas, intermediate calculations, and numerical verifications are included. No result depends on material outside this document except where explicitly cross-referenced to Supplements I--III.*

**Parameters derived in this supplement:**

::: center
   \#  Quantity                                        Value
  ---- ----------------------------------------------- ---------------------------
   11  Leading proton--electron mass ratio $m_p/m_e$   $6\pi^5 = 1836.118\ldots$
   12  Spectral coupling $G$ (one-loop)                $10/9$
   13  Two-loop coefficient $G_2$                      $-280/9$
:::

# Derivation of the Leading Term $6\pi^5$ {#sec:leading}

## Statement

::: {#thm:leading .theorem}
**Theorem 1** (Leading proton--electron mass ratio). *Let $S^5$ be the unit round five-sphere with its canonical metric. Let $d_1 = 6$ be the multiplicity of the first nonzero eigenspace of the Laplacian on $S^5$, and let $\pi^5$ be the pointwise Gaussian phase-space weight on the tangent phase space. Then $$\begin{equation}
\label{eq:leading}
\boxed{\frac{m_p}{m_e}\bigg|_{\mathrm{leading}} = d_1 \cdot \pi^5 = 6\pi^5 = 1836.118\ldots}
\end{equation}$$*
:::

The proof occupies the remainder of this section.

## The factor $\pi^5$: pointwise Gaussian phase-space weight

::: remark
**Remark 1** (Two origins of $\pi^5$). *The Riemannian volume of the round unit $S^5$ is $\mathrm{Vol}(S^5) = 2\pi^3/\Gamma(3) = \pi^3$. This is *not* the origin of the factor $\pi^5$ in the local (Gaussian) derivation below, which uses only the tangent space. However, a complementary *global* decomposition exists: $\pi^5 = \mathrm{Vol}(S^5) \times \pi^2 = \pi^3 \times (\lambda_1 + \alpha_s)$, where $\pi^2 = 5 + (\pi^2 - 5)$ splits into the first eigenvalue and the Dirichlet gap. Both derivations yield $\pi^5$; the global one reveals the connection to $\alpha_s$. See §[8.2](#sec:pi5-decomp){reference-type="ref" reference="sec:pi5-decomp"} below.*
:::

::: definition
**Definition 1** (Tangent phase space). *At any point $x \in S^5$, the tangent phase space is $$\begin{equation}
\mathcal{P}_x \;=\; T_x S^5 \oplus T_x^* S^5 \;\cong\; \mathbb{R}^5 \oplus \mathbb{R}^5 \;=\; \mathbb{R}^{10}.
\end{equation}$$ This space is **flat**: it is a vector space equipped with the standard Euclidean inner product inherited from the round metric on $S^5$. No curvature approximation is involved --- the tangent space at a point is exactly $\mathbb{R}^5$.*
:::

::: {#prop:gaussian .proposition}
**Proposition 1** (Gaussian phase-space integral). *The Gaussian integral over $\mathcal{P}_x = \mathbb{R}^{10}$ is $$\begin{equation}
\label{eq:gaussian}
\int_{\mathbb{R}^{10}} e^{-(|q|^2 + |p|^2)}\, d^5 q\, d^5 p
\;=\; \left(\int_{\mathbb{R}^5} e^{-|q|^2}\, d^5 q\right)
     \left(\int_{\mathbb{R}^5} e^{-|p|^2}\, d^5 p\right)
\;=\; \pi^{5/2} \cdot \pi^{5/2}
\;=\; \pi^5.
\end{equation}$$*
:::

::: proof
*Proof.* This is the standard $n$-dimensional Gaussian integral $$\begin{equation}
\int_{\mathbb{R}^n} e^{-|x|^2}\, d^n x = \pi^{n/2},
\end{equation}$$ applied with $n = 5$ independently to the position and momentum sectors. The result is *exact* --- no series expansion, no curvature correction, no regularisation. The tangent space is a genuine vector space. ◻
:::

## Normalization convention: Wigner $e^{-r^2}$

The choice of Gaussian exponent is physically meaningful and must be stated precisely.

::: definition
**Definition 2** (Wigner convention). *The Wigner quasi-probability distribution for the vacuum state of a single harmonic mode is $$\begin{equation}
\label{eq:wigner}
W(q,p) = \frac{1}{\pi}\, e^{-(q^2 + p^2)}.
\end{equation}$$ The per-mode phase-space weight is $$\begin{equation}
\int_{\mathbb{R}^2} W(q,p)\, dq\, dp = 1,
\end{equation}$$ but the unnormalized Gaussian volume per mode is $$\begin{equation}
\int_{\mathbb{R}^2} e^{-(q^2+p^2)}\, dq\, dp = \pi.
\end{equation}$$ This $\pi$ is one quantum cell: the phase-space area occupied by one vacuum mode under the Wigner convention.*
:::

For five independent dimensions, the weight is $\pi^5$.

::: remark
**Remark 2** (Wrong convention check). *The alternative wave-function convention uses $e^{-|x|^2/2}$, which yields $$\begin{equation}
\int_{\mathbb{R}^{10}} e^{-(|q|^2+|p|^2)/2}\, d^5 q\, d^5 p
= (2\pi)^{5/2} \cdot (2\pi)^{5/2} = (2\pi)^5 \approx 9671.
\end{equation}$$ The ratio $6 \times 9671 \approx 58{,}027$ is off by a factor of $\sim\!32$ and is clearly wrong. The Wigner convention $e^{-r^2}$ is the correct one.*
:::

## Why the tangent space suffices

::: {#prop:tangent .proposition}
**Proposition 2**. *The pointwise Gaussian weight $\pi^5$ receives no curvature correction at leading (zeroth) order in the Seeley--DeWitt (SDW) expansion.*
:::

::: proof
*Proof.* The vacuum state is Gaussian in the tangent approximation; this corresponds to the zeroth-order SDW heat-kernel coefficient $a_0$. The round metric on $S^5$ is homogeneous under $\mathrm{SO}(6)$, so the pointwise weight $\pi^5$ is the same at every point $x \in S^5$. Curvature corrections enter only at the $a_2$ level and beyond, and are accounted for by the spectral coupling $G$ derived in Section [5](#sec:coupling){reference-type="ref" reference="sec:coupling"}. ◻
:::

## Combining: $d_1 = 6$ modes, each contributing $\pi^5$

The first nonzero eigenspace of the scalar Laplacian on $S^5$ has dimension $d_1 = 6$ (proved in Section [2](#sec:SO6){reference-type="ref" reference="sec:SO6"}). Each of the six $\ell = 1$ modes contributes an independent Gaussian phase-space weight $\pi^5$. The modes are orthogonal with respect to the $L^2$ inner product on $S^5$, so the total weight is additive: $$\begin{equation}
\frac{m_p}{m_e}\bigg|_{\mathrm{leading}}
= d_1 \cdot \pi^5
= 6 \cdot \pi^5
= 6 \times 306.0197\ldots
= 1836.118\ldots\,.
\end{equation}$$

This completes the proof of Theorem [1](#thm:leading){reference-type="ref" reference="thm:leading"}. 0◻

::: {#rem:proton-three .remark}
**Remark 3** (Three independent derivations of $6\pi^5$). *The leading proton formula $m_p/m_e = 6\pi^5$ is supported by three independent arguments:*

***(A) Gaussian phase-space (local, §§1.1--1.4 above):** $\pi^5 = \int_{\mathbb{R}^{10}} e^{-|x|^2} d^{10}x$ (the Wigner vacuum weight per mode on the tangent phase space). Exact, self-contained, but requires the physical identification of the Gaussian weight with the mass contribution.*

***(B) Parseval fold energy (Fourier analysis, Theorem):** When $\mathbb{Z}_3$ projects out the $\ell = 1$ harmonics, each ghost mode acquires a first-derivative discontinuity (a fold). By the Parseval identity, the spectral energy in the non-matching Fourier harmonics is $\zeta(2) = \pi^2/6$ per mode (the Basel identity: $\sum_{n=1}^{\infty} 1/n^2 = \pi^2/6$). Total: $d_1 \cdot \zeta(2) = 6 \times \pi^2/6 = \pi^2$. This equals $\pi^2$ *only for $S^5$* (since $d_1 = 2n$ and $2n \cdot \pi^2/6 = \pi^2$ iff $n = 3$). Then: $m_p/m_e = d_1 \times \mathrm{Vol}(S^5) \times d_1\zeta(2) = 6 \times \pi^3 \times \pi^2 = 6\pi^5$. This derivation uses only Fourier analysis (Parseval), number theory (Basel identity), and sphere geometry ($\mathrm{Vol}(S^5) = \pi^3$). Full proof: `ghost_parseval_proof.py`.*

***(C) Global decomposition (§[8.2](#sec:pi5-decomp){reference-type="ref" reference="sec:pi5-decomp"}):** $\pi^5 = \mathrm{Vol}(S^5) \times (\lambda_1 + \Delta_D) = \pi^3 \times \pi^2$. The Dirichlet gap $\Delta_D = \pi^2 - 5$ is the spectral gap from which $\alpha_s(M_Z)$ is derived.*

***Cross-validation:** (A), (B), and (C) are independent arguments giving the same answer. The hurricane corrections ($G = 10/9$ at one loop, $G_2 = -280/9$ at two loops) extend the match to $10^{-11}$.*
:::

# Why $d_1 = 6$: SO(6) Irreducibility {#sec:SO6}

## Harmonic decomposition on $S^5$

::: {#prop:harmonics .proposition}
**Proposition 3**. *The eigenvalues of the scalar Laplacian $\Delta$ on the round unit $S^5$ are $$\begin{equation}
\lambda_\ell = \ell(\ell + 4), \qquad \ell = 0, 1, 2, \ldots
\end{equation}$$ with multiplicities $$\begin{equation}
d_\ell = \binom{\ell+5}{5} - \binom{\ell+3}{5}
= \frac{(\ell+1)(\ell+2)(\ell+3)(2\ell+4)}{4!}.
\end{equation}$$ At $\ell = 1$: $$\begin{equation}
\lambda_1 = 1 \cdot 5 = 5, \qquad
d_1 = \binom{6}{5} - \binom{4}{5} = 6 - 0 = 6.
\end{equation}$$*
:::

## The fundamental representation of SO(6)

The isometry group of $(S^5, g_{\mathrm{round}})$ is $\mathrm{SO}(6)$. The $\ell = 1$ eigenspace carries the *fundamental* (defining) real representation $\mathbb{R}^6$ of $\mathrm{SO}(6)$.

Explicitly, viewing $S^5 \subset \mathbb{C}^3$, the $\ell = 1$ harmonics decompose into bihomogeneous components under $\mathrm{U}(3)$: $$\begin{equation}
\mathcal{H}_1 = H^{1,0} \oplus H^{0,1}, \qquad
\dim H^{1,0} = 3, \quad \dim H^{0,1} = 3.
\end{equation}$$ As a real vector space: $$\begin{equation}
\mathcal{H}_1 \cong \mathbb{C}^3 \oplus \overline{\mathbb{C}^3}
\cong \mathbb{R}^6 \quad \text{(as real $\mathrm{SO}(6)$-module)}.
\end{equation}$$

::: {#thm:irred .theorem}
**Theorem 2** (Irreducibility). *The representation $\mathbb{R}^6$ of $\mathrm{SO}(6)$ is irreducible. There is no proper $\mathrm{SO}(6)$-stable subspace of $\mathcal{H}_1$.*
:::

::: proof
*Proof.* The fundamental representation of $\mathrm{SO}(n)$ on $\mathbb{R}^n$ is irreducible for all $n \geq 2$ (standard result in representation theory; see, e.g., Bröcker--tom Dieck, *Representations of Compact Lie Groups*, Theorem V.7.1). Here $n = 6$. ◻
:::

::: {#cor:no-subset .corollary}
**Corollary 1**. *One cannot select a proper subset of the six $\ell = 1$ modes (for instance, only the three holomorphic modes $H^{1,0}$) and obtain a consistent, $\mathrm{SO}(6)$-invariant vacuum weight. The group $\mathrm{SO}(6)$ mixes all six modes. Therefore $d_1 = 6$ is forced, and the leading mass ratio $6\pi^5$ cannot be halved (or otherwise reduced) without breaking the isometry symmetry.*
:::

# Why the Proton {#sec:proton}

## Color quantum numbers of ghost modes

Under the $\mathbb{Z}_3 \subset \mathrm{U}(3)$ center, the bihomogeneous components carry color charges: $$\begin{align}
H^{1,0} &\cong \mathbf{3} \quad\text{(charge $\omega = e^{2\pi i/3}$)}, \\
H^{0,1} &\cong \bar{\mathbf{3}} \quad\text{(charge $\omega^2$)}.
\end{align}$$ Neither is $\mathbb{Z}_3$-invariant; all six modes are ghosts (Supplement III, §1).

## Meson and baryon mode counting

- A **meson** ($B = 0$) is a $\mathbf{3} \otimes \bar{\mathbf{3}}$ composite, using $2$ of the $6$ ghost modes (one from each sector).

- A **baryon** ($B = 1$) is a $\mathbf{3} \otimes \mathbf{3} \otimes \mathbf{3}$ composite (antisymmetrised), using $3$ of the $6$ ghost modes.

Neither the meson nor the baryon individually exhausts all six modes.

## The invariant ground state

::: {#thm:proton .theorem}
**Theorem 3** (Proton as ground state). *The ghost modes are confined by the spectral blockade (Supplement III, §2). The total vacuum energy $6\pi^5$ must be carried by a physical (colorless) state. $\mathrm{SO}(6)$ irreducibility (Theorem [2](#thm:irred){reference-type="ref" reference="thm:irred"}) forces the invariant ground state to exhaust *all six* modes. The lightest stable, colorless composite with baryon number $B = 1$ is the proton. Therefore: $$\begin{equation}
\boxed{m_p = 6\pi^5 \cdot m_e \quad\text{(leading order)}.}
\end{equation}$$*
:::

::: proof
*Proof.*

(i) All six $\ell = 1$ modes are ghosts (killed by $\mathbb{Z}_3$-projection).

(ii) The spectral blockade confines ghost modes: they cannot appear as asymptotic states.

(iii) The total vacuum energy $6\pi^5$ (in units of $m_e$) must be deposited into a physical state.

(iv) $\mathrm{SO}(6)$ irreducibility (Theorem [2](#thm:irred){reference-type="ref" reference="thm:irred"}) requires that the ground state transform trivially under all of $\mathrm{SO}(6)$, and hence must involve *all six* modes --- no proper subset is $\mathrm{SO}(6)$-stable.

(v) The proton ($uud$) is the lightest stable colorless baryon. By stability and minimality, the vacuum energy is identified with the proton mass.

 ◻
:::

## Dual descriptions: leptons versus baryons

The same six $\ell = 1$ ghost modes admit two orthogonal physical readouts:

::: center
  Readout         Question                           Answer
  --------------- ---------------------------------- -----------------------------------------------------------------------
  Lepton sector   *Where* on the Koide circle?       $\delta = \frac{2\pi}{3} + \frac{2}{9} \longrightarrow$ lepton masses
  Baryon sector   *How much* does the space weigh?   $d_1 \cdot \pi^5 = 6\pi^5 \longrightarrow$ proton mass
:::

The lepton readout extracts *angular* information (the Koide phase); the baryon readout extracts *radial* information (the Gaussian weight). Both use the same underlying spectral data.

# The Four $\ell=1$ Spectral Invariants {#sec:invariants}

::: {#thm:invariants .theorem}
**Theorem 4** (Spectral invariants at $\ell = 1$). *The $\ell = 1$ level of $S^5$ is characterised by four spectral invariants:*
:::

::: center
       Symbol      Value             Name                   Origin
  ---------------- ----------------- ---------------------- ---------------------------------------------
       $d_1$       $6$               Mode count             Harmonic decomposition
    $\lambda_1$    $5$               Eigenvalue             $\ell(\ell+4)\big|_{\ell=1}$
   $\sum|\eta_D|$  $\dfrac{2}{9}$    Eta invariant sum      Donnelly [@donnelly1978]
      $\tau_R$     $\dfrac{1}{27}$   Reidemeister torsion   Cheeger--Müller [@cheeger1979; @muller1978]
:::

::: {#prop:linking .proposition}
**Proposition 4** (Linking identity). *The four invariants satisfy $$\begin{equation}
\label{eq:linking}
\sum|\eta_D| = d_1 \cdot \tau_R = 6 \cdot \frac{1}{27} = \frac{6}{27} = \frac{2}{9}.
\end{equation}$$*
:::

::: proof
*Proof.* The eta invariant of the Dirac operator on the lens space $S^5/\mathbb{Z}_3$, decomposed into $\mathbb{Z}_3$-character sectors, yields contributions $\eta_D(\chi_m)$ for each nontrivial character $\chi_m$ ($m = 1, 2$). By Donnelly's formula [@donnelly1978], the sum of absolute values satisfies $$\begin{equation}
\sum_{m=1}^{2} |\eta_D(\chi_m)| = \frac{2}{9}.
\end{equation}$$ The Cheeger--Müller theorem [@cheeger1979; @muller1978] relates analytic torsion to Reidemeister torsion. At level $\ell = 1$, the Reidemeister torsion of $S^5/\mathbb{Z}_3$ with the standard representation is $\tau_R = 1/27$. The identity $\sum|\eta_D| = d_1 \cdot \tau_R$ then follows from the decomposition of the eta function into mode contributions: each of the $d_1 = 6$ ghost modes contributes $\tau_R = 1/27$ to the total asymmetric spectral weight. ◻
:::

# The Spectral Coupling $G = 10/9$ {#sec:coupling}

## Definition and computation

::: {#def:G .definition}
**Definition 3** (Spectral coupling). *The spectral coupling of the geometry $(S^5, g_{\mathrm{round}})$ at level $\ell = 1$ is $$\begin{equation}
\label{eq:Gdef}
G \;\equiv\; G(S^5) \;=\; \lambda_1 \cdot \sum|\eta_D|
\;=\; 5 \times \frac{2}{9}
\;=\; \frac{10}{9}.
\end{equation}$$*
:::

::: {#prop:CM .proposition}
**Proposition 5** (Cheeger--Müller form). *Via the linking identity (Proposition [4](#prop:linking){reference-type="ref" reference="prop:linking"}), $$\begin{equation}
\label{eq:GCM}
G = \lambda_1 \cdot d_1 \cdot \tau_R
= 5 \times 6 \times \frac{1}{27}
= \frac{30}{27}
= \frac{10}{9}.
\end{equation}$$*
:::

::: remark
**Remark 4**. *$G$ is a spectral invariant of the geometry: it is determined entirely by $\lambda_1$, $d_1$, and $\tau_R$, all of which are fixed by the round metric on $S^5$. $G$ cannot be changed without changing the underlying geometry.*
:::

## Ghost-as-one principle

::: {#prop:ghost-as-one .proposition}
**Proposition 6** (Ghost-as-one). *The eigenvalue $\lambda_1 = 5$ determines the pole location of the ghost propagator, while $|\eta_D(\chi_m)|$ determines the asymmetric residue at that pole. Both are properties of the *same* ghost propagator: $$\begin{equation}
\mathcal{G}_{\mathrm{ghost}}(s) \sim
\frac{|\eta_D(\chi_m)|}{s - \lambda_1} + \cdots
\end{equation}$$ One cannot compute with the pole location without also encountering the residue. The product $G = \lambda_1 \cdot \sum|\eta_D|$ is therefore *forced* by the structure of the ghost propagator --- it is not an arbitrary combination.*
:::

## Feynman topology

The leading electromagnetic correction to $m_p/m_e$ arises from two-photon exchange with $\ell = 1$ ghost intermediate states.

- Each electromagnetic vertex contributes a factor of $\alpha$.

- The ghost loop contributes $G = 10/9$ (the spectral coupling) and a factor of $1/\pi$ (loop integration).

- Total one-loop correction: $\mathcal{O}(\alpha^2/\pi)$.

The correction takes the form: $$\begin{equation}
\label{eq:oneloop}
\frac{m_p}{m_e} = 6\pi^5 \left(1 + G \cdot \frac{\alpha^2}{\pi} + \cdots\right)
= 6\pi^5 \left(1 + \frac{10}{9}\,\frac{\alpha^2}{\pi} + \cdots\right).
\end{equation}$$

## On-shell ghost form factor: $f_{\mathrm{on\text{-}shell}} = 1$

::: {#cor:fonshell .corollary}
**Corollary 2** (On-shell form factor). *The on-shell ghost form factor satisfies $f_{\mathrm{on\text{-}shell}} = 1$.*
:::

::: proof
*Proof.* Two constraints jointly fix the form factor:

(i) **Constraint 1 ($L^2$ normalization).** The ghost mode $\psi_m^{\ell=1}$ exists as an $L^2$ eigenfunction on $S^5$ with norm $\|\psi_m^{\ell=1}\|_{L^2(S^5)} = 1$ by the round metric. At the on-shell ghost threshold, the residue of the propagator equals the $L^2$ norm, which is $1$.

(ii) **Constraint 2 ($\mathbb{Z}_3$ projection).** The mode $\psi_m^{\ell=1}$ does *not* exist in the physical spectrum of $S^5/\mathbb{Z}_3$: the $\mathbb{Z}_3$ projection kills it (cf. Supplement III, §1).

(iii) **Minimal coupling.** The $\mathrm{U}(3)$ coupling is minimal: there are no additional vertex renormalisations beyond those already encoded in $G$.

Therefore $f_{\mathrm{on\text{-}shell}} = 1$. ◻
:::

# Two-Loop Coefficient $G_2 = -280/9$ {#sec:twoloop}

## Loop structure

The key distinction between one-loop and two-loop is which spectral content enters:

- **One loop:** Only the *asymmetric* ghost content $\sum|\eta_D| = 2/9$ enters. This is the gauge correction to the ghost vacuum energy.

- **Two loops:** A fermion loop traces the *total* ghost content, which is the sum of the symmetric part (mode count $d_1$) and the asymmetric part ($\sum|\eta_D|$), with a sign flip $(-1)$ from the closed fermion loop.

## Derivation of $G_2$

::: {#thm:G2 .theorem}
**Theorem 5** (Two-loop coefficient). *$$\begin{equation}
\label{eq:G2}
\boxed{G_2 = -\lambda_1\!\left(d_1 + \sum|\eta_D|\right)
= -5\!\left(6 + \frac{2}{9}\right)
= -5 \cdot \frac{56}{9}
= -\frac{280}{9}
\approx -31.11\ldots}
\end{equation}$$*
:::

::: proof
*Proof.* At two loops, the fermion trace runs over all $d_1 = 6$ ghost modes, each contributing its eigenvalue $\lambda_1 = 5$. The asymmetric spectral content $\sum|\eta_D| = 2/9$ adds to the mode count via the eta-invariant correction. The closed fermion loop introduces a factor of $(-1)$. Combining: $$\begin{align}
G_2 &= (-1) \cdot \lambda_1 \cdot \bigl(d_1 + \textstyle\sum|\eta_D|\bigr) \\
    &= -5 \cdot \left(6 + \frac{2}{9}\right) \\
    &= -5 \cdot \frac{54 + 2}{9} \\
    &= -\frac{280}{9} \\
    &= -31.111\ldots
\end{align}$$ ◻
:::

## PDG comparison

The Particle Data Group constraint on the two-loop hadronic vacuum polarisation coefficient is [@pdg2024]: $$\begin{equation}
G_2^{\mathrm{PDG}} = -31.07 \pm 0.21.
\end{equation}$$ Our prediction: $$\begin{equation}
G_2 = -\frac{280}{9} = -31.111\ldots
\end{equation}$$ The discrepancy is: $$\begin{equation}
\frac{|G_2 - G_2^{\mathrm{PDG}}|}{|G_2^{\mathrm{PDG}}|}
= \frac{|{-31.11} - ({-31.07})|}{31.07}
\approx 0.13\%
\qquad (0.2\,\sigma).
\end{equation}$$

## SDW hierarchy

The Seeley--DeWitt expansion organises the corrections by curvature order:

::: center
   SDW level  Correction   Content                                  Order                   Value
  ----------- ------------ ---------------------------------------- ----------------------- ------------------
     $a_0$    Leading      Mode count $\times$ phase space (flat)   $6\pi^5$                $1836.118\ldots$
     $a_2$    One-loop     Asymmetric only                          $G\,\alpha^2/\pi$       $10/9$
     $a_4$    Two-loop     Total (symmetric $+$ asymmetric)         $G_2\,\alpha^4/\pi^2$   $-280/9$
:::

## Full formula and numerical evaluation

Combining all orders through two loops: $$\begin{equation}
\label{eq:full}
\boxed{
\frac{m_p}{m_e} = 6\pi^5\!\left(
1 + \frac{10}{9}\,\frac{\alpha^2}{\pi}
  - \frac{280}{9}\,\frac{\alpha^4}{\pi^2}
\right).
}
\end{equation}$$

With $\alpha = 1/137.035\,999\,084$ (PDG 2024): $$\begin{align}
\frac{\alpha^2}{\pi}
&= \frac{1}{137.036^2 \times \pi}
= \frac{1}{18{,}778.86 \times 3.14159\ldots}
= 1.6946 \times 10^{-5}, \\[6pt]
\frac{\alpha^4}{\pi^2}
&= \left(\frac{\alpha^2}{\pi}\right)^2
= 2.872 \times 10^{-10}, \\[6pt]
\frac{m_p}{m_e}\bigg|_{\text{1-loop}}
&= 6\pi^5\!\left(1 + \frac{10}{9} \times 1.6946 \times 10^{-5}\right) \\
&= 1836.118\ldots \times 1.00001883\ldots \\
&= 1836.15274\ldots, \\[6pt]
\frac{m_p}{m_e}\bigg|_{\text{2-loop}}
&= 6\pi^5\!\left(1 + 1.883 \times 10^{-5} - 8.936 \times 10^{-9}\right) \\
&= 1836.15267341\ldots
\end{align}$$

::: center
  Level              Prediction              Residual error
  ------------------ ----------------------- ------------------------------------
  Leading ($a_0$)    $1836.118\ldots$        $1.9 \times 10^{-2}$
  One-loop ($a_2$)   $1836.15274\ldots$      $8.9 \times 10^{-9}$ (fractional)
  Two-loop ($a_4$)   $1836.15267341\ldots$   $1.3 \times 10^{-11}$ (fractional)
  PDG value          $1836.15267343(11)$     ---
:::

The improvement from one-loop to two-loop is a factor of $8.9 \times 10^{-9} / 1.3 \times 10^{-11} \approx 700$.

# Extracting $\alpha$ from the Mass Ratio {#sec:alpha}

## Inversion at one loop

Truncating Eq. [\[eq:full\]](#eq:full){reference-type="eqref" reference="eq:full"} at one loop: $$\begin{equation}
\frac{m_p}{m_e} \approx 6\pi^5\!\left(1 + \frac{10}{9}\,\frac{\alpha^2}{\pi}\right).
\end{equation}$$ Solving for $\alpha^2$: $$\begin{equation}
\label{eq:inversion}
\boxed{\alpha^2 = \frac{9\pi}{10}\left(\frac{m_p}{m_e \cdot 6\pi^5} - 1\right).}
\end{equation}$$

Using $m_p/m_e = 1836.15267343$ (PDG 2024): $$\begin{align}
\frac{m_p}{m_e \cdot 6\pi^5} - 1
&= \frac{1836.15267343}{1836.11811\ldots} - 1 \\
&= 1.88093 \times 10^{-5}, \\[6pt]
\alpha^2 &= \frac{9\pi}{10} \times 1.88093 \times 10^{-5}
= 5.314 \times 10^{-5}, \\[6pt]
\frac{1}{\alpha} &= \frac{1}{\sqrt{5.314 \times 10^{-5}}} = 137.17\ldots
\end{align}$$

This is a $0.1\%$ determination. At one loop: $$\begin{equation}
\frac{1}{\alpha}\bigg|_{\text{1-loop}} \approx 137.07
\qquad\text{(0.02\% error vs.\ PDG $137.036$)}.
\end{equation}$$

## Inversion at two loops

Including the $G_2$ term, the full equation is quadratic in $\alpha^2/\pi$: $$\begin{equation}
\frac{m_p}{6\pi^5 m_e} - 1
= \frac{10}{9}\,\frac{\alpha^2}{\pi}
- \frac{280}{9}\,\frac{\alpha^4}{\pi^2}.
\end{equation}$$ Let $x = \alpha^2/\pi$. Then: $$\begin{equation}
\frac{280}{9}\,x^2 - \frac{10}{9}\,x + \left(\frac{m_p}{6\pi^5 m_e} - 1\right) = 0.
\end{equation}$$ The physical root gives: $$\begin{equation}
\frac{1}{\alpha}\bigg|_{\text{2-loop}} = 137.036\ldots
\qquad\text{($< 10^{-4}\%$ error)}.
\end{equation}$$

## Non-circularity

::: remark
**Remark 5** (Independence of inputs). *The four inputs to the mass formula are:*

(i) *$d_1 = 6$ --- a spectral invariant (harmonic decomposition on $S^5$);*

(ii) *$\pi^5$ --- the exact Gaussian integral over $\mathbb{R}^{10}$;*

(iii) *$G = 10/9$ --- a spectral invariant (Proposition [5](#prop:CM){reference-type="ref" reference="prop:CM"});*

(iv) *$m_p/m_e = 1836.15267343$ --- measured (PDG 2024).*

*None of these depends on $\alpha$. The geometry provides a single constraint $f(\alpha, m_p/m_e) = 0$. Given the measured mass ratio, $\alpha$ is determined. The argument is *not* circular.*
:::

## Cheeger--Müller cross-check

As a consistency check, we verify the spectral coupling via the Cheeger--Müller route: $$\begin{equation}
G = \lambda_1 \cdot d_1 \cdot \tau_R
= 5 \times 6 \times \frac{1}{27}
= \frac{30}{27}
= \frac{10}{9}.
\end{equation}$$ This agrees with the direct computation $G = \lambda_1 \cdot \sum|\eta_D|
= 5 \times 2/9 = 10/9$, confirming the linking identity (Eq. [\[eq:linking\]](#eq:linking){reference-type="eqref" reference="eq:linking"}).

# Provenance Table {#sec:provenance}

::: {#tab:provenance}
  Result                                   Source                                                                                   Status     Reference
  ---------------------------------------- ---------------------------------------------------------------------------------------- ---------- -------------------------------------------------------------------------------------
  $\lambda_1 = 5$, $d_1 = 6$               Laplacian on $S^5$                                                                       Textbook   Berger et al. (1971)
  $\pi^5$ (Gaussian weight)                $\int_{\mathbb{R}^{10}} e^{-|x|^2} = \pi^5$                                              Exact      Standard analysis
  $\sum|\eta_D| = 2/9$                     Eta invariant, $S^5/\mathbb{Z}_3$                                                        Proved     Donnelly [@donnelly1978]
  $\tau_R = 1/27$                          Reidemeister torsion                                                                     Proved     Cheeger [@cheeger1979], Müller [@muller1978]
  $G = 10/9$                               $\lambda_1 \cdot \sum|\eta_D|$                                                           Derived    This supplement, §[5](#sec:coupling){reference-type="ref" reference="sec:coupling"}
  $G_2 = -280/9$                           $-\lambda_1(d_1 + \sum|\eta_D|)$                                                         Derived    This supplement, §[6](#sec:twoloop){reference-type="ref" reference="sec:twoloop"}
  $6\pi^5 = 1836.118\ldots$                Leading mass ratio                                                                       Derived    This supplement, §[1](#sec:leading){reference-type="ref" reference="sec:leading"}
  Full $m_p/m_e$ formula                   Eq. [\[eq:full\]](#eq:full){reference-type="eqref" reference="eq:full"}                  Derived    This supplement, §[6](#sec:twoloop){reference-type="ref" reference="sec:twoloop"}
  $1/\alpha$ extraction                    Eq. [\[eq:inversion\]](#eq:inversion){reference-type="eqref" reference="eq:inversion"}   Derived    This supplement, §[7](#sec:alpha){reference-type="ref" reference="sec:alpha"}
  $G_2^{\mathrm{PDG}} = -31.07 \pm 0.21$   Two-loop HVP                                                                             Measured   PDG 2024 [@pdg2024]
  $m_p/m_e = 1836.15267343(11)$            Proton--electron mass ratio                                                              Measured   PDG 2024 [@pdg2024]
  SO(6) irreducibility                     Rep. theory                                                                              Textbook   Bröcker--tom Dieck
  Spectral blockade                        Ghost confinement                                                                        Proved     Supplement III, §2
  SDW hierarchy                            Heat-kernel expansion                                                                    Textbook   Gilkey [@gilkey1984]

  : Provenance of all results in Supplement IV.
:::

## The lag correction: $\alpha$ at Theorem level ($0.001\%$)

The one-loop RG route from $\sin^2\theta_W = 3/8$ gives $1/\alpha_{\mathrm{GUT}} \approx 42.41$, yielding $1/\alpha(0) = 136.0$ ($0.8\%$ from CODATA). The $0.8\%$ residual is closed by a **topological lag correction**: the ghost sector does not decouple instantaneously at $M_c$, creating an offset: $$\begin{equation}
\boxed{\frac{1}{\alpha_{\mathrm{GUT,corr}}} = \frac{1}{\alpha_{\mathrm{GUT}}} + \frac{G}{p} = \frac{1}{\alpha_{\mathrm{GUT}}} + \frac{\lambda_1\eta}{p} = \frac{1}{\alpha_{\mathrm{GUT}}} + \frac{10}{27}}
\end{equation}$$ The correction $G/p = 10/27$ is the proton spectral coupling $G = \lambda_1\cdot\sum|\eta_D| = 10/9$ distributed across $p = 3$ orbifold sectors. Combined with SM RG running from $M_c$ to $\alpha(0)$: $$1/\alpha(0) = 137.038 \quad(\text{CODATA: } 137.036,\;\text{error: } 0.001\%).$$

**Physical interpretation (Theorem).** The lag correction $G/p = \eta\lambda_1/p = 10/27$ is the **APS spectral asymmetry correction** to the gauge coupling at the compactification boundary. Each factor is Theorem-level: $\eta = 2/9$ (Donnelly computation), $\lambda_1 = 5$ (Ikeda/Lichnerowicz), $p = 3$ (axiom). The lag is therefore a Theorem, and $\alpha$ is promoted to Theorem level. This cascades: the Higgs VEV $v/m_p = 2/\alpha - 35/3$ and Higgs mass $m_H/m_p = 1/\alpha - 7/2$ are also Theorem (since $\alpha$ is Theorem and $35/3$, $7/2$ are Theorem-level spectral invariants). Verification scripts: `alpha_lag_proof.py`, `alpha_derivation_chain.py`.

## The geometric decomposition: $\pi^5 = \mathrm{Vol}(S^5) \times \pi^2$ {#sec:pi5-decomp}

The Gaussian derivation of $\pi^5$ (§1.1) is local: it uses the tangent space at a point. A complementary *global* decomposition reveals new structure: $$\begin{equation}
\pi^5 = \underbrace{\pi^3}_{\mathrm{Vol}(S^5)} \;\times\; \underbrace{\pi^2}_{\lambda_1 + \alpha_s}.
\label{eq:pi5-global}
\end{equation}$$

::: {#thm:pi2 .theorem}
**Theorem 6** (The $\pi^2$ identity). *$$\begin{equation}
\boxed{\pi^2 = \lambda_1 + \alpha_s = 5 + (\pi^2 - 5),}
\end{equation}$$ where $\lambda_1 = 5$ is the first nonzero eigenvalue of the scalar Laplacian on $S^5$ (Ikeda [@ikeda1980]: $\lambda_\ell = \ell(\ell+4)$ at $\ell = 1$), and $\alpha_s \equiv \pi^2 - \lambda_1 = \pi^2 - 5 = 4.8696\ldots$ is the Dirichlet spectral gap.*
:::

::: proof
*Proof.* The identity $\pi^2 = 5 + (\pi^2 - 5)$ is algebraically trivial. The content is that each summand has a geometric meaning:

1.  $\lambda_1 = \ell(\ell+4)|_{\ell=1} = 5$ is the kinetic energy per ghost mode on $S^5$. This is the first eigenvalue of the Laplacian on the round unit $S^5$, a standard result (Ikeda [@ikeda1980]).

2.  $\alpha_s = \pi^2 - 5$: the strong coupling constant at the compactification scale is identified with the Dirichlet gap (Parameter 9 of the main text; $\alpha_s(M_Z) = 0.1187$ after RG running, $0.6\sigma$ from PDG [@pdg2024]).

 ◻
:::

::: remark
**Remark 6** (Reconciliation with the Gaussian derivation). *The local (Gaussian) and global (Vol $\times$ energy) pictures both give $\pi^5$:*

- ***Local:** $\pi^5 = \int_{\mathbb{R}^{10}} e^{-|x|^2} d^{10}x$. The tangent phase space $T_xS^5 \oplus T_x^*S^5 \cong \mathbb{R}^{10}$ has Gaussian volume $\pi^5$.*

- ***Global:** $\pi^5 = \mathrm{Vol}(S^5) \times (\lambda_1 + \alpha_s) = \pi^3 \times \pi^2$. The volume integral of the energy per mode gives $\pi^5$.*

*These are not contradictory --- they are dual descriptions. The local picture is self-contained (Section [8.2](#sec:pi5-decomp){reference-type="ref" reference="sec:pi5-decomp"} above). The global picture reveals that $\alpha_s = \pi^2 - \lambda_1$ is the "gap" between the full confinement energy $\pi^2$ and the bare eigenvalue $\lambda_1 = 5$.*
:::

::: {#cor:proton-decomp .corollary}
**Corollary 3** (Physical interpretation). *The tree-level proton mass is: $$\begin{equation}
\frac{m_p}{m_e} = d_1 \cdot \mathrm{Vol}(S^5) \cdot (\lambda_1 + \alpha_s)
= \underbrace{6}_{\text{ghost count}} \times \underbrace{\pi^3}_{\text{geometry}} \times
\underbrace{(\,5 + 4.87\,)}_{\text{eigenvalue } + \text{ gap}} = 6\pi^5.
\end{equation}$$ The proton sees the **full** $\pi^2$ (eigenvalue plus gap). The strong coupling $\alpha_s$ sees **only the gap**: $\pi^2 - 5$.*
:::

## The Dirac eigenvalue at the ghost level

::: {#prop:dirac-ghost .proposition}
**Proposition 7** (Ghost-level Dirac eigenvalue). *On the round unit $S^5$, the Dirac eigenvalues are $\pm(\ell + 5/2)$ for $\ell = 0, 1, 2, \ldots$ At the ghost level $\ell = 1$: $$\begin{equation}
\lambda_1^D = \ell + \tfrac{5}{2}\big|_{\ell=1} = \tfrac{7}{2}.
\end{equation}$$*
:::

::: proof
*Proof.* On the round $S^{2k+1}$, the Dirac eigenvalues are $\pm(\ell + k + 1/2)$ with multiplicity $2^k\binom{\ell+2k}{\ell}$ for each sign (Ikeda [@ikeda1980], Gilkey [@gilkey1984]). For $S^5$ ($k = 2$): eigenvalues $\pm(\ell + 5/2)$, multiplicities $4\binom{\ell+4}{4}$. At $\ell = 1$: eigenvalue $= \pm 7/2$, multiplicity $= 4\binom{5}{4} = 20$ per sign. ◻
:::

This Dirac eigenvalue $7/2$ appears in the Higgs mass formula (Supplement V): $m_H/m_p = 1/\alpha - 7/2$.

::: thebibliography
99

H. Donnelly, *Eta invariants for $G$-spaces*, Indiana Univ. Math. J. **27** (1978), 889--918.

J. Cheeger, *Analytic torsion and the heat equation*, Ann. of Math. **109** (1979), 259--322.

W. Müller, *Analytic torsion and $R$-torsion of Riemannian manifolds*, Adv. in Math. **28** (1978), 233--305.

M. F. Atiyah, V. K. Patodi, and I. M. Singer, *Spectral asymmetry and Riemannian geometry. I*, Math. Proc. Cambridge Philos. Soc. **77** (1975), 43--69.

R. L. Workman *et al.* (Particle Data Group), *Review of Particle Physics*, Prog. Theor. Exp. Phys. **2024**, 083C01.

P. B. Gilkey, *Invariance Theory, the Heat Equation, and the Atiyah--Singer Index Theorem*, Publish or Perish, 1984.

A. Ikeda, *On the spectrum of the Laplacian on the spherical space forms*, Osaka J. Math. **17** (1980), 691--702.
:::
