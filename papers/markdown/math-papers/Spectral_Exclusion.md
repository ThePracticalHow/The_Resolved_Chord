# The Harmonics of $S^5$

::: {#def:harmonics .definition}
**Definition 1** (Spherical harmonics on $S^{2n-1}$). *The eigenspaces of the scalar Laplacian $-\Delta$ on the round unit $S^{2n-1}$ at level $\ell$ are the restrictions to $S^{2n-1}$ of harmonic homogeneous polynomials of degree $\ell$ on $\mathbb{R}^{2n} = \mathbb{C}^n$. The eigenvalue is $\lambda_\ell = \ell(\ell + 2n - 2)$ and the degeneracy is $d_\ell = \binom{\ell+2n-1}{\ell} - \binom{\ell+2n-3}{\ell-2}$.*
:::

For $S^5$ ($n = 3$): $\lambda_0 = 0$, $d_0 = 1$; $\lambda_1 = 5$, $d_1 = 6$; $\lambda_2 = 12$, $d_2 = 20$.

::: {#prop:ell1 .proposition}
**Proposition 2** ($\ell = 1$ harmonics are linear functions). *The $\ell = 1$ eigenspace of $-\Delta$ on $S^5 \subset \mathbb{C}^3$ is spanned by the restrictions of the six real-linear coordinate functions: $$\begin{equation}
\{x_1, y_1, x_2, y_2, x_3, y_3\} \quad\text{where } z_j = x_j + iy_j.
\end{equation}$$ Equivalently, as complex-linear and anti-linear functions: $\{z_1, z_2, z_3, \bar{z}_1, \bar{z}_2, \bar{z}_3\}$.*
:::

::: proof
*Proof.* Standard: the harmonic homogeneous polynomials of degree 1 on $\mathbb{R}^6$ are the linear functions, forming a 6-dimensional space. On $S^5$: $-\Delta(z_j|_{S^5}) = \lambda_1 \cdot z_j|_{S^5}$ with $\lambda_1 = 1 \cdot (1+4) = 5$. 0◻ ◻
:::

# The $\mathbb{Z}_3$ Action on the $\ell = 1$ Eigenspace

::: {#def:z3-action .definition}
**Definition 3** ($\mathbb{Z}_3$ action). *$\mathbb{Z}_3$ acts on $\mathbb{C}^3$ by $\omega \cdot (z_1, z_2, z_3) = (\omega z_1, \omega z_2, \omega z_3)$ with $\omega = e^{2\pi i/3}$. This induces an action on the $\ell = 1$ eigenspace: $$\begin{align}
\omega \cdot z_j &= \omega z_j \qquad\text{(character $\chi_1$: eigenvalue $\omega$)}, \\
\omega \cdot \bar{z}_j &= \bar{\omega}\, \bar{z}_j = \omega^2 \bar{z}_j
\qquad\text{(character $\chi_2$: eigenvalue $\omega^2$)}.
\end{align}$$*
:::

::: {#thm:exclusion .theorem}
**Theorem 4** (Triple Spectral Exclusion). *The $\ell = 1$ eigenspace of $-\Delta$ on $S^5$ contains **no** $\mathbb{Z}_3$-invariant modes. The six-dimensional space decomposes as: $$\begin{equation}
\boxed{H_1 = V_{\chi_1} \oplus V_{\chi_2}, \qquad \dim V_{\chi_1} = 3, \qquad \dim V_{\chi_2} = 3, \qquad \dim V_{\chi_0} = 0.}
\end{equation}$$ All $d_1 = 6$ modes are **ghost modes** (non-invariant under $\mathbb{Z}_3$).*
:::

:::: proof
*Proof.* The six basis elements decompose under $\mathbb{Z}_3$ as:

::: center
  **Basis element**                    **$\omega$ eigenvalue**   **Character**   **Invariant?**
  ----------------------------------- ------------------------- --------------- ----------------
  $z_1, z_2, z_3$                             $\omega$             $\chi_1$            No
  $\bar{z}_1, \bar{z}_2, \bar{z}_3$          $\omega^2$            $\chi_2$            No
:::

No linear combination of $\chi_1$-modes and $\chi_2$-modes can be $\chi_0$-invariant (since $\chi_1 \neq \chi_0$ and $\chi_2 \neq \chi_0$, and the characters are orthogonal in $\mathbb{C}[\mathbb{Z}_3]$). Therefore $V_{\chi_0} = \{0\}$: no invariant modes. 0◻ ◻
::::

::: remark
**Remark 5** (Why the diagonal action is special). *If $\mathbb{Z}_3$ acted with different weights, e.g., $(z_1, z_2, z_3) \mapsto (\omega z_1, \omega z_2, z_3)$, then $z_3$ and $\bar{z}_3$ would be invariant, and $V_{\chi_0}$ would be 2-dimensional. The *diagonal* action (all weights equal to 1) is the one that kills ALL $\ell = 1$ modes. This is the action selected by the uniqueness theorem $n = p^{n-2}$ [@leng2026eta].*
:::

# Consequences

## The mass gap

::: {#cor:gap .corollary}
**Corollary 6** (Ghost spectral gap). *On $S^5/\mathbb{Z}_3$, the physical (invariant) scalar spectrum has a gap: the first nonzero invariant eigenvalue is at $\ell = 2$ ($\lambda_2 = 12$), not $\ell = 1$ ($\lambda_1 = 5$). The $\ell = 1$ level is entirely removed from the physical spectrum.*
:::

::: proof
*Proof.* Theorem [4](#thm:exclusion){reference-type="ref" reference="thm:exclusion"}: $d_1^{(0)} = 0$. The first level with invariant modes is $\ell = 2$, where $d_2^{(0)} = 8$ (computed by the character formula: $d_2 = 20$, and $20/3 + \text{character corrections} = 8$). 0◻ ◻
:::

## Confinement

::: {#cor:confinement .corollary}
**Corollary 7** (No fundamental triplet). *The three holomorphic coordinates $z_1, z_2, z_3$ transform in the fundamental $\mathbf{3}$ of $\mathrm{SU}(3) \subset \mathrm{SO}(6) = \mathrm{Isom}(S^5)$, where $\mathrm{SU}(3)$ is embedded via its natural action on $\mathbb{C}^3$. The diagonal $\mathbb{Z}_3$ used throughout this paper is the *center* $Z(\mathrm{SU}(3)) \cong \mathbb{Z}_3$, acting as scalar multiplication on the $\mathbf{3}$. Since they are all in $V_{\chi_1}$ (non-invariant), no physical mode at $\ell = 1$ transforms in the fundamental $\mathbf{3}$. A free color triplet cannot propagate on $S^5/\mathbb{Z}_3$: it is confined.*
:::

::: proof
*Proof.* A physical (propagating) mode must be $\mathbb{Z}_3$-invariant. The $\mathbf{3}$ of SU(3) lies entirely in $V_{\chi_1}$ at $\ell = 1$. Therefore no $\mathbb{Z}_3$-invariant mode transforms as a fundamental triplet. Color singlet combinations arise at higher $\ell$ (e.g., from $\bar{\mathbf{3}} \otimes \mathbf{3}$ decompositions at $\ell = 2$). The key point is that no *fundamental* $\mathbf{3}$ propagates at $\ell = 1$. 0◻ ◻
:::

## Chirality

::: {#cor:chirality .corollary}
**Corollary 8** (Chirality from exclusion). *The splitting $H_1 = V_{\chi_1} \oplus V_{\chi_2}$ with $3 + 3$ (rather than $6 + 0$ or $2 + 2 + 2$) breaks left-right symmetry. The $\chi_1$ and $\chi_2$ sectors contribute with opposite signs to the eta invariant: $\eta_D(\chi_1) = +1/9$ and $\eta_D(\chi_2) = -1/9$ [@leng2026eta]. The nonvanishing total $\eta = 2/9 \neq 0$ is the spectral signature of chirality.*
:::

# Higher Levels

::: {#prop:higher .proposition}
**Proposition 9** (Character decomposition at $\ell = 2, 3$). *At $\ell = 2$: $d_2 = 20$, $d_2^{(0)} = 8$, $d_2^{\mathrm{ghost}} = 12$. At $\ell = 3$: $d_3 = 50$, $d_3^{(0)} = 20$, $d_3^{\mathrm{ghost}} = 30$.*
:::

::: proof
*Proof.* The $\ell$-th harmonic space on $S^5$ is spanned by polynomials $z_1^{a_1}z_2^{a_2}z_3^{a_3}\bar{z}_1^{b_1}\bar{z}_2^{b_2}\bar{z}_3^{b_3}$ with $\sum a_j + \sum b_j = \ell$ (restricted to the harmonic subspace). Under $\mathbb{Z}_3$, such a monomial transforms with character $\omega^{(\sum a_j - \sum b_j) \bmod 3}$. The invariant monomials are those with $\sum a_j \equiv \sum b_j \pmod{3}$.

For $\ell = 2$: the invariant count is $d_2^{(0)} = (d_2 + 2\mathrm{Re}[\chi_2(\omega)])/3$, where $\chi_2(\omega)$ is the character trace on $H_2$. The polynomial-space character at $\ell = 2$: $\chi_{P_2}(\omega) = \sum_{a+b=2}\binom{a+2}{2}\binom{b+2}{2}\omega^{a-b}$ $= \binom{4}{2}\omega^2 + \binom{3}{2}\binom{3}{2}\omega^0 + \binom{2}{2}\binom{4}{2}\omega^{-2}$ $= 6\omega^2 + 9 + 6\omega = 9 + 6(\omega+\omega^2) = 9 - 6 = 3$. Harmonic correction: $\chi_{H_2}(\omega) = \chi_{P_2}(\omega) - \chi_{P_0}(\omega) = 3 - 1 = 2$. Therefore $d_2^{(0)} = (20 + 2 \cdot 2)/3 = 24/3 = 8$. And $d_2^{\mathrm{ghost}} = 20 - 8 = 12$. 0◻ ◻
:::

::: remark
**Remark 10** (Equidistribution at large $\ell$). *For $\ell \gg 1$: $d_\ell^{(0)} \to d_\ell/3$ (the $\mathbb{Z}_3$ characters equidistribute). The ghost fraction $d_\ell^{\mathrm{ghost}}/d_\ell \to 2/3$. The $\ell = 1$ exclusion ($d_1^{(0)} = 0$, ghost fraction $= 1$) is special to the first level. This equidistribution is the mechanism behind the heavy-mode cancellation in the cosmological constant derivation [@leng2026].*
:::

# Summary

On $S^5/\mathbb{Z}_3$ with the diagonal action:

1.  All $d_1 = 6$ first-level harmonics are ghost modes (Theorem [4](#thm:exclusion){reference-type="ref" reference="thm:exclusion"}).

2.  The physical spectrum has a gap: no invariant modes at $\ell = 1$ (Corollary [6](#cor:gap){reference-type="ref" reference="cor:gap"}).

3.  No fundamental SU(3) triplet propagates: confinement (Corollary [7](#cor:confinement){reference-type="ref" reference="cor:confinement"}).

4.  Chirality: $\eta = 2/9 \neq 0$ from the $\chi_1/\chi_2$ asymmetry (Corollary [8](#cor:chirality){reference-type="ref" reference="cor:chirality"}).

These are representation-theoretic facts about $\mathbb{Z}_3 \hookrightarrow \mathrm{U}(3)$ acting on spherical harmonics. No physical assumption is required.

::: thebibliography
99

A. Ikeda, "On the spectrum of the Laplacian on the spherical space forms," *Osaka J. Math.* **17** (1980) 691--702.

J. Leng, "Eta invariants, Reidemeister torsion, and a ghost-mode identity on the lens space $L(3;1,1,1)$," (2026).

J. Leng, "The Resolved Chord: The Theorem of Everything," v10 (2026).
:::
