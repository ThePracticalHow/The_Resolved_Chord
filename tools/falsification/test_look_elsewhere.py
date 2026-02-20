"""
Look-Elsewhere Correction (LEE)
==============================
Quantitative p_LEE under null: "best-of-N" Monte Carlo.
Blocks: "You scanned 96 models and reported the best hit."

Match score: S(δ) = max(|m_μ_pred/m_μ_PDG - 1|, |m_τ_pred/m_τ_PDG - 1|)
             (mass prediction accuracy, not Koide K which is always 2/3 for r=√2).

Note: K = (1 + r²/2)/3 = 2/3 identically for all δ when r = √2 (Brannen's algebraic
identity). The LEE test therefore measures the physically meaningful quantity: whether
the predicted mass RATIOS m_μ/m_e and m_τ/m_e match PDG, which is a strong constraint
on δ. The theoretical δ = 2π/3 + 2/9 predicts both to < 0.01%, while random δ values
from [0, 2π] almost never achieve this simultaneously.

Null: random δ ~ Uniform[0, 2π] (broadest possible prior; conservative for the claim).
"""

import numpy as np
import pytest

from leng_test_utils import (
    brannen_masses,
    mu_from_electron,
    koide_ratio,
    R,
    M_E_PDG,
    M_MU_PDG,
    M_TAU_PDG,
    K_THEORY,
)

N_CANDIDATES = 96  # n∈[1,8], p∈[2,13]
N_MONTECARLO = 100_000  # trials for empirical p_LEE (scaled for referee robustness)
DELTA_RANGE = (0.0, 2 * np.pi)  # broadest prior; conservative against the framework
RANDOM_SEED = 42  # Reproducibility; documented for Supplement C


def match_score_koide(masses: np.ndarray) -> float:
    """Koide deviation: |K - 2/3|. Always ≈ 0 for r=√2 Brannen masses (Brannen identity).
    Not used for LEE; retained for other tests."""
    return abs(koide_ratio(masses) - K_THEORY)


def match_score_mass(delta: float) -> float:
    """Mass prediction accuracy: max relative error of (m_μ, m_τ) vs PDG.
    Uses m_e as calibration input (scale parameter). Lower = better.
    This is the meaningful LEE score: K=2/3 is trivial, mass ratios are not."""
    mu = mu_from_electron(M_E_PDG, R, delta)
    masses = np.sort(brannen_masses(mu, R, delta))
    if len(masses) < 3 or np.any(masses <= 0):
        return 1.0  # unphysical
    rel_err_mu = abs(masses[1] - M_MU_PDG) / M_MU_PDG
    rel_err_tau = abs(masses[2] - M_TAU_PDG) / M_TAU_PDG
    return float(max(rel_err_mu, rel_err_tau))


def match_score_residual(m_e: float, m_mu: float, m_tau: float) -> float:
    """Relative residual norm vs PDG. Lower = better."""
    pred = np.array([m_e, m_mu, m_tau])
    obs = np.array([M_E_PDG, M_MU_PDG, M_TAU_PDG])
    return float(np.sqrt(np.sum(((pred - obs) / obs) ** 2)))


def best_match_under_null_mass(n_trials: int, delta_lo: float, delta_hi: float, seed: int = RANDOM_SEED) -> np.ndarray:
    """For each trial: draw random δ from [0, 2π], return mass prediction score."""
    rng = np.random.default_rng(seed)
    scores = np.zeros(n_trials)
    for t in range(n_trials):
        delta = rng.uniform(delta_lo, delta_hi)
        scores[t] = match_score_mass(delta)
    return scores


def full_pipeline_null_scores(
    n_trials: int,
    n_max: int = 8,
    p_max: int = 13,
    seed: int = RANDOM_SEED,
) -> np.ndarray:
    """
    Full-pipeline null: sample (n,p) from candidate set and δ from [0,2π].
    For (n,p)=(3,3) only does the framework predict 3 lepton masses; else score=1.
    Used for referee-robust 'full pipeline null' (ChatGPT recommendation).
    """
    rng = np.random.default_rng(seed)
    scores = np.ones(n_trials)
    for t in range(n_trials):
        n = rng.integers(1, n_max + 1)
        p = rng.integers(2, p_max + 1)
        delta = rng.uniform(DELTA_RANGE[0], DELTA_RANGE[1])
        if n == 3 and p == 3:
            scores[t] = match_score_mass(delta)
        # else: score stays 1 (worst)
    return scores


def best_of_N_empirical(null_scores: np.ndarray, observed_score: float, N: int) -> float:
    """
    Empirical p_LEE: fraction of trials where best-of-N under null ≤ observed.
    Uses conservative approximation: each trial gives one 'best' draw from null.
    """
    # Each null trial is one random geometry's score. Best-of-N means we take
    # the minimum over N such draws. We simulate: for each of n_trials 'scans',
    # take min of N random null scores; count how often that min ≤ observed.
    rng = np.random.default_rng(43)
    n_trials = len(null_scores)
    best_of_N = np.minimum.reduceat(
        rng.choice(null_scores, size=n_trials * N),
        np.arange(0, n_trials * N, N)
    )[:n_trials]
    return np.mean(best_of_N <= observed_score)


def p_LEE_analytic(p_single: float, N: int) -> float:
    """p_LEE ≈ 1 - (1 - p_single)^N (simplified LEE formula)."""
    return 1.0 - (1.0 - p_single) ** N


def wilson_ci(hits: int, n: int, z: float = 1.96) -> tuple:
    """Wilson score interval for binomial proportion. z=1.96 gives 95% CI."""
    p = hits / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denom
    margin = z * np.sqrt((p * (1 - p) + z**2 / (4 * n)) / n) / denom
    return (max(0, center - margin), min(1, center + margin))


class TestLookElsewhereCorrection:
    """Quantitative LEE for Supplement C."""

    def test_observed_mass_score(self):
        """Observed mass prediction score: max relative error vs PDG."""
        from leng_test_utils import DELTA
        obs_score = match_score_mass(DELTA)
        assert obs_score < 0.001  # < 0.1% max relative error
        assert obs_score < 0.0001  # < 0.01% (LENG actually achieves ~0.007%)

    def test_observed_residual_score(self):
        """Observed residual norm vs PDG."""
        from leng_test_utils import predict_masses_from_e
        m_e, m_mu, m_tau = predict_masses_from_e(M_E_PDG)
        score = match_score_residual(m_e, m_mu, m_tau)
        assert score < 0.01  # < 1% RMS relative error

    def test_null_mass_hits_rare(self):
        """Under uniform δ in [0, 2π], hitting mass prediction score < 0.01% is rare."""
        scores = best_match_under_null_mass(
            N_MONTECARLO, DELTA_RANGE[0], DELTA_RANGE[1]
        )
        from leng_test_utils import DELTA
        obs_score = match_score_mass(DELTA)
        hit_threshold = max(obs_score * 2, 0.0001)  # 2× observed or 0.01% floor
        hits = np.sum(scores < hit_threshold)
        hit_rate = hits / N_MONTECARLO
        lo, hi = wilson_ci(int(hits), N_MONTECARLO)
        assert hit_rate < 0.05  # < 5% of random δ predict masses this accurately
        assert hi < 0.10  # 95% CI upper bound < 10%

    def test_effect_size_null_vs_observed(self):
        """Median null score >> observed score (effect size)."""
        scores = best_match_under_null_mass(
            N_MONTECARLO, DELTA_RANGE[0], DELTA_RANGE[1]
        )
        from leng_test_utils import DELTA
        obs_score = match_score_mass(DELTA)
        median_null = np.median(scores)
        ratio = median_null / max(obs_score, 1e-12)
        assert ratio > 10  # Observed is orders of magnitude better than typical null

    def test_p_single_mass(self):
        """Single-candidate p under null: fraction hitting observed tolerance."""
        scores = best_match_under_null_mass(
            N_MONTECARLO, DELTA_RANGE[0], DELTA_RANGE[1]
        )
        from leng_test_utils import DELTA
        obs_score = match_score_mass(DELTA)
        p_single = np.mean(scores <= obs_score)
        assert p_single < 0.1  # < 10% under broad prior

    def test_analytic_lee_formula(self):
        """p_LEE = 1 - (1 - p_single)^N."""
        p_single = 0.01  # example
        p_lee = p_LEE_analytic(p_single, N_CANDIDATES)
        assert p_lee > p_single
        assert p_lee == pytest.approx(1 - (0.99 ** 96), rel=0.01)

    def test_wilson_ci_computable(self):
        """Wilson CI for hit rate: observed h, estimated p_single = h/M, 95% CI."""
        scores = best_match_under_null_mass(
            N_MONTECARLO, DELTA_RANGE[0], DELTA_RANGE[1]
        )
        from leng_test_utils import DELTA
        obs_score = match_score_mass(DELTA)
        hits = np.sum(scores <= obs_score)
        p_single = hits / N_MONTECARLO
        lo, hi = wilson_ci(int(hits), N_MONTECARLO)
        assert lo <= p_single <= hi
        assert 0 <= lo <= 1 and 0 <= hi <= 1

    def test_full_pipeline_null_hits_rare(self):
        """Full-pipeline null (sample n,p and δ): hits at LENG accuracy are rare."""
        scores = full_pipeline_null_scores(N_MONTECARLO)
        from leng_test_utils import DELTA
        obs_score = match_score_mass(DELTA)
        hits = int(np.sum(scores <= obs_score))
        hit_rate = hits / N_MONTECARLO
        # Full pipeline: P(hit) = (1/96)*p_single; with p_single ~ 0, expect hit_rate < 0.01
        assert hit_rate < 0.02, f"Full-pipeline null: hit_rate={hit_rate:.4f} should be < 2%"


class TestOmnibusLandscapeLEE:
    """Unified test: landscape (n,p) selection + LEE δ-scan in one pipeline.
    Combines test_universe_landscape ablation with LEE Monte Carlo."""

    def test_landscape_then_lee(self):
        """Full pipeline: (n,p) sieve → δ scan on survivor → p_LEE."""
        from leng_test_utils import (
            is_self_consistent, all_positive_universe, is_prime, koide_p, DELTA,
        )
        # Step 1: Landscape sieve (deterministic)
        survivors = []
        for n in range(1, 9):
            for p in range(2, 14):
                if is_self_consistent(n, p) and all_positive_universe(n, p):
                    if is_prime(p) and abs(koide_p(p) - 1.0) > 1e-9:
                        survivors.append((n, p))
        assert survivors == [(3, 3)], f"Landscape sieve should give unique survivor, got {survivors}"

        # Step 2: LEE on the survivor's δ
        obs_score = match_score_mass(DELTA)
        assert obs_score < 0.0001  # < 0.01% max relative error

        # Step 3: Monte Carlo null (random δ for the survivor geometry)
        scores = best_match_under_null_mass(10_000, 0.0, 2 * np.pi, seed=99)
        p_single = np.mean(scores <= obs_score)
        assert p_single < 0.05  # < 5% of random δ match this well

    def test_permutation_scramble(self):
        """Negative control: scramble mass assignments and check score degrades."""
        from leng_test_utils import DELTA
        rng = np.random.default_rng(77)
        obs_score = match_score_mass(DELTA)

        # Scramble: assign predicted masses to wrong leptons
        mu = mu_from_electron(M_E_PDG, R, DELTA)
        masses = np.sort(brannen_masses(mu, R, DELTA))
        n_worse = 0
        for _ in range(1000):
            perm = rng.permutation(masses)
            rel_mu = abs(perm[1] - M_MU_PDG) / M_MU_PDG
            rel_tau = abs(perm[2] - M_TAU_PDG) / M_TAU_PDG
            scramble_score = max(rel_mu, rel_tau)
            if scramble_score > obs_score:
                n_worse += 1
        # Most permutations should be worse than the correct assignment
        assert n_worse / 1000 > 0.6


def lee_report() -> dict:
    """
    Full statistical report for Supplement C. Returns dict with:
    - hits: k hits out of M draws (mass prediction accuracy ≤ observed)
    - M: number of Monte Carlo trials
    - p_single: estimated null hit rate
    - ci_95: (lo, hi) Wilson 95% CI on p_single
    - seed: random seed used
    - delta_range: (lo, hi) rad (broad prior [0, 2π])
    - p_LEE: corrected p-value for N=96
    - observed_score: LENG mass prediction score S(δ)
    - score_note: explains the score metric
    """
    from leng_test_utils import DELTA
    scores = best_match_under_null_mass(N_MONTECARLO, DELTA_RANGE[0], DELTA_RANGE[1])
    obs_score = match_score_mass(DELTA)
    hits = int(np.sum(scores <= obs_score))
    p_single = hits / N_MONTECARLO
    lo, hi = wilson_ci(hits, N_MONTECARLO)
    p_lee = p_LEE_analytic(p_single, N_CANDIDATES)
    # Conservative: use Wilson upper bound for p_single when reporting p_LEE
    p_lee_upper = p_LEE_analytic(hi, N_CANDIDATES)
    return {
        "hits": hits,
        "M": N_MONTECARLO,
        "p_single": p_single,
        "ci_95": (lo, hi),
        "seed": RANDOM_SEED,
        "delta_range": DELTA_RANGE,
        "delta_justification": "Broad prior [0, 2π]; most conservative against the framework",
        "N": N_CANDIDATES,
        "p_LEE": p_lee,
        "p_LEE_upper": p_lee_upper,
        "observed_score": obs_score,
        "score_note": "max relative error of (m_mu, m_tau) vs PDG; K=2/3 always for r=sqrt(2)",
    }
