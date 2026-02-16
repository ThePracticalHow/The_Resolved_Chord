# The Resolved Chord

**The Spectral Geometry of Everything: Standard Model, Gravity, and Cosmology from S⁵/Z₃**

*Jixiang Leng, February 2026*

## What this is

A single compact manifold — the lens space S⁵/Z₃ — determines all 26 parameters of the Standard Model, the gravitational coupling (0.10%), the cosmological constant (1.4%), and the broad features of the cosmological history, with zero free parameters.

The electron mass is the unit of measurement, not an input ([Duff, Okun, Veneziano 2002](https://arxiv.org/abs/physics/0110060)). Every dimensionless ratio is fixed by the spectral geometry.

## Key results

| Prediction | Precision | Status | Verification script |
|-----------|-----------|--------|---------------------|
| Koide ratio K = 2/3 | Exact | Theorem | `leng_replication.py` |
| Generations N_g = 3 | Exact | Theorem | `lotus_aps_generation.py` |
| Proton mass 6π⁵ | 10⁻¹¹ | Theorem (Parseval) | `ghost_parseval_proof.py` |
| Fine-structure constant | 0.001% | Derived | `alpha_from_spectral_geometry.py` |
| Higgs VEV | 0.005% | Derived | `vev_overlap.py` |
| Higgs mass | 0.034% | Derived | `spectral_action_derivation.py` |
| All 6 quark masses | RMS 0.111% | Theorem | `downtype_spectral_ordering.py` |
| Planck mass | 0.10% | Theorem (5-lock) | `gravity_theorem_proof.py` |
| Cosmological constant | 1.4% | Derived | `cc_aps_proof.py` |
| Strong CP θ̄ = 0 | Exact | Theorem | — |

30 predictions. 11 at Theorem level. Zero free parameters.

## Quick verification

```bash
pip install numpy scipy mpmath
python verification/EtaInvariant.py          # eta = 2/9
python verification/ghost_parseval_proof.py   # proton mass to 50 digits
python verification/gravity_theorem_proof.py  # Planck mass (16/16 checks)
python verification/cc_aps_proof.py           # cosmological constant
python -m pytest falsification/ -v            # full test suite
```

## Paper

- **paper/The_Resolved_Chord_v9.tex** — Main paper (18 sections, 1615 lines)
- **supplements/** — Supplements I–XI (complete proof chains)
- **math-papers/** — Four standalone math papers

## How to read this

1. Start with the abstract in **paper/The_Resolved_Chord_v9.tex**
2. If you doubt a claim, check **supplements/Supplement_XI_DerivationStatus.tex** — it lists every claim with its status, proof location, and the response to the strongest skeptical objection
3. Run the verification script for any prediction you want to check
4. Run `pytest falsification/` for the full adversarial test suite

## Citation

```bibtex
@article{leng2026resolved,
  author = {Leng, Jixiang},
  title = {The Resolved Chord: The Spectral Geometry of Everything},
  year = {2026},
  note = {arXiv:XXXX.XXXXX}
}
```

## License

MIT
