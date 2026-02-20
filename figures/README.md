# Figures

Generate figures for The Resolved Chord. Run from the repo root (`public-release/`):

```bash
python figures/prediction_scatter.py   # all 87 predictions vs measurement (live from LOTUS)
python figures/roadmap_figure.py       # From Tr(f(D²)) to all 87 predictions
python figures/equations_derivation.py # 7 classical equations from one trace
python figures/lotus_potential.py      # LOTUS potential V(φ) landscape
python figures/lotus_emblem.py         # The LOTUS emblem
python figures/bare_orbifold.py        # S⁵/Z₃ orbifold cross-section
python figures/five_invariants.py      # Five spectral invariants on the orbifold
python figures/spectral_cascade.py     # Spectral cascade derivation levels
python figures/lepton_simplex.py       # Lepton mass simplex
```

Outputs: `prediction_scatter.png`, `roadmap_figure.png`, `equations_derivation.png`, `lotus_potential_figure.png`, `lotus_emblem.png`, `bare_orbifold.png`, `five_invariants.png`, `spectral_cascade.png`, `lepton_simplex.png`

These `.png` files are also copied into `papers/figures/` so the LaTeX paper can reference them directly when uploaded to Overleaf.

Requires: `numpy`, `matplotlib`, `lotus` (for prediction_scatter.py)
