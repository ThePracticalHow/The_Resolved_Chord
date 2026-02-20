# Falsification Test Suite

Adversarial pytest tests that stress the LOTUS framework. Run from project root:

```bash
python tools/falsification/run_falsification.py -v
# or
pytest tools/falsification/ -v
```

## Key Files

| File | Purpose |
|------|---------|
| **run_falsification.py** | Entry point; runs pytest and writes results |
| **leng_test_utils.py** | Shared constants, Koide formulas, PDG values |
| **dictionary_spec.py** | Frozen geometry-to-physics mapping spec |
| **conftest.py** | Pytest config; adds project paths |

## Test Files (by topic)

- **Core**: `test_universe.py`, `test_falsification.py`, `test_geometry.py`
- **Masses**: `test_masses.py`, `test_koide.py`, `test_eta.py`, `test_proton.py`
- **Neutrinos**: `test_neutrinos.py`
- **Couplings**: `test_alpha_s.py`
- **Controls**: `test_negative_controls.py`, `test_ghost_modes.py`, `test_look_elsewhere.py`
- **Scalar**: `test_scalar_95.py`, `test_permutation_scramble.py`
- **Spec**: `test_dictionary.py`, `test_replication.py`, `test_precision.py`
- **Integration**: `test_paper_claims.py`, `test_verification_scripts.py`, `test_universe_landscape.py`
