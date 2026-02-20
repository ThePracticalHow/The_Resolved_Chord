# Papers

LaTeX sources and generated outputs for The Resolved Chord.

## Structure

| Path | Contents |
|------|----------|
| **The_Resolved_Chord_v11.tex** | Main paper (72 predictions, zero free parameters) |
| **Supplement_\*.tex** | Supplements I–XIV (14 files: Geometry, Leptons, Gauge, Baryons, Higgs, Quark Flavor, Neutrinos, Adversarial Defense, Strange Castles, Spectral Action, Derivation Status, Companion Guide, Perfect LOTUS, Why FMA) |
| **math-papers/** | Standalone math papers (Eta/Ghost Identity, Spectral Exclusion, N1 Bridge, APS Three Generations) |
| **figures/** | Figures used by the paper |
| **pdf/** | Pre-compiled PDFs |
| **markdown/** | Markdown versions (generated from .tex via `tools/convert_tex_pdf_md.py`) |

> [!NOTE]
> Previous versions (v10) are archived in `05_Project_LENG/_archive/`.

## Generating Outputs

```bash
# From project root:
python tools/convert_tex_pdf_md.py   # .tex → .md + .pdf
```

## Key Files

- **The_Resolved_Chord_v11.tex** — Main paper (72 predictions, zero free parameters)
- **Supplement_X_SpectralAction** — Math-to-physics derivation chains
- **Supplement_XI_DerivationStatus** — Derivation status table
- **Supplement_XII_CompanionGuide** — Topic-based reference
