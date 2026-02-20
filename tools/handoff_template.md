# AI Handoff Template

Use this template for Claude, Composer, Copilot, Gemini, Grok, Perplexity, ChatGPT, and other coding agents.

---

## 1) Files Changed

- Path: short reason
- Path: short reason

## 2) What Changed (Summary)

- One-line summary of the functional/doc impact.
- Note any policy/count/path updates.

## 3) Validation Run

Commands executed:

```bash
python tools/release_gate.py
python -m pytest tools/falsification/ -m "not slow" -q
```

Results:

- `release_gate`: PASS/FAIL (+ one-line issue summary if fail)
- `pytest`: PASS/FAIL (+ counts)

## 4) Open Decisions / Risks

- Decision needed: (owner, options, recommendation)
- Risk: (scope, impact, mitigation)

## 5) Next Owner Action

- Single highest-priority next step.

## 6) Structural Docs Touched?

- [ ] `V11_INVENTORY.md`
- [ ] `MASTER_CODE_INDEX.md`
- [ ] `CASTLE_CONSTRUCTION_PLANS.md`
- [ ] `public-release/README.md`
- [ ] `public-release/CONTRIBUTING.md`

If checked, confirm whether changes were append-only where required.
