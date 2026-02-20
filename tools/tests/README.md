# Basic Tests

Quick smoke tests for LOTUS. Run from project root.

| Script | Purpose |
|--------|---------|
| **test_new_functionality.py** | LOTUS API smoke test: import, Universe, derive(), equations |
| **test_docker.py** | Docker container validation (run inside container via `make docker-test`) |

```bash
python tools/tests/test_new_functionality.py
```
