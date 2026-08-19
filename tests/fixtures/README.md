# KGO test fixtures

Binary reference outputs for regression tests. Regenerate after intentional physics or data changes:

```bash
python3 scripts/generate_kgo_reference.py
```

Then run `pytest tests/` and commit updated `.npz` files if outputs changed as expected.

Fixtures include **non-flat** spectra: energy power laws (E⁻⁷, E⁻²) and rigidity power law (R⁻⁷) per species.
