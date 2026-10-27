# Verification protocol

Run from a clean source checkout:

```bash
python -m pip install -r requirements-analysis.txt
python -m pytest -q tests
python -m compileall -q analysis tools tests
node --check repo_additions/static/cfps/model-free/model_free_injector.js
python tools/build_full_model_free_colab.py
python -m pytest -q tests/test_notebook_parity.py
```

The unit tests use synthetic records (not participant-level data). The suite validates canonical exact-recipe identity, external ZIP parsing, exact single-factor contrasts, deterministic candidate selection, unannotated PocketBase-source normalization, mode eligibility, and byte-for-byte embedded notebook/module parity.

For a locally authorized regression fixture, run the published export once with the CLI, store the original file SHA-256 and selected ID lists privately, and compare an independently executed Colab run. A PocketBase-only run may contain a later or different eligible candidate universe and **must not** be expected to reproduce a historical annotated export's portfolio.

Live PocketBase authentication and protected Keycloak deployment require authorized credentials and application-level integration testing; local notebook/CLI tests do not certify those services.


## Decision Lab source integration

Run `PYTHONPATH=. python -m pytest -q tests/test_decision_framework.py` and `node --check repo_additions/static/cfps/decision/decision.js`. The self-contained Colab notebook is `notebooks/cfps_round1_decision_framework_colab.ipynb`; it must be regenerated after mirrored-module changes. The source-layout evidence, scientific math and deployment restrictions are in `docs/DECISION_FRAMEWORK.md`.
