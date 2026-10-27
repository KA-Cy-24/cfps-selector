# Validation record — publication-source candidate v1.1.0

**Evaluation date:** 27 October 2026. Tests were run with local authorized snapshots and synthetic fixtures; no production PocketBase or Keycloak credentials were used.

## Automated results

- **43 synthetic/regression tests passed** (`python -m pytest -q tests`). These include canonical formulation identity, mode separation, complete/partial source validation, external exact overlap, and notebook/module source parity.
- All Python modules and test files compile; the additive JavaScript injector passes `node --check`.
- The canonical Colab notebook parses under `nbformat` and embeds `pairwise_correlation_clustering.py`, `selection_framework.py`, `input_sources.py`, and `model_free_pipeline.py` **verbatim** in executable cells.
- A headless execution of the self-contained notebook against a synthetic offline PocketBase snapshot completed source normalization, both permitted model-free objective variants, private visualization JSON generation, and a validated private ZIP export. No GP-informed baseline was fabricated.
- A full frozen-source CLI run reproduced the original benchmark: **318 sfGFP, 284 PETase and 179 Reteplase** unique GP-informed baseline formulations and **320 per target** for each hybrid and model-free objective variant (working quota assumption).
- A separate headless self-contained notebook run on that same frozen source, with external overlap ZIP, reproduced **identical selected IDs for all 15 protein/strategy combinations** against the CLI. This is an exact-result parity check, not independent GP retraining.

## Scientific and technical caveats

- Input-based regression checks are specific to the frozen files identified by SHA-256 in the private manifest; newer live database content is not assumed to be the same data universe.
- Live PocketBase authorization via bearer token or email/password has not been tested in production. Keycloak agent tokens are not assumed to work as PocketBase JWTs.
- The additional SvelteKit source is additive; **no deployment to the protected application, end-to-end Keycloak test, or remote repository modification has been performed**.
- Source-level validation does not establish accuracy of the original Gaussian-process forecasts, GP covariance, virtual Latin-hypercube sensitivity model, official quota, sample eligibility, or biological inferential value.
- The final portfolio remains provisional pending manual and rules-based controls/replicate constraints.

Full input/provenance conventions: `INPUT_DATA_AND_REPRODUCIBILITY.md`. Research attribution: `CONTRIBUTIONS.md`.
