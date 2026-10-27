## Plate Decision Lab — new mixed-strategy decision framework

The authoritative integration instructions are in **[`COMMIT_GUIDE.txt`](COMMIT_GUIDE.txt)** and the complete mathematics, provenance, limitations and quota specifications in **[`docs/DECISION_FRAMEWORK.md`](docs/DECISION_FRAMEWORK.md)**. The new, self-contained **[`notebooks/cfps_round1_decision_framework_colab.ipynb`](notebooks/cfps_round1_decision_framework_colab.ipynb)** prepares private dataset inputs and a downloadable interactive viewer. The Python implementation is `analysis/cfps_selection/decision_framework.py`; the additive authenticated website is `/cfps/decision` with a protected `/cfps/decision-data` endpoint and plain-JS offline preview in `viewer/CFPS_Decision_Lab_Standalone.html`.

**Critical:** Ginkgo example plate files reserve columns 1–4 of 24 and specify four technical replicates: **320 physical experimental positions represent 80 distinct formulations per plate**, not 320 distinct submissions. Replication and GP share are configurable. The original GP `clusters_designs.json` has published analytical group/sweep lists, **not participant `batch_id` metadata**; PocketBase normalization now preserves explicit participant batch metadata for source-accurate warnings. Independent correlation clustering is *not equivalent* to GP shared-error groups; direct GP published designations are a first-class option.

# CFPS Round 1 — Reproducible Candidate-Space Analysis

**Release type:** source-only computational research package. This project analyzes **existing submitted cell-free protein synthesis (CFPS) formulations** for sfGFP, PETase, and Reteplase. It supports an independent model-free design-space coverage analysis and a read-only, GP-informed baseline comparison when the original annotated export is available. Its outputs are **provisional experiment portfolios**, not assayed response measurements or approved final submissions.

## Canonical entry points

- `notebooks/cfps_round1_full_model_free_colab.ipynb` — **self-contained Google Colab notebook** with source selection, read-only database fetch, input normalization, reproducible analysis, and private ZIP export; does **not** import the mirrored modules.
- `analysis/cfps_selection/input_sources.py` — validated conversion of PocketBase `cfps_designs` records into a formulation-only candidate pool.
- `analysis/cfps_selection/selection_framework.py` — audit of the published GP-informed design export, including check and sweep comparisons.
- `analysis/cfps_selection/model_free_pipeline.py` — exact duplicate/overlap screening, single-factor sweep detection, normalized-volume representation, Maximin/facility-location selection, and private visualization export.
- `analysis/cfps_correlation/pairwise_correlation_clustering.py` — independent, descriptive reagent-vector correlation clustering (not GP output-error covariance).
- `tools/cfps_model_free_cli.py` — deterministic offline equivalent of the computational pipeline.
- `repo_additions/` — exclusively additive authenticated visualization route and scoped assets for `RonanChance/htgaa-automation`.
- `viewer/CFPS_Model_Free_Standalone_Viewer.html` — standalone browser viewer; opens an authorized, private result file locally.

- `notebooks/cfps_gp_export_pocketbase_reconciliation_colab.ipynb` — read-only GP-export versus PocketBase integrity comparison; produces private discrepancy report and a verified GP copy only on exact experimental-pool parity. See `docs/GP_EXPORT_RECONCILIATION.md`.

## Select the analysis data source

Change **one setting** in the notebook's configuration cell:

```python
INPUT_SOURCE = "published_snapshot"    # default: GP-informed clusters_designs.json
# INPUT_SOURCE = "pocketbase_snapshot" # a downloaded read-only pocketbase_snapshot.json
# INPUT_SOURCE = "pocketbase_live"     # authorized live fetch via Colab Secrets
```

The default reads a frozen authorized `clusters_designs.json` export. Live PocketBase uses `PB_URL=https://db.rcdonovan.com` by default and requests authorized database credentials; An authorized `PB_TOKEN` (PocketBase bearer token) may be provided through Colab Secrets; otherwise `PB_EMAIL`/`PB_PASSWORD` are requested securely. Keycloak web or agent tokens are **not** assumed to be PocketBase authentication tokens. No credentials are required for either saved-file mode.

**Scientific separation:** The published annotated snapshot can be analyzed in `gp_model_baseline`, `hybrid`, and `model_free` modes. PocketBase-only snapshots lack GP-predicted responses, model-derived groups, and baseline run-role assignments: they are analyzed in `model_free` mode **only**. The pipeline never invents those fields or silently substitutes a historical baseline for live records.

Optional `cfps-plates.zip` is used exclusively for external exact-duplicate screening, **not** as eligible participant submissions or predictive training data. The optional `clusters.json.gz` supports a display-link diagnostic on the published snapshot; it does not train a GP.

## Offline CLI

```bash
python -m pip install -r requirements-analysis.txt
python -m pytest -q tests
PYTHONPATH=. python tools/cfps_model_free_cli.py \
  --source published_snapshot --designs /private/clusters_designs.json \
  --external-zip /private/cfps-plates.zip \
  --outdir /private/cfps-selection-results --quota 320 --external-policy flag
```

To analyze a downloaded PocketBase export instead:

```bash
PYTHONPATH=. python tools/cfps_model_free_cli.py \
  --source pocketbase_snapshot --designs /private/pocketbase_snapshot.json \
  --outdir /private/cfps-selection-live --quota 320
```

All method assumptions, frozen provenance hashes, output schemas, validation conditions, and known limitations are documented in `docs/INPUT_DATA_AND_REPRODUCIBILITY.md` and `docs/METHODOLOGY_AND_REPRODUCIBILITY.md`.

## Privacy and deployment

Source distributions contain **no participant-level recipes, identifiers, access tokens, private visualization JSON, or external plate records**. Source-only additions do not modify existing GP cluster route files. Protected outputs must live in an ignored server-side directory, not under `static/`; website data access requires server-side authorization. Validate against the current application checkout before deploying.

## Research contributions

Authorship and contributions are documented once in `docs/CONTRIBUTIONS.md` and in the notebook's attribution section. Public method names, APIs, routes, and visualization labels use scientific terminology only.
