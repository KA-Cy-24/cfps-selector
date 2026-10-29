# CFPS Selector: A mutli-strategy decision-making framework for selecting human-designed cell-free protein synthesis (CFPS) reaction formulations

As a part of How To Biomanufacture Almost Anything 2026 ([HTBAA 2026](https://2026b.htgaa.org/)), a global online course led by MIT Media Lab's Community Biotechnology Initiative, a *coopetition* is in progress (as of October 2026) between the community's group of human researchers (totaling to around 4 teams) and a group of LLM agents used by OpenAI to optimize reaction compositions for *Escherichia coli*-based cell-free systems. The overarching goal of the coopetition to is achieve the optimal yield/cost balance (i.e. the **specific cost** or $/g) for *E. coli* for three protein targets and across different scales. Nebula, an Autonomous Laboratory based on reconfigurable automation units and instruments developed by Ginkgo Bioworks, serves as the main testing platform. This coopetition is to be run in three rounds, with each comprising of 384 well-plates. For each of the rounds, not all of the submissions made by human teams could be included due to the availability of only 320 samples (or 80 unique compositions excluding the technical replicates) per protein target. For more details and supplementary information, please refer to the main repository hosting the coopetition ([RonanChance/htgaa-automation](https://github.com/RonanChance/htgaa-automation)).

The HTBAA community alone, which constitutes one of the human research groups, has produced compositions almost double number of free wells available for experimental validation in Round 1. Hence, the community has recognized the need for rational yet fair inclusion of submissions made by teaching assistants and students of the course. The present repository, consisting of both a Jupyter/Python-based framework as well as a graphical web application, constitutes part of the effort to address this need and a contribution towards the preparation of the best set of test samples for Round 1 and beyond.

## Integration Instructions/Metadata and Future Additions 

The authoritative integration instructions are in **[`COMMIT_GUIDE.txt`](COMMIT_GUIDE.txt)** and the complete mathematics, provenance, limitations and quota specifications in **[`docs/DECISION_FRAMEWORK.md`](docs/DECISION_FRAMEWORK.md)**. The new, self-contained **[`notebooks/cfps_round1_decision_framework_colab.ipynb`](notebooks/cfps_round1_decision_framework_colab.ipynb)** prepares private dataset inputs and a downloadable interactive viewer. The Python implementation is `analysis/cfps_selection/decision_framework.py`; the additive authenticated website is `/cfps/decision` with a protected `/cfps/decision-data` endpoint and plain-JS offline preview in `viewer/CFPS_Decision_Lab_Standalone.html`.

**Critical:** Ginkgo example plate files reserve columns 1 to 4 of 24 and specify four technical replicates: **320 physical experimental positions represent 80 distinct formulations per plate**, not 320 distinct submissions. Replication and GP share are configurable. The original GP `clusters_designs.json` has published analytical group/sweep lists, **not participant `batch_id` metadata**; PocketBase normalization now preserves explicit participant batch metadata for source-accurate warnings. Independent correlation clustering is *not equivalent* to GP shared-error groups; direct GP published designations are a first-class option.

Features that are to be added in the near future are:
1. Individual fairness rules and batch-inclusion weights integrated into the graphical user interface and the codebase
2. Support for multiple models (other than the Gaussian Process model)
3. Support for Latin Hypercube Sampling as another method of model-free selection of samples
4. Better integration of cluster-based and rule-based selection of samples
5. Support for more than two selection methods as modularly-added weighted terms (e.g. GP-informed sampling, Random Forest-informed sampling, model-free sampling using LHS, etc.)

## CFPS Round 1: Reproducible Candidate-Space Analysis

This project analyzes **existing submitted cell-free protein synthesis (CFPS) formulations** for sfGFP, PETase, and Reteplase. It supports an independent model-free design-space coverage analysis and a read-only, GP-informed baseline comparison when the original annotated export is available. Its outputs are **provisional experiment portfolios**, not assayed response measurements or approved final submissions.

### Canonical entry points

- `notebooks/cfps_round1_full_model_free_colab.ipynb` — **self-contained Google Colab notebook** with source selection, read-only database fetch, input normalization, reproducible analysis, and private ZIP export; does **not** import the mirrored modules.
- `analysis/cfps_selection/input_sources.py` — validated conversion of PocketBase `cfps_designs` records into a formulation-only candidate pool.
- `analysis/cfps_selection/selection_framework.py` — audit of the published GP-informed design export, including check and sweep comparisons.
- `analysis/cfps_selection/model_free_pipeline.py` — exact duplicate/overlap screening, single-factor sweep detection, normalized-volume representation, Maximin/facility-location selection, and private visualization export.
- `analysis/cfps_correlation/pairwise_correlation_clustering.py` — independent, descriptive reagent-vector correlation clustering (not GP output-error covariance).
- `tools/cfps_model_free_cli.py` — deterministic offline equivalent of the computational pipeline.
- `repo_additions/` — exclusively additive authenticated visualization route and scoped assets for `RonanChance/htgaa-automation`.
- `viewer/CFPS_Model_Free_Standalone_Viewer.html` — standalone browser viewer; opens an authorized, private result file locally.

- `notebooks/cfps_gp_export_pocketbase_reconciliation_colab.ipynb` — read-only GP-export versus PocketBase integrity comparison; produces private discrepancy report and a verified GP copy only on exact experimental-pool parity. See `docs/GP_EXPORT_RECONCILIATION.md`.

### Select the analysis data source

Change **one setting** in the notebook's configuration cell:

```python
INPUT_SOURCE = "published_snapshot"    # default: GP-informed clusters_designs.json
# INPUT_SOURCE = "pocketbase_snapshot" # a downloaded read-only pocketbase_snapshot.json
# INPUT_SOURCE = "pocketbase_live"     # authorized live fetch via Colab Secrets
```

The default reads a frozen authorized `clusters_designs.json` export. Live PocketBase uses `PB_URL=https://db.rcdonovan.com` by default and requests authorized database credentials; An authorized `PB_TOKEN` (PocketBase bearer token) may be provided through Colab Secrets; otherwise `PB_EMAIL`/`PB_PASSWORD` are requested securely. Keycloak web or agent tokens are **not** assumed to be PocketBase authentication tokens. No credentials are required for either saved-file mode.

**Scientific separation:** The published annotated snapshot can be analyzed in `gp_model_baseline`, `hybrid`, and `model_free` modes. PocketBase-only snapshots lack GP-predicted responses, model-derived groups, and baseline run-role assignments: they are analyzed in `model_free` mode **only**. The pipeline never invents those fields or silently substitutes a historical baseline for live records.

Optional `cfps-plates.zip` is used exclusively for external exact-duplicate screening, **not** as eligible participant submissions or predictive training data. The optional `clusters.json.gz` supports a display-link diagnostic on the published snapshot; it does not train a GP.

### Offline CLI

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

### Privacy and deployment

Source distributions contain **no participant-level recipes, identifiers, access tokens, private visualization JSON, or external plate records**. Source-only additions do not modify existing GP cluster route files. Protected outputs must live in an ignored server-side directory, not under `static/`; website data access requires server-side authorization. Validate against the current application checkout before deploying.

### Research contributions

Authorship and contributions are documented once in `docs/CONTRIBUTIONS.md` and in the notebook's attribution section. Public method names, APIs, routes, and visualization labels use scientific terminology only.
