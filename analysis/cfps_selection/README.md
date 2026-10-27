# CFPS Selection Analysis — Python API

This package provides independent, read-only audits and provisional experiment-selection strategies for **existing CFPS submissions**.

## Modules

- `input_sources.py` — validate a published GP-informed export or normalize a complete PocketBase `cfps_designs` snapshot into compatible experimental formulation records.
- `selection_framework.py` — reproducible audit of an annotated GP-informed export: canonical recipe duplicates, control checks, matched comparisons, sweep completeness, model-link matrix diagnostics.
- `model_free_pipeline.py` — candidate-pool geometry, exact independent/external deduplication, exact single-factor contrast/sweep detection, facility-location/Maximin selection, protected ID-level JSON.

## Source and mode contracts

- **`published_snapshot`**: dated `clusters_designs.json`, baseline annotations already present. All three modes available: `gp_model_baseline`, `hybrid`, and `model_free`.
- **`pocketbase_snapshot`**: fetched/read-only PocketBase `cfps_designs` records, normalized from `design_json` and `reagent_list`. Only `model_free` is available unless the GP-annotated data are independently supplied and aligned; the adapter never fabricates model outputs.

`canonical_recipe` excludes balancing water for exact dispensing-recipe identity. A zero and absent ingredient count as equal. Candidate geometry uses normalized nL dispensing volumes and optionally ingredient presence, never assumed molarity. Each distinct candidate formulation receives equal coverage weight, regardless of duplicate submitted records.

## Commands

```bash
python tools/cfps_model_free_cli.py --source published_snapshot \
  --designs /private/clusters_designs.json --outdir /private/results
python tools/cfps_model_free_cli.py --source pocketbase_snapshot \
  --designs /private/pocketbase_snapshot.json --outdir /private/live_results
```

The CLI never performs network access. The self-contained Colab notebook optionally fetches the PocketBase snapshot in read-only mode. See `docs/INPUT_DATA_AND_REPRODUCIBILITY.md` for precise provenance, limitations, and comparison rules.
