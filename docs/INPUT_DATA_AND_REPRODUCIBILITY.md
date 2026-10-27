# Input source contract, provenance, and reproducibility

## Input sources: exactly what is analyzed

| Mode | Raw source | Data used for selection | Available strategies | Authentication |
|---|---|---|---|---|
| `published_snapshot` | Authorized, frozen `clusters_designs.json` | 1,600 submitted participant records in the dated analysis snapshot: target, `id`, `recipeNl` (stock dispensing volume), available annotations | `gp_model_baseline`, `hybrid`, `model_free` | None for a locally held authorized file |
| `pocketbase_snapshot` | Authorized `pocketbase_snapshot.json` with `collections.cfps_designs` | Records of `sample_type="experimental"`, target from `design_json.metadata.protein_target`, volumes from `design_json.reagent_list`, immutable PocketBase record ID | `model_free` only | None after authorized export |
| `pocketbase_live` | Live read-only PocketBase query | Same as above, after a complete paginated fetch | `model_free` only | Authorized PocketBase account via optional PocketBase `PB_TOKEN` secret, otherwise authorized email/password credentials |

**Do not confuse a model annotation with a measurement.** The published export contains pre-run GP forecasts, LHS-based sensitivity annotations, run roles and group memberships. The PocketBase source adapter does not reconstruct any of those properties. It also does not combine records across snapshots by similar names or nearest neighbors.

The default notebook setup is:

```python
INPUT_SOURCE = "published_snapshot"
PUBLISHED_SNAPSHOT_PATH = "/content/clusters_designs.json"
POCKETBASE_SNAPSHOT_PATH = "/content/pocketbase_snapshot.json"
EXTERNAL_PLATES_ZIP = "/content/cfps-plates.zip"
AUTO_DOWNLOAD_PRIVATE_RESULTS = True
```

The database fetch is optional and uses `PB_URL=https://db.rcdonovan.com` by default, according to the October 27 site update. Credentials are never placed into output manifests. A Keycloak web-login/agent token is **not guaranteed** to be a valid PocketBase `_superusers` password login. Live authentication cannot be fully validated without authorized credentials.

## Input normalization

`analysis/cfps_selection/input_sources.py` validates identity, target, record schema, nonnegative finite reagent volumes, 20,000 nL well volume, and data completeness indicators. It explicitly counts non-experimental records excluded from the participant pool. The normalized snapshot is written privately with its own SHA-256; original source files remain immutable. Any source mismatch that prevents reliable target assignment causes an error, not heuristic recoding. The original GP-informed JSON passes through without invented modifications.

External ZIP: `cfps-plates.zip` comprises **other teams' plate recipes**. It is a comparator for exact matches of reagent key/nL volume/protein only. It is neither merged into the candidate pool nor used to train a predictive model. Excluding an external overlap is explicitly optional (`external_policy=exclude`), otherwise results are flagged. Absolute chemical equivalence is not established without verified stock metadata.

`clusters.json.gz` contains published UI group/link information. It is optional for the **published-source** audit and never used for geometry-based selection.

## Reproducible output bundle

The final notebook cell produces `cfps_private_run_YYYYMMDDTHHMMSSZ.zip` with:

- `source_manifest.json`: source mode, immutable SHA-256 checksums, input counts, adapter/algorithm versions, external comparison mode and quota assumption.
- `normalized_candidate_pool_PRIVATE.json`: input-adapted recipes with original immutable IDs. **Private.**
- `model_free_analysis.json`: ID-level visualizations and all supported selection strategies. **Private.**
- `selected_ids_PRIVATE.json`: exactly which original IDs are selected by each mode. **Private.**
- `strategy_comparison.csv`: per-target strategy-level coverage metrics.

A separate GP-informed baseline audit is produced only when a GP-annotated input source is used. All files must remain behind authentication; do not put private outputs into a public source repository, its `static/` directory, or public notebook outputs.

## Interpretation boundary

- The dataset describes **submitted recipes**. There are no new, hypothetical, LHS-generated candidate submissions.
- Model-free distances use reagent-wise scaled **dispensing volumes and ingredient-presence**, not inferred molarities or biological distance.
- Facility-location reduces average distance to the selected recipes; Maximin prioritizes maximal separation. Neither is expected information gain or activity optimization.
- Comparisons between snapshots are not direct algorithm comparisons unless candidate IDs, source dates, and recipes are aligned and recorded.
- `320` per protein is a working quota assumption; controls, official eligibility, final format, and external-data policy must be verified before submission.

## Self-contained notebook parity

`tools/build_full_model_free_colab.py` embeds the input adapter, published-snapshot auditor, correlation clusterer, and model-free analysis source verbatim. The notebook does not import those `.py` modules at runtime. `tests/test_notebook_parity.py` verifies equivalence. The Python CLI uses the same modules directly. Changing a module requires rebuilding and rerunning the notebook tests before publication.
