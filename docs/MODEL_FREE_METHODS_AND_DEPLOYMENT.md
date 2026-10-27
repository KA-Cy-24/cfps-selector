> **Current release (2026-10-27):** The canonical data-source contract is in [`INPUT_DATA_AND_REPRODUCIBILITY.md`](INPUT_DATA_AND_REPRODUCIBILITY.md). The notebook supports a frozen GP-informed export or read-only PocketBase input. A PocketBase-only input does **not** reproduce the GP-informed baseline or enable hybrid mode. For contributor attribution, see [`CONTRIBUTIONS.md`](CONTRIBUTIONS.md). Older notes below are historical methodology context and must not override this source contract.

# CFPS pre-run selection: independent model-free pipeline (v1.1.0)

## Scientific boundary and scope

This computational package combines a **read-only GP-informed experiment-selection baseline** with **independent model-free coverage optimization**. The former imports an existing proposal of first/check/sweep/independent experiments and associated synthetic model sensitivities. The latter applies exact deduplication, chemistry-vector representations based on dispensing volumes, matched one-factor comparisons, and facility-location/Maximin selection to **actual** submitted formulations.

Neither framework incorporates new Round 1 outcomes or experimentally measured efficacy. The original fitted Gaussian process and its virtual-LHS samples are not reconstructed here. For named research contributions see [`CONTRIBUTIONS.md`](CONTRIBUTIONS.md).

## Input data and assumed experimental unit

Input source: an authorized, frozen GP-informed `clusters_designs.json` (one record per submitted design, original ID and `recipeNl` in nL), **or** an independently fetched/read-only PocketBase `cfps_designs` snapshot after strict normalization. PocketBase-only snapshots support **model-free mode only** and do not contain historical GP baseline role assignments. Optional `clusters.json.gz` only for validation of the **display** expected-link matrix. External source: ZIP containing plate-level `samples` with `sample_type="experimental"` and `reagent_list`. Supports both single plates and `{"plates":[...]}` batches. Protein is extracted from plate's `protein_target`, `plate_id`, or filename; unrecognizable target raises an error. External designs are **not appended** to the participant pool.

For each protein, canonical dispensing signature is

\[
\sigma(x)=\operatorname{sorted}\{(j,\mathrm{round}(v_j,9)):j\ne\mathrm{water},\ v_j>0\}.
\]

This is an exact *stock dispensing* equivalence test, **not verified chemical equivalence** across stock concentrations and lot numbers. Zero and missing nL are treated as the same volume; the missing/zero semantics must be changed if future exports distinguish unknown values. The fluid-balance reagent `nuclease_free_water` is ignored. Duplicated students' recipe records are grouped as one *physical formulation*, but all source IDs and author attribution are preserved in private output.

**Working quota 320 per target** is an assumption, not an independently verified final rule. There are no new idealized LHS designs in this pipeline.

## External deduplication

1. Hash each experimental external recipe signature for its protein and count unique values, regardless of repeated plates or replicate factors.
2. Flag exact signatures that appear in participant records. Report source plate, sample ID privately.
3. Default `external_policy="flag"`: preserve the original eligible selection and show overlap.
4. Explicit `external_policy="exclude"`: remove exact matches from the selectable model-free pool. Hybrid fails with a clear error if an externally matching recipe was already locked by the GP baseline; this prevents silently rewriting the original GP proposal.
5. **Optional descriptive near-overlap distance**, NOT a deduplication criterion: append external unique recipes to the same per-target feature scale, calculate \(\min_{e\in E_p}d(x,e)\), and flag \(\le0.20\) in *normalized feature units*. This threshold is purely heuristic and only for interactive inspection. No automatic filtering, chemical-equivalence claims, or statistical testing based on it.

## Mathematical feature construction

Within each protein, for reagent \(j\), let \(v_{ij}\) denote stock-dispensing nL, and remove invariant reagent columns (including fixed template and lysate if constant). For the default normalization

\[
s_j=Q_{0.95}\{v_{ij}:v_{ij}>0\},\qquad
z_{ij}=\min(v_{ij}/s_j,1),\qquad
b_{ij}=\mathbf{1}[v_{ij}>0].
\]

The data-driven 95th percentile clips unusually large values to reduce leverage; this is a methodological choice, not a physical chemistry metric. An alternative uses per-reagent maximum. Three modes are supported: `scaled_volume`, `presence`, `scaled_plus_presence` (default):

\[
\phi(x_i)=(z_{i1},\dots,z_{ip},\sqrt{\alpha}\,b_{i1},\dots,\sqrt{\alpha}\,b_{ip}),\quad \alpha=0.35.
\]

Distance is \(d(x_i,x_k)=\lVert\phi(x_i)-\phi(x_k)\rVert_2\). No final molar concentrations are inferred: these require verified stock molarities. The distance reflects **composition geometry**, not predicted output covariance. Pairwise Pearson/Spearman *across reagent features* in the bundled separate correlation module is an optional sensitivity analysis; it is never conflated with GP prediction-error covariance.

### Exact contrasts and sweep reservation

A matched pair \((x_i,x_k)\) for factor \(j\) satisfies equality in every other non-water factor, and \(v_{ij}\ne v_{kj}\). A sweep family requires at least 3 **distinct** levels of \(j\) in a common unchanged non-water background. The algorithm inventories all such families, ranks by number of levels, then span and deterministic lexical keys. It prioritizes endpoints and distributed intermediate levels (up to 5 per family), capping families per factor at 5.

New sweep slots are bounded by BOTH \(\lfloor 0.20n\rfloor\) and \(\lfloor 0.50(n-|S_0|)\rfloor\) where \(S_0\) is the locked baseline/mandatory seeds. This prevents tiny hybrid gaps being completely consumed by sweep additions, leaving room for geometric coverage. These are **configurable heuristics**, not optimal budget allocations; one physical recipe can contribute to multiple pairwise contrasts.

### Geometric selection

Let \(X\) be the distinct participant formulations (never external rows) and \(S\subseteq X\) the currently selected physical recipes. Coverage is evaluated against all \(X\), weighting **each distinct formulation once**, rather than giving more influence to designs repeatedly submitted by a student.

- Greedy facility location: choose the candidate with largest incremental reduction of \(\sum_{x\in X}\min_{s\in S}d(x,s)\).
- Greedy Maximin: choose the candidate with largest \(\min_{s\in S}d(x,s)\). Without seeds, facility begins from the discrete medoid; Maximin begins from a maximally separated extreme, deterministically.

Both optimize over original submitted recipes only, excluding already selected canonical signatures. Ties break by lexical design ID. Neither optimizes yield nor chooses a causal experimental effect. Facility location tends to represent the observed candidate distribution; Maximin may favor edge designs.

### Diagnostics, not direct optimization scores

For each ingredient use observed **distinct positive** candidate nL levels to construct five quantile bins (collapsing tied edges) and represent zero as a separate bin if present. Report proportion of bins represented in \(S\). Report factor-pair bin occupancy for the eight factors with most occupied bins; denominator is **actually observed bin combinations**, never the impossible full Cartesian product. Display mean/90th percentile/maximum distance from each unique candidate to the selected set, count exact selected one-factor contrast pairs, representation of distinct author names, externally matching recipes and selected sweeps. A deterministic PCA scatter is a **visual projection only**, not the optimizing space.

Do not rank methods using one universally optimal threshold or present these metrics as independently learned biological effects.

## Three selection modes

| Mode | Locks GP proposal? | Uses model forecasts, links, or GP sensitivity to **choose** recipes? | Uses model-free data? |
|---|---|---|---|
| `gp_model_baseline` | Returns existing distinct proposed runs, **without padding** | Yes, inherited from original strategy, not re-run here | For audit/reporting only |
| `hybrid` | Yes, all GP-selected physical recipes | Only through that **fixed baseline**, not in geometric additions | Yes, exact sweeps + coverage |
| `model_free` | No | **No**. Model annotations may appear in viewer for post-hoc comparison only. | Yes, exact sweeps + coverage |

`ModelFreeConfig` controls `objective`, `feature_mode`, `distance_scale`, `presence_weight`, `max_contrast_fraction`, `max_sweep_share_of_new_slots`, `factor_bins`, `external_policy` and maximum sweep sizes. Explicit `seeds` allow retaining manual controls even in model-free mode. Any scientific control requirements/replicate counts should be set before final selection; the provisional comparison run has not yet finalized them.

## How to reproduce this release

Python dependency versions and checksums are captured in source manifests. Recommended: Python 3.11+, NumPy >=1.24, SciPy >=1.10, pandas >=2, nbformat >=5.7. For exact archival comparability, record precise `pip freeze`/`python --version` on the executing host.

```bash
# In a *fresh* local checkout of RonanChance/htgaa-automation, after adding the source files:
python -m pip install -r requirements-analysis.txt
PYTHONPATH=. python -m pytest -q tests

# Follow Ronan's authenticated workflow to obtain data/cfps-clusters/ (do NOT commit).
# E.g. npm run clusters:fetch -- --token=<agent token> when supported in that checkout.
# Keep agent token out of shell history, logs, notebooks and tracked files.

OPENBLAS_NUM_THREADS=1 PYTHONPATH=. python tools/cfps_model_free_cli.py \
  --designs data/cfps-clusters/clusters_designs.json \
  --links data/cfps-clusters/clusters.json.gz \
  --external-zip /private/path/cfps-plates.zip \
  --outdir data/cfps-model-free --quota 320 --external-policy flag
```

If Ronan's download command is not yet on your branch, obtain the authorized JSON through his supported process and supply its **local file path**. The pipeline does **no HTTP fetch** and never stores tokens. `--external-zip` and `--links` are optional; only the `--designs` file is required.

Standalone Colab: choose `INPUT_SOURCE = "published_snapshot"`, `"pocketbase_snapshot"`, or `"pocketbase_live"` in the first configuration cell. Upload the appropriate private JSON when not using live mode. Run all cells and download the generated **private results ZIP**. The notebook embeds verbatim sources; it does **not** import the mirrored Python modules at runtime. Default output under `/content/cfps_private_snapshots/model_free_analysis/`.

### Expected snapshot regression anchors

For the original uploaded participant JSON: design records sfGFP 549, PETase 500, Reteplase 551. Unique physical formulations: 484 / 459 / 488. Original GP-run baseline unique proposed recipes: 318 / 284 / 179. External ZIP: 2,880 experimental row records per protein across multiple plates, unique external formulations 717 / 717 / 714, exact internal cross-pool overlaps 0 / 0 / 1. These are **input-specific** anchor values, not universal properties of the competition.

## Deploy **only additive code** to the public repository

The repo's existing `src/routes/cfps/clusters/` remains intact. `repo_additions/` holds only NEW source paths:

```text
src/routes/cfps/model-free/+page.server.js
src/routes/cfps/model-free/+page.svelte
src/routes/cfps/model-free-data/+server.js
static/cfps/model-free/model_free_injector.js
static/cfps/model-free/model_free.css
data/cfps-model-free/.gitignore
```

Paths are relative to the **root of the current checkout**. `src/routes/cfps/model-free/+page.server.js` checks the same Keycloak cookie as the existing cluster route. The separate data endpoint checks `requestUser(event)` on **every** fetch (cookie or valid agent token), then reads only `data/cfps-model-free/model_free_analysis.json`. The server sends `Cache-Control: private, no-store`. The private file belongs in the ignored directory, NEVER under `static/` and NEVER in Git; no participant-level inputs or outputs are shipped in the source bundle. Original `/cfps/clusters` and source assets are not touched. The site references the published GP-informed public `clusters.css` plus narrowly scoped supplementary CSS.

**Important:** The public site has a rewritten Git history. Use a clean fresh clone and copy only these additive files. Do NOT run `git reset --hard` on a working directory containing uncommitted work. Inspect `git status --short` and explicitly stage only these source paths; verify `git check-ignore data/cfps-model-free/model_free_analysis.json` and `git diff --cached --name-only` before publishing. Do not share Keycloak credentials or agent tokens in code or logs.

For a local browser test without deployment, open `viewer/model_free_offline_preview.html` and use its file picker to load the private JSON; this reads the file locally, with no upload. An integration smoke test on Ronan's actual SvelteKit checkout is required before claiming successful production deployment.

## Remaining limitations (publication transparency)

- GP training, Ginkgo plate inference, actual covariance and LHS model sensitivities are not independently recreated.
- No verified stock concentration metadata -> the feature space is dispensing nL, not molarity, osmolarity, pH or biochemical activity; clipping and feature weights are sensitivity choices.
- External recipe equality requires identical reagent keys and volumes. External records may contain unknown stock provenance; exact signature only flags an **apparent** match.
- Site route and JS syntax are validated from inspected public source but have not been integrated into an authenticated deployed SvelteKit build; server authentication must be independently tested.
- Quotas, mandatory controls/replicates, plate constraints, allowed designs, author fairness and final submission validity remain for the separate final decision-framework phase.
