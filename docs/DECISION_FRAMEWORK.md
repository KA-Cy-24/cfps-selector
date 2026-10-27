# CFPS Experimental Plate Decision Framework — methodology and deployment specification

**Status:** pre-experimental decision-support; **not** an execution-ready Ginkgo plate file. The implementation chooses only submitted recipes and uses source-provided GP annotations, independent dispensing-volume geometry, or a stated mixture. The deployment must restrict participant-level data to authorized users.

## Inputs and provenance

1. `clusters_designs.json` (immutable published GP-annotated export): contains source `proteins[].designs[]`, original submission IDs, reagent dispensing volumes in nL, original `round1Role` annotations, GP forecasts and original `groups[].designs`/`firstRun`/`checkRuns` and `sweeps[].sweeps[].toRun`. **It does not contain `batch_id`, `batch_size`, or `batch_index` fields.** Accordingly, original participant-defined batching *cannot be validated from this export alone*. Its GP groups are analytical families, not submission batches.
2. PocketBase `cfps_designs` snapshot (optional): authorized fetch or static export. The input adapter preserves explicit `design_json.metadata.batch_id`, `batch_size`, `batch_index` in normalized records; this permits genuine participant-batch checks. It does **not** include GP forecasts, groups, or roles. The GP quota is therefore fixed to 0 unless a fully validated GP-annotated snapshot is independently available; none is silently manufactured.
3. The example Ginkgo result/plate files in `cfps-plates.zip` are **not part of the eligible submission pool**. Their plate metadata is used **solely to validate plate geometry**. Each examined 16×24 example declares `reserved_columns=[1,2,3,4]`, `replicate_factor=4` and 320 experimental sample positions, with 80 distinct `sample_id` values repeated four times in examined representative plates. This is consistent with the public [Ginkgo CFPS plate model](https://github.com/ginkgobioworks/ginkgo-automation-cfps). Actual run format and any mixed controls require approval from the experiment operator.

Inputs are hashed. The publishing ZIP contains source and tests only, **no private competition recipes or participant IDs**. The Colab exports a separate PRIVATE ZIP, deliberately unsuitable for Git.

## Operational constraints

The 384 physical wells are rows A–P and columns 1–24. The 64 wells in columns 1–4 are reserved; columns 5–24 furnish 320 available experimental positions. For configurable replicate factor \(r\) dividing 320, the maximal number of independent conditions is \(K=320/r\). At \(r=4\), **K=80 distinct submitted recipes per protein plate**. This is **not** the same as the earlier study's capacity of 320 unique selected recipes per target, which was a screening/ranking budget rather than a certified four-replicate physical plate plan.

The GP quota is either an exact multiple-of-r physical-well count or a percentage of the K unique experimental conditions, using half-up rounding. The model-free quota is the remaining conditions. For example at r=4, 50% yields 40 GP and 40 model-free submissions per plate; 160 GP physical wells plus 160 model-free physical wells. A GP-only plate uses a 100% share, and a model-free-only plate uses a 0% share. Over-capacity, duplicate manual assignments, forbidden identifiers and missing GP annotations fail closed. The website visually reassigns whole replicate groups, not single technical replicates.

**Important:** The displayed deterministic well order is a planning preview. Ginkgo's experimental execution software may randomize sample placement to mitigate spatial position effects; do not treat this preview as a certified Ginkgo run file.

## GP-informed scoring and original designations

Two explicitly different user-facing policies:

- `published_designations` (recommended for scientific fidelity): take original source `round1Role` priorities in order `first run`, `check run`, `sweep run`, `independent run`, with forecast-score tie breaking. When GP quota exceeds the available source-designated roles, extend it using forecast-ranked additional submissions and record all such IDs as `gp_forecast_extensions`. This is a **new decision policy using** the published labels, **not an exact reproduction of the original plate proposal**.
- `forecast_ranked`: rank every eligible candidate by a user-weighted normalized sum of ordinal ranks of `pBeat` (higher), `lnSpecCost` (lower), and `sd` (higher, exploration). All values are **model predictions**, not experimentally observed yields or causal effects. The score is \(G(x)=\frac{w_y R(pBeat_x)+w_c R(-\log C_x)+w_u R(sd_x)}{w_y+w_c+w_u}\). Equal-valued attributes contribute zero rather than invent rank separation. When all weights are zero, deterministic ID tie-breaking applies.

The `published_designations` route **does not recompute** GP link scores, the original clusters, or GP regression. It directly reads the original source's annotations. This is essential because the independent reagent-pattern correlation algorithm is **not** the GP shared-output-error cluster model.

## Model-free criteria

The unique candidate pool is deduplicated by *exact dispensing recipe*, ignoring balancing water and preserving source record aliases. Let \(v_{ij}\) be the submitted nL of reagent j. Each reagent's volume feature is \(z_{ij}=\min(1,v_{ij}/Q_{.95}(v_{\cdot j}>0))\), with a separate binary presence feature weighted by \(\sqrt{\lambda}\). Euclidean pairwise distance is normalized to the largest observed pairwise distance for the target. It is **not a molarity-calibrated or empirically validated biological response metric**.

For a currently selected set S and candidate pool X:

- `maximin`: select the remaining \(x\) maximizing \(m(x;S)=\min_{s\in S}d(x,s)\).
- `facility_location`: select \(x\) maximizing the coverage gain \(F(x;S)=|X|^{-1}\sum_{u\in X}\max(0,d(u,S)-d(u,x))\).
- `balanced`: use \(\alpha m(x;S)+(1-\alpha)F(x;S)\), with \(\alpha\) exposed in the interface.

If S is empty, the normalized initial nearest-selected distance is 1. All criteria use the same deterministic lexical-ID tie breaking; the exact GP and model-free quotas are maintained. The previous full model-free pipeline also investigates exact contrasts and sweeps; **this interactive greedy geometric tool is a narrower configurable decision policy**, and does not claim to reproduce that prior multi-criterion recipe-selection sequence.

## Source groups, candidate batches and incompleteness

The strict published-label integrity audit checks every original cluster's members and its `firstRun` and `checkRuns` against the source candidate IDs, and verifies per-design `group.letter` matches membership. On the frozen dataset, these checks pass for 86 sfGFP, 79 PETase and 54 Reteplase groups. **This is internal source consistency, not equivalence to an independent reconstruction.** Because the independently recreated correlation algorithm measures a different statistical object, exact equivalence **is not established** and should never be asserted.

Warnings are separated by mathematical meaning:

- **High:** a participant-defined PocketBase `metadata.batch_id` was partly selected, or the source declares more original batch members than present (when `batch_size` is provided). A warning does not automatically force a partial batch into the plate; that could violate quota.
- **Medium:** a published GP sweep family's `toRun` set is partly selected; interpretation may require the omitted dose levels.
- **Information:** an analytical GP cluster's `designs` set is partly selected. This is *expected* in subset sampling and not automatically a protocol violation.

Both machine exports and the UI present warnings. They cannot be removed by hiding a viewer panel; they remain included in the exported decision JSON. Batch completeness requires an authorized PocketBase snapshot because the frozen enriched export lacks original batch metadata. Nonexistent IDs are never fabricated.

## Outputs and privacy

The offline Colab and CLI generate (outside Git): `cfps_decision_input_private.json` (source recipes for interactive recomputation), `decision_plans_private.json` (exact selected IDs, quotas, assignments, warnings), `decision_well_assignments_private.csv` (960 experimental positions for three 384-well plates), and the private input checksums. The notebook also generates `cfps_decision_live_preview_PRIVATE.html` and a downloadable private ZIP. The HTML source-only viewer does not embed participants' recipes and prompts for local JSON. The SvelteKit `/cfps/decision-data` route reads `data/cfps-decision/cfps_decision_input_private.json` **from the server filesystem** only after existing authenticated-session checks; it sends `Cache-Control: private, no-store`.

No Keycloak confidential-client secret is supplied, stored, logged, embedded, or required by the decision engine. Live PocketBase credentials are requested at runtime in Colab, using Colab Secrets and `getpass` if a token is not available. If the deployment requires a Keycloak confidential-client credential, the existing server-side auth layer must obtain it separately from its secrets manager or a secure interactive provisioner; **do not add it to the viewer, notebook files, or repository**.

## Tests and caveats

- Python CLI, source adapters and standalone notebook are built from mirrored source; JS mirrors the same deterministic GP/model-free greedy policy.
- The sample plate metadata verifies reserved control columns and fourfold replication but **does not approve experimental running conditions, reagent stocks, fixed control compositions, robot ordering, or fulfillment quotas for this particular competition**.
- GP predictions are hypotheses. Geometric representation is a recipe-space proxy, not a biological performance guarantee.
- Neither JSON exports nor visualization should be loaded in anonymous static routes.
