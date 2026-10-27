> **Current release (2026-10-27):** The canonical data-source contract is in [`INPUT_DATA_AND_REPRODUCIBILITY.md`](INPUT_DATA_AND_REPRODUCIBILITY.md). The notebook supports a frozen GP-informed export or read-only PocketBase input. A PocketBase-only input does **not** reproduce the GP-informed baseline or enable hybrid mode. For contributor attribution, see [`CONTRIBUTIONS.md`](CONTRIBUTIONS.md). Older notes below are historical methodology context and must not override this source contract.

# HTGAA CFPS — pre-run selection and reproducibility contract

Date: 2026-10-27. Version 1.1.0. Source: the original **pre-run GP-informed**
`clusters_designs.json` and display asset `clusters.json.gz` (uploaded by the
researcher). No class Round 1 measurements appear in these files.

## A. Clear comparison of mathematical approaches

### GP-informed modeling: GP-based predicted response + cost/risks

A standard GP for protein `p` can be represented as

\[
 f_p(x)\sim\mathcal{GP}(m_p(x),k_p(x,x')),
 \quad y_i=f_p(x_i)+\epsilon_i.
\]

For target ingredients represented in Ginkgo's control plates, the fitted model
predicts \(\mu_p(x)\) and \(\sigma_p(x)\). The model has separate fits for the
96- and 384-well sources, according to the published GP-informed forum explanation. The GP's
predicted **output** is not a measurement from a student's formulation.

The priced recipe is \(C(x)=\sum_j p_j v_j\), with the caveat that "exact"
means *exact under the current supplied price table*, not a measured production cost.
Predicted specific cost involves \(C(x)/\hat Y_p(x)\), not an independently
observed yield. A reported `lnSpecCost` can be interpreted as a predicted log
specific cost on the published GP-informed scale; the complete GP/price implementation is needed
to reconstruct the precise transformation or expectation under uncertainty.

For reagents not varied in Ginkgo's controls (e.g., Reteplase DsbC), **the GP
cannot learn a dose-response effect from those controls**. the fitted GP-informed framework separately
assumes risk/sensitivity by ingredient type. A large cost sensitivity of a
`priceOnly=True` ingredient is primarily a **cost derivative**, *not a learned
protein yield derivative*.

The published GP-informed clustering links model-expected co-movement of prediction errors,
which may arise from shared extrapolation/untested changes as well as recipe
similarity. These links are **not** empirical Pearson correlations among protein
outputs, and are not equivalent to the submission-vector correlation module in
this distribution. Group check runs need not be one-factor experiments.

### LHS: virtual model interrogation, not recipe generation for lab execution

For each feature \(j\) and `n` virtual samples, LHS stratifies a factor's
specified range into `n` equiprobable intervals, samples one level per interval,
and independently permutes the assignments among dimensions. Inside each group,
the GP-informed framework evaluates >=160 virtual recipes **through the original GP**, then regresses
the model's cost-per-gram forecasts on ingredient values. Conceptually, with
standardized virtual inputs \(z^{(v)}\) and response \(g(x^{(v)})\), a simple
surrogate fits

\[
 \hat\beta=\arg\min_\beta\sum_v\bigl(g(x^{(v)})-\beta_0-\sum_j\beta_j z_j^{(v)}\bigr)^2.
\]

The resulting coefficient magnitude is a *local sensitivity of the model*,
not measured biological importance. In particular, the surrogate regression
\(R^2\) values close to 1 reflect the internal synthetic response and can
be dominated by a highly priced additive. the GP-informed framework then picks **actual student
submissions** as checks varying influential ingredients, using predicted odds.
No published source currently establishes a "virtual point -> nearest real
submission" matching step. Do not implement or claim one without evidence.

### Independent model-free space/contrast analysis

Let \(X=\{x_1,\dots,x_N\}\) be existing, submitted reagent vectors, with
water excluded, fixed reagents excluded, and all unique dispensing formulations
given **equal weight**. Until stock concentrations are verified, the coordinates
are *reagent-specific dispensing nL*, not molar concentrations. The default normalization divides each ingredient by its 95th percentile of **positive** observed nL and clips values above 1; an alternative uses each dimension's maximum; append optional presence indicators for sparse mixtures.

- **Recipe-pattern correlation** (optional explanatory comparator):
  \(r_{ij}=\mathrm{corr}(x_i,x_j)\) across reagent features;
  hierarchical clustering can use \(1-r_{ij}\), excluding Ward linkage.
  Sparse zeros and compositional constraints limit interpretation.
- **Greedy maximin**, given selected set \(S\): choose
  \(x^*=\arg\max_{x\in X\setminus S}\min_{s\in S}d(x,s)\).
  Promotes separation but may prefer extremes/outliers.
- **Facility-location**, given \(S\): choose the existing candidate \(x\)
  minimizing \(\sum_{u\in X}\min(d(u,x),\min_{s\in S}d(u,s))\).
  Improves mean representation of the submitted candidate population, not
  necessarily the *feasible unobserved chemical design space*.
- **Exact contrast** between designs `i`,`j`: exactly one non-water reagent
  changes; compute \(\sum_k\mathbf1\{x_{ik}\neq x_{jk}\}=1\).
- **Exact one-factor sweep:** ≥3 submitted concentrations for one reagent
  on an identical background for all other non-water reagents.
- **Approximate controls or biochemical interpretability:** **not** assessed
  by the code; require additional design constraints and scientific review.

### Three strategy modes

| Mode | Uses predicted yields / the GP-informed model role assignments | What it actually does |
|---|---|---|
| `gp_model_baseline` | Uses model-derived **roles**, no re-fit | Deduplicate and reproduce the published pre-run baseline, without filling |
| `hybrid` | Uses roles to lock baseline; geometry to fill | Keep baseline unique recipes, then fill the quota using distance optimization |
| `model_free` | Does **not** use GP predictions/roles in selection | Select actual recipes via distance objective alone (provisional benchmark) |

A full "model-based independent reoptimization" **is not implemented** because
we have neither the verified training data nor complete GP/selection code.
The `gp_model_baseline` mode is the *published proposal*, not a reconstructed model fit.
Neither the `hybrid` nor the `model_free` preview is submission-ready;
matching controls, sweeps, student attribution, final quota, batch feasibility,
and alternative scaling must be incorporated before a final decision.

## B. Published link matrix audit

The **displayed rounded** `tips[target].C` matrices were found to be symmetric,
but their minimum eigenvalues are negative:

- sfGFP: approximately -0.05316
- PETase: approximately -0.05579
- Reteplase: approximately -0.18777

Thus these *stored matrices* are not positive-semidefinite at numerical tolerance
1e-6 and must **not** be inserted directly into GP information-gain formulae.
It remains possible that underlying *unrounded* model covariance is PSD, but
that was not verified and should not be inferred from this output.

For the published GP-informed 217 group checks, the published `newInformation` percentage agrees
to integer rounding with \(100(1-r^2)\), where `r` is the check's displayed
link to the first run (maximum discrepancy about 0.50 percentage points).
This is a **residual-variance heuristic**, not a derivation of expected mutual
information or experiment utility. The audit logs the exact deviation.

## C. Frozen-input regression record

Snapshot timestamp: `2026-10-25T22:38:54+00:00`.

| Source | SHA-256 |
|---|---|
| `clusters_designs.json` | `0c312d2a4c2112e5319d988017d69c10a19bf8c5688b01ef190707ffb6ca48f6` |
| `clusters.json.gz` | `ecd2a59286e64ff3c6bacc40195be3274c16b9fd820349236688dc59e8ccc065` |

Reference results for precisely these hashes:

| Protein | records | unique dispensing formulations | unique proposed | exact single-factor group checks | unselected exact check counterparts |
|---|---:|---:|---:|---:|---:|
| sfGFP | 549 | 484 | 318 | 18 | 5 |
| PETase | 500 | 459 | 284 | 25 | 10 |
| Reteplase | 551 | 488 | 179 | 18 | 3 |

`320` is a **working** quota, not a certified competition rule.
The snapshot is immutable: fresh PocketBase exports need a new SHA and a new
report. The notebook tests the pinned regression figures *only when input SHA
matches*. The general unit suite uses synthetic data and never contacts Ginkgo
or PocketBase.

## D. Build, test, and source parity

From the distribution root:

```bash
python -m pip install -r requirements-analysis.txt
python -m unittest discover -s tests -v
python tools/build_full_model_free_colab.py
python tools/audit_gp_model_baseline.py --designs /private/clusters_designs.json \
  --links /private/clusters.json.gz --out /private/run001 --mode audit
```

The notebook contains the importable modules' **verbatim source in executable
cells**. `tools/build_full_model_free_colab.py` regenerates the canonical notebook and stores
SHA-256 hashes of both sources under notebook metadata `analysis_provenance`.
At run time, it **does not import** those mirror modules, and the GP-informed model-only mode
does not need a PocketBase login or Git checkout. Python package/dependency
versions and hashes are recorded in the private output `manifest.json`.

## E. Limits to a publication-standard analysis

**Reproduced:** input checksums, exact recipe equivalence, volume totals,
cluster check comparisons, sweep completeness, model-link *display* diagnostics,
unique baselines, deterministic geometrical selection previews and test suite.

**Not yet reproduced independently:** the published GP-informed GP training and validation,
LHS random seeds and virtual data, construction of uncertainty links, true
measurement-based validation, reference-price schedule, official quota/eligibility
rules, reagent stock-molarity registry, controls/batch constraints, scientific
expected information gain, and the final approved selection.

These should be identified as *missing dependencies*, not silently replaced
by assumptions. Do not include unapproved row-level competition exports in Git.
