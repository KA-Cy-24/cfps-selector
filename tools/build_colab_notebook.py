#!/usr/bin/env python3
"""Build a self-contained Colab notebook from the supplied original and module sources.

This is a *build-time* tool only. The generated notebook does not import the mirrored
modules at runtime. The module source code is embedded verbatim as executable cells.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import textwrap
import nbformat as nbf

ROOT = Path(__file__).resolve().parent.parent
BASE = ROOT / "HTGAA_CFPS_Fetcher_Colab.ipynb"
DEST = ROOT / "notebooks" / "cfps_database_snapshot_and_exploration.ipynb"
CORR = ROOT / "analysis/cfps_correlation/pairwise_correlation_clustering.py"
SELECT = ROOT / "analysis/cfps_selection/selection_framework.py"


def md(s): return nbf.v4.new_markdown_cell(textwrap.dedent(s).strip())
def py(s): return nbf.v4.new_code_cell(textwrap.dedent(s).strip())

def guarded(cell, flag):
    return py(f"# Optional section (enable {flag} in configuration cell).\nif {flag}:\n" + textwrap.indent(cell.source, "    "))


def main():
    original = nbf.read(BASE, as_version=4)
    cells = []
    first = original.cells[0].source
    first = first.replace("# HTGAA CFPS Coopetition — database snapshot and exploratory visualizations", "# HTGAA CFPS — reproducible clustering and pre-run selection audit")
    cells.append(md(first + "\n\n**Two independent data paths:** optional read-only PocketBase snapshot/recipe-pattern correlation (Sections 1–12); and the published GP-informed uploaded `clusters_designs.json` pre-run audit (Sections 13–16). The second path requires no credentials, live database, or installed local `.py` modules. All row-level results are private by default. The compiled notebook embeds the exact module source; it does not import or create the mirrored `.py` modules during its analysis."))
    cells.append(py('''
    # Run configuration. Default: audit the published GP-informed already exported design data, without PocketBase access.
    RUN_POCKETBASE_FETCH = False      # opt-in: credentials and live database snapshot
    RUN_CORRELATION_CLUSTERING = False # requires RUN_POCKETBASE_FETCH=True
    RUN_GP_BASELINE_AUDIT = True
    RUN_PORTFOLIO_PREVIEWS = False    # opt-in: PROVISIONAL selections, never automatic submissions
    GP_DESIGNS_PATH = "/content/clusters_designs.json"
    GP_LINKS_PATH = "/content/clusters.json.gz"  # optional; skip if not uploaded
    WORKING_QUOTA_NOTEBOOK = 320     # working assumption; verify official rules before use
    GEOMETRY_OBJECTIVE = "facility" # facility | maximin
    print("PocketBase:", RUN_POCKETBASE_FETCH, "GP-informed audit:", RUN_GP_BASELINE_AUDIT)
    '''))
    cells.append(py('''
    # Environment (Colab). GP-export-only analysis does not clone GitHub or request credentials.
    from pathlib import Path
    import os, sys, subprocess
    REPO_URL = "https://github.com/RonanChance/htgaa-automation.git"
    REPO_ROOT = Path("/content/htgaa-automation")
    PRIVATE_ROOT = Path("/content/cfps_private_snapshots")
    PRIVATE_ROOT.mkdir(parents=True, exist_ok=True)
    if RUN_POCKETBASE_FETCH:
        if not (REPO_ROOT / ".git").exists():
            if REPO_ROOT.exists() and any(REPO_ROOT.iterdir()):
                raise RuntimeError(f"{REPO_ROOT} exists but is not a Git checkout")
            subprocess.run(["git", "clone", "--depth", "1", REPO_URL, str(REPO_ROOT)], check=True)
        REPO_SHA = subprocess.check_output(["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"], text=True).strip()
        REPORT_DIR = PRIVATE_ROOT / "aggregate_review"
        REPORT_DIR.mkdir(parents=True, exist_ok=True)
    else:
        REPO_SHA = None
        REPORT_DIR = PRIVATE_ROOT / "aggregate_review"
        REPORT_DIR.mkdir(parents=True, exist_ok=True)
    subprocess.run([sys.executable, "-m", "pip", "-q", "install", "numpy", "pandas", "scipy", "requests", "matplotlib"], check=True)
    print("Read-only repository reference:", REPO_URL, "(cloned only in PocketBase mode)")
    print("Private analysis folder:", PRIVATE_ROOT)
    '''))
    # Existing snapshot, exploratory figures, and QC. Wrapped so GP-export-only run is independent of credentials.
    for c in original.cells[2:13]:
        if c.cell_type=="code":
            cells.append(guarded(c,"RUN_POCKETBASE_FETCH"))
        else:
            cells.append(c)
    cells.append(md(r'''
    ## 12. Recipe-vector pairwise correlation (independent of the published GP-informed model links)

    Define each submitted recipe as a reagent-feature vector \(x_i\). The Spearman/Pearson
    similarity \(\mathrm{corr}(x_i,x_j)\) correlates their patterns **across reagent features**, not
    their observed yields and not the covariance of the published GP-informed Gaussian-process predictions.
    A dissimilarity \(d_{ij}=1-r_{ij}\) with average-linkage clustering is only a
    descriptive, model-free alternative. Features with zero variation within each protein are
    dropped. Spearman on sparse recipes can be dominated by tied zeros and is a sensitivity
    comparison—not the recommended chemical-space optimization distance.

    The entire correlation class source is defined in the next cell, with no file writing.
    '''))
    corr_source = CORR.read_text(encoding="utf-8")
    cells.append(nbf.v4.new_code_cell(corr_source))
    # Old module fitting, output and sensitivity adapted to optional snapshot, without import from analysis.
    for i in [19,20,22]:
        s = original.cells[i].source
        s = s.replace("from analysis.cfps_correlation import CFPSSubmissionCorrelationAdapter", "# CFPSSubmissionCorrelationAdapter defined above (self-contained).")
        cells.append(guarded(py(s),"RUN_CORRELATION_CLUSTERING and RUN_POCKETBASE_FETCH"))
    cells.append(md('''
    ## 13. GP-informed design audit: exact submitted recipes, no GP retraining

    **Source:** the published GP-informed published `clusters_designs.json`. The notebook reads exactly this
    file from `GP_DESIGNS_PATH` (upload it using Colab's left-side Files pane). The
    compressed `clusters.json.gz` contains presentation data and rounded link matrices;
    upload it optionally to obtain a separate link-matrix audit.

    **Audit subroutine:** `audit_gp_model_snapshot(data, working_quota=...)`.

    - Preserve stable original IDs, targets, roles, and pre-run forecasts.
    - Verify nominal 20,000 nL volumes; identify duplicate dispensing formulations
      ignoring compensating water; count *unique* first/check/independent/sweep runs.
    - For every cluster check, compare it against its group's first run and count all
      changed reagents (not merely the labeled change).
    - Find existing exact matched counterparts and audit ≥3-level single-factor sweeps.
    - Emit a manifest with SHA-256 of input, package versions, parameters, and private CSVs.

    **Caveat:** dispensing volumes are not automatically concentrations. No inferred stock
    molarities. These are pre-experimental selection diagnostics, not measured effects.
    '''))
    selection_source = SELECT.read_text(encoding="utf-8")
    cells.append(nbf.v4.new_code_cell(selection_source))
    cells.append(py('''
    # 14. Run pre-experimental audit from a local export without requiring cloned modules.
    from pathlib import Path
    from datetime import datetime, timezone
    import json
    import pandas as pd
    if RUN_GP_BASELINE_AUDIT:
        gp_snapshot_path = Path(GP_DESIGNS_PATH).expanduser()
        if not gp_snapshot_path.exists():
            raise FileNotFoundError(
                f"Upload clusters_designs.json into Colab's Files panel and set GP_DESIGNS_PATH. "
                f"Current path: {gp_snapshot_path}"
            )
        gp_data = read_gp_model_file(gp_snapshot_path)
        audit = audit_gp_model_snapshot(gp_data, working_quota=WORKING_QUOTA_NOTEBOOK)
        audit_out = PRIVATE_ROOT / ("gp_model_audit_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
        gp_links_path = Path(GP_LINKS_PATH).expanduser()
        audit_manifest = audit.save(audit_out, input_paths=[gp_snapshot_path]+([gp_links_path] if gp_links_path.exists() else []))
        if gp_links_path.exists():
            link_summary = audit_gp_model_link_matrices(gp_links_path, design_audit=audit)
            (audit_out/"model_link_matrix_audit.json").write_text(json.dumps(link_summary,indent=2),encoding="utf-8")
            print("Published link-matrix audit:", {k:v["min_eigenvalue_of_symmetrized_DISPLAY_matrix"] for k,v in link_summary["targets"].items()})
        audit_summary_df = pd.DataFrame([
            {"protein":name, "candidate_count":x["candidateDesigns"],
             "distinct_formulations":x["distinctFormulations"],
             "gp_baseline_unique_runs":x["selectedUniqueFormulations"],
             "quota_gap":x["gapToWorkingQuota"],
             "exact_one_factor_checks":x["exactOneFactorChecksVsFirst"],
             "additional_exact_matches":x["newExactCounterpartsForChecks"]}
            for name,x in sorted(audit.summary["proteins"].items())
        ])
        display(audit_summary_df)
        print("SHA256", audit_manifest["input_sha256"])
        print("PRIVATE output:", audit_out)
        if "snapshot" in globals() and isinstance(snapshot,dict):
            live_ids = {str(r.get("id")) for r in snapshot["collections"].get("cfps_designs",[])}
            gp_ids = {rid for block in audit.records.values() for rid in block}
            crosswalk = {"shared_design_ids":len(live_ids & gp_ids),
                         "live_only_ids":len(live_ids-gp_ids),
                         "gp_export_only_ids":len(gp_ids-live_ids)}
            print("Snapshot crosswalk (potentially different snapshot times):",crosswalk)
    '''))
    cells.append(py('''
    # 15. Diagnostic selection modes on ACTUAL recipe IDs (opt-in; not final selection).
    if RUN_GP_BASELINE_AUDIT and RUN_PORTFOLIO_PREVIEWS:
        from copy import deepcopy
        preview = {}
        for target in sorted(audit.records):
            rows = [audit.records[target][rid] for rid in sorted(audit.records[target])]
            preview[target] = {}
            for mode in ("gp_model_baseline", "hybrid", "model_free"):
                alternative = suggest_portfolio(
                    rows, quota=WORKING_QUOTA_NOTEBOOK, mode=mode, objective=GEOMETRY_OBJECTIVE,
                )
                diagnostics = coverage_metrics(rows, set(alternative["ids"]))
                preview[target][mode] = alternative
                print(target,mode,"n:", alternative["count"],
                      "mean coverage distance:", round(diagnostics["mean_candidate_to_selected_distance"],5))
        write_selection_overlay(audit, preview, audit_out/"selection_overlay_private_summary.json")
        (audit_out/"selection_ids_private.json").write_text(json.dumps(
            {target:{mode:item["ids"] for mode,item in modes.items()}
             for target,modes in preview.items()}, indent=2), encoding="utf-8")
        print("Provisional candidate IDs saved privately; do not submit without controls/contrast review.")
    '''))
    cells.append(md(r'''
    ## 16. Exact mathematical interpretation and decision framework

    **the published GP-informed model-based approach.** A protein-specific Gaussian-process (GP)
    regression is fitted on Ginkgo's measured control wells, schematically
    \(f_p(x)\sim\mathcal{GP}(m_p(x),k_p(x,x'))\). For each submitted recipe \(x_i\),
    the GP yields predicted output and uncertainty for the *plate-tested* features.
    The price \(C(x_i)\) comes from the reagent price schedule. His model projects
    cost per gram \(C(x_i)/\hat Y_p(x_i)\), with additional **assumed**, not
    measured, risks for untested ingredients. His expected co-movement links combine
    recipe similarity and shared untested risks; those links are not empirical
    correlations or guaranteed covariance matrices.

    **What LHS actually does.** Within a cluster, ≥160 **virtual** formulations
    sample factor ranges with Latin hypercube stratification; each is passed to the
    *existing model*, and a second regression approximates the model's predictions
    over that region, identifying influential factors. GP-informed then picks **actual
    submitted recipes** as checks that differ on those factors, considering predicted
    odds. This is *not* established to be a nearest-neighbor projection of virtual
    recipes to submitted ones, nor a claim that LHS generates lab-ready candidates.
    It tests the model's sensitivities rather than establishing causal effects.

    **Your model-free complement.** For existing recipes \(x_i\), omit water and
    constant reagent dimensions, scale each variable stock-volume dimension by its
    maximum observed value (to be replaced with chemically meaningful molar limits
    when verified), and optionally append an ingredient-presence indicator. Then:

    - Greedy **maximin** fill: \(x^*=\arg\max_{x\in X\\S}\min_{s\in S} d(x,s)\).
    - Greedy **facility-location** fill: choose \(x^*\) minimizing
      \(\sum_{u\in X}\min(d(u,x),\min_{s\in S}d(u,s))\).
    - *Exact matched contrast* exists when two recipes vary in exactly one
      non-water ingredient; a clean sweep has a common non-swept background
      and at least three distinct levels.

    Three configurable strategies are included: `gp_model_baseline` (model-driven baseline),
    `hybrid` (baseline locked, geometry fills quota), `model_free` (no GP roles or
    forecasts used). **These strategies are provisional computational baselines.**
    The final portfolio requires independent review of controls, design validity,
    sweep coverage, author representation, quota rules, and biochemical interactions.

    Critically, facility location covers the **observed candidate distribution**, not
    all feasible CFPS chemistry; maximin can emphasize outliers. A distance-based
    output is not a formal expected information-gain calculation.

    ### Reproducibility requirements

    Keep frozen `clusters_designs.json`, input hash, code snapshot, raw Ginkgo model
    dependencies if attempting to reproduce the published GP-informed GP, deterministic tie-breaking,
    audit version, and dependency versions. Raw ID-level outputs remain outside Git.
    The current notebook does **not** reproduce the published GP-informed fitted GP, its error
    covariance, virtual LHS seed, or model-risk hyperparameters from first principles;
    those would require his original training data and full source.
    '''))
    cells.append(py('''
    # 16.1 Snapshot-specific regression checks. Adapt explicitly if the source data change.
    if RUN_GP_BASELINE_AUDIT:
        EXPECTED_FILE_SHA256 = "0c312d2a4c2112e5319d988017d69c10a19bf8c5688b01ef190707ffb6ca48f6"
        if input_sha256(gp_snapshot_path) == EXPECTED_FILE_SHA256:
            expected = {
                "sfGFP": (549,484,318,18,5),
                "PETase": (500,459,284,25,10),
                "Reteplase": (551,488,179,18,3),
            }
            for protein,(n,distinct,chosen,one_factor,new_pairs) in expected.items():
                p = audit.summary["proteins"][protein]
                observed=(p["candidateDesigns"],p["distinctFormulations"],p["selectedUniqueFormulations"],
                          p["exactOneFactorChecksVsFirst"],p["newExactCounterpartsForChecks"])
                assert observed==(n,distinct,chosen,one_factor,new_pairs), (protein,observed)
            print("PASS: known-input regression checks")
        else:
            print("Source differs from the pinned data SHA256: review revised dataset and version outputs.")
    '''))
    cells.append(md('''
    ## 17. Deliverables, permissions and sources

    This notebook is **standalone**: the correlation and selection functions are
    *defined in executable cells*, so the notebook never imports mirrored project `.py`
    modules. The repository distribution includes importable modules and synthetic unit
    tests for headless analyses. The reproducibility manifest records the data hash,
    runtime dependencies, version and normalization choices.

    - Notebook source: `notebooks/cfps_database_snapshot_and_exploration.ipynb`
    - Modules: `analysis/cfps_correlation/pairwise_correlation_clustering.py`,
      `analysis/cfps_selection/selection_framework.py`
    - CLI: `tools/audit_gp_model_baseline.py`
    - Tests: `python -m unittest discover -s tests -v`

    No Git operations, private data publishing, or competition submissions occur automatically.
    The `selection_overlay_private_summary.json` is a preliminary scaffold for later
    integration into the published GP-informed website; his `index.md`, loader, and stylesheet should
    only be edited after reviewing the actual visualization source and permissions.
    '''))
    result = nbf.v4.new_notebook(cells=cells, metadata=original.metadata.copy())
    result.metadata.setdefault("colab",{})["name"] = DEST.name
    result.metadata["analysis_provenance"] = {
        "generation": "tools/build_colab_notebook.py",
        "correlation_source_sha256": hashlib.sha256(corr_source.encode()).hexdigest(),
        "selection_source_sha256": hashlib.sha256(selection_source.encode()).hexdigest(),
        "embedded_modules": True,
        "requires_external_mirror_modules": False,
    }
    # Stable IDs: original bundle uses legacy cells without IDs; normalize for current nbformat.
    for index, cell in enumerate(result.cells):
        if not cell.get("id"):
            cell["id"] = hashlib.sha1((str(index)+"|"+cell.source).encode()).hexdigest()[:12]
    nbf.validate(result)
    DEST.parent.mkdir(parents=True,exist_ok=True)
    nbf.write(result,DEST)
    print(DEST, len(cells),"cells")
    print(json.dumps(result.metadata["analysis_provenance"],indent=2))


if __name__ == "__main__":
    main()
