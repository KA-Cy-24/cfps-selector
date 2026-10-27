"""Reproducible, pre-experimental CFPS candidate audit and portfolio diagnostics.

Inputs: the dated GP-informed clusters_designs.json export (no measurements required).
All output is a description of submitted recipes or PRE-RUN model predictions.
No code here retrains the published GP or estimates empirical response correlation.

The same source is embedded *verbatim* into the independently executable Colab notebook.
Version: 1.1.0 (2026-10-27).
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import gzip
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

VERSION = "1.1.0"
WORKING_QUOTA = 320  # NOT verified final competition rules
REACTION_NL = 20_000.0
EXCLUDED_BALANCERS = frozenset({"nuclease_free_water"})
BASELINE_ROLES = frozenset({"first run", "check run", "independent run", "sweep run"})
TARGETS = ("sfGFP", "PETase", "Reteplase")


def canonical_recipe(recipe: Mapping[str, Any], *, exclude_water: bool = True) -> tuple:
    """Immutable recipe identity; zeros and absent keys are equivalent.

    This identifies exact *dispensing recipes*, NOT chemical equivalence when stock
    concentrations or identities differ. Does not compare across protein targets.
    """
    if not isinstance(recipe, Mapping):
        raise TypeError("recipe must be a mapping")
    out = []
    for key, raw in recipe.items():
        if exclude_water and key in EXCLUDED_BALANCERS:
            continue
        if isinstance(raw, bool):
            raise ValueError(f"Boolean volume for {key}")
        try:
            vol = float(raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Non-numeric volume for {key}: {raw!r}") from exc
        if not np.isfinite(vol) or vol < 0:
            raise ValueError(f"Invalid volume for {key}: {vol!r}")
        if vol != 0:
            out.append((str(key), round(vol, 9)))
    return tuple(sorted(out))


def signature_except(recipe: Mapping[str, Any], factor: str) -> tuple:
    return tuple((k, v) for k, v in canonical_recipe(recipe) if k != factor)


def recipe_factor_value(recipe: Mapping[str, Any], factor: str) -> float:
    return float(recipe.get(factor, 0))


def changed_reagents(recipe_a: Mapping[str, Any], recipe_b: Mapping[str, Any]) -> list[str]:
    a, b = dict(canonical_recipe(recipe_a)), dict(canonical_recipe(recipe_b))
    return [key for key in sorted(a.keys() | b.keys()) if not np.isclose(a.get(key, 0), b.get(key, 0), rtol=0, atol=1e-8)]


def input_sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _factor_key_from_display(ingredient: str, record: Mapping[str, Any]) -> str:
    """Resolve check's display label to stored reagent key using check's own data.

    Prefer exact `group.checks.ingredient` mapped through displayed names from
    vsGinkgoPlates or sweeps; handle common aliases for non-tested reagents.
    Fail closed when unknown, instead of guessing from textual similarity.
    """
    common = {
        "DsbC (E. coli, 25.5 kDa)": "dsbc_ecoli",
        "PEG-8000": "peg_8000",
        "Magnesium Glutamate": "magnesium_glutamate",
        "HEPES-KOH pH 7.5": "hepes_koh",
        "Coenzyme A": "coa",
        "PEP mono (0.1M stock)": "pep_mono",
        "17 Amino Acid Mix": "amino_acid_mix_17",
        "Nicotinamide": "nicotinamide",
        "Glucose": "glucose",
        "Folinic Acid": "folinic_acid",
        "Potassium Glutamate": "potassium_glutamate",
        "Potassium Oxalate": "potassium_oxalate",
        "Ammonium Glutamate": "ammonium_glutamate",
        "Reduced Glutathione (GSH)": "reduced_glutathione",
        "Oxidized Glutathione (GSSG)": "oxidized_glutathione",
        "Maltodextrin 17": "maltodextrin_17",
        "Spermidine": "spermidine",
        "Tyrosine pH 12": "tyrosine",
        "Potassium phosphate 1.6:1 dibasic:monobasic (molar ratio)": "kpo_dibasic_mix",
        "Potassium phosphate 1.6:1 monobasic:dibasic (molar ratio)": "kpo_monobasic_mix",
        "Putrescine Dihydrochloride": "putrescine",
    }
    for category in ("neverOnAPlate", "pastThePlatesRange"):
        for d in record.get("vsGinkgoPlates", {}).get(category, []):
            if d.get("name") == ingredient:
                return str(d["reagent"])
    if ingredient in common:
        return common[ingredient]
    slug = ingredient.lower().replace("-", "_").replace(" ", "_")
    if slug in record.get("recipeNl", {}):
        return slug
    for reagent in record.get("recipeNl", {}):
        if str(reagent).casefold() == ingredient.casefold():
            return reagent
    return ""  # Caller must flag unresolved; no inferential fallback.


def _validated_records(data: Mapping[str, Any]) -> tuple[dict[str, dict[str, dict]], list[dict]]:
    if not isinstance(data, Mapping) or not isinstance(data.get("proteins"), Mapping):
        raise ValueError("Expected published clusters_designs.json with proteins mapping")
    by_target = {}
    rows = []
    global_seen = set()
    for target in sorted(data["proteins"]):
        block = data["proteins"][target]
        designs = block.get("designs")
        if not isinstance(designs, list):
            raise ValueError(f"Missing designs list for {target}")
        lookup = {}
        for design in designs:
            rid = str(design["id"])
            if rid in lookup or rid in global_seen:
                raise ValueError(f"Nonunique design ID: {rid}")
            global_seen.add(rid)
            if design.get("protein") != target:
                raise ValueError(f"Design {rid} target mismatch")
            recipe = design.get("recipeNl")
            key = canonical_recipe(recipe)
            total = sum(float(v) for v in recipe.values())
            lookup[rid] = design
            rows.append({
                "protein": target, "id": rid,
                "sampleId": design.get("sampleId"), "author": design.get("author"),
                "round1Role": design.get("round1Role"),
                "selectedProposed": bool(design.get("round1Role") in BASELINE_ROLES),
                "sameRecipeAs": design.get("sameRecipeAs"),
                "group": (design.get("group") or {}).get("letter"),
                "recipe_key": key,
                "recipe_hash": hashlib.sha256(json.dumps(key).encode()).hexdigest(),
                "total_nL": total,
                "total_valid": abs(total - REACTION_NL) < 1e-5,
                "n_active_reagents": len(key),
                "pBeat_model": (design.get("forecast") or {}).get("pBeat"),
                "lnCostPerGram_model": (design.get("forecast") or {}).get("lnSpecCost"),
            })
        by_target[target] = lookup
    return by_target, rows


def _sweep_rows(data: Mapping[str, Any], by_target: dict, selected_ids: dict) -> list[dict]:
    result = []
    for target in sorted(by_target):
        lookup = by_target[target]
        for ingredient_block in data["proteins"][target].get("sweeps", []):
            factor = ingredient_block["reagent"]
            for family_index, sweep in enumerate(ingredient_block.get("sweeps", [])):
                to_run = sweep.get("toRun") or []
                if not to_run:
                    continue
                ids = [x.get("id") if isinstance(x, Mapping) else x for x in to_run]
                ids = [str(x) for x in ids]
                missing = [rid for rid in ids if rid not in lookup]
                if missing:
                    raise ValueError(f"Unknown IDs in sweep {target}/{factor}: {missing}")
                unique_selected = {}
                background = set()
                for rid in ids:
                    if rid not in selected_ids[target]:
                        continue
                    design = lookup[rid]
                    unique_selected[canonical_recipe(design["recipeNl"])] = design
                    background.add(signature_except(design["recipeNl"], factor))
                levels = sorted({float(r["recipeNl"].get(factor, 0)) for r in unique_selected.values()})
                result.append({
                    "protein": target, "factor": factor, "family_index": family_index,
                    "recommendedByModel": bool(ingredient_block.get("recommended")),
                    "n_toRun_ids": len(ids), "n_toRun_unique": len(unique_selected),
                    "n_selected_levels": len(levels), "selected_levels_nL": ";".join(map(lambda x: str(int(x)) if x.is_integer() else str(x), levels)),
                    "background_signature_count": len(background),
                    "clean_single_factor_sweep": len(background) == 1,
                    "at_least_three_levels": len(levels) >= 3,
                    "run_ids": ";".join(ids),
                })
    return result


def _cluster_checks(data: Mapping[str, Any], by_target: dict, selected_ids: dict) -> tuple[list[dict], list[dict]]:
    checks, additions = [], []
    for target in sorted(by_target):
        lookup = by_target[target]
        for group in data["proteins"][target].get("groups", []):
            base_id = str(group["firstRun"])
            if base_id not in lookup:
                raise ValueError(f"Unknown first run: {base_id}")
            base = lookup[base_id]
            for check in group.get("checkRuns", []):
                check_id = str(check["id"])
                if check_id not in lookup:
                    raise ValueError(f"Unknown check run: {check_id}")
                design = lookup[check_id]
                factor = _factor_key_from_display(check["ingredient"], design)
                if not factor:
                    # Try first-run fields, as named ingredient might be absent in check
                    factor = _factor_key_from_display(check["ingredient"], base)
                changed = changed_reagents(base["recipeNl"], design["recipeNl"])
                named_changed = factor in changed if factor else False
                exact_selected = []
                exact_unselected = []
                if factor:
                    sig = signature_except(design["recipeNl"], factor)
                    for other_id, other in lookup.items():
                        if other_id == check_id:
                            continue
                        if signature_except(other["recipeNl"], factor) != sig:
                            continue
                        if np.isclose(recipe_factor_value(other["recipeNl"], factor), recipe_factor_value(design["recipeNl"], factor)):
                            continue
                        (exact_selected if other_id in selected_ids[target] else exact_unselected).append(other_id)
                exact_selected = sorted(exact_selected)
                exact_unselected = sorted(exact_unselected)
                checks.append({
                    "protein": target, "group": group["letter"],
                    "firstRunId": base_id, "checkRunId": check_id,
                    "namedFactorLabel": check["ingredient"], "namedFactor": factor,
                    "named_factor_changed": named_changed,
                    "changed_factor_count": len(changed),
                    "changed_factors": ";".join(changed),
                    "exact_one_factor_with_first_run": len(changed) == 1 and named_changed,
                    "any_selected_exact_counterpart": bool(exact_selected),
                    "any_unselected_exact_counterpart": bool(exact_unselected),
                    "model_link_with_first": (design.get("group") or {}).get("linkWithFirstRun"),
                    "model_new_information_pct": check.get("newInformation"),
                    "model_price_only": check.get("priceOnly"),
                })
                if not exact_selected and exact_unselected and factor:
                    # Stable deterministic candidate: prefer the largest dose contrast, then lowest ID.
                    own_level = recipe_factor_value(design["recipeNl"], factor)
                    chosen = sorted(exact_unselected, key=lambda rid: (-abs(recipe_factor_value(lookup[rid]["recipeNl"], factor) - own_level), rid))[0]
                    additions.append({
                        "protein": target, "group": group["letter"],
                        "checkRunId": check_id, "namedFactor": factor,
                        "checkFactorNl": own_level, "candidateId": chosen,
                        "candidateFactorNl": recipe_factor_value(lookup[chosen]["recipeNl"], factor),
                        "is_exact_one_factor_against_check": True,
                    })
    return checks, additions


@dataclass
class AuditResult:
    summary: dict
    candidates: pd.DataFrame
    duplicates: pd.DataFrame
    check_runs: pd.DataFrame
    sweep_families: pd.DataFrame
    additional_contrasts: pd.DataFrame
    records: dict[str, dict[str, dict]]

    def save(self, out: str | Path, *, input_paths: Sequence[str | Path] = ()) -> dict:
        """Write private row-level CSVs and provenance manifest.

        The `publish_row_level` flag requires explicit opt-in; caller controls actual
        directory permissions. No code uploads or pushes files to Git automatically.
        """
        out = Path(out)
        out.mkdir(parents=True, exist_ok=True)
        exports = {
            "candidate_audit.csv": self.candidates.drop(columns=["recipe_key"], errors="ignore"),
            "duplicate_formulations.csv": self.duplicates,
            "cluster_check_audit.csv": self.check_runs,
            "sweep_family_audit.csv": self.sweep_families,
            "additional_exact_contrasts.csv": self.additional_contrasts,
        }
        for filename, frame in exports.items():
            frame.to_csv(out / filename, index=False)
        with open(out / "audit_summary.json", "w") as fh:
            json.dump(self.summary, fh, indent=2, sort_keys=True, allow_nan=False)
        manifest = {
            "schema": "cfps-selection-audit/v1", "analysis_version": VERSION,
            "input_sha256": {Path(p).name: input_sha256(p) for p in input_paths},
            "working_quota": self.summary["working_quota_per_protein"],
            "python": __import__("sys").version.split()[0],
            "numpy": np.__version__, "pandas": pd.__version__,
            "scipy": __import__("scipy").__version__,
            "parameters": {
                "reaction_nL": REACTION_NL,
                "baseline_roles": sorted(BASELINE_ROLES),
                "duplicate_recipe_excludes": sorted(EXCLUDED_BALANCERS),
                "recipe_volume_unit": "nL (stock dispensing volumes, not molarity)",
                "selected_unique_rule": "one exact dispensing formulation per target",
            },
            "files": sorted(exports) + ["audit_summary.json"],
            "classification": "row-level competition data; private until approved",
            "publication_review_required": True,
        }
        (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return manifest


def audit_gp_model_snapshot(data: Mapping[str, Any], *, working_quota: int = WORKING_QUOTA) -> AuditResult:
    """Deterministic end-to-end audit; never modifies raw data or model forecasts."""
    if int(working_quota) < 1:
        raise ValueError("working_quota must be positive")
    by_target, rows = _validated_records(data)
    candidates = pd.DataFrame(rows).sort_values(["protein", "id"], kind="stable").reset_index(drop=True)
    # Internal exact-key maps; no cross-target merging.
    selected_ids = {}
    duplicates = []
    summary_proteins = {}
    for target in sorted(by_target):
        part = candidates.loc[candidates["protein"] == target].copy()
        selected = part.loc[part["selectedProposed"]]
        all_groups = defaultdict(list)
        for r in part.itertuples(index=False):
            all_groups[r.recipe_key].append(r.id)
        selected_groups = defaultdict(list)
        for r in selected.itertuples(index=False):
            selected_groups[r.recipe_key].append(r.id)
        selected_ids[target] = set(selected["id"])
        for signature, ids in all_groups.items():
            if len(ids) > 1:
                duplicates.append({
                    "protein": target, "count": len(ids), "ids": ";".join(sorted(ids)),
                    "canonical_signature_sha256": hashlib.sha256(repr(signature).encode()).hexdigest(),
                    "selected_proposed_count": sum(rid in selected_ids[target] for rid in ids),
                })
        block = data["proteins"][target]
        untested_groups = sum(bool(g.get("linkFromSharedUntested", 0) > 0) for g in block.get("groups", []))
        summary_proteins[target] = {
            "candidateDesigns": int(len(part)),
            "distinctFormulations": int(len(all_groups)),
            "duplicateExcess": int(len(part) - len(all_groups)),
            "selectedRoleRecords": int(len(selected)),
            "selectedUniqueFormulations": int(len(selected_groups)),
            "selectedDuplicateExcess": int(len(selected) - len(selected_groups)),
            "gapToWorkingQuota": int(working_quota - len(selected_groups)),
            "roleBreakdownUniqueDesignIDs": {str(k): int(v) for k,v in selected["round1Role"].value_counts().sort_index().items()},
            "groups": len(block.get("groups", [])),
            "groupsWithSharedUntested": int(untested_groups),
            "authorsAll": int(part["author"].nunique(dropna=True)),
            "authorsSelected": int(selected["author"].nunique(dropna=True)),
            "volumeValidRecords": int(part["total_valid"].sum()),
            "sourceTotals": block.get("totals", {}),
        }
    checks, add = _cluster_checks(data, by_target, selected_ids)
    sweeps = _sweep_rows(data, by_target, selected_ids)
    check_cols = ["protein", "group", "firstRunId", "checkRunId", "namedFactorLabel", "namedFactor", "named_factor_changed", "changed_factor_count", "changed_factors", "exact_one_factor_with_first_run", "any_selected_exact_counterpart", "any_unselected_exact_counterpart", "model_link_with_first", "model_new_information_pct", "model_price_only"]
    sweep_cols = ["protein", "factor", "family_index", "recommendedByModel", "n_toRun_ids", "n_toRun_unique", "n_selected_levels", "selected_levels_nL", "background_signature_count", "clean_single_factor_sweep", "at_least_three_levels", "run_ids"]
    add_cols = ["protein", "group", "checkRunId", "namedFactor", "checkFactorNl", "candidateId", "candidateFactorNl", "is_exact_one_factor_against_check"]
    checks_df = pd.DataFrame(checks,columns=check_cols).sort_values(["protein", "group", "checkRunId"]).reset_index(drop=True)
    sweeps_df = pd.DataFrame(sweeps,columns=sweep_cols).sort_values(["protein", "factor", "family_index"]).reset_index(drop=True)
    add_df = pd.DataFrame(add,columns=add_cols).sort_values(["protein", "group", "checkRunId"]).reset_index(drop=True)
    for target in summary_proteins:
        c = checks_df[checks_df["protein"] == target]
        s = sweeps_df[sweeps_df["protein"] == target]
        a = add_df[add_df["protein"] == target]
        summary_proteins[target].update({
            "clusterChecks": int(len(c)),
            "exactOneFactorChecksVsFirst": int(c["exact_one_factor_with_first_run"].sum()),
            "multifactorChecksVsFirst": int((c["changed_factor_count"] > 1).sum()),
            "unresolvedNamedFactorChecks": int((c["namedFactor"] == "").sum()),
            "newExactCounterpartsForChecks": int(len(a)),
            "listedSweepFamiliesWithRuns": int(len(s)),
            "cleanSweepFamilies": int(s["clean_single_factor_sweep"].sum()),
            "sweepFamiliesWithThreeOrMoreLevels": int(s["at_least_three_levels"].sum()),
        })
    duplicate_df = pd.DataFrame(duplicates).sort_values(["protein", "canonical_signature_sha256"]).reset_index(drop=True)
    summary = {
        "analysis_version": VERSION,
        "source_built_utc": data.get("built"),
        "scope": "audit of published GP-informed proposed runs, not measured results or final selection",
        "working_quota_per_protein": int(working_quota),
        "proteins": summary_proteins,
        "interpretation": {
            "model_based": "forecast, expected links and LHS sensitivities are pre-run model products",
            "model_free": "recipe identities, exact contrasts, sweep completeness, distance/coverage",
            "price_per_gram": "model forecast; not validated observed cost per gram",
            "experimental_info_gain": "NOT calculated: model newInformation is a heuristic",
        },
    }
    return AuditResult(summary, candidates, duplicate_df, checks_df, sweeps_df, add_df, by_target)


def read_gp_model_file(path: str | Path) -> dict:
    with open(path, "rb") as fh:
        return json.load(fh)


def _variable_recipe_matrix(rows: list[dict], *, feature_mode: str = "scaled_volume", presence_weight: float = 0.35) -> tuple[np.ndarray, list[str]]:
    """Candidate rows -> scaled recipe matrix; no model predictions involved.

    Excludes water and columns invariant in this protein. For concentrations,
    use stock nL only as a per-reagent within-stock proxy, NOT chemical molarity.
    Normalize by max observed nL (>=1) to make each ingredient dimension comparable.
    """
    if feature_mode not in {"scaled_volume", "presence", "scaled_plus_presence"}:
        raise ValueError("Invalid feature_mode")
    all_recipes = [r["recipeNl"] for r in rows]
    keys = sorted({k for r in all_recipes for k in r if k not in EXCLUDED_BALANCERS})
    A = np.asarray([[float(rec.get(k, 0)) for k in keys] for rec in all_recipes], dtype=float)
    var = np.ptp(A, axis=0) > 1e-9
    A = A[:, var]
    keys = [k for k, keep in zip(keys, var) if keep]
    if not keys:
        raise ValueError("No variable reagent dimensions")
    scales = np.maximum(np.max(A, axis=0), 1)
    Z = A / scales
    if feature_mode == "presence":
        return (A > 0).astype(float), keys
    if feature_mode == "scaled_plus_presence":
        Z = np.column_stack([Z, (A > 0).astype(float)*math.sqrt(presence_weight)])
        return Z, keys + ["presence:"+k for k in keys]
    return Z, keys


def coverage_metrics(candidate_rows: list[dict], selected_ids: set[str], *, metric: str = "euclidean") -> dict:
    """Model-free dispersion/coverage diagnostics of an existing-recipe portfolio.

    A candidate-population metric: no claim of uniformity in unobserved chemical space.
    """
    selected_indices = [i for i,r in enumerate(candidate_rows) if r["id"] in selected_ids]
    if not selected_indices:
        raise ValueError("No selected IDs in candidate population")
    Z, keys = _variable_recipe_matrix(candidate_rows, feature_mode="scaled_plus_presence")
    ds = cdist(Z, Z[selected_indices], metric=metric)
    nearest = ds.min(axis=1)
    # Unique physical formulations, not student submission records, are the
    # scientific sampling units. Equal-weight each distinct recipe.
    signatures = [canonical_recipe(r["recipeNl"]) for r in candidate_rows]
    representative_indices = list({sig:i for i,sig in reversed(list(enumerate(signatures)))}.values())
    unique_nearest = nearest[representative_indices]
    factors = {}
    for k in sorted({k for r in candidate_rows for k in r["recipeNl"] if k not in EXCLUDED_BALANCERS}):
        all_values = np.array([float(r["recipeNl"].get(k,0)) for r in candidate_rows])
        sel_values = all_values[selected_indices]
        unique_positive = sorted(set(all_values[all_values>0]))
        # Bins drawn from candidate quantiles: a coverage *diagnostic*, not a perfect LHS.
        if len(unique_positive) >= 2:
            bounds = np.unique(np.quantile(unique_positive, np.linspace(0,1,6)))
            pos_all = np.digitize(all_values[all_values>0], bounds[1:-1])
            pos_sel = np.digitize(sel_values[sel_values>0], bounds[1:-1])
            all_bins = set(pos_all)
            sel_bins = set(pos_sel)
            factor_coverage = len(sel_bins & all_bins)/len(all_bins) if all_bins else 1.0
        else:
            factor_coverage = 1.0 if (len(unique_positive)==0 or (sel_values>0).any()) else 0.0
        factors[k] = {"candidate_bin_coverage": round(float(factor_coverage),6),
                      "has_absent": bool((sel_values==0).any()),
                      "has_present": bool((sel_values>0).any())}
    return {
        "n_selected_ids": len(selected_indices),
        "n_distinct_candidate_formulations": len(representative_indices),
        "feature_mode": "scaled_plus_presence", "n_features": Z.shape[1],
        "mean_candidate_to_selected_distance": float(unique_nearest.mean()),
        "p90_candidate_to_selected_distance": float(np.quantile(unique_nearest, 0.90)),
        "max_candidate_to_selected_distance": float(unique_nearest.max()),
        "factor_coverage": factors,
    }


def suggest_portfolio(
    candidate_rows: list[dict], *, quota: int = WORKING_QUOTA,
    mode: str = "hybrid", objective: str = "facility",
    seed_ids: Sequence[str] = (),
    exact_unique: bool = True,
) -> dict:
    """Deterministic proposed IDs. No new recipes; no GP calls.

    - 'gp_model_baseline': return the original GP-informed proposed selection *without padding*.
    - 'hybrid': lock the published GP-informed baseline plus explicitly given seeds, fill to quota.
    - 'model_free': ignore model-based roles/predictions, use only seeds and geometry.

    A *diagnostic selection* only. Does not guarantee experimental contrast quality,
    student fairness, batch feasibility, or a scientifically optimal design.
    """
    if mode not in {"gp_model_baseline", "hybrid", "model_free"}:
        raise ValueError("mode must be gp_model_baseline/hybrid/model_free")
    if objective not in {"facility", "maximin"}:
        raise ValueError("objective must be facility/maximin")
    if quota < 1:
        raise ValueError("quota must be positive")
    ordered = sorted(candidate_rows, key=lambda r:r["id"])
    idmap = {r["id"]:i for i,r in enumerate(ordered)}
    if len(idmap) != len(ordered):
        raise ValueError("Non-unique IDs")
    # Published export labels include multiple IDs for the *same* dispensing formulation.
    # Physical baseline counts are unique chemical recipes; preserve provenance via
    # duplicate-group output, and prefer a first run, check, sweep, independent.
    priority = {"first run": 0, "check run": 1, "sweep run": 2, "independent run": 3}
    reps = {}
    for r in sorted(ordered, key=lambda x:(priority.get(x["round1Role"],99),x["id"])):
        if r["round1Role"] in BASELINE_ROLES:
            reps.setdefault(canonical_recipe(r["recipeNl"]), r["id"])
    baseline = set(reps.values())
    if mode == "gp_model_baseline":
        return {"mode":mode,"objective": "baseline", "ids": sorted(baseline),
                "count":len(baseline),"quota_gap":quota-len(baseline)}
    selected = set(seed_ids)
    if unknown := selected-set(idmap):
        raise ValueError(f"Unknown seeds: {sorted(unknown)}")
    if mode == "hybrid":
        already_represented = {canonical_recipe(ordered[idmap[rid]]["recipeNl"]) for rid in selected}
        selected.update(rid for signature,rid in reps.items() if signature not in already_represented)
    if len(selected) > quota:
        raise ValueError(f"Selected seeds/baseline exceed quota {quota}")
    # Identity buckets: even when records have separate IDs, avoid paying twice for same recipe.
    keys = [canonical_recipe(r["recipeNl"]) for r in ordered]
    represented = {keys[idmap[rid]] for rid in selected}
    if exact_unique and len(represented) < len(selected):
        raise ValueError("Seed/baseline contains chemically identical formulations")
    allowed = [i for i, r in enumerate(ordered) if (not exact_unique or keys[i] not in represented)]
    Z, feature_names = _variable_recipe_matrix(ordered, feature_mode="scaled_plus_presence")
    dist = cdist(Z, Z, metric="euclidean")
    # Equal weight to distinct formulations; otherwise students submitting
    # identical recipes many times would bias coverage optimization.
    signature_count = Counter(keys)
    weights = np.asarray([1.0 / signature_count[k] for k in keys])
    weights = weights / weights.sum()
    # Initial seed-free selection: medoid minimizes total candidate distance for facility,
    # maximum total distance for maximin (farthest extremum; documented heuristic).
    best_distance = np.min(dist[:,[idmap[i] for i in selected]], axis=1) if selected else None
    while len(selected)<quota and allowed:
        if best_distance is None:
            scores = weights @ dist
            pick = min(allowed, key=lambda i:(scores[i],ordered[i]["id"])) if objective=="facility" else max(allowed,key=lambda i:(scores[i],ordered[i]["id"]))
        elif objective=="maximin":
            pick = max(allowed, key=lambda i:(float(best_distance[i]),ordered[i]["id"]))
        else:
            # Weighted candidate coverage: minimize sum of nearest distances.
            pick = min(allowed,key=lambda i:(float(weights @ np.minimum(best_distance,dist[:,i])),ordered[i]["id"]))
        chosen_id = ordered[pick]["id"]
        selected.add(chosen_id)
        represented.add(keys[pick])
        best_distance = dist[:,pick].copy() if best_distance is None else np.minimum(best_distance,dist[:,pick])
        allowed = [i for i in allowed if (i != pick and (not exact_unique or keys[i] not in represented))]
    return {"mode":mode,"objective":objective,"ids":sorted(selected),
            "count":len(selected), "quota_gap":quota-len(selected),
            "n_features":len(feature_names),
            "coverage_mean":float(weights @ best_distance) if best_distance is not None else None,
            "selection_status":"PROVISIONAL: geometry-only fill; controls/contrasts/fairness require separate review"}


def write_selection_overlay(
    audit: AuditResult, portfolios: Mapping[str, Mapping[str, dict]], path: str | Path,
    *, include_private_ids: bool = False,
) -> None:
    """Precomputed JSON for an opt-in browser visualization; author names omitted.

    ID-level overlay remains private by default; explicit approval before publishing.
    """
    result = {"schema":"cfps-selection-overlay/v1", "audit_version":VERSION,
              "snapshot_built":audit.summary["source_built_utc"], "targets":{}}
    for target in sorted(audit.records):
        port = portfolios.get(target, {})
        pitems = {}
        for label, d in port.items():
            pitems[label] = {k:v for k,v in d.items() if k not in {"ids"}}
            if include_private_ids:
                pitems[label]["ids"] = d.get("ids",[])
        result["targets"][target] = {
            "audit": {k:v for k,v in audit.summary["proteins"][target].items() if k!="sourceTotals"},
            "portfolios":pitems,
        }
    Path(path).write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")


def audit_gp_model_link_matrices(
    links: Mapping[str, Any] | str | Path, *, design_audit: AuditResult | None = None
) -> dict:
    """Audit *published/displayed* link matrices (not estimated GP covariance).

    For cross-snapshot validity, check matching row counts separately; no IDs are
    inferred from the rendered `names` strings. A rounded correlation-like matrix
    with negative eigenvalues should NOT be used as a GP covariance without an
    explicitly justified PSD repair and model verification.
    """
    if isinstance(links, (str, Path)):
        p = Path(links)
        with (gzip.open(p, "rt", encoding="utf-8") if p.suffix == ".gz" else p.open("r", encoding="utf-8")) as f:
            data = json.load(f)
    else:
        data = links
    if not isinstance(data, Mapping) or "tips" not in data:
        raise ValueError("Expected published clusters.json.gz with 'tips'")
    result = {"kind":"DISPLAY link-matrix audit, not biological outcome correlation",
              "built":data.get("built"), "targets":{}}
    names = {"sfgfp":"sfGFP", "petase":"PETase", "reteplase":"Reteplase"}
    for key, tip in sorted(data["tips"].items()):
        A = np.asarray(tip["C"],dtype=float)
        if A.ndim != 2 or A.shape[0] != A.shape[1] or not np.isfinite(A).all():
            raise ValueError(f"Invalid link matrix for {key}")
        n=A.shape[0]
        symmetric_max=float(np.max(np.abs(A-A.T)))
        B=(A+A.T)/2.0
        eig=np.linalg.eigvalsh(B)
        target=names.get(key,key)
        expected = (design_audit.summary["proteins"].get(target) or {}).get("candidateDesigns") if design_audit else None
        result["targets"][target] = {
            "n_rows":int(n),"n_names":len(tip.get("names",[])),
            "candidate_count_matches":None if expected is None else bool(n==expected),
            "entry_min":float(A.min()),"entry_max":float(A.max()),
            "diagonal_max_abs_error_from_one":float(np.max(np.abs(np.diag(A)-1))),
            "max_absolute_asymmetry":symmetric_max,
            "min_eigenvalue_of_symmetrized_DISPLAY_matrix":float(eig[0]),
            "negative_eigenvalue_count_below_minus_1e_6":int((eig < -1e-6).sum()),
            "PSD_as_stored_tolerance_1e_6": bool(eig[0]>=-1e-6),
            "interpretation":"Display values rounded; not independently verified GP covariance",
        }
    if design_audit is not None:
        check=design_audit.check_runs.dropna(subset=["model_link_with_first","model_new_information_pct"])
        if len(check):
            v=np.asarray(check["model_link_with_first"],dtype=float)
            scores=np.asarray(check["model_new_information_pct"],dtype=float)
            deviations=np.abs(scores-100*(1-v**2))
            result["newInformation_heuristic_check"] = {
                "n_checked":int(len(check)),
                "max_absolute_pct_point_difference":float(np.max(deviations)),
                "median_absolute_pct_point_difference":float(np.median(deviations)),
                "tested_formula":"100 * (1 - displayed_link_to_first**2)",
                "not_expected_information_gain":True,
            }
    return result
