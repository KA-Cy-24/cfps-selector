"""Validated participant-pool input adapters for CFPS design selection.

Accepted inputs:
- ``published_snapshot``: the dated GP-informed ``clusters_designs.json`` export;
- ``pocketbase_snapshot``: an authorized, read-only PocketBase snapshot containing
  ``collections.cfps_designs`` or the corresponding paginated records payload.

The PocketBase adapter deliberately does not infer GP forecasts, cluster membership,
run roles, experimental output, or chemical stock concentrations. Consequently,
GP-informed baseline and hybrid modes must not run on PocketBase-only input.

This module is embedded as an independent code cell in the Colab notebook.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Mapping

SOURCE_ADAPTER_VERSION = "1.1.0"
SOURCE_TYPES = ("published_snapshot", "pocketbase_snapshot")
CANONICAL_PROTEINS = ("sfGFP", "PETase", "Reteplase")


def _json_object(value: Any, label: str) -> dict:
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except (ValueError, TypeError) as exc:
            raise ValueError(f"{label} is not valid JSON") from exc
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a JSON object")
    return dict(value)


def _canonical_target(value: Any) -> str:
    text = str(value or "").strip().casefold().replace(" ", "").replace("-", "")
    lookup = {"sfgfp": "sfGFP", "petase": "PETase", "reteplase": "Reteplase"}
    if text not in lookup:
        raise ValueError(f"Unknown/missing protein target: {value!r}; no target guessed from recipe")
    return lookup[text]


def normalize_pocketbase_records(raw_snapshot: Mapping[str, Any], *,
                                 require_complete: bool = True) -> tuple[dict, dict]:
    """Create an analysis-compatible recipe pool from an authorized raw snapshot.

    The output contains only fields verifiable in PocketBase; GP-derived fields are
    *absent or null*. Invalid experimental records fail closed (no silent deletion).
    Nonexperimental records are counted and excluded from the eligible pool.
    """
    if not isinstance(raw_snapshot, Mapping):
        raise TypeError("PocketBase input must be an object")
    meta = raw_snapshot.get("metadata") or {}
    if meta and not isinstance(meta, Mapping):
        raise ValueError("Snapshot metadata must be an object")
    if require_complete and meta.get("designs_possibly_truncated") is True:
        raise ValueError("PocketBase snapshot explicitly marked truncated; use a complete fetch")
    if "collections" in raw_snapshot:
        collections = _json_object(raw_snapshot["collections"], "collections")
        records = collections.get("cfps_designs")
    elif "cfps_designs" in raw_snapshot:
        records = raw_snapshot["cfps_designs"]
    elif "records" in raw_snapshot:
        records = raw_snapshot["records"]
    elif "items" in raw_snapshot:
        records = raw_snapshot["items"]
    else:
        raise ValueError("No cfps_designs collection or records/items array")
    if not isinstance(records, list):
        raise ValueError("cfps_designs must be a list")
    if require_complete and (meta.get("designs_advertised_total_items") is not None):
        expected = int(meta["designs_advertised_total_items"])
        if len(records) < expected:
            raise ValueError(f"Incomplete PocketBase response: {len(records)}/{expected} records")
    designs = {p: [] for p in CANONICAL_PROTEINS}
    ids = set()
    excluded_types = Counter()
    missing_author_count = 0
    for index, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise ValueError(f"Record #{index} is not an object")
        payload = _json_object(record.get("design_json"), f"record[{index}].design_json")
        kind = str(payload.get("sample_type", "")).strip().lower()
        if kind != "experimental":
            excluded_types[kind or "missing"] += 1
            continue
        metadata = _json_object(payload.get("metadata") or {}, f"record[{index}].metadata")
        target = _canonical_target(metadata.get("protein_target"))
        recipe = _json_object(payload.get("reagent_list"), f"record[{index}].reagent_list")
        if not recipe:
            raise ValueError(f"Experimental record #{index} has no reagent_list")
        converted = {}
        total = 0.0
        for k, value in recipe.items():
            if isinstance(value, bool):
                raise ValueError(f"Boolean volume {k} in record #{index}")
            try:
                volume = float(value)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Non-numeric reagent volume {k} in record #{index}") from exc
            import math
            if not math.isfinite(volume) or volume < 0:
                raise ValueError(f"Invalid reagent volume {k} in record #{index}")
            total += volume
            if volume > 0:
                converted[str(k)] = volume
        if abs(total - 20000.0) > 1e-5:
            raise ValueError(f"Experimental record #{index} total volume {total} nL, expected 20000 nL")
        rid = str(record.get("id") or "").strip()
        if not rid or rid in ids:
            raise ValueError(f"Missing or duplicate immutable record ID: {rid!r}")
        ids.add(rid)
        author = (record.get("author") or record.get("author_sub") or
                  record.get("author_name") or "")
        if not author:
            missing_author_count += 1
        designs[target].append({
            "id": rid,
            "sampleId": str(payload.get("sample_id") or record.get("title") or rid),
            "protein": target,
            "created": record.get("created"),
            "author": str(author),
            "recipeNl": converted,
            # Preserve explicit participant batch metadata for honest completeness checks.
            # Never infer submitted batches from GP-defined clusters or sample ID prefixes.
            "metadata": {key: metadata[key] for key in ("batch_id", "batch_size", "batch_index") if key in metadata},
            # No GP-based properties fabricated; these are neutral placeholders.
            "round1Role": None,
            "sameRecipeAs": None,
            "group": None,
            "forecast": None,
            "vsGinkgoPlates": None,
            "dataOrigin": "PocketBase/cfps_designs",
        })
    if not ids:
        raise ValueError("No valid experimental CFPS records in PocketBase snapshot")
    out = {
        "schema": "cfps-selection-input/v1",
        "built": meta.get("created_utc") or None,
        "input_source": "pocketbase_snapshot",
        "model_based_annotations_available": False,
        "proteins": {p: {"designs": designs[p], "groups": [], "sweeps": [],
                         "totals": {"designs": len(designs[p])}}
                     for p in CANONICAL_PROTEINS},
    }
    summary = {
        "source_type": "pocketbase_snapshot",
        "adapter_version": SOURCE_ADAPTER_VERSION,
        "record_count": len(records),
        "experimental_eligible": len(ids),
        "counts_by_protein": {p: len(designs[p]) for p in CANONICAL_PROTEINS},
        "nonexperimental_excluded": dict(sorted(excluded_types.items())),
        "author_field_missing_count": missing_author_count,
        "gp_baseline_available": False,
        "quantity_semantics": "PocketBase reagent_list stored dispensing volumes (nL)",
        "note": "No GP forecasts/roles were inferred; baseline and hybrid modes unavailable",
    }
    return out, summary


def normalize_analysis_source(path: str | Path, source_type: str) -> tuple[dict, dict]:
    """Read an immutable source file and return compatible data with provenance.

    The caller hashes both source bytes and any normalized cache file separately.
    """
    if source_type not in SOURCE_TYPES:
        raise ValueError(f"source_type must be one of {SOURCE_TYPES}; got {source_type!r}")
    source = Path(path)
    raw = json.loads(source.read_text(encoding="utf-8"))
    if source_type == "pocketbase_snapshot":
        return normalize_pocketbase_records(raw)
    if not isinstance(raw, dict) or not isinstance(raw.get("proteins"), dict):
        raise ValueError("Published snapshot requires top-level proteins mapping")
    for protein in CANONICAL_PROTEINS:
        if protein not in raw["proteins"] or not isinstance(raw["proteins"][protein].get("designs"), list):
            raise ValueError(f"Published snapshot is missing designs for {protein}")
    return raw, {
        "source_type": "published_snapshot", "adapter_version": SOURCE_ADAPTER_VERSION,
        "counts_by_protein": {p: len(raw["proteins"][p]["designs"]) for p in CANONICAL_PROTEINS},
        "gp_baseline_available": True,
        "quantity_semantics": "published recipeNl dispensing volumes (nL)",
        "note": "GP predictions and roles are source-provided, not regenerated",
    }
