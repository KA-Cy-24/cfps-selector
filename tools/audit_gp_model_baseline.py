#!/usr/bin/env python3
"""CLI audit for the published GP-informed CFPS pre-run selection; no network/API usage.

Run at repository root:
python tools/audit_gp_model_baseline.py --designs /private/clusters_designs.json --out /private/audit --mode audit
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

# Executable from checkout without pip-installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from analysis.cfps_selection import (
    audit_gp_model_snapshot, audit_gp_model_link_matrices, read_gp_model_file, coverage_metrics,
    suggest_portfolio, write_selection_overlay,
)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--designs", type=Path, required=True, help="GP-informed clusters_designs.json")
    p.add_argument("--out", type=Path, required=True, help="Private output directory (outside Git)")
    p.add_argument("--links", type=Path, help="Optional GP-informed clusters.json.gz link matrix source")
    p.add_argument("--quota", type=int, default=320, help="Working quota, NOT automatically verified")
    p.add_argument("--mode", choices=["audit", "gp_model_baseline", "hybrid", "model_free", "all"], default="audit")
    p.add_argument("--objective", choices=["facility", "maximin"], default="facility")
    a = p.parse_args(argv)
    if a.out.resolve().is_relative_to(Path(__file__).resolve().parent.parent):
        raise SystemExit("Refusing row-level outputs inside a Git checkout; choose a PRIVATE --out directory")
    data = read_gp_model_file(a.designs)
    audit = audit_gp_model_snapshot(data, working_quota=a.quota)
    manifest = audit.save(a.out, input_paths=[a.designs]+([a.links] if a.links is not None else []))
    if a.links is not None:
        link_report = audit_gp_model_link_matrices(a.links, design_audit=audit)
        (a.out/"model_link_matrix_audit.json").write_text(json.dumps(link_report,indent=2),encoding="utf-8")
    portfolios = {}
    chosen_modes = ["gp_model_baseline", "hybrid", "model_free"] if a.mode == "all" else ([a.mode] if a.mode != "audit" else [])
    for target in sorted(audit.records):
        rows = [audit.records[target][rid] for rid in sorted(audit.records[target])]
        portfolios[target] = {}
        for mode in chosen_modes:
            result = suggest_portfolio(rows, quota=a.quota, mode=mode, objective=a.objective)
            portfolios[target][mode] = result
            result["coverage"] = coverage_metrics(rows, set(result["ids"]))
        if chosen_modes:
            print(target, {mode:(v["count"],v["quota_gap"]) for mode,v in portfolios[target].items()})
        else:
            print(target, {"candidate":audit.summary["proteins"][target]["candidateDesigns"],
                           "proposed_unique":audit.summary["proteins"][target]["selectedUniqueFormulations"],
                           "gap":audit.summary["proteins"][target]["gapToWorkingQuota"]})
    write_selection_overlay(audit, portfolios, a.out/"selection_overlay_private_summary.json", include_private_ids=False)
    if chosen_modes:
        # never publish or write IDs to the source checkout without deliberate review.
        (a.out/"portfolio_ids_private.json").write_text(json.dumps({t:{mode:v["ids"] for mode,v in pp.items()} for t,pp in portfolios.items()},indent=2),encoding="utf-8")
    print("Private data directory:",a.out)
    print("Manifest input SHA256:", manifest["input_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
