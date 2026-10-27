#!/usr/bin/env python3
"""Run fully reproducible pre-run CFPS selection audit against existing submissions.

Defaults keep all row-level content under gitignored data/cfps-model-free.
No login, no API, no submission, no mutation of GP source files.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from dataclasses import asdict
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analysis.cfps_selection.selection_framework import audit_gp_model_snapshot, audit_gp_model_link_matrices, input_sha256
from analysis.cfps_selection.input_sources import normalize_analysis_source, SOURCE_TYPES
from analysis.cfps_selection.model_free_pipeline import model_free_run,write_private_analysis_payload,MODEL_FREE_VERSION


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--designs',type=Path,required=True,help='PRIVATE input: GP-informed clusters_designs.json OR read-only PocketBase snapshot')
    p.add_argument('--source',choices=SOURCE_TYPES,default='published_snapshot',help='Input schema; PocketBase-only mode has no GP baseline')
    p.add_argument('--external-zip',type=Path,help='OPTIONAL external team cfps-plates.zip')
    p.add_argument('--links',type=Path,help='OPTIONAL original clusters.json.gz')
    p.add_argument('--outdir',type=Path,default=ROOT/'data/cfps-model-free',help='Private, gitignored output directory')
    p.add_argument('--quota',type=int,default=320,help='Working selection quota, VERIFY final rules')
    p.add_argument('--external-policy',choices=['flag','exclude'],default='flag')
    p.add_argument('--seed-ids-json',type=Path,help='OPTIONAL private mapping protein -> mandatory seed IDs')
    p.add_argument('--skip-baseline-audit',action='store_true')
    args=p.parse_args(argv)
    if not args.designs.is_file():p.error('Missing --designs file')
    if args.external_zip and not args.external_zip.is_file():p.error('Missing --external-zip file')
    args.outdir.mkdir(parents=True,exist_ok=True)
    seed=json.loads(args.seed_ids_json.read_text()) if args.seed_ids_json else None
    normalized, source_summary = normalize_analysis_source(args.designs,args.source)
    if args.source == 'pocketbase_snapshot':
        normalized_path=args.outdir/'normalized_candidate_pool_PRIVATE.json'
        normalized_path.write_text(json.dumps(normalized,sort_keys=True,indent=2)+'\n',encoding='utf-8')
        normalized_path.chmod(0o600)
        modes=('model_free',)
    else:
        normalized_path=args.designs
        modes=('gp_model_baseline','hybrid','model_free')
    run=model_free_run(normalized_path,args.external_zip,quota=args.quota,
                       external_policy=args.external_policy,seed_ids_by_target=seed,modes=modes)
    paths=[args.designs]+([args.external_zip] if args.external_zip else [])
    provenance={'method':'pre-experimental existing-recipe portfolio coverage',
                'source_type':args.source, 'source_adapter':source_summary,
                'working_quota_not_verified':args.quota,'external_policy':args.external_policy,
                'source_files':{p.name:input_sha256(p) for p in paths}}
    v=write_private_analysis_payload(run['results'],args.outdir/'model_free_analysis.json',
                      input_paths=paths,external_summary=run['external_summary'],provenance=provenance)
    if not args.skip_baseline_audit and args.source=='published_snapshot':
        audit=audit_gp_model_snapshot(run['data'],working_quota=args.quota)
        audit_dir=args.outdir/'baseline_audit';audit.save(audit_dir,input_paths=paths)
        if args.links:
            link=audit_gp_model_link_matrices(args.links,design_audit=audit)
            (audit_dir/'model_link_matrix_audit.json').write_text(json.dumps(link,indent=2))
    def safe_round(value, precision=6):
        return round(float(value), precision) if value is not None else None
    public_summary={'schema':'cfps-model-free-aggregate/v1', 'algorithm_version':MODEL_FREE_VERSION,
                    'privacy':'aggregate only; no recipe or participant IDs',
                    'external_overlap_summary':run['external_summary'],
                    'input_source':args.source,
                    'proteins':{target:{k:{'selected_unique':v['n_selected'],
                           'quota_gap':v['gap_to_working_quota'],
                           'mean_nearest_distance':round(v['distance']['mean'],6),
                           'p90_nearest_distance':round(v['distance']['p90'],6),
                           'factor_bin_coverage':safe_round(v['factor_bin_coverage_mean'],6),
                           'pair_bin_coverage':safe_round(v['pair_bin_coverage_mean'],6),
                           'exact_sweep_reserved_new_slots':v['sweep_reserved_count'],
                           'authors_represented':v['authors_represented'],
                           'external_overlaps_selected':v['external_overlaps_selected'],
                           'near_external_selected_heuristic':v['near_external_selected_count']}
                           for k,v in strategies.items()}
                           for target,strategies in run['results'].items()},
                    'provenance':provenance}
    (args.outdir/'aggregate_summary.json').write_text(json.dumps(public_summary,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps({'private_payload':v,**public_summary},indent=2,sort_keys=True))
    return 0

if __name__=='__main__':raise SystemExit(main())
