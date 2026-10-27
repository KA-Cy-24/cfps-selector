#!/usr/bin/env python3
"""Private-input, read-only decision planning CLI (never a robot-ready plan)."""
from __future__ import annotations
import argparse
from dataclasses import replace
from pathlib import Path
import json
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analysis.cfps_selection.decision_framework import DecisionConfig,TARGETS,write_private_input,select_plate,candidate_payload,verify_published_groups
from analysis.cfps_selection.input_sources import normalize_analysis_source
from analysis.cfps_selection.selection_framework import input_sha256


def run(source_path:str, out_dir:str, *, source_type:str='published_snapshot',
        replicate_factor:int=4, gp_share:float=.5, gp_wells:int|None=None,
        gp_policy:str='published_designations', objective:str='maximin')->dict:
    cfg=DecisionConfig(replicate_factor=replicate_factor, gp_share=gp_share,
                       gp_wells=gp_wells,gp_policy=gp_policy,
                       model_free_objective=objective,
                       source_mode=source_type)
    cfg.validate()
    data,summary=normalize_analysis_source(source_path,source_type)
    dest=Path(out_dir);dest.mkdir(parents=True,exist_ok=True)
    (dest/'.gitignore').write_text('*\n!.gitignore\n')
    provenance={'input_source_type':source_type,'source_sha256':input_sha256(source_path),
                'normalization':summary,'warning':'PRIVATE — never commit participant-level output'}
    generated=write_private_input(data,dest/'cfps_decision_input_private.json',provenance=provenance)
    result={'schema':'cfps-decision-plan-collection/v1','privacy':'PRIVATE DO NOT COMMIT',
            'provenance':provenance,'input':generated,
            'cluster_parity':verify_published_groups(data) if source_type=='published_snapshot' else None,
            'plans':{}}
    for target in TARGETS:
        payload=candidate_payload(data,target,presence_weight=cfg.presence_weight)
        result['plans'][target]=select_plate(payload,cfg)
    (dest/'decision_plans_private.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))
    warnings=[dict(target=t,**w) for t,p in result['plans'].items() for w in p['warnings']]
    (dest/'decision_warnings_private.json').write_text(json.dumps(warnings,ensure_ascii=False,indent=2,allow_nan=False))
    with (dest/'decision_well_assignments_private.csv').open('w') as f:
        f.write('protein,well,submission_id,selection_method,replicate,provisional\n')
        for t,p in result['plans'].items():
            for a in p['well_assignment']:
                f.write(f'{t},{a["well"]},{a["id"]},{a["method"]},{a["replicate"]},true\n')
    print(json.dumps({'input':str(generated['path']), 'unique_conditions_per_plate':cfg.unique_capacity,
                      'GP_conditions_per_plate':cfg.gp_unique_quota,
                      'warnings':{t:len(result['plans'][t]['warnings']) for t in TARGETS},
                      'privacy':'PRIVATE; not committed'},indent=2))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',choices=['published_snapshot','pocketbase_snapshot'],default='published_snapshot')
    p.add_argument('--input',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--replicates',type=int,default=4)
    p.add_argument('--gp-share',type=float,default=.5)
    p.add_argument('--gp-wells',type=int)
    p.add_argument('--gp-policy',choices=['published_designations','forecast_ranked'],default='published_designations')
    p.add_argument('--objective',choices=['maximin','facility_location','balanced'],default='maximin')
    a=p.parse_args()
    run(a.input,a.output,source_type=a.source,replicate_factor=a.replicates,gp_share=a.gp_share,
        gp_wells=a.gp_wells,gp_policy=a.gp_policy,objective=a.objective)
