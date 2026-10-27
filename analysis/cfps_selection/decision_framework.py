"""Auditable, deterministic mixed-strategy CFPS plate planning.

Submitted formulations only; descriptive proxy objectives, NOT measured effects or
validated GP predictive optima. First four columns of a 16x24 plate are reserved,
per inspected Ginkgo CFPS example plate metadata. The 320 available physical wells
are divided into independent formulations via a configurable replicate factor.

No secrets, network access, or private-data output paths are embedded here.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping
import hashlib
import json
import math
import numpy as np
from scipy.spatial.distance import cdist

from .selection_framework import canonical_recipe, input_sha256

TARGETS = ('sfGFP', 'PETase', 'Reteplase')
ROLE_ORDER = {'first run': 0, 'check run': 1, 'sweep run': 2, 'independent run': 3}
RESERVED_COLUMNS = (1, 2, 3, 4)
EXPERIMENTAL_WELLS = [f'{chr(65+r)}{c:02d}' for c in range(5,25) for r in range(16)]
ALL_WELLS = [f'{chr(65+r)}{c:02d}' for r in range(16) for c in range(1,25)]

@dataclass(frozen=True)
class DecisionConfig:
    replicate_factor: int = 4
    gp_share: float = 0.50      # percentage normalized 0..1
    gp_wells: int | None = None # override fraction; physical wells, multiple of replicate_factor
    gp_policy: str = 'published_designations'  # or 'forecast_ranked'
    model_free_objective: str = 'maximin'        # or 'facility_location' or 'balanced'
    geometry_weight: float = 0.5  # balanced mix maximin vs facility gain
    presence_weight: float = 0.35
    gp_yield_weight: float = 0.60
    gp_cost_weight: float = 0.25
    gp_uncertainty_weight: float = 0.15
    source_mode: str = 'published_snapshot'

    def validate(self):
        if self.replicate_factor <= 0 or 320 % self.replicate_factor:
            raise ValueError('Replicate factor must be a positive divisor of 320')
        if not 0 <= self.gp_share <= 1:
            raise ValueError('gp_share must lie in [0, 1]')
        if self.gp_wells is not None and (not isinstance(self.gp_wells, int) or
                 self.gp_wells < 0 or self.gp_wells > 320 or self.gp_wells % self.replicate_factor):
            raise ValueError('gp_wells must be 0..320 and divisible by replicate_factor')
        if self.gp_policy not in {'published_designations','forecast_ranked'}:
            raise ValueError('Unknown GP policy')
        if self.model_free_objective not in {'maximin','facility_location','balanced'}:
            raise ValueError('Unknown geometry objective')
        for v in (self.geometry_weight,self.presence_weight,self.gp_yield_weight,self.gp_cost_weight,self.gp_uncertainty_weight):
            if not math.isfinite(v) or v < 0:raise ValueError('Weights must be nonnegative and finite')
        if self.geometry_weight > 1:raise ValueError('geometry_weight must be <=1')
        if self.source_mode not in {'published_snapshot','pocketbase_snapshot','pocketbase_live'}:
            raise ValueError('Unknown source mode')
        return self

    @property
    def unique_capacity(self):
        return 320//self.replicate_factor

    @property
    def gp_unique_quota(self):
        return (self.gp_wells//self.replicate_factor if self.gp_wells is not None
                else math.floor(self.unique_capacity * self.gp_share + 0.5))


def batch_id_for_design(d):
    """Only explicitly declared batch IDs qualify as submission batches."""
    meta = d.get('metadata') if isinstance(d.get('metadata'),dict) else {}
    value = d.get('batch_id') or meta.get('batch_id')
    return str(value) if value else None


def candidate_payload(data: Mapping[str,Any], target:str, *, presence_weight:float=.35) -> dict:
    """Compact authenticated web input from a frozen GP-annotated export or normalized PB data.

    Includes participant-level identifiers: NEVER write to static/ or Git.
    """
    if target not in TARGETS:raise ValueError('Unknown protein')
    block = data['proteins'][target]
    groups=block.get('groups',[])
    sweeps=block.get('sweeps',[])
    raw=block['designs']
    canonical=defaultdict(list)
    for d in raw:
        canonical[canonical_recipe(d['recipeNl'])].append(d)
    reps=[]
    for _,dupes in canonical.items():
        dupes.sort(key=lambda d:(ROLE_ORDER.get(d.get('round1Role'),99),str(d['id'])))
        representative=dupes[0]
        rep=dict(representative)
        rep['aliases']=sorted(str(d['id']) for d in dupes)
        rep['batch_ids']=sorted({x for d in dupes if (x:=batch_id_for_design(d))})
        rep['roles']=sorted({d['round1Role'] for d in dupes if d.get('round1Role') in ROLE_ORDER}, key=lambda r:ROLE_ORDER[r])
        reps.append(rep)
    reps.sort(key=lambda d:str(d['id']))
    factors=sorted({k for r in reps for k,v in r['recipeNl'].items() if k!='nuclease_free_water' and float(v)!=0})
    raw_matrix=np.array([[float(r['recipeNl'].get(k,0)) for k in factors] for r in reps],dtype=float)
    positive=np.array([np.quantile(raw_matrix[raw_matrix[:,i]>0,i],.95) if np.any(raw_matrix[:,i]>0) else 1. for i in range(len(factors))])
    scaled=np.minimum(1,raw_matrix/np.maximum(1e-9,positive))
    presence=(raw_matrix>0).astype(float)
    features=np.hstack([scaled, math.sqrt(presence_weight)*presence])
    fmax=cdist(features,features).max() if len(reps)>1 else 1.
    if fmax==0:fmax=1
    # Published groups, sweeps and batch IDs are stored separately: they are not interchangeable.
    id_alias={alias:r['id'] for r in reps for alias in r['aliases']}
    published_groups=[]
    for g in groups:
        ids=sorted({id_alias[i] for i in g.get('designs',[]) if i in id_alias})
        published_groups.append({'name':str(g['letter']),'ids':ids,'firstRun':id_alias.get(g.get('firstRun')),
                                 'checkRuns':[id_alias.get(c['id']) for c in g.get('checkRuns',[]) if c.get('id') in id_alias]})
    sweep_families=[]
    for group in sweeps:
        for idx,s in enumerate(group.get('sweeps',[])):
            ids=sorted({id_alias[x.get('id')] for x in s.get('toRun',[]) if isinstance(x,dict) and x.get('id') in id_alias})
            if len(ids)>1:sweep_families.append({'name':f'{group["reagent"]} #{idx+1}','ids':ids})
    explicit=defaultdict(list)
    expected_batch_sizes=defaultdict(set)
    for d in raw:
        batch=batch_id_for_design(d)
        meta=d.get('metadata') if isinstance(d.get('metadata'),dict) else {}
        declared=meta.get('batch_size') or d.get('batch_size')
        if batch and declared is not None:
            try: value=int(declared)
            except (TypeError,ValueError):value=0
            if value>0:expected_batch_sizes[batch].add(value)
    for r in reps:
        for batch in r['batch_ids']:explicit[batch].append(r['id'])
    has_gp=bool(groups) and any(r.get('forecast') for r in reps)
    return {'target':target,'source_built':data.get('built'),'n_original_submissions':len(raw),
            'n_distinct':len(reps),'source_has_gp_annotations':has_gp,
            'metric':{'name':'p95-scaled dispensing-volume Euclidean distance with reagent presence',
                      'presence_weight':presence_weight,'normalizer':float(fmax),'factor_keys':factors,
                      'scale_p95_nL':positive.tolist(), 'not_molarity':True},
            'candidates':[{'id':r['id'],'aliases':r['aliases'],'sampleId':r.get('sampleId'),'author':r.get('author'),
                    'recipeNl':r['recipeNl'],'features':[round(float(x),7) for x in features[i]],
                    'forecast':r.get('forecast'),'role':r.get('round1Role'),'roles':r['roles'],
                    'groups':sorted({str(d['group']['letter']) for d in raw if d['id'] in r['aliases'] and d.get('group')}),
                    'batch_ids':r['batch_ids']} for i,r in enumerate(reps)],
            'published_groups':published_groups,'sweep_families':sweep_families,
            'submission_batches':[{'name':b,'ids':sorted(set(ids)),
                        'reported_size':next(iter(expected_batch_sizes[b])) if len(expected_batch_sizes[b])==1 else None,
                        'reported_sizes_disagree':len(expected_batch_sizes[b])>1,
                        'possible_missing_source_members':(max(0,next(iter(expected_batch_sizes[b]))-len(ids)) if len(expected_batch_sizes[b])==1 else None)}
                        for b,ids in sorted(explicit.items()) if len(set(ids))>1 or expected_batch_sizes[b]]}


def _gp_scores(candidates,cfg):
    fields=[('pBeat',cfg.gp_yield_weight,1),('lnSpecCost',cfg.gp_cost_weight,-1),('sd',cfg.gp_uncertainty_weight,1)]
    score=np.zeros(len(candidates),float)
    denom=sum(w for _,w,_ in fields)
    if denom==0:return score
    for field,w,sign in fields:
        vals=np.array([(c.get('forecast') or {}).get(field) for c in candidates],dtype=object)
        nums=np.array([float(v) if v is not None and math.isfinite(float(v)) else np.nan for v in vals])
        valid=np.isfinite(nums)
        if valid.sum()<2 or np.nanmax(nums)==np.nanmin(nums):continue
        order=np.argsort(nums[valid],kind='stable')
        ranks=np.empty(order.size)
        ranks[order]=np.arange(order.size)/(order.size-1)
        part=np.full(len(candidates),.0)
        part[valid]=ranks if sign>0 else (1-ranks)
        score+=w/denom*part
    return score


def _distance_matrix(candidates):
    X=np.array([c['features'] for c in candidates],float)
    d=cdist(X,X)
    return d / max(float(d.max()),1e-12)


def detect_partial_sets(payload,selected):
    """Explicit submission batches are authoritative; GP groups/sweeps advisory."""
    selected=set(selected);warnings=[]
    for typ,field,severity in [('submitted_batch','submission_batches','HIGH'),
                               ('published_gp_group','published_groups','INFO'),
                               ('model_sweep','sweep_families','WARNING')]:
        for g in payload.get(field,[]):
            ids=set(g['ids'])
            present=sorted(ids & selected)
            missing=sorted(ids - selected)
            if present and missing:
                warnings.append({'type':typ,'severity':severity,'name':g['name'],
                                 'selected_count':len(present),'expected_count':len(ids),
                                 'selected_ids':present,'missing_ids':missing,
                                 'all_or_none_required':typ=='submitted_batch'})
    for group in payload.get('submission_batches',[]):
        if not (set(group.get('ids',[])) & selected):continue
        if group.get('reported_sizes_disagree') or (group.get('possible_missing_source_members') or 0)>0:
            warnings.append({'type':'source_batch_incomplete','severity':'HIGH','name':group['name'],
                'selected_count':len(set(group['ids']) & selected),'expected_count':group.get('reported_size'),
                'selected_ids':sorted(set(group['ids']) & selected),'missing_ids':[],
                'missing_source_count':group.get('possible_missing_source_members'),
                'all_or_none_required':True})
    return warnings


def select_plate(payload:dict, cfg:DecisionConfig=DecisionConfig(),
                 *, locked_gp=(),locked_model_free=(),excluded=()):
    """Exact GP/MF quota with locked manual choices. No implicit fallback across methods.

    Published GP roles are used *as published*. A GP quota beyond those roles
    is extended by forecast ranking and recorded explicitly. Unannotated PB data
    is rejected for GP share > 0. Returns a suggestion, never execution-ready.
    """
    cfg.validate(); candidates=payload['candidates'];n=len(candidates)
    capacity=cfg.unique_capacity
    if n<capacity:raise ValueError(f'Only {n} distinct recipes for capacity {capacity}')
    gp_capacity=cfg.gp_unique_quota;mf_capacity=capacity-gp_capacity
    if gp_capacity and not payload.get('source_has_gp_annotations'):raise ValueError('GP quota requires annotated published snapshot')
    ids={c['id']:i for i,c in enumerate(candidates)}
    locked_gp=set(locked_gp);locked_model_free=set(locked_model_free);excluded=set(excluded)
    invalid=(locked_gp | locked_model_free | excluded)-ids.keys()
    if invalid:raise ValueError(f'Unknown recipe ID(s): {sorted(invalid)[:5]}')
    if (locked_gp & locked_model_free) or ((locked_gp|locked_model_free)&excluded):raise ValueError('Conflicting manual locks')
    if len(locked_gp)>gp_capacity or len(locked_model_free)>mf_capacity:raise ValueError('Manual locks exceed method quota')
    if n-len(excluded)<capacity:raise ValueError('Too few eligible recipes')
    gp_scores=_gp_scores(candidates,cfg)
    picked_gp=list(sorted(locked_gp))
    chosen=set(locked_gp)|set(locked_model_free)
    gp_extension=[]
    available=[i for i,c in enumerate(candidates) if c['id'] not in chosen|excluded]
    if cfg.gp_policy=='published_designations':
        available.sort(key=lambda i:(min((ROLE_ORDER.get(role,99) for role in candidates[i].get('roles',[])),default=99),
                                      -gp_scores[i],candidates[i]['id']))
    else:available.sort(key=lambda i:(-gp_scores[i],candidates[i]['id']))
    for i in available:
        if len(picked_gp)>=gp_capacity:break
        c=candidates[i]
        picked_gp.append(c['id']);chosen.add(c['id'])
        if not c.get('roles'):gp_extension.append(c['id'])
    if len(picked_gp)!=gp_capacity:raise ValueError('Unable to satisfy GP quota')
    # Greedy normalized chemical-distance coverage. Seeded by GP+manual selections.
    D=_distance_matrix(candidates)
    idx=np.array([ids[i] for i in chosen],int)
    nearest=np.min(D[:,idx],axis=1) if len(idx) else np.full(n,1.,float)
    picked_free=sorted(locked_model_free)
    while len(picked_free)<mf_capacity:
        eligible=[i for i,c in enumerate(candidates) if c['id'] not in chosen|excluded]
        if not eligible:raise ValueError('Insufficient distinct eligible candidate recipes')
        if cfg.model_free_objective=='maximin':scores=nearest[eligible]
        else:
            facilities=np.maximum(0.,nearest[None,:]-D[eligible,:]).mean(axis=1)
            if cfg.model_free_objective=='facility_location':scores=facilities
            else:scores=cfg.geometry_weight*nearest[eligible]+(1-cfg.geometry_weight)*facilities
        best=eligible[min(range(len(eligible)),key=lambda j:(-float(scores[j]),candidates[eligible[j]]['id']))]
        rid=candidates[best]['id'];picked_free.append(rid);chosen.add(rid)
        nearest=np.minimum(nearest,D[:,best])
    # Column-major physical order (as in supplied 384-well plate examples).
    assignment=[];ordered=[('gp',i) for i in picked_gp]+[('model_free',i) for i in picked_free]
    physical=[x for category,rid in ordered for x in range(cfg.replicate_factor)]
    if len(physical)!=320:raise AssertionError('320 available wells must be mapped')
    for i,(method,rid) in enumerate(ordered):
        for replicate in range(cfg.replicate_factor):
            assignment.append({'well':EXPERIMENTAL_WELLS[i*cfg.replicate_factor+replicate],
                               'id':rid,'method':method,'replicate':replicate+1})
    warnings=detect_partial_sets(payload,chosen)
    return {'schema':'cfps-decision-plan/v1','target':payload['target'],'status':'PROVISIONAL, requires verification of plate rules, batches, physical layout and approval',
            'config':asdict(cfg),'unique_capacity':capacity,'gp_unique_quota':gp_capacity,
            'model_free_unique_quota':mf_capacity,'gp_well_quota':gp_capacity*cfg.replicate_factor,
            'model_free_well_quota':mf_capacity*cfg.replicate_factor,
            'selected_ids':sorted(chosen),'by_method':{'gp':picked_gp,'model_free':picked_free},
            'gp_forecast_extensions':gp_extension,
            'reserved_columns':list(RESERVED_COLUMNS),'experimental_wells':EXPERIMENTAL_WELLS,
            'well_assignment':assignment,'warnings':warnings,
            'partial_submitted_batches':sum(w['type'] in ('submitted_batch','source_batch_incomplete') for w in warnings),
            'coverage_mean_nearest_distance':float(np.mean(nearest)),
            'distance_not_empirical_biology':True}


def verify_published_groups(data:Mapping[str,Any]) -> dict:
    """Internal GP export integrity + explicit non-equivalence of independent clustering.

    Independent raw-volume Pearson clusters are NOT reconstructing GP shared-error
    links, hence equality cannot be asserted without the original GP algorithm.
    """
    out={'source':'original GP-annotated clusters_designs.json',
         'independent_cluster_equivalence':'NOT ESTABLISHED — algorithms represent different quantities',
         'published_designations_directly_available':True,'proteins':{}}
    for target in TARGETS:
        block=data['proteins'][target];all_ids={d['id'] for d in block['designs']}
        errs=[];count=0
        for g in block['groups']:
            members=set(g['designs']);count+=len(members)
            if len(members)!=len(g['designs']):errs.append(f'duplicate IDs group {g["letter"]}')
            if members-all_ids:errs.append(f'unknown design IDs group {g["letter"]}')
            if g.get('firstRun') not in members:errs.append(f'first run absent from {g["letter"]}')
            for c in g.get('checkRuns',[]):
                if c['id'] not in members:errs.append(f'check run absent from {g["letter"]}')
        source_group={rid:g['letter'] for g in block['groups'] for rid in g['designs']}
        for design in block['designs']:
            annotated=(design.get('group') or {}).get('letter')
            if annotated != source_group.get(design['id']):
                errs.append(f'group membership annotation disagreement for {design["id"]}')
        out['proteins'][target]={'n_published_groups':len(block['groups']),
                                 'n_group_memberships':count,'errors':errs,'integrity_ok':not errs}
    return out


def write_private_input(data:Mapping[str,Any],path:str|Path,*,provenance:dict|None=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    obj={'schema':'cfps-decision-input/v1','privacy':'PRIVATE; DO NOT COMMIT TO GIT OR STATIC',
         'provenance':provenance or {},'layout':{'n_rows':16,'n_columns':24,'reserved_columns':[1,2,3,4],
                    'experimental_columns':list(range(5,25)),'array_order':'column',
                    'source':'inspected cfps-plates.zip example plate metadata'},
         'cluster_parity':verify_published_groups(data) if all(data.get('proteins',{}).get(t,{}).get('groups') for t in TARGETS) else {'independent_cluster_equivalence':'NOT ESTABLISHED'},
         'targets':{t:candidate_payload(data,t) for t in TARGETS}}
    path.write_text(json.dumps(obj,ensure_ascii=False,separators=(',',':'),allow_nan=False),encoding='utf8')
    return {'path':str(path),'sha256':input_sha256(path),'bytes':path.stat().st_size}
