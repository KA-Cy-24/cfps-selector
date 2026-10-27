"""Model-free CFPS portfolio selection, exact external overlap, and private visualization export.

Scientific scope: existing submitted recipe IDs, no output GP and no assay outcomes.
Scientific methods and author contributions are documented in docs/CONTRIBUTIONS.md.

Module API mirrored verbatim inside self-contained Colab notebook (not imported there).
Private input/output: never commit recipe-level snapshots, IDs, or overlay JSON publicly.
"""
from __future__ import annotations
from collections import defaultdict, Counter
from dataclasses import dataclass, asdict
import hashlib
import io
import json
import math
import os
from pathlib import Path
from typing import Any, Mapping, Sequence
import zipfile

import numpy as np
from scipy.spatial.distance import cdist
# canonical_recipe, changed_reagents, signature_except and input_sha256 are supplied
# by selection_framework.py in the modular distribution; in Colab notebook,
# their definitions are embedded in a preceding code cell.
if __name__ != '__main__':
    # Import only in the modular package. Colab executes this exact source as
    # a __main__ notebook cell, using the already embedded audit definitions.
    from .selection_framework import canonical_recipe, changed_reagents, signature_except, input_sha256, BASELINE_ROLES

MODEL_FREE_VERSION = '1.1.0'
TARGETS_CANONICAL = ('sfGFP', 'PETase', 'Reteplase')
WATER_FIELDS = frozenset({'nuclease_free_water'})


def protein_target(text: str) -> str:
    """Canonicalize plate's protein name or its template identifier."""
    lower = str(text).casefold()
    if 'reteplase' in lower: return 'Reteplase'
    if 'petase' in lower: return 'PETase'
    if 'sfgfp' in lower or 'gfp' in lower: return 'sfGFP'
    raise ValueError(f'Cannot infer protein target: {text!r}')


def parse_external_plates_zip(path: str | Path) -> dict:
    """Read external team plate submissions, *not* import them into candidate pool.

    Handles individual plate schema and batch schema with nested plates. Experimental
    sample types only. Deduplication is per protein and exact non-water nL signature;
    external replicate factors are not unfolded into fake independent formulations.
    Target inferred from protein_target or recognized plate/sample/file labels.
    Avoids storing external authors in public summaries.
    """
    path = Path(path)
    by_target = {p: defaultdict(list) for p in TARGETS_CANONICAL}
    records_per_target = Counter()
    files = []
    with zipfile.ZipFile(path) as z:
        for info in sorted(z.infolist(),key=lambda x:x.filename):
            name=info.filename
            if info.is_dir() or not name.lower().endswith('.json') or name.startswith('__MACOSX/'):
                continue
            data=json.loads(z.read(info))
            if not isinstance(data,dict): raise ValueError(f'External JSON not a mapping: {name}')
            plates=data.get('plates',[data]); file_count=0
            if not isinstance(plates,list): raise ValueError(f'Unexpected plates object: {name}')
            for plate in plates:
                if not isinstance(plate,dict): raise ValueError(f'Invalid plate: {name}')
                samples=plate.get('samples',[])
                if not isinstance(samples,list): raise ValueError(f'Invalid samples: {name}')
                for sample in samples:
                    if sample.get('sample_type')!='experimental':continue
                    try:
                        target=protein_target(plate.get('protein_target') or plate.get('plate_id') or sample.get('sample_id') or name)
                    except ValueError:
                        target=protein_target(name)
                    recipe=sample.get('reagent_list')
                    if not isinstance(recipe,dict): raise ValueError(f'External sample missing reagent_list: {name}')
                    sig=canonical_recipe(recipe)
                    by_target[target][sig].append({'source_file':Path(name).name,'plate_id':str(plate.get('plate_id','')),'sample_id':str(sample.get('sample_id',''))})
                    file_count+=1
                    records_per_target[target]+=1
            files.append({'file':Path(name).name,'experimental_records':file_count})
    return {'signatures':{p:dict(by_target[p]) for p in TARGETS_CANONICAL},
            'summary':{'experimental_records':{p:int(records_per_target[p]) for p in TARGETS_CANONICAL},
                       'distinct_formulations':{p:len(by_target[p]) for p in TARGETS_CANONICAL},
                       'source_json_files':len(files),'file_counts':files,
                       'overlap_definition':'exact non-water dispensing nL by reagent key and protein, not confirmed chemical equivalence'},
            'sha256':input_sha256(path)}


def external_nearest_distance(data:Mapping[str,Any], external:Mapping[str,Any],
                              target:str,cfg:ModelFreeConfig|None=None)->dict[str,float]:
    """Descriptive near-overlap distance only, NEVER used for automatic exclusion.

    WARNING: assumes shared reagent-key/stock-volume semantics. This is a
    screening metric, not an assertion that different recipes are bioequivalent.
    """
    if external is None:return {}
    rows,_=unique_candidate_pool(data,target)
    others=[{'id':f'external_{i}','recipeNl':dict(sig)}
            for i,sig in enumerate(external['signatures'][target])]
    if not others:return {}
    cfg=cfg or ModelFreeConfig()
    z,_,_=feature_matrix(rows+others,cfg)
    distances=cdist(z[:len(rows)],z[len(rows):],'euclidean')
    return {r['id']:float(v) for r,v in zip(rows,distances.min(axis=1))}


def exact_external_overlap(data: Mapping[str,Any], external: Mapping[str,Any]) -> tuple[dict,dict]:
    """Return private ID-level overlaps + aggregate summaries, by target."""
    if external is None: return ({p:{} for p in TARGETS_CANONICAL}, {'status':'not_supplied'})
    matches={}; summary={}
    for target in TARGETS_CANONICAL:
        look=external['signatures'].get(target,{})
        block=data['proteins'][target]
        results={}
        selected_count=0
        distinct=set()
        for design in block['designs']:
            sig=canonical_recipe(design['recipeNl'])
            if sig in look:
                results[design['id']]=look[sig]
                distinct.add(sig)
                if design.get('round1Role') in BASELINE_ROLES: selected_count+=1
        matches[target]=results
        summary[target]={'participant_records_matching':len(results),
                         'distinct_participant_formulations_matching':len(distinct),
                         'proposed_gp_baseline_records_matching':selected_count,
                         'external_distinct_formulations':len(look)}
    return matches,summary


@dataclass(frozen=True)
class ModelFreeConfig:
    """Frozen configuration; no tuning against predicted yield.

    Quota: working capacity assumption, not official competition validation.
    max_contrast_fraction: maximum *slots* allocated to new complete 3-level
      one-factor sweeps (remaining slots chosen by geometric coverage).
    external_policy: 'flag' or 'exclude' (never auto-purge the underlying pool).
    hybrid: locked GP baseline unique recipes, with optional manual seed IDs.
    model_free: no GP roles or forecasts consulted (except in *separate* comparisons).
    """
    quota: int=320
    external_policy: str='flag'
    objective: str='facility' # facility or maximin
    feature_mode: str='scaled_plus_presence'  # scaled_volume, presence, scaled_plus_presence
    presence_weight: float=0.35
    distance_scale: str='p95'   # p95 or max
    max_contrast_fraction: float=0.20
    max_sweep_share_of_new_slots: float=0.50
    max_families_per_factor: int=5
    max_levels_per_sweep: int=5
    min_levels_per_sweep: int=3
    factor_bins: int=5
    near_external_threshold: float=0.20  # scaled-distance HEURISTIC; annotate only
    n_projection_components: int=2

    def validate(self):
        if self.quota<1: raise ValueError('quota must be positive')
        if self.external_policy not in {'flag','exclude'}: raise ValueError('external_policy must be flag|exclude')
        if self.objective not in {'facility','maximin'}: raise ValueError('objective must be facility|maximin')
        if self.feature_mode not in {'scaled_volume','presence','scaled_plus_presence'}: raise ValueError('Bad feature_mode')
        if self.distance_scale not in {'p95','max'}: raise ValueError('Bad distance_scale')
        if not 0<=self.max_contrast_fraction<=1: raise ValueError('Bad contrast fraction')
        if not 0<=self.max_sweep_share_of_new_slots<=1: raise ValueError('Bad new-slot sweep share')
        if self.presence_weight<0: raise ValueError('Bad presence weight')
        if self.max_families_per_factor<0 or self.max_levels_per_sweep<2: raise ValueError('Bad sweep caps')
        if self.min_levels_per_sweep<2 or self.min_levels_per_sweep>self.max_levels_per_sweep: raise ValueError('Bad min sweep levels')
        if self.factor_bins<2: raise ValueError('Bad factor bins')
        if self.near_external_threshold<0: raise ValueError('Bad near-external threshold')


def unique_candidate_pool(data: Mapping[str,Any],target:str, *, forbidden_ids: set[str] | None=None) -> tuple[list[dict],dict]:
    """Distinct input recipes with ID provenance, independent of GP-derived grouping.

    Representative ID is lexical minimum; the origin's round1Role is stored for
    comparison only and is never consulted for the *model_free* selection.
    """
    forbidden_ids=forbidden_ids or set()
    unique=defaultdict(list)
    for raw in data['proteins'][target]['designs']:
        key=canonical_recipe(raw['recipeNl']);unique[key].append(raw)
    rows=[]; all_records=[]
    for key, originals in sorted(unique.items(),key=lambda x: min(r['id'] for r in x[1])):
        ids=sorted(r['id'] for r in originals)
        row={'id':ids[0], 'ids':ids, 'recipeNl':dict(key),
             'author_count':len(set(str(r.get('author','')) for r in originals)),
             'external_overlap':bool(any(rid in forbidden_ids for rid in ids)),
             'baseline_ids':[r['id'] for r in originals if r.get('round1Role') in BASELINE_ROLES],
             'group_labels':sorted({r['group']['letter'] for r in originals if r.get('group')}),
             'protein':target}
        rows.append(row)
        all_records.extend(ids)
    return rows, {'records':len(all_records),'distinct':len(rows),'duplicate_excess':len(all_records)-len(rows)}


def feature_matrix(rows:list[dict],cfg:ModelFreeConfig) -> tuple[np.ndarray,list[str],dict]:
    """Explicit recipe volume feature map; no stock molarities inferred.

    Per-reagent scaling by positive p95 (or observed maximum), capped at 1 for
    larger values to limit outliers. Presence indicator adds orthogonal weight.
    Water is excluded, and features constant across target are removed.
    """
    cfg.validate()
    factors=sorted({k for r in rows for k in r['recipeNl'] if k not in WATER_FIELDS})
    if not factors:raise ValueError('No chemistry factors')
    A=np.array([[float(r['recipeNl'].get(k,0)) for k in factors] for r in rows],float)
    variable=np.ptp(A,axis=0)>1e-9
    A=A[:,variable]; factors=[f for f,b in zip(factors,variable) if b]
    if not len(factors):raise ValueError('No variable chemistry factors')
    scale=np.array([np.quantile(A[A[:,j]>0,j],0.95) if cfg.distance_scale=='p95' else max(A[:,j]) for j in range(A.shape[1])],float)
    scale=np.maximum(scale,1e-12)
    volume=np.minimum(A/scale,1.0)
    presence=(A>0).astype(float)
    if cfg.feature_mode=='scaled_volume':Z=volume; names=factors
    elif cfg.feature_mode=='presence':Z=presence; names=['presence:'+k for k in factors]
    else:
        Z=np.hstack([volume,math.sqrt(cfg.presence_weight)*presence])
        names=factors+['presence:'+k for k in factors]
    meta={'factors':factors,'scaling':dict(zip(factors,scale.tolist())),
          'distance_unit':'scaled nL per reagent (within-stock), NOT molarity',
          'n_features':int(Z.shape[1]),'exclude_balancers':sorted(WATER_FIELDS)}
    return Z,names,meta


def discover_exact_sweeps(rows:list[dict],cfg:ModelFreeConfig)->list[dict]:
    """One-factor sweeps discovered on *candidate pool* (not modeled trajectories).

    Signature groups are identical on all non-water, non-factor reagents.
    A family requires >= min_levels_per_sweep distinct factor values; it may
    overlap with another family and so cannot be counted as an independent DoE.
    """
    cfg.validate()
    all_factors=sorted({k for r in rows for k in r['recipeNl'] if k not in WATER_FIELDS})
    families=[]
    for factor in all_factors:
        grouped=defaultdict(list)
        for i,row in enumerate(rows):
            grouped[signature_except(row['recipeNl'],factor)].append(i)
        for background,idxs in grouped.items():
            levels=sorted({float(rows[i]['recipeNl'].get(factor,0)) for i in idxs})
            if len(levels)<cfg.min_levels_per_sweep: continue
            by_level=defaultdict(list)
            for i in idxs:by_level[float(rows[i]['recipeNl'].get(factor,0))].append(rows[i]['id'])
            chosen_ids=[min(by_level[level]) for level in levels]
            families.append({'factor':factor,'levels':levels,'ids':chosen_ids,
                             'n_levels':len(levels),'span_nL':levels[-1]-levels[0],
                             'background_hash':hashlib.sha256(repr(background).encode()).hexdigest()[:16]})
    return sorted(families,key=lambda f:(-f['n_levels'],-f['span_nL'],f['factor'],f['background_hash']))


def _sweep_subset(family:dict, max_levels:int) ->list[str]:
    """Pick range extremes and spread intermediate indices deterministically."""
    n=len(family['ids']);k=min(n,max_levels)
    if k==n:return family['ids'][:]
    # Index fractions might collide under rounding: fill from left if so.
    indices=sorted(set(int(round(i*(n-1)/(k-1))) for i in range(k)))
    for i in range(n):
        if len(indices)>=k:break
        if i not in indices:indices.append(i)
    return [family['ids'][i] for i in sorted(indices)]


def _select_sweeps(rows:list[dict],start:set[str],capacity:int,cfg:ModelFreeConfig,forbidden:set[str])->tuple[set[str],list[dict]]:
    """Reserve limited preexisting exact 3+level sweeps, no model statistics."""
    families=discover_exact_sweeps(rows,cfg)
    allowed={r['id'] for r in rows}-forbidden
    selected=set(start); taken=[];factor_count=Counter()
    free_slots=max(0,capacity-len(start))
    sweep_extra_limit=min(free_slots,int(math.floor(cfg.quota*cfg.max_contrast_fraction)),
                          int(math.floor(free_slots*cfg.max_sweep_share_of_new_slots)))
    spent=0
    # Ordering spans factor families; per-factor cap prevents a single dominant ingredient.
    for fam in families:
        factor=fam['factor']
        if factor_count[factor]>=cfg.max_families_per_factor: continue
        ids=_sweep_subset(fam,cfg.max_levels_per_sweep)
        if any(rid not in allowed for rid in ids):continue
        extra=[rid for rid in ids if rid not in selected]
        if len(extra)+spent>sweep_extra_limit or len(selected)+len(extra)>capacity:continue
        if not extra: continue
        selected.update(extra); spent+=len(extra)
        factor_count[factor]+=1
        taken.append({'factor':factor,'ids':ids,'levels_nL':fam['levels'],
                      'new_ids':extra,'n_new':len(extra),'background_hash':fam['background_hash']})
    return selected,taken


def _distance_and_select(rows:list[dict], selected_ids:set[str],cfg:ModelFreeConfig,quota:int,allowed_indices:set[int]) ->tuple[set[str],np.ndarray,np.ndarray,dict]:
    Z,names,meta=feature_matrix(rows,cfg)
    D=cdist(Z,Z,'euclidean')
    id_idx={r['id']:i for i,r in enumerate(rows)}
    selected=set(selected_ids)
    indices=[id_idx[i] for i in selected]
    nearest=D[:,indices].min(axis=1) if indices else np.full(len(rows),float('inf'))
    W=np.full(len(rows),1/len(rows),float)
    while len(selected)<quota and allowed_indices:
        ix=np.array(sorted(allowed_indices),dtype=int)
        if not indices:
            total=W@D
            pick=int(ix[np.argmin(total[ix])]) if cfg.objective=='facility' else int(ix[np.argmax(total[ix])])
        elif cfg.objective=='maximin':
            scores=nearest[ix]
            pick=int(ix[np.argmax(scores)])
        else:
            # marginal reduction in weighted sum of distance to nearest representative
            improvements=np.maximum(nearest[:,None]-D[:,ix],0)
            gains=W@improvements
            pick=int(ix[np.argmax(gains)])
        selected.add(rows[pick]['id']);indices.append(pick)
        nearest=np.minimum(nearest,D[:,pick]);allowed_indices.remove(pick)
    return selected,Z,D,meta


def _factor_bins(rows:list[dict], selected_ids:set[str],n_bins:int)->tuple[dict,dict]:
    factors=sorted({k for r in rows for k in r['recipeNl'] if k not in WATER_FIELDS})
    cov={}; enc={}
    for f in factors:
        values=np.array([float(r['recipeNl'].get(f,0)) for r in rows])
        positive=np.unique(values[values>0])
        if len(positive)==0: continue
        bounds=np.unique(np.quantile(positive,np.linspace(0,1,n_bins+1)))
        # Bin 0 is absence; positive bins are +1.. ; ties robustly collapsed.
        assigned=np.where(values==0,0,1+np.searchsorted(bounds[1:-1],values,side='right'))
        enc[f]=assigned
        active={int(assigned[i]) for i,r in enumerate(rows) if r['id'] in selected_ids}
        possible={int(x) for x in assigned}
        cov[f]={'represented_bins':len(active), 'possible_bins':len(possible),
                'coverage_fraction':float(len(active)/len(possible)),
                'missing_bins':sorted(possible-active),
                'has_zero_in_candidates':0 in possible}
    return cov,enc


def _pair_bin_coverage(encoded:dict,rows:list[dict],chosen:set[str],count_top:int=8)->dict:
    variability=sorted(encoded,key=lambda f:(-len(set(encoded[f])),f))[:count_top]
    selidx={i for i,r in enumerate(rows) if r['id'] in chosen}
    output={}
    for j,a in enumerate(variability):
        for b in variability[j+1:]:
            combos=set(zip(encoded[a].tolist(),encoded[b].tolist()))
            have={(int(encoded[a][i]),int(encoded[b][i])) for i in selidx}
            output[f'{a} × {b}']={'occupied_candidate_bins':len(combos),
                                  'represented_bins':len(have),
                                  'coverage_fraction':float(len(have)/len(combos))}
    return output


def exact_selected_contrast_summary(rows:list[dict], chosen:set[str],max_examples:int=36)->dict:
    """Find exact one-factor pairs among the chosen physical formulations.

    This counts pairs within existing recipe backgrounds, not independent degrees
    of freedom. Reports sampled examples only to bound browser payload size.
    """
    subset=[r for r in rows if r['id'] in chosen]
    factors=sorted({k for r in subset for k in r['recipeNl'] if k not in WATER_FIELDS})
    pairs=set();per_factor=Counter(); examples=[]
    for factor in factors:
        groups=defaultdict(list)
        for r in subset:
            groups[signature_except(r['recipeNl'],factor)].append(r)
        for members in groups.values():
            for i in range(len(members)):
                for j in range(i+1,len(members)):
                    a,b=members[i],members[j]
                    if math.isclose(float(a['recipeNl'].get(factor,0)),float(b['recipeNl'].get(factor,0)),abs_tol=1e-8):continue
                    pair=tuple(sorted([a['id'],b['id']]))
                    if pair in pairs:continue
                    pairs.add(pair);per_factor[factor]+=1
                    if len(examples)<max_examples:
                        examples.append({'a':pair[0],'b':pair[1],'factor':factor,
                                'levels_nL':sorted([float(a['recipeNl'].get(factor,0)),float(b['recipeNl'].get(factor,0))])})
    return {'unique_exact_pairs':len(pairs),'factor_counts':dict(sorted(per_factor.items())),
            'examples':examples,
            'note':'Pairs sharing families or ingredients are not statistically independent contrasts'}


def _pca_projection(Z:np.ndarray)->np.ndarray:
    """Deterministic PCA coordinates as a display ONLY; distances use full feature matrix."""
    C=Z-Z.mean(axis=0)
    if min(C.shape)<2:
        result=np.zeros((len(Z),2));result[:,0]=C[:,0] if C.shape[1] else 0
        return result
    _u,s,vt=np.linalg.svd(C,full_matrices=False)
    arr=C@vt[:2].T
    for j in range(2):
        pivot=int(np.argmax(abs(vt[j])))
        if vt[j,pivot]<0: arr[:,j]*=-1
    return arr


def _quantiles(nearest:np.ndarray)->dict:
    return {'mean':float(np.mean(nearest)),'p90':float(np.quantile(nearest,.9)),
            'max':float(np.max(nearest))}


def analyze_portfolio(data:Mapping[str,Any],target:str,*,config:ModelFreeConfig=ModelFreeConfig(),
                      mode:str='model_free',external_matches:dict|None=None,
                      seeds:Sequence[str]=(), include_geometry:bool=True,
                      external_nearest_by_id:Mapping[str,float]|None=None)->dict:
    """Full deterministic model-free/hybrid portfolio for existing submissions.

    * gp_model_baseline = read-only reconstruction of the published GP-informed selection.
    * hybrid = lock that baseline and use model-free sweeps/coverage to fill gaps.
    * model_free = no forecast/GP link roles consulted during selection.
    External exact overlaps may be flagged or excluded; locked external overlaps
    in hybrid are an error when exclude chosen (never silently violate baseline).
    """
    config.validate()
    if mode not in {'gp_model_baseline','hybrid','model_free'}: raise ValueError(f'Invalid mode {mode}')
    if target not in data['proteins']:raise ValueError('Unknown target')
    external_matches=external_matches or {}
    external_nearest_by_id=external_nearest_by_id or {}
    rows,pool_stats=unique_candidate_pool(data,target,forbidden_ids=set(external_matches))
    if config.quota>len(rows):raise ValueError(f'Quota {config.quota} exceeds distinct pool {len(rows)}')
    map_alias={ident:r['id'] for r in rows for ident in r['ids']}
    forbid={r['id'] for r in rows if r['external_overlap']} if config.external_policy=='exclude' else set()
    locked={map_alias[i] for i in seeds}
    if mode in {'gp_model_baseline','hybrid'}:
        locked.update(r['id'] for r in rows if r['baseline_ids'])
    if unknown := set(seeds)-set(map_alias):raise ValueError(f'Unknown seed IDs: {sorted(unknown)[:5]}')
    if len(locked)>config.quota:raise ValueError('Locked baseline larger than quota')
    if conflicted := locked & forbid:
        raise ValueError(f'GP baseline/manual seeds overlap external data under exclude policy ({len(conflicted)}). Use flag mode or resolve explicitly.')
    selected=set(locked);chosen_sweeps=[]
    if mode!='gp_model_baseline':
        selected,chosen_sweeps=_select_sweeps(rows,selected,config.quota,config,forbid)
    allowed={i for i,r in enumerate(rows) if r['id'] not in selected|forbid}
    # GP baseline is reconstructed only, never secretly padded.
    if mode=='gp_model_baseline':
        Z,fn,meta=feature_matrix(rows,config);D=cdist(Z,Z)
    else:
        selected,Z,D,meta=_distance_and_select(rows,selected,config,config.quota,allowed)
    id_idx={r['id']:i for i,r in enumerate(rows)}
    idx=[id_idx[i] for i in sorted(selected)]
    nearest=D[:,idx].min(axis=1) if idx else np.full(len(rows),float('nan'))
    factor_cov,bins=_factor_bins(rows,selected,config.factor_bins)
    pair_cov=_pair_bin_coverage(bins,rows,selected)
    projection=_pca_projection(Z) if include_geometry else np.zeros((len(rows),2))
    contrasts=exact_selected_contrast_summary(rows,selected)
    n_external_selected=sum(r['external_overlap'] and r['id'] in selected for r in rows)
    authors=set(); selected_authors=set()
    for d in data['proteins'][target]['designs']:
        authors.add(str(d.get('author','')))
        if map_alias[d['id']] in selected: selected_authors.add(str(d.get('author','')))
    return {
        'schema':'cfps-model-free-portfolio/v1','algorithm_version':MODEL_FREE_VERSION,
        'target':target,'mode':mode,'selection_status':'PROVISIONAL: pending manual controls, rules, lab feasibility and group decision',
        'config':asdict(config),'selected_ids':sorted(selected),
        'selected_all_submission_ids':sorted(d['id'] for d in data['proteins'][target]['designs'] if map_alias[d['id']] in selected),
        'n_selected':len(selected),'gap_to_working_quota':config.quota-len(selected),
        'locked_count':len(locked),'sweep_reserved_count':sum(x['n_new'] for x in chosen_sweeps),
        'selected_sweeps':chosen_sweeps,'external_overlaps_selected':n_external_selected,
        'external_overlaps_pool':sum(r['external_overlap'] for r in rows),
        'near_external_pool_count':sum(external_nearest_by_id.get(r['id'],float('inf'))<=config.near_external_threshold for r in rows),
        'near_external_selected_count':sum(external_nearest_by_id.get(r['id'],float('inf'))<=config.near_external_threshold for r in rows if r['id'] in selected),
        'pool':pool_stats,'authors_represented':len(selected_authors),'authors_total':len(authors),
        'exact_selected_contrasts':contrasts,
        'distance':_quantiles(nearest),'dimension_count':meta['n_features'],
        'feature_meta':meta,
        'factor_bin_coverage_mean':float(np.mean([v['coverage_fraction'] for v in factor_cov.values()])) if factor_cov else None,
        'factor_bins':factor_cov,'pair_bin_coverage_mean':float(np.mean([v['coverage_fraction'] for v in pair_cov.values()])) if pair_cov else None,
        'pair_bins':pair_cov,
        'visual_rows':[{'id':r['id'],'ids':r['ids'], 'selected':r['id'] in selected,
                        'externalOverlap':r['external_overlap'],
                        'externalNearestDistance':round(float(external_nearest_by_id[r['id']]),6) if r['id'] in external_nearest_by_id else None,
                        'point':[round(float(projection[i,j]),6) for j in range(2)],
                        'groups':r['group_labels'],
                        'nReagents':len(r['recipeNl'])}
                       for i,r in enumerate(rows)]
    }


def write_private_analysis_payload(results:dict, out:str|Path, *, input_paths:Sequence[str|Path],
                                   external_summary:dict|None=None, provenance:dict|None=None)->dict:
    """Write *private* JSON; never place in static/ or Git. Contains participant IDs."""
    path=Path(out);path.parent.mkdir(parents=True,exist_ok=True)
    body={'schema':'cfps-model-free-visualization/v1','privacy':'PRIVATE: authenticated participants only',
          'algorithm_version':MODEL_FREE_VERSION,'input_hashes':{Path(p).name:input_sha256(p) for p in input_paths},
          'targets':results,'external_summary':external_summary,
          'provenance':provenance or {},
          'caveat':'PCA projection is visual only; model-free coverage is not biological efficacy'}
    path.write_text(json.dumps(body,sort_keys=True,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    try:os.chmod(path,0o600)
    except OSError:pass
    return {'path':str(path),'sha256':input_sha256(path),'size_bytes':path.stat().st_size}


def model_free_run(data_path:str|Path, external_zip:str|Path|None=None, *,
                   quota:int=320, external_policy:str='flag',seed_ids_by_target:dict|None=None,
                   objectives:Sequence[str]=('facility','maximin'),
                   modes:Sequence[str]=('gp_model_baseline','hybrid','model_free')) ->dict:
    """Top-level reproducible pipeline: read inputs, external audit, run strategy grid.

    No network, no credential use, no row-level public output.
    """
    data_path=Path(data_path)
    data=json.loads(data_path.read_text(encoding='utf-8'))
    if data.get('model_based_annotations_available') is False and any(m != 'model_free' for m in modes):
        raise ValueError('PocketBase-only input has no GP-informed baseline or run-role annotations. ' +
                         'Use modes=("model_free",), or analyze a published GP-annotated snapshot separately.')
    external=parse_external_plates_zip(external_zip) if external_zip else None
    matches,summary=exact_external_overlap(data,external)
    results={}
    for target in TARGETS_CANONICAL:
        t_results={}
        ext_dist=external_nearest_distance(data,external,target) if external else {}
        for mode in modes:
            for obj in ('facility',) if mode=='gp_model_baseline' else objectives:
                config=ModelFreeConfig(quota=quota,external_policy=external_policy,objective=obj)
                output=analyze_portfolio(data,target,config=config,mode=mode,
                                         external_matches=matches.get(target,{}),
                                         seeds=(seed_ids_by_target or {}).get(target,()),
                                         external_nearest_by_id=ext_dist)
                t_results[f'{mode}:{obj}']=output
        results[target]=t_results
    return {'data':data,'external':external,'matches':matches,'external_summary':summary,'results':results}
