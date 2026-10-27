"""Tests for source schema normalization and strategy eligibility (synthetic data only)."""
from __future__ import annotations
import copy
import json
import tempfile
from pathlib import Path
import pytest
from analysis.cfps_selection.input_sources import normalize_analysis_source, normalize_pocketbase_records
from analysis.cfps_selection.model_free_pipeline import model_free_run


def sample_source():
    def rec(rid, protein, a, kind='experimental'):
        return {'id':rid, 'author':'synthetic', 'design_json':{
            'sample_type':kind,'sample_id':rid,
            'metadata':{'protein_target':protein},
            'reagent_list':{'a':a,'nuclease_free_water':20000-a}}}
    rows=[]
    for p in ('sfGFP','PETase','Reteplase'):
        for i,v in enumerate((100.,250.,400.)):
            rows.append(rec(f'{p}_{i}',p,v))
    rows.append(rec('control','sfGFP',100,'control'))
    return {'metadata':{'designs_possibly_truncated':False,'designs_advertised_total_items':10},
            'collections':{'cfps_designs':rows}}


def test_pocketbase_adapter_three_targets_and_provenance():
    converted, meta=normalize_pocketbase_records(sample_source())
    assert meta['experimental_eligible']==9
    assert meta['nonexperimental_excluded']=={'control':1}
    assert meta['gp_baseline_available'] is False
    assert [len(converted['proteins'][p]['designs']) for p in ('sfGFP','PETase','Reteplase')]==[3,3,3]
    assert all(r['forecast'] is None and r['group'] is None and r['round1Role'] is None for block in converted['proteins'].values() for r in block['designs'])


def test_adapters_use_json_serialized_designs_and_preserve_ids():
    src=sample_source()
    for rec in src['collections']['cfps_designs']:
        rec['design_json']=json.dumps(rec['design_json'])
    normalized,_=normalize_pocketbase_records(src)
    assert normalized['proteins']['sfGFP']['designs'][0]['id']=='sfGFP_0'


def test_source_adapter_rejects_partial_and_invalid_targets():
    s=sample_source();s['metadata']['designs_possibly_truncated']=True
    with pytest.raises(ValueError,match='truncated'):normalize_pocketbase_records(s)
    s=sample_source();s['collections']['cfps_designs'][0]['design_json']['metadata']['protein_target']='unknown'
    with pytest.raises(ValueError,match='Unknown'):normalize_pocketbase_records(s)


def test_source_adapter_rejects_bad_volume_and_missing_id():
    s=sample_source();s['collections']['cfps_designs'][0]['design_json']['reagent_list']['a']=-1
    with pytest.raises(ValueError,match='Invalid'):normalize_pocketbase_records(s)
    s=sample_source();s['collections']['cfps_designs'][0]['id']=''
    with pytest.raises(ValueError,match='Missing or duplicate'):normalize_pocketbase_records(s)


def test_pocketbase_run_model_free_only():
    src=sample_source()
    with tempfile.TemporaryDirectory() as tmp:
        path=Path(tmp)/'pocketbase_snapshot.json';path.write_text(json.dumps(src))
        normalized, meta=normalize_analysis_source(path,'pocketbase_snapshot')
        adapted=Path(tmp)/'adapted.json';adapted.write_text(json.dumps(normalized))
        result=model_free_run(adapted,quota=2,modes=('model_free',),objectives=('facility',))
        assert all(list(res)==['model_free:facility'] for res in result['results'].values())
        assert all(res['model_free:facility']['n_selected']==2 for res in result['results'].values())
        with pytest.raises(ValueError,match='no GP-informed baseline'):
            model_free_run(adapted,quota=2,modes=('gp_model_baseline',))


def test_published_snapshot_shape_fail_closed():
    with tempfile.TemporaryDirectory() as tmp:
        f=Path(tmp)/'source.json';f.write_text(json.dumps({'proteins':{'sfGFP':{'designs':[]}}}))
        with pytest.raises(ValueError,match='missing designs'):
            normalize_analysis_source(f,'published_snapshot')


def test_source_rejects_invalid_source_enum():
    with pytest.raises(ValueError,match='source_type'):
        normalize_analysis_source('missing.json','not_a_source')


def test_source_rejects_wrong_reaction_volume():
    src=sample_source()
    src['collections']['cfps_designs'][0]['design_json']['reagent_list']['a']=110
    with pytest.raises(ValueError,match='total volume'):
        normalize_pocketbase_records(src)
