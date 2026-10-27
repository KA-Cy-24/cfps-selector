import json
from pathlib import Path
import pytest
from analysis.cfps_selection.decision_framework import (
    DecisionConfig,candidate_payload,select_plate,detect_partial_sets,
    verify_published_groups,EXPERIMENTAL_WELLS,RESERVED_COLUMNS)


def fake_data(n=12,batch=True):
    raw=[]
    for i in range(n):
        d={'id':f'id{i:03}','protein':'sfGFP','sampleId':f's{i}',
           'recipeNl':{'base_buffer':10000,'magnesium_glutamate':100+i*50,'nuclease_free_water':9900-i*50},
           'forecast':{'pBeat':i/n,'sd':.2+i*.01,'lnSpecCost':i/2},
           'round1Role':('first run' if i<3 else None), 'group':({'letter':'A'} if i<3 else None)}
        if batch and i in [1,2]:d['metadata']={'batch_id':'batch_A'}
        raw.append(d)
    return {'proteins':{'sfGFP':{'designs':raw,'groups':[{'letter':'A','designs':['id000','id001','id002'],
          'firstRun':'id000','checkRuns':[{'id':'id001'}]}], 'sweeps':[]},
          'PETase':{'designs':[],'groups':[],'sweeps':[]},
          'Reteplase':{'designs':[],'groups':[],'sweeps':[]}}}


def test_384_layout():
    assert RESERVED_COLUMNS==(1,2,3,4)
    assert len(EXPERIMENTAL_WELLS)==len(set(EXPERIMENTAL_WELLS))==320
    assert EXPERIMENTAL_WELLS[:3]==['A05','B05','C05']
    assert EXPERIMENTAL_WELLS[-1]=='P24'

@pytest.mark.parametrize('rep,gp,wells',[(1,0,0),(1,1,320),(2,.5,160),(4,.5,160),(8,0,0),(16,.5,160)])
def test_exact_quota_replication(rep,gp,wells):
    cfg=DecisionConfig(replicate_factor=rep,gp_share=gp)
    assert cfg.gp_unique_quota*rep==wells
    assert cfg.gp_unique_quota<=cfg.unique_capacity

@pytest.mark.parametrize('rep,gpw',[(4,3),(3,60),(4,324),(-4,None)])
def test_invalid_capacity(rep,gpw):
    with pytest.raises(ValueError):DecisionConfig(replicate_factor=rep,gp_wells=gpw).validate()


def test_published_integrity():
    parity=verify_published_groups(fake_data())
    assert parity['proteins']['sfGFP']['integrity_ok']
    assert parity['independent_cluster_equivalence'].startswith('NOT ESTABLISHED')


def test_partial_batch_with_gp_group():
    payload=candidate_payload(fake_data(),'sfGFP')
    warnings=detect_partial_sets(payload,{'id001'})
    assert any(w['type']=='submitted_batch' and w['severity']=='HIGH' for w in warnings)
    assert any(w['type']=='published_gp_group' and w['severity']=='INFO' for w in warnings)


def test_lock_and_published_only():
    payload=candidate_payload(fake_data(),'sfGFP')
    cfg=DecisionConfig(replicate_factor=32,gp_share=.5)  # capacity 10, 5+5
    result=select_plate(payload,cfg,locked_gp=['id000'],locked_model_free=['id011'])
    assert len(result['well_assignment'])==320
    assert result['by_method']['gp'].count('id000')==1
    assert 'id011' in result['by_method']['model_free']
    assert len(set(result['selected_ids']))==10
    assert result['gp_well_quota']==160
    assert result['model_free_well_quota']==160


def test_no_gp_forecasts_means_gp_forbidden():
    data=fake_data()
    for d in data['proteins']['sfGFP']['designs']:d['forecast']=None
    data['proteins']['sfGFP']['groups']=[]
    payload=candidate_payload(data,'sfGFP')
    with pytest.raises(ValueError,match='annotated'):
        select_plate(payload,DecisionConfig(replicate_factor=32,gp_share=.2))
    result=select_plate(payload,DecisionConfig(replicate_factor=32,gp_share=0))
    assert result['gp_unique_quota']==0


def test_forecast_extension_fully_labeled():
    payload=candidate_payload(fake_data(),'sfGFP')
    result=select_plate(payload,DecisionConfig(replicate_factor=32,gp_share=1))
    assert result['gp_unique_quota']==10
    assert len(result['gp_forecast_extensions'])==7


def test_duplicate_recipe_collapsed():
    data=fake_data()
    data['proteins']['sfGFP']['designs'].append({**data['proteins']['sfGFP']['designs'][0],'id':'new_dup'})
    p=candidate_payload(data,'sfGFP')
    assert p['n_distinct']==12
    assert 'new_dup' in p['candidates'][0]['aliases']


def test_batch_without_provenance_not_fabricated():
    p=candidate_payload(fake_data(batch=False),'sfGFP')
    assert p['submission_batches']==[]


def test_source_batch_expected_size_flagged():
    data=fake_data(batch=False)
    data['proteins']['sfGFP']['designs'][1]['metadata']={'batch_id':'planned_three','batch_size':3,'batch_index':0}
    p=candidate_payload(data,'sfGFP')
    batches=p['submission_batches']
    assert len(batches)==1 and batches[0]['possible_missing_source_members']==2
    warnings=detect_partial_sets(p,{'id001'})
    assert any(w['type']=='source_batch_incomplete' and w['severity']=='HIGH' for w in warnings)


def test_bad_group_membership_annotation_detected():
    data=fake_data()
    data['proteins']['sfGFP']['designs'][0]['group']={'letter':'B'}
    report=verify_published_groups(data)
    assert not report['proteins']['sfGFP']['integrity_ok']
    assert any('annotation disagreement' in x for x in report['proteins']['sfGFP']['errors'])


def test_js_python_selected_ids_match_with_synthetic_data(tmp_path):
    import shutil,subprocess
    if not shutil.which('node'):
        pytest.skip('Node is not available')
    from analysis.cfps_selection.decision_framework import DecisionConfig
    src=Path(__file__).resolve().parents[1]
    p=candidate_payload(fake_data(),'sfGFP')
    cfg=DecisionConfig(replicate_factor=32,gp_share=.5)
    expected=select_plate(p,cfg)['selected_ids']
    path=tmp_path/'payload.json'
    path.write_text(json.dumps(p))
    js=src/'repo_additions/static/cfps/decision/decision.js'
    code=f'''require({str(js)!r});const fs=require('fs');const p=JSON.parse(fs.readFileSync({str(path)!r}));
       const c=globalThis.CFPSDecisionLab.defaultConfig();c.reps=32;c.gpPercent=50;
       const got=globalThis.CFPSDecisionLab.solve({{targets:{{sfGFP:p}}}},'sfGFP',c,{{sfGFP:{{slots:{{}},excluded:[]}}}},{{}});
       console.log(JSON.stringify(got.selected))'''
    got=json.loads(subprocess.check_output(['node','-e',code],text=True))
    assert got==expected
