"""Model-free pre-run selection tests: synthetic data; no live site or private snapshots."""
import copy
import json
from pathlib import Path
import zipfile

import pytest
from analysis.cfps_selection.model_free_pipeline import (
    ModelFreeConfig,analyze_portfolio,canonical_recipe,discover_exact_sweeps,
    exact_external_overlap,model_free_run,parse_external_plates_zip,
    protein_target,unique_candidate_pool,feature_matrix,exact_selected_contrast_summary,
)


def sample_data():
    records=[]
    for i in range(8):
        r={'mg':100+i*25,'peg':0 if i<3 else (i-2)*100,
           'nuclease_free_water':20000-100-i*25-(0 if i<3 else (i-2)*100)}
        records.append({'id':f'id{i:02d}','protein':'sfGFP','recipeNl':r,
                        'round1Role':'first run' if i==0 else 'follows its group',
                        'forecast':{'pBeat':i/100,'lnSpecCost':100*i},'group':None,'author':f'person{i%3}'})
    # A deliberate exact duplicate with different author ID and balancing water.
    dup=copy.deepcopy(records[1]);dup.update(id='dupe',author='another',round1Role='sweep run');records.append(dup)
    return {'proteins':{'sfGFP':{'designs':records},'PETase':{'designs':[]},'Reteplase':{'designs':[]}}}


def test_target_aliases():
    assert protein_target('sfgfp-jewett-m9134733')=='sfGFP'
    assert protein_target('reteplase-m9307400')=='Reteplase'
    assert protein_target('R1-petase-p1')=='PETase'
    with pytest.raises(ValueError):protein_target('unrelated')


def test_unique_counts_and_water():
    rows,stats=unique_candidate_pool(sample_data(),'sfGFP')
    assert stats['records']==9 and stats['distinct']==8
    assert canonical_recipe({'A':100,'nuclease_free_water':19900})==canonical_recipe({'A':100,'nuclease_free_water':19600})
    assert any(len(r['ids'])==2 for r in rows)


def test_exact_three_level_sweep_and_pairs():
    rows,_=unique_candidate_pool(sample_data(),'sfGFP')
    cfg=ModelFreeConfig(quota=5)
    families=discover_exact_sweeps(rows,cfg)
    assert any(f['factor']=='mg' and f['n_levels']>=3 for f in families)
    summary=exact_selected_contrast_summary(rows,{r['id'] for r in rows})
    assert summary['unique_exact_pairs']>0


def test_mode_model_free_independent_of_forecasts_and_roles():
    data=sample_data();cfg=ModelFreeConfig(quota=5,max_contrast_fraction=0.4)
    a=analyze_portfolio(data,'sfGFP',config=cfg,mode='model_free')
    changed=copy.deepcopy(data)
    for row in changed['proteins']['sfGFP']['designs']:
        row['forecast']={'pBeat':-1,'lnSpecCost':-99999}
        row['group']={'letter':'CHANGE'}
        row['round1Role']='independent run'
    b=analyze_portfolio(changed,'sfGFP',config=cfg,mode='model_free')
    assert a['selected_ids']==b['selected_ids']
    assert a['distance']==b['distance']


def test_fixed_seed_honored_and_quota():
    data=sample_data();cfg=ModelFreeConfig(quota=5)
    a=analyze_portfolio(data,'sfGFP',config=cfg,mode='model_free',seeds=['id07'])
    assert a['n_selected']==5 and 'id07' in a['selected_ids']
    assert a['gap_to_working_quota']==0


def test_baseline_not_padded():
    data=sample_data();cfg=ModelFreeConfig(quota=5)
    result=analyze_portfolio(data,'sfGFP',config=cfg,mode='gp_model_baseline')
    assert result['n_selected']==2 and result['gap_to_working_quota']==3


def test_hybrid_preserves_baseline():
    data=sample_data();cfg=ModelFreeConfig(quota=5)
    baseline=analyze_portfolio(data,'sfGFP',config=cfg,mode='gp_model_baseline')
    hybrid=analyze_portfolio(data,'sfGFP',config=cfg,mode='hybrid')
    assert set(baseline['selected_ids']).issubset(hybrid['selected_ids'])
    assert hybrid['n_selected']==5


def test_external_flag_and_exclude():
    data=sample_data();match={'id07':[{'file':'abc'}]}
    a=analyze_portfolio(data,'sfGFP',config=ModelFreeConfig(quota=5,external_policy='exclude'),
                        mode='model_free',external_matches=match)
    assert 'id07' not in a['selected_ids']
    with pytest.raises(ValueError,match='external'):
        analyze_portfolio(data,'sfGFP',config=ModelFreeConfig(quota=5,external_policy='exclude'),
                          mode='model_free',external_matches=match,seeds=['id07'])


def test_external_zip_batch_and_single(tmp_path):
    data=sample_data()
    match=data['proteins']['sfGFP']['designs'][0]
    recipe=match['recipeNl']
    z=tmp_path/'other_teams.zip'
    with zipfile.ZipFile(z,'w') as archive:
        archive.writestr('cfps-plates/team_plate.json',json.dumps({'protein_target':'sfgfp-jewett-m9134733','samples':[{'sample_type':'experimental','sample_id':'xyz','reagent_list':recipe}]}))
        archive.writestr('cfps-plates/more.json',json.dumps({'plates':[{'plate_id':'R1-petase-p2','samples':[{'sample_type':'experimental','sample_id':'other','reagent_list':{'a':10}}]}]}))
        archive.writestr('__MACOSX/cfps-plates/._foo.json','garbage')
    other=parse_external_plates_zip(z)
    found,stats=exact_external_overlap(data,other)
    assert stats['sfGFP']['participant_records_matching']==1
    assert stats['PETase']['participant_records_matching']==0
    assert other['summary']['experimental_records']=={'sfGFP':1,'PETase':1,'Reteplase':0}


def test_features_bounded_and_finite():
    rows,_=unique_candidate_pool(sample_data(),'sfGFP')
    z,names,meta=feature_matrix(rows,ModelFreeConfig(quota=5))
    assert z.shape[0]==8 and len(names)==z.shape[1]
    assert (z>=0).all() and (z<=1).all()
    assert 'nuclease_free_water' not in names
    assert meta['distance_unit'].startswith('scaled nL')


def test_deterministic():
    data=sample_data();cfg=ModelFreeConfig(quota=5,objective='facility')
    a=analyze_portfolio(data,'sfGFP',config=cfg,mode='model_free')
    b=analyze_portfolio(data,'sfGFP',config=cfg,mode='model_free')
    assert a['selected_ids']==b['selected_ids']
    assert a['factor_bins']==b['factor_bins']


def test_model_free_does_not_invent_recipes():
    data=sample_data();cfg=ModelFreeConfig(quota=5)
    result=analyze_portfolio(data,'sfGFP',config=cfg,mode='model_free')
    assert set(result['selected_all_submission_ids']).issubset({r['id'] for r in data['proteins']['sfGFP']['designs']})
    assert set(result['selected_ids']).issubset({r['id'] for r in data['proteins']['sfGFP']['designs']})


def test_config_validations():
    for config in (ModelFreeConfig(quota=0),ModelFreeConfig(external_policy='bad'),
                   ModelFreeConfig(max_contrast_fraction=2),ModelFreeConfig(presence_weight=-1)):
        with pytest.raises(ValueError):config.validate()


def test_model_free_visualization_fields_present():
    result=analyze_portfolio(sample_data(),'sfGFP',config=ModelFreeConfig(quota=5),mode='model_free')
    assert result['schema']=='cfps-model-free-portfolio/v1'
    assert len(result['visual_rows'])==result['pool']['distinct']
    assert all('point' in item and len(item['point'])==2 for item in result['visual_rows'])
