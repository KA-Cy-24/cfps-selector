#!/usr/bin/env python3
"""Generate the self-contained decision Colab from the exact mirrored source files."""
from pathlib import Path
import nbformat as nbf
r=Path(__file__).resolve().parents[1]
N=nbf.v4.new_notebook()
M=lambda s:nbf.v4.new_markdown_cell(s)
C=lambda s:nbf.v4.new_code_cell(s)
N.cells=[M('''# CFPS · Publication-ready Plate Decision Lab

**Purpose:** Combine the original GP-informed submitted-recipe designations and a model-free reagent-volume distance strategy under exact physical-well and replicate constraints, then build auditable **provisional** maps for sfGFP, PETase and Reteplase. 

**Plate evidence:** inspected Ginkgo CFPS example plate JSONs use 16×24 wells, reserve columns 1–4, expose columns 5–24 (320 experimental wells), and set `replicate_factor=4`, i.e. **80 independent formulations with four technical replicates**. These assumptions are configurable, and actual laboratory approval is still required. The three maps are visualization/planning artifacts, **not robot execution instructions**. Ginkgo may randomize sample placement in execution.

**Sources:** `published_snapshot` loads the original GP-annotated `clusters_designs.json` (immutable input). `pocketbase_snapshot` loads a private authorized raw PocketBase JSON snapshot; `pocketbase_live` fetches `cfps_designs` through a supplied authorized PocketBase token or password. PocketBase-only data has *no* published GP roles or forecasts: its allowed GP share is zero. No GP forecast or cluster is fabricated. The notebook records SHA-256 hashes.

**Contributions:** Original GP-informed forecasts/group designations: Clinton Monteiro; model-free geometry/decision framework: Kourosh Afshinjoo; application infrastructure: Ronan Chance. Methodological terminology in the implementation remains objective.

**Privacy:** all participant-level outputs, the HTML with private data loader and private JSON, are kept in Colab's private runtime and downloaded as a ZIP. **Do not commit private ZIP contents or generated `decision_input_private.json` to Git.**'''),
C('''%pip -q install numpy scipy requests nbformat pytest
from pathlib import Path
import os, json, hashlib, datetime, getpass, zipfile
PRIVATE = Path('/content/cfps_decision_private')
PRIVATE.mkdir(parents=True, exist_ok=True)
ROOT = Path('/content/cfps_decision_runtime')
ROOT.mkdir(parents=True, exist_ok=True)
'''),
M('''## 1 · Change inputs and strategy parameters here

`gp_wells=None` means the percentage slider controls the final quota. `gp_wells=160` means *exactly* 160 physical wells reserved for GP-informed selections. The remaining experimental wells are filled model-free. Replicate factor must divide 320, and explicit GP wells must divide the replicate factor.

The **published-designations** policy uses the original `round1Role` ordering for GP candidates; if its quota requires more candidates than annotated original roles, any additional selections are ranked by GP forecasts and **explicitly flagged**. The `forecast_ranked` policy ignores role priority and uses rank-normalized forecast criteria. These criteria are not calibrated biological utilities.'''),
C('''#@title Input source and experimental design
INPUT_SOURCE = 'published_snapshot'  # 'published_snapshot' | 'pocketbase_snapshot' | 'pocketbase_live'
INPUT_PATH = '/content/clusters_designs.json'  # or PocketBase JSON path
PB_URL = 'https://db.rcdonovan.com'  # only used for pocketbase_live
PB_AUTH_COLLECTION = '_superusers'  # must match YOUR authorized deployment
REPLICATE_FACTOR = 4  # inspected Ginkgo example default (80 unique / protein plate)
GP_SHARE_PERCENT = 50  # set 0 for PocketBase-only input
GP_WELLS = None  # integer physical wells overrides GP_SHARE_PERCENT (e.g., 160)
GP_POLICY = 'published_designations'  # or 'forecast_ranked'
MODEL_FREE_CRITERION = 'maximin'  # or 'facility_location' / 'balanced'
BALANCED_MAXIMIN_WEIGHT = 0.5
GP_YIELD_WEIGHT = 0.60
GP_COST_WEIGHT = 0.25
GP_UNCERTAINTY_WEIGHT = 0.15
'''),
M('''## 2 · Self-contained implementation

The code below writes unmodified source modules bundled into this notebook to a temporary Colab package and imports them for this session. **No GitHub checkout or external source download is required.** Do not edit generated code in the output cells; edit the mirrored source module and rebuild the notebook in the repository.''')]
for rel in ['analysis/__init__.py','analysis/cfps_selection/__init__.py','analysis/cfps_selection/selection_framework.py','analysis/cfps_selection/input_sources.py','analysis/cfps_selection/decision_framework.py','tools/cfps_decision_cli.py']:
 src=(r/rel).read_text()
 N.cells.append(C(f"rel={rel!r}\npath=ROOT/rel\npath.parent.mkdir(parents=True,exist_ok=True)\npath.write_text({src!r},encoding='utf-8')\nprint('Embedded module:',rel,'bytes:',path.stat().st_size)"))
N.cells += [C('''import sys
sys.path.insert(0,str(ROOT))
from analysis.cfps_selection.input_sources import normalize_analysis_source
from analysis.cfps_selection.decision_framework import (
    DecisionConfig,TARGETS,candidate_payload,select_plate,
    verify_published_groups,write_private_input)
from analysis.cfps_selection.selection_framework import input_sha256
'''),
M('''## 3 · Acquire an authorized, immutable snapshot

For live PocketBase, provide a token through Colab Secrets (`PB_TOKEN`) **or** authenticate interactively using email/password. This cell never prints credentials or writes tokens to disk. `PB_TOKEN` is not a Keycloak token; do not use a Keycloak confidential-client secret for PocketBase authentication.'''),
C('''from google.colab import files
if INPUT_SOURCE not in {'published_snapshot','pocketbase_snapshot','pocketbase_live'}:
    raise ValueError('Unsupported source mode')
if INPUT_SOURCE == 'pocketbase_live':
    import requests
    try:
        from google.colab import userdata
    except ImportError:
        userdata = None
    def private_value(name):
        x=os.environ.get(name,'')
        if not x and userdata is not None:
            try: x=userdata.get(name) or ''
            except Exception: pass
        return x
    token=private_value('PB_TOKEN')
    if not token:
        email=private_value('PB_EMAIL') or input('PocketBase authorized email: ')
        pw=private_value('PB_PASSWORD') or getpass.getpass('PocketBase password: ')
        reply=requests.post(f'{PB_URL.rstrip("/")}/api/collections/{PB_AUTH_COLLECTION}/auth-with-password',
                 json={'identity':email,'password':pw},timeout=30)
        reply.raise_for_status()
        token=reply.json()['token']
        del pw
    http=requests.Session()
    http.headers.update({'Authorization':token,'Accept':'application/json'})
    records=[];page=1;expected=None
    while True:
        reply=http.get(f'{PB_URL.rstrip("/")}/api/collections/cfps_designs/records',
                       params={'page':page,'perPage':200,'sort':'id'},timeout=60)
        reply.raise_for_status()
        body=reply.json()
        if expected is None: expected=int(body['totalItems'])
        batch=body['items']; records.extend(batch)
        if not batch or len(records)>=expected:break
        page+=1
    assert len(records)==expected, 'Incomplete paginated fetch: FAIL CLOSED'
    assert len({str(r['id']) for r in records})==len(records),'Duplicate immutable record IDs'
    source_path=PRIVATE/'pocketbase_snapshot_private.json'
    source_path.write_text(json.dumps({'metadata':{'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'designs_advertised_total_items':expected,'designs_possibly_truncated':False,'source':'PocketBase/cfps_designs'},
      'collections':{'cfps_designs':records}},ensure_ascii=False),encoding='utf8')
    del token
    source_type='pocketbase_snapshot'
else:
    source_type=INPUT_SOURCE
    source_path=Path(INPUT_PATH)
    if not source_path.exists():
        print('Upload your authorized, private source JSON file:')
        up=files.upload()
        if not up:raise ValueError('No source uploaded')
        source_path=PRIVATE/next(iter(up))
        source_path.write_bytes(next(iter(up.values())))
print('Input mode:',INPUT_SOURCE,'input SHA-256:',input_sha256(source_path))
'''),
M('''## 4 · Validate input, original groups, quotas, and source provenance'''),
C('''data,adapter_summary=normalize_analysis_source(source_path,source_type)
config=DecisionConfig(replicate_factor=REPLICATE_FACTOR,
    gp_share=GP_SHARE_PERCENT/100,gp_wells=GP_WELLS,gp_policy=GP_POLICY,
    model_free_objective=MODEL_FREE_CRITERION,geometry_weight=BALANCED_MAXIMIN_WEIGHT,
    gp_yield_weight=GP_YIELD_WEIGHT,gp_cost_weight=GP_COST_WEIGHT,
    gp_uncertainty_weight=GP_UNCERTAINTY_WEIGHT,source_mode=INPUT_SOURCE).validate()
print('Input normalization:',json.dumps(adapter_summary,indent=2))
if INPUT_SOURCE=='published_snapshot':
    integrity=verify_published_groups(data)
    print('Published group consistency:',json.dumps(integrity,indent=2))
    assert all(s['integrity_ok'] for s in integrity['proteins'].values()),'Published group references inconsistent'
else:
    integrity={'independent_cluster_equivalence':'GP groups are absent; no equivalence asserted'}
    assert config.gp_unique_quota==0, 'PocketBase-only mode CANNOT allocate GP quota'
print('Capacity:',config.unique_capacity,'unique formulations per plate;',config.replicate_factor,'replicates per formulation')
print('Requested allocations: GP',config.gp_unique_quota,'model-free',config.unique_capacity-config.gp_unique_quota)
'''),
M('''## 5 · Produce provisional decisions and all 384-well visualizations'''),
C('''from dataclasses import asdict
provenance={'source_mode':INPUT_SOURCE,'source_sha256':input_sha256(source_path),
            'adapter_summary':adapter_summary,'published_group_integrity':integrity,
            'warning':'Participant-level private input; do not publish'}
source_result=write_private_input(data,PRIVATE/'cfps_decision_input_private.json',provenance=provenance)
plans={t:select_plate(candidate_payload(data,t,presence_weight=config.presence_weight),config) for t in TARGETS}
output={'schema':'cfps-decision-plan-collection/v1','privacy':'PRIVATE; DO NOT COMMIT',
        'provenance':provenance,'plans':plans}
(PRIVATE/'decision_plans_private.json').write_text(json.dumps(output,ensure_ascii=False,indent=2,allow_nan=False))
rows=['protein,well,submission_id,method,replicate,provisional']
for target,plan in plans.items():
    print(target,'GP:',plan['gp_unique_quota'],'Model-free:',plan['model_free_unique_quota'],
          'GP forecast-only extensions:',len(plan['gp_forecast_extensions']),
          'Partial explicit batches:',plan['partial_submitted_batches'],
          'All partial analytical families:',len(plan['warnings']))
    for a in plan['well_assignment']:
        rows.append(','.join(map(str,[target,a['well'],a['id'],a['method'],a['replicate'],'true'])))
(PRIVATE/'decision_well_assignments_private.csv').write_text('\\n'.join(rows)+'\\n')
(PRIVATE/'decision_warnings_private.json').write_text(json.dumps({t:p['warnings'] for t,p in plans.items()},indent=2,ensure_ascii=False))
'''),
M('''## 6 · Offline interactive dashboard

The HTML dashboard is generated with the *same* JavaScript and CSS as the authenticated SvelteKit route. Load `cfps_decision_input_private.json` using its on-page file picker, or use the embedded preload version below (contains private submissions; never commit it). The dashboard allows manual well assignments and download of edited plans.''')]
for rel,var in [('repo_additions/static/cfps/decision/decision.js','viewer_js'),('repo_additions/static/cfps/decision/decision.css','viewer_css')]:
 src=(r/rel).read_text()
 N.cells.append(C(f'{var}={src!r}\nprint("Embedded client asset:",{rel!r},"length",len({var}))'))
N.cells += [C('''from IPython.display import HTML,display
input_json=(PRIVATE/'cfps_decision_input_private.json').read_text()
# JSON safely embedded as <script type=application/json>; encoded angle bracket avoids script-tag injection.
input_json=input_json.replace('<','\\\\u003c').replace('>','\\\\u003e').replace('&','\\\\u0026')
page=f"""<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CFPS Decision Lab · Private</title><style>{viewer_css}</style></head>
<body style="margin:0"><div id="cfps-root"></div><script type="application/json" id="private-data">{input_json}</script>
<script>{viewer_js}</script><script>window.CFPSDecisionLab.mount(document.getElementById('cfps-root'),{{data:JSON.parse(document.getElementById('private-data').textContent)}});</script></body></html>"""
(PRIVATE/'cfps_decision_live_preview_PRIVATE.html').write_text(page,encoding='utf8')
print('Preview saved:',PRIVATE/'cfps_decision_live_preview_PRIVATE.html')
print('Open the downloaded HTML locally for manual editing, then download edited plans from its Export controls.')
'''),
M('''## 7 · Download private results (do NOT commit)

Includes the normalized decision input, machine-readable selected IDs, a 960-row experimental-well CSV across the three plates, full warnings, input hash, and a self-contained interactive HTML viewer. **No live Keycloak secret is needed for local notebook execution.**'''),
C('''archive=Path('/content/cfps_decision_results_PRIVATE.zip')
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for path in sorted(PRIVATE.iterdir()):
        if path.is_file():z.write(path,arcname=path.name)
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None,'ZIP verification failed'
print('Saved private results:',archive,'bytes:',archive.stat().st_size)
files.download(str(archive))
''')]
N.metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python'}}
p=r/'notebooks/cfps_round1_decision_framework_colab.ipynb'
nbf.write(N,p)
print('notebook',p,'cells',len(N.cells),'bytes',p.stat().st_size)
