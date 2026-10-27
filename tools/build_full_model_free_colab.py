#!/usr/bin/env python3
"""Compile the independently executable CFPS Colab notebook from mirrored modules.

The notebook embeds sources verbatim and never imports the mirrored analysis modules.
Source data are not bundled; live database credentials are never written to disk.
"""
from __future__ import annotations
import hashlib
from pathlib import Path
import subprocess
import sys
import textwrap
import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'notebooks/cfps_round1_full_model_free_colab.ipynb'
MODULE = ROOT / 'analysis/cfps_selection/model_free_pipeline.py'
ADAPTER = ROOT / 'analysis/cfps_selection/input_sources.py'


def main():
    subprocess.run([sys.executable, str(ROOT / 'tools/build_colab_notebook.py')], check=True, cwd=ROOT)
    nb = nbf.read(ROOT / 'notebooks/cfps_database_snapshot_and_exploration.ipynb', as_version=4)
    nb.cells[1].source = '''# Source selection: change ONE variable before running the notebook.
# "published_snapshot"   : upload clusters_designs.json; includes GP-informed annotations.
# "pocketbase_snapshot"  : upload a previously fetched pocketbase_snapshot.json.
# "pocketbase_live"      : fetch the current read-only snapshot using authorized credentials.
INPUT_SOURCE = "published_snapshot"
PUBLISHED_SNAPSHOT_PATH = "/content/clusters_designs.json"
POCKETBASE_SNAPSHOT_PATH = "/content/pocketbase_snapshot.json"
GP_LINKS_PATH = "/content/clusters.json.gz"  # optional; relevant to published snapshot only
EXTERNAL_PLATES_ZIP = "/content/cfps-plates.zip"  # optional
EXTERNAL_POLICY = "flag"  # "flag" or "exclude"; never silently discard an ID
WORKING_QUOTA_NOTEBOOK = 320  # unverified working assumption
GEOMETRY_OBJECTIVE = "facility"  # facility | maximin; both modes evaluated by default
AUTO_DOWNLOAD_PRIVATE_RESULTS = True

RUN_POCKETBASE_FETCH = INPUT_SOURCE == "pocketbase_live"
RUN_CORRELATION_CLUSTERING = False  # optional independent descriptive analysis
RUN_GP_BASELINE_AUDIT = INPUT_SOURCE == "published_snapshot"
RUN_PORTFOLIO_PREVIEWS = False
RUN_MODEL_FREE_PIPELINE = True
GP_DESIGNS_PATH = PUBLISHED_SNAPSHOT_PATH  # backward-compatible base-audit input
print("Input source:", INPUT_SOURCE)
print("GP-informed comparison available:", RUN_GP_BASELINE_AUDIT)
'''
    nb.cells[0].source += '''\n\n**Research attribution (methods only):** The GP-informed forecasts, shared-error clustering, virtual-LHS sensitivity analysis and original proposed runs were developed by **Clinton Monteiro**. The model-free recipe-space auditing, deduplication, exact-contrast analysis, facility-location and Maximin selection, and comparison visualization were developed by **Kourosh Afshinjoo**. The `htgaa-automation` platform, authentication and deployment infrastructure are maintained by **Ronan Chance**. Scientific methods are otherwise named by their mathematical purpose, not their contributors. See `docs/CONTRIBUTIONS.md`.\n\n**Run contract:** Choose `INPUT_SOURCE` in the configuration cell. Each run records the original input SHA-256 and a normalized snapshot SHA-256. PocketBase-only runs cannot reconstruct the GP-informed comparison without the independently supplied, version-matched annotated export. Private ID-level exports remain outside Git.\n'''
    nb.cells.append(nbf.v4.new_markdown_cell('''
## 18. Source normalization and provenance

Choose `INPUT_SOURCE` in the first configuration cell. `published_snapshot` uses the
GP-informed `clusters_designs.json` provided with the analysis; `pocketbase_snapshot`
uses a previously authorized `pocketbase_snapshot.json`; `pocketbase_live` first fetches
a read-only PocketBase snapshot using Colab Secrets. Both PocketBase modes normalize
`cfps_designs.design_json` into an analysis-compatible *experimental-recipe pool*.

**No forecasts or selection annotations are fabricated.** PocketBase-only runs enable
`model_free` selection; the historical GP-informed baseline and hybrid are unavailable
unless a validated GP-annotated export is analyzed separately. Non-experimental records
are excluded and counted, while malformed experimental records cause an error.
Inputs and normalized copies are private; nominally 20 µL volumes are checked by the
audit, not chemically inferred from stock names.
'''))
    adapter_source = ADAPTER.read_text(encoding='utf-8')
    nb.cells.append(nbf.v4.new_code_cell(adapter_source))
    nb.cells.append(nbf.v4.new_markdown_cell('''
## 19. Reproducible pre-experimental portfolio analysis

Methods: direct recipe identity, exact external overlap, controlled single-factor
contrasts and sweeps, normalized dispensing-volume/presence features, deterministic
facility-location and Maximin selection. Distances are *not* output measurements,
yield predictions, empirical response correlations, or statistical information gain.

Modes: `gp_model_baseline` reproduces the source's GP-informed roles; `hybrid` keeps
those roles while adding geometric coverage; `model_free` does not consult those roles
or GP forecasts. Only the third mode is available from an unannotated PocketBase pool.

The following implementation is **byte-for-byte identical** to the accompanying
`analysis/cfps_selection/model_free_pipeline.py`, but defined inside this notebook.
'''))
    mf_source = MODULE.read_text(encoding='utf-8')
    nb.cells.append(nbf.v4.new_code_cell(mf_source))
    nb.cells.append(nbf.v4.new_code_cell(textwrap.dedent('''
# 20. Resolve the explicitly selected input, normalize where necessary, and run.
from pathlib import Path
import json, hashlib
from IPython.display import display
import pandas as pd
if RUN_MODEL_FREE_PIPELINE:
    if INPUT_SOURCE == 'published_snapshot':
        raw_source_path = Path(PUBLISHED_SNAPSHOT_PATH).expanduser()
        normalized_kind = 'published_snapshot'
    elif INPUT_SOURCE == 'pocketbase_snapshot':
        raw_source_path = Path(POCKETBASE_SNAPSHOT_PATH).expanduser()
        normalized_kind = 'pocketbase_snapshot'
    elif INPUT_SOURCE == 'pocketbase_live':
        if not RUN_POCKETBASE_FETCH or 'RAW_SNAPSHOT_PATH' not in globals():
            raise RuntimeError('Live PocketBase fetch did not produce RAW_SNAPSHOT_PATH')
        raw_source_path = Path(RAW_SNAPSHOT_PATH)
        normalized_kind = 'pocketbase_snapshot'
    else:
        raise ValueError('INPUT_SOURCE must be published_snapshot, pocketbase_snapshot or pocketbase_live')
    if not raw_source_path.is_file():
        raise FileNotFoundError(f'Missing private source file: {raw_source_path}')

    data, source_provenance = normalize_analysis_source(raw_source_path, normalized_kind)
    private_folder = PRIVATE_ROOT / 'portfolio_analysis'
    private_folder.mkdir(parents=True, exist_ok=True)
    normalized_input = private_folder / 'normalized_candidate_pool_PRIVATE.json'
    normalized_input.write_text(json.dumps(data, sort_keys=True, indent=2) + '\\n', encoding='utf-8')
    normalized_input.chmod(0o600)
    print('Source summary:', source_provenance)
    if normalized_kind == 'published_snapshot':
        analysis_input = raw_source_path
        analysis_modes = ('gp_model_baseline', 'hybrid', 'model_free')
    else:
        analysis_input = normalized_input
        analysis_modes = ('model_free',)
    ext = Path(EXTERNAL_PLATES_ZIP).expanduser()
    ext = ext if ext.is_file() else None
    if ext is None:
        print('External plates not supplied: external exact-overlap audit disabled.')

    run = model_free_run(analysis_input, ext, quota=WORKING_QUOTA_NOTEBOOK,
                         external_policy=EXTERNAL_POLICY, modes=analysis_modes)
    original_paths = [raw_source_path] + ([ext] if ext else [])
    source_manifest = {
        'input_source': INPUT_SOURCE,
        'adapter_version': SOURCE_ADAPTER_VERSION,
        'original_input_sha256': input_sha256(raw_source_path),
        'normalized_input_sha256': input_sha256(normalized_input),
        'input_counts': source_provenance['counts_by_protein'],
        'gp_baseline_available': source_provenance['gp_baseline_available'],
        'external_input_sha256': input_sha256(ext) if ext else None,
        'working_quota_not_verified': WORKING_QUOTA_NOTEBOOK,
        'external_policy': EXTERNAL_POLICY,
        'algorithm_version': MODEL_FREE_VERSION,
    }
    (private_folder / 'source_manifest.json').write_text(json.dumps(source_manifest, indent=2, sort_keys=True) + '\\n')
    private_payload = write_private_analysis_payload(
        run['results'], private_folder / 'model_free_analysis.json',
        input_paths=original_paths, external_summary=run['external_summary'],
        provenance=source_manifest)
    table = pd.DataFrame([{
       'protein': t, 'strategy': key, 'selected': entry['n_selected'],
       'mean_nearest_distance': round(entry['distance']['mean'], 6),
       'factor_bin_coverage': round(entry['factor_bin_coverage_mean'], 4) if entry['factor_bin_coverage_mean'] is not None else None,
       'exact_contrast_pairs': entry['exact_selected_contrasts']['unique_exact_pairs'],
       'external_exact_matches': entry['external_overlaps_selected']
    } for t, strategies in run['results'].items() for key, entry in strategies.items()])
    display(table)
    table.to_csv(private_folder / 'strategy_comparison.csv', index=False)
    (private_folder / 'selected_ids_PRIVATE.json').write_text(json.dumps({
       t: {strategy: portfolio['selected_ids'] for strategy, portfolio in group.items()}
       for t, group in run['results'].items()}, indent=2, sort_keys=True) + '\\n')
    print('Protected visualization JSON:', private_payload['path'])
''')))
    nb.cells.append(nbf.v4.new_code_cell(textwrap.dedent('''
# 21. Export the reproducibility bundle for comparison with another run.
# Contains participant IDs, external-overlap flags and potentially sensitive
# predictions: do not commit or post the archive publicly.
if RUN_MODEL_FREE_PIPELINE:
    from zipfile import ZipFile, ZIP_DEFLATED
    from datetime import datetime, timezone
    export_path = PRIVATE_ROOT / ('cfps_private_run_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '.zip')
    export_files = sorted(p for p in private_folder.rglob('*') if p.is_file())
    with ZipFile(export_path, 'w', ZIP_DEFLATED) as z:
        for item in export_files:
            z.write(item, arcname=item.relative_to(private_folder))
    with ZipFile(export_path) as z:
        assert z.testzip() is None, 'Corrupt results ZIP'
    print('Download this PRIVATE results archive for cross-run comparison:', export_path)
    print('ZIP SHA256:', input_sha256(export_path))
    if AUTO_DOWNLOAD_PRIVATE_RESULTS:
        try:
            from google.colab import files as colab_files
            colab_files.download(str(export_path))
        except ImportError:
            print('Outside Colab: the private ZIP remains at the path printed above.')
''')))
    nb.cells.append(nbf.v4.new_markdown_cell('''
## 22. Deployment and evidence boundaries

- The source distribution contains no private participant records or tokens.
- Transfer `model_free_analysis.json` exclusively to a Keycloak-gated, Git-ignored
  deployment directory. Do **not** place it under `static/` or commit it.
- The optional `clusters.json.gz` is a *published-link diagnostic*, not an empirical
  response covariance matrix; its absence does not affect model-free selection.
- PocketBase requires valid authorized credentials in the live mode and may expose a
  different, later participant snapshot than the published GP-informed export.
- Run identifiers, exact input checksums and population counts are saved for comparison.
- All selection results are provisional until final quota, controls and eligibility
  constraints are independently verified.
'''))
    # Remove contributor names outside explicit attribution. Methodological narrative
    # and public labels use neutral scientific terminology.
    # Narrative cells use neutral method names. Contributor names occur only in
    # the explicit research-attribution section above, not in API identifiers.
    # Preserve embedded module *exactly*: transform only narrative/code cells other
    # than the three byte-for-byte source mirrors.
    exact_sources = {
        (ROOT/'analysis/cfps_correlation/pairwise_correlation_clustering.py').read_text(encoding='utf-8'),
        (ROOT/'analysis/cfps_selection/selection_framework.py').read_text(encoding='utf-8'),
        adapter_source, mf_source,
    }
    for i, cell in enumerate(nb.cells):
        if cell.source not in exact_sources and i != 0:
            cell.source = cell.source.replace('**the published GP-informed', '**GP-informed')
        cell['id'] = hashlib.sha1((str(i)+'|'+cell.source).encode('utf-8')).hexdigest()[:12]
        if cell.cell_type == 'code':
            cell['execution_count'] = None
            cell['outputs'] = []
    nb.metadata['analysis_provenance'] = {
        'generator': 'tools/build_full_model_free_colab.py',
        'input_adapter_sha256': hashlib.sha256(adapter_source.encode()).hexdigest(),
        'source_correlation_sha256': hashlib.sha256((ROOT/'analysis/cfps_correlation/pairwise_correlation_clustering.py').read_bytes()).hexdigest(),
        'source_baseline_audit_sha256': hashlib.sha256((ROOT/'analysis/cfps_selection/selection_framework.py').read_bytes()).hexdigest(),
        'source_model_free_sha256': hashlib.sha256(mf_source.encode()).hexdigest(),
        'module_source_embedded_verbatim': True,
        'runtime_imports_mirrored_python_modules': False,
    }
    nb.metadata.setdefault('colab', {})['name'] = DEST.name
    nbf.validate(nb)
    nbf.write(nb, DEST)
    print('Built', DEST, 'cells', len(nb.cells))


if __name__ == '__main__':
    main()
