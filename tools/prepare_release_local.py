#!/usr/bin/env python3
"""Create a source-only, checksum-verified CFPS publication release archive.

The archive contains no private input snapshots, submission records or result JSON.
Its manifest hashes the *source* files; manifests are excluded from their own hash.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import zipfile

ROOT = Path(__file__).resolve().parents[1]
INCLUDE_EXTENSIONS={'.py','.ipynb','.md','.js','.svelte','.css','.html','.txt'}
EXCLUDED_DIRS={'data','__pycache__','.pytest_cache','.git','.ipynb_checkpoints','node_modules'}
MANIFEST_NAME='SOURCE_RELEASE_MANIFEST.json'


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def source_files():
    selected=[]
    for p in ROOT.rglob('*'):
        if not p.is_file():
            continue
        rel=p.relative_to(ROOT)
        if any(k in EXCLUDED_DIRS for k in rel.parts):
            continue
        if p.name in {MANIFEST_NAME,'RELEASE_MANIFEST.json'}:
            continue
        if p.name == '.gitignore' or p.suffix.lower() in INCLUDE_EXTENSIONS:
            selected.append(p)
    return sorted(selected,key=lambda p:str(p.relative_to(ROOT)))


def create_release(destination: Path):
    selected=source_files()
    source_hashes={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()) for p in selected}
    manifest={
        'schema':'cfps-selection-source-release/v1',
        'version':'1.1.0',
        'source_only':True,
        'private_participant_data_included':False,
        'canonical_notebook':'notebooks/cfps_round1_full_model_free_colab.ipynb',
        'contribution_statement':'docs/CONTRIBUTIONS.md',
        'input_source_contract':'docs/INPUT_DATA_AND_REPRODUCIBILITY.md',
        'source_files_sha256':source_hashes,
        'python_build_version':platform.python_version(),
        'scientific_status':'PRE-EXPERIMENTAL; unvalidated working quota; no measured outcomes',
        'deployment_status':'SOURCE ONLY; live PocketBase/Keycloak integration not authenticated here',
    }
    manifest_bytes=(json.dumps(manifest,indent=2,sort_keys=True)+'\n').encode('utf-8')
    destination.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(destination,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for p in selected:
            z.write(p,arcname='cfps_selection_research/'+p.relative_to(ROOT).as_posix())
        z.writestr('cfps_selection_research/'+MANIFEST_NAME,manifest_bytes)
    with zipfile.ZipFile(destination) as z:
        assert z.testzip() is None
        for path,digest in source_hashes.items():
            assert sha256(z.read('cfps_selection_research/'+path)) == digest, path
    (ROOT/MANIFEST_NAME).write_bytes(manifest_bytes)
    print(json.dumps({'archive':str(destination),'files':len(selected),'sha256':sha256(destination.read_bytes()),
                      'source_integrity':'passed','private_data':'excluded'},indent=2))
    return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path('/mnt/data/HTGAA_CFPS_Publication_Ready_Framework.zip'))
    args=p.parse_args()
    create_release(args.output)
