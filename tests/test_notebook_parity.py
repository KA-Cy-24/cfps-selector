"""Notebook self-containment and source parity regression (no private inputs)."""
from pathlib import Path
import nbformat
ROOT=Path(__file__).resolve().parents[1]

def test_exact_embedded_module_source_parity():
    nb=nbformat.read(ROOT/'notebooks/cfps_round1_full_model_free_colab.ipynb',as_version=4)
    source_cells=[c.source for c in nb.cells if c.cell_type=='code']
    for rel in ['analysis/cfps_correlation/pairwise_correlation_clustering.py',
                'analysis/cfps_selection/input_sources.py',
                'analysis/cfps_selection/selection_framework.py',
                'analysis/cfps_selection/model_free_pipeline.py']:
        assert (ROOT/rel).read_text(encoding='utf-8') in source_cells, rel
    assert nb.metadata['analysis_provenance']['runtime_imports_mirrored_python_modules'] is False
    assert 'RUN_MODEL_FREE_PIPELINE = True' in nb.cells[1].source


def test_all_notebook_code_cells_compile():
    nb=nbformat.read(ROOT/'notebooks/cfps_round1_full_model_free_colab.ipynb',as_version=4)
    for index,cell in enumerate(nb.cells):
        if cell.cell_type=='code':
            compile(cell.source,f'notebook cell {index}','exec')


def test_source_switch_and_private_download_cells_present():
    nb=nbformat.read(ROOT/'notebooks/cfps_round1_full_model_free_colab.ipynb',as_version=4)
    combined='\n'.join(c.source for c in nb.cells)
    assert 'INPUT_SOURCE = "published_snapshot"' in nb.cells[1].source
    for source in ('pocketbase_snapshot','pocketbase_live','published_snapshot'):
        assert source in nb.cells[1].source
    assert 'colab_files.download(str(export_path))' in combined
    assert 'normalize_analysis_source(raw_source_path, normalized_kind)' in combined
