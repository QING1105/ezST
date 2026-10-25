import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy import sparse
from ezst import visium as v
from ezst.annotate_cells import _to_counts


@pytest.fixture
def sample(tmp_path):
    rng = np.random.default_rng(121)
    a = ad.AnnData(sparse.csr_matrix(rng.poisson(4, (60, 20)).astype(float)))
    a.var_names = ['Actb'] + [f'Gene{i}' for i in range(19)]
    a.obs['kind'] = pd.Categorical(['A'] * 30 + ['B'] * 30)
    a.obsm['spatial'] = rng.uniform(0, 200, (60, 2))
    path = tmp_path / 'input.h5ad'
    a.write_h5ad(path)
    return a, path


def test_domain_graph_changes_only_when_spatial_weight_enabled(tmp_path, sample):
    a, path = sample
    graphs = {}
    for weight in [0, .5]:
        for permuted in [False, True]:
            current = a.copy()
            if permuted:
                current.obsm['spatial'] = current.obsm['spatial'][::-1].copy()
            current.write_h5ad(path)
            output = tmp_path / f'domains-{weight}-{permuted}.h5ad'
            v.identify_spatial_domains(str(path), str(output), spatial_weight=weight, n_pcs=5, resolution=.3)
            result = ad.read_h5ad(output)
            key = result.uns['domain_neighbors']['connectivities_key']
            graphs[weight, permuted] = result.obsp[key].toarray()
            assert result.uns['domain_clustering']['spatial_weight'] == weight
    np.testing.assert_allclose(graphs[0,False], graphs[0,True])
    assert not np.allclose(graphs[.5,False], graphs[.5,True])


@pytest.mark.parametrize('source', ['counts', 'raw'])
def test_count_models_recover_preserved_counts(sample, source):
    a, _ = sample
    counts = a.X.copy()
    if source == 'counts':
        a.layers['counts'] = counts
    else:
        a.raw = a.copy()
    a.X = a.X.astype(float)
    a.X.data = np.log1p(a.X.data * .37)
    a.uns['log1p'] = {'base':None}
    np.testing.assert_array_equal(_to_counts(a).X.toarray(), counts.toarray())


def test_count_model_rejects_normalized_values_without_counts(sample):
    a, _ = sample
    a.X.data = np.log1p(a.X.data * .37)
    a.uns['log1p'] = {'base':None}
    with pytest.raises(ValueError, match='raw counts'):
        _to_counts(a)


def test_spotlight_prefers_counts_layer_over_normalized_raw(sample, tmp_path, monkeypatch):
    from scipy.io import mmread
    a, path = sample
    a.layers['counts'] = a.X.copy()
    a.X.data = np.log1p(a.X.data * .37)
    a.raw = a.copy()
    a.uns['log1p'] = {'base':None}
    a.write_h5ad(path)
    captured = []
    def inspect_r(executable, script_path, output_dir):
        for prefix in ['ref', 'st']:
            data = mmread(Path(output_dir)/'_bridge'/f'{prefix}_matrix.mtx').T.toarray()
            np.testing.assert_array_equal(data, a.layers['counts'].toarray())
            captured.append(prefix)
        return SimpleNamespace(returncode=1, stdout='', stderr='Bridge input inspection complete')
    monkeypatch.setattr(v, '_run_r_script', inspect_r)
    v.deconvolve_spatial_spotlight(str(path), str(path), 'kind', str(tmp_path/'spotlight'))
    assert captured == ['ref','st']


@pytest.mark.parametrize('mode', ['expression', 'spatial'])
def test_cellchat_bridge_supplies_real_mode_and_raw_counts(tmp_path, sample, monkeypatch, mode):
    from scipy.io import mmread
    a, path = sample
    a.layers['counts'] = a.X.copy()
    a.X.data = np.log1p(a.X.data * .37)
    a.uns['log1p'] = {'base':None}
    a.write_h5ad(path)
    out = tmp_path / "CellChat user's data"
    seen = {}
    def inspect_r(executable, script_path, output_dir):
        script = Path(script_path).read_text(encoding='utf-8')
        seen['script'] = script
        pd.DataFrame({'source':['A'], 'target':['B']}).to_csv(Path(output_dir)/'cellchat_communication.csv', index=False)
        return SimpleNamespace(returncode=0, stdout='', stderr='')
    monkeypatch.setattr(v, '_run_r_script', inspect_r)
    kwargs = dict(microns_per_coordinate=.5, spatial_tolerance=5, contact_range=20) if mode == 'spatial' else {}
    v.infer_spatial_cell_communication(str(path), 'kind', str(out), species='mouse', mode=mode, **kwargs)
    np.testing.assert_allclose(mmread(out/'_bridge/st_matrix.mtx').T.toarray(), a.layers['counts'].toarray())
    assert ('datatype = "spatial"' in seen['script']) == (mode == 'spatial')
    if mode == 'spatial':
        assert 'ratio = 0.5' in seen['script'] and 'tol = 5.0' in seen['script']
        assert 'distance.use = TRUE' in seen['script'] and 'contact.range = 20.0' in seen['script']
    assert json.loads((out/'run_settings.json').read_text())['mode'] == mode


@pytest.mark.parametrize('kwargs', [{}, {'microns_per_coordinate':-1, 'spatial_tolerance':5},
                                   {'microns_per_coordinate':1, 'spatial_tolerance':-1}])
def test_spatial_cellchat_requires_valid_physical_scale(sample, tmp_path, kwargs):
    _, path = sample
    with pytest.raises(ValueError):
        v.infer_spatial_cell_communication(str(path), 'kind', str(tmp_path/'bad'), species='mouse', mode='spatial', **kwargs)


def test_plugin_runtime_matches_package():
    root = Path(__file__).resolve().parents[1]
    bundle = root/'plugins/spatial-transcriptomics/skills/spatial-visium-hd'
    for source in (root/'src/ezst').glob('*.py'):
        assert (bundle/source.name).read_bytes() == source.read_bytes(), source.name
