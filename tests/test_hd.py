import gzip
import pickle
import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy import sparse
from scipy.io import mmwrite
from ezst._geometry import (filter_nuclei_instances, load_nuclei_polygons,
                            save_nuclei_polygons, read_proseg_counts, nuclei_xy_polygons)
from ezst.visium_new_platforms import bin_to_cell


def polygons():
    # StarDist stores row (y), col (x); nuclei centers are (1, 1), (3, 1).
    return {'coord': np.array([[[0, 0, 2, 2], [0, 2, 2, 0]],
                               [[0, 0, 2, 2], [2, 4, 4, 2]]], dtype=float),
            'prob': np.array([.7, .9]), 'points': np.array([[1, 1], [1, 3]])}


def test_touching_instances_preserve_identity():
    labels = np.array([[1, 1, 2, 2], [1, 1, 2, 2]])
    result, poly = filter_nuclei_instances(labels, polygons(), min_area=4)
    np.testing.assert_array_equal(result, labels)
    np.testing.assert_array_equal(poly['label_ids'], [1, 2])


def test_small_instance_filter_remaps_geometry_and_labels(tmp_path):
    labels = np.array([[1, 0, 2, 2], [0, 0, 2, 2]])
    result, poly = filter_nuclei_instances(labels, polygons(), min_area=2)
    np.testing.assert_array_equal(result, [[0, 0, 1, 1], [0, 0, 1, 1]])
    np.testing.assert_array_equal(poly['coord'], polygons()['coord'][1:])
    np.testing.assert_array_equal(poly['prob'], [.9])
    path = tmp_path / 'nuclei.npz'
    save_nuclei_polygons(path, poly)
    np.testing.assert_array_equal(load_nuclei_polygons(path)['coord'], poly['coord'])
    np.testing.assert_array_equal(nuclei_xy_polygons(path, object_ids=[1])[0], [[2,0],[4,0],[4,2],[2,2]])


def test_legacy_geometry_requires_explicit_opt_in(tmp_path):
    path = tmp_path / 'nuclei.pckl'
    path.write_bytes(pickle.dumps(polygons()))
    with pytest.raises(ValueError, match='allow_legacy_pickle'):
        load_nuclei_polygons(path)
    np.testing.assert_array_equal(load_nuclei_polygons(path, allow_legacy_pickle=True)['coord'], polygons()['coord'])


def test_windows_model_cache_loading_preserves_checksum(tmp_path, monkeypatch):
    import sys
    from types import ModuleType
    from ezst.visium_new_platforms import _load_stardist_pretrained
    directory = tmp_path / 'model_extracted'
    directory.mkdir()
    (directory / 'config.json').write_text('{}')
    class Model:
        @staticmethod
        def from_pretrained(name):
            error = OSError('symlink privilege')
            error.winerror = 1314
            raise error
        def __init__(self, config, name, basedir):
            assert basedir / name == directory
            self.basedir = basedir
    bridge = ModuleType('csbdeep.models.pretrained')
    bridge.get_model_details = lambda cls, name: ('model', None, {'url':'https://example.invalid/model.zip', 'hash':'verified_hash'})
    def get_file(**kwargs):
        assert kwargs['file_hash'] == 'verified_hash' and kwargs['extract']
        return str(directory)
    bridge.get_file = get_file
    monkeypatch.setitem(sys.modules, 'csbdeep.models.pretrained', bridge)
    assert _load_stardist_pretrained(Model, 'model').basedir is None


@pytest.mark.parametrize('dense', [False, True])
def test_overlap_unique_assignment_and_count_conservation(tmp_path, dense):
    counts = np.array([[10, 2], [20, 3], [30, 4], [99, 99]])
    a = ad.AnnData(counts if dense else sparse.csr_matrix(counts))
    a.var_names = ['GeneA', 'GeneB']
    a.obs_names = ['tie', 'right', 'left', 'outside']
    a.obsm['spatial'] = np.array([[2, 1], [2.8, 1], [1, 1], [100, 100]])
    source, geometry, output = [tmp_path / name for name in ('bins.h5ad', 'poly.npz', 'cells.h5ad')]
    a.write_h5ad(source)
    save_nuclei_polygons(geometry, polygons())
    mapping = tmp_path / 'mapping.csv'
    bin_to_cell(str(source), str(geometry), output_path=str(output), expand_px=1, barcode_map_path=str(mapping))
    cells = ad.read_h5ad(output)
    actual = cells.X.toarray() if sparse.issparse(cells.X) else cells.X
    np.testing.assert_array_equal(actual, [[40, 6], [20, 3]])
    np.testing.assert_array_equal(actual.sum(axis=0), counts[:3].sum(axis=0))
    np.testing.assert_array_equal(cells.obsm['spatial'], [[1.5, 1], [2.8, 1]])
    assert list(cells.var_names) == ['GeneA', 'GeneB']
    frame = pd.read_csv(mapping)
    assert frame['index'].is_unique and len(frame) == 3
    assert cells.uns['bin_to_cell']['n_bins_unassigned'] == 1


def proseg_files(tmp_path):
    matrix = tmp_path/'counts.mtx.gz'
    with gzip.open(matrix, 'wb') as handle:
        mmwrite(handle, sparse.csr_matrix([[2, 3], [4, 5]]))
    cells, genes = tmp_path/'cells.csv.gz', tmp_path/'genes.csv.gz'
    pd.DataFrame({'cell':[1,0], 'centroid_x':[8.,2.], 'centroid_y':[10.,4.]}).to_csv(cells, index=False)
    pd.DataFrame({'gene':['Gapdh','Actb'], 'total_count':[6,8]}).to_csv(genes, index=False)
    return matrix, cells, genes


def test_proseg_metadata_restores_names_row_order_and_scale(tmp_path):
    a = read_proseg_counts(*proseg_files(tmp_path), microns_per_pixel=.5)
    assert list(a.var_names) == ['Gapdh', 'Actb']
    assert list(a.obs_names) == ['cell_0', 'cell_1']
    np.testing.assert_array_equal(a.X.toarray(), [[2,3],[4,5]])
    np.testing.assert_array_equal(a.obsm['spatial'], [[4,8],[16,20]])
    assert a.uns['cellpose_scale'] == .5


@pytest.mark.parametrize('fault', ['dimensions', 'duplicate_gene', 'cell_ids'])
def test_proseg_metadata_rejects_misalignment(tmp_path, fault):
    matrix, cells, genes = proseg_files(tmp_path)
    if fault == 'dimensions':
        pd.DataFrame({'gene':['Gapdh']}).to_csv(genes, index=False)
    elif fault == 'duplicate_gene':
        pd.DataFrame({'gene':['Gapdh','Gapdh']}).to_csv(genes, index=False)
    else:
        data = pd.read_csv(cells); data['cell'] = [0,0]; data.to_csv(cells, index=False)
    with pytest.raises(ValueError):
        read_proseg_counts(matrix, cells, genes, microns_per_pixel=.5)
