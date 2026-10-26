import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
import zarr
from ezst.atera import load_atera_zarr, plot_atera_segmentation, plot_atera_annotations


def official_matrix(tmp_path):
    # Unfiltered row, negative controls and reordered gene IDs exercise the
    # distinctions that would change biological counts or spatial identity.
    counts = np.array([[2, 7, 3], [1, 9, 0], [4, 2, 5], [6, 3, 8]], dtype=np.int32)
    obs = pd.DataFrame({
        'cell_id': [4294967300, 4294967301, 4294967302, 4294967303],
        'filtered': [True, False, True, True], 'centroid_column': [10., 20., 30., 40.],
        'centroid_row': [1., 2., 3., 4.], 'cell_area': [20.]*4, 'nucleus_area': [5.]*4,
        'transcript_counts': counts[:, [0, 2]].sum(1), 'total_counts': counts.sum(1),
        'control_probe_counts': counts[:, 1], 'control_codeword_counts': [0]*4,
        'azimuth_coarse': pd.Categorical(['A', 'B', 'A', 'B']),
    }, index=pd.Index(['c4', 'c5', 'c6', 'c7'], name='barcode'))
    var = pd.DataFrame({'feature_name': ['G2', 'Negative', 'G1'],
                        'feature_type': pd.Categorical(['Gene Expression', 'Negative Control Probe', 'Gene Expression']),
                        'feature_id': ['id2', 'ctrl', 'id1']}, index=['id2', 'ctrl', 'id1'])
    a = ad.AnnData(sparse.csr_matrix(counts), obs=obs, var=var)
    a.obsm['X_umap'] = np.arange(8).reshape(4, 2).astype(float)
    path = tmp_path/'cell_feature_matrix.zarr'
    a.write_zarr(path)
    (tmp_path/'experiment.spatial').write_text(json.dumps({'analysis_sw_version':'atera-test',
        'roi_uuid':'roi1', 'explorer_files':{'cell_feature_zarr_filepath':path.name}}))
    return path, a


def test_controls_filter_coordinates_and_counts_survive_conversion(tmp_path):
    source, original = official_matrix(tmp_path)
    target = tmp_path/'out.h5ad'
    summary = json.loads(load_atera_zarr(tmp_path, target, chunk_rows=2))
    a = ad.read_h5ad(target)
    np.testing.assert_array_equal(a.X.toarray(), [[2,3],[4,5],[6,8]])
    np.testing.assert_array_equal(a.layers['counts'].toarray(), a.X.toarray())
    np.testing.assert_array_equal(a.obsm['spatial'], [[10,1],[30,3],[40,4]])
    assert list(a.obs_names) == ['c4','c6','c7']
    assert list(a.var_names) == ['G2','G1']
    assert list(a.var.feature_id) == ['id2','id1']
    assert list(a.obs.control_probe_counts) == [7,2,3]
    assert a.uns['spatial_meta']['coordinate_unit'] == 'micrometers'
    np.testing.assert_array_equal(a.obsm['X_umap_10x'], original.obsm['X_umap'][[0,2,3]])
    a.X.data[:] = 0
    assert a.layers['counts'].sum() == 28
    assert summary['official_retained_cells'] == 3
    assert summary['gene_counts_sum'] == 29


def test_roi_is_physical_and_half_open_and_filter_can_be_disabled(tmp_path):
    source, _ = official_matrix(tmp_path)
    target = tmp_path/'roi.h5ad'
    load_atera_zarr(source, target, roi=(20, 0, 40, 5), use_official_filter=False, chunk_rows=1)
    a = ad.read_h5ad(target)
    assert list(a.obs_names) == ['c5', 'c6']
    np.testing.assert_array_equal(a.X.toarray(), [[1,0],[4,5]])
    assert list(a.obs.official_keep) == [False, True]


@pytest.mark.parametrize('label_key', ['reviewed_type', None])
def test_native_input_uses_actual_sample_and_annotation_field(tmp_path, label_key):
    from ezst.visium_new_platforms import load_atera
    source, original = official_matrix(tmp_path)
    original.obs.rename(columns={'azimuth_coarse': 'reviewed_type'}, inplace=True)
    original.write_zarr(source)
    target = tmp_path/'another_sample.h5ad'
    load_atera(str(tmp_path/'experiment.spatial'), sample_id='another_sample',
               output_path=str(target), cell_type_key=label_key)
    result = ad.read_h5ad(target)
    assert set(result.obs.sample_id) == {'another_sample'}
    assert result.uns['spatial_meta']['cell_type_source'] == (label_key or '')
    if label_key is None:
        assert 'cell_type' not in result.obs
    else:
        assert result.obs.cell_type.tolist() == ['A', 'A', 'B']
    np.testing.assert_array_equal(result.obsm['spatial'], [[10, 1], [30, 3], [40, 4]])


def test_h5ad_reload_preserves_physical_metadata_counts_and_selected_labels(tmp_path):
    from ezst.visium_new_platforms import load_atera
    source, _ = official_matrix(tmp_path)
    initial = tmp_path/'initial.h5ad'
    load_atera_zarr(source, initial)
    a = ad.read_h5ad(initial)
    a.obs['reviewed_type'] = pd.Categorical(['Type 1', 'Type 2', 'Type 2'])
    a.write_h5ad(initial)
    target = tmp_path/'reloaded.h5ad'
    load_atera(str(initial), sample_id='new_sample', output_path=str(target), cell_type_key='reviewed_type')
    result = ad.read_h5ad(target)
    assert result.uns['spatial_meta']['coordinate_unit'] == 'micrometers'
    assert result.uns['spatial_meta']['cell_type_source'] == 'reviewed_type'
    assert result.obs['cell_type'].tolist() == ['Type 1', 'Type 2', 'Type 2']
    assert set(result.obs.sample_id) == {'new_sample'}
    np.testing.assert_array_equal(result.obsm['spatial'], a.obsm['spatial'])
    assert (result.layers['counts'] != a.layers['counts']).nnz == 0


def test_bad_official_counts_do_not_replace_existing_output(tmp_path):
    source, _ = official_matrix(tmp_path)
    g = zarr.open_group(str(source), mode='a')
    g['obs/transcript_counts'][0] = 999
    target = tmp_path/'out.h5ad'
    target.write_bytes(b'existing user output')
    with pytest.raises(ValueError, match='disagree'):
        load_atera_zarr(source, target, chunk_rows=2)
    assert target.read_bytes() == b'existing user output'


def test_empty_roi_and_visualization_matrix_fail_clearly(tmp_path):
    source, _ = official_matrix(tmp_path)
    with pytest.raises(ValueError, match='no cells'):
        load_atera_zarr(source, tmp_path/'empty.h5ad', roi=(100,100,101,101))
    wrong = tmp_path/'csc_cell_feature_matrix.zarr.zip'
    wrong.touch()
    with pytest.raises(ValueError, match='includes cell metadata'):
        load_atera_zarr(wrong, tmp_path/'wrong.h5ad')


def test_segmentation_matches_ids_even_when_source_order_differs(tmp_path, monkeypatch):
    source, _ = official_matrix(tmp_path)
    target = tmp_path/'roi.h5ad'
    load_atera_zarr(source, target)
    cells = zarr.open_group(str(tmp_path/'cells.zarr'), mode='w', zarr_format=2)
    cells.create_array('cell_id', data=np.array([[7,1],[4,1],[6,1],[5,1]],dtype=np.uint32))
    summary = cells.create_array('cell_summary', data=np.array([[40,4],[10,1],[30,3],[20,2]],dtype=float))
    summary.attrs['column_names'] = ['cell_centroid_x','cell_centroid_y']
    cells.create_array('masks/homogeneous_transform', data=np.eye(4))
    cells.create_array('polygon_sets/1/cell_index', data=np.arange(4,dtype=np.uint32))
    cells.create_array('polygon_sets/1/num_vertices', data=np.full(4,4,dtype=np.int32))
    cells.create_array('polygon_sets/1/vertices', data=np.array([[x-1,y-1,x+1,y-1,x,y+1,x-1,y-1] for x,y in [[40,4],[10,1],[30,3],[20,2]]],dtype=np.float32))
    result = plot_atera_segmentation(target, tmp_path/'cells.zarr', tmp_path/'segmentation.png')
    assert result['matched_cells'] == 3
    assert result['max_centroid_error_um'] == 0
    # Capture the colors actually handed to the polygon renderer: segmentation
    # order is c7, c4, c6, while expression order is c4, c6, c7.
    from matplotlib.collections import PolyCollection
    rendered = []
    original_init = PolyCollection.__init__

    def capture(self, vertices, *args, **kwargs):
        rendered.append(list(kwargs['facecolors']))
        original_init(self, vertices, *args, **kwargs)

    monkeypatch.setattr(PolyCollection, '__init__', capture)
    output = tmp_path/'annotation.png'
    plot_atera_annotations(target, output, cells_path=tmp_path/'cells.zarr', source_matrix_path=source)
    mapping = pd.read_csv(output.with_suffix('.cells.csv'))
    assert mapping.barcode.tolist() == ['c4', 'c6', 'c7']
    assert mapping.official_label.tolist() == ['A', 'A', 'B']
    assert mapping.polygon_row.tolist() == [1, 2, 0]
    assert mapping.color.iloc[0] == mapping.color.iloc[1]
    assert mapping.color.iloc[0] != mapping.color.iloc[2]
    assert all(colors == mapping.color.iloc[[2, 0, 1]].tolist() for colors in rendered)
    assert len(rendered) == 2  # Main ROI and magnified detail use the same mapping.
    report = json.loads(output.with_suffix('.json').read_text())
    assert report['plotted_polygons'] == 3
    assert report['roi_label_counts'] == {'A': 2, 'B': 1}
    assert Path(report['full_tissue_figure']).is_file()
