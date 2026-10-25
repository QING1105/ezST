import json
import sys
from pathlib import Path
from types import SimpleNamespace
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
from ezst import visium
from ezst import run_hd_downstream
from ezst.visium_new_platforms import run_proseg


def test_downstream_keeps_annotations_domains_and_custom_key(tmp_path, monkeypatch):
    a = ad.AnnData(sparse.csr_matrix(np.ones((20, 3))))
    a.obsm['spatial'] = np.arange(40).reshape(20, 2)
    source = tmp_path/'cells.h5ad'
    a.write_h5ad(source)
    out = tmp_path/'results'
    seen = []
    def annotate(command, **kwargs):
        result = ad.read_h5ad(source)
        result.obs['cell_type'] = pd.Categorical(['A']*10 + ['B']*10)
        result.obs['cell_type_prob'] = 1.0
        result.write_h5ad(out/'annotated_cells.h5ad')
    def domains(path, output_path, **kwargs):
        result = ad.read_h5ad(path)
        assert 'custom_type' in result.obs and 'cell_type_prob' in result.obs
        result.obs['spatial_domain'] = pd.Categorical(['0']*10 + ['1']*10)
        result.write_h5ad(output_path)
        seen.append('domains')
        return 'domains complete'
    def downstream(path, cell_type_key, **kwargs):
        result = ad.read_h5ad(path)
        assert cell_type_key == 'custom_type'
        assert 'spatial_domain' in result.obs and 'custom_type' in result.obs
        seen.append('downstream')
        return 'complete'
    monkeypatch.setattr(run_hd_downstream.subprocess, 'run', annotate)
    monkeypatch.setattr(visium, 'identify_spatial_domains', domains)
    monkeypatch.setattr(visium, 'spatial_neighborhood_enrichment', downstream)
    monkeypatch.setattr(visium, 'infer_spatial_cell_communication', downstream)
    monkeypatch.setattr(sys, 'argv', ['run_hd_downstream', '--cells', str(source), '--ref', 'reference.h5ad',
                                    '--out-dir', str(out), '--cell-type-key', 'custom_type'])
    run_hd_downstream.main()
    assert seen == ['domains', 'downstream', 'downstream']


def test_proseg_command_exports_gene_metadata_and_imports_by_name(tmp_path, monkeypatch):
    import subprocess
    import gzip
    directory = tmp_path/'sample/binned_outputs/square_002um'
    directory.mkdir(parents=True)
    with h5py.File(directory/'filtered_feature_bc_matrix.h5', 'w') as handle:
        handle.create_dataset('matrix/data', data=np.array([2, 3, 4, 5]))
    labels_path = tmp_path/'labels.npz'
    np.savez(labels_path, labels=np.array([[1,1,0],[0,2,2]], dtype=np.uint32))
    output = tmp_path/'proseg'
    def export(command, **kwargs):
        assert command[0] == 'chosen-proseg'
        assert command[command.index('--recorded-samples') + 1] == '1'
        gene_path = command[command.index('--output-gene-metadata') + 1]
        count_path = command[command.index('--output-counts') + 1]
        cell_path = command[command.index('--output-cell-metadata') + 1]
        with gzip.open(count_path, 'wb') as handle:
            mmwrite(handle, sparse.csr_matrix([[2,3],[4,5]]))
        pd.DataFrame({'gene':['Gapdh','Actb']}).to_csv(gene_path, index=False)
        pd.DataFrame({'cell':[0,1], 'centroid_x':[2.,8.], 'centroid_y':[4.,10.]}).to_csv(cell_path, index=False)
        return SimpleNamespace(returncode=0, stdout='completed fixture', stderr='')
    monkeypatch.setattr(subprocess, 'run', export)
    run_proseg(str(tmp_path/'sample'), str(labels_path), str(output), cellpose_scale=.5,
               samples=1, burnin_samples=1, proseg_executable='chosen-proseg')
    result = ad.read_h5ad(output/'proseg_cells.h5ad')
    assert list(result.var_names) == ['Gapdh', 'Actb']
    np.testing.assert_array_equal(result.obsm['spatial'], [[4,8],[16,20]])
