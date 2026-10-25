import gzip, sys, types
from unittest.mock import patch
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmwrite
import pytest
import ezst.visium as v

@pytest.fixture(autouse=True)
def classic_data(tmp_path):
    global base, raw, dest
    dest = tmp_path
    rng = np.random.default_rng(42)
    counts = rng.negative_binomial(3, rng.uniform(.15,.9,600), size=(120,600)).astype(np.float32)
    counts[:60,20:45] += 30
    counts[60:,45:70] += 30
    base = sc.AnnData(sp.csr_matrix(counts))
    base.var_names = ['mt-Co1'] + [f'Gene{i}' for i in range(1,600)]
    base.obs_names = [f'spot{i}' for i in range(120)]
    base.obsm['spatial'] = np.column_stack([np.arange(120)%12,np.arange(120)//12]).astype(float)
    base.obs['cell_type'] = pd.Categorical(['A']*60 + ['B']*60)
    raw = dest/'raw.h5ad'
    base.write_h5ad(raw)

def test_gz10x():
    d=dest/'tenx';d.mkdir(exist_ok=True)
    with gzip.open(d/'matrix.mtx.gz','wb') as f: mmwrite(f,base.X.T)
    pd.DataFrame({'id':base.var_names,'symbol':base.var_names,'type':'Gene Expression'}).to_csv(d/'features.tsv.gz',sep='\t',header=False,index=False)
    pd.Series(base.obs_names).to_csv(d/'barcodes.tsv.gz',header=False,index=False)
    out=dest/'loaded.h5ad';v.load_visium_data(str(d),output_path=str(out))
    a=sc.read_h5ad(out); assert a.shape==base.shape and (a.X!=base.X).nnz==0
def test_mouseqc():
    a=base.copy();x=a.X.tolil();x[0,0]=100000;a.X=x.tocsr(); p=dest/'mouse.h5ad';a.write_h5ad(p)
    out=dest/'qc.h5ad';v.filter_visium_spots(str(p),str(out),min_counts=0,min_genes=0,pct_mt=20)
    b=sc.read_h5ad(out);assert 'spot0' not in b.obs_names,'High-mito mouse spot survived QC'
def test_hvg():
    out=dest/'normalized.h5ad';log=v.normalize_visium(str(raw),str(out),n_top_genes=100)
    actual=sc.read_h5ad(out);expected=base.copy();sc.pp.highly_variable_genes(expected,n_top_genes=100,flavor='seurat_v3')
    assert np.array_equal(actual.var.highly_variable,expected.var.highly_variable),'Selected genes differ from seurat_v3 despite log claim'
    assert (actual.raw.X!=base.X).nnz==0
    assert 'seurat_v3' in log
def test_cluster_default():
    v.normalize_visium(str(raw),str(dest/'normalized.h5ad'),n_top_genes=100)
    out=dest/'clustered.h5ad';v.cluster_spatial_data(str(dest/'normalized.h5ad'),str(out),n_pcs=10)
    a=sc.read_h5ad(out);assert 'leiden' in a.obs
def test_dense_domains():
    a=base[:,:50].copy();a.X=a.X.toarray();p=dest/'dense.h5ad';a.write_h5ad(p)
    out=dest/'domains.h5ad';v.identify_spatial_domains(str(p),str(out),resolution=.5,n_pcs=10)
    assert 'spatial_domain' in sc.read_h5ad(out).obs
def test_svg_after_domain():
    a=base.copy();sc.pp.normalize_total(a);sc.pp.log1p(a);a.obs['spatial_domain']=a.obs.cell_type
    p=dest/'svg_input.h5ad';a.write_h5ad(p)
    out=dest/'svg.csv';v.find_spatially_variable_genes(str(p),str(out),n_top_genes=5)
    df=pd.read_csv(out,index_col=0);assert len(df)==5
def test_svg_all_genes():
    import squidpy as sq
    seen={}
    def spy(a,**kw):
        seen['genes']=list(kw['genes'])
        a.uns['moranI']=pd.DataFrame({'I':np.linspace(1,0,len(kw['genes']))},index=kw['genes'])
    with patch.object(sq.gr,'spatial_autocorr',spy):v.find_spatially_variable_genes(str(raw),str(dest/'allgenes.csv'),n_top_genes=5)
    assert len(seen['genes'])==base.n_vars,f"Only {len(seen['genes'])}/{base.n_vars} genes evaluated"
def test_destvi_counts():
    v.normalize_visium(str(raw),str(dest/'normalized.h5ad'),n_top_genes=100)
    normalized=sc.read_h5ad(dest/'normalized.h5ad');seen={}
    class Cond:
        @staticmethod
        def setup_anndata(a,**kw):seen['ref']=a.X.copy()
        def __init__(self,*a,**kw):pass
        def train(self,**kw):pass
    class Dest:
        @staticmethod
        def setup_anndata(a,**kw):seen['spatial']=a.layers['counts'].copy()
        @staticmethod
        def from_rna_model(a,b):return Dest()
        def train(self,**kw):pass
        def get_proportions(self):return pd.DataFrame(np.full((base.n_obs,2),.5),columns=['A','B'])
    fake=types.ModuleType('scvi.model');fake.CondSCVI=Cond;fake.DestVI=Dest
    with patch.dict(sys.modules,{'scvi.model':fake}):
        v.deconvolve_spatial_destvi(str(dest/'normalized.h5ad'),str(dest/'normalized.h5ad'),'cell_type',str(dest/'destvi'),max_epochs=1,destvi_max_epochs=1)
    for label in ['ref','spatial']:
        assert np.allclose(seen[label].toarray(),base.X.toarray()),f'{label} model input differs from original counts'
