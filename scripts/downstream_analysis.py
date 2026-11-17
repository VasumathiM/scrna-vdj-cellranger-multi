import scanpy as sc
import scirpy as ir

base = "T1k_run/outs/per_sample_outs/T1k_run"
adata = sc.read_10x_mtx(f"{base}/count/sample_filtered_feature_bc_matrix")
adata_tcr = ir.io.read_10x_vdj(f"{base}/vdj_t/filtered_contig_annotations.csv")

# Gene expression
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs["pct_counts_mt"] < 10].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# TCR: merge by barcode, define clonotypes
ir.pp.merge_airr(adata, adata_tcr)      # older scirpy: ir.pp.merge_with_ir
ir.pp.ir_dist(adata)
ir.tl.define_clonotypes(adata)
ir.tl.clonotype_network(adata)
ir.pl.clonal_expansion(adata)
