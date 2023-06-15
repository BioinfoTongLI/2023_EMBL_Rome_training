#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import fire
import pandas as pd
import scanpy as sc


def compute_embeddings(adata):
    sc.pp.neighbors(adata, n_neighbors=300, n_pcs=10)
    sc.tl.umap(adata, random_state = 777)
    sc.tl.leiden(adata, resolution = 0.25, random_state = 777)
    sc.pl.umap(adata, color=['total_counts', "leiden"])
    sc.tl.pca(adata, svd_solver='arpack', n_comps=30,
              use_highly_variable=False)

def to_h5ad(countTable, centroids, stem, n_gene_min=4):
    adata = pd.read_csv(countTable, index_col=0)
    adata = adata.loc[:, ~adata.columns.isin(['background', 'infeasible'])]
    coord = pd.read_csv(centroids, index_col=0)
    adata = sc.AnnData(adata)
    adata.obsm['spatial'] = coord[['x', 'y']].values
    adata.obs['sample'] = stem

    # Remove cells with no mRNA
    adata.obs['total_counts'] = adata.X.sum(1)
    adata.obs['n_genes_by_counts'] = (adata.X > 0).sum(1)
    adata = adata[adata.obs.n_genes_by_counts > 0, :]

    adata.var['total_counts'] = adata.X.sum(0)
    adata.var['n_cells_by_counts'] = (adata.X > 0).sum(0)

    adata.raw = adata
    # compute_embeddings(adata)
    adata.write_h5ad(f"{stem}.h5ad")

    adata_copy = adata[adata.obs['total_counts'] >= n_gene_min]
    # compute_embeddings(adata_copy)
    adata_copy.write_h5ad(f"{stem}_n_gene_min_{n_gene_min}.h5ad")



if __name__ == "__main__":
    fire.Fire(to_h5ad)
