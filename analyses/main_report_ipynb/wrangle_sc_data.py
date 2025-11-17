import logging
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_adata_samples(samples, logger, min_genes=200, min_cells=3):
    adatas = {}

    for sampleName, path in samples.items():
        sample_adata = sc.read_10x_mtx(
            path,
            cache=True
        )
        sample_adata.var_names_make_unique()

        # mitochondrial genes, "MT-" for human, "Mt-" for mouse
        sample_adata.var["mt"] = sample_adata.var_names.str.startswith("MT-")
        # ribosomal genes
        sample_adata.var["ribo"] = sample_adata.var_names.str.startswith(("RPS", "RPL"))
        # hemoglobin genes
        sample_adata.var["hb"] = sample_adata.var_names.str.contains("^HB[^(P)]")

        ## Filtering on a per sample basis
        sc.pp.filter_cells(sample_adata, min_genes=min_genes)

        sc.pp.filter_genes(sample_adata, min_cells=min_cells)

        sc.pp.calculate_qc_metrics(
            sample_adata, 
            qc_vars=["mt", "ribo", "hb"], 
            inplace=True, 
            log1p=True
        )

        sc.pl.violin(
            sample_adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
            jitter=0.4,
            multi_panel=False,
            show=False
        )

        plt.savefig(f"../../../../html_local/report_figures/{sampleName}_qc_violin.pdf")
        plt.close()

        adatas[sampleName] = sample_adata
        logger.info(f"Loaded {sampleName}: {sample_adata.n_obs} cells, {sample_adata.n_vars} genes")
    return logger, adatas



