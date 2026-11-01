"""Write QC-passing cells and common genes as raw counts for cycleVI.

Usage: python checks/prepare_cyclevi.py <rna.h5ad> <out.h5ad>
"""

import sys

import anndata as ad
import numpy as np

rna = ad.read_h5ad(sys.argv[1])
rna = rna[~rna.obs["qc_fail"].values].copy()
rna.layers.clear()
detected = np.asarray((rna.X > 0).mean(0)).ravel()
rna = rna[:, detected >= 0.01].copy()  # rare genes add noise at about 1,450 counts per cell
rna.obs = rna.obs[["Well", "condition", "guide", "target_gene", "rna_total_counts"]].astype({"Well": str})
print(rna.shape)
rna.write_h5ad(sys.argv[2])
