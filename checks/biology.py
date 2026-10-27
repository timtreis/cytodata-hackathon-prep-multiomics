"""Sanity-check guide and condition labels against expected biology (pseudobulk per well).

Usage: python checks/biology.py <tables.zarr>
"""

import sys

import numpy as np
import pandas as pd
import spatialdata as sd
from scipy import stats

s = sd.read_zarr(sys.argv[1])
rna = s.tables["global-table-ThreePrimeUntargeted"]
img = s.tables["global-table-imaging"]
obs = rna.obs
keep = (obs["pass_qc"]).values
rna, img, obs = rna[keep], img[keep], obs[keep]

# pseudobulk log-CPM per (well, target_gene)
groups = obs.groupby(["Well", "condition", "target_gene"], observed=True).indices
keys = list(groups)
counts = np.vstack([np.asarray(rna.X[idx].sum(0)).ravel() for idx in groups.values()])
pb = pd.DataFrame(np.log1p(counts / counts.sum(1, keepdims=True) * 1e6), columns=rna.var_names)
pb.index = pd.MultiIndex.from_tuples(keys, names=["Well", "condition", "target_gene"])


def compare(genes, a, b, level, within=None):
    """Paired-by-well logFC (log-CPM) of group a vs b, with a paired t-test over wells."""
    sel = pb if within is None else pb.xs(within[1], level=within[0], drop_level=False)
    ga = sel.xs(a, level=level).groupby("Well").mean()[genes]
    gb = sel.xs(b, level=level).groupby("Well").mean()[genes]
    if level == "condition":  # TNFa and baseline are different wells: unpaired
        t, p = stats.ttest_ind(ga, gb)
    else:
        ga, gb = ga.loc[ga.index.intersection(gb.index)], gb.loc[ga.index.intersection(gb.index)]
        t, p = stats.ttest_rel(ga, gb)
    return pd.DataFrame({"log_fc": (ga.mean() - gb.mean()).round(3), "p": p}, index=genes)


print("== on-target transcript, KO vs NTC (all wells, paired by well)")
for g in ["TGFBR2", "SMAD7", "CTNNB1", "GSK3B", "IL1R1", "OR1I1"]:
    print(g, compare([g], g, "NTC", "target_gene").iloc[0].to_dict())

print("\n== TGFBR2 KO vs NTC: SMAD targets (expect down)")
print(compare(["PMEPA1", "TGFBI", "CDH2", "IGFBP7", "SERPINE1", "SKIL"], "TGFBR2", "NTC", "target_gene"))

print("\n== GSK3B / CTNNB1 KO vs NTC: Wnt targets (GSK3B KO expect up, CTNNB1 KO expect down)")
for ko in ["GSK3B", "CTNNB1"]:
    print(ko); print(compare(["AXIN2", "NKD1", "LGR5", "SP5"], ko, "NTC", "target_gene"))

print("\n== TNFa vs baseline, NTC cells: NF-kB early targets (expect up)")
print(compare(["NFKBIA", "TNFAIP3", "CXCL8", "IER3", "NFKBIZ", "ZFP36", "JUNB", "CCL2"], "TNFa", "baseline", "condition", within=("target_gene", "NTC")))

print("\n== TNFa vs baseline, all cells: protein (mean per well)")
prot = [v for v in img.var_names if v in ("Phospho-p38", "Phospho-ERK1_Phospho-ERK2", "ERK1_ERK2", "Phospho-HSP27")]
pw = pd.DataFrame(img[:, prot].X, columns=prot, index=obs.index).groupby([obs["Well"], obs["condition"]], observed=True).mean()
ta, tb = pw.xs("TNFa", level="condition"), pw.xs("baseline", level="condition")
print(pd.DataFrame({"log2_fc": np.log2(ta.mean() / tb.mean()).round(3), "p": stats.ttest_ind(ta, tb).pvalue}))
