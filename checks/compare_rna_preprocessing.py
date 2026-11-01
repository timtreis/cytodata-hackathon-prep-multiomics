"""Compare RNA preprocessing for notebook 01's baselines on the same cells, over several samples.

Usage: python checks/compare_rna_preprocessing.py <dir with h5ads> <out.csv>
"""

import sys
from pathlib import Path

import anndata as ad
import mantispy as mt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.model_selection import cross_val_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler

N_CELLS = 60_000
SEEDS = [0, 1, 2, 3]
KNOCKOUTS = ["TGFBR2", "SMAD7", "IL1R1", "GSK3B", "CTNNB1", "OR1I1"]
N_PCS = {"rna": 30, "protein": 10, "paint": 30}


def rna_block(rna: ad.AnnData, new: bool) -> np.ndarray:
    a = rna.copy()
    if new:
        a = a[:, np.asarray((a.X > 0).mean(0)).ravel() >= 0.01].copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000)
    x = StandardScaler().fit_transform(a[:, a.var["highly_variable"]].X.toarray())
    return x.clip(-10, 10) if new else x


def knn(x, y):
    return cross_val_score(KNeighborsClassifier(15), x, y, cv=5, scoring="balanced_accuracy").mean()


def score(pcs: dict, obs: pd.DataFrame) -> list[dict]:
    target = np.asarray(obs["target_gene"], dtype=str)
    condition = np.asarray(obs["condition"], dtype=str)
    well, guide = obs["Well"].astype(str).values, obs["guide"].astype(str).values
    meta = obs[["Well", "guide", "target_gene", "condition"]].astype(str).groupby([well, guide], observed=True).first().reset_index(drop=True)
    meta = meta.rename(columns={"Well": "Metadata_Well", "guide": "Metadata_Guide", "target_gene": "Metadata_Perturbation", "condition": "Metadata_Condition"})
    meta["Metadata_Plate"] = "Cytos-MPS_TNFa"
    means = {k: pd.DataFrame(v).groupby([well, guide]).mean().to_numpy(np.float32) for k, v in pcs.items()}
    agg = ad.AnnData(X=means["joint"], obs=meta, obsm={f"X_{k}": v for k, v in means.items()})
    mt.io.stamp(agg, resolution="aggregate", grouped_by=["Metadata_Well", "Metadata_Guide"], feature_kind="embedding")
    mt.pp.annotate_controls(agg, negcon=["NTC"])
    rows = []
    for name, x in pcs.items():
        rows.append({"metric": "knn", "representation": name, "target": "TNFa", "value": knn(x, condition)})
        for ko in KNOCKOUTS:
            per = [knn(x[m], target[m]) for m in (((condition == c) & np.isin(target, [ko, "NTC"])) for c in ("baseline", "TNFa"))]
            rows.append({"metric": "knn", "representation": name, "target": ko, "value": float(np.mean(per))})
        mt.tl.map(agg, mode="activity", batch="Metadata_Condition", use_rep=f"X_{name}", null_size=1000, key_added=f"map_{name}")
        table = agg.uns["mantispy"][f"map_{name}"].groupby("Metadata_Perturbation")["score"].mean()
        rows += [{"metric": "map", "representation": name, "target": k, "value": float(v)} for k, v in table.items()]
    return rows


def main(data: str, out: str) -> None:
    data = Path(data)
    rna_all = ad.read_h5ad(data / "rna.h5ad", backed="r")
    img_all = ad.read_h5ad(data / "imaging.h5ad", backed="r")
    passing = np.flatnonzero(~rna_all.obs["qc_fail"].values)
    selected, modality = img_all.var["selected"].values, img_all.var["modality"].astype(str).values
    rows = []
    for seed in SEEDS:
        cells = np.sort(np.random.default_rng(seed).choice(passing, N_CELLS, replace=False))
        rna = rna_all[cells].to_memory()
        rna.layers.clear()
        img = img_all[cells].to_memory()
        x = np.nan_to_num(np.clip(np.asarray(img.X), -10, 10))
        fixed = {"protein": x[:, selected & (modality == "protein")], "paint": x[:, selected & np.isin(modality, ["cell_paint", "shape"])]}
        fixed = {k: PCA(N_PCS[k], random_state=0).fit_transform(v) for k, v in fixed.items()}
        for label, new in {"old": False, "new": True}.items():
            pcs = {"rna": PCA(N_PCS["rna"], random_state=0).fit_transform(rna_block(rna, new)), **fixed}
            pcs["joint"] = np.hstack([p / np.sqrt((p**2).sum(1).mean()) for p in pcs.values()])
            rows += [{**r, "seed": seed, "preprocessing": label} for r in score(pcs, rna.obs)]
        pd.DataFrame(rows).to_csv(out, index=False)
        print("done seed", seed, flush=True)

    df = pd.DataFrame(rows)
    pd.set_option("display.width", 220, "display.max_columns", 30)
    for metric in ["map", "knn"]:
        d = df[(df["metric"] == metric)]
        print(f"\n=== {metric}: mean over seeds (sd)")
        mean = d.groupby(["preprocessing", "representation", "target"])["value"].mean().unstack().round(3)
        sd = d.groupby(["preprocessing", "representation", "target"])["value"].std().unstack().round(3)
        print(mean.to_string()); print("sd\n" + sd.to_string())
    m = df[df["metric"] == "map"].pivot_table(index=["seed", "preprocessing", "target"], columns="representation", values="value").reset_index()
    print("\n=== margins per seed")
    rows = []
    for (seed, label), g in m.groupby(["seed", "preprocessing"]):
        g = g.set_index("target")
        for ko in ["TGFBR2", "CTNNB1"]:
            rows.append({"seed": seed, "preprocessing": label, "knockout": ko, "rna - OR1I1": g.loc[ko, "rna"] - g.loc["OR1I1", "rna"], "rna - joint": g.loc[ko, "rna"] - g.loc[ko, "joint"], "paint - joint": g.loc[ko, "paint"] - g.loc[ko, "joint"]})
    margins = pd.DataFrame(rows)
    print(margins.round(3).to_string(index=False))
    print(margins.groupby(["knockout", "preprocessing"])[["rna - OR1I1", "rna - joint", "paint - joint"]].agg(["mean", "std"]).round(3).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:])
