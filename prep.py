"""Annotate the Element Cytos CRISPR-screen SpatialData and write hackathon tiers.

Usage: python prep.py <unzipped .zarr> <out dir>

- annotates the full zarr in place (QC flags + guide/condition columns on every cell table,
  cell line corrected to A549)
- <out>/tables.zarr: tables only, the input to export_h5ad.py
- <out>/subset_A1_A2.zarr: one baseline + one TNFa well with images and labels
- <out>/qc/: summary plots
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spatialdata as sd

CELL_TABLES = ["ThreePrimeUntargeted", "imaging", "B01-Cellecta", "B01-QC_Panel"]
GUIDE_CONF = 0.5  # reproduces Element's reported 87% assignment rate
MIN_RNA_COUNTS = 200  # ~1st percentile


def annotations(sdata: sd.SpatialData) -> pd.DataFrame:
    t = {k: sdata.tables[f"global-table-{k}"] for k in CELL_TABLES}
    ref = t["imaging"].obs[["Well", "cell_id_global"]]
    for a in t.values():
        assert a.obs[["Well", "cell_id_global"]].equals(ref), "cell tables are not row-aligned"

    rna, img = t["ThreePrimeUntargeted"], t["imaging"]
    ga = t["B01-Cellecta"].obsm["target_assignment"]
    mt = rna.var_names.str.startswith("MT-")
    counts = np.asarray(rna.X.sum(1)).ravel()

    df = pd.DataFrame(index=ref.index)
    df["condition"] = np.where(img.obs["WellLabel"].astype(str).str.endswith("TNFa"), "TNFa", "baseline")
    df["guide"] = ga["MostAbundantTarget"].astype(str).values
    df["target_gene"] = df["guide"].str.rsplit("_", n=1).str[0].replace({"NT_Human_CRko_BB": "NTC"})
    df["guide_confidence"] = ga["Confidence"].values
    df["guide_assigned"] = (ga["Confidence"].values > GUIDE_CONF) & (df["guide"] != "Undetermined")
    df["rna_total_counts"] = counts
    df["rna_n_genes"] = np.asarray((rna.X > 0).sum(1)).ravel()
    df["rna_pct_mt"] = np.asarray(rna.X[:, mt].sum(1)).ravel() / np.maximum(counts, 1) * 100
    df["qc_low_rna"] = counts < MIN_RNA_COUNTS
    df["qc_guide_unassigned"] = ~df["guide_assigned"]
    df["qc_spatial_artifact"] = rna.obsm["spatial_artifacts"]["any"].astype(bool).values
    # *_NuclearRatio is NaN for up to 14% of otherwise fine cells (ratio artifact), so it doesn't count
    df["qc_imaging_nan"] = np.isnan(img.X[:, ~img.var_names.str.endswith("_NuclearRatio")]).any(1)
    df["qc_no_nucleus"] = img.obsm["cell_geometry"]["NuclearArea"].values == 0
    df["pass_qc"] = ~df.filter(like="qc_").any(axis=1)
    df["condition"] = df["condition"].astype("category")
    df["target_gene"] = df["target_gene"].astype("category")
    return df


def annotate(sdata: sd.SpatialData, ann: pd.DataFrame) -> None:
    for k in CELL_TABLES:
        a = sdata.tables[f"global-table-{k}"]
        a.obs = a.obs.drop(columns=ann.columns, errors="ignore").join(ann)
    wells = sdata.tables["global-table-metadata"].uns["wells"]
    if "cell_type_original" not in wells:
        wells["cell_type_original"] = wells["cell_type"]
    # ponytail: Element ships HELA here; expression (Y genes, NRF2 program, CDKN2A loss) shows A549, see README
    wells["cell_type"] = "A549"


def subset(sdata: sd.SpatialData, wells: list[str]) -> sd.SpatialData:
    keep = [n for _, n, _ in sdata.gen_spatial_elements() if n.split("-")[0] in wells]
    tables = {}
    for name, a in sdata.tables.items():
        if a.n_obs == 0:
            tables[name] = a
            continue
        a = a[a.obs["Well"].isin(wells)].copy()
        a.obs["region_key"] = a.obs["region_key"].astype(str).astype("category")
        a.uns["spatialdata_attrs"]["region"] = sorted(a.obs["region_key"].cat.categories)
        tables[name] = a
    elements = {n: sdata[n] for n in keep}
    return sd.SpatialData.init_from_elements(elements | tables, attrs=sdata.attrs)


def qc_plots(ann: pd.DataFrame, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    ann["rna_total_counts"].clip(upper=8000).plot.hist(bins=100, ax=ax[0, 0], title="RNA counts/cell")
    ax[0, 0].axvline(MIN_RNA_COUNTS, color="k", ls="--")
    ann["guide_confidence"].plot.hist(bins=50, ax=ax[0, 1], title="guide confidence")
    ax[0, 1].axvline(GUIDE_CONF, color="k", ls="--")
    pd.crosstab(ann.loc[ann.guide_assigned, "target_gene"], ann.loc[ann.guide_assigned, "condition"]).plot.barh(ax=ax[1, 0], title="assigned cells per target")
    ann.filter(like="qc_").mean().mul(100).plot.barh(ax=ax[1, 1], title="% cells flagged")
    fig.tight_layout()
    fig.savefig(out / "qc_summary.png", dpi=120)


def main(zarr: str, out: str) -> None:
    out = Path(out)
    sdata = sd.read_zarr(zarr)
    ann = annotations(sdata)
    annotate(sdata, ann)
    for name in [*(f"global-table-{k}" for k in CELL_TABLES), "global-table-metadata"]:
        # spatialdata can't overwrite in place; tables are in memory and the .zip is the backup
        sdata.delete_element_from_disk(name)
        sdata.write_element(name)

    sd.SpatialData(tables=dict(sdata.tables), attrs=sdata.attrs).write(out / "tables.zarr", overwrite=True)
    subset(sdata, ["A1", "A2"]).write(out / "subset_A1_A2.zarr", overwrite=True)

    qc_plots(ann, out / "qc")
    summary = ann.groupby("condition", observed=True)[["pass_qc", *ann.filter(like="qc_").columns]].mean()
    summary.loc["all"] = ann[summary.columns].mean()
    print(f"cells: {len(ann)}\n{summary.round(4)}")
    print(ann[ann.pass_qc].groupby(["target_gene", "condition"], observed=True).size().unstack())


if __name__ == "__main__":
    main(*sys.argv[1:])
