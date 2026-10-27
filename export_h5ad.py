"""Write the tables tier as plain h5ad files, one per modality, sharing the same rows and obs.

Usage: python export_h5ad.py <tables.zarr> <out dir>

Writes rna.h5ad, imaging.h5ad and guides.h5ad. Needs the `prep` dependency group (mantispy).
The imaging features go through mantispy's single-cell workflow: QC metrics, normalisation against
the NTC cells of the plate, feature selection and outlier flags. mantispy's qc_pass combines every QC flag.
"""

import sys
from pathlib import Path

import anndata as ad
import mantispy as mt
from mantispy._core.features import load_blocklist  # ponytail: private, no public accessor yet
import numpy as np
import pandas as pd
import zarr

PROTEIN = ["Phospho-ERK1_Phospho-ERK2", "ERK1_ERK2", "Phospho-HSP27", "Phospho-p38"]
PAINT = ["Mitochondria", "Cell-Membrane", "Nucleus", "Actin"]
CONTROL = ["Mouse_IgG", "Rabbit_IgG"]
OBS = [
    "Well", "Cell", "tile_id", "condition", "guide", "target_gene", "guide_confidence", "guide_assigned",
    "rna_total_counts", "rna_n_genes", "rna_pct_mt",
    "qc_low_rna", "qc_guide_unassigned", "qc_spatial_artifact", "qc_imaging_nan", "qc_no_nucleus",
]
OUR_FLAGS = [c for c in OBS if c.startswith("qc_")]


def read_table(root: Path, name: str) -> ad.AnnData:
    a = ad.io.read_elem(zarr.open(root / "tables" / f"global-table-{name}"))
    a.uns.pop("spatialdata_attrs", None)
    a.obs = a.obs[OBS]
    a.obs_names = a.obs["Well"].astype(str) + "_" + a.obs["Cell"].astype(str)
    return a


def channel_of(name: str) -> str | None:
    # longest match first, so "Phospho-ERK1_Phospho-ERK2" doesn't match as "ERK1_ERK2"
    for ch in sorted(PROTEIN + PAINT + CONTROL, key=len, reverse=True):
        if ch in name:
            return ch
    return None


def imaging(root: Path) -> ad.AnnData:
    img = read_table(root, "imaging")
    shape = img.obsm.pop("cell_geometry").drop(columns=["Xum_tile", "Yum_tile"])
    keep = ~img.var_names.str.startswith("Location_")  # pixel coordinates, not features
    X = np.hstack([img.X[:, keep], shape.to_numpy(np.float32)])
    var = pd.DataFrame(index=pd.Index([*img.var_names[keep], *shape.columns]))
    var["channel"] = [channel_of(n) for n in var.index]
    var["modality"] = np.select(
        [var.index.isin(shape.columns), var.index.str.startswith("Bkg_"), var["channel"].isin(CONTROL),
         var["channel"].isin(PROTEIN), var["channel"].isin(PAINT)],
        ["shape", "background", "control", "protein", "cell_paint"], default="other",
    )
    assert not (var["modality"] == "other").any(), var.index[var["modality"] == "other"][:5]
    # ponytail: mantispy can't parse Element's names (no object prefix, "_" in channel names), so we only
    # tell it which columns are areas; that's all calculate_qc_metrics needs
    var["feature"] = np.where(var.index.isin(["Area", "NuclearArea"]), "Area", None)
    return ad.AnnData(X=X, obs=img.obs, var=var, obsm={"spatial": img.obsm["spatial"].to_numpy()})


def run_mantispy(img: ad.AnnData) -> None:
    obs = img.obs
    obs["Metadata_Plate"] = "Cytos-MPS_TNFa"
    obs["Metadata_Well"] = obs["Well"].astype(str)
    obs["Metadata_ImageID"] = obs["Well"].astype(str) + "_" + obs["tile_id"].astype(str)
    obs["Metadata_ObjectType"] = "cell"
    obs["Metadata_ObjectNumber"] = obs["Cell"].astype(int)
    obs["Metadata_Perturbation"] = obs["target_gene"].astype(str)
    mt.io.stamp(img, resolution="object", feature_kind="measurement")
    mt.pp.annotate_controls(img, negcon=["NTC"])
    mt.pp.calculate_qc_metrics(img)
    # one plate, so TNFa stays in X; per-well normalisation would erase it, since condition is confounded with well
    mt.pp.normalize(img, by="Metadata_Plate", reference="negcon", keep_raw=True)
    # IgG isotype controls and background aren't biology, so keep them out of var["selected"]
    not_biology = img.var_names[img.var["modality"].isin(["control", "background"])]
    mt.pp.feature_select(img, blocklist=[*load_blocklist(), *not_biology])
    mt.pp.outliers(img, contamination=0.01, by="Metadata_Well")
    obs = img.obs
    obs["qc_pass"] = obs["qc_pass"] & ~obs["qc_outlier"] & ~obs[OUR_FLAGS].any(axis=1)


def main(zarr_path: str, out: str) -> None:
    root, out = Path(zarr_path), Path(out)
    out.mkdir(parents=True, exist_ok=True)

    img = imaging(root)
    run_mantispy(img)
    shared = img.obs  # every file gets the same obs, so rows and QC stay aligned

    rna = read_table(root, "ThreePrimeUntargeted")
    rna.obs = shared
    rna.var = rna.var[["gene_id", "gene_type", "seqname", "start", "end", "strand", "hgnc_id"]]
    rna.obsm = {"spatial": rna.obsm["spatial"].to_numpy()}

    guides = read_table(root, "B01-Cellecta")
    guides.obs = shared
    guides.var = guides.var[["target_gene", "target_sequence"]].astype(str)
    guides.obsm = {"target_assignment": guides.obsm["target_assignment"]}

    for name, a in {"rna": rna, "guides": guides}.items():
        a.uns["mantispy"] = img.uns["mantispy"]  # stamp, so mantispy accepts every file
        a.write_h5ad(out / f"{name}.h5ad", compression="gzip")
    # mt.io.write can't compress, so validate the same way it does and write gzipped
    mt.io.validate(img, raise_on_error=True)
    img.write_h5ad(out / "imaging.h5ad", compression="gzip")

    for f in sorted(out.glob("*.h5ad")):
        print(f"{f.name}  {f.stat().st_size / 1e9:.2f} GB")
    flags = [c for c in shared if c.startswith("qc_") and shared[c].dtype == bool and c != "qc_pass"]
    print(shared[[*flags, "qc_pass"]].mean().round(4).to_string())
    print("selected imaging features:", int(img.var["selected"].sum()), "of", img.n_vars)


if __name__ == "__main__":
    main(*sys.argv[1:])
