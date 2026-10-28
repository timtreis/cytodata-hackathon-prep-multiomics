# SBI2 CytoData hackathon: multiomics dataset

The dataset is a pooled CRISPR screen of 258,896 A549 cells, measured with Element Biosciences' Cytos on the AVITI24. Each cell has cell painting images, four protein stains, its CRISPR guide and a 3' transcriptome, and because all of these measurements come from the same cell, the tables line up row for row.

Element Biosciences published the [original data](https://www.elementbiosciences.com/resources/dataset-explore-multimodal-pooled-screening-on-aviti24).

## Experimental design

- The cells are A549 cells that express Cas9. Element's metadata says HeLa, but the expression data show A549 (see Known issues).
- Each cell carries one of 28 guides. There are 4 guides each against CTNNB1, GSK3B, IL1R1, SMAD7 and TGFBR2, 4 against OR1I1, which A549 cells don't express and which controls for the effect of cutting, and 4 non-targeting guides (`NTC`).
- Six wells are untreated (A1, B1, C1, D2, E2, F2), and six wells were treated with TNFα for 30 minutes (A2, B2, C2, D1, E1, F1).

## Download

| File | Contents | Size |
|---|---|---|
| `rna.h5ad` | RNA counts for 54,949 genes | 0.7 GB |
| `imaging.h5ad` | 1,131 cell painting, protein and shape features | 1.6 GB |
| `guides.h5ad` | Guide counts for 28 guides | 0.04 GB |
| `subset_A1_A2.zarr.zip` | Wells A1 and A2 with images, masks and transcript locations (SpatialData) | 17 GB |
| `full_<well>.zarr.zip` | One file per well (`full_A1` to `full_F2`) with images, masks and transcript locations (SpatialData) | 8.8 to 9.6 GB each |

Most analyses only need the three `.h5ad` files. You only need the zarr files if you want to work with the images.

```bash
URL=https://exampledata.scverse.org/sbi2-cytodata/multiomics
for f in rna.h5ad imaging.h5ad guides.h5ad; do curl -O $URL/$f; done
shasum -a 256 -c <<EOF
983b74f8a7f9a70b4dfdc1cc4cbc54168dd2c8b36099724bfb740b331f66d984  rna.h5ad
8be68212c1f710095e7f0400c28a76abf4a5465c6715d6d9856c9ebf76db8cc1  imaging.h5ad
1dd75ee3c80931b3e35698ffddad419915aa14275f9fac91838e3bf56c26d093  guides.h5ad
EOF
```

The zarr files are at the same URL, e.g. `curl -O $URL/full_B2.zarr.zip` for well B2. Their sha256 checksums are:

```
bfa0442ec3183fc379f711a89f65487c10d96119340a92fb184180c539280317  subset_A1_A2.zarr.zip
16ff5051e595b9d96891b9efd05d8a5e1942b090d5b8710ca87226417adff1da  full_A1.zarr.zip
9a32b79ea46a7db2caeb66940a12ed909aed62d59f1c785c0db1b4460f7b75b0  full_A2.zarr.zip
287d91d0fd7a00fe8ba92395cea902332fbba1bb96889659b45a46f18a24f80e  full_B1.zarr.zip
001c757a4946eda97bede539c616e8241853cdbde82889653591ca2db94707be  full_B2.zarr.zip
3839173bc68be210869c4179f090f62c4da986b3d81d05ea63d123b80ccc74d9  full_C1.zarr.zip
93a2bedf5a4d5a295798be9b932df39924a454dd25374c5291c2586b537e0567  full_C2.zarr.zip
21cc889de09b6d56dfeac8bf98bbbae3e0bb7ccff927676c0bb9a31dddf63a39  full_D1.zarr.zip
3ab0fd284ce8d70f075913995b142694648a9f817ff8a9a8a0e5bbf392523b50  full_D2.zarr.zip
52ee000ba873b506b4ce35685300626e00d13b1865d5ee900d51fb64629c12c4  full_E1.zarr.zip
add56e9d6a7ee2b8b8cd7da34c1ad7aa872462a1f58eb492efc14d1cc1b5ccea  full_E2.zarr.zip
87c2ac200b341ef6fe38f94e3d208292a726671d65fd4419502ab021b6848f85  full_F1.zarr.zip
e06676e7b8d054e2afb132cf579d4251e4abf7cebc34ae3838bdf5fdbba77fff  full_F2.zarr.zip
```

```python
import anndata as ad

rna = ad.read_h5ad("rna.h5ad")
imaging = ad.read_h5ad("imaging.h5ad")
guides = ad.read_h5ad("guides.h5ad")
cells = rna.obs["qc_pass"].values  # the same rows in all three files
```

The notebook `notebooks/01_starter.ipynb` downloads the h5ad files, joins the modalities, builds a first joint embedding, and scores each representation on how well it separates TNFα-treated cells and each knockout from the non-targeting cells. It runs on a laptop in about 3 minutes and needs less than 4 GB of memory. You can open the zarr files with `spatialdata.read_zarr` after unzipping them.

To run the notebook, install the packages it uses:

```bash
pip install anndata scanpy scikit-learn matplotlib "mantispy[map]>=0.1.0a1"
```

## Contents of the h5ad files

The three files share the same rows and the same `obs`, and we didn't remove any cells. The main `obs` columns are:

- `condition` (`baseline` or `TNFa`), `Well`, `guide`, `target_gene` (`NTC` for non-targeting guides, `Undetermined` when no guide was called) and `guide_confidence`
- `guide_assigned`, which is true when the guide confidence is above 0.5. With that cutoff, 87% of cells have a guide, which matches Element's number.
- The QC flags, which are true for cells with a problem: `qc_low_rna` (fewer than 200 RNA counts), `qc_guide_unassigned`, `qc_spatial_artifact`, `qc_imaging_nan`, `qc_no_nucleus`, and mantispy's `qc_area_outlier` and `qc_outlier`
- `qc_pass`, which is true when none of the QC flags is set. 79.6% of cells pass.

We processed `imaging.h5ad` with [mantispy](https://github.com/scverse/mantispy). `X` holds robust z-scores against all non-targeting cells, the raw values are in `layers["raw"]`, and `var["selected"]` marks the 369 features that mantispy's feature selection kept, leaving out the IgG controls and the background. `var["modality"]` says whether a feature comes from cell painting, protein, cell shape, the IgG controls or the background. The zarr files hold Element's tables with an earlier QC flag, `pass_qc`, which doesn't exclude mantispy's outliers.

## Checks on the biology

We ran `checks/biology.py` to make sure the guide and treatment labels are correct, and the results match the known biology:

- Each knockout lowers the mRNA of its target gene, except for SMAD7 and OR1I1.
- The TGFBR2 knockout lowers TGF-β target genes such as PMEPA1 and TGFBI.
- TNFα raises NF-κB target genes such as CXCL8, NFKBIA and TNFAIP3, and it raises phosphorylated HSP27.

Wnt signalling is low in A549 cells, so the CTNNB1 and GSK3B knockouts have little effect.

## Known issues

- Each well is either untreated or treated, so the condition is confounded with the well. Use `Well` as a covariate rather than normalising it away, because normalising per well removes the TNFα response.
- The `*_NuclearRatio` imaging features are missing for up to 14% of cells.
- Element's metadata lists the cells as HeLa, but the data show they are A549. The cells express Y chromosome genes, they have high NRF2 target expression, and they don't express CDKN2A, and all of these fit A549 but not HeLa.

## Reproduce the files

```bash
uv sync --group prep
python prep.py Cytos-MPS_TNFa-CRISPR-Screen.zarr out
pip install "mantispy @ git+https://github.com/scverse/mantispy.git@3552713"  # the version the h5ad files were made with
python export_h5ad.py out/tables.zarr out/h5ad
python split_wells.py out/full.zarr out/wells
python checks/biology.py out/tables.zarr
```
