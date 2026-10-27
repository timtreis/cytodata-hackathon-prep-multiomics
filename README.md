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
| `full.zarr.zip` | All 12 wells with images, masks and transcript locations (SpatialData) | 104 GB |

Most analyses only need the three `.h5ad` files. You only need the zarr files if you want to work with the images.

```bash
URL=https://exampledata.scverse.org/sbi2-cytodata/multiomics
for f in rna.h5ad imaging.h5ad guides.h5ad; do curl -O $URL/$f; done
shasum -a 256 -c <<EOF
6bcaf34969e6277c3b1ba334520dbc4a627da1b291d8ee64b304d43b53e59e7a  rna.h5ad
0d922ae40fad54226f26af4a4cd3f9d62adff56461c338adac4a4e0cca7d48f0  imaging.h5ad
5c0b0bdaefee317fe52541aa9213b0c93a96cf8678dc538448ca536460556feb  guides.h5ad
EOF
```

The zarr files are at the same URL, e.g. `$URL/subset_A1_A2.zarr.zip`.

```python
import anndata as ad

rna = ad.read_h5ad("rna.h5ad")
imaging = ad.read_h5ad("imaging.h5ad")
guides = ad.read_h5ad("guides.h5ad")
cells = rna.obs["qc_pass"].values  # the same rows in all three files
```

The notebook `notebooks/01_starter.ipynb` downloads the h5ad files, joins the modalities, and builds a first joint embedding. It runs on a laptop in about 3 minutes and needs less than 4 GB of memory. You can open the zarr files with `spatialdata.read_zarr` after unzipping them.

## Contents of the h5ad files

The three files share the same rows and the same `obs`, and we didn't remove any cells. The main `obs` columns are:

- `condition` (`baseline` or `TNFa`), `Well`, `guide`, `target_gene` (`NTC` for non-targeting guides, `Undetermined` when no guide was called) and `guide_confidence`
- `guide_assigned`, which is true when the guide confidence is above 0.5. With that cutoff, 87% of cells have a guide, which matches Element's number.
- The QC flags, which are true for cells with a problem: `qc_low_rna` (fewer than 200 RNA counts), `qc_guide_unassigned`, `qc_spatial_artifact`, `qc_imaging_nan`, `qc_no_nucleus`, and mantispy's `qc_area_outlier` and `qc_outlier`
- `qc_pass`, which is true when none of the QC flags is set. 79.6% of cells pass.

We processed `imaging.h5ad` with [mantispy](https://github.com/scverse/mantispy). `X` holds robust z-scores against all non-targeting cells, the raw values are in `layers["raw"]`, and `var["selected"]` marks the 465 features that mantispy's feature selection kept. `var["modality"]` says whether a feature comes from cell painting, protein, cell shape, the IgG controls or the background. The zarr files hold Element's tables with an earlier QC flag, `pass_qc`, which doesn't exclude mantispy's outliers.

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
python export_h5ad.py out/tables.zarr out/h5ad
python checks/biology.py out/tables.zarr
```
