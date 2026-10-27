# SBI2 CytoData hackathon: multiomics dataset

The dataset is a pooled CRISPR screen in A549 cells, measured with Element Biosciences' Cytos on the AVITI24. Each cell has cell painting images, four protein stains, its guide and a 3' transcriptome, and because all four measurements come from the same cell, the tables line up row for row.

Element Biosciences published the [original data](https://www.elementbiosciences.com/resources/dataset-explore-multimodal-pooled-screening-on-aviti24).

## Experimental design

- The cells are A549 cells that express Cas9, and they carry one of 28 guides. There are 4 guides each against CTNNB1, GSK3B, IL1R1, SMAD7, TGFBR2 and OR1I1, plus 4 non-targeting guides. A549 cells don't express OR1I1, so its guides control for the effect of cutting.
- Six wells are untreated (A1, B1, C1, D2, E2, F2), and six wells were treated with TNFα for 30 minutes (A2, B2, C2, D1, E1, F1).
- The dataset has 258,896 cells, and 82% of them pass QC.

## Download

| File | Contents | Size |
|---|---|---|
| `tables.zarr` | Per-cell tables only | 1.5 GB |
| `subset_A1_A2.zarr` | Wells A1 and A2 with images, masks and transcript locations | 17 GB |
| `full.zarr` | All 12 wells | 104 GB |

Link: TBA

You can open each file with `spatialdata.read_zarr(path)`. The notebook `notebooks/01_starter.ipynb` shows how to load and explore the tables.

## Contents

The images have 14 channels:

- Cell painting stains for mitochondria, membrane, nucleus and actin
- Antibody stains for pERK1/2, ERK1/2, pHSP27 and p-p38, plus two IgG controls
- Four background channels

The `labels` hold the cell and nuclear masks, and the `points` hold each RNA and guide read with its position.

There are four per-cell tables, and they share the same rows. The `Cell` column in each table matches the label value in the cell mask.

| Table | Contents |
|---|---|
| `global-table-ThreePrimeUntargeted` | RNA counts for 54,949 genes |
| `global-table-imaging` | 1,150 intensity and morphology features |
| `global-table-B01-Cellecta` | Guide counts |
| `global-table-B01-QC_Panel` | Counts from GAPDH and ACTB control probes |

We added the following columns to `obs` in each table, and we didn't remove any cells:

- `condition` (`baseline` or `TNFa`), `guide`, `target_gene` and `guide_confidence`
- `guide_assigned`, which is true when the guide confidence is above 0.5. With that cutoff, 87% of cells have a guide, which matches Element's number.
- The QC flags `qc_low_rna` (fewer than 200 RNA counts), `qc_guide_unassigned`, `qc_spatial_artifact`, `qc_imaging_nan` and `qc_no_nucleus`
- `pass_qc`, which is true when none of the QC flags is set

## Checks on the biology

We ran `checks/biology.py` to make sure the guide and treatment labels are correct, and the results match the known biology:

- Each knockout lowers the mRNA of its target gene, except for SMAD7 and OR1I1.
- The TGFBR2 knockout lowers TGF-β target genes such as PMEPA1 and TGFBI.
- TNFα raises NF-κB target genes such as CXCL8, NFKBIA and TNFAIP3, and it raises phosphorylated HSP27.

Wnt signalling is low in A549 cells, so the CTNNB1 and GSK3B knockouts have little effect.

## Known issues

- The `*_NuclearRatio` features are missing for up to 14% of cells.
- The `Bkg_*` features measure the background, not the cells.
- Wells differ in sequencing depth and staining, so normalise each well separately.
- Element's metadata lists the cells as HeLa, but the data show they are A549. The cells express Y chromosome genes, they have high NRF2 target expression, and they don't express CDKN2A, and all of these fit A549 but not HeLa. We changed `cell_type` to A549 and kept Element's value in `cell_type_original`.

## Reproduce the files

```bash
uv sync
python prep.py Cytos-MPS_TNFa-CRISPR-Screen.zarr out
python checks/biology.py out/tables.zarr
```
