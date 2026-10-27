# CytoData hackathon @ SBI2: Element Cytos multiomics CRISPR screen

Pooled CRISPR screen in A549 cells. Cell painting, a protein panel, CRISPR guide identity and the 3' transcriptome are all read out from **the same cells** in situ on Element's AVITI24 (Direct In Sample Sequencing).

Source: [Element Biosciences, "Explore Multimodal Pooled Screening with Cytos"](https://www.elementbiosciences.com/resources/dataset-explore-multimodal-pooled-screening-on-aviti24). Public S3 download: `Cytos-MPS_TNFa-CRISPR-Screen.zarr.zip` (111 GB, SpatialData).

## Design

- A549 Cas9 (Cellecta) with a 28-guide library: CTNNB1, GSK3B, IL1R1, SMAD7, TGFBR2 and OR1I1 (4 guides each, OR1I1 is a cutting control with no expected phenotype), plus 4 non-targeting guides (`NTC`).
- 12 wells: 6 baseline (A1, B1, C1, D2, E2, F2) and 6 with 30 min TNFα (A2, B2, C2, D1, E1, F1). About 21k cells per well, 258,896 cells in total.
- Modalities:
  - cell paint: mitochondria, membrane, nucleus, actin
  - protein IF: pERK1/2, ERK1/2, pHSP27, p38, plus IgG controls
  - 3' whole transcriptome
  - guide reads

## Tiers

| Tier | Contents | Size |
|---|---|---|
| `full` | Everything: images, cell + nuclear masks, transcript points, tables | 104 GB |
| `subset_A1_A2` | One baseline (A1) and one TNFα well (A2), all elements | 17 GB |
| `tables` | Tables only, no images | 1.5 GB |

Open any tier with `spatialdata.read_zarr(...)`.

## Tables

The four cell tables have the same rows, in the same order. Cells are linked to the `<well>-labels-Cell` masks through `region_key` / `Cell`.

| Table | `X` |
|---|---|
| `global-table-ThreePrimeUntargeted` | RNA counts (54,949 genes; nuclear counts in `layers["nuclear"]`) |
| `global-table-imaging` | 1,150 intensity/morphology features (cell paint + protein). `*_NuclearRatio` features are NaN for up to 14% of cells: drop or impute them. |
| `global-table-B01-Cellecta` | guide counts (28 guides) |
| `global-table-B01-QC_Panel` | targeted GAPDH/ACTB control probes |
| `global-table-metadata` | run-level stats in `uns` |

Columns added to `obs` of every cell table:

- `condition`: `baseline` or `TNFa`.
- `guide`, `target_gene`, `guide_confidence`: Element's most abundant guide per cell and its share of guide reads.
- `guide_assigned`: confidence > 0.5. This reproduces Element's reported 87% assignment rate.
- `rna_total_counts`, `rna_n_genes`, `rna_pct_mt`.
- QC flags (nothing is removed):
  - `qc_low_rna`: fewer than 200 RNA counts.
  - `qc_guide_unassigned`
  - `qc_spatial_artifact`: Element's per-cell artifact call across any modality.
  - `qc_imaging_nan`
  - `qc_no_nucleus`
- `pass_qc`: no QC flag set.

## Cell line note

Element's well metadata lists `cell_type = HELA`, but the dataset page says A549. The expression data agree with A549:

- **Male:** Y-chromosome genes are expressed (RPS4Y1, EIF1AY, DDX3Y) and XIST is absent. HeLa is female.
- **KEAP1-mutant NRF2 program:** very high NQO1, AKR1C1/2/3, AKR1B10, ALDH3A1 and ALDH1A1.
- **CDKN2A absent:** deleted in A549. HeLa expresses p16 strongly because of HPV E7.

We set `cell_type` to `A549` and kept Element's value in `cell_type_original`.

## Reproduce

```bash
uv sync
python prep.py <unzipped .zarr> <out dir>
```

Processing reference: [Elembio/cytoprofiling](https://github.com/Elembio/cytoprofiling).
