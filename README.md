# CytoData hackathon @ SBI2: Element Cytos multiomics CRISPR screen

A pooled CRISPR knockout screen in A549 lung adenocarcinoma cells, with and without TNFα. Four modalities are read out from **the same cells**, in situ, on one instrument (Element AVITI24 with Direct In Sample Sequencing):

1. **Cell painting:** 4 stains, full-resolution images plus segmentation masks.
2. **Protein immunofluorescence:** 4 MAPK/stress markers plus 2 IgG controls.
3. **CRISPR guide identity:** sequenced in situ, so no barcode imaging or decoding is needed.
4. **3' whole transcriptome:** about 1,450 counts per cell (median), with every transcript located in space.

There's no cell matching, no separate assays and no integration guesswork: one row is one cell in every table.

Source: [Element Biosciences, "Explore Multimodal Pooled Screening with Cytos"](https://www.elementbiosciences.com/resources/dataset-explore-multimodal-pooled-screening-on-aviti24). The data here is Element's release with QC columns added (see [Changes from the original](#changes-from-the-original)).

## Experimental design

| | |
|---|---|
| Cells | A549 Cas9 (Cellecta), fixed |
| Library | 28 guides: 4 each against CTNNB1, GSK3B, IL1R1, SMAD7, TGFBR2, OR1I1, plus 4 non-targeting (`NTC`) |
| Stimulus | 6 wells baseline (A1, B1, C1, D2, E2, F2); 6 wells 30 min TNFα (A2, B2, C2, D1, E1, F1) |
| Size | 258,896 in total, about 21k per well; 212,857 (82%) pass QC |

The targets cover three pathways:
- **TGF-β:** TGFBR2 (receptor) and SMAD7 (negative regulator).
- **Wnt/β-catenin:** CTNNB1, and GSK3B, which degrades β-catenin.
- **IL-1:** IL1R1.

OR1I1 is an olfactory receptor that A549 doesn't express. It's a cutting control: Cas9 cuts DNA, but no phenotype is expected.

Cells passing QC per target:

| | CTNNB1 | GSK3B | IL1R1 | NTC | OR1I1 | SMAD7 | TGFBR2 |
|---|---|---|---|---|---|---|---|
| baseline | 24,022 | 17,318 | 7,323 | 14,125 | 36,217 | 5,300 | 3,672 |
| TNFα | 22,790 | 17,674 | 6,816 | 13,677 | 35,129 | 5,252 | 3,542 |

## What's known to be in there

We checked these with pseudobulk per well on QC-passing cells. Script: `checks/biology.py`.

- **Knockouts work.** The targeted gene's own transcript is lower than in NTC cells for TGFBR2, CTNNB1, GSK3B and IL1R1 (log-CPM −0.4 to −0.75, p < 1e-5). OR1I1 shows no change, as expected.
- **TGFBR2 KO lowers TGF-β/SMAD targets:** PMEPA1, TGFBI, CDH2, IGFBP7 and SERPINE1 (log-CPM −0.3 to −0.7, p < 1e-6). Element also reports a matching morphology phenotype.
- **TNFα triggers a strong early NF-κB response:** CXCL8 +3.5, ZFP36 +2.4, TNFAIP3 +2.1, NFKBIA +2.0, NFKBIZ +1.5 (log-CPM, NTC cells, p < 1e-7).
- **TNFα activates stress MAPKs at the protein level:** pHSP27 up 2.7× and pERK1/2 up 1.3×. Phospho-p38 itself barely moves (1.05×).
- **Wnt signalling is close to silent in A549.** GSK3B KO raises AXIN2 only slightly (+0.11) and CTNNB1 KO has no detectable effect on Wnt targets. Don't expect a strong Wnt phenotype.

Open questions to work on:
- Which knockouts change morphology or protein without a transcriptional signature, or the other way round?
- Does knocking out TGF-β pathway genes change the TNFα response?
- How well does one modality predict another?
- How should guide-assignment uncertainty be handled?

## Data tiers

| Tier | Contents | Size | Good for |
|---|---|---|---|
| `tables.zarr` | All per-cell tables, no images | 1.5 GB | Most analyses: perturbation effects, multimodal integration, prediction |
| `subset_A1_A2.zarr` | One baseline (A1) and one TNFα well (A2), everything | 17 GB | Image workflows on a laptop, segmentation, viewers |
| full | All 12 wells, everything | 104 GB | Image-based models at scale |

Download: scverse S3 (link TBA). Open any tier with `spatialdata.read_zarr(path)`. Start with `notebooks/01_starter.ipynb`.

## What's inside a SpatialData object

| Element | Per well | Description |
|---|---|---|
| `images/<well>-images-cell_paint` | 1 | uint16, about 14.7k × 16.5k px, 6 pyramid levels, 14 channels (listed below) |
| `labels/<well>-labels-Cell`, `<well>-labels-Nuclear` | 2 | Cell and nuclear segmentation masks; label value = `Cell` in the tables |
| `shapes/` | 2 | Cell and nuclear boundaries |
| `points/<well>-points-B01` … `B06` | 6 | Individual decoded reads with x/y. B01 holds guides and GAPDH/ACTB control probes; B02–B06 hold the 3' transcriptome, sequenced in 5 batches |
| `tables/` | shared | Per-cell tables, below |

Image channels:
- **Cell paint:** `Mitochondria`, `Cell-Membrane`, `Nucleus`, `Actin`.
- **Protein:** `Phospho-ERK1_Phospho-ERK2`, `ERK1_ERK2`, `Phospho-HSP27`, `Phospho-p38`, plus `Mouse_IgG` and `Rabbit_IgG` as negative controls.
- **Background:** four channels (`Bkg_*`).

### Tables

The four cell tables have the same rows, in the same order, and identical `obs`. `region_key` (e.g. `A1-labels-Cell`) plus `Cell` link each row to its mask.

| Table | `X` | Notes |
|---|---|---|
| `global-table-ThreePrimeUntargeted` | RNA counts, 54,949 genes (GENCODE, 20k protein-coding) | Nuclear-only counts in `layers["nuclear"]`. Per-batch totals in `obsm["total_counts"]` |
| `global-table-imaging` | 1,150 features | Per-channel intensity: whole cell, nuclear, cytoplasmic, above-threshold. Also about 900 CellProfiler-style texture, granularity and intensity features. Shape features are in `obsm["cell_geometry"]` |
| `global-table-B01-Cellecta` | Guide read counts, 28 guides | Element's per-cell call in `obsm["target_assignment"]` |
| `global-table-B01-QC_Panel` | Targeted GAPDH/ACTB probes | Technical control |
| `global-table-metadata` | Empty | Run- and well-level stats in `uns` (`wells`, `batches`, `c2s_runstats_*`, …) |

Every cell table shares these `obs` columns. The first five are original; the rest were added by `prep.py`.

| Column | Meaning |
|---|---|
| `Well`, `WellLabel`, `Cell`, `region_key`, `tile_id` | Original Element identifiers |
| `condition` | `baseline` or `TNFa` |
| `guide`, `target_gene` | Most abundant guide per cell and its target. Non-targeting guides are labelled `NTC`; cells without a call are `Undetermined` |
| `guide_confidence` | Share of the cell's guide reads that come from its top guide |
| `guide_assigned` | `guide_confidence > 0.5`. Reproduces Element's 87% assignment rate |
| `rna_total_counts`, `rna_n_genes`, `rna_pct_mt` | RNA QC metrics |
| `qc_low_rna` | Fewer than 200 RNA counts (0.7% of cells) |
| `qc_guide_unassigned` | Not `guide_assigned` (12.9%) |
| `qc_spatial_artifact` | Element's per-cell artifact call in any modality (7.1%) |
| `qc_imaging_nan` | NaN in any imaging feature other than `*_NuclearRatio` (0.3%) |
| `qc_no_nucleus` | No segmented nucleus (0.2%) |
| `pass_qc` | No `qc_*` flag set (82.2%) |

**Nothing has been filtered.** Use `pass_qc`, or build your own filter from the flags.

## Gotchas

- `*_NuclearRatio` imaging features are NaN in up to 14% of otherwise good cells. Drop or impute them.
- `Bkg_*` features measure background, not biology. Leave them out of embeddings, or use them for correction.
- Guide calls are probabilistic. If an analysis depends on clean labels, use a stricter `guide_confidence` threshold.
- Wells differ in depth and intensity. Normalise per well, or include `Well` as a covariate.
- Guide abundance is uneven: OR1I1 has 10× more cells than TGFBR2.

## Changes from the original

- Added the `obs` columns above to every cell table.
- Element's well metadata lists `cell_type = HELA`, but the cells are A549. We set `uns["wells"]["cell_type"] = "A549"` and kept the original in `cell_type_original`. The expression data support A549:
  - **Male:** Y-chromosome genes are expressed (RPS4Y1, EIF1AY, DDX3Y) and XIST is absent. HeLa is female.
  - **KEAP1-mutant NRF2 program:** very high NQO1, AKR1C1/2/3, AKR1B10, ALDH3A1 and ALDH1A1.
  - **CDKN2A absent:** deleted in A549. HeLa expresses p16 strongly because of HPV E7.

## Reproduce

```bash
uv sync
# download and unzip Cytos-MPS_TNFa-CRISPR-Screen.zarr.zip (111 GB) from Element's dataset page
python prep.py Cytos-MPS_TNFa-CRISPR-Screen.zarr out    # annotates in place, writes tables + subset tiers
python checks/biology.py out/tables.zarr                 # sanity checks above
pytest test_prep.py
```

Element's own processing tools: [Elembio/cytoprofiling](https://github.com/Elembio/cytoprofiling).
