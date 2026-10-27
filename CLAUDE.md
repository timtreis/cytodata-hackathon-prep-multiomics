# AI Agent Instructions

Data prep for the CytoData hackathon at SBI2: Element Biosciences Cytos multiomics (cell paint + 3' transcriptome + protein + CRISPR guide, same cells), QC'd and packaged for participants.

## Dataset

- Cytos Multimodal Pooled Screening, A549 CRISPR screen, ±TNFα, ~258k cells.
- Source: `https://element-public-data.s3.us-west-2.amazonaws.com/Cytos/elembio-multiomics/Cytos-MPS_TNFa-CRISPR-Screen/0c47c400756647deb9b3a37f502458de/spatialdata/Cytos-MPS_TNFa-CRISPR-Screen.zarr.zip` (111 GB, SpatialData zarr).
- Processing reference: https://github.com/Elembio/cytoprofiling

## Rules

- Terse output. Lead with results. Say plainly when something failed.
- Data never goes in git. Too big for laptops, so it lives on the cluster.
- Cluster hostnames, paths, accounts and login details live only in `CLAUDE.local.md` (gitignored). Never write them into tracked files, commits, PRs or issues.
- Use scverse tooling (spatialdata, anndata, scanpy). No new dependencies without asking.
- Never commit or push without being asked. Stage only files you touched.
- `tasks/` is local scratch and gitignored.
