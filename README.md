# JSM-ED-H2 — CEMIP in erectile dysfunction

Reproducibility materials for:

**CEMIP expression in human cavernous endothelium: an exploratory transcriptomic study of erectile dysfunction**

Authors: Misael Alejandro Juárez Hernández, Yuli Hernández González, Ricardo Alonso Castillejos Molina.

## Repository status

**GITHUB-READY RELEASE CANDIDATE — keep private until the final repository/Zenodo release step.**

The frozen analysis code, manuscript-facing derived results, runtime records, and the two canonical derived GSE206528 pseudobulk inputs are included. The GSE144803 source file is not duplicated because it is publicly available from NCBI GEO; its expected SHA-256 is recorded for exact input verification.

## Data sources

The study is a secondary analysis of publicly available transcriptomic datasets:

- **GSE206528** — human corpus cavernosum single-cell transcriptomics; donor-level endothelial pseudobulk was used for the primary ED analysis.
- **GSE144803** — paired primary HUVEC expression before and after TNF-alpha exposure; used as orthogonal endothelial perturbational evidence.

Original source datasets remain available from NCBI GEO under their respective accession numbers. This repository is a reproducibility resource and does not replace GEO.

## Canonical analysis

For GSE206528, the manuscript-facing genome-wide endothelial differential-expression analysis uses edgeR quasi-likelihood with donor as the biological replicate, design-aware `filterByExpr`, TMM normalization, robust dispersion estimation, quasi-likelihood inference, and Benjamini-Hochberg correction across retained genes.

For GSE144803, the frozen targeted analysis evaluates probe `ILMN_1813704` (KIAA1199/CEMIP) in 39 paired HUVEC cultures.

## Repository structure

- `R/` — frozen R analysis scripts.
- `RESULTS/` — manuscript-facing frozen derived tables and supplementary result tables.
- `RUNTIME/` — R/session/package-version records.
- `METHODS/` — frozen statistical-method description and manuscript-facing numerical values.
- `PROVENANCE/` — input SHA-256 identifiers, repository manifest/checksums, and method-finalization provenance.
- `AUDIT/` — final audit records retained for reproducibility/provenance.
- `DATA/derived/` — the two canonical frozen GSE206528 pseudobulk inputs required for re-execution.
- `DATA/README.md` — source/input instructions, including the external GSE144803 file.

## Required inputs for re-execution

The frozen scripts require:

1. `DATA/derived/GSE206528_target_pseudobulk_counts_v1.0.tsv.gz`
2. `DATA/derived/GSE206528_target_pseudobulk_schema_v1.1.tsv`
3. `GSE144803_non-normalized.txt.gz`

The first two are included in this repository. The third is an original publicly deposited GEO file and should be obtained from GSE144803 rather than duplicated here. Exact expected SHA-256 values are recorded in `PROVENANCE/INPUT_SHA256SUMS.txt`.

## Re-execution

From the repository root, after obtaining the GSE144803 source file:

```bash
Rscript R/run_all.R \
  DATA/derived/GSE206528_target_pseudobulk_counts_v1.0.tsv.gz \
  DATA/derived/GSE206528_target_pseudobulk_schema_v1.1.tsv \
  /path/to/GSE144803_non-normalized.txt.gz \
  FINAL_R_RUN
```

A faithful rerun should return `READY_FOR_MANUSCRIPT_DRAFTING` and reproduce the frozen manuscript-facing outputs.

## Software environment

The frozen execution record used R 4.6.1 and edgeR 4.10.5. See `RUNTIME/sessionInfo.txt` and `RUNTIME/PACKAGE_VERSIONS.tsv` for the recorded environment.

## Input integrity

The included GSE206528 derived inputs have been verified against the frozen expected checksums:

- `GSE206528_target_pseudobulk_counts_v1.0.tsv.gz` — SHA-256 `98b0c7d37506999a4c0c39e702f399a3141069afa2fbb7bfbba7823b6d3f82fc`
- `GSE206528_target_pseudobulk_schema_v1.1.tsv` — SHA-256 `dd43cbed41f6fc9523c0ea61f820af220ab545c906ca410ed4463e2dd0cc43f0`

## Final release step

Before making the GitHub repository public or publishing the Zenodo deposit:

- create/reserve the final Zenodo DOI;
- add the final repository/DOI citation to the release metadata and manuscript placeholders;
- perform one last public-file audit;
- then make the GitHub repository public and publish the archived Zenodo release.
