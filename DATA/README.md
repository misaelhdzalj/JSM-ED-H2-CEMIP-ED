# Frozen analysis inputs

## Included derived GSE206528 inputs

The repository includes the two canonical frozen pseudobulk inputs used for the GSE206528 analyses:

- `derived/GSE206528_target_pseudobulk_counts_v1.0.tsv.gz`
  - SHA-256: `98b0c7d37506999a4c0c39e702f399a3141069afa2fbb7bfbba7823b6d3f82fc`
  - size: 1,043,706 bytes

- `derived/GSE206528_target_pseudobulk_schema_v1.1.tsv`
  - SHA-256: `dd43cbed41f6fc9523c0ea61f820af220ab545c906ca410ed4463e2dd0cc43f0`
  - size: 4,430 bytes

These are frozen derived inputs from the publicly available GSE206528 source dataset.

## External GSE144803 source input

The original source file used by the frozen orthogonal TNF analysis is:

- `GSE144803_non-normalized.txt.gz`
  - SHA-256: `c1a8dacd6ad93574299f30886d5832758bdd3291af3f2e26de3f3d16925d7417`
  - expected size: 23,351,928 bytes

Because this file is publicly deposited in GEO, it is not duplicated in this repository. Obtain it from GSE144803 and verify the checksum before re-execution.
