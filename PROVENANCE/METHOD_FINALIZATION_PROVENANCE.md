# Method-finalization provenance resolution

## Why the first final-R audit stopped

The first frozen edgeR-QL execution produced genome-wide values that differed from numerical values inherited from the earlier Phase-01 discovery implementation. Under the original audit rules, this correctly triggered a material-discrepancy hard stop.

## Evidence recovered after the stop

The original packages were subsequently recovered:

- `JSM_ED_H2_PHASE01_DISCOVERY_FINAL_v1.0.zip`
- `JSM_ED_H2_PHASE02_CEMIP_TRIAGE_FINAL_v1.0.zip`

The Phase-01 report documents that its discovery ranking used a donor-level count-model approximation because R/DESeq2 was not available in that execution environment, and explicitly states that before manuscript submission the selected analysis should be rerun once in a canonical DESeq2/edgeR environment using the same donor-level counts.

The Phase-02 decision independently reiterates that before manuscript submission the human count model should be reproduced in DESeq2 or edgeR using the same donor-level EC pseudobulk and disclosed as the definitive count-model implementation.

Therefore the canonical edgeR-QL analysis in this package is a planned method-finalization step, not a post-result switch to a more favorable method.

## Historical versus final values

Historical Phase-01 discovery implementation, CEMIP EC:

- log2FC: approximately 4.9231
- nominal p: approximately 0.00061136
- expression-filtered BH-FDR: approximately 0.96633
- no genes at FDR<0.10 under that discovery implementation

Final canonical edgeR-QL implementation:

- log2FC: 5.11848760800753
- nominal p: 0.000161192568249643
- BH-FDR: 0.0964187076198762
- 26 genes at FDR<0.10
- 0 genes at FDR<0.05

Both are preserved in the project record. Only the canonical final edgeR-QL values are manuscript-facing inferential estimates.

## Audit consequence

The four former `MATERIAL DISCREPANCY` rows are reclassified as **RESOLVED METHOD-FINALIZATION DIFFERENCE** because their provenance is now established. They are not deleted or relabeled as exact matches. The original stop-generating audit is retained verbatim under `AUDIT/ORIGINAL_CANONICAL_RUN_STOP/`.
