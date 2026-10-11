# JSM-ED-H2 — FINAL R AUDIT REPORT

## Final disposition

**Unresolved material discrepancy: NO**

# READY FOR MANUSCRIPT DRAFTING

The canonical R execution completed and all manuscript-facing analyses are frozen. The initial hard stop was caused by four genome-wide differences between the recovered historical Phase-01 discovery implementation and the canonical final edgeR-QL implementation.

After recovery of the original Phase-01 and Phase-02 records, those four differences are classified as:

**RESOLVED METHOD-FINALIZATION DIFFERENCE**

They are not hidden, deleted, or reclassified as exact matches. The original stop report and numerical audit are preserved under `AUDIT/ORIGINAL_CANONICAL_RUN_STOP/`, and the resolved comparison is preserved in `AUDIT/HISTORICAL_COMPARISON_RESOLVED.tsv`.

## Resolution basis

The recovered Phase-01 record explicitly required a later canonical DESeq2/edgeR rerun before manuscript submission, and Phase-02 independently repeated that requirement. The final edgeR-QL execution therefore fulfills the documented method-finalization step. No alternative DE engines were compared after seeing the final p/FDR values.

## Manuscript-facing rule

For genome-wide GSE206528 inference, cite only the frozen canonical edgeR-QL estimates in `RESULTS/TABLE_GSE206528_GENOMEWIDE_DE_SUMMARY.tsv` and the complete table in `RESULTS/SUPPLEMENT_GSE206528_GENOMEWIDE_EDGER_QL.tsv`.

The earlier Phase-01 numerical values remain provenance records and are not manuscript-facing inferential estimates.

## Other analyses

Donor-level CEMIP robustness, exact permutations, leave-one-donor-out analyses, the prespecified upstream triad, and GSE144803 paired TNF analysis reproduced the expected structures and conclusions. Small TMM-related numerical shifts relative to historical tables are preserved in the historical comparison and do not change inferential conclusions.
