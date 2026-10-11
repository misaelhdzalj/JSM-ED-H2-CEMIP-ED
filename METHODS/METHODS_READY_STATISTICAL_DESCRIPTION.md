# Methods-ready statistical description

## GSE206528: human erectile-dysfunction endothelial pseudobulk

The biological replicate is the **donor**. Single cells are not treated as independent replicates. The primary analysis uses the eight endothelial-cell (EC) donor pseudobulks from GSE206528 (3 controls, 3 non-diabetic ED, 2 diabetic ED).

For donor-level CEMIP display and robustness analyses, TMM normalization is estimated across the complete EC pseudobulk count matrix using **edgeR**. CEMIP is summarized as TMM-normalized counts per million (CPM) and, because all eight donor values are positive, as log2(TMM-CPM) without a pseudocount. Individual donor values, the proportion of ED donors exceeding the maximum control value, an exact donor-label permutation test, Hedges g with a 95% confidence interval obtained by inversion of the noncentral t distribution, a leave-one-donor-out analysis with TMM re-estimated after each omission, and a prespecified non-diabetic ED sensitivity analysis are reported. The exact ED-versus-control permutation enumerates all C(8,3)=56 allocations; the non-diabetic ED sensitivity enumerates all C(6,3)=20 allocations. The permutation statistic is the ED-minus-control difference in mean donor log2(TMM-CPM). These candidate-level statistics are post-selection robustness measures and do not replace genome-wide multiplicity control.

The final reproducible genome-wide exploratory differential-expression analysis is performed once under a frozen **edgeR quasi-likelihood** workflow. Genes are filtered with `filterByExpr` using the 3-versus-5 donor design, retained libraries are TMM-normalized, dispersions are estimated with robust empirical-Bayes estimation, and a quasi-likelihood F test is fit for ED versus control. The complete gene universe is retained in the output table, with filtered genes marked and inferential columns left missing. Benjamini-Hochberg false-discovery rates are calculated over the tested genes. If this final R rerun differs numerically from historical Phase-01 values, the final R values are canonical and the discrepancy is documented rather than selecting the method that yields a more favorable p value.

Diabetic ED (n=2) is shown descriptively and is not used for a separate small-subgroup inferential claim.

## Prespecified upstream triad

Three fixed EC expression scores are used without gene substitution: hypoxia (ADM, ANGPTL4, EGLN3, PGF, VEGFA), inflammation/TNF–NF-κB (ICAM1, VCAM1, SELE, CCL2, TNFAIP3), and glycolysis (SLC2A1, HK2, PFKFB3, PDK1, LDHA). TMM-normalized donor-level expression is z-standardized separately for each gene across eligible donors, and each score is the equal-weight mean of the five gene z scores.

Association between donor-level CEMIP expression and each score is quantified using Spearman rho. Exact one-sided positive-association p values are obtained by enumeration of all 8! rank permutations, and the three primary p values are corrected using Holm's method. Robustness analyses include leave-one-donor-out recalculation in which TMM normalization and score standardization are repeated after each omission, and a non-diabetic sensitivity analysis restricted to the three controls and three non-diabetic ED donors with TMM and scores re-estimated in that six-donor subset. These analyses assess transcriptomic compatibility with an upstream state and do not establish cytokine exposure, pathway activation, or causal induction of CEMIP.

## GSE144803: paired TNF-α HUVEC perturbation

A targeted analysis uses probe **ILMN_1813704 (KIAA1199/CEMIP)** from the deposited non-normalized Illumina HumanHT-12 v4 dataset for the 39 paired HUVEC cultures retained after source-study quality control. For each pair, untreated and TNF-α raw signals and detection p values are extracted. Detectability is defined as detection p<0.05.

The analysis reports individual paired values, the number of pairs in which TNF-α signal exceeds untreated signal, detectability transitions, the paired log2 raw-signal change, its median and mean, a two-sided exact Wilcoxon signed-rank test, an exact two-sided sign test, and an exact paired detectability test based on discordant pairs. Paired standardized effect size is expressed as Cohen dz and its small-sample Hedges-corrected counterpart. This is a targeted pre-filter interrogation of the deposited probe and is not presented as part of the source study's filtered genome-wide differential-expression analysis.
