# Manuscript-facing final values

Only the values below (or values directly contained in the frozen final tables) should be used for manuscript drafting.

## GSE206528 — genome-wide EC pseudobulk differential expression

Canonical method: edgeR quasi-likelihood; donor as replicate; 3 controls vs 5 ED; filterByExpr; TMM; robust dispersion/QL; BH correction over tested genes.

- genes input: 33,694
- genes tested: 14,947
- CEMIP log2FC, ED − control: 5.11848760800753
- CEMIP nominal p: 0.000161192568249643
- CEMIP BH-FDR: 0.0964187076198762
- genes BH-FDR <0.10: 26
- minimum BH-FDR: 0.0856666804660712
- genes BH-FDR <0.05: 0

Interpretation: CEMIP meets an exploratory FDR<0.10 threshold but not the conventional FDR<0.05 threshold. Do not call it conventionally genome-wide significant.

## GSE206528 — donor-level CEMIP robustness

- mean log2 TMM-CPM difference, ED − control: 5.33627412695182
- Cohen d: 4.01575251883641
- Hedges g: 3.49195871203166
- Hedges g 95% CI: 1.16519985396976 to 5.74853169250463
- exact one-sided donor-label permutation: 1/56 = 0.0178571428571429
- ED > maximum control: 5/5
- all 8 leave-one-donor-out effects positive: yes

### Non-diabetic ED sensitivity

- mean log2 TMM-CPM difference: 5.92176783000785
- Hedges g: 3.47780661679928
- Hedges g 95% CI: 0.80876723443636 to 6.09097123363009
- exact one-sided permutation: 1/20 = 0.05
- non-DMED > maximum control: 3/3

## Frozen upstream triad

### Hypoxia
- Spearman rho: 0.166666666666667
- exact one-sided p: 0.351661706349206
- Holm p: 0.703323412698413
- LOO positive: 6/8
- non-DMED rho: -0.0285714285714286

### Inflammation/TNF-NF-kB
- Spearman rho: 0.666666666666667
- exact one-sided p: 0.0415426587301587
- Holm p: 0.124627976190476
- LOO positive: 8/8
- non-DMED rho: 0.771428571428571
- non-DMED exact one-sided p: 0.0513888888888889

Interpretation: supportive exploratory association; not significant after Holm correction.

### Glycolysis
- Spearman rho: -0.142857142857143
- exact one-sided p: 0.648338293650794
- Holm p: 0.703323412698413
- LOO positive: 1/8
- non-DMED rho: -0.142857142857143

## GSE144803 — targeted paired TNF HUVEC perturbation

- probe: ILMN_1813704 (KIAA1199/CEMIP)
- pairs: 39
- TNF > untreated: 30/39
- untreated detected: 0/39
- TNF detected: 23/39
- not-detected → detected: 23/39
- median delta log2: 0.983179438696014
- mean delta log2: 1.58829493611596
- mean delta log2 95% CI: 1.02935628764448 to 2.14723358458744
- geometric mean fold-change: 3.00693761399727
- paired Hedges gz: 0.902847816401271
- exact Wilcoxon paired two-sided p: 8.54459358379246e-07
- exact sign-test two-sided p: 0.0010650196345523
- exact paired detectability p: 2.38418579101562e-07
