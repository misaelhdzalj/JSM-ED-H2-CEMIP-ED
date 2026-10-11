# 04 — Final reproducibility audit against the canonical manuscript freeze.
# Historical Phase-01 method-finalization differences are preserved separately in
# AUDIT/HISTORICAL_COMPARISON_RESOLVED.tsv and do not constitute unresolved discrepancies.

run_final_audit <- function(outdir, audit_outdir) {
  safe_dir(audit_outdir)
  read1 <- function(f) read.delim(file.path(outdir,f), check.names=FALSE, stringsAsFactors=FALSE)

  p <- read1("TABLE_GSE206528_CEMIP_PRIMARY_ROBUSTNESS.tsv")
  n <- read1("TABLE_GSE206528_CEMIP_NONDM_SENSITIVITY.tsv")
  de <- read1("TABLE_GSE206528_GENOMEWIDE_DE_SUMMARY.tsv")
  tri <- read1("TABLE_GSE206528_UPSTREAM_TRIAD_FINAL_SUMMARY.tsv")
  tnf <- read1("TABLE_GSE144803_CEMIP_TNF_FINAL_SUMMARY.tsv")
  loo <- read1("SUPPLEMENT_GSE206528_CEMIP_LOO_REESTIMATED_TMM.tsv")

  rows <- list()
  add_num <- function(section, metric, final, expected, material_rule=FALSE, note="") {
    rows[[length(rows)+1L]] <<- data.frame(
      section=section, metric=metric, final_value=as.character(signif(final,15)),
      frozen_expected=as.character(signif(expected,15)),
      classification=classify_numeric(final, expected, material=material_rule),
      note=note, stringsAsFactors=FALSE
    )
  }
  add_log <- function(section, metric, final, expected, note="") {
    rows[[length(rows)+1L]] <<- data.frame(
      section=section, metric=metric, final_value=as.character(final),
      frozen_expected=as.character(expected),
      classification=classify_logical(final, expected),
      note=note, stringsAsFactors=FALSE
    )
  }

  # Canonical final expected values from the executed R 4.6.1 / edgeR 4.10.5 freeze.
  add_num("GSE206528 CEMIP donor robustness","delta mean log2",p$delta_mean_log2_ED_minus_control,5.33627412695182,FALSE)
  add_num("GSE206528 CEMIP donor robustness","Hedges g",p$hedges_g,3.49195871203166,FALSE)
  add_num("GSE206528 CEMIP donor robustness","Hedges g CI low",p$hedges_g_ci95_low,1.16519985396976,FALSE)
  add_num("GSE206528 CEMIP donor robustness","Hedges g CI high",p$hedges_g_ci95_high,5.74853169250463,FALSE)
  add_num("GSE206528 CEMIP donor robustness","exact permutation total",p$permutation_total,56,TRUE)
  add_num("GSE206528 CEMIP donor robustness","exact permutation p",p$exact_one_sided_permutation_p,1/56,FALSE)
  add_log("GSE206528 CEMIP donor robustness","5/5 ED above maximum control",p$all_ED_above_max_control,TRUE)

  add_num("GSE206528 CEMIP non-DMED","exact permutation total",n$permutation_total,20,TRUE)
  add_num("GSE206528 CEMIP non-DMED","exact permutation p",n$exact_one_sided_permutation_p,0.05,FALSE)
  add_log("GSE206528 CEMIP non-DMED","3/3 non-DMED above maximum control",n$all_nonDM_above_max_control,TRUE)

  add_num("GSE206528 genome-wide DE","CEMIP log2FC",de$CEMIP_log2FC_ED_minus_control,5.11848760800753,TRUE)
  add_num("GSE206528 genome-wide DE","CEMIP nominal p",de$CEMIP_nominal_p,0.000161192568249643,TRUE)
  add_num("GSE206528 genome-wide DE","CEMIP BH-FDR",de$CEMIP_BH_FDR,0.0964187076198762,TRUE)
  add_num("GSE206528 genome-wide DE","number genes BH-FDR <0.10",de$n_genes_FDR_lt_0_10,26,TRUE)
  add_num("GSE206528 genome-wide DE","minimum BH-FDR",de$min_BH_FDR,0.0856666804660712,TRUE)

  add_log("GSE206528 CEMIP LOO","all 8 LOO effects positive",all(loo$positive_direction),TRUE)

  exp_tri <- data.frame(
    axis=c("HYPOXIA","INFLAMMATION_TNF_NFKB","GLYCOLYSIS"),
    rho=c(0.166666666666667,0.666666666666667,-0.142857142857143),
    p=c(0.351661706349206,0.0415426587301587,0.648338293650794),
    holm=c(0.703323412698413,0.124627976190476,0.703323412698413),
    loo=c(6,8,1),
    nondm=c(-0.0285714285714286,0.771428571428571,-0.142857142857143),
    stringsAsFactors=FALSE
  )
  for (i in seq_len(nrow(exp_tri))) {
    a <- exp_tri$axis[i]; f <- tri[tri$axis==a,,drop=FALSE]
    add_num("GSE206528 upstream triad",paste(a,"rho"),f$rho_primary,exp_tri$rho[i],TRUE)
    add_num("GSE206528 upstream triad",paste(a,"exact p"),f$exact_one_sided_p_primary,exp_tri$p[i],TRUE)
    add_num("GSE206528 upstream triad",paste(a,"Holm p"),f$Holm_p_3primary,exp_tri$holm[i],TRUE)
    add_num("GSE206528 upstream triad",paste(a,"LOO positive n"),f$LOO_positive_n,exp_tri$loo[i],TRUE)
    add_num("GSE206528 upstream triad",paste(a,"non-DMED rho"),f$rho_nonDM,exp_tri$nondm[i],TRUE)
  }

  add_num("GSE144803 TNF","n pairs",tnf$n_pairs,39,TRUE)
  add_num("GSE144803 TNF","TNF > untreated pairs",tnf$TNF_gt_untreated_n,30,TRUE)
  add_num("GSE144803 TNF","untreated detected",tnf$untreated_detected_n,0,TRUE)
  add_num("GSE144803 TNF","TNF detected",tnf$TNF_detected_n,23,TRUE)
  add_num("GSE144803 TNF","median delta log2",tnf$median_delta_log2,0.983179438696014,TRUE)
  add_num("GSE144803 TNF","mean delta log2",tnf$mean_delta_log2,1.58829493611596,TRUE)
  add_num("GSE144803 TNF","Wilcoxon two-sided p",tnf$wilcoxon_signed_rank_two_sided_p,8.54459358379246e-07,TRUE)
  add_num("GSE144803 TNF","exact sign-test p",tnf$exact_sign_test_two_sided_p,0.0010650196345523,TRUE)
  add_num("GSE144803 TNF","paired detectability exact p",tnf$exact_paired_detectability_two_sided_p,2.38418579101562e-07,TRUE)
  add_num("GSE144803 TNF","paired Hedges gz",tnf$paired_hedges_gz,0.902847816401271,TRUE)

  audit <- do.call(rbind,rows)
  write_tsv(audit, file.path(audit_outdir,"FINAL_REPRODUCIBILITY_AUDIT.tsv"))

  material <- any(audit$classification=="MATERIAL DISCREPANCY")
  lines <- c(
    "# JSM-ED-H2 — FINAL REPRODUCIBILITY AUDIT",
    "",
    paste0("Unresolved material discrepancy: **",ifelse(material,"YES","NO"),"**"),
    "",
    if (material) "**STOP:** the rerun does not reproduce the frozen manuscript values."
    else "# READY FOR MANUSCRIPT DRAFTING",
    "",
    "Historical Phase-01 method-finalization differences are provenance records and are not used as the reproducibility target for the final canonical edgeR-QL model."
  )
  writeLines(lines,file.path(audit_outdir,"FINAL_AUDIT_REPORT.md"))
  invisible(list(audit=audit,material=material))
}
