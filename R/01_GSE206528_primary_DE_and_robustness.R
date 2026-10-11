# 01 — GSE206528 final donor-level EC genome-wide DE + CEMIP robustness
# Frozen inferential engine: edgeR quasi-likelihood, donor-level pseudobulk, EC only.
# No cell is treated as an independent replicate.

run_gse206528_primary <- function(counts_file, schema_file, outdir) {
  safe_dir(outdir)
  obj <- read_gse206528(counts_file, schema_file)
  mat <- obj$mat; ec <- obj$ec; sc <- obj$sc; idx <- obj$idx

  # ---------------------------------------------------------------------------
  # A. Descriptive/candidate-facing donor normalization.
  # TMM is estimated on the COMPLETE 33,694-gene EC pseudobulk matrix.
  # This is kept separate from the design-aware filtered genome-wide DE model.
  # ---------------------------------------------------------------------------
  full <- tmm_full_ec(mat)
  if (!("CEMIP" %in% rownames(full$cpm))) stop("CEMIP not present.")
  cemip_cpm <- as.numeric(full$cpm["CEMIP",])
  if (any(cemip_cpm <= 0)) stop("CEMIP contains non-positive donor CPM; frozen no-pseudocount log2 convention cannot be used.")
  cemip_log2 <- log2(cemip_cpm)

  donor <- data.frame(
    pseudobulk_id=ec[[sc$profile]],
    donor_id=ec[[sc$donor]],
    group=ec[[sc$group]],
    broad_group=ifelse(seq_len(nrow(ec)) %in% idx$control, "Control", "ED"),
    raw_CEMIP=as.integer(mat["CEMIP",]),
    raw_library_size=colSums(mat),
    TMM_norm_factor=full$norm_factors,
    effective_library_size=full$effective_library_size,
    CEMIP_TMM_CPM=cemip_cpm,
    log2_CEMIP_TMM_CPM=cemip_log2,
    stringsAsFactors=FALSE
  )
  write_tsv(donor, file.path(outdir, "TABLE_GSE206528_CEMIP_DONOR_VALUES.tsv"))

  # Frozen 3-vs-5 robustness.
  hg <- hedges_g_independent(cemip_log2[idx$ed], cemip_log2[idx$control])
  perm <- exact_group_permutation(cemip_log2, idx$control)
  sep_all <- all(cemip_cpm[idx$ed] > max(cemip_cpm[idx$control]))

  primary <- data.frame(
    endpoint="CEMIP donor-level descriptive robustness: ED vs control",
    n_control=length(idx$control),
    n_ED=length(idx$ed),
    delta_mean_log2_ED_minus_control=unname(perm$observed),
    cohen_d=unname(hg["cohen_d"]),
    hedges_g=unname(hg["hedges_g"]),
    hedges_g_ci95_low=unname(hg["ci95_low"]),
    hedges_g_ci95_high=unname(hg["ci95_high"]),
    permutation_hits=perm$hits,
    permutation_total=perm$total,
    exact_one_sided_permutation_p=perm$p_one_sided,
    ED_above_max_control_n=sum(cemip_cpm[idx$ed] > max(cemip_cpm[idx$control])),
    ED_above_max_control_den=length(idx$ed),
    all_ED_above_max_control=sep_all,
    stringsAsFactors=FALSE
  )
  write_tsv(primary, file.path(outdir, "TABLE_GSE206528_CEMIP_PRIMARY_ROBUSTNESS.tsv"))

  # Non-DMED frozen sensitivity, re-estimating TMM within 3 control + 3 non-DMED donors.
  k6 <- c(idx$control, idx$nondm)
  mat6 <- mat[,k6,drop=FALSE]
  full6 <- tmm_full_ec(mat6)
  c6 <- as.numeric(full6$cpm["CEMIP",])
  l6 <- log2(c6)
  local_ctl <- seq_along(k6)[k6 %in% idx$control]
  local_ed <- seq_along(k6)[k6 %in% idx$nondm]
  perm6 <- exact_group_permutation(l6, local_ctl)
  hg6 <- hedges_g_independent(l6[local_ed], l6[local_ctl])
  nondm <- data.frame(
    endpoint="CEMIP non-DMED sensitivity: non-diabetic ED vs control",
    n_control=3L, n_nonDM_ED=3L,
    delta_mean_log2_nonDM_minus_control=unname(perm6$observed),
    cohen_d=unname(hg6["cohen_d"]),
    hedges_g=unname(hg6["hedges_g"]),
    hedges_g_ci95_low=unname(hg6["ci95_low"]),
    hedges_g_ci95_high=unname(hg6["ci95_high"]),
    permutation_hits=perm6$hits,
    permutation_total=perm6$total,
    exact_one_sided_permutation_p=perm6$p_one_sided,
    nonDM_above_max_control_n=sum(c6[local_ed] > max(c6[local_ctl])),
    nonDM_above_max_control_den=length(local_ed),
    all_nonDM_above_max_control=all(c6[local_ed] > max(c6[local_ctl])),
    stringsAsFactors=FALSE
  )
  write_tsv(nondm, file.path(outdir, "TABLE_GSE206528_CEMIP_NONDM_SENSITIVITY.tsv"))

  # Descriptive 3-level condition table. DMED n=2 remains descriptive only.
  lev <- list(Control=idx$control, NonDiabetic_ED=idx$nondm, DMED=idx$dmed)
  group_desc <- do.call(rbind, lapply(names(lev), function(nm) {
    z <- lev[[nm]]
    data.frame(
      group=nm, n=length(z),
      mean_TMM_CPM=mean(cemip_cpm[z]),
      median_TMM_CPM=median(cemip_cpm[z]),
      mean_log2_TMM_CPM=mean(cemip_log2[z]),
      median_log2_TMM_CPM=median(cemip_log2[z]),
      min_TMM_CPM=min(cemip_cpm[z]), max_TMM_CPM=max(cemip_cpm[z])
    )
  }))
  write_tsv(group_desc, file.path(outdir, "TABLE_GSE206528_CEMIP_THREE_GROUP_DESCRIPTIVE.tsv"))

  # LOO: omit one donor, RE-ESTIMATE TMM from complete remaining EC matrix.
  loo <- do.call(rbind, lapply(seq_len(ncol(mat)), function(i) {
    keep <- setdiff(seq_len(ncol(mat)), i)
    sub <- tmm_full_ec(mat[,keep,drop=FALSE])
    v <- as.numeric(sub$cpm["CEMIP",])
    l <- log2(v)
    control_local <- which(keep %in% idx$control)
    ed_local <- which(keep %in% idx$ed)
    data.frame(
      omitted_donor=ec[[sc$donor]][i],
      omitted_group=ec[[sc$group]][i],
      n_remaining=length(keep),
      delta_mean_log2_ED_minus_control=mean(l[ed_local])-mean(l[control_local]),
      positive_direction=(mean(l[ed_local])-mean(l[control_local])>0),
      stringsAsFactors=FALSE
    )
  }))
  write_tsv(loo, file.path(outdir, "SUPPLEMENT_GSE206528_CEMIP_LOO_REESTIMATED_TMM.tsv"))

  # ---------------------------------------------------------------------------
  # B. Complete genome-wide exploratory EC differential-expression rerun.
  # Frozen engine: edgeR QL. This supersedes historical Phase-01 p/FDR if numeric
  # values differ. No alternative method is tested/selected by significance.
  # ---------------------------------------------------------------------------
  condition <- factor(ifelse(seq_len(ncol(mat)) %in% idx$control, "Control", "ED"),
                      levels=c("Control","ED"))
  design <- model.matrix(~condition)
  y0 <- edgeR::DGEList(counts=mat, group=condition)
  keep_gene <- edgeR::filterByExpr(y0, design=design)
  y <- y0[keep_gene,,keep.lib.sizes=FALSE]
  y <- edgeR::calcNormFactors(y, method="TMM")
  y <- edgeR::estimateDisp(y, design=design, robust=TRUE)
  fit <- edgeR::glmQLFit(y, design=design, robust=TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef="conditionED")
  tab <- edgeR::topTags(qlf, n=Inf, sort.by="none")$table

  allgenes <- data.frame(
    gene_symbol=rownames(mat),
    retained_by_filter=rownames(mat) %in% rownames(tab),
    logFC_ED_minus_control=NA_real_,
    logCPM=NA_real_,
    F=NA_real_,
    PValue=NA_real_,
    FDR_BH=NA_real_,
    stringsAsFactors=FALSE
  )
  m <- match(allgenes$gene_symbol, rownames(tab))
  ok <- !is.na(m)
  allgenes$logFC_ED_minus_control[ok] <- tab$logFC[m[ok]]
  if ("logCPM" %in% names(tab)) allgenes$logCPM[ok] <- tab$logCPM[m[ok]]
  if ("F" %in% names(tab)) allgenes$F[ok] <- tab$F[m[ok]]
  allgenes$PValue[ok] <- tab$PValue[m[ok]]
  allgenes$FDR_BH[ok] <- tab$FDR[m[ok]]
  write_tsv(allgenes, file.path(outdir, "SUPPLEMENT_GSE206528_GENOMEWIDE_EDGER_QL.tsv"))

  cemip_de <- allgenes[allgenes$gene_symbol=="CEMIP",,drop=FALSE]
  if (nrow(cemip_de)!=1L) stop("CEMIP genome-wide DE row not unique.")
  de_summary <- data.frame(
    method="edgeR quasi-likelihood; filterByExpr; TMM; robust dispersion/QL",
    n_donors=8L, n_control=3L, n_ED=5L,
    genes_input=nrow(mat),
    genes_retained=sum(keep_gene),
    CEMIP_retained=cemip_de$retained_by_filter,
    CEMIP_log2FC_ED_minus_control=cemip_de$logFC_ED_minus_control,
    CEMIP_nominal_p=cemip_de$PValue,
    CEMIP_BH_FDR=cemip_de$FDR_BH,
    n_genes_FDR_lt_0_10=sum(allgenes$FDR_BH < 0.10, na.rm=TRUE),
    min_BH_FDR=min(allgenes$FDR_BH, na.rm=TRUE),
    stringsAsFactors=FALSE
  )
  write_tsv(de_summary, file.path(outdir, "TABLE_GSE206528_GENOMEWIDE_DE_SUMMARY.tsv"))

  invisible(list(donor=donor, primary=primary, nondm=nondm, group_desc=group_desc,
                 loo=loo, allgenes=allgenes, de_summary=de_summary))
}
