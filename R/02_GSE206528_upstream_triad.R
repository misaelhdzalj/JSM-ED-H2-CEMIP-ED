# 02 — GSE206528 prespecified upstream triad
# Exactly three frozen panels; no substitutions/additions.

run_gse206528_upstream_triad <- function(counts_file, schema_file, outdir) {
  safe_dir(outdir)
  obj <- read_gse206528(counts_file, schema_file)
  mat <- obj$mat; ec <- obj$ec; sc <- obj$sc; idx <- obj$idx

  panels <- list(
    HYPOXIA=c("ADM","ANGPTL4","EGLN3","PGF","VEGFA"),
    INFLAMMATION_TNF_NFKB=c("ICAM1","VCAM1","SELE","CCL2","TNFAIP3"),
    GLYCOLYSIS=c("SLC2A1","HK2","PFKFB3","PDK1","LDHA")
  )
  required <- c("CEMIP", unlist(panels, use.names=FALSE))
  if (!all(required %in% rownames(mat))) stop("At least one frozen triad gene is absent.")

  score_from_mat <- function(mm) {
    z <- tmm_full_ec(mm)
    ex <- z$cpm
    score <- sapply(panels, function(genes) {
      # t(ex[genes,]) => donors x genes; scale() z-standardizes each gene across donors.
      rowMeans(scale(t(ex[genes,,drop=FALSE]), center=TRUE, scale=TRUE))
    })
    if (is.null(dim(score))) score <- matrix(score, ncol=1L)
    rownames(score) <- colnames(mm)
    list(cpm=ex, score=score, cemip=as.numeric(ex["CEMIP",]))
  }

  full <- score_from_mat(mat)
  score <- full$score
  cemip <- full$cemip

  primary <- do.call(rbind, lapply(colnames(score), function(axis) {
    q <- spearman_exact_one_sided(cemip, score[,axis])
    data.frame(
      axis=axis, rho=unname(q["rho"]),
      permutation_hits=as.integer(q["hits"]),
      permutation_total=as.integer(q["total"]),
      exact_one_sided_p=unname(q["p_one_sided"]),
      stringsAsFactors=FALSE
    )
  }))
  primary$Holm_p_3primary <- p.adjust(primary$exact_one_sided_p, method="holm")
  write_tsv(primary, file.path(outdir, "TABLE_GSE206528_UPSTREAM_TRIAD_PRIMARY.tsv"))

  donor_scores <- data.frame(
    donor_id=ec[[sc$donor]],
    group=ec[[sc$group]],
    CEMIP_TMM_CPM=cemip,
    score,
    check.names=FALSE,
    stringsAsFactors=FALSE
  )
  write_tsv(donor_scores, file.path(outdir, "SUPPLEMENT_GSE206528_UPSTREAM_TRIAD_DONOR_SCORES.tsv"))

  # Disease-facing descriptive context retained solely to audit historical Phase05N.
  disease_context <- do.call(rbind, lapply(colnames(score), function(axis) {
    s <- score[,axis]
    perm <- exact_group_permutation(s, idx$control)
    hg <- hedges_g_independent(s[idx$ed], s[idx$control])
    data.frame(
      axis=axis,
      ED_minus_control_score_delta=unname(perm$observed),
      hedges_g=unname(hg["hedges_g"]),
      exact_ED_higher_hits=perm$hits,
      exact_ED_higher_total=perm$total,
      exact_ED_higher_p=perm$p_one_sided,
      stringsAsFactors=FALSE
    )
  }))
  write_tsv(disease_context, file.path(outdir, "SUPPLEMENT_GSE206528_UPSTREAM_TRIAD_DISEASE_CONTEXT.tsv"))

  # Mandatory LOO: omit donor, re-estimate TMM, recompute z-scores, score and rho.
  loo <- do.call(rbind, lapply(seq_len(ncol(mat)), function(i) {
    keep <- setdiff(seq_len(ncol(mat)),i)
    tmp <- score_from_mat(mat[,keep,drop=FALSE])
    do.call(rbind, lapply(colnames(tmp$score), function(axis) {
      rho <- suppressWarnings(cor(tmp$cemip, tmp$score[,axis], method="spearman"))
      data.frame(
        omitted_donor=ec[[sc$donor]][i],
        omitted_group=ec[[sc$group]][i],
        axis=axis,
        rho=rho,
        positive=(rho>0),
        stringsAsFactors=FALSE
      )
    }))
  }))
  write_tsv(loo, file.path(outdir, "SUPPLEMENT_GSE206528_UPSTREAM_TRIAD_LOO.tsv"))

  # Mandatory non-DMED sensitivity: n=6, re-estimate TMM and scores.
  k6 <- c(idx$control, idx$nondm)
  tmp6 <- score_from_mat(mat[,k6,drop=FALSE])
  sens <- do.call(rbind, lapply(colnames(tmp6$score), function(axis) {
    q <- spearman_exact_one_sided(tmp6$cemip, tmp6$score[,axis])
    data.frame(
      axis=axis, n=6L, rho=unname(q["rho"]),
      permutation_hits=as.integer(q["hits"]),
      permutation_total=as.integer(q["total"]),
      exact_one_sided_p=unname(q["p_one_sided"]),
      stringsAsFactors=FALSE
    )
  }))
  write_tsv(sens, file.path(outdir, "TABLE_GSE206528_UPSTREAM_TRIAD_NONDM_SENSITIVITY.tsv"))

  # Final classification under the previously frozen directional rubric.
  loo_pos <- aggregate(positive ~ axis, data=loo, FUN=sum)
  names(loo_pos)[2] <- "LOO_positive_n"
  merged <- merge(primary, disease_context, by="axis", sort=FALSE)
  merged <- merge(merged, loo_pos, by="axis", sort=FALSE)
  merged <- merge(merged, sens[,c("axis","rho","exact_one_sided_p")], by="axis",
                  suffixes=c("_primary","_nonDM"), sort=FALSE)
  merged$classification <- ifelse(
    merged$rho_primary <= 0, "C — NOT SUPPORTIVE",
    ifelse(merged$LOO_positive_n >= 6 &
             merged$ED_minus_control_score_delta > 0 &
             merged$rho_nonDM > 0,
           "A — SUPPORTIVE", "B — MIXED")
  )
  write_tsv(merged, file.path(outdir, "TABLE_GSE206528_UPSTREAM_TRIAD_FINAL_SUMMARY.tsv"))

  invisible(list(primary=primary, donor_scores=donor_scores, disease_context=disease_context,
                 loo=loo, nonDM=sens, final=merged))
}
