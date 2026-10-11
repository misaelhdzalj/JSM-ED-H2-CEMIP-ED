# 03 — GSE144803 targeted paired HUVEC TNF-alpha analysis for ILMN_1813704 (KIAA1199/CEMIP)

run_gse144803_tnf_cemip <- function(raw_file, outdir) {
  safe_dir(outdir)
  x <- read.delim(gzfile(raw_file), check.names=FALSE, stringsAsFactors=FALSE)
  probe_col <- intersect(c("PROBE_ID","ProbeID","ID","ID_REF"), names(x))
  if (length(probe_col)!=1L) stop("Could not uniquely identify probe column.")
  row <- x[x[[probe_col]]=="ILMN_1813704",,drop=FALSE]
  if (nrow(row)!=1L) stop("ILMN_1813704 not found uniquely.")

  # 39 source-study retained paired HUVEC donors; donor 35 is not retained.
  donors <- c(0:34,36:39)

  find_col <- function(d, condition, measure) {
    nm <- names(row)
    z <- tolower(gsub("[[:space:]_.-]+","",nm))
    dpat <- paste0("huvec",d)
    cond_pat <- if (tolower(condition)=="untreated") "untreated" else "tnf"
    meas_pat <- if (tolower(measure)=="signal") "avgsignal" else "detectionpval"
    # Anchor donor+condition at the start so donor 1 cannot match donor 10/11/etc.
    prefix <- paste0("^", dpat, cond_pat)
    hit <- which(grepl(prefix,z) & grepl(meas_pat,z,fixed=TRUE))
    if (length(hit)!=1L) stop("Could not uniquely map donor ",d," ",condition," ",measure,"; hits=",length(hit))
    nm[hit]
  }

  out <- do.call(rbind, lapply(donors, function(d) {
    u_sig <- find_col(d,"untreated","signal")
    u_p <- find_col(d,"untreated","detection")
    t_sig <- find_col(d,"tnf","signal")
    t_p <- find_col(d,"tnf","detection")
    data.frame(
      donor=d,
      untreated_raw_signal=as.numeric(row[[u_sig]]),
      tnf_raw_signal=as.numeric(row[[t_sig]]),
      untreated_detection_p=as.numeric(row[[u_p]]),
      tnf_detection_p=as.numeric(row[[t_p]]),
      stringsAsFactors=FALSE
    )
  }))
  if (nrow(out)!=39L) stop("Expected 39 paired donors.")

  out$TNF_gt_untreated <- out$tnf_raw_signal > out$untreated_raw_signal
  out$untreated_detected_p_lt_0_05 <- out$untreated_detection_p < 0.05
  out$tnf_detected_p_lt_0_05 <- out$tnf_detection_p < 0.05
  out$detectability_transition <- ifelse(
    !out$untreated_detected_p_lt_0_05 & out$tnf_detected_p_lt_0_05, "not_detected_to_detected",
    ifelse(out$untreated_detected_p_lt_0_05 & !out$tnf_detected_p_lt_0_05, "detected_to_not_detected",
           ifelse(out$untreated_detected_p_lt_0_05 & out$tnf_detected_p_lt_0_05, "detected_both","not_detected_both"))
  )
  out$log2_untreated_raw <- log2(out$untreated_raw_signal)
  out$log2_tnf_raw <- log2(out$tnf_raw_signal)
  out$delta_log2_TNF_minus_untreated <- out$log2_tnf_raw - out$log2_untreated_raw
  out$fold_change_TNF_over_untreated <- out$tnf_raw_signal/out$untreated_raw_signal
  write_tsv(out, file.path(outdir, "SUPPLEMENT_GSE144803_CEMIP_39_PAIRED_VALUES.tsv"))

  delta <- out$delta_log2_TNF_minus_untreated
  # Exact signed-rank: no zero differences and no ties in absolute differences in this dataset.
  wil <- wilcox.test(delta, mu=0, alternative="two.sided", exact=TRUE, correct=FALSE)
  sgn <- binom.test(sum(delta>0), length(delta), p=0.5, alternative="two.sided")
  b <- sum(!out$untreated_detected_p_lt_0_05 & out$tnf_detected_p_lt_0_05)
  c <- sum(out$untreated_detected_p_lt_0_05 & !out$tnf_detected_p_lt_0_05)
  det_p <- if ((b+c)>0L) binom.test(b, b+c, p=0.5, alternative="two.sided")$p.value else NA_real_
  eff <- paired_hedges_dz(delta)
  ci <- mean_ci_t(delta)

  summ <- data.frame(
    probe="ILMN_1813704",
    annotation="KIAA1199/CEMIP",
    n_pairs=nrow(out),
    TNF_gt_untreated_n=sum(out$TNF_gt_untreated),
    TNF_gt_untreated_den=nrow(out),
    untreated_detected_n=sum(out$untreated_detected_p_lt_0_05),
    TNF_detected_n=sum(out$tnf_detected_p_lt_0_05),
    not_detected_to_detected_n=b,
    detected_to_not_detected_n=c,
    median_untreated_raw_signal=median(out$untreated_raw_signal),
    median_TNF_raw_signal=median(out$tnf_raw_signal),
    median_delta_log2=median(delta),
    mean_delta_log2=mean(delta),
    mean_delta_log2_ci95_low=unname(ci["low"]),
    mean_delta_log2_ci95_high=unname(ci["high"]),
    geometric_mean_fold_change=2^mean(delta),
    paired_cohen_dz=unname(eff["cohen_dz"]),
    paired_hedges_gz=unname(eff["hedges_gz"]),
    wilcoxon_signed_rank_two_sided_p=wil$p.value,
    exact_sign_test_two_sided_p=sgn$p.value,
    exact_paired_detectability_two_sided_p=det_p,
    stringsAsFactors=FALSE
  )
  write_tsv(summ, file.path(outdir, "TABLE_GSE144803_CEMIP_TNF_FINAL_SUMMARY.tsv"))

  invisible(list(pairs=out, summary=summ))
}
