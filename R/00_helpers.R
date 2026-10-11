# JSM-ED-H2 FINAL R FREEZE — shared helpers
# Scope: donor-level EC pseudobulk only for GSE206528; paired targeted HUVEC analysis for GSE144803.

write_tsv <- function(x, file) {
  write.table(x, file=file, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
}

safe_dir <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive=TRUE, showWarnings=FALSE)
  normalizePath(x, winslash="/", mustWork=FALSE)
}

sha256_file <- function(path) {
  if (!requireNamespace("tools", quietly=TRUE)) stop("Base R tools package unavailable.")
  # tools::md5sum is not SHA256; SHA256 is therefore verified outside R via provenance manifest.
  # This helper intentionally does not silently substitute MD5.
  NA_character_
}

detect_schema_columns <- function(schema) {
  pick <- function(candidates) {
    hit <- intersect(candidates, names(schema))
    if (length(hit) != 1L) stop("Could not uniquely resolve schema column among: ", paste(candidates, collapse=", "))
    hit[[1]]
  }
  list(
    profile=pick(c("pseudobulk_id","profile","sample","column")),
    donor=pick(c("donor_id","donor")),
    cell=pick(c("cell_type","celltype","celltype_final")),
    group=pick(c("group","disease_group","condition"))
  )
}

read_gse206528 <- function(counts_file, schema_file) {
  counts <- read.delim(gzfile(counts_file), check.names=FALSE, stringsAsFactors=FALSE)
  if (ncol(counts) < 2L) stop("Counts file has fewer than 2 columns.")
  gene_col <- names(counts)[1]
  if (anyDuplicated(counts[[gene_col]])) stop("Duplicate gene symbols detected.")
  rownames(counts) <- counts[[gene_col]]
  counts[[gene_col]] <- NULL

  schema <- read.delim(schema_file, check.names=FALSE, stringsAsFactors=FALSE)
  sc <- detect_schema_columns(schema)
  ec <- schema[toupper(schema[[sc$cell]]) %in% c("EC","ENDOTHELIAL","ENDOTHELIAL CELLS"), , drop=FALSE]
  if (nrow(ec) != 8L) stop("Expected exactly 8 EC donor pseudobulks; found ", nrow(ec))
  if (!all(ec[[sc$profile]] %in% colnames(counts))) stop("At least one EC profile is absent from counts matrix.")

  # Preserve schema order; it is donor-defined and not outcome-optimized.
  mat <- as.matrix(counts[, ec[[sc$profile]], drop=FALSE])
  storage.mode(mat) <- "integer"
  colnames(mat) <- ec[[sc$profile]]

  # Frozen group map.
  g <- tolower(ec[[sc$group]])
  ctl <- which(g %in% c("normal","control","ctrl"))
  nondm <- which(g %in% c("non_diabetic_ed","non-diabetic_ed","non-diabetic ed","non_diabetic ed"))
  dmed <- which(g %in% c("diabetic_ed","diabetic-ed","diabetic ed","dmed"))
  ed <- sort(c(nondm, dmed))
  if (length(ctl)!=3L || length(nondm)!=3L || length(dmed)!=2L || length(ed)!=5L) {
    stop("Frozen group sizes not recovered: expected Control=3, non-DMED=3, DMED=2.")
  }

  list(counts=counts, schema=schema, sc=sc, ec=ec, mat=mat,
       idx=list(control=ctl, nondm=nondm, dmed=dmed, ed=ed))
}

tmm_full_ec <- function(mat) {
  if (!requireNamespace("edgeR", quietly=TRUE)) stop("edgeR is required.")
  y <- edgeR::DGEList(counts=mat)
  y <- edgeR::calcNormFactors(y, method="TMM")
  cpm <- edgeR::cpm(y, normalized.lib.sizes=TRUE, log=FALSE)
  list(y=y, cpm=cpm,
       norm_factors=y$samples$norm.factors,
       library_size=y$samples$lib.size,
       effective_library_size=y$samples$lib.size*y$samples$norm.factors)
}

exact_group_permutation <- function(values, control_idx) {
  n <- length(values)
  k <- length(control_idx)
  combos <- combn(seq_len(n), k)
  stat_fun <- function(cidx) mean(values[-cidx]) - mean(values[cidx])
  obs <- stat_fun(control_idx)
  perm <- apply(combos, 2L, stat_fun)
  tol <- 1e-12
  hits <- sum(perm >= obs - tol)
  list(observed=obs, total=ncol(combos), hits=hits, p_one_sided=hits/ncol(combos),
       perm_stats=perm)
}

all_permutations <- function(v) {
  if (length(v) == 1L) return(matrix(v, nrow=1L))
  do.call(rbind, lapply(seq_along(v), function(i) {
    cbind(v[i], all_permutations(v[-i]))
  }))
}

spearman_exact_one_sided <- function(x, y) {
  if (length(x) != length(y)) stop("x and y length mismatch.")
  n <- length(x)
  if (n > 9L) stop("Exact rank enumeration intentionally limited to n<=9.")
  rx <- rank(x, ties.method="average")
  ry <- rank(y, ties.method="average")
  obs <- suppressWarnings(cor(rx, ry, method="pearson"))
  perms <- all_permutations(seq_len(n))
  null_rho <- apply(perms, 1L, function(idx) suppressWarnings(cor(rx, ry[idx], method="pearson")))
  tol <- 1e-12
  hits <- sum(null_rho >= obs - tol)
  c(rho=obs, hits=hits, total=nrow(perms), p_one_sided=hits/nrow(perms))
}

hedges_g_independent <- function(x_treatment, x_control, conf.level=0.95) {
  x1 <- as.numeric(x_treatment); x0 <- as.numeric(x_control)
  n1 <- length(x1); n0 <- length(x0)
  if (n1 < 2L || n0 < 2L) stop("Need >=2 observations per group.")
  df <- n1+n0-2
  sp <- sqrt(((n1-1)*var(x1)+(n0-1)*var(x0))/df)
  d <- (mean(x1)-mean(x0))/sp
  J <- 1 - 3/(4*df - 1)
  g <- J*d
  se_scale <- sqrt(1/n1 + 1/n0)
  t_obs <- d/se_scale

  alpha <- 1-conf.level
  root_ncp <- function(target) {
    f <- function(ncp) stats::pt(t_obs, df=df, ncp=ncp) - target
    lo <- -100; hi <- 100
    flo <- f(lo); fhi <- f(hi)
    if (is.na(flo) || is.na(fhi) || flo*fhi > 0) stop("Could not bracket noncentral-t CI root.")
    stats::uniroot(f, interval=c(lo,hi), tol=1e-12)$root
  }
  # For a fixed positive/negative t, CDF decreases monotonically with ncp.
  ncp_lo <- root_ncp(1-alpha/2)
  ncp_hi <- root_ncp(alpha/2)
  d_lo <- ncp_lo*se_scale
  d_hi <- ncp_hi*se_scale
  g_ci <- sort(J*c(d_lo,d_hi))

  c(n_treatment=n1, n_control=n0, mean_treatment=mean(x1), mean_control=mean(x0),
    cohen_d=d, hedges_g=g, ci95_low=g_ci[1], ci95_high=g_ci[2], df=df, t_equiv=t_obs)
}

paired_hedges_dz <- function(delta) {
  delta <- as.numeric(delta)
  n <- length(delta)
  dz <- mean(delta)/sd(delta)
  # Small-sample correction for standardized mean change based on n paired differences.
  J <- 1 - 3/(4*(n-1)-1)
  c(n=n, cohen_dz=dz, hedges_gz=J*dz, correction_J=J)
}

mean_ci_t <- function(x, conf.level=0.95) {
  x <- as.numeric(x); n <- length(x); a <- 1-conf.level
  se <- sd(x)/sqrt(n)
  crit <- qt(1-a/2, df=n-1)
  c(mean=mean(x), low=mean(x)-crit*se, high=mean(x)+crit*se)
}

classify_numeric <- function(final, expected, material=FALSE, exact_tol=1e-12,
                             equiv_abs=1e-8, equiv_rel=1e-6) {
  if (is.na(final) || is.na(expected)) return("MATERIAL DISCREPANCY")
  if (abs(final-expected) <= exact_tol) return("EXACT MATCH")
  if (abs(final-expected) <= max(equiv_abs, equiv_rel*max(1,abs(expected)))) return("NUMERICALLY EQUIVALENT")
  if (material) return("MATERIAL DISCREPANCY")
  "MINOR DISCREPANCY"
}

classify_logical <- function(final, expected) {
  if (isTRUE(identical(as.logical(final), as.logical(expected)))) "EXACT MATCH" else "MATERIAL DISCREPANCY"
}
