# 05 — Archive exact R runtime and package versions.

write_session_info <- function(outdir) {
  safe_dir(outdir)
  capture.output(sessionInfo(), file=file.path(outdir,"sessionInfo.txt"))
  pkgs <- c("edgeR")
  pv <- do.call(rbind, lapply(pkgs, function(p) {
    data.frame(package=p,
               installed=requireNamespace(p,quietly=TRUE),
               version=if (requireNamespace(p,quietly=TRUE)) as.character(utils::packageVersion(p)) else NA_character_,
               stringsAsFactors=FALSE)
  }))
  base <- data.frame(package="R", installed=TRUE, version=R.version.string, stringsAsFactors=FALSE)
  write_tsv(rbind(base,pv), file.path(outdir,"PACKAGE_VERSIONS.tsv"))
  invisible(pv)
}
