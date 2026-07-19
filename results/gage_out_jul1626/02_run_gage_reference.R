#!/usr/bin/env Rscript
# Run the REAL gage R package on its own demo data (gse16873 + kegg.gs) and
# export (a) the prepared per-gene fold-change matrix, (b) the gene sets, and
# (c) gage's greater/less result tables for t-test, z-test, and Fisher meta.
# pygage will be fed the IDENTICAL prepared matrix so the comparison isolates
# the gene-set test + meta step (no prep differences).
options(stringsAsFactors = FALSE)
outdir <- "gage_out"; dir.create(outdir, showWarnings = FALSE)

# --- obtain gage: installed package preferred, else source the cloned repo ---
if (requireNamespace("gage", quietly = TRUE)) {
  library(gage)
  data(gse16873, package = "gage"); data(kegg.gs, package = "gage")
  message("Using installed gage package: ", as.character(packageVersion("gage")))
} else if (dir.exists("gage_R")) {
  for (f in c("gage-internal.R","invalid.R","odd.R","gagePrep.R",
              "gs.tTest.R","gs.zTest.R","gs.KSTest.R","gageSum.R","gage.R"))
    source(file.path("gage_R","R",f))
  load("gage_R/data/gse16873.rda"); load("gage_R/data/kegg.gs.rda")
  message("Using gage source from ./gage_R (core path: base R + stats only)")
} else stop("Need either the gage package or a ./gage_R clone (run 01_setup.sh)")

# --- the standard vignette comparison: 6 HN vs 6 DCIS, paired ---
cn   <- colnames(gse16873)
hn   <- grep("HN",   cn, ignore.case = TRUE)
dcis <- grep("DCIS", cn, ignore.case = TRUE)
message(sprintf("gse16873: %d genes x %d samples | |HN|=%d |DCIS|=%d",
                nrow(gse16873), ncol(gse16873), length(hn), length(dcis)))

# prepared per-gene fold-change matrix (this is the shared input for pygage)
prep <- gagePrep(gse16873, ref = hn, samp = dcis, compare = "paired", use.fold = TRUE)
write.csv(data.frame(gene_id = rownames(prep), prep, check.names = FALSE),
          file.path(outdir, "prepared_matrix.csv"), row.names = FALSE)

# gage runs: default (t-test + Stouffer), z-test, and Fisher/gamma meta
res_t      <- gage(gse16873, gsets = kegg.gs, ref = hn, samp = dcis, compare = "paired")
res_z      <- gage(gse16873, gsets = kegg.gs, ref = hn, samp = dcis, compare = "paired", saaTest = gs.zTest)
res_fisher <- gage(gse16873, gsets = kegg.gs, ref = hn, samp = dcis, compare = "paired", use.stouffer = FALSE)

wr <- function(x, f) write.csv(data.frame(gene_set = rownames(x), x, check.names = FALSE),
                               file.path(outdir, f), row.names = FALSE)
wr(res_t$greater,      "gage_tTest_greater.csv")
wr(res_t$less,         "gage_tTest_less.csv")
wr(res_z$greater,      "gage_zTest_greater.csv")
wr(res_fisher$greater, "gage_fisher_greater.csv")

# gene sets -> JSON (so pygage reads the exact same set memberships)
con <- file(file.path(outdir, "kegg_gs.json"), "w"); cat("{", file = con); nm <- names(kegg.gs)
for (i in seq_along(kegg.gs)) {
  cat(paste0('"', gsub('"','\\\\"', nm[i]), '":[',
             paste0('"', kegg.gs[[i]], '"', collapse = ","), ']',
             if (i < length(kegg.gs)) "," else ""), file = con)
}
cat("}", file = con); close(con)

message("Wrote gage outputs to ", normalizePath(outdir))
