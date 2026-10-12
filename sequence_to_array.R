#!/usr/bin/env Rscript
# Emulate a 450K / EPICv1 / EPICv2 microarray from a WGBS bsseq object (hg38),
# by filtering the bsseq CpGs down to the CpGs profiled by the chosen array.
# Optionally estimates epigenetic age (methylCIPHER core clocks) and pace of
# aging (DunedinPACE) from the emulated array data.
#
# Usage:
#   Rscript sequence_to_array.R my_bsseq_object.rda [--array EPICv1]
#           [--runClocks TRUE] [--methylOut emulated_data.csv]
#           [--ageOut age_estimates.csv] [--covMin 5]
#
# Coordinate handling: array probe positions are taken from hg38 Illumina
# manifests. For 450K/EPICv1 the default annotation columns are legacy hg19,
# so the embedded CHR_hg38/Start_hg38 columns are used instead (Start_hg38 is
# 0-based, so +1 gives the 1-based CpG position); EPICv2's hg38 manifest is
# hg38-native in its default chr/pos columns already. The bsseq object's own
# coordinate convention (0- vs 1-based) is auto-detected by comparing overlap
# counts against the array probes as-is vs. shifted by +1, and the better-
# matching convention is used.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--array", type = "character", default = "EPICv1",
              help = "Array to emulate: 450K, EPICv1, or EPICv2 [default %default]"),
  make_option("--runClocks", type = "character", default = "TRUE",
              help = "Run methylCIPHER + DunedinPACE epigenetic age/pace estimates [default %default]"),
  make_option("--methylOut", type = "character", default = "emulated_data.csv",
              help = "Output CSV path for the emulated array beta matrix [default %default]"),
  make_option("--ageOut", type = "character", default = "age_estimates.csv",
              help = "Output CSV path for epigenetic age/pace estimates [default %default]"),
  make_option("--covMin", type = "integer", default = 5,
              help = "Minimum read coverage to call a beta value; below this the value is NA [default %default]")
)

parser <- OptionParser(usage = "%prog <bsseq_object.rda> [options]", option_list = option_list)
parsed <- parse_args(parser, positional_arguments = 1)
opt <- parsed$options
bsseq_path <- parsed$args[1]

array <- switch(toupper(opt$array),
                 "450K" = "450K", "450" = "450K",
                 "EPICV1" = "EPICv1", "EPIC" = "EPICv1", "EPIC1" = "EPICv1",
                 "EPICV2" = "EPICv2", "EPIC2" = "EPICv2",
                 stop("Unknown --array '", opt$array, "'; must be one of 450K, EPICv1, EPICv2"))
runClocks <- as.logical(opt$runClocks)
if (is.na(runClocks)) stop("--runClocks must be TRUE or FALSE")
COV_MIN <- opt$covMin

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(minfi)
})

## ---- 1. Load the bsseq object ---------------------------------------------

cat("Loading bsseq object from", bsseq_path, "...\n")
e <- new.env()
loaded <- load(bsseq_path, envir = e)
is_bsseq <- vapply(loaded, function(n) is(get(n, envir = e), "BSseq"), logical(1))
if (!any(is_bsseq)) stop("No BSseq object found in ", bsseq_path)
obj_name <- loaded[is_bsseq][1]
bs <- get(obj_name, envir = e)
rm(e)
cat("Loaded object '", obj_name, "', dim: ", paste(dim(bs), collapse = " x "),
    " (CpGs x samples)\n", sep = "")

## ---- 2. Build the hg38 array probe panel ----------------------------------

cat("Building", array, "hg38 probe panel...\n")
if (array %in% c("450K", "EPICv1")) {
  suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38))
  anno <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38))
  if (array == "450K") anno <- anno[!is.na(anno$Methyl450_Loci) & anno$Methyl450_Loci == "TRUE", ]
  anno$Start_hg38 <- as.numeric(anno$Start_hg38)
  anno <- anno[!is.na(anno$Start_hg38) & !is.na(anno$CHR_hg38) & anno$CHR_hg38 != "", ]
  probe_gr <- GRanges(anno$CHR_hg38, IRanges(anno$Start_hg38 + 1, width = 1))
  names(probe_gr) <- rownames(anno)
} else {
  suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
  anno <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
  anno$pos <- as.numeric(anno$pos)
  anno <- anno[!is.na(anno$pos) & !is.na(anno$chr) & anno$chr != "", ]
  probe_gr <- GRanges(anno$chr, IRanges(anno$pos, width = 1))
  names(probe_gr) <- rownames(anno)
}
cat(array, "probe panel size:", length(probe_gr), "\n")

## ---- 3. Align chromosome naming, then detect 0- vs 1-based bsseq coords ---

bs_gr <- granges(bs)

if (length(intersect(seqlevels(bs_gr), seqlevels(probe_gr))) == 0) {
  if (all(grepl("^chr", seqlevels(probe_gr))) && !any(grepl("^chr", seqlevels(bs_gr)))) {
    seqlevels(bs_gr) <- paste0("chr", seqlevels(bs_gr))
  } else if (!any(grepl("^chr", seqlevels(probe_gr))) && all(grepl("^chr", seqlevels(bs_gr)))) {
    seqlevels(bs_gr) <- sub("^chr", "", seqlevels(bs_gr))
  }
}

ov_asis  <- findOverlaps(probe_gr, bs_gr, type = "equal")
shift_gr <- GRanges(seqnames(bs_gr), IRanges(start(bs_gr) + 1, width = 1))
ov_shift <- findOverlaps(probe_gr, shift_gr, type = "equal")

n_asis  <- length(unique(queryHits(ov_asis)))
n_shift <- length(unique(queryHits(ov_shift)))
cat(sprintf("Coordinate check: %d/%d probes match as-is; %d/%d match with a +1 shift.\n",
            n_asis, length(probe_gr), n_shift, length(probe_gr)))

if (n_shift > n_asis) {
  cat("-> bsseq object appears 0-based; shifting +1 to align with 1-based array probe positions.\n")
  ov <- ov_shift
  coord_mode <- "0-based (shifted +1)"
} else {
  cat("-> bsseq object appears standard 1-based; using positions as-is.\n")
  ov <- ov_asis
  coord_mode <- "1-based"
}

ov_df <- as.data.frame(ov)
ov_df <- ov_df[!duplicated(ov_df$queryHits), ]
cat(sprintf("Matched %d / %d %s probes to CpGs in the bsseq object (%s).\n",
            nrow(ov_df), length(probe_gr), array, coord_mode))
if (nrow(ov_df) == 0) {
  stop("No probes matched. Check that the bsseq object is hg38 and that chromosome naming ",
       "(e.g. 'chr1' vs '1') is compatible with the array manifest.")
}

## ---- 4. Extract the emulated array beta matrix (CpGs x samples) ----------

cat("Extracting beta values (min coverage =", COV_MIN, ")...\n")
M   <- as.matrix(getCoverage(bs, type = "M")[ov_df$subjectHits, , drop = FALSE])
Cov <- as.matrix(getCoverage(bs, type = "Cov")[ov_df$subjectHits, , drop = FALSE])
beta <- M / Cov
beta[Cov < COV_MIN] <- NA
rownames(beta) <- names(probe_gr)[ov_df$queryHits]

sample_ids <- sampleNames(bs)
if (is.null(sample_ids) || anyNA(sample_ids) || any(sample_ids == "") || any(duplicated(sample_ids))) {
  cat("Note: bsseq object has no valid/unique sampleNames(); using generic sample IDs (sample1, sample2, ...).\n")
  sample_ids <- sprintf("sample%d", seq_len(ncol(beta)))
}
colnames(beta) <- sample_ids
beta <- beta[order(rownames(beta)), , drop = FALSE]

cat(sprintf("Emulated %s array: %d CpGs x %d samples. Fraction NA (Cov<%d): %.2f%%\n",
            array, nrow(beta), ncol(beta), COV_MIN, 100 * mean(is.na(beta))))

cat("Writing emulated array data to", opt$methylOut, "...\n")
write.csv(beta, opt$methylOut, row.names = TRUE)

## ---- 5. Optionally run epigenetic clocks ----------------------------------

if (runClocks) {
  cat("\n=== Running DunedinPACE ===\n")
  suppressPackageStartupMessages(library(DunedinPACE))
  pace_res <- PACEProjector(beta, proportionOfProbesRequired = 0.8)
  pace <- pace_res$DunedinPACE
  cat("DunedinPACE computed for", sum(!is.na(pace)), "/", ncol(beta), "samples.\n")

  cat("\n=== Running all methylCIPHER clocks ===\n")
  cat("Please remember to cite the individual clocks you use; see methylCIPHER's README.md.\n")
  suppressPackageStartupMessages(library(methylCIPHER))
  DNAm <- as.data.frame(t(beta))

  clock_fns <- c(
    "calcHannum", "calcHorvath1", "calcHorvath2", "calcPhenoAge", "calcEpiTOC",
    "calcEpiTOC2", "calcDunedinPoAm38", "calcAlcoholMcCartney", "calcBMIMcCartney",
    "calcBocklandt", "calcBohlin", "calcDNAmClockCortical", "calcDNAmTL",
    "calcGaragnani", "calcHRSInChPhenoAge", "calcHypoClock", "calcKnight",
    "calcLeeControl", "calcLeeRefinedRobust", "calcLeeRobust", "calcLin",
    "calcMayne", "calcMiAge", "calcPEDBE", "calcSmokingMcCartney",
    "calcVidalBralo", "calcWeidner", "calcZhang", "calcZhang2019"
  )

  pheno <- data.frame(row.names = rownames(DNAm))
  for (fn in clock_fns) {
    pheno <- tryCatch(get(fn)(DNAm, pheno, imputation = FALSE),
                       error = function(e) {
                         cat("  [skipped]", fn, "-", conditionMessage(e), "\n")
                         pheno
                       })
  }

  age_est <- pheno
  age_est$DunedinPACE <- pace[rownames(age_est)]
  age_est <- data.frame(SampleID = rownames(age_est), age_est, row.names = NULL,
                         check.names = FALSE, stringsAsFactors = FALSE)

  pd <- tryCatch(pData(bs), error = function(e) NULL)
  if (!is.null(pd)) {
    pd_df <- as.data.frame(pd)
    if (ncol(pd_df) > 0) {
      cat("Appending pData() (", ncol(pd_df), " columns) to age estimates.\n", sep = "")
      rownames(pd_df) <- sample_ids
      age_est <- cbind(age_est, pd_df[age_est$SampleID, , drop = FALSE])
    } else {
      cat("pData() has 0 columns; skipping.\n")
    }
  } else {
    cat("pData() is NULL; skipping.\n")
  }

  cat("Writing age/pace estimates to", opt$ageOut, "...\n")
  write.csv(age_est, opt$ageOut, row.names = FALSE)
}

cat("\nDONE.\n")
