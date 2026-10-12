#!/usr/bin/env Rscript
# Infer candidate allele-specific methylation (ASM) regions from Illumina
# array-style beta values (e.g. the output of sequence_to_array.R).
#
# Rationale: at a locus with ASM one allele is (near-)fully methylated and the
# other (near-)unmethylated, so the bulk beta value sits near 0.5. A single
# intermediate CpG is weak evidence (noise, cell mixture, low coverage), but a
# run of >= --minCpGs neighbouring CpGs that are all intermediate is a much
# stronger ("soft") indication of ASM. Imprinted control regions (ICRs) are
# the textbook example, so calls are validated against the Human Imprintome
# ICR set and a curated imprinted-gene list.
#
# Usage:
#   Rscript epic_infer_asm.R emulated_data.csv [--low 0.4] [--high 0.6]
#           [--minCpGs 3] [--maxGap 1000] [--maxMiss 0] [--minSamples 2]
#           [--array auto] [--outPrefix asm] [--knownRegions my_asm.bed]
#
# Input: CSV with probe IDs in the first column and one beta column per sample.
# If the CSV already has `chr` and `pos` columns (hg38, 1-based) those are
# used; otherwise probe IDs are mapped to hg38 with the Illumina annotation
# (450K/EPICv1: CHR_hg38 / Start_hg38 + 1, since the default chr/pos columns
# of the ilm10b5.hg38 package are legacy hg19; EPICv2: native hg38).
#
# Outputs:
#   <outPrefix>_sample_regions.csv     one row per ASM run per sample
#   <outPrefix>_consensus_regions.csv  per-sample runs merged across samples,
#                                      with sample support + ICR/gene overlap
#   <outPrefix>_probe_calls.csv        per-probe flags for probes in any run
#   <outPrefix>_validation.txt         enrichment / recall vs known imprinting

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--low", type = "double", default = 0.4,
              help = "Lower bound of the 'intermediate' beta window [default %default]"),
  make_option("--high", type = "double", default = 0.6,
              help = "Upper bound of the 'intermediate' beta window [default %default]"),
  make_option("--minCpGs", type = "integer", default = 3,
              help = "Minimum number of consecutive intermediate CpGs to call a region [default %default]"),
  make_option("--maxGap", type = "integer", default = 1000,
              help = "Max distance (bp) between neighbouring CpGs within a run [default %default]"),
  make_option("--maxMiss", type = "integer", default = 0,
              help = "Non-intermediate CpGs tolerated inside a run (never at its ends) [default %default]"),
  make_option("--skipNA", type = "character", default = "TRUE",
              help = "Skip NA probes rather than letting them break a run [default %default]"),
  make_option("--minSamples", type = "integer", default = 2,
              help = "Min samples supporting a consensus region for it to be flagged as 'consensus' [default %default]"),
  make_option("--dropSexChr", type = "character", default = "FALSE",
              help = "Drop chrX/chrY (X-inactivation gives ASM-like signal in females) [default %default]"),
  make_option("--array", type = "character", default = "auto",
              help = "Probe annotation: auto, 450K, EPICv1, EPICv2 [default %default]"),
  make_option("--imprintome", type = "character",
              default = "/media/data/WGBS/Genomes/Human_Genome/variants/blacklisted/imprinted_regions_hg38.rda",
              help = "RDA with a GRanges of known ICRs (hg38) for validation; 'none' to skip [default %default]"),
  make_option("--imprintedGenes", type = "character",
              default = "/media/data/WGBS/Genomes/Human_Genome/variants/blacklisted/imprinted_genes.txt",
              help = "TSV of imprinted genes (Chromosome/Start/End/Gene name); 'none' to skip [default %default]"),
  make_option("--knownRegions", type = "character", default = NULL,
              help = "Optional extra BED/CSV (chr,start,end[,name]) of known ASM regions (e.g. monocyte ASM), hg38"),
  make_option("--outPrefix", type = "character", default = "asm",
              help = "Prefix for output files [default %default]")
)

parser <- OptionParser(usage = "%prog <betas.csv> [options]", option_list = option_list)
parsed <- parse_args(parser, positional_arguments = 1)
opt <- parsed$options
in_csv <- parsed$args[1]
skipNA <- as.logical(opt$skipNA)
dropSex <- as.logical(opt$dropSexChr)
if (is.na(skipNA) || is.na(dropSex)) stop("--skipNA / --dropSexChr must be TRUE or FALSE")
if (!(opt$low < opt$high)) stop("--low must be < --high")
if (opt$minCpGs < 2) stop("--minCpGs must be >= 2")

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

## ---- 1. Read betas ---------------------------------------------------------

cat("Reading", in_csv, "...\n")
dt <- fread(in_csv, header = TRUE)
setnames(dt, 1, "probe")
if (dt$probe[1] == "" || is.na(dt$probe[1])) stop("First column must hold probe IDs")
has_coords <- all(c("chr", "pos") %in% names(dt))
sample_cols <- setdiff(names(dt), c("probe", "chr", "pos", "strand"))
cat(nrow(dt), "probes x", length(sample_cols), "samples\n")

## ---- 2. Probe coordinates (hg38) ------------------------------------------

if (!has_coords) {
  arr <- toupper(opt$array)
  if (arr == "AUTO") arr <- if (any(grepl("_(TC|BC|TO|BO)[0-9]+$", head(dt$probe, 1000)))) "EPICV2" else "EPICV1"
  cat("Annotating probes with", arr, "hg38 manifest...\n")
  if (arr %in% c("450K", "EPICV1", "EPIC")) {
    suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38))
    anno <- as.data.frame(minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38))
    anno <- data.table(probe = rownames(anno), chr = anno$CHR_hg38,
                       pos = suppressWarnings(as.numeric(anno$Start_hg38)) + 1)
  } else if (arr %in% c("EPICV2", "EPIC2")) {
    suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
    anno <- as.data.frame(minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
    anno <- data.table(probe = rownames(anno), chr = anno$chr, pos = as.numeric(anno$pos))
  } else stop("Unknown --array ", opt$array)
  anno <- anno[!is.na(pos) & !is.na(chr) & chr != ""]
  dt <- merge(dt, anno, by = "probe")
  cat(nrow(dt), "probes with hg38 coordinates\n")
}
dt <- dt[grepl("^chr([0-9]+|X|Y)$", chr)]
if (dropSex) dt <- dt[!chr %in% c("chrX", "chrY")]
chr_levels <- paste0("chr", c(1:22, "X", "Y"))
dt[, chr := factor(chr, levels = chr_levels)]
setorder(dt, chr, pos)
dt[, chr := as.character(chr)]
# Collapse duplicate positions (e.g. EPICv2 replicate probes) by mean beta
if (anyDuplicated(dt, by = c("chr", "pos"))) {
  dt <- dt[, c(list(probe = paste(probe, collapse = ";")),
               lapply(.SD, mean, na.rm = TRUE)), by = .(chr, pos), .SDcols = sample_cols]
  for (s in sample_cols) set(dt, which(is.nan(dt[[s]])), s, NA_real_)
}

## ---- 3. Run detection ------------------------------------------------------

# Walk probes in genomic order; a run grows while CpGs are intermediate,
# tolerating up to maxMiss non-intermediate CpGs (only kept if the run later
# resumes) and never spanning a gap > maxGap or a chromosome change.
call_runs <- function(chr, pos, beta, low, high, minCpGs, maxGap, maxMiss, skipNA) {
  keep <- if (skipNA) !is.na(beta) else rep(TRUE, length(beta))
  idx <- which(keep)
  chr <- chr[idx]; pos <- pos[idx]; b <- beta[idx]
  inter <- !is.na(b) & b >= low & b <= high
  n <- length(b)
  # break points: chromosome change or gap too large
  brk <- c(TRUE, chr[-1] != chr[-n] | (pos[-1] - pos[-n]) > maxGap)
  out <- list(); k <- 0L
  start <- NA_integer_; last_hit <- NA_integer_; miss <- 0L; nhit <- 0L
  flush <- function() {
    if (!is.na(start) && nhit >= minCpGs) {
      k <<- k + 1L
      out[[k]] <<- c(idx[start], idx[last_hit], nhit)
    }
  }
  for (i in seq_len(n)) {
    if (brk[i]) { flush(); start <- NA_integer_; miss <- 0L; nhit <- 0L }
    if (inter[i]) {
      if (is.na(start)) { start <- i; nhit <- 0L; miss <- 0L }
      last_hit <- i; nhit <- nhit + 1L; miss <- 0L
    } else if (!is.na(start)) {
      miss <- miss + 1L
      if (miss > maxMiss) { flush(); start <- NA_integer_; miss <- 0L; nhit <- 0L }
    }
  }
  flush()
  if (k == 0) return(NULL)
  m <- do.call(rbind, out)
  data.table(i_start = m[, 1], i_end = m[, 2], n_intermediate = m[, 3])
}

cat(sprintf("Calling runs: beta in [%.2f, %.2f], >= %d CpGs, gap <= %d bp, maxMiss %d\n",
            opt$low, opt$high, opt$minCpGs, opt$maxGap, opt$maxMiss))
sample_regions <- rbindlist(lapply(sample_cols, function(s) {
  b <- dt[[s]]
  r <- call_runs(dt$chr, dt$pos, b, opt$low, opt$high, opt$minCpGs,
                 opt$maxGap, opt$maxMiss, skipNA)
  if (is.null(r)) return(NULL)
  r[, `:=`(sample = s, chr = dt$chr[i_start], start = dt$pos[i_start], end = dt$pos[i_end])]
  r[, c("n_probes", "mean_beta", "mean_abs_dev_0.5", "probes") := {
    res <- lapply(seq_len(.N), function(j) {
      ii <- i_start[j]:i_end[j]
      bb <- b[ii]; bb <- bb[!is.na(bb)]
      list(length(bb), mean(bb), mean(abs(bb - 0.5)), paste(dt$probe[ii], collapse = ","))
    })
    lapply(1:4, function(q) sapply(res, `[[`, q))
  }]
  cat(sprintf("  %-70s %6d regions\n", s, nrow(r)))
  r
}))
if (nrow(sample_regions) == 0) stop("No ASM runs found with these settings")
# Soft score in (0,1): rises with the number of intermediate CpGs (saturating)
# and with how tightly the betas sit around 0.5 (0.5 weight at the window edge)
half_w <- max(abs(c(opt$high, opt$low) - 0.5))
sample_regions[, asm_score := pmax(0, 1 - 0.5 * mean_abs_dev_0.5 / half_w) *
                 (1 - exp(-n_intermediate / opt$minCpGs))]
sample_regions[, width := end - start + 1]
setcolorder(sample_regions, c("sample", "chr", "start", "end", "width", "n_probes",
                              "n_intermediate", "mean_beta", "mean_abs_dev_0.5", "asm_score", "probes"))
sample_regions[, c("i_start", "i_end") := NULL]

## ---- 4. Consensus across samples ------------------------------------------

sr_gr <- GRanges(sample_regions$chr, IRanges(sample_regions$start, sample_regions$end))
cons_gr <- reduce(sr_gr)
hits <- findOverlaps(cons_gr, sr_gr)
cons <- data.table(chr = as.character(seqnames(cons_gr)), start = start(cons_gr), end = end(cons_gr))
agg <- data.table(cid = queryHits(hits), sample = sample_regions$sample[subjectHits(hits)],
                  score = sample_regions$asm_score[subjectHits(hits)])
agg <- agg[, .(n_samples = uniqueN(sample), samples = paste(unique(sample), collapse = ";"),
               max_score = max(score)), by = cid]
cons[agg$cid, `:=`(n_samples = agg$n_samples, max_score = agg$max_score, samples = agg$samples)]
# Mean beta across all samples for the probes inside each consensus region
probe_gr <- GRanges(dt$chr, IRanges(dt$pos, width = 1))
ph <- findOverlaps(cons_gr, probe_gr)
bm <- as.matrix(dt[, ..sample_cols])
pm <- data.table(cid = queryHits(ph), pid = subjectHits(ph))
cons_stats <- pm[, .(n_probes = .N,
                     mean_beta_all = mean(bm[pid, ], na.rm = TRUE),
                     frac_intermediate_all = mean(bm[pid, ] >= opt$low & bm[pid, ] <= opt$high, na.rm = TRUE),
                     probes = paste(dt$probe[pid], collapse = ",")), by = cid]
cons[cons_stats$cid, `:=`(n_probes = cons_stats$n_probes, mean_beta_all = cons_stats$mean_beta_all,
                          frac_intermediate_all = cons_stats$frac_intermediate_all,
                          probes = cons_stats$probes)]
cons[, consensus := n_samples >= min(opt$minSamples, length(sample_cols))]
cons[, sex_chr := chr %in% c("chrX", "chrY")]

## ---- 5. Annotate with known imprinting / ASM ------------------------------

annotate_overlap <- function(target_gr, ref_gr, label_col) {
  h <- findOverlaps(target_gr, ref_gr)
  lab <- rep(NA_character_, length(target_gr))
  if (length(h)) {
    s <- tapply(as.character(mcols(ref_gr)[[label_col]][subjectHits(h)]), queryHits(h),
                function(v) paste(unique(v), collapse = ";"))
    lab[as.integer(names(s))] <- s
  }
  lab
}

ref_sets <- list()
if (!identical(opt$imprintome, "none") && file.exists(opt$imprintome)) {
  e <- new.env(); nm <- load(opt$imprintome, envir = e)
  icr <- get(nm[vapply(nm, function(n) is(get(n, envir = e), "GRanges"), logical(1))][1], envir = e)
  mcols(icr)$label <- paste0(icr$ID,
                             ifelse(!is.null(icr$Parental.Origin.of.Methylation) &
                                      icr$Parental.Origin.of.Methylation != "",
                                    paste0("(", icr$Parental.Origin.of.Methylation, ")"), ""),
                             ":", icr$Nearest.Transcript)
  ref_sets$ICR <- icr
}
if (!identical(opt$imprintedGenes, "none") && file.exists(opt$imprintedGenes)) {
  g <- fread(opt$imprintedGenes)
  ref_sets$imprinted_gene <- GRanges(g$Chromosome, IRanges(g$Start, g$End),
                                     label = paste0(g$`Gene name`, "(", g$`Expressed allele`, ")"))
}
if (!is.null(opt$knownRegions)) {
  k <- fread(opt$knownRegions, header = "auto")
  kn <- names(k)
  lbl <- if (ncol(k) >= 4) as.character(k[[4]]) else paste0("known_", seq_len(nrow(k)))
  ref_sets$known_ASM <- GRanges(k[[1]], IRanges(as.numeric(k[[2]]) + 1L * grepl("\\.bed$", opt$knownRegions),
                                                as.numeric(k[[3]])), label = lbl)
}
for (nm in names(ref_sets)) {
  cons[, (nm) := annotate_overlap(cons_gr, ref_sets[[nm]], "label")]
  sample_regions[, (nm) := annotate_overlap(sr_gr, ref_sets[[nm]], "label")]
}

## ---- 6. Validation ---------------------------------------------------------

val <- character()
say <- function(...) { line <- sprintf(...); val <<- c(val, line); cat(line, "\n") }
say("ASM inference: %s", in_csv)
say("Settings: beta [%.2f, %.2f], minCpGs %d, maxGap %d bp, maxMiss %d, skipNA %s, minSamples %d",
    opt$low, opt$high, opt$minCpGs, opt$maxGap, opt$maxMiss, skipNA, opt$minSamples)
say("Probes tested: %d", nrow(dt))
say("")
say("Per-sample regions:")
ps <- sample_regions[, .(regions = .N, cpgs = sum(n_intermediate)), by = sample]
for (i in seq_len(nrow(ps))) say("  %-70s %6d regions  %7d CpGs", ps$sample[i], ps$regions[i], ps$cpgs[i])
say("Consensus regions: %d total, %d supported by >= %d samples, %d in all %d samples",
    nrow(cons), sum(cons$consensus), opt$minSamples, sum(cons$n_samples == length(sample_cols)),
    length(sample_cols))

auto <- !(dt$chr %in% c("chrX", "chrY"))
in_cons <- overlapsAny(probe_gr, cons_gr[cons$consensus])
for (nm in names(ref_sets)) {
  ref <- ref_sets[[nm]]
  in_ref <- overlapsAny(probe_gr, ref)
  say("")
  say("== Validation vs %s (%d regions; autosomal probes only) ==", nm, length(ref))
  tab <- table(called = in_cons[auto], ref = in_ref[auto])
  if (all(dim(tab) == 2)) {
    ft <- fisher.test(tab)
    say("  Probes in %s: %d / %d (%.2f%%) of all tested; %d / %d (%.2f%%) of consensus-ASM probes",
        nm, sum(in_ref[auto]), sum(auto), 100 * mean(in_ref[auto]),
        sum(in_ref[auto] & in_cons[auto]), sum(in_cons[auto]), 100 * mean(in_ref[in_cons & auto]))
    say("  Enrichment odds ratio %.2f (95%% CI %.2f-%.2f), Fisher p = %.3g",
        ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value)
  }
  # Recall: reference regions with >= minCpGs probes that contain a consensus call
  rp <- countOverlaps(ref, probe_gr)
  testable <- ref[rp >= opt$minCpGs & !(as.character(seqnames(ref)) %in% c("chrX", "chrY"))]
  rec <- overlapsAny(testable, cons_gr[cons$consensus])
  say("  Recall: %d / %d %s regions with >= %d array probes overlap a consensus ASM call (%.1f%%)",
      sum(rec), length(testable), nm, opt$minCpGs, 100 * mean(rec))
  # Precision-style: consensus regions that hit the reference
  hitc <- overlapsAny(cons_gr[cons$consensus & !cons$sex_chr], ref)
  say("  Consensus autosomal ASM regions overlapping %s: %d / %d (%.1f%%)",
      nm, sum(hitc), length(hitc), 100 * mean(hitc))
}
if ("ICR" %in% names(ref_sets)) {
  say("")
  say("Canonical imprinted loci (ICRs with >= 1 probe): called in n samples")
  canon <- "(^|\\|)(H19|IGF2|KCNQ1|SNRPN|SNURF|PEG3|ZIM2|MEST|PLAGL1|GNAS|MEG3|DLK1|PEG10|NNAT|ZDBF2|GRB10|IGF2R|AIRN|L3MBTL1|NAP1L5|INPP5F|PEG13|FAM50B|ZNF331|MCTS2|DIRAS3|ZNF597|NAA60)(\\||$)"
  icr <- ref_sets$ICR
  cl <- icr[grepl(canon, icr$Nearest.Transcript) & countOverlaps(icr, probe_gr) > 0]
  h <- findOverlaps(cl, sr_gr)
  ns <- integer(length(cl))
  if (length(h)) {
    tt <- tapply(sample_regions$sample[subjectHits(h)], queryHits(h), uniqueN)
    ns[as.integer(names(tt))] <- tt
  }
  mb <- sapply(seq_along(cl), function(i) mean(bm[overlapsAny(probe_gr, cl[i]), ], na.rm = TRUE))
  np <- countOverlaps(cl, probe_gr)
  for (i in order(as.character(seqnames(cl)), start(cl)))
    say("  %-11s %-26s %-28s probes %3d  mean beta %.2f  called in %d/%d",
        cl$ID[i], paste0(seqnames(cl)[i], ":", start(cl)[i], "-", end(cl)[i]),
        substr(cl$Nearest.Transcript[i], 1, 28), np[i], mb[i], ns[i], length(sample_cols))
}

## ---- 7. Write outputs ------------------------------------------------------

setorder(cons, -n_samples, -max_score)
cons[, region_id := sprintf("ASM_%05d", seq_len(.N))]
setcolorder(cons, c("region_id", "chr", "start", "end", "n_samples", "consensus", "max_score",
                    "n_probes", "mean_beta_all", "frac_intermediate_all", "sex_chr"))
fwrite(cons, paste0(opt$outPrefix, "_consensus_regions.csv"))
fwrite(sample_regions, paste0(opt$outPrefix, "_sample_regions.csv"))
inter_probes <- unique(unlist(strsplit(sample_regions$probes, ",")))
pc <- dt[probe %in% inter_probes]
fwrite(pc, paste0(opt$outPrefix, "_probe_calls.csv"))
writeLines(val, paste0(opt$outPrefix, "_validation.txt"))
cat("\nWrote", paste0(opt$outPrefix, c("_consensus_regions.csv", "_sample_regions.csv",
                                        "_probe_calls.csv", "_validation.txt"), collapse = ", "), "\n")
