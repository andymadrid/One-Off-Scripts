suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer)
  library(GenomeInfoDb);  library(data.table)
})

# Feature priority, highest first. Matches ChIPseeker's ordering so that
# results stay directly comparable to an existing annotatePeak() run.
TIERS <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron",
           "Downstream", "Distal Intergenic")

# Strip the version suffix but keep the _PAR_Y tag, otherwise the X and Y
# copies of a pseudoautosomal gene collapse onto the same id.
.strip_ver <- function(x) sub("\\.[0-9]+(_PAR_Y)?$", "\\1", x)

# rtracklayer collapses GENCODE's repeated `tag` attribute to a single value, so
# the "basic" flag is lost on import. Read the attribute column directly to get
# the set of full-length, 5'-complete transcripts.
.basic_tx <- function(gtf) {
  cmd <- sprintf("gzip -dcf %s | awk -F'\\t' '$3==\"transcript\"{print $9}'",
                 shQuote(path.expand(gtf)))
  a <- data.table::fread(cmd = cmd, sep = "\n", header = FALSE, quote = "")[[1]]
  ok <- grepl('tag "basic"', a, fixed = TRUE) &
        !grepl("mRNA_start_NF", a, fixed = TRUE)
  sub('.*transcript_id "([^"]+)".*', "\\1", a[ok])
}

.cat <- function(x, sep = ";") {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x)) paste(sort(x), collapse = sep) else NA_character_
}

# ---- input -----------------------------------------------------------------
#' Coerce GRanges / BED path / data.frame to GRanges.
#' zero_based = NA auto-detects: TRUE for a file path (BED is 0-based
#' half-open by spec), FALSE for a data.frame (dmrseq/DSS/bsseq are 1-based).
as_regions <- function(x, zero_based = NA) {
  if (is(x, "GRanges")) return(x)
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop("No such file: ", x)
    if (is.na(zero_based)) zero_based <- TRUE
    x <- data.table::fread(x, header = FALSE, select = 1:3,
                           col.names = c("chr", "start", "end"))
  }
  if (is.na(zero_based)) zero_based <- FALSE
  x  <- as.data.frame(x)
  cn <- tolower(names(x))
  get <- function(o) { i <- which(cn %in% o)[1]; if (is.na(i)) NULL else x[[i]] }
  chr <- get(c("chr", "chrom", "seqnames", "chromosome"))
  st  <- get(c("start", "chromstart", "begin"))
  en  <- get(c("end", "chromend", "stop"))
  if (is.null(chr) || is.null(st) || is.null(en))
    stop("Need chr/start/end columns; found: ", paste(names(x), collapse = ", "))
  GRanges(as.character(chr),
          IRanges(as.integer(st) + if (zero_based) 1L else 0L, as.integer(en)))
}

# ---- reference -------------------------------------------------------------
#' Path to a GENCODE GTF, downloading it on first use.
#' release: "50" (human) or "M36" (mouse).
gencode_gtf <- function(release = "50", species = c("human", "mouse"),
                        dir = "ref") {
  species <- match.arg(species)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(dir, sprintf("gencode.v%s.annotation.gtf.gz", release))
  if (!file.exists(f)) {
    u <- sprintf(paste0("https://ftp.ebi.ac.uk/pub/databases/gencode/",
                        "Gencode_%s/release_%s/gencode.v%s.annotation.gtf.gz"),
                 species, release, release)
    message("Downloading ", u)
    utils::download.file(u, f, mode = "wb", quiet = TRUE)
  }
  f
}

#' Parse a GTF into the feature sets used for annotation. ~60 s and ~3.5 GB the
#' first time; cached to <gtf>.annoref.rds and instant thereafter.
build_annotation <- function(gtf = gencode_gtf(),
                             cache = sub("\\.gtf(\\.gz)?$", ".annoref.rds", gtf),
                             keep_seqlevels = paste0("chr", c(1:22, "X", "Y", "M")),
                             basic_only = TRUE,
                             force = FALSE, verbose = TRUE) {

  if (!force && !is.null(cache) && file.exists(cache)) {
    if (verbose) message("Using cached reference: ", cache)
    return(readRDS(cache))
  }
  if (verbose) message("Parsing ", basename(gtf), " (one-time) ...")

  # colnames= keeps memory sane: GENCODE carries ~20 attributes per row and
  # only these five are ever used.
  gr <- rtracklayer::import(
    gtf, format = "gtf",
    colnames = c("type", "gene_id", "gene_name", "gene_type", "transcript_id"),
    feature.type = c("gene", "transcript", "exon", "CDS",
                     "UTR", "five_prime_UTR", "three_prime_UTR"))

  # Restricting to primary chromosomes is what stops genes from becoming
  # multi-locus (and therefore droppable) in the first place.
  gr <- GenomeInfoDb::keepSeqlevels(gr, intersect(keep_seqlevels, seqlevels(gr)),
                                    pruning.mode = "coarse")
  gr$gene_id   <- .strip_ver(gr$gene_id)
  gr$gene_name <- ifelse(is.na(gr$gene_name) | !nzchar(gr$gene_name),
                         gr$gene_id, gr$gene_name)
  ty <- as.character(gr$type)

  keep3 <- function(x) {
    mcols(x) <- DataFrame(gene_id = x$gene_id, gene_name = x$gene_name,
                          gene_type = x$gene_type)
    x
  }
  genes <- keep3(gr[ty == "gene"])

  # Restrict transcript-derived features to GENCODE's `basic` set: full-length,
  # 5'-complete, representative models -- what a genome browser shows by
  # default. Without this, a retained_intron or start-not-found fragment
  # contributes a spurious TSS inside the gene body (and spurious exons over
  # real introns), so regions get called Promoter when they are intronic.
  if (basic_only) {
    keep <- gr$transcript_id %in% .basic_tx(gtf)
    gr <- gr[ty == "gene" | keep]
    ty <- as.character(gr$type)
  }

  tx    <- gr[ty == "transcript"]
  tss   <- keep3(GenomicRanges::resize(tx, 1L, fix = "start"))   # strand-aware
  exons <- keep3(gr[ty == "exon"])

  # 5'/3' UTRs: use explicit records when the GTF has them, otherwise split the
  # generic "UTR" records by position relative to that transcript's CDS.
  utr  <- gr[ty %in% c("UTR", "five_prime_UTR", "three_prime_UTR")]
  side <- c(five_prime_UTR = "utr5", three_prime_UTR = "utr3")[as.character(utr$type)]
  need <- which(is.na(side))
  if (length(need)) {
    cds <- gr[ty == "CDS"]
    cb  <- data.table(tx = cds$transcript_id, s = start(cds), e = end(cds))[
             , .(cs = min(s), ce = max(e)), by = tx]
    u <- cb[data.table(tx  = utr$transcript_id[need],
                       s   = start(utr)[need], e = end(utr)[need],
                       str = as.character(strand(utr))[need]), on = "tx"]
    side[need] <- fifelse(
      is.na(u$cs), NA_character_,
      fifelse(u$str == "-",
              fifelse(u$s > u$ce, "utr5", fifelse(u$e < u$cs, "utr3", NA_character_)),
              fifelse(u$e < u$cs, "utr5", fifelse(u$s > u$ce, "utr3", NA_character_))))
  }
  utr5 <- keep3(utr[which(side == "utr5")])
  utr3 <- keep3(utr[which(side == "utr3")])

  # Gene-level introns: gene body minus every exon of that gene, i.e. sequence
  # that is intronic in *all* isoforms.
  exl <- GenomicRanges::reduce(split(GenomicRanges::granges(exons), exons$gene_id))
  gbl <- split(GenomicRanges::granges(genes), genes$gene_id)
  ids <- intersect(names(gbl), names(exl))
  # range() is a base group generic with S4 methods and is not masked by the
  # tidyverse, so it is deliberately left unqualified.
  il  <- GenomicRanges::psetdiff(
           unlist(range(gbl[ids]), use.names = FALSE), exl[ids])
  introns <- unlist(il, use.names = FALSE)
  iid <- rep(ids, S4Vectors::elementNROWS(il))
  m   <- match(iid, genes$gene_id)
  mcols(introns) <- DataFrame(gene_id = iid, gene_name = genes$gene_name[m],
                              gene_type = genes$gene_type[m])

  out <- list(genes = genes, tss = tss, exons = exons, utr5 = utr5,
              utr3 = utr3, introns = introns,
              source = basename(gtf), basic_only = basic_only,
              built = Sys.time())
  if (verbose)
    message(sprintf("  %s genes | %s transcripts | %s exons",
                    format(length(genes), big.mark = ","),
                    format(length(tss),   big.mark = ","),
                    format(length(exons), big.mark = ",")))
  if (!is.null(cache)) {
    saveRDS(out, cache)
    if (verbose) message("Cached -> ", cache)
  }
  out
}

# ---- annotation ------------------------------------------------------------
#' Annotate regions against a reference from build_annotation().
#'
#' @param regions    GRanges, path to a .bed, or data.frame with chr/start/end
#' @param anno       list from build_annotation()
#' @param tss_region promoter window relative to the TSS, e.g. c(-5000, 0)
#' @param downstream bp past the gene end still called "Downstream"
#' @param flank_kb   width of the genes_within_<n>kb column; 0 disables it
#' @param gene_types restrict to e.g. "protein_coding"; NULL keeps all
#' @param zero_based override BED/1-based auto-detection (see as_regions)
#' @return data.frame with one row per input region, in the input order, so it
#'         is always safe to cbind() back onto the original object.
annotate_regions <- function(regions, anno,
                             tss_region = c(-5000, 0), downstream = 3000,
                             flank_kb = 100, gene_types = NULL,
                             zero_based = NA, sep = ";", verbose = TRUE) {

  gr <- as_regions(regions, zero_based)
  n  <- length(gr)
  # tolerate "1" vs "chr1" between the regions and the reference
  try(seqlevelsStyle(gr) <- seqlevelsStyle(anno$genes)[1], silent = TRUE)

  f <- anno[c("genes", "tss", "exons", "utr5", "utr3", "introns")]
  if (!is.null(gene_types)) f <- lapply(f, function(x) x[x$gene_type %in% gene_types])

  # promoter = TSS + tss_region[1] .. TSS + tss_region[2], inclusive of the TSS
  prom <- GenomicRanges::promoters(f$tss,
            upstream   = max(0L, as.integer(-tss_region[1])),
            downstream = max(0L, as.integer(tss_region[2])) + 1L)
  dstr <- GenomicRanges::flank(f$genes, width = downstream, start = FALSE)

  hit <- function(x, tier) {
    if (!length(x)) return(NULL)
    h <- findOverlaps(gr, x, ignore.strand = TRUE)
    if (!length(h)) return(NULL)
    data.table(i = queryHits(h), gene_id = x$gene_id[subjectHits(h)], tier = tier)
  }
  ht <- rbindlist(list(hit(prom, 1L), hit(f$utr5, 2L), hit(f$utr3, 3L),
                       hit(f$exons, 4L), hit(f$introns, 5L), hit(dstr, 6L)))

  gm  <- data.table(gene_id = f$genes$gene_id, gene_name = f$genes$gene_name,
                    gene_type = f$genes$gene_type)
  mid <- (start(gr) + end(gr)) / 2
  res <- data.table(idx = seq_len(n))

  feat_cols <- c("promoter_genes", "utr5_genes", "utr3_genes",
                 "exonic_genes", "intronic_genes", "downstream_genes")

  if (nrow(ht)) {
    # every (region, gene, feature) membership -> the per-feature columns, so a
    # gene that is both a promoter hit and an exon hit appears in both
    ht_all <- gm[unique(ht, by = c("i", "gene_id", "tier")), on = "gene_id"]
    # collapsed to the best tier per gene -> the single primary call
    ht <- gm[unique(ht[order(i, tier)], by = c("i", "gene_id")), on = "gene_id"]
    ht[, mid := mid[i]]

    # signed distance from the region midpoint to that gene's closest TSS;
    # negative means the region sits upstream of the TSS
    td <- data.table(gene_id = f$tss$gene_id, tpos = start(f$tss),
                     tstr = as.character(strand(f$tss)))
    j <- td[ht, on = "gene_id", allow.cartesian = TRUE]
    j[, d := fifelse(tstr == "-", tpos - mid, mid - tpos)][, ad := abs(d)]
    j <- unique(j[order(i, gene_id, ad)], by = c("i", "gene_id"))

    # primary gene: best tier, then protein_coding, then closest TSS
    j[, pc := as.integer(gene_type != "protein_coding")]
    p <- unique(j[order(i, tier, pc, ad)], by = "i")
    res[p, on = c(idx = "i"), `:=`(
      annotation      = TIERS[p$tier],   gene_symbol = p$gene_name,
      gene_id         = p$gene_id,       gene_type   = p$gene_type,
      distance_to_TSS = as.integer(round(p$d)))]

    for (t in seq_along(feat_cols)) {
      L <- ht_all[tier == t, .(v = .cat(gene_name, sep)), by = i]
      if (nrow(L)) res[L, on = c(idx = "i"), (feat_cols[t]) := L$v]
    }
  }

  # every gene whose body overlaps the region -- the part annotatePeak drops
  ho <- findOverlaps(gr, f$genes, ignore.strand = TRUE)
  if (length(ho)) {
    o <- data.table(i = queryHits(ho), gn = f$genes$gene_name[subjectHits(ho)],
                    gt = f$genes$gene_type[subjectHits(ho)])[
           , .(g = .cat(gn, sep), t = .cat(gt, sep), n = uniqueN(gn)), by = i]
    res[o, on = c(idx = "i"), `:=`(overlapping_genes = o$g,
                                   overlapping_gene_types = o$t,
                                   n_overlapping_genes = o$n)]
  }

  # nearest gene regardless of overlap -- the useful column for intergenic hits
  dn <- distanceToNearest(gr, f$genes, ignore.strand = TRUE)
  if (length(dn))
    res[queryHits(dn), `:=`(nearest_gene = f$genes$gene_name[subjectHits(dn)],
                            nearest_gene_type = f$genes$gene_type[subjectHits(dn)],
                            nearest_gene_distance = mcols(dn)$distance)]

  if (flank_kb > 0) {
    hf <- findOverlaps(suppressWarnings(GenomicRanges::trim(gr + as.integer(flank_kb * 1000))),
                       f$genes, ignore.strand = TRUE)
    if (length(hf)) {
      fl <- data.table(i = queryHits(hf),
                       gn = f$genes$gene_name[subjectHits(hf)])[
              , .(v = .cat(gn, sep)), by = i]
      res[fl, on = c(idx = "i"), (sprintf("genes_within_%gkb", flank_kb)) := fl$v]
    }
  }

  if (!"annotation" %in% names(res)) res[, annotation := NA_character_]
  res[is.na(annotation), annotation := "Distal Intergenic"]
  if (!"n_overlapping_genes" %in% names(res)) res[, n_overlapping_genes := 0L]
  res[is.na(n_overlapping_genes), n_overlapping_genes := 0L]

  res[, annotation_detail := annotation]
  res[annotation == "Promoter", annotation_detail := paste0(
        "Promoter (", fcase(abs(distance_to_TSS) <= 1000, "<=1kb",
                            abs(distance_to_TSS) <= 2000, "1-2kb",
                            abs(distance_to_TSS) <= 3000, "2-3kb",
                            default = ">3kb"), ")")]

  out <- data.frame(seqnames = as.character(seqnames(gr)), start = start(gr),
                    end = end(gr), width = width(gr),
                    as.data.frame(res[, !"idx"]),
                    stringsAsFactors = FALSE, check.names = FALSE)
  want <- c("annotation", "annotation_detail", "gene_symbol", "gene_id",
            "gene_type", "distance_to_TSS", "overlapping_genes",
            "overlapping_gene_types", "n_overlapping_genes", feat_cols,
            "nearest_gene", "nearest_gene_type", "nearest_gene_distance")
  for (w in want) if (!w %in% names(out)) out[[w]] <- NA
  head4 <- c("seqnames", "start", "end", "width")
  out <- out[, c(head4, want, setdiff(names(out), c(head4, want)))]
  rownames(out) <- NULL

  if (verbose) {
    message("Annotated ", n, " regions:")
    print(sort(table(out$annotation), decreasing = TRUE))
  }
  out
}

# ---- helpers ---------------------------------------------------------------
#' One row per region-gene pair instead of ";"-collapsed strings.
region_gene_pairs <- function(ann, column = "overlapping_genes", sep = ";") {
  dt <- as.data.table(ann)[, .rid := .I]
  s  <- strsplit(dt[[column]], sep, fixed = TRUE)
  out <- dt[rep(.rid, pmax(lengths(s), 1L))]
  out[, gene := unlist(lapply(s, function(z) if (length(z)) z else NA_character_))]
  out[!is.na(gene)][, .rid := NULL][]
}

#' Flat, de-duplicated gene vector for enrichment analysis.
region_genes <- function(ann, column = "overlapping_genes",
                         exclude_intergenic = TRUE, named_only = TRUE, sep = ";") {
  v <- ann[[column]]
  if (exclude_intergenic) v <- v[ann$annotation != "Distal Intergenic"]
  g <- unique(unlist(strsplit(stats::na.omit(v), sep, fixed = TRUE)))
  if (named_only) g <- g[!grepl("^ENS[A-Z]*G[0-9]", g)]
  sort(g)
}

# ---- command line ----------------------------------------------------------
if (sys.nframe() == 0L && !interactive()) {
  a <- commandArgs(trailingOnly = TRUE)
  if (!length(a)) {
    cat("Usage: Rscript annotate_regions.R regions.bed [out.csv]\n")
  } else {
    res <- annotate_regions(a[1], build_annotation())
    o   <- if (length(a) > 1) a[2] else sub("\\.[^.]*$", "_annotated.csv", a[1])
    write.csv(res, o, row.names = FALSE)
    message("Wrote ", o)
  }
}
