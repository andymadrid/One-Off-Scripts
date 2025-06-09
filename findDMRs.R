# A function to find DMRs that is better than the callDMR() function from DSS
findDMRs <- function(dml, pCutoff = 1e-3, delta = 0, pt.sig = 0.5, minCG = 3, maxGap = 150) {

  # Identify which CpGs are significant
  dml <- dml %>%
    mutate(isSig = (pval < pCutoff) & (abs(diff) >= delta)) %>%
    arrange(chr, pos)

  # Filter to CpGs that pass either filter to be considered for DMRs
  candidates <- dml %>% filter(pval < pCutoff)

  if (nrow(candidates) == 0) return(data.frame())

  # Initialize
  dmrs <- list()
  current <- candidates[1, , drop=FALSE]

  for (i in 2:nrow(candidates)) {
    cpg <- candidates[i, , drop=FALSE]

    # Same chr & within maxGap then extend region
    if (cpg$chr == current$chr[nrow(current)] && cpg$pos - current$pos[nrow(current)] <= maxGap) {
      current <- bind_rows(current, cpg)
    } else {
      # Only keep region if enough CpGs and high pct.sig
      nSig <- sum(current$isSig)
      pct.sig <- nSig / nrow(current)
      if (nrow(current) >= minCG && pct.sig >= pt.sig) {
        dmrs[[length(dmrs) + 1]] <- current
      }
      current <- cpg
    }
  }

  # Check last region
  nSig <- sum(current$isSig)
  pct.sig <- nSig / nrow(current)
  if (nrow(current) >= minCG && pct.sig >= pt.sig) {
    dmrs[[length(dmrs) + 1]] <- current
  }

  # Format output
  dmrResults <- do.call(rbind, lapply(dmrs, function(region) {
    data.frame(
      chr = region$chr[1],
      start = min(region$pos),
      end = max(region$pos),
      length = max(region$pos) - min(region$pos),
      nCG = nrow(region),
      pct.sig = mean(region$isSig),
      diff.Methy = mean(region$diff),
      min.pval = min(region$pval)
    )
  }))

  return(dmrResults)
}
