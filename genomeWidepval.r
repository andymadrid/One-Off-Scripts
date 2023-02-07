# get the genome-wide p-value cutoff from WGBS data
gwPval <- 0.05/nrow(dmlTest)

# plot value
plotPval <- -log10(0.05/nrow(dmlTest))
