### The following assumes you have normalized, identified variable genes, scaled, and run PCA on a single-cell RNAseq sample
### The process can be run on multiple samples using lapply()

run_doubletfinder <- function(sample, preprocess = FALSE, resol = 0.1) {

  # Preprocess (normalize, find variable genes, scale) if needed
  if (preprocess == TRUE) {
  sample <- NormalizeData(seu_sample_subset)
  sample <- FindVariableFeatures(sample)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample)
  }
  
  # Find optimal number of significant PCs
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)

  # Finish pre-processing with min_pc
  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(object = sample, dims = 1:min_pc)              
  sample <- FindClusters(object = sample, resolution = resol)

  # combinations of the pN and pK
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats)
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))

  # Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)

  # Generate multiplet rate
  multiplet_rate <- NULL
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(sample@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

  sample <- doubletFinder(seu = sample, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk,
                          nExp = nExp.poi.adj)
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"

  # Return object with singlet/doublet annotations
  return(sample)

}
