nano2epic <- function(bs, offset = 1, arrayType = NULL) {

    # load packages
    message("[nano2epic]: beginning launch sequence")
    suppressPackageStartupMessages({
        library(bsseq)
        library(GenomicRanges)
    })

    # load in epic array CpGs
    if (arrayType == "both" || is.null(arrayType)) {
        url <- "https://github.com/andymadrid/One-Off-Scripts/raw/refs/heads/main/EPIC_Array_v1_v2_CpGs_hg38.both.rda"
        message("[nano2epic]: will overlap CpGs profiled on v1 and/or v2 of the array")
    } else if (arrayType == "v1") {
        url <- "https://github.com/andymadrid/One-Off-Scripts/raw/refs/heads/main/EPIC_Array_v1_v2_CpGs_hg38.v1.rda"
        message("[nano2epic]: will overlap CpGs profiled on v1 of the array")
    } else if (arrayType == "v2") {
        url <- "https://github.com/andymadrid/One-Off-Scripts/raw/refs/heads/main/EPIC_Array_v1_v2_CpGs_hg38.v2.rda"
        message("[nano2epic]: will overlap CpGs profiled on v2 of the array")
    } else {
        stop("Option arrayType must be one of both, v1, or v2.")
    }

    temp <- tempfile(fileext = ".rda")
    download.file(url, destfile = temp, quiet = TRUE, mode = "wb")
    load(temp)
    if (arrayType == "both" || is.null(arrayType)) {
        epicArrayCpGs <- epicArrayCpGs.both
    } else if (arrayType == "v1") {
        epicArrayCpGs <- epicArrayCpGs.v1
    } else if (arrayType == "v2") {
        epicArrayCpGs <- epicArrayCpGs.v2
    }

    # create granges object from bsseq object
    message("[nano2epic]: overlapping sequence data with array coordinates")
    nano <- GRanges(seqnames(bs), IRanges(start(bs) + offset, start(bs) + offset))

    # find overlapping CpGs between platforms
    o <- as.data.frame(findOverlaps(nano, epicArrayCpGs))

    # subset granges objects for overlapping CpGs
    bs.subset <- bs[o[,1], ]
    array.subset <- as.data.frame(epicArrayCpGs[o[,2], ])

    # get methylation levels from subset CpGs
    message("[nano2epic]: generating matrices")
    meth <- bsseq::getMeth(bs.subset, type = "raw")
    totCov <- getCoverage(bs.subset, type = "Cov")
    rownames(meth) <- paste0(array.subset$CpG_ID, "_", array.subset$seqnames, "_", array.subset$start, "_", array.subset$Gene)
    rownames(totCov) <- paste0(array.subset$CpG_ID, "_", array.subset$seqnames, "_", array.subset$start, "_", array.subset$Gene)
#    methCov <- getCoverage(bs.subset, type = "M")
#    unmethCov <- as.data.frame(totCov - methCov)
#   colnames(methCov) <- "Methylated_Read_Count"
#    colnames(unmethCov) <- "Unmethylated_Read_Count"
#    colnames(totCov) <- "Total_Read_Count"
#    colnames(meth) <- "Methylation_Level"

    # bringing it all back home
    converted <- list(meth, totCov)
    names(converted) <- c("Methylation_Levels", "Total_Coverage")
    return(converted)
}
