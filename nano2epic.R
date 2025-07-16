nano2epic <- function(bs, offset = 1) {
    # load packages
    library(bsseq)
    library(GenomicRanges)

    # load in epic array CpGs
    url <- "https://github.com/andymadrid/One-Off-Scripts/raw/refs/heads/main/EPIC_Array_v1_v2_CpGs_hg38.rda"
    temp <- tempfile(fileext = ".rda")
    download.file(url, destfile = temp, mode = "wb")
    load(temp)

    # create granges object from bsseq object
    nano <- GRanges(seqnames(bs), IRanges(start(bs) + offset, start(bs) + offset))

    # find overlapping CpGs between platforms
    o <- as.data.frame(findOverlaps(nano, epicArrayCpGs))

    # subset granges objects for overlapping CpGs
    bs.subset <- bs[o[,1], ]
    array.subset <- as.data.frame(epicArrayCpGs[o[,2], ])

    # get methylation levels from subset CpGs
    meth <- bsseq::getMeth(bs.subset, type = "raw")
    totCov <- getCoverage(bs.subset, type = "Cov")
    methCov <- getCoverage(bs.subset, type = "M")
    unmethCov <- as.data.frame(totCov - methCov)
    colnames(methCov) <- "Methylated_Read_Count"
    colnames(unmethCov) <- "Unmethylated_Read_Count"
    colnames(totCov) <- "Total_Read_Count"
    colnames(meth) <- "Methylation_Level"

    # bringing it all back home
    converted <- data.frame(CpG_ID = array.subset$CpG_ID,
		Methylation_Level = meth,
		Methylated_Read_Count = methCov,
		Unmethylated_Read_Count = unmethCov,
		Total_Read_Count = totCov,
		Chr = array.subset$seqnames,
		Position = array.subset$start,
		Gene = array.subset$Gene,
		Annotation = array.subset$Annotation,
		Relation_to_Island = array.subset$Relation_to_Island)

    return(converted)
}
