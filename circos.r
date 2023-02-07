### a circos plot...but sexier
library(RCircos)

# we’ll be trying out rat data…
data(UCSC.Baylor.3.4.Rat.cytoBandIdeogram)
#data(UCSC.HG19.Human.CytoBandIdeogram)
chr.exclude <- NULL
cyto.info <- UCSC.Baylor.3.4.Rat.cytoBandIdeogram
tracks.inside <- 10
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info,chr.exclude,tracks.inside,tracks.outside)
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.cyto <- RCircos.Get.Plot.Ideogram()
rcircos.position <- RCircos.Get.Plot.Positions()
RCircos.List.Plot.Parameters()
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

# load in gene data for labels
data(RCircos.Gene.Label.Data)
# this produces a data.frame that has four columns (Chromosome,chromStart,chromEnd,Gene)
# chromStart and End are just the start and end of the gene and Gene is the symbol
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data,track.num,side)
track.num <- 2
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data,name.col,track.num,side)

# plot vertical lines
lines <- DMRs.data
lines$PlotColor <- "red" # set colors to as wanted/needed
RCircos.Vertical.Line.Plot(line.data=lines,track.num=5,side,is.sorted=FALSE)
