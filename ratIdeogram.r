# Ideogram for rat genome

library(karyoploteR)

setwd("~/Desktop/Postdoc\ Research/Retinal\ Methylation/")

hyper <- read.table("dmrs.ir.inr.hyper.bed",header=T)
hypo <- read.table("dmrs.ir.inr.hypo.bed",header=T)

kp <- plotKaryotype(genome="rn6",plot.type=1)

kpPlotRegions(kp, hyper, col="red")
kpPlotRegions(kp, hypo, col="blue")
