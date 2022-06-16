phenoData = read.table("/workdir/aalba-trimmed/bup/targetFile.txt", header = T)
sampleNames <- as.character(phenoData$sample.name)
ID <- as.character(phenoData$sample.id)
phenoData$treatment = c(0,0,0,1,1,1)
treatment <- phenoData$treatment

file.list= as.list(sampleNames)
ID = as.list(ID)


myobj = methRead(file.list,
                 sample.id = ID,
                 pipeline = "bismarkCoverage",
                 assembly = "Aalba_v_1.1",
                 treatment = treatment, mincov = 1) # control 0
sink(file = "MethylationSummary.log")
for (i in 1:length(myobj)) {
  print(paste0(phenoData$sample.name[i], "_MethylationSummary:"))
  getMethylationStats(myobj[[i]], plot = F, both.strands = F)
}
sink()

suppressPackageStartupMessages(library("graphics"))
system("mkdir MethylationStatsPlots")
for (i in 1:length(myobj)) {
  jpeg(filename = paste0(phenoData$sample.name[i],".methylationStats.jpg"), width = 1080, height = 1080, res = 155)
  getMethylationStats(myobj[[i]], plot = T, both.strands = F)
  dev.off()
}
system("mv *.methylationStats.jpg MethylationStatsPlots")

system("mkdir CoverageStatsPlots")
for (i in 1:length(myobj)) {
  jpeg(filename = paste0(phenoData$sample.name[i],".coverageStats.jpg"), width = 1080, height = 1080, res = 155)
  getCoverageStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
  dev.off()
}
system("mv *.coverageStats.jpg CoverageStatsPlots")



# sample corelation

meth=unite(myobj, destrand=FALSE)


outputPrefix <- "methylKit"

system("mkdir DiagnosticPlots")
jpeg(filename = paste0(outputPrefix,".SampleCorelationPlot" ,".jpg"), width = 1080, height = 1080, res = 155)
getCorrelation(meth, plot = T)
dev.off()
# cluster based on methyl profile
jpeg(filename = paste0(outputPrefix,".ClusterPlot" ,".jpg"), width = 1080, height = 1080, res = 155)
clusterSamples(meth, dist = "correlation", method = "ward", plot = TRUE)
dev.off()
hc <- clusterSamples(meth, dist = "correlation", method = "ward", plot = FALSE)
hc
# PCA
jpeg(filename = paste0(outputPrefix, ".PCA.plot" ,".jpg"), width = 1080, height = 1080, res = 155)
PCASamples(meth, screeplot = TRUE)
dev.off()
jpeg(filename = paste0(outputPrefix, ".PCA2.plot" ,".jpg"), width = 1080, height = 1080, res = 155)
PCASamples(meth)
dev.off()

system("mv *.jpg DiagnosticPlots")


myDiff <- calculateDiffMeth(meth)

rawDiff <- getData(myDiff)

write.csv(rawDiff, paste0(outputPrefix, "_Final_DiffCalls.csv"), quote = F)

library(dplyr)
filtered.final.OUT <- rawDiff %>%
  filter(abs(meth.diff) > 20) %>%
  filter(pvalue <= 0.10)

write.csv(filtered.final.OUT, paste0(outputPrefix, "_Filtered_DiffCalls.csv"), quote = F)

