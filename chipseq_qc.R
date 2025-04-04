
# ChIP-Seq quality accessment  ancillary R script for ChIP-Seq Peaks QC
# Author Zhen Y 07222023
#

options(stringsAsFactors=F)
suppressMessages(library(ChIPQC))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

args = commandArgs(T)
samplesheet = args[1] # *_prep_meta.csv
genome = args[2]  # hg38, mm10 ...
output = args[3]

samples = read.table(samplesheet, header=T, sep="\t")
chipqc = suppressMessages(ChIPQCsample(reads=samples$bamReads, peaks=samples$Peak, annotation=genome))


setwd(output)
fd = regi(chipqc); names(fd) = gsub("All", "", names(fd)); names(fd) = gsub("Long", "", names(fd));
pdf("feature_distribution.pdf", w=6, h=5)
par(mar=c(5.1, 13 ,4.1 ,2.1))
barplot(height=fd, names=names(fd), horiz=T, las=1)
dev.off()

QCs = QCmetrics(chipqc)
fileConn <- "ChIPQC_output.txt"
df = data.frame(group = c("Fragment length", "FragCC", "RelCC", "RiP", "RiBL"), value = c(QCs["FragL"],
round(FragmentLengthCrossCoverage(chipqc),3), round(RelativeCrossCoverage(chipqc), 3), round(frip(chipqc),3), round(QCs["RiBL%"],3)))
write.table(df, fileConn, sep="\t", quote=F, col.names=F, row.names=F)
saveRDS(chipqc, file = "ChIPQC_results.rds")
ChIPQCreport(chipqc, reportName="ChIP_QC_report", reportFolder=output)

