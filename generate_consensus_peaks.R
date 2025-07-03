suppressMessages(suppressWarnings(library(optigrab)))

meta = opt_get(c('meta', 'm'), required = T, description = "metasheet table with headers contains peaks files and bam files: Sample header: sample ID, 
	Peaks header: files of peak bed; Bam header: bam files of corresponding sample")
cutoff = opt_get(c('cutoff', 'c'), required = F, default = 2, description = "sample depth cutoff for merging peaks, default is 2")
prefix = opt_get(c('prefix', 'p'), required = T, description = "output files prefix")
gap = opt_get(c('gap', 'g'), required = F, default = 100, description = "gap size cutoff for merging split peaks, default is 100")
workdir = opt_get(c('dir', 'd'), required = F, default = ".", description = "working directory, default is current directory")
#taxonomy = opt_get(c('tax', 't'), required = F, default = 9606, description = "Taxonomy ID, default is human 9606, ps. mouse is ")
opt_help(c('help', 'h'))



suppressMessages(suppressWarnings(library(rtracklayer)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(DiffBind)))

options(stringsAsFactors=F)

print(paste("INFO: feeding with meta sheet file:", meta))
print(paste("INFO: output files with prefix:", prefix))
print(paste("INFO: sample depth cutoff is:", cutoff))
print(paste("INFO: merging split peaks with gap size smaller than:", gap))



Handle_peaks = function(x, cutoff = 2, prefix, gap = 100) {
 	peak_granges <- lapply(x, import)
	peak_grangeslist <- GRangesList(peak_granges)
	peak_coverage <- coverage(peak_grangeslist)
	covered_ranges <- IRanges::slice(peak_coverage, lower = cutoff, rangesOnly=T)
	covered_granges <- GRanges(covered_ranges)
	export(covered_granges, paste(prefix,"consensus_nomerge.bed",sep="_"))
	covered_granges1 = reduce(covered_granges, min.gapwidth = gap)
	export(covered_granges, paste(prefix,"consensus_merged.bed",sep="_"))
 }


print("INFO: reading meta sheet file...")
info = read.table(meta, sep="\t", header=T)




listfile = as.list(info$Peaks)
print("INFO: Number of input peaks files:", length(listfile))

setwd(workdir)
print("INFO: generating consensus peaks...")
Handle_peaks(listfile, cutoff = cutoff, prefix, gap)


