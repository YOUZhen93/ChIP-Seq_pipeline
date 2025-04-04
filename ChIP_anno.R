# Rscript to handle super enhancers from ROSE and do consensus merging and annotations


library(optparse)

Ranno_list <- list( 
    make_option(c("-p", "--peaks"), type="character", default=FALSE,
        help="ROSE output peaks *_ROSE_input_enhancer_AllStitched.table.txt, multiple file separated by comma"),
    make_option(c("-o", "--outdir"), type="character", default=FALSE, 
        help="Indicate a output directory here"),
    make_option(c("-f", "--prefix"), type="character", default=FALSE, 
        help="Tell me your output file prefix"),
    make_option(c("-t", "--thredhold"), type="character", default=FALSE, 
        help="sample depth cutoff for merging peaks, default is 1"),
    make_option(c("-g", "--gap"), type="character", default=FALSE, 
        help="gap size cutoff for merging split peaks, default is 300"),
    make_option(c("-c", "--chip1"), type="character", default="both", 
        help = "ChIPpeakAnno parameter output, default=\"both\", go read ChIPpeakAnno doc if don't know"),
    make_option(c("-d", "--chip2"), type="character", default="TSS", 
        help="ChIPpeakAnno parameter FeatureLocForDistance, default=\"TSS\", go read ChIPpeakAnno doc if don't know")
    )

opt <- parse_args(OptionParser(option_list = Ranno_list))

if (any(is.na(c(opt$peaks, opt$outdir, opt$prefix)))) {
  stop("Error: peaks input bed, output dir, and output file prefix are all required!")
}

chip1="both"
if (!is.na(opt$chip1)) {
    chip1 = opt$chip1
}

chip2="TSS"
if (!is.na(opt$chip2)) {
    chip2 = opt$chip2
}

suppressPackageStartupMessages(library(ChIPpeakAnno))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DiffBind))


peaks = opt$peaks
out_prefix = opt$prefix
outdir = opt$outdir
thredhold = opt$thredhold
gap = opt$gap



read_txt_to_GRanges = function(x) {
  x = read.table(x, header = T)
  x = makeGRangesFromDataFrame(x, keep.extra.columns=T)
  return(x)
}

# generate consensus super enhancer peaks
Handle_peaks = function(x, thredhold = 1, prefix, gap = 12500) {
 peak_granges <- lapply(x, read_txt_to_GRanges)
 peak_grangeslist <- GRangesList(peak_granges)
 peak_coverage <- coverage(peak_grangeslist)
 covered_ranges <- IRanges::slice(peak_coverage, lower = thredhold, rangesOnly=T)
 covered_granges <- GRanges(covered_ranges)
 covered_granges$id = paste("SP", seq(length(covered_granges)), sep = "_")
 covered_granges_df = as.data.frame(covered_granges)
 write.table(covered_granges_df, file = paste(prefix,"consensus_nomerge.bed",sep="_"), sep = "\t", quote = F, col.names = F, row.names = F)
 covered_granges1 = reduce(covered_granges, min.gapwidth = gap)
 covered_granges1$id = paste("SP", seq(length(covered_granges1)), sep = "_")
 covered_granges1_df = as.data.frame(covered_granges1)
 write.table(covered_granges1_df, file = paste(prefix,"consensus_merged.bed",sep="_"), sep = "\t", quote = F, col.names = F, row.names = F)
 return(covered_granges1)
}


listfile = as.list(unlist(strsplit(peaks, ",")))
print("INFO: reading meta sheet file...")

converged_peaks = Handle_peaks(listfile, thredhold = thredhold, out_prefix, gap)
#converged_peaks = Handle_peaks(listfile, thredhold = thredhold, "DM_all_StitchedPeaks", 300)


print("INFO: annotating converged peaks file ...")


ucsc.hg38.knownGene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
macs.anno <- annotatePeakInBatch(converged_peaks, 
                                 AnnotationData=ucsc.hg38.knownGene, output=chip1, FeatureLocForDistance=chip2)

macs.anno <- addGeneIDs(annotatedPeak=macs.anno, 
                        orgAnn="org.Hs.eg.db", 
                        feature_id_type="entrez_id",
                        IDs2Add="symbol")
macs = data.frame(macs.anno)
write.table(macs, paste0(outdir, "/", out_prefix, "_annotation.txt"), sep = "\t", quote=F, col.names = T, row.names = F)
macs1 = macs[,c(1:3,6,12,16)]
macs1 = macs1 %>% group_by(id) %>% mutate_at(c(5,6), paste0, collapse = ", ")
macs1 = data.frame(macs1); macs1 = macs1[!duplicated(macs1$id), ]
write.table(macs1, paste0(outdir, "/", out_prefix, "_short_annotation.txt"), sep = "\t", quote=F, col.names = T, row.names = F)

pdf(paste0(outdir, "/", out_prefix, "_Spiral_annotation.pdf"), w=8, h=6)
mypa = brewer.pal(11, "Paired")
out <- genomicElementDistribution(converged_peaks, 
                                  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                   labels = list(geneLevel = c(promoter = "Promoter", geneDownstream = "Downstream",
    geneBody = "Gene body", distalIntergenic = "Distal Intergenic"), ExonIntron = c(exon
    = "Exon", intron = "Intron", intergenic = "Intergenic"), Exons = c(utr5 = "5' UTR",
    utr3 = "3' UTR", CDS = "CDS", otherExon = "Other exon"), group = c(geneLevel =
    "Gene Level", promoterLevel = "Promoter Level", Exons = "Exon level", ExonIntron =
    "Exon/Intron/Intergenic")),
  labelColors = c(promoter = mypa[1], geneDownstream = mypa[2], geneBody =
    mypa[3], distalIntergenic = mypa[4], exon = mypa[5], intron = mypa[6],
    intergenic = mypa[7], utr5 = mypa[8], utr3 = mypa[9], CDS = mypa[10],
    otherExon = mypa[11]),
                                  promoterRegion=c(upstream=2000, downstream=500),
                                  geneDownstream=c(upstream=0, downstream=2000),
                                  promoterLevel=list(
                                    # from 5' -> 3', fixed precedence 3' -> 5'
                                    breaks = c(-2000, 0, 500),
                                    labels = c("upstream 1-2Kb", 
                                               "TSS - 500b"),
                                    colors = c("#FFE5CC",  
                                               "#FF8E32")))
dev.off()

