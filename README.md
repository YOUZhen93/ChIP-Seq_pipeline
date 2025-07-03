# ChIP-Seq analysis pipeline (work on linux)
## Author: Zhen Y

### main script: ChIP-Seq_pipeline.sh
example: ChIP-Seq_pipeline.sh -1 1.fastq.gz -2 2.fastq.gz -M input|peak -o outdir -a illumina -r ref_genome_bowtie_index -s sample_id -t 10 -b path_to_blacklist_bed -P primary_assembly_bed -C chrom.size \
-G gtf_annotation -B TRUE|FALSE -c 0.01 -q 20 -S hs|mm -I path_to_input_bam

Notes: -M: analysis mode: input or peak. input: inputs are input fastq and only generate input bam file and bigwig track; peak: input files are peaks fastq.
       -a: adapters has 4 options: stranded_illumina or illumina or none or BGI; users can provide their own adapter sequences, if two ends have different adapter seqs, use comma to separate them, e.g., ATGC,GTAC
       -B: TRUE for broadpeak and FALSE for narrow one.
       -P (optional): primary assembly bed file. Only keep primary assembly e.g.:

|----------|----------|----------|
| chr1| 1| 248956422|
| chr2| 1| 242193529|
| chr3| 1| 198295559|
...
### config file: ChIP_default.config
Notes used for QC script ChIP_Seq_QC_report_v1.sh, modified as needed;

### QC script: ChIP_Seq_QC_report_v1.sh
Notes: this would generated .tsv QC report along with .markdown and .html QC report files; keep generated figures with the .html file; require qualimap tool and chipseq_qc.R script;

### QC R script: chipseq_qc.R
Notes: require ChIPQC R package, currently support hg38 only (modify R script with mm10 UCSC TxDb package to adapt mouse data)


