#!/bin/bash

#  Fri Dec 16 14:23:34 2022
#  Author: Zhen Y
#  prerequisite: bowtie2 & bowtie2 index; samtools; macs2; bedtools; picard; deeptools bamCoverage; bedGraphToBigWig
#  only support PE; MACS2 for human. and mouse
set -e

trap 'error_handler $LINENO' ERR
error_handler() {
    local line_num=$1
    echo "Error: pipeline failed at line $line_num"
    exit 1
}

helpFunc()
{
    echo "Usage: ./ChIP-Seq_pipeline.sh
    -1 <input abs path fastq1>
    -2 <input abs path fastq2>
    -M <analyze mode: input: analyze input data; peak: analyze peaks data, default: peak>
    -o <output directory>
    -a <adapter pass to trim_galore>
    -r <ref index path>
    -s <sample ID string>
    -t <thread default 4>
    -b <blacklist regions>
    -P <primary assembly bed file>
    -C <chromosome size file>
    -G <gtf or gff annotation file>
    -B <bool: broadpeak: TRUE>
    -c <macs2 qval cutoff; default=1e-2>
    -q <bowtie2 mininum MAPQ, default is 20>
  -S <species, support human and mouse, hs for human, mm for mouse, default: hs>
    -I <input/IgG bam file>
    -h <showing this message>

    ________     \\\        //
           /  /   \\\      //
          /  /     \\\    //
         /  /       \\\  //
        /  /         \\\//
       /  /           ||
      /  /            ||
     /  /             ||
    /  /              ||
    ==========        ||

    "
    echo ""
 }



mapq=20
sample="default"
output="."
cutoff=1e-2
thread=4
species="hs"
MODE="peak"
DATE=`date +"%d-%m-%y-%T" | tr ":" "-"`
while getopts "1:2:o:M:a:r:s:t:b:P:C:G:B:c:q:S:I:h" opt; do
  case $opt in
    1 ) fastq1="$OPTARG"
    ;;
    2 ) fastq2="$OPTARG"
    ;;
    o ) output="$OPTARG"
    ;;
    M ) MODE="$OPTARG"
    ;;
    a ) adapter="$OPTARG"
    ;;
    r ) refindex="$OPTARG"
    ;;
    s ) sample="$OPTARG"
    ;;
    t ) thread="$OPTARG"
    ;;
    b ) blacklist="$OPTARG"
    ;;
    P ) primary="$OPTARG"
    ;;
    C ) chromsize="$OPTARG"
    ;;
    G ) GTF="$OPTARG"
    ;;
    B ) broadpeak="$OPTARG"
    ;;
    c ) cutoff="$OPTARG"
    ;;
    q ) mapq="$OPTARG"
    ;;
    S ) species="$OPTARG"
    ;;
    I ) input="$OPTARG"
    ;;
    h ) helpFunc ; exit 0
    ;;
    \? )
     echo "Invalid Option"
     exit 1
    ;;
  esac
done


#majorgeno=/public/users/youz/database/hg38/hg38.primary_assembly.bed
BOWTIE2=`which bowtie2 || true`
SAMTOOLS=`which samtools || true`
BEDTOOLS=`which bedtools || true`
PICARD=`which picard.jar || true`
BAMCOVERAGE=`which bamCoverage || true`
MACS2=`which macs3 || true`
bgTobw=`which bedGraphToBigWig || true`
timega=`which trim_galore || true`

if [[ -z "$BOWTIE2" ]]
then
    echo -e "Error: Where is your bowtie2? \n"
    exit 1
elif [[ -z  "$SAMTOOLS" ]]
then
    echo -e "Error: Where is your samtooles? \n"
    exit 1
elif [[ -z "$BEDTOOLS" ]]
then
    echo -e "Error: Where is your bedtools? \n"
    exit 1
elif [[ -z "$PICARD" ]]
then
    echo -e "Error: Where is your picard tools? \n"
    exit 1
elif [[ -z "$BAMCOVERAGE" ]]
then
    echo -e "Error: Where is your deeptools bamCoverage? \n"
    exit 1
elif [[ -z "$MACS2" ]]
then
    echo -e "Error: Where is your MACS2? \n"
    exit 1
elif [[ -z "$bgTobw" ]]
then
    echo -e "Error: Where is your bedGraphToBigWig? \n"
    exit 1
fi


majorgeno=${primary}

if [[ "$species" == "hs" ]]; then
    gnomsize=2913022398
    macs2mode="hs"
elif [[ "$species" == "mm" ]]; then
    gnomsize=2652783500
    macs2mode="mm"
else
    echo "
    Error: invalid species! Only human (hs) and mouse (mm) are supported ...
    "
    helpFunc
    exit 2
fi


if [[ $(($# / 2)) -lt 6 ]]; then
    echo ""
    echo "only $(($# / 2)) arguments listed; Too few arguments"
    helpFunc
    exit 2
fi


if [ -d ${fastq1} ] || [ -d ${fastq2} ] || [ -d ${blacklist} ]; then
    echo ""
    echo "Please input absolute path for fastq files and blacklist file; not just folder"
    helpFunc
    exit 2
fi


if [[ ! -d ${output} ]]; then
  echo ""
  echo "output folder: ${output} doesn't exist"
  echo "
  Creating new output folder:
  ${output}

  "
  mkdir -p ${output}
fi

if test -z "$output"; then
    echo "no output dir provide, set to current directory ..."
    output="./"
fi

if test -z "$refindex"; then
    echo "no refindex specify, Plz specify a bowtie2 index prefix ..."
    exit 2
fi


if test -z "$input"; then
    echo "no input/IgG bam specify, Running without input ..."
    RUN_MACS2_1="${MACS2} callpeak -t ${output}/${sample}.sorted.mdup.filted.bam"
else
    echo "running with input ..."
    RUN_MACS2_1="${MACS2} callpeak -t ${output}/${sample}.sorted.mdup.filted.bam -c ${input}"
fi

if test -z "$broadpeak" || [[ "$broadpeak" == [Ff]* ]]; then
    echo "RUNING NARROWPEAK MODE ..."
    broadpeak="FALSE"
    RUN_MACS2="${RUN_MACS2_1} --name ${sample} -f BAMPE -g ${macs2mode} --nomodel -B -q ${cutoff} --outdir ${output}/${DATE}_${sample}"
else
    echo "RUNING BROADPEAK MODE ..."
    broadpeak="TURE"
    RUN_MACS2="${RUN_MACS2_1} --name ${sample} -f BAMPE -g ${macs2mode} --nomodel --broad -B -q ${cutoff} --broad-cutoff 1e-4 --outdir ${output}/${DATE}_${sample}"
fi


echo "

START ANALYZING CHIP-SEQ PIPELINE ....

"

echo "
input fastq files are :
${fastq1}
${fastq2}

output dir is:
${output}

reference index is:
${refindex}

sample ID is:
${sample}

blacklist is:
${blacklist}

primary assembly is:
${majorgeno}

chromosome size is:
${chromsize}

GTF/GFF annotation is:
${GTF}

MACS2 q value cutoff is:
${cutoff}

Adater is:
${adapter}
"
echo "doing trim_galore ..."
mkdir -p ${output}/clean_reads

if [[ ${thread} -gt 8 ]]; then
cores4trim=8
echo "n_thread greater than 8. using cores 8 for trim_galore"
else
cores4trim=${thread}
fi

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

ln -sf ${fastq1} ${output}/clean_reads/${sample}_R1.fastq.gz
ln -sf ${fastq2} ${output}/clean_reads/${sample}_R2.fastq.gz
export fastq1=${output}/clean_reads/${sample}_R1.fastq.gz
export fastq2=${output}/clean_reads/${sample}_R2.fastq.gz

# QC first you change the -q --stringency --length
if [[ ${adapter} != "none" && ${adapter} == *","* ]]; then
    adapter1=$( echo ${adapter} | awk -F"," '{print $1}' | tr -d " " )
    adapter2=$( echo ${adapter} | awk -F"," '{print $2}' | tr -d " " )
    echo -e "Using adapter: ${adapter1} and adapter2: ${adapter2} ..."
    RUN_TRIM="${timega} -q 20 -j ${cores4trim} --stringency 3 --gzip --paired --length 20 -a ${adapter1} -a2 ${adapter2} --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}"

elif [[ ${adapter} =~ ^[A|T|C|G]{3,} && ${adapter} != *","* ]]; then
    echo -e "Using universal adapter: ${adapter} ..."
    RUN_TRIM="${timega} -q 20 -j ${cores4trim} --stringency 3 --gzip --paired --length 20 -a ${adapter} --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}"

elif [[ ${adapter} != "none" && ${adapter} == "stranded_illumina" ]]; then
    echo -e "Using stranded_illumina default adapters ..."
    RUN_TRIM="${timega} -q 20 -j ${cores4trim} --stringency 3 --gzip --paired --length 20 --stranded_illumina --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}"

elif [[ ${adapter} == "illumina" ]]; then
    echo -e "Using illumina default adapters ..."
    RUN_TRIM="${timega} -q 20 -j ${cores4trim} --stringency 3 --gzip --paired --length 20 --illumina --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}"

elif [[ ${adapter} != "none" || ${adapter} == "BGI" ]]; then
    echo -e "Using BGI default adapters ..."
    RUN_TRIM="${timega} -q 20 -j ${cores4trim} --stringency 3 --paired -a AAGTCGGAGGCCAAGCGG -a2 AAGTCGGATCGTAGCCATG --gzip --length 20 --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}"

elif [[ ${adapter} == "none" ]]; then
    echo -e "Using trim_galore searched adapters ..."
    RUN_TRIM="${timega} -q 20 -j ${cores4trim} --stringency 3 --paired --gzip --length 20 --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}"

else
echo -e "Error: Please provide adapter sequences (separated by comma) or sequencing platform \n"
exit 1
fi

${RUN_TRIM} && echo -e "Trim_galore is done for ${sample}! \n Now changing the input fastq for alignment"

#trim_galore -q 20 -j ${thread} --illumina --stringency 3 --gzip --length 20 --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}

export fastq1=$( echo `basename ${fastq1}` | sed 's/.fq.gz/_val_1.fq.gz/g' )
export fastq2=$( echo `basename ${fastq2}` | sed 's/.fq.gz/_val_2.fq.gz/g' )
export fastq1=${output}/clean_reads/${fastq1}
export fastq2=${output}/clean_reads/${fastq2}

echo "
trim_galore is done;
newly generated clean reads are:
R1: ${fastq1}
R2: ${fastq2}
"

export BOWTIE2_INDEXES=${refindex}

echo "
###
INFO: bowtie2 aligning ...
###
"

if [[ "$MODE" == "input" ]]; then
    echo -e "\n    input mode was detected .....      \n
    Running in input mode ...........  \n"
    cd ${output}
    ${BOWTIE2} -p ${thread} -q --local -1 ${fastq1} -2 ${fastq2} -x ${refindex} -S ${output}/${sample}.sam 2>bowtie2_report.log &&
    ${SAMTOOLS} view -Sbq ${mapq} -@ ${thread} ${output}/${sample}.sam | ${SAMTOOLS} sort -@ ${thread} > ${output}/${sample}.sorted.bam &&
    java -XX:ParallelGCThreads=${thread} -jar ${PICARD} MarkDuplicates \
          I=${output}/${sample}.sorted.bam \
          O=${output}/${sample}.sorted.mdup.bam \
          TMP_DIR=${output} \
          M=${output}/${sample}_marked_dup_metrics.txt &&
    rm -rf ${output}/${sample}.sam ${output}/${sample}.sorted.bam;

    ${SAMTOOLS} view -b -L ${majorgeno} ${output}/${sample}.sorted.mdup.bam > ${output}/${sample}.sorted.mdup.maj.bam &&

    ${BEDTOOLS} intersect -abam ${output}/${sample}.sorted.mdup.maj.bam -b ${blacklist} -v > ${output}/${sample}.sorted.mdup.filted.bam &&
    ${SAMTOOLS} index ${output}/${sample}.sorted.mdup.filted.bam &&

    rm -rf ${output}/${sample}.sorted.mdup.bam;

    echo "
    ###
    INFO: now making bam to bw ...
    ###
    "

    ${BAMCOVERAGE} -b ${output}/${sample}.sorted.mdup.filted.bam -o ${sample}_10bp.CPM.bedgraph -of bedgraph --samFlagExclude 1024 --minMappingQuality 20 --normalizeUsing CPM -bs 10 -p ${thread} --effectiveGenomeSize ${gnomsize}
    awk '($1 ~ /^chr[0-9XYM]+/) {print}' ${sample}_10bp.CPM.bedgraph > ${sample}_temp; bedSort ${sample}_temp ${sample}_10bp.majchr.CPM.bedgraph;
    ${bgTobw} ${sample}_10bp.majchr.CPM.bedgraph ${chromsize} ${sample}_10bp.CPM.bw &&
    rm -rf ${sample}_temp ${sample}_10bp.CPM.bedgraph


    echo "
    ###
    INFO: bowtie2 alignment is done!
    INFO: cleaning up is done!
    INFO: input mode is done! input bam files were generated...
    ###
    "
elif [[ "$MODE" == "peak" ]]; then
    echo -e "\n    INFO: peak mode was detected .....      \n
    INFO: Running in peak mode ...........  \n"
    cd ${output}
    ${BOWTIE2} -p ${thread} -q --local -1 ${fastq1} -2 ${fastq2} -x ${refindex} --un-conc-gz ${output}/${sample}.unmapped.fq.gz -S ${output}/${sample}.sam 2>bowtie2_report.log &&
    ${SAMTOOLS} view -Sbq ${mapq} -@ ${thread} ${output}/${sample}.sam | ${SAMTOOLS} sort -@ ${thread} > ${output}/${sample}.sorted.bam &&
    java -XX:ParallelGCThreads=${thread} -jar ${PICARD} MarkDuplicates \
          I=${output}/${sample}.sorted.bam \
          O=${output}/${sample}.sorted.mdup.bam \
          TMP_DIR=${output} \
          M=${output}/${sample}_marked_dup_metrics.txt &&
    rm -rf ${output}/${sample}.sam ${output}/${sample}.sorted.bam;

    ${SAMTOOLS} view -b -L ${majorgeno} ${output}/${sample}.sorted.mdup.bam > ${output}/${sample}.sorted.mdup.maj.bam &&

    ${BEDTOOLS} intersect -abam ${output}/${sample}.sorted.mdup.maj.bam -b ${blacklist} -v > ${output}/${sample}.sorted.mdup.filted.bam &&
    ${SAMTOOLS} index ${output}/${sample}.sorted.mdup.filted.bam &&

    rm -rf ${output}/${sample}.sorted.mdup.bam;

    echo "
    ###
    INFO: bowtie2 alignment is done!
    INFO: cleaning up is done!
    ###
    "
    echo "
    ###
    INFO: Running MACS2 now ...
    ###
    "

    mkdir -p ${output}/${DATE}_${sample}
    ${RUN_MACS2} &&
    echo "
    ###
    INFO: MACS2 callpeak is done!
    ###
    "

    echo "
    ###
    INFO: now making bam to bw ...
    ###
    "

    ${BAMCOVERAGE} -b ${output}/${sample}.sorted.mdup.filted.bam -o ${sample}_10bp.CPM.bedgraph -of bedgraph --samFlagExclude 1024 --minMappingQuality 20 --normalizeUsing CPM -bs 10 -p ${thread} --effectiveGenomeSize ${gnomsize}
    awk '($1 ~ /^chr[0-9XYM]+/) {print}' ${sample}_10bp.CPM.bedgraph > ${sample}_temp; bedSort ${sample}_temp ${sample}_10bp.majchr.CPM.bedgraph;
    ${bgTobw} ${sample}_10bp.majchr.CPM.bedgraph ${chromsize} ${sample}_10bp.CPM.bw &&
    rm -rf ${sample}_temp ${sample}_10bp.CPM.bedgraph

    echo "
    ###
    INFO: Now running QC step ...
    ###
    "
    # Create config file for QC
    config_default=${SCRIPT_DIR}/ChIP_default.config
    cp ${config_default} ${output}
    config_default=${output}/$( basename ${config_default} )
    cd ${output}
    echo -e "ID="${sample} >> ${config_default}
    echo -e "IP=H3K27ac" >> ${config_default}
    echo -e "Tissue=COLO320DM" >> ${config_default}
    echo -e "Replicates=1" >> ${config_default}
    echo -e "Condition=None" >> ${config_default}
    echo -e "Peaks="$( find ~+ -name *narrowPeak ) >> ${config_default}
    echo -e "trim_galore_report1="$( find ~+ -name *1.fastq.gz_trimming_report.txt ) >> ${config_default}
    echo -e "trim_galore_report2="$( find ~+ -name *2.fastq.gz_trimming_report.txt ) >> ${config_default}
    echo -e "BAM=${output}/${sample}.sorted.mdup.filted.bam" >> ${config_default}
    echo -e "rawBAM=${output}/${sample}.sorted.mdup.maj.bam" >> ${config_default}
    echo -e "bowtie2_report=${output}/bowtie2_report.log" >> ${config_default}
    echo -e "blacklist=${blacklist}" >> ${config_default}
    echo -e "GTF=${GTF}" >> ${config_default}
    
    ${SCRIPT_DIR}/ChIP_Seq_QC_report_v1.sh --conf ${config_default} --out_dir ${output}/res_summary
    echo "
    ###
    INFO: peak mode ChIP-Seq pipeline analysis was completely finished!
    ###
    "
fi




