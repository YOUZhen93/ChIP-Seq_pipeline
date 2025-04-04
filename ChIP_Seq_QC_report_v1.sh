
#!/bin/bash

# ChIP-Seq quality accessment
# Author Zhen Y 07222023
# additional prerequisite: qualimap
module load apps/java/1.8.0_202
set -e

trap 'error_handler $LINENO' ERR
error_handler() {
    local line_num=$1
    echo "Error: pipeline failed at line $line_num"
    exit 1
}



function usage
{
    echo -e "usage: bash ChIP_Seq_QC_report_v1.sh
    --conf ${conf} --out_dir ${out_dir}
    --help <showing this message>

   ________     \\        //
          /  /   \\      //
         /  /     \\    //
        /  /       \\  //
       /  /         \\//
      /  /           ||
     /  /            ||
    /  /             ||
   /  /              ||
   ==========        ||

    "
}

while [ "$1" != "" ]; do
    case $1 in
        -c | --conf )           shift
                                conf=$1
                                ;;
        -o | --out_dir )        shift
                                out_dir=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done


if [[ ! -d ${out_dir} ]]; then
  echo ""
  echo "output folder: ${out_dir} doesn't exist"
  echo "
  Creating new output folder:
  ${out_dir}

  "
  mkdir -p ${out_dir}
fi


source ${conf}
cd ${out_dir}
out_file=${ID}.chipseq_QC.tsv
out_file_md=${ID}.chipseq_QC.markdown
dates=`date`
echo -e "# ChIP-Seq QC summary table" >> ${out_file_md}
echo -e '
<br>
<br>
<table border="1" style="border:1px solid white; border-collapse: collapse; font-size:15px; width:60%; height=40%" class="table">
  <thead>
    <tr bgcolor="#f0a1a8" style="text-align: left;">
      <th style="padding: 10px">Title</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 5px">Author</td>
      <td>Zhen Y</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Pipeline</td>
      <td>ChIP-Seq_pipeline</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Date</td>
      <td>'`date`'</td>
    </tr>
  </tbody>
</table>
<br>
<br>' >> ${out_file_md}


echo -e "ChIP-Seq QC summary table" >> ${out_file}
echo -e "Author\tZhen Y" >> ${out_file}
echo -e "Pipeline\tdefault ChIP-Seq_pipeline.sh" >> ${out_file}
echo -e "Date\t"$( echo `date` ) >> ${out_file}

echo -e "## Basic summary" >> ${out_file_md}
echo -e '
<br>
<br>
<table border="1" style="border:1px solid white; border-collapse: collapse; font-size:15px; width:60%; height=40%" class="table">
  <thead>
    <tr bgcolor="#f0a1a8" style="text-align: left;">
      <th style="padding: 10px">Items</th>
      <th>Stats</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 5px">Sample_ID</td>
      <td>'${ID}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Library_type</td>
      <td>'${Library_type}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Tissue</td>
      <td>'${Tissue}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Condition</td>
      <td>'${Condition}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Replicates</td>
      <td>'${Replicates}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">IP factor</td>
      <td>'${IP}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Platform</td>
      <td>'${Platform}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Reference_Genome</td>
      <td>'${Reference_Genome}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Alignment_software</td>
      <td>'${Alignment_software}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Peaks caller</td>
      <td>'${Peaks_caller}'</td>
    </tr>
  </tbody>
</table>
<br>
<br>' >> ${out_file_md}


# echo -e "## Basic summary" >> ${out_file_md}
# echo -e "> > **Sample_ID:**&emsp;"${ID}"<br>" >> ${out_file_md}
# echo -e "> > **Library_type:**&emsp;"${Library_type}"<br>" >> ${out_file_md}
# echo -e "> > **Tissue:**&emsp;"${Tissue}"<br>" >> ${out_file_md}
# echo -e "> > **Condition:**&emsp;"${Condition}"<br>" >> ${out_file_md}
# echo -e "> > **Replicates:**&emsp;"${Replicates}"<br>" >> ${out_file_md}
# echo -e "> > **IP factor:**&emsp;"${IP}"<br>" >> ${out_file_md}
# echo -e "> > **Platform:**&emsp;"${Platform}"<br>" >> ${out_file_md}
# echo -e "> > **Reference_Genome:**&emsp;"${Reference_Genome}"<br>" >> ${out_file_md}
# echo -e "> > **Alignment_software:**&emsp;"${Alignment_software}"<br>" >> ${out_file_md}
# echo -e "> > **Peaks_caller:**&emsp;"${Peaks_caller} >> ${out_file_md}

echo -e "Basic summary" >> ${out_file}
echo -e "Sample_ID\t"${ID} >> ${out_file}
echo -e "Library_type\t"${Library_type} >> ${out_file}
echo -e "Tissue\t"${Tissue} >> ${out_file}
echo -e "Condition\t"${Condition} >> ${out_file}
echo -e "Condition\t"${Replicates} >> ${out_file}
echo -e "IP factor\t"${IP} >> ${out_file}
echo -e "Platform\t"${Platform} >> ${out_file}
echo -e "Reference_Genome\t"${Reference_Genome} >> ${out_file}
echo -e "Alignment_software\t"${Alignment_software} >> ${out_file}
echo -e "Peaks_caller\t"${Peaks_caller} >> ${out_file}


echo -e "## Reads QC summary" >> ${out_file_md}
echo -e "Reads QC summary" >> ${out_file}

R1totalbases=$( grep "Total written (filtered):" ${trim_galore_report1} | awk '{print $4}' | sed 's/[(|)]//g' )
R2totalbases=$( grep "Total written (filtered):" ${trim_galore_report2} | awk '{print $4}' | sed 's/[(|)]//g' )
TotalbasespostQC=$( expr $(echo ${R1totalbases} | sed 's/,//g') + $(echo ${R2totalbases} | sed 's/,//g') )

echo -e '
<br>
<br>
<table border="1" style="border:1px solid white; border-collapse: collapse; font-size:15px; width:60%; height=40%" class="table">
  <thead>
    <tr bgcolor="#f0a1a8" style="text-align: left;">
      <th style="padding: 10px">Items</th>
      <th>Stats</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 5px">Total reads processed</td>
      <td>'$( grep "Total reads processed:" ${trim_galore_report1} | awk '{print $NF}' )'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Reads with adapters from R1</td>
      <td>'$( grep "Reads with adapters:" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Reads with adapters from R2</td>
      <td>'$( grep "Reads with adapters:" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g')'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Passed Reads from R1</td>
      <td>'$( grep "Reads written (passing filters):" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Passed Reads from R2</td>
      <td>'$( grep "Reads written (passing filters):" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Total bases processed from R1</td>
      <td>'$( grep "Total basepairs processed:" ${trim_galore_report1} | awk -F" " '{print $4}' )'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Total bases processed from R2</td>
      <td>'$( grep "Total basepairs processed:" ${trim_galore_report2} | awk -F" " '{print $4}' )'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Bases trimmed from R1</td>
      <td>'$( grep "Quality-trimmed:" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Bases trimmed from R2</td>
      <td>'$( grep "Quality-trimmed:" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Passed bases from R1</td>
      <td>'$( grep "Total written (filtered):" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Passed bases from R2</td>
      <td>'$( grep "Total written (filtered):" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
  </tbody>
</table>
<br>
<br>' >> ${out_file_md}

echo -e "Total reads processed\t"$( grep "Total reads processed:" ${trim_galore_report1} | awk '{print $NF}' ) >> ${out_file}
# echo -e "> > **Total reads processed:**&emsp;"$( grep "Total reads processed:" ${trim_galore_report1} | awk '{print $NF}' )"<br>" >> ${out_file_md}

echo -e "Reads with adapters from R1:\t"$( grep "Reads with adapters:" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Reads with adapters from R1:**&emsp;"$( grep "Reads with adapters:" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )"<br>" >> ${out_file_md}

echo -e "Reads with adapters from R2:\t"$( grep "Reads with adapters:" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g') >> ${out_file}
# echo -e "> > **Reads with adapters from R2:**&emsp;"$( grep "Reads with adapters:" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g')"<br>" >> ${out_file_md}

echo -e "Passed Reads from R1:\t"$( grep "Reads written (passing filters):" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Passed Reads from R1:**&emsp;"$( grep "Reads written (passing filters):" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )"<br>" >> ${out_file_md}

echo -e "Passed Reads from R2:\t"$( grep "Reads written (passing filters):" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Passed Reads from R2:**&emsp;"$( grep "Reads written (passing filters):" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' )"<br>" >> ${out_file_md}

echo -e "Total bases processed from R1:\t"$( grep "Total basepairs processed:" ${trim_galore_report1} | awk -F" " '{print $4}' ) >> ${out_file}
# echo -e "> > **Total bases processed from R1:**&emsp;"$( grep "Total basepairs processed:" ${trim_galore_report1} | awk -F" " '{print $4}' )"<br>" >> ${out_file_md}

echo -e "Total bases processed from R2:\t"$( grep "Total basepairs processed:" ${trim_galore_report2} | awk -F" " '{print $4}' ) >> ${out_file}
# echo -e "> > **Total bases processed from R2:**&emsp;"$( grep "Total basepairs processed:" ${trim_galore_report2} | awk -F" " '{print $4}' )"<br>" >> ${out_file_md}

echo -e "Bases trimmed from R1:\t"$( grep "Quality-trimmed:" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Bases trimmed from R1:**&emsp;"$( grep "Quality-trimmed:" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )"<br>" >> ${out_file_md}

echo -e "Bases trimmed from R2:\t"$( grep "Quality-trimmed:" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Bases trimmed from R2:**&emsp;"$( grep "Quality-trimmed:" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' )"<br>" >> ${out_file_md}

echo -e "Passed bases from R1:\t"$( grep "Total written (filtered):" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Passed bases from R1:**&emsp;"$( grep "Total written (filtered):" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )"<br>" >> ${out_file_md}

echo -e "Passed bases from R2:\t"$( grep "Total written (filtered):" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Passed bases from R2:**&emsp;"$( grep "Total written (filtered):" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file_md}


echo -e "## Alignment summary" >> ${out_file_md}
echo -e "Alignment summary" >> ${out_file}

echo "
Running QC summary on bam file
"
qualimap bamqc --java-mem-size=${RAM} -bam ${BAM} -c -nw 200 -hm 3 -dl 100 -gff ${GTF} -nt ${cores} -outdir ${out_dir} -outfile ${ID}_bamquali_report && echo "Bam qc is done"

reportbam=$( find -name "genome_results.txt" )
TotalReads=$( grep -oP "(?<=number of reads = ).*" ${reportbam}  )
mappedbases=$( grep -oP "(?<=number of mapped bases = ).*" ${reportbam} | tr -d " bp" )
dupreads=$( grep -oP "(?<=number of duplicated reads \(flagged\) = ).*" ${reportbam}  )
pairedreads=$( grep -oP "(?<=number of mapped paired reads \(both in pair\) = ).*" ${reportbam} | head -n1 )
singletons=$( grep -oP "(?<=number of mapped paired reads \(singletons\) = ).*" ${reportbam} | head -n1 )
insertsize=$( grep -oP "(?<=median insert size = ).*" ${reportbam} )
meanmapq=$( grep -oP "(?<=mean mapping quality = ).*" ${reportbam} )
meancov=$( grep -oP "(?<=mean coverageData = ).*" ${reportbam} )
mitoreads=$( grep chrM ${reportbam} | awk '{print $3}' )
dupratio=$(echo "$(echo ${dupreads} | tr -d -c 0-9) / $(echo ${TotalReads} | tr -d -c 0-9)" | bc -l | xargs printf "%.3f")
mitoratio=$(echo "$(echo ${mitoreads} | tr -d -c 0-9) / $(echo ${TotalReads} | tr -d -c 0-9)" | bc -l | xargs printf "%.3f")

# $(echo ${mappedbases} | tr -d -c 0-9) /
echo -e '
<br>
<br>
<table border="1" style="border:1px solid white; border-collapse: collapse; font-size:15px; width:60%; height=40%" class="table">
  <thead>
    <tr bgcolor="#f0a1a8" style="text-align: left;">
      <th style="padding: 10px">Items</th>
      <th>Stats</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 5px">Overall alignment rate</td>
      <td>'$( grep overal ${bowtie2_report} | awk '{print $1}' )'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Duplicate reads</td>
      <td>'${dupreads}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Duplicate reads ratio</td>
      <td>'${dupratio}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Paired mapped reads (MAPQ>20)</td>
      <td>'${pairedreads}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Singleton</td>
      <td>'${singletons}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Insert size</td>
      <td>'${insertsize}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Average MAPQ</td>
      <td>'${meanmapq}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Average coverage</td>
      <td>'${meancov}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Mitochondrial reads</td>
      <td>'$(echo ${mitoreads} | xargs printf "%'.f\n")'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Mitochondrial reads ratio</td>
      <td>'${mitoratio}'</td>
    </tr>
  </tbody>
</table>
<br>
<br>' >> ${out_file_md}

echo -e "Overall alignment rate\t"$( grep overal ${bowtie2_report} | awk '{print $1}' ) >> ${out_file}
# echo -e "> > **Overall alignment rate:**&emsp;"$( grep overal ${bowtie2_report} | awk '{print $1}' )"<br>" >> ${out_file_md}
echo -e "Duplicate reads\t"${dupreads} >> ${out_file}
# echo -e "> > **Duplicate reads:**&emsp;"${dupreads}"<br>" >> ${out_file_md}
echo -e "Duplicate reads ratio\t"${dupratio} >> ${out_file}
# echo -e "> > **Duplicate reads ratio:**&emsp;"${dupratio}"<br>" >> ${out_file_md}
echo -e "Paired mapped reads\t"${pairedreads} >> ${out_file}
# echo -e "> > **Paired mapped reads:**&emsp;"${pairedreads}"<br>" >> ${out_file_md}
echo -e "Singleton\t"${singletons} >> ${out_file}
# echo -e "> > **Singleton:**&emsp;"${singletons}"<br>" >> ${out_file_md}
echo -e "Insert size\t"${insertsize} >> ${out_file}
# echo -e "> > **Insert size:**&emsp;"${insertsize}"<br>" >> ${out_file_md}
echo -e "Average MAPQ\t"${meanmapq} >> ${out_file}
# echo -e "> > **Average MAPQ:**&emsp;"${meanmapq}"<br>" >> ${out_file_md}
echo -e "Average coverage\t"${meancov} >> ${out_file}
# echo -e "> > **Average coverage:**&emsp;"${meancov}"<br>" >> ${out_file_md}
echo -e "Mitochondrial reads\t"$(echo ${mitoreads} | xargs printf "%'.f\n") >> ${out_file}
# echo -e "> > **Mitochondrial reads:**&emsp;"$(echo ${mitoreads} | xargs printf "%'.f\n")"<br>" >> ${out_file_md}
echo -e "Mitochondrial reads ratio\t"${mitoratio} >> ${out_file}
# echo -e "> > **Mitochondrial reads ratio:**&emsp;"${mitoratio} >> ${out_file_md}




echo -e "## Peaks calling summary" >> ${out_file_md}
echo -e "Peaks calling summary" >> ${out_file}

echo -e "SampleID\tTissue\tFactor\tReplicate\tbamReads\tPeaks\tPeakCaller\tBlacklist\n${ID}\t${Tissue}\t${IP}\t${Replicates}\t${rawBAM}\t${Peaks}\tnarrow\t${blacklist}" > ${ID}_prep_meta.csv

if test -f ${ID}_prep_meta.csv; then
  echo "INFO: meta file exists."
else
  echo "ERROR: meta file did not generated!"
fi


samtools index ${rawBAM}
Rscript ${Chipseq_QC} ${ID}_prep_meta.csv ${Reference_Genome} ${out_dir} && echo "ChIPQC is done"

reportf=${out_dir}/ChIPQC_output.txt
FragLen=$( grep "Fragment" ${reportf} | cut -f2 )
FragCC=$( grep "FragCC" ${reportf} | cut -f2 )
RelCC=$( grep "RelCC" ${reportf} | cut -f2 )
#RiP=$( grep "RiP" ${reportf} | cut -f2 )
RiBL=$( grep "RiBL" ${reportf} | cut -f2 )
RiP1=$( samtools view -c -F 3488 -L ${Peaks} ${BAM} ) 
TiP1=$( samtools view -h -F 3488 ${BAM} | wc -l )
RiP=$( echo ${RiP1} / ${TiP1} | bc -l | xargs printf "%.3f\n" )
PeakN=$( find -name *narrowPeak | xargs wc -l | awk '{print $1}' )

echo -e '
<br>
<br>
<table border="1" style="border:1px solid white; border-collapse: collapse; font-size:15px; width:60%; height=40%" class="table">
  <thead>
    <tr bgcolor="#f0a1a8" style="text-align: left;">
      <th style="padding: 10px">Items</th>
      <th>Stats</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 5px">Fragment length</td>
      <td>'${FragLen}'</td>
    </tr>
    <tr>
      <td style="padding: 5px">Peak number</td>
      <td>'${PeakN}'</td>
    </tr>    
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">FragCC</td>
      <td>'${FragCC}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">RelCC</td>
      <td>'${RelCC}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">RiP</td>
      <td>'${RiP}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">RiBL</td>
      <td>'${RiBL}'</td>
    </tr>
  </tbody>
</table>
<br>
<br>' >> ${out_file_md}

echo -e "Fragment length\t"${FragLen} >> ${out_file}
# echo -e "> > **Fragment length:**&emsp;"${FragLen}"<br>" >> ${out_file_md}
echo -e "FragCC\t"${FragCC} >> ${out_file}
# echo -e "> > **FragCC:**&emsp;"${FragCC}"<br>" >> ${out_file_md}
echo -e "RelCC\t"${RelCC} >> ${out_file}
# echo -e "> > **RelCC:**&emsp;"${RelCC}"<br>" >> ${out_file_md}
echo -e "RiP\t"${RiP} >> ${out_file}
# echo -e "> > **RiP:**&emsp;"${RiP} >> ${out_file_md}
echo -e "RiBL\t"${RiBL} >> ${out_file}
echo -e "Peak number\t"${PeakN} >> ${out_file}
echo -e "<br>
<br>
<br>" >> ${out_file_md}

picard=`which picard.jar || true`
module unload apps/java/1.8.0_202
module load apps/java/21.0.3 
java -jar ${picard} CollectInsertSizeMetrics -I ${BAM} -O ${ID}_insertsize.txt -H insertsize.pdf

echo -e '
### Insert size histogram
<br>
<embed src="insertsize.pdf" alt="drawing" width="600" height="700"/>
<br>
<br>' >> ${out_file_md}
echo -e '
### Genomic features overlapping with peaks
<br>
<embed src="feature_distribution.pdf" alt="drawing" width="600" height="700"/>
<br>
<br>' >> ${out_file_md}
echo -e '
<br>
<img src="CCPlot.png" alt="drawing" width="600"/>
<br>' >> ${out_file_md}

python3 -m markdown ${out_file_md} > ${out_file_md}.html
