date='0221'
organ='Brain'
tox='human'
SRA='8AD10NAD.txt'
workplace=/home/liusai/altersplice/data/workplace/${date}${organ}
prj_path_root=/home/liusai/altersplice/data/${date}${organ}/
prj_path_starindex=/home/liusai/index/STAR/${tox}/ensem/index
prj_path_featureCounts_gtf=/home/liusai/index/reference/${tox}/ensem/Homo_sapiens.GRCh38.110.gtf
prj_path_RSEMindex=/home/liusai/index/RSEM/${tox}/ensem/RSEM_index/
prj_path_salmonindex=/home/liusai/index/Salmon/${tox}/ensem/transcripts_index
mkdir ../../${date}${organ}
mkdir ../../${date}${organ}/rawdata
mkdir ../../${date}${organ}/fastq
mkdir ../../${date}${organ}/clean
mkdir ../../${date}${organ}/aligned_bam
mkdir ../../${date}${organ}/trans_bam
mkdir ../../${date}${organ}/STARmix
mkdir ../../${date}${organ}/RSEMmix
mkdir ../../${date}${organ}/RSEM_result
mkdir $prj_path_root/salmon_result
cat ${SRA} | while read id
  do
    source /home/liusai/miniconda3/bin/activate synapse
    cd ${prj_path_root}/rawdata
    prefetch ${id}
    mkdir ${prj_path_root}/fastq/${id}
    cd ${prj_path_root}/fastq/${id}
    fasterq-dump --split-3 $prj_path_root/rawdata/${id}.sra
    mkdir ${prj_path_root}/clean/${id}
    cd ${prj_path_root}/clean/${id}
    source /home/liusai/miniconda3/bin/activate rna
    fastp -i $prj_path_root/fastq/${id}/${id}_1.fastq -o $prj_path_root/clean/${id}/${id}_1.clean.fastq -I $prj_path_root/fastq/${id}/${id}_2.fastq -O $prj_path_root/clean/${id}/${id}_2.clean.fastq
    mkdir ${prj_path_root}/STARmix/${id}
    STAR --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --runThreadN 8 --genomeDir $prj_path_starindex\
         --alignIntronMin 20 --alignIntronMax 50000 --outSAMtype BAM SortedByCoordinate --sjdbOverhang 149\
         --outSAMattrRGline ID:sample SM:sample PL:ILLUMINA --outFilterMismatchNmax 2 \
         --outSJfilterReads Unique --outSAMmultNmax 1 --outFileNamePrefix $prj_path_root/STARmix/${id}. \
         --outSAMmapqUnique 60 --readFilesIn $prj_path_root/clean/${id}/${id}_1.clean.fastq $prj_path_root/clean/${id}/${id}_2.clean.fastq
    cd $prj_path_root/STARmix
    mv ${id}.Aligned.sortedByCoord.out.bam $prj_path_root/aligned_bam
    mv ${id}.Aligned.toTranscriptome.out.bam $prj_path_root/trans_bam
done
cd $workplace
featureCounts -p -t exon -g gene_id -a $prj_path_featureCounts_gtf -o COUNTS.txt $prj_path_root/aligned_bam/*
cp -r $prj_path_RSEMindex/* ${prj_path_root}/trans_bam
cat ${SRA} | while read id
 do
  rsem-calculate-expression --paired-end -no-bam-output --alignments -p 16 $prj_path_root/trans_bam/${id}.Aligned.toTranscriptome.out.bam $prj_path_RSEMindex/RSEM_index $prj_path_root/RSEMmix/${id}
  mv $prj_path_root/RSEMmix/${id}.isoforms.results $prj_path_root/RSEM_result
done

source /home/liusai/miniconda3/bin/activate py36
touch ${date}${organ}rmats.txt
rm rmats$disease.txt
rm rmatsCON.txt
cat ${MOD} | while read id
    do
      echo -n "$prj_path_root/aligned_bam/${id}.Aligned.sortedByCoord.out.bam,">>rmats$disease.txt
done
sed -i '$ s/.$//' rmats$disease.txt
cat CON.txt | while read id
    do
      echo -n "$prj_path_root/aligned_bam/${id}.Aligned.sortedByCoord.out.bam,">>rmatsCON.txt
done
sed -i '$ s/.$//' rmatsCON.txt
source /home/liusai/miniconda3/bin/activate py36
echo rmats$disease.txt
rmats.py  \
 --b1 rmats$disease.txt --b2 rmatsCON.txt \
 --gtf $prj_path_featureCounts_gtf\
 --od  $workplace/result\
 --tmp $workplace/rmatstmp\
 -t paired \
 --readLength 149 \
 --variable-read-length \
 --cstat 0.05 \
 --novelSS \
 --libType fr-unstranded
