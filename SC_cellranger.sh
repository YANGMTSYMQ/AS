a=0
info='info.txt'
while read info
    do
       ((a+=1))
       if [ $a -eq 1 ];
       then date=$info
       elif [ $a -eq 2 ];
       then organ=$info
       elif [ $a -eq 3 ];
       then tox=$info
       elif [ $a -eq 4 ];
       then disease=$info
       else length=$(($info+1))
       fi
done<$info
CON='CONsc.txt'
MOD=${disease}sc.txt
rm corresponds.txt
if [ "$tox" = "human" ] ; then
gtf=Homo_sapiens.GRCh38.110.gtf
else
gtf=Mus_musculus.GRCm39.110.gtf 
fi
prj_path_workplace=/home/liusai/altersplice/Singlecell/workplace/${date}${organ}
prj_path_root=/home/liusai/altersplice/Singlecell/data/${date}${organ}/
prj_path_starsolo=/home/liusai/index/STARsolo/${tox}
prj_path_featureCounts_gtf=/home/liusai/index/reference/${tox}/ensem/${gtf}
prj_path_RSEMindex=/home/liusai/index/RSEM/${tox}/ensem/RSEM_index/
prj_path_Cellindex=/home/liusai/index/cellranger/human/refdata-gex-GRCh38-2020-A
prj_path_whitelist=~/SC/cellranger-7.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt
i=0
j=0
k=1
l=1
while read id
  do
     ((j+=1))
     if [ ${j} -lt ${length} ];then
     echo ${id}=${disease}sc_${k}_${j}>>corresponds.txt
     else
     echo ${disease}sc_${k}>>ALL1.txt
     ((k+=1))
     j=1
     echo ${id}=${disease}sc_${k}_${j}>>corresponds.txt
     fi
done<$MOD
while read id
   do
     ((i+=1))
     if [ ${i} -lt $length ];then
     echo ${id}=CONsc_${l}_${i}>>corresponds.txt
     else
     echo CONsc_$l>>ALL1.txt
     ((l+=1)) 
     i=1
     echo ${id}=CONsc_${l}_${i}>>corresponds.txt
     fi
done<$CON
mv ALL1.txt ${k}${disease}${l}CON.txt
i=0
j=0
SRA=corresponds.txt
mkdir $prj_path_root
mkdir $prj_path_root/rawdata
mkdir $prj_path_root/fastq
mkdir $prj_path_root/clean
mkdir $prj_path_root/STARmix
mkdir $prj_path_root/Cellranger
mkdir $prj_path_root/integrate
cat ${SRA} | while read A
do
    source /home/liusai/miniconda3/bin/activate synapse
    id=${A#*=}sc
    prefetch ${A%=*}
    mv /home/liusai/ncbi/public/sra/${A%=*}.sra $prj_path_root/rawdata/${A%=*}.sra
    mkdir ${prj_path_root}/fastq/
    cd ${prj_path_root}/fastq/
    fasterq-dump --split-files $prj_path_root/rawdata/${A%=*}.sra --include-technical
#    mv /home/liusai/ncbi/public/sra/${id}.sra $prj_path_root/rawdata
#    source /home/liusai/miniconda3/bin/activate rna
#    ((i+=1))
#    echo ${i}
#    if [ $((${i}%2)) -eq 1 ] ; then 
#	    ((j+=1))
#	    mkdir $prj_path_root/integrate/${id%_*}_${j}
#      mv $prj_path_root/fastq/${id}/${id}_1.fastq $prj_path_root/integrate/${id%_*}_${j}/${id%_*}_${j}_S1_L001_R1_001.fastq.gz
#      mv $prj_path_root/fastq/${id}/${id}_2.fastq $prj_path_root/integrate/${id%_*}_${j}/${id%_*}_${j}_S1_L001_R2_001.fastq.gz
#   # elif [ $((${i}%4)) -eq 2 ] ; then
#    else
#      mv $prj_path_root/fastq/${id}/${id}_1.fastq $prj_path_root/integrate/${id%_*}_${j}/${id%_*}_${j}_S1_L002_R1_001.fastq.gz
#      mv $prj_path_root/fastq/${id}/${id}_2.fastq $prj_path_root/integrate/${id%_*}_${j}/${id%_*}_${j}_S1_L002_R2_001.fastq.gz
#   # elif [ $((${i}%4)) -eq 3 ] ; then
#    #  mv $prj_path_root/fastq/${id}/${id}_1.fastq $prj_path_root/integrate/${id%_*}_${j}/${id%_*}_${j}_S1_L003_R1_001.fastq.gz
#     # mv $prj_path_root/fastq/${id}/${id}_2.fastq $prj_path_root/integrate/${id%_*}_${j}/${id%_*}_${j}_S1_L003_R2_001.fastq.gz
#   # else 
#    #  mv $prj_path_root/fastq/${id}/${id}_1.fastq $prj_path_root/integrate/${id%_*}_${j}/${id%_*}_${j}_S1_L004_R1_001.fastq.gz
#     # mv $prj_path_root/fastq/${id}/${id}_2.fastq $prj_path_root/integrate/${id%_*}_${j}/${id%_*}_${j}_S1_L004_R2_001.fastq.gz
#      cd $prj_path_workplace
#     export PATH=/home/liusai/SC/cellranger-7.1.0:$PATH
#     cellranger count --id=${id%_*}_${j}counts \
#    	         --fastqs=$prj_path_root/integrate/${id%_*}_${j} \
#       	 --sample=${id%_*}_${j} \
#                --transcriptome=$prj_path_Cellindex
#     source /home/liusai/miniconda3/bin/activate rna   
#     cd  $prj_path_workplace/${id%_*}_${j}counts
#     STAR --runThreadN 16 \
#            --genomeDir $prj_path_starsolo \
#            --soloType CB_UMI_Simple \
#            --readFilesIn $prj_path_workplace/${id%_*}_${j}counts/outs/possorted_genome_bam.bam \
#            --readFilesCommand samtools view -F 0x100 \
#            --readFilesType SAM SE \
#            --soloInputSAMattrBarcodeSeq CR UR \
#            --soloInputSAMattrBarcodeQual CY UY \
#           --soloCBwhitelist ${prj_path_whitelist} \
#            --soloFeatures Gene SJ \
#            --soloBarcodeReadLength 28     
#     mv $prj_path_workplace/${id%_*}_${j}counts/SJ.out.tab $prj_path_workplace/${id%_*}_${j}counts/Solo.out
#     mkdir $prj_path_workplace/results
#     cd $prj_path_workplace/${id%_*}_${j}counts
#     zip -r -y ${id%_*} Solo.out 
#     mv ${id%_*}.zip $prj_path_workplace/results   
#    fi      

done     
