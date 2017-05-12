#!/usr/bin/env zsh
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
     -gtf|--annotation) 
     geneGTF=$2
     shift
     ;; 
     -thread|--ThreadN) 
     Thread=$2
     shift
     ;; 
     -index|--STAR_RSEMindex) 
     SRindex=$2
     shift
     ;; 
     -read1|--read1) 
     read1=$2
     shift
     ;; 
     -read2|--read2) 
     read2=$2
     shift
     ;;
     -intra|--cicular) 
     CIRCULAR=$2
     shift
     ;;
     -inter|--fusion)
     FUSION=$2
     shift
     ;;
     -sce|--SCE)
     SCE=$2
     shift
     ;;
     *)

esac
shift
done



if [[ -z "$geneGTF" ]]; then
   echo "Usage:"
   echo "./NCLcomparator.sh -gtf [annotation gtf file] -thread [number of thread] -index [STAR_RSEM_index folder] -read1 [read1 fastq.gz file] -read2 [read2 fastq.gz file]  -intra [circular results folder] -inter [fusion results folder] -sce [sce bed file] "
   echo ""
   echo "ERROR: no annotation gtf file !"
   echo ""
   exit
fi

if [[ -z "$SRindex" ]]; then
   echo "Usage:"
   echo "./NCLcomparator.sh -gtf [annotation gtf file] -thread [number of thread] -index [STAR_RSEM_index folder] -read1 [read1 fastq.gz file] -read2 [read2 fastq.gz file]  -intra [circular results folder] -inter [fusion results folder] -sce [sce bed file] "
   echo ""
   echo "ERROR: no STAR_RSEM_index folder !"
   echo ""
   exit
fi

if [[ -z "$read1" ]]; then
   echo "Usage:"
   echo "./NCLcomparator.sh -gtf [annotation gtf file] -thread [number of thread] -index [STAR_RSEM_index folder] -read1 [read1 fastq.gz file] -read2 [read2 fastq.gz file]  -intra [circular results folder] -inter [fusion results folder] -sce [sce bed file] "
   echo ""
   echo "ERROR: no read1 file !"
   echo ""
   exit
fi

if [[ -z "$read2" ]]; then
   echo "Usage:"
   echo "./NCLcomparator.sh -gtf [annotation gtf file] -thread [number of thread] -index [STAR_RSEM_index folder] -read1 [read1 fastq.gz file] -read2 [read2 fastq.gz file]  -intra [circular results folder] -inter [fusion results folder] -sce [sce bed file] "
   echo ""
   echo "ERROR: no read2 file !"
   echo ""
   exit
fi

if [[-z "$Thread"]]; then
   Thread=6
fi


BASEDIR=$(dirname $(readlink -f "$0"))

rm -r -f $BASEDIR\/STAR_RSEM_out
mkdir $BASEDIR\/STAR_RSEM_out
cd $BASEDIR\/STAR_RSEM_out
STAR --chimSegmentMin 10 --runThreadN $Thread --genomeDir $SRindex --readFilesIn $read1 $read2 --readFilesCommand zcat --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM
rm -r -f Aligned.out.bam
rsem-calculate-expression -p $Thread --paired-end --bam Aligned.toTranscriptome.out.bam $SRindex\/RSEM RSEMout

totalRead=$(cat $BASEDIR\/STAR_RSEM_out\/Log.final.out | grep 'Number of input reads' | awk '{print $6}')
UMRead=$(cat $BASEDIR\/STAR_RSEM_out\/Log.final.out | grep 'Uniquely mapped reads number'| awk '{print $6}')

rm -r -f $BASEDIR\/comparison
mkdir $BASEDIR\/comparison

## build exon boundary ##
echo ""
echo "Step: building exon boundary"
cd $BASEDIR\/comparison

chrCheck=$(cat $geneGTF | grep -v "#" | head -1 | awk '{print $1}' | grep 'chr')

if [[ -n "$chrCheck" ]]
then
    cat $geneGTF | awk '{for(i=1;i<=NF;i++) if($i=="exon_id"){print  $(i+1) "\t" $0} }' | awk '{for(i=1;i<=NF;i++) if($i=="gene_name"){print  $(i+1) "\t" $0} }' | awk '{for(i=1;i<=NF;i++) if($i=="gene_id"){print  $(i+1) "\t" $0} }' | awk '{print $3 "\t" $2 "\t" $1 "\t" $4 "\t" $7 "\t" $8 "\t" $10}' | sed 's/;//g' | sed 's/"//g' | awk '{print $4 "\t" $5 "\t" $6 "\t" $7 "\t" $2"_"$3}' | sort -k1,1 -k2,2n | uniq > exons.txt
    cat exons.txt | awk '{print $1 "\t" $2-6 "\t" $2+5 "\t" $1"_"$2"_"$4"_"$5}' | sort -k1,1 -k2,2n | uniq > exons_boundary.tmp1
    cat exons.txt | awk '{print $1 "\t" $3-6 "\t" $3+5 "\t" $1"_"$3"_"$4"_"$5}' | sort -k1,1 -k2,2n | uniq > exons_boundary.tmp2
    cat exons_boundary.tmp1 exons_boundary.tmp2 | sort -k1,1 -k2,2n | uniq  > exons_boundary.bed
else
    cat $geneGTF | awk '{for(i=1;i<=NF;i++) if($i=="exon_id"){print  $(i+1) "\t" $0} }' | awk '{for(i=1;i<=NF;i++) if($i=="gene_name"){print  $(i+1) "\t" $0} }' | awk '{for(i=1;i<=NF;i++) if($i=="gene_id"){print  $(i+1) "\t" $0} }' | awk '{print $3 "\t" $2 "\t" $1 "\t" "chr"$4 "\t" $7 "\t" $8 "\t" $10}' | sed 's/;//g' | sed 's/"//g' | awk '{print $4 "\t" $5 "\t" $6 "\t" $7 "\t" $2"_"$3}' | sort -k1,1 -k2,2n | uniq > exons.txt
    cat exons.txt | awk '{print $1 "\t" $2-6 "\t" $2+5 "\t" $1"_"$2"_"$4"_"$5}' | sort -k1,1 -k2,2n | uniq > exons_boundary.tmp1   
    cat exons.txt | awk '{print $1 "\t" $3-6 "\t" $3+5 "\t" $1"_"$3"_"$4"_"$5}' | sort -k1,1 -k2,2n | uniq > exons_boundary.tmp2
    cat exons_boundary.tmp1 exons_boundary.tmp2 | sort -k1,1 -k2,2n | uniq  > exons_boundary.bed
fi

cat $geneGTF | awk '{for(i=1;i<=NF;i++) if($i=="gene_name"){print  $(i+1) "\t" $0} }'  | awk '{for(i=1;i<=NF;i++) if($i=="gene_id") {print $(i+1) "\t" $0}}' | awk '{print $1 "\t" $2}' | sed 's/;//g' | sed 's/"//g' | sort  | uniq > ENSG_GeneID_convert.txt

chrCheck=$(cat $BASEDIR\/STAR_RSEM_out\/SJ.out.tab | head -10 | grep 'chr')
if [[ -n "$chrCheck" ]]
then
      cat $BASEDIR\/STAR_RSEM_out\/SJ.out.tab | awk '{print $1 "\t" $2-1 "\t" $7}' | sort -k1,1 -k2,2n | bedtools groupby -g 1,2 -c 3 -o sum > SJ_exon.pre1
      cat $BASEDIR\/STAR_RSEM_out\/SJ.out.tab | awk '{print $1 "\t" $3+1 "\t" $7}' | sort -k1,1 -k2,2n | bedtools groupby -g 1,2 -c 3 -o sum > SJ_exon.pre2
      cat SJ_exon.pre1 SJ_exon.pre2 | awk '{print $1"_"$2 "\t" $3}' | sort -k1,1 > SJ_exon.count
else
      cat $BASEDIR\/STAR_RSEM_out\/SJ.out.tab | awk '{print "chr"$1 "\t" $2-1 "\t" $7}' | sort -k1,1 -k2,2n | bedtools groupby -g 1,2 -c 3 -o sum > SJ_exon.pre1
      cat $BASEDIR\/STAR_RSEM_out\/SJ.out.tab | awk '{print "chr"$1 "\t" $3+1 "\t" $7}' | sort -k1,1 -k2,2n | bedtools groupby -g 1,2 -c 3 -o sum > SJ_exon.pre2
      cat SJ_exon.pre1 SJ_exon.pre2 | awk '{print $1"_"$2 "\t" $3}' | sort -k1,1 > SJ_exon.count      
fi 

cat $BASEDIR\/STAR_RSEM_out\/RSEMout.genes.results | awk '{print $1 "\t" $6 "\t" $7}' | sed -e '1d' > genes_RSEM.txt
join ENSG_GeneID_convert.txt genes_RSEM.txt > ENSG_GeneID_RSEM.txt  
## ENSG versus GeneID is multiple to one, choose choose large FPKM as represent of GeneID ## 
cat ENSG_GeneID_RSEM.txt | awk '{print $2 "\t" $3 "\t" $4}' |  sort -rnk2 | awk '!x[$1]++' | sort -k1,1 > ENSG_GeneID_RSEM.txt.tmp 
    

if [[ -n "$CIRCULAR" ]]; then
   
   echo ""
  
   date "+%H:%M:%S   %d/%m/%y ...... intra-Comparation analysis"

   rm -r -f  $BASEDIR\/comparison\/intra
   mkdir $BASEDIR\/comparison\/intra

   cd $BASEDIR\/comparison
   ## intra results ##
   echo ""
   echo "Step: adjusting intra results to exon boundary"
   ls $CIRCULAR > sampleintra.tmp 
   sampleintra_num=$(cat sampleintra.tmp | wc -l)
   echo "intra sample number: "$sampleintra_num

   for i in $(seq 1 "$sampleintra_num")
   do 
     getOne=$(cat sampleintra.tmp | awk 'NR==k {print $1}' k=$i)
     echo $getOne
     getOneName=$(echo $getOne | sed 's/[.]/\t/g' | awk '{print $1}') 
     cat $CIRCULAR/$getOne | sort -k1,1 -k2,2n | awk '$1~/^chr[0-9XY]*$/ {print $0}' | awk '$3~/^chr[0-9XY]*$/ {print $0}' | awk '$2~/^[0-9]*$/ {print $0}'| awk '$4~/^[0-9]*$/ {print $0}' | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $1"_"$2"_"$3"_"$4"_"$5}' | sort -k1,1 -k2,2n | uniq  > intra_tmp1.bed
     bedtools intersect -a intra_tmp1.bed -b exons_boundary.bed -wa -wb | awk '{print $4 "\t" $8}' | sed 's/_/\t/g' | awk '{print $3 "\t" $4-1 "\t" $4 "\t" $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_"$9}' | sort -k1,1 -k2,2n | uniq > intra_tmp2.bed
     bedtools intersect -a intra_tmp2.bed -b exons_boundary.bed  -wa -wb | awk '{print $4 "\t" $8}' | sed 's/_/\t/g' | awk '{print $6 "\t" $7 "\t" $8 "\t" $10 "\t" $11 "\t" $12 "\t" $5 "\t" $9 "\t" $13}'| sort -k1,1 -k2,2n | uniq | awk '$3==$6 && $8==$9 {print $0}'  > intra_tmp.result
     cat intra_tmp.result | sed 's/+/sense/g' | sed 's/-/anti/g' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8}' | sort -k1,1 -k2,2n | uniq | bedtools groupby -g 1,2,3,4,5,6,7 -c 8 -o collapse | sed 's/sense/+/g' | sed 's/anti/-/g' | awk '{print $0 "\t" $8}'| sort -k1,1 -k2,2n | uniq > intra_tmp1.result
     cat intra_tmp1.result | awk '{if($3=="+" && $2 > $5) print $0; else if($3=="+" && $2 < $5) print $1 "\t" $5 "\t" $3 "\t" $4 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9; else if($3=="-" && $2 < $5) print $0; else if($3=="-"&& $2> $5) print $1 "\t" $5 "\t" $3 "\t" $4 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9}'| sort -k1,1 -k2,2n | uniq > intra_tmp2.result
     cat intra_tmp2.result | sed 's/+/sense/g' | sed 's/-/anti/g' | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$8 "\t" $7}' | sort | bedtools groupby -g 1 -c 2 -o sum > intra_tmp3.result  
     cat intra_tmp3.result | sed 's/_/\t/g'| sed 's/sense/+/g' | sed 's/anti/-/g' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $8 "\t" $7 "\t" $7 }'> intra\/$getOneName.circEB 
 
   done

   ## merge intra ##
    echo ""
    echo "Step: merging intra results"
    cd  $BASEDIR
    cat comparison\/intra\/*.circEB | sort -k1,1 -k2,2n | uniq  > all.circEB
    cat all.circEB | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" $8}' | sort -k1 | uniq > all.circEB.tmp 

    ls comparison\/intra > sampleintra 
    sampleintra_num=$(cat sampleintra | wc -l)
    echo "intra sample number: "$sampleintra_num

    header=$(echo "Chr" "\t" "Pos" "\t" "Strand" "\t" "Chr" "\t" "Pos" "\t" "Strand" "\t" "Genename")  
    sampleintra_num=$(echo $sampleintra_num)
    for i in $(seq 1 $sampleintra_num)
    do 
       getOne=$(cat sampleintra | awk 'NR==k {print $1}' k=$i)
       echo $getOne
       getOneName=$(echo $getOne | sed 's/[.]/\t/g' | awk '{print $1}')
       header=$(echo $header "\t" $getOneName)
       cat comparison\/intra\/$getOne | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" 1}' | sort -k1 > $getOneName.circEB.tmp
       join -o 1.1 1.2 2.2 all.circEB.tmp  $getOneName.circEB.tmp -a1 -e0 | awk '{print $1 "\t" $2":"$3}' | sort -k1 | uniq > all.circEB_2.tmp
       cat all.circEB_2.tmp > all.circEB.tmp
    done

    cat all.circEB.tmp  | sed 's/:/\t/g' | tr ' ' \\t > intraMerged.tmp
    cat <(echo $header) intraMerged.tmp > intraMerged.result
    rm -r -f *.tmp
    rm -r -f all.circEB
    rm -r -f sampleintra
    $BASEDIR\/bin\/graphic.R intraMerged.result intra

   echo "Step - Characteristics of intra-NCL events"
   rm -r -f $BASEDIR\/comparison\/intra_Characteristics
   mkdir $BASEDIR\/comparison\/intra_Characteristics
   cd $BASEDIR\/comparison
   ls intra > sampleCircEB.tmp 
   sampleCircEB_num=$(cat sampleCircEB.tmp | wc -l)
   echo "intra sample number: "$sampleCircEB_num

   for i in $(seq 1 "$sampleCircEB_num")
   do 
       getOne=$(cat sampleCircEB.tmp | awk 'NR==k {print $1}' k=$i)
       echo $getOne
       getOneName=$(echo $getOne | sed 's/[.]/\t/g' | awk '{print $1}') 
       #echo "Step - to calcuate LR, RPM_MappedReads, RPM_totalReads, CF, and NCLratio"
       cat intra/$getOne > NCLsample.intra
       cat NCLsample.intra | awk '{print $1"_"$2 "\t" $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_"$9}' | sort -k1,1 > NCLsample.intra.pre1
       join -o 1.1 1.2 2.2 NCLsample.intra.pre1 SJ_exon.count -a1 -e0 | awk '{print $2 "\t" $3}' | sed 's/_/\t/g' | awk '{print $4"_"$5 "\t" $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10}' | sort -k1,1 > NCLsample.intra.pre2
       join -o 1.1 1.2 2.2 NCLsample.intra.pre2 SJ_exon.count -a1 -e0 | awk '{print $2 "\t" $3}' | sed 's/_/\t/g' > NCLsample.count
       cat NCLsample.count | awk '{print $0 "\t" ($7/M)*1000000 "\t" ($7/T)*1000000 "\t" $7/($10+$11+1) "\t" (2*$7)/(2*$7+$10+$11)}' T=$totalRead  M=$UMRead > NCLsample.events

       #echo "Step - to merge FPKM and TPM"
       ## remove junction side overlap in more than one geneID ##
       cat NCLsample.events | awk '$8!~/[,]/ {print $0}' | awk '{print $8 "\t" $0}' | sort  -k1,1 > NCLsample.events.tmp 
       join NCLsample.events.tmp ENSG_GeneID_RSEM.txt.tmp  | cut -f 2- -d' ' > NCLsample_RSEM.txt

       #echo "Step - to calculate SJpj using CalSJpj.sh"
       ## generate SJ_pj.txt and gene_pmed.txt ##
       $BASEDIR/bin/CalSJpj.sh intra/$getOne $geneGTF $BASEDIR\/STAR_RSEM_out\/SJ.out.tab

       #echo "Step - to calculate pj and pmed"
       ### merge SJ_pj ###
       cat SJ_pj.txt | awk '{print $1"_"$2 "\t" $4}' | sort -k1,1 | uniq > SJ_pj.tmp 
       cat gene_pmed.txt | sort -k1,1 | uniq > gene_pmed.tmp
       ## one geneID ##
       cat NCLsample_RSEM.txt | awk '{print $1"_"$2"_"$8 "\t" $0}' | sort -k1,1 > NCLsample_RSEM.txt.tmp1
       join NCLsample_RSEM.txt.tmp1 SJ_pj.tmp | cut -f 2- -d ' ' | sort -k1,1 -k2,2n | uniq  > NCLsample_RSEM_pj.tmp1
       join NCLsample_RSEM.txt.tmp1 SJ_pj.tmp -v1 | cut -f 2- -d ' ' | sort -k1,1 -k2,2n | uniq | awk '{print $0 " " "0"}' > NCLsample_RSEM_pj.tmp2
       cat NCLsample_RSEM_pj.tmp1 NCLsample_RSEM_pj.tmp2 | awk '{print $1"_"$5"_"$8 "\t" $0}' | sort -k1,1 > NCLsample_RSEM_pj.tmp3
       join NCLsample_RSEM_pj.tmp3 SJ_pj.tmp | cut -f 2- -d ' ' | sort -k1,1 -k2,2n | uniq > NCLsample_RSEM_pj.tmp4
       join NCLsample_RSEM_pj.tmp3 SJ_pj.tmp -v1 | cut -f 2- -d ' ' | sort -k1,1 -k2,2n | uniq | awk '{print $0 " " "0"}' > NCLsample_RSEM_pj.tmp5
       cat NCLsample_RSEM_pj.tmp4 NCLsample_RSEM_pj.tmp5 | sort -k1,1 -k2,2n | uniq > NCLsample_RSEM_pj.txt        

       cat NCLsample_RSEM_pj.txt | awk '{print $8 "\t" $0}' | sort -k1,1 | uniq > NCLsample_RSEM_pj.txt.tmp1
       join NCLsample_RSEM_pj.txt.tmp1 gene_pmed.tmp | cut -f 2- -d ' '| sort -k1,1 -k2,2n | uniq > NCLsample_RSEM_pj_pmed.tmp1
       join NCLsample_RSEM_pj.txt.tmp1 gene_pmed.tmp -v1 | cut -f 2- -d ' '| sort -k1,1 -k2,2n | uniq | awk '{print $0 " " "0"}' > NCLsample_RSEM_pj_pmed.tmp2
       cat NCLsample_RSEM_pj_pmed.tmp1 NCLsample_RSEM_pj_pmed.tmp2 | sort -k1,1 -k2,2n | uniq > NCLsample_RSEM_pj_pmed.txt

       #echo "Step - to examine out of circular using OutOfCircular.sh"
       ## generate OC.final ##
       $BASEDIR/bin/OutOfCircular.sh intra/$getOne $BASEDIR\/STAR_RSEM_out\/Chimeric.out.junction $BASEDIR\/STAR_RSEM_out\/Chimeric.out.sam $BASEDIR\/STAR_RSEM_out\/Log.final.out
       cat OC.final | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6 "\t" $8}' | sort -k1 > OC.final.tmp
       cat NCLsample_RSEM_pj_pmed.txt | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6 "\t" $0}' | sort -k1 > NCLsample_RSEM_pj_pmed_OC.txt.tmp1
       join NCLsample_RSEM_pj_pmed_OC.txt.tmp1 OC.final.tmp | cut -f 2- -d ' ' | sort -k1,1 -k2,2n | uniq > NCLsample_RSEM_pj_pmed_OC.txt.tmp2 
       join NCLsample_RSEM_pj_pmed_OC.txt.tmp1 OC.final.tmp -v1 | cut -f 2- -d ' ' | sort -k1,1 -k2,2n | uniq | awk '{print $0 " " "0"}' > NCLsample_RSEM_pj_pmed_OC.txt.tmp3
       cat NCLsample_RSEM_pj_pmed_OC.txt.tmp2 NCLsample_RSEM_pj_pmed_OC.txt.tmp3 | sort -k1,1 -k2,2n | uniq > NCLsample_RSEM_pj_pmed_OC.txt

       if [[ -n "$SCE" ]]; then 
          #echo "Step - to determine donor/acceptor side within SCE"
          paste <(cat NCLsample_RSEM_pj_pmed_OC.txt | awk '{print $1 "\t" $2-1 "\t" $2}')  <(cat NCLsample_RSEM_pj_pmed_OC.txt | tr ' ' '_')  > NCLsample_RSEM_pj_pmed_OC_SCE.tmp1
          bedtools intersect -a NCLsample_RSEM_pj_pmed_OC_SCE.tmp1 -b ../example/SCE/SCE_hg38.bed -wa | sort | uniq | awk '{print $4"_1"}' > NCLsample_RSEM_pj_pmed_OC_SCE.tmp2
          bedtools intersect -a NCLsample_RSEM_pj_pmed_OC_SCE.tmp1 -b ../example/SCE/SCE_hg38.bed -v | sort | uniq | awk '{print $4"_0"}' > NCLsample_RSEM_pj_pmed_OC_SCE.tmp3
          cat NCLsample_RSEM_pj_pmed_OC_SCE.tmp2 NCLsample_RSEM_pj_pmed_OC_SCE.tmp3 | tr '_' ' ' > NCLsample_RSEM_pj_pmed_OC_SCE.tmp4
          paste <(cat NCLsample_RSEM_pj_pmed_OC_SCE.tmp4 | awk '{print $4 "\t" $5-1 "\t" $5}')  <(cat NCLsample_RSEM_pj_pmed_OC_SCE.tmp4 | tr ' ' '_')  > NCLsample_RSEM_pj_pmed_OC_SCE.tmp5
          bedtools intersect -a NCLsample_RSEM_pj_pmed_OC_SCE.tmp5 -b ../example/SCE/SCE_hg38.bed -wa  | sort | uniq | awk '{print $4"_1"}' > NCLsample_RSEM_pj_pmed_OC_SCE.tmp6
          bedtools intersect -a NCLsample_RSEM_pj_pmed_OC_SCE.tmp5 -b ../example/SCE/SCE_hg38.bed -v  | sort | uniq| awk '{print $4"_0"}' > NCLsample_RSEM_pj_pmed_OC_SCE.tmp7
          cat NCLsample_RSEM_pj_pmed_OC_SCE.tmp6 NCLsample_RSEM_pj_pmed_OC_SCE.tmp7 | tr '_' ' ' > NCLsample_RSEM_pj_pmed_OC_SCE.txt
          echo "Chr_donor Pos_donor Strand_donor Chr_acceptot Pos_acceptor Strand_acceptor JunctionRead Gene_donor Gene_acceptor LinearRead_donor LinearRead_acceptor RPM_totalRead RPM_mappedRead CF NCLratio FPKM TPM pj_dornor pj_acceptor pmedian OC SCE_donor SCE_acceptor" > HeaderPost.txt
          cat HeaderPost.txt NCLsample_RSEM_pj_pmed_OC_SCE.txt | tr ' ' \\t > intra_Characteristics/$getOneName.final
       fi

       if [[ -z "$SCE" ]]; then 
          echo "Chr_donor Pos_donor Strand_donor Chr_acceptot Pos_acceptor Strand_acceptor JunctionRead Gene_donor Gene_acceptor LinearRead_donor LinearRead_acceptor RPM_totalRead RPM_mappedRead CF NCLratio FPKM TPM pj_dornor pj_acceptor pmedian OC" > HeaderPost.txt
          cat HeaderPost.txt NCLsample_RSEM_pj_pmed_OC.txt | tr ' ' \\t > intra_Characteristics/$getOneName.final
       fi
   done

fi

if [[ -n "$FUSION" ]]; then

    echo ""
    date "+%H:%M:%S   %d/%m/%y ...... intra-Comparation analysis"

    rm -r -f  $BASEDIR\/comparison\/inter
    mkdir $BASEDIR\/comparison\/inter
    cd $BASEDIR\/comparison
    ## inter results ##
    echo ""
    echo "Step: adjusting inter results to exon boundary"
    ls $FUSION > sampleinter.tmp 
    sampleinter_num=$(cat sampleinter.tmp | wc -l)
    echo "inter sample number: "$sampleinter_num
    for i in $(seq 1 "$sampleinter_num")
    do 
       getOne=$(cat sampleinter.tmp | awk 'NR==k {print $1}' k=$i)
        echo $getOne
        getOneName=$(echo $getOne | sed 's/[.]/\t/g' | awk '{print $1}') 
        cat $FUSION/$getOne | awk '$1~/^chr[0-9XY]*$/ {print $0}' | awk '$3~/^chr[0-9XY]*$/ {print $0}' | awk '$2~/^[0-9]*$/ {print $0}'| awk '$4~/^[0-9]*$/ {print $0}'| awk '{print $1 "\t" $2-1 "\t" $2 "\t" $1"_"$2"_"$3"_"$4"_"$5}' | sort -k1,1 -k2,2n | uniq  > inter_tmp1.bed
        bedtools intersect -a inter_tmp1.bed -b exons_boundary.bed -wa -wb | awk '{print $4 "\t" $8}' | sed 's/_/\t/g' | awk '{print $3 "\t" $4-1 "\t" $4 "\t" $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_"$9}' | sort -k1,1 -k2,2n | uniq > inter_tmp2.bed
        bedtools intersect -a inter_tmp2.bed -b exons_boundary.bed  -wa -wb | awk '{print $4 "\t" $8}' | sed 's/_/\t/g' | awk '{print $6 "\t" $7 "\t" $8 "\t" $10 "\t" $11 "\t" $12 "\t" $5 "\t" $9 "\t" $13}' | sort -k1,1 -k2,2n | uniq | awk '$8!=$9 {print $0}' > inter_tmp1.result
        cat inter_tmp1.result | awk '{if($2 >= $5) print $0; else if($2 < $5) print $4 "\t" $5 "\t" $6 "\t" $1 "\t" $2 "\t" $3 "\t" $7 "\t" $9 "\t" $8}' | sed 's/+/sense/g' | sed 's/-/anti/g'| sort -k1,1 -k2,2n | uniq > inter_tmp2.result
        cat inter_tmp2.result | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6 "\t" $7}' | sort | bedtools groupby -g 1 -c 2 -o sum > inter_tmp3.result
        cat inter_tmp2.result | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6 "\t" $8}' | sort | uniq | bedtools groupby -g 1 -c 2 -o collapse > inter_tmp4.result
        cat inter_tmp2.result | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6 "\t" $9}' | sort | uniq | bedtools groupby -g 1 -c 2 -o collapse > inter_tmp5.result
        join inter_tmp3.result inter_tmp4.result | sort > inter_tmp6.result
        join inter_tmp6.result inter_tmp5.result| sed 's/_/\t/g' | sed 's/sense/+/g' | sed 's/anti/-/g' | tr ' ' \\t | sort -k1,1 -k2,2n | uniq > inter\/$getOneName.fusionEB
     done

    ## merge inter ##
    echo ""
    echo "Step: merging inter results"
    cd  $BASEDIR
    cat comparison/inter/*.fusionEB | sort -k1,1 -k2,2n | uniq  > all.fusionEB
    cat all.fusionEB | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" $8":"$9}' | sort -k1 | uniq > all.fusionEB.tmp 

    ls comparison/inter > sampleinter 
    sampleinter_num=$(cat sampleinter | wc -l)
    echo "inter sample number: "$sampleinter_num
    header=$(echo "Chr" "\t" "Pos" "\t" "Strand" "\t" "Chr" "\t" "Pos" "\t" "Strand" "\t" "Genename" "\t" "Genename")
    sampleinter_num=$(echo $sampleinter_num)
    for i in $(seq 1 $sampleinter_num)
    do 
        getOne=$(cat sampleinter | awk 'NR==k {print $1}' k=$i)
        echo $getOne
        getOneName=$(echo $getOne | sed 's/[.]/\t/g' | awk '{print $1}')
        header=$(echo $header "\t" $getOneName)
        cat comparison/inter/$getOne | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" 1}' | sort -k1 > $getOneName.fusionEB.tmp
        join -o 1.1 1.2 2.2 all.fusionEB.tmp  $getOneName.fusionEB.tmp -a1 -e0 | awk '{print $1 "\t" $2":"$3}' | sort -k1 | uniq > all.fusionEB_2.tmp
        cat all.fusionEB_2.tmp > all.fusionEB.tmp
    done

    cat all.fusionEB.tmp  | sed 's/:/\t/g' | tr ' ' \\t > interMerged.tmp
    cat <(echo $header) interMerged.tmp > interMerged.result
    cut  -d$'\t' -f 7 --complement interMerged.result > interMerged.result.tmp
    $BASEDIR/bin/graphic.R interMerged.result.tmp inter
    rm -r -f *.tmp
    rm -r -f all.fusionEB
    rm -r -f sampleinter
    
   
    echo "Step - Characteristics of inter-NCL events"
    cd $BASEDIR\/comparison
    rm -r -f $BASEDIR\/comparison\/inter_Characteristics
    mkdir $BASEDIR\/comparison\/inter_Characteristics

    ls inter > sampleFusionEB.tmp 
    sampleFusionEB_num=$(cat sampleFusionEB.tmp | wc -l)
    echo "inter sample number: "$sampleFusionEB_num

    for i in $(seq 1 "$sampleFusionEB_num")
    do 
        getOne=$(cat sampleFusionEB.tmp | awk 'NR==k {print $1}' k=$i)
        echo $getOne
        getOneName=$(echo $getOne | sed 's/[.]/\t/g' | awk '{print $1}') 
        #echo "Step - to RPM_MappedReads and RPM_totalReads"
        cat inter/$getOne > NCLsample.inter
        cat NCLsample.inter | awk '{print $1"_"$2 "\t" $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_"$9}' | sort -k1,1 > NCLsample.inter.pre1
        join -o 1.1 1.2 2.2 NCLsample.inter.pre1 SJ_exon.count -a1 -e0 | awk '{print $2 "\t" $3}' | sed 's/_/\t/g' | awk '{print $4"_"$5 "\t" $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10}' | sort -k1,1 > NCLsample.inter.pre2
        join -o 1.1 1.2 2.2 NCLsample.inter.pre2 SJ_exon.count -a1 -e0 | awk '{print $2 "\t" $3}' | sed 's/_/\t/g' > NCLsampleinter.count
        cat NCLsampleinter.count | awk '{print $0 "\t" ($7/M)*1000000 "\t" ($7/T)*1000000 }' T=$totalRead  M=$UMRead > NCLsampleinter.events

        #echo "Step - to merge FPKM and TPM"
        ## one geneID ##
        cat NCLsampleinter.events | awk '$8!~/[,]/ {print $0}' | awk '{print $8 "\t" $0}' | sort  -k1,1 > NCLsampleinter.events.tmp1 
        join NCLsampleinter.events.tmp1 ENSG_GeneID_RSEM.txt.tmp  | cut -f 2- -d' ' > NCLsampleinter_RSEM.tmp1
        cat NCLsampleinter_RSEM.tmp1 | awk '$9!~/[,]/ {print $0}' | awk '{print $9 "\t" $0}' | sort  -k1,1 > NCLsampleinter.events.tmp2 
        join NCLsampleinter.events.tmp2 ENSG_GeneID_RSEM.txt.tmp  | cut -f 2- -d' ' > NCLsampleinter_RSEM.txt

        #echo "Step - to calculate SJpj using CalSJpj.sh"
        ## generate SJ_pj.txt and gene_pmed.txt ##
        $BASEDIR/bin/CalSJpj.sh inter/$getOne $geneGTF $BASEDIR\/STAR_RSEM_out\/SJ.out.tab


        #echo "Step - to calculate pj and pmed"
        ## merge SJ_pj ##
        cat SJ_pj.txt | awk '{print $1"_"$2 "\t" $4}' | sort -k1,1 | uniq > SJ_pj.tmp 
        cat gene_pmed.txt | sort -k1,1 | uniq > gene_pmed.tmp

        cat NCLsampleinter_RSEM.txt | awk '{print $1"_"$2"_"$8 "\t" $0}' | sort -k1,1 > NCLsampleinter_RSEM.txt.tmp1
        join NCLsampleinter_RSEM.txt.tmp1 SJ_pj.tmp | cut -f 2- -d ' ' | sort -k1,1 -k2,2n | uniq  > NCLsampleinter_RSEM_pj.tmp1
        join NCLsampleinter_RSEM.txt.tmp1 SJ_pj.tmp -v1 | cut -f 2- -d ' ' | sort -k1,1 -k2,2n | uniq | awk '{print $0 " " "0"}' > NCLsampleinter_RSEM_pj.tmp2
        cat NCLsampleinter_RSEM_pj.tmp1 NCLsampleinter_RSEM_pj.tmp2 | awk '{print $4"_"$5"_"$9 "\t" $0}' | sort -k1,1 > NCLsampleinter_RSEM_pj.tmp3
        join NCLsampleinter_RSEM_pj.tmp3 SJ_pj.tmp | cut -f 2- -d ' ' | sort -k1,1 -k2,2n | uniq > NCLsampleinter_RSEM_pj.tmp4
        join NCLsampleinter_RSEM_pj.tmp3 SJ_pj.tmp -v1 | cut -f 2- -d ' ' | sort -k1,1 -k2,2n | uniq | awk '{print $0 " " "0"}' > NCLsampleinter_RSEM_pj.tmp5
        cat NCLsampleinter_RSEM_pj.tmp4 NCLsampleinter_RSEM_pj.tmp5 | sort -k1,1 -k2,2n | uniq > NCLsampleinter_RSEM_pj.txt        

        cat NCLsampleinter_RSEM_pj.txt | awk '{print $8 "\t" $0}' | sort -k1,1 | uniq > NCLsampleinter_RSEM_pj.txt.tmp1
        join  NCLsampleinter_RSEM_pj.txt.tmp1 gene_pmed.tmp | cut -f 2- -d ' '| sort -k1,1 -k2,2n | uniq > NCLsampleinter_RSEM_pj_pmed.tmp1
        join  NCLsampleinter_RSEM_pj.txt.tmp1 gene_pmed.tmp -v1 | cut -f 2- -d ' '| sort -k1,1 -k2,2n | uniq | awk '{print $0 " " "0"}' > NCLsampleinter_RSEM_pj_pmed.tmp2
        cat NCLsampleinter_RSEM_pj_pmed.tmp1 NCLsampleinter_RSEM_pj_pmed.tmp2 | sort -k1,1 -k2,2n | uniq > NCLsampleinter_RSEM_pj_pmed.tmp3
        cat NCLsampleinter_RSEM_pj_pmed.tmp3 | awk '{print $9 "\t" $0}' | sort -k1,1 | uniq > NCLsampleinter_RSEM_pj_pmed.tmp4
        join  NCLsampleinter_RSEM_pj_pmed.tmp4 gene_pmed.tmp | cut -f 2- -d ' '| sort -k1,1 -k2,2n | uniq > NCLsampleinter_RSEM_pj_pmed.tmp5
        join  NCLsampleinter_RSEM_pj_pmed.tmp4 gene_pmed.tmp -v1 | cut -f 2- -d ' '| sort -k1,1 -k2,2n | uniq | awk '{print $0 " " "0"}' > NCLsampleinter_RSEM_pj_pmed.tmp6
        cat NCLsampleinter_RSEM_pj_pmed.tmp5 NCLsampleinter_RSEM_pj_pmed.tmp6 | sort -k1,1 -k2,2n | uniq > NCLsampleinter_RSEM_pj_pmed.txt

        if [[ -n "$SCE" ]]; then 
           #echo "Step - to determine donor/acceptor side within SCE"
           paste <(cat NCLsampleinter_RSEM_pj_pmed.txt | awk '{print $1 "\t" $2-1 "\t" $2}')  <(cat NCLsampleinter_RSEM_pj_pmed.txt | tr ' ' '_')  > NCLsampleinter_RSEM_pj_pmed_SCE.tmp1
           bedtools intersect -a NCLsampleinter_RSEM_pj_pmed_SCE.tmp1 -b $SCE -wa | sort | uniq | awk '{print $4"_1"}' > NCLsampleinter_RSEM_pj_pmed_SCE.tmp2
           bedtools intersect -a NCLsampleinter_RSEM_pj_pmed_SCE.tmp1 -b $SCE -v | sort | uniq | awk '{print $4"_0"}' > NCLsampleinter_RSEM_pj_pmed_SCE.tmp3
           cat NCLsampleinter_RSEM_pj_pmed_SCE.tmp2 NCLsampleinter_RSEM_pj_pmed_SCE.tmp3 | tr '_' ' ' > NCLsampleinter_RSEM_pj_pmed_SCE.tmp4
           paste <(cat NCLsampleinter_RSEM_pj_pmed_SCE.tmp4 | awk '{print $4 "\t" $5-1 "\t" $5}')  <(cat NCLsampleinter_RSEM_pj_pmed_SCE.tmp4 | tr ' ' '_')  > NCLsampleinter_RSEM_pj_pmed_SCE.tmp5
           bedtools intersect -a NCLsampleinter_RSEM_pj_pmed_SCE.tmp5 -b $SCE -wa  | sort | uniq | awk '{print $4"_1"}' > NCLsampleinter_RSEM_pj_pmed_SCE.tmp6
           bedtools intersect -a NCLsampleinter_RSEM_pj_pmed_SCE.tmp5 -b $SCE -v  | sort | uniq| awk '{print $4"_0"}' > NCLsampleinter_RSEM_pj_pmed_SCE.tmp7
           cat NCLsampleinter_RSEM_pj_pmed_SCE.tmp6 NCLsampleinter_RSEM_pj_pmed_SCE.tmp7 | tr '_' ' ' > NCLsampleinter_RSEM_pj_pmed_SCE.txt
           echo "Chr_donor Pos_donor Strand_donor Chr_acceptor Pos_acceptor Strand_acceptor JunctionRead Gene_donor Gene_acceptor LinearRead_donor LinearRead_acceptor RPM_totalRaed RPM_mappedRead FPKM_donor TPM_donor FPKM_acceptor TPM_acceptor pj_dornor pj_acceptor pmedian_donor pmedian_acceptor SCE_donor SCE_acceptor" > HeaderPostInter.txt
           cat HeaderPostInter.txt NCLsampleinter_RSEM_pj_pmed_SCE.txt | tr ' ' \\t > inter_Characteristics/$getOneName.final
       fi

       if [[ -z "$SCE" ]]; then 
          echo "Chr_donor Pos_donor Strand_donor Chr_acceptor Pos_acceptor Strand_acceptor JunctionRead Gene_donor Gene_acceptor LinearRead_donor LinearRead_acceptor RPM_totalRaed RPM_mappedRead FPKM_donor TPM_donor FPKM_acceptor TPM_acceptor pj_dornor pj_acceptor pmedian_donor pmedian_acceptor" > HeaderPostInter.txt
          cat HeaderPostInter.txt NCLsampleinter_RSEM_pj_pmed.txt | tr ' ' \\t > inter_Characteristics/$getOneName.final
       fi
   done

fi


rm -r -f $BASEDIR\/comparison\/exons_boundary.bed
rm -r -f $BASEDIR\/comparison\/exons_boundary.tmp1
rm -r -f $BASEDIR\/comparison\/exons_boundary.tmp2
rm -r -f $BASEDIR\/comparison\/exons.txt
rm -r -f $BASEDIR\/comparison\/intra_tmp1.bed
rm -r -f $BASEDIR\/comparison\/intra_tmp2.bed
rm -r -f $BASEDIR\/comparison\/inter_tmp1.bed
rm -r -f $BASEDIR\/comparison\/inter_tmp2.bed 
rm -r -f $BASEDIR\/comparison\/intra_tmp.result
rm -r -f $BASEDIR\/comparison\/inter_tmp.result
rm -r -f $BASEDIR\/comparison\/sampleinter.tmp
rm -r -f $BASEDIR\/comparison\/sampleintra.tmp


rm -r -f NCLsample*
rm -r -f *tmp*
rm -r -f OC.final
rm -r -f *pre*
rm -r -f HeaderPost.txt
rm -r -f HeaderPostInter.txt


