#!/usr/bin/env zsh
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
     -gtf|--annotation) 
     geneGTF=$(readlink -f $2)
     shift
     ;; 
     -genome|--genome) 
     genome=$(readlink -f $2)
     shift
     ;; 
     -trpts|--trpts) 
     trpts=$(readlink -f $2)
     shift
     ;; 
     -thread|--ThreadN) 
     Thread=$2
     shift
     ;; 
     -index|--STAR_RSEMindex) 
     SRindex=$(readlink -f $2)
     shift
     ;; 
     -read1|--read1) 
     read1=$(readlink -f $2)
     shift
     ;; 
     -read2|--read2) 
     read2=$(readlink -f $2)
     shift
     ;;
     -intra|--cicular) 
     CIRCULAR=$(readlink -f $2)
     shift
     ;;
     -inter|--fusion)
     FUSION=$(readlink -f $2)
     shift
     ;;
     -sce|--SCE)
     SCE=$(readlink -f $2)
     shift
     ;;
     -o|--outputName)
     OPN=$2
     shift
     ;;
     *)

esac
shift
done



if [[ -z "$geneGTF" ]]; then
   echo ""
   echo "Usage:"
   echo "./NCLcomparator.sh -gtf [annotation gtf file] -genome [genome annotation] -trpts [transcripts folder] -thread [number of thread] -index [STAR_RSEM_index folder] -read1 [read1 fastq.gz file] -read2 [read2 fastq.gz file]  -intra [circular results folder] -inter [fusion results folder] -sce [sce bed file] -o [ output Name] "
   echo ""
   echo "ERROR: no annotation gtf file !"
   echo ""
   exit
fi

if [[ -z "$genome" ]]; then
   echo ""
   echo "Usage:"
   echo "./NCLcomparator.sh -gtf [annotation gtf file] -genome [genome annotation] -trpts [transcripts folder] -thread [number of thread] -index [STAR_RSEM_index folder] -read1 [read1 fastq.gz file] -read2 [read2 fastq.gz file]  -intra [circular results folder] -inter [fusion results folder] -sce [sce bed file] -o [ output Name] "
   echo ""
   echo "ERROR: no genome file !"
   echo ""
   exit
fi

if [[ -z "$trpts" ]]; then
   echo ""
   echo "Usage:"
   echo "./NCLcomparator.sh -gtf [annotation gtf file] -genome [genome annotation] -trpts [transcripts folder] -thread [number of thread] -index [STAR_RSEM_index folder] -read1 [read1 fastq.gz file] -read2 [read2 fastq.gz file]  -intra [circular results folder] -inter [fusion results folder] -sce [sce bed file] -o [ output Name] "
   echo ""
   echo "ERROR: no trpts folder !"
   echo ""
   exit
fi

if [[ -z "$SRindex" ]]; then
   echo ""
   echo "Usage:"
   echo "./NCLcomparator.sh -gtf [annotation gtf file] -genome [genome annotation] -trpts [transcripts folder] -thread [number of thread] -index [STAR_RSEM_index folder] -read1 [read1 fastq.gz file] -read2 [read2 fastq.gz file]  -intra [circular results folder] -inter [fusion results folder] -sce [sce bed file] -o [ output Name] "
   echo ""
   echo "ERROR: no STAR_RSEM_index folder !"
   echo ""
   exit
fi

if [[ -z "$read1" ]]; then
   echo ""
   echo "Usage:"
   echo "./NCLcomparator.sh -gtf [annotation gtf file] -thread [number of thread] -index [STAR_RSEM_index folder] -read1 [read1 fastq.gz file] -read2 [read2 fastq.gz file]  -intra [circular results folder] -inter [fusion results folder] -sce [sce bed file] -o [ output Name] "
   echo ""
   echo "ERROR: no read1 fastq.gz file !"
   echo ""
   exit
fi

if [[ -z "$read2" ]]; then
   echo ""
   echo "Usage:"
   echo "./NCLcomparator.sh -gtf [annotation gtf file] -genome [genome annotation] -trpts [transcripts folder] -thread [number of thread] -index [STAR_RSEM_index folder] -read1 [read1 fastq.gz file] -read2 [read2 fastq.gz file]  -intra [circular results folder] -inter [fusion results folder] -sce [sce bed file] -o [ output Name] "  
   echo ""
   echo "ERROR: no read2 fastq.gz file !"
   echo ""
   exit
fi


if [[ -z "$Thread" ]]; then
   Thread=6
fi

if [[ -z "$OPN" ]]; then
   OPN=$(echo "out")
fi

BASEDIR=$(dirname $(readlink -f "$0"))
BASEDIR1=$(pwd)
rm -r -f $BASEDIR1\/$OPN
mkdir $BASEDIR1\/$OPN
mkdir $BASEDIR1\/$OPN\/STAR_RSEM_out
cd $BASEDIR1\/$OPN\/STAR_RSEM_out
STAR --chimSegmentMin 10 --runThreadN $Thread --genomeDir $SRindex --readFilesIn $read1 $read2 --readFilesCommand zcat --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM
rm -r -f Aligned.out.bam
rsem-calculate-expression -p $Thread --paired-end --bam Aligned.toTranscriptome.out.bam $SRindex\/RSEM RSEMout
echo ""
date "+%d/%m/%y %H:%M:%S ...... begining NCLcomparator analysis!"

totalRead=$(cat $BASEDIR1\/$OPN\/STAR_RSEM_out\/Log.final.out | grep 'Number of input reads' | awk '{print $6}')
UMRead=$(cat $BASEDIR1\/$OPN\/STAR_RSEM_out\/Log.final.out | grep 'Uniquely mapped reads number'| awk '{print $6}')

rm -r -f $BASEDIR1\/$OPN\/comparison
mkdir $BASEDIR1\/$OPN\/comparison

## build exon boundary ##
echo ""
date "+%d/%m/%y %H:%M:%S ...... building exon boundaries"
cd $BASEDIR1\/$OPN\/comparison

chrCheck=$(cat $geneGTF | grep -v "#" | head -1 | awk '{print $1}' | grep 'chr')

if [[ -n "$chrCheck" ]]
then
    cat $geneGTF | awk '{for(i=1;i<=NF;i++) if($i=="exon_id"){print  $(i+1) "\t" $0} }' | awk '{for(i=1;i<=NF;i++) if($i=="gene_name"){print  $(i+1) "\t" $0} }' | awk '{for(i=1;i<=NF;i++) if($i=="gene_id"){print  $(i+1) "\t" $0} }' | awk '{print $3 "\t" $2 "\t" $1 "\t" $4 "\t" $7 "\t" $8 "\t" $10}' | sed 's/;//g' | sed 's/"//g' | awk '{print $4 "\t" $5 "\t" $6 "\t" $7 "\t" $2":"$3}' | sort -k1,1 -k2,2n | uniq > exons.txt
    cat exons.txt | awk '{print $1 "\t" $2-6 "\t" $2+5 "\t" $1":"$2":"$4":"$5}' | sort -k1,1 -k2,2n | uniq > exons_boundary.tmp1
    cat exons.txt | awk '{print $1 "\t" $3-6 "\t" $3+5 "\t" $1":"$3":"$4":"$5}' | sort -k1,1 -k2,2n | uniq > exons_boundary.tmp2
    cat exons_boundary.tmp1 exons_boundary.tmp2 | sort -k1,1 -k2,2n | uniq  > exons_boundary.bed
else
    cat $geneGTF | awk '{for(i=1;i<=NF;i++) if($i=="exon_id"){print  $(i+1) "\t" $0} }' | awk '{for(i=1;i<=NF;i++) if($i=="gene_name"){print  $(i+1) "\t" $0} }' | awk '{for(i=1;i<=NF;i++) if($i=="gene_id"){print  $(i+1) "\t" $0} }' | awk '{print $3 "\t" $2 "\t" $1 "\t" "chr"$4 "\t" $7 "\t" $8 "\t" $10}' | sed 's/;//g' | sed 's/"//g' | awk '{print $4 "\t" $5 "\t" $6 "\t" $7 "\t" $2":"$3}' | sort -k1,1 -k2,2n | uniq > exons.txt
    cat exons.txt | awk '{print $1 "\t" $2-6 "\t" $2+5 "\t" $1":"$2":"$4":"$5}' | sort -k1,1 -k2,2n | uniq > exons_boundary.tmp1   
    cat exons.txt | awk '{print $1 "\t" $3-6 "\t" $3+5 "\t" $1":"$3":"$4":"$5}' | sort -k1,1 -k2,2n | uniq > exons_boundary.tmp2
    cat exons_boundary.tmp1 exons_boundary.tmp2 | sort -k1,1 -k2,2n | uniq  > exons_boundary.bed
fi

cat $geneGTF | awk '{for(i=1;i<=NF;i++) if($i=="gene_name"){print  $(i+1) "\t" $0} }'  | awk '{for(i=1;i<=NF;i++) if($i=="gene_id") {print $(i+1) "\t" $0}}' | awk '{print $1 "\t" $2}' | sed 's/;//g' | sed 's/"//g' | sort -k1,1 | uniq > ENSG_GeneID_convert.txt
cat $geneGTF | awk '{for(i=1;i<=NF;i++) if($i=="exon_number") {print $(i+1) "\t" $0} }'| awk '{for(i=1;i<=NF;i++) if($i=="gene_name"){print  $(i+1) "\t" $0} }' | awk '{print $1 "\t" $2}' | sed 's/;//g' | sed 's/"//g'| sort -k1,1 -k2,2nr | awk '!x[$1]++' | sort -k1,1 > GeneID_exonNum.txt

 
cat $BASEDIR1\/$OPN\/STAR_RSEM_out\/RSEMout.genes.results | awk '{print $1 "\t" $6 "\t" $7}' | sed -e '1d' > genes_RSEM.txt
join ENSG_GeneID_convert.txt genes_RSEM.txt > ENSG_GeneID_RSEM.txt  
## ENSG versus GeneID is multiple to one, choose choose large FPKM as represent of GeneID ## 
cat ENSG_GeneID_RSEM.txt | awk '{print $2 "\t" $3 "\t" $4}' |  sort -k2,2nr | awk '!x[$1]++' | sort -k1,1 > ENSG_GeneID_RSEM.txt.tmp 
 

if [[ -n "$CIRCULAR" ]]; then
   
   echo ""
  
   date "+%d/%m/%y %H:%M:%S ...... intra-Comparation analysis"

   rm -r -f  $BASEDIR1\/$OPN\/comparison\/intra
   mkdir $BASEDIR1\/$OPN\/comparison\/intra

   cd $BASEDIR1\/$OPN\/comparison
   ## intra results ##
   echo "Step: to adjust the junction coordinates of intra results into exon boundaries"
   ls $CIRCULAR > sampleintra.tmp 
   sampleintra_num=$(cat sampleintra.tmp | wc -l)
   
   for i in $(seq 1 "$sampleintra_num")
   do 
     getOne=$(cat sampleintra.tmp | awk 'NR==k {print $1}' k=$i)
     #echo $getOne
     getOneName=$(echo $getOne | sed 's/[.]/\t/g' | awk '{print $1}') 
     cat $CIRCULAR/$getOne | sort -k1,1 -k2,2n | awk '$1~/^chr[0-9XY]*$/ {print $0}' | awk '$3~/^chr[0-9XY]*$/ {print $0}' | awk '$2~/^[0-9]*$/ {print $0}'| awk '$4~/^[0-9]*$/ {print $0}' | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $1":"$2":"$3":"$4":"$5}' | sort -k1,1 -k2,2n | uniq  > intra_tmp1.bed
     bedtools intersect -a intra_tmp1.bed -b exons_boundary.bed -wa -wb | awk '{print $4 "\t" $8}' | sed 's/:/\t/g' | awk '{print $3 "\t" $4-1 "\t" $4 "\t" $1":"$2":"$3":"$4":"$5":"$6":"$7":"$8":"$9}' | sort -k1,1 -k2,2n | uniq > intra_tmp2.bed
     bedtools intersect -a intra_tmp2.bed -b exons_boundary.bed  -wa -wb | awk '{print $4 "\t" $8}' | sed 's/:/\t/g' | awk '{print $6 "\t" $7 "\t" $8 "\t" $10 "\t" $11 "\t" $12 "\t" $5 "\t" $9 "\t" $13}'| sort -k1,1 -k2,2n | uniq | awk '$3==$6 && $8==$9 {print $0}'  > intra_tmp.result
     cat intra_tmp.result | sed 's/+/sense/g' | sed 's/-/anti/g' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8}' | sort -k1,1 -k2,2n | uniq | bedtools groupby -g 1,2,3,4,5,6,7 -c 8 -o collapse | sed 's/sense/+/g' | sed 's/anti/-/g' | awk '{print $0 "\t" $8}'| sort -k1,1 -k2,2n | uniq > intra_tmp1.result
     cat intra_tmp1.result | awk '{if($3=="+" && $2 > $5) print $0; else if($3=="+" && $2 < $5) print $1 "\t" $5 "\t" $3 "\t" $4 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9; else if($3=="-" && $2 < $5) print $0; else if($3=="-"&& $2> $5) print $1 "\t" $5 "\t" $3 "\t" $4 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9}'| sort -k1,1 -k2,2n | uniq > intra_tmp2.result
     cat intra_tmp2.result | sed 's/+/sense/g' | sed 's/-/anti/g' | awk '{print $1":"$2":"$3":"$4":"$5":"$6":"$8 "\t" $7}' | sort -k1,1 | bedtools groupby -g 1 -c 2 -o sum > intra_tmp3.result  
     cat intra_tmp3.result | sed 's/:/\t/g'| sed 's/sense/+/g' | sed 's/anti/-/g' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $8 "\t" $7 "\t" $7 }'> intra\/$getOneName.circEB 
 
   done

   ## merge intra ##
    echo "Step: to merge intra results to intraMerged.result"
    cd  $BASEDIR1\/$OPN


    cat comparison\/intra\/*.circEB | sort -k1,1 -k2,2n | uniq  > all.circEB
    cat all.circEB | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" $8}' | sort -k1,1 | uniq > all.circEB.tmp 

    ls comparison\/intra > sampleintra 
    sampleintra_num=$(cat sampleintra | wc -l)
    echo "Number of intra tools: "$sampleintra_num

    header=$(echo "Chr" "\t" "Pos" "\t" "Strand" "\t" "Chr" "\t" "Pos" "\t" "Strand" "\t" "Genename")  
    sampleintra_num=$(echo $sampleintra_num)
    for i in $(seq 1 $sampleintra_num)
    do 
       getOne=$(cat sampleintra | awk 'NR==k {print $1}' k=$i)
       echo $i"." $getOne
       getOneName=$(echo $getOne | sed 's/[.]/\t/g' | awk '{print $1}')
       header=$(echo $header "\t" $getOneName)
       cat comparison\/intra\/$getOne | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" $7}' | sort -k1,1 > $getOneName.circEB.tmp
       join -o 1.1 1.2 2.2 all.circEB.tmp  $getOneName.circEB.tmp -a1 -e0 | awk '{print $1 "\t" $2":"$3}' | sort -k1,1 | uniq > all.circEB_2.tmp
       cat all.circEB_2.tmp > all.circEB.tmp
    done

    cat all.circEB.tmp  | sed 's/:/\t/g' | tr ' ' \\t > $BASEDIR1\/$OPN\/comparison\/intraMerged.out
    cat <(echo $header) $BASEDIR1\/$OPN\/comparison\/intraMerged.out > intraMerged.result
    
    echo "Step: to graph intraMerged.result"
    # run graphic_intra.R
    # output file : intra.pdf and intraMerged_junction.result
    $BASEDIR\/bin\/graphic_intra.R intraMerged.result intra
    rm -r -f *.tmp
    rm -r -f all.circEB
    rm -r -f sampleintra
    
    echo ""
    date "+%d/%m/%y %H:%M:%S ...... checking intra ambiguous alignment"
    # run ambiguous alignment: checkAA.sh
    # output file: AA_intraMerged_junction.result
    mkdir $BASEDIR1\/$OPN\/checkAA_intra_out
    mv $BASEDIR1\/$OPN\/intraMerged_junction.result $BASEDIR1\/$OPN\/checkAA_intra_out
    cd $BASEDIR1\/$OPN\/checkAA_intra_out
    $BASEDIR\/bin\/checkAA.sh -NCLevents intraMerged_junction.result -thread $Thread -genome $genome -trpts $trpts -tools $BASEDIR\/bin
    mv AA_intraMerged_junction.result $BASEDIR1\/$OPN
   
   cd $BASEDIR1\/$OPN\/comparison
   cat intraMerged.out | awk '{print $7 "\t" $1 ":" $2 ":" $3 ":" $4 ":" $5 ":" $6 }' | sort -k1,1 > intraMerged.tmp1
   echo "Step: to add the exon number of its host gene"
   ## to add exon number ## 
   join -o 1.1 1.2 2.2 intraMerged.tmp1 GeneID_exonNum.txt -a1 -e nd  | sort -k1,1 | uniq | tr ' ' \\t  > intraMerged.tmp.exonNum
   echo "Step: to add  TPM and FPKM of its host gene"
   ## to add TPM and FPKM ##
   join -o 1.1 1.2 2.2 2.3 intraMerged.tmp1 ENSG_GeneID_RSEM.txt.tmp -a1 -e nd | sort -k1,1 | uniq | tr ' ' \\t  > intraMerged.tmp.TPM_FPKM
   
   
   echo "Step: to calculate SJpj using CalSJpj.sh"
   ### generate SJ_pj.txt and gene_pmed.txt ##
   $BASEDIR/bin/CalSJpj.sh intraMerged.out $geneGTF $BASEDIR1\/$OPN\/STAR_RSEM_out\/SJ.out.tab
   cat SJ_pj.txt > intra_SJ_pj.txt
   cat gene_pmed.txt > intra_gene_pmed.txt
   
   echo "Step: to examine out of circular using OutOfCircular.sh"
   ## generate OC.final ##
   $BASEDIR/bin/OutOfCircular.sh intraMerged.out $BASEDIR1\/$OPN\/STAR_RSEM_out\/Chimeric.out.junction $BASEDIR1\/$OPN\/STAR_RSEM_out\/Chimeric.out.sam $BASEDIR1\/$OPN\/STAR_RSEM_out\/Log.final.out
   
   
   join -o 1.1 1.2 1.3 2.3 2.4 intraMerged.tmp.exonNum intraMerged.tmp.TPM_FPKM -a1 -e nd | sort | uniq | tr ' ' \\t | awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\t" $5}' | tr ':' \\t > intraMerged.tmp2
   paste <(cat intraMerged.tmp2 | awk '{print $1"_"$2}') <(cat intraMerged.tmp2 | tr '\t' ':') | sort -k1,1  > intraMerged.tmp3.1
   join -o 1.2 2.2 intraMerged.tmp3.1 SJ_exon.count -a1 -e nd | tr ' ' \\t | tr ':' \\t | sort | uniq > intraMerged.tmp3.2
   paste <(cat intraMerged.tmp3.2 | awk '{print $4"_"$5}') <(cat intraMerged.tmp3.2 | tr '\t' ':') | sort -k1,1 | uniq > intraMerged.tmp3.3  
   join -o 1.2 2.2 intraMerged.tmp3.3 SJ_exon.count -a1 -e nd | tr ' ' \\t | tr ':' \\t | sort  | uniq > intraMerged.tmp3
   
   cat intra_SJ_pj.txt | awk '{print $1"_"$2 "\t" $4}' | sort | uniq > intra_SJ_pj.tmp1
   paste <(cat intraMerged.tmp3 | awk '{print $1"_"$2"_"$7}') <(cat intraMerged.tmp3 | tr '\t' ':') | sort -k1,1 -k2,2n | uniq > intraMerged.tmp4.1
   join -o 1.2 2.2 intraMerged.tmp4.1 intra_SJ_pj.tmp1 -a1 -e nd | sort -k1,1 -rnk2 | awk '!x[$1]++' | tr ' ' \\t | tr ':' \\t > intraMerged.tmp4.2
   paste <(cat intraMerged.tmp4.2 | awk '{print $4"_"$5"_"$7 }') <(cat intraMerged.tmp4.2 | tr '\t' ':' ) | sort -k1,1 -k2,2n | uniq  > intraMerged.tmp4.3    
   join -o 1.2 2.2 intraMerged.tmp4.3 intra_SJ_pj.tmp1 -a1 -e nd | tr ':' \\t | tr ' ' \\t | sort | uniq > intraMerged.tmp4.4
   paste <(cat intraMerged.tmp4.4 | awk '{print $7}') <(cat intraMerged.tmp4.4 | tr '\t' ':') | sort -k1,1  | uniq > intraMerged.tmp4.5
   join -o 1.2 2.2 intraMerged.tmp4.5 intra_gene_pmed.txt -a1 -e nd | tr ' ' \\t | tr ':' \\t | sort  | uniq > intraMerged.tmp4
   paste <(cat intraMerged.tmp4 | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6}') <(cat intraMerged.tmp4 | tr '\t' ':') | sort -k1,1 | uniq > intraMerged.tmp5.1
   cat OC.final | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6 "\t" $8}' | sort -k1,1 | uniq > OC.final.tmp1
   join -o 1.2 2.2 intraMerged.tmp5.1 OC.final.tmp1 -a1 -e nd | tr ' ' \\t | tr ':' \\t > intraMerged.tmp5

   if [[ -z "$SCE" ]]; then 
       echo "Chr_donor Pos_donor Strand_donor Chr_acceptot Pos_acceptor Strand_acceptor Gene_name Number_of_Exon TPM FPKM LinearRead_donor LinearRead_acceptor P_D P_A P_median Out_of_Circle" | tr ' ' \\t > Header.tmp1  
    fi


   if [[ -n "$SCE" ]]; then 
      echo "Step: to determine donor/acceptor side within SCE"
      paste <(cat intraMerged.out | awk '{print $1 "\t" $2-1 "\t" $2}')  <(cat intraMerged.out | awk '{print $1":"$2":"$3":"$4":"$5":"$6}')  > intraMerged_SCE.tmp1
      bedtools intersect -a intraMerged_SCE.tmp1 -b $SCE -wa | sort | uniq | awk '{print $4":1"}' > intraMerged_SCE.tmp2
      bedtools intersect -a intraMerged_SCE.tmp1 -b $SCE -v | sort | uniq | awk '{print $4":0"}' > intraMerged_SCE.tmp3
      cat intraMerged_SCE.tmp2 intraMerged_SCE.tmp3 | tr ':' '\t' > intraMerged_SCE.tmp4
      paste <(cat intraMerged_SCE.tmp4 | awk '{print $4 "\t" $5-1 "\t" $5}')  <(cat intraMerged_SCE.tmp4 | tr '\t' ':')  > intraMerged_SCE.tmp5
      bedtools intersect -a intraMerged_SCE.tmp5 -b $SCE -wa  | sort | uniq | awk '{print $4":1"}' > intraMerged_SCE.tmp6
      bedtools intersect -a intraMerged_SCE.tmp5 -b $SCE -v  | sort | uniq| awk '{print $4":0"}' > intraMerged_SCE.tmp7
      cat intraMerged_SCE.tmp6 intraMerged_SCE.tmp7 | tr ':' '\t'  | sort -k1,1 -k2,2n > intraMerged_tmp_SCE.txt
      paste <(cat intraMerged.tmp5 | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6}') <(cat intraMerged.tmp5 | tr '\t' ':') | sort -k1,1 | uniq > intraMerged.tmp6.1
      cat intraMerged_tmp_SCE.txt | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6 "\t" $7":"$8 }' | sort -k1,1 | uniq > intraMerged_SCE.txt.tmp1
      join -o 1.2 2.2 intraMerged.tmp6.1 intraMerged_SCE.txt.tmp1 -a1 -e nd | tr ' ' \\t | tr ':' \\t > intraMerged.tmp5
      echo "Chr_donor Pos_donor Strand_donor Chr_acceptot Pos_acceptor Strand_acceptor Gene_name Number_of_Exon TPM FPKM LinearRead_donor LinearRead_acceptor P_D P_A P_median Out_of_Circle SCE_D SCE_A" | tr ' ' \\t > Header.tmp1  
    fi

   echo "Step: to calculate RPM, CF and RNCL in each tool "
   
   cat intraMerged.tmp3 | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" $11 "\t" $12}' | sort -k1,1 > intraMerged_LR.tmp
   cat intraMerged.tmp5 > intraMerged.tmp5.1
   cat Header.tmp1 > Header.tmp1.1
   for i in $(seq 1 $sampleintra_num)
   do 
      cat $BASEDIR1\/$OPN\/intraMerged.result | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $(7+k)}' k=$i > intraMerged_toolOne.tmp1
      toolOneName=$(cat intraMerged_toolOne.tmp1 | head -n1  | awk '{print $7}') 
      #echo $toolOneName
      cat intraMerged_toolOne.tmp1 | sed -e '1d' | awk '$7 > 0 {print $1":"$2":"$3":"$4":"$5":"$6 "\t" $7}' | sort -k1,1 | uniq > intraMerged_toolOne.tmp2
      cat intraMerged_toolOne.tmp2 | awk '{print $1 "\t" ($2/T)*1000000 "\t" ($2/M)*1000000}' T=$totalRead  M=$UMRead | sort -k1,1 > intraMerged_toolOne.tmp3
      join -o 1.1 1.2 2.2 2.3 intraMerged_toolOne.tmp2 intraMerged_LR.tmp -a1 -e nd | sort -k1,1 | uniq > intraMerged_toolOne.tmp4.1
      cat intraMerged_toolOne.tmp4.1 | sed 's/nd/0/g' | awk '{print $1 "\t" $2/($2+$3+$4+1) "\t" (2*$2)/(2*$2+$3+$3)}' | sort -k1,1  | uniq  > intraMerged_toolOne.tmp4
      
      join -o 1.1 1.2 2.2 2.3 intraMerged_toolOne.tmp2 intraMerged_toolOne.tmp3 -a1 -e nd | tr ' ' \\t > intraMerged_toolOne.tmp5.1
      join -o 1.1 1.2 1.3 1.4 2.2 2.3 intraMerged_toolOne.tmp5.1 intraMerged_toolOne.tmp4 -a1 -e nd |tr ' ' \\t > intraMerged_toolOne.tmp5
      
      paste <(cat intraMerged.tmp5.1 | awk '{print $1":"$2":"$3":"$4":"$5":"$6}') <(cat intraMerged.tmp5.1 | tr '\t' ':') | sort -k1,1 | uniq > intraMerged.tmp6 
      cat intraMerged_toolOne.tmp5 | awk '{print $1 "\t" $2":"$3":"$4":"$5":"$6}' | sort -k1,1 | uniq  > intraMerged_toolOne.tmp6.1
      join -o 1.2 2.2 intraMerged.tmp6 intraMerged_toolOne.tmp6.1 -a1 -e nd:nd:nd:nd:nd | sort -k1,1 | tr ' ' \\t | tr ':' \\t >  intraMerged.tmp7
      echo $toolOneName"_JunctionRead" $toolOneName"_RPM_total" $toolOneName"_RPM_mapped" $toolOneName"_CF" $toolOneName"_RNCL" | tr ' ' \\t > toolHeader.tmp
      paste Header.tmp1.1 toolHeader.tmp | tr ' ' \\t > Header.tmp2
      cat Header.tmp2 > Header.tmp1.1
      cat intraMerged.tmp7 > intraMerged.tmp5.1
  done 
  cat Header.tmp2 intraMerged.tmp7 > $BASEDIR1\/$OPN\/intraMerged_characteristic.result

fi

if [[ -n "$FUSION" ]]; then

    echo ""
    date "+%d/%m/%y %H:%M:%S ...... inter-Comparation analysis"

    rm -r -f  $BASEDIR1\/$OPN\/comparison\/inter
    mkdir $BASEDIR1\/$OPN\/comparison\/inter

    cd $BASEDIR1\/$OPN\/comparison
    ## inter results ##
    echo "Step: to adjust the junction coordinates of inter results into exon boundaries"
    ls $FUSION > sampleinter.tmp 
    sampleinter_num=$(cat sampleinter.tmp | wc -l)
    for i in $(seq 1 "$sampleinter_num")
    do 
     getOne=$(cat sampleinter.tmp | awk 'NR==k {print $1}' k=$i)
     #echo $getOne
     getOneName=$(echo $getOne | sed 's/[.]/\t/g' | awk '{print $1}') 
     cat $FUSION/$getOne | sort -k1,1 -k2,2n | awk '$1~/^chr[0-9XY]*$/ {print $0}' | awk '$3~/^chr[0-9XY]*$/ {print $0}' | awk '$2~/^[0-9]*$/ {print $0}'| awk '$4~/^[0-9]*$/ {print $0}' | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $1":"$2":"$3":"$4":"$5}' | sort -k1,1 -k2,2n | uniq  > inter_tmp1.bed
     bedtools intersect -a inter_tmp1.bed -b exons_boundary.bed -wa -wb | awk '{print $4 "\t" $8}' | sed 's/:/\t/g' | awk '{print $3 "\t" $4-1 "\t" $4 "\t" $1":"$2":"$3":"$4":"$5":"$6":"$7":"$8":"$9}' | sort -k1,1 -k2,2n | uniq > inter_tmp2.bed
     bedtools intersect -a inter_tmp2.bed -b exons_boundary.bed  -wa -wb | awk '{print $4 "\t" $8}' | sed 's/:/\t/g' | awk '{print $6 "\t" $7 "\t" $8 "\t" $10 "\t" $11 "\t" $12 "\t" $5 "\t" $9 "\t" $13}'| sort -k1,1 -k2,2n | uniq | awk '$8!=$9 {print $0}' > inter_tmp.result
     cat inter_tmp.result | sed 's/+/sense/g' | sed 's/-/anti/g' | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" $7}' | sort -k1,1 -k2,2n | uniq | bedtools groupby -g 1 -c 2 -o sum | sort | uniq > inter_tmp1.result
     cat inter_tmp.result | sed 's/+/sense/g' | sed 's/-/anti/g' | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" $8}' | sort -k1,1 -k2,2n | uniq | bedtools groupby -g 1 -c 2 -o collapse | sort | uniq > inter_tmp2.result
     cat inter_tmp.result | sed 's/+/sense/g' | sed 's/-/anti/g' | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" $9}' | sort -k1,1 -k2,2n | uniq | bedtools groupby -g 1 -c 2 -o collapse | sort | uniq > inter_tmp3.result
     cat inter_tmp.result | sed 's/+/sense/g' | sed 's/-/anti/g' | awk '{print $1":"$2":"$3":"$4":"$5":"$6 }' | sort | uniq > inter_tmp4.result
     join -o 1.1 2.2 inter_tmp4.result inter_tmp1.result -a1 -e nd | tr ' ' \\t | sort  | uniq > inter_tmp5.result
     join -o 1.1 1.2 2.2 inter_tmp5.result inter_tmp2.result -a1 -e nd | tr ' ' \\t | sort  | uniq > inter_tmp6.result
     join -o 1.1 1.2 1.3 2.2 inter_tmp6.result inter_tmp3.result -a1 -e nd | tr ' ' \\t | sort  | sed 's/sense/+/g' | sed 's/anti/-/g' | sed 's/:/\t/g' > inter\/$getOneName.fusionEB  
    done

    ## merge inter ##
    echo "Step: to merge inter results to interMerged.result"
    cd  $BASEDIR1\/$OPN
    
    cat comparison/inter/*.fusionEB | sort -k1,1 -k2,2n | uniq  > all.fusionEB
    cat all.fusionEB | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" $8":"$9}' | sort -k1,1 | uniq > all.fusionEB.tmp 
    
    ls comparison/inter > sampleinter 
    sampleinter_num=$(cat sampleinter | wc -l)
    echo "Number of inter tools: "$sampleinter_num
    header=$(echo "Chr" "\t" "Pos" "\t" "Strand" "\t" "Chr" "\t" "Pos" "\t" "Strand" "\t" "Genename" "\t" "Genename")
    sampleinter_num=$(echo $sampleinter_num)
    for i in $(seq 1 $sampleinter_num)
    do 
        getOne=$(cat sampleinter | awk 'NR==k {print $1}' k=$i)
        echo $i"." $getOne
        getOneName=$(echo $getOne | sed 's/[.]/\t/g' | awk '{print $1}')
        header=$(echo $header "\t" $getOneName)
        cat comparison/inter/$getOne | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" $7}' | sort -k1,1 > $getOneName.fusionEB.tmp
        join -o 1.1 1.2 2.2 all.fusionEB.tmp  $getOneName.fusionEB.tmp -a1 -e0 | awk '{print $1 "\t" $2":"$3}' | sort -k1,1 | uniq > all.fusionEB_2.tmp
        cat all.fusionEB_2.tmp > all.fusionEB.tmp
    done
    
    cat all.fusionEB.tmp  | sed 's/:/\t/g' | tr ' ' \\t > $BASEDIR1\/$OPN\/comparison\/interMerged.out
    cat <(echo $header) $BASEDIR1\/$OPN\/comparison\/interMerged.out > interMerged.result
    
    echo "Step: to graph interMerged.result"
    # run graphic_intra.R
    # output file : intra.pdf and intraMerged_junction.result
    $BASEDIR\/bin\/graphic_inter.R interMerged.result inter
    rm -r -f *.tmp
    rm -r -f all.fusionEB
    rm -r -f sampleinter
     
    echo ""
    date "+%d/%m/%y %H:%M:%S ...... checking inter ambiguous alignment"
    # run ambiguous alignment: checkAA.sh
    # output file: AA_interMerged_junction.result
    mkdir $BASEDIR1\/$OPN\/checkAA_inter_out
    mv $BASEDIR1\/$OPN\/interMerged_junction.result $BASEDIR1\/$OPN\/checkAA_inter_out
    cd $BASEDIR1\/$OPN\/checkAA_inter_out
    $BASEDIR\/bin\/checkAA.sh -NCLevents interMerged_junction.result -thread $Thread -genome $genome -trpts $trpts -tools $BASEDIR\/bin
    mv AA_interMerged_junction.result $BASEDIR1\/$OPN
   


    cd $BASEDIR1\/$OPN\/comparison
    cat interMerged.out | awk '{print $7 "\t" $1":"$2":"$3":"$4":"$5":"$6":"$7":"$8 }' | sort -k1,1 | uniq > interMerged.tmp1
    echo "Step: to add  TPM and FPKM of its host gene"
    ## to add TPM and FPKM ##
    join -o 1.2 2.2 2.3 interMerged.tmp1 ENSG_GeneID_RSEM.txt.tmp -a1 -e nd | sort | uniq | tr ' ' \\t  > interMerged.tmp1.TPM_FPKM
    paste <(cat interMerged.tmp1.TPM_FPKM | tr ':' \\t | awk '{print $8}') <(cat interMerged.tmp1.TPM_FPKM | tr '\t' ':')  | sort  -k1,1 > interMerged.tmp2.TPM_FPKM
    join -o 1.2 2.2 2.3 interMerged.tmp2.TPM_FPKM ENSG_GeneID_RSEM.txt.tmp -a1 -e nd | sort | uniq | tr ' ' \\t  | tr ':' \\t > interMerged.tmp.TPM_FPKM
    
    echo "Step: to calculate SJpj using CalSJpj.sh"
    ## generate SJ_pj.txt and gene_pmed.txt ##
    $BASEDIR/bin/CalSJpj.sh interMerged.out $geneGTF $BASEDIR1\/$OPN\/STAR_RSEM_out\/SJ.out.tab
    cat SJ_pj.txt > inter_SJ_pj.txt
    cat gene_pmed.txt > inter_gene_pmed.txt
     
    paste <(cat interMerged.tmp.TPM_FPKM | awk '{print $1"_"$2}') <(cat interMerged.tmp.TPM_FPKM | tr '\t' ':') | sort -k1,1 | uniq > interMerged.tmp2.1
    join -o 1.2 2.2 interMerged.tmp2.1 SJ_exon.count -a1 -e nd | tr ' ' \\t | tr ':' \\t | sort | uniq > interMerged.tmp2.2
    paste <(cat interMerged.tmp2.2 | awk '{print $4"_"$5}') <(cat interMerged.tmp2.2 | tr '\t' ':') | sort -k1,1 | uniq > interMerged.tmp2.3
    join -o 1.2 2.2 interMerged.tmp2.3 SJ_exon.count -a1 -e nd | tr ' ' \\t | tr ':' \\t | sort | uniq > interMerged.tmp2
    
    cat inter_SJ_pj.txt | awk '{print $1"_"$2 "\t" $4}' | sort | uniq > inter_SJ_pj.tmp1
    paste <(cat interMerged.tmp2 | awk '{print $1"_"$2"_"$7}') <(cat interMerged.tmp2 | tr '\t' ':') | sort -k1,1 > interMerged.tmp3.1 
    join -o 1.2 2.2 interMerged.tmp3.1 inter_SJ_pj.tmp1 -a1 -e nd | sort -k1,1 -k2,2nr | awk '!x[$1]++' | tr ' ' \\t | tr ':' \\t | sort | uniq > interMerged.tmp3.2
    paste <(cat interMerged.tmp3.2 | awk '{print $4"_"$5"_"$8}') <(cat interMerged.tmp3.2 | tr '\t' ':') | sort -k1,1 | uniq > interMerged.tmp3.3
    join -o 1.2 2.2 interMerged.tmp3.3 inter_SJ_pj.tmp1 -a1 -e nd | sort -k1,1 -k2,2nr | awk '!x[$1]++' | tr ' ' \\t | tr ':' \\t | sort | uniq > interMerged.tmp3
    paste <(cat interMerged.tmp3 | awk '{print $7}') <(cat interMerged.tmp3 | tr '\t' ':') | sort -k1,1 | uniq > interMerged.tmp4.1
    join -o 1.2 2.2 interMerged.tmp4.1 inter_gene_pmed.txt -a1 -e nd  | tr ' ' \\t | tr ':' \\t |  sort -k1,1 | uniq > interMerged.tmp4.2
    paste <(cat interMerged.tmp4.2 | awk '{print $8}') <(cat interMerged.tmp4.2 | tr '\t' ':') | sort -k1,1 | uniq > interMerged.tmp4.3
    join -o 1.2 2.2 interMerged.tmp4.3 inter_gene_pmed.txt -a1 -e nd | tr ' ' \\t | tr ':' '\t' | sort -k1,1 | uniq > interMerged.tmp4

     if [[ -z "$SCE" ]]; then 
       echo "Chr_donor Pos_donor Strand_donor Chr_acceptot Pos_acceptor Strand_acceptor GeneName_donor GeneName_acceptor TPM_donor FPKM_donor TPM_acceptor FPKM_acceptor LinearRead_donor LinearRead_acceptor P_D P_A Pmedian_donor Pmedian_acceptor" | tr ' ' \\t > Header.tmp1  
    fi


   if [[ -n "$SCE" ]]; then 
      echo "Step: to determine donor/acceptor side within SCE"
      paste <(cat interMerged.out | awk '{print $1 "\t" $2-1 "\t" $2}')  <(cat interMerged.out | awk '{print $1":"$2":"$3":"$4":"$5":"$6}')  > interMerged_SCE.tmp1
      bedtools intersect -a interMerged_SCE.tmp1 -b $SCE -wa | sort | uniq | awk '{print $4":1"}' > interMerged_SCE.tmp2
      bedtools intersect -a interMerged_SCE.tmp1 -b $SCE -v | sort | uniq | awk '{print $4":0"}' > interMerged_SCE.tmp3
      cat interMerged_SCE.tmp2 interMerged_SCE.tmp3 | tr ':' '\t' > interMerged_SCE.tmp4
      paste <(cat interMerged_SCE.tmp4 | awk '{print $4 "\t" $5-1 "\t" $5}')  <(cat interMerged_SCE.tmp4 | tr '\t' ':')  > interMerged_SCE.tmp5
      bedtools intersect -a interMerged_SCE.tmp5 -b $SCE -wa  | sort | uniq | awk '{print $4":1"}' > interMerged_SCE.tmp6
      bedtools intersect -a interMerged_SCE.tmp5 -b $SCE -v  | sort | uniq| awk '{print $4":0"}' > interMerged_SCE.tmp7
      cat interMerged_SCE.tmp6 interMerged_SCE.tmp7 | tr ':' '\t'  | sort -k1,1 -k2,2n > interMerged_tmp_SCE.txt
      paste <(cat interMerged.tmp4 | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6}') <(cat interMerged.tmp4 | tr '\t' ':') | sort -k1,1 | uniq > interMerged.tmp5.1
      cat interMerged_tmp_SCE.txt | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6 "\t" $7":"$8 }' | sort -k1,1 | uniq > interMerged_SCE.txt.tmp1
      join -o 1.2 2.2 interMerged.tmp5.1 interMerged_SCE.txt.tmp1 -a1 -e nd | tr ' ' \\t | tr ':' \\t > interMerged.tmp4
      echo "Chr_donor Pos_donor Strand_donor Chr_acceptot Pos_acceptor Strand_acceptor GeneName_donor GeneName_acceptor TPM_donor FPKM_donor TPM_acceptor FPKM_acceptor LinearRead_donor LinearRead_acceptor P_D P_A Pmedian_donor Pmedian_acceptor SCE_D SCE_A" | tr ' ' \\t > Header.tmp1  
   fi

   echo "Step: to calculate RPM in each tool " 
   cat interMerged.tmp4 > interMerged.tmp4.1
   cat Header.tmp1 > Header.tmp1.1
   for i in $(seq 1 $sampleinter_num)
   do 
      cat $BASEDIR1\/$OPN\/interMerged.result | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $(8+k)}' k=$i > interMerged_toolOne.tmp1
      toolOneName=$(cat interMerged_toolOne.tmp1 | head -n1  | awk '{print $7}') 
      echo $toolOneName
      cat interMerged_toolOne.tmp1 | sed -e '1d' | awk '$7 > 0 {print $1":"$2":"$3":"$4":"$5":"$6 "\t" $7}' | sort -k1,1 | uniq > interMerged_toolOne.tmp2
      cat interMerged_toolOne.tmp2 | awk '{print $1 "\t" $2 "\t" ($2/T)*1000000 "\t" ($2/M)*1000000}' T=$totalRead  M=$UMRead | sort -k1,1 > interMerged_toolOne.tmp3
      paste <(cat interMerged.tmp4.1 | awk '{print $1":"$2":"$3":"$4":"$5":"$6}') <(cat interMerged.tmp4.1 | tr '\t' ':') | sort -k1,1 | uniq > interMerged.tmp5 
      join -o 1.2 2.2 2.3 2.4 interMerged.tmp5 interMerged_toolOne.tmp3 -a1 -e nd | sort -k1,1 | tr ' ' \\t | tr ':' \\t >  interMerged.tmp6
      echo $toolOneName"_JunctionRead" $toolOneName"_RPM_total" $toolOneName"_RPM_mapped" | tr ' ' \\t > toolHeader.tmp
      paste Header.tmp1.1 toolHeader.tmp | tr ' ' \\t > Header.tmp2
      cat Header.tmp2 > Header.tmp1.1
      cat interMerged.tmp6 > interMerged.tmp4.1
  done 
  cat Header.tmp2 interMerged.tmp6 > $BASEDIR1\/$OPN\/interMerged_characteristic.result

fi


cd $BASEDIR1\/$OPN\/comparison
rm  -r -f *tmp*        
rm -r -f SJ_pj.txt
rm -r -f gene_pmed.txt
rm -r -f *pre* 
cd $BASEDIR1\/$OPN\/STAR_RSEM_out
rm -r -f Aligned.toTranscriptome.out.bam
rm -r -f RSEMout.transcript.bam

echo ""
date "+%d/%m/%y %H:%M:%S ...... NCLcomparator analysis completed !"

