#!/usr/bin/env zsh
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
     -NCLevents) 
     NCLevents=$(readlink -f $2)
     shift
     ;; 
     -thread) 
     thread=$2
     shift
     ;; 
     -genome) 
     genome=$(readlink -f $2)
     shift
     ;; 
     -trpts) 
     trpts=$(readlink -f $2)
     shift
     ;; 
     -tools)
     tools_bin=$(readlink -f $2)
     shift
     ;;
     
     *)

esac
shift
done

if [[ -z "$NCLevents" ]]; then
   echo ""
   echo "Usage:"
   echo "./checkAA.sh -NCLevents [circRNAs.txt] -thread [number of thread] -genome [genome.fa] -trpts [transcriptome folder] -tools [bin path]"
   echo ""
   exit
fi

nameout=$(echo $NCLevents | sed 's/\//\n/g' | tail -n 1)
echo "Step1: fetch flanking +- 100 read"
cat $NCLevents | cut -d$'\t' -f 7- --complement | tr \\t ':' > Merged_pos.tmp
cat $NCLevents | cut -d$'\t' -f 7- | tr \\t ':' > Merged_juncs.tmp
paste Merged_pos.tmp Merged_juncs.tmp > Merged_pos_juncs.tmp

cat Merged_pos.tmp |sed -e '1d' | sed 's/:/\t/g' | awk '$3=="+" {print $1":"$2":"$3":"$4":"$5":"$6".1" "\t" $1 "\t" $2-99 "\t" $2 "\t" $3}' > NCL_Dplus.tmp
cat Merged_pos.tmp |sed -e '1d' | sed 's/:/\t/g' | awk '$3=="-"{print $1":"$2":"$3":"$4":"$5":"$6".1" "\t" $1 "\t" $2 "\t" $2+99 "\t" $3}' > NCL_Dminus.tmp
cat Merged_pos.tmp |sed -e '1d' | sed 's/:/\t/g' | awk '$6=="+" {print $1":"$2":"$3":"$4":"$5":"$6".2" "\t" $4 "\t" $5 "\t" $5+99 "\t" $6}' > NCL_Aplus.tmp
cat Merged_pos.tmp |sed -e '1d' | sed 's/:/\t/g' | awk '$6=="-" {print $1":"$2":"$3":"$4":"$5":"$6".2" "\t" $4 "\t" $5-99 "\t" $5 "\t" $6}' > NCL_Aminus.tmp

cat NCL_Dplus.tmp NCL_Aplus.tmp NCL_Dminus.tmp NCL_Aminus.tmp |sort -k1,1 > NCL_200_interval.txt
cat NCL_200_interval.txt |awk '{print $2 "\t" $3-1 "\t" $4 "\t" $1 "\t" "1" "\t" $5}' > NCL_200_interval.bed

chrCheck=$(cat $genome | head -n 1 | awk '{print $1}' | grep 'chr')
if [[ -n "$chrCheck" ]]
then
     cat $genome > genome.fa
else
     cat $genome | sed 's/>\([0-9XYMT]\)/>chr\1/' | sed 's/dna.*//'  > genome.fa
fi

bedtools getfasta -fi genome.fa -bed NCL_200_interval.bed -s -name -fo NCL_200_interval.fa
$tools_bin/merge_paired_sequences.py NCL_200_interval.fa  NCL_200.merged.fa


echo ""
echo  "Step2: runs blat"
mkdir BLAT_tmp
echo "chrM" > chrM.list
cat $genome | $tools_bin/SeqOut chrM.list 1 > chrM.fa
cat <(echo ">RepChrM") <(cat chrM.fa | sed -e '1d' ) <(cat chrM.fa | sed -e '1d') > RepChrM.fa  
cat $trpts/*.fa RepChrM.fa > others.fa

$tools_bin/mp_blat.py $genome NCL_200.merged.fa NCL.rG.1.psl -p $thread --blat_bin blat --tmp_path BLAT_tmp
$tools_bin/mp_blat.py $genome NCL_200.merged.fa NCL.rG.2.psl -p $thread --blat_bin blat --blat_opt "-tileSize=9 -stepSize=9 -repMatch=32768" --tmp_path BLAT_tmp
$tools_bin/mp_blat.py others.fa NCL_200.merged.fa NCL.rO.1.psl -p $thread --blat_bin blat --tmp_path BLAT_tmp
$tools_bin/mp_blat.py others.fa NCL_200.merged.fa NCL.rO.2.psl -p $thread --blat_bin blat --blat_opt "-tileSize=9 -stepSize=9 -repMatch=32768" --tmp_path BLAT_tmp

echo ""
echo "Step3: evaluate blat results"
cat NCL.rG.1.psl | sed -e '1,5d' | awk 'substr($14,0,4)!~/^chr[0-9XY]*/' | awk '($1+$3)/$11 > 0.8' | awk '{print $10}' | sort  | uniq  > GL.rG.1.list
cat NCL.rG.1.psl | sed -e '1,5d' | awk '($1+$3)/$11 > 0.8 {print $10}' | sort  | uniq  > colinear.rG.1.list
cat NCL.rG.1.psl | sed -e '1,5d' | awk 'substr($14,0,4) ~/^chr[0-9XY]*/' | $tools_bin/PslChimeraFilter 30 3 > NCL.chi0.1.bed
cat NCL.rG.1.psl | sed -e '1,5d' | $tools_bin/RemoveInList 10 NCL.chi0.1.bed 4  | $tools_bin/RemoveInList 10 colinear.rG.1.list 1 | $tools_bin/RemoveInList 10 GL.rG.1.list 1 > NCL.rG.1_ambigous.psl 
cat NCL.rG.1_ambigous.psl | awk '{print $10}' | sort  | uniq -c | awk '$1=="1" {print $2}' > badBlat.rG.1.list
cat NCL.rG.1_ambigous.psl | $tools_bin/RemoveInList 10 badBlat.rG.1.list 1 | awk '{print $10}' | sort  | uniq > multipleHit.rG.1.list

cat NCL.rG.2.psl | sed -e '1,5d' | awk 'substr($14,0,4)!~/^chr[0-9XY]*/' | awk '($1+$3)/$11 > 0.8' | awk '{print $10}' | sort  | uniq  > GL.rG.2.list
cat NCL.rG.2.psl | sed -e '1,5d' | awk '($1+$3)/$11 > 0.8 {print $10}' | sort  | uniq > colinear.rG.2.list
cat NCL.rG.2.psl | sed -e '1,5d' | awk 'substr($14,0,4) ~/^chr[0-9XY]*/' | $tools_bin/PslChimeraFilter 30 3 > NCL.chi0.2.bed
cat NCL.rG.2.psl | sed -e '1,5d' | $tools_bin/RemoveInList 10 NCL.chi0.2.bed 4  | $tools_bin/RemoveInList 10 colinear.rG.2.list 1 | $tools_bin/RemoveInList 10 GL.rG.2.list 1 > NCL.rG.2_ambigous.psl 
cat NCL.rG.2_ambigous.psl | awk '{print $10}' | sort  | uniq -c | awk '$1=="1" {print $2}' > badBlat.rG.2.list
cat NCL.rG.2_ambigous.psl | $tools_bin/RemoveInList 10 badBlat.rG.2.list 1 | awk '{print $10}' | sort  | uniq > multipleHit.rG.2.list

cat NCL.rO.1.psl | sed -e '1,5d' | awk '($1+$3)/$11 > 0.8 {print $10}' | awk '{print $10}' | sort  | uniq  > colinear.rO.1.list
cat NCL.rO.2.psl | sed -e '1,5d' | awk '($1+$3)/$11 > 0.8 {print $10}' | sort  | uniq  > colinear.rO.2.list
cat GL.rG.1.list GL.rG.2.list | sort | uniq > GL.rG.list
cat colinear.rG.1.list colinear.rG.2.list colinear.rO.1.list colinear.rO.2.list GL.rG.list | sort  | uniq | sed '/^$/d'  > colinear.rGO.list
cat multipleHit.rG.1.list multipleHit.rG.2.list | sort  | uniq  > multipleHit.rG.list
cat colinear.rGO.list | sed '/^$/d'> colinear.clear.list
cat multipleHit.rG.list | sort  | uniq > multipleHit.rG.GL.list 
join multipleHit.rG.GL.list colinear.clear.list -v1 | sed '/^$/d' > multipleHit.clear.list
cat Merged_pos_juncs.tmp | sed -e '1d' | sort -k1,1 > Merged_pos_juncs_noheader.tmp
cat Merged_pos_juncs.tmp | head -n 1 | awk 'BEGIN{FS="\t"}{print $1 "\t" $2 "\t" "AA_colinear" "\t" "AA_multipleHit" }'> header_Merged_pos_juncs.tmp
cat colinear.clear.list  | awk '{print $0 "\t" "1"}' | sort  -k1,1 > Merged_colinear.tmp
cat multipleHit.clear.list  | awk '{print $0 "\t" "1"}' | sort  -k1,1 > Merged_multipleHit.tmp
join -o 1.1 1.2 2.2 Merged_pos_juncs_noheader.tmp Merged_colinear.tmp -a1 -e 0 | sort | uniq | tr ' ' \\t   > Merged_pos_juncs_noheader_colinear.tmp
join -o 1.1 1.2 1.3 2.2 Merged_pos_juncs_noheader_colinear.tmp Merged_multipleHit.tmp -a1 -e 0 | sort | uniq | tr ' ' \\t  > Merged_pos_juncs_noheader_colinear_mutipleHit.tmp

cat header_Merged_pos_juncs.tmp Merged_pos_juncs_noheader_colinear_mutipleHit.tmp | tr ':' \\t > AA_$nameout   

rm -r -f *.tmp*
rm -r -f others.fa
rm -r -f genome.fa
rm -r -f genome.fa.fai
rm -r -f RepChrM.fa
