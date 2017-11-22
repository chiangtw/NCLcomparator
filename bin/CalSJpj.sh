#!/usr/bin/env zsh

if [[ -z "$1" ]]; then
   echo ""
   echo "Usage:"
   echo "./CalSJpj.sh [NCLin] [GTF] [SJ.out.tab]"
   echo "Output files: SJ_pj.txt and genes_pmed.txt"
   echo ""
   echo "ERROR: no NCL input file !"
   echo ""
   exit
fi

if [[ -z "$2" ]]; then
   echo ""
   echo "Usage:"
   echo "./CalSJpj.sh [NCLin] [GTF] [SJ.out.tab]"
   echo "Output files: SJ_pj.txt and genes_pmed.txt"
   echo ""
   echo "ERROR: no annotation gtf file !"
   echo ""
   exit
fi


if [[ -z "$3" ]]; then
   echo ""
   echo "Usage:"
   echo "./CalSJpj.sh [NCLin] [GTF] [SJ.out.tab]"
   echo "Output files: SJ_pj.txt and genes_pmed.txt"
   echo ""
   echo "ERROR: no SJ.out.tab !"
   echo ""
   exit
fi


NCLresult=$1
geneGTF=$2
starSJ=$3

cat $NCLresult |  awk '{print $1"_"$2 "\t" $7}' >SJ_h.tmp1 
cat $NCLresult |  awk '{print $4"_"$5 "\t" $7}' >SJ_h.tmp2 
cat SJ_h.tmp1 SJ_h.tmp2 | sort -k1,1 | grep ',' | sed 's/,/\t/g' | awk '{for(i=2; i<=NF;i++){print $1 "\t" $i}}' | sort -k1,1 | uniq > SJ_h.tmp3
cat SJ_h.tmp1 SJ_h.tmp2 | sort -k1,1 | grep -v ',' > SJ_h.tmp4
cat SJ_h.tmp3 SJ_h.tmp4 | sort -k1,1 | uniq > SJ_h.tmp 
cat SJ_h.tmp | awk '{print $2}' | sort  | uniq > SJ_h_Genename.tmp 

chrCheck=$(cat $starSJ | grep 'chr')
if [[ -n "$chrCheck" ]]
then
      cat $starSJ | awk '{print $1 "\t" $2-1 "\t" $7}' | sort -k1,1 -k2,2n | bedtools groupby -g 1,2 -c 3 -o sum > SJ_exon.pre1
      cat $starSJ | awk '{print $1 "\t" $3+1 "\t" $7}' | sort -k1,1 -k2,2n | bedtools groupby -g 1,2 -c 3 -o sum > SJ_exon.pre2
      cat SJ_exon.pre1 SJ_exon.pre2 | awk '{print $1"_"$2 "\t" $3}' | sort -k1,1 -k2,2n | bedtools groupby -g 1 -c 2 -o sum > SJ_exon.count      

else
      cat $starSJ | awk '{print "chr"$1 "\t" $2-1 "\t" $7}' | sort -k1,1 -k2,2n | bedtools groupby -g 1,2 -c 3 -o sum > SJ_exon.pre1
      cat $starSJ | awk '{print "chr"$1 "\t" $3+1 "\t" $7}' | sort -k1,1 -k2,2n | bedtools groupby -g 1,2 -c 3 -o sum > SJ_exon.pre2
      cat SJ_exon.pre1 SJ_exon.pre2 | awk '{print $1"_"$2 "\t" $3}' | sort -k1,1 -k2,2n | bedtools groupby -g 1 -c 2 -o sum > SJ_exon.count      
     
fi 


echo -n > SJ_pj.txt
echo -n > gene_pmed.txt

split -l 1000 SJ_h_Genename.tmp SJ_h_Genename.split.
for CHUNK in SJ_h_Genename.split.*;do
 
row_num=$(cat $CHUNK | wc -l)
echo $row_num
cat $geneGTF | grep -f $CHUNK > SJ_gtf.tmp
     
    for i in $(seq 1 $row_num)
    do
       SJ_pjCheckOne=$(cat "$CHUNK"| awk 'NR==k {print $0}' k=$i) 
       cat SJ_h.tmp | grep "$SJ_pjCheckOne$" > SJ_One.tmp
       cat SJ_gtf.tmp | grep "\"$SJ_pjCheckOne\"\;" | awk '{print $1"_"$4 "\t" $1"_"$5 "\t" $10}' | sed 's/;//g' | sed 's/"//g'| sort  -k1,1 | uniq  > getw_accdon.tmp
       cat getw_accdon.tmp | awk '{print $1 "\t" $3}' > getw_acc.tmp
       cat getw_accdon.tmp | awk '{print $2 "\t" $3}' > getw_don.tmp

       chrCheck2=$(cat getw_acc.tmp | grep 'chr')
       if [[ -n "$chrCheck2" ]]
       then
            cat getw_acc.tmp getw_don.tmp | sort -k1,1 | uniq > getOnew_SJ.tmp
       else
            cat getw_acc.tmp getw_don.tmp | awk '{print "chr"$0}' | sort -k1,1 | uniq > getOnew_SJ.tmp
       fi
  
       join -o 1.1 1.2 2.2 getOnew_SJ.tmp SJ_exon.count > getOnew_SJ_count.tmp
       #cat getOnew_SJ_count.tmp | awk '$3>0 {print $0}' > getOnew_SJ_countB0.tmp
       #cat getOnew_SJ_count.tmp | awk '$3==0 {print $1 "\t" $2 "\t" $3+1}' > getOnew_SJ_countZ.tmp
       #cat getOnew_SJ_countB0.tmp getOnew_SJ_countZ.tmp > getOnew_SJ_countAdj.tmp
       getOneSizew=$(cat getOnew_SJ_count.tmp | wc -l)
         if [ $getOneSizew -gt 0 ]    
          then  
            getOneSumw=$(cat getOnew_SJ_count.tmp | awk '{sum +=$3}END {print sum}')
            if [ $getOneSumw -gt 0 ]
             then  
               pmedtmp=$(cat getOnew_SJ_count.tmp | awk '{print $3/s}' s=$getOneSumw | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')
               cat getOnew_SJ_count.tmp | awk '{print $1 "\t" $2 "\t" $3/s}' s=$getOneSumw | sort -k1,1 > pjw.tmp
            fi
         fi
       #pjwSize=$(cat pjw.tmp | wc -l)
       
       if [[ -e pjw.tmp ]]  
       then       
           join SJ_One.tmp pjw.tmp >> SJ_pj.txt
           echo $SJ_pjCheckOne $pmedtmp >> gene_pmed.txt
       fi  
    done
done
rm -r -f SJ_h_Genename.split.*
#rm -r -f *.tmp*
#rm -r -f *.pre*
