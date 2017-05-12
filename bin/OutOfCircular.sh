#!/usr/bin/env zsh 

NCLIn=$1
FChimericOutJunc=$2
FSAM=$3
starlog=$4

echo -n > output.OR
echo -n > NCL.SRR

cat $NCLIn | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' | sort -k1,1 -k2,2n > NCL_6col.tmp

cat NCL_6col.tmp | grep '+' | awk '{print $1 "\t" $2+1 "\t" $3 "\t" $4 "\t" $5-1 "\t" $6}' > NCL.pos.pre1
cat NCL_6col.tmp | grep '-' | awk '{print $1 "\t" $2-1 "\t" $3 "\t" $4 "\t" $5+1 "\t" $6}' > NCL.pos.pre2
cat NCL.pos.pre1 NCL.pos.pre2 | sort -k1,1 -k2,2n > NCL_atIntronPos

RL=$(cat $starlog | grep 'Average input read length'  | awk '{print $6/2}')

split -l 100 NCL_atIntronPos NCL_atIntronPos.split.
for CHUNK in NCL_atIntronPos.split.*;do
    cat $FChimericOutJunc | grep -f <(cat "$CHUNK" | awk '{print $2"\t"}') | grep -f <(cat "$CHUNK" | awk '{print $5"\t"}')| grep -f <(cat "$CHUNK" | awk '{print $1"\t"}')>> NCL.SRR
done
rm -r -f NCL_atIntronPos.split.*


sort -k1,1 $FSAM > tmp_sorted.sam
#join tmp_sorted.sam <(cat NCL.SRR | awk '{print $10}' | sort -k1,1 )  | awk '$2=="163" || $2=="99" {print $0}' > check_All.sam
#join tmp_sorted.sam <(cat NCL.SRR | awk '{print $10}' | sort -k1,1 )  > check_All.sam
join tmp_sorted.sam <(cat NCL.SRR | awk '{print $10}' | sort -k1,1 )  | awk '$2=="83" || $2=="163" || $2=="337" || $2=="99" || $2=="147" || $2=="401" || $2=="385" {print $0}' > check_All.sam

declare -i row_num=$(cat NCL_6col.tmp | wc -l)
for  i in $(seq 1 $row_num)
do
   cat NCL_atIntronPos | awk 'NR==k {print $0}' k=$i > checkOne.pos 
   check_chr=$(cat checkOne.pos | awk '{print $1}')
   check_donor=$(cat checkOne.pos | awk '{print $2}')
   check_acceptor=$(cat checkOne.pos | awk '{print $5}')
   cat NCL.SRR | awk '$1==var1 {print $0}' var1=$check_chr | awk '$4==var1 {print $0}' var1=$check_chr | awk '$2==var2 || $5==var2 {print $0}' var2=$check_donor | awk '$2==var3 || $5==var3 {print $0}' var3=$check_acceptor | awk '{print $10}'| sort -k1,1 > checkOne.SRR 
   
   join check_All.sam checkOne.SRR > check_OR.sam
   JS_num=$(cat checkOne.SRR | wc -l)
 
   if [[ $check_donor -gt $check_acceptor ]]
    then
     OR1_num=$(cat check_OR.sam | awk '$3==var1 {print $0}' var1=$check_chr | awk '$4 >=var2 {print $0}' var2=$check_donor var4=$RL | wc -l)
     OR2_num=$(cat check_OR.sam | awk '$3==var1 {print $0}' var1=$check_chr | awk '($4+var4) <=var3 {print $0}' var3=$check_acceptor var4=$RL | wc -l)      
     OR1_SRR=$(cat check_OR.sam | awk '$3==var1 {print $0}' var1=$check_chr | awk '$4 >=var2 {print $0}' var2=$check_donor var4=$RL | awk '{print $1}' | sort -k1,1 | tr "\n\r" ",")     
     OR2_SRR=$(cat check_OR.sam | awk '$3==var1 {print $0}' var1=$check_chr | awk '($4+var4) <=var3 {print $0}' var3=$check_acceptor var4=$RL | awk '{print $1}' | sort -k1,1| tr "\n\r" ",")    
 
    else
     OR1_num=$(cat check_OR.sam | awk '$3==var1 {print $0}' var1=$check_chr | awk '($4+var4) <=var2 {print $0}' var2=$check_donor var4=$RL | wc -l)
     OR2_num=$(cat check_OR.sam | awk '$3==var1 {print $0}' var1=$check_chr | awk '$4 >=var3 {print $0}' var3=$check_acceptor var4=$RL | wc -l)
     OR1_SRR=$(cat check_OR.sam | awk '$3==var1 {print $0}' var1=$check_chr | awk '($4+var4) <=var2 {print $0}' var2=$check_donor var4=$RL | awk '{print $1}' | sort -k1,1 | tr "\n\r" ",")     
     OR2_SRR=$(cat check_OR.sam | awk '$3==var1 {print $0}' var1=$check_chr | awk '$4 >=var3 {print $0}' var3=$check_acceptor var4=$RL | awk '{print $1}' | sort -k1,1 | tr "\n\r" ",")     
   fi  

     OR_num=$(($OR1_num + $OR2_num))
     OR_SRR="$OR1_SRR$OR2_SRR"
     OR="$(cat checkOne.pos) $JS_num $OR_num $OR_SRR"
     
     echo $OR | sed 's/\r//g' >> output.OR   
done

cat output.OR | grep '+' | awk '{print $1 "\t" $2-1 "\t" $3 "\t" $4 "\t" $5+1 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' > output.OR1
cat output.OR | grep '-' | awk '{print $1 "\t" $2+1 "\t" $3 "\t" $4 "\t" $5-1 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' > output.OR2
cat output.OR1 output.OR2 | sort -k1,1 -k2,2n > OC.final


rm -r -f checkOne.pos
rm -r -f checkOne.SRR
rm -r -f check_OR.sam
rm -r -f NCL.pos.pre1
rm -r -f NCL.pos.pre2
rm -r -f NCL.SRR
rm -r -f NCL_atIntronPos
rm -r -f output.OR
rm -r -f output.OR1
rm -r -f output.OR2
rm -r -f tmp_sorted.sam 
rm -r -f check_All.sam
rm -r -f NCL_6col.tmp

