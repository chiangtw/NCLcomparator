######################
### building index ###
######################

#mkdir STAR_RSEM_index
#cd STAR_RSEM_index
#rsem-prepare-reference --gtf Homo_sapiens.GRCh38.87.gtf --star -p 6 Homo_sapiens.GRCh38.dna.primary_assembly.fa RSEM

##########################
### run NCLcomparator  ###
##########################

#./NCLcomparator.sh -gtf Homo_sapiens.GRCh38.87.gtf -thread 6 -index STAR_RSEM_index -read1 SRR1637089_1.fastq.gz -read2 SRR1637089_2.fastq.gz -intra example/5col/intra -inter example/5col/inter -sce example/SCE/SCE_hg38.bed -o HeLa  

./NCLcomparator_new.sh -gtf Homo_sapiens.GRCh38.87.gtf \
 -genome GRCh38.p10.genome.fa \
 -trpts /home/cswu/Desktop/NCLcompartor_revise/NCLcomparator/NCLcomparator_v0.01/trpts \
 -thread 6 \
 -index /home/cswu/Desktop/NCLcomparator/STAR_RSEM_index \
 -read1 /home/cswu/Desktop/NCLcomparator/HeLa_reads/SRR1637089_1.fastq.gz \
 -read2 /home/cswu/Desktop/NCLcomparator/HeLa_reads/SRR1637089_2.fastq.gz \
 -intra /home/cswu/Desktop/NCLcompartor_revise/NCLcomparator/NCLcomparator_v0.01/example/5col/intra \
 -inter /home/cswu/Desktop/NCLcompartor_revise/NCLcomparator/NCLcomparator_v0.01/example/5col/inter \
 -sce/home/cswu/Desktop/NCLcompartor_revise/NCLcomparator/NCLcomparator_v0.01/example/SCE \
 -o HeLa
 
