######################
### building index ###
######################

mkdir STAR_RSEM_index
cd STAR_RSEM_index
rsem-prepare-reference --gtf /home/sctech/Desktop/test_NCLcomparator/Homo_sapiens.GRCh38.87.gtf --star -p 20 /home/sctech/Desktop/test_NCLcomparator/Homo_sapiens.GRCh38.dna.primary_assembly.fa RSEM

##########################
### run NCLcomparator  ###
##########################

./NCLcomparator.sh -gtf '/media/cswu/My Passport/toolsDev/NCLcomparator/Homo_sapiens.GRCh38.87.gtf' -thread 6 -index STAR_RSEM_index -read1 '/media/cswu/My Passport/toolsDev/HeLa_reads/SRR1637089_1.fastq.gz' -read2 '/media/cswu/My Passport/toolsDev/HeLa_reads/SRR1637089_2.fastq.gz' -intra '/media/cswu/My Passport/toolsDev/test_NCLcomparator/example/5col/intra' -inter '/media/cswu/My Passport/toolsDev/test_NCLcomparator/example/5col/inter' -sce '/media/cswu/My Passport/toolsDev/test_NCLcomparator/example/SCE/SCE_hg38.bed' -o HeLa  

