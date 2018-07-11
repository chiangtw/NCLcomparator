# Manual of NCLcomparator
--------------------------
#### Version: 0.0.1
#### NCLcomparator, a comprehensive tool for analyzing non-co-linear (NCL) transcripts (fusion, trans-splicing, and circular RNA), enables to combine several NCL results from different detection tools, and to provide the characteristics of NCL events. 
------------------------------
The NCLcomparator program, document, and test set can be downloaded from our FTP site: ftp://treeslab1.genomics.sinica.edu.tw/NCLcomparator or GitHub: https://github.com/TreesLab/NCLcomparator.

### 1. System requirements
  NCLcomparator runs under Linux-like environment (i.e. Bio-Linux, also see http://environmentalomics.org/bio-linux/) with at least 30 GB RAM.
       
### 2. Installation
```sh
$ tar zxvf NCLcomparator.tar.gz
$ cd NCLcomparator
$ chmod +x NCLcomparator.sh
$ chmod +x bin/*
```
### 3. Installation external tools
   (1) bedtools (http://bedtools.readthedocs.io/en/latest/)
   
   (2) BLAT (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/)
   
   (3) STAR (https://github.com/alexdobin/STAR)
   
   (4) RSEM (https://github.com/deweylab/RSEM)
   
   (5) R (https://www.r-project.org/)
   
Bedtools was applied to manipulate genomic coordinate data (BED file). STAR was applied to align RNA-seq reads aganist reference genome for retrieving the number of mapped reads and the reads spanning NCL junctions to calculate the statistics of detected NCL events, such as RPM, RNCL, CF , PD , PA  and Pmedian. RSEM was applied to calculate the expression of  NCL events’ corresponding co-linear host genes (TPM and FPKM).  R was applied to graph the coverage of identified NCL events among the compared tools. BLAT was applied to identify false positive NCL events caused by ambiguous alignment originating from repetitive sequences or paralog genes.    
   

Get latest bedtools source from releases and install it 
```sh
$ wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools2.25.0.tar.gz
$ tar -zxvf bedtools-2.25.0.tar.gz
$ cd bedtools2
$ make
$ sudo cp ./bin/ /usr/local/bin
```
Get executable BLAT prgram
```sh
$ wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat
$ sudo cp blat /usr/local/bin
```
Get latest STAR source from releases and install it 
```sh
$ wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
$ tar -xzf 2.5.3a.tar.gz
$ cd STAR-2.5.3a
$ sudo cp bin/Linux_x86_64_static/STAR /usr/local/bin
```
Get latest RSEM source code and install it
```sh
$ wget https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz
$ tar -xzf v1.3.0.tar.gz
$ cd RSEM-1.3.0
$ make
$ make ebseq
$ make install
```
Get R in Ubuntu environment
```sh
$ sudo apt-get update
$ sudo apt-get install r-base
```
### 4. Preparation
 (1) Genome, transcritome and its annotation, which can be download from the GENCODE website (http://www.gencodegenes.org/) or ensembl FTP   (http://www.ensembl.org/info/data/ftp/index.html). Given Human as an example, go to ensembl FTP (http://www.ensembl.org/info/data/ftp/index.html) to download human genome and annotation.
 
```sh 
$ wget ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
$ wget ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz
$ gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
$ gunzip Homo_sapiens.GRCh38.87.gtf.gz

# download types of transctips files into the trpts folder
$ mkdir trpts
$ cd trpts
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.pc_transcripts.fa.gz
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.lncRNA_transcripts.fa.gz
$ gunzip gencode.v28.pc_transcripts.fa.gz
$ gunzip gencode.v28.lncRNA_transcripts.fa.gz
``` 
 (2) (Optional) Synonymous Constraint elements (SCE), which can be download from (http://compbio.mit.edu/SCE/)

 (3) Building STAR and RSEM index  
 
The following steps are to generate the index files before running NCLcomparator tool. 
```sh
$ mkdir STAR_RSEM_index
$ cd STAR_RSEM_index
$ rsem-prepare-reference --gtf /path/to/Homo_sapiens.GRCh38.87.gtf \
  --star -p 10 /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa RSEM
```
The STAR and RSEM index of the genome hg38 and the annotation ensemble 87 is prepared and can be downloaded in our FTP (ftp://treeslab1.genomics.sinica.edu.tw/NCLcomparator/STAR_RSEM_index.tar.gz). 

### 5. Execution of NCLcomparator
    
Usage:
 ```sh
 $ ./NCLcomparator.sh -gtf [annotation GTF file] \
  -genome [genome fasta file] \
  -trpts [the folder with transcripts] \
  -index [STAR_RSEM index folder] \
  -thread [number of thread] \
  -read1 [fastq_1.gz] \
  -read2 [fastq_2.gz] \
  -intra [circular result folder] \
  -inter [fusion result folder] \
  -sce [SCE bed file] \
  -o prefix of output folder
```
An example: 
```sh 
$ ./NCLcomparator.sh -gtf Homo_sapiens.GRCh38.87.gtf \
  -genome Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -trpts /path/to/trpts \
  -index /path/to/STAR_RSEM_index \
  -thread 6\
  -read1 HeLa_1.fastq.gz \
  -read2 HeLa_2.fastq.gz \
  -intra /path/to/example/5col/intra \
  -inter /path/to/example/5col/inter \
  -sce /path/to/example/SCE_hg38.bed \
  -o HeLa  
```
The basic options to run a job as follow:
- -trpts /path/to/the folder with the transcripts downloaded from different types or sources.
- -intra /path/to/NCL-intra result folder
- -inter (optional) /path/to/NCL-inter result folder
- -sce (optional) /path/to/SCE.bed

The output file of running HeLa rRNA depleted RNA-seq data are provided in the example_output/HeLa folder.

### 6. Input format of NCL events
The input format of NCL detection results are required to modify as 5-col format, which includes the positons of donor/acceptor sides and the number of supporting NCL-junction reads. Intra-NCL modified 5-col format results are gathered into a folder, and inter-NCL modified 5-col format results are gathered into another folder. An example of NCL detectors' result on HeLa rRNA depleted RNA-seq data are provided at the foler (example_output/HeLa) 
  
Given a 5-col result, CIRCexplorer2.5col (tab-delimited text file), as an example,  

| Chromosome | Coordination | Chromosome | Coordination | Total number of supporting NCL-junction reads | 
|------------|--------------|------------|--------------| -------------------------------|
| chr1	| 945176 | chr1	| 945517 |	1 |
| chr1	| 955922 | chr1	| 956013 |	1 |
| chr1	| 955922 | chr1	| 957273 |  5 |
     
### 7. Output files
After executing NCLcomparator program, two merged circular RNA and fusion RNA tools' results (AA_intraMerged_junction.result and AA_interMerged_junction.result) accompanied by its graphic report (intra.pdf and inter.pdf) are generated, and two newly folder is created (comparison and STAR_RSEM_out). The genomic positions of NCL events (within 5 bp of franking exon boundaries) in each NCL detection tool are adjusticed to the exact exonic boundaries. The other two characteristic files of NCL events (intraMerged_characteristic.result and interMerged_characteristic.result) are produced.  

In STAR_RSEM_out folder, several output files are generated by running STAR and RSEM programs, these files are given as the input files of characteristic process in NCLcomaprator. In comparison folder, the two folders (intra and inter) are generated, where NCL events of the positions in each tool with <= 5 bp of franking exon boundaries are adjusted to the exact positions of exon boundaries.



The column formats of outputs are described as follow:

#### AA_intraMerged_junction.result (tab-delimited text file)

|No. of column | Description |
|---------------|-------------|
| (1) | Chromosome name of the donor side (5' ss) |
| (2) | Junction coordination of the donor side |
| (3) | Strand of the donor side |
| (4) | Chromosome name of the acceptor side (3' ss) |
| (5) | Junction coordinate of the acceptor side |
| (6) | Strand of the acceptor side |
| (7) | Gene name |
| (8) | Circular tool No.1 (Yes: number of supporting intragenic NCL-junction reads; No: 0) |
| (9) | Circular tool No.2	(Yes: number of supporting intragenic NCL-junction reads; No: 0) |
|(10) | ... |
|( the last 6th ) | Number of supported circular tools in the event |
|( the last 5th ) | Median of junction reads |
|( the last 4th ) | Tau: evaluating variantion of identified junction reads among circular tools |
|( the last 3th) | NCLscore: confident level of the event |
|( the last 2th) | AA_colinear : ambiguous alignemnt diagnosis of NCL event, colinear effect: Yes(1); No(0) |
|( the last )| AA_multipleHit : ambiguous alignement diagnosis of NCL event. multiple hit effect: Yes(1); No(0) |  


#### AA_interMerged_junction.result (tab-delimited text file)

|No. of column | Description |
|---------------|-------------|
| (1) | Chromosome name of the donor side (5' ss) |
| (2) | Junction coordination of the donor side |
| (3) | Strand of the donor side |
| (4) | Chromosome name of the acceptor side (3' ss) |
| (5) | Junction coordinate of the acceptor side |
| (6) | Strand of the acceptor side |
| (7) | Gene name of the donor side |
| (8) | Gene nqamae of the acceptor side |
| (9) | Fusion tool No.1 (Yes: number of supporting intergenic NCL-junction reads; No: 0) |
| (10) | Fusion tool No.2	(Yes: number of supporting intergenic NCL-junction reads; No: 0) |
|(11) | ...|
|( the last 6th ) | Number of supported fusion fools in the event |
|( the last 5th ) | Median of junction reads |
|( the last 4th ) | Tau: evaluating variation of identified junction reads among fusion tools |
|( the last 3th) | NCLscore: confident level of the event |
|( the last 2th) | AA_colinear : ambiguous alignemnt diagnosis of NCL event, colinear effect: Yes(1); No(0) |
|( the last )| AA_multipleHit : ambiguous alignement diagnosis of NCL event. multiple hit effect: Yes(1); No(0) |  


The other two ouputs, intraMerged_characteristic.result and interMerged_characteristic.result, provide the useful features of NCL events.


#### intraMerged_characteristic.result (tab-delimited text file)

| No. of column | Description |
|---------------|-------------|
| (1) | Chromosome name of the donor side (5' ss) |
| (2) | Junction coordination of the donor site |
| (3) | Strand of the donor site |
| (4) | Chromosome name of the acceptor site (3' ss) |
| (5) | Junction coordinate of the acceptor site |
| (6) | Strand of the acceptor site |
| (7) | Gene name|
| (8) | Total number of exons in the gene |
| (9) | TPM of the gene |
| (10) | FPKM of the gene |
| (11) | Number of reads spanning the co-linearly spliced juncations at NCL donor site|
| (12) | Number of reads spanning the co-linearly splice junctions at NCL acceptor site |
| (13) | Usage of the co-linear junctions at NCL donor splice site (P_D) | 
| (14) | Usage of the co-linear junctions at NCL acceptor spliced site (P_A) | 
| (15) | Median frequency of occurrance of all well-annotated splice sites (co-linear) in the host gene (P_median)|
| (16) | Out of circle (Yes: 1; No: 0) |
| (17) | Circular tool No.1 : number of supporting NCL-junction reads |
| (18) | Circular tool No.1: RPM based on total RNA-seq reads |
| (19) | Circular tool No.1: RPM based on uniquely mapping reads |
| (20) | Circular tool No.1 : circular fraction (CF) |
| (21) | Circular tool No.1: non-co-linear ratio (R_NCL) |
| (22) | Circular tool No.2 : number of supporting NCL-junction reads |
| (23) | Circular tool No.2: RPM based on total RNA-seq reads |
| (24) | Circular tool No.2: RPM based on uniquely mapping reads |
| (25) | Circular tool No.2 : circular fraction (CF) |
| (26) | Circular tool No.2: non-co-linear ratio (R_NCL) |
|...| ...|
|(Optional)| The donor size within SCE (1) or outside SCE (0) |	
|(Optional)| The acceptor size within SCE (1) or outside SCE (0) |


#### interMerged_characteristic.result (tab-delimited text file)

| No. of column | Description |      
|---------------|-------------|
| (1) | Chromosome name of the donor side (5' ss) |
| (2) | Junction coordination of the donor side |
| (3) | Strand of the donor side |
| (4) | Chromosome name of the acceptor side (3' ss) |
| (5) | Junction coordinate of the acceptor side |
| (6) | Strand of the acceptor side |
| (7) | Gene name of the donor side |
| (8) | Gene name of the acceptor side |
| (9)| TPM of the donor gene |
| (10)| FPKM of the donor gene |
| (11)| TPM of the acceptor gene |
| (12)| FPKM of the acceptor gene |
| (13)| Number of reads spanning the co-linearly spliced junctions at NCL donor site |	
| (14)| Number of reads spanning the co-linearly spliced junctions at NCL acceptor site |
| (15)| Usage of co-linear juncations at NCL donor splice site in the donor gene (P_D) |
| (16)| Usage of co-linear junctions at NCL acceptor splice site in the acceptor gene (P_A)| 
| (17)| Median frequency of occurrence of well-annotated splice sites (co-linear) in the donor gene (P_median_D) |
| (18)| Median frequency of occurrence of well-annoataed splice sites (co-linear) in the acceptor gene (P_median_A) |
| (19) | Fusion tool No.1 : number of supporting intergenic NCL-junction reads |
| (20) | Fusion tool No.1 : RPM based on raw reads |
| (21) | Fusion tool No.1 : RPM based on uniquely mapped reads |
| (22) | Fusion tool No.2 : number of supporting intergenic NCL-junction reads |
| (23) | Fusion tool No.2 : RPM based on raw reads |
| (24) | Fusion tool No.2 : RPM based on uniquely mapped reads |
|...|...|
| (Optional) | The donor side within SCE (1) or outside SCE (0) |
| (Optional) | The acceptor size within SCE (1) or outside SCE (0) |

In comparison folder, the intra or inter folder are created for adjusting the positions of NCL detection tools to the positions of exon boundaries.


#### toolName.circEB or toolName.fusionEB (tab-delimited text file)

| No. of column | Description |
|------------|--------------|
| (1) | Chromsome name of the donor side (5'ss) |
| (2) | Exonic junction coordination of the donor side |
| (3)| Strand of the donor side |
| (4) | Chromsome name of the acceptor side (3'ss) |
| (5) | Exonic junction coordination of the acceptor side |
| (6) | Strand of the acceptor side |
| (7) | Total number of supporting NCL-junction reads | 
| (8) | Gene name of the donor side |
| (9) | Gene name of the acceptor side |

