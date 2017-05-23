## Manual of NCLcomparator
###### Version: 1.0

##### NCLcomparator, a comprehensive tool for analyzing non-co-linear (NCL) transcripts (fusion, trans-splicing, and circular RNA), enables to combine several NCL results from different detection tools, and to provide the characteristics of NCL events.
--------------
####

The NCLcomparator program, document, and test set can be downloaded from our FTP site: ftp://treeslab1.genomics.sinica.edu.tw/NCLcomparator or GitHub: https://github.com/TreesLab/NCLcomparator.

### 1. System requirements
  NCLcomparator runs under Linux-like environment (i.e. Bio-Linux, also see http://environmentalomics.org/bio-linux/) with at least 
       
     
### 2. Installation
     $ tar zxvf NCLcomparator.tar.gz
     $ cd NCLcomparator
     $ chmod +x NCLcomparator.sh
     $ chmod +x /bin/*

### 3. External tools
   (1) bedtools (http://bedtools.readthedocs.io/en/latest/) <br>
   (2) STAR (https://github.com/alexdobin/STAR) <br>
   (3) RSEM (https://github.com/deweylab/RSEM) <br>
   (4) R (https://www.r-project.org/)<br>

##### Get latest bedtools source from releases and install it 
     $ wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools2.25.0.tar.gz
     $ tar -zxvf bedtools-2.25.0.tar.gz
     $ cd bedtools2
     $ make
     $ cp ./bin/ /usr/local/bin
##### Get latest STAR source from releases and install it 
     $ wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
     $ tar -xzf 2.5.3a.tar.gz
     $ cd STAR-2.5.3a
     $ cp bin/Linux_x86_64_static/STAR /usr/local/bin
##### Get latest RSEM source code and install it #
     $ wget https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz
     $ tar -xzf v1.3.0.tar.gz
     $ cd RSEM-1.3.0
     $ make
     $ make ebseq
##### Get R in Ubuntu environment
     $ sudo apt-get update
     $ sudo apt-get install r-base

### 4. Preparation
 (1)	Genome and its annotation, which can be download from the GENCODE website (http://www.gencodegenes.org/) or ensembl FTP   (http://www.ensembl.org/info/data/ftp/index.html). Given Human as an example, go to ensembl FTP (http://www.ensembl.org/info/data/ftp/index.html) to download human genome and annotation.
      
     $ wget ftp://ftp.ensembl.org/pub/release88/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
     $ wget ftp://ftp.ensembl.org/pub/release88/gtf/homo_sapiens/Homo_sapiens.GRCh38.88.gtf.gz
     $ gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
     $ gunzip Homo_sapiens.GRCh38.88.gtf.gz
 
 (2)	(optional) Synonymous Constraint elements (SCE), which can be download from (http://compbio.mit.edu/SCE/)

 (3) Building STAR and RSEM index  
 
The following steps are to generate the index files before running NCLcomparator tool. 

    $ mkdir STAR_RSEM_index
    $ cd STAR_RSEM_index
    $ STAR --runMode genomeGenerate --runThreadN 10 --genomeDir . --genomeFastaFiles /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa
    $ rsem-prepare-reference --gtf /path/to/Homo_sapiens.GRCh38.85.gtf --star -p 10 /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa RSEM

The STAR and RSEM index of the genome hg38 and the annotation ensemble 85 is prepared and can be downloaded in our FTP. 

  (4) Installation of R packages in R console
  
      $ R
      $ install.packages("ggplot2")
      $ install.packages("grodExtra")
      $ Would you like to use a personal library instead?  (y/n) y
      $ Would you like to create a personal library~/R/x86_64-pc-linux-gnu-library/3.2 to install packages into?  (y/n) y
      $ q()
      $ Save workspace image? [y/n/c]: n
          
### 5. Execution of NCLcomparator
    
 ##### Usage:
 
     $ ./NCLcomparator.sh -gtf [annotation GTF file] -thread [number of thread] -read1[read1 fastq.gz] -read2 [read2 fastq.gz] -index [STAR_RSEM index folder] -intra [circular result folder] -inter [fusion result folder] -sce [SCE bed file]

 ##### An example: 
 
     $ ./NCLcomparator.sh -gtf Homo_sapiens.GRCh38.85.gtf -intra /path/to/intra -inter /path/to/inter -sce SCE_hg38.bed -read1 GM12878_1.fastq.gz -read2 GM12878_2.fastq.gz -index /path/to/STAR_RSEM_index

##### The basic options to run a job as follow:
##### -intra/--circular /path/to/NCL-intra folder
##### -inter/--fusion (optional) /path/to/NCL-inter folder
##### -sce/--SCE (optional) /path/to/SCE.bed

### 6. Input format of NCL events
The outputs of NCL detection tools are required to modify as 5-col format, which includes the positons of donor/acceptor sides and the number of junction reads. Intra-NCL tools' modified 5-col format results are gathered into a folder, and inter-NCL tools' modified 5-col format results are gathered into another folder.  
  
Given a 5-col result, CIRCexplorer2.5col (tab-delimited text file), as an example,  

| Chromosome | Coordination | Chromosome | Coordination | Total number of junction reads | 
|------------|--------------|------------|--------------| -------------------------------|
| chr1	| 945176 | chr1	| 945517 |	1 |
| chr1	| 955922 | chr1	| 956013 |	1 |
| chr1	| 955922 | chr1	| 957273 |  5 |
     
### 7. Output files
After executing NCLcomparator program, two merged circular RNA and fusion RNA tools' results (intraMerged.result and interMerged.result) accompanied by its graphic report (intra.pdf and inter.pdf) are generated, and two newly folder is created (comparison and STAR_RSEM_out).

In STAR_RSEM_out folder, several output files are generated by running STAR and RSEM programs, these files are given as the input files of characteristic process in NCLcomaprator. In comparison folder, the two folders (intra and inter) are generated, where NCL events of the positions in each tool with <= 5 bp of franking exon boundary are adjusted to the exact positions of exon boundary, and the other two folders intra_Characteristic and inter_Characteristic) are produced, where many genomic feathers of NCL events in each tool are given.  

The column formats of outputs are described as follow:
##### intraMerged.result (tab-delimited text file)

|No. of columun | Description |
|---------------|-------------|
| (1) | Chromosome name of the donor side (5' ss) |
| (2) | Junction coordination of the donor side |
| (3) | Strand of the donor side |
| (4) | Chromosome name of the acceptor side (3' ss) |
| (5) | Junction coordinate of the acceptor side |
| (6) | Strand of the acceptor side |
| (7) | Gene name |
| (8) | Circular tool1 (Yes: 1; No: 0) |
| (9) | Circular tool2	(Yes: 1; No: 0) |
|(10) | ..|

##### interMerged.result (tab-delimited text file)

|No. of columun | Description |
|---------------|-------------|
| (1) | Chromosome name of the donor side (5' ss) |
| (2) | Junction coordination of the donor side |
| (3) | Strand of the donor side |
| (4) | Chromosome name of the acceptor side (3' ss) |
| (5) | Junction coordinate of the acceptor side |
| (6) | Strand of the acceptor side |
| (7) | Gene name of the donor side |
| (8) | Gene nqamae of the acceptor side |
| (9) | Fusion tool1 (Yes: 1; No: 0) |
| (10) | Fusion tool2	(Yes: 1; No: 0) |
|(11) | ..|


In comparison folder, the intra or inter folder are created for adjusting the positions of NCL detection tools to the positions of exon boundary.
##### toolName.circEB or toolName.fusionEB (tab-delimited text file)

| No. of column | Description |
|---------------|-------------|
| (1) | Chromosome name of the donor side (5' ss) |
| (2) | Junction coordination of the donor side |
| (3) | Strand of the donor side |
| (4) | Chromosome name of the acceptor side (3' ss) |
| (5) | Junction coordinate of the acceptor side |
| (6) | Strand of the acceptor side |
| (7) | Total number of junction reads |
| (8) | Gene name of the donor side |
| (9) | Gene name of the acceptor side |	


The other two folders, intra_Characteristic and inter_Characteristic, provide the genomic feathers of NCL events with respect to each tool's result.
##### intra-NCLToolName.final (tab-delimited text file)

| No. of column | Description |
|---------------|-------------|
| (1) | Chromosome name of the donor side (5' ss) |
| (2) | Junction coordination of the donor side |
| (3) | Strand of the donor side |
| (4) | Chromosome name of the acceptor side (3' ss) |
| (5) | Junction coordinate of the acceptor side |
| (6) | Strand of the acceptor side |
| (7) | Total number of junction read |
| (8) | Gene name |
| (9) | Total number of linear junction reads at the donor side |
| (10) | RPM based on total RNA-seq reads |
| (11) | RPM based on uniquely mapping reads |
| (12) | FPKM of the gene |
| (13) | TPM of the gene |
| (14) | Frequency of occurrence of linear junction reads at the donor side of the gene |
| (15) | Frequency of occurrence of linear junction reads at the acceptor side of the gene |
| (16) | Median of frequency of occurrence of linear junction reads among all junctions in the gene |
| (17) | Out of circle (Yes: 1; No: 0) |
| (18) | (Optional) the donor size within SCE (1) or outside SCE (0) |	
| (19) | (Optional) the acceptor size within SCE (1) or outside SCE (0) |

##### intra-NCLToolName.final (tab-delimited text file)
     
| No. of column | Description |      
|---------------|-------------|
| (1) | Chromosome name of the donor side (5' ss) |
| (2) | Junction coordination of the donor side |
| (3) | Strand of the donor side |
| (4) | Chromosome name of the acceptor side (3' ss) |
| (5) | Junction coordinate of the acceptor side |
| (6) | Strand of the acceptor side |
| (7) | Total number of junction read |
| (8) | Gene name of the donor side |
| (9) | Gene name of the acceptor side |	
| (10)| Number of linear junction reads at the donor side |	
| (11)| Number of linear junction reads at the acceptor side |
| (12)| RPM based on total reads |
| (13)| RPM based on uniquely mapping reads |
| (14)| FPKM of the donor gene |
| (15)| TPM of the donor gene |
| (16)| FPKM of the acceptor gene |
| (17)| TPM of the acceptor gene |
| (18)| Frequency of occurrence of linear junction reads at the donor side of the donor gene |
| (19)| Frequency of occurrence of linear junction reads at the acceptor side of the acceptor gene |
| (20)| Median of frequency of occurrence of linear junction reads among all junctions in the donor gene |
| (21)| Median of frequency if occurrence of linear junction reads among all junctions on the acceptor gene |
| (22)| (Optional) the donor side within SCE (1) or outside SCE (0) |
| (23)| (Optional) the acceptor size within SCE (1) or outside SCE (0) |
