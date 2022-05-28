# Genome assembly Bond Lab protocol 
Here we are mapping low-coverage pair-end 10x genomic reads to a reference genome. We will end up with **Sequence Alignment Map (SAM)** files- one per specimen. <br>

## Table of contents
1. [Clean reads](#cleaning)
2. [OPTIONAL Index the genome](#index) 
3. [Map reads to genome](#mapping)  
    1. [OPTIONAL selecting the 'best genome assembly'](#choosegenome)
    2. [OPTIONAL installing the programs](#install)
    3. 


# STEP 1: Clean raw reads using illumiprocessor                                       

### Set up your contig file <a name="mapping"></a>
Your contig file has 4 sections 
1. [adapters] - Universal Adapters - these do not change 
2. [tag sequences] - barcode for each column and row
3. [tag map] - unique barcode pair combination for each well
4. [names] - link your barcode pair combos to your specimen ID

```
[adapters]
i7:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT*GTGTAGATCTCGGTGGTCGCCGTATCATT
```
Tag sequence are your unique barcodes (i7 and i5) for each well
```
[tag sequences]
```
tag map
```

```

### Run illumprocessor 
illumiprocessor can take a while if you are working with low coverage 10X data. After illumiprocessor finishes you will have a folder called **clean-reads** which contains cleaned R1.fasta.gz and R2.fasta.gz for all paired ends. 

```
#!/bin/bash       

[INSERT YOUR SLURM SCRIPT BATCH BLOCK HERE - you can copy from earlier script]

module load phyluce/1.7.1

source activate phyluce

cd /directory/to/your/working/folder

illumiprocessor \
    --input raw-reads/ \
    --output clean-reads \
    --config illumiprocessor.conf \
    --cores 20
```
# STEP 2: OPTIONAL select your genome assembly (optional) <a name="choosegenome"></a>

If your have multiple whole genome assemblies in your folder, you'll need to select the 'best' genome to map your reads to. Look at the metrics text file. This will have statistics for both your contigs and your scaffolds.

contigs = number of contigs<br>
largest contig = longest/most bp<br>
total legnth = entire length of genome assembly<br>
GC (%) = quality score?<br>
N50 = the length of the smallest contig at 50% <br>
N75 = the length of the smallest contig at 75% <br>
L50 =<br>
L75 = <br>
N's per 100 kbp = number of ambiguious bases per 100 kbp - the higher the number the more ambiguities you have <br>
```
Assembly                qqAptStep1.NCBI.p_ctg   qqAptStep1.NCBI.a_ctg
 contigs               3036                    17333
Largest contig          25992589                17784358
Total length            3634786791              3919114782
GC (%)                  40.55                   40.66
N50                     2589198                 781664
N75                     1298795                 356880
L50                     380                     1305
L75                     876                     3125
 scaffolds             1198                    13461
Largest scaffold        58246652                84834143
Total length            3634971263              3919502458
GC (%)                  40.55                   40.66
N50                     12871321                3263039
N75                     4219297                 878507
L50                     67                      227
L75                     183                     812
 N's per 100 kbp       5.07                    9.89       
```

# STEP 3: Map reads to a reference genome <a name="mapping"></a>

### Installing the programs (optional) <a name="install"></a>
_Note: You only need to download these if you are working on locally your computer. You do NOT need to install these if you are runnning the analysis on the UC Davis farm cluster_

Install Bowtie2<br>
_type this command into your terminal and hit enter_
```
conda install -c bioconda bowtie2<br>
```
Install samtools
_type this command into your terminal and hit enter_
```
conda install -c bioconda samtools
``~

help and list of commands
```
bowtie2 --help
bowtie2-build --help
```

### Index the reference genome <a name="index"></a>

Before we map our reads, we need to index our reference file into a set of files that bowtie will understand. After running this step you will end up with a directory with a set of files The last "bowtie" in the command below is simply a prefix in your name. Make sure you genome is in the directory when you execute the commmand.

```
bowtie2-build your_reference_genome.fasta bowtie2
```
### Map to reference

-x path to your index files. You must include the prefix that you indicated when you generated your prefix files. In this case we used the prefix "bowtie2"<br>
-l READ1 reads<br>
-2 READ2 reads<br>
-S output samfile<br>

_single read pair_

```
bowtie2 --very-fast-local -x /folder/with/index/files/bowtie2 -1 /clean/read1/files_R1.fastq -2 /clean/read1/files_R2.fastq -S /directory/to/samfile/samfile.sam
```
_to loop the command across all reads use this command_

```
cd /dirctory/to/cleanreads/

for sample in `ls /media/sample/fastqfiles/*READ1.fastq.gz`
do
dir="/directory/to/cleanreads/$i/split-adapter-quality-trimmed"
base=$(basename $sample "-READ1.fastq.gz")
bowtie2 -x path_to_my_index -1 ${dir}/${base}-READ1.fastq.gz -2 ${dir}/${base}-READ2.fastq.gz -S ${dir}/${base}.sam
done
```

take a look at your output
```
cat /folder/to/samfile/samfile.sam|less
```

most analyses use BAM files not SAM files- SAM files are human readable. Must convert into Binary Alingment Map (BAM)- NO PROGRAMS WILL WORK WITH SAM FILES
```
samtools view -S -b /folder/to/samfile/samfile.sam > /folder/to/samfile/bamfile.bam
```
-b specifies you want bAMfile output

```
cat /folder/to/samfile/bamfile.bam|less
cat /folder/to/samfile/samfile.sam|less
```



