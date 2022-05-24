# Genome assembly Bond Lab protocol 



## Clean raw reads using illumiprocessor                                       

### Set up your contig file 
Your contig file has 3 sections 
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


```
#!/bin/bash

#note: to enter phyluce environment: module load phyluce, then source activate $
#to exit: source deactivate
# doesn't like # within command blocks

#SBATCH --job-name=P002-name
#SBATCH --nodes=1
#SBATCH --ntasks=2                               
#SBATCH --mem=64G                               
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --cpus-per-task=16
#SBATCH --time=2-20:00                          
#SBATCH --partition=high2                           
#SBATCH --reservation=                           
#SBATCH --output=phyluce-%N-%j.out               
#SBATCH --error=phyluce-%N-%j.err                
#SBATCH --mail-type=ALL                          
#SBATCH --mail-user=lchamberland@ucdavis.edu        


#might need to include this if running not with module, not sure
#OMPI_MCA_opal_cuda_support=true


# To activate this environment, use
#
#  $ conda activate phyluce-1.7.0
#
# To deactivate an active environment, use
#
#     $ conda deactivate


#module load phyluce/1.6.8
#module load phyluce
module load phyluce/1.7.1

source activate phyluce
#conda deactivate

cd yourdirectory

illumiprocessor \
    --input raw-reads/ \
    --output clean-reads \
    --config illumiprocessor.conf \
    --cores 20
```

Illumina/454/IonTorrent paired-end reads longer than ~70bp:

```
bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
```

If your have multiple whole genome assemblies in your folder, look at the metrics text file. This will have statistics for both your contigs and your scaffolds.

contigs = number of contigs
largest contig = longest/most bp
total legnth = entire length of genome assembly
GC (%) = quality score?
N50 = the length of the smallest contig at 50% 
N75 = the length of the smallest contig at 75% 
L50 =
L75 = 
N's per 100 kbp = number of ambiguious bases per 100 kbp - the higher the number the more ambiguities you have 


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

# Installing the programs
_Note: You only need to download these if you are working on locally your computer. You do NOT need to install these if you are runnning the analysis on the UC Davis farm cluster_

### Install Bowtie2 
type this command into your terminal and hit enter
```
conda install -c bioconda bowtie2
```
### Install samtools
```
conda install -c bioconda samtools
```

# Bowtie 2 - map reads to a reference genome
help and list of commands
```
bowtie2 --help
bowtie2-build --help
```
### Index the reference genome
Index our reference file into a file that bowtie will understand. The last "bowtie" in the command below is simply a prefix in your name. Make sure you genome is in the directory when you execute the commmand.

```
bowtie2-build your_reference_genome.fasta bowtie2
```
### map to reference 
Here you only need to specify the prefix (-x), not each individual files. -1 are your R1 reads, -2 are your R2 reads -S is the output samfile 
```
bowtie2 --very-fast-local -x /folder/with/index/files/bowtie2 -1 /clean/read1/files_R1.fastq -2 /clean/read1/files_R2.fastq -S /directory/to/samfile/samfile.sam
```
to loop the command across all reads use 
```
for sample in `ls /media/sample/fastqfiles/*R1.fastq.gz`
do
dir="/media/sample/fastqfiles"
base=$(basename $sample "_R1.fastq")
bowtie2 -x path_to_my_index -1 ${dir}/${base}_R1.fastq.gz -2 ${dir}/${base}_R2.fastq.gz -S ${dir}/${base}.sam
done
```

```


 module load GCC/9.3.0 Bowtie2/2.4.1 SAMtools/1.11
cd /directory/with/samfile/

  for f in *_R1_paired.fastq.gz
  do
 n=${f%%_R1_paired.fastq.gz} 

 bowtie2 –-threads 5 -x $INPUT_DIRECTORY/ARS-UCD1.2 -1 $INPUT_DIRECTORY/${n}_R1_paired.fastq.gz -2 $INPUT_DIRECTORY/${n}_R2_paired.fastq.gz -S $INPUT_DIRECTORY/${n}_R1_host_mapped_and_unmapped.sam
 done
 ```


```
echo "bowtie2 -x path_to_my_index -1 ${dir}/${base}_R1.fastq -2 ${dir}/${base}_R2.fastq -S ${dir}/${base}.sam"
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



