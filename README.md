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
for sample in `ls /media/sample/fastqfiles/*R1.fastq`
do
dir="/media/sample/fastqfiles"
base=$(basename $sample "_R1.fastq")
bowtie2 -x path_to_my_index -1 ${dir}/${base}_R1.fastq -2 ${dir}/${base}_R2.fastq -S ${dir}/${base}.sam
done
```

```
echo "bowtie2 -x path_to_my_index -1 ${dir}/${base}_R1.fastq -2 ${dir}/${base}_R2.fastq -S ${dir}/${base}.sam"
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

