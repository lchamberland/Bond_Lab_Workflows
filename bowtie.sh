#!/bin/bash

#note: to enter phyluce environment: module load phyluce, then source activate phyluce
#to exit: source deactivate
# doesn't like # within command blocks

#SBATCH --job-name=bowtie
#SBATCH --nodes=1
#SBATCH --ntasks=2                               
#SBATCH --mem=24G                               
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

###################################################
module load bowtie2

cd /directory/to/your/genome

### Index the reference genome
#Index our reference file into a file that bowtie will understand. The last "bowtie" in the command below is simply a prefix in your name. Make sure you genome is in the directory when you execute the commmand.

bowtie2-build your_reference_genome.fasta bowtie2

### running bowtie2
bowtie2 --very-fast-local -x /folder/with/index/files/bowtie2 -1 /clean/read1/files_R1.fastq -2 /clean/read1/files_R2.fastq -S /directory/to/samfile/samfile.sam



