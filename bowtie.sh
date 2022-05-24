#!/bin/bash

#note: to enter phyluce environment: module load phyluce, then source activate phyluce
#to exit: source deactivate
# doesn't like # within command blocks

#SBATCH --job-name=bowtie
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

###################################################
module load bowtie2

source activate bowtie2
###################################################
### Index the reference genome
#Index our reference file into a file that bowtie will understand. The last "bowtie" in the command below is simply a prefix in your name. Make sure you genome is in the directory when you execute the commmand.

bowtie2-build your_reference_genome.fasta bowtie2
###################################################
#map to reference- Here you only need to specify the prefix (-x), not each individual files. -1 are your R1 reads, -2 are your R2 reads -S is the output samfile 

bowtie2 --very-fast-local -x /folder/with/index/files/bowtie2 -1 /clean/read1/files_R1.fastq -2 /clean/read1/files_R2.fastq -S /directory/to/samfile/samfile.sam

#to loop the command across all reads use 
#for sample in `ls /media/sample/fastqfiles/*R1.fastq.gz`
#do
#dir="/media/sample/fastqfiles"
#base=$(basename $sample "_R1.fastq")
#bowtie2 -x path_to_my_index -1 ${dir}/${base}_R1.fastq.gz -2 ${dir}/${base}_R2.fastq.gz -S ${dir}/${base}.sam
#done

INPUT_DIRECTORY=/mnt/home/vascokar/mastitis_study/trimmed

module load SAMtools

source activate samtools

cd $INPUT_DIRECTORY

for f in *_R1_paired.fastq.gz
do
n=${f%%_R1_paired.fastq.gz} 
bowtie2 –-threads 5 -x $INPUT_DIRECTORY/ARS-UCD1.2 -1 $INPUT_DIRECTORY/${n}_R1_paired.fastq.gz -2 $INPUT_DIRECTORY/${n}_R2_paired.fastq.gz -S $INPUT_DIRECTORY/${n}_R1_host_mapped_and_unmapped.sam
done
 

