###################################################
#map to reference- Here you only need to specify the prefix (-x), not each individual files. -1 are your R1 reads, -2 are your R2 reads -S is the output samfile 

module load SAMtools

source activate samtools
###################################################

#to loop the command across all reads use 
#for sample in `ls /media/sample/fastqfiles/*R1.fastq.gz`
#do
#dir="/media/sample/fastqfiles"
#base=$(basename $sample "_R1.fastq")
#bowtie2 -x path_to_my_index -1 ${dir}/${base}_R1.fastq.gz -2 ${dir}/${base}_R2.fastq.gz -S ${dir}/${base}.sam
#done

INPUT_DIRECTORY=/mnt/home/vascokar/mastitis_study/trimmed

cd $INPUT_DIRECTORY

for f in *_R1_paired.fastq.gz
do
n=${f%%_R1_paired.fastq.gz} 
bowtie2 –-threads 5 -x $INPUT_DIRECTORY/ARS-UCD1.2 -1 $INPUT_DIRECTORY/${n}_R1_paired.fastq.gz -2 $INPUT_DIRECTORY/${n}_R2_paired.fastq.gz -S $INPUT_DIRECTORY/${n}_R1_host_mapped_and_unmapped.sam
done
 
