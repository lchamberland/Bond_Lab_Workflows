# CCGP_aptostichus
10x genomics assembly


```
ragtag.py scaffold -f 1200 -t 10 ref/YBT-1520.fasta output_dir/ZZQ-130.fasta.
```
DOI : https://doi.org/10.1128/mra.00887-21


cats :https://doi.org/10.1093/jhered/esaa057

Universal Adapters - these do not change 

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
