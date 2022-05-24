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

module load bowtie2

source activate bowtie2
