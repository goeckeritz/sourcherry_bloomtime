#!/bin/sh --login
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --time=3:59:59
#SBATCH --mem=36G
#SBATCH -J bwa_mem_cherry_index

module load GNU/6.4.0-2.28 OpenMPI/2.1.1 
module load BWA/0.7.17

#To use the simpler alignment method, evidently I had to have more memory. Grrrrr. Looks like the uncompressed fasta file was ~2G, so the memory needed was 5.37(~2). Over 10G, in other words. So I had to submit this as a job -.-
bwa index -p /mnt/home/goeckeri/genome_files/mont/genome_submission/montAv5_p99_400k_final/Pcer_montAv5_p99_400k_scaffolded_assembly -a is /mnt/home/goeckeri/genome_files/mont/genome_submission/montAv5_p99_400k_final/Pcer_montAv5_p99_400k_scaffolded_assembly.fasta