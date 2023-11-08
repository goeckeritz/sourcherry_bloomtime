#!/bin/sh --login
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH -J full_genome_orthofinder
#SBATCH --time 14:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16g
#SBATCH -o /mnt/gs21/scratch/goeckeri/orthofinder_montAv5_scaffolded/arabidopsis_orthofinder_%jar

#troubleshooting the weird error that orhofinder is giving me
module purge
module load GCC/8.3.0 OpenMPI/3.1.4
module load OrthoFinder/2.5.4-Python-2.7.16 #not sure why the python version of orthofinder/2.5.4 wasn't working anymore; somehow the python dependencies got fucked up. 
module load FastTree/2.1.11


cd /mnt/gs21/scratch/goeckeri/orthofinder_montAv5_scaffolded/Genomes/

orthofinder -t 24 -a 5 -M msa -A mafft -z -oa -p /mnt/gs21/scratch/goeckeri/orthofinder_montAv5_scaffolded/Genomes/temp/ -n orthofinder_sourcherry_arabidopsis_full -f proteins/

#well fuck. this version of orthofinder (2.4.1) doesn't recognize the z flag. that prevents automatic trimming.
#however, if I don't need to back-convert the sequences to cds for tree-making, that might be fine for these purposes. 

#don't have to do trees -- looking for orthogroups our candidate genes are in only; I also don't even really need alignments, honestly. So I'll use -og to stop it at inferring orthogroups.

#for subgenome A: awk '/^Pcer_chr[1-8]A\t/' mRNAs | grep -o 'ID=Pcer_.*-R[A,B];Parent' | sed 's/ID=//' | sed 's/;Parent//'
#repeat with others

# awk '/^Pcer_chr[1-8]A__\t/' mRNAs | grep -o 'ID=Pcer_.*-R[A,B];Parent' | sed 's/ID=//' | sed 's/;Parent//'
# awk '/^Pcer_chr[1-8]B\t/' mRNAs | grep -o 'ID=Pcer_.*-R[A,B];Parent' | sed 's/ID=//' | sed 's/;Parent//'




