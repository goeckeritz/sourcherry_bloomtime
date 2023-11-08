#!/bin/sh --login
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH -J lastal
#SBATCH --time 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12g
#SBATCH -o /mnt/scratch/goeckeri/MCScan_%j

module purge
module load icc/2017.4.196-GCC-6.4.0-2.28  impi/2017.3.196
module load Python/3.7.2
module load texlive/20210316
module load LAST/914

cd /mnt/home/goeckeri/software/jcvi

source jcvi/bin/activate

cd /mnt/scratch/goeckeri/REGINA_v_sourcherry/

python -m jcvi.formats.gff bed --type=mRNA --key=ID /mnt/scratch/goeckeri/REGINA_v_sourcherry/chr4A.gff3 -o subA.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID /mnt/scratch/goeckeri/REGINA_v_sourcherry/PAV04_REGINA.gff3 -o regina.bed

python -m jcvi.formats.fasta format /mnt/scratch/goeckeri/REGINA_v_sourcherry/Mont_chr4A_cds.fasta subA.cds
python -m jcvi.formats.fasta format /mnt/scratch/goeckeri/REGINA_v_sourcherry/regina_chr4_cds.fasta regina.cds

python -m jcvi.compara.catalog ortholog subA regina --no_strip_names
