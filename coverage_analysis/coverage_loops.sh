#!/bin/sh --login
#SBATCH -J cov_loop
#SBATCH --mail-user=goeckeri@msu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=48g
#SBATCH --time=10:00:00
#SBATCH -o /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/coverage_loop_%j

cd /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/coverage/

#chr_list.txt is a simple text file with chromosomes 1[A,A'B] - 8 [A,A',B] listed, one chromosome per line  
#be sure there is a tab at the end of each line in chr_list.txt --> sed -i 's/$/\t/' chr_list.txt
#it just provided an anchor for me to get fgrep to fucking work.

#this loop goes through the depth files, creating a directory for every sample. 
#then, an internal loop subsets every chromosome's read depth for each sample, then goes through each chromosome depth file to average the depth at every 500 kb window.
#first I separate the depth files into their 3 component columns, average column 2 (the base position of the chromosome)
#then I do the same for column 3, which is the read depth at each position. c2 and c3 should have the same number of lines after this
#then the first column is created, which is simply the name of the chromosome repeated however many lines that c2 and c3 are now. (important for plotting). The 3 columns then get pasted together. 
#make sure you're in the directory where all your depth files are. These depth files should already be sorted and have had their duplicates removed. 
for i in $(ls -1 *depth.post.dups | sed 's/.depth.post.dups//'); do
mkdir $(echo $i)
    for j in $(less chr_list.txt); do
    fgrep "$(echo $j)	" $(echo $i).depth.post.dups > ./$(echo $i)/$(echo $j)_coverage
    awk '{print $2}' ./$(echo $i)/$(echo $j)_coverage | awk '{sum+=$1} (NR%500000)==0{print sum/500000; sum=0;}' > ./$(echo $i)/$(echo $j)_500000_c2
    awk '{print $3}' ./$(echo $i)/$(echo $j)_coverage | awk '{sum+=$1} (NR%500000)==0{print sum/500000; sum=0;}' > ./$(echo $i)/$(echo $j)_500000_c3
    yes $(echo $j) | head -`wc -l ./$(echo $i)/$(echo $j)_500000_c3 | awk '{print $1}'`> ./$(echo $i)/$(echo $j)_500000_c1
    paste ./$(echo $i)/$(echo $j)_500000_c1 ./$(echo $i)/$(echo $j)_500000_c2 ./$(echo $i)/$(echo $j)_500000_c3 > ./$(echo $i)/$(echo $j)_500000
done
    	done
    	
#this loop then goes into each sample's directory to combine all the chromosome coverage files into one big fat file for the sample.    	
for i in $(ls -1 *depth.post.dups | sed 's/.depth.post.dups//'); do
        cd /mnt/scratch/goeckeri/pop4_variant_calling/bwa_mem/MarkDups/coverage/$(echo $i)
        cat Pcer_chr*_500000 > $(echo $i)_500000
        cp $(echo $i)_500000 /mnt/home/goeckeri/pop4_coverage/
done

#lastly, we simply add headers to each sample's big fat combination file.
cd /mnt/home/goeckeri/pop4_coverage/

for i in $(ls -1 /mnt/home/goeckeri/pop4_coverage/*500000 | sed 's/\/mnt\/home\/goeckeri\/pop4_coverage\///'); do
$(echo $i)
echo -e "chr\tpos\tcov" | cat - $(echo $i) > $(echo $i)_headers.tsv
	done


#These coverage files are then read into R for coverage plotting. 
#Make sure the headers get added correctly.

