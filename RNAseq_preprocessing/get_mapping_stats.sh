#!/bin/sh

for i in $(ls -1 *relaxedLog.final.out); do
	echo $i >> mapping_stats.txt
	less $i >> mapping_stats.txt
	echo "############################################################################" >> mapping_stats.txt
done

	