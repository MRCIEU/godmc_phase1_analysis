#!/bin/bash

# This script takes about 20 minutes to run with 74 chunks and 62 million hits

# 1. Start with 74 association lists. These are split by CPG names.
# 2. Sort each list to be alphabetical wrt column 3 (CPG name)
# 3. Combine into single file
# 4. Split into 500 files
# 5. For each split file, check that the last CPG is not the same as the first CPG in the next file
# 6. If it is, remove those CPGs from the file and add to the start of the next file.
# 7. Get probe list and gzip assoclist


# assoclist="testassoc.txt.gz"
outname="assoclist"
intermediate="temp"


# 1-3

for f in /panfs/panasas01/shared-godmc/counts/combined/*cg*ge1.2*.txt.gz
do
	zcat ${f} | sort -t" " -k3,3 
done > ${intermediate}

# 4

awk 'NR%125000==1 { file = FILENAME "_" sprintf("%d", NR/125000+1) } { print > file }' ${intermediate}
rm ${intermediate}


# 5-6

n=`ls ${intermediate}_[0-9]* | wc -w`
num=$(($n - 1))

for i in $(seq 1 ${num})
do
	ii=$(($i + 1))
	echo "${i} ${ii}"

	lastcpg=`tail -n 1 ${intermediate}_${i} | cut -d " " -f 3`
	firstcpg=`head -n 1 ${intermediate}_${ii} | cut -d " " -f 3`

	if [ "$lastcpg" == "$firstcpg" ]; then
		echo "${lastcpg} ${firstcpg} same"
		grep -w ${lastcpg} ${intermediate}_${i} > temp
		grep -wv ${lastcpg} ${intermediate}_${i} > temp2
		mv temp2 ${intermediate}_${i}
		cat temp ${intermediate}_${ii} > temp2
		mv temp2 ${intermediate}_${ii}
		rm temp
	else
		echo "${lastcpg} ${firstcpg} different"
	fi
done


# 7

for i in $(seq 1 ${n})
do
	echo ${i}
	mv ${intermediate}_${i} ${outname}_${i}
	cut -d " " -f 3 ${outname}_${i} | uniq > ${outname}_${i}.probes
	gzip ${outname}_${i}
done

