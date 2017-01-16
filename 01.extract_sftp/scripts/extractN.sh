#!/bin/bash

MYDIR="/panfs/panasas01/sscm/epzjlm"
filename="$MYDIR/repo/godmc_phase1_analysis/01.extract_sftp/data/cohorts.txt"

nocohorts=`cat $MYDIR/repo/godmc_phase1_analysis/01.extract_sftp/data/cohorts.txt |wc -l`

echo "user" "cohort" "N" > $MYDIR/repo/godmc_phase1_analysis/01.extract_sftp/data/cohort_samplesizes.txt

for i in `seq 1 $nocohorts`;do

user=`head -n ${i} ${filename} | tail -n 1 | cut -d " " -f 1`
cohort=`head -n ${i} ${filename} | tail -n 1 | cut -d " " -f 2`

if [ $cohort = "Phase1_SCZ" ]; then

N=`grep "people pass filters and QC." /panfs/panasas01/shared-godmc/sftp/GoDMC/${user}/${cohort}/results/04/positive_control_cg07959070.log | cut -d " " -f 4`

elif [ $cohort = "TwinsUK" ]; then

N="833"	

else

N=`grep "people pass filters and QC." /panfs/panasas01/shared-godmc/sftp/GoDMC/${user}/${cohort}/results/04/positive_control_pcadjusted_cg07959070.log | cut -d " " -f 4`

fi

echo $i
echo $user
echo $cohort
echo $N

echo $user $cohort $N >> $MYDIR/repo/godmc_phase1_analysis/01.extract_sftp/data/cohort_samplesizes.txt
done

