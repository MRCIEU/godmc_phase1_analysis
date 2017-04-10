#!/bin/bash

#PBS -N counts
#PBS -o /panfs/panasas01/shared-godmc/job_report_2017/counts_file_number-output
#PBS -e /panfs/panasas01/shared-godmc/job_report_2017/counts_file_number-error
#PBS -l walltime=00:30:00
#PBS -t 1-23
# PBS -t 19
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash

set -e
if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

echo "Running on ${HOSTNAME}"

# Check that all counts have been performed for all CpG subsets and P Value thresholds

# Number of files expected:

## 100 meQTL chunks
### .gz of 100 chunks * 74 subsets for both cis and trans = (9*100)*2 = 1800
### combined counts for 74 subsets * pvalue thresholds for both cis and trans = (74*9)*2 = 1332
### combined counts for pvalue thresholds for both cis and trans = 9*2 = 18
#### 1800+18+1332 = 3150

## 500 meQTL chunks
### .gz of 100 chunks * 74 subsets for both cis and trans = (9*100)*2 = 9000
### combined counts for 74 subsets * pvalue thresholds for both cis and trans = (74*9)*2 = 1332
### combined counts for pvalue thresholds for both cis and trans = 9*2 = 18
#### 9000+18+1332 = 10,350


## 1000 meQTL chunks
### .gz of 100 chunks * 74 subsets for both cis and trans = (9*1000)*2 = 18000
### combined counts for 74 subsets * pvalue thresholds for both cis and trans = (74*9)*2 = 1332
### combined counts for pvalue thresholds for both cis and trans = 9*2 = 18
#### 18000+18+1332 = 19,350

echo "Running file check"

MYDIR="/panfs/panasas01/sscm/epwkb/GoDMC_Analysis/CpG_SNP_FINAL"
filename="$MYDIR/repo/godmc_phase1_analysis/01.extract_sftp/data/cohorts.txt"

user=`head -n ${PBS_ARRAYID} ${filename} | tail -n 1 | cut -d " " -f 1`
cohort=`head -n ${PBS_ARRAYID} ${filename} | tail -n 1 | cut -d " " -f 2`

echo $user
echo $cohort

cd /panfs/panasas01/shared-godmc/counts_2017/${user}_${cohort}

echo "if 100 meQTL chunks, then 3150 files"
echo "if 500 meQTL chunks, then 10,350 files"
echo "if 1000 meQTL chucks, then 19,350 files"

ls -1 | wc -l
