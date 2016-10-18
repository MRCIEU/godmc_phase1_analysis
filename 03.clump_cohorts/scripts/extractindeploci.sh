#!/bin/bash

#PBS -N indeploci
#PBS -o /panfs/panasas01/shared-godmc/job_report/indeploci-output
#PBS -e /panfs/panasas01/shared-godmc/job_report/indeploci-error
#PBS -l walltime=12:00:00
#PBS -t 1-15
# PBS -t 3
#PBS -l nodes=1:ppn=8
#PBS -S /bin/bash

set -e
if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

echo "Running on ${HOSTNAME}"


filename="/panfs/panasas01/shared-godmc/scripts/cohorts.txt"

user=`head -n ${PBS_ARRAYID} ${filename} | tail -n 1 | cut -d " " -f 1`
cohort=`head -n ${PBS_ARRAYID} ${filename} | tail -n 1 | cut -d " " -f 2`

echo $user
echo $cohort

###extract cis/trans SNP.CpG pairs from *Robj and generate count summaries
 
cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/clump_cohorts/data


Rscript /panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/clump_cohorts/scripts/extractindeploci.R $user $cohort