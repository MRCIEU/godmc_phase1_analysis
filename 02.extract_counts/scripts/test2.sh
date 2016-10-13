#!/bin/bash

#PBS -N counts
#PBS -o /panfs/panasas01/shared-godmc/job_report/counts-output
#PBS -e /panfs/panasas01/shared-godmc/job_report/counts-error
#PBS -l walltime=12:00:00
#PBS -t 1-17
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

mkdir -p /panfs/panasas01/shared-godmc/meta-analysis/inputfiles/${user}_${cohort}

#mv /panfs/panasas01/shared-godmc/meta-analysis/inputfiles/${user}_${cohort}.*.gwama.formatted.txt /panfs/panasas01/shared-godmc/meta-analysis/inputfiles/${user}_${cohort}
mv /panfs/panasas01/shared-godmc/meta-analysis/inputfiles/${user}_${cohort}/${user}_${cohort}.*.gwama.formatted.txt.gz /panfs/panasas01/shared-godmc/meta-analysis/inputfiles
#mv /panfs/panasas01/shared-godmc/meta-analysis/inputfiles/${user}\_${cohort}.probes /panfs/panasas01/shared-godmc/meta-analysis/inputfiles/${user}_${cohort}
