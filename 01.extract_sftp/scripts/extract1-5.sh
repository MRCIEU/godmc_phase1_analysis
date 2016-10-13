#!/bin/bash

#PBS -N extract1-4
#PBS -o /panfs/panasas01/shared-godmc/job_report/extract1-5-output
#PBS -e /panfs/panasas01/shared-godmc/job_report/extract1-5-error
#PBS -l walltime=12:00:00
#PBS -t 3
#PBS -l nodes=1:ppn=8
#PBS -S /bin/bash

set -e
if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

echo "Running on ${HOSTNAME}"


filename="/panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/01.extract_sftp/data/cohorts.txt"

user=`head -n ${PBS_ARRAYID} ${filename} | tail -n 1 | cut -d " " -f 1`
cohort=`head -n ${PBS_ARRAYID} ${filename} | tail -n 1 | cut -d " " -f 2`

echo $user
echo $cohort

mkdir  -p /panfs/panasas01/shared-godmc/sftp/GoDMC/$user
cd /panfs/panasas01/shared-godmc/sftp/GoDMC/$user
mkdir -p $cohort

if [ -f $cohort\_01.tgz ]; then

   mv $cohort\_01.tgz ./$cohort
   cd $cohort
   tar -zxvf $cohort\_01.tgz
   cd ..
else
  echo "$cohort\_01.tgz does not exist"
fi

if [ -f $cohort\_02.tgz ]; then

   mv $cohort\_02.tgz ./$cohort
   cd $cohort
   tar -zxvf $cohort\_02.tgz
   cd ..
else
  echo "$cohort\_02.tgz does not exist"
fi

if [ -f $cohort\_03.tgz ]; then

   mv $cohort\_03.tgz ./$cohort
   cd $cohort
   tar -zxvf $cohort\_03.tgz
   cd .. 
else
  echo "$cohort\_03.tgz does not exist"
fi

if [ -f $cohort\_04.tgz ]; then

   mv $cohort\_04.tgz ./$cohort
   cd $cohort
   tar -zxvf $cohort\_04.tgz
   cd ..
else
  echo "$cohort\_04.tgz does not exist"
fi

if [ -f $cohort\_05.tgz ]; then

   mv $cohort\_05.tgz ./$cohort
   cd $cohort
   tar -zxvf $cohort\_05.tgz
   cd ..
else
  echo "$cohort\_05.tgz does not exist"
fi


#mv $cohort\_01.tgz ./$cohort
#mv $cohort\_02.tgz ./$cohort
#mv $cohort\_03.tgz ./$cohort
#mv $cohort\_04.tgz ./$cohort
#cd $cohort
#tar -zxvf $cohort\_01.tgz
#tar -zxvf $cohort\_02.tgz
#tar -zxvf $cohort\_03.tgz
#tar -zxvf $cohort\_04.tgz

