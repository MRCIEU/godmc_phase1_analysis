#!/bin/bash

#PBS -N counts
#PBS -o /panfs/panasas01/shared-godmc/job_report/counts14-output
#PBS -e /panfs/panasas01/shared-godmc/job_report/counts14-error
#PBS -l walltime=12:00:00
#PBS -t 1
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

#cd /panfs/panasas01/shared-godmc/sftp/GoDMC/$user
#mkdir -p $cohort
#mv $cohort\_05.tgz ./$cohort
#cd $cohort
#tar -zxvf $cohort\_05.tgz

###extract cis/trans SNP.CpG pairs from *Robj and generate count summaries
 
cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/extract_counts/scripts
mkdir -p /panfs/panasas01/shared-godmc/counts/14/${user}_${cohort}

Rscript mqtl.count.R $user $cohort

#74 subgroups

pvals=("1e-05" "1e-06" "1e-07" "1e-08" "1e-09" "1e-10" "1e-11" "1e-12" "1e-13")
cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01[0-4]" "cg01[5-9]" "cg02[0-4]" "cg02[5-9]" "cg03[0-4]" "cg03[5-9]" "cg04[0-4]" "cg04[5-9]" "cg05[0-4]" "cg05[5-9]" "cg06[0-4]" "cg06[5-9]" "cg07[0-4]" "cg07[5-9]" "cg08[0-4]" "cg08[5-9]" "cg09[0-4]" "cg09[5-9]" "cg10[0-4]" "cg10[5-9]" "cg11[0-4]" "cg11[5-9]" "cg12[0-4]" "cg12[5-9]" "cg13[0-4]" "cg13[5-9]" "cg14[0-4]" "cg14[5-9]" "cg15[0-4]" "cg15[5-9]" "cg16[0-4]" "cg16[5-9]" "cg17[0-4]" "cg17[5-9]" "cg18[0-4]" "cg18[5-9]" "cg19[0-4]" "cg19[5-9]" "cg20[0-4]" "cg20[5-9]" "cg21[0-4]" "cg21[5-9]" "cg22[0-4]" "cg22[5-9]" "cg23[0-4]" "cg23[5-9]" "cg24[0-4]" "cg24[5-9]" "cg25[0-4]" "cg25[5-9]" "cg26[0-4]" "cg26[5-9]" "cg27[0-4]" "cg27[5-9]" "_ch")


#cpgs=("ch1" "ch2" "ch3" "ch4" "ch5" "ch6" "ch7" "ch8" "ch9" "chX")
#cpgs=("_ch")
#cpgs=("cg007" "cg008" "cg009")

cd /panfs/panasas01/shared-godmc/counts/14/${user}_${cohort}

for i in ${pvals[@]}; do
echo $i
#rm cis.assoc.*.${user}_${cohort}.$i.txt.gz
#rm trans.assoc.*.${user}_${cohort}.$i.txt.gz
cat cis.assoc.*.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.txt
cat trans.assoc.*.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.txt
gzip cis.assoc.*.${user}_${cohort}.$i.txt
gzip trans.assoc.*.${user}_${cohort}.$i.txt

#the line below doesn't work on the node
#mv *gz /projects/MRC-IEU/groups/godmc/counts/

cd /panfs/panasas01/shared-godmc/counts/14/${user}_${cohort}

for j in ${cpgs[@]}; do
echo $j
grep $j cis.assoc.${user}\_${cohort}.$i.txt > cis.assoc.${user}\_${cohort}.${i}.${j}.txt
grep $j trans.assoc.${user}\_${cohort}.$i.txt > trans.assoc.${user}\_${cohort}.${i}.${j}.txt
done
done


