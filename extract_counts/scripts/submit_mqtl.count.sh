#!/bin/bash

#PBS -N counts
#PBS -o /panfs/panasas01/shared-godmc/job_report/counts-output
#PBS -e /panfs/panasas01/shared-godmc/job_report/counts-error
#PBS -l walltime=12:00:00
#PBS -t 2,13,14
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
 
cd /panfs/panasas01/shared-godmc/scripts/

Rscript mqtl.count.R $user $cohort

pvals=("1e-05" "1e-06" "1e-07" "1e-08" "1e-09" "1e-10" "1e-11" "1e-12" "1e-13")
cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01" "cg02" "cg03" "cg04" "cg05" "cg06" "cg07" "cg08" "cg09" "cg10" "cg11" "cg12" "cg13" "cg14" "cg15" "cg16" "cg17" "cg18" "cg19" "cg20" "cg21" "cg22" "cg23" "cg24" "cg25" "cg26" "cg27" "_ch")
#cpgs=("ch1" "ch2" "ch3" "ch4" "ch5" "ch6" "ch7" "ch8" "ch9" "chX")
#cpgs=("_ch")
#cpgs=("cg007" "cg008" "cg009")

cd /panfs/panasas01/shared-godmc/counts

for i in ${pvals[@]}; do
echo $i
cat cis.assoc.*.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.txt
cat trans.assoc.*.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.txt
gzip cis.assoc.*.${user}_${cohort}.$i.txt
gzip trans.assoc.*.${user}_${cohort}.$i.txt

#the line below doesn't work on the node
#mv *gz /projects/MRC-IEU/groups/godmc/counts/

cd /panfs/panasas01/shared-godmc/counts

for j in ${cpgs[@]}; do
echo $j

grep $j cis.assoc.${user}\_${cohort}.$i.txt > cis.assoc.${user}\_${cohort}.${i}.${j}.txt
grep $j trans.assoc.${user}\_${cohort}.$i.txt > trans.assoc.${user}\_${cohort}.${i}.${j}.txt
done
done

#grep "cg0000[0-9]" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg0000_0_9.txt
#grep "cg0001" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg0001.txt
#grep "cg0002" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg0002.txt
#grep "cg0003" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg0003.txt
#grep "cg0004" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg0004.txt
#grep "cg0005" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg0005.txt
#grep "cg0006" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg0006.txt
#grep "cg0007" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg0007.txt
#grep "cg0008" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg0008.txt
#grep "cg0009" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg0009.txt
#grep "cg001" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg001.txt
#grep "cg002" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg002.txt
#grep "cg003" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg003.txt
#grep "cg004" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg004.txt
#grep "cg005" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg005.txt
#grep "cg006" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg006.txt
#grep "cg007" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg007.txt
#grep "cg008" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg008.txt
#grep "cg009" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg009.txt
#grep "cg01" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg01.txt
#grep "cg02" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg02.txt
#grep "cg03" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg03.txt
#grep "cg04" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg04.txt
#grep "cg05" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg05.txt
#grep "cg06" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg06.txt
#grep "cg07" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg07.txt
#grep "cg08" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg08.txt
#grep "cg09" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg09.txt
#grep "cg10" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg10.txt
#grep "cg11" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg11.txt
#grep "cg12" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg12.txt
#grep "cg13" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg13.txt
#grep "cg14" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg14.txt
#grep "cg15" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg15.txt
#grep "cg16" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg16.txt
#grep "cg17" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg17.txt
#grep "cg18" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg18.txt
#grep "cg19" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg19.txt
#grep "cg20" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg20.txt
#grep "cg21" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg21.txt
#grep "cg22" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg22.txt
#grep "cg23" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg23.txt
#grep "cg24" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg24.txt
#grep "cg25" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg25.txt
#grep "cg26" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg26.txt
#grep "cg27" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg27.txt
#grep "ch" cis.assoc.${user}_${cohort}.$i.txt > cis.assoc.${user}_${cohort}.$i.cg28.txt

#grep "cg0000[0-9]" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg0000_0_9.txt
#grep "cg0001" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg0001.txt
#grep "cg0002" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg0002.txt
#grep "cg0003" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg0003.txt
#grep "cg0004" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg0004.txt
#grep "cg0005" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg0005.txt
#grep "cg0006" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg0006.txt
#grep "cg0007" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg0007.txt
#grep "cg0008" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg0008.txt
#grep "cg0009" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg0009.txt
#grep "cg001" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg001.txt
#grep "cg002" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg002.txt
#grep "cg003" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg003.txt
#grep "cg004" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg004.txt
#grep "cg005" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg005.txt
#grep "cg006" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg006.txt
#grep "cg007" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg007.txt
#grep "cg008" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg008.txt
#grep "cg009" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg009.txt
#grep "cg01" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg01.txt
#grep "cg02" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg02.txt
#grep "cg03" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg03.txt
#grep "cg04" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg04.txt
#grep "cg05" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg05.txt
#grep "cg06" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg06.txt
#grep "cg07" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg07.txt
#grep "cg08" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg08.txt
#grep "cg09" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg09.txt
#grep "cg10" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg10.txt
#grep "cg11" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg11.txt
#grep "cg12" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg12.txt
#grep "cg13" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg13.txt
#grep "cg14" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg14.txt
#grep "cg15" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg15.txt
#grep "cg16" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg16.txt
#grep "cg17" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg17.txt
#grep "cg18" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg18.txt
#grep "cg19" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg19.txt
#grep "cg20" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg20.txt
#grep "cg21" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg21.txt
#grep "cg22" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg22.txt
#grep "cg23" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg23.txt
#grep "cg24" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg24.txt
#grep "cg25" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg25.txt
#grep "cg26" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg26.txt
#grep "cg27" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg27.txt
#grep "ch" trans.assoc.${user}_${cohort}.$i.txt > trans.assoc.${user}_${cohort}.$i.cg28.txt

#done
#done

