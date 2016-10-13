#!/bin/bash

#PBS -N split
#PBS -o /panfs/panasas01/shared-godmc/job_report/split-o
#PBS -e /panfs/panasas01/shared-godmc/job_report/split-e
#PBS -l walltime=12:00:00
#PBS -t 1-12
# PBS -t 3
#PBS -l nodes=1:ppn=8
#PBS -S /bin/bash

set -e
if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

echo "Running on ${HOSTNAME}"

TMPDIR=/panfs/panasas01/shared-godmc/tmp

filename="/panfs/panasas01/shared-godmc/scripts/cohorts.txt"

user=`head -n ${PBS_ARRAYID} ${filename} | tail -n 1 | cut -d " " -f 1`
cohort=`head -n ${PBS_ARRAYID} ${filename} | tail -n 1 | cut -d " " -f 2`

echo $user
echo $cohort

#copy positive control
#cp /panfs/panasas01/shared-godmc/sftp/GoDMC/$user/$cohort/results/04/positive_control_pcadjusted_cg07959070_qqplot.png /panfs/panasas01/shared-godmc/methQTL_poscon/$user\_$cohort.positive_control_pcadjusted_cg07959070_qqplot.png


#merge matrixQTL to AFs

cd /panfs/panasas01/shared-godmc/meta-analysis/inputfiles
#echo SNP CpG ID BETA SE P |perl -pe 's/ /\t/g' >${user}\_${cohort}.gwama.formatted.txt
#cat ${user}\_${cohort}.*.gwama.formatted.txt >>${user}\_${cohort}.gwama.formatted.txt
#awk 'NR>1 {print $2}' <${user}\_${cohort}.gwama.formatted.txt |sort -u >${user}\_${cohort}.probes
#gzip ${user}\_${cohort}.*.gwama.formatted.txt

#split -l $(( $( wc -l < ${user}\_${cohort}.probes ) / 50 + 1 )) ${user}\_${cohort}.probes probe.${user}\_${cohort}
#ls probe.${user}\_${cohort}* >probe.${user}\_${cohort}.files
zcat /panfs/panasas01/shared-godmc/sftp/GoDMC/$user/$cohort/results/02/data.frq.gz | sed -e 's/[[:space:]]\+/ /g' |perl -pe 's/^ //g'|perl -pe 's/ /\t/g'|awk -v OFS='\t' '{ if(NR>1) print $1,$2,$3,$4,$5,$6/2; else print $0;}'|perl -pe 's/A1/EA/g' |perl -pe 's/A2/NEA/g' |perl -pe 's/MAF/EAF/g'|perl -pe 's/NCHROBS/N/g' |perl -pe 's/ /\t/g'>$cohort.frq.tmp
perl /panfs/panasas01/shared-godmc/scripts/join_file.pl -i "${user}_${cohort}.gwama.formatted.txt,TAB,0 $cohort.frq.tmp,TAB,1" -o ${user}_${cohort}.gwama.formatted.txt2 -a 1
awk -F'\t' '{ if(NR>1) print $0; else print "POS","CpG","ID","BETA","SE","P","CHR","SNP","EA","NEA","EAF","N";}' OFS='\t'<${user}\_${cohort}.gwama.formatted.txt2 | sed 's/\:SNP//1'|sed 's/\:INDEL//1' |sed 's/[^:]*://1' >${user}\_${cohort}.gwama.formatted.txt3
rm ${user}\_${cohort}.gwama.formatted.txt2

#cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01" "cg02" "cg03" "cg04" "cg05" "cg06" "cg07" "cg08" "cg09" "cg10" "cg11" "cg12" "cg13" "cg14" "cg15" "cg16" "cg17" "cg18" "cg19" "cg20" "cg21" "cg22" "cg23" "cg24" "cg25" "cg26" "cg27" "_ch")

cpgs=("_ch")
for j in ${cpgs[@]}; do
echo $j
grep $j ${user}_${cohort}.gwama.formatted.txt3 > ${user}_${cohort}.gwama.formatted.$j.txt

k=$j
if [ $k = "_ch" ]; then
k="ch"
fi

grep $k ${user}\_${cohort}.probes > ${user}\_${cohort}.$j.probes
done

ls ${user}\_${cohort}.*.probes >probe.${user}\_${cohort}.files


#DON'T USE THIS
#awk 'BEGIN {FS = "\t"}; {print > ("'$user'_'${cohort}.'"$2".txt")}' <${user}\_${cohort}.gwama.formatted.txt
#gzip ${user}\_${cohort}.gwama.formatted.txt
#mkdir -p $cohort
#cd $cohort
#awk 'BEGIN {FS = "\t"}; {print > ("'$user'_'${cohort}.'"$2".txt")}' </panfs/panasas01/shared-godmc/meta-analysis/inputfiles/${user}\_${cohort}.gwama.formatted.txt

#split -l $(( $( wc -l < ${user}\_${cohort}.probes ) / 100 + 1 )) ${user}\_${cohort}.probes probe.${user}\_${cohort}
#ls probe.${user}\_${cohort}* >probe.${user}\_${cohort}.files
