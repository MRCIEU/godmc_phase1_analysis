#!/bin/bash

#PBS -N combinecounts
#PBS -o /panfs/panasas01/shared-godmc/job_report/combinecounts-output
#PBS -e /panfs/panasas01/shared-godmc/job_report/combinecounts-error
#PBS -l walltime=12:00:00
# PBS -t 1-10
#PBS -l nodes=1:ppn=16
#PBS -S /bin/bash

TMPDIR=/panfs/panasas01/shared-godmc/tmp

filename="/panfs/panasas01/shared-godmc/scripts/cohorts.txt"

cd /panfs/panasas01/shared-godmc/counts

pvals=("1e-05" "1e-06" "1e-07" "1e-08" "1e-09" "1e-10" "1e-11" "1e-12" "1e-13")
cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01" "cg02" "cg03" "cg04" "cg05" "cg06" "cg07" "cg08" "cg09" "cg10" "cg11" "cg12" "cg13" "cg14" "cg15" "cg16" "cg17" "cg18" "cg19" "cg20" "cg21" "cg22" "cg23" "cg24" "cg25" "cg26" "cg27" "_ch")

#for i in ${pvals[@]}; do
#echo $i

#for j in ${cpgs[@]}; do
#echo $j
    
#    name=`perl -pe 's/ /\_/g'<../scripts/cohorts.txt |sed -e 's/^/cis.assoc./' |sed -e 's/$/.'${i}'.'${j}'.txt/'`
#    echo $name
    
#    #cat $name |sort |uniq -c |sort -n -r > cis.${i}\_${j}.allcohorts.txt
#    echo cis.${i}\_${j}.allcohorts.txt
    
#    name=`perl -pe 's/ /\_/g'<../scripts/cohorts.txt|sed -e 's/^/trans.assoc./' |sed -e 's/$/.'${i}'.'${j}'.txt/'`
#    echo $name   

#    #cat $name |sort |uniq -c |sort -n -r > trans.${i}\_${j}.allcohorts.txt
#    echo trans.${i}\_${j}.allcohorts.txt

#done

#cat cis.${i}*.allcohorts.txt > cis.$i\_allcohorts.txt
#cat trans.${i}*.allcohorts.txt > trans.$i\_allcohorts.txt
#done

nocoh=`cat ../scripts/cohorts.txt |wc -l`
echo $nocoh
nocoh="1"

#for (( no=1; no<=$nocoh; c++ ))
#do
#	echo $no

START=1
END=$nocoh
## save $START, just in case if we need it later ##
no=$START
while [[ $no -le $END ]]
do
    echo "$no"

touch counts.allcohorts.combined.$no.txt

for i in ${pvals[@]}; do
echo $i

cis=`awk '$1=='$no' {print $0}' cis.${i}\_allcohorts.txt | grep -v chr23 |wc -l`
trans=`awk '$1=='$no' {print $0}' trans.${i}\_allcohorts.txt |grep -v chr23 |wc -l`
echo $cis
echo $trans

cis_all=`wc -l cis.${i}\_allcohorts.txt |awk '{print $1}'`
trans_all=`wc -l trans.${i}\_allcohorts.txt|awk '{print $1}'`
cperc=$((100*$cis / $cis_all))
tperc=$((100*$trans / $trans_all))

echo "Pvalue" ${i} >>counts.allcohorts.combined.$no.txt
echo "Overlapping_cis_SNP-CpG_pairs(N)" $cis >>counts.allcohorts.combined.$no.txt
echo "Overlapping_trans_SNP-CpG_pairs(N)" $trans >>counts.allcohorts.combined.$no.txt
echo "All_cis_SNP-CpG_pairs(N)" $cis_all >>counts.allcohorts.combined.$no.txt
echo "All_trans_SNP-CpG_pairs(N)" $trans_all >>counts.allcohorts.combined.$no.txt
echo "Overlapping_cis_SNP-CpG_pairs(%)" $cperc >>counts.allcohorts.combined.$no.txt
echo "Overlapping_trans_SNP-CpG_pairs(%)" $tperc >>counts.allcohorts.combined.$no.txt

done
((no = no + 1))

done

