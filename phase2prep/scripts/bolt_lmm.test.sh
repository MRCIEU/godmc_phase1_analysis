#!/bin/bash

#PBS -N combinecounts
#PBS -o /panfs/panasas01/shared-godmc/job_report/boltlmm-output
#PBS -e /panfs/panasas01/shared-godmc/job_report/boltlmm-error
#PBS -l walltime=12:00:00
# PBS -t 1-10
#PBS -l nodes=1:ppn=16
#PBS -S /bin/bash


#A minimal BOLT-LMM invocation looks like:
#./bolt --bfile=geno --phenoFile=pheno.txt --phenoCol=phenoName 
#     --lmm --LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz 
#     --statsFile=stats.tab

module add apps/bolt-lmm-2.2
#bolt

mydir="/panfs/panasas01/shared-godmc/counts"

cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01[0-4]" "cg01[5-9]" "cg02[0-4]" "cg02[5-9]" "cg03[0-4]" "cg03[5-9]" "cg04[0-4]" "cg04[5-9]" "cg05[0-4]" "cg05[5-9]" "cg06[0-4]" "cg06[5-9]" "cg07[0-4]" "cg07[5-9]" "cg08[0-4]" "cg08[5-9]" "cg09[0-4]" "cg09[5-9]" "cg10[0-4]" "cg10[5-9]" "cg11[0-4]" "cg11[5-9]" "cg12[0-4]" "cg12[5-9]" "cg13[0-4]" "cg13[5-9]" "cg14[0-4]" "cg14[5-9]" "cg15[0-4]" "cg15[5-9]" "cg16[0-4]" "cg16[5-9]" "cg17[0-4]" "cg17[5-9]" "cg18[0-4]" "cg18[5-9]" "cg19[0-4]" "cg19[5-9]" "cg20[0-4]" "cg20[5-9]" "cg21[0-4]" "cg21[5-9]" "cg22[0-4]" "cg22[5-9]" "cg23[0-4]" "cg23[5-9]" "cg24[0-4]" "cg24[5-9]" "cg25[0-4]" "cg25[5-9]" "cg26[0-4]" "cg26[5-9]" "cg27[0-4]" "cg27[5-9]" "_ch")

cpgs=("cg0000[0-9]")
pvals=("1e-13")
no=1

for i in ${pvals[@]}; do
echo $i

for j in ${cpgs[@]}; do
echo $j

cat $mydir/combined/trans.$i\_${j}.allcohorts.txt $mydir/combined/cis.$i\_${j}.allcohorts.txt | awk '$1>='$no' {print $0}' |perl -pe 's/\_/ /g'|awk '{print $2,$3}'  >$mydir/combined/cis_trans.$i\_${j}.gt${no}.allcohorts.txt
awk '{print $2}'<$mydir/combined/cis_trans.$i\_${j}.gt${no}.allcohorts.txt | sort -u >$mydir/combined/cis_trans.$i\_${j}.gt${no}.allcohorts.probes

filename=$mydir/combined/cis_trans.$i\_${j}.gt${no}.allcohorts.probes

#prepare probes

Rscript resources/methylation/extractprobesubsets.R \
		${betas} \
		${covariates_combined}.txt \
        $mydir/combined/cis_trans.$i\_${j}.gt${no}.allcohorts.probes \
		${bfile}.fam \
		${methylation_subset}.$i\_${j}.gt${no}.txt
        ${covariates_combined}.bolt_lmm



while read -r line
do
    probe="$line"
    echo $probe

    fgrep $probe $mydir/combined/cis_trans.$i\_${j}.gt${no}.allcohorts.txt |awk '{print $1}'> $mydir/combined/cis_trans.$i\_${j}.gt${no}.allcohorts.snps
    grep -v -w $mydir/combined/cis_trans.$i\_${j}.gt${no}.allcohorts.snps ${bfile}.bim |awk '{print $2}' > $mydir/combined/snp.exclude.$i\_${j}.gt${no}
    nocols=`awk -F' ' '{print NF; exit}'< ${covariates_combined}.bolt_lmm`
    echo $nocols
    qcov=`grep -o QCOV ${covariates_combined}.bolt_lmm |wc -l`
    catcov=`grep -o CATCOV ${covariates_combined}.bolt_lmm |wc -l`

 bolt \
    --bfile=${bfile} \ 
    --phenoFile=${methylation_subset}.$i\_${j}.gt${no}.txt \
    --phenoCol=$probe \
    --exclude $mydir/combined/snp.exclude.$i\_${j}.gt${no} \
    --covarFile ${covariates_combined}.bolt_lmm \
    --covarCol {1:catcov}\
    --qCovarCol {1:qcov}\
    --lmmForceNonInf \
    --LDscoresFile=tables/LDSCORE.1000G_EUR.tab.gz \
    --statsFile=stats.tab.$i\_${j}.gt${no}
   

   #LC_ALL=C awk -v mytrait="$probe" 'BEGIN{FS=OFS="\t"} $2==mytrait {print $1}'<$prefix.gwama.formatted.txt > $myout/$prefix/$prefix.$probe.snps
    
    #time (awk -v mytrait="$probe" 'BEGIN{FS=OFS="\t"} $2==mytrait {print $0}'<$prefix.gwama.formatted.txt >> $myout/$prefix/$prefix.$probe.in)
    #real  1m5.080s
    #user  1m0.498s
    #sys 0m3.807s
    
    #time (LC_ALL=C fgrep $probe $prefix.gwama.formatted.txt >> $myout/$prefix/$prefix.$probe.in)
    #real    0m11.420s
    #user    0m3.911s
    #sys 0m2.744s
    
    #time (fgrep $probe $prefix.gwama.formatted.txt >> $myout/$prefix/$prefix.$probe.in)
    #real 0m10.975s
    #user    0m3.841s
    #sys 0m2.628s
    
    #time(fgrep $probe $mydir/$prefix.gwama.formatted.${j}.txt >> $myout/$prefix/$prefix.$probe.assoc)

    #real    0m0.095s
    #user    0m0.052s
    #sys 0m0.029s  
      fi
done < "$filename"
    
done
done