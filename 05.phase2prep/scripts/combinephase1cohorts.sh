#!/bin/bash

#PBS -N combinecounts_noRAINE
#PBS -o /panfs/panasas01/shared-godmc/job_report_2017/combinecounts_noRAINE-output
#PBS -e /panfs/panasas01/shared-godmc/job_report_2017/combinecounts_noRAINE-error
#PBS -l walltime=12:00:00
# PBS -t 1-10
#PBS -l nodes=1:ppn=12
#PBS -S /bin/bash


#I had an error something like, "/tmp/sortA3aLjF: No space left on device"
#The shared-godmc area had enough space and but perhaps it is temporary space on the node that is running out when executing the 'sort' command.
#Sort creates temporary files when working on large input files.  These are typically placed in /tmp.  However, space in /tmp can run out if the inputs are large enough.  You can change where sort places these temporary files.  You can do this by adding, say,
#TMPDIR=/panfs/panasas01/shared-godmc/tmp
#on a line after '#!/bin/bash' 

TMPDIR=/panfs/panasas01/shared-godmc/tmp

#filename="/panfs/panasas01/sscm/epwkb/GoDMC_Analysis/CpG_SNP_FINAL/repo/godmc_phase1_analysis/01.extract_sftp/data/cohorts_noRAINE.txt"
#filename="/panfs/panasas01/sscm/epwkb/GoDMC_Analysis/CpG_SNP_FINAL/repo/godmc_phase1_analysis/01.extract_sftp/data/cohorts.txt"
filename="/panfs/panasas01/shared-godmc/scripts/cohorts_noRAINE.txt"

mydir="/panfs/panasas01/shared-godmc/counts_2017"

mkdir -p /panfs/panasas01/shared-godmc/counts_2017/combined


pvals=("1e-05" "1e-06" "1e-07" "1e-08" "1e-09" "1e-10" "1e-11" "1e-12" "1e-13")
#pvals=("1e-13")

cpgs=("cg0000[0-9]" "cg0001" "cg0002" "cg0003" "cg0004" "cg0005" "cg0006" "cg0007" "cg0008" "cg0009" "cg001" "cg002" "cg003" "cg004" "cg005" "cg006" "cg007" "cg008" "cg009" "cg01[0-4]" "cg01[5-9]" "cg02[0-4]" "cg02[5-9]" "cg03[0-4]" "cg03[5-9]" "cg04[0-4]" "cg04[5-9]" "cg05[0-4]" "cg05[5-9]" "cg06[0-4]" "cg06[5-9]" "cg07[0-4]" "cg07[5-9]" "cg08[0-4]" "cg08[5-9]" "cg09[0-4]" "cg09[5-9]" "cg10[0-4]" "cg10[5-9]" "cg11[0-4]" "cg11[5-9]" "cg12[0-4]" "cg12[5-9]" "cg13[0-4]" "cg13[5-9]" "cg14[0-4]" "cg14[5-9]" "cg15[0-4]" "cg15[5-9]" "cg16[0-4]" "cg16[5-9]" "cg17[0-4]" "cg17[5-9]" "cg18[0-4]" "cg18[5-9]" "cg19[0-4]" "cg19[5-9]" "cg20[0-4]" "cg20[5-9]" "cg21[0-4]" "cg21[5-9]" "cg22[0-4]" "cg22[5-9]" "cg23[0-4]" "cg23[5-9]" "cg24[0-4]" "cg24[5-9]" "cg25[0-4]" "cg25[5-9]" "cg26[0-4]" "cg26[5-9]" "cg27[0-4]" "cg27[5-9]" "_ch")
#cpgs=("cg0000[0-9]")

for i in ${pvals[@]}; do
echo $i

for j in ${cpgs[@]}; do
echo $j

    cd $mydir
    
    date
    #escape brackets
    k=$(sed -e 's/[][?*]/\\&/g' <<< "$j")
    echo $k

    #generate cat list
    
    #list="$(find $mydir -mindepth 2 -type f -name "cis.assoc.*.${i}."${k}".txt")"
	list="$(find $mydir -mindepth 2 -type f -name "cis.assoc.*.${i}."${k}".txt" ! -name "cis.assoc.cpennell_Raine.${i}."${k}".txt")"
    echo $list

    #concatenate and count
    cat $list | awk '{c[$0]++}END{for(l in c){print c[l], l}}' | sort -n -r > $mydir/combined/cis.${i}\_${j}.allcohorts.txt
    echo cis.${i}\_${j}.allcohorts.txt
    date
    
    ciscount=`wc -l $mydir/combined/cis.${i}\_${j}.allcohorts.txt`
    echo $ciscount

    date
    
    #list="$(find $mydir -mindepth 2 -type f -name "cis.assoc.*.${i}."${k}".txt")"
	list="$(find $mydir -mindepth 2 -type f -name "trans.assoc.*.${i}."${k}".txt" ! -name "trans.assoc.cpennell_Raine.${i}."${k}".txt")"
    echo $list
    
    cat $list | awk '{c[$0]++}END{for(l in c){print c[l], l}}' | sort -n -r > $mydir/combined/trans.${i}\_${j}.allcohorts.txt
    echo trans.${i}\_${j}.allcohorts.txt
    date
    
    transcount=`wc -l $mydir/combined/trans.${i}\_${j}.allcohorts.txt`
    echo $transcount
  
    
    #find $mydir -mindepth 2 -type f -name "cis.assoc.*.${i}."${j}".txt" -exec cat {} \; |perl -ne 'for (split /\s+/, $_){ $w{$_}++ }; END{ for $key (sort keys %w) { print "$key $w{$key}\n"}}' | sort -n -r > $mydir/combined/cis.${i}\_${j}.allcohorts2.txt
    #date
    #ciscount=`wc -l $mydir/combined/cis.${i}\_${j}.allcohorts2.txt`
    #echo $ciscount
   
    
    ####PROBLEMS with sort eg. not enough tmp space even after setting TMPDIR
    #date    
    #find $mydir -mindepth 2 -type f -name "cis.assoc.*.${i}."${j}".txt" -exec cat {} \;  |sort |uniq -c |sort -n -r > $mydir/combined/cis.${i}\_${j}.allcohorts3.txt
    #date
    #echo cis.${i}\_${j}.allcohorts.txt
    #ciscount=`wc -l $mydir/combined/cis.${i}\_${j}.allcohorts3.txt`
    #echo $ciscount
   
done

cd $mydir/combined
cat cis.${i}*.allcohorts.txt > cis.$i\_allcohorts.txt
cat trans.${i}*.allcohorts.txt > trans.$i\_allcohorts.txt
done

cd $mydir/combined

nocoh=`cat $filename |wc -l`
echo $nocoh

START=1
END=$nocoh


#START=7
#END=7
#filename="/panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/extract_sftp/data/cohorts.txt"
#mydir="/panfs/panasas01/shared-godmc/counts"
#pvals=("1e-05" "1e-06" "1e-07" "1e-08" "1e-09" "1e-10" "1e-11" "1e-12" "1e-13")


## save $START, just in case if we need it later ##
no=$START
echo $no
while [[ $no -le $END ]]
do
    echo "$no"
rm counts.allcohorts.combined.$no.txt
touch counts.allcohorts.combined.$no.txt

for i in ${pvals[@]}; do
echo $i

cis=`awk '$1=='$no' {print $0}' cis.${i}\_allcohorts.txt | wc -l`
trans=`awk '$1=='$no' {print $0}' trans.${i}\_allcohorts.txt |wc -l`
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
((no = $no + 1))

done

#pval cut-off
i="1e-05"
#cpgs=("cg0000[0-9]")

#found in number of cohorts
no="2"

for j in ${cpgs[@]}; do
echo $j

#cat $mydir/combined/cis.${i}\_${j}.allcohorts.txt $mydir/combined/trans.${i}\_${j}.allcohorts.txt | perl -pe 's/_/ /g' | awk '$1>='$no' {print $0}' > $mydir/combined/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt
#awk '{print $3}' <$mydir/combined/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt | sort -u > $mydir/combined/cis_trans.${i}\_${j}.ge${no}.allcohorts.probes
#rsync -av $mydir/combined/cis_trans.${i}\_${j}.ge${no}.allcohorts.probes /panfs/panasas01/sscm/epzjlm/repo/godmc/processed_data/methylation_data/
#rsync -av $mydir/combined/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt /panfs/panasas01/sscm/epzjlm/repo/godmc/processed_data/methylation_data/

#no cis =1, no trans =2

nocis="1"
perl -pe 's/_/ /g' <$mydir/combined/cis.${i}\_${j}.allcohorts.txt | awk '$1>='$nocis' {print $0}' > $mydir/combined/cis_trans.${i}\_${j}.ge${nocis}.allcohorts.cis.txt
notrans="2"
perl -pe 's/_/ /g' <$mydir/combined/trans.${i}\_${j}.allcohorts.txt| awk '$1>='$notrans' {print $0}' > $mydir/combined/cis_trans.${i}\_${j}.ge${notrans}.allcohorts.trans.txt

no="1.2"

#k=$(echo $j | perl -pe 's/\[/\\[/g'|perl -pe 's/\]/\\]/g')
#echo $k

cat $mydir/combined/cis_trans.${i}\_${j}.ge${nocis}.allcohorts.cis.txt $mydir/combined/cis_trans.${i}\_${j}.ge${notrans}.allcohorts.trans.txt > $mydir/combined/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt
awk '{print $3}' <$mydir/combined/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt | sort -u > $mydir/combined/cis_trans.${i}\_${j}.ge${no}.allcohorts.probes
gzip $mydir/combined/cis_trans.${i}\_${j}.ge${no}.allcohorts.txt

done