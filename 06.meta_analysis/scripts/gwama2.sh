MYDIR=/panfs/panasas01/shared-godmc/meta-analysis/inputfiles
outdir="/panfs/panasas01/shared-godmc/meta-analysis/output"
#outdir="/projects/MRC-IEU/groups/godmc/meta-analysis/output"
#traits='$MYDIR/traits_all.txt' 

#i=0
#while [ $i -lt "${#traits[*]}" ]; do

while read line
do
    trait=$line
    echo "$trait"
outfile="${trait}.3cohorts"
echo $outfile

if [ ! -f $outdir/$outfile.out ]; then
    

awk 'NR==1 {print $0}' <$MYDIR/epzjlm_ARIES.gwama.formatted.txt3 >$MYDIR/epzjlm_ARIES.gwama.in
awk -v mytrait="$trait" 'BEGIN{FS=OFS="\t"} $2==mytrait {print $0}'<$MYDIR/epzjlm_ARIES.gwama.formatted.txt3 >>$MYDIR/epzjlm_ARIES.gwama.in

awk 'NR==1 {print $0}' <$MYDIR/rluijk_Leiden_Longevity_Study.gwama.formatted.txt3 >$MYDIR/rluijk_Leiden_Longevity_Study.gwama.in
awk -v mytrait="$trait" 'BEGIN{FS=OFS="\t"} $2==mytrait {print $0}'<$MYDIR/rluijk_Leiden_Longevity_Study.gwama.formatted.txt3 >>$MYDIR/rluijk_Leiden_Longevity_Study.gwama.in

awk 'NR==1 {print $0}' <$MYDIR/ehannon_Phase1_SCZ.gwama.formatted.txt3 >$MYDIR/ehannon_Phase1_SCZ.gwama.in
awk -v mytrait="$trait" 'BEGIN{FS=OFS="\t"} $2==mytrait {print $0}'<$MYDIR/ehannon_Phase1_SCZ.gwama.formatted.txt3 >>$MYDIR/ehannon_Phase1_SCZ.gwama.in


touch $MYDIR/meta-analysis.files.in
files=`ls $MYDIR/*gwama.in`
echo $files > $MYDIR/meta-analysis.files.in
perl -pe 's/ /\n/g' <$MYDIR/meta-analysis.files.in >$MYDIR/meta-analysis.files.in2
mv $MYDIR/meta-analysis.files.in2 $MYDIR/meta-analysis.files.in




GWAMA -i $MYDIR/meta-analysis.files.in -o $outdir/$outfile -qt --indel_alleles --name_marker MARKER
rm $MYDIR/meta-analysis.files.in

#let i++
#done
fi

done < $MYDIR/traits_all.txt

#done < $MYDIR/test.txt
