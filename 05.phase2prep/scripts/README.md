### Count across cohorts

## combinephase1cohorts.sh

- This script combines the cohort lists of SNP-CpG pairs in a list with counts of SNP-CpG pair.

- For example:
```
head /panfs/panasas01/shared-godmc/counts/combined/trans.1e-05_cg27[5-9].allcohorts.txt
17 chr1:226598314:SNP_cg27539482
```

SNP CpG pair `chr1:226598314:SNP_cg27539482` is found in 17 cohorts.

###INPUT

`/panfs/panasas01/shared-godmc/counts/${user}_${cohort}/trans.assoc.${user}_${cohort}.{pvalthreshold}.${cpgsubset}.txt`
`/panfs/panasas01/shared-godmc/counts/${user}_${cohort}/cis.assoc.${user}_${cohort}.{pvalthreshold}.${cpgsubset}.txt`
###OUTPUT

`/panfs/panasas01/shared-godmc/counts/combined/trans.{pvalthreshold}_${cpgsubset}.allcohorts.txt`
`/panfs/panasas01/shared-godmc/counts/combined/cis.{pvalthreshold}_${cpgsubset}.allcohorts.txt`

## countsacrosscohorts.R

This script generates plots such as the ones in 

`/panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/05.phase2prep/data`