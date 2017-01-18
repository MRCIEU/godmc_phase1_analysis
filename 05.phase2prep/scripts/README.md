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


## sort_and_split.sh

This script takes the 74 `*ge1.2*` association lists and splits it into ~497 roughly equal lists.

- Each list has about 125000 associations
- The files are numbered from 1-497, ordered by CpG alphabetical order
- There are no CpGs present in more than one file
- Also creates probe list for each file

Takes about 20 minutes to run. 

## annotate.R

Takes output from sort_and_split.sh and adds cis/trans annotation


## prepare_upload.sh

This will tarball the necessary files and create md5 checksum. The .tar and .md5sum files need to be copied to

```
/srv/sftponly/GoDMC/resources/phase2
```

on the GoDMC SFTP server `sscmv-filetran.epi.bris.ac.uk`.


## countsacrosscohorts.R

This script generates plots such as the ones in 

`/panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/05.phase2prep/data`
