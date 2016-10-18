# Counting the number of hits at different thresholds (cis/trans)

## mqtl.count.R 

Each PBS_ARRAYID is a cohort number in /panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/01.extract_sftp/data/cohorts.txt

###INPUT
- reads in i matrixqtl files

###OUTPUT
- generates `/panfs/panasas01/shared-godmc/meta-analysis/inputfiles/{user}.{cohort}.{i}.gwama.formatted.txt` which contains snp, cpg, beta,se, p
- generates `/panfs/panasas01/shared-godmc/counts/{user}_{cohort}/cis.assoc.{i}.{user}_{cohort}.{pvaluethreshold}.txt` which contains a list of SNP_CpG at a particular pval threshold (1e-5 till 1e-13)
- generates `/panfs/panasas01/shared-godmc/counts/{user}_{cohort}/trans.assoc.{i}.{user}_{cohort}.{pvaluethreshold}.txt`which contains a list of SNP_CpG at a particular pval threshold (1e-5 till 1e-13)
- cis has been defined as 1Mb between SNP and CpG
- it generates `/panfs/panasas01/shared-godmc/cohort_summary/mqtlcount.{user}_{cohort}.txt` with counts for each i file for each pvalue threshold
- it generates `/panfs/panasas01/shared-godmc/cohort_summary/mqtlcount.{user}_{cohort}.summary.txt` which generates summary counts for each pvalue threshold


## submit_mqtl.count.sh
- runs `Rscript mqtl.count.R $user $cohort`
- combines `{cis/trans}.assoc.*.${user}_${cohort}.{pvaluethreshold}.txt` to `{cis/trans}.assoc.${user}_${cohort}.{pvaluethreshold}.txt` eg. combines all chunks to one file
- split files in j CpG subsets `{cis/trans}.assoc.${user}\_${cohort}.{pvaluethreshold}.${j}.txt` eg. 74 subsets based on CpG name

- merge `${user}_${cohort}.gwama.formatted.txt` to `/panfs/panasas01/shared-godmc/sftp/GoDMC/$user/$cohort/results/02/data.frq.gz` and generate `${user}\_${cohort}.gwama.formatted.txt` which now contains POS,CpG,ID,BETA,SE,P,CHR,SNP,EA,NEA,EAF,N
- splits files in j CpG subsets `${user}_${cohort}.gwama.formatted.txt` -> `/panfs/panasas01/shared-godmc/meta-analysis/inputfiles/${user}_${cohort}.gwama.formatted.$j.txt` which can be used for clumping, meta-analysis and all further phase 1 analysis.

## countscohortscombined.R
- generates plots from  `/panfs/panasas01/shared-godmc/cohort_summary/mqtlcount.{user}_{cohort}.summary.txt` in `/panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/02.extract_counts/data`

##join_file.pl
- You should use it as below:

```
join_file.pl –i “infile1.txt,SPACE,1 infile2,COMMA,3” –a 1 –o outfile.txt
``` 

- It means: the join is based on infile1’s 2nd column (count start from 0), which is SPACE separated; the join is based on infile2’s 4th column, which is COMMA separated.
“-a 1 “ means to include all rows from infile1, even for those having no match in infile2.


