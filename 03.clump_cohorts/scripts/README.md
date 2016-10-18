# Generating the number of independent hits at pval 1e-05

## clump_cohorts[0-8] 

Each PBS_ARRAYID is a row in `/panfs/panasas01/shared-godmc/meta-analysis/inputfiles/allprobefiles.txt` (15 rows * 47 = 705 jobs)

Each cohort is splitted in 47 subsets. This is old and should be revised to 74 probesets which are used later. Also data from an older release are used for TwinsUK and ARIES.
Clumping has been done on 15 cohorts only and is extremely time consuming although the jobs are optimised and profiled extensively. If the job is killed due to walltime you can resubmit using clump_cohorts[0-8].resubmit

```
qsub clump_cohorts1
```

###INPUT

It uses /panfs/panasas01/shared-godmc/meta-analysis/inputfiles/{$user}_{$cohort}.gwama.formatted.${cpgprobeset}.txt as input file

###output
/panfs/panasas01/shared-godmc/clump/{$user_$cohort}/{user}_{cohort}.indexSNP.{cpgprobeset}

- Please note that with a new round of clumping you should use 74 probesets rather than 47 probesets.

## extractindeploci.R

```
Rscript /panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/clump_cohorts/scripts/extractindeploci.R $user $cohort

```

- It extracts the number of independent SNPs from `/panfs/panasas01/shared-godmc/clump/{$user_$cohort}/{user}_{cohort}.indexSNP.{cpgprobeset}` and generates a `/panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/clump_cohorts/data/$user_$cohort.numberofindependentloci.Robj` file which has counts of number of independent snps, independent cis and trans snps.

## plotindeploci.R

- Plots the number of independent loci by cohort ordered by cohort samplesize. THe black line in the plot represents the number of expected hits.
