### BETA comparison


### betacomp2.R

- This script checks whether there are any strand issues. It alignes betas and effect alleles to 1000G for chromosome 20.

### INPUT

It takes `/panfs/panasas01/shared-godmc/meta-analysis/betacomparison/*.gwama.formatted.chr20.txt` as input. Please note these are extracted from old files and should be extracted from `/panfs/panasas01/shared-godmc/meta-analysis/inputfiles/{$user}_{$cohort}/${user}_${cohort}.gwama.formatted.txt` 

### OUTPUT

`/panfs/panasas01/shared-godmc/meta-analysis/betacomparison/betacomp.Robj`

- This is a list where coh[[1]] are the results for the first cohort.

### betacompplot.R 

- This script plots cohort effect sizes against each other

`R CMD BATCH /panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/04.strand_effectsize/scripts/betacompplot.R /panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/04.strand_effectsize/scripts/betacompplot.Rout`

- The current plot is saved here:

`/panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/strand_effectsize/data`

- I have done this using pvalues <1e-05 but it should be done with a more stringent pvalue threshold.