#!/bin/bash

#PBS -N betacomp
#PBS -o /panfs/panasas01/shared-godmc/job_report/betacomp-output
#PBS -e /panfs/panasas01/shared-godmc/job_report/betacomp-error
#PBS -l walltime=12:00:00
# PBS -t 1-3
#PBS -l nodes=1:ppn=16
#PBS -S /bin/bash

set -e

echo "Running on ${HOSTNAME}"

#cd /panfs/panasas01/shared-godmc/meta-analysis/betacomparison/
#R CMD BATCH betacomparison_3cohorts.R /panfs/panasas01/shared-godmc/Rout/betacomparison_3cohorts.Rout
#R CMD BATCH betacomparison_3cohorts2.R /panfs/panasas01/shared-godmc/Rout/betacomparison_3cohorts2.Rout
#R CMD BATCH betacomparison_3cohorts3.R /panfs/panasas01/shared-godmc/Rout/betacomparison_3cohorts3.Rout

#R CMD BATCH /panfs/panasas01/shared-godmc/scripts/betacomp2.R /panfs/panasas01/shared-godmc/Rout/betacomp2.Rout
#R CMD BATCH /panfs/panasas01/shared-godmc/scripts/betacompplot.R /panfs/panasas01/shared-godmc/Rout/betacomp2.Rout
cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/strand_effectsize/data
R CMD BATCH /panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/strand_effectsize/scripts/betacompplot.R /panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/strand_effectsize/scripts/betacompplot.Rout

