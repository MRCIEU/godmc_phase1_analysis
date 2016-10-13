# mirror results from sftp 

1. Mirror the phase 1 results to the RDSF

```
ssh epzjlm@sscmv-filetran.epi.bris.ac.uk
cd  /srv/sftponly
sudo su
rsync -Ov GoDMC epzjlm@bluecrystalp3.bris.ac.uk://projects/MRC-IEU/groups/godmc/sftp/
```

2. Extract data from cohort tgz files

All cohorts from `/panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/01.extract_sftp/data/cohorts.txt` can be extracted. 

Data needs to be copied to shared-godmc:

rsync -av //projects/MRC-IEU/groups/godmc/sftp/GoDMC/ /panfs/panasas01/shared-godmc/sftp/GoDMC/

To extract data from the first cohort:

```
rsync -av //projects/MRC-IEU/groups/godmc/sftp/GoDMC/epzjlm/ARIES_0[1-5].tgz /panfs/panasas01/shared-godmc/sftp/GoDMC/epzjlm/
sh extract1-5.sh 1
```
Data will be written to: `/panfs/panasas01/shared-godmc/sftp/GoDMC/$user/$cohort`

3. Extract samplesizes:

Please note that ./results/01/cohort_descriptives.RData gives incorrect sample sizes.

Please run:

```
extractN.sh
```

It will give you a file called:

`/panfs/panasas01/sscm/epzjlm/repo/godmc_phase1_analysis/01.extract_sftp/data/cohort_samplesizes.txt`

4. I run the following checks:

- check plots in `/panfs/panasas01/shared-godmc/sftp/GoDMC/$user/$cohort/results/01`
- check cohort_descriptives.RData in `/panfs/panasas01/shared-godmc/sftp/GoDMC/$user/$cohort/results/01`
- check pcaplot in `/panfs/panasas01/shared-godmc/sftp/GoDMC/$user/$cohort/results/02/pcaplot.pdf`


