# mirror results from sftp to shared-godmc

1. Mirror the phase 1 results

```
ssh epzjlm@sscmv-filetran.epi.bris.ac.uk
cd  /srv/sftponly
sudo su
rsync -av GoDMC epzjlm@bluecrystalp3.bris.ac.uk:/panfs/panasas01/shared-godmc/sftp/
```

2. Extract data from cohort tgz files

All cohorts from `/panfs/panasas01/shared-godmc/scripts/cohorts.txt` can be extracted

```
sh extract1-5.sh 1
```

