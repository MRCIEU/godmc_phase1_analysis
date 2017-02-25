#!/bin/bash

Rscript make.R
cd ../lists
tar czvf ../lists_17.tgz *
cd ../
md5sum lists_17.tgz > lists_17.tgz.md5sum
