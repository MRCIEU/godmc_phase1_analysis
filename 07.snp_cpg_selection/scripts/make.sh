#!/bin/bash

rm ../lists/*
Rscript make.R
cd ../lists
tar cv * | gzip -n > ../lists_17.tgz
cd ../
md5sum lists_17.tgz > lists_17.tgz.md5sum
