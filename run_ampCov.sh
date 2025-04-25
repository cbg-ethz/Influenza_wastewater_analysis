#!/bin/sh

#This is an example set-up on how run the amplicon read count calculation using cov_amplicon_insert.py.
#Input parameters will need to be adjusted for every dataset and per segment.
#The cov_amplicon_insert.py was based on original code from David Dreifuss, see https://github.com/dr-david/amplicon_cov

#### Amplicon coverage for time series data set ####

bed="work-IA_H1/v-pipe/references/primers/HA_scheme.insert_mod.bed"
samples="work-IA_H1/v-pipe/samples_timeseries.tsv"

#target_dir needs to be the 
target_dir="work-IA_H1/v-pipe/results_timeseries"
out_dir="work-IA_H1/v-pipe/results_timeseries"

python cov_amplicon_insert.py -r $bed -s $samples -p $target_dir -o $out_dir
