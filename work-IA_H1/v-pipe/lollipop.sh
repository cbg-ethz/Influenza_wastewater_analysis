#!/bin/sh

mkdir results_timeseries/Lollipop/

lollipop generate-mutlist --output results_timeseries/Lollipop/mutlist.tsv  -- vocs/5a2a1_mutations_full.yaml vocs/5a2a_mutations_full.yaml

#location and data of each sample documented under Plot_scripts/sample_mapping
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode01.mut.tsv --location "GE" --date "2023-12-09" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode01/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode02.mut.tsv --location "GE" --date "2023-12-12" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode02/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode03.mut.tsv --location "GE" --date "2023-12-19" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode03/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode04.mut.tsv --location "GE" --date "2023-12-25" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode04/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode05.mut.tsv --location "GE" --date "2024-01-04" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode05/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode06.mut.tsv --location "GE" --date "2024-01-10" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode06/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode07.mut.tsv --location "GE" --date "2024-01-15" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode07/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode08.mut.tsv --location "GE" --date "2024-01-23" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode08/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode09.mut.tsv --location "GE" --date "2024-01-30" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode09/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode10.mut.tsv --location "GE" --date "2024-02-06" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode10/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode11.mut.tsv --location "GE" --date "2024-02-15" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode11/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode12.mut.tsv --location "GE" --date "2024-02-20" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode12/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode13.mut.tsv --location "GE" --date "2024-02-27" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode13/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode14.mut.tsv --location "GE" --date "2024-03-06" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode14/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode15.mut.tsv --location "GE" --date "2024-03-12" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode15/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode16.mut.tsv --location "GE" --date "2024-03-20" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode16/batch1/alignments/basecnt.tsv.gz
lollipop getmutations from-basecount --based 1 --output results_timeseries//Lollipop/barcode17.mut.tsv --location "GE" --date "2024-03-28" -m results_timeseries//Lollipop/mutlist.tsv -- results_timeseries//barcode17/batch1/alignments/basecnt.tsv.gz


xsv cat rows --output results_timeseries//Lollipop/tallymut.tsv results_timeseries//Lollipop/barcode*.mut.tsv

lollipop deconvolute --output=results_timeseries/Lollipop/deconvoluted.tsv --out-json=results_timeseries/Lollipop/deconvoluted_upload.json --var=variants_conf.yaml --vd=variants_dates.yaml --dec=deconv_bootstrap_cowwid.yaml --seed=42 --location "GE" --n-cores=8 -- results_timeseries/Lollipop/tallymut.tsv