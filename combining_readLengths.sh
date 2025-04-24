
#!/bin/sh

#1. excecute readLength.sh on samples of interest

#$1 = path to results directory (v-pipe output) of interest 
target_dir=$1


find $target_dir -type f -name "read_length.txt" | while read file; do
cat "$file" 
done | gzip > $target_dir"/readLengths_combined.tsv.gz"

