
#!/bin/sh

#$1 = path to results directory (v-pipe output) of interest 
target_dir=$1
[[ "$target_dir=" =~ (H(1|3)|N(1|2)|M) ]] && segment="${BASH_REMATCH[0]}"

find $target_dir -type f -name "coverage.tsv.gz" | while read file; do
    zcat "$file" | awk 'BEGIN {OFS="\t"} NR==1 {header=$3; print $0, header} NR>1 {print $0, header}' | tail -n +2
done | gzip > $target_dir"/"$segment"_coverages_combined.tsv.gz"

