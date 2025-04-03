#!/bin/sh

BASE_DIR=xxx #path to raw data directory, formatted as described under README 1.2 ii)


for dir in "$BASE_DIR"/barcode*/ ; do
  
  dirname=$(basename "$dir")
  TARGET_DIR="$BASE_DIR/$dirname/batch1/raw_data"
  fastq="$BASE_DIR/$dirname/batch1/raw_data/$dirname.fastq.gz"

  echo $fastq

  zcat $fastq | awk '{if(NR%4==2) print length($1)}' | sort -n  > "$TARGET_DIR/read_length.txt"
  
  echo "$TARGET_DIR/read_length.txt"
#   mkdir -p "$TARGET_DIR"
  
#   cat "$dir"/* > "$TARGET_DIR/$dirname.fastq.gz"

done
