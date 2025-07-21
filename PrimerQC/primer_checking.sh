#!/bin/bash

#environment specified by environment.yaml
source primerQC

#Made with consensus_buliding.ipynb
declare -A refs=(
  [H1]=HA_H1N1_consensus.fasta
  [H3]=HA_H3N2_consensus.fasta
  [N1]=NA_H1N1_consensus.fasta
  [N2]=NA_H3N2_consensus.fasta
  [M_H1N1]=M_H1N1_consensus.fasta
  [M_H3N2]=M_H3N2_consensus.fasta
)

# Based on Source Data 4
declare -A inputs=(
  [H1]=H1.fasta
  [H3]=H3.fasta
  [N1]=N1.fasta
  [N2]=N2.fasta
  [M_H1N1]=M_H1N1.fasta
  [M_H3N2]=M_H3N2.fasta
)

mkdir -p results


for seg in H1 H3 N1 N2 M_H1N1 M_H3N2; do
  ref=${refs[$seg]}
  input=${inputs[$seg]}
  bam_out=${seg}.consensus_23TO24.bam
  vcf_out=${seg}.consensus_23TO24.vcf

  echo "Processing $seg..."

  bwa index "$ref"

  bwa mem -t 4 "$ref" "$input" | samtools sort -@ 4 -o "$bam_out"

  samtools index "$bam_out"

  lofreq call -f "$ref" "$bam_out" -o "$vcf_out"

  echo "$seg processing complete."

done

