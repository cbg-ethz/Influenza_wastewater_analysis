library(data.table)
library(tidyverse)
library(seqinr)
library(ggmsa)
library(RColorBrewer)


# setwd("path_to_Figure_5/")
source("helper_functions.R")

dir <- '../../work-IA_M/v-pipe/results_timeseries'

fasta_files = list.files(path = dir, 
              full.names = TRUE, 
              recursive = TRUE,
              pattern="*support.fas")

sequences = map(fasta_files,fasta_seq) %>% do.call(c,.)

#M2e is based on splicing event (2 parts), ORF was based on sources below
#1. https://www.mdpi.com/2076-393X/3/1/105
#2. https://www.ncbi.nlm.nih.gov/nuccore/KU933488.1

#default settings of filter_seq() --> posterior=0.99,ave_reads=7
#Only high confident haplotypes (VILOCA output) are considered

vac_side_1 = subseq(sequences, start=23, end=47)

#sequences with "-" cannot be converted to AA
#only one sequence with "-" so neglectable
filtered_sequences <- vac_side_1[!grepl("-", as.character(sequences))]
M2e_side1_AA = filter_seq(filtered_sequences,"1-201",aa_conversion = TRUE)
M2e_side1_AA_cons = aa_consensus(M2e_side1_AA)

vac_side_2 = subseq(sequences, start=65, end=109)
M2e_side2_AA = filter_seq(vac_side_2 ,"671-871",aa_conversion = TRUE)
M2e_side2_AA_cons = aa_consensus(M2e_side2_AA)

#sequences were always the same so can just join spliced part
M2e_ww_consensus =paste0(M2e_side1_AA_cons,M2e_side2_AA_cons)


#### Clinical consensus sequence #####

#Function to generate consensus sequence from gisaid fasta and metadata files.
#EPI IDs of final sequences (after filtering by timeframe and location) are 
#recorded in SourceData_3.txt

# x = making_M_consensus("Intermediate_DTs/Gisaid_AbundanceAnalysis/sequences.fasta",
#                        "Intermediate_DTs/Gisaid_AbundanceAnalysis/metadata.tsv.zst")

# writeXStringSet(DNAStringSet(x), 'M_clinical_consensus.fa')

consensus_sequence = readDNAStringSet('M_clinical_consensus.fa')

Me2_1 = Biostrings::translate(subseq(consensus_sequence, start=4, end=27))
Me2_2 = Biostrings::translate(subseq(consensus_sequence, start=733, end=777))

Me2_consensus_CH = paste0(Me2_1,Me2_2)



#### Literature consensus sequence #####

#https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-019-0572-3

#Fig1, lineage 3 in Mezhenskaya et al.
Me2_consensus_lit1 = "SLLTEVETPTRSEWECRCSGSSD"
#Fig1, lineage 5 in Mezhenskaya et al.
Me2_consensus_lit2 = "SLLTEVETPIRNGWECKCNDSSD" 


#### Comparison ww, clinical and literature consensus ####

M2e_comparison = AAStringSet(c(M2e_ww_consensus,Me2_consensus_CH, Me2_consensus_lit1,Me2_consensus_lit2))
names(M2e_comparison) = c("M2e Wastewater","M2e Clinical","M2e Literature 1","M2e Literature 2")

my_custom <- data.frame(names = c(LETTERS[1:25], "-"), 
                        color =  c(
                          "#E0FFFF", "#A6CEE3", "#EE799F", "#5CACEE", "#B4EEB4",
                          "#E0FFFF", "#FFD700", "#E31A1C", "#FFB6C1", "#FDBF6F",
                          "#43CD80", "#FFEC8B", "#6A3D9A", "#FF6A6A", "#FF0000",
                          "#CD96CD", "#E0FFFF", "#B0E2FF", "#40E0D0", "#FFC125",
                          "#E0FFFF", "#99C2FF", "#EED5B7", "#D9E5F1", "#8DD3C7",
                          "#E6B3D4"
                        ), 
                        stringsAsFactors = FALSE)

my_custom_2 <- data.frame(names = c(LETTERS[1:25], "-"), 
                          colors <- c(
                            "#FF6F61", "#6B5B95", "#43CD80", "#F7CAC9", "#92A8D1", "#F0E68C", 
                            "#FFB6C1", "#D1C6B1", "#F1D6A1", "#FFD700", "#D3B897", 
                            "#98FB98", "#40E0D0", "#EE82EE", "#DAA520", "#ADFF2F", 
                            "#20B2AA", "#F08080", "#B0E0E6", "#FFD700", "#D2691E", "#8A2BE2", 
                            "#00FA9A", "#00BFFF", "#FF1493", "#A52A2A"
                          ), 
                        stringsAsFactors = FALSE)

letters = c("S", "L", "T", "E", "V", "P", "I", "R", "N", "G", "W", "C", "D","K")
# colors <- brewer.pal(12, "Paired")
# colors[12] <- "#00BFFF"
extended_colors <- c(brewer.pal(12, "Set3"), "#FF6F61", "#6B5B95")

my_custom_3 <- data.frame(names = letters, 
                          colors = extended_colors, 
                          stringsAsFactors = FALSE)


p <- ggmsa(M2e_comparison, seq_name = TRUE,color = "LETTER",custom_color = my_custom_3 ) + 
  # geom_seqlogo(adaptive = FALSE) +
  xlab("AA position") +
  theme(
    axis.text = element_text(size = 18),
    axis.title.x = element_text(size = 22,margin = margin(t=10))
      
  )

p

ggsave("/Users/anjohn/Desktop/manuscripts/Influenza/Figures/FIG5_B.png",
       plot = p, width = 10, height = 5, units = "in", dpi = 600)
