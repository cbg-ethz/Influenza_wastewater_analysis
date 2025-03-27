library(tidyverse)
library(data.table)

#TODO: needs cleaning
setwd("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_5/")

source("/Users/anjohn/Desktop/euler/wastewater/Influenza/AA_mutationAnalysis_IAV/mutationTranslation.R")
source("/Users/anjohn/Desktop/euler/wastewater/Influenza/AA_mutationAnalysis_IAV/readingFunctions.R")
# source("/Users/anjohn/projects/wastewater/influenza/readingFunctions.R")


dir_euler <- '/Users/anjohn/Desktop/euler/wastewater/Influenza/work-IA_H1/v-pipe/'
fasta_file <- paste0(dir_euler,"references/nextstrain/H1.fasta")
results_aire <-  "results_Aire/"


start = 72
end = 1718

samples_aire <- read.table('../sample_mapping/samples_Info_timeseries.csv',
                           sep = ",", header = TRUE)
samples_aire$Barcode <- paste0("barcode",str_trim(samples_aire$Barcode))
samples_dir <- paste0(samples_aire$Barcode,"/","batch1")
var_files_aire <- paste0(dir_euler,results_aire,samples_dir,"/variants_lofreq/SNVs/snvs.vcf")

#Total drop-out could not compute anything there
filtered_var_files_aire <- var_files_aire[!grepl("barcode13", var_files_aire)]




var_list_aire <- filtered_var_files_aire %>% 
  map(function(x) read_and_mark_vcf(x, regex_sample = "barcode\\d{2}"))

var_list_aire <- var_list_aire[!sapply(var_list_aire, is.null)]


HA_mutAA_aire <- var_list_aire %>% 
  map(function(x) identify_AA(fasta_file,start,end,x)) %>% 
  bind_rows() %>%
  merge(.,samples_aire[c('Barcode','Season','Place','Date')],
        by.x=c('sample_name'),
        by.y=c('Barcode'))


setDT(HA_mutAA_aire)

#don't care about synonymous mutations
dt_aire <- HA_mutAA_aire[new_AA != org_AA]
dt_aire[, mut_lable := paste0(org_AA,gene_codonIDX,new_AA)]
dt_aire[, Date := as.Date(Date, format="%Y_%m_%d")]
dt_aire[, loc_lable := sapply(str_split(Place, pattern = " "), `[`, 2)]
dt_aire[["week"]] = strftime(as.Date(dt_aire$Date, format="%Y_%m_%d"), 
                             format = "%V")
dt_aire[["year"]] = strftime(as.Date(dt_aire$Date, format="%Y_%m_%d"), 
                             format = "%Y")

dt_aire[["year"]] = factor(dt_aire[["year"]],levels=c("2024","2023"))

#very low coverage after that
dt_clean_H1 = dt_aire %>% filter(POS <= 1032)
