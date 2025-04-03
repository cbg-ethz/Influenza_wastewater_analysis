library(data.table)
library(tidyverse)
library(boot)
library(ggplot2)

setwd("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_3/")

#Metadata as provided by Gisaid
meta_dt = fread("zstdcat metadata.tsv.zst")


#### Calculating clinical abundance estimates #####

#Filtering criteria which were applied. Find the IDs of the filtered metadata
#in Source Data 2

meta_dt[["Collection_Date"]] = strftime(as.Date(meta_dt$Collection_Date, 
                                                format="%Y-%m-%d"))

# Filter for dates between Dec 23 and Mar 24 => time-series sample in this time
meta_filt_dt = meta_dt[Collection_Date >= as.Date('2023-12-01') & 
                         Collection_Date <= as.Date('2024-03-31') & 
                         grepl("Europe", Location)]

IDs_dt = meta_filt_dt[, c("Isolate_Id")]

meta_filt_dt = meta_filt_dt[, c("Collection_Date", "Isolate_Id","Clade")]

s = readDNAStringSet("sequences.fasta")
seq_name = names(s)

seq_dt <- as.data.table(data.frame(seq_name))

seq_dt[, c("EPI_ID") := str_extract(seq_name,"EPI_ISL_[0-9]+") ]
#unique segment names found, not sure if they are all correct, but needed ones are there
seq_dt[, c("segment") := str_extract(seq_name,"NA|PB2|PB1|MP|NP|NS|PA|HA|HE|P3") ]

#only HA needed for typing
dt_typing = seq_dt[segment == 'HA', ]

dt_typing = merge(dt_typing,meta_filt_dt, by.x="EPI_ID",by.y="Isolate_Id")
dt_typing[["week"]] = strftime(as.Date(dt_typing$Collection_Date, format="%Y-%m-%d"), 
                               format = "%V")

dt_typing[["week"]] = factor(dt_typing[["week"]],
                             levels=c("48","49","50","51","52",
                                      "01","02","03","04",
                                      "05","06","07","08",
                                      "09","10","11","12","13"))

dt_typing[["year"]] = strftime(as.Date(dt_typing$Collection_Date, format="%Y-%m-%d"), 
                               format = "%Y")

#Lollipop results are based on alignment against H1N1 ref, hence need to select only 
#those in gisaid comparison
H1N1_subtypes = c('6B.1A.5a.2a','6B.1A.5a.2a.1',"6B.1", "6B.1A.5a.2")


dt_gisaid = dt_typing %>% 
  group_by(year,week,Clade) %>% 
  filter(Clade %in% H1N1_subtypes) %>% 
  summarise(Clade_count = n()) %>% 
  ungroup() %>% 
  group_by(year, week) %>%  
  mutate(total_count = sum(Clade_count),  
         relative_abundance_gisaid =Clade_count / total_count) %>% 
  ungroup()

#Acession IDs in in Source Data 2
#fwrite(dt_gisaid, file = "Clinical_abundance.tsv")

#### Calculating clinical abundance estimates #####

#Output of deconvoultion, see 2. in the README
deconv = fread('../../work-IA_H1/v-pipe/results_timeseries/Lollipop/deconvoluted.tsv',
               header = TRUE)
deconv[, clade_ww := ifelse(variant == "5a2a1", "6B.1A.5a.2a.1","6B.1A.5a.2a")]
deconv = deconv[variant != 'undetermined', ]
setnames(deconv, old = c("proportion"), new = c("relative_abundance_ww_deconv_lin"))
setnames(deconv, old = c("date"), new = c("Date"))
deconv[["Date"]] = as.Date(deconv$Date, format="%Y-%m-%d")
deconv[["week"]] = strftime(as.Date(deconv$Date, format="%Y-%m-%d"), 
                            format = "%V")

dt_all = merge(dt_gisaid,
                deconv[,c("week", "relative_abundance_ww_deconv_lin","clade_ww")],
                by.x=c('week','Clade'),
                by.y=c('week','clade_ww'),
                all.y = TRUE)

fwrite(dt_all, file = "dt_Fig3B.tsv")


