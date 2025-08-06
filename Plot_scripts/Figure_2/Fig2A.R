library(tidyverse)
library(scales)
library(ggplot2)


# setwd("path_to_Figure_2")
source("helper_functions.R")

out_dir = "path_of_your_choice"

##ampCoverage.tsv FILES WERE GENERATED using run_ampCov.sh ####

H1_ampCov <- fread("H1_comb_ampCoverage.tsv")
H3_ampCov <- fread("H3_comb_ampCoverage.tsv")

N1_ampCov <- fread("N1_comb_ampCoverage.tsv")
N2_ampCov <- fread("N2_comb_ampCoverage.tsv")

HA_ampCov <- rbind(H1_ampCov,H3_ampCov)
NA_ampCov <- rbind(N1_ampCov,N2_ampCov)
M_ampCov <- fread("M_comb_ampCoverage.tsv")


HA_plot = amplicon_plotting(HA_ampCov)
NA_plot = amplicon_plotting(NA_ampCov)
M_plot = amplicon_plotting(M_ampCov)



ggsave(paste0(out_dir,"FIG2_A_HA.png"),
         plot = HA_plot, width = 10.5, 
       height = 6, units = "in", dpi = 600)

ggsave(paste0(out_dir,"FIG2_A_NA.png"),
       plot = NA_plot, width = 10.5, 
       height = 6, units = "in", dpi = 600)

ggsave(paste0(out_dir,"FIG2_A_M.png"),
       plot = M_plot, width = 10.5, 
       height = 6, units = "in", dpi = 600)




##### Mean HA read count incld. time series data #####

H1_ampCov_time <- fread("H1_timeseries_ampCoverage.tsv", header = TRUE)
H1_mean_time = mean(colMeans(H1_ampCov_time[,-1]))

H3_ampCov_time <- fread("H3_timeseries_ampCoverage.csv", header = TRUE)
H3_mean_time = mean(colMeans(H3_ampCov_time[,-1]))

HA_mean_readcount = mean(H1_mean_time,H3_mean_time, mean(HA_ampCov$read_count))
