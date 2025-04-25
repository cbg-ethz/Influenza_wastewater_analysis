library(tidyverse)


# setwd("path_to_Figure_2")
source("helper_functions.R")


out_dir = "/Users/anjohn/Desktop/manuscripts/Influenza/Figures/"


#### GENERATING PLOTS #####


#XX_coverages_combined.tsv.gz can be generated using the bash script combining_coverages.sh
# example command:
# ./combining_coverages.sh IAV_wastewater/work-IA_H1/v-pipe/results_timeseries/

H1_cov_ex1 = "../../work-IA_H1/v-pipe/results_experiment_1/H1_coverages_combined.tsv.gz"
H1_cov_ex2 = "../../work-IA_H1/v-pipe/results_experiment_2/H1_coverages_combined.tsv.gz"

N1_cov_ex1 = "../../work-IA_N1/v-pipe/results_experiment_1/N1_coverages_combined.tsv.gz"
N1_cov_ex2 = "../../work-IA_N1/v-pipe/results_experiment_2/N1_coverages_combined.tsv.gz"

M_cov_ex1 = "../../work-IA_M/v-pipe/results_experiment_1/M_coverages_combined.tsv.gz"
M_cov_ex2 = "../../work-IA_M/v-pipe/results_experiment_2/M_coverages_combined.tsv.gz"



ex1_trans_file = "../sample_mapping/samples_Info_experiment_1.csv"
ex2_trans_file = "../sample_mapping/samples_Info_experiment_2.csv"


#combining results from both experiments
dt_H1 = bind_rows(formatting_dt(H1_cov_ex1,ex1_trans_file,"_ex1"),
                  formatting_dt(H1_cov_ex2,ex2_trans_file,"_ex2"))

dt_N1 = bind_rows(formatting_dt(N1_cov_ex1,ex1_trans_file,"_ex1"),
                  formatting_dt(N1_cov_ex2,ex2_trans_file,"_ex2"))

dt_M = bind_rows(formatting_dt(M_cov_ex1,ex1_trans_file,"_ex1"),
                 formatting_dt(M_cov_ex2,ex2_trans_file,"_ex2"))



H1_coverage = coverage_plotting(dt_H1)
H1_coverage$plot

N1_coverage = coverage_plotting(dt_N1)
N1_coverage$plot

M_coverage = coverage_plotting(dt_M)
M_coverage$plot


ggsave(paste0(out_dir,"FIG2_B_H1.png"), 
       plot = H1_coverage$plot, width = 10.5, 
       height = 11, units = "in", dpi = 600)

ggsave(paste0(out_dir,"FIG2_B_N1.png"), 
       plot = N1_coverage$plot, width = 10.5, 
       height = 11, units = "in", dpi = 600)

ggsave(paste0(out_dir,"FIG2_B_M.png"), 
       plot = M_coverage$plot, width = 10.5, 
       height = 11, units = "in", dpi = 600)



