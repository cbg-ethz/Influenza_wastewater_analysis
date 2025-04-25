library(tidyverse)
library(ggplot2)
library(gridExtra)
library(Biostrings)
library(data.table)
library(stringr)

# setwd("path_to_Figure_S4/")


formatting_dt <- function(file,
                          trans_file,
                          extension=""){
  
  dt = fread(file, col.names = c("Ref","Position","Coverage","sample"))
  
  dt[, c("name", "batch") := tstrsplit(sample, "/", 2)]
  dt[,sample:=NULL]
  dt$name = paste0(str_trim(dt$name),extension)
  
  location_translation = fread(trans_file, header = TRUE)
  location_translation$Barcode = paste0("barcode",
                                        str_trim(location_translation$Barcode),
                                        extension) 
  
  dt = merge(dt,location_translation,
             by.x=c('name'),
             by.y=c('Barcode'))
  
  dt$Place = str_trim(dt$Place)
  dt$Date = str_trim(dt$Date)
  dt$Season = str_trim(dt$Season)
  
  
  return(dt)
  
  
}

coverage_plotting = function(dt){

  dt$Place <- str_trim(gsub("^[A-Z]+", "", dt$Place))
  dt$Place_f = factor(dt$Place, levels=c('Chur','Werdhoelzli','Lugano','Aire'))
  
  
  p <- dt %>%
    ggplot(aes(x=Position,y=Coverage),group = Date)+
    geom_line(stat = 'identity',aes(colour = Date), alpha=0.8,linewidth = 0.75)+
    scale_y_continuous('Read Depth',
                       trans='log10',
                       labels = scales::scientific
    )+
    #NA only
    scale_x_continuous("Position")+
    theme(legend.position="none",
          title = element_text(size = 16),
          axis.title.y = element_text(margin = margin(r=10), size = 20),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.title.x = element_text(margin = margin(t=10), size = 20),
          strip.text = element_text(size = 18))+
    facet_grid(vars(Place_f),vars(Season))
  
  return(list(plot = p, plot_data = dt))
  
  
}


#### GENERATING PLOTS #####


#XX_coverages_combined.tsv.gz can be generated using the bash script combining_coverages.sh
# example command:
# ./combining_coverages.sh IAV_wastewater/work-IA_H1/v-pipe/results_timeseries/

H3_cov_ex1 = "../../work-IA_H3/v-pipe/results_experiment_1/H3_coverages_combined.tsv.gz"
H3_cov_ex2 = "../../work-IA_H3/v-pipe/results_experiment_2/H3_coverages_combined.tsv.gz"

N2_cov_ex1 = "../../work-IA_N2/v-pipe/results_experiment_1/N2_coverages_combined.tsv.gz"
N2_cov_ex2 = "../../work-IA_N2/v-pipe/results_experiment_2/N2_coverages_combined.tsv.gz"


ex1_trans_file = "../sample_mapping/samples_Info_experiment_1.csv"
ex2_trans_file = "../sample_mapping/samples_Info_experiment_2.csv"


#combining results from both experiments
dt_H3 = bind_rows(formatting_dt(H3_cov_ex1,ex1_trans_file,"_ex1"),
                  formatting_dt(H3_cov_ex2,ex2_trans_file,"_ex2"))

dt_N2 = bind_rows(formatting_dt(N2_cov_ex1,ex1_trans_file,"_ex1"),
                  formatting_dt(N2_cov_ex2,ex2_trans_file,"_ex2"))



H3_coverage = coverage_plotting(dt_H3)
H3_coverage$plot

N2_coverage = coverage_plotting(dt_N2)
N2_coverage$plot


