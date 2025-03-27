library(tidyverse)
library(ggplot2)
library(gridExtra)
library(Biostrings)
library(data.table)
library(ggh4x)
library(stringr)

#TODO: remove before making public
setwd("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_S5_S6")
source("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_2/helper_functions.R")
source("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_4/Data_cleaning_N1.R")
source("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_5/Data_cleaning_H1.R")


timeseries_trans_file = "../sample_mapping/samples_Info_timeseries.csv"


#' Function to produce Supplementary plots S5 and S6 plots
#'
#' @param file combined coverage (over all time series samples), made with combining_coverages.sh
#' @param trans_file csv file to translate sample information (here samples_Info_timeseries.csv)
#' @param legend boolean value, if to include a legend or not
#'
#' @return plot like seeen in either S5 or S6
#' @export
#'
#' @examples
#' S5 = coverage_plotting("../../work-IA_N1/v-pipe/results_timeseries/N1_coverages_combined.tsv.gz")


coverage_plotting = function(file, 
                             trans_file = "../sample_mapping/samples_Info_timeseries.csv",
                             legend = FALSE){

  location_translation = fread(trans_file, header = TRUE)
  location_translation$Barcode = paste0("barcode",str_trim(location_translation$Barcode)) 
  
  dt = fread(file)
  dt[, c("name", "batch") := tstrsplit(V4, "/", 2)]
  dt[,V4:=NULL]
  
  dt = merge(dt,location_translation,
             by.x=c('name'),
             by.y=c('Barcode'))
  
  legend= TRUE
  
  if(legend == FALSE){
    
    p = dt %>% 
      ggplot(aes(x=V2,y=V3))+
      geom_line(stat = 'identity',aes(colour = Date))+
      scale_y_continuous('Read Depth',
                         trans='log10')+
      scale_x_continuous("Position")+ 
      theme(legend.position='none')+
      facet_wrap(vars(Place))
    
    }else{
      
      p = dt %>% 
        ggplot(aes(x=V2,y=V3))+
        geom_line(stat = 'identity',aes(colour = Date))+
        scale_y_continuous('Read Depth',
                           trans='log10')+
        scale_x_continuous("Position")+
        facet_wrap(vars(Place))
  }

  
  return(p)
  
  
}


#### GENERATING PLOTS #####


#XX_coverages_combined.tsv.gz can be generated using the bash script combining_coverages.sh
# example command:
# ./combining_coverages.sh IAV_wastewater/work-IA_H1/v-pipe/results_timeseries/



#Supplementary Figure 5A (N1 genome-wide coverage)
S5A = coverage_plotting("../../work-IA_N1/v-pipe/results_timeseries/N1_coverages_combined.tsv.gz")

#Supplementary Figure 5A (N1 coverage of target muations)
N1_cov_timeseries = "../../work-IA_N1/v-pipe/results_timeseries/N1_coverages_combined.tsv.gz"
dt_N1 = formatting_dt(N1_cov_timeseries,timeseries_trans_file)

dt_pos = dt_clean_N1[,c('POS','mut_lable','sample_name','gene_codonIDX')]
dt_pos = dt_pos %>% distinct(POS, .keep_all = TRUE)
dt_plot = merge(dt_N1 %>% filter(Position %in% dt_pos$POS), 
                dt_pos, by.x = "Position", by.y = "POS", all.y = TRUE)


dt_plot[, Date := as.Date(Date, format="%Y_%m_%d")]
dt_plot[["week"]] = strftime(as.Date(dt_plot$Date, format="%Y_%m_%d"), 
                             format = "%V")
dt_plot[["year"]] = strftime(as.Date(dt_plot$Date, format="%Y_%m_%d"), 
                             format = "%Y")

dt_plot[["year"]] = factor(dt_plot[["year"]],levels=c("2024","2023"))



S5B = ggplot(dt_plot, aes(x = as.factor(Position), 
                          y = as.factor(Date), 
                          fill = Coverage)) +
        geom_tile()+
        facet_grid(vars(year), switch = "y", scales = "free") +
        force_panelsizes(rows = c(0.8,0.35))+
        scale_fill_gradient(low = "lightblue", high = "darkblue",
                            trans = "log", breaks = c(1, 10, 100, 1000))+
        scale_x_discrete(breaks=dt_plot$Position, labels=dt_plot$mut_lable)+
        scale_y_discrete(name = "Week",
                         labels = setNames(dt_plot$week, dt_plot$Date))+
        theme_light()+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_text(size = 18,margin = margin(r=10)),
              axis.text.y = element_text(size = 16),
              axis.text.x = element_text(size = 16,angle = 90, hjust = 1,vjust = 0.3),
              legend.title = element_text(vjust = .8, size = 16),
              legend.text = element_text(size = 14),
              legend.position = "bottom",
              legend.key.width = unit(1, "cm"),
              strip.text = element_text(size = 16))

ggsave(paste0("/Users/anjohn/Desktop/manuscripts/Influenza/Figures/Supplements/S5A.png"),
       plot = S5A, width = 8, height = 6, units = "in", dpi = 600)


ggsave(paste0("/Users/anjohn/Desktop/manuscripts/Influenza/Figures/Supplements/S5B.png"),
       plot = S5B, width = 12, height = 8, units = "in", dpi = 600)





S6A = coverage_plotting("../../work-IA_H1/v-pipe/results_timeseries/H1_coverages_combined.tsv.gz")

H1_cov_timeseries = "../../work-IA_H1/v-pipe/results_timeseries/H1_coverages_combined.tsv.gz"
dt_H1 = formatting_dt(H1_cov_timeseries,timeseries_trans_file)

dt_pos = dt_clean_H1[,c('POS','mut_lable','sample_name','gene_codonIDX')]
dt_pos = dt_pos %>% distinct(POS, .keep_all = TRUE)
dt_plot = merge(dt_H1 %>% filter(Position %in% dt_pos$POS), 
                dt_pos, by.x = "Position", by.y = "POS", all.y = TRUE)


dt_plot[, Date := as.Date(Date, format="%Y_%m_%d")]
dt_plot[["week"]] = strftime(as.Date(dt_plot$Date, format="%Y_%m_%d"), 
                             format = "%V")
dt_plot[["year"]] = strftime(as.Date(dt_plot$Date, format="%Y_%m_%d"), 
                             format = "%Y")

dt_plot[["year"]] = factor(dt_plot[["year"]],levels=c("2024","2023"))



S6B = ggplot(dt_plot, aes(x = as.factor(Position), 
                          y = as.factor(Date), 
                          fill = Coverage)) +
  geom_tile()+
  facet_grid(vars(year), switch = "y", scales = "free") +
  force_panelsizes(rows = c(0.8,0.35))+
  scale_fill_gradient(low = "lightblue", high = "darkblue",
                      trans = "log", breaks = c(1, 10, 100, 1000))+
  scale_x_discrete(breaks=dt_plot$Position, labels=dt_plot$mut_lable)+
  scale_y_discrete(name = "Week",
                   labels = setNames(dt_plot$week, dt_plot$Date))+
  theme_light()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18,margin = margin(r=10)),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16,angle = 90, hjust = 1,vjust = 0.3),
        legend.title = element_text(vjust = .8, size = 16),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        strip.text = element_text(size = 16))


ggsave(paste0("/Users/anjohn/Desktop/manuscripts/Influenza/Figures/Supplements/S6A.png"),
       plot = S6A, width = 8, height = 6, units = "in", dpi = 600)


ggsave(paste0("/Users/anjohn/Desktop/manuscripts/Influenza/Figures/Supplements/S6B.png"),
       plot = S6B, width = 12, height = 8, units = "in", dpi = 600)



### mean coverage as listed in Tabel S1 ####

timeseries_cov = "../../work-IA_M/v-pipe/results_timeseries/M_coverages_combined.tsv.gz"
trans_file = "../sample_mapping/samples_Info_timeseries.csv"



dt_cov = formatting_dt(timeseries_cov,trans_file,"_timeseries")

dt = dt_cov %>% 
  group_by(name,Season,Date,Place) %>% 
  summarise(mean_pos_coverag = mean(Coverage)) %>% 
  ungroup() %>% 
  setDT()

dt

