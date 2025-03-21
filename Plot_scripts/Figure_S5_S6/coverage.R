library(tidyverse)
library(ggplot2)
library(gridExtra)
library(Biostrings)
library(data.table)
library(stringr)

#TODO: remove before making public
setwd("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_S5_S6")
source("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_2/helper_functions.R")

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

S5 = coverage_plotting("../../work-IA_N1/v-pipe/results_timeseries/N1_coverages_combined.tsv.gz")
S5

S6 = coverage_plotting("../../work-IA_H1/v-pipe/results_timeseries/H1_coverages_combined.tsv.gz")
S6


### mean coverage as listed in S1 ####

timeseries_cov = "../../work-IA_M/v-pipe/results_timeseries/M_coverages_combined.tsv.gz"
trans_file = "../sample_mapping/samples_Info_timeseries.csv"



dt_cov = formatting_dt(timeseries_cov,trans_file,"_timeseries")

dt = dt_cov %>% 
  group_by(name,Season,Date,Place) %>% 
  summarise(mean_pos_coverag = mean(Coverage)) %>% 
  ungroup() %>% 
  setDT()

dt

