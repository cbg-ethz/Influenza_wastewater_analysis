library(tidyverse)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggforce)
library(ggh4x)

setwd("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_4/")
source("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_4/Data_cleaning.R")

dt_plot = dt_clean_N1

p <- ggplot(dt_plot, aes(x = as.factor(gene_codonIDX), y = as.factor(Date), fill = AF)) +
  geom_tile()+
  facet_grid(vars(year), switch = "y", scales = "free") +
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  force_panelsizes(rows = c(0.8,0.35))+
  scale_x_discrete(breaks=dt_plot$gene_codonIDX, labels=dt_plot$mut_lable)+
  scale_y_discrete(name = "Week",
                   labels = setNames(dt_plot$week, dt_plot$Date))+
  # geom_text(aes(label = round(frac, 2)), color = "white", size = 4) +
  theme_light()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18,margin = margin(r=10)),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16,angle = 90, hjust = 1,vjust = 0.3),
        legend.title = element_text(vjust = .8, size = 16),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        strip.text = element_text(size = 16))+
  labs(fill = "Mutation Frequency")

p

# ggsave(paste0("xxx"), 
#        plot = q, width = 12, height = 8, units = "in", dpi = 300)





#needs to be separate supplementary file

source("~/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_2/helper_functions.R")
timeseries_trans_file = "../sample_mapping/samples_Info_timeseries.csv"

N1_cov_timeseries = "../../work-IA_N1/v-pipe/results_timeseries/N1_coverages_combined.tsv.gz"

dt_N1 = formatting_dt(N1_cov_timeseries,timeseries_trans_file)



dt_pos = dt_clean_N1[,c('POS','mut_lable','sample_name','gene_codonIDX')]
dt_pos = dt_pos %>% distinct(POS, .keep_all = TRUE)
z = dt_N1 %>% filter(Position %in% dt_pos$POS)
dt_plot = merge(z, dt_pos, by.x = "Position", by.y = "POS", all.y = TRUE)


dt_plot[, Date := as.Date(Date, format="%Y_%m_%d")]
dt_plot[["week"]] = strftime(as.Date(dt_plot$Date, format="%Y_%m_%d"), 
                             format = "%V")
dt_plot[["year"]] = strftime(as.Date(dt_plot$Date, format="%Y_%m_%d"), 
                             format = "%Y")

dt_plot[["year"]] = factor(dt_plot[["year"]],levels=c("2024","2023"))








ggplot(dt_plot, aes(x = as.factor(Position), y = as.factor(Date), fill = Coverage)) +
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






