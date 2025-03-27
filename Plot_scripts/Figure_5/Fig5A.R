library(tidyverse)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggforce)
library(ggh4x)

setwd("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_5/")
source("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_5/Data_cleaning_H1.R")

dt_plot = dt_clean_H1

p <- ggplot(dt_plot , 
            aes(x = as.factor(gene_codonIDX), y = as.factor(Date), fill = AF)) +
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
