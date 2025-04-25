library(data.table)
library(tidyverse)
library(ggplot2)

#setwd("path_to_Figure_S6")

dt = fread("../Figure_3/dt_Fig3B.tsv")

##### Small amount of cleaning ####

dt[["year"]] = factor(dt[["year"]],
                      levels=c("2023","2024"))

dt[["week"]] = factor(as.character(dt[["week"]]),
                      levels= c("49","50","51","52",
                                "1","2","3","4",
                                "5","6","7","8",
                                "9","10","11","12"))

### Plotting #####

p1 <- ggplot(dt)+  
  geom_point(aes(week, relative_abundance_ww_deconv_lin, colour = "ww-estimates"), alpha = 0.7,size=3)+
  geom_point(aes(week, relative_abundance_gisaid, colour = "gisaid-estimates"), alpha = 0.7,size=3)+
  geom_smooth(aes(x = as.numeric(week), y = relative_abundance_ww_deconv_lin, colour = "ww-estimates"),
              method = "loess")+
  geom_smooth(aes(x = as.numeric(week), y = relative_abundance_gisaid, colour = "gisaid-estimates"),
              method = "loess") +
  scale_color_discrete(labels=c('Clinical Data Europe', 'Wastewater Data Geneva'))+
  facet_grid(rows=vars(Clade))+
  ylab("Abundance estimates")+
  xlab("Week (2023/2024)")+
  theme(text = element_text(size = 20),
        legend.title =element_blank(),
        legend.position="bottom",
        legend.text =  element_text(margin = margin(r = 15, l=5)),
        axis.title.x = element_text(margin = margin(t = 20)),  
        axis.title.y = element_text(margin = margin(r = 20)))

p1
