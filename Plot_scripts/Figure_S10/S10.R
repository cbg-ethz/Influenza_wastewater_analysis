library(data.table)
library(tidyverse)
library(boot)
library(ggplot2)
library(RColorBrewer)

setwd("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_SX2")

dt_CH = fread("dt_Fig3B_CH.tsv")
dt_EU = fread("../Figure_3/dt_Fig3B.tsv")

out_dir = "PATH/TO/OUTDIR"

prepare_plot_data <- function(dt, clinical_level) {
 
  dt[["year"]] = factor(dt[["year"]], levels=c("2023","2024"))
  
  dt[["week"]] = factor(as.character(dt[["week"]]),
                        levels= c("49","50","51","52",
                                  "1","2","3","4",
                                  "5","6","7","8",
                                  "9","10","11","12"))
  
  # Format for proportional area chart
  dt_plot_long <- melt(dt, 
                       id.vars = c("week", "Clade", "year"), 
                       measure.vars = c("relative_abundance_ww_deconv_lin", 
                                        "relative_abundance_gisaid"), 
                       variable.name = "data_source", 
                       value.name = "relative_abundance")
  
  dt_plot_long$data_source <- factor(dt_plot_long$data_source, 
                                     levels = c("relative_abundance_ww_deconv_lin", 
                                                "relative_abundance_gisaid"), 
                                     labels = c("Wastewater Data Geneva", clinical_level))
  
  # Calculate unassigned fraction
  dt_unassaigned = dt_plot_long %>% 
    group_by(week, data_source) %>% 
    summarise(relative_abundance = round(1-sum(relative_abundance),3)) %>% 
    setDT()
  
  dt_unassaigned[["Clade"]] = "unassigned"
  
  # Calculate statistics for wastewater
  median_unassigned = median(dt_unassaigned[data_source == 'Wastewater Data Geneva'][['relative_abundance']]) 
  mean_unassigned = mean(dt_unassaigned[data_source == 'Wastewater Data Geneva'][['relative_abundance']])
  
  # Combine and format final dataset
  dt_plot = bind_rows(dt_plot_long, dt_unassaigned) %>% 
    filter(week != 52) #drop-out
  
  dt_plot[["Clade"]] = factor(dt_plot[["Clade"]],
                              levels= c("unassigned", "6B.1A.5a.2a","6B.1A.5a.2a.1"))
  
  # Return both the plot data and the statistics
  return(dt_plot)
}


dt_plot_CH = prepare_plot_data(dt_CH,"Clinical Data Switzerland")

dt_plot_EU = prepare_plot_data(dt_EU,"Clinical Data Europe")

dt_plot = rbind(dt_plot_EU[grepl("^Clin", data_source)],
                dt_plot_CH[grepl("^Clin", data_source)])

colour_palette<-c(
  "unassigned"="grey",
  "6B.1A.5a.2a.1"="#1F78B4",
  "6B.1A.5a.2a"="#A6CEE3")

dt_missing = dt_plot %>%
  group_by(week, data_source) %>%
  mutate(relative_abundance = if(any(is.na(relative_abundance))) NA else relative_abundance) %>%
  ungroup()

dt_zeroAbund = dt_plot %>%
  mutate(relative_abundance = replace_na(relative_abundance, 0))



SX2 <- ggplot(dt_missing, aes(group=Clade))+
  geom_area(aes(x = week, 
                y = relative_abundance,
                fill = Clade),
            stat="identity", position="stack",
            colour="white")+
  ylab("Abundance estimates")+
  xlab("Week (2023/2024)")+
  facet_grid(cols=vars(data_source))+
  scale_fill_manual(values = colour_palette,
                    breaks=c("6B.1A.5a.2a",
                             "6B.1A.5a.2a.1",
                             "unassigned"))+
  theme(text = element_text(size = 20),
        legend.title =element_blank(),
        legend.position="bottom",
        # strip.text = element_text(size = 22),
        legend.text =  element_text(margin = margin(r = 15, l=5)),
        axis.title.x = element_text(margin = margin(t = 20)),  
        axis.title.y = element_text(margin = margin(r = 20)),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linewidth = 0.5),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "lightgrey"), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "lightgrey", 
                                    fill=NA, 
                                    linewidth=1))


SX2_zero <- ggplot(dt_zeroAbund, aes(group=Clade))+
  geom_area(aes(x = week, 
                y = relative_abundance,
                fill = Clade),
            stat="identity", position="stack",
            colour="white")+
  ylab("Abundance estimates")+
  xlab("Week (2023/2024)")+
  facet_grid(cols=vars(data_source))+
  scale_fill_manual(values = colour_palette,
                    breaks=c("6B.1A.5a.2a",
                             "6B.1A.5a.2a.1",
                             "unassigned"))+
  theme(text = element_text(size = 20),
        legend.title =element_blank(),
        legend.position="bottom",
        # strip.text = element_text(size = 22),
        legend.text =  element_text(margin = margin(r = 15, l=5)),
        axis.title.x = element_text(margin = margin(t = 20)),  
        axis.title.y = element_text(margin = margin(r = 20)),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linewidth = 0.5),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "lightgrey"), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "lightgrey", 
                                    fill=NA, 
                                    linewidth=1))



SX2_zero

# ggsave(paste0(out_dir,"FIG_SX2.png"),
#        plot = SX2, width = 10, height = 5, units = "in", dpi = 600)
# 
# ggsave(paste0(out_dir,"FIG_SX2_zeroAbund.png"),
#        plot = SX2_zero, width = 10, height = 5, units = "in", dpi = 600)
